// File: warplda_kernels.cpp
// Package: craftgrn
// Author: Yaoxiang Li
//
// Native OpenMP topic-model kernels for CraftGRN Module 3.
//
// The public R contract intentionally mirrors the previous text2vec WarpLDA
// output contract used by CraftGRN: theta is document-topic probability,
// phi is topic-term probability, and metrics include token count,
// perplexity, and approximate log-likelihood.
//
// The implementation is inspired by the WarpLDA/text2vec Rcpp interface:
//   Copyright (C) 2015 - 2017 Dmitriy Selivanov and contributors
//   text2vec is distributed under GPL (>= 2)
//
// CraftGRN modifications:
//   - package-local Rcpp entrypoint
//   - OpenMP document-parallel sampling
//   - text2vec-compatible transform semantics for the default sampler
//   - no runtime dependency on text2vec
//
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

using namespace Rcpp;

namespace {

struct SparseDtmTokens {
  int n_doc;
  int n_word;
  std::vector<int> token_doc;
  std::vector<int> token_word;
  std::vector<int> doc_ptr;
  std::vector<int> doc_token_index;
  std::vector<int> word_ptr;
  std::vector<int> nz_doc;
  std::vector<int> nz_word;
  std::vector<int> nz_count;
};


static inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

static inline double token_unif(std::uint64_t seed,
                                std::uint64_t phase,
                                std::uint64_t iter,
                                std::uint64_t token) {
  const std::uint64_t x = splitmix64(seed ^ (phase * 0xbf58476d1ce4e5b9ULL) ^
                                     (iter * 0x94d049bb133111ebULL) ^ token);
  return static_cast<double>(x >> 11) * (1.0 / 9007199254740992.0);
}

static inline int sample_topic_from_unit(const std::vector<double>& prob,
                                         double total,
                                         double u) {
  if (!(total > 0.0) || !R_finite(total)) {
    const std::uint64_t idx = static_cast<std::uint64_t>(std::floor(u * static_cast<double>(prob.size())));
    return static_cast<int>(std::min<std::uint64_t>(idx, static_cast<std::uint64_t>(prob.size() - 1)));
  }
  const double target = u * total;
  double acc = 0.0;
  for (std::size_t k = 0; k < prob.size(); ++k) {
    acc += prob[k];
    if (target <= acc) return static_cast<int>(k);
  }
  return static_cast<int>(prob.size() - 1);
}

SparseDtmTokens read_dgC_tokens(const S4& dtm) {
  IntegerVector dims = dtm.slot("Dim");
  IntegerVector i_slot = dtm.slot("i");
  IntegerVector p_slot = dtm.slot("p");
  NumericVector x_slot = dtm.slot("x");

  if (dims.size() != 2) stop("DTM must have two dimensions.");
  const int n_doc = dims[0];
  const int n_word = dims[1];
  if (n_doc <= 0 || n_word <= 0) stop("DTM must have positive dimensions.");
  if (p_slot.size() != n_word + 1) stop("DTM column pointer has invalid length.");
  if (i_slot.size() != x_slot.size()) stop("DTM index and value slots have inconsistent lengths.");

  SparseDtmTokens out;
  out.n_doc = n_doc;
  out.n_word = n_word;
  out.doc_ptr.assign(static_cast<std::size_t>(n_doc) + 1U, 0);
  out.word_ptr.assign(static_cast<std::size_t>(n_word) + 1U, 0);

  double token_total_d = 0.0;
  for (int w = 0; w < n_word; ++w) {
    const int bgn = p_slot[w];
    const int end = p_slot[w + 1];
    if (bgn < 0 || end < bgn || end > i_slot.size()) stop("DTM column pointer is invalid.");
    for (int idx = bgn; idx < end; ++idx) {
      const int d = i_slot[idx];
      const double x = x_slot[idx];
      if (d < 0 || d >= n_doc) stop("DTM row index is out of range.");
      if (!R_finite(x) || x <= 0.0) continue;
      const double rounded = std::floor(x + 0.5);
      if (std::fabs(x - rounded) > 1e-7) {
        stop("WarpLDA DTM counts must be positive integer-like values.");
      }
      if (rounded > static_cast<double>(std::numeric_limits<int>::max())) {
        stop("WarpLDA DTM contains a count larger than supported integer range.");
      }
      const int cnt = static_cast<int>(rounded);
      token_total_d += static_cast<double>(cnt);
      out.doc_ptr[static_cast<std::size_t>(d) + 1U] += cnt;
      out.nz_doc.push_back(d);
      out.nz_word.push_back(w);
      out.nz_count.push_back(cnt);
    }
  }

  if (!(token_total_d > 0.0)) stop("DTM has no positive token counts.");
  if (token_total_d > static_cast<double>(std::numeric_limits<int>::max())) {
    stop("Expanded token count exceeds supported integer range.");
  }

  for (int d = 0; d < n_doc; ++d) {
    out.doc_ptr[d + 1] += out.doc_ptr[d];
  }

  const int n_token = static_cast<int>(token_total_d);
  out.token_doc.reserve(n_token);
  out.token_word.reserve(n_token);
  out.word_ptr[0] = 0;

  for (int w = 0; w < n_word; ++w) {
    out.word_ptr[w] = static_cast<int>(out.token_doc.size());
    const int bgn = p_slot[w];
    const int end = p_slot[w + 1];
    for (int idx = bgn; idx < end; ++idx) {
      const int d = i_slot[idx];
      const double x = x_slot[idx];
      if (!R_finite(x) || x <= 0.0) continue;
      const int cnt = static_cast<int>(std::floor(x + 0.5));
      for (int c = 0; c < cnt; ++c) {
        out.token_doc.push_back(d);
        out.token_word.push_back(w);
      }
    }
  }
  out.word_ptr[n_word] = static_cast<int>(out.token_doc.size());

  out.doc_token_index.assign(out.token_doc.size(), 0);
  std::vector<int> cursor = out.doc_ptr;
  for (int t = 0; t < n_token; ++t) {
    const int d = out.token_doc[t];
    out.doc_token_index[cursor[d]++] = t;
  }

  return out;
}

double compute_loglik(const SparseDtmTokens& dat,
                      const std::vector<int>& doc_topic,
                      const std::vector<int>& word_topic,
                      const std::vector<int>& topic_count,
                      double alpha,
                      double beta,
                      int K) {
  const double beta_sum = static_cast<double>(dat.n_word) * beta;
  const double alpha_sum = static_cast<double>(K) * alpha;
  double ll = 0.0;
  for (std::size_t idx = 0; idx < dat.nz_count.size(); ++idx) {
    const int d = dat.nz_doc[idx];
    const int w = dat.nz_word[idx];
    const int cnt = dat.nz_count[idx];
    const int doc_len = dat.doc_ptr[d + 1] - dat.doc_ptr[d];
    double p = 0.0;
    for (int k = 0; k < K; ++k) {
      const double theta = (static_cast<double>(doc_topic[d * K + k]) + alpha) /
        (static_cast<double>(doc_len) + alpha_sum);
      const double phi = (static_cast<double>(word_topic[w * K + k]) + beta) /
        (static_cast<double>(topic_count[k]) + beta_sum);
      p += theta * phi;
    }
    if (p <= 0.0 || !R_finite(p)) p = std::numeric_limits<double>::min();
    ll += static_cast<double>(cnt) * std::log(p);
  }
  return ll;
}

void rebuild_training_counts(const SparseDtmTokens& dat,
                             const std::vector<int>& topic,
                             std::vector<int>& doc_topic,
                             std::vector<int>& word_topic,
                             std::vector<int>& topic_count,
                             int K,
                             int n_threads) {
  std::fill(doc_topic.begin(), doc_topic.end(), 0);
  std::fill(word_topic.begin(), word_topic.end(), 0);
  std::fill(topic_count.begin(), topic_count.end(), 0);

  if (n_threads <= 1) {
    const int n_token = static_cast<int>(topic.size());
    for (int tok = 0; tok < n_token; ++tok) {
      const int k = topic[tok];
      const int d = dat.token_doc[tok];
      const int w = dat.token_word[tok];
      doc_topic[d * K + k] += 1;
      word_topic[w * K + k] += 1;
      topic_count[k] += 1;
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int d = 0; d < dat.n_doc; ++d) {
    int* doc_row = &doc_topic[static_cast<std::size_t>(d) * static_cast<std::size_t>(K)];
    for (int ptr = dat.doc_ptr[d]; ptr < dat.doc_ptr[d + 1]; ++ptr) {
      const int tok = dat.doc_token_index[ptr];
      doc_row[topic[tok]] += 1;
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int w = 0; w < dat.n_word; ++w) {
    int* word_row = &word_topic[static_cast<std::size_t>(w) * static_cast<std::size_t>(K)];
    for (int tok = dat.word_ptr[w]; tok < dat.word_ptr[w + 1]; ++tok) {
      word_row[topic[tok]] += 1;
    }
  }

  for (int k = 0; k < K; ++k) {
    int total = 0;
    for (int w = 0; w < dat.n_word; ++w) {
      total += word_topic[w * K + k];
    }
    topic_count[k] = total;
  }
}

void rebuild_doc_counts(const SparseDtmTokens& dat,
                        const std::vector<int>& topic,
                        std::vector<int>& doc_topic,
                        int K,
                        int n_threads) {
  std::fill(doc_topic.begin(), doc_topic.end(), 0);

  if (n_threads <= 1) {
    const int n_token = static_cast<int>(topic.size());
    for (int tok = 0; tok < n_token; ++tok) {
      const int k = topic[tok];
      const int d = dat.token_doc[tok];
      doc_topic[d * K + k] += 1;
    }
    return;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
  for (int d = 0; d < dat.n_doc; ++d) {
    int* doc_row = &doc_topic[static_cast<std::size_t>(d) * static_cast<std::size_t>(K)];
    for (int ptr = dat.doc_ptr[d]; ptr < dat.doc_ptr[d + 1]; ++ptr) {
      const int tok = dat.doc_token_index[ptr];
      doc_row[topic[tok]] += 1;
    }
  }
}

void sample_training_iteration(const SparseDtmTokens& dat,
                               std::vector<int>& topic,
                               std::vector<int>& doc_topic,
                               std::vector<int>& word_topic,
                               std::vector<int>& topic_count,
                               double alpha,
                               double beta,
                               int K,
                               int n_threads,
                               std::uint64_t seed,
                               int iter,
                               std::vector<int>& next_topic) {
  // Sample from the previous iteration counts, then rebuild counts by
  // independent doc/word partitions. This keeps output invariant to
  // OpenMP thread count and scheduling.
  const int n_token = static_cast<int>(topic.size());
  next_topic.resize(static_cast<std::size_t>(n_token));
  const double beta_sum = static_cast<double>(dat.n_word) * beta;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    std::vector<double> prob(static_cast<std::size_t>(K));
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int tok = 0; tok < n_token; ++tok) {
      const int d = dat.token_doc[tok];
      const int w = dat.token_word[tok];
      const int old_topic = topic[tok];
      double total = 0.0;
      for (int k = 0; k < K; ++k) {
        const int old = (old_topic == k) ? 1 : 0;
        const double dt = static_cast<double>(doc_topic[d * K + k] - old) + alpha;
        const double wt = static_cast<double>(word_topic[w * K + k] - old) + beta;
        double denom = static_cast<double>(topic_count[k] - old) + beta_sum;
        if (!(denom > 0.0)) denom = beta_sum;
        const double val = dt * wt / denom;
        prob[k] = val;
        total += val;
      }
      next_topic[tok] = sample_topic_from_unit(
        prob, total, token_unif(seed, 1U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok))
      );
    }
  }

  topic.swap(next_topic);
  rebuild_training_counts(dat, topic, doc_topic, word_topic, topic_count, K, n_threads);
}

void sample_inference_iteration(const SparseDtmTokens& dat,
                                std::vector<int>& topic,
                                std::vector<int>& doc_topic,
                                const std::vector<int>& word_topic,
                                const std::vector<int>& topic_count,
                                double alpha,
                                double beta,
                                int K,
                                int n_threads,
                                std::uint64_t seed,
                                int iter,
                                std::vector<int>& next_topic) {
  // Keep inference thread-invariant for the same reason as training.
  const int n_token = static_cast<int>(topic.size());
  next_topic.resize(static_cast<std::size_t>(n_token));
  const double beta_sum = static_cast<double>(dat.n_word) * beta;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    std::vector<double> prob(static_cast<std::size_t>(K));
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int tok = 0; tok < n_token; ++tok) {
      const int d = dat.token_doc[tok];
      const int w = dat.token_word[tok];
      const int old_topic = topic[tok];
      double total = 0.0;
      for (int k = 0; k < K; ++k) {
        const int old = (old_topic == k) ? 1 : 0;
        const double wt = static_cast<double>(word_topic[w * K + k]) + beta;
        const double dt = static_cast<double>(doc_topic[d * K + k] - old) + alpha;
        const double denom = static_cast<double>(topic_count[k]) + beta_sum;
        const double val = dt * wt / denom;
        prob[k] = val;
        total += val;
      }
      next_topic[tok] = sample_topic_from_unit(
        prob, total, token_unif(seed, 2U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok))
      );
    }
  }

  topic.swap(next_topic);
  rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
}

void sample_warp_doc_pass(const SparseDtmTokens& dat,
                          std::vector<int>& topic,
                          std::vector<int>& proposal_topic,
                          const std::vector<int>& doc_topic,
                          const std::vector<int>& topic_count,
                          double alpha,
                          double beta,
                          int K,
                          int n_threads,
                          std::uint64_t seed,
                          int iter,
                          std::uint64_t phase_offset,
                          std::vector<int>& next_topic) {
  const double beta_bar = static_cast<double>(K) * beta;
  next_topic.resize(topic.size());

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int d = 0; d < dat.n_doc; ++d) {
      const int bgn = dat.doc_ptr[d];
      const int end = dat.doc_ptr[d + 1];
      const int len = end - bgn;
      if (len <= 0) continue;

      for (int ptr = bgn; ptr < end; ++ptr) {
        const int tok = dat.doc_token_index[ptr];
        const int old_topic = topic[tok];
        const int new_topic = proposal_topic[tok];
        int accepted_topic = old_topic;
        if (new_topic != old_topic) {
          const double num_doc = static_cast<double>(doc_topic[d * K + new_topic]) + alpha;
          const double den_doc = static_cast<double>(doc_topic[d * K + old_topic]) + alpha;
          const double num_all = static_cast<double>(topic_count[old_topic]) + beta_bar;
          const double den_all = static_cast<double>(topic_count[new_topic]) + beta_bar;
          const double accept_prob = (num_doc / den_doc) * (num_all / den_all);
          const double u_accept = token_unif(seed, phase_offset + 1U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok));
          if (u_accept < accept_prob) accepted_topic = new_topic;
        }
        next_topic[tok] = accepted_topic;
      }

      const double copy_prob = static_cast<double>(len) /
        (static_cast<double>(len) + static_cast<double>(K) * alpha);
      for (int ptr = bgn; ptr < end; ++ptr) {
        const int tok = dat.doc_token_index[ptr];
        const double u_kind = token_unif(seed, phase_offset + 2U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok));
        if (u_kind < copy_prob) {
          const double u_pos = token_unif(seed, phase_offset + 3U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok));
          int offset = static_cast<int>(std::floor(u_pos * static_cast<double>(len)));
          if (offset >= len) offset = len - 1;
          const int source_tok = dat.doc_token_index[bgn + offset];
          proposal_topic[tok] = next_topic[source_tok];
        } else {
          const double u_topic = token_unif(seed, phase_offset + 4U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok));
          int sampled_topic = static_cast<int>(std::floor(u_topic * static_cast<double>(K)));
          if (sampled_topic >= K) sampled_topic = K - 1;
          proposal_topic[tok] = sampled_topic;
        }
      }
    }
  }

  topic.swap(next_topic);
}

void sample_warp_word_pass(const SparseDtmTokens& dat,
                           std::vector<int>& topic,
                           std::vector<int>& proposal_topic,
                           const std::vector<int>& word_topic,
                           const std::vector<int>& topic_count,
                           double beta,
                           int K,
                           int n_threads,
                           std::uint64_t seed,
                           int iter,
                           std::uint64_t phase_offset,
                           bool update_word_counts_within_pass,
                           std::vector<int>& next_topic) {
  const double beta_bar = static_cast<double>(K) * beta;
  next_topic.resize(topic.size());

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
    std::vector<int> local_count(static_cast<std::size_t>(K));
    std::vector<double> prob(static_cast<std::size_t>(K));
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int w = 0; w < dat.n_word; ++w) {
      const int bgn = dat.word_ptr[w];
      const int end = dat.word_ptr[w + 1];
      if (bgn >= end) continue;

      for (int k = 0; k < K; ++k) {
        local_count[k] = word_topic[w * K + k];
      }

      for (int tok = bgn; tok < end; ++tok) {
        const int old_topic = topic[tok];
        const int new_topic = proposal_topic[tok];
        int accepted_topic = old_topic;
        if (new_topic != old_topic) {
          const double num_word = static_cast<double>(local_count[new_topic]) + beta;
          const double den_word = static_cast<double>(local_count[old_topic]) + beta;
          const double num_all = static_cast<double>(topic_count[old_topic]) + beta_bar;
          const double den_all = static_cast<double>(topic_count[new_topic]) + beta_bar;
          const double accept_prob = (num_word / den_word) * (num_all / den_all);
          const double u_accept = token_unif(seed, phase_offset + 1U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok));
          if (u_accept < accept_prob) accepted_topic = new_topic;
        }

        if (update_word_counts_within_pass && accepted_topic != old_topic) {
          local_count[accepted_topic] += 1;
          local_count[old_topic] -= 1;
        }
        next_topic[tok] = accepted_topic;
      }

      double total = 0.0;
      for (int k = 0; k < K; ++k) {
        prob[k] = static_cast<double>(local_count[k]) + beta;
        total += prob[k];
      }
      for (int tok = bgn; tok < end; ++tok) {
        proposal_topic[tok] = sample_topic_from_unit(
          prob, total, token_unif(seed, phase_offset + 2U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tok))
        );
      }
    }
  }

  topic.swap(next_topic);
}

void sample_warp_training_iteration(const SparseDtmTokens& dat,
                                    std::vector<int>& topic,
                                    std::vector<int>& proposal_topic,
                                    std::vector<int>& doc_topic,
                                    std::vector<int>& word_topic,
                                    std::vector<int>& topic_count,
                                    double alpha,
                                    double beta,
                                    int K,
                                    int n_threads,
                                    std::uint64_t seed,
                                    int iter,
                                    std::vector<int>& next_topic) {
  sample_warp_doc_pass(dat, topic, proposal_topic, doc_topic, topic_count,
                       alpha, beta, K, n_threads, seed, iter, 10U, next_topic);
  rebuild_training_counts(dat, topic, doc_topic, word_topic, topic_count, K, n_threads);
  sample_warp_word_pass(dat, topic, proposal_topic, word_topic, topic_count,
                        beta, K, n_threads, seed, iter, 20U, true, next_topic);
  rebuild_training_counts(dat, topic, doc_topic, word_topic, topic_count, K, n_threads);
}

void sample_warp_inference_iteration(const SparseDtmTokens& dat,
                                     std::vector<int>& topic,
                                     std::vector<int>& proposal_topic,
                                     std::vector<int>& doc_topic,
                                     const std::vector<int>& word_topic,
                                     const std::vector<int>& topic_count,
                                     double alpha,
                                     double beta,
                                     int K,
                                     int n_threads,
                                     std::uint64_t seed,
                                     int iter,
                                     std::vector<int>& next_topic) {
  sample_warp_doc_pass(dat, topic, proposal_topic, doc_topic, topic_count,
                       alpha, beta, K, n_threads, seed, iter, 30U, next_topic);
  rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
  sample_warp_word_pass(dat, topic, proposal_topic, word_topic, topic_count,
                        beta, K, n_threads, seed, iter, 40U, false, next_topic);
  rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
}

void sample_text2vec_compat_inference_iteration(const SparseDtmTokens& dat,
                                                std::vector<int>& topic,
                                                std::vector<int>& proposal_topic,
                                                std::vector<int>& doc_topic,
                                                const std::vector<int>& word_topic,
                                                const std::vector<int>& topic_count,
                                                double alpha,
                                                double beta,
                                                int K,
                                                int n_threads,
                                                std::uint64_t seed,
                                                int iter,
                                                std::vector<int>& next_topic) {
  // text2vec transform() keeps word-topic counts fixed and reports
  // document-topic counts after the document pass. The word pass still updates
  // latent token topics for the next iteration, but it does not refresh the
  // reported document-topic counts.
  rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
  sample_warp_doc_pass(dat, topic, proposal_topic, doc_topic, topic_count,
                       alpha, beta, K, n_threads, seed, iter, 50U, next_topic);
  rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
  sample_warp_word_pass(dat, topic, proposal_topic, word_topic, topic_count,
                        beta, K, n_threads, seed, iter, 60U, false, next_topic);
}

} // namespace

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export(.craftgrn_warplda_fit_cpp)]]
List craftgrn_warplda_fit_cpp(const S4& dtm,
                              int K,
                              int iterations = 1000,
                              double alpha = 0.1,
                              double beta = 0.1,
                              int seed = 1,
                              double convergence_tol = 1e-3,
                              int n_check_convergence = 10,
                              int n_iter_inference = 10,
                              int n_threads = 1,
                              std::string sampler = "text2vec_compat") {
  if (K <= 1) stop("K must be greater than 1.");
  if (iterations < 1) stop("iterations must be positive.");
  if (!R_finite(alpha) || alpha <= 0.0) stop("alpha must be positive.");
  if (!R_finite(beta) || beta <= 0.0) stop("beta must be positive.");
  if (n_check_convergence < 0) stop("n_check_convergence must be non-negative.");
  if (n_iter_inference < 0) stop("n_iter_inference must be non-negative.");
  if (sampler != "gibbs_sync" && sampler != "warp_mh" && sampler != "text2vec_compat") {
    stop("sampler must be 'text2vec_compat', 'warp_mh', or 'gibbs_sync'.");
  }
  if (n_threads < 1) n_threads = 1;
#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  if (n_threads > max_threads) n_threads = max_threads;
#else
  n_threads = 1;
#endif

  SparseDtmTokens dat = read_dgC_tokens(dtm);
  const double delta_cells_per_thread = static_cast<double>(dat.n_word) * static_cast<double>(K);
  const double max_delta_cells = 100000000.0;
  if (delta_cells_per_thread > 0.0) {
    const int memory_safe_threads = std::max(1, static_cast<int>(std::floor(max_delta_cells / delta_cells_per_thread)));
    if (n_threads > memory_safe_threads) n_threads = memory_safe_threads;
  }
  const int n_token = static_cast<int>(dat.token_word.size());
  std::vector<int> topic(static_cast<std::size_t>(n_token), 0);
  std::vector<int> doc_topic(static_cast<std::size_t>(dat.n_doc) * static_cast<std::size_t>(K), 0);
  std::vector<int> word_topic(static_cast<std::size_t>(dat.n_word) * static_cast<std::size_t>(K), 0);
  std::vector<int> topic_count(static_cast<std::size_t>(K), 0);
  std::vector<int> next_topic(static_cast<std::size_t>(n_token), 0);
  std::vector<int> proposal_topic(static_cast<std::size_t>(n_token), 0);

  std::mt19937_64 init_rng(static_cast<std::uint64_t>(seed));
  for (int tok = 0; tok < n_token; ++tok) {
    const int k = static_cast<int>(init_rng() % static_cast<std::uint64_t>(K));
    const int proposal = static_cast<int>(init_rng() % static_cast<std::uint64_t>(K));
    const int d = dat.token_doc[tok];
    const int w = dat.token_word[tok];
    topic[tok] = k;
    proposal_topic[tok] = proposal;
    doc_topic[d * K + k] += 1;
    word_topic[w * K + k] += 1;
    topic_count[k] += 1;
  }

  std::vector<int> trace_iter;
  std::vector<double> trace_loglik;
  double previous_loglik = NA_REAL;
  int actual_iterations = iterations;
  for (int iter = 1; iter <= iterations; ++iter) {
    if (sampler == "warp_mh" || sampler == "text2vec_compat") {
      sample_warp_training_iteration(dat, topic, proposal_topic, doc_topic, word_topic, topic_count,
                                     alpha, beta, K, n_threads,
                                     static_cast<std::uint64_t>(seed), iter, next_topic);
    } else {
      sample_training_iteration(dat, topic, doc_topic, word_topic, topic_count,
                                alpha, beta, K, n_threads,
                                static_cast<std::uint64_t>(seed), iter, next_topic);
    }
    if (n_check_convergence > 0 && iter % n_check_convergence == 0) {
      const double ll = compute_loglik(dat, doc_topic, word_topic, topic_count, alpha, beta, K);
      trace_iter.push_back(iter);
      trace_loglik.push_back(ll);
      if (R_finite(previous_loglik) && convergence_tol >= 0.0) {
        const double denom = std::max(1.0, std::fabs(previous_loglik));
        const double rel = std::fabs(ll - previous_loglik) / denom;
        if (rel < convergence_tol) {
          actual_iterations = iter;
          break;
        }
      }
      previous_loglik = ll;
    }
    if (iter % 10 == 0) Rcpp::checkUserInterrupt();
  }

  std::vector<int> trained_word_topic;
  std::vector<int> trained_topic_count;
  if (sampler == "warp_mh" || sampler == "text2vec_compat") {
    trained_word_topic = word_topic;
    trained_topic_count = topic_count;

    std::fill(topic.begin(), topic.end(), 0);
    std::fill(proposal_topic.begin(), proposal_topic.end(), 0);
    std::mt19937_64 infer_rng(static_cast<std::uint64_t>(seed));
    for (int tok = 0; tok < n_token; ++tok) {
      topic[tok] = static_cast<int>(infer_rng() % static_cast<std::uint64_t>(K));
      proposal_topic[tok] = static_cast<int>(infer_rng() % static_cast<std::uint64_t>(K));
    }
    if (sampler == "text2vec_compat") {
      rebuild_training_counts(dat, topic, doc_topic, word_topic, topic_count, K, n_threads);
      word_topic = trained_word_topic;
    } else {
      word_topic = trained_word_topic;
      topic_count = trained_topic_count;
      rebuild_doc_counts(dat, topic, doc_topic, K, n_threads);
    }
  }

  const int inference_burnin = iterations;
  for (int iter = 1; iter <= inference_burnin; ++iter) {
    if (sampler == "text2vec_compat") {
      sample_text2vec_compat_inference_iteration(dat, topic, proposal_topic, doc_topic, word_topic, topic_count,
                                                 alpha, beta, K, n_threads,
                                                 static_cast<std::uint64_t>(seed), iter, next_topic);
    } else if (sampler == "warp_mh") {
      sample_warp_inference_iteration(dat, topic, proposal_topic, doc_topic, word_topic, topic_count,
                                      alpha, beta, K, n_threads,
                                      static_cast<std::uint64_t>(seed), iter, next_topic);
    } else {
      sample_inference_iteration(dat, topic, doc_topic, word_topic, topic_count,
                                 alpha, beta, K, n_threads,
                                 static_cast<std::uint64_t>(seed), iter, next_topic);
    }
    if (iter % 10 == 0) Rcpp::checkUserInterrupt();
  }

  std::vector<double> doc_topic_accum(static_cast<std::size_t>(dat.n_doc) * static_cast<std::size_t>(K), 0.0);
  const int n_accum = std::max(1, n_iter_inference);
  if (n_iter_inference == 0) {
    for (std::size_t idx = 0; idx < doc_topic.size(); ++idx) {
      doc_topic_accum[idx] = static_cast<double>(doc_topic[idx]);
    }
  } else {
    for (int iter = 1; iter <= n_iter_inference; ++iter) {
      if (sampler == "text2vec_compat") {
        sample_text2vec_compat_inference_iteration(dat, topic, proposal_topic, doc_topic, word_topic, topic_count,
                                                   alpha, beta, K, n_threads,
                                                   static_cast<std::uint64_t>(seed + 104729), iter, next_topic);
      } else if (sampler == "warp_mh") {
        sample_warp_inference_iteration(dat, topic, proposal_topic, doc_topic, word_topic, topic_count,
                                        alpha, beta, K, n_threads,
                                        static_cast<std::uint64_t>(seed + 104729), iter, next_topic);
      } else {
        sample_inference_iteration(dat, topic, doc_topic, word_topic, topic_count,
                                   alpha, beta, K, n_threads,
                                   static_cast<std::uint64_t>(seed + 104729), iter, next_topic);
      }
      for (std::size_t idx = 0; idx < doc_topic.size(); ++idx) {
        doc_topic_accum[idx] += static_cast<double>(doc_topic[idx]);
      }
    }
  }

  NumericMatrix theta(dat.n_doc, K);
  const double alpha_sum = static_cast<double>(K) * alpha;
  for (int d = 0; d < dat.n_doc; ++d) {
    const int doc_len = dat.doc_ptr[d + 1] - dat.doc_ptr[d];
    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const double v = doc_topic_accum[d * K + k] / static_cast<double>(n_accum) + alpha;
      theta(d, k) = v;
      row_sum += v;
    }
    if (!(row_sum > 0.0)) row_sum = static_cast<double>(doc_len) + alpha_sum;
    for (int k = 0; k < K; ++k) theta(d, k) /= row_sum;
  }

  NumericMatrix phi(K, dat.n_word);
  const double beta_sum = static_cast<double>(dat.n_word) * beta;
  std::vector<int> phi_topic_count(static_cast<std::size_t>(K), 0);
  for (int k = 0; k < K; ++k) {
    int total = 0;
    for (int w = 0; w < dat.n_word; ++w) {
      total += word_topic[w * K + k];
    }
    phi_topic_count[k] = total;
  }
  for (int k = 0; k < K; ++k) {
    const double denom = static_cast<double>(phi_topic_count[k]) + beta_sum;
    for (int w = 0; w < dat.n_word; ++w) {
      phi(k, w) = (static_cast<double>(word_topic[w * K + k]) + beta) / denom;
    }
  }

  double ll = 0.0;
  for (std::size_t idx = 0; idx < dat.nz_count.size(); ++idx) {
    const int d = dat.nz_doc[idx];
    const int w = dat.nz_word[idx];
    const int cnt = dat.nz_count[idx];
    double p = 0.0;
    for (int k = 0; k < K; ++k) p += theta(d, k) * phi(k, w);
    if (p <= 0.0 || !R_finite(p)) p = std::numeric_limits<double>::min();
    ll += static_cast<double>(cnt) * std::log(p);
  }
  const double perplexity = std::exp(-ll / static_cast<double>(n_token));

  return List::create(
    _["theta"] = theta,
    _["phi"] = phi,
    _["n_tokens"] = static_cast<double>(n_token),
    _["perplexity"] = perplexity,
    _["loglik"] = ll,
    _["iterations"] = actual_iterations,
    _["threads"] = n_threads,
    _["sampler"] = sampler,
    _["loglik_trace"] = DataFrame::create(
      _["iter"] = wrap(trace_iter),
      _["loglikelihood"] = wrap(trace_loglik),
      _["stringsAsFactors"] = false
    )
  );
}
