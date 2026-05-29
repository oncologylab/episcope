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

static inline double rng_unif(std::mt19937_64& rng) {
  return (static_cast<double>(rng()) + 0.5) /
    (static_cast<double>(std::numeric_limits<std::uint64_t>::max()) + 1.0);
}

static inline std::uint64_t mix_seed(std::uint64_t seed,
                                     std::uint64_t phase,
                                     std::uint64_t iter,
                                     std::uint64_t thread_id) {
  std::uint64_t x = seed + 0x9e3779b97f4a7c15ULL;
  x ^= phase + 0xbf58476d1ce4e5b9ULL + (x << 6) + (x >> 2);
  x ^= iter + 0x94d049bb133111ebULL + (x << 6) + (x >> 2);
  x ^= thread_id + 0xda942042e4dd58b5ULL + (x << 6) + (x >> 2);
  return x;
}

static inline int sample_topic(const std::vector<double>& prob,
                               double total,
                               std::mt19937_64& rng) {
  if (!(total > 0.0) || !R_finite(total)) {
    return static_cast<int>(rng() % static_cast<std::uint64_t>(prob.size()));
  }
  const double target = rng_unif(rng) * total;
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
                               int iter) {
  const std::size_t word_k = static_cast<std::size_t>(dat.n_word) * static_cast<std::size_t>(K);
  std::vector< std::vector<int> > local_word_delta(static_cast<std::size_t>(n_threads));
  std::vector< std::vector<int> > touched(static_cast<std::size_t>(n_threads));
  std::vector<int> local_topic_delta(static_cast<std::size_t>(n_threads) * static_cast<std::size_t>(K), 0);
  for (int t = 0; t < n_threads; ++t) {
    local_word_delta[t].assign(word_k, 0);
    touched[t].reserve(1024);
  }
  const double beta_sum = static_cast<double>(dat.n_word) * beta;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    std::mt19937_64 rng(mix_seed(seed, 1U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tid)));
    std::vector<double> prob(static_cast<std::size_t>(K));
    std::vector<int>& delta = local_word_delta[tid];
    std::vector<int>& touched_t = touched[tid];
    int* topic_delta = &local_topic_delta[static_cast<std::size_t>(tid) * static_cast<std::size_t>(K)];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int d = 0; d < dat.n_doc; ++d) {
      for (int ptr = dat.doc_ptr[d]; ptr < dat.doc_ptr[d + 1]; ++ptr) {
        const int tok = dat.doc_token_index[ptr];
        const int w = dat.token_word[tok];
        const int old_topic = topic[tok];

        doc_topic[d * K + old_topic] -= 1;
        const int old_idx = w * K + old_topic;
        if (delta[old_idx] == 0) touched_t.push_back(old_idx);
        delta[old_idx] -= 1;
        topic_delta[old_topic] -= 1;

        double total = 0.0;
        for (int k = 0; k < K; ++k) {
          const int wk = w * K + k;
          const double wt = static_cast<double>(word_topic[wk] + delta[wk]) + beta;
          const double dt = static_cast<double>(doc_topic[d * K + k]) + alpha;
          double denom = static_cast<double>(topic_count[k] + topic_delta[k]) + beta_sum;
          if (!(denom > 0.0)) denom = beta_sum;
          const double val = dt * wt / denom;
          prob[k] = val;
          total += val;
        }

        const int new_topic = sample_topic(prob, total, rng);
        topic[tok] = new_topic;
        doc_topic[d * K + new_topic] += 1;
        const int new_idx = w * K + new_topic;
        if (delta[new_idx] == 0) touched_t.push_back(new_idx);
        delta[new_idx] += 1;
        topic_delta[new_topic] += 1;
      }
    }
  }

  for (int tid = 0; tid < n_threads; ++tid) {
    std::vector<int>& touched_t = touched[tid];
    if (!touched_t.empty()) {
      std::sort(touched_t.begin(), touched_t.end());
      touched_t.erase(std::unique(touched_t.begin(), touched_t.end()), touched_t.end());
      std::vector<int>& delta = local_word_delta[tid];
      for (std::size_t i = 0; i < touched_t.size(); ++i) {
        const int idx = touched_t[i];
        word_topic[idx] += delta[idx];
      }
    }
    for (int k = 0; k < K; ++k) {
      topic_count[k] += local_topic_delta[static_cast<std::size_t>(tid) * static_cast<std::size_t>(K) + k];
    }
  }
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
                                int iter) {
  const double beta_sum = static_cast<double>(dat.n_word) * beta;

#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    std::mt19937_64 rng(mix_seed(seed, 2U, static_cast<std::uint64_t>(iter), static_cast<std::uint64_t>(tid)));
    std::vector<double> prob(static_cast<std::size_t>(K));

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int d = 0; d < dat.n_doc; ++d) {
      for (int ptr = dat.doc_ptr[d]; ptr < dat.doc_ptr[d + 1]; ++ptr) {
        const int tok = dat.doc_token_index[ptr];
        const int w = dat.token_word[tok];
        const int old_topic = topic[tok];
        doc_topic[d * K + old_topic] -= 1;

        double total = 0.0;
        for (int k = 0; k < K; ++k) {
          const double wt = static_cast<double>(word_topic[w * K + k]) + beta;
          const double dt = static_cast<double>(doc_topic[d * K + k]) + alpha;
          const double denom = static_cast<double>(topic_count[k]) + beta_sum;
          const double val = dt * wt / denom;
          prob[k] = val;
          total += val;
        }

        const int new_topic = sample_topic(prob, total, rng);
        topic[tok] = new_topic;
        doc_topic[d * K + new_topic] += 1;
      }
    }
  }
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
                              int n_threads = 1) {
  if (K <= 1) stop("K must be greater than 1.");
  if (iterations < 1) stop("iterations must be positive.");
  if (!R_finite(alpha) || alpha <= 0.0) stop("alpha must be positive.");
  if (!R_finite(beta) || beta <= 0.0) stop("beta must be positive.");
  if (n_check_convergence < 0) stop("n_check_convergence must be non-negative.");
  if (n_iter_inference < 0) stop("n_iter_inference must be non-negative.");
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

  std::mt19937_64 init_rng(static_cast<std::uint64_t>(seed));
  for (int tok = 0; tok < n_token; ++tok) {
    const int k = static_cast<int>(init_rng() % static_cast<std::uint64_t>(K));
    const int d = dat.token_doc[tok];
    const int w = dat.token_word[tok];
    topic[tok] = k;
    doc_topic[d * K + k] += 1;
    word_topic[w * K + k] += 1;
    topic_count[k] += 1;
  }

  std::vector<int> trace_iter;
  std::vector<double> trace_loglik;
  double previous_loglik = NA_REAL;
  int actual_iterations = iterations;
  for (int iter = 1; iter <= iterations; ++iter) {
    sample_training_iteration(dat, topic, doc_topic, word_topic, topic_count,
                              alpha, beta, K, n_threads,
                              static_cast<std::uint64_t>(seed), iter);
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

  const int inference_burnin = iterations;
  for (int iter = 1; iter <= inference_burnin; ++iter) {
    sample_inference_iteration(dat, topic, doc_topic, word_topic, topic_count,
                               alpha, beta, K, n_threads,
                               static_cast<std::uint64_t>(seed), iter);
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
      sample_inference_iteration(dat, topic, doc_topic, word_topic, topic_count,
                                 alpha, beta, K, n_threads,
                                 static_cast<std::uint64_t>(seed + 104729), iter);
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
  for (int k = 0; k < K; ++k) {
    const double denom = static_cast<double>(topic_count[k]) + beta_sum;
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
    _["loglik_trace"] = DataFrame::create(
      _["iter"] = wrap(trace_iter),
      _["loglikelihood"] = wrap(trace_loglik),
      _["stringsAsFactors"] = false
    )
  );
}
