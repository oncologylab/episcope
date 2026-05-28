// File: corr_kernels.cpp
// Package: craftgrn
// Author: Yaoxiang Li
//
// Shared C++ correlation kernels for CraftGRN.
//
// These exported kernels are intentionally module-neutral:
// - sparse_pair_correlations_cpp is used by Module 1 motif-supported tests
//   and Module 2 TF-target / FP-target restricted correlation tests.
// - dense_prediction_stats_cpp is used by Module 1 all-TF prediction
//   chunks after canonical-bound FP filtering.
//
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

static inline double cor_complete(const NumericMatrix& x, const NumericMatrix& y, int i, int j) {
  const int n = x.ncol();
  double mx = 0.0, my = 0.0;
  for (int k = 0; k < n; ++k) {
    mx += x(i, k);
    my += y(j, k);
  }
  mx /= n;
  my /= n;
  double num = 0.0, sx = 0.0, sy = 0.0;
  for (int k = 0; k < n; ++k) {
    const double dx = x(i, k) - mx;
    const double dy = y(j, k) - my;
    num += dx * dy;
    sx += dx * dx;
    sy += dy * dy;
  }
  const double den = std::sqrt(sx * sy);
  if (den == 0.0) return NA_REAL;
  double r = num / den;
  if (r > 1.0) r = 1.0;
  if (r < -1.0) r = -1.0;
  return r;
}

static inline double cor_p_value(double r, int n) {
  if (!R_finite(r) || n <= 2) return NA_REAL;
  if (std::fabs(r) >= 1.0) return 0.0;
  const double denom = std::max(1.0 - r * r, DBL_EPSILON);
  const double stat = r * std::sqrt((n - 2.0) / denom);
  return 2.0 * R::pt(std::fabs(stat), n - 2.0, 0, 0);
}

static std::vector<double> bh_adjust_column(const NumericMatrix& p, int col) {
  const int n = p.nrow();
  std::vector<std::pair<double, int> > vals;
  vals.reserve(n);
  for (int i = 0; i < n; ++i) {
    const double v = p(i, col);
    if (R_finite(v)) vals.push_back(std::make_pair(v, i));
  }
  std::sort(vals.begin(), vals.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
    return a.first < b.first;
  });
  std::vector<double> out(n, NA_REAL);
  const int m = static_cast<int>(vals.size());
  double prev = 1.0;
  for (int rank = m; rank >= 1; --rank) {
    const double raw = vals[rank - 1].first * static_cast<double>(m) / static_cast<double>(rank);
    const double adj = std::min(prev, std::min(1.0, raw));
    out[vals[rank - 1].second] = adj;
    prev = adj;
  }
  return out;
}


static inline double cor_indexed_complete(const NumericMatrix& x, const NumericMatrix& y, int xi, int yi) {
  return cor_complete(x, y, xi, yi);
}

// Sparse pairwise kernel shared by Module 1 and Module 2. It avoids
// materializing a full FP x TF or TF x gene matrix when only restricted
// candidate pairs are needed.
// [[Rcpp::export(.sparse_pair_correlations_cpp)]]
DataFrame sparse_pair_correlations_cpp(NumericMatrix fp,
                                       NumericMatrix tf,
                                       NumericMatrix fp_rank,
                                       NumericMatrix tf_rank,
                                       IntegerVector fp_index,
                                       IntegerVector tf_index,
                                       int n_threads = 1) {
  const int npair = fp_index.size();
  const int ncond = fp.ncol();
  if (tf.ncol() != ncond || fp_rank.nrow() != fp.nrow() || fp_rank.ncol() != ncond ||
      tf_rank.nrow() != tf.nrow() || tf_rank.ncol() != ncond) {
    stop("Input matrices have inconsistent dimensions.");
  }
  NumericVector pearson_r(npair), pearson_p(npair), spearman_r(npair), spearman_p(npair);
#ifdef _OPENMP
  if (n_threads < 1) n_threads = 1;
  omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < npair; ++k) {
    const int i = fp_index[k];
    const int j = tf_index[k];
    if (i < 0 || i >= fp.nrow() || j < 0 || j >= tf.nrow()) {
      pearson_r[k] = NA_REAL;
      pearson_p[k] = NA_REAL;
      spearman_r[k] = NA_REAL;
      spearman_p[k] = NA_REAL;
      continue;
    }
    const double pr = cor_indexed_complete(fp, tf, i, j);
    const double sr = cor_indexed_complete(fp_rank, tf_rank, i, j);
    pearson_r[k] = pr;
    pearson_p[k] = cor_p_value(pr, ncond);
    spearman_r[k] = sr;
    spearman_p[k] = cor_p_value(sr, ncond);
  }
  return DataFrame::create(
    _["pearson_r"] = pearson_r,
    _["pearson_p"] = pearson_p,
    _["spearman_r"] = spearman_r,
    _["spearman_p"] = spearman_p,
    _["stringsAsFactors"] = false
  );
}
// Dense all-pairs kernel used by Module 1 prediction chunks after
// canonical-bound FP filtering. It emits only passing FP-TF predictions
// to keep R-side memory usage bounded.
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export(.dense_prediction_stats_cpp)]]
DataFrame dense_prediction_stats_cpp(NumericMatrix fp,
                                     NumericMatrix tf,
                                     NumericMatrix fp_rank,
                                     NumericMatrix tf_rank,
                                     CharacterVector fp_id,
                                     CharacterVector atac_peak,
                                     CharacterVector tf_name,
                                     double r_cutoff,
                                     double p_cutoff,
                                     double fdr_cutoff,
                                     int n_threads = 1) {
  const int nfp = fp.nrow();
  const int ntf = tf.nrow();
  const int ncond = fp.ncol();
  if (tf.ncol() != ncond || fp_rank.nrow() != nfp || fp_rank.ncol() != ncond ||
      tf_rank.nrow() != ntf || tf_rank.ncol() != ncond) {
    stop("Input matrices have inconsistent dimensions.");
  }
  NumericMatrix pearson(nfp, ntf), spearman(nfp, ntf), pearson_p(nfp, ntf), spearman_p(nfp, ntf);
#ifdef _OPENMP
  if (n_threads < 1) n_threads = 1;
  omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < nfp; ++i) {
    for (int j = 0; j < ntf; ++j) {
      const double pr = cor_complete(fp, tf, i, j);
      const double sr = cor_complete(fp_rank, tf_rank, i, j);
      pearson(i, j) = pr;
      spearman(i, j) = sr;
      pearson_p(i, j) = cor_p_value(pr, ncond);
      spearman_p(i, j) = cor_p_value(sr, ncond);
    }
  }

  std::vector<std::vector<double> > pearson_adj(ntf), spearman_adj(ntf);
  for (int j = 0; j < ntf; ++j) {
    pearson_adj[j] = bh_adjust_column(pearson_p, j);
    spearman_adj[j] = bh_adjust_column(spearman_p, j);
  }

  std::size_t n_pass = 0;
  for (int j = 0; j < ntf; ++j) {
    for (int i = 0; i < nfp; ++i) {
      double best = pearson(i, j);
      bool use_s = false;
      const double sr = spearman(i, j);
      if (!R_finite(best) || (R_finite(sr) && sr > best)) {
        best = sr;
        use_s = true;
      }
      const double bp = use_s ? spearman_p(i, j) : pearson_p(i, j);
      const double bf = use_s ? spearman_adj[j][i] : pearson_adj[j][i];
      const bool pass_r = R_finite(best) && best >= r_cutoff;
      const bool pass_p = !R_finite(p_cutoff) || (R_finite(bp) && bp <= p_cutoff);
      const bool pass_f = !R_finite(fdr_cutoff) || (R_finite(bf) && bf <= fdr_cutoff);
      if (pass_r && pass_p && pass_f) ++n_pass;
    }
  }

  CharacterVector out_fp(n_pass), out_peak(n_pass), out_tf(n_pass), out_method(n_pass);
  NumericVector out_pr(n_pass), out_pp(n_pass), out_ppa(n_pass), out_sr(n_pass), out_sp(n_pass), out_spa(n_pass);
  NumericVector out_best(n_pass), out_bp(n_pass), out_bf(n_pass);
  LogicalVector out_pass_r(n_pass), out_pass_p(n_pass), out_pass_f(n_pass), out_pass(n_pass);
  std::size_t k = 0;
  for (int j = 0; j < ntf; ++j) {
    for (int i = 0; i < nfp; ++i) {
      double best = pearson(i, j);
      bool use_s = false;
      const double sr = spearman(i, j);
      if (!R_finite(best) || (R_finite(sr) && sr > best)) {
        best = sr;
        use_s = true;
      }
      const double bp = use_s ? spearman_p(i, j) : pearson_p(i, j);
      const double bf = use_s ? spearman_adj[j][i] : pearson_adj[j][i];
      const bool pass_r = R_finite(best) && best >= r_cutoff;
      const bool pass_p = !R_finite(p_cutoff) || (R_finite(bp) && bp <= p_cutoff);
      const bool pass_f = !R_finite(fdr_cutoff) || (R_finite(bf) && bf <= fdr_cutoff);
      if (!(pass_r && pass_p && pass_f)) continue;
      out_fp[k] = fp_id[i];
      out_peak[k] = atac_peak[i];
      out_tf[k] = tf_name[j];
      out_pr[k] = pearson(i, j);
      out_pp[k] = pearson_p(i, j);
      out_ppa[k] = pearson_adj[j][i];
      out_sr[k] = spearman(i, j);
      out_sp[k] = spearman_p(i, j);
      out_spa[k] = spearman_adj[j][i];
      out_best[k] = best;
      out_method[k] = use_s ? "spearman" : "pearson";
      out_bp[k] = bp;
      out_bf[k] = bf;
      out_pass_r[k] = pass_r;
      out_pass_p[k] = pass_p;
      out_pass_f[k] = pass_f;
      out_pass[k] = true;
      ++k;
    }
  }
  return DataFrame::create(
    _["fp_id"] = out_fp,
    _["atac_peak"] = out_peak,
    _["tf"] = out_tf,
    _["pearson_r"] = out_pr,
    _["pearson_p"] = out_pp,
    _["pearson_p_adj"] = out_ppa,
    _["spearman_r"] = out_sr,
    _["spearman_p"] = out_sp,
    _["spearman_p_adj"] = out_spa,
    _["best_r"] = out_best,
    _["best_method"] = out_method,
    _["best_p"] = out_bp,
    _["best_fdr"] = out_bf,
    _["pass_r"] = out_pass_r,
    _["pass_p"] = out_pass_p,
    _["pass_fdr"] = out_pass_f,
    _["pass"] = out_pass,
    _["stringsAsFactors"] = false
  );
}
