// Shared-nearest-neighbor edges for experimental Module 3 benchmarks.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export(.topic_snn_edges_cpp)]]
DataFrame topic_snn_edges_cpp(const IntegerMatrix& neighbor_index) {
  const int n = neighbor_index.nrow();
  const int k = neighbor_index.ncol();
  std::vector<int> from;
  std::vector<int> to;
  std::vector<int> shared;
  std::vector<double> weight;
  std::vector<int> marks(n, -1);

  from.reserve(static_cast<std::size_t>(n) * k / 3);
  to.reserve(static_cast<std::size_t>(n) * k / 3);
  shared.reserve(static_cast<std::size_t>(n) * k / 3);
  weight.reserve(static_cast<std::size_t>(n) * k / 3);

  for (int i = 0; i < n; ++i) {
    marks[i] = i;
    for (int a = 0; a < k; ++a) {
      const int neighbor = neighbor_index(i, a) - 1;
      if (neighbor >= 0 && neighbor < n) marks[neighbor] = i;
    }
    for (int a = 0; a < k; ++a) {
      const int j = neighbor_index(i, a) - 1;
      if (j <= i || j >= n) continue;

      bool mutual = false;
      for (int b = 0; b < k; ++b) {
        if (neighbor_index(j, b) - 1 == i) {
          mutual = true;
          break;
        }
      }
      if (!mutual) continue;

      int overlap = 1;
      for (int b = 0; b < k; ++b) {
        const int candidate = neighbor_index(j, b) - 1;
        if (candidate >= 0 && candidate < n && marks[candidate] == i) {
          ++overlap;
        }
      }
      const int set_union = 2 * (k + 1) - overlap;
      const double jaccard = set_union > 0
        ? static_cast<double>(overlap) / static_cast<double>(set_union)
        : 0.0;
      from.push_back(i + 1);
      to.push_back(j + 1);
      shared.push_back(overlap);
      weight.push_back(jaccard);
    }
  }

  return DataFrame::create(
    _["from"] = from,
    _["to"] = to,
    _["shared_neighbors"] = shared,
    _["weight"] = weight
  );
}

// [[Rcpp::export(.topic_finalize_sparse_counts_cpp)]]
List topic_finalize_sparse_counts_cpp(const IntegerVector& row_pointer,
                                      const IntegerVector& column_index,
                                      const NumericVector& source_count,
                                      const LogicalVector& gene_term,
                                      const NumericVector& idf_multiplier,
                                      const double peak_gene_ratio = NA_REAL) {
  const int n_rows = row_pointer.size() - 1;
  const int n_terms = gene_term.size();
  if (column_index.size() != source_count.size()) {
    stop("Sparse column and count vectors must have equal length.");
  }
  if (idf_multiplier.size() != n_terms) {
    stop("IDF multipliers must have one value per term.");
  }
  const bool balance_modalities = R_finite(peak_gene_ratio);
  if (balance_modalities && peak_gene_ratio <= 0) {
    stop("Peak/Gene token ratio must be positive.");
  }
  NumericVector output(source_count.size());
  double source_total = 0;
  double final_total = 0;
  double gene_total = 0;
  double peak_total = 0;

  auto allocate_group = [&](const int begin,
                            const int end,
                            const std::vector<double>& score,
                            const std::vector<unsigned char>& selected,
                            const long long target) {
    int n_selected = 0;
    double score_total = 0;
    for (int offset = 0; offset < end - begin; ++offset) {
      if (!selected[offset]) continue;
      ++n_selected;
      score_total += score[offset];
    }
    if (target < n_selected) {
      stop("Token target cannot retain every sparse term.");
    }
    const long long remaining = target - n_selected;
    double cumulative = 0;
    long long previous_floor = 0;
    long long allocated_total = 0;
    int maximum_position = -1;
    double maximum_count = -1;
    for (int offset = 0; offset < end - begin; ++offset) {
      if (!selected[offset]) continue;
      const double share = score_total > 0
        ? score[offset] / score_total
        : 1.0 / static_cast<double>(n_selected);
      const double quota = share * static_cast<double>(remaining);
      cumulative += quota - std::floor(quota);
      const long long cumulative_floor = static_cast<long long>(
        std::floor(cumulative + 1e-8)
      );
      const long long residual = cumulative_floor - previous_floor;
      previous_floor = cumulative_floor;
      const double count = 1.0 + std::floor(quota) + residual;
      output[begin + offset] = count;
      allocated_total += static_cast<long long>(count);
      if (count > maximum_count) {
        maximum_count = count;
        maximum_position = begin + offset;
      }
    }
    const long long correction = target - allocated_total;
    if (correction != 0 && maximum_position >= 0) {
      output[maximum_position] += static_cast<double>(correction);
    }
  };

  for (int row = 0; row < n_rows; ++row) {
    const int begin = row_pointer[row];
    const int end = row_pointer[row + 1];
    const int n = end - begin;
    std::vector<double> score(n, 0);
    std::vector<unsigned char> all(n, 1);
    std::vector<unsigned char> gene(n, 0);
    std::vector<unsigned char> peak(n, 0);
    long long row_total = 0;
    int n_gene = 0;
    int n_peak = 0;
    for (int offset = 0; offset < n; ++offset) {
      const int position = begin + offset;
      const int term = column_index[position];
      if (term < 0 || term >= n_terms) stop("Sparse term index is invalid.");
      const double value = source_count[position];
      if (!R_finite(value) || value <= 0) {
        stop("Sparse model counts must be positive and finite.");
      }
      row_total += static_cast<long long>(std::llround(value));
      score[offset] = value * idf_multiplier[term];
      if (gene_term[term] == TRUE) {
        gene[offset] = 1;
        ++n_gene;
      } else {
        peak[offset] = 1;
        ++n_peak;
      }
    }
    source_total += static_cast<double>(row_total);
    if (!balance_modalities) {
      allocate_group(begin, end, score, all, row_total);
    } else {
      if (n_gene == 0 || n_peak == 0) {
        stop("Every document must contain Gene and Peak terms.");
      }
      long long target_peak = static_cast<long long>(std::llround(
        static_cast<double>(row_total) * peak_gene_ratio /
          (1.0 + peak_gene_ratio)
      ));
      target_peak = std::max<long long>(
        n_peak,
        std::min<long long>(row_total - n_gene, target_peak)
      );
      allocate_group(begin, end, score, gene, row_total - target_peak);
      allocate_group(begin, end, score, peak, target_peak);
    }
    for (int position = begin; position < end; ++position) {
      const int term = column_index[position];
      final_total += output[position];
      if (gene_term[term] == TRUE) gene_total += output[position];
      else peak_total += output[position];
    }
  }

  return List::create(
    _["counts"] = output,
    _["source_tokens"] = source_total,
    _["final_tokens"] = final_total,
    _["gene_tokens"] = gene_total,
    _["peak_tokens"] = peak_total
  );
}
