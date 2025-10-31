#' @import data.table
#' @importFrom methods as
#' @importFrom stats aggregate ave end na.omit setNames start
#' @importFrom utils data head modifyList
NULL

utils::globalVariables(c(
  ":", ":=", ".", ".N", ".data",
  "..cond_cols", "..have", "..keep_cols", "..sample_cols",
  "any_bound", "ATAC_bound", "atac_peak", "ATAC_score",
  "bed_end", "bed_start", "bound_peaks", "chr", "comp_id",
  "condition", "corr_fp_tf_p", "corr_fp_tf_p_adj", "corr_fp_tf_r",
  "e_exp", "expressed_genes", "focal_point", "fp_peak", "fp_peak_bak",
  "group_size", "i.atac_peak", "i.group_size", "i.new_peak_ID",
  "metric", "mid", "motifs", "N", "new_end", "new_peak_ID",
  "new_start", "peak_ATAC", "peak_chr", "peak_end", "peak_ID",
  "peak_ID_old", "peak_start", "pearson_corr_fp_tf_r", "s_exp",
  "sample_id", "score", "sig_id", "spearman_corr_fp_tf_r",
  "TFBS_chr", "TFBS_end", "TFBS_name", "TFBS_start", "tfs",
  "type", "value"
))
