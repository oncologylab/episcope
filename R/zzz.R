#' @import data.table
#' @useDynLib craftgrn, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics mtext symbols
#' @importFrom methods as
#' @importFrom stats aggregate ave end na.omit setNames start
#' @importFrom utils data head modifyList
NULL

utils::globalVariables(c(
  ":", ":=", ".", ".N", ".data",
  "..cond_cols", "..have", "..keep_cols", "..sample_cols",
  "any_bound", "ATAC_bound", "atac_peak", "ATAC_score",
  "bed_end", "bed_start", "bound_peaks", "chr", "comp_id",
  "condition", "condition_support", "corr_fp_tf_p", "corr_fp_tf_p_adj", "corr_fp_tf_r",
  "e_exp", "expressed_genes", "focal_point", "fp_peak", "fp_peak_bak",
  "group_size", "i.atac_peak", "i.group_size", "i.new_peak_ID",
  "metric", "mid", "motif_use", "motifs", "N", "new_end", "new_peak_ID",
  "new_start", "peak_ATAC", "peak_chr", "peak_end", "peak_ID",
  "peak_ID_old", "peak_start", "pearson_corr_fp_tf_r", "s_exp",
  "sample_id", "score", "sig_id", "spearman_corr_fp_tf_r",
  "TFBS_chr", "TFBS_end", "TFBS_name", "TFBS_start", "tfs", "tfs_clean",
  "type", "value"
))

utils::globalVariables(c(
  "..from_cols", "..keep", "..topic_cols", "case_id", "category",
  "cellline", "cluster_label", "cluster_size", "col_index",
  "comparison_label", "condition_label", "ctrl_id", "d2_loglik",
  "delta_abs", "delta_fp", "delta_fp_bed_score", "delta_fp_score",
  "delta_gene", "delta_gene_expr", "direction_label", "doc_id",
  "doc_idx", "doc_weight_total", "dot_top_n_per_topic", "fail",
  "fc_fp_bed_score", "fc_fp_score", "fc_gene_expr", "fc_mag_fp",
  "fc_mag_gene", "fc_mag_tf", "fill_group", "fit_file",
  "fp_bound_case", "fp_bound_ctrl", "fp_dir_sum", "fp_id",
  "fp_score_case", "fp_score_condition", "fp_score_ctrl",
  "fp_weight", "fp_weight_raw", "gene", "gene_expr_case",
  "gene_expr_condition", "gene_expr_ctrl", "gene_expr_flag_case",
  "gene_expr_flag_ctrl", "gene_gamma_cutoff", "gene_idx", "gene_name",
  "gene_pass", "gene_prob", "gene_prob_pass", "gene_scale",
  "gene_score", "gene_symbol", "gene_term", "gene_term_id",
  "gene_weight", "gs_name", "has_matching_peak", "HGNC",
  "in_marker_heatmap", "in_topic", "is_sig", "is_topic_top_pathway",
  "K", "leading_edge", "link_efdr_p", "link_efdr_q", "link_id",
  "link_pass", "link_score", "link_score_prob", "log2_fc_z_c",
  "log2fc_fp", "log2FC_fp_bed_score", "log2FC_fp_score",
  "log2fc_gene", "log2fc_tf", "log2FC_tf_expr", "loglik", "logp",
  "marker_row_order", "max_doc_score", "max_score", "min_topic",
  "modality", "n_down_total", "N_gene", "N_peak", "n_tokens",
  "n_topics", "n_up_total", "NES", "other_gamma_cutoff",
  "owner_topic", "padj", "page_label", "panel", "pass",
  "pathway_make_dotplot", "pathway_order", "peak", "peak_id",
  "peak_idx", "peak_key", "peak_pass", "peak_scale", "peak_score",
  "peak_term", "peaks_gamma_cutoff", "perplexity",
  "plot_tf_network_delta", "prob", "pseudo_count", "pseudo_count_bin",
  "pseudo_count_log", "pval", "rid___", "row_index", "row_order",
  "run_id", "score_num", "score_val", "size_val", "sub_cluster_name",
  "target_total", "term_base", "term_id", "term_pass", "term_prefix",
  "term_str", "tf", "tf_doc", "tf_expr_case", "tf_expr_condition",
  "tf_expr_ctrl", "tf_expr_flag_case", "tf_expr_flag_ctrl", "tf_idx",
  "tf_key", "tf_sum_abs_delta", "tf_sum_delta", "tf_term",
  "tf_weight", "tf_weight_norm", "topic", "topic_label", "topic_max",
  "topic_mean_abs_delta", "topic_n_links", "topic_n_target_genes",
  "topic_n_TFs", "topic_num", "topic_rank", "total", "value_diff",
  "value_max", "value_min", "value_raw", "value_tol", "w", "weight_raw"
))
