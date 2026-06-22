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
  "e_exp", "expressed_genes", "focal_point", "fp_peak", "source_fp_peak",
  "group_size", "i.atac_peak", "i.group_size", "i.new_peak_ID",
  "metric", "mid", "motif_use", "motifs", "N", "new_end", "new_peak_ID",
  "new_start", "peak_ATAC", "peak_chr", "peak_end", "peak_ID",
  "peak_ID_old", "peak_start", "pearson_corr_fp_tf_r", "s_exp",
  "sample_id", "score", "sig_id", "spearman_corr_fp_tf_r",
  "TFBS_chr", "TFBS_end", "TFBS_name", "TFBS_start", "tfbs_id", "tfs", "tfs_clean",
  "type", "value"
))

utils::globalVariables(c(
  "..from_cols", "..keep", "..topic_cols", "case_id", "category",
  "cellline", "cluster_label", "cluster_size", "col_index",
  "comparison_id", "comparison_label", "condition_group", "condition_label", "context_type",
  "cond1_id", "cond2_id", "cond1_label", "cond2_label", "cond1_matrix_id",
  "cond2_matrix_id", "ctrl_id", "d2_loglik",
  "active_any", "active_both", "active_cond1", "active_cond2",
  "abs_net_target_gene_count",
  "comparison_group",
  "delta_abs", "delta_fp", "delta_fp_bed_score", "delta_fp_score",
  "delta_link_score",
  "delta_gene", "delta_gene_expr", "de_source", "de_test_id",
  "direction_label", "display_label", "doc_id", "dominant_direction", "down",
  "doc_idx", "doc_weight_total", "dot_top_n_per_topic", "fail",
  "family", "fc_fp_bed_score", "fc_fp_score", "fc_gene_expr", "fc_mag_fp",
  "fc_mag_gene", "fc_mag_tf", "fill_group", "fit_file",
  "fp_bound_case", "fp_bound_cond1", "fp_bound_cond2", "fp_bound_ctrl",
  "fp_dir_sum", "fp_id",
  "fp_score_case", "fp_score_cond1", "fp_score_condition", "fp_score_cond2",
  "fp_score_ctrl",
  "fp_weight", "fp_weight_raw", "gene", "gene_expr_case",
  "gene_expr_cond1", "gene_expr_condition", "gene_expr_cond2",
  "gene_expr_ctrl", "gene_expr_flag_case", "gene_expr_flag_cond1",
  "gene_expr_flag_cond2", "gene_expr_flag_ctrl", "gene_gamma_cutoff", "gene_idx", "gene_name",
  "gene_in", "gene_out", "gene_pass", "gene_prob", "gene_prob_pass", "gene_scale",
  "gene_score", "gene_symbol", "gene_term", "gene_term_id",
  "gene_total", "gene_total_universe", "gene_total_universe_key",
  "gene_weight", "genes", "group_label", "gs_name", "has_matching_peak", "HGNC",
  "in_marker_heatmap", "in_topic", "inner_distance", "inter_distance", "is_sig", "is_topic_top_pathway",
  "item",
  "K", "leading_edge", "link_efdr_p", "link_efdr_q", "link_id",
  "candidate_id", "candidate_source", "distance_to_tss", "fp_center", "point_end", "point_start", "prior_id", "prior_score", "prior_source", "prior_status", "prior_supported", "target_chr", "target_gene", "target_strand", "target_tss", "window_end", "window_start", "within_tss_window",
  "link_pass", "module2_link_pass", "link_score", "link_score_prob", "log2_fc_z_c",
  "log2fc_fp", "log2FC_fp_bed_score", "log2FC_fp_score",
  "log2fc_gene", "log2fc_rna", "log2fc_tf", "log2FC_gene_expr",
  "log2FC_gene_expr_direct", "log2FC_tf_expr", "log2FC_tf_expr_direct", "loglik", "logp",
  "marker_row_order", "max_doc_score", "max_score", "median_log2FC_tf_expr", "method", "method_order",
  "method_setup", "metric_group", "min_topic",
  "model_label", "modality", "n_down_total", "N_gene", "N_peak", "n_target_genes_down",
  "n_target_genes_up", "n_tokens",
  "n_topics", "n_up_total", "net_target_gene_count", "NES", "other_gamma_cutoff",
  "owner_topic", "padj", "padj_rna", "page_label", "panel", "pass", "path",
  "pathway_make_dotplot", "pathway_order", "peak", "peak_id",
  "peak_idx", "peak_key", "peak_pass", "peak_scale", "peak_score",
  "peak_term", "peaks_gamma_cutoff", "perplexity",
  "plot_label", "plot_tf_network_delta", "prob", "pseudo_count", "pseudo_count_bin",
  "p_adj_rna_de_gene", "p_adj_rna_de_tf", "pseudo_count_log", "pval",
  "pvalue_rna", "p_rna_de_gene", "p_rna_de_tf", "rid___", "row_index", "row_order",
  "run_id", "sampler", "score_num", "score_val", "size_val", "sub_cluster_name",
  "target_total", "term_base", "term_id", "term_pass", "term_prefix",
  "term_str", "tf", "tf_doc", "tf_expr_case", "tf_expr_cond1",
  "tf_expr_condition", "tf_expr_cond2", "tf_expr_ctrl", "tf_expr_flag_case",
  "tf_expr_flag_cond1", "tf_expr_flag_cond2", "tf_expr_flag_ctrl", "tf_idx",
  "selected_k", "setup", "setup_label", "shape_value", "split_panel_label",
  "tf_key", "tf_sum_abs_delta", "tf_sum_delta", "tf_term",
  "tf_weight", "tf_weight_norm", "topic", "topic_label", "topic_max",
  "topic_mean_abs_delta", "topic_n_links", "topic_n_target_genes",
  "theta_condition_label_score", "theta_condition_separation_score", "theta_mean",
  "abs_edge_score", "best_distance_to_tss", "best_peak_ID",
  "color", "color_value", "edge_rank", "edge_score", "edge_score_row", "edge_r",
  "from", "max_abs_r_gene", "max_abs_r_rna_gene",
  "label_mode", "mds_label", "median_log2FC_gene_expr", "MDS1", "MDS2",
  "n_edges", "n_genes", "nearest_scaled_dist",
  "n_supporting_links", "n_tfs", "r_gene", "r_rna_gene",
  "scaled_x", "scaled_y", "tf_edge_score_sum", "tf_rank", "tf_target_count", "to",
  "topic_n_TFs", "topic_num", "topic_rank", "total", "value_diff",
  "up", "value_max", "value_min", "value_raw", "value_tol", "w", "weight_raw"
))

utils::globalVariables(c(
  "bar_group", "base_panel_label", "comparison_display",
  "comparison_display_from_edges", "comparison_label_manifest", "count",
  "count_basis", "design_display_label", "doc_design", "doc_display_label",
  "doc_display_label_from_edges", "fraction", "gammafit_scope_label",
  "k_label", "method_key", "method_label", "method_setup_label", "model_k",
  "model_k_short", "n_items", "panel_label", "src", "status", "unit"
))
