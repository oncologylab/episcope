# File: utils_step1_qc.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 QC plots for the compact Step 1 workflow.
#
# Inputs:
# - prepared compact Step 1 tables
# - output directories, thresholds, and labeling information
#
# Outputs:
# - footprint merge, normalization/bound, and gene-expression QC PDFs
#
# Notes:
# - Keep plotting and output reporting separate from core data transformations.

#' Module 1 QC helpers
#'
#' @noRd
NULL

plot_fp_merge_summary <- function(
    fp_aligned,
    out_dir,
    db,
    verbose = TRUE
) {
  if (!is.list(fp_aligned) || !is.data.frame(fp_aligned$id_map)) {
    .log_abort("`fp_aligned$id_map` is missing or invalid.")
  }
  if (!is.character(out_dir) || !nzchar(out_dir)) {
    .log_abort("`out_dir` must be a non-empty path.")
  }
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  id_map <- fp_aligned$id_map
  if (!"group_size" %in% names(id_map)) {
    id_map$group_size <- 1L
  }
  id_map$group_size <- as.numeric(id_map$group_size)

  counts_df <- data.frame(
    metric = c("Raw footprints", "Aligned footprints", "ATAC peaks"),
    value = c(
      length(unique(id_map$source_fp_peak)),
      length(unique(id_map$peak_ID)),
      length(unique(id_map$atac_peak))
    ),
    stringsAsFactors = FALSE
  )
  counts_df$metric <- factor(counts_df$metric, levels = counts_df$metric)

  base_theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, hjust = 1)
    )

  p_counts <- ggplot2::ggplot(counts_df, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_col(fill = "#3182bd") +
    ggplot2::labs(
      title = "Footprint merge: peak counts (raw vs aligned vs ATAC)",
      x = "Entity",
      y = "Count",
      caption = "Counts derived from fp_aligned$id_map."
    ) +
    base_theme

  p_hist <- ggplot2::ggplot(id_map, ggplot2::aes(x = group_size)) +
    ggplot2::geom_histogram(bins = 50, fill = "#6baed6", color = "white") +
    ggplot2::labs(
      title = "Footprint merge: group size distribution",
      x = "Merged group size (raw footprints per aligned peak)",
      y = "Aligned peaks",
      caption = "group_size from fp_aligned$id_map."
    ) +
    base_theme

  p_ecdf <- ggplot2::ggplot(id_map, ggplot2::aes(x = group_size)) +
    ggplot2::stat_ecdf(geom = "step", linewidth = 0.8, color = "#08519c") +
    ggplot2::labs(
      title = "Footprint merge: ECDF of group size",
      x = "Merged group size (raw footprints per aligned peak)",
      y = "ECDF",
      caption = "ECDF of group_size from fp_aligned$id_map."
    ) +
    base_theme

  pdf_path <- file.path(out_dir, "02_fp_merge_summary.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 6)
  print(p_counts)
  print(p_hist)
  print(p_ecdf)
  grDevices::dev.off()

  if (isTRUE(verbose)) .log_inform("Footprint merge summary saved: {pdf_path}")
  invisible(pdf_path)
}

plot_fp_norm_bound_qc <- function(
    omics_data = NULL,
    grn_set = NULL,
    out_dir,
    db,
    threshold_fp_score,
    max_points = 100000L,
    verbose = TRUE
) {
  if (is.null(omics_data)) omics_data <- grn_set
  if (!is.list(omics_data) || !is.data.frame(omics_data$fp_score_condition)) {
    .log_abort("`omics_data$fp_score_condition` is missing or invalid.")
  }
  if (!is.data.frame(omics_data$fp_score_condition_qn)) {
    .log_abort("`omics_data$fp_score_condition_qn` is missing; run grn_add_fp_score_qn() first.")
  }
  if (!is.data.frame(omics_data$fp_bound_condition)) {
    .log_abort("`omics_data$fp_bound_condition` is missing; run grn_add_fp_bound_condition() first.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fp_raw <- omics_data$fp_score_condition
  fp_qn <- omics_data$fp_score_condition_qn
  n_peaks <- nrow(fp_qn)
  if (!n_peaks) .log_abort("`fp_score_condition_qn` has zero rows.")

  keep_n <- min(as.integer(max_points), n_peaks)
  keep_idx <- if (keep_n >= n_peaks) {
    seq_len(n_peaks)
  } else {
    unique(pmax(1L, pmin(n_peaks, round(seq(1, n_peaks, length.out = keep_n)))))
  }

  fp_raw_sub <- fp_raw[keep_idx, , drop = FALSE]
  fp_qn_sub <- fp_qn[keep_idx, , drop = FALSE]

  fp_raw_long <- tidyr::pivot_longer(fp_raw_sub, cols = -peak_ID, names_to = "condition", values_to = "score")
  fp_raw_long$type <- "raw"
  fp_qn_long <- tidyr::pivot_longer(fp_qn_sub, cols = -peak_ID, names_to = "condition", values_to = "score")
  fp_qn_long$type <- "quantile_normalized"
  fp_long <- rbind(fp_raw_long, fp_qn_long)

  bound_tbl <- omics_data$fp_bound_condition
  bound_mat <- as.matrix(bound_tbl[, setdiff(names(bound_tbl), "peak_ID"), drop = FALSE])
  bound_counts <- colSums(bound_mat > 0L, na.rm = TRUE)
  bound_df <- data.frame(
    condition = names(bound_counts),
    bound_peaks = as.integer(bound_counts),
    stringsAsFactors = FALSE
  )

  base_theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, hjust = 1)
    )

  p_violin <- ggplot2::ggplot(fp_long, ggplot2::aes(x = condition, y = score, fill = type)) +
    ggplot2::geom_violin(
      position = ggplot2::position_dodge(width = 0.8),
      trim = TRUE,
      scale = "width",
      alpha = 0.85,
      color = "grey30",
      linewidth = 0.25
    ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.12,
      outlier.size = 0.1,
      linewidth = 0.2,
      alpha = 0.45
    ) +
    ggplot2::labs(
      title = "Footprint score distributions by condition (raw vs quantile-normalized)",
      x = "Condition",
      y = "Footprint score",
      caption = paste0("Showing ", keep_n, " of ", n_peaks, " peaks (sampled).")
    ) +
    base_theme +
    ggplot2::coord_flip()

  p_bound <- ggplot2::ggplot(bound_df, ggplot2::aes(x = condition, y = bound_peaks)) +
    ggplot2::geom_col(fill = "#3182bd") +
    ggplot2::labs(
      title = "Total bound footprints per condition",
      x = "Condition",
      y = "Bound footprints",
      caption = paste0("Bound: fp_score >= ", threshold_fp_score, ".")
    ) +
    base_theme +
    ggplot2::coord_flip()

  pdf_path <- file.path(out_dir, "03_fp_norm_bound_summary.pdf")
  grDevices::pdf(pdf_path, width = 7, height = 10)
  print(p_violin)
  print(p_bound)
  grDevices::dev.off()

  if (isTRUE(verbose)) .log_inform("FP normalization/bound summary saved: {pdf_path}")
  invisible(pdf_path)
}

plot_gene_expr_qc <- function(
    omics_data = NULL,
    grn_set = NULL,
    out_dir,
    db,
    threshold_gene_expr,
    verbose = TRUE
) {
  if (is.null(omics_data)) omics_data <- grn_set
  if (!is.list(omics_data) || !is.data.frame(omics_data$rna_expressed)) {
    .log_abort("`omics_data$rna_expressed` is missing or invalid.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expr_tbl <- omics_data$rna_expressed
  expr_mat <- as.matrix(expr_tbl[, setdiff(names(expr_tbl), c("ensembl_gene_id", "HGNC")), drop = FALSE])
  expr_counts <- colSums(expr_mat > 0L, na.rm = TRUE)
  expr_df <- data.frame(
    condition = names(expr_counts),
    expressed_genes = as.integer(expr_counts),
    stringsAsFactors = FALSE
  )

  base_theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, hjust = 1)
    )

  p_expr <- ggplot2::ggplot(expr_df, ggplot2::aes(x = condition, y = expressed_genes)) +
    ggplot2::geom_col(fill = "#31a354") +
    ggplot2::labs(
      title = "Total expressed genes per condition",
      x = "Condition",
      y = "Expressed genes",
      caption = paste0("Expressed: RNA >= ", threshold_gene_expr, ".")
    ) +
    base_theme +
    ggplot2::coord_flip()

  pdf_path <- file.path(out_dir, "05_gene_expr_flag_summary.pdf")
  grDevices::pdf(pdf_path, width = 7, height = 9)
  print(p_expr)
  grDevices::dev.off()

  if (isTRUE(verbose)) .log_inform("Gene expression QC saved: {pdf_path}")
  invisible(pdf_path)
}
