# File: utils_step1_qc.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 QC plots and numbered Step 1 output writers.
#
# Inputs:
# - prepared Step 1 tables and correlation outputs
# - output directories, thresholds, and labeling information
#
# Outputs:
# - numbered Step 1 CSV outputs
# - PDF QC plots
# - Step 1 summary tables and overview files
# - grouped Module 1 motif-clustering outputs under 04_tf_motif_clustering/
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
      length(unique(id_map$fp_peak_bak)),
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

write_grn_outputs <- function(grn_set, out_dir, db, qn_base_dir = NULL, verbose = TRUE) {
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.data.frame(grn_set$fp_score)) .log_abort("`grn_set$fp_score` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_annotation)) .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  if (!is.data.frame(grn_set$rna_expressed)) .log_abort("`grn_set$rna_expressed` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_bound_condition)) .log_abort("`grn_set$fp_bound_condition` is missing or invalid.")

  fp_score_out <- grn_set$fp_score
  if ("peak_ID" %in% names(fp_score_out)) {
    fp_score_out <- dplyr::distinct(fp_score_out, .data$peak_ID, .keep_all = TRUE)
  }
  readr::write_csv(fp_score_out, file.path(out_dir, sprintf("02_fp_score_raw_%s.csv", db)))
  readr::write_csv(grn_set$fp_annotation, file.path(out_dir, sprintf("03_fp_annotation_%s.csv", db)))
  readr::write_csv(grn_set$rna_expressed, file.path(out_dir, sprintf("05_gene_expr_flag_%s.csv", db)))

  fp_bound_out <- grn_set$fp_bound_condition
  if ("peak_ID" %in% names(fp_bound_out)) {
    fp_bound_out <- dplyr::distinct(fp_bound_out, .data$peak_ID, .keep_all = TRUE)
  }

  if (is.data.frame(grn_set$fp_score_condition_qn)) {
    fp_qn_out <- grn_set$fp_score_condition_qn
    if ("peak_ID" %in% names(fp_qn_out)) {
      fp_qn_out <- dplyr::distinct(fp_qn_out, .data$peak_ID, .keep_all = TRUE)
    }
    if ("peak_ID" %in% names(fp_qn_out) && "peak_ID" %in% names(fp_bound_out)) {
      idx_bound <- match(fp_qn_out$peak_ID, fp_bound_out$peak_ID)
      if (anyNA(idx_bound)) {
        .log_warn("Some fp_bound_condition peaks are missing in fp_score_qn; dropping unmatched rows.")
        keep <- which(!is.na(idx_bound))
        fp_qn_out <- fp_qn_out[keep, , drop = FALSE]
        idx_bound <- idx_bound[keep]
      }
      fp_bound_out <- fp_bound_out[idx_bound, , drop = FALSE]
    }
    readr::write_csv(fp_bound_out, file.path(out_dir, sprintf("03_fp_bound_condition_%s.csv", db)))
    readr::write_csv(fp_qn_out, file.path(out_dir, sprintf("03_fp_score_qn_%s.csv", db)))
    if (!is.null(qn_base_dir) && is.character(qn_base_dir) && nzchar(qn_base_dir)) {
      dir.create(qn_base_dir, recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(fp_qn_out, file.path(qn_base_dir, sprintf("fp_scores_qn_%s.csv", db)))
    }
  }

  if (isTRUE(verbose)) {
    .log_inform("Module 1 numbered prep outputs saved in {.path {out_dir}}")
  }
  invisible(out_dir)
}

#' Write numbered Module 1 preprocessing outputs
#'
#' Write the standard numbered Module 1 preprocessing outputs from a prepared
#' multi-omic object without rerunning preprocessing. This includes the raw
#' footprint score table, normalized footprint outputs, gene-expression flag
#' table, and the associated QC PDFs.
#'
#' @param omics_data Prepared Module 1 multi-omic object.
#' @param out_dir Output directory, typically `predict_tf_binding_sites/`.
#' @param db Database label used in output file names.
#' @param threshold_fp_score Footprint-score threshold used when creating the
#'   condition-specific bound matrix.
#' @param threshold_gene_expr Expression threshold used when creating the
#'   condition-specific gene-expression flag matrix.
#' @param fp_cache_dir Optional aligned-footprint cache directory used to
#'   regenerate `02_fp_merge_summary.pdf` when it is missing.
#' @param fp_cache_tag Optional cache tag used with `fp_cache_dir`.
#' @param output_mode Cached aligned-footprint output mode used with
#'   `load_fp_aligned_from_cache()`.
#' @param cluster_dir Optional source directory for motif-clustering outputs.
#'   When supplied, `04_*` files are copied into
#'   `file.path(out_dir, "04_tf_motif_clustering")`.
#' @param verbose Logical; if `TRUE`, emit concise progress messages.
#'
#' @return Output directory path (invisible).
#' @keywords internal
write_multiomic_prep_outputs <- function(
    omics_data,
    out_dir,
    db,
    threshold_fp_score,
    threshold_gene_expr,
    fp_cache_dir = NULL,
    fp_cache_tag = NULL,
    output_mode = "distinct",
    cluster_dir = NULL,
    verbose = TRUE
) {
  if (!is.list(omics_data)) {
    .log_abort("`omics_data` must be a list.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_subdir <- file.path(out_dir, "04_tf_motif_clustering")
  required_files <- c(
    file.path(out_dir, sprintf("02_fp_score_raw_%s.csv", db)),
    file.path(out_dir, "02_fp_merge_summary.pdf"),
    file.path(out_dir, sprintf("03_fp_annotation_%s.csv", db)),
    file.path(out_dir, sprintf("03_fp_bound_condition_%s.csv", db)),
    file.path(out_dir, sprintf("03_fp_score_qn_%s.csv", db)),
    file.path(out_dir, "03_fp_norm_bound_summary.pdf"),
    file.path(out_dir, sprintf("05_gene_expr_flag_%s.csv", db)),
    file.path(out_dir, "05_gene_expr_flag_summary.pdf")
  )
  has_04_outputs <- dir.exists(cluster_subdir) &&
    length(list.files(cluster_subdir, pattern = "^04_", full.names = TRUE)) > 0L
  if (all(file.exists(required_files)) && has_04_outputs) {
    if (isTRUE(verbose)) {
      .log_inform("Module 1 numbered prep outputs already present in {.path {out_dir}}")
    }
    return(invisible(out_dir))
  }

  write_grn_outputs(
    grn_set = omics_data,
    out_dir = out_dir,
    db = db,
    verbose = verbose
  )

  merge_summary_path <- file.path(out_dir, "02_fp_merge_summary.pdf")
  if (!file.exists(merge_summary_path) &&
      is.character(fp_cache_dir) && nzchar(fp_cache_dir) &&
      is.character(fp_cache_tag) && nzchar(fp_cache_tag) &&
      exists("load_fp_aligned_from_cache", mode = "function")) {
    fp_aligned <- load_fp_aligned_from_cache(
      cache_dir = fp_cache_dir,
      cache_tag = fp_cache_tag,
      output_mode = output_mode,
      verbose = FALSE
    )
    if (is.list(fp_aligned) && is.data.frame(fp_aligned$id_map)) {
      plot_fp_merge_summary(
        fp_aligned = fp_aligned,
        out_dir = out_dir,
        db = db,
        verbose = verbose
      )
    }
  }

  plot_fp_norm_bound_qc(
    omics_data = omics_data,
    out_dir = out_dir,
    db = db,
    threshold_fp_score = threshold_fp_score,
    verbose = verbose
  )
  plot_gene_expr_qc(
    omics_data = omics_data,
    out_dir = out_dir,
    db = db,
    threshold_gene_expr = threshold_gene_expr,
    verbose = verbose
  )

  if (is.character(cluster_dir) && nzchar(cluster_dir) && dir.exists(cluster_dir)) {
    cluster_files <- list.files(cluster_dir, pattern = "^04_", full.names = TRUE)
    if (length(cluster_files)) {
      dir.create(cluster_subdir, recursive = TRUE, showWarnings = FALSE)
      file.copy(cluster_files, cluster_subdir, overwrite = TRUE)
    }
  }

  invisible(out_dir)
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
  seed_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  set.seed(1L)
  keep_idx <- sort(sample.int(n_peaks, keep_n))
  if (!is.null(seed_state)) .Random.seed <<- seed_state

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

plot_tf_corr_stats_pdf <- function(
    ann_pearson,
    ann_spearman,
    out_dir,
    db,
    mode = "canonical",
    r_thr = NULL,
    p_thr = NULL,
    ann_pearson_all = NULL,
    ann_spearman_all = NULL,
    max_tfs = NULL,
    verbose = TRUE
) {
  if (!is.data.frame(ann_pearson) && !is.data.frame(ann_spearman)) {
    .log_abort("Both Pearson and Spearman annotation tables are missing.")
  }
  if (!is.character(out_dir) || !nzchar(out_dir)) {
    .log_abort("`out_dir` must be a non-empty path.")
  }
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tf_col <- function(dt) {
    if (!is.data.frame(dt)) return(NA_character_)
    if ("tfs" %in% names(dt)) return("tfs")
    if ("TF" %in% names(dt)) return("TF")
    NA_character_
  }
  r_col <- "corr_fp_tf_r"

  if (!is.data.frame(ann_pearson_all)) ann_pearson_all <- ann_pearson
  if (!is.data.frame(ann_spearman_all)) ann_spearman_all <- ann_spearman

  tf_col_p <- tf_col(ann_pearson)
  tf_col_s <- tf_col(ann_spearman)
  tf_col_p_all <- tf_col(ann_pearson_all)
  tf_col_s_all <- tf_col(ann_spearman_all)
  if (is.data.frame(ann_pearson) && is.na(tf_col_p)) .log_abort("Pearson annotation missing TF column (tfs/TF).")
  if (is.data.frame(ann_spearman) && is.na(tf_col_s)) .log_abort("Spearman annotation missing TF column (tfs/TF).")
  if (is.data.frame(ann_pearson_all) && is.na(tf_col_p_all)) .log_abort("Pearson (all) annotation missing TF column (tfs/TF).")
  if (is.data.frame(ann_spearman_all) && is.na(tf_col_s_all)) .log_abort("Spearman (all) annotation missing TF column (tfs/TF).")

  get_r_vals <- function(dt, tf_col_name, tf = NULL) {
    if (!is.data.frame(dt) || !r_col %in% names(dt)) return(numeric(0))
    if (!is.null(tf) && !is.na(tf_col_name)) {
      dt <- dt[dt[[tf_col_name]] == tf, , drop = FALSE]
    }
    v <- dt[[r_col]]
    v <- v[is.finite(v)]
    v
  }

  r_all_p <- get_r_vals(ann_pearson, tf_col_p)
  r_all_s <- get_r_vals(ann_spearman, tf_col_s)
  r_all_p_all <- get_r_vals(ann_pearson_all, tf_col_p_all)
  r_all_s_all <- get_r_vals(ann_spearman_all, tf_col_s_all)

  tfs_p <- if (is.data.frame(ann_pearson)) unique(ann_pearson[[tf_col_p]]) else character(0)
  tfs_s <- if (is.data.frame(ann_spearman)) unique(ann_spearman[[tf_col_s]]) else character(0)
  tfs_p_all <- if (is.data.frame(ann_pearson_all)) unique(ann_pearson_all[[tf_col_p_all]]) else character(0)
  tfs_s_all <- if (is.data.frame(ann_spearman_all)) unique(ann_spearman_all[[tf_col_s_all]]) else character(0)
  tfs <- sort(unique(c(tfs_p, tfs_s, tfs_p_all, tfs_s_all)))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  if (!is.null(max_tfs) && is.finite(max_tfs) && max_tfs > 0L) {
    tfs <- tfs[seq_len(min(length(tfs), as.integer(max_tfs)))]
  }

  col_p <- "#3182bd"
  col_s <- "#e6550d"
  col_p_all <- grDevices::adjustcolor("grey60", alpha.f = 0.35)
  col_s_all <- grDevices::adjustcolor("grey60", alpha.f = 0.35)

  panel_hist_overlay <- function(vals_all, vals_filt, title, col, col_all) {
    vals_all <- suppressWarnings(as.numeric(vals_all))
    vals_all <- vals_all[is.finite(vals_all)]
    vals_filt <- suppressWarnings(as.numeric(vals_filt))
    vals_filt <- vals_filt[is.finite(vals_filt)]
    if (!length(vals_all) && !length(vals_filt)) return(FALSE)
    vals_all <- pmin(pmax(vals_all, -1), 1)
    vals_filt <- pmin(pmax(vals_filt, -1), 1)

    hist_all <- NULL
    if (length(vals_all)) {
      hist_all <- graphics::hist(vals_all, breaks = seq(-1, 1, length.out = 41), plot = FALSE)
    }

    hist_filt <- NULL
    if (length(vals_filt)) {
      hist_filt <- graphics::hist(vals_filt, breaks = seq(-1, 1, length.out = 41), plot = FALSE)
    }

    ylim_vals <- c(
      if (!is.null(hist_all)) hist_all$counts else numeric(0),
      if (!is.null(hist_filt)) hist_filt$counts else numeric(0)
    )
    ylim_max <- max(c(ylim_vals, 0), na.rm = TRUE)
    if (!is.finite(ylim_max) || ylim_max <= 0) ylim_max <- 1

    graphics::plot.new()
    graphics::plot.window(xlim = c(-1, 1), ylim = c(0, ylim_max * 1.05))
    graphics::title(main = title, font.main = 2)
    graphics::axis(1)
    graphics::axis(2)
    graphics::box()
    graphics::mtext("R value", side = 1, line = 2.2, font = 2)
    graphics::mtext("Count", side = 2, line = 2.4, font = 2)

    if (!is.null(hist_all)) {
      graphics::rect(
        xleft = hist_all$breaks[-length(hist_all$breaks)],
        ybottom = 0,
        xright = hist_all$breaks[-1],
        ytop = hist_all$counts,
        col = col_all,
        border = NA
      )
    }
    if (!is.null(hist_filt)) {
      graphics::rect(
        xleft = hist_filt$breaks[-length(hist_filt$breaks)],
        ybottom = 0,
        xright = hist_filt$breaks[-1],
        ytop = hist_filt$counts,
        col = grDevices::adjustcolor(col, alpha.f = 0.5),
        border = NA
      )
    }
    TRUE
  }

  .mk_panel <- function(fun, vals_all, vals_filt, title, col, col_all) {
    local({
      vals_all <- vals_all
      vals_filt <- vals_filt
      title <- title
      col <- col
      col_all <- col_all
      function() fun(vals_all, vals_filt, title, col, col_all)
    })
  }

  panels <- list(
    .mk_panel(panel_hist_overlay, r_all_p_all, r_all_p, "Overall Pearson R", col_p, col_p_all),
    .mk_panel(panel_hist_overlay, r_all_s_all, r_all_s, "Overall Spearman R", col_s, col_s_all)
  )

  for (tf in tfs) {
    r_tf_p <- get_r_vals(ann_pearson, tf_col_p, tf = tf)
    r_tf_s <- get_r_vals(ann_spearman, tf_col_s, tf = tf)
    r_tf_p_all <- get_r_vals(ann_pearson_all, tf_col_p_all, tf = tf)
    r_tf_s_all <- get_r_vals(ann_spearman_all, tf_col_s_all, tf = tf)
    if (length(r_tf_p) || length(r_tf_p_all)) {
      panels <- c(panels, list(
        .mk_panel(panel_hist_overlay, r_tf_p_all, r_tf_p, paste0(tf, " Pearson R"), col_p, col_p_all)
      ))
    }
    if (length(r_tf_s) || length(r_tf_s_all)) {
      panels <- c(panels, list(
        .mk_panel(panel_hist_overlay, r_tf_s_all, r_tf_s, paste0(tf, " Spearman R"), col_s, col_s_all)
      ))
    }
  }

  pdf_path <- file.path(out_dir, sprintf("06_tf_corr_stats_%s%s.pdf", db, if (identical(mode, "all")) "" else paste0("_", mode)))
  grDevices::pdf(pdf_path, width = 11, height = 8.5)

  n_panels <- length(panels)
  idx <- 1L
  while (idx <= n_panels) {
    graphics::par(mfrow = c(3, 4), mar = c(4.2, 4.4, 2.5, 0.8), oma = c(2, 1, 2, 1), font.lab = 2, font.main = 2)
    for (k in 0:11) {
      if (idx + k <= n_panels) {
        panels[[idx + k]]()
      } else {
        graphics::plot.new()
      }
    }
    caption_parts <- character(0)
    if (!is.null(r_thr) && is.finite(r_thr)) {
      caption_parts <- c(caption_parts, paste0("r >", r_thr))
    }
    if (!is.null(p_thr) && is.finite(p_thr) && p_thr < 1) {
      caption_parts <- c(caption_parts, paste0("FDR <", p_thr))
    }
    if (length(caption_parts)) {
      graphics::mtext(
        paste0("Filtered: ", paste(caption_parts, collapse = ", ")),
        side = 1,
        outer = TRUE,
        line = 0.5,
        cex = 0.7
      )
    }
    graphics::mtext(paste0("TF correlation stats (", mode, ")"), side = 3, outer = TRUE, line = 0.2, cex = 0.9, font = 2)
    idx <- idx + 12L
  }

  grDevices::dev.off()
  if (isTRUE(verbose)) .log_inform("TF correlation stats PDF saved: {pdf_path}")
  invisible(pdf_path)
}

.step1_hist_counts <- function(x, breaks) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  x <- pmin(pmax(x, -1), 1)
  out <- integer(length(breaks) - 1L)
  if (!length(x)) return(out)
  h <- graphics::hist(x, breaks = breaks, plot = FALSE)
  as.integer(h$counts)
}

.step1_panel_hist_counts <- function(counts_all, counts_filt, breaks, title, col, col_all) {
  counts_all <- as.numeric(counts_all)
  counts_filt <- as.numeric(counts_filt)
  ylim_max <- max(c(counts_all, counts_filt, 0), na.rm = TRUE)
  if (!is.finite(ylim_max) || ylim_max <= 0) ylim_max <- 1

  graphics::plot.new()
  graphics::plot.window(xlim = c(-1, 1), ylim = c(0, ylim_max * 1.05))
  graphics::title(main = title, font.main = 2)
  graphics::axis(1)
  graphics::axis(2)
  graphics::box()
  graphics::mtext("R value", side = 1, line = 2.2, font = 2)
  graphics::mtext("Count", side = 2, line = 2.4, font = 2)
  graphics::rect(
    xleft = breaks[-length(breaks)],
    ybottom = 0,
    xright = breaks[-1],
    ytop = counts_all,
    col = col_all,
    border = NA
  )
  graphics::rect(
    xleft = breaks[-length(breaks)],
    ybottom = 0,
    xright = breaks[-1],
    ytop = counts_filt,
    col = grDevices::adjustcolor(col, alpha.f = 0.5),
    border = NA
  )
  invisible(TRUE)
}

plot_tf_corr_stats_pdf_from_cache <- function(
    cache_dir,
    out_dir,
    db,
    mode = "all",
    r_thr = NULL,
    p_thr = NULL,
    max_tfs = NULL,
    verbose = TRUE
) {
  if (!is.character(cache_dir) || !nzchar(cache_dir)) {
    .log_abort("`cache_dir` must be a non-empty path.")
  }
  if (!dir.exists(cache_dir)) {
    .log_abort("`cache_dir` does not exist: {cache_dir}")
  }
  if (!is.character(out_dir) || !nzchar(out_dir)) {
    .log_abort("`out_dir` must be a non-empty path.")
  }
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pear_dir <- file.path(cache_dir, "pearson_by_tf")
  spear_dir <- file.path(cache_dir, "spearman_by_tf")
  pear_files <- if (dir.exists(pear_dir)) sort(list.files(pear_dir, pattern = "\\.rds$", full.names = TRUE)) else character(0)
  spear_files <- if (dir.exists(spear_dir)) sort(list.files(spear_dir, pattern = "\\.rds$", full.names = TRUE)) else character(0)
  tf_names <- sort(unique(c(
    sub("\\.rds$", "", basename(pear_files)),
    sub("\\.rds$", "", basename(spear_files))
  )))
  if (!length(tf_names)) {
    .log_abort("No per-TF correlation cache files were found in {cache_dir}.")
  }
  if (!is.null(max_tfs) && is.finite(max_tfs) && max_tfs > 0L) {
    tf_names <- tf_names[seq_len(min(length(tf_names), as.integer(max_tfs)))]
  }

  breaks <- seq(-1, 1, length.out = 41)

  load_one <- function(dir_path, tf_name) {
    path <- file.path(dir_path, paste0(tf_name, ".rds"))
    if (!file.exists(path)) return(NULL)
    readRDS(path)
  }

  get_r <- function(tbl) {
    if (!is.data.frame(tbl) || !"corr_fp_tf_r" %in% names(tbl)) return(numeric(0))
    tbl$corr_fp_tf_r
  }

  col_p <- "#3182bd"
  col_s <- "#e6550d"
  col_p_all <- grDevices::adjustcolor("grey60", alpha.f = 0.35)
  col_s_all <- grDevices::adjustcolor("grey60", alpha.f = 0.35)

  pdf_path <- file.path(out_dir, sprintf("06_tf_corr_stats_%s%s.pdf", db, if (identical(mode, "all")) "" else paste0("_", mode)))
  grDevices::pdf(pdf_path, width = 11, height = 8.5)
  panels <- list()

  append_tf_panels <- function(tf_name) {
    pear_tbl <- load_one(pear_dir, tf_name)
    spear_tbl <- load_one(spear_dir, tf_name)
    pear_r <- get_r(pear_tbl)
    spear_r <- get_r(spear_tbl)
    out <- list()
    if (length(pear_r)) {
      out[[length(out) + 1L]] <- local({
        counts_all <- .step1_hist_counts(pear_r, breaks)
        counts_filt <- .step1_hist_counts(pear_r[is.finite(pear_r) & pear_r > r_thr], breaks)
        function() .step1_panel_hist_counts(counts_all, counts_filt, breaks, paste0(tf_name, " Pearson R"), col_p, col_p_all)
      })
    }
    if (length(spear_r)) {
      out[[length(out) + 1L]] <- local({
        counts_all <- .step1_hist_counts(spear_r, breaks)
        counts_filt <- .step1_hist_counts(spear_r[is.finite(spear_r) & spear_r > r_thr], breaks)
        function() .step1_panel_hist_counts(counts_all, counts_filt, breaks, paste0(tf_name, " Spearman R"), col_s, col_s_all)
      })
    }
    out
  }

  current_panels <- panels
  tf_idx <- 1L
  while (length(current_panels) || tf_idx <= length(tf_names)) {
    while (length(current_panels) < 12L && tf_idx <= length(tf_names)) {
      current_panels <- c(current_panels, append_tf_panels(tf_names[[tf_idx]]))
      tf_idx <- tf_idx + 1L
    }
    graphics::par(mfrow = c(3, 4), mar = c(4.2, 4.4, 2.5, 0.8), oma = c(2, 1, 2, 1), font.lab = 2, font.main = 2)
    for (k in seq_len(12L)) {
      if (k <= length(current_panels)) {
        current_panels[[k]]()
      } else {
        graphics::plot.new()
      }
    }
    caption_parts <- character(0)
    if (!is.null(r_thr) && is.finite(r_thr)) {
      caption_parts <- c(caption_parts, paste0("r >", r_thr))
    }
    if (!is.null(p_thr) && is.finite(p_thr) && p_thr < 1) {
      caption_parts <- c(caption_parts, paste0("FDR <", p_thr))
    }
    if (length(caption_parts)) {
      graphics::mtext(
        paste0("Filtered: ", paste(caption_parts, collapse = ", ")),
        side = 1,
        outer = TRUE,
        line = 0.5,
        cex = 0.7
      )
    }
    graphics::mtext(paste0("TF correlation stats (", mode, ")"), side = 3, outer = TRUE, line = 0.2, cex = 0.9, font = 2)
    if (length(current_panels) > 12L) {
      current_panels <- current_panels[-seq_len(12L)]
    } else {
      current_panels <- list()
    }
  }

  grDevices::dev.off()
  if (isTRUE(verbose)) .log_inform("TF correlation stats PDF saved: {pdf_path}")
  invisible(pdf_path)
}
