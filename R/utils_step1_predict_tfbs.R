# File: utils_step1_predict_tfbs.R
# Purpose: Public Module 1 TFBS prediction wrapper.

.module1_condition_columns <- function(omics_data) {
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.data.frame(omics_data$fp_score_condition_qn)) {
    .log_abort("`omics_data$fp_score_condition_qn` is required.")
  }
  if (!is.data.frame(omics_data$rna_condition)) {
    .log_abort("`omics_data$rna_condition` is required.")
  }
  intersect(
    setdiff(names(omics_data$fp_score_condition_qn), "peak_ID"),
    setdiff(names(omics_data$rna_condition), c("ensembl_gene_id", "HGNC"))
  )
}

.module1_prepare_predict_omics <- function(omics_data, label_col = NULL, verbose = TRUE) {
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a prepared Module 1 list.")
  if (!is.data.frame(omics_data$fp_annotation)) {
    .log_abort("`omics_data$fp_annotation` is required.")
  }
  if (!"tfs" %in% names(omics_data$fp_annotation)) {
    omics_data <- grn_add_fp_tfs(omics_data, force = TRUE, verbose = verbose)
  }
  if (!is.data.frame(omics_data$rna_condition)) {
    if (is.null(label_col) || !nzchar(label_col)) {
      .log_abort("`label_col` is required when `omics_data$rna_condition` is missing.")
    }
    omics_data <- grn_add_rna_condition(omics_data, label_col = label_col, verbose = verbose)
  }
  if (!is.data.frame(omics_data$fp_score_condition_qn)) {
    if (!is.data.frame(omics_data$fp_score_condition)) {
      if (is.null(label_col) || !nzchar(label_col)) {
        .log_abort("`label_col` is required when footprint condition matrices are missing.")
      }
      omics_data <- grn_add_fp_score_condition(omics_data, label_col = label_col, verbose = verbose)
    }
    omics_data <- grn_add_fp_score_qn(omics_data, id_col = "peak_ID", verbose = verbose)
  }
  omics_data
}

.module1_expressed_tfs <- function(omics_data, tf_subset = NULL) {
  rna_condition <- omics_data$rna_condition
  tfs <- unique(as.character(rna_condition$HGNC))
  tfs <- tfs[!is.na(tfs) & nzchar(tfs)]
  if (!is.null(omics_data$tf_list)) {
    tf_list <- unique(as.character(omics_data$tf_list))
    tf_list <- tf_list[!is.na(tf_list) & nzchar(tf_list)]
    if (length(tf_list)) {
      tfs <- intersect(tfs, tf_list)
    }
  }
  if (is.data.frame(omics_data$rna_expressed) && "HGNC" %in% names(omics_data$rna_expressed)) {
    cond_cols <- intersect(
      setdiff(names(omics_data$rna_expressed), c("ensembl_gene_id", "HGNC")),
      .module1_condition_columns(omics_data)
    )
    if (length(cond_cols)) {
      expr_mat <- as.matrix(omics_data$rna_expressed[, cond_cols, drop = FALSE])
      expressed <- rowSums(expr_mat > 0, na.rm = TRUE) > 0
      tfs <- intersect(tfs, as.character(omics_data$rna_expressed$HGNC[expressed]))
    }
  }
  if (!is.null(tf_subset)) {
    tfs <- intersect(tfs, unique(as.character(tf_subset)))
  }
  sort(tfs)
}

.module1_motif_supported_pairs <- function(omics_data, tf_subset = NULL) {
  ann <- omics_data$fp_annotation
  if (!is.data.frame(ann)) .log_abort("`omics_data$fp_annotation` is required.")
  if (!all(c("fp_peak", "atac_peak", "tfs") %in% names(ann))) {
    .log_abort("`omics_data$fp_annotation` must include fp_peak, atac_peak, and tfs.")
  }
  tf_subset <- if (is.null(tf_subset)) NULL else unique(as.character(tf_subset))
  tf_subset <- tf_subset[!is.na(tf_subset) & nzchar(tf_subset)]
  if (length(tf_subset)) {
    tfs_index <- paste0(",", gsub("[[:space:]]+", "", as.character(ann$tfs)), ",")
    keep <- rep(FALSE, length(tfs_index))
    for (tf in tf_subset) {
      keep <- keep | grepl(paste0(",", tf, ","), tfs_index, fixed = TRUE)
    }
    ann <- ann[keep, , drop = FALSE]
  }
  if (!nrow(ann)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  motif_values <- if ("motifs" %in% names(ann)) as.character(ann$motifs) else rep(NA_character_, nrow(ann))
  tf_parts <- strsplit(gsub("[[:space:]]+", "", as.character(ann$tfs)), ",", fixed = FALSE)
  tf_parts <- lapply(tf_parts, function(x) {
    x <- unique(x[!is.na(x) & nzchar(x)])
    if (!is.null(tf_subset)) x <- intersect(x, tf_subset)
    x
  })
  n_tf <- lengths(tf_parts)
  keep_rows <- n_tf > 0L
  if (!any(keep_rows)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  row_idx <- rep(which(keep_rows), n_tf[keep_rows])
  out <- data.table::data.table(
    fp_id = as.character(ann$fp_peak[row_idx]),
    atac_peak = as.character(ann$atac_peak[row_idx]),
    motifs = motif_values[row_idx],
    tf = unlist(tf_parts[keep_rows], use.names = FALSE)
  )
  if (!nrow(out)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  out <- out[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf)]
  tibble::as_tibble(unique(out))
}

.module1_all_prediction_pairs <- function(omics_data, fp_ids, tf_subset = NULL) {
  ann <- omics_data$fp_annotation
  if (!is.data.frame(ann)) .log_abort("`omics_data$fp_annotation` is required.")
  fp_ids <- unique(as.character(fp_ids))
  fp_ids <- fp_ids[!is.na(fp_ids) & nzchar(fp_ids)]
  tfs <- .module1_expressed_tfs(omics_data, tf_subset = tf_subset)
  if (!length(fp_ids) || !length(tfs)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  base <- ann[match(fp_ids, ann$fp_peak), c("fp_peak", "atac_peak", "motifs"), drop = FALSE]
  base <- base[!is.na(base$fp_peak), , drop = FALSE]
  base <- base[!duplicated(base$fp_peak), , drop = FALSE]
  pairs <- tidyr::crossing(
    fp_id = as.character(base$fp_peak),
    tf = tfs
  )
  dplyr::left_join(
    pairs,
    tibble::tibble(
      fp_id = as.character(base$fp_peak),
      atac_peak = as.character(base$atac_peak),
      motifs = as.character(base$motifs)
    ),
    by = "fp_id"
  )
}

.module1_compute_pair_correlations <- function(omics_data, pairs, min_non_na = 3L) {
  if (!is.data.frame(pairs) || !nrow(pairs)) {
    empty <- tibble::tibble(
      fp_id = character(),
      atac_peak = character(),
      motifs = character(),
      tf = character(),
      pearson_r = numeric(),
      pearson_p = numeric(),
      pearson_p_adj = numeric(),
      spearman_r = numeric(),
      spearman_p = numeric(),
      spearman_p_adj = numeric()
    )
    return(empty)
  }
  cond_cols <- .module1_condition_columns(omics_data)
  min_non_na <- as.integer(min_non_na)[[1L]]
  if (!length(cond_cols) || length(cond_cols) < min_non_na) {
    .log_abort("Not enough shared condition columns for Module 1 correlations.")
  }

  fp_tbl <- omics_data$fp_score_condition_qn
  rna_tbl <- omics_data$rna_condition
  fp_mat <- data.matrix(fp_tbl[, cond_cols, drop = FALSE])
  rownames(fp_mat) <- as.character(fp_tbl$peak_ID)
  rna_tbl <- rna_tbl[!duplicated(rna_tbl$HGNC), , drop = FALSE]
  rna_mat <- data.matrix(rna_tbl[, cond_cols, drop = FALSE])
  rownames(rna_mat) <- as.character(rna_tbl$HGNC)

  row_rank <- function(x, ok) {
    ranked <- x
    for (i in seq_len(nrow(ranked))) {
      values <- ranked[i, ]
      values[!ok[i, ]] <- NA_real_
      ranked[i, ] <- rank(values, na.last = "keep", ties.method = "average")
    }
    ranked
  }

  row_cor <- function(x, y, method) {
    ok <- is.finite(x) & is.finite(y)
    if (identical(method, "spearman")) {
      x <- row_rank(x, ok)
      y <- row_rank(y, ok)
      ok <- is.finite(x) & is.finite(y)
    }
    n <- rowSums(ok)
    x[!ok] <- NA_real_
    y[!ok] <- NA_real_
    x_mean <- rowSums(x, na.rm = TRUE) / pmax(n, 1L)
    y_mean <- rowSums(y, na.rm = TRUE) / pmax(n, 1L)
    x_centered <- sweep(x, 1L, x_mean, "-")
    y_centered <- sweep(y, 1L, y_mean, "-")
    x_centered[!ok] <- 0
    y_centered[!ok] <- 0
    numerator <- rowSums(x_centered * y_centered)
    denominator <- sqrt(rowSums(x_centered^2) * rowSums(y_centered^2))
    r <- numerator / denominator
    r[n < min_non_na | denominator == 0] <- NA_real_
    r <- pmax(pmin(r, 1), -1)
    p <- rep(NA_real_, length(r))
    p_ok <- is.finite(r) & n > 2L
    perfect <- p_ok & abs(r) >= 1
    p[perfect] <- 0
    ordinary <- p_ok & !perfect
    statistic <- r[ordinary] * sqrt((n[ordinary] - 2) / pmax(1 - r[ordinary]^2, .Machine$double.eps))
    p[ordinary] <- 2 * stats::pt(abs(statistic), df = n[ordinary] - 2, lower.tail = FALSE)
    list(r = r, p = p)
  }

  out <- tibble::as_tibble(pairs)
  out$pearson_r <- NA_real_
  out$pearson_p <- NA_real_
  out$spearman_r <- NA_real_
  out$spearman_p <- NA_real_
  fp_idx <- match(out$fp_id, rownames(fp_mat))
  tf_idx <- match(out$tf, rownames(rna_mat))
  valid <- !is.na(fp_idx) & !is.na(tf_idx)
  if (any(valid)) {
    x <- fp_mat[fp_idx[valid], , drop = FALSE]
    y <- rna_mat[tf_idx[valid], , drop = FALSE]
    pearson <- row_cor(x, y, "pearson")
    spearman <- row_cor(x, y, "spearman")
    out$pearson_r[valid] <- pearson$r
    out$pearson_p[valid] <- pearson$p
    out$spearman_r[valid] <- spearman$r
    out$spearman_p[valid] <- spearman$p
  }
  out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
  out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
  out
}

.module1_rank_matrix_rows <- function(x) {
  ranked <- x
  for (i in seq_len(nrow(ranked))) {
    ranked[i, ] <- rank(ranked[i, ], na.last = "keep", ties.method = "average")
  }
  ranked
}

.module1_cor_matrix_vector <- function(x, y, min_non_na = 3L) {
  y <- as.numeric(y)
  ok_y <- is.finite(y)
  ok <- is.finite(x)
  ok[, !ok_y] <- FALSE
  n <- rowSums(ok)

  x_work <- x
  x_work[!ok] <- NA_real_
  y_work <- matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  y_work[!ok] <- NA_real_

  x_mean <- rowSums(x_work, na.rm = TRUE) / pmax(n, 1L)
  y_mean <- rowSums(y_work, na.rm = TRUE) / pmax(n, 1L)
  x_centered <- sweep(x_work, 1L, x_mean, "-")
  y_centered <- sweep(y_work, 1L, y_mean, "-")
  x_centered[!ok] <- 0
  y_centered[!ok] <- 0

  numerator <- rowSums(x_centered * y_centered)
  denominator <- sqrt(rowSums(x_centered^2) * rowSums(y_centered^2))
  r <- numerator / denominator
  r[n < as.integer(min_non_na)[[1L]] | denominator == 0] <- NA_real_
  r <- pmax(pmin(r, 1), -1)

  p <- rep(NA_real_, length(r))
  p_ok <- is.finite(r) & n > 2L
  perfect <- p_ok & abs(r) >= 1
  p[perfect] <- 0
  ordinary <- p_ok & !perfect
  statistic <- r[ordinary] * sqrt((n[ordinary] - 2) / pmax(1 - r[ordinary]^2, .Machine$double.eps))
  p[ordinary] <- 2 * stats::pt(abs(statistic), df = n[ordinary] - 2, lower.tail = FALSE)

  list(r = r, p = p)
}

.module1_predict_links_streamed <- function(omics_data,
                                            high_confidence_footprints,
                                            r_cutoff = 0.3,
                                            tf_subset = NULL,
                                            min_non_na = 3L) {
  if (!is.data.frame(high_confidence_footprints)) {
    .log_abort("`high_confidence_footprints` must be a data.frame.")
  }
  need_high <- c("fp_id", "chr", "start", "end", "atac_peak")
  if (!all(need_high %in% names(high_confidence_footprints))) {
    .log_abort("`high_confidence_footprints` is missing required columns.")
  }

  cond_cols <- .module1_condition_columns(omics_data)
  if (!length(cond_cols) || length(cond_cols) < as.integer(min_non_na)[[1L]]) {
    .log_abort("Not enough shared condition columns for Module 1 correlations.")
  }

  high <- high_confidence_footprints[!duplicated(high_confidence_footprints$fp_id), need_high, drop = FALSE]
  tfs <- .module1_expressed_tfs(omics_data, tf_subset = tf_subset)
  pair_count <- as.double(nrow(high)) * as.double(length(tfs))
  if (!nrow(high) || !length(tfs)) {
    empty_links <- .module1_build_tfbs_links(
      tfbs_stats = tibble::tibble(
        fp_id = character(),
        atac_peak = character(),
        tf = character(),
        best_r = numeric(),
        best_method = character(),
        pass = logical()
      ),
      high_confidence_footprints = high,
      omics_data = omics_data
    )
    return(list(tfbs_links = empty_links, prediction_stats = empty_links[0, ], prediction_pair_count = pair_count))
  }

  fp_tbl <- omics_data$fp_score_condition_qn
  rna_tbl <- omics_data$rna_condition
  fp_idx <- match(high$fp_id, fp_tbl$peak_ID)
  keep_fp <- !is.na(fp_idx)
  high <- high[keep_fp, , drop = FALSE]
  fp_idx <- fp_idx[keep_fp]
  fp_mat <- data.matrix(fp_tbl[fp_idx, cond_cols, drop = FALSE])

  rna_tbl <- rna_tbl[!duplicated(rna_tbl$HGNC), , drop = FALSE]
  tf_idx <- match(tfs, rna_tbl$HGNC)
  keep_tf <- !is.na(tf_idx)
  tfs <- tfs[keep_tf]
  tf_idx <- tf_idx[keep_tf]
  rna_mat <- data.matrix(rna_tbl[tf_idx, cond_cols, drop = FALSE])
  rownames(rna_mat) <- tfs

  if (!nrow(high) || !length(tfs)) {
    empty_links <- .module1_build_tfbs_links(
      tfbs_stats = tibble::tibble(
        fp_id = character(),
        atac_peak = character(),
        tf = character(),
        best_r = numeric(),
        best_method = character(),
        pass = logical()
      ),
      high_confidence_footprints = high,
      omics_data = omics_data
    )
    return(list(tfbs_links = empty_links, prediction_stats = empty_links[0, ], prediction_pair_count = pair_count))
  }

  fp_rank <- .module1_rank_matrix_rows(fp_mat)
  rna_rank <- .module1_rank_matrix_rows(rna_mat)
  stats_list <- vector("list", length(tfs))
  r_cutoff <- as.numeric(r_cutoff)[[1L]]

  for (i in seq_along(tfs)) {
    pearson <- .module1_cor_matrix_vector(fp_mat, rna_mat[i, ], min_non_na = min_non_na)
    spearman <- .module1_cor_matrix_vector(fp_rank, rna_rank[i, ], min_non_na = min_non_na)
    best <- .module1_best_corr(pearson$r, spearman$r)
    pass <- is.finite(best$best_r) & best$best_r >= r_cutoff
    if (!any(pass)) next
    stats_list[[i]] <- tibble::tibble(
      fp_id = high$fp_id[pass],
      atac_peak = high$atac_peak[pass],
      tf = tfs[[i]],
      pearson_r = pearson$r[pass],
      pearson_p = pearson$p[pass],
      pearson_p_adj = stats::p.adjust(pearson$p[pass], method = "BH"),
      spearman_r = spearman$r[pass],
      spearman_p = spearman$p[pass],
      spearman_p_adj = stats::p.adjust(spearman$p[pass], method = "BH"),
      best_r = best$best_r[pass],
      best_method = best$best_method[pass],
      pass = TRUE
    )
  }

  prediction_stats <- dplyr::bind_rows(stats_list)
  if (!nrow(prediction_stats)) {
    prediction_stats <- tibble::tibble(
      fp_id = character(),
      atac_peak = character(),
      tf = character(),
      pearson_r = numeric(),
      pearson_p = numeric(),
      pearson_p_adj = numeric(),
      spearman_r = numeric(),
      spearman_p = numeric(),
      spearman_p_adj = numeric(),
      best_r = numeric(),
      best_method = character(),
      pass = logical()
    )
  }
  tfbs_links <- .module1_build_tfbs_links(
    tfbs_stats = prediction_stats,
    high_confidence_footprints = high,
    omics_data = omics_data
  )

  list(
    tfbs_links = tfbs_links,
    prediction_stats = prediction_stats,
    prediction_pair_count = pair_count
  )
}

#' Predict transcription factor binding sites from matched footprint and RNA data
#'
#' Run the Module 1 TFBS workflow as one user-facing operation. The function
#' first uses motif-supported FP-TF correlations to define high-confidence
#' footprints, then predicts sparse FP-TF binding links for expressed TFs.
#'
#' @param omics_data Prepared Module 1 multiomic object.
#' @param out_dir Output directory.
#' @param db Motif database label used in output metadata.
#' @param label_col Metadata column used to build condition-level matrices when
#'   missing from `omics_data`.
#' @param r_cutoff Minimum positive correlation used for motif-supported and
#'   prediction calls.
#' @param tf_subset Optional TF subset.
#' @param write_outputs Write compact Module 1 output files.
#' @param write_stats Retain and write full FP-TF correlation statistics.
#' @param write_bed Write optional BED-like browser files. Not used by default.
#' @param min_non_na Minimum finite condition pairs required for correlation.
#' @param verbose Emit concise progress messages.
#'
#' @return A list containing `omics_data`, `high_confidence_footprints`,
#'   `motif_supported_correlations`, `tfbs_links`, `tfbs_stats`, `reports`, and
#'   `parameters`.
#' @export
predict_tfbs <- function(omics_data,
                         out_dir = "predict_tf_binding_sites",
                         db = "JASPAR2024",
                         label_col = NULL,
                         r_cutoff = 0.3,
                         tf_subset = NULL,
                         write_outputs = TRUE,
                         write_stats = FALSE,
                         write_bed = FALSE,
                         min_non_na = 3L,
                         verbose = TRUE) {
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a prepared Module 1 list.")
  if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    .log_abort("`out_dir` must be a non-empty path.")
  }
  if (!is.character(db) || length(db) != 1L || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  r_cutoff <- as.numeric(r_cutoff)[[1L]]
  if (!is.finite(r_cutoff)) .log_abort("`r_cutoff` must be finite.")
  omics_data <- .module1_prepare_predict_omics(
    omics_data = omics_data,
    label_col = label_col,
    verbose = verbose
  )

  if (isTRUE(verbose)) .log_inform("Module 1 TFBS prediction: motif-supported FP-TF correlations.")
  motif_supported_pairs <- .module1_motif_supported_pairs(omics_data, tf_subset = tf_subset)
  if (isTRUE(verbose)) {
    .log_inform("Module 1 TFBS prediction: {nrow(motif_supported_pairs)} motif-supported FP-TF pair(s).")
  }
  motif_supported_correlations_raw <- .module1_compute_pair_correlations(
    omics_data = omics_data,
    pairs = motif_supported_pairs,
    min_non_na = min_non_na
  )
  motif_supported_correlations <- .module1_merge_tfbs_stats(
    pearson_stats = motif_supported_correlations_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
    spearman_stats = motif_supported_correlations_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
    r_cutoff = r_cutoff
  )
  high_confidence_footprints <- .module1_select_high_confidence_footprints(
    motif_supported_correlations = motif_supported_correlations,
    r_cutoff = r_cutoff
  )
  if (isTRUE(verbose)) {
    .log_inform("Module 1 TFBS prediction: {nrow(high_confidence_footprints)} high-confidence footprint(s).")
  }

  if (isTRUE(verbose)) .log_inform("Module 1 TFBS prediction: prediction correlations.")
  if (isTRUE(write_stats)) {
    prediction_pairs <- .module1_all_prediction_pairs(
      omics_data = omics_data,
      fp_ids = high_confidence_footprints$fp_id,
      tf_subset = tf_subset
    )
    if (isTRUE(verbose)) {
      .log_inform("Module 1 TFBS prediction: {nrow(prediction_pairs)} prediction FP-TF pair(s).")
    }
    prediction_stats_raw <- .module1_compute_pair_correlations(
      omics_data = omics_data,
      pairs = prediction_pairs,
      min_non_na = min_non_na
    )
    prediction_stats <- .module1_merge_tfbs_stats(
      pearson_stats = prediction_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
      spearman_stats = prediction_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
      r_cutoff = r_cutoff
    )
    tfbs_links <- .module1_build_tfbs_links(
      tfbs_stats = prediction_stats,
      high_confidence_footprints = high_confidence_footprints,
      omics_data = omics_data
    )
  } else {
    streamed <- .module1_predict_links_streamed(
      omics_data = omics_data,
      high_confidence_footprints = high_confidence_footprints,
      r_cutoff = r_cutoff,
      tf_subset = tf_subset,
      min_non_na = min_non_na
    )
    prediction_stats <- streamed$prediction_stats
    tfbs_links <- streamed$tfbs_links
    if (isTRUE(verbose)) {
      .log_inform("Module 1 TFBS prediction: {streamed$prediction_pair_count} prediction FP-TF pair(s) evaluated without materializing all pairs.")
    }
  }

  reports <- list()
  if (isTRUE(write_outputs)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    high_path <- file.path(out_dir, "module1_high_confidence_footprints.csv")
    links_path <- file.path(out_dir, "module1_tfbs_links.csv.gz")
    params_path <- file.path(out_dir, "module1_parameters.yml")
    readr::write_csv(high_confidence_footprints, high_path)
    readr::write_csv(tfbs_links, links_path)
    yaml::write_yaml(
      list(db = db, r_cutoff = r_cutoff, min_non_na = as.integer(min_non_na), write_stats = isTRUE(write_stats)),
      params_path
    )
    reports <- c(reports, list(
      high_confidence_footprints = high_path,
      tfbs_links = links_path,
      parameters = params_path
    ))
    if (isTRUE(write_stats)) {
      stats_path <- file.path(out_dir, "module1_tfbs_stats.csv.gz")
      readr::write_csv(prediction_stats, stats_path)
      reports$tfbs_stats <- stats_path
    }
  }
  if (isTRUE(write_bed) && isTRUE(verbose)) {
    .log_warn("Optional Module 1 BED output is not implemented in `predict_tfbs()` yet.")
  }

  list(
    omics_data = omics_data,
    high_confidence_footprints = high_confidence_footprints,
    motif_supported_correlations = motif_supported_correlations,
    tfbs_links = tfbs_links,
    tfbs_stats = if (isTRUE(write_stats)) prediction_stats else NULL,
    reports = reports,
    parameters = list(
      db = db,
      label_col = label_col,
      r_cutoff = r_cutoff,
      min_non_na = as.integer(min_non_na),
      tf_subset = tf_subset,
      write_outputs = isTRUE(write_outputs),
      write_stats = isTRUE(write_stats),
      write_bed = isTRUE(write_bed)
    )
  )
}
