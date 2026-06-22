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

.module1_motif_gene_col <- function(motif_db) {
  if (!is.data.frame(motif_db)) return(NULL)
  if ("gene_symbol" %in% names(motif_db)) return("gene_symbol")
  if ("HGNC" %in% names(motif_db)) return("HGNC")
  NULL
}

.module1_add_fp_tfs <- function(grn_set, force = FALSE, verbose = TRUE) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_tfs") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_annotation)) {
    .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  }

  annot <- grn_set$fp_annotation
  if (!"fp_peak" %in% names(annot)) .log_abort("`fp_annotation` needs fp_peak.")
  if (!"motifs" %in% names(annot)) .log_abort("`fp_annotation` needs motifs.")

  gene_col <- .module1_motif_gene_col(grn_set$motif_db)
  has_db <- is.data.frame(grn_set$motif_db) &&
    !is.null(gene_col) &&
    all(c("motif", gene_col) %in% names(grn_set$motif_db))

  if (isTRUE(has_db)) {
    motif_join_tbl <- grn_set$motif_db |>
      dplyr::select("motif", gene_symbol = dplyr::all_of(gene_col))
    if ("motifs" %in% names(annot)) {
      annot$motifs <- sub("^_+", "", annot$motifs)
    }
    if ("motif" %in% names(motif_join_tbl)) {
      motif_join_tbl$motif <- sub("^_+", "", motif_join_tbl$motif)
    }
    annot <- annot |>
      dplyr::left_join(
        motif_join_tbl,
        by = c("motifs" = "motif")
      )
    annot$tfs <- ifelse(
      !is.na(annot$gene_symbol) & nzchar(annot$gene_symbol),
      gsub("\\s*::\\s*", ",", annot$gene_symbol),
      toupper(sub("\\..*$", "", sub("^_+", "", annot$motifs)))
    )
    annot$gene_symbol <- NULL
  } else {
    annot$tfs <- toupper(sub("\\..*$", "", sub("^_+", "", annot$motifs)))
  }

  grn_set$fp_annotation <- annot
  fp_tfs <- annot |>
    dplyr::select("fp_peak", "atac_peak", "tfs") |>
    tidyr::separate_rows("tfs", sep = "\\s*,\\s*|\\s*::\\s*") |>
    dplyr::filter(!is.na(.data$tfs), .data$tfs != "") |>
    dplyr::distinct(.data$fp_peak, .data$atac_peak, .data$tfs) |>
    dplyr::group_by(.data$fp_peak, .data$atac_peak) |>
    dplyr::summarise(tfs = paste(unique(.data$tfs), collapse = ","), .groups = "drop")

  grn_set$fp_tfs <- fp_tfs
  grn_status_set(grn_set, "fp_tfs")
}

.module1_prepare_predict_omics <- function(omics_data, label_col = NULL, verbose = TRUE) {
  if (!is_multiomic_object(omics_data)) {
    .log_abort("`omics_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  }
  omics_data <- .as_module1_analysis_data(omics_data, verbose = verbose)
  if (!is.data.frame(omics_data$fp_annotation)) {
    .log_abort("`omics_data$fp_annotation` is required.")
  }
  if (!"tfs" %in% names(omics_data$fp_annotation)) {
    omics_data <- .module1_add_fp_tfs(omics_data, force = TRUE, verbose = verbose)
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

  keep_cols <- intersect(c("fp_peak", "atac_peak", "motifs", "tfs"), names(ann))
  dt <- data.table::as.data.table(ann[, keep_cols, drop = FALSE])
  if (!"motifs" %in% names(dt)) dt[, motifs := NA_character_]
  dt[, tfs_clean := gsub("[[:space:]]+", "", as.character(tfs))]
  dt <- dt[!is.na(tfs_clean) & nzchar(tfs_clean)]
  if (!nrow(dt)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  dt[, tf := strsplit(tfs_clean, ",", fixed = TRUE)]
  out <- dt[, .(tf = unlist(tf, use.names = FALSE)), by = .(
    fp_id = fp_peak,
    atac_peak = atac_peak,
    motifs = motifs
  )]
  if (length(tf_subset)) out <- out[tf %chin% tf_subset]

  out <- unique(out[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf)])
  if (!nrow(out)) {
    return(tibble::tibble(fp_id = character(), atac_peak = character(), motifs = character(), tf = character()))
  }
  tibble::as_tibble(out)
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

.module1_compute_pair_correlations_cpp <- function(out,
                                           fp_mat,
                                           rna_mat,
                                           fp_idx,
                                           tf_idx,
                                           min_non_na = 3L,
                                           cores = NULL) {
  cpp_fun <- get0(".sparse_pair_correlations_cpp", envir = asNamespace("craftgrn"), mode = "function")
  if (!is.function(cpp_fun)) return(NULL)
  valid <- !is.na(fp_idx) & !is.na(tf_idx)
  if (!any(valid)) return(NULL)
  min_non_na <- as.integer(min_non_na)[[1L]]
  if (ncol(fp_mat) < min_non_na || ncol(fp_mat) != ncol(rna_mat)) return(NULL)
  if (any(!is.finite(fp_mat)) || any(!is.finite(rna_mat))) return(NULL)

  unique_fp_idx <- unique(fp_idx[valid])
  unique_tf_idx <- unique(tf_idx[valid])
  fp_sub <- fp_mat[unique_fp_idx, , drop = FALSE]
  rna_sub <- rna_mat[unique_tf_idx, , drop = FALSE]
  fp_rank <- .module1_rank_matrix_rows(fp_sub)
  rna_rank <- .module1_rank_matrix_rows(rna_sub)
  pair_fp_idx <- match(fp_idx[valid], unique_fp_idx) - 1L
  pair_tf_idx <- match(tf_idx[valid], unique_tf_idx) - 1L

  cores <- .module1_default_cores(cores)
  stats <- tibble::as_tibble(cpp_fun(
    fp = fp_sub,
    tf = rna_sub,
    fp_rank = fp_rank,
    tf_rank = rna_rank,
    fp_index = as.integer(pair_fp_idx),
    tf_index = as.integer(pair_tf_idx),
    n_threads = cores
  ))
  out$pearson_r[valid] <- stats$pearson_r
  out$pearson_p[valid] <- stats$pearson_p
  out$spearman_r[valid] <- stats$spearman_r
  out$spearman_p[valid] <- stats$spearman_p
  out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
  out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
  out
}

.module1_compute_pair_correlations <- function(omics_data, pairs, min_non_na = 3L, cores = NULL) {
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

  out <- tibble::as_tibble(pairs)
  out$pearson_r <- NA_real_
  out$pearson_p <- NA_real_
  out$spearman_r <- NA_real_
  out$spearman_p <- NA_real_
  fp_idx <- match(out$fp_id, rownames(fp_mat))
  tf_idx <- match(out$tf, rownames(rna_mat))
  valid <- !is.na(fp_idx) & !is.na(tf_idx)
  cpp_out <- tryCatch(
    .module1_compute_pair_correlations_cpp(
      out = out,
      fp_mat = fp_mat,
      rna_mat = rna_mat,
      fp_idx = fp_idx,
      tf_idx = tf_idx,
      min_non_na = min_non_na,
      cores = cores
    ),
    error = function(e) NULL
  )
  if (!is.null(cpp_out)) return(cpp_out)

  if (any(valid)) {
    unique_fp_idx <- unique(fp_idx[valid])
    unique_tf_idx <- unique(tf_idx[valid])
    fp_sub <- fp_mat[unique_fp_idx, , drop = FALSE]
    rna_sub <- rna_mat[unique_tf_idx, , drop = FALSE]
    pearson_all <- .module1_cor_matrix_matrix(fp_sub, rna_sub, min_non_na = min_non_na)
    if (!is.null(pearson_all)) {
      fp_rank <- .module1_rank_matrix_rows(fp_sub)
      rna_rank <- .module1_rank_matrix_rows(rna_sub)
      spearman_all <- .module1_cor_matrix_matrix(fp_rank, rna_rank, min_non_na = min_non_na)
      if (!is.null(spearman_all)) {
        pair_idx <- cbind(match(fp_idx[valid], unique_fp_idx), match(tf_idx[valid], unique_tf_idx))
        out$pearson_r[valid] <- pearson_all$r[pair_idx]
        out$pearson_p[valid] <- pearson_all$p[pair_idx]
        out$spearman_r[valid] <- spearman_all$r[pair_idx]
        out$spearman_p[valid] <- spearman_all$p[pair_idx]
        out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
        out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
        return(out)
      }
    }
  }

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

.module1_cor_matrix_matrix <- function(x, y, min_non_na = 3L) {
  min_non_na <- as.integer(min_non_na)[[1L]]
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!nrow(x) || !nrow(y) || !ncol(x) || ncol(x) != ncol(y)) {
    return(NULL)
  }
  if (ncol(x) < min_non_na || any(!is.finite(x)) || any(!is.finite(y))) {
    return(NULL)
  }

  x_centered <- sweep(x, 1L, rowMeans(x), "-")
  y_centered <- sweep(y, 1L, rowMeans(y), "-")
  x_ss <- rowSums(x_centered^2)
  y_ss <- rowSums(y_centered^2)
  denominator <- sqrt(outer(x_ss, y_ss, "*"))
  r <- tcrossprod(x_centered, y_centered) / denominator
  r[denominator == 0] <- NA_real_
  r <- pmax(pmin(r, 1), -1)

  p <- matrix(NA_real_, nrow = nrow(r), ncol = ncol(r))
  p_ok <- is.finite(r) & ncol(x) > 2L
  perfect <- p_ok & abs(r) >= 1
  p[perfect] <- 0
  ordinary <- p_ok & !perfect
  statistic <- r[ordinary] * sqrt((ncol(x) - 2) / pmax(1 - r[ordinary]^2, .Machine$double.eps))
  p[ordinary] <- 2 * stats::pt(abs(statistic), df = ncol(x) - 2, lower.tail = FALSE)

  list(r = r, p = p)
}

.module1_default_cores <- function(cores = NULL) {
  if (!is.null(cores)) {
    cores <- as.integer(cores)[[1L]]
    if (!is.finite(cores) || cores < 1L) .log_abort("`cores` must be a positive integer.")
    return(cores)
  }
  cores <- if (requireNamespace("parallelly", quietly = TRUE)) {
    parallelly::availableCores()
  } else {
    parallel::detectCores(logical = TRUE)
  }
  cores <- suppressWarnings(as.integer(cores)[[1L]])
  if (!is.finite(cores) || cores < 1L) 1L else cores
}

.module1_output_format <- function(output_format = c("csv", "parquet", "auto")) {
  output_format <- match.arg(output_format)
  if (identical(output_format, "auto")) {
    if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "csv"
  } else {
    output_format
  }
}

.module1_write_prediction_stats_chunk <- function(x, chunk_id, out_dir, output_format = c("csv", "parquet", "auto")) {
  output_format <- .module1_output_format(output_format)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- sprintf("module1_prediction_stats_chunk_%04d", as.integer(chunk_id))
  if (identical(output_format, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
    path <- file.path(out_dir, paste0(stem, ".parquet"))
    .write_parquet_table(x, path)
  } else {
    output_format <- "csv"
    path <- file.path(out_dir, paste0(stem, ".csv.gz"))
    readr::write_csv(x, path)
  }
  tibble::tibble(
    chunk_id = as.integer(chunk_id),
    path = path,
    format = output_format,
    n_rows = nrow(x)
  )
}

.module1_write_prediction_stats_manifest <- function(manifest, out_dir) {
  path <- file.path(out_dir, "module1_prediction_stats_manifest.csv")
  readr::write_csv(manifest, path)
  path
}

.module1_write_canonical_support_cache <- function(high_confidence_footprints, qc_summary, out_dir) {
  cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  keep <- intersect(
    c(
      "fp_id", "chr", "start", "end", "atac_peak",
      "n_canonical_bound_tfs", "canonical_bound_tfs",
      "canonical_bound_motifs", "has_canonical_bound"
    ),
    names(high_confidence_footprints)
  )
  canonical_fps <- if (is.data.frame(high_confidence_footprints) && nrow(high_confidence_footprints) && length(keep)) {
    high_confidence_footprints[, keep, drop = FALSE]
  } else {
    tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      n_canonical_bound_tfs = integer(),
      canonical_bound_tfs = character(),
      canonical_bound_motifs = character(),
      has_canonical_bound = logical()
    )
  }
  if (!"has_canonical_bound" %in% names(canonical_fps)) canonical_fps$has_canonical_bound <- TRUE
  canonical_fps <- canonical_fps[!duplicated(as.character(canonical_fps$fp_id)), , drop = FALSE]
  fps_path <- file.path(cache_dir, "module1_canonical_bound_fps.csv.gz")
  summary_path <- file.path(cache_dir, "module1_canonical_support_summary.csv")
  readr::write_csv(canonical_fps, fps_path)
  summary_tbl <- tibble::tibble(
    metric = c(
      "n_canonical_bound_fps",
      "n_canonical_pairs_pass",
      "n_motif_supported_pairs",
      "n_prediction_fps",
      "cache_schema"
    ),
    value = c(
      as.character(nrow(canonical_fps)),
      as.character(qc_summary[["n_canonical_pairs_pass"]] %||% NA_real_),
      as.character(qc_summary[["n_motif_supported_pairs"]] %||% NA_real_),
      as.character(qc_summary[["n_prediction_fps"]] %||% NA_real_),
      "module1_canonical_support_cache_v1"
    )
  )
  readr::write_csv(summary_tbl, summary_path)
  list(
    canonical_bound_fps_cache = fps_path,
    canonical_support_summary_cache = summary_path
  )
}


.module1_write_bed_outputs <- function(high_confidence_footprints, prediction_stats, out_dir) {
  if (!is.data.frame(high_confidence_footprints)) .log_abort("high_confidence_footprints must be a data.frame.")
  need_high <- c("chr", "start", "end", "fp_id")
  if (!all(need_high %in% names(high_confidence_footprints))) .log_abort("high_confidence_footprints is missing BED columns.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  high_keep <- intersect(c("chr", "start", "end", "fp_id", "n_canonical_bound_tfs"), names(high_confidence_footprints))
  high_bed <- high_confidence_footprints[, high_keep, drop = FALSE]
  names(high_bed)[seq_len(min(4L, ncol(high_bed)))] <- c("chrom", "chromStart", "chromEnd", "name")[seq_len(min(4L, ncol(high_bed)))]
  high_path <- file.path(out_dir, "module1_high_confidence_footprints.bed")
  data.table::fwrite(data.table::as.data.table(high_bed), high_path, sep = "\t", col.names = FALSE)
  reports <- list(high_confidence_footprints_bed = high_path)
  need_links <- c("chr", "start", "end", "tf", "fp_id", "best_r")
  if (is.data.frame(prediction_stats) && nrow(prediction_stats) && all(need_links %in% names(prediction_stats))) {
    score <- pmax(0L, pmin(1000L, as.integer(round(as.numeric(prediction_stats$best_r) * 1000))))
    bed <- data.frame(chrom = as.character(prediction_stats$chr), chromStart = as.integer(prediction_stats$start), chromEnd = as.integer(prediction_stats$end), name = paste(as.character(prediction_stats$tf), as.character(prediction_stats$fp_id), sep = "|"), score = score, strand = ".", stringsAsFactors = FALSE)
    links_path <- file.path(out_dir, "module1_prediction_stats.bed")
    data.table::fwrite(data.table::as.data.table(bed), links_path, sep = "\t", col.names = FALSE)
    reports$prediction_stats_bed <- links_path
  }
  reports
}
.module1_compute_prediction_stats_by_tf <- function(high,
                                                     tfs,
                                                     fp_mat,
                                                     rna_mat,
                                                     fp_rank,
                                                     rna_rank,
                                                     cutoffs,
                                                     min_non_na = 3L,
                                                     cores = 1L) {
  worker <- function(i) {
    pearson <- .module1_cor_matrix_vector(fp_mat, rna_mat[i, ], min_non_na = min_non_na)
    spearman <- .module1_cor_matrix_vector(fp_rank, rna_rank[i, ], min_non_na = min_non_na)
    best <- .module1_best_corr(pearson$r, spearman$r)
    stats_i <- tibble::tibble(
      fp_id = high$fp_id,
      atac_peak = high$atac_peak,
      tf = tfs[[i]],
      pearson_r = pearson$r,
      pearson_p = pearson$p,
      pearson_p_adj = stats::p.adjust(pearson$p, method = "BH"),
      spearman_r = spearman$r,
      spearman_p = spearman$p,
      spearman_p_adj = stats::p.adjust(spearman$p, method = "BH"),
      best_r = best$best_r,
      best_method = best$best_method
    )
    stats_i <- .module1_apply_tfbs_cutoffs(stats_i, cutoffs)
    stats_i[stats_i$pass %in% TRUE, , drop = FALSE]
  }
  idx <- seq_along(tfs)
  if (.Platform$OS.type != "windows" && cores > 1L && length(idx) > 1L) {
    dplyr::bind_rows(parallel::mclapply(idx, worker, mc.cores = min(cores, length(idx))))
  } else {
    dplyr::bind_rows(lapply(idx, worker))
  }
}

.module1_compute_prediction_stats_cpp <- function(high,
                                                   tfs,
                                                   fp_mat,
                                                   rna_mat,
                                                   fp_rank,
                                                   rna_rank,
                                                   cutoffs,
                                                   cores = 1L,
                                                   tf_chunk_size = NULL) {
  cpp_fun <- get0(".dense_prediction_stats_cpp", envir = asNamespace("craftgrn"), mode = "function")
  if (!is.function(cpp_fun)) return(NULL)
  if (!nrow(fp_mat) || !nrow(rna_mat) || ncol(fp_mat) != ncol(rna_mat)) return(NULL)
  if (any(!is.finite(fp_mat)) || any(!is.finite(rna_mat)) ||
      any(!is.finite(fp_rank)) || any(!is.finite(rna_rank))) {
    return(NULL)
  }

  if (is.null(tf_chunk_size)) {
    tf_chunk_size <- getOption("craftgrn.module1_cpp_tf_chunk_size", 50L)
  }
  tf_chunk_size <- suppressWarnings(as.integer(tf_chunk_size)[[1L]])
  if (!is.finite(tf_chunk_size) || tf_chunk_size < 1L) tf_chunk_size <- 50L
  idx <- seq_along(tfs)
  chunks <- split(idx, ceiling(idx / tf_chunk_size))
  out <- lapply(chunks, function(ii) {
    tibble::as_tibble(cpp_fun(
      fp = fp_mat,
      tf = rna_mat[ii, , drop = FALSE],
      fp_rank = fp_rank,
      tf_rank = rna_rank[ii, , drop = FALSE],
      fp_id = as.character(high$fp_id),
      atac_peak = as.character(high$atac_peak),
      tf_name = as.character(tfs[ii]),
      r_cutoff = cutoffs$r,
      p_cutoff = cutoffs$p,
      fdr_cutoff = cutoffs$fdr,
      n_threads = min(as.integer(cores)[[1L]], length(ii))
    ))
  })
  dplyr::bind_rows(out)
}

.module1_predict_tfbs_streamed <- function(omics_data,
                                            high_confidence_footprints,
                                            r_cutoff = 0.3,
                                            p_cutoff = NULL,
                                            fdr_cutoff = NULL,
                                            tf_subset = NULL,
                                            min_non_na = 3L,
                                            cores = NULL,
                                            return_prediction_stats = TRUE,
                                            prediction_stats_out_dir = NULL,
                                            output_format = c("csv", "parquet", "auto"),
                                            verbose = TRUE) {
  output_format <- .module1_output_format(output_format)
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
  stream_prediction_stats <- isFALSE(return_prediction_stats) && !is.null(prediction_stats_out_dir)

  make_empty <- function() {
    empty_stats <- .module1_build_prediction_stats(
      correlation_stats = tibble::tibble(
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
    list(
      prediction_stats = empty_stats,
      prediction_pair_count = pair_count,
      prediction_stats_manifest = NULL,
      prediction_stats_manifest_path = NULL,
      n_prediction_stats = 0L
    )
  }
  if (!nrow(high) || !length(tfs)) return(make_empty())

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
  pair_count <- as.double(nrow(high)) * as.double(length(tfs))
  if (!nrow(high) || !length(tfs)) return(make_empty())

  fp_rank <- .module1_rank_matrix_rows(fp_mat)
  rna_rank <- .module1_rank_matrix_rows(rna_mat)
  cutoffs <- .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
  cores <- .module1_default_cores(cores)

  tf_chunk_size <- getOption("craftgrn.module1_prediction_tf_chunk_size", 25L)
  tf_chunk_size <- suppressWarnings(as.integer(tf_chunk_size)[[1L]])
  if (!is.finite(tf_chunk_size) || tf_chunk_size < 1L) tf_chunk_size <- 50L
  chunks <- split(seq_along(tfs), ceiling(seq_along(tfs) / tf_chunk_size))
  prediction_stats_chunks <- if (stream_prediction_stats) NULL else vector("list", length(chunks))
  manifest_chunks <- vector("list", length(chunks))

  for (chunk_i in seq_along(chunks)) {
    ii <- chunks[[chunk_i]]
    if (isTRUE(verbose) && length(chunks) > 1L) {
      .log_inform("Module 1 prediction chunk {chunk_i}/{length(chunks)}: {length(ii)} TF(s).")
    }
    prediction_stats_i <- tryCatch(
      .module1_compute_prediction_stats_cpp(
        high = high,
        tfs = tfs[ii],
        fp_mat = fp_mat,
        rna_mat = rna_mat[ii, , drop = FALSE],
        fp_rank = fp_rank,
        rna_rank = rna_rank[ii, , drop = FALSE],
        cutoffs = cutoffs,
        cores = cores,
        tf_chunk_size = length(ii)
      ),
      error = function(e) NULL
    )
    if (is.null(prediction_stats_i)) {
      prediction_stats_i <- .module1_compute_prediction_stats_by_tf(
        high = high,
        tfs = tfs[ii],
        fp_mat = fp_mat,
        rna_mat = rna_mat[ii, , drop = FALSE],
        fp_rank = fp_rank,
        rna_rank = rna_rank[ii, , drop = FALSE],
        cutoffs = cutoffs,
        min_non_na = min_non_na,
        cores = cores
      )
    }
    prediction_stats_i <- .module1_build_prediction_stats(
      correlation_stats = prediction_stats_i,
      high_confidence_footprints = high,
      omics_data = omics_data
    )
    if (stream_prediction_stats) {
      manifest_chunks[[chunk_i]] <- .module1_write_prediction_stats_chunk(
        x = prediction_stats_i,
        chunk_id = chunk_i,
        out_dir = prediction_stats_out_dir,
        output_format = output_format
      )
    } else {
      prediction_stats_chunks[[chunk_i]] <- prediction_stats_i
    }
    rm(prediction_stats_i)
  }

  if (stream_prediction_stats) {
    prediction_stats_manifest <- dplyr::bind_rows(manifest_chunks)
    prediction_stats_manifest_path <- .module1_write_prediction_stats_manifest(prediction_stats_manifest, prediction_stats_out_dir)
    prediction_stats <- tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      best_r = numeric(),
      best_method = character(),
      condition_support = integer()
    )
    n_prediction_stats <- sum(prediction_stats_manifest$n_rows)
  } else {
    prediction_stats <- dplyr::bind_rows(prediction_stats_chunks)
    prediction_stats_manifest <- NULL
    prediction_stats_manifest_path <- NULL
    n_prediction_stats <- nrow(prediction_stats)
  }

  list(
    prediction_stats = prediction_stats,
    prediction_pair_count = pair_count,
    prediction_stats_manifest = prediction_stats_manifest,
    prediction_stats_manifest_path = prediction_stats_manifest_path,
    n_prediction_stats = n_prediction_stats
  )
}
#' Predict transcription factor binding sites from matched footprint and RNA data
#'
#' Run the Module 1 TFBS workflow as one user-facing operation. The function
#' first uses motif-supported FP-TF correlations to define high-confidence
#' footprints, then predicts sparse FP-TF binding events for expressed TFs.
#'
#' @param omics_data CraftGRN multiomic object returned by `load_prep_multiomic_data()`.
#' @param out_dir Optional output directory. Required when `write_outputs = TRUE`.
#' @param db Motif database label used in output metadata.
#' @param label_col Metadata column used to build condition-level matrices when
#'   missing from `omics_data`.
#' @param r_cutoff Minimum positive correlation used for motif-supported and
#'   prediction calls.
#' @param p_cutoff Optional best-method p-value cutoff. If `NULL`, p-value
#'   filtering is disabled.
#' @param fdr_cutoff Optional best-method adjusted p-value cutoff. If `NULL`,
#'   FDR filtering is disabled.
#' @param filter_to_canonical_bound Logical; if `TRUE`, only footprints with
#'   at least one motif-supported TF passing the cutoffs are used for the
#'   all-expressed-TF prediction stage.
#' @param tf_subset Optional TF subset.
#' @param write_outputs Write Module 1 output files. Defaults to `FALSE` so
#'   calls do not write into the working directory unless an output directory is
#'   explicitly supplied.
#' @param write_stats Retain and write full FP-TF correlation statistics.
#' @param write_bed Write optional BED-like browser files for high-confidence footprints and in-memory TFBS prediction statistics.
#' @param write_qc_report Write a Module 1 HTML QC report when outputs are written.
#' @param qc_report_scan Scan predicted TFBS chunks for top-TF summaries in the QC report.
#' @param output_format Output format for large streamed TFBS prediction statistic chunks.
#' @param return_prediction_stats Return the TFBS prediction statistic table in memory. If `NULL`,
#'   large output-writing runs are streamed to disk and return a manifest.
#' @param prediction_return_limit Maximum number of predicted events to keep in
#'   memory when `return_prediction_stats = NULL` and `write_outputs = TRUE`.
#' @param min_non_na Minimum finite condition pairs required for correlation.
#' @param cores Number of worker cores for the dense prediction correlation step.
#'   If `NULL`, use available cores.
#' @param verbose Emit concise progress messages.
#' @param time_log Logical; if TRUE, emit elapsed-time messages.
#'
#' @return A list containing `omics_data`, `high_confidence_footprints`,
#'   `motif_supported_correlations`, `prediction_stats`, `predicted_tfbs`,
#'   `prediction_stats_manifest`, `reports`, and `parameters`.
#' @export
predict_tfbs <- function(omics_data,
                         out_dir = NULL,
                         db = "JASPAR2024",
                         label_col = NULL,
                         r_cutoff = 0.3,
                         p_cutoff = NULL,
                         fdr_cutoff = NULL,
                         filter_to_canonical_bound = TRUE,
                         tf_subset = NULL,
                         write_outputs = FALSE,
                         write_stats = FALSE,
                         write_bed = FALSE,
                         write_qc_report = TRUE,
                         qc_report_scan = FALSE,
                         output_format = c("csv", "parquet", "auto"),
                         return_prediction_stats = NULL,
                         prediction_return_limit = getOption("craftgrn.module1_prediction_return_limit", 5000000),
                         min_non_na = 3L,
                         cores = NULL,
                         verbose = TRUE,
                         time_log = verbose) {
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a prepared Module 1 list.")
  stopifnot(is.logical(write_outputs), length(write_outputs) == 1L, !is.na(write_outputs))
  if (isTRUE(write_outputs) && (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir))) {
    .log_abort("`out_dir` must be a non-empty path when `write_outputs = TRUE`.")
  }
  if (!is.null(out_dir) && (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir))) {
    .log_abort("`out_dir` must be NULL or a non-empty path.")
  }
  if (!is.character(db) || length(db) != 1L || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  output_format <- .module1_output_format(output_format)
  stopifnot(is.logical(write_qc_report), length(write_qc_report) == 1L, !is.na(write_qc_report))
  stopifnot(is.logical(qc_report_scan), length(qc_report_scan) == 1L, !is.na(qc_report_scan))
  if (!is.null(return_prediction_stats)) {
    stopifnot(is.logical(return_prediction_stats), length(return_prediction_stats) == 1L, !is.na(return_prediction_stats))
  }
  prediction_return_limit <- suppressWarnings(as.numeric(prediction_return_limit)[[1L]])
  if (!is.finite(prediction_return_limit) || prediction_return_limit < 0) {
    .log_abort("`prediction_return_limit` must be a non-negative number.")
  }
  cutoffs <- .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
  r_cutoff <- cutoffs$r
  p_cutoff <- if (is.infinite(cutoffs$p)) NULL else cutoffs$p
  fdr_cutoff <- if (is.infinite(cutoffs$fdr)) NULL else cutoffs$fdr
  stopifnot(is.logical(filter_to_canonical_bound), length(filter_to_canonical_bound) == 1L, !is.na(filter_to_canonical_bound))
  omics_data <- .module1_prepare_predict_omics(
    omics_data = omics_data,
    label_col = label_col,
    verbose = verbose
  )

  cond_cols <- .module1_condition_columns(omics_data)
  expressed_tfs <- .module1_expressed_tfs(omics_data, tf_subset = tf_subset)
  fp_universe <- unique(as.character(omics_data$fp_score_condition_qn$peak_ID))
  fp_universe <- fp_universe[!is.na(fp_universe) & nzchar(fp_universe)]
  bound_fp_ids <- fp_universe
  if (is.data.frame(omics_data$fp_bound_condition) && "peak_ID" %in% names(omics_data$fp_bound_condition)) {
    bound_cols <- intersect(setdiff(names(omics_data$fp_bound_condition), "peak_ID"), cond_cols)
    if (length(bound_cols)) {
      bound_mat <- as.matrix(omics_data$fp_bound_condition[, bound_cols, drop = FALSE])
      bound_mat[is.na(bound_mat)] <- 0L
      bound_fp_ids <- as.character(omics_data$fp_bound_condition$peak_ID[rowSums(bound_mat > 0, na.rm = TRUE) > 0])
    }
  }
  if (isTRUE(verbose)) {
    .log_inform("Module 1 preprocessing gates: {length(fp_universe)} FP(s), {length(bound_fp_ids)} bound/accessibility-supported FP(s), {length(expressed_tfs)} expressed TF(s).")
  }
  fp_universe <- intersect(fp_universe, bound_fp_ids)

  if (isTRUE(verbose)) .log_inform("Module 1 TFBS prediction: motif-supported FP-TF correlations.")
  motif_supported_pairs <- .module1_motif_supported_pairs(omics_data, tf_subset = expressed_tfs)
  motif_supported_pairs <- motif_supported_pairs[motif_supported_pairs$fp_id %in% fp_universe, , drop = FALSE]
  if (isTRUE(verbose)) {
    .log_inform("Module 1 TFBS prediction: {nrow(motif_supported_pairs)} motif-supported FP-TF pair(s).")
  }
  motif_supported_correlations_raw <- .module1_compute_pair_correlations(
    omics_data = omics_data,
    pairs = motif_supported_pairs,
    min_non_na = min_non_na,
    cores = cores
  )
  motif_supported_correlations <- .module1_merge_correlation_stats(
    pearson_stats = motif_supported_correlations_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
    spearman_stats = motif_supported_correlations_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff
  )
  high_confidence_footprints <- .module1_select_high_confidence_footprints(
    motif_supported_correlations = motif_supported_correlations,
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff
  )
  canonical_fp_ids <- unique(high_confidence_footprints$fp_id)
  if (isTRUE(filter_to_canonical_bound)) {
    prediction_footprints <- high_confidence_footprints
  } else {
    ann_base <- omics_data$fp_annotation[!duplicated(omics_data$fp_annotation$fp_peak), c("fp_peak", "atac_peak"), drop = FALSE]
    ann_base <- ann_base[ann_base$fp_peak %in% fp_universe, , drop = FALSE]
    coords <- .module1_parse_fp_coordinates(ann_base$fp_peak)
    prediction_footprints <- tibble::as_tibble(cbind(
      tibble::tibble(fp_id = as.character(ann_base$fp_peak)),
      coords,
      tibble::tibble(atac_peak = as.character(ann_base$atac_peak))
    ))
    prediction_footprints <- dplyr::left_join(
      prediction_footprints,
      high_confidence_footprints[, c("fp_id", "n_canonical_bound_tfs", "canonical_bound_tfs", "canonical_bound_motifs"), drop = FALSE],
      by = "fp_id"
    )
    prediction_footprints$n_canonical_bound_tfs[is.na(prediction_footprints$n_canonical_bound_tfs)] <- 0L
    prediction_footprints$canonical_bound_tfs[is.na(prediction_footprints$canonical_bound_tfs)] <- ""
    prediction_footprints$canonical_bound_motifs[is.na(prediction_footprints$canonical_bound_motifs)] <- ""
    prediction_footprints$has_canonical_bound <- prediction_footprints$fp_id %in% canonical_fp_ids
  }
  if (isTRUE(verbose)) {
    .log_inform("Module 1 canonical-bound FP-TF pairs passing cutoffs: {sum(motif_supported_correlations$pass %in% TRUE)}.")
    .log_inform("Module 1 FPs with at least one canonical-bound TF: {length(canonical_fp_ids)} / {length(fp_universe)}.")
    if (isTRUE(filter_to_canonical_bound)) {
      .log_inform("Module 1 canonical-bound FP filter removed {length(setdiff(fp_universe, canonical_fp_ids))} FP(s).")
    } else {
      .log_inform("Module 1 canonical-bound FP filter disabled; retaining {nrow(prediction_footprints)} FP(s) for prediction.")
    }
  }

  if (isTRUE(verbose)) .log_inform("Module 1 TFBS prediction: prediction correlations.")
  prediction_pair_count <- 0
  keep_prediction_stats <- if (is.null(return_prediction_stats)) {
    !isTRUE(write_outputs) || (as.double(nrow(prediction_footprints)) * as.double(length(expressed_tfs))) <= prediction_return_limit
  } else {
    isTRUE(return_prediction_stats)
  }
  prediction_stats_chunk_dir <- if (isTRUE(write_outputs) && !isTRUE(keep_prediction_stats) && !isTRUE(write_stats)) {
    file.path(out_dir, "module1_prediction_stats_chunks")
  } else {
    NULL
  }
  if (isTRUE(verbose) && !isTRUE(keep_prediction_stats) && !isTRUE(write_stats)) {
    .log_inform("Module 1 TFBS prediction statistics will be streamed to disk in {output_format} chunks.")
  }

  if (isTRUE(write_stats)) {
    prediction_pairs <- .module1_all_prediction_pairs(
      omics_data = omics_data,
      fp_ids = prediction_footprints$fp_id,
      tf_subset = expressed_tfs
    )
    prediction_pair_count <- nrow(prediction_pairs)
    if (isTRUE(verbose)) {
      .log_inform("Module 1 TFBS prediction: {prediction_pair_count} prediction FP-TF pair(s).")
    }
    prediction_stats_raw <- .module1_compute_pair_correlations(
      omics_data = omics_data,
      pairs = prediction_pairs,
      min_non_na = min_non_na,
      cores = cores
    )
    prediction_stats <- .module1_merge_correlation_stats(
      pearson_stats = prediction_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
      spearman_stats = prediction_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
      r_cutoff = r_cutoff,
      p_cutoff = p_cutoff,
      fdr_cutoff = fdr_cutoff
    )
    prediction_stats <- .module1_build_prediction_stats(
      correlation_stats = prediction_stats,
      high_confidence_footprints = prediction_footprints,
      omics_data = omics_data
    )
    prediction_stats_manifest <- NULL
    prediction_stats_manifest_path <- NULL
    n_prediction_stats <- nrow(prediction_stats)
  } else {
    streamed <- .module1_predict_tfbs_streamed(
      omics_data = omics_data,
      high_confidence_footprints = prediction_footprints,
      r_cutoff = r_cutoff,
      p_cutoff = p_cutoff,
      fdr_cutoff = fdr_cutoff,
      tf_subset = expressed_tfs,
      min_non_na = min_non_na,
      cores = cores,
      return_prediction_stats = keep_prediction_stats,
      prediction_stats_out_dir = prediction_stats_chunk_dir,
      output_format = output_format,
      verbose = verbose
    )
    prediction_stats <- streamed$prediction_stats
    prediction_pair_count <- streamed$prediction_pair_count
    prediction_stats_manifest <- streamed$prediction_stats_manifest
    prediction_stats_manifest_path <- streamed$prediction_stats_manifest_path
    n_prediction_stats <- streamed$n_prediction_stats
    if (isTRUE(verbose)) {
      .log_inform("Module 1 TFBS prediction: {streamed$prediction_pair_count} prediction FP-TF pair(s) evaluated without materializing all pairs.")
    }
  }

  predicted_tfbs_paths <- NULL
  if (!is.null(prediction_stats_manifest_path) && !isTRUE(keep_prediction_stats)) {
    if (isTRUE(verbose)) .log_inform("Module 1 predicted TFBS handoff: writing chunks.")
    predicted_tfbs_paths <- .write_predicted_tfbs_from_prediction_stats_manifest(prediction_stats_manifest, out_dir = out_dir, output_format = output_format)
    predicted_tfbs <- tibble::tibble(
      tfbs_id = character(),
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      condition_support = integer()
    )
    n_predicted_tfbs <- predicted_tfbs_paths$n_rows
  } else {
    predicted_tfbs <- build_predicted_tfbs(prediction_stats)
    n_predicted_tfbs <- nrow(predicted_tfbs)
  }

  qc_summary <- list(
    n_fp_input = length(fp_universe),
    n_fp_bound_accessible = length(bound_fp_ids),
    n_expressed_tfs = length(expressed_tfs),
    n_motif_supported_pairs = nrow(motif_supported_pairs),
    n_canonical_pairs_pass = sum(motif_supported_correlations$pass %in% TRUE),
    n_canonical_bound_fps = length(canonical_fp_ids),
    n_prediction_fps = nrow(prediction_footprints),
    n_prediction_pairs = as.numeric(prediction_pair_count),
    n_prediction_stats = as.numeric(n_prediction_stats),
    n_predicted_tfbs = as.numeric(n_predicted_tfbs)
  )

  reports <- list()
  if (isTRUE(write_outputs)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    high_path <- file.path(out_dir, "module1_high_confidence_footprints.csv")
    canonical_stats_path <- file.path(out_dir, "module1_canonical_prediction_stats.csv.gz")
    qc_summary_path <- file.path(out_dir, "module1_qc_summary.csv")
    links_path <- file.path(out_dir, "module1_prediction_stats.csv.gz")
    if (is.null(predicted_tfbs_paths)) predicted_tfbs_paths <- .write_predicted_tfbs_table(predicted_tfbs, out_dir = out_dir, output_format = output_format)
    readr::write_csv(high_confidence_footprints, high_path)
    readr::write_csv(motif_supported_correlations, canonical_stats_path)
    readr::write_csv(tibble::tibble(metric = names(qc_summary), value = as.numeric(unlist(qc_summary, use.names = FALSE))), qc_summary_path)
    cache_reports <- .module1_write_canonical_support_cache(
      high_confidence_footprints = high_confidence_footprints,
      qc_summary = qc_summary,
      out_dir = out_dir
    )
    if (!is.null(prediction_stats_manifest_path)) {
      reports <- c(reports, list(
        high_confidence_footprints = high_path,
        canonical_correlation_stats = canonical_stats_path,
        qc_summary = qc_summary_path,
        predicted_tfbs = predicted_tfbs_paths$path,
        predicted_tfbs_manifest = predicted_tfbs_paths$manifest,
        prediction_stats_manifest = prediction_stats_manifest_path,
        prediction_stats_chunks = prediction_stats_chunk_dir
      ), cache_reports)
    } else {
      readr::write_csv(prediction_stats, links_path)
      reports <- c(reports, list(
        high_confidence_footprints = high_path,
        canonical_correlation_stats = canonical_stats_path,
        qc_summary = qc_summary_path,
        predicted_tfbs = predicted_tfbs_paths$path,
        predicted_tfbs_manifest = predicted_tfbs_paths$manifest,
        correlation_stats = links_path
      ), cache_reports)
    }
    if (isTRUE(write_stats)) {
      stats_path <- file.path(out_dir, "module1_prediction_stats.csv.gz")
      readr::write_csv(prediction_stats, stats_path)
      reports$prediction_stats <- stats_path
    }
  }
  if (isTRUE(write_outputs) && isTRUE(write_bed)) {
    bed_reports <- .module1_write_bed_outputs(
      high_confidence_footprints = high_confidence_footprints,
      prediction_stats = prediction_stats,
      out_dir = out_dir
    )
    reports <- c(reports, bed_reports)
  }
  if (isTRUE(write_outputs) && isTRUE(write_qc_report)) {
    qc_report <- build_module1_qc_report(
      module1 = list(
        omics_data = omics_data,
        high_confidence_footprints = high_confidence_footprints,
        motif_supported_correlations = motif_supported_correlations,
        prediction_stats = prediction_stats,
        predicted_tfbs = predicted_tfbs,
        prediction_stats_manifest = prediction_stats_manifest,
        reports = reports,
        parameters = list(qc_summary = qc_summary)
      ),
      output_dir = file.path(out_dir, "reports"),
      scan_predicted_tfbs = isTRUE(qc_report_scan),
      verbose = verbose
    )
    reports$qc_html <- qc_report
  }

  list(
    omics_data = omics_data,
    high_confidence_footprints = high_confidence_footprints,
    motif_supported_correlations = motif_supported_correlations,
    prediction_stats = prediction_stats,
    predicted_tfbs = predicted_tfbs,
    prediction_stats_manifest = prediction_stats_manifest,
    reports = reports,
    parameters = list(
      db = db,
      label_col = label_col,
      r_cutoff = r_cutoff,
      min_non_na = as.integer(min_non_na),
      cores = .module1_default_cores(cores),
      tf_subset = tf_subset,
      expressed_tfs = expressed_tfs,
      p_cutoff = p_cutoff,
      fdr_cutoff = fdr_cutoff,
      filter_to_canonical_bound = isTRUE(filter_to_canonical_bound),
      write_outputs = isTRUE(write_outputs),
      write_stats = isTRUE(write_stats),
      write_bed = isTRUE(write_bed),
      write_qc_report = isTRUE(write_qc_report),
      qc_report_scan = isTRUE(qc_report_scan),
      output_format = output_format,
      return_prediction_stats = isTRUE(keep_prediction_stats),
      qc_summary = qc_summary
    )
  )
}
NA
NA
NA
