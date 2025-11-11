# File: utils_step1_tf_correlation.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 TF-to-footprint correlation, filtering, and TFBS call
# helpers.
#
# Inputs:
# - prepared Step 1 multi-omic object
# - TF subsets
# - correlation settings and output directories
#
# Outputs:
# - Pearson and Spearman correlation tables
# - all-mode TFBS correlation summaries
# - downstream direct-bound compatible TFBS call summaries
# - per-TF overview tables with one global TFBS call column (`any_bound`) plus
#   condition-specific bound columns
#
# Notes:
# - The rebuilt path keeps the all-mode correlation pass as the primary
#   correlation engine.
# - Direct-bound is a downstream classification layer derived from the
#   all-mode correlation output plus direct-bound support information.
# - Final TFBS calls should be based on the maximum of Pearson and Spearman
#   correlation scores.
# - The overview column `any_bound` means the TFBS passes the global
#   correlation threshold, before applying condition-specific fp/gene filters.

#' Module 1 TF correlation helpers
#'
#' @noRd
NULL

.step1_mode_suffix <- function(mode) {
  if (is.null(mode)) mode <- ""
  mode <- as.character(mode)
  if (!nzchar(mode) || identical(mode, "all")) return("")
  paste0("_", mode)
}

.format_gib <- function(bytes) {
  if (!is.finite(bytes) || is.na(bytes) || bytes < 0) return(NA_character_)
  sprintf("%.1f GiB", as.numeric(bytes) / (1024^3))
}

.get_mem_available_bytes <- function() {
  meminfo <- "/proc/meminfo"
  if (!file.exists(meminfo)) return(NA_real_)
  lines <- tryCatch(readLines(meminfo, warn = FALSE), error = function(e) character(0))
  if (!length(lines)) return(NA_real_)
  hit <- grep("^MemAvailable:", lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", hit[[1]])))
  if (!is.finite(kb)) return(NA_real_)
  kb * 1024
}

.estimate_tf_corr_workload <- function(grn_set, mode, tf_subset = NULL, chunk_size = 5000L, cores = 1L) {
  mode <- match.arg(mode, c("canonical", "all"))
  tf_subset <- if (is.null(tf_subset)) NULL else unique(as.character(tf_subset))

  fp_score_qn <- grn_set$fp_score_condition_qn
  rna_cond <- grn_set$rna_condition
  fp_annotation <- grn_set$fp_annotation

  n_samples <- length(intersect(
    setdiff(names(fp_score_qn), "peak_ID"),
    setdiff(names(rna_cond), c("ensembl_gene_id", "HGNC"))
  ))

  if (identical(mode, "canonical")) {
    ann_exp <- .expand_fp_annotation_tfs(fp_annotation, tfs_col = "tfs")
    if (!is.null(tf_subset) && length(tf_subset)) {
      ann_exp <- ann_exp[ann_exp$tfs %in% tf_subset, , drop = FALSE]
    }
    n_pairs <- nrow(ann_exp)
    n_peaks <- length(unique(ann_exp$fp_peak))
    n_tfs <- length(unique(ann_exp$tfs))
  } else {
    tf_all <- if (is.character(grn_set$tf_list) && length(grn_set$tf_list)) {
      intersect(unique(rna_cond$HGNC), unique(as.character(grn_set$tf_list)))
    } else if (is.data.frame(grn_set$fp_annotation) && "tfs" %in% names(grn_set$fp_annotation)) {
      motif_tfs <- unique(unlist(strsplit(paste(grn_set$fp_annotation$tfs, collapse = ","), ",", fixed = TRUE)))
      motif_tfs <- trimws(as.character(motif_tfs))
      motif_tfs <- motif_tfs[nzchar(motif_tfs)]
      intersect(unique(rna_cond$HGNC), unique(motif_tfs))
    } else {
      unique(rna_cond$HGNC)
    }
    tf_all <- tf_all[!is.na(tf_all) & nzchar(tf_all)]
    if (!is.null(tf_subset) && length(tf_subset)) {
      tf_all <- intersect(tf_all, tf_subset)
    }
    base <- fp_annotation[, c("fp_peak", "motifs"), drop = FALSE]
    base <- base[!duplicated(base[, "fp_peak"]), , drop = FALSE]
    base <- base[!is.na(base$motifs) & nzchar(base$motifs), , drop = FALSE]
    n_peaks <- nrow(base)
    n_tfs <- length(tf_all)
    n_pairs <- as.numeric(n_peaks) * as.numeric(n_tfs)
  }

  chunk_size <- max(1L, as.integer(chunk_size))
  n_chunks <- if (n_pairs > 0) ceiling(n_pairs / chunk_size) else 0L
  cores <- max(1L, as.integer(cores))
  if (identical(mode, "all")) {
    matrix_bytes <- (as.numeric(n_peaks) * as.numeric(n_samples) * 8) +
      (as.numeric(n_tfs) * as.numeric(n_samples) * 8)
    per_tf_pairs <- as.numeric(n_peaks)
    ann_bytes <- per_tf_pairs * 48
    pair_index_bytes <- per_tf_pairs * 24
    corr_bytes <- per_tf_pairs * 24
    est_bytes <- (matrix_bytes + ann_bytes + pair_index_bytes + corr_bytes) * 1.2
  } else {
    matrix_bytes <- (as.numeric(n_peaks) * as.numeric(n_samples) * 8) +
      (as.numeric(n_tfs) * as.numeric(n_samples) * 8)
    ann_bytes <- as.numeric(n_pairs) * 48
    pair_index_bytes <- as.numeric(n_pairs) * 24
    corr_bytes <- as.numeric(n_pairs) * 24
    est_bytes <- (matrix_bytes + ann_bytes + pair_index_bytes + corr_bytes) * 1.2
  }

  list(
    mode = mode,
    n_peaks = as.numeric(n_peaks),
    n_tfs = as.numeric(n_tfs),
    n_pairs = as.numeric(n_pairs),
    n_samples = as.numeric(n_samples),
    chunk_size = chunk_size,
    n_chunks = as.numeric(n_chunks),
    cores = cores,
    est_bytes = as.numeric(est_bytes),
    mem_available_bytes = .get_mem_available_bytes()
  )
}

.warn_tf_corr_workload <- function(workload, verbose = TRUE) {
  if (!isTRUE(verbose) || !is.list(workload)) return(invisible(workload))

  .log_inform(
    paste0(
      "Step 1 TF correlation workload: mode=", workload$mode,
      ", peaks=", format(workload$n_peaks, big.mark = ",", scientific = FALSE),
      ", TFs=", format(workload$n_tfs, big.mark = ",", scientific = FALSE),
      ", pairs=", format(workload$n_pairs, big.mark = ",", scientific = FALSE),
      ", chunks=", format(workload$n_chunks, big.mark = ",", scientific = FALSE),
      ", chunk_size=", format(workload$chunk_size, big.mark = ",", scientific = FALSE),
      ", cores=", workload$cores,
      ", est_RAM=", .format_gib(workload$est_bytes)
    )
  )

  avail <- workload$mem_available_bytes
  if (!is.finite(avail) || is.na(avail) || avail <= 0) return(invisible(workload))

  frac <- workload$est_bytes / avail
  if (is.finite(frac) && frac >= 0.8) {
    .log_warn(
      paste0(
        "Requested TF correlation may be too aggressive: estimated RAM ",
        .format_gib(workload$est_bytes),
        " vs available ",
        .format_gib(avail),
        ". Consider fewer TFs, fewer cores, or running from cache."
      )
    )
  } else if (is.finite(frac) && frac >= 0.5) {
    .log_warn(
      paste0(
        "Requested TF correlation is moderately heavy: estimated RAM ",
        .format_gib(workload$est_bytes),
        " vs available ",
        .format_gib(avail),
        ". The run should be feasible, but watch memory if other jobs are active."
      )
    )
  }

  invisible(workload)
}

.tf_cache_slug <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  if (!nzchar(out)) out <- "tf"
  out
}

.motif_gene_col <- function(motif_db) {
  if (!is.data.frame(motif_db)) return(NULL)
  if ("gene_symbol" %in% names(motif_db)) return("gene_symbol")
  if ("HGNC" %in% names(motif_db)) return("HGNC")
  NULL
}

grn_add_fp_tfs <- function(grn_set, force = FALSE, verbose = TRUE) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_tfs") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_annotation)) {
    .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  }

  annot <- grn_set$fp_annotation
  if (!"fp_peak" %in% names(annot)) .log_abort("`fp_annotation` needs fp_peak.")
  if (!"motifs" %in% names(annot)) .log_abort("`fp_annotation` needs motifs.")

  gene_col <- .motif_gene_col(grn_set$motif_db)
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

.expand_fp_annotation_tfs <- function(fp_annotation, tfs_col = "tfs") {
  if (!is.data.frame(fp_annotation)) .log_abort("`fp_annotation` must be a data.frame.")
  if (!tfs_col %in% names(fp_annotation)) {
    .log_abort("`fp_annotation` is missing column {.val {tfs_col}}.")
  }
  tfs_sym <- rlang::sym(tfs_col)
  fp_annotation |>
    tidyr::separate_rows(!!tfs_sym, sep = "\\s*,\\s*|\\s*::\\s*") |>
    dplyr::filter(!is.na(!!tfs_sym), !!tfs_sym != "")
}

.fp_tf_corr_add_p_adj <- function(fp_annotation, p_col = "corr_fp_tf_p") {
  if (!is.data.frame(fp_annotation)) .log_abort("`fp_annotation` must be a data.frame.")
  need <- c("fp_peak", "tfs", p_col)
  if (!all(need %in% names(fp_annotation))) {
    .log_abort("`fp_annotation` must include columns: {paste(need, collapse = ', ')}.")
  }
  annot <- fp_annotation |>
    dplyr::distinct(.data$fp_peak, .data$tfs, .keep_all = TRUE)
  p_vals <- annot[[p_col]]
  annot$corr_fp_tf_p_adj <- NA_real_
  good_p <- is.finite(p_vals)
  if (any(good_p)) {
    adj <- stats::p.adjust(p_vals[good_p], method = "BH")
    annot$corr_fp_tf_p_adj[good_p] <- pmax(adj, p_vals[good_p])
  }
  annot
}

make_fp_annotation_corr <- function(
    grn_set,
    method = c("pearson", "spearman", "kendall"),
    cores = 10L,
    chunk_size = 5000L,
    min_non_na = 3L
) {
  method <- match.arg(method)
  ann <- grn_set$fp_annotation
  fp <- grn_set$fp_score
  rna <- grn_set$rna

  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id", "HGNC"))
  )
  min_non_na <- as.integer(min_non_na)
  if (length(sample_cols) < min_non_na) {
    .log_abort("Not enough shared samples between fp_score and rna.")
  }

  fp_m <- as.matrix(fp[, sample_cols, drop = FALSE])
  rownames(fp_m) <- fp$peak_ID

  rna_m <- rna |>
    dplyr::filter(!is.na(.data$HGNC), .data$HGNC != "") |>
    dplyr::select("HGNC", dplyr::all_of(sample_cols)) |>
    dplyr::distinct(.data$HGNC, .keep_all = TRUE) |>
    as.data.frame()
  rownames(rna_m) <- rna_m$HGNC
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  fp_id <- match(ann$fp_peak, rownames(fp_m))
  tf_id <- match(ann$tfs, rownames(rna_m))

  pairs_u <- unique(data.frame(fp_id = fp_id, tf_id = tf_id))
  pairs_u <- pairs_u[!is.na(pairs_u$fp_id) & !is.na(pairs_u$tf_id), , drop = FALSE]
  if (!nrow(pairs_u)) {
    out <- ann
    out$corr_fp_tf_r <- NA_real_
    out$corr_fp_tf_p <- NA_real_
    out$corr_fp_tf_p_adj <- NA_real_
    return(out)
  }

  idx <- split(seq_len(nrow(pairs_u)), ceiling(seq_len(nrow(pairs_u)) / chunk_size))
  worker <- function(ii) {
    sub <- pairs_u[ii, , drop = FALSE]
    r <- numeric(nrow(sub))
    p <- numeric(nrow(sub))
    for (i in seq_len(nrow(sub))) {
      x <- fp_m[sub$fp_id[i], ]
      y <- rna_m[sub$tf_id[i], ]
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      if (n < min_non_na) {
        r[i] <- NA_real_
        p[i] <- NA_real_
        next
      }
      ct <- if (identical(method, "pearson")) {
        suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
      } else {
        suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
      }
      r[i] <- as.numeric(unname(ct$estimate))
      p[i] <- as.numeric(ct$p.value)
    }
    data.frame(fp_id = sub$fp_id, tf_id = sub$tf_id, r = r, p = p)
  }

  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)
  corr_tbl <- do.call(rbind, res_list)

  key_u <- paste(corr_tbl$fp_id, corr_tbl$tf_id, sep = "_")
  key_a <- paste(fp_id, tf_id, sep = "_")
  m <- match(key_a, key_u)

  out <- ann
  out$corr_fp_tf_r <- corr_tbl$r[m]
  out$corr_fp_tf_p <- corr_tbl$p[m]
  out$corr_fp_tf_p_adj <- stats::p.adjust(out$corr_fp_tf_p, method = "BH")
  out
}

.build_corr_context <- function(grn_set, min_non_na = 3L) {
  ann <- grn_set$fp_annotation
  fp <- grn_set$fp_score
  rna <- grn_set$rna

  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id", "HGNC"))
  )
  min_non_na <- as.integer(min_non_na)
  if (length(sample_cols) < min_non_na) {
    .log_abort("Not enough shared samples between fp_score and rna.")
  }

  fp_m <- as.matrix(fp[, sample_cols, drop = FALSE])
  rownames(fp_m) <- fp$peak_ID

  rna_m <- rna |>
    dplyr::filter(!is.na(.data$HGNC), .data$HGNC != "") |>
    dplyr::select("HGNC", dplyr::all_of(sample_cols)) |>
    dplyr::distinct(.data$HGNC, .keep_all = TRUE) |>
    as.data.frame()
  rownames(rna_m) <- rna_m$HGNC
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  list(
    sample_cols = sample_cols,
    fp_m = fp_m,
    rna_m = rna_m
  )
}

.make_fp_annotation_corr_one_tf <- function(
    base_ann,
    tf,
    corr_ctx,
    method = c("pearson", "spearman", "kendall"),
    cores = 10L,
    chunk_size = 5000L,
    min_non_na = 3L
) {
  method <- match.arg(method)
  if (!is.data.frame(base_ann) || !nrow(base_ann)) {
    out <- base_ann
    out$tfs <- character(0)
    out$corr_fp_tf_r <- numeric(0)
    out$corr_fp_tf_p <- numeric(0)
    out$corr_fp_tf_p_adj <- numeric(0)
    return(out)
  }

  tf <- as.character(tf)[1]
  tf_idx <- match(tf, rownames(corr_ctx$rna_m))
  out <- base_ann
  out$tfs <- tf

  if (is.na(tf_idx)) {
    out$corr_fp_tf_r <- NA_real_
    out$corr_fp_tf_p <- NA_real_
    out$corr_fp_tf_p_adj <- NA_real_
    return(out)
  }

  fp_id <- match(base_ann$fp_peak, rownames(corr_ctx$fp_m))
  keep <- which(!is.na(fp_id))
  out$corr_fp_tf_r <- NA_real_
  out$corr_fp_tf_p <- NA_real_
  out$corr_fp_tf_p_adj <- NA_real_
  if (!length(keep)) return(out)

  fp_id_use <- fp_id[keep]
  idx <- split(seq_along(fp_id_use), ceiling(seq_along(fp_id_use) / max(1L, as.integer(chunk_size))))
  y_all <- corr_ctx$rna_m[tf_idx, ]

  worker <- function(ii) {
    sub_idx <- fp_id_use[ii]
    r <- numeric(length(sub_idx))
    p <- numeric(length(sub_idx))
    for (i in seq_along(sub_idx)) {
      x <- corr_ctx$fp_m[sub_idx[i], ]
      ok <- is.finite(x) & is.finite(y_all)
      n <- sum(ok)
      if (n < min_non_na) {
        r[i] <- NA_real_
        p[i] <- NA_real_
        next
      }
      ct <- if (identical(method, "pearson")) {
        suppressWarnings(stats::cor.test(x[ok], y_all[ok], method = method))
      } else {
        suppressWarnings(stats::cor.test(x[ok], y_all[ok], method = method, exact = FALSE))
      }
      r[i] <- as.numeric(unname(ct$estimate))
      p[i] <- as.numeric(ct$p.value)
    }
    data.frame(idx = ii, r = r, p = p)
  }

  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)
  corr_tbl <- do.call(rbind, res_list)
  out$corr_fp_tf_r[keep[corr_tbl$idx]] <- corr_tbl$r
  out$corr_fp_tf_p[keep[corr_tbl$idx]] <- corr_tbl$p
  out$corr_fp_tf_p_adj <- stats::p.adjust(out$corr_fp_tf_p, method = "BH")
  out
}

grn_add_fp_tf_corr <- function(
    grn_set,
    method = c("pearson", "spearman", "kendall"),
    mode = c("canonical", "all"),
    tf_subset = NULL,
    cores = 10L,
    chunk_size = 5000L,
    min_non_na = 5L,
    cache_dir = NULL,
    use_cache = TRUE,
    stream_only = FALSE,
    force = FALSE,
    verbose = TRUE
) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  grn_set <- grn_status_init(grn_set)
  status_key <- paste0("fp_tf_corr_", method, "_", mode)
  if (grn_status_is(grn_set, status_key) && !isTRUE(force)) return(grn_set)

  if (!is.data.frame(grn_set$fp_score_condition_qn)) {
    .log_abort("`grn_set$fp_score_condition_qn` is missing; run grn_add_fp_score_qn() first.")
  }
  if (!is.data.frame(grn_set$rna_condition)) {
    .log_abort("`grn_set$rna_condition` is missing; run grn_add_rna_condition() first.")
  }
  if (!is.data.frame(grn_set$fp_annotation)) {
    .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  }
  if (!"tfs" %in% names(grn_set$fp_annotation)) {
    grn_set <- grn_add_fp_tfs(grn_set, force = TRUE, verbose = verbose)
  }

  tf_subset <- if (is.null(tf_subset)) NULL else unique(as.character(tf_subset))
  if (identical(mode, "canonical")) {
    ann_exp <- .expand_fp_annotation_tfs(grn_set$fp_annotation, tfs_col = "tfs")
    if (!is.null(tf_subset) && length(tf_subset)) {
      ann_exp <- ann_exp[ann_exp$tfs %in% tf_subset, , drop = FALSE]
    }
  } else {
    rna_cond <- grn_set$rna_condition
    tf_all <- if (is.character(grn_set$tf_list) && length(grn_set$tf_list)) {
      intersect(unique(rna_cond$HGNC), unique(as.character(grn_set$tf_list)))
    } else if (is.data.frame(grn_set$fp_annotation) && "tfs" %in% names(grn_set$fp_annotation)) {
      motif_tfs <- unique(unlist(strsplit(paste(grn_set$fp_annotation$tfs, collapse = ","), ",", fixed = TRUE)))
      motif_tfs <- trimws(as.character(motif_tfs))
      motif_tfs <- motif_tfs[nzchar(motif_tfs)]
      intersect(unique(rna_cond$HGNC), unique(motif_tfs))
    } else {
      unique(rna_cond$HGNC)
    }
    tf_all <- tf_all[!is.na(tf_all) & nzchar(tf_all)]
    if (!is.null(tf_subset) && length(tf_subset)) {
      tf_all <- intersect(tf_all, tf_subset)
    }
    if (!length(tf_all)) {
      ann_name <- paste0("fp_annotation_", method)
      ann_exp <- grn_set$fp_annotation[0, , drop = FALSE]
      if (!"tfs" %in% names(ann_exp)) ann_exp$tfs <- character(0)
      grn_set[[ann_name]] <- ann_exp
      grn_status_set(grn_set, status_key)
      return(grn_set)
    }
    base <- grn_set$fp_annotation[, c("fp_peak", "atac_peak", "motifs"), drop = FALSE]
    base <- base[!duplicated(base[, c("fp_peak")]), , drop = FALSE]
    base <- base[!is.na(base$motifs) & nzchar(base$motifs), , drop = FALSE]
    ann_exp <- NULL
  }

  if (identical(mode, "canonical") && !nrow(ann_exp)) {
    ann_name <- paste0("fp_annotation_", method)
    grn_set[[ann_name]] <- ann_exp
    grn_status_set(grn_set, status_key)
    return(grn_set)
  }

  fp_score_qn <- grn_set$fp_score_condition_qn
  if (anyDuplicated(fp_score_qn$peak_ID)) {
    fp_score_qn <- fp_score_qn[!duplicated(fp_score_qn$peak_ID), , drop = FALSE]
  }
  rna_cond <- grn_set$rna_condition
  if (anyDuplicated(rna_cond$HGNC)) {
    rna_cond <- rna_cond[!duplicated(rna_cond$HGNC), , drop = FALSE]
  }

  if (identical(mode, "all")) {
    tmp_ctx <- list(fp_annotation = base, fp_score = fp_score_qn, rna = rna_cond)
    corr_ctx <- .build_corr_context(tmp_ctx, min_non_na = min_non_na)
    tf_cache_dir <- NULL
    if (!is.null(cache_dir) && nzchar(cache_dir)) {
      tf_cache_dir <- file.path(cache_dir, paste0(method, "_by_tf"))
      dir.create(tf_cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
    ann_corr_list <- if (isTRUE(stream_only)) NULL else vector("list", length(tf_all))
    for (i in seq_along(tf_all)) {
      tf_i <- tf_all[[i]]
      cache_path <- NULL
      if (!is.null(tf_cache_dir)) {
        cache_path <- file.path(tf_cache_dir, paste0(.tf_cache_slug(tf_i), ".rds"))
      }
      ann_tf <- NULL
      if (isTRUE(use_cache) && is.character(cache_path) && file.exists(cache_path)) {
        ann_tf <- readRDS(cache_path)
      } else {
        ann_tf <- .make_fp_annotation_corr_one_tf(
          base_ann = base,
          tf = tf_i,
          corr_ctx = corr_ctx,
          method = method,
          cores = cores,
          chunk_size = chunk_size,
          min_non_na = min_non_na
        )
        if (is.character(cache_path) && nzchar(cache_path)) {
          saveRDS(ann_tf, cache_path)
        }
      }
      if (!isTRUE(stream_only)) ann_corr_list[[i]] <- ann_tf
      if (isTRUE(verbose) && (i == 1L || i == length(tf_all) || (i %% 10L) == 0L)) {
        .log_inform("Module 1 {.val {method}} correlation: processed {i}/{length(tf_all)} TFs")
      }
      if ((i %% 10L) == 0L) gc(verbose = FALSE)
    }
    ann_corr <- if (isTRUE(stream_only)) NULL else dplyr::bind_rows(ann_corr_list)
  } else {
    tmp <- list(fp_annotation = ann_exp, fp_score = fp_score_qn, rna = rna_cond)
    ann_corr <- make_fp_annotation_corr(
      tmp,
      method = method,
      cores = cores,
      chunk_size = chunk_size,
      min_non_na = min_non_na
    )
    ann_corr <- .fp_tf_corr_add_p_adj(ann_corr, p_col = "corr_fp_tf_p")
  }

  ann_name <- paste0("fp_annotation_", method)
  if (!isTRUE(stream_only)) grn_set[[ann_name]] <- ann_corr
  grn_status_set(grn_set, status_key)
}

grn_filter_fp_tf_corr <- function(
    grn_set,
    method = c("pearson", "spearman", "kendall"),
    mode = c("canonical", "all"),
    r_thr = 0.3,
    p_thr = 0.05,
    set_active = TRUE,
    output_bed = NULL,
    output_bed_condition = NULL,
    label_col = NULL,
    verbose = TRUE
) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  ann_name <- paste0("fp_annotation_", method)
  if (!is.data.frame(grn_set[[ann_name]])) {
    .log_abort("`grn_set` is missing {ann_name}; run grn_add_fp_tf_corr() first.")
  }

  ann_all <- grn_set[[ann_name]]
  if (!"corr_fp_tf_p_adj" %in% names(ann_all)) {
    ann_all <- .fp_tf_corr_add_p_adj(ann_all, p_col = "corr_fp_tf_p")
  }
  grn_set[[ann_name]] <- ann_all

  keep <- is.finite(ann_all$corr_fp_tf_r) & (ann_all$corr_fp_tf_r > r_thr)
  ann_filt <- ann_all[keep, , drop = FALSE]
  grn_set[[paste0(ann_name, "_filtered")]] <- ann_filt
  grn_status_set(grn_set, paste0("fp_tf_corr_filtered_", method))

  if (isTRUE(set_active)) {
    grn_set$fp_annotation <- ann_filt
    if (!nrow(ann_filt)) {
      grn_set$fp_score <- grn_set$fp_score[0, , drop = FALSE]
      grn_set$fp_bound <- grn_set$fp_bound[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_score_condition)) grn_set$fp_score_condition <- grn_set$fp_score_condition[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_bound_condition)) grn_set$fp_bound_condition <- grn_set$fp_bound_condition[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_score_condition_qn)) grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[0, , drop = FALSE]
    } else {
      peaks_keep <- unique(ann_filt$fp_peak)
      grn_set <- filter_fp_rows(
        grn_set = grn_set,
        peaks = peaks_keep,
        score_key = "peak_ID",
        bound_key = "peak_ID",
        annot_key = "fp_peak",
        verbose = FALSE,
        warn_on_mismatch = FALSE
      )
      if (is.data.frame(grn_set$fp_score_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition$peak_ID)
        if (anyNA(ord)) .log_abort("fp_score_condition missing peaks after TF-corr filter.")
        grn_set$fp_score_condition <- grn_set$fp_score_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_bound_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_bound_condition$peak_ID)
        if (anyNA(ord)) .log_abort("fp_bound_condition missing peaks after TF-corr filter.")
        grn_set$fp_bound_condition <- grn_set$fp_bound_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_score_condition_qn)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition_qn$peak_ID)
        if (anyNA(ord)) .log_abort("fp_score_condition_qn missing peaks after TF-corr filter.")
        grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[ord, , drop = FALSE]
      }
    }
  }
  grn_set
}

.merge_fp_tf_corr_methods <- function(ann_pearson, ann_spearman, r_thr = 0.3) {
  keys <- c("fp_peak", "atac_peak", "tfs", "motifs")
  ann_p <- ann_pearson[, unique(c(keys, "corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")), drop = FALSE]
  ann_s <- ann_spearman[, unique(c(keys, "corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")), drop = FALSE]
  names(ann_p)[names(ann_p) %in% c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")] <- c("pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj")
  names(ann_s)[names(ann_s) %in% c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")] <- c("spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj")
  merged <- dplyr::full_join(ann_p, ann_s, by = keys)
  merged$pearson_corr_fp_tf_r[!is.finite(merged$pearson_corr_fp_tf_r)] <- NA_real_
  merged$spearman_corr_fp_tf_r[!is.finite(merged$spearman_corr_fp_tf_r)] <- NA_real_
  merged$corr_fp_tf_r <- pmax(merged$pearson_corr_fp_tf_r, merged$spearman_corr_fp_tf_r, na.rm = TRUE)
  merged$corr_fp_tf_r[!is.finite(merged$corr_fp_tf_r)] <- NA_real_
  merged$corr_fp_tf_p <- ifelse(
    is.finite(merged$pearson_corr_fp_tf_r) & (merged$pearson_corr_fp_tf_r >= merged$spearman_corr_fp_tf_r | !is.finite(merged$spearman_corr_fp_tf_r)),
    merged$pearson_corr_fp_tf_p,
    merged$spearman_corr_fp_tf_p
  )
  merged$corr_fp_tf_p_adj <- ifelse(
    is.finite(merged$pearson_corr_fp_tf_r) & (merged$pearson_corr_fp_tf_r >= merged$spearman_corr_fp_tf_r | !is.finite(merged$spearman_corr_fp_tf_r)),
    merged$pearson_corr_fp_tf_p_adj,
    merged$spearman_corr_fp_tf_p_adj
  )
  merged[is.finite(merged$corr_fp_tf_r) & (merged$corr_fp_tf_r > r_thr), , drop = FALSE]
}

write_grn_tf_corr_outputs <- function(grn_set, out_dir, db, mode = c("canonical", "all")) {
  mode <- match.arg(mode)
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- .step1_mode_suffix(mode)
  readr::write_csv(grn_set$fp_score, file.path(cache_dir, sprintf("fp_score_strict_tf_filtered_corr_%s%s.csv", db, suffix)))
  readr::write_csv(grn_set$fp_bound, file.path(cache_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s%s.csv", db, suffix)))
  readr::write_csv(grn_set$fp_annotation, file.path(cache_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s%s.csv", db, suffix)))
  invisible(out_dir)
}

write_tf_tfbs_overviews <- function(
    omics_data = NULL,
    grn_set = NULL,
    ann_pearson,
    ann_spearman,
    out_dir,
    db,
    mode = c("canonical", "all"),
    label_col = NULL,
    r_thr = 0.3,
    p_thr = 0.05,
    overwrite_tf_overview = FALSE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  if (is.null(omics_data)) omics_data <- grn_set
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  grn_set <- omics_data
  if (!is.data.frame(grn_set$fp_bound_condition)) .log_abort("`omics_data$fp_bound_condition` is missing.")
  if (!is.data.frame(grn_set$rna_expressed)) .log_abort("`omics_data$rna_expressed` is missing.")

  prep_ann <- function(dt, prefix) {
    if (!is.data.frame(dt) || !nrow(dt)) return(data.table::data.table())
    if (!"corr_fp_tf_p_adj" %in% names(dt)) dt <- .fp_tf_corr_add_p_adj(dt, p_col = "corr_fp_tf_p")
    fp_col <- if ("fp_peak" %in% names(dt)) "fp_peak" else "peak_ID"
    tf_col <- if ("tfs" %in% names(dt)) "tfs" else "TF"
    dt <- data.table::as.data.table(dt)
    out <- dt[, .(
      fp_peak = as.character(get(fp_col)),
      atac_peak = as.character(atac_peak),
      motifs = as.character(motifs),
      tfs = as.character(get(tf_col)),
      corr_fp_tf_r = corr_fp_tf_r,
      corr_fp_tf_p = corr_fp_tf_p,
      corr_fp_tf_p_adj = corr_fp_tf_p_adj
    )]
    data.table::setnames(out, c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj"), paste0(prefix, c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")))
    out
  }

  ann_p <- prep_ann(ann_pearson, "pearson_")
  ann_s <- prep_ann(ann_spearman, "spearman_")
  tfs_vec <- sort(unique(c(ann_p$tfs, ann_s$tfs)))
  tfs_vec <- tfs_vec[!is.na(tfs_vec) & nzchar(tfs_vec)]

  cond_cols <- intersect(
    setdiff(names(grn_set$fp_bound_condition), "peak_ID"),
    setdiff(names(grn_set$rna_expressed), c("ensembl_gene_id", "HGNC"))
  )
  fp_bound_tbl <- data.table::as.data.table(grn_set$fp_bound_condition[, c("peak_ID", cond_cols), drop = FALSE])
  fp_bound_tbl <- unique(fp_bound_tbl, by = "peak_ID")
  tf_expr_tbl <- grn_set$rna_expressed[, c("HGNC", cond_cols), drop = FALSE]
  tf_expr_tbl$HGNC <- toupper(tf_expr_tbl$HGNC)
  tf_expr_tbl <- tf_expr_tbl[!duplicated(tf_expr_tbl$HGNC), , drop = FALSE]
  tf_expr_mat <- as.matrix(tf_expr_tbl[, cond_cols, drop = FALSE])
  storage.mode(tf_expr_mat) <- "integer"
  rownames(tf_expr_mat) <- tf_expr_tbl$HGNC

  summary_rows <- list()
  .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
  for (tf in tfs_vec) {
    overview_path <- file.path(out_dir, paste0(.slug(tf), "_overview.txt"))
    if (!isTRUE(overwrite_tf_overview) && file.exists(overview_path)) {
      next
    }
    sub_p <- ann_p[tfs == tf]
    sub_s <- ann_s[tfs == tf]
    sub_all <- if (nrow(sub_p) && nrow(sub_s)) {
      data.table::merge.data.table(sub_p, sub_s, by = c("fp_peak", "atac_peak", "motifs", "tfs"), all = TRUE, sort = FALSE)
    } else if (nrow(sub_p)) {
      data.table::copy(sub_p)
    } else {
      data.table::copy(sub_s)
    }
    if (!nrow(sub_all)) next
    for (cn in c("pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj", "spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj")) {
      if (!cn %in% names(sub_all)) sub_all[[cn]] <- NA_real_
    }
    sub_all[, any_bound := as.integer((is.finite(pearson_corr_fp_tf_r) & pearson_corr_fp_tf_r > r_thr) | (is.finite(spearman_corr_fp_tf_r) & spearman_corr_fp_tf_r > r_thr))]
    fp_s <- data.table::tstrsplit(sub_all[["fp_peak"]], "[:-]", perl = TRUE)
    at_s <- data.table::tstrsplit(sub_all[["atac_peak"]], "[:-]", perl = TRUE)
    sub_all[, `:=`(
      TFBS_chr = fp_s[[1]],
      TFBS_start = suppressWarnings(as.integer(fp_s[[2]])),
      TFBS_end = suppressWarnings(as.integer(fp_s[[3]])),
      TFBS_name = motifs,
      peak_chr = at_s[[1]],
      peak_start = suppressWarnings(as.integer(at_s[[2]])),
      peak_end = suppressWarnings(as.integer(at_s[[3]])),
      TF = tfs
    )]
    tf_idx <- match(toupper(tf), rownames(tf_expr_mat))
    fp_idx <- match(sub_all$fp_peak, fp_bound_tbl$peak_ID)
    fp_ok <- matrix(
      0L,
      nrow = nrow(sub_all),
      ncol = length(cond_cols),
      dimnames = list(NULL, cond_cols)
    )
    good_fp_idx <- which(!is.na(fp_idx))
    if (length(good_fp_idx)) {
      fp_ok[good_fp_idx, ] <- as.matrix(fp_bound_tbl[fp_idx[good_fp_idx], ..cond_cols])
    }
    storage.mode(fp_ok) <- "integer"
    fp_ok[is.na(fp_ok)] <- 0L
    if (is.na(tf_idx)) {
      fp_ok[,] <- 0L
    } else {
      fp_ok <- sweep(fp_ok, 2, as.integer(tf_expr_mat[tf_idx, , drop = TRUE]), `*`)
    }
    fp_ok <- sweep(fp_ok, 1, as.integer(sub_all$any_bound > 0L), `*`)
    colnames(fp_ok) <- paste0(cond_cols, "_bound")
    fp_ok_dt <- data.table::as.data.table(fp_ok)
    data.table::setnames(fp_ok_dt, paste0(cond_cols, "_bound"))
    sub_all[, (names(fp_ok_dt)) := fp_ok_dt]
    keep_cols <- c(
      "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
      "peak_chr", "peak_start", "peak_end", "TF",
      "pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj",
      "spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj",
      "any_bound", paste0(cond_cols, "_bound")
    )
    data.table::fwrite(sub_all[, ..keep_cols], file.path(out_dir, paste0(.slug(tf), "_overview.txt")), sep = "\t", col.names = TRUE, quote = FALSE, na = "NA")
    summary_rows[[length(summary_rows) + 1L]] <- data.table::as.data.table(c(list(TF = tf), stats::setNames(as.list(colSums(fp_ok > 0L, na.rm = TRUE)), cond_cols)))
  }

  if (length(summary_rows)) {
    summary_tbl <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
    summary_path <- file.path(dirname(out_dir), sprintf("06_tf_binding_site_counts_%s%s.csv", db, .step1_mode_suffix(mode)))
    if (!isTRUE(overwrite_tf_overview) && file.exists(summary_path)) {
      old_tbl <- tryCatch(data.table::fread(summary_path, sep = ",", showProgress = FALSE), error = function(e) NULL)
      if (is.data.frame(old_tbl) && nrow(old_tbl)) {
        summary_tbl <- data.table::rbindlist(list(data.table::as.data.table(old_tbl), summary_tbl), use.names = TRUE, fill = TRUE)
        if ("TF" %in% names(summary_tbl)) summary_tbl <- unique(summary_tbl, by = "TF", fromLast = TRUE)
      }
    }
    data.table::fwrite(summary_tbl, summary_path, sep = ",", col.names = TRUE)
  }
  invisible(out_dir)
}

.step1_tf_cache_files <- function(cache_dir, method) {
  dir_path <- file.path(cache_dir, paste0(method, "_by_tf"))
  if (!dir.exists(dir_path)) return(character(0))
  sort(list.files(dir_path, pattern = "\\.rds$", full.names = TRUE))
}

.step1_load_tf_cache <- function(cache_dir, method, tf_slug) {
  path <- file.path(cache_dir, paste0(method, "_by_tf"), paste0(tf_slug, ".rds"))
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

write_tf_tfbs_overviews_from_cache <- function(
    omics_data,
    cache_dir,
    out_dir,
    db,
    mode = c("canonical", "all"),
    label_col = NULL,
    r_thr = 0.3,
    p_thr = 0.05,
    overwrite_tf_overview = FALSE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.character(cache_dir) || !nzchar(cache_dir)) .log_abort("`cache_dir` must be a non-empty path.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.data.frame(omics_data$fp_bound_condition)) .log_abort("`omics_data$fp_bound_condition` is missing.")
  if (!is.data.frame(omics_data$rna_expressed)) .log_abort("`omics_data$rna_expressed` is missing.")

  pear_files <- .step1_tf_cache_files(cache_dir, "pearson")
  spear_files <- .step1_tf_cache_files(cache_dir, "spearman")
  tf_slugs <- sort(unique(c(
    sub("\\.rds$", "", basename(pear_files)),
    sub("\\.rds$", "", basename(spear_files))
  )))
  tf_slugs <- tf_slugs[nzchar(tf_slugs)]
  if (!length(tf_slugs)) {
    .log_abort("No per-TF cache files were found in {cache_dir}.")
  }


  summary_path <- file.path(dirname(out_dir), sprintf("06_tf_binding_site_counts_%s%s.csv", db, .step1_mode_suffix(mode)))
  overview_paths <- file.path(out_dir, paste0(tf_slugs, "_overview.txt"))
  if (!isTRUE(overwrite_tf_overview) && file.exists(summary_path) && all(file.exists(overview_paths))) {
    if (isTRUE(verbose)) .log_inform("Module 1 TFBS overviews already present; skipping cached TF overview rewrite.")
    return(invisible(out_dir))
  
  }
  cond_cols <- intersect(
    setdiff(names(omics_data$fp_bound_condition), "peak_ID"),
    setdiff(names(omics_data$rna_expressed), c("ensembl_gene_id", "HGNC"))
  )
  fp_bound_tbl <- data.table::as.data.table(omics_data$fp_bound_condition[, c("peak_ID", cond_cols), drop = FALSE])
  fp_bound_tbl <- unique(fp_bound_tbl, by = "peak_ID")
  tf_expr_tbl <- omics_data$rna_expressed[, c("HGNC", cond_cols), drop = FALSE]
  tf_expr_tbl$HGNC <- toupper(tf_expr_tbl$HGNC)
  tf_expr_tbl <- tf_expr_tbl[!duplicated(tf_expr_tbl$HGNC), , drop = FALSE]
  tf_expr_mat <- as.matrix(tf_expr_tbl[, cond_cols, drop = FALSE])
  storage.mode(tf_expr_mat) <- "integer"
  rownames(tf_expr_mat) <- tf_expr_tbl$HGNC

  summary_rows <- vector("list", length(tf_slugs))
  .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))

  for (i in seq_along(tf_slugs)) {
    tf_slug <- tf_slugs[[i]]
    overview_path <- file.path(out_dir, paste0(tf_slug, "_overview.txt"))
    if (!isTRUE(overwrite_tf_overview) && file.exists(overview_path)) {
      next
    }
    sub_p <- .step1_load_tf_cache(cache_dir, "pearson", tf_slug)
    sub_s <- .step1_load_tf_cache(cache_dir, "spearman", tf_slug)
    if (!is.data.frame(sub_p) && !is.data.frame(sub_s)) next

    sub_p <- if (is.data.frame(sub_p)) data.table::as.data.table(sub_p) else data.table::data.table()
    sub_s <- if (is.data.frame(sub_s)) data.table::as.data.table(sub_s) else data.table::data.table()

    if (nrow(sub_p) && !"corr_fp_tf_p_adj" %in% names(sub_p)) sub_p <- data.table::as.data.table(.fp_tf_corr_add_p_adj(as.data.frame(sub_p), p_col = "corr_fp_tf_p"))
    if (nrow(sub_s) && !"corr_fp_tf_p_adj" %in% names(sub_s)) sub_s <- data.table::as.data.table(.fp_tf_corr_add_p_adj(as.data.frame(sub_s), p_col = "corr_fp_tf_p"))

    if (nrow(sub_p)) {
      keep_p <- c("fp_peak", "atac_peak", "motifs", "tfs", "corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")
      sub_p <- sub_p[, keep_p, with = FALSE]
      data.table::setnames(sub_p, c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj"), c("pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj"))
    }
    if (nrow(sub_s)) {
      keep_s <- c("fp_peak", "atac_peak", "motifs", "tfs", "corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj")
      sub_s <- sub_s[, keep_s, with = FALSE]
      data.table::setnames(sub_s, c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj"), c("spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj"))
    }

    sub_all <- if (nrow(sub_p) && nrow(sub_s)) {
      data.table::merge.data.table(sub_p, sub_s, by = c("fp_peak", "atac_peak", "motifs", "tfs"), all = TRUE, sort = FALSE)
    } else if (nrow(sub_p)) {
      data.table::copy(sub_p)
    } else {
      data.table::copy(sub_s)
    }
    if (!nrow(sub_all)) next

    for (cn in c("pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj", "spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj")) {
      if (!cn %in% names(sub_all)) sub_all[[cn]] <- NA_real_
    }
    sub_all[, any_bound := as.integer((is.finite(pearson_corr_fp_tf_r) & pearson_corr_fp_tf_r > r_thr) | (is.finite(spearman_corr_fp_tf_r) & spearman_corr_fp_tf_r > r_thr))]
    fp_s <- data.table::tstrsplit(sub_all[["fp_peak"]], "[:-]", perl = TRUE)
    at_s <- data.table::tstrsplit(sub_all[["atac_peak"]], "[:-]", perl = TRUE)
    sub_all[, `:=`(
      TFBS_chr = fp_s[[1]],
      TFBS_start = suppressWarnings(as.integer(fp_s[[2]])),
      TFBS_end = suppressWarnings(as.integer(fp_s[[3]])),
      TFBS_name = motifs,
      peak_chr = at_s[[1]],
      peak_start = suppressWarnings(as.integer(at_s[[2]])),
      peak_end = suppressWarnings(as.integer(at_s[[3]])),
      TF = tfs
    )]

    tf_name <- unique(sub_all$tfs)
    tf_name <- tf_name[!is.na(tf_name) & nzchar(tf_name)]
    tf_name <- if (length(tf_name)) tf_name[[1]] else tf_slug
    tf_idx <- match(toupper(tf_name), rownames(tf_expr_mat))
    fp_idx <- match(sub_all$fp_peak, fp_bound_tbl$peak_ID)
    fp_ok <- matrix(
      0L,
      nrow = nrow(sub_all),
      ncol = length(cond_cols),
      dimnames = list(NULL, cond_cols)
    )
    good_fp_idx <- which(!is.na(fp_idx))
    if (length(good_fp_idx)) {
      fp_ok[good_fp_idx, ] <- as.matrix(fp_bound_tbl[fp_idx[good_fp_idx], ..cond_cols])
    }
    storage.mode(fp_ok) <- "integer"
    fp_ok[is.na(fp_ok)] <- 0L
    if (is.na(tf_idx)) {
      fp_ok[,] <- 0L
    } else {
      fp_ok <- sweep(fp_ok, 2, as.integer(tf_expr_mat[tf_idx, , drop = TRUE]), `*`)
    }
    fp_ok <- sweep(fp_ok, 1, as.integer(sub_all$any_bound > 0L), `*`)
    colnames(fp_ok) <- paste0(cond_cols, "_bound")
    fp_ok_dt <- data.table::as.data.table(fp_ok)
    data.table::setnames(fp_ok_dt, paste0(cond_cols, "_bound"))
    sub_all[, (names(fp_ok_dt)) := fp_ok_dt]

    keep_cols <- c(
      "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
      "peak_chr", "peak_start", "peak_end", "TF",
      "pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj",
      "spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj",
      "any_bound", paste0(cond_cols, "_bound")
    )
    data.table::fwrite(sub_all[, ..keep_cols], file.path(out_dir, paste0(.slug(tf_name), "_overview.txt")), sep = "\t", col.names = TRUE, quote = FALSE, na = "NA")
    summary_rows[[i]] <- data.table::as.data.table(c(list(TF = tf_name), stats::setNames(as.list(colSums(fp_ok > 0L, na.rm = TRUE)), cond_cols)))
    if (isTRUE(verbose) && (i == 1L || i == length(tf_slugs) || (i %% 10L) == 0L)) {
      .log_inform("Module 1 TFBS overview writing: processed {i}/{length(tf_slugs)} TFs")
    }
    if ((i %% 10L) == 0L) gc(verbose = FALSE)
  }

  summary_rows <- Filter(Negate(is.null), summary_rows)
  if (length(summary_rows)) {
    summary_tbl <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
    summary_path <- file.path(dirname(out_dir), sprintf("06_tf_binding_site_counts_%s%s.csv", db, .step1_mode_suffix(mode)))
    if (!isTRUE(overwrite_tf_overview) && file.exists(summary_path)) {
      old_tbl <- tryCatch(data.table::fread(summary_path, sep = ",", showProgress = FALSE), error = function(e) NULL)
      if (is.data.frame(old_tbl) && nrow(old_tbl)) {
        summary_tbl <- data.table::rbindlist(list(data.table::as.data.table(old_tbl), summary_tbl), use.names = TRUE, fill = TRUE)
        if ("TF" %in% names(summary_tbl)) summary_tbl <- unique(summary_tbl, by = "TF", fromLast = TRUE)
      }
    }
    data.table::fwrite(summary_tbl, summary_path, sep = ",", col.names = TRUE)
  }
  invisible(out_dir)
}

#' Correlate TF expression to footprint scores for Module 1 TFBS prediction
#'
#' Run Pearson and Spearman correlations between condition-level TF expression
#' and quantile-normalized footprint scores, then keep TFBS calls using the
#' rebuilt Step 1 rule based on the maximum of the two correlation coefficients.
#' The default rebuilt output structure writes directly under
#' `predict_tf_binding_sites/` without a separate `all/` subfolder.
#'
#' @param omics_data Optional prepared Module 1 multi-omic object.
#' @param grn_set Deprecated synonym for `omics_data`.
#' @param config Optional YAML config path used when `omics_data` is not
#'   supplied.
#' @param genome Optional genome override.
#' @param gene_symbol_col Gene-symbol column in the RNA table.
#' @param fp_aligned Optional aligned footprint object used when rebuilding
#'   `omics_data` inline.
#' @param do_preprocess Logical; if `TRUE`, rebuild Step 1 preprocessing before
#'   correlation.
#' @param do_motif_clustering Logical; if `TRUE`, allow motif clustering during
#'   inline preprocessing.
#' @param fp_root_dir,fp_cache_dir,fp_cache_tag Optional footprint source and
#'   cache settings for inline preprocessing.
#' @param mode Correlation mode. The rebuilt Step 1 path uses `"all"` as the
#'   default and primary mode.
#' @param label_col Metadata column used to aggregate matched conditions.
#' @param out_dir Output directory for Step 1 TFBS files.
#' @param db Motif database label used in output file names.
#' @param r_thr Correlation threshold applied to `max(Pearson_r, Spearman_r)`.
#' @param p_thr Retained for compatibility but not used in the final rebuilt
#'   TFBS call rule.
#' @param tf_subset Optional TF subset for quick downstream correlation checks.
#' @param cores_pearson,cores_spearman Optional worker counts for the two
#'   correlation methods.
#' @param chunk_size Correlation chunk size.
#' @param all_mode_tf_chunk_size,canonical_mode_tf_chunk_size Retained for
#'   compatibility with the legacy interface.
#' @param min_non_na Minimum number of shared finite values required to compute
#'   a correlation.
#' @param qc Logical; if `TRUE`, write correlation QC and overview outputs.
#' @param write_bed Logical; retained for compatibility. The rebuilt Step 1 path
#'   does not yet write BED outputs.
#' @param overwrite_tf_overview Logical; if TRUE, rewrite existing per-TF
#'   overview files. By default existing overview files are kept and skipped.
#' @param write_direct_bound Logical; if `TRUE`, write condition-specific
#'   direct-bound BED files as part of the standard Module 1 workflow.
#' @param direct_bound_conditions Optional condition subset for direct-bound BED
#'   writing. Defaults to all conditions in `fp_bound_condition`.
#' @param write_outputs Logical; if `TRUE`, save filtered score, bound, and
#'   annotation outputs.
#' @param use_cache Logical; if `TRUE`, reuse cached correlation results when
#'   present.
#' @param cache_dir Optional cache directory for per-method correlation RDS
#'   files.
#' @param verbose Logical; if `TRUE`, emit concise progress messages.
#'
#' @return The updated Module 1 multi-omic object with correlation outputs.
#' @examples
#' \dontrun{
#' omics_data <- correlate_tf_to_fp(
#'   omics_data = omics_data,
#'   mode = "all",
#'   out_dir = file.path("predict_tf_binding_sites"),
#'   label_col = "strict_match_rna",
#'   db = "JASPAR2024",
#'   tf_subset = c("HNF1A", "SOX9"),
#'   use_cache = TRUE,
#'   verbose = TRUE
#' )
#' }
#' @export
correlate_tf_to_fp <- function(
    omics_data = NULL,
    grn_set = NULL,
    config = NULL,
    genome = NULL,
    gene_symbol_col = "HGNC",
    fp_aligned = NULL,
    do_preprocess = FALSE,
    do_motif_clustering = FALSE,
    fp_root_dir = NULL,
    fp_cache_dir = NULL,
    fp_cache_tag = NULL,
    mode = c("canonical", "all"),
    label_col,
    out_dir = NULL,
    db = NULL,
    r_thr = 0.3,
    p_thr = 0.05,
    tf_subset = NULL,
    cores_pearson = NULL,
    cores_spearman = NULL,
    chunk_size = 5000L,
    all_mode_tf_chunk_size = 5L,
    canonical_mode_tf_chunk_size = 50L,
    min_non_na = 5L,
    qc = TRUE,
    write_bed = FALSE,
    overwrite_tf_overview = FALSE,
    write_direct_bound = TRUE,
    direct_bound_conditions = NULL,
    write_outputs = TRUE,
    use_cache = TRUE,
    cache_dir = NULL,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  default_cores <- max(1L, as.integer(floor(parallel::detectCores(logical = TRUE) * 0.8)))
  cores_pearson <- if (is.null(cores_pearson) || !is.finite(as.numeric(cores_pearson))) default_cores else max(1L, as.integer(cores_pearson))
  cores_spearman <- if (is.null(cores_spearman) || !is.finite(as.numeric(cores_spearman))) default_cores else max(1L, as.integer(cores_spearman))
  if (is.null(omics_data)) omics_data <- grn_set

  cfg_val <- function(name, default = NULL) .cfg_get(name, default = default)
  if (is.null(omics_data) && is.null(config) && is.null(fp_aligned) && !isTRUE(do_preprocess)) {
    .log_abort("`omics_data` not found. Provide it or set `config` or `do_preprocess`.")
  }
  if (is.null(omics_data) && (!is.null(config) || !is.null(fp_aligned) || isTRUE(do_preprocess))) {
    if (!is.null(config)) {
      if (is.character(config) && length(config) == 1L && file.exists(config)) {
        load_config(config)
      } else {
        .log_abort("`config` must be a path to a YAML file.")
      }
    }
    if (!is.null(genome) && nzchar(genome)) .cfg_set("ref_genome", genome)
    if (is.null(fp_cache_dir)) {
      base_dir_cfg <- cfg_val("base_dir")
      if (is.character(base_dir_cfg) && nzchar(base_dir_cfg)) fp_cache_dir <- file.path(base_dir_cfg, "cache")
    }
    if (is.null(fp_cache_tag)) fp_cache_tag <- cfg_val("db")
    omics_data <- load_prep_multiomic_data(
      config = config,
      genome = genome,
      gene_symbol_col = gene_symbol_col,
      fp_aligned = fp_aligned,
      do_preprocess = do_preprocess,
      do_motif_clustering = do_motif_clustering,
      fp_root_dir = fp_root_dir,
      fp_cache_dir = fp_cache_dir,
      fp_cache_tag = fp_cache_tag,
      label_col = label_col,
      expected_n = cfg_val("expected_n"),
      tf_list = cfg_val("tf_list"),
      motif_db = cfg_val("motif_db"),
      threshold_gene_expr = cfg_val("threshold_gene_expr"),
      threshold_fp_score = cfg_val("threshold_fp_score"),
      use_parallel = TRUE,
      verbose = verbose
    )
  }

  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.character(label_col) || !nzchar(label_col)) .log_abort("`label_col` must be a non-empty string.")
  if (is.null(out_dir) || !is.character(out_dir) || !nzchar(out_dir)) {
    base_dir_cfg <- cfg_val("base_dir")
    if (!is.character(base_dir_cfg) || !nzchar(base_dir_cfg)) .log_abort("`out_dir` is not set and `base_dir` is unavailable.")
    out_dir <- file.path(base_dir_cfg, "predict_tf_binding_sites")
  }
  if (!is.character(db) || !nzchar(db)) db <- cfg_val("db")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cache_dir)) cache_dir <- file.path(out_dir, "cache", "fp_tf_corr")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  grn_set <- omics_data
  grn_set <- grn_add_rna_condition(grn_set, label_col = label_col, verbose = verbose)
  grn_set <- grn_add_fp_tfs(grn_set, verbose = verbose)
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID", verbose = verbose)

  workload <- .estimate_tf_corr_workload(
    grn_set = grn_set,
    mode = mode,
    tf_subset = tf_subset,
    chunk_size = chunk_size,
    cores = max(cores_pearson, cores_spearman)
  )
  .warn_tf_corr_workload(workload, verbose = verbose)

  .load_cache <- function(method) {
    cache_path <- file.path(cache_dir, sprintf("fp_tf_corr_%s%s.rds", method, .step1_mode_suffix(mode)))
    if (isTRUE(use_cache) && file.exists(cache_path)) readRDS(cache_path) else NULL
  }
  .save_cache <- function(method, tbl) {
    cache_path <- file.path(cache_dir, sprintf("fp_tf_corr_%s%s.rds", method, .step1_mode_suffix(mode)))
    if (isTRUE(use_cache) && is.data.frame(tbl)) saveRDS(tbl, cache_path)
  }

  use_streaming_all <- identical(mode, "all") && isTRUE(qc) && isTRUE(use_cache)
  pear_dir <- file.path(cache_dir, "pearson_by_tf")
  spear_dir <- file.path(cache_dir, "spearman_by_tf")
  cache_ready <- FALSE
  if (use_streaming_all && dir.exists(pear_dir) && dir.exists(spear_dir)) {
    n_pear <- length(list.files(pear_dir, pattern = "\\.rds$", full.names = FALSE))
    n_spear <- length(list.files(spear_dir, pattern = "\\.rds$", full.names = FALSE))
    cache_ready <- (n_pear >= workload$n_tfs) && (n_spear >= workload$n_tfs)
  }

  ann_p <- if (use_streaming_all) NULL else .load_cache("pearson")
  if (is.data.frame(ann_p)) {
    grn_set$fp_annotation_pearson <- ann_p
  } else if (use_streaming_all && isTRUE(cache_ready)) {
    if (isTRUE(verbose)) .log_inform("Module 1 pearson cache is complete; skipping all-TF in-memory reload.")
  } else {
    grn_set <- grn_add_fp_tf_corr(grn_set, method = "pearson", mode = mode, tf_subset = tf_subset, cores = cores_pearson, chunk_size = chunk_size, min_non_na = min_non_na, cache_dir = cache_dir, use_cache = use_cache, stream_only = use_streaming_all, verbose = verbose)
    if (!use_streaming_all) .save_cache("pearson", grn_set$fp_annotation_pearson)
  }
  ann_s <- if (use_streaming_all) NULL else .load_cache("spearman")
  if (is.data.frame(ann_s)) {
    grn_set$fp_annotation_spearman <- ann_s
  } else if (use_streaming_all && isTRUE(cache_ready)) {
    if (isTRUE(verbose)) .log_inform("Module 1 spearman cache is complete; skipping all-TF in-memory reload.")
  } else {
    grn_set <- grn_add_fp_tf_corr(grn_set, method = "spearman", mode = mode, tf_subset = tf_subset, cores = cores_spearman, chunk_size = chunk_size, min_non_na = min_non_na, cache_dir = cache_dir, use_cache = use_cache, stream_only = use_streaming_all, verbose = verbose)
    if (!use_streaming_all) .save_cache("spearman", grn_set$fp_annotation_spearman)
  }

  if (!use_streaming_all) {
    grn_set <- grn_filter_fp_tf_corr(grn_set, method = "pearson", mode = mode, r_thr = r_thr, p_thr = p_thr, set_active = FALSE, verbose = verbose)
    grn_set <- grn_filter_fp_tf_corr(grn_set, method = "spearman", mode = mode, r_thr = r_thr, p_thr = p_thr, set_active = FALSE, verbose = verbose)
    combined_filtered <- .merge_fp_tf_corr_methods(grn_set$fp_annotation_pearson, grn_set$fp_annotation_spearman, r_thr = r_thr)
    grn_set$fp_annotation_combined <- combined_filtered
    grn_set$fp_annotation <- combined_filtered

    if (!nrow(combined_filtered)) {
      grn_set$fp_score <- grn_set$fp_score[0, , drop = FALSE]
      grn_set$fp_bound <- grn_set$fp_bound[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_score_condition)) grn_set$fp_score_condition <- grn_set$fp_score_condition[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_bound_condition)) grn_set$fp_bound_condition <- grn_set$fp_bound_condition[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_score_condition_qn)) grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[0, , drop = FALSE]
    } else {
      peaks_keep <- unique(combined_filtered$fp_peak)
      grn_set <- filter_fp_rows(grn_set = grn_set, peaks = peaks_keep, score_key = "peak_ID", bound_key = "peak_ID", annot_key = "fp_peak", verbose = FALSE, warn_on_mismatch = FALSE)
      if (is.data.frame(grn_set$fp_score_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition$peak_ID)
        grn_set$fp_score_condition <- grn_set$fp_score_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_bound_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_bound_condition$peak_ID)
        grn_set$fp_bound_condition <- grn_set$fp_bound_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_score_condition_qn)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition_qn$peak_ID)
        grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[ord, , drop = FALSE]
      }
    }
  }

  if (isTRUE(write_bed)) {
    .log_warn("`write_bed = TRUE` is not implemented yet in the rebuilt Step 1 path.")
  }
  if (isTRUE(write_outputs) && !use_streaming_all) {
    write_grn_tf_corr_outputs(grn_set, out_dir = out_dir, db = db, mode = mode)
  }
  if (isTRUE(qc)) {
    pdf_out_path <- file.path(out_dir, sprintf("06_tf_corr_stats_%s%s.pdf", db, if (identical(mode, "all")) "" else paste0("_", mode)))
    if (use_streaming_all && !isTRUE(overwrite_tf_overview) && file.exists(pdf_out_path)) {
      if (isTRUE(verbose)) .log_inform("Module 1 TF correlation stats PDF already present; skipping PDF rewrite.")
    } else if (use_streaming_all && exists("plot_tf_corr_stats_pdf_from_cache", mode = "function")) {
      plot_tf_corr_stats_pdf_from_cache(
        cache_dir = cache_dir,
        out_dir = out_dir,
        db = db,
        mode = mode,
        r_thr = r_thr,
        p_thr = p_thr,
        verbose = verbose
      )
    } else if (exists("plot_tf_corr_stats_pdf", mode = "function")) {
      plot_tf_corr_stats_pdf(
        ann_pearson = grn_set$fp_annotation_pearson_filtered,
        ann_spearman = grn_set$fp_annotation_spearman_filtered,
        out_dir = out_dir,
        db = db,
        mode = mode,
        r_thr = r_thr,
        p_thr = p_thr,
        ann_pearson_all = grn_set$fp_annotation_pearson,
        ann_spearman_all = grn_set$fp_annotation_spearman,
        verbose = verbose
      )
    }
    overview_dir <- file.path(out_dir, sprintf("06_fp_predicted_tfbs_%s%s", db, .step1_mode_suffix(mode)))
    if (use_streaming_all) {
      write_tf_tfbs_overviews_from_cache(
        omics_data = grn_set,
        cache_dir = cache_dir,
        out_dir = overview_dir,
        db = db,
        mode = mode,
        label_col = label_col,
        r_thr = r_thr,
        p_thr = p_thr,
        verbose = verbose
      )
    } else {
      write_tf_tfbs_overviews(
        omics_data = grn_set,
        ann_pearson = grn_set$fp_annotation_pearson,
        ann_spearman = grn_set$fp_annotation_spearman,
        out_dir = overview_dir,
        db = db,
        mode = mode,
        label_col = label_col,
        r_thr = r_thr,
        p_thr = p_thr,
        overwrite_tf_overview = overwrite_tf_overview,
        verbose = verbose
      )
    }
  }
  if (isTRUE(write_direct_bound)) {
    if (!identical(mode, "all")) {
      if (isTRUE(verbose)) {
        .log_warn("Direct-bound BED writing is only supported for mode = 'all'; skipping.")
      }
    } else if (!exists("write_direct_bound_filter_bed", mode = "function")) {
      .log_warn("Direct-bound BED writer is unavailable; skipping.")
    } else {
      genome_use <- cfg_val("ref_genome")
      if (!is.character(genome_use) || !nzchar(genome_use)) {
        .log_abort("`ref_genome` must be available to write direct-bound BED files.")
      }
      if (is.null(direct_bound_conditions)) {
        if (!is.data.frame(grn_set$fp_bound_condition)) {
          .log_abort("`fp_bound_condition` is required to derive direct-bound conditions.")
        }
        direct_bound_conditions <- setdiff(names(grn_set$fp_bound_condition), "peak_ID")
        direct_bound_conditions <- sub("_bound$", "", direct_bound_conditions)
      }
      write_direct_bound_filter_bed(
        step1_out_dir = out_dir,
        conditions = direct_bound_conditions,
        db = db,
        genome = genome_use,
        r_cutoff = r_thr,
        tf_subset = tf_subset,
        verbose = verbose
      )
    }
  }
  grn_set
}

output_per_condition_tfbs_bed_files <- function(
    omics_data,
    out_dir,
    db,
    label_col,
    mode = c("canonical", "all"),
    r_thr = 0.3,
    p_thr = 0.05,
    cores_pearson = NULL,
    cores_spearman = NULL,
    chunk_size = 5000L,
    min_non_na = 5L,
    use_cache = TRUE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  correlate_tf_to_fp(
    omics_data = omics_data,
    mode = mode,
    out_dir = out_dir,
    label_col = label_col,
    r_thr = r_thr,
    p_thr = p_thr,
    db = db,
    cores_pearson = cores_pearson,
    cores_spearman = cores_spearman,
    chunk_size = chunk_size,
    min_non_na = min_non_na,
    qc = FALSE,
    write_bed = TRUE,
    write_outputs = FALSE,
    use_cache = use_cache,
    verbose = verbose
  )
}

experimental_benchmarking_of_tfbs_predictions <- function(
    omics_data,
    out_dir,
    db,
    label_col,
    r_thr = 0.3,
    p_thr = 1,
    cores_pearson = NULL,
    cores_spearman = NULL,
    chunk_size = 5000L,
    min_non_na = 5L,
    use_cache = TRUE,
    verbose = TRUE
) {
  correlate_tf_to_fp(
    omics_data = omics_data,
    mode = "all",
    out_dir = out_dir,
    label_col = label_col,
    r_thr = r_thr,
    p_thr = p_thr,
    db = db,
    cores_pearson = cores_pearson,
    cores_spearman = cores_spearman,
    chunk_size = chunk_size,
    min_non_na = min_non_na,
    qc = TRUE,
    write_bed = FALSE,
    write_outputs = TRUE,
    use_cache = use_cache,
    verbose = verbose
  )
}
