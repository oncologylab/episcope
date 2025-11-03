# File: utils_step1_condition_matrices.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 helpers for building condition-level RNA, ATAC, and
# footprint matrices and flags.
#
# Inputs:
# - Step 1 multi-omic object
# - thresholds and condition label settings
#
# Outputs:
# - condition matrices
# - bound or expressed flag matrices
# - normalized footprint score matrices
#
# Notes:
# - Keep thresholding logic and matrix construction localized here.

#' Module 1 condition matrix helpers
#'
#' @noRd
NULL

make_rna_expressed <- function(
    rna_tbl,
    metadata,
    label_col,
    threshold_gene_expr = 10,
    min_samples = 1L
) {
  if (!is.data.frame(rna_tbl)) .log_abort("`rna_tbl` must be a data.frame.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl))) {
    .log_abort("`rna_tbl` must include ensembl_gene_id and HGNC.")
  }
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  if (!is.numeric(threshold_gene_expr) || length(threshold_gene_expr) != 1L) {
    .log_abort("`threshold_gene_expr` must be a single numeric value.")
  }
  min_samples <- as.integer(min_samples)
  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 1L) {
    .log_abort("`min_samples` must be a positive integer.")
  }

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(rna_tbl))
  if (!length(ids_present)) .log_abort("No overlap between metadata ids and RNA columns.")
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  out <- tibble::tibble(
    ensembl_gene_id = rna_tbl$ensembl_gene_id,
    HGNC = rna_tbl$HGNC
  )

  for (cond in cond_levels) {
    ids_cond <- meta_use$id[meta_use[[label_col]] == cond]
    ids_cond <- intersect(ids_cond, names(rna_tbl))
    if (!length(ids_cond)) {
      .log_warn("No RNA columns found for condition {.val {cond}}.")
      out[[cond]] <- integer(nrow(rna_tbl))
      next
    }
    mat <- as.matrix(rna_tbl[, ids_cond, drop = FALSE])
    flag <- rowSums(mat >= threshold_gene_expr, na.rm = TRUE) >= min_samples
    out[[cond]] <- as.integer(flag)
  }
  out
}

make_rna_condition <- function(
    rna_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(rna_tbl)) .log_abort("`rna_tbl` must be a data.frame.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl))) {
    .log_abort("`rna_tbl` must include ensembl_gene_id and HGNC.")
  }
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  agg_fun <- match.arg(agg_fun)

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(rna_tbl))
  if (!length(ids_present)) .log_abort("No overlap between metadata ids and RNA columns.")
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  out <- tibble::tibble(
    ensembl_gene_id = rna_tbl$ensembl_gene_id,
    HGNC = rna_tbl$HGNC
  )

  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], names(rna_tbl))
      if (!length(ids_cond)) {
        .log_warn("No RNA columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      out[[cond]] <- rna_tbl[[ids_cond[[1L]]]]
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], names(rna_tbl))
      if (!length(ids_cond)) {
        .log_warn("No RNA columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      mat <- as.matrix(rna_tbl[, ids_cond, drop = FALSE])
      if (identical(agg_fun, "median")) {
        out[[cond]] <- apply(mat, 1L, stats::median, na.rm = TRUE)
      } else {
        out[[cond]] <- rowMeans(mat, na.rm = TRUE)
      }
    }
  }
  out
}

make_fp_score_condition <- function(
    fp_score_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(fp_score_tbl)) .log_abort("`fp_score_tbl` must be a data.frame.")
  if (!"peak_ID" %in% names(fp_score_tbl)) .log_abort("`fp_score_tbl` needs a peak_ID column.")
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  agg_fun <- match.arg(agg_fun)

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(fp_score_tbl))
  if (!length(ids_present)) .log_abort("No overlap between metadata ids and fp_score columns.")
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  out <- tibble::tibble(peak_ID = fp_score_tbl$peak_ID)
  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], names(fp_score_tbl))
      if (!length(ids_cond)) {
        .log_warn("No fp_score columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      out[[cond]] <- fp_score_tbl[[ids_cond[[1L]]]]
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], names(fp_score_tbl))
      if (!length(ids_cond)) {
        .log_warn("No fp_score columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      mat <- as.matrix(fp_score_tbl[, ids_cond, drop = FALSE])
      if (identical(agg_fun, "median")) {
        out[[cond]] <- apply(mat, 1L, stats::median, na.rm = TRUE)
      } else {
        out[[cond]] <- rowMeans(mat, na.rm = TRUE)
      }
    }
  }
  out
}

make_atac_score_condition <- function(
    atac_score_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(atac_score_tbl)) .log_abort("`atac_score_tbl` must be a data.frame.")
  if (!"atac_peak" %in% names(atac_score_tbl)) .log_abort("`atac_score_tbl` needs an atac_peak column.")
  atac_tmp <- dplyr::rename(atac_score_tbl, peak_ID = .data$atac_peak)
  out <- make_fp_score_condition(
    fp_score_tbl = atac_tmp,
    metadata = metadata,
    label_col = label_col,
    agg_fun = agg_fun
  )
  dplyr::rename(out, atac_peak = .data$peak_ID)
}

make_fp_bound_condition <- function(
    fp_bound_tbl,
    fp_score_tbl,
    atac_overlap_tbl,
    fp_annotation_tbl,
    metadata,
    label_col,
    threshold_fp_score = 2,
    min_samples = 1L
) {
  if (!"peak_ID" %in% names(fp_bound_tbl)) .log_abort("`fp_bound_tbl` needs peak_ID.")
  if (!"peak_ID" %in% names(fp_score_tbl)) .log_abort("`fp_score_tbl` needs peak_ID.")
  if (!"atac_peak" %in% names(atac_overlap_tbl)) .log_abort("`atac_overlap_tbl` needs atac_peak.")
  if (!all(c("fp_peak", "atac_peak") %in% names(fp_annotation_tbl))) {
    .log_abort("`fp_annotation_tbl` needs fp_peak and atac_peak.")
  }
  if (!is.data.frame(metadata) || !"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  min_samples <- as.integer(min_samples)
  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 1L) {
    .log_abort("`min_samples` must be a positive integer.")
  }

  ids_bound <- setdiff(names(fp_bound_tbl), "peak_ID")
  ids_score <- setdiff(names(fp_score_tbl), "peak_ID")
  ids_atac <- setdiff(names(atac_overlap_tbl), "atac_peak")
  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_use <- intersect(meta_use$id, intersect(ids_bound, intersect(ids_score, ids_atac)))
  if (!length(ids_use)) {
    .log_abort("No overlapping ids across fp_bound, fp_score, atac_overlap, and metadata.")
  }
  meta_use <- meta_use[meta_use$id %in% ids_use, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")
  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  if (!identical(fp_bound_tbl$peak_ID, fp_score_tbl$peak_ID)) {
    idx <- match(fp_bound_tbl$peak_ID, fp_score_tbl$peak_ID)
    if (anyNA(idx)) .log_abort("fp_score rows do not align with fp_bound by peak_ID.")
    score_tbl <- fp_score_tbl[idx, , drop = FALSE]
  } else {
    score_tbl <- fp_score_tbl
  }

  bound_mat <- as.matrix(fp_bound_tbl[, ids_use, drop = FALSE])
  score_mat <- as.matrix(score_tbl[, ids_use, drop = FALSE])
  bound_mat[is.na(bound_mat)] <- 0L
  score_mat[is.na(score_mat)] <- NA_real_
  bound_mat <- bound_mat > 0 & score_mat >= threshold_fp_score
  storage.mode(bound_mat) <- "integer"

  map_fp_atac <- dplyr::distinct(fp_annotation_tbl, .data$fp_peak, .data$atac_peak)
  dup <- map_fp_atac |>
    dplyr::count(.data$fp_peak, name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup)) .log_abort("Each fp_peak must map to a single atac_peak.")

  idx_fp <- match(fp_bound_tbl$peak_ID, map_fp_atac$fp_peak)
  if (anyNA(idx_fp)) .log_abort("Some fp_bound peak_ID are not found in fp_annotation fp_peak.")
  atac_peaks <- map_fp_atac$atac_peak[idx_fp]
  idx_atac <- match(atac_peaks, atac_overlap_tbl$atac_peak)
  if (anyNA(idx_atac)) .log_abort("Mapped atac_peak is missing in atac_overlap.")
  atac_tbl <- atac_overlap_tbl[idx_atac, ids_use, drop = FALSE]
  atac_mat <- as.matrix(atac_tbl)
  atac_mat[is.na(atac_mat)] <- 0L
  atac_mat <- atac_mat > 0

  out <- tibble::tibble(peak_ID = fp_bound_tbl$peak_ID)
  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], ids_use)
      if (!length(ids_cond)) {
        .log_warn("No fp_bound columns found for condition {.val {cond}}.")
        out[[cond]] <- integer(nrow(fp_bound_tbl))
        next
      }
      id <- ids_cond[[1L]]
      out[[cond]] <- as.integer(bound_mat[, id, drop = FALSE] > 0 & atac_mat[, id, drop = FALSE] > 0)
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- intersect(cond_groups[[cond]], ids_use)
      if (!length(ids_cond)) {
        .log_warn("No fp_bound columns found for condition {.val {cond}}.")
        out[[cond]] <- integer(nrow(fp_bound_tbl))
        next
      }
      b_ok <- rowSums(bound_mat[, ids_cond, drop = FALSE], na.rm = TRUE) >= min_samples
      a_ok <- rowSums(atac_mat[, ids_cond, drop = FALSE], na.rm = TRUE) > 0L
      out[[cond]] <- as.integer(b_ok & a_ok)
    }
  }
  out
}

get_min_bound_peaks <- function(
    fp_bound_tbl,
    min_bound = 1L,
    samples = NULL,
    na_as_unbound = TRUE,
    use_parallel = FALSE,
    n_workers = getOption("mc.cores", 1L)
) {
  if (!is.data.frame(fp_bound_tbl)) .log_abort("`fp_bound_tbl` must be a data.frame.")
  if (!"peak_ID" %in% names(fp_bound_tbl)) .log_abort("`fp_bound_tbl` needs a peak_ID column.")
  min_bound <- as.integer(min_bound)
  if (length(min_bound) != 1L || is.na(min_bound) || min_bound < 1L) {
    .log_abort("`min_bound` must be a positive integer.")
  }

  all_ids <- setdiff(names(fp_bound_tbl), "peak_ID")
  use_ids <- if (is.null(samples)) all_ids else samples
  if (!is.null(samples)) {
    miss <- setdiff(samples, all_ids)
    if (length(miss)) {
      .log_abort("Some `samples` are not found in fp_bound: {paste(miss, collapse = ', ')}.")
    }
  }
  if (!length(use_ids)) .log_abort("No sample columns selected to evaluate binding.")

  mat_raw <- as.matrix(fp_bound_tbl[, use_ids, drop = FALSE])
  if (isTRUE(na_as_unbound)) {
    mat_raw[is.na(mat_raw)] <- 0L
  }
  mat <- mat_raw > 0

  if (isTRUE(use_parallel) && .Platform$OS.type == "windows") {
    .log_warn("Multicore is not supported on Windows; using single-core rowSums().")
    use_parallel <- FALSE
  }

  if (isTRUE(use_parallel) && n_workers > 1L) {
    n_rows <- nrow(mat)
    if (!n_rows) return(character(0))
    chunk_size <- ceiling(n_rows / n_workers)
    chunks <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))
    sums <- unlist(
      parallel::mclapply(
        chunks,
        function(idx) rowSums(mat[idx, , drop = FALSE]),
        mc.cores = n_workers
      ),
      use.names = FALSE
    )
  } else {
    sums <- rowSums(mat)
  }

  keep <- if (!isTRUE(na_as_unbound)) {
    any_na <- apply(is.na(mat_raw), 1L, any)
    !any_na & sums >= min_bound
  } else {
    sums >= min_bound
  }

  fp_bound_tbl$peak_ID[keep]
}

qn_footprints <- function(fp_tbl, id_col = "peak_ID", tie_method = c("average", "first")) {
  tie_method <- match.arg(tie_method)
  stopifnot(id_col %in% names(fp_tbl))

  samp <- setdiff(names(fp_tbl), id_col)
  if (!length(samp)) return(fp_tbl)

  ids_all <- fp_tbl[[id_col]]
  uniq_first <- !duplicated(ids_all)
  ids_unique <- ids_all[uniq_first]
  grp <- match(ids_all, ids_unique)

  x <- as.matrix(fp_tbl[uniq_first, samp, drop = FALSE])
  storage.mode(x) <- "double"
  n <- nrow(x)
  p <- ncol(x)
  if (n == 0L || p == 0L) return(fp_tbl)

  ref_sum <- numeric(n)
  ref_count <- integer(n)
  ord_list <- vector("list", p)
  m_list <- integer(p)

  for (j in seq_len(p)) {
    v <- x[, j]
    o <- order(v, na.last = TRUE)
    m <- sum(!is.na(v))
    ord_list[[j]] <- o
    m_list[j] <- m
    if (m > 0L) {
      ref_sum[seq_len(m)] <- ref_sum[seq_len(m)] + v[o][seq_len(m)]
      ref_count[seq_len(m)] <- ref_count[seq_len(m)] + 1L
    }
  }

  ref <- ref_sum
  nz <- ref_count > 0L
  ref[nz] <- ref_sum[nz] / ref_count[nz]
  ref[!nz] <- NA_real_

  x_qn <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    v <- x[, j]
    o <- ord_list[[j]]
    m <- m_list[j]
    if (m == 0L) next

    if (identical(tie_method, "first")) {
      x_qn[o[seq_len(m)], j] <- ref[seq_len(m)]
    } else {
      s <- v[o[seq_len(m)]]
      r <- rle(s)
      ends <- cumsum(r$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)
      for (k in seq_along(ends)) {
        idx <- starts[k]:ends[k]
        x_qn[o[idx], j] <- mean(ref[idx], na.rm = TRUE)
      }
    }
  }

  out <- fp_tbl
  for (j in seq_len(p)) {
    out[[samp[j]]] <- x_qn[, j][grp]
  }
  out
}

grn_add_rna_expressed <- function(
    grn_set,
    label_col,
    threshold_gene_expr = 10,
    min_samples = 1L,
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "rna_expressed") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$rna)) .log_abort("`grn_set$rna` is missing or invalid.")
  grn_set$rna_expressed <- make_rna_expressed(
    rna_tbl = grn_set$rna,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    threshold_gene_expr = threshold_gene_expr,
    min_samples = min_samples
  )
  grn_status_set(grn_set, "rna_expressed")
}

grn_add_rna_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "rna_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$rna)) .log_abort("`grn_set$rna` is missing or invalid.")
  grn_set$rna_condition <- make_rna_condition(
    rna_tbl = grn_set$rna,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  grn_status_set(grn_set, "rna_condition")
}

grn_add_fp_score_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_score_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_score)) .log_abort("`grn_set$fp_score` is missing or invalid.")
  grn_set$fp_score_condition <- make_fp_score_condition(
    fp_score_tbl = grn_set$fp_score,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  grn_status_set(grn_set, "fp_score_condition")
}

grn_add_atac_score_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "atac_score_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$atac_score)) .log_abort("`grn_set$atac_score` is missing or invalid.")
  grn_set$atac_score_condition <- make_atac_score_condition(
    atac_score_tbl = grn_set$atac_score,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  grn_status_set(grn_set, "atac_score_condition")
}

grn_add_fp_bound_condition <- function(
    grn_set,
    label_col,
    threshold_fp_score = 2,
    min_samples = 1L,
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_bound_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_bound)) .log_abort("`grn_set$fp_bound` is missing or invalid.")
  grn_set$fp_bound_condition <- make_fp_bound_condition(
    fp_bound_tbl = grn_set$fp_bound,
    fp_score_tbl = grn_set$fp_score,
    atac_overlap_tbl = grn_set$atac_overlap,
    fp_annotation_tbl = grn_set$fp_annotation,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    threshold_fp_score = threshold_fp_score,
    min_samples = min_samples
  )
  grn_status_set(grn_set, "fp_bound_condition")
}

grn_filter_fp_bound_condition <- function(
    grn_set,
    min_bound = 1L,
    use_parallel = TRUE,
    force = FALSE,
    verbose = FALSE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_filtered") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_bound_condition)) {
    .log_abort("`grn_set$fp_bound_condition` is missing; run grn_add_fp_bound_condition() first.")
  }
  keep_peaks <- get_min_bound_peaks(
    grn_set$fp_bound_condition,
    min_bound = min_bound,
    use_parallel = use_parallel
  )
  grn_set <- filter_fp_rows(
    grn_set = grn_set,
    peaks = keep_peaks,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose = FALSE
  )
  if (is.data.frame(grn_set$fp_score_condition)) {
    ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition$peak_ID)
    if (anyNA(ord)) .log_abort("fp_score_condition is missing peaks after filtering.")
    grn_set$fp_score_condition <- grn_set$fp_score_condition[ord, , drop = FALSE]
  }
  if (is.data.frame(grn_set$fp_bound_condition)) {
    ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_bound_condition$peak_ID)
    if (anyNA(ord)) .log_abort("fp_bound_condition is missing peaks after filtering.")
    grn_set$fp_bound_condition <- grn_set$fp_bound_condition[ord, , drop = FALSE]
  }
  grn_status_set(grn_set, "fp_filtered")
}

grn_add_fp_score_qn <- function(
    grn_set,
    id_col = "peak_ID",
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_score_qn") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_score_condition)) {
    .log_abort("`grn_set$fp_score_condition` is missing; run grn_add_fp_score_condition() first.")
  }
  grn_set$fp_score_condition_qn <- qn_footprints(grn_set$fp_score_condition, id_col = id_col)
  grn_status_set(grn_set, "fp_score_qn")
}
