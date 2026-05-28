# File: utils_step1_tfbs_links.R
# Purpose: Internal Module 1 helpers for compact TFBS prediction tables.

.module1_parse_fp_coordinates <- function(fp_id) {
  fp_id <- as.character(fp_id)
  chr <- sub(":.*$", "", fp_id, perl = TRUE)
  pos <- sub("^[^:]+:", "", fp_id, perl = TRUE)
  start <- suppressWarnings(as.integer(sub("-.*$", "", pos, perl = TRUE)))
  end <- suppressWarnings(as.integer(sub("^.*-", "", pos, perl = TRUE)))
  tibble::tibble(chr = chr, start = start, end = end)
}

.module1_best_corr <- function(pearson_r, spearman_r) {
  pearson_r <- as.numeric(pearson_r)
  spearman_r <- as.numeric(spearman_r)
  p_use <- is.finite(pearson_r)
  s_use <- is.finite(spearman_r)
  best_r <- rep(NA_real_, length(pearson_r))
  best_method <- rep(NA_character_, length(pearson_r))

  use_p <- p_use & (!s_use | pearson_r >= spearman_r)
  use_s <- s_use & (!p_use | spearman_r > pearson_r)
  best_r[use_p] <- pearson_r[use_p]
  best_method[use_p] <- "pearson"
  best_r[use_s] <- spearman_r[use_s]
  best_method[use_s] <- "spearman"

  list(best_r = best_r, best_method = best_method)
}

.module1_tfbs_cutoffs <- function(r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  r_cutoff <- as.numeric(r_cutoff)[[1L]]
  if (!is.finite(r_cutoff)) .log_abort("`r_cutoff` must be finite.")
  p_cutoff <- if (is.null(p_cutoff) || length(p_cutoff) == 0L || is.na(p_cutoff[[1L]])) Inf else as.numeric(p_cutoff)[[1L]]
  fdr_cutoff <- if (is.null(fdr_cutoff) || length(fdr_cutoff) == 0L || is.na(fdr_cutoff[[1L]])) Inf else as.numeric(fdr_cutoff)[[1L]]
  if (!is.finite(p_cutoff) && !identical(p_cutoff, Inf)) .log_abort("`p_cutoff` must be finite or NULL.")
  if (!is.finite(fdr_cutoff) && !identical(fdr_cutoff, Inf)) .log_abort("`fdr_cutoff` must be finite or NULL.")
  list(r = r_cutoff, p = p_cutoff, fdr = fdr_cutoff)
}

.module1_apply_tfbs_cutoffs <- function(x, cutoffs) {
  if (!is.data.frame(x)) .log_abort("`x` must be a data.frame.")
  if (!all(c("best_r", "best_method") %in% names(x))) {
    .log_abort("`x` must include best_r and best_method.")
  }
  pearson_best <- rep(FALSE, nrow(x))
  if ("pearson_p" %in% names(x)) {
    pearson_best <- is.finite(x$pearson_r) & (x$best_method == "pearson")
  }
  spearman_best <- rep(FALSE, nrow(x))
  if ("spearman_p" %in% names(x)) {
    spearman_best <- is.finite(x$spearman_r) & (x$best_method == "spearman")
  }
  best_p <- rep(NA_real_, nrow(x))
  best_fdr <- rep(NA_real_, nrow(x))
  if ("pearson_p" %in% names(x)) best_p[pearson_best] <- x$pearson_p[pearson_best]
  if ("spearman_p" %in% names(x)) best_p[spearman_best] <- x$spearman_p[spearman_best]
  if ("pearson_p_adj" %in% names(x)) best_fdr[pearson_best] <- x$pearson_p_adj[pearson_best]
  if ("spearman_p_adj" %in% names(x)) best_fdr[spearman_best] <- x$spearman_p_adj[spearman_best]
  x$best_p <- best_p
  x$best_fdr <- best_fdr
  x$pass_r <- is.finite(x$best_r) & x$best_r >= cutoffs$r
  x$pass_p <- is.infinite(cutoffs$p) | (is.finite(x$best_p) & x$best_p <= cutoffs$p)
  x$pass_fdr <- is.infinite(cutoffs$fdr) | (is.finite(x$best_fdr) & x$best_fdr <= cutoffs$fdr)
  x$pass <- x$pass_r & x$pass_p & x$pass_fdr
  tibble::as_tibble(x)
}

.module1_merge_tfbs_stats <- function(pearson_stats, spearman_stats, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  keys <- c("fp_id", "atac_peak", "tf")
  if (!is.data.frame(pearson_stats)) .log_abort("`pearson_stats` must be a data.frame.")
  if (!is.data.frame(spearman_stats)) .log_abort("`spearman_stats` must be a data.frame.")
  if (!all(keys %in% names(pearson_stats))) .log_abort("`pearson_stats` is missing required key columns.")
  if (!all(keys %in% names(spearman_stats))) .log_abort("`spearman_stats` is missing required key columns.")

  pearson_keep <- unique(c(keys, "pearson_r", "pearson_p", "pearson_p_adj", "motifs"))
  spearman_keep <- unique(c(keys, "spearman_r", "spearman_p", "spearman_p_adj", "motifs"))
  pearson_keep <- intersect(pearson_keep, names(pearson_stats))
  spearman_keep <- intersect(spearman_keep, names(spearman_stats))

  merged <- dplyr::full_join(
    pearson_stats[, pearson_keep, drop = FALSE],
    spearman_stats[, spearman_keep, drop = FALSE],
    by = keys,
    suffix = c(".pearson", ".spearman")
  )

  if (!"pearson_r" %in% names(merged)) merged$pearson_r <- NA_real_
  if (!"spearman_r" %in% names(merged)) merged$spearman_r <- NA_real_
  best <- .module1_best_corr(merged$pearson_r, merged$spearman_r)
  merged$best_r <- best$best_r
  merged$best_method <- best$best_method
  .module1_apply_tfbs_cutoffs(
    merged,
    .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
  )
}

.module1_select_high_confidence_footprints <- function(motif_supported_correlations, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  if (!is.data.frame(motif_supported_correlations)) {
    .log_abort("`motif_supported_correlations` must be a data.frame.")
  }
  need <- c("fp_id", "atac_peak")
  if (!all(need %in% names(motif_supported_correlations))) {
    .log_abort("`motif_supported_correlations` must include fp_id and atac_peak.")
  }
  if (!"best_r" %in% names(motif_supported_correlations)) {
    if (!all(c("pearson_r", "spearman_r") %in% names(motif_supported_correlations))) {
      .log_abort("`motif_supported_correlations` must include best_r or Pearson/Spearman columns.")
    }
    best <- .module1_best_corr(
      motif_supported_correlations$pearson_r,
      motif_supported_correlations$spearman_r
    )
    motif_supported_correlations$best_r <- best$best_r
    motif_supported_correlations$best_method <- best$best_method
  }

  if (!"pass" %in% names(motif_supported_correlations)) {
    motif_supported_correlations <- .module1_apply_tfbs_cutoffs(
      motif_supported_correlations,
      .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
    )
  }
  keep <- motif_supported_correlations$pass %in% TRUE
  passed <- motif_supported_correlations[keep, , drop = FALSE]
  if (!nrow(passed)) {
    return(tibble::tibble(
      fp_id = character(), chr = character(), start = integer(), end = integer(),
      atac_peak = character(), n_canonical_bound_tfs = integer(),
      canonical_bound_tfs = character(), canonical_bound_motifs = character(),
      has_canonical_bound = logical()
    ))
  }
  motif_cols <- intersect(c("motifs", "motifs.pearson", "motifs.spearman"), names(passed))
  if (length(motif_cols)) {
    passed$motif_use <- passed[[motif_cols[[1L]]]]
  } else {
    passed$motif_use <- NA_character_
  }
  dt <- data.table::as.data.table(passed)
  out <- dt[, .(
    atac_peak = as.character(atac_peak[[1L]]),
    n_canonical_bound_tfs = data.table::uniqueN(tf),
    canonical_bound_tfs = paste(sort(unique(as.character(tf))), collapse = ";"),
    canonical_bound_motifs = paste(sort(unique(as.character(motif_use[!is.na(motif_use) & nzchar(motif_use)]))), collapse = ";"),
    has_canonical_bound = TRUE
  ), by = .(fp_id)]
  coords <- .module1_parse_fp_coordinates(out$fp_id)
  out_tbl <- tibble::as_tibble(out)
  tibble::as_tibble(cbind(out_tbl["fp_id"], coords, out_tbl[c("atac_peak", "n_canonical_bound_tfs", "canonical_bound_tfs", "canonical_bound_motifs", "has_canonical_bound")]))
}

.module1_condition_support <- function(fp_id, tf, omics_data) {
  if (!is.list(omics_data)) return(rep(NA_integer_, length(fp_id)))
  if (!is.data.frame(omics_data$fp_bound_condition) || !is.data.frame(omics_data$rna_expressed)) {
    return(rep(NA_integer_, length(fp_id)))
  }
  fp_bound <- omics_data$fp_bound_condition
  rna_expr <- omics_data$rna_expressed
  cond_cols <- intersect(
    setdiff(names(fp_bound), "peak_ID"),
    setdiff(names(rna_expr), c("ensembl_gene_id", "HGNC"))
  )
  if (!length(cond_cols)) return(rep(NA_integer_, length(fp_id)))

  fp_idx <- match(fp_id, fp_bound$peak_ID)
  rna_expr <- rna_expr[!duplicated(toupper(rna_expr$HGNC)), , drop = FALSE]
  tf_idx <- match(toupper(tf), toupper(rna_expr$HGNC))
  out <- rep(0L, length(fp_id))
  valid <- !is.na(fp_idx) & !is.na(tf_idx)
  if (!any(valid)) return(out)

  fp_mat <- as.matrix(fp_bound[, cond_cols, drop = FALSE]) > 0L
  tf_mat <- as.matrix(rna_expr[, cond_cols, drop = FALSE]) > 0L
  fp_mat[is.na(fp_mat)] <- FALSE
  tf_mat[is.na(tf_mat)] <- FALSE

  valid_idx <- which(valid)
  chunk_size <- 1000000L
  chunks <- split(valid_idx, ceiling(seq_along(valid_idx) / chunk_size))
  for (idx in chunks) {
    out[idx] <- rowSums(fp_mat[fp_idx[idx], , drop = FALSE] & tf_mat[tf_idx[idx], , drop = FALSE])
  }
  out
}

.module1_build_tfbs_links <- function(tfbs_stats, high_confidence_footprints, omics_data = NULL) {
  if (!is.data.frame(tfbs_stats)) .log_abort("`tfbs_stats` must be a data.frame.")
  if (!is.data.frame(high_confidence_footprints)) .log_abort("`high_confidence_footprints` must be a data.frame.")
  need_stats <- c("fp_id", "atac_peak", "tf", "best_r", "best_method", "pass")
  need_high <- c("fp_id", "chr", "start", "end", "atac_peak")
  if (!all(need_stats %in% names(tfbs_stats))) .log_abort("`tfbs_stats` is missing required columns.")
  if (!all(need_high %in% names(high_confidence_footprints))) .log_abort("`high_confidence_footprints` is missing required columns.")

  high <- high_confidence_footprints[!duplicated(high_confidence_footprints$fp_id), need_high, drop = FALSE]
  stats_pass <- tfbs_stats[tfbs_stats$pass %in% TRUE, , drop = FALSE]
  if (!nrow(stats_pass)) {
    return(tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      best_r = numeric(),
      best_method = character(),
      condition_support = integer()
    ))
  }

  high_idx <- match(stats_pass$fp_id, high$fp_id)
  keep <- !is.na(high_idx)
  stats_pass <- stats_pass[keep, , drop = FALSE]
  high_idx <- high_idx[keep]
  if (!nrow(stats_pass)) {
    return(tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      best_r = numeric(),
      best_method = character(),
      condition_support = integer()
    ))
  }

  out <- data.table::as.data.table(stats_pass)
  out[, chr := as.character(high$chr[high_idx])]
  out[, start := as.integer(high$start[high_idx])]
  out[, end := as.integer(high$end[high_idx])]
  out[, condition_support := .module1_condition_support(fp_id, tf, omics_data)]
  keep_cols <- c(
    "fp_id", "chr", "start", "end", "atac_peak", "tf",
    intersect(c("motifs", "pearson_r", "spearman_r", "pearson_p", "spearman_p"), names(out)),
    "best_r", "best_method", "condition_support"
  )
  out <- out[, ..keep_cols]
  data.table::setorder(out, fp_id, tf)
  tibble::as_tibble(out)
}
