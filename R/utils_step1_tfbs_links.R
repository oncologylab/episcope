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

.module1_merge_tfbs_stats <- function(pearson_stats, spearman_stats, r_cutoff = 0.3) {
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
  merged$pass <- is.finite(merged$best_r) & merged$best_r >= as.numeric(r_cutoff)
  tibble::as_tibble(merged)
}

.module1_select_high_confidence_footprints <- function(motif_supported_correlations, r_cutoff = 0.3) {
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
  }

  keep <- is.finite(motif_supported_correlations$best_r) &
    motif_supported_correlations$best_r >= as.numeric(r_cutoff)
  out <- motif_supported_correlations[keep, need, drop = FALSE]
  out <- out[!duplicated(out$fp_id), , drop = FALSE]
  coords <- .module1_parse_fp_coordinates(out$fp_id)
  tibble::as_tibble(cbind(out["fp_id"], coords, out["atac_peak"]))
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
  tf_idx <- match(toupper(tf), toupper(rna_expr$HGNC))
  out <- rep(0L, length(fp_id))
  for (i in seq_along(fp_id)) {
    if (is.na(fp_idx[[i]]) || is.na(tf_idx[[i]])) next
    fp_vec <- as.integer(unlist(fp_bound[fp_idx[[i]], cond_cols, drop = FALSE], use.names = FALSE))
    tf_vec <- as.integer(unlist(rna_expr[tf_idx[[i]], cond_cols, drop = FALSE], use.names = FALSE))
    fp_vec[is.na(fp_vec)] <- 0L
    tf_vec[is.na(tf_vec)] <- 0L
    out[[i]] <- sum(fp_vec > 0L & tf_vec > 0L)
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
  stats_pass <- tfbs_stats[isTRUE(tfbs_stats$pass) | tfbs_stats$pass %in% TRUE, , drop = FALSE]
  stats_pass <- stats_pass[stats_pass$fp_id %in% high$fp_id, , drop = FALSE]
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

  out <- dplyr::left_join(stats_pass, high, by = c("fp_id", "atac_peak"))
  out$condition_support <- .module1_condition_support(out$fp_id, out$tf, omics_data)
  keep <- c(
    "fp_id", "chr", "start", "end", "atac_peak", "tf",
    intersect(c("motifs", "pearson_r", "spearman_r", "pearson_p", "spearman_p"), names(out)),
    "best_r", "best_method", "condition_support"
  )
  out <- out[, keep, drop = FALSE]
  out <- out[order(out$fp_id, out$tf), , drop = FALSE]
  tibble::as_tibble(out)
}
