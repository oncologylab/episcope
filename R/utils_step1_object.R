# File: utils_step1_object.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide helpers for constructing and validating the core Module 1 multi-omic
# object.
#
# Inputs:
# - loaded ATAC, RNA, metadata, and footprint tables
# - motif database information
#
# Outputs:
# - a Step 1 multi-omic object ready for downstream condition-level processing
#
# Notes:
# - Keep object construction separate from filtering and QC.

#' Module 1 object helpers
#'
#' @noRd
NULL

grn_status_init <- function(grn_set) {
  if (!is.list(grn_set)) {
    .log_abort("`grn_set` must be a list.")
  }
  if (is.null(grn_set$status) || !is.list(grn_set$status)) {
    grn_set$status <- list()
  }
  grn_set
}

grn_status_is <- function(grn_set, key) {
  is.list(grn_set) && is.list(grn_set$status) && isTRUE(grn_set$status[[key]])
}

grn_status_set <- function(grn_set, key, value = TRUE) {
  grn_set <- grn_status_init(grn_set)
  grn_set$status[[key]] <- isTRUE(value)
  grn_set
}

build_grn_set <- function(
    fp_score,
    fp_bound,
    fp_annotation,
    atac_score,
    atac_overlap,
    rna,
    metadata,
    tf_list,
    motif_db = NULL,
    label_col = NULL,
    expected_n = NULL,
    filter_atac_by_fp_annotation = TRUE
) {
  if (!"peak_ID" %in% names(fp_score)) .log_abort("`fp_score` needs 'peak_ID'.")
  if (!"peak_ID" %in% names(fp_bound)) .log_abort("`fp_bound` needs 'peak_ID'.")
  if (!"atac_peak" %in% names(atac_score)) .log_abort("`atac_score` needs 'atac_peak'.")
  if (!"atac_peak" %in% names(atac_overlap)) .log_abort("`atac_overlap` needs 'atac_peak'.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna))) {
    .log_abort("`rna` needs 'ensembl_gene_id' and 'HGNC'.")
  }
  if (!"id" %in% names(metadata)) .log_abort("`metadata` needs an 'id' column.")
  if (!is.null(label_col) && !label_col %in% names(metadata)) {
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  }
  if (!is.character(tf_list) && !is.null(tf_list)) {
    .log_abort("`tf_list` must be a character vector or NULL.")
  }

  meta_use <- metadata[!is.na(metadata$id), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  if (!nrow(meta_use)) .log_abort("No usable sample ids found in `metadata$id`.")

  id_col_fp <- if ("id_fp" %in% names(meta_use)) "id_fp" else "id"
  id_col_rna <- if ("id_rna" %in% names(meta_use)) "id_rna" else "id"
  id_col_atac <- if ("id_atac" %in% names(meta_use)) "id_atac" else id_col_fp

  fp_ids <- as.character(meta_use[[id_col_fp]])
  rna_ids <- as.character(meta_use[[id_col_rna]])
  atac_ids <- as.character(meta_use[[id_col_atac]])

  miss_fp_score <- meta_use$id[is.na(fp_ids) | !(fp_ids %in% names(fp_score))]
  miss_fp_bound <- meta_use$id[is.na(fp_ids) | !(fp_ids %in% names(fp_bound))]
  miss_rna <- meta_use$id[is.na(rna_ids) | !(rna_ids %in% names(rna))]
  miss_atac_sc <- meta_use$id[is.na(atac_ids) | !(atac_ids %in% names(atac_score))]
  miss_atac_ol <- meta_use$id[is.na(atac_ids) | !(atac_ids %in% names(atac_overlap))]

  if (length(miss_fp_score)) .log_warn("Dropping {length(miss_fp_score)} id(s) missing in `fp_score` via {id_col_fp}.")
  if (length(miss_fp_bound)) .log_warn("Dropping {length(miss_fp_bound)} id(s) missing in `fp_bound` via {id_col_fp}.")
  if (length(miss_rna)) .log_warn("Dropping {length(miss_rna)} id(s) missing in `rna` via {id_col_rna}.")
  if (length(miss_atac_sc)) .log_warn("`atac_score` is missing {length(miss_atac_sc)} id(s) via {id_col_atac}.")
  if (length(miss_atac_ol)) .log_warn("`atac_overlap` is missing {length(miss_atac_ol)} id(s) via {id_col_atac}.")

  keep_idx <- which(
    !is.na(fp_ids) & (fp_ids %in% names(fp_score)) &
      (fp_ids %in% names(fp_bound)) &
      !is.na(rna_ids) & (rna_ids %in% names(rna))
  )
  if (!length(keep_idx)) {
    .log_abort("No overlapping ids across fp_score, fp_bound, and rna after mapping checks.")
  }

  keep_ids <- as.character(meta_use$id[keep_idx])
  keep_fp_ids <- fp_ids[keep_idx]
  keep_rna_ids <- rna_ids[keep_idx]
  keep_atac_ids <- atac_ids[keep_idx]

  if (!is.null(expected_n) && length(keep_ids) != expected_n) {
    .log_warn("Aligned sample count is {length(keep_ids)} but expected {expected_n}.")
  }

  if (isTRUE(filter_atac_by_fp_annotation)) {
    peaks_keep_atac <- unique(fp_annotation$atac_peak)
    peaks_keep_atac <- peaks_keep_atac[!is.na(peaks_keep_atac)]
    if (!length(peaks_keep_atac)) {
      .log_abort("No usable 'atac_peak' values found in `fp_annotation` after filtering.")
    }
    atac_score <- dplyr::semi_join(
      atac_score,
      tibble::tibble(atac_peak = peaks_keep_atac),
      by = "atac_peak"
    )
    atac_overlap <- dplyr::semi_join(
      atac_overlap,
      tibble::tibble(atac_peak = peaks_keep_atac),
      by = "atac_peak"
    )
  }

  fp_score_sub <- fp_score[, c("peak_ID", keep_fp_ids), drop = FALSE]
  names(fp_score_sub)[-1] <- keep_ids
  fp_score_sub <- tibble::as_tibble(fp_score_sub)

  fp_bound_sub <- fp_bound[, c("peak_ID", keep_fp_ids), drop = FALSE]
  names(fp_bound_sub)[-1] <- keep_ids
  fp_bound_sub <- tibble::as_tibble(fp_bound_sub)

  atac_keep_ok <- !is.na(keep_atac_ids) &
    (keep_atac_ids %in% names(atac_score)) &
    (keep_atac_ids %in% names(atac_overlap))
  if (!all(atac_keep_ok)) {
    bad <- keep_ids[!atac_keep_ok]
    if (length(bad)) {
      .log_warn("Dropping {length(bad)} id(s) missing ATAC source columns after mapping.")
    }
    keep_ids <- keep_ids[atac_keep_ok]
    keep_fp_ids <- keep_fp_ids[atac_keep_ok]
    keep_rna_ids <- keep_rna_ids[atac_keep_ok]
    keep_atac_ids <- keep_atac_ids[atac_keep_ok]
    fp_score_sub <- fp_score_sub[, c("peak_ID", keep_ids), drop = FALSE]
    fp_bound_sub <- fp_bound_sub[, c("peak_ID", keep_ids), drop = FALSE]
  }
  if (!length(keep_ids)) {
    .log_abort("No overlapping ids across fp, rna, and atac after mapping checks.")
  }

  atac_score_sub <- atac_score[, c("atac_peak", keep_atac_ids), drop = FALSE]
  names(atac_score_sub)[-1] <- keep_ids
  atac_score_sub <- tibble::as_tibble(atac_score_sub)

  atac_overlap_sub <- atac_overlap[, c("atac_peak", keep_atac_ids), drop = FALSE]
  names(atac_overlap_sub)[-1] <- keep_ids
  atac_overlap_sub <- tibble::as_tibble(atac_overlap_sub)

  rna_sub <- rna[, c("ensembl_gene_id", "HGNC", keep_rna_ids), drop = FALSE]
  names(rna_sub)[-(1:2)] <- keep_ids
  rna_sub <- tibble::as_tibble(rna_sub)

  meta_used <- meta_use[match(keep_ids, meta_use$id), , drop = FALSE]
  sample_labels <- if (!is.null(label_col)) meta_used[[label_col]] else NULL

  grn_set <- list(
    fp_score = fp_score_sub,
    fp_bound = fp_bound_sub,
    atac_score = atac_score_sub,
    atac_overlap = atac_overlap_sub,
    rna = rna_sub,
    fp_annotation = fp_annotation,
    tf_list = tf_list,
    motif_db = motif_db,
    sample_metadata_used = meta_used,
    sample_labels = sample_labels,
    dropped_ids = setdiff(as.character(meta_use$id), keep_ids)
  )

  grn_set <- grn_status_init(grn_set)
  grn_status_set(grn_set, "built")
}

filter_fp_rows <- function(
    grn_set,
    peaks = NULL,
    predicate = NULL,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose = TRUE,
    warn_on_mismatch = TRUE
) {
  stopifnot(
    is.list(grn_set),
    all(c("fp_score", "fp_bound", "fp_annotation") %in% names(grn_set))
  )

  fp_score <- grn_set$fp_score
  fp_bound <- grn_set$fp_bound
  fp_annot <- grn_set$fp_annotation

  s_ids <- fp_score[[score_key]]
  b_ids <- fp_bound[[bound_key]]
  a_ids <- fp_annot[[annot_key]]

  if (!is.character(s_ids) || !is.character(b_ids) || !is.character(a_ids)) {
    .log_abort("Keys must be character; coerce before calling.")
  }

  if (length(s_ids) == 0L && length(b_ids) == 0L && length(a_ids) == 0L) {
    if (!is.null(predicate)) {
      keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
      if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
        .log_abort("`predicate` must return a logical vector of length nrow(fp_score).")
      }
    }
    if (!is.null(peaks) && length(peaks) > 0L) {
      .log_warn("Requested {length(peaks)} peak ID(s), but the current set is empty; returning empty triplet.")
    }
    grn_set$fp_score <- fp_score[0, , drop = FALSE]
    grn_set$fp_bound <- fp_bound[0, , drop = FALSE]
    grn_set$fp_annotation <- fp_annot[0, , drop = FALSE]
    return(grn_set)
  }

  common_ids <- Reduce(intersect, list(unique(s_ids), unique(b_ids), unique(a_ids)))
  if (!length(common_ids)) {
    .log_abort("No common keys across fp_score, fp_bound, and fp_annotation.")
  }

  if (!(setequal(s_ids, b_ids) && setequal(s_ids, a_ids))) {
    n_drop <- c(
      score_only = sum(!(s_ids %in% common_ids)),
      bound_only = sum(!(b_ids %in% common_ids)),
      annot_only = sum(!(a_ids %in% common_ids))
    )
    if (isTRUE(warn_on_mismatch)) {
      .log_warn(
        c(
          "!" = "Key sets differ across tables; reconciling to the intersection ({length(common_ids)} rows).",
          "i" = "Dropped score-only={n_drop['score_only']}, bound-only={n_drop['bound_only']}, annot-only={n_drop['annot_only']}."
        )
      )
    }

    a_keep <- fp_annot[[annot_key]] %in% common_ids
    a_ord <- fp_annot[[annot_key]][a_keep]

    s_keep <- fp_score[[score_key]] %in% common_ids
    b_keep <- fp_bound[[bound_key]] %in% common_ids

    fp_annot <- fp_annot[a_keep, , drop = FALSE]

    fp_score <- fp_score[s_keep, , drop = FALSE]
    idx_s <- match(a_ord, fp_score[[score_key]])
    fp_score <- fp_score[idx_s, , drop = FALSE]

    fp_bound <- fp_bound[b_keep, , drop = FALSE]
    idx_b <- match(a_ord, fp_bound[[bound_key]])
    fp_bound <- fp_bound[idx_b, , drop = FALSE]

    s_ids <- fp_score[[score_key]]
    b_ids <- fp_bound[[bound_key]]
    a_ids <- fp_annot[[annot_key]]
  }

  if (!is.null(predicate)) {
    keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
    if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
      .log_abort("`predicate` must return a logical vector of length nrow(fp_score).")
    }
    keep_ids <- s_ids[keep_lgl %in% TRUE]
  } else if (!is.null(peaks)) {
    if (!all(peaks %in% s_ids)) {
      miss <- setdiff(peaks, s_ids)
      if (length(miss) && isTRUE(warn_on_mismatch)) {
        .log_warn("Some requested `peaks` are not present after reconciliation; dropping {length(miss)} missing id(s).")
      }
    }
    peaks_uniq <- peaks[!duplicated(peaks)]
    keep_ids <- peaks_uniq[peaks_uniq %in% s_ids]
  } else {
    grn_set$fp_score <- fp_score
    grn_set$fp_bound <- fp_bound
    grn_set$fp_annotation <- fp_annot
    return(grn_set)
  }

  if (!length(keep_ids)) {
    fp_score2 <- fp_score[0, , drop = FALSE]
    fp_bound2 <- fp_bound[0, , drop = FALSE]
    fp_annot2 <- fp_annot[0, , drop = FALSE]
  } else {
    keep_s <- fp_score[[score_key]] %in% keep_ids
    keep_b <- fp_bound[[bound_key]] %in% keep_ids
    keep_a <- fp_annot[[annot_key]] %in% keep_ids

    fp_score2 <- fp_score[keep_s, , drop = FALSE]
    fp_bound2 <- fp_bound[keep_b, , drop = FALSE]
    fp_annot2 <- fp_annot[keep_a, , drop = FALSE]

    ord_a <- fp_annot2[[annot_key]]
    if (!is.null(peaks)) {
      target <- keep_ids
      ord_idx <- match(ord_a, target)
      o <- order(ord_idx, seq_along(ord_idx), na.last = TRUE)
      fp_annot2 <- fp_annot2[o, , drop = FALSE]
      m_s <- match(fp_annot2[[annot_key]], fp_score2[[score_key]])
      m_b <- match(fp_annot2[[annot_key]], fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    } else {
      m_s <- match(ord_a, fp_score2[[score_key]])
      m_b <- match(ord_a, fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    }
  }

  s2 <- fp_score2[[score_key]]
  b2 <- fp_bound2[[bound_key]]
  a2 <- fp_annot2[[annot_key]]
  if (!(identical(s2, a2) && identical(b2, a2))) {
    .log_abort("Filtered key sets are not identical across fp_score, fp_bound, and fp_annotation.")
  }

  grn_set$fp_score <- fp_score2
  grn_set$fp_bound <- fp_bound2
  grn_set$fp_annotation <- fp_annot2
  grn_set
}
