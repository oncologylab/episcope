#' Prepare a GRN set for condition-aware linking
#'
#' Ensures the `grn_set` carries the right slots for `light_by_condition()`
#' (normalized signal tables and labeled metadata) and disambiguates any
#' duplicate sample labels by appending the sample id.
#'
#' @param grn_set List returned by `build_grn_set()` plus downstream slots.
#' @param label_col Character scalar naming the metadata column whose values
#'   should label conditions (default `"strict_match_rna"`).
#' @return `grn_set` clone ready for `light_by_condition()`.
#' @export
prepare_grn_set_for_light_by_condition <- function(grn_set, label_col = "strict_match_rna") {
  stopifnot(is.list(grn_set), !is.null(grn_set$sample_metadata_used))
  ds <- grn_set
  if (is.data.frame(grn_set$fp_score_condition_qn)) {
    ds$fp_score <- grn_set$fp_score_condition_qn
  }
  if (is.data.frame(grn_set$rna_condition)) {
    ds$rna <- grn_set$rna_condition
  }
  ds$sample_metadata_used <- ensure_grn_sample_metadata_label(grn_set, label_col)
  ds
}

ensure_grn_sample_metadata_label <- function(grn_set, label_col) {
  meta <- grn_set$sample_metadata_used
  stopifnot(is.data.frame(meta), "id" %in% names(meta))
  pick <- label_col
  if (!(pick %in% names(meta))) {
    pick <- "id"
  }
  label <- meta[[pick]]
  missing <- which(is.na(label) | label == "")
  if (length(missing)) {
    label[missing] <- meta$id[missing]
  }
  if (anyDuplicated(label)) {
    dup_idx <- which(duplicated(label) | duplicated(label, fromLast = TRUE))
    for (i in dup_idx) {
      label[i] <- paste0(label[i], "_", meta$id[i])
    }
  }
  # Keep sample id stable; only update the label column.
  meta[[pick]] <- label
  meta
}
