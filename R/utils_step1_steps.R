# File: utils_step1_steps.R
# Purpose: Public step-by-step helpers for Module 1 TFBS prediction.

#' Prepare Module 1 TFBS prediction inputs
#'
#' @param omics_data CraftGRN multiomic object returned by `load_prep_multiomic_data()`.
#' @param label_col Optional metadata column used to rebuild condition matrices.
#' @param tf_subset Optional TF symbols to keep.
#' @param verbose Emit concise progress messages.
#' @return A list containing prepared data, condition columns, TFs, and footprint universe.
#' @export
module1_prepare_tfbs_inputs <- function(omics_data, label_col = NULL, tf_subset = NULL, verbose = TRUE) {
  omics_data <- .module1_prepare_predict_omics(omics_data, label_col = label_col, verbose = verbose)
  cond_cols <- .module1_condition_columns(omics_data)
  expressed_tfs <- .module1_expressed_tfs(omics_data, tf_subset = tf_subset)
  fp_input <- unique(as.character(omics_data$fp_score_condition_qn$peak_ID))
  fp_input <- fp_input[!is.na(fp_input) & nzchar(fp_input)]
  bound_fp_ids <- fp_input
  if (is.data.frame(omics_data$fp_bound_condition) && "peak_ID" %in% names(omics_data$fp_bound_condition)) {
    bound_cols <- intersect(setdiff(names(omics_data$fp_bound_condition), "peak_ID"), cond_cols)
    if (length(bound_cols)) {
      bound_mat <- as.matrix(omics_data$fp_bound_condition[, bound_cols, drop = FALSE])
      bound_mat[is.na(bound_mat)] <- 0L
      bound_fp_ids <- as.character(omics_data$fp_bound_condition$peak_ID[rowSums(bound_mat > 0, na.rm = TRUE) > 0])
    }
  }
  fp_universe <- intersect(fp_input, bound_fp_ids)
  if (isTRUE(verbose)) .log_inform("Module 1 inputs: {length(fp_universe)} usable FP(s), {length(bound_fp_ids)} bound FP(s), {length(expressed_tfs)} expressed TF(s).")
  list(omics_data = omics_data, condition_cols = cond_cols, expressed_tfs = expressed_tfs, fp_input = fp_input, bound_fp_ids = bound_fp_ids, fp_universe = fp_universe)
}

#' Correlate motif-supported FP-TF pairs
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional best-method p-value cutoff.
#' @param fdr_cutoff Optional best-method FDR cutoff.
#' @param min_non_na Minimum finite condition pairs required.
#' @param cores Number of worker cores; NULL uses all available cores.
#' @param verbose Emit concise progress messages.
#' @return A tibble with Pearson, Spearman, best-method statistics, and pass flags.
#' @noRd
module1_correlate_motif_supported_tfbs <- function(module1_inputs, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, min_non_na = 3L, cores = NULL, verbose = TRUE) {
  omics_data <- module1_inputs$omics_data
  pairs <- .module1_motif_supported_pairs(omics_data, tf_subset = module1_inputs$expressed_tfs)
  pairs <- pairs[pairs$fp_id %in% module1_inputs$fp_universe, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 1 motif-supported correlation: testing {nrow(pairs)} FP-TF pair(s).")
  raw <- .module1_compute_pair_correlations(omics_data = omics_data, pairs = pairs, min_non_na = min_non_na, cores = cores)
  .module1_merge_correlation_stats(raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE], raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE], r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
}

#' Correlate TFs to their canonical TFBS
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional best-method p-value cutoff.
#' @param fdr_cutoff Optional best-method FDR cutoff.
#' @param min_non_na Minimum finite condition pairs required.
#' @param cores Number of worker cores; NULL uses all available cores.
#' @param verbose Emit concise progress messages.
#' @return A tibble with Pearson, Spearman, best-method statistics, and pass flags.
#' @export
module1_correlate_TF_to_canonical_tfbs <- function(module1_inputs, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, min_non_na = 3L, cores = NULL, verbose = TRUE) {
  module1_correlate_motif_supported_tfbs(
    module1_inputs = module1_inputs,
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff,
    min_non_na = min_non_na,
    cores = cores,
    verbose = verbose
  )
}

#' Select footprints for Module 1 TFBS prediction
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param motif_supported_correlations Output from module1_correlate_TF_to_canonical_tfbs.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional p-value cutoff.
#' @param fdr_cutoff Optional FDR cutoff.
#' @param filter_to_canonical_bound Keep only footprints with a passing motif-supported TF.
#' @param verbose Emit concise progress messages.
#' @return A list with high-confidence and prediction footprint tables.
#' @noRd
module1_select_prediction_footprints <- function(module1_inputs, motif_supported_correlations, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, filter_to_canonical_bound = TRUE, verbose = TRUE) {
  high <- .module1_select_high_confidence_footprints(motif_supported_correlations, r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
  canonical_fp_ids <- unique(high$fp_id)
  if (isTRUE(filter_to_canonical_bound)) {
    pred <- high
  } else {
    ann <- module1_inputs$omics_data$fp_annotation[!duplicated(module1_inputs$omics_data$fp_annotation$fp_peak), c("fp_peak", "atac_peak"), drop = FALSE]
    ann <- ann[ann$fp_peak %in% module1_inputs$fp_universe, , drop = FALSE]
    coords <- .module1_parse_fp_coordinates(ann$fp_peak)
    pred <- tibble::as_tibble(cbind(tibble::tibble(fp_id = as.character(ann$fp_peak)), coords, tibble::tibble(atac_peak = as.character(ann$atac_peak))))
    pred <- dplyr::left_join(pred, high[, c("fp_id", "n_canonical_bound_tfs", "canonical_bound_tfs", "canonical_bound_motifs"), drop = FALSE], by = "fp_id")
    pred$n_canonical_bound_tfs[is.na(pred$n_canonical_bound_tfs)] <- 0L
    pred$canonical_bound_tfs[is.na(pred$canonical_bound_tfs)] <- ""
    pred$canonical_bound_motifs[is.na(pred$canonical_bound_motifs)] <- ""
    pred$has_canonical_bound <- pred$fp_id %in% canonical_fp_ids
  }
  if (isTRUE(verbose)) .log_inform("Module 1 footprint selection: {length(canonical_fp_ids)} high-confidence FP(s), {nrow(pred)} prediction FP(s).")
  list(high_confidence_footprints = high, prediction_footprints = pred, canonical_fp_ids = canonical_fp_ids, n_removed = length(setdiff(module1_inputs$fp_universe, canonical_fp_ids)))
}

#' Filter footprints with canonical binding for full TFBS prediction
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param motif_supported_correlations Output from module1_correlate_TF_to_canonical_tfbs.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional p-value cutoff.
#' @param fdr_cutoff Optional FDR cutoff.
#' @param filter_to_canonical_bound Keep only footprints with a passing motif-supported TF.
#' @param verbose Emit concise progress messages.
#' @return A list with canonical-bound and prediction footprint tables.
#' @export
module1_filter_canonical_bound_tfbs <- function(module1_inputs, motif_supported_correlations, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, filter_to_canonical_bound = TRUE, verbose = TRUE) {
  module1_select_prediction_footprints(
    module1_inputs = module1_inputs,
    motif_supported_correlations = motif_supported_correlations,
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff,
    filter_to_canonical_bound = filter_to_canonical_bound,
    verbose = verbose
  )
}

#' Predict TFBS from selected Module 1 footprints
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param prediction_footprints Footprint table from module1_filter_canonical_bound_tfbs.
#' @param out_dir Optional output directory. Required when `write_outputs = TRUE`.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional best-method p-value cutoff.
#' @param fdr_cutoff Optional best-method FDR cutoff.
#' @param min_non_na Minimum finite condition pairs required.
#' @param cores Number of worker cores; NULL uses all available cores.
#' @param write_outputs Write predicted TFBS outputs.
#' @param output_format One of csv, parquet, or auto.
#' @param return_prediction_stats Return full prediction statistics in memory.
#' @param verbose Emit concise progress messages.
#' @return A list with prediction statistics or manifests and predicted_tfbs.
#' @noRd
module1_predict_tfbs_from_correlations <- function(module1_inputs, prediction_footprints, out_dir = NULL, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, min_non_na = 3L, cores = NULL, write_outputs = FALSE, output_format = c("csv", "parquet", "auto"), return_prediction_stats = NULL, verbose = TRUE) {
  stopifnot(is.logical(write_outputs), length(write_outputs) == 1L, !is.na(write_outputs))
  if (isTRUE(write_outputs) && (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir))) {
    .log_abort("`out_dir` must be a non-empty path when `write_outputs = TRUE`.")
  }
  output_format <- .module1_output_format(output_format)
  pair_count <- as.double(nrow(prediction_footprints)) * as.double(length(module1_inputs$expressed_tfs))
  keep_stats <- if (is.null(return_prediction_stats)) !isTRUE(write_outputs) || pair_count <= getOption("craftgrn.module1_prediction_return_limit", 5000000) else isTRUE(return_prediction_stats)
  stats_dir <- if (isTRUE(write_outputs) && !isTRUE(keep_stats)) file.path(out_dir, "module1_prediction_stats_chunks") else NULL
  streamed <- .module1_predict_tfbs_streamed(module1_inputs$omics_data, prediction_footprints, r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff, tf_subset = module1_inputs$expressed_tfs, min_non_na = min_non_na, cores = cores, return_prediction_stats = keep_stats, prediction_stats_out_dir = stats_dir, output_format = output_format, verbose = verbose)
  if (!is.null(streamed$prediction_stats_manifest_path) && !isTRUE(keep_stats)) {
    pred_paths <- .write_predicted_tfbs_from_prediction_stats_manifest(streamed$prediction_stats_manifest, out_dir = out_dir, output_format = output_format)
    predicted_tfbs <- tibble::tibble(tfbs_id = character(), fp_id = character(), chr = character(), start = integer(), end = integer(), atac_peak = character(), tf = character(), condition_support = integer())
  } else {
    predicted_tfbs <- build_predicted_tfbs(streamed$prediction_stats)
    pred_paths <- if (isTRUE(write_outputs)) .write_predicted_tfbs_table(predicted_tfbs, out_dir = out_dir, output_format = output_format) else NULL
  }
  list(prediction_stats = streamed$prediction_stats, prediction_stats_manifest = streamed$prediction_stats_manifest, prediction_stats_manifest_path = streamed$prediction_stats_manifest_path, predicted_tfbs = predicted_tfbs, predicted_tfbs_paths = pred_paths, prediction_pair_count = streamed$prediction_pair_count, n_prediction_stats = streamed$n_prediction_stats)
}

#' Predict full TFBS for all expressed TFs
#'
#' @param module1_inputs Output from module1_prepare_tfbs_inputs.
#' @param prediction_footprints Footprint table from module1_filter_canonical_bound_tfbs.
#' @param out_dir Optional output directory. Required when `write_outputs = TRUE`.
#' @param r_cutoff Minimum positive best correlation.
#' @param p_cutoff Optional best-method p-value cutoff.
#' @param fdr_cutoff Optional best-method FDR cutoff.
#' @param min_non_na Minimum finite condition pairs required.
#' @param cores Number of worker cores; NULL uses all available cores.
#' @param write_outputs Write predicted TFBS outputs. Defaults to `FALSE` so
#'   calls do not write into the working directory unless an output directory is
#'   explicitly supplied.
#' @param output_format One of csv, parquet, or auto.
#' @param return_prediction_stats Return full prediction statistics in memory.
#' @param verbose Emit concise progress messages.
#' @return A list with prediction statistics or manifests and predicted TFBS outputs.
#' @export
module1_predict_full_tfbs <- function(module1_inputs, prediction_footprints, out_dir = NULL, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL, min_non_na = 3L, cores = NULL, write_outputs = FALSE, output_format = c("csv", "parquet", "auto"), return_prediction_stats = NULL, verbose = TRUE) {
  module1_predict_tfbs_from_correlations(
    module1_inputs = module1_inputs,
    prediction_footprints = prediction_footprints,
    out_dir = out_dir,
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff,
    fdr_cutoff = fdr_cutoff,
    min_non_na = min_non_na,
    cores = cores,
    write_outputs = write_outputs,
    output_format = output_format,
    return_prediction_stats = return_prediction_stats,
    verbose = verbose
  )
}
