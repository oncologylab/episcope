
# Step-by-step Module 2 public helpers.
#' Prepare Module 2 linking inputs
#'
#' @param multiomic_data CraftGRN multiomic object.
#' @param predicted_tfbs Predicted TFBS table or path from Module 1.
#' @param gene_tss Optional gene TSS table or path.
#' @param regulatory_prior Optional generic FP-target prior.
#' @param project_config Optional project config path or list.
#' @param max_distance_bp Maximum signed distance to TSS.
#' @param verbose Emit concise progress messages.
#' @return A list of normalized Module 2 inputs.
#' @noRd
module2_prepare_link_inputs <- function(multiomic_data, predicted_tfbs, gene_tss = NULL, regulatory_prior = NULL, project_config = NULL, max_distance_bp = NULL, verbose = TRUE) {
  cfg <- .module2_cfg(project_config)
  if (is.null(max_distance_bp)) max_distance_bp <- as.numeric(.module2_cfg_value(cfg, "max_distance_bp", .module2_cfg_value(cfg, "link_window_bp", 100000)))[[1L]]
  if (!is_multiomic_object(multiomic_data)) .log_abort("multiomic_data must be a CraftGRN multiomic object.")
  validate_multiomic_object(multiomic_data)
  if (is.character(predicted_tfbs) && length(predicted_tfbs) == 1L && file.exists(predicted_tfbs)) predicted_tfbs <- load_predicted_tfbs(predicted_tfbs)
  predicted_tfbs <- build_predicted_tfbs(predicted_tfbs)
  gene_tss <- .module2_resolve_gene_tss(gene_tss, cfg = cfg, multiomic_data = multiomic_data, verbose = verbose)
  mats <- multiomic_data$matrices
  expressed_genes <- rownames(mats$gene_on)[rowSums(mats$gene_on > 0, na.rm = TRUE) > 0]
  bound_fps <- rownames(mats$fp_bound)[rowSums(mats$fp_bound > 0, na.rm = TRUE) > 0]
  predicted_tfbs <- predicted_tfbs[predicted_tfbs$fp_id %in% bound_fps & predicted_tfbs$tf %in% expressed_genes, , drop = FALSE]
  target_genes <- intersect(expressed_genes, gene_tss$target_gene)
  tfs <- sort(unique(as.character(predicted_tfbs$tf)))
  if (isTRUE(verbose)) .log_inform("Module 2 inputs: {length(tfs)} TF(s), {length(target_genes)} target gene(s), {nrow(predicted_tfbs)} predicted TFBS row(s).")
  list(multiomic_data = multiomic_data, predicted_tfbs = predicted_tfbs, gene_tss = gene_tss, regulatory_prior = regulatory_prior, cfg = cfg, max_distance_bp = max_distance_bp, expressed_genes = expressed_genes, bound_fps = bound_fps, target_genes = target_genes, tfs = tfs)
}
#' Correlate TF expression with target gene expression
#'
#' @param module2_inputs Output from module2_identify_candidate_links.
#' @param n_cores Number of worker cores; NULL uses all available cores.
#' @param verbose Emit concise progress messages.
#' @return A TF-target correlation table with pass flags.
#' @export
module2_correlate_tf_targets <- function(module2_inputs, n_cores = NULL, verbose = TRUE) {
  tf_pairs <- tidyr::crossing(tf = module2_inputs$tfs, target_gene = module2_inputs$target_genes)
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: testing {nrow(tf_pairs)} pair(s).")
  cfg <- module2_inputs$cfg
  cutoffs <- .module2_corr_cutoffs(cfg, "tf_target", r_default = .module2_cfg_value(cfg, "threshold_rna_gene_corr_r", 0.3))
  if (is.null(cutoffs$p)) cutoffs$p <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_p", NULL)
  if (is.null(cutoffs$fdr)) cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_fdr", NULL)
  out <- .module2_pair_correlations(module2_inputs$multiomic_data$matrices$gene_expr, module2_inputs$multiomic_data$matrices$gene_expr, tf_pairs, "tf", "target_gene", cutoffs, n_cores = n_cores)
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: {sum(out$pass %in% TRUE)} pair(s) passed.")
  out
}
#' Build restricted FP-target candidates
#'
#' @param module2_inputs Output from module2_identify_candidate_links.
#' @param tf_target_corr Output from module2_correlate_tf_targets.
#' @param verbose Emit concise progress messages.
#' @return A candidate table restricted by TF-target pass calls and genomic priors.
#' @noRd
module2_build_fp_target_candidates <- function(module2_inputs, tf_target_corr, verbose = TRUE) {
  tf_pass <- tf_target_corr[tf_target_corr$pass %in% TRUE, , drop = FALSE]
  cand <- .module2_build_candidates(module2_inputs$predicted_tfbs, tf_pass, module2_inputs$gene_tss, regulatory_prior = module2_inputs$regulatory_prior, max_distance_bp = module2_inputs$max_distance_bp)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target candidates: {nrow(cand)} row(s).")
  cand
}

#' Link TFs to potential target genes based on TFBS-TSS proximity or 3D interaction data
#'
#' @param multiomic_data CraftGRN multiomic object.
#' @param predicted_tfbs Predicted TFBS table or path from Module 1.
#' @param gene_tss Optional gene TSS table or path.
#' @param regulatory_prior Optional generic FP-target prior.
#' @param project_config Optional project config path or list.
#' @param max_distance_bp Maximum signed distance to TSS.
#' @param verbose Emit concise progress messages.
#' @return A list of normalized Module 2 inputs used by downstream step functions.
#' @export
module2_identify_candidate_links <- function(multiomic_data, predicted_tfbs, gene_tss = NULL, regulatory_prior = NULL, project_config = NULL, max_distance_bp = NULL, verbose = TRUE) {
  module2_prepare_link_inputs(
    multiomic_data = multiomic_data,
    predicted_tfbs = predicted_tfbs,
    gene_tss = gene_tss,
    regulatory_prior = regulatory_prior,
    project_config = project_config,
    max_distance_bp = max_distance_bp,
    verbose = verbose
  )
}

#' Build restricted candidate FP-target links
#'
#' @param module2_inputs Output from internal Module 2 input preparation.
#' @param tf_target_corr Output from module2_correlate_tf_targets.
#' @param verbose Emit concise progress messages.
#' @return A candidate table restricted by TF-target pass calls and genomic priors.
#' @export
module2_link_fp_targets <- function(module2_inputs, tf_target_corr, verbose = TRUE) {
  module2_build_fp_target_candidates(module2_inputs, tf_target_corr, verbose = verbose)
}
#' Correlate FP score with target gene expression
#'
#' @param module2_inputs Output from module2_identify_candidate_links.
#' @param candidates Output from module2_link_fp_targets.
#' @param n_cores Number of worker cores; NULL uses all available cores.
#' @param verbose Emit concise progress messages.
#' @return An FP-target correlation table with pass flags.
#' @export
module2_correlate_fp_targets <- function(module2_inputs, candidates, n_cores = NULL, verbose = TRUE) {
  fp_pairs <- unique(candidates[, c("fp_id", "target_gene"), drop = FALSE])
  cfg <- module2_inputs$cfg
  cutoffs <- .module2_corr_cutoffs(cfg, "fp_target", r_default = .module2_cfg_value(cfg, "threshold_fp_gene_corr_r", 0.3))
  if (is.null(cutoffs$p)) cutoffs$p <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_p", NULL)
  if (is.null(cutoffs$fdr)) cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_fdr", NULL)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: testing {nrow(fp_pairs)} restricted pair(s).")
  out <- .module2_pair_correlations(module2_inputs$multiomic_data$matrices$fp_score, module2_inputs$multiomic_data$matrices$gene_expr, fp_pairs, "fp_id", "target_gene", cutoffs, n_cores = n_cores)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: {sum(out$pass %in% TRUE)} pair(s) passed.")
  out
}
#' Assemble final Module 2 TF-FP-target links
#'
#' @param module2_inputs Output from module2_identify_candidate_links.
#' @param candidates Candidate table.
#' @param tf_target_corr TF-target correlation table.
#' @param fp_target_corr FP-target correlation table.
#' @param output_dir Optional output directory.
#' @param output_format One of auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#' @return A Module 2 result list.
#' @noRd
module2_assemble_links <- function(module2_inputs, candidates, tf_target_corr, fp_target_corr, output_dir = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE) {
  output_format <- match.arg(output_format)
  pred_dt <- unique(data.table::as.data.table(module2_inputs$predicted_tfbs[, c("fp_id", "tf"), drop = FALSE]))
  cand_dt <- unique(data.table::as.data.table(candidates[, c("candidate_id", "fp_id", "target_gene"), drop = FALSE]))
  tf_pass_dt <- unique(data.table::as.data.table(tf_target_corr[tf_target_corr$pass %in% TRUE, c("tf", "target_gene"), drop = FALSE]))
  fp_pass_dt <- unique(data.table::as.data.table(fp_target_corr[fp_target_corr$pass %in% TRUE, c("fp_id", "target_gene"), drop = FALSE]))
  links <- pred_dt[cand_dt, on = "fp_id", allow.cartesian = TRUE, nomatch = 0L]
  links <- tf_pass_dt[links, on = c("tf", "target_gene"), nomatch = 0L]
  links <- fp_pass_dt[links, on = c("fp_id", "target_gene"), nomatch = 0L]
  links <- unique(links); links$link_id <- sprintf("link_%08d", seq_len(nrow(links))); links$tf_target_pass <- TRUE; links$fp_target_pass <- TRUE; links$module2_link_pass <- TRUE
  links <- tibble::as_tibble(links[, c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "tf_target_pass", "fp_target_pass", "module2_link_pass"), drop = FALSE])
  activity <- .module2_condition_activity(links, module2_inputs$predicted_tfbs, module2_inputs$multiomic_data)
  qc_summary <- tibble::tibble(metric = c("n_predicted_tfbs", "n_tfs", "n_target_genes", "n_tf_target_pairs_tested", "n_tf_target_pairs_pass", "n_fp_target_candidates", "n_fp_target_pairs_tested", "n_fp_target_pairs_pass", "n_module2_links", "n_active_link_conditions"), value = c(nrow(module2_inputs$predicted_tfbs), length(module2_inputs$tfs), length(module2_inputs$target_genes), nrow(tf_target_corr), sum(tf_target_corr$pass %in% TRUE), nrow(candidates), nrow(unique(candidates[, c("fp_id", "target_gene"), drop = FALSE])), sum(fp_target_corr$pass %in% TRUE), nrow(links), sum(activity$active %in% TRUE)))
  out <- list(predicted_tfbs = module2_inputs$predicted_tfbs, candidates = candidates, tf_target_corr = tf_target_corr, fp_target_corr = fp_target_corr, links = links, condition_activity = activity, qc_summary = qc_summary, manifest = tibble::tibble(), reports = list(), parameters = list(max_distance_bp = module2_inputs$max_distance_bp, n_cores = .module2_default_cores(NULL), output_format = output_format))
  if (!is.null(output_dir) && nzchar(output_dir)) {
    if (isTRUE(verbose)) .log_inform("Module 2 writing step-by-step outputs to {output_dir}.")
    data_dir <- file.path(output_dir, "data")
    manifest <- dplyr::bind_rows(
      .module2_manifest_table(.module2_write_table(module2_inputs$predicted_tfbs, file.path(data_dir, "predicted_tfbs"), "predicted_tfbs", output_format), "module1_predicted_tfbs"),
      .module2_manifest_table(.module2_write_table(candidates, file.path(data_dir, "candidates"), "fp_target_candidates", output_format), "module2_fp_target_candidates"),
      .module2_manifest_table(.module2_write_table(tf_target_corr, file.path(data_dir, "correlations"), "tf_target_corr", output_format), "module2_tf_target_corr"),
      .module2_manifest_table(.module2_write_table(fp_target_corr, file.path(data_dir, "correlations"), "fp_target_corr", output_format), "module2_fp_target_corr"),
      .module2_manifest_table(.module2_write_table(links, file.path(data_dir, "links"), "module2_links", output_format), "module2_links"),
      .module2_manifest_table(.module2_write_table(activity, file.path(data_dir, "condition_activity"), "condition_activity", output_format), "module2_condition_activity"),
      .module2_write_run_summary(qc_summary, output_dir)
    )
    manifest_path <- file.path(output_dir, "module2_manifest.csv")
    readr::write_csv(manifest, manifest_path)
    out$manifest <- manifest
    out$reports$manifest <- manifest_path
  }
  class(out) <- c("craftgrn_module2", "list"); out
}

#' Assemble, filter, and output final predicted TF-FP-target links
#'
#' @param module2_inputs Output from internal Module 2 input preparation.
#' @param candidates Candidate table.
#' @param tf_target_corr TF-target correlation table.
#' @param fp_target_corr FP-target correlation table.
#' @param output_dir Optional output directory.
#' @param output_format One of auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#' @return A Module 2 result list.
#' @noRd
output_predicted_links <- function(module2_inputs, candidates, tf_target_corr, fp_target_corr, output_dir = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE) {
  module2_assemble_links(
    module2_inputs = module2_inputs,
    candidates = candidates,
    tf_target_corr = tf_target_corr,
    fp_target_corr = fp_target_corr,
    output_dir = output_dir,
    output_format = output_format,
    verbose = verbose
  )
}

#' Assemble, filter, and output final predicted TF-FP-target links
#'
#' @param module2_inputs Output from [module2_identify_candidate_links()].
#' @param candidates Candidate table from [module2_link_fp_targets()].
#' @param tf_target_corr TF-target correlation table from
#'   [module2_correlate_tf_targets()].
#' @param fp_target_corr FP-target correlation table from
#'   [module2_correlate_fp_targets()].
#' @param output_dir Optional output directory.
#' @param output_format One of auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#'
#' @return A Module 2 result list.
#' @export
module2_output_predicted_links <- function(module2_inputs,
                                           candidates,
                                           tf_target_corr,
                                           fp_target_corr,
                                           output_dir = NULL,
                                           output_format = c("auto", "parquet", "csv"),
                                           verbose = TRUE) {
  module2_assemble_links(
    module2_inputs = module2_inputs,
    candidates = candidates,
    tf_target_corr = tf_target_corr,
    fp_target_corr = fp_target_corr,
    output_dir = output_dir,
    output_format = output_format,
    verbose = verbose
  )
}
