# File: utils_config.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide package-level configuration state, validation, and path normalization.
#
# Inputs:
# - YAML config files
# - required config keys and path keys
# - optional target environment for validation or normalization
#
# Outputs:
# - normalized in-memory config state used by package helpers
# - validated config values with clear errors on missing or invalid settings
#
# Notes:
# - Keep this file generic and reusable across modules.
# - Step-specific interpretation of config values belongs in module files.

# Internal package config state used by package helpers.
.craftgrn_state <- new.env(parent = emptyenv())

.craftgrn_get_option <- function(name, default = NULL) {
  getOption(paste0("craftgrn.", name), default = default)
}

.craftgrn_getenv <- function(name, unset = "") {
  Sys.getenv(paste0("CRAFTGRN_", name), unset = unset)
}

.cfg_get <- function(name, default = NULL) {
  if (exists(name, envir = .craftgrn_state, inherits = FALSE)) {
    return(get(name, envir = .craftgrn_state, inherits = FALSE))
  }
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    return(get(name, envir = .GlobalEnv, inherits = FALSE))
  }
  default
}

.cfg_set <- function(name, value) {
  assign(name, value, envir = .craftgrn_state)
  invisible(value)
}

.path_is_absolute <- function(path) {
  path <- path.expand(path)
  grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
}

.resolve_config_path <- function(path, config_dir = NULL) {
  if (!is.character(path) || !length(path) || !nzchar(path[1])) {
    return(path)
  }
  path <- path.expand(path[1])
  if (.path_is_absolute(path) || is.null(config_dir) || !nzchar(config_dir)) {
    return(path)
  }
  normalizePath(file.path(config_dir, path), winslash = "/", mustWork = FALSE)
}

.cfg_apply_aliases <- function(env = .craftgrn_state) {
  if (exists("threshold_expr", envir = env, inherits = FALSE)) {
    thr_expr <- get("threshold_expr", envir = env, inherits = FALSE)
    if (!exists("threshold_gene_expr", envir = env, inherits = FALSE)) {
      assign("threshold_gene_expr", thr_expr, envir = env)
    }
    if (!exists("threshold_tf_expr", envir = env, inherits = FALSE)) {
      assign("threshold_tf_expr", thr_expr, envir = env)
    }
  }

  if (!exists("link_window_bp", envir = env, inherits = FALSE)) {
    assign("link_window_bp", 200000L, envir = env)
  }

  invisible(TRUE)
}

.project_config_allowed_keys <- function() {
  c(
    "all_mode_tf_chunk_size", "atac_master", "atac_score_threshold",
    "basal_links_dir", "base_dir", "comparison_metadata", "condition_colors",
    "db", "demo_chromosome", "diff_dir_name", "filter_to_canonical_bound",
    "footprint_sample_scope", "force_rebuild_fp_manifest", "fp_delta_cutoff",
    "fp_input_format", "fp_tools_compact_dir",
    "fp_filter_mode", "fp_log2fc_cutoff", "fp_root_dir", "fp_score_threshold",
    "gene_de_padj_cutoff", "gene_log2fc_cutoff", "gene_tss",
    "genehancer_path", "link_gene_id_col", "link_score_threshold",
    "link_window_bp", "metadata_path", "module1_label_col", "module2",
    "module2_run_label", "module3", "module3_cell_subset",
    "module3_comparison_workers", "module3_condition_multiomic",
    "module3_label_col", "motif_db_path", "multiomic_object", "multiomic_object_raw",
    "pathway_backend", "pathway_databases", "pathway_organism",
    "pathway_species", "ref_genome", "report", "report_state",
    "rna_de_results", "rna_mapped", "rna_path", "sample_metadata",
    "batch_correction", "resources",
    "source_project", "step1_output_dir_name", "tf_binding_mode",
    "tf_target_link_mode", "threshold_atac_gene_corr_abs_r",
    "threshold_atac_gene_corr_p", "threshold_expr",
    "threshold_fp_gene_corr_abs_r", "threshold_fp_gene_corr_fdr",
    "threshold_fp_gene_corr_p", "threshold_fp_gene_corr_r",
    "threshold_fp_score", "threshold_fp_tf_corr_fdr",
    "threshold_fp_tf_corr_p", "threshold_fp_tf_corr_r",
    "threshold_fp_target_corr_fdr", "threshold_fp_target_corr_p",
    "threshold_fp_target_corr_r",
    "threshold_gene_expr", "threshold_link_score",
    "threshold_rna_gene_corr_abs_r", "threshold_rna_gene_corr_fdr",
    "threshold_rna_gene_corr_p", "threshold_rna_gene_corr_r",
    "threshold_tf_target_corr_fdr", "threshold_tf_target_corr_p",
    "threshold_tf_target_corr_r",
    "threshold_tf_expr", "topic_benchmark_enabled", "topic_benchmark_k_grid",
    "topic_benchmark_methods", "topic_binarize_method", "topic_count_input",
    "topic_count_method", "topic_count_scale", "topic_dir_name", "topic_gammafit_min_terms",
    "topic_gammafit_scope", "topic_gammafit_scopes", "topic_gammafit_thrP",
    "topic_k", "topic_k_grid", "topic_link_method", "topic_link_output",
    "topic_link_prob_cutoff", "topic_method", "topic_score_method",
    "topic_term_assignment_method",
    "topic_optimize_topics", "topic_merge_min_genes",
    "topic_merge_min_links", "topic_merge_similarity_threshold",
    "topic_assignment_qc", "topic_qc_umap_links_per_condition",
    "topic_qc_top_tfs", "topic_qc_reference_condition",
    "topic_condition_gene_weighting", "topic_condition_gene_expression_file",
    "topic_condition_peak_weighting",
    "topic_condition_specificity_temperature", "topic_condition_specificity_floor",
    "topic_condition_specificity_expression_min",
    "topic_condition_gene_expression_min",
    "topic_qc_upregulated_log2fc_min",
    "topic_qc_upregulated_pseudocount", "topic_qc_seed",
    "topic_raw_theta_document_heatmap", "run_raw_theta_document_heatmap",
    "topic_tf_membership_cutoff", "topic_tf_primary_margin_cutoff",
    "topic_vae_batch_size", "topic_vae_device", "warplda_iterations",
    "topic_warplda_sampler", "topic_warplda_beta", "topic_warplda_seed",
    "waterfall_min_abs_net",
    paste0("module3_", c(
      "pathway_backend", "pathway_databases", "pathway_organism",
      "pathway_species", "topic_benchmark_enabled", "topic_benchmark_k_grid",
      "topic_benchmark_methods", "topic_count_input", "topic_count_method",
      "topic_count_scale",
      "topic_gammafit_min_terms", "topic_gammafit_scope",
      "topic_gammafit_scopes", "topic_gammafit_thrP", "topic_k",
      "topic_k_grid", "topic_link_method", "topic_link_output",
      "topic_link_prob_cutoff", "topic_method", "topic_score_method",
      "topic_term_assignment_method",
      "topic_optimize_topics", "topic_merge_min_genes",
      "topic_merge_min_links", "topic_merge_similarity_threshold",
      "topic_assignment_qc", "topic_qc_umap_links_per_condition",
      "topic_qc_top_tfs", "topic_qc_reference_condition",
      "topic_condition_gene_weighting", "topic_condition_gene_expression_file",
      "topic_condition_peak_weighting",
      "topic_condition_specificity_temperature", "topic_condition_specificity_floor",
      "topic_condition_specificity_expression_min",
      "topic_condition_gene_expression_min",
      "topic_qc_upregulated_log2fc_min",
      "topic_qc_upregulated_pseudocount", "topic_qc_seed",
      "topic_raw_theta_document_heatmap",
      "topic_tf_membership_cutoff", "topic_tf_primary_margin_cutoff",
      "topic_vae_batch_size", "topic_vae_device", "warplda_iterations",
      "topic_warplda_sampler", "topic_warplda_beta", "topic_warplda_seed"
    )),
    "topic_warplda_iterations",
    paste0("module3_", c(
      "optimize_topics", "merge_min_genes", "merge_min_links",
      "merge_similarity_threshold", "run_topic_assignment_qc",
      "qc_umap_links_per_condition", "qc_top_tfs",
      "qc_reference_condition", "qc_upregulated_log2fc_min",
      "qc_upregulated_pseudocount", "qc_seed"
    ))
  )
}

.project_config_suggest_key <- function(key, allowed) {
  hit <- utils::adist(key, allowed)
  if (!length(hit)) return("")
  best <- allowed[[which.min(hit)]]
  if (min(hit) > max(2L, floor(nchar(key) / 3))) return("")
  paste0(" Did you mean `", best, "`?")
}

.project_config_assert_mapping <- function(x, label) {
  if (is.list(x) && !length(x)) return(invisible(TRUE))
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x)))) {
    cli::cli_abort("`{label}` must be a named YAML mapping.")
  }
  if (anyDuplicated(names(x))) {
    cli::cli_abort("`{label}` contains duplicate entries: {paste(unique(names(x)[duplicated(names(x))]), collapse = ', ')}.")
  }
  invisible(TRUE)
}

.project_config_scalar <- function(x, type, nullable = FALSE) {
  if (is.null(x)) return(isTRUE(nullable))
  length(x) == 1L && !is.na(x) && switch(
    type,
    character = is.character(x) && nzchar(x),
    logical = is.logical(x),
    numeric = is.numeric(x) && is.finite(x),
    integer = is.numeric(x) && is.finite(x) && x == as.integer(x),
    FALSE
  )
}

.validate_report_config <- function(report) {
  if (is.null(report)) return(invisible(TRUE))
  .project_config_assert_mapping(report, "report")
  allowed <- c("condition_colors", "condition_order", "defaults")
  unknown <- setdiff(names(report), allowed)
  if (length(unknown)) {
    cli::cli_abort("Unknown `report` entries: {paste(unknown, collapse = ', ')}.")
  }

  if (!is.null(report$condition_order)) {
    order <- report$condition_order
    if (!is.character(order) || anyNA(order) || any(!nzchar(order))) {
      cli::cli_abort("`report.condition_order` must contain non-empty condition IDs.")
    }
    if (anyDuplicated(order)) {
      cli::cli_abort("`report.condition_order` must not contain duplicate condition IDs.")
    }
  }
  if (!is.null(report$condition_colors)) {
    colors <- report$condition_colors
    .project_config_assert_mapping(colors, "report.condition_colors")
    vals <- unlist(colors, use.names = FALSE)
    if (!is.character(vals) || any(!grepl("^#[0-9A-Fa-f]{6}$", vals))) {
      cli::cli_abort("`report.condition_colors` values must be six-digit hex colors such as `#4E79A7`.")
    }
  }
  if (!is.null(report$defaults)) {
    defaults <- report$defaults
    .project_config_assert_mapping(defaults, "report.defaults")
    allowed_defaults <- c(
      "condition_1", "condition_2", "condition_1_suffix",
      "condition_2_suffix", "pathway", "short_condition_names", "tf", "topic",
      "topic_space", "trained_k"
    )
    unknown <- setdiff(names(defaults), allowed_defaults)
    if (length(unknown)) {
      cli::cli_abort("Unknown `report.defaults` entries: {paste(unknown, collapse = ', ')}.")
    }
    char_keys <- intersect(
      c("condition_1", "condition_2", "condition_1_suffix", "condition_2_suffix", "pathway", "tf"),
      names(defaults)
    )
    bad_char <- char_keys[!vapply(defaults[char_keys], .project_config_scalar, logical(1L), type = "character", nullable = TRUE)]
    if (length(bad_char)) {
      cli::cli_abort("Report defaults must be non-empty strings or null: {paste(bad_char, collapse = ', ')}.")
    }
    if (!is.null(defaults[["topic"]]) &&
        (!.project_config_scalar(defaults[["topic"]], "integer") ||
          defaults[["topic"]] < 1L)) {
      cli::cli_abort("`report.defaults.topic` must be a positive integer.")
    }
    if (!is.null(defaults[["trained_k"]]) &&
        (!.project_config_scalar(defaults[["trained_k"]], "integer") ||
          defaults[["trained_k"]] < 2L)) {
      cli::cli_abort("`report.defaults.trained_k` must be an integer greater than 1.")
    }
    if (!is.null(defaults[["topic_space"]]) &&
        (!.project_config_scalar(defaults[["topic_space"]], "character") ||
          !defaults[["topic_space"]] %in% c("raw", "combined"))) {
      cli::cli_abort("`report.defaults.topic_space` must be raw or combined.")
    }
    if (!is.null(defaults[["short_condition_names"]]) &&
        !.project_config_scalar(defaults[["short_condition_names"]], "logical")) {
      cli::cli_abort("`report.defaults.short_condition_names` must be true or false.")
    }
  }
  invisible(TRUE)
}

.validate_resource_config <- function(resources) {
  if (is.null(resources)) return(invisible(TRUE))
  .project_config_assert_mapping(resources, "resources")
  allowed <- c("max_memory_fraction", "reserve_memory_gb", "max_workers")
  unknown <- setdiff(names(resources), allowed)
  if (length(unknown)) {
    cli::cli_abort("Unknown `resources` entries: {paste(unknown, collapse = ', ')}.")
  }
  if (!is.null(resources$max_memory_fraction) &&
      (!.project_config_scalar(resources$max_memory_fraction, "numeric") ||
       resources$max_memory_fraction < 0.1 || resources$max_memory_fraction > 0.9)) {
    cli::cli_abort("`resources.max_memory_fraction` must be between 0.1 and 0.9.")
  }
  if (!is.null(resources$reserve_memory_gb) &&
      (!.project_config_scalar(resources$reserve_memory_gb, "numeric") || resources$reserve_memory_gb < 0)) {
    cli::cli_abort("`resources.reserve_memory_gb` must be non-negative.")
  }
  if (!is.null(resources$max_workers) &&
      (!.project_config_scalar(resources$max_workers, "integer") || resources$max_workers < 1L)) {
    cli::cli_abort("`resources.max_workers` must be a positive integer.")
  }
  invisible(TRUE)
}

.validate_batch_correction_config <- function(batch_correction) {
  if (is.null(batch_correction)) return(invisible(TRUE))
  .project_config_assert_mapping(batch_correction, "batch_correction")
  allowed <- c(
    "enabled", "method", "batch_column", "preserve_column", "k_candidates",
    "condition_distance_spearman_min", "effect_spearman_min",
    "direction_concordance_min", "between_condition_variance_min"
  )
  unknown <- setdiff(names(batch_correction), allowed)
  if (length(unknown)) {
    cli::cli_abort("Unknown `batch_correction` entries: {paste(unknown, collapse = ', ')}.")
  }
  if (!is.null(batch_correction$enabled) &&
      !.project_config_scalar(batch_correction$enabled, "logical")) {
    cli::cli_abort("`batch_correction.enabled` must be true or false.")
  }
  if (!is.null(batch_correction$method) &&
      (!.project_config_scalar(batch_correction$method, "character") ||
       !batch_correction$method %in% c("none", "ruvr"))) {
    cli::cli_abort("`batch_correction.method` must be `none` or `ruvr`.")
  }
  char_keys <- intersect(c("batch_column", "preserve_column"), names(batch_correction))
  bad_char <- char_keys[!vapply(batch_correction[char_keys], .project_config_scalar, logical(1L), type = "character")]
  if (length(bad_char)) {
    cli::cli_abort("Batch-correction column entries must be non-empty strings: {paste(bad_char, collapse = ', ')}.")
  }
  if (!is.null(batch_correction$k_candidates)) {
    values <- batch_correction$k_candidates
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
        any(values < 0) || any(values != as.integer(values))) {
      cli::cli_abort("`batch_correction.k_candidates` must contain non-negative integers.")
    }
  }
  metric_keys <- intersect(c(
    "condition_distance_spearman_min", "effect_spearman_min",
    "direction_concordance_min", "between_condition_variance_min"
  ), names(batch_correction))
  bad_metric <- metric_keys[!vapply(batch_correction[metric_keys], function(x) {
    .project_config_scalar(x, "numeric") && x >= 0 && x <= 1
  }, logical(1L))]
  if (length(bad_metric)) {
    cli::cli_abort("Batch-correction acceptance metrics must be between 0 and 1: {paste(bad_metric, collapse = ', ')}.")
  }
  invisible(TRUE)
}

.validate_module2_config <- function(module2) {
  if (is.null(module2)) return(invisible(TRUE))
  .project_config_assert_mapping(module2, "module2")
  allowed <- c(
    "gene_tss", "gene_tss_path", "link_window_bp", "max_distance_bp",
    "ref_genome", "threshold_fp_gene_corr_fdr", "threshold_fp_gene_corr_p",
    "threshold_fp_gene_corr_r", "threshold_rna_gene_corr_fdr",
    "threshold_rna_gene_corr_p", "threshold_rna_gene_corr_r",
    "threshold_tf_target_corr_fdr", "threshold_tf_target_corr_p",
    "threshold_tf_target_corr_r", "threshold_fp_target_corr_fdr",
    "threshold_fp_target_corr_p", "threshold_fp_target_corr_r"
  )
  unknown <- setdiff(names(module2), allowed)
  if (length(unknown)) {
    cli::cli_abort("Unknown `module2` entries: {paste(unknown, collapse = ', ')}.")
  }
  for (key in intersect(c("link_window_bp", "max_distance_bp"), names(module2))) {
    if (!.project_config_scalar(module2[[key]], "numeric") || module2[[key]] <= 0) {
      cli::cli_abort("`module2.{key}` must be a positive number.")
    }
  }
  for (key in intersect(c(
    "threshold_fp_gene_corr_p", "threshold_fp_gene_corr_fdr",
    "threshold_rna_gene_corr_p", "threshold_rna_gene_corr_fdr",
    "threshold_tf_target_corr_p", "threshold_tf_target_corr_fdr",
    "threshold_fp_target_corr_p", "threshold_fp_target_corr_fdr"
  ), names(module2))) {
    value <- module2[[key]]
    if (!is.null(value) && (!.project_config_scalar(value, "numeric") || value < 0 || value > 1)) {
      cli::cli_abort("`module2.{key}` must be between 0 and 1 or null.")
    }
  }
  for (key in intersect(c(
    "threshold_fp_gene_corr_r", "threshold_rna_gene_corr_r",
    "threshold_tf_target_corr_r", "threshold_fp_target_corr_r"
  ), names(module2))) {
    value <- module2[[key]]
    if (!.project_config_scalar(value, "numeric") || value < -1 || value > 1) {
      cli::cli_abort("`module2.{key}` must be between -1 and 1.")
    }
  }
  for (key in intersect(c("gene_tss", "gene_tss_path", "ref_genome"), names(module2))) {
    if (!.project_config_scalar(module2[[key]], "character")) {
      cli::cli_abort("`module2.{key}` must be a non-empty string.")
    }
  }
  invisible(TRUE)
}

.validate_module3_config <- function(module3) {
  if (is.null(module3)) return(invisible(TRUE))
  .project_config_assert_mapping(module3, "module3")
  allowed <- sub(
    "^module3_",
    "",
    grep("^module3_", .project_config_allowed_keys(), value = TRUE)
  )
  unknown <- setdiff(names(module3), allowed)
  if (length(unknown)) {
    cli::cli_abort("Unknown `module3` entries: {paste(unknown, collapse = ', ')}.")
  }
  flattened <- module3
  names(flattened) <- paste0("module3_", names(flattened))
  .validate_project_config(flattened, source = "module3")
  invisible(TRUE)
}

.validate_project_config <- function(cfg, source = "project config") {
  .project_config_assert_mapping(cfg, source)
  allowed <- .project_config_allowed_keys()
  unknown <- setdiff(names(cfg), allowed)
  if (length(unknown)) {
    hints <- vapply(unknown, .project_config_suggest_key, character(1L), allowed = allowed)
    cli::cli_abort(paste0(
      "Unknown project config entr", if (length(unknown) == 1L) "y" else "ies",
      ": ", paste0("`", unknown, "`", hints, collapse = ", "), "."
    ))
  }

  logical_keys <- intersect(c(
    "filter_to_canonical_bound", "force_rebuild_fp_manifest",
    "topic_benchmark_enabled", "module3_topic_benchmark_enabled",
    "topic_optimize_topics", "module3_topic_optimize_topics",
    "topic_assignment_qc", "module3_topic_assignment_qc",
    "topic_raw_theta_document_heatmap",
    "module3_topic_raw_theta_document_heatmap",
    "run_raw_theta_document_heatmap",
    "module3_optimize_topics", "module3_run_topic_assignment_qc"
  ), names(cfg))
  bad_logical <- logical_keys[!vapply(cfg[logical_keys], .project_config_scalar, logical(1L), type = "logical")]
  if (length(bad_logical)) {
    cli::cli_abort("Project config entries must be true or false: {paste(bad_logical, collapse = ', ')}.")
  }

  positive_integer_keys <- intersect(c(
    "all_mode_tf_chunk_size", "link_window_bp", "module3_comparison_workers",
    "topic_gammafit_min_terms", "module3_topic_gammafit_min_terms",
    "topic_vae_batch_size", "module3_topic_vae_batch_size",
    "topic_merge_min_genes", "module3_topic_merge_min_genes",
    "topic_merge_min_links", "module3_topic_merge_min_links",
    "topic_qc_umap_links_per_condition",
    "module3_topic_qc_umap_links_per_condition",
    "topic_qc_top_tfs", "module3_topic_qc_top_tfs",
    "topic_qc_seed", "module3_topic_qc_seed",
    "module3_merge_min_genes", "module3_merge_min_links",
    "module3_qc_umap_links_per_condition", "module3_qc_top_tfs",
    "module3_qc_seed",
    "warplda_iterations", "topic_warplda_iterations", "module3_warplda_iterations",
    "topic_warplda_seed", "module3_topic_warplda_seed"
  ), names(cfg))
  bad_integer <- positive_integer_keys[!vapply(cfg[positive_integer_keys], function(x) {
    .project_config_scalar(x, "integer") && x >= 1L
  }, logical(1L))]
  if (length(bad_integer)) {
    cli::cli_abort("Project config entries must be positive integers: {paste(bad_integer, collapse = ', ')}.")
  }

  positive_numeric_keys <- intersect(c(
    "topic_count_scale", "module3_topic_count_scale",
    "topic_warplda_beta", "module3_topic_warplda_beta",
    "topic_condition_specificity_temperature",
    "module3_topic_condition_specificity_temperature"
  ), names(cfg))
  bad_numeric <- positive_numeric_keys[!vapply(cfg[positive_numeric_keys], function(x) {
    is.null(x) || (.project_config_scalar(x, "numeric") && x > 0)
  }, logical(1L))]
  if (length(bad_numeric)) {
    cli::cli_abort("Project config entries must be null or positive numbers: {paste(bad_numeric, collapse = ', ')}.")
  }

  sampler_keys <- intersect(c(
    "topic_warplda_sampler", "module3_topic_warplda_sampler"
  ), names(cfg))
  allowed_samplers <- c("warp_omp", "warp_ref", "warp_mh", "gibbs_sync")
  bad_sampler <- sampler_keys[!vapply(cfg[sampler_keys], function(x) {
    .project_config_scalar(x, "character") && x %in% allowed_samplers
  }, logical(1L))]
  if (length(bad_sampler)) {
    cli::cli_abort("Topic WarpLDA samplers must be one of: {paste(allowed_samplers, collapse = ', ')}.")
  }

  reference_keys <- intersect(c(
    "topic_qc_reference_condition", "module3_topic_qc_reference_condition",
    "module3_qc_reference_condition"
  ), names(cfg))
  bad_reference <- reference_keys[!vapply(cfg[reference_keys], function(x) {
    is.null(x) || (.project_config_scalar(x, "character") && nzchar(trimws(x)))
  }, logical(1L))]
  if (length(bad_reference)) {
    cli::cli_abort(
      "Topic QC reference conditions must be one non-empty condition ID or null: {paste(bad_reference, collapse = ', ')}."
    )
  }
  expression_file_keys <- intersect(c(
    "topic_condition_gene_expression_file",
    "module3_topic_condition_gene_expression_file"
  ), names(cfg))
  bad_expression_file <- expression_file_keys[!vapply(cfg[expression_file_keys], function(x) {
    is.null(x) || (.project_config_scalar(x, "character") && nzchar(trimws(x)))
  }, logical(1L))]
  if (length(bad_expression_file)) {
    cli::cli_abort(
      "Condition gene expression files must be one non-empty path or null: {paste(bad_expression_file, collapse = ', ')}."
    )
  }

  nonnegative_qc_keys <- intersect(c(
    "topic_condition_gene_expression_min",
    "module3_topic_condition_gene_expression_min",
    "topic_condition_specificity_expression_min",
    "module3_topic_condition_specificity_expression_min",
    "topic_qc_upregulated_log2fc_min",
    "module3_topic_qc_upregulated_log2fc_min",
    "module3_qc_upregulated_log2fc_min"
  ), names(cfg))
  bad_nonnegative_qc <- nonnegative_qc_keys[
    !vapply(cfg[nonnegative_qc_keys], function(x) {
      .project_config_scalar(x, "numeric") && x >= 0
    }, logical(1L))
  ]
  if (length(bad_nonnegative_qc)) {
    cli::cli_abort(
      "Module 3 expression and QC log2FC cutoffs must be non-negative: {paste(bad_nonnegative_qc, collapse = ', ')}."
    )
  }

  positive_qc_keys <- intersect(c(
    "topic_qc_upregulated_pseudocount",
    "module3_topic_qc_upregulated_pseudocount",
    "module3_qc_upregulated_pseudocount"
  ), names(cfg))
  bad_positive_qc <- positive_qc_keys[!vapply(cfg[positive_qc_keys], function(x) {
    .project_config_scalar(x, "numeric") && x > 0
  }, logical(1L))]
  if (length(bad_positive_qc)) {
    cli::cli_abort(
      "Topic QC upregulated pseudocounts must be positive: {paste(bad_positive_qc, collapse = ', ')}."
    )
  }

  gammafit_probability_keys <- intersect(c(
    "topic_gammafit_thrP", "module3_topic_gammafit_thrP"
  ), names(cfg))
  valid_gammafit_probability <- function(x) {
    if (.project_config_scalar(x, "numeric")) return(x > 0 && x < 1)
    values <- if (is.list(x)) unlist(x, use.names = TRUE) else x
    value_names <- names(values)
    allowed_names <- c("lda", "multivi", "vae_mlp", "default")
    normalized_names <- tolower(gsub("-", "_", value_names, fixed = TRUE))
    normalized_names[normalized_names == "warplda"] <- "lda"
    normalized_names[normalized_names == "vae"] <- "vae_mlp"
    length(values) > 0L && !is.null(value_names) && all(nzchar(value_names)) &&
      all(normalized_names %in% allowed_names) && !anyDuplicated(normalized_names) &&
      all(vapply(as.list(values), function(value) {
        .project_config_scalar(value, "numeric") && value > 0 && value < 1
      }, logical(1L)))
  }
  bad_gammafit_probability <- gammafit_probability_keys[
    !vapply(cfg[gammafit_probability_keys], valid_gammafit_probability, logical(1L))
  ]
  if (length(bad_gammafit_probability)) {
    cli::cli_abort(paste0(
      "GammaFit probability config entries must be one probability strictly between 0 and 1, ",
      "or a named lda/multivi/vae_mlp/default mapping: ",
      paste(bad_gammafit_probability, collapse = ", "), "."
    ))
  }

  similarity_keys <- intersect(c(
    "topic_merge_similarity_threshold",
    "module3_topic_merge_similarity_threshold",
    "module3_merge_similarity_threshold"
  ), names(cfg))
  bad_similarity <- similarity_keys[!vapply(cfg[similarity_keys], function(x) {
    .project_config_scalar(x, "numeric") && x > 0 && x <= 1
  }, logical(1L))]
  if (length(bad_similarity)) {
    cli::cli_abort(
      "Topic merge similarity thresholds must be greater than 0 and at most 1: {paste(bad_similarity, collapse = ', ')}."
    )
  }

  probability_keys <- intersect(c(
    "gene_de_padj_cutoff", "threshold_fp_gene_corr_fdr",
    "threshold_fp_gene_corr_p", "threshold_fp_tf_corr_fdr",
    "threshold_fp_tf_corr_p", "threshold_fp_target_corr_fdr",
    "threshold_fp_target_corr_p", "threshold_rna_gene_corr_fdr",
    "threshold_rna_gene_corr_p", "threshold_tf_target_corr_fdr",
    "threshold_tf_target_corr_p", "topic_tf_membership_cutoff",
    "module3_topic_tf_membership_cutoff", "topic_tf_primary_margin_cutoff",
    "module3_topic_tf_primary_margin_cutoff"
  ), names(cfg))
  bad_probability <- probability_keys[!vapply(cfg[probability_keys], function(x) {
    is.null(x) || (.project_config_scalar(x, "numeric") && x >= 0 && x <= 1)
  }, logical(1L))]
  if (length(bad_probability)) {
    cli::cli_abort("Probability config entries must be between 0 and 1 or null: {paste(bad_probability, collapse = ', ')}.")
  }
  specificity_floor_keys <- intersect(c(
    "topic_condition_specificity_floor",
    "module3_topic_condition_specificity_floor"
  ), names(cfg))
  bad_specificity_floor <- specificity_floor_keys[
    !vapply(cfg[specificity_floor_keys], function(x) {
      .project_config_scalar(x, "numeric") && x > 0 && x < 1
    }, logical(1L))
  ]
  if (length(bad_specificity_floor)) {
    cli::cli_abort(
      "Condition specificity floors must be strictly between 0 and 1: {paste(bad_specificity_floor, collapse = ', ')}."
    )
  }
  correlation_keys <- grep("^threshold_.*_corr_(abs_)?r$", names(cfg), value = TRUE)
  bad_correlation <- correlation_keys[!vapply(correlation_keys, function(key) {
    value <- cfg[[key]]
    lower <- if (grepl("_abs_r$", key)) 0 else -1
    .project_config_scalar(value, "numeric") && value >= lower && value <= 1
  }, logical(1L))]
  if (length(bad_correlation)) {
    cli::cli_abort("Correlation cutoffs must be finite values between -1 and 1; absolute-R cutoffs must be non-negative: {paste(bad_correlation, collapse = ', ')}.")
  }
  for (key in intersect(c("topic_link_prob_cutoff", "module3_topic_link_prob_cutoff"), names(cfg))) {
    value <- cfg[[key]]
    valid <- identical(value, "max") ||
      (.project_config_scalar(value, "numeric") && value >= 0 && value <= 1)
    if (!valid) {
      cli::cli_abort("`{key}` must be between 0 and 1 or `max`.")
    }
  }

  enum_values <- list(
    fp_filter_mode = c("delta", "log2fc", "both"),
    fp_input_format = c("overview", "fp_tools_compact"),
    pathway_backend = c("enrichly", "enrichr"),
    module3_pathway_backend = c("enrichly", "enrichr"),
    pathway_species = c("human", "mouse", "human_mouse_best"),
    module3_pathway_species = c("human", "mouse", "human_mouse_best"),
    tf_binding_mode = c("all", "canonical_bound"),
    topic_binarize_method = c("gammafit", "probability"),
    topic_count_method = c("bin", "log"),
    module3_topic_count_method = c("bin", "log"),
    topic_gammafit_scope = c("topic_term_group", "global_term_group"),
    module3_topic_gammafit_scope = c("topic_term_group", "global_term_group"),
    topic_link_output = c("pass", "full", "both", "none"),
    module3_topic_link_output = c("pass", "full", "both", "none"),
    topic_score_method = c("normtop_specificity", "rowmax_phi"),
    module3_topic_score_method = c("normtop_specificity", "rowmax_phi"),
    topic_term_assignment_method = c("gammafit_maxprob", "max_phi", "gammafit"),
    module3_topic_term_assignment_method = c("gammafit_maxprob", "max_phi", "gammafit"),
    topic_condition_gene_weighting = c("none", "specificity"),
    module3_topic_condition_gene_weighting = c("none", "specificity"),
    topic_condition_peak_weighting = c("none", "tf_expression"),
    module3_topic_condition_peak_weighting = c("none", "tf_expression"),
    topic_vae_device = c("auto", "cpu", "cuda"),
    module3_topic_vae_device = c("auto", "cpu", "cuda")
  )
  for (key in intersect(names(enum_values), names(cfg))) {
    value <- cfg[[key]]
    if (!.project_config_scalar(value, "character") || !value %in% enum_values[[key]]) {
      cli::cli_abort("`{key}` must be one of: {paste(enum_values[[key]], collapse = ', ')}.")
    }
  }

  if (!is.null(cfg[["report"]]) && !is.null(cfg[["report_state"]])) {
    cli::cli_abort("Use `report`, not both `report` and the legacy `report_state` entry.")
  }
  .validate_report_config(cfg[["report"]] %||% cfg[["report_state"]])
  .validate_resource_config(cfg[["resources"]])
  .validate_batch_correction_config(cfg[["batch_correction"]])
  .validate_module2_config(cfg[["module2"]])
  .validate_module3_config(cfg[["module3"]])
  invisible(TRUE)
}

#' Validate config values
#'
#' Ensures required config keys (e.g. thresholds and db) exist in the chosen
#' environment before running pipelines.
#'
#' @param required Character vector of required variable names.
#' @param numeric_required Character vector of required numeric variable names.
#' @param env Environment to check. Defaults to the internal CraftGRN config state.
#' @param project_config Optional parsed project config or YAML path to validate
#'   against the package project schema. Unknown or malformed entries fail.
#'
#' @return \code{TRUE} invisibly when validation passes.
#' @export
validate_config <- function(
    required = c(
      "db",
      "ref_genome",
      "threshold_expr",
      "threshold_fp_score",
      "threshold_fp_tf_corr_r",
      "link_window_bp",
      "threshold_rna_gene_corr_r",
      "threshold_fp_gene_corr_r"
    ),
    numeric_required = c(
      "threshold_expr",
      "threshold_fp_score",
      "threshold_fp_tf_corr_r",
      "link_window_bp",
      "threshold_rna_gene_corr_r",
      "threshold_fp_gene_corr_r"
    ),
    env = .craftgrn_state,
    project_config = NULL
) {
  if (!is.null(project_config)) {
    cfg <- if (is.character(project_config) && length(project_config) == 1L) {
      if (!file.exists(project_config)) {
        cli::cli_abort("Config file not found: {project_config}")
      }
      yaml::read_yaml(project_config)
    } else if (is.list(project_config)) {
      project_config
    } else {
      cli::cli_abort("`project_config` must be a YAML path or named list.")
    }
    .validate_project_config(cfg)
  }

  missing <- required[!vapply(required, function(nm) exists(nm, envir = env, inherits = FALSE), logical(1))]
  if (length(missing)) {
    cli::cli_abort("Missing config values: {paste(missing, collapse = ', ')}.")
  }

  if (length(numeric_required)) {
    bad_numeric <- numeric_required[!vapply(numeric_required, function(nm) {
      val <- get(nm, envir = env)
      is.numeric(val) && length(val) == 1L && is.finite(val)
    }, logical(1))]
    if (length(bad_numeric)) {
      cli::cli_abort("Config values must be finite numeric scalars: {paste(bad_numeric, collapse = ', ')}.")
    }
  }

  invisible(TRUE)
}


#' Normalize configured path variables
#'
#' Expands \code{~} and resolves relative paths against the config directory.
#'
#' @param keys Character vector of config keys to normalize.
#' @param env Environment to update. Defaults to the internal CraftGRN config state.
#'
#' @return Invisibly returns the updated values (named character vector).
#' @noRd
normalize_config_paths <- function(
    keys = c(
      "base_dir",
      "fp_root_dir",
      "fp_manifest_path",
      "benchmark_tfbs_dir",
      "comparison_metadata",
      "pathway_dir",
      "genehancer_path",
      "gene_tss",
      "loop_path",
      "multiomic_object",
      "multiomic_object_raw",
      "motif_db_path",
      "preprocess_dir",
      "sample_metadata",
      "metadata_path",
      "sample_metadata_path",
      "atac_master",
      "atac_data_path",
      "atac_master_path",
      "rna_mapped",
      "rna_de_results",
      "rna_mapped_path",
      "rna_path",
      "basal_links_dir",
      "data_dir",
      "sample_metadata_out",
      "strict_metadata_out",
      "lenient_metadata_out",
      "atac_master_out",
      "rna_filtered_out",
      "strict_rna_out",
      "lenient_rna_out",
      "rna_mapped_out"
    ),
    env = .craftgrn_state
) {
  existing <- keys[vapply(keys, function(nm) exists(nm, envir = env, inherits = FALSE), logical(1))]
  if (!length(existing)) return(invisible(character(0)))

  vals <- list()
  for (nm in existing) {
    val <- get(nm, envir = env)
    if (!is.character(val) || !length(val) || !nzchar(val[1])) {
      next
    }
    config_dir <- if (exists(".config_dir", envir = env, inherits = FALSE)) {
      get(".config_dir", envir = env, inherits = FALSE)
    } else {
      NULL
    }
    vals[[nm]] <- .resolve_config_path(val[1], config_dir = config_dir)
  }

  if (!length(vals)) return(invisible(character(0)))

  list2env(vals, envir = env)
  invisible(unlist(vals, use.names = TRUE))
}

#' Load a CraftGRN YAML config into an environment
#'
#' Reads a YAML file and assigns each top-level key as a variable in
#' the target environment (e.g., `db`, `threshold_tf_expr`, etc.).
#' Also runs standard config initialization helpers when available.
#'
#' @param path Character path to a YAML file.
#' @param env Environment to populate. Defaults to the internal CraftGRN config state.
#' @return (Invisibly) the parsed list.
#' @examples
#' \donttest{
#' config_path <- tempfile(fileext = ".yaml")
#' writeLines(c(
#'   "db: JASPAR2024",
#'   "ref_genome: hg38",
#'   "threshold_expr: 1",
#'   "threshold_fp_score: 0",
#'   "threshold_fp_tf_corr_r: 0.3",
#'   "threshold_rna_gene_corr_r: 0.3",
#'   "threshold_fp_gene_corr_r: 0.3"
#' ), config_path)
#' load_config(config_path)
#' # Config values are now available to CraftGRN helper functions.
#' }
#' @export
load_config <- function(path, env = .craftgrn_state) {
  if (!file.exists(path)) {
    cli::cli_abort("Config file not found: {path}")
  }
  config_dir <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  cfg <- yaml::read_yaml(path)
  .validate_project_config(cfg, source = basename(path))
  list2env(cfg, envir = env)
  assign(".config_dir", config_dir, envir = env)
  .cfg_apply_aliases(env = env)
  if (exists("validate_config", mode = "function")) {
    validate_config(env = env)
  }
  if (exists("normalize_config_paths", mode = "function")) {
    normalize_config_paths(env = env)
  }
  if (exists("init_motif_db", mode = "function") && exists("db", envir = env, inherits = FALSE)) {
    ref_genome_val <- if (exists("ref_genome", envir = env, inherits = FALSE)) get("ref_genome", envir = env, inherits = FALSE) else NULL
    init_formals <- names(formals(init_motif_db))
    if (!is.null(init_formals) && "ref_genome" %in% init_formals) {
      motif_init <- init_motif_db(get("db", envir = env, inherits = FALSE), ref_genome = ref_genome_val)
    } else {
      motif_init <- init_motif_db(get("db", envir = env, inherits = FALSE))
    }
    assign("motif_init", motif_init, envir = env)
    assign("motif_db", motif_init$motif_db, envir = env)
    assign("tf_list_all", motif_init$tf_list, envir = env)
    assign("tf_list", sort(unique(motif_init$tf_list)), envir = env)
  }
  invisible(cfg)
}
