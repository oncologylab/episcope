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

#' Validate config values
#'
#' Ensures required config keys (e.g. thresholds and db) exist in the chosen
#' environment before running pipelines.
#'
#' @param required Character vector of required variable names.
#' @param numeric_required Character vector of required numeric variable names.
#' @param env Environment to check. Defaults to the internal CraftGRN config state.
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
    env = .craftgrn_state
) {
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
      "pathway_dir",
      "genehancer_path",
      "loop_path",
      "preprocess_dir",
      "sample_metadata",
      "metadata_path",
      "sample_metadata_path",
      "atac_master",
      "atac_data_path",
      "atac_master_path",
      "rna_mapped",
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
