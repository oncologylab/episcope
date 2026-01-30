#' Validate config values
#'
#' Ensures required config keys (e.g. thresholds and db) exist in the chosen
#' environment before running pipelines.
#'
#' @param required Character vector of required variable names.
#' @param numeric_required Character vector of required numeric variable names.
#' @param env Environment to check. Default \code{.GlobalEnv}.
#'
#' @return \code{TRUE} invisibly when validation passes.
#' @export
validate_config <- function(
    required = c(
      "db",
      "ref_genome",
      "threshold_tf_expr",
      "threshold_gene_expr",
      "threshold_fp_score",
      "threshold_link_score",
      "threshold_fp_tf_corr_p",
      "threshold_fp_tf_corr_r",
      "threshold_fp_gene_corr_p",
      "threshold_fp_gene_corr_abs_r",
      "threshold_atac_gene_corr_p",
      "threshold_atac_gene_corr_abs_r",
      "link_score_threshold",
      "fp_score_threshold"
    ),
    numeric_required = c(
      "threshold_tf_expr",
      "threshold_gene_expr",
      "threshold_fp_score",
      "threshold_link_score",
      "threshold_fp_tf_corr_p",
      "threshold_fp_tf_corr_r",
      "threshold_fp_gene_corr_p",
      "threshold_fp_gene_corr_abs_r",
      "threshold_atac_gene_corr_p",
      "threshold_atac_gene_corr_abs_r",
      "link_score_threshold",
      "fp_score_threshold"
    ),
    env = .GlobalEnv
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
#' Expands \code{~} and environment variables in configured paths.
#'
#' @param keys Character vector of config keys to normalize.
#' @param env Environment to update. Default \code{.GlobalEnv}.
#'
#' @return Invisibly returns the updated values (named character vector).
#' @export
normalize_config_paths <- function(
    keys = c("base_dir", "fp_root_dir", "benchmark_tfbs_dir", "pathway_dir"),
    env = .GlobalEnv
) {
  existing <- keys[vapply(keys, function(nm) exists(nm, envir = env, inherits = FALSE), logical(1))]
  if (!length(existing)) return(invisible(character(0)))

  vals <- vapply(existing, function(nm) {
    val <- get(nm, envir = env)
    if (!is.character(val) || !length(val) || !nzchar(val[1])) return(val)
    path.expand(val[1])
  }, character(1))

  list2env(as.list(vals), envir = env)
  invisible(vals)
}

#' Load an episcope YAML config into the global environment
#'
#' Reads a YAML file and assigns each top-level key as a variable in
#' the global environment (e.g., `db`, `threshold_tf_expr`, etc.).
#' Also runs standard config initialization helpers when available.
#'
#' @param path Character path to a YAML file.
#' @return (Invisibly) the parsed list.
#' @examples
#' \dontrun{
#' load_config("episcope_grn.yaml")
#' # Now `db`, `threshold_tf_expr`, `link_score_threshold`, ... are available globally.
#' }
#' @export
load_config <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("Config file not found: {path}")
  }
  cfg <- yaml::read_yaml(path)
  # Assign everything at top level directly into the global environment
  list2env(cfg, envir = .GlobalEnv)
  if (exists("validate_config", mode = "function")) {
    validate_config()
  }
  if (exists("normalize_config_paths", mode = "function")) {
    normalize_config_paths()
  }
  if (exists("init_motif_db", mode = "function") && exists("db", envir = .GlobalEnv)) {
    ref_genome_val <- if (exists("ref_genome", envir = .GlobalEnv)) get("ref_genome", envir = .GlobalEnv) else NULL
    init_formals <- names(formals(init_motif_db))
    if (!is.null(init_formals) && "ref_genome" %in% init_formals) {
      motif_init <- init_motif_db(get("db", envir = .GlobalEnv), ref_genome = ref_genome_val)
    } else {
      motif_init <- init_motif_db(get("db", envir = .GlobalEnv))
    }
    assign("motif_init", motif_init, envir = .GlobalEnv)
    assign("motif_db", motif_init$motif_db, envir = .GlobalEnv)
    assign("tf_list_all", motif_init$tf_list, envir = .GlobalEnv)
    assign("tf_list", sort(unique(motif_init$tf_list)), envir = .GlobalEnv)
  }
  invisible(cfg)
}
