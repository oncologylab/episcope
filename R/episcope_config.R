#' Load an episcope YAML config into the global environment
#'
#' Reads a YAML file and assigns each top-level key as a variable in
#' the global environment (e.g., `db`, `threshold_tf_expr`, etc.).
#'
#' @param path Character path to a YAML file.
#' @return (Invisibly) the parsed list.
#' @examples
#' \dontrun{
#' load_episcope_config("episcope_grn.yaml")
#' # Now `db`, `threshold_tf_expr`, `link_score_threshold`, ... are available globally.
#' }
#' @export
load_episcope_config <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("Config file not found: {path}")
  }
  cfg <- yaml::read_yaml(path)
  # Assign everything at top level directly into the global environment
  list2env(cfg, envir = .GlobalEnv)
  invisible(cfg)
}
