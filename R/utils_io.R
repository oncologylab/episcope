#' Save a multi-omic data object to disk
#'
#' @param omics_data A multi-omic data list (e.g., output of load_prep_multiomic_data()).
#' @param file Optional full path to an RDS file. If NULL, uses out_dir/db/prefix.
#' @param out_dir Output directory used when file is NULL.
#' @param db Optional database tag appended to the filename when file is NULL.
#' @param prefix Filename prefix used when file is NULL.
#' @param compress Compression passed to saveRDS().
#' @param verbose Emit status messages.
#'
#' @return Path to the written file (invisible).
#' @export
save_omics_data <- function(
    omics_data,
    file = NULL,
    out_dir = NULL,
    db = NULL,
    prefix = "omics_data",
    compress = "xz",
    verbose = TRUE
) {
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")

  if (is.null(file)) {
    if (!is.character(out_dir) || !nzchar(out_dir)) {
      .log_abort("`out_dir` must be provided when `file` is NULL.")
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    suffix <- if (is.character(db) && nzchar(db)) paste0("_", db) else ""
    file <- file.path(out_dir, paste0(prefix, suffix, ".rds"))
  }

  saveRDS(omics_data, file = file, compress = compress)
  if (isTRUE(verbose)) .log_inform("Saved omics data: {file}")
  invisible(file)
}

#' Load a multi-omic data object from disk
#'
#' @param file Path to an RDS file produced by save_omics_data().
#' @param verbose Emit status messages.
#'
#' @return The loaded multi-omic data list.
#' @export
load_omics_data <- function(file, verbose = TRUE) {
  if (!is.character(file) || !nzchar(file)) .log_abort("`file` must be a non-empty path.")
  if (!file.exists(file)) .log_abort("Omics data file not found: {file}")
  out <- readRDS(file)
  if (isTRUE(verbose)) .log_inform("Loaded omics data: {file}")
  out
}
