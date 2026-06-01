# File: utils_demo_data.R
# Author: Yaoxiang Li
# Created: 2026-05-29
# Updated: 2026-06-01
#
# Purpose:
# Provide small helpers for listing and downloading external CraftGRN demo data.
#
# Notes:
# - Demo data are intentionally kept out of the source package and git history.
# - Network access is explicit and user-triggered.

#' CraftGRN demo data helpers
#'
#' @noRd
NULL

.craftgrn_demo_data_registry <- function() {
  list()
}

.craftgrn_demo_data_empty <- function() {
  tibble::tibble(
    name = character(),
    title = character(),
    version = character(),
    release_tag = character(),
    file = character(),
    project_dir = character(),
    url = character(),
    md5 = character(),
    size = character(),
    description = character()
  )
}

#' Return metadata for configured external CraftGRN demo data
#'
#' @param demo Optional demo bundle name. No external demo bundle is currently
#'   configured.
#'
#' @return A data frame with the bundle URL, checksum, archive file name, and
#'   extracted project directory name. When no demo bundle is configured, the
#'   returned data frame has zero rows.
#' @export
craftgrn_demo_data_info <- function(demo = NULL) {
  registry <- .craftgrn_demo_data_registry()
  if (!length(registry)) {
    return(.craftgrn_demo_data_empty())
  }
  choices <- names(registry)
  if (is.null(demo)) {
    demo <- choices[[1L]]
  }
  if (!is.character(demo) || length(demo) != 1L || is.na(demo) || !nzchar(demo)) {
    .log_abort("`demo` must be a non-empty character scalar.")
  }
  if (!demo %in% choices) {
    .log_abort("Unknown CraftGRN demo data bundle: {demo}.")
  }
  x <- registry[[demo]]
  tibble::tibble(
    name = x$name,
    title = x$title,
    version = x$version,
    release_tag = x$release_tag,
    file = x$file,
    project_dir = x$project_dir,
    url = x$url,
    md5 = x$md5,
    size = x$size,
    description = x$description
  )
}

#' Download and unpack configured external CraftGRN demo data
#'
#' Downloads a processed demo data archive from a configured external source,
#' verifies its MD5 checksum by default, extracts it, and returns the extracted
#' project directory. Demo bundles are external to the R package so package
#' installation remains small and CRAN-friendly. No external demo bundle is
#' currently configured.
#'
#' If the download fails, inspect `craftgrn_demo_data_info()` and download the
#' configured asset manually. If checksum verification fails, rerun with
#' `overwrite = TRUE` to replace a stale or partial archive. The extracted
#' project uses `base_dir: "."`, so pass the returned directory or its
#' project config path directly to package functions after moving the folder.
#'
#' @param destdir Directory where the archive should be downloaded and unpacked.
#' @param demo Optional demo bundle name. No external demo bundle is currently
#'   configured.
#' @param overwrite Logical; if `TRUE`, download the archive again and replace an
#'   existing extracted project directory.
#' @param checksum Logical; if `TRUE`, verify the downloaded archive MD5.
#' @param verbose Logical; if `TRUE`, emit concise status messages.
#'
#' @return The normalized path to the extracted demo project directory.
#' @examples
#' craftgrn_demo_data_info()
#' @export
download_craftgrn_demo_data <- function(destdir = ".",
                                        demo = NULL,
                                        overwrite = FALSE,
                                        checksum = TRUE,
                                        verbose = TRUE) {
  info <- craftgrn_demo_data_info(demo)
  if (!nrow(info)) {
    .log_abort("No external CraftGRN demo data bundle is currently configured.")
  }

  stopifnot(is.logical(overwrite), length(overwrite) == 1L, !is.na(overwrite))
  stopifnot(is.logical(checksum), length(checksum) == 1L, !is.na(checksum))
  stopifnot(is.logical(verbose), length(verbose) == 1L, !is.na(verbose))

  if (!is.character(destdir) || length(destdir) != 1L || !nzchar(destdir)) {
    .log_abort("`destdir` must be a non-empty character scalar.")
  }
  destdir <- path.expand(destdir)
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)
  destdir <- normalizePath(destdir, winslash = "/", mustWork = TRUE)

  archive <- file.path(destdir, info$file[[1L]])
  project_dir <- file.path(destdir, info$project_dir[[1L]])

  if (dir.exists(project_dir) && !isTRUE(overwrite)) {
    if (isTRUE(verbose)) {
      .log_inform("Using existing CraftGRN demo data at {.path {project_dir}}.")
    }
    return(normalizePath(project_dir, winslash = "/", mustWork = TRUE))
  }

  if (!file.exists(archive) || isTRUE(overwrite)) {
    if (isTRUE(verbose)) {
      .log_inform("Downloading CraftGRN demo data: {info$title[[1L]]}.")
    }
    status <- utils::download.file(
      url = info$url[[1L]],
      destfile = archive,
      mode = "wb",
      quiet = !isTRUE(verbose)
    )
    if (!identical(status, 0L)) {
      .log_abort("Failed to download CraftGRN demo data from {info$url[[1L]]}.")
    }
  }

  if (isTRUE(checksum)) {
    observed <- unname(as.character(tools::md5sum(archive)))
    if (!identical(tolower(observed), tolower(info$md5[[1L]]))) {
      .log_abort(c(
        "Downloaded CraftGRN demo data checksum does not match.",
        i = paste0("Expected: ", info$md5[[1L]]),
        i = paste0("Observed: ", observed)
      ))
    }
  }

  if (dir.exists(project_dir) && isTRUE(overwrite)) {
    unlink(project_dir, recursive = TRUE, force = TRUE)
  }
  utils::untar(archive, exdir = destdir)

  if (!dir.exists(project_dir)) {
    .log_abort("Demo archive did not extract expected project directory: {.path {project_dir}}.")
  }
  if (isTRUE(verbose)) {
    .log_inform("CraftGRN demo data ready: {.path {project_dir}}.")
  }
  normalizePath(project_dir, winslash = "/", mustWork = TRUE)
}
