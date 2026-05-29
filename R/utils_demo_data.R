# File: utils_demo_data.R
# Author: Yaoxiang Li
# Created: 2026-05-29
# Updated: 2026-05-29
#
# Purpose:
# Provide small helpers for downloading external CraftGRN demo data bundles.
#
# Notes:
# - Demo data are intentionally kept out of the source package and git history.
# - Network access is explicit and user-triggered.

#' CraftGRN demo data helpers
#'
#' @noRd
NULL

.craftgrn_demo_data_registry <- function() {
  list(
    nutrient_stress_chr22 = list(
      name = "nutrient_stress_chr22",
      title = "Nutrient stress chr22 processed demo",
      version = "0.1.0",
      release_tag = "demo-data-v0.1.0",
      file = "nutrient_stress_strict_JASPAR2024_chr22_demo_inputs.tar.gz",
      project_dir = "nutrient_stress_strict_JASPAR2024_chr22_demo",
      url = paste0(
        "https://github.com/oncologylab/craftgrn/releases/download/",
        "demo-data-v0.1.0/",
        "nutrient_stress_strict_JASPAR2024_chr22_demo_inputs.tar.gz"
      ),
      md5 = "fa754783186b0e5119387b0405652331",
      size = "5.2 MB",
      description = paste(
        "Input-only chr22 subset with matched ATAC/RNA tables and compact",
        "aligned JASPAR2024 footprint cache. Raw TOBIAS files and generated",
        "analysis outputs are not included."
      )
    )
  )
}

#' Return metadata for an external CraftGRN demo data bundle
#'
#' @param demo Demo bundle name. Currently `"nutrient_stress_chr22"`.
#'
#' @return A data frame with the bundle URL, checksum, archive file name, and
#'   extracted project directory name.
#' @export
craftgrn_demo_data_info <- function(demo = c("nutrient_stress_chr22")) {
  demo <- match.arg(demo)
  x <- .craftgrn_demo_data_registry()[[demo]]
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

#' Download and unpack an external CraftGRN demo data bundle
#'
#' Downloads a processed demo data archive from the CraftGRN GitHub Releases
#' page, verifies its MD5 checksum by default, extracts it, and returns the
#' extracted project directory. The demo bundle is external to the R package so
#' package installation remains small and CRAN-friendly.
#'
#' @param destdir Directory where the archive should be downloaded and unpacked.
#' @param demo Demo bundle name. Currently `"nutrient_stress_chr22"`.
#' @param overwrite Logical; if `TRUE`, download the archive again and replace an
#'   existing extracted project directory.
#' @param checksum Logical; if `TRUE`, verify the downloaded archive MD5.
#' @param verbose Logical; if `TRUE`, emit concise status messages.
#'
#' @return The normalized path to the extracted demo project directory.
#' @examples
#' \dontrun{
#' demo_dir <- download_craftgrn_demo_data(destdir = tempdir())
#' config <- file.path(demo_dir, "project.yaml")
#' omics <- load_prep_multiomic_data(
#'   config = config,
#'   label_col = "strict_match_rna",
#'   do_preprocess = FALSE,
#'   verbose = TRUE
#' )
#' }
#' @export
download_craftgrn_demo_data <- function(destdir = ".",
                                        demo = c("nutrient_stress_chr22"),
                                        overwrite = FALSE,
                                        checksum = TRUE,
                                        verbose = TRUE) {
  demo <- match.arg(demo)
  info <- craftgrn_demo_data_info(demo)

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
