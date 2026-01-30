# utils_tflinks.R
# Utilities for working with TFLink interaction tables stored as
# chunked .rds files in the external episcope-data repo.
#
# Layout of episcope-data (already on GitHub):
#   episcope-data/
#     TFLink/
#       Homo_sapiens/
#         TFLink_Homo_sapiens_chunk_001.rds
#         TFLink_Homo_sapiens_chunk_002.rds
#         ...
#       Mus_musculus/
#         TFLink_Mus_musculus_chunk_001.rds
#         ...
#
# Default behaviour:
#   * On first use, episcope will download the required species chunks from
#     GitHub into a per-user cache directory and reuse them on later calls.
#   * The cache directory is tools::R_user_dir("episcope","cache")/TFLink.
#   * Advanced users can point to a local clone of episcope-data with
#     tflink_set_data_dir() or options(episcope.tflink_data_dir = "...").

.tflink_cache <- new.env(parent = emptyenv())

# -------------------------------------------------------------------------
# Internal helpers (not exported)
# -------------------------------------------------------------------------

.tflink_normalise_species <- function(species) {
  alias <- c(
    hs = "Homo_sapiens",
    mm = "Mus_musculus",
    rn = "Rattus_norvegicus",
    ce = "Caenorhabditis_elegans",
    dr = "Danio_rerio",
    dm = "Drosophila_melanogaster",
    sc = "Saccharomyces_cerevisiae"
  )

  if (!is.character(species) || length(species) != 1L) {
    cli::cli_abort("`species` must be a length-1 character vector.")
  }

  if (species %in% unname(alias)) {
    return(species)
  }

  if (species %in% names(alias)) {
    return(unname(alias[species]))
  }

  cli::cli_abort("Unknown species {.val {species}} for TFLink.")
}

.tflink_default_cache_dir <- function() {
  base <- tools::R_user_dir("episcope", which = "cache")
  dir  <- file.path(base, "TFLink")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  dir
}

.tflink_base_dir <- function() {
  opt <- getOption("episcope.tflink_data_dir", default = NULL)
  if (!is.null(opt)) {
    if (!dir.exists(opt)) {
      cli::cli_abort(
        c(
          "TFLink data directory from option {.code episcope.tflink_data_dir} does not exist.",
          "x" = "Current value: {.path {opt}}",
          "i" = "Unset the option or call {.code tflink_set_data_dir()} with a valid path."
        )
      )
    }
    return(normalizePath(opt, winslash = "/", mustWork = TRUE))
  }
  .tflink_default_cache_dir()
}

.tflink_species_dir <- function(species_norm) {
  file.path(.tflink_base_dir(), species_norm)
}

.tflink_chunk_pattern <- "^TFLink_.*_chunk_\\d{3}\\.rds$"

.tflink_list_chunk_files <- function(species_norm) {
  species_dir <- .tflink_species_dir(species_norm)
  if (!dir.exists(species_dir)) {
    return(character())
  }
  list.files(
    species_dir,
    pattern = .tflink_chunk_pattern,
    full.names = TRUE
  )
}

# Robust downloader: stop at first missing chunk, no noisy warnings for 404s
.tflink_download_species <- function(species_norm,
                                     max_chunks = 999L,
                                     overwrite = FALSE,
                                     quiet = TRUE) {
  base_dir   <- .tflink_base_dir()
  species_dir <- file.path(base_dir, species_norm)
  if (!dir.exists(species_dir)) {
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
  }

  existing <- .tflink_list_chunk_files(species_norm)
  if (length(existing) && !overwrite) {
    return(invisible(existing))
  }

  gh_base <- "https://raw.githubusercontent.com/oncologylab/episcope-data/main/TFLink"
  downloaded <- character()

  for (i in seq_len(max_chunks)) {
    fname <- sprintf("TFLink_%s_chunk_%03d.rds", species_norm, i)
    url   <- sprintf("%s/%s/%s", gh_base, species_norm, fname)
    dest  <- file.path(species_dir, fname)

    status <- suppressWarnings(
      tryCatch(
        utils::download.file(
          url      = url,
          destfile = dest,
          quiet    = quiet,
          mode     = "wb"
        ),
        error = function(e) 1L
      )
    )

    ok <- identical(status, 0L) &&
      file.exists(dest) &&
      is.finite(info <- file.info(dest)$size) &&
      info > 0

    if (!ok) {
      if (file.exists(dest)) {
        unlink(dest)
      }
      # First non-existent or zero-length chunk -> assume end of series
      break
    }

    downloaded <- c(downloaded, dest)
  }

  files <- .tflink_list_chunk_files(species_norm)

  if (!length(files)) {
    cli::cli_abort(
      c(
        "Failed to download TFLink chunks for species {.val {species_norm}}.",
        "i" = "Tried URL prefix: {.url {gh_base}}."
      )
    )
  }

  invisible(files)
}

# -------------------------------------------------------------------------
# Public helpers
# -------------------------------------------------------------------------

#' Configure TFLink data directory
#'
#' Set the base directory where chunked TFLink `.rds` files live.
#' This is only needed if you prefer a local clone of the episcope-data
#' repository. If not set, episcope uses a per-user cache directory and
#' downloads chunks from GitHub on first use.
#'
#' @param path Path to the TFLink data directory (the folder that contains
#'   species subdirectories like `Homo_sapiens`, `Mus_musculus`, etc.).
#'
#' @return Invisibly returns the normalised path.
#' @export
tflink_set_data_dir <- function(path) {
  if (!dir.exists(path)) {
    cli::cli_abort("Directory {.path {path}} does not exist.")
  }
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  options(episcope.tflink_data_dir = path_norm)
  invisible(path_norm)
}

#' Show the current TFLink data directory
#'
#' @return A character string with the resolved TFLink base directory.
#' @export
tflink_data_dir <- function() {
  .tflink_base_dir()
}

#' List available TFLink species
#'
#' Lists the species for which chunked TFLink data is available under
#' the current TFLink data directory or cache.
#'
#' @return A character vector of species folder names.
#' @export
tflink_list_species <- function() {
  base <- .tflink_base_dir()
  dirs <- list.dirs(base, full.names = FALSE, recursive = FALSE)
  dirs[nzchar(dirs)]
}

#' Clear in-memory TFLink cache
#'
#' Clears any TFLink interaction tables cached in memory by [tflink_load()].
#' This does not delete on-disk `.rds` files.
#'
#' @return Invisibly returns `TRUE`.
#' @export
tflink_clear_cache <- function() {
  rm(list = ls(envir = .tflink_cache, all.names = TRUE), envir = .tflink_cache)
  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Preprocessing raw TFLink *.tsv(.gz) into chunked RDS files
# (used when *creating* episcope-data; most users never need this)
# -------------------------------------------------------------------------

#' Pre-process a raw TFLink TSV into chunked RDS files
#'
#' Offline helper for building the episcope-data repository. Converts the raw
#' `TFLink_*_interactions_All_simpleFormat_v1.0.tsv(.gz)` file into a series of
#' compressed `.rds` chunks that are suitable for GitHub storage.
#'
#' @param in_file Path to the raw TFLink TSV (plain or gzipped).
#' @param species Species name or alias (e.g. `"Homo_sapiens"` or `"hs"`).
#' @param out_dir Output directory for this species. It will be created
#'   if missing.
#' @param chunk_rows Number of data rows per chunk. Default `250000L`.
#' @param compress Compression method for [base::saveRDS()]. Default `"xz"`.
#' @param verbose Logical; if `TRUE` prints basic progress.
#'
#' @return Invisibly returns a character vector of output file paths.
#' @export
tflink_preprocess_tsv <- function(in_file,
                                  species,
                                  out_dir,
                                  chunk_rows = 250000L,
                                  compress = "xz",
                                  verbose = TRUE) {
  if (!file.exists(in_file)) {
    cli::cli_abort("Input file {.path {in_file}} does not exist.")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  species_norm <- .tflink_normalise_species(species)
  out_paths <- character()

  cb_fun <- function(x, pos) {
    idx <- ceiling(pos / chunk_rows)
    fname <- sprintf("TFLink_%s_chunk_%03d.rds", species_norm, idx)
    out_path <- file.path(out_dir, fname)

    if (verbose) {
      cli::cli_inform("Writing TFLink chunk {.path {fname}} (rows starting at {pos})")
    }

    base::saveRDS(x, file = out_path, compress = compress)
    out_paths <<- c(out_paths, out_path)
    invisible(NULL)
  }

  readr::read_tsv_chunked(
    file       = in_file,
    callback   = readr::DataFrameCallback$new(cb_fun),
    chunk_size = chunk_rows,
    progress   = verbose,
    col_types  = readr::cols(.default = readr::col_character())
  )

  invisible(out_paths)
}

# -------------------------------------------------------------------------
# Loading + caching
# -------------------------------------------------------------------------

#' Load TFLink interactions for a species (with optional subsetting)
#'
#' Loads chunked TFLink interaction tables for a given species. On first use,
#' if no local chunks are found, episcope will download them into a per-user cache and
#' reuse them for subsequent calls.
#'
#' You can optionally filter by TF name and/or target gene name while reading
#' the chunks to reduce memory usage.
#'
#' @param species Species name or alias. Supported aliases:
#'   \describe{
#'     \item{\code{"hs"}}{Homo sapiens (\code{"Homo_sapiens"})}
#'     \item{\code{"mm"}}{Mus musculus (\code{"Mus_musculus"})}
#'     \item{\code{"rn"}}{Rattus norvegicus (\code{"Rattus_norvegicus"})}
#'     \item{\code{"ce"}}{Caenorhabditis elegans (\code{"Caenorhabditis_elegans"})}
#'     \item{\code{"dr"}}{Danio rerio (\code{"Danio_rerio"})}
#'     \item{\code{"dm"}}{Drosophila melanogaster (\code{"Drosophila_melanogaster"})}
#'     \item{\code{"sc"}}{Saccharomyces cerevisiae (\code{"Saccharomyces_cerevisiae"})}
#'   }
#'   You may also pass the full species identifier (e.g. \code{"Homo_sapiens"}).
#' @param filter_tf Optional character vector of TF names (`Name.TF`) to keep.
#' @param filter_target Optional character vector of target names
#'   (`Name.Target`) to keep.
#' @param select_cols Optional character vector of columns to keep. If `NULL`,
#'   all columns are returned.
#' @param use_cache Logical; if `TRUE` (default) the fully loaded table for
#'   that species is cached in memory and reused on subsequent calls.
#'
#' @return A tibble with TFLink interactions.
#' @export
tflink_load <- function(species,
                        filter_tf = NULL,
                        filter_target = NULL,
                        select_cols = NULL,
                        use_cache = TRUE) {
  species_norm <- .tflink_normalise_species(species)

  cache_key   <- paste0("TFLink_", species_norm)
  has_filters <- !is.null(filter_tf) || !is.null(filter_target) || !is.null(select_cols)

  if (use_cache && !has_filters && exists(cache_key, envir = .tflink_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .tflink_cache, inherits = FALSE))
  }

  files <- .tflink_list_chunk_files(species_norm)

  if (!length(files)) {
    .tflink_download_species(species_norm)
    files <- .tflink_list_chunk_files(species_norm)
  }

  if (!length(files)) {
    cli::cli_abort(
      c(
        "No chunked TFLink {.code .rds} files found for species {.val {species_norm}}.",
        "i" = "Ensure episcope-data is reachable or call {.code tflink_preprocess_tsv()} for a local copy."
      )
    )
  }

  files <- sort(files)
  out_list <- vector("list", length(files))

  for (i in seq_along(files)) {
    dat <- base::readRDS(files[[i]])

    if (!is.null(filter_tf)) {
      if (!"Name.TF" %in% names(dat)) {
        cli::cli_abort("Column {.field Name.TF} not found in TFLink chunk {.path {files[[i]]}}.")
      }
      dat <- dat[dat[["Name.TF"]] %in% filter_tf, , drop = FALSE]
    }

    if (!is.null(filter_target)) {
      if (!"Name.Target" %in% names(dat)) {
        cli::cli_abort("Column {.field Name.Target} not found in TFLink chunk {.path {files[[i]]}}.")
      }
      dat <- dat[dat[["Name.Target"]] %in% filter_target, , drop = FALSE]
    }

    if (!is.null(select_cols)) {
      keep <- intersect(select_cols, names(dat))
      if (!length(keep)) {
        cli::cli_abort(
          c(
            "None of the requested columns present in TFLink data.",
            "x" = "Requested: {.val {select_cols}}"
          )
        )
      }
      dat <- dat[, keep, drop = FALSE]
    }

    out_list[[i]] <- dat
  }

  res <- dplyr::bind_rows(out_list)

  if (use_cache && !has_filters) {
    assign(cache_key, res, envir = .tflink_cache)
  }

  res
}





# ## Base directory where the original TFLink TSV/TSV.GZ live
# base_dir <- "~/tflink"
# ## species -> file name
# tflink_files <- c(
#   "Caenorhabditis_elegans"      = "TFLink_Caenorhabditis_elegans_interactions_All_simpleFormat_v1.0.tsv",
#   "Danio_rerio"                 = "TFLink_Danio_rerio_interactions_All_simpleFormat_v1.0.tsv",
#   "Drosophila_melanogaster"     = "TFLink_Drosophila_melanogaster_interactions_All_simpleFormat_v1.0.tsv",
#   "Homo_sapiens"                = "TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz",
#   "Mus_musculus"                = "TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv.gz",
#   "Rattus_norvegicus"           = "TFLink_Rattus_norvegicus_interactions_All_simpleFormat_v1.0.tsv",
#   "Saccharomyces_cerevisiae"    = "TFLink_Saccharomyces_cerevisiae_interactions_All_simpleFormat_v1.0.tsv"
# )
# ## Loop over species and preprocess each one
# for (species in names(tflink_files)) {
#   in_file <- file.path(base_dir, tflink_files[[species]])
#   # chunks will be written to
#   out_dir <- file.path(base_dir, species)
#   cli::cli_inform("Processing {species} from {.path {basename(in_file)}}")
#   tflink_preprocess_tsv(
#     in_file    = in_file,
#     species    = species,
#     out_dir    = out_dir,
#     chunk_rows = 250000L,
#     compress   = "xz",
#     verbose    = TRUE
#   )
# }
# cli::cli_inform("All TFLink species processed.")

# hs_hnf <- tflink_load("hs",
#   filter_tf  = c("HNF1A", "HNF4A"),
#   select_cols = c("Name.TF", "Name.Target"))
hs_hnf <- tflink_load("hs", filter_tf  = c("HNF1A"))
