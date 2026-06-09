# File: utils_helpers.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide small internal helpers that do not fit a more specific shared utility
# file.
#
# Inputs:
# - helper-specific arguments such as external service names or logging hooks
#
# Outputs:
# - internal setup state or logical flags for optional integrations
#
# Notes:
# - Keep this file generic and reusable across modules.
# - Move code into a more specific utility or module file once ownership is
#   clear.

#' Internal miscellaneous helpers
#'
#' @noRd
NULL

#' @noRd
.assert_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .log_abort("Package {.pkg {pkg}} is required but not installed.")
  }
  invisible(TRUE)
}

.normalize_enrichr_sleep_time <- function(sleep_time = 0) {
  sleep_time <- suppressWarnings(as.numeric(sleep_time[[1L]]))
  if (!is.finite(sleep_time) || sleep_time < 0) sleep_time <- 0
  sleep_time
}

.enrichr_cache_key <- function(genes, dbs, site = "Enrichr") {
  .assert_pkg("digest")
  genes <- sort(unique(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  dbs <- sort(unique(as.character(dbs)))
  dbs <- dbs[!is.na(dbs) & nzchar(dbs)]
  digest::digest(
    list(site = site, dbs = dbs, genes = genes),
    algo = "xxhash64",
    serialize = TRUE
  )
}

.enrichr_cache_path <- function(cache_dir, genes, dbs, site = "Enrichr") {
  if (is.null(cache_dir) || !nzchar(as.character(cache_dir)[[1L]])) {
    return(NULL)
  }
  file.path(as.character(cache_dir)[[1L]], paste0(.enrichr_cache_key(genes, dbs, site = site), ".rds"))
}

.module3_default_enrichr_cache_dir <- function(out_dir) {
  out_dir <- as.character(out_dir)[[1L]]
  parent <- dirname(out_dir)
  if (grepl("^K[0-9]+$", basename(parent))) {
    return(file.path(dirname(parent), "cache", "enrichr"))
  }
  file.path(parent, "cache", "enrichr")
}

.run_enrichr_cached <- function(genes,
                                dbs,
                                sleep_time = 0,
                                cache_dir = NULL,
                                site = "Enrichr") {
  sleep_time <- .normalize_enrichr_sleep_time(sleep_time)
  cache_path <- .enrichr_cache_path(cache_dir, genes, dbs, site = site)
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.list(cached)) return(cached)
  }
  res <- enrichR::enrichr(genes, dbs, sleepTime = sleep_time)
  if (!is.null(cache_path)) {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(cache_path, ".tmp.", Sys.getpid())
    tryCatch(
      {
        saveRDS(res, tmp)
        file.rename(tmp, cache_path)
      },
      error = function(e) {
        if (file.exists(tmp)) unlink(tmp, force = TRUE)
      }
    )
  }
  res
}

# Internal Enrichr setup helper (no library(enrichR) needed in package code).
.ensure_enrichr_ready <- function(site = "Enrichr", verbose = TRUE, log_fun = NULL) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    if (isTRUE(verbose) && is.function(log_fun)) log_fun("enrichR not installed.")
    return(FALSE)
  }

  .log_local <- function(msg) {
    if (is.function(log_fun)) log_fun(msg)
    invisible(NULL)
  }

  base_addr <- getOption("enrichR.sites.base.address")
  if (is.null(base_addr) || !nzchar(base_addr)) {
    base_addr <- "https://maayanlab.cloud/"
    options(enrichR.sites.base.address = base_addr)
  }
  sites <- getOption("enrichR.sites")
  if (is.null(sites) || !length(sites)) {
    options(enrichR.sites = site)
  }
  base <- getOption("enrichR.base.address")
  if (is.null(base) || !nzchar(base)) {
    options(enrichR.base.address = paste0(base_addr, site, "/"))
  }
  live <- getOption("enrichR.live")
  if (is.null(live) || !is.logical(live) || !length(live)) {
    options(enrichR.live = FALSE)
  }
  quiet <- getOption("enrichR.quiet")
  if (is.null(quiet) || !is.logical(quiet) || !length(quiet)) {
    options(enrichR.quiet = !isTRUE(verbose))
  }

  ok <- TRUE
  tryCatch(
    {
      utils::capture.output(
        suppressMessages(enrichR::setEnrichrSite(site)),
        type = "output"
      )
      .log_local(sprintf("Enrichr site set to '%s'.", site))
    },
    error = function(e) {
      ok <<- FALSE
      options(
        enrichR.sites.base.address = base_addr,
        enrichR.sites = site,
        enrichR.base.address = paste0(base_addr, site, "/")
      )
      .log_local(sprintf("Enrichr site init fallback applied after error: %s", conditionMessage(e)))
    }
  )

  # If setEnrichrSite failed, keep going with fallback options.
  TRUE
}
