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
