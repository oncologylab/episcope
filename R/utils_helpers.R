#' helpers 
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
NULL

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
