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

.pathway_backend <- function(backend = NULL) {
  if (is.null(backend) || !length(backend) || is.na(backend[[1L]]) || !nzchar(as.character(backend)[[1L]])) {
    opt <- getOption("craftgrn.pathway_backend", NULL)
    env <- Sys.getenv("CRAFTGRN_PATHWAY_BACKEND", unset = "")
    backend <- if (!is.null(opt) && length(opt) && nzchar(as.character(opt)[[1L]])) opt else env
  }
  if (is.null(backend) || !length(backend) || is.na(backend[[1L]]) || !nzchar(as.character(backend)[[1L]])) {
    backend <- "enrichly"
  }
  backend <- tolower(as.character(backend)[[1L]])
  backend <- gsub("[^a-z0-9]+", "_", backend)
  if (backend %in% c("local", "local_first", "enrichly_local")) backend <- "enrichly"
  if (backend %in% c("web", "api", "enrichr_api")) backend <- "enrichr"
  if (!(backend %in% c("enrichly", "enrichr"))) {
    .log_abort("`pathway_backend` must be either 'enrichly' or 'enrichr'.")
  }
  backend
}

.pathway_backend_available <- function(backend = NULL) {
  backend <- .pathway_backend(backend)
  if (identical(backend, "enrichly") && .optional_namespace_available("enrichly")) {
    return(TRUE)
  }
  requireNamespace("enrichR", quietly = TRUE)
}

.optional_namespace_available <- function(pkg) {
  pkg <- as.character(pkg)[[1L]]
  tryCatch(
    requireNamespace(pkg, quietly = TRUE),
    error = function(e) FALSE
  )
}

.enrichr_cache_key <- function(genes, dbs, site = "Enrichr", backend = NULL) {
  .assert_pkg("digest")
  genes <- sort(unique(as.character(genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  dbs <- sort(unique(as.character(dbs)))
  dbs <- dbs[!is.na(dbs) & nzchar(dbs)]
  backend <- .pathway_backend(backend)
  backend_version <- if (identical(backend, "enrichly")) "enrichly_enrichr_standard_bg20000_v1" else "enrichr_api"
  digest::digest(
    list(site = site, backend = backend, backend_version = backend_version, dbs = dbs, genes = genes),
    algo = "xxhash64",
    serialize = TRUE
  )
}

.enrichr_cache_path <- function(cache_dir, genes, dbs, site = "Enrichr", backend = NULL) {
  if (is.null(cache_dir) || !nzchar(as.character(cache_dir)[[1L]])) {
    return(NULL)
  }
  file.path(as.character(cache_dir)[[1L]], paste0(.enrichr_cache_key(genes, dbs, site = site, backend = backend), ".rds"))
}

.module3_default_enrichr_cache_dir <- function(out_dir) {
  out_dir <- as.character(out_dir)[[1L]]
  parent <- dirname(out_dir)
  if (grepl("^K[0-9]+$", basename(parent))) {
    return(file.path(dirname(parent), "cache", "enrichr"))
  }
  file.path(parent, "cache", "enrichr")
}

.normalize_enrichr_n_cores <- function(n_cores = 1L) {
  if (is.null(n_cores) || !length(n_cores) || is.na(n_cores[[1L]])) {
    env <- Sys.getenv("CRAFTGRN_ENRICHR_MAX_CORES", unset = "")
    n_cores <- if (nzchar(env)) suppressWarnings(as.integer(env)) else 1L
  }
  n_cores <- suppressWarnings(as.integer(n_cores[[1L]]))
  if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L
  max(1L, min(n_cores, .available_cores(logical = TRUE)))
}

.enrichly_db_cache_dir <- function(cache_dir = NULL) {
  opt <- getOption("craftgrn.enrichly.db_cache", NULL)
  env <- Sys.getenv("CRAFTGRN_ENRICHLY_DB_CACHE", unset = "")
  if (!is.null(opt) && length(opt) && nzchar(as.character(opt)[[1L]])) {
    return(as.character(opt)[[1L]])
  }
  if (nzchar(env)) return(env)
  file.path(tools::R_user_dir("craftgrn", which = "cache"), "enrichly_databases")
}

.enrichly_result_to_enrichr_list <- function(x) {
  .assert_pkg("data.table")
  dt <- data.table::as.data.table(x)
  required <- c("database", "term", "p_value", "adjusted_p_value")
  if (!all(required %in% names(dt))) {
    .log_abort("enrichly result is missing required columns.")
  }
  if (!nrow(dt)) {
    return(list())
  }
  if (!("overlap" %in% names(dt))) data.table::set(dt, j = "overlap", value = NA_character_)
  if (!("overlap_genes" %in% names(dt))) data.table::set(dt, j = "overlap_genes", value = NA_character_)
  if (!("odds_ratio" %in% names(dt))) data.table::set(dt, j = "odds_ratio", value = NA_real_)
  if (!("combined_score" %in% names(dt))) data.table::set(dt, j = "combined_score", value = NA_real_)
  data.table::setorderv(dt, c("database", "adjusted_p_value", "p_value", "term"))
  out <- lapply(split(dt, dt[["database"]]), function(z) {
    data.frame(
      Term = as.character(z[["term"]]),
      Overlap = as.character(z[["overlap"]]),
      P.value = as.numeric(z[["p_value"]]),
      Adjusted.P.value = as.numeric(z[["adjusted_p_value"]]),
      Odds.Ratio = as.numeric(z[["odds_ratio"]]),
      Combined.Score = as.numeric(z[["combined_score"]]),
      Genes = as.character(z[["overlap_genes"]]),
      stringsAsFactors = FALSE
    )
  })
  out[sort(names(out))]
}

.run_enrichly_local <- function(genes, dbs, cache_dir = NULL, universe = NULL) {
  if (!.optional_namespace_available("enrichly")) {
    return(NULL)
  }
  db_cache_dir <- .enrichly_db_cache_dir(cache_dir)
  enrichly_download <- getExportedValue("enrichly", "enrichly_download")
  enrichly_load <- getExportedValue("enrichly", "enrichly_load")
  enrichly_enrich <- getExportedValue("enrichly", "enrichly_enrich")
  manifest <- enrichly_download(
    databases = dbs,
    cache_dir = db_cache_dir,
    overwrite = FALSE,
    verbose = FALSE
  )
  db <- enrichly_load(manifest$path, databases = dbs)
  res <- enrichly_enrich(
    genes = genes,
    db = db,
    query_id = "query",
    universe = universe,
    background_size = 20000L
  )
  .enrichly_result_to_enrichr_list(res)
}

.run_enrichr_cached <- function(genes,
                                dbs,
                                sleep_time = 0,
                                cache_dir = NULL,
                                site = "Enrichr",
                                backend = NULL,
                                universe = NULL) {
  sleep_time <- .normalize_enrichr_sleep_time(sleep_time)
  backend <- .pathway_backend(backend)
  cache_path <- .enrichr_cache_path(cache_dir, genes, dbs, site = site, backend = backend)
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (is.list(cached)) return(cached)
  }
  res <- NULL
  if (identical(backend, "enrichly")) {
    res <- tryCatch(
      .run_enrichly_local(genes = genes, dbs = dbs, cache_dir = cache_dir, universe = universe),
      error = function(e) NULL
    )
  }
  if (!is.list(res)) {
    if (!requireNamespace("enrichR", quietly = TRUE)) {
      .log_abort("Pathway enrichment requires either {.pkg enrichly} for local analysis or {.pkg enrichR} for web API analysis.")
    }
    .ensure_enrichr_ready(site = site, verbose = FALSE)
    res <- enrichR::enrichr(genes, dbs, sleepTime = sleep_time)
  }
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
