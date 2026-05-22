# File: utils_logging.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide shared package logging helpers with a consistent timestamp prefix.
#
# Inputs:
# - message text and optional cli interpolation arguments
#
# Outputs:
# - concise user-facing log, warning, and error messages
#
# Notes:
# - Keep this file generic and reusable across modules.
# - Do not add Module-specific workflow logic here.

#' Internal logging helpers with timestamp prefix
#'
#' @noRd
.log_stamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

#' @noRd
.log_inform <- function(msg, ..., .envir = parent.frame()) {
  if (is.character(msg)) {
    msg[1] <- paste0("[", .log_stamp(), "] ", msg[1])
  }
  cli::cli_inform(msg, ..., .envir = .envir)
}

#' @noRd
.log_warn <- function(msg, ..., .envir = parent.frame()) {
  if (is.character(msg)) {
    msg[1] <- paste0("[", .log_stamp(), "] ", msg[1])
  }
  cli::cli_warn(msg, ..., .envir = .envir)
}

#' @noRd
.log_abort <- function(msg, ..., .envir = parent.frame()) {
  if (is.character(msg)) {
    msg[1] <- paste0("[", .log_stamp(), "] ", msg[1])
  }
  cli::cli_abort(msg, ..., .envir = .envir)
}

.progress_env <- new.env(parent = emptyenv())
.progress_env$id <- NULL
.progress_env$verbose <- FALSE

#' @noRd
.progress_start <- function(total = NULL, message = NULL, verbose = TRUE) {
  .progress_env$verbose <- isTRUE(verbose)
  .progress_env$id <- NULL
  if (!isTRUE(verbose)) return(invisible(NULL))
  total <- if (is.null(total)) NA_integer_ else as.integer(total)[[1L]]
  if (!is.finite(total) || total < 0L) total <- NA_integer_
  message <- if (is.null(message) || !nzchar(message)) "progress" else as.character(message)[[1L]]
  .progress_env$id <- cli::cli_progress_bar(
    name = message,
    total = total,
    clear = TRUE
  )
  invisible(.progress_env$id)
}

#' @noRd
.progress_update <- function(increment = 1L) {
  if (!isTRUE(.progress_env$verbose) || is.null(.progress_env$id)) return(invisible(NULL))
  increment <- as.integer(increment)[[1L]]
  if (!is.finite(increment) || increment < 0L) increment <- 1L
  cli::cli_progress_update(id = .progress_env$id, inc = increment)
  invisible(NULL)
}

#' @noRd
.progress_done <- function() {
  if (isTRUE(.progress_env$verbose) && !is.null(.progress_env$id)) {
    cli::cli_progress_done(id = .progress_env$id)
  }
  .progress_env$id <- NULL
  .progress_env$verbose <- FALSE
  invisible(NULL)
}

#' @noRd
.with_progress <- function(total = NULL, message = NULL, verbose = TRUE, code) {
  .progress_start(total = total, message = message, verbose = verbose)
  on.exit(.progress_done(), add = TRUE)
  force(code)
}
