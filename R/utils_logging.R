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
