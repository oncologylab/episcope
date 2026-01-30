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
