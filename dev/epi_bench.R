set.seed(2025)
n <- 10000
p <- 100
X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y <- X %*% stats::rnorm(p) + stats::rnorm(n)
check_for_equal_coefs <- function(values) {
  tol <- 1e-12
  max_error <- max(
    c(
      abs(values[[1]] - values[[2]]),
      abs(values[[2]] - values[[3]]),
      abs(values[[1]] - values[[3]])
    )
  )
  max_error < tol
}
epi_lm_solver <- function(t = 10,
                          report_every_sec = 300) {
  start_time   <- Sys.time()
  end_time     <- start_time + t * 3600
  last_report  <- start_time
  iter         <- 0L
  while (Sys.time() < end_time) {
    b1 <- stats::lm(y ~ X + 0)$coef
    XtX <- base::crossprod(X)        # t(X) %*% X
    Xty <- base::crossprod(X, y)     # t(X) %*% y
    b2 <- base::solve(XtX) %*% Xty
    b3 <- base::solve(XtX, Xty)
    ok <- check_for_equal_coefs(list(b1, b2, b3))
    if (!ok) {
      stop("Coefficient mismatch between methods; aborting stress test.")
    }

    iter <- iter + 1L

    now <- Sys.time()
    if (as.numeric(difftime(now, last_report, units = "secs")) >= report_every_sec) {
      elapsed_sec <- as.numeric(difftime(now, start_time, units = "secs"))
      cat(
        sprintf(
          "[%s] iterations: %d | elapsed: %.2f h\n",
          format(now, "%Y-%m-%d %H:%M:%S"),
          iter,
          elapsed_sec / 3600
        )
      )
      last_report <- now
    }
  }

  total_elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(
    sprintf(
      "Done. Total iterations: %d | total elapsed: %.2f h\n",
      iter,
      total_elapsed_sec / 3600
    )
  )

  invisible(list(iterations = iter, hours = total_elapsed_sec / 3600))
}


epi_lm_solver(10)
