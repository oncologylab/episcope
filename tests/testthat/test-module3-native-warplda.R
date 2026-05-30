make_native_warplda_dtm <- function() {
  skip_if_not_installed("Matrix")
  counts <- matrix(1L, nrow = 8L, ncol = 12L)
  counts[1:4, 1:4] <- 25L
  counts[1:4, 5:8] <- 2L
  counts[5:8, 1:4] <- 2L
  counts[5:8, 5:8] <- 25L
  counts[, 9:12] <- 1L
  rownames(counts) <- paste0("doc", seq_len(nrow(counts)))
  colnames(counts) <- paste0("term", seq_len(ncol(counts)))
  methods::as(Matrix::Matrix(counts, sparse = TRUE), "dgCMatrix")
}

expect_probability_matrix <- function(x, n_row, n_col) {
  expect_equal(dim(x), c(n_row, n_col))
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0))
  expect_equal(as.numeric(rowSums(x)), rep(1, n_row), tolerance = 1e-8)
}

test_that("native WarpLDA returns normalized topic matrices", {
  dtm <- make_native_warplda_dtm()

  fit <- fit_warplda_one(
    dtm,
    K = 2L,
    iterations = 20L,
    alpha = 0.5,
    beta = 0.1,
    seed = 11L,
    n_check_convergence = 0L,
    n_iter_inference = 2L,
    n_threads = 1L,
    progressbar = FALSE
  )

  expect_identical(fit$model$backend, "craftgrn_native_warplda")
  expect_probability_matrix(fit$theta, nrow(dtm), 2L)
  expect_probability_matrix(fit$phi, 2L, ncol(dtm))
  expect_equal(rownames(fit$theta), rownames(dtm))
  expect_equal(colnames(fit$phi), colnames(dtm))
  expect_true(is.finite(fit$metrics$perplexity))
  expect_true(fit$metrics$n_tokens > 0)
})

test_that("native WarpLDA is deterministic with one thread and fixed seed", {
  dtm <- make_native_warplda_dtm()

  fit1 <- fit_warplda_one(dtm, K = 2L, iterations = 12L, alpha = 0.5, beta = 0.1, seed = 42L, n_check_convergence = 0L, n_iter_inference = 1L, n_threads = 1L, progressbar = FALSE)
  fit2 <- fit_warplda_one(dtm, K = 2L, iterations = 12L, alpha = 0.5, beta = 0.1, seed = 42L, n_check_convergence = 0L, n_iter_inference = 1L, n_threads = 1L, progressbar = FALSE)

  expect_equal(fit1$theta, fit2$theta, tolerance = 1e-12)
  expect_equal(fit1$phi, fit2$phi, tolerance = 1e-12)
})

test_that("native WarpLDA supports multicore fitting and recovers marker terms", {
  dtm <- make_native_warplda_dtm()
  n_threads <- min(2L, .warplda_default_threads())

  fit <- fit_warplda_one(
    dtm,
    K = 2L,
    iterations = 25L,
    alpha = 0.5,
    beta = 0.1,
    seed = 123L,
    n_check_convergence = 0L,
    n_iter_inference = 2L,
    n_threads = n_threads,
    progressbar = FALSE
  )

  expect_gte(fit$metrics$threads, 1L)
  expect_lte(fit$metrics$threads, n_threads)
  top_terms <- apply(fit$phi, 1L, function(x) names(sort(x, decreasing = TRUE))[1:4])
  marker_a <- paste0("term", 1:4)
  marker_b <- paste0("term", 5:8)
  overlap_a <- max(apply(top_terms, 2L, function(z) sum(z %in% marker_a)))
  overlap_b <- max(apply(top_terms, 2L, function(z) sum(z %in% marker_b)))
  expect_gte(overlap_a, 3L)
  expect_gte(overlap_b, 3L)
})


test_that("native WarpLDA is thread-count invariant for fixed seeds", {
  dtm <- make_native_warplda_dtm()
  n_threads <- min(4L, .warplda_default_threads())
  skip_if(n_threads < 2L, "OpenMP multicore test needs at least two threads.")

  fit1 <- fit_warplda_one(dtm, K = 2L, iterations = 12L, alpha = 0.5, beta = 0.1, seed = 101L, n_check_convergence = 0L, n_iter_inference = 2L, n_threads = 1L, progressbar = FALSE)
  fitn <- fit_warplda_one(dtm, K = 2L, iterations = 12L, alpha = 0.5, beta = 0.1, seed = 101L, n_check_convergence = 0L, n_iter_inference = 2L, n_threads = n_threads, progressbar = FALSE)

  expect_identical(fit1$theta, fitn$theta)
  expect_identical(fit1$phi, fitn$phi)
  expect_identical(fit1$metrics$perplexity, fitn$metrics$perplexity)
  expect_identical(fit1$metrics$loglik_approx, fitn$metrics$loglik_approx)
})

test_that("native WarpLDA validates positive integer-like counts", {
  skip_if_not_installed("Matrix")
  dtm <- Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1.5, 2),
    dims = c(2L, 2L),
    dimnames = list(c("doc1", "doc2"), c("term1", "term2"))
  )

  expect_error(
    fit_warplda_one(dtm, K = 2L, iterations = 2L, n_threads = 1L, progressbar = FALSE),
    "integer-like"
  )
})

test_that("run_warplda_models uses native fit files and cache reuse", {
  skip_if_not_installed("data.table")
  dtm <- make_native_warplda_dtm()
  tmp <- withr::local_tempdir()
  metrics_file <- file.path(tmp, "model_metrics.csv")

  out1 <- run_warplda_models(
    dtm,
    K_grid = c(2L, 3L),
    iterations = 8L,
    alpha_by_topic = FALSE,
    alpha = 0.5,
    beta = 0.1,
    seed = 9L,
    save_tmp_dir = file.path(tmp, "tmp_models"),
    workers = 1L,
    threads_per_model = 2L,
    metrics_file = metrics_file,
    verbose = FALSE
  )
  expect_equal(out1$metrics$K, c(2L, 3L))
  expect_true(all(file.exists(out1$fit_files)))
  fit_k2 <- readRDS(out1$fit_files[["2"]])
  expect_identical(fit_k2$model$backend, "craftgrn_native_warplda")
  expect_gte(fit_k2$metrics$threads, 1L)
  expect_lte(fit_k2$metrics$threads, 2L)

  out2 <- run_warplda_models(
    dtm,
    K_grid = c(2L, 3L),
    iterations = 8L,
    alpha_by_topic = FALSE,
    alpha = 0.5,
    beta = 0.1,
    seed = 9L,
    save_tmp_dir = file.path(tmp, "tmp_models"),
    workers = 1L,
    threads_per_model = 2L,
    metrics_file = metrics_file,
    verbose = FALSE
  )
  expect_equal(out2$metrics$K, c(2L, 3L))
  expect_equal(out1$fit_files, out2$fit_files)
})
