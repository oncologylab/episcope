test_that("RUVr candidates are audited and preserve matrix dimensions", {
  skip_if_not_installed("RUVSeq")
  skip_if_not_installed("limma")
  skip_if_not_installed("matrixStats")
  set.seed(11)
  sample_id <- paste0("S", seq_len(12L))
  condition_id <- rep(c("A", "B"), each = 6L)
  study_id <- rep(rep(c("X", "Y"), each = 3L), 2L)
  matrix <- matrix(stats::rpois(1200L, lambda = 20), nrow = 100L)
  rownames(matrix) <- paste0("G", seq_len(nrow(matrix)))
  colnames(matrix) <- sample_id
  matrix[seq_len(20L), condition_id == "B"] <- matrix[seq_len(20L), condition_id == "B"] + 20
  matrix[, study_id == "Y"] <- matrix[, study_id == "Y"] + 3
  metadata <- data.frame(sample_id, condition_id, study_id)

  result <- .evaluate_ruvr_batch_correction(
    matrix,
    metadata,
    k_candidates = c(0L, 1L, 2L),
    top_n = 80L,
    verbose = FALSE
  )
  expect_equal(dim(result$corrected_log), dim(matrix))
  expect_equal(colnames(result$corrected_log), sample_id)
  expect_true(all(c("k", "accepted", "effect_spearman") %in% names(result$audit)))
  expect_true(result$selected_k %in% result$audit$k)
})
