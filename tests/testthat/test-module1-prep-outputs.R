test_that("Module 1 preparation writes quantile-normalized FP score CSV", {
  out_dir <- tempfile("craftgrn-module1-prep-output-")
  grn_set <- list(
    fp_score_condition_qn = tibble::tibble(
      peak_ID = c("fp1", "fp2"),
      sample_a = c(1.5, 2.5),
      sample_b = c(3.5, 4.5)
    )
  )

  out_path <- .write_fp_score_qn_csv(
    grn_set = grn_set,
    out_dir = out_dir,
    db = "JASPAR2024",
    verbose = FALSE
  )

  expect_equal(basename(out_path), "01_fp_scores_qn_JASPAR2024.csv")
  expect_true(file.exists(out_path))
  written <- data.table::fread(out_path, showProgress = FALSE)
  expect_equal(as.data.frame(written), as.data.frame(grn_set$fp_score_condition_qn))
})
