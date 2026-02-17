topic_distance_script <- normalizePath(file.path(
  Sys.getenv(
    "CRAFTGRN_REPO_DIR",
    unset = normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
  ),
  "dev/benchmark/03_benchmark_topic_theta_phi_distance_correlation.R"
), mustWork = FALSE)

test_that("theta-phi distance panels use shared setup and method order", {
  testthat::skip_if_not(file.exists(topic_distance_script), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(topic_distance_script))))
  on.exit(setwd(old_wd), add = TRUE)
  source(topic_distance_script, local = TRUE)

  summary <- data.table::data.table(
    setup = rep(c(
      "std_tf_cond_fp_uniq",
      "std_tf_diff_fp_uniq",
      "std_tf_cond_fp_aggr",
      "std_tf_diff_fp_aggr",
      "std_tf_cond_fp_aggr_weight",
      "std_tf_diff_fp_aggr_weight"
    ), each = 2L),
    setup_label = rep(c(
      "cond fp uniq",
      "diff fp uniq",
      "cond fp aggr",
      "diff fp aggr",
      "cond fp aggr weight",
      "diff fp aggr weight"
    ), each = 2L),
    model_label = rep(c("LDA", "MultiVI"), 6L)
  )

  out <- .theta_phi_panel_order(summary)

  expect_equal(
    out,
    c(
      "cond fp uniq\nMultiVI",
      "cond fp uniq\nLDA",
      "diff fp uniq\nMultiVI",
      "diff fp uniq\nLDA",
      "cond fp aggr\nMultiVI",
      "cond fp aggr\nLDA",
      "diff fp aggr\nMultiVI",
      "diff fp aggr\nLDA",
      "cond fp aggr weight\nMultiVI",
      "cond fp aggr weight\nLDA",
      "diff fp aggr weight\nMultiVI",
      "diff fp aggr weight\nLDA"
    )
  )
})
