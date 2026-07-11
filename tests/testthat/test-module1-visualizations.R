test_that("TFBS condition comparison uses active sites and summed FP scores", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  result <- predict_tfbs(compact, r_cutoff = 0.8, verbose = FALSE)

  plot <- plot_tfbs_condition_comparison(
    result,
    cond1 = "cond_a",
    cond2 = "cond_d",
    pseudocount = 1,
    verbose = FALSE
  )
  summary <- attr(plot, "tfbs_summary")
  tf_a <- summary[summary$tf == "TF_A", , drop = FALSE]

  expect_s3_class(plot, "ggplot")
  expect_equal(tf_a$n_tfbs_cond1, 1L)
  expect_equal(tf_a$n_tfbs_cond2, 1L)
  expect_equal(tf_a$delta_n_predicted_tfbs, 0L)
  expect_equal(tf_a$log2fc_fp_score, log2((1 + 1) / (9 + 1)))
  expect_equal(tf_a$max_tf_expression, 9)
  expect_equal(tf_a$log2fc_tf_expression, log2((1 + 1) / (9 + 1)))
  expect_equal(plot$labels$title, "cond_a v cond_d")
})

test_that("TF-TF co-binding heatmap supports counts and JSD", {
  predicted <- tibble::tibble(
    tf = c("TF1", "TF1", "TF2", "TF2", "TF3"),
    fp_id = c("fp1", "fp2", "fp2", "fp3", "fp4")
  )

  count_plot <- plot_tf_tf_cobinding_heatmap(predicted, metric = "absolute", verbose = FALSE)
  count_matrix <- attr(count_plot, "cobinding_matrix")
  jsd_plot <- plot_tf_tf_cobinding_heatmap(predicted, metric = "jsd", verbose = FALSE)
  jsd <- attr(jsd_plot, "cobinding_matrix")

  expect_s3_class(count_plot, "ggplot")
  expect_equal(count_matrix["TF1", "TF2"], 1)
  expect_equal(diag(count_matrix), c(TF1 = 2, TF2 = 2, TF3 = 1))
  expect_equal(jsd, t(jsd))
  expect_equal(unname(diag(jsd)), rep(0, 3))
  expect_true(all(jsd >= 0 & jsd <= 1))
})

test_that("TFBS UMAP report provides interactive cluster choices", {
  skip_if_not_installed("uwot")
  predicted <- tibble::tibble(
    tf = rep(paste0("TF", 1:6), each = 3),
    fp_id = c(
      "p1", "p2", "p3", "p1", "p2", "p4", "p2", "p3", "p5",
      "p3", "p4", "p6", "p4", "p5", "p7", "p5", "p6", "p8"
    )
  )
  out <- tempfile(fileext = ".html")

  path <- build_tfbs_umap_report(
    predicted,
    output_file = out,
    top_variable_tfbs = 8L,
    cluster_range = 2:4,
    default_clusters = 3L,
    seed = 11L,
    verbose = FALSE
  )
  page <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(file.exists(path))
  expect_match(page, "Number of clusters", fixed = TRUE)
  expect_match(page, "option value=\"3\" selected", fixed = TRUE)
  expect_match(page, "TF UMAP from highly variable predicted TFBS", fixed = TRUE)
  expect_match(page, "predicted TFBS", fixed = TRUE)
})

test_that("Module 1 config provenance keeps relevant values", {
  cfg <- craftgrn:::.module1_relevant_config(list(
    db = "JASPAR2024",
    ref_genome = "hg38",
    threshold_fp_tf_corr_r = 0.4,
    module3_backend = "lda",
    api_key = "do-not-keep"
  ))

  expect_equal(cfg$db, "JASPAR2024")
  expect_equal(cfg$threshold_fp_tf_corr_r, 0.4)
  expect_false("module3_backend" %in% names(cfg))
  expect_false("api_key" %in% names(cfg))
})
