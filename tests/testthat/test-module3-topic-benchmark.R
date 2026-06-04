test_that("Module 3 topic benchmark scores existing models and writes review reports", {
  root <- tempfile("module3-topic-benchmark-")
  dir.create(root, recursive = TRUE)

  plan <- .module3_topic_method_plan(
    methods = "condition_aggr_weight_lda",
    k_grid = 2L
  )
  expect_equal(plan$method_setup, "cond fp aggr weight | LDA")

  model_dir <- plan$model_dir[[1L]]
  model_dir <- file.path(root, model_dir)
  vae_dir <- file.path(model_dir, "vae_models")
  dir.create(vae_dir, recursive = TRUE)

  theta <- data.table::data.table(
    doc_id = c("CondA_rep1::TF1", "CondA_rep2::TF2", "CondB_rep1::TF1", "CondB_rep2::TF2"),
    Topic1 = c(0.90, 0.82, 0.12, 0.18),
    Topic2 = c(0.10, 0.18, 0.88, 0.82)
  )
  phi <- data.table::data.table(
    term_id = c("GENE:G1", "GENE:G2", "PEAK:P1"),
    Topic1 = c(0.7, 0.2, 0.1),
    Topic2 = c(0.1, 0.7, 0.2)
  )
  data.table::fwrite(theta, file.path(vae_dir, "theta_K2.csv"))
  data.table::fwrite(phi, file.path(vae_dir, "phi_K2.csv"))

  extraction_dir <- file.path(
    root,
    plan$setup[[1L]],
    "02_topic_extraction",
    "condition_aggr_weight_lda_K2",
    plan$combo_id[[1L]]
  )
  dir.create(extraction_dir, recursive = TRUE)
  topic_links <- data.table::data.table(
    doc_id = theta$doc_id,
    topic_num = c(1L, 1L, 2L, 2L),
    tf = c("TF1", "TF2", "TF1", "TF2"),
    peak_id = c("P1", "P2", "P3", "P4"),
    gene_key = c("G1", "G2", "G3", "G4"),
    peak_score = c(1, 1, 1, 1),
    gene_score = c(1, 1, 1, 1),
    peak_pass = TRUE,
    gene_pass = TRUE,
    link_pass = TRUE
  )
  data.table::fwrite(topic_links, file.path(extraction_dir, "topic_links.csv"))

  comparisons <- data.table::data.table(
    condition_label = c("CondA_rep1", "CondA_rep2", "CondB_rep1", "CondB_rep2"),
    condition_group = c("CondA", "CondA", "CondB", "CondB")
  )

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = comparisons,
    methods = "condition_aggr_weight_lda",
    k_grid = 2L,
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = TRUE,
    reuse_if_exists = TRUE,
    verbose = FALSE
  )

  csv_dir <- file.path(root, "review_topic_experiments", "csv")
  html_dir <- file.path(root, "review_topic_experiments", "html")

  expect_s3_class(res$method_plan, "data.table")
  expect_true(file.exists(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv")))
  expect_true(file.exists(file.path(csv_dir, "theta_condition_separation_score_long.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_pass_state_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_shared_topic_counts.csv")))
  expect_true(file.exists(file.path(root, "review_topic_experiments", "tf_std_six_setups_pass_state_counts.pdf")))
  expect_true(file.exists(file.path(root, "review_topic_experiments", "tf_std_six_setups_shared_topic_counts.pdf")))
  expect_true(file.exists(file.path(html_dir, "theta_phi_and_group_mds.html")))
  expect_true(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report.html")))
  expect_true(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report_global_term_group.html")))

  score_mat <- data.table::fread(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv"))
  expect_equal(score_mat$method_setup, "cond fp aggr weight | LDA")
  expect_true("K2" %in% names(score_mat))
  expect_true(is.finite(score_mat$K2[[1L]]))
})
