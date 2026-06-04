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
  theta3 <- data.table::copy(theta)
  theta3[, Topic3 := c(0.02, 0.03, 0.04, 0.05)]
  theta3[, `:=`(
    Topic1 = Topic1 - Topic3 / 2,
    Topic2 = Topic2 - Topic3 / 2
  )]
  phi3 <- data.table::copy(phi)
  phi3[, Topic3 := c(0.2, 0.2, 0.6)]
  data.table::fwrite(theta3, file.path(vae_dir, "theta_K3.csv"))
  data.table::fwrite(phi3, file.path(vae_dir, "phi_K3.csv"))

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
    k_grid = c(2L, 3L),
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
  expect_true("K3" %in% names(score_mat))
  expect_true(is.finite(score_mat$K2[[1L]]))
  expect_true(is.finite(score_mat$K3[[1L]]))

  pass_counts <- data.table::fread(file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  expect_equal(unique(pass_counts$selected_k), 2L)
})

test_that("Module 3 topic benchmark uses a clean standard layout for one selected method", {
  root <- tempfile("module3-topic-standard-layout-")
  dir.create(root, recursive = TRUE)

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      condition_label = c("CondA", "CondB"),
      condition_group = c("CondA", "CondB")
    ),
    methods = "condition_aggr_weight_lda",
    k_grid = 2L,
    output_layout = "standard",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = FALSE,
    verbose = FALSE
  )

  expect_equal(res$output_layout, "standard")
  expect_equal(res$method_plan$run_id, "selected")
  expect_equal(res$method_plan$topic_documents_dir, file.path(root, "topic_documents"))
  expect_equal(res$method_plan$topic_models_dir, file.path(root, "topic_models", "lda"))
  expect_equal(res$method_plan$topic_extraction_dir, file.path(root, "topic_extraction"))
  expect_equal(res$review_dir, file.path(root, "review"))
  expect_false(grepl("std_tf_", res$method_plan$topic_models_dir, fixed = TRUE))
})

test_that("Module 3 topic benchmark uses shallow run folders for method grids", {
  root <- tempfile("module3-topic-shallow-benchmark-")
  dir.create(root, recursive = TRUE)

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      condition_label = c("CondA", "CondB"),
      condition_group = c("CondA", "CondB")
    ),
    methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    k_grid = 2L,
    output_layout = "benchmark",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = FALSE,
    verbose = FALSE
  )

  expect_equal(res$output_layout, "benchmark")
  expect_equal(res$review_dir, file.path(root, "review"))
  expect_equal(res$method_plan$run_id, c("run_001", "run_002"))
  expect_true(all(grepl("^run_[0-9]{3}_", basename(res$method_plan$run_dir))))
  expect_equal(res$method_plan$topic_documents_dir, file.path(res$method_plan$run_dir, "topic_documents"))
  expect_equal(res$method_plan$topic_models_dir, file.path(res$method_plan$run_dir, "topic_models"))
  expect_equal(res$method_plan$topic_extraction_dir, file.path(res$method_plan$run_dir, "topic_extraction"))
  expect_true(file.exists(file.path(root, "runs.csv")))
})

test_that("Module 3 topic benchmark writes typed empty topic-link summaries", {
  root <- tempfile("module3-topic-benchmark-empty-links-")
  dir.create(root, recursive = TRUE)

  plan <- .module3_topic_method_plan(
    methods = "condition_aggr_weight_lda",
    k_grid = 2L
  )
  vae_dir <- file.path(root, plan$model_dir[[1L]], "vae_models")
  dir.create(vae_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(doc_id = c("CondA::TF1", "CondB::TF1"), Topic1 = c(0.9, 0.1), Topic2 = c(0.1, 0.9)),
    file.path(vae_dir, "theta_K2.csv")
  )
  data.table::fwrite(
    data.table::data.table(term_id = c("GENE:G1", "PEAK:P1"), Topic1 = c(0.8, 0.2), Topic2 = c(0.2, 0.8)),
    file.path(vae_dir, "phi_K2.csv")
  )

  comparisons <- data.table::data.table(
    condition_label = c("CondA", "CondB"),
    condition_group = c("CondA", "CondB")
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
    verbose = FALSE
  )

  expect_named(
    res$topic_link_summary$pass,
    c("method_order", "method_setup", "setup", "model_label", "selected_k", "status", "count")
  )
  expect_equal(nrow(res$topic_link_summary$pass), 0L)
  expect_equal(nrow(res$topic_link_summary$shared), 0L)
})

test_that("Module 3 topic benchmark accepts explicit context-aware design rows", {
  design <- data.table::data.table(
    context_type = c("condition", "comparison"),
    comparison_label = c("CondA_rep1", "CondA_vs_CondB::Target-Up"),
    display_label = c("CondA rep1", "CondA vs CondB Up"),
    metric_group = c("CondA", "CondA_vs_CondB")
  )
  out <- .m3tb_design_table(design)
  expect_equal(out$context_type, c("condition", "comparison"))
  expect_equal(out$metric_group, c("CondA", "CondA_vs_CondB"))
})

test_that("Module 3 topic benchmark scores replicate-resolved condition and comparison documents", {
  root <- tempfile("module3-topic-benchmark-replicates-")
  dir.create(root, recursive = TRUE)

  plan <- .module3_topic_method_plan(
    methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    k_grid = 2L
  )
  for (i in seq_len(nrow(plan))) {
    vae_dir <- file.path(root, plan$model_dir[[i]], "vae_models")
    dir.create(vae_dir, recursive = TRUE)
    if (identical(plan$context_type[[i]], "condition")) {
      theta <- data.table::data.table(
        doc_id = c("CondA_1::TF1", "CondA_2::TF1", "CondB_1::TF1", "CondB_2::TF1"),
        Topic1 = c(0.90, 0.86, 0.14, 0.10),
        Topic2 = c(0.10, 0.14, 0.86, 0.90)
      )
    } else {
      theta <- data.table::data.table(
        doc_id = c(
          "CondA_vs_CondB__CondA_1::TF1::Target-Up",
          "CondA_vs_CondB__CondA_2::TF1::Target-Up",
          "CondA_vs_CondB__CondB_1::TF1::Target-Down",
          "CondA_vs_CondB__CondB_2::TF1::Target-Down"
        ),
        Topic1 = c(0.88, 0.82, 0.16, 0.12),
        Topic2 = c(0.12, 0.18, 0.84, 0.88)
      )
    }
    phi <- data.table::data.table(
      term_id = c("GENE:G1", "GENE:G2"),
      Topic1 = c(0.8, 0.2),
      Topic2 = c(0.2, 0.8)
    )
    data.table::fwrite(theta, file.path(vae_dir, "theta_K2.csv"))
    data.table::fwrite(phi, file.path(vae_dir, "phi_K2.csv"))
  }

  multiomic_data <- list(
    matrices = list(
      fp_score = matrix(
        1,
        nrow = 1,
        ncol = 4,
        dimnames = list("fp1", c("CondA_1", "CondA_2", "CondB_1", "CondB_2"))
      )
    )
  )
  comparisons <- data.table::data.table(
    cond1_label = "CondA",
    cond2_label = "CondB"
  )

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    multiomic_data = multiomic_data,
    comparisons = comparisons,
    output_dir = root,
    methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    k_grid = 2L,
    replicate_documents = TRUE,
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = TRUE,
    verbose = FALSE
  )

  csv_dir <- file.path(root, "review_topic_experiments", "csv")
  expect_true(file.exists(file.path(csv_dir, "theta_condition_replicate_separation_score_heatmap_values_matrix.csv")))
  expect_false(file.exists(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv")))
  expect_equal(res$score$score_prefix, "theta_condition_replicate_separation")

  per_label <- data.table::fread(file.path(csv_dir, "theta_condition_replicate_separation_per_label.csv"))
  expect_equal(unique(per_label[context_type == "condition", group_size]), 2L)
  expect_equal(unique(per_label[context_type == "comparison", group_size]), 2L)
  expect_setequal(
    unique(per_label[context_type == "condition", metric_group]),
    c("CondA", "CondB")
  )
  expect_setequal(
    unique(per_label[context_type == "comparison", metric_group]),
    c("CondA_vs_CondB::Target-Up", "CondA_vs_CondB::Target-Down")
  )
})
