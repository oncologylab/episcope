test_that("Module 3 report API exposes publication-facing helpers", {
  exports <- getNamespaceExports("craftgrn")
  expect_true("build_module3_qc_report" %in% exports)
  expect_true("visualize_topic_modeling_results" %in% exports)
  expect_true("visualize_differential_grns" %in% exports)
  expect_false("report_differential_pathways" %in% exports)
  expect_false("report_master_tfs" %in% exports)
  expect_false("default_diff_grn_pathway_databases" %in% exports)
  expect_false("run_diff_grn_pathway_analysis" %in% exports)
  expect_false("run_diff_grn_pathway_enrichment" %in% exports)
  expect_false("read_diff_grn_pathway_enrichment" %in% exports)
  expect_false("select_diff_grn_pathway_terms" %in% exports)
  expect_false("plot_diff_grn_pathway_dotplot" %in% exports)
  expect_false("plot_diff_grn_pathway_network" %in% exports)
  expect_false("run_diff_grn_master_tf_summary" %in% exports)
  expect_false("summarize_diff_grn_master_tf_links" %in% exports)
  expect_false("plot_diff_grn_master_tf_summary" %in% exports)
})

test_that("Module 3 QC report writes differential TF summaries", {
  root <- tempfile("module3-qc-topic-")
  diff_dir <- tempfile("module3-qc-diff-")
  dir.create(file.path(root, "review", "csv"), recursive = TRUE)
  dir.create(diff_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(method = "condition_aggr_weight_lda", method_setup = "cond fp aggr weight | LDA"),
    file.path(root, "review", "csv", "module3_topic_method_plan.csv")
  )
  links <- data.table::data.table(
    comparison_id = "CondA_vs_CondB",
    comparison_group = "Group1",
    cond1_id = "CondA",
    cond2_id = "CondB",
    direction_group = c("up", "up", "down"),
    tf = c("TF1", "TF1", "TF2"),
    gene_key = c("GeneA", "GeneB", "GeneC"),
    peak_id = paste0("peak", 1:3),
    delta_fp_score = c(1.2, 0.8, -1.1),
    log2FC_tf_expr = c(1.1, 1.1, -0.7)
  )
  data.table::fwrite(links[direction_group == "up"], file.path(diff_dir, "CondA_vs_CondB_filtered_links_up.csv"))
  data.table::fwrite(links[direction_group == "down"], file.path(diff_dir, "CondA_vs_CondB_filtered_links_down.csv"))

  report <- build_module3_qc_report(root, differential_links_dir = diff_dir, verbose = FALSE)

  expect_true(file.exists(report))
  summary_dir <- file.path(dirname(report), "differential_tf_summary")
  expect_true(file.exists(file.path(summary_dir, "module3_top_differential_tfs.csv")))
  expect_true(file.exists(file.path(summary_dir, "per_comparison", "CondA_vs_CondB_differential_tfs.csv")))
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_true(grepl("Differential GRN Summary", html, fixed = TRUE))
  expect_true(grepl("Top Differential TFs", html, fixed = TRUE))
})

test_that("Module 3 visualization helpers write self-contained index pages", {
  topic_dir <- tempfile("module3-topic-vis-")
  diff_dir <- tempfile("module3-diff-vis-")
  dir.create(file.path(topic_dir, "review", "html"), recursive = TRUE)
  dir.create(diff_dir, recursive = TRUE)
  writeLines("<!doctype html><html><body>Topic MDS</body></html>", file.path(topic_dir, "review", "html", "topic_method_k_topic_mds_report.html"))
  data.table::fwrite(
    data.table::data.table(
      comparison_id = "CondA_vs_CondB",
      comparison_group = "Group1",
      cond1_id = "CondA",
      cond2_id = "CondB",
      direction_group = "up",
      tf = "TF1",
      gene_key = "GeneA",
      peak_id = "peak1",
      delta_fp_score = 1,
      log2FC_tf_expr = 1
    ),
    file.path(diff_dir, "CondA_vs_CondB_filtered_links_up.csv")
  )

  topic_html <- visualize_topic_modeling_results(topic_dir, verbose = FALSE)
  diff_html <- visualize_differential_grns(diff_dir, verbose = FALSE)

  expect_true(file.exists(topic_html))
  expect_true(file.exists(file.path(dirname(topic_html), "topic_modeling_results_manifest.csv")))
  expect_true(file.exists(diff_html))
  expect_true(file.exists(file.path(dirname(diff_html), "differential_grn_summary.csv")))
})
