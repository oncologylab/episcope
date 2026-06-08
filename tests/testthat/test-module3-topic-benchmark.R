test_that("Module 3 topic benchmark scores existing models and writes review reports", {
  root <- tempfile("module3-topic-benchmark-")
  dir.create(root, recursive = TRUE)

  plan <- .module3_topic_method_plan(
    methods = "condition_aggr_weight_lda",
    k_grid = 2L
  )
  expect_equal(plan$method_setup, "cond fp aggr weight | LDA")

  model_dir <- file.path(root, "topic_models", "lda")
  vae_dir <- file.path(model_dir, "vae_models")
  dir.create(vae_dir, recursive = TRUE)

  theta <- data.table::data.table(
    doc_id = c("CondA_rep1::TF1", "CondA_rep2::TF2", "CondB_rep1::TF1", "CondB_rep2::TF2"),
    Topic1 = c(0.90, 0.82, 0.12, 0.18),
    Topic2 = c(0.10, 0.18, 0.88, 0.82)
  )
  phi <- data.table::data.table(
    term_id = c("Topic1", "Topic2"),
    `GENE:G1` = c(0.7, 0.1),
    `GENE:G2` = c(0.2, 0.7),
    `PEAK:P1` = c(0.1, 0.2)
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
  phi3 <- data.table::rbindlist(list(
    phi3,
    data.table::data.table(term_id = "Topic3", `GENE:G1` = 0.2, `GENE:G2` = 0.2, `PEAK:P1` = 0.6)
  ), use.names = TRUE)
  data.table::fwrite(theta3, file.path(vae_dir, "theta_K3.csv"))
  data.table::fwrite(phi3, file.path(vae_dir, "phi_K3.csv"))

  extraction_dir <- file.path(
    root,
    "topic_extraction",
    "K2"
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
  pathway_rows <- data.table::data.table(
    topic = c(1L, 1L, 2L),
    pathway = c("Hallmark: Signal A", "Reactome: Signal B", "Hallmark: Signal A"),
    padj = c(0.001, 0.002, 0.003),
    overlap = c("2/3", "1/3", "1/3"),
    overlap_hits = c(2L, 1L, 1L),
    genes = c("G1;G2", "G2", "G1")
  )
  data.table::fwrite(
    pathway_rows,
    file.path(extraction_dir, "topic_pathway_enrichment_peak_and_gene_dotplot.csv")
  )
  dir.create(file.path(model_dir, "rds"), recursive = TRUE)
  saveRDS(
    data.table::data.table(term_id = c("GENE:G1", "GENE:G2", "GENE:G3", "PEAK:P1")),
    file.path(model_dir, "rds", "doc_term.rds")
  )
  data.table::fwrite(
    data.table::data.table(
      topic = 0L,
      pathway = c("Hallmark: Signal A", "Reactome: Signal B"),
      pathway_key = c("Hallmark: Signal A", "Reactome: Signal B"),
      padj = c(0.01, 0.02),
      overlap = c("10/100", "5/100"),
      overlap_hits = c(10L, 5L),
      genes = c("G1;G2;G3", "G2;G3")
    ),
    file.path(model_dir, "topic_pathway_enrichment_gene_universe_all.csv")
  )
  per_group_dir <- file.path(extraction_dir, "per_comparison_pathway_peak_pass_gene_pass")
  dir.create(per_group_dir, recursive = TRUE)
  data.table::fwrite(
    pathway_rows[topic == 1L],
    file.path(per_group_dir, "CondA_rep1_All_dotplot.csv")
  )

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
    output_layout = "standard",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = TRUE,
    reuse_if_exists = TRUE,
    verbose = FALSE
  )

  csv_dir <- file.path(root, "review", "tables")
  html_dir <- file.path(root, "review")

  expect_s3_class(res$method_plan, "data.table")
  expect_true(file.exists(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv")))
  expect_true(file.exists(file.path(csv_dir, "theta_condition_separation_score_long.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_pass_state_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_shared_topic_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "theta_group_mds_points.csv")))
  expect_true(file.exists(file.path(root, "review", "tf_std_six_setups_pass_state_counts.pdf")))
  expect_true(file.exists(file.path(root, "review", "tf_std_six_setups_shared_topic_counts.pdf")))
  expect_true(file.exists(file.path(html_dir, "theta_phi_and_group_mds.html")))
  expect_true(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report.html")))
  expect_true(file.exists(file.path(html_dir, "topic_method_k_condition_mds_report.html")))
  expect_false(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report_global_term_group.html")))
  expect_false(file.exists(file.path(html_dir, "topic_method_k_condition_mds_report_global_term_group.html")))

  topic_html <- paste(readLines(file.path(html_dir, "topic_method_k_topic_mds_report.html"), warn = FALSE), collapse = "\n")
  condition_html <- paste(readLines(file.path(html_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  theta_html <- paste(readLines(file.path(html_dir, "theta_phi_and_group_mds.html"), warn = FALSE), collapse = "\n")
  expect_match(topic_html, "Intertopic Distance Map", fixed = TRUE)
  expect_match(topic_html, "Condition Waterfall", fixed = TRUE)
  expect_match(topic_html, "Pathways", fixed = TRUE)
  expect_match(topic_html, "Export SVG", fixed = TRUE)
  expect_match(topic_html, "Genes in full document gene-universe enrichment", fixed = TRUE)
  expect_match(topic_html, "topic genes", fixed = TRUE)
  expect_match(topic_html, "universe remainder", fixed = TRUE)
  expect_match(topic_html, "gene_total_universe", fixed = TRUE)
  expect_match(topic_html, "mdsLeader", fixed = TRUE)
  expect_match(topic_html, "paletteSelect", fixed = TRUE)
  expect_match(topic_html, "mdsLayer", fixed = TRUE)
  expect_match(topic_html, "waterfallLayer", fixed = TRUE)
  expect_match(topic_html, "pathLayer", fixed = TRUE)
  expect_match(topic_html, "mdsTooltip", fixed = TRUE)
  expect_match(topic_html, "wfTooltip", fixed = TRUE)
  expect_match(topic_html, "pathTooltip", fixed = TRUE)
  expect_match(topic_html, "mdsStats", fixed = TRUE)
  expect_match(topic_html, "waterfallStats", fixed = TRUE)
  expect_match(topic_html, "pathStats", fixed = TRUE)
  expect_match(topic_html, "Mean document-to-topic probability", fixed = TRUE)
  expect_match(topic_html, "marker", fixed = TRUE)
  expect_match(condition_html, "Condition/Comparison MDS", fixed = TRUE)
  expect_match(condition_html, "Topic Waterfall", fixed = TRUE)
  expect_match(condition_html, "mdsImage", fixed = TRUE)
  expect_match(condition_html, "mdsHotspotLayer", fixed = TRUE)
  expect_match(condition_html, "function drawMdsHotspots()", fixed = TRUE)
  expect_match(condition_html, "mdsHotspot", fixed = TRUE)
  expect_match(condition_html, "groupColor", fixed = TRUE)
  expect_match(condition_html, "selectGroup", fixed = TRUE)
  expect_match(condition_html, "pathLabelTopicSpecific", fixed = TRUE)
  expect_match(condition_html, "pathLabelGroupSpecific", fixed = TRUE)
  expect_match(condition_html, "pathLabelBothSpecific", fixed = TRUE)
  expect_match(condition_html, "Pathway name colors", fixed = TRUE)
  expect_match(condition_html, "GROUP_MDS", fixed = TRUE)
  expect_match(condition_html, "selectedGroupRows", fixed = TRUE)
  expect_match(condition_html, "Document-to-topic probability", fixed = TRUE)
  expect_match(condition_html, "waterfallLayer", fixed = TRUE)
  expect_match(condition_html, "waterfallStats", fixed = TRUE)
  expect_match(theta_html, "theta_group_mds_k2.png", fixed = TRUE)

  condition_svg <- file.path(root, "review", "assets", "condition_mds_K2.svg")
  expect_true(file.exists(condition_svg))

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

test_that("Module 3 topic benchmark exposes project-agnostic sample subsetting", {
  expect_false("celllines" %in% names(formals(run_module3_topic_benchmark)))
  expect_true("sample_subset" %in% names(formals(run_module3_topic_benchmark)))
  expect_true("analysis_label" %in% names(formals(run_module3_topic_benchmark)))
  expect_false("celllines" %in% names(formals(train_topic_models)))
  expect_true("sample_subset" %in% names(formals(train_topic_models)))
  expect_true("analysis_label" %in% names(formals(train_topic_models)))
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
  vae_dir <- file.path(root, "topic_models", "lda", "fixture_model", "vae_models")
  dir.create(vae_dir, recursive = TRUE, showWarnings = FALSE)
  expect_true(dir.exists(vae_dir))
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
    output_layout = "standard",
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

  plan <- .m3tb_apply_output_layout(
    .module3_topic_method_plan(
      methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
      k_grid = 2L
    ),
    root,
    "benchmark"
  )
  for (i in seq_len(nrow(plan))) {
    vae_dir <- file.path(plan$topic_models_dir[[i]], "fixture_model", "vae_models")
    dir.create(vae_dir, recursive = TRUE, showWarnings = FALSE)
    expect_true(dir.exists(vae_dir))
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
    output_layout = "benchmark",
    replicate_documents = TRUE,
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = TRUE,
    verbose = FALSE
  )

  csv_dir <- file.path(root, "review", "tables")
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

test_that("Module 3 topic links default to compact pass-only output", {
  edges <- data.table::data.table(
    doc_id = c("D1", "D1"),
    tf = c("TF1", "TF1"),
    peak_id = c("P1", "P2"),
    gene_key = c("G1", "G2")
  )
  score_mat <- matrix(
    c(0.9, 0.1, 0.8, 0.2, 0.2, 0.7, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE
  )
  rownames(score_mat) <- c("Topic1", "Topic2")
  colnames(score_mat) <- c("PEAK:P1", "PEAK:P2", "GENE:G1", "GENE:G2")
  out_dir <- tempfile("module3-topic-links-")
  dir.create(out_dir, recursive = TRUE)
  res <- compute_topic_links(
    edges,
    score_mat,
    link_method = "gene_prob",
    link_prob_cutoff = "max",
    out_file = file.path(out_dir, "topic_links.csv"),
    overwrite = TRUE
  )
  expect_s3_class(res, "data.table")
  expect_false(file.exists(file.path(out_dir, "topic_links.csv")))
  expect_true(file.exists(file.path(out_dir, "topic_links_pass.csv")))
  expect_true(file.exists(file.path(out_dir, "topic_link_summary.csv")))
  pass <- data.table::fread(file.path(out_dir, "topic_links_pass.csv"))
  expect_true(nrow(pass) <= nrow(res))
  expect_true(all(pass$gene_prob_pass))
})

test_that("Module 3 prepares reusable topic input caches", {
  root <- tempfile("module3-topic-inputs-")
  filtered_dir <- file.path(root, "filtered")
  out_dir <- file.path(root, "topic_documents")
  dir.create(filtered_dir, recursive = TRUE)
  links <- data.table::data.table(
    comparison_id = "Cond1_vs_Cond2",
    cond1_id = "Cond1",
    cond2_id = "Cond2",
    tf = c("TF1", "TF1"),
    gene_key = c("G1", "G2"),
    peak_id = c("P1", "P2"),
    log2FC_fp_score = c(2, -2),
    delta_fp_score = c(2, -2),
    log2FC_gene_expr = c(2, -2),
    log2FC_tf_expr = c(1, 1),
    fp_bound_cond1 = TRUE,
    fp_bound_cond2 = TRUE,
    tf_expr_cond1 = 10,
    tf_expr_cond2 = 10,
    gene_expr_cond1 = 10,
    gene_expr_cond2 = 10,
    fp_score_cond1 = 10,
    fp_score_cond2 = 10
  )
  active_path <- file.path(filtered_dir, "Cond1_vs_Cond2_delta_links_filtered_up.csv")
  stale_path <- file.path(filtered_dir, "Stale_Contrast_delta_links_filtered_up.csv")
  data.table::fwrite(links, active_path)
  stale_links <- data.table::copy(links)
  stale_links[, comparison_id := "Stale_Contrast"]
  data.table::fwrite(stale_links, stale_path)
  data.table::fwrite(
    data.table::data.table(
      comparison_id = "Cond1_vs_Cond2",
      up_path = active_path,
      down_path = NA_character_
    ),
    file.path(filtered_dir, "filtered_links_manifest.csv")
  )
  res <- module3_construct_docs(
    filtered_dir = filtered_dir,
    output_dir = out_dir,
    tf_cluster_map = c(TF1 = "K01"),
    doc_mode = "tf",
    doc_design = "comparison",
    fp_term_mode = "aggregate_weight",
    min_df = 1,
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    threshold_tf_expr = 0,
    direction_consistency = "none",
    require_tf_expr_either = FALSE,
    require_gene_expr_either = FALSE,
    require_fp_bound_either = FALSE,
    overwrite = TRUE,
    verbose = FALSE
  )
  expect_true(file.exists(file.path(out_dir, "rds", "edges_docs.rds")))
  expect_true(file.exists(file.path(out_dir, "rds", "doc_term.rds")))
  expect_true(file.exists(file.path(out_dir, "rds", "dtm.rds")))
  expect_true(file.exists(file.path(out_dir, "topic_input_summary.csv")))
  expect_equal(res$summary$n_link_rows_loaded[[1L]], nrow(links))
  expect_gt(res$summary$n_documents[[1L]], 0)
  reused <- module3_construct_docs(
    filtered_dir = filtered_dir,
    output_dir = out_dir,
    tf_cluster_map = c(TF1 = "K01"),
    doc_mode = "tf",
    doc_design = "comparison",
    fp_term_mode = "aggregate_weight",
    overwrite = FALSE,
    verbose = FALSE
  )
  expect_true(reused$reused)
})

test_that("Module 3 production wrapper exposes compact defaults and QC report", {
  expect_true("run_topic_modeling" %in% getNamespaceExports("craftgrn"))
  expect_true("module3_construct_docs" %in% getNamespaceExports("craftgrn"))
  expect_true("build_module3_qc_report" %in% getNamespaceExports("craftgrn"))
  expect_true("warplda_iterations" %in% names(formals(run_topic_modeling)))
  root <- tempfile("module3-qc-")
  dir.create(file.path(root, "review", "tables"), recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(method = "condition_aggr_weight_lda", method_setup = "cond fp aggr weight | LDA"),
    file.path(root, "review", "tables", "module3_topic_method_plan.csv")
  )
  report <- build_module3_qc_report(root, verbose = FALSE)
  expect_true(file.exists(report))
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_true(grepl("Module 3 QC report", html, fixed = TRUE))
  expect_true(grepl("Method Plan", html, fixed = TRUE))
})
