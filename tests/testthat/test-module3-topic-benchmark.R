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
  data.table::fwrite(
    data.table::data.table(
      link_method = "gene_prob",
      output_mode = "pass",
      n_scored_rows = 8L,
      n_pass_rows = 4L,
      full_file = NA_character_,
      pass_file = "topic_links_pass.csv"
    ),
    file.path(extraction_dir, "topic_link_summary.csv")
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

  pass_counts <- data.table::fread(file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  shared_counts <- data.table::fread(file.path(csv_dir, "topic_setup_shared_topic_counts.csv"))
  expect_setequal(pass_counts$status, c("Pass", "Fail"))
  expect_equal(sum(pass_counts[status == "Pass", count]), 4)
  expect_equal(sum(pass_counts[status == "Fail", count]), 0)
  expect_equal(unique(pass_counts$count_basis), "Standard TF-FP-gene links")
  expect_true("Pathways" %in% shared_counts$unit)

  topic_html <- paste(readLines(file.path(html_dir, "topic_method_k_topic_mds_report.html"), warn = FALSE), collapse = "\n")
  condition_html <- paste(readLines(file.path(html_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  theta_html <- paste(readLines(file.path(html_dir, "theta_phi_and_group_mds.html"), warn = FALSE), collapse = "\n")
  expect_match(topic_html, "Topic Method/K Reports", fixed = TRUE)
  expect_match(topic_html, "Method <select", fixed = TRUE)
  expect_match(topic_html, "K <select", fixed = TRUE)
  expect_match(topic_html, "topic_reports/", fixed = TRUE)
  expect_match(topic_html, "embedSrc(hit.src)", fixed = TRUE)
  expect_match(topic_html, "embed=1", fixed = TRUE)
  expect_no_match(topic_html, "srcdoc", fixed = TRUE)
  expect_match(condition_html, "Topic Method/K Condition Reports", fixed = TRUE)
  expect_match(condition_html, "Method <select", fixed = TRUE)
  expect_match(condition_html, "K <select", fixed = TRUE)
  expect_match(condition_html, "condition_topic_reports/", fixed = TRUE)
  expect_match(condition_html, "embedSrc(hit.src)", fixed = TRUE)
  expect_match(condition_html, "embed=1", fixed = TRUE)
  expect_no_match(condition_html, "srcdoc", fixed = TRUE)

  topic_report_files <- list.files(
    file.path(html_dir, "topic_reports"),
    pattern = "_topic_report[.]html$",
    full.names = TRUE
  )
  condition_report_files <- list.files(
    file.path(html_dir, "condition_topic_reports"),
    pattern = "_condition_topic_report[.]html$",
    full.names = TRUE
  )
  topic_page_files <- list.files(
    file.path(html_dir, "topic_reports", "pages"),
    pattern = "_topic_report[.]html$",
    full.names = TRUE
  )
  condition_page_files <- list.files(
    file.path(html_dir, "condition_topic_reports", "pages"),
    pattern = "_condition_topic_report[.]html$",
    full.names = TRUE
  )
  expect_true(length(topic_report_files) >= 1L)
  expect_true(length(condition_report_files) >= 1L)
  expect_true(length(topic_page_files) >= 1L)
  expect_true(length(condition_page_files) >= 1L)
  topic_wrapper_html <- paste(readLines(topic_report_files[[1L]], warn = FALSE), collapse = "\n")
  condition_wrapper_html <- paste(readLines(condition_report_files[[1L]], warn = FALSE), collapse = "\n")
  expect_match(topic_wrapper_html, "K <select", fixed = TRUE)
  expect_match(topic_wrapper_html, "pages/", fixed = TRUE)
  expect_match(condition_wrapper_html, "K <select", fixed = TRUE)
  expect_match(condition_wrapper_html, "pages/", fixed = TRUE)
  topic_detail_html <- paste(readLines(topic_page_files[[1L]], warn = FALSE), collapse = "\n")
  condition_detail_html <- paste(readLines(condition_page_files[[1L]], warn = FALSE), collapse = "\n")
  expect_match(topic_detail_html, "Intertopic Distance Map", fixed = TRUE)
  expect_match(topic_detail_html, "Condition Waterfall", fixed = TRUE)
  expect_match(topic_detail_html, "Pathways", fixed = TRUE)
  expect_match(topic_detail_html, "Export SVG", fixed = TRUE)
  expect_match(topic_detail_html, "image/svg+xml", fixed = TRUE)
  expect_match(topic_detail_html, "XMLSerializer", fixed = TRUE)
  expect_match(topic_detail_html, "exportSvg", fixed = TRUE)
  expect_no_match(topic_detail_html, "document.documentElement.outerHTML", fixed = TRUE)
  expect_no_match(topic_detail_html, "type:'text/html'", fixed = TRUE)
  expect_match(topic_detail_html, "Genes in full document gene-universe enrichment", fixed = TRUE)
  expect_match(topic_detail_html, "topic genes", fixed = TRUE)
  expect_match(topic_detail_html, "universe remainder", fixed = TRUE)
  expect_match(topic_detail_html, "gene_total_universe", fixed = TRUE)
  expect_match(topic_detail_html, "mdsLeader", fixed = TRUE)
  expect_match(topic_detail_html, "paletteSelect", fixed = TRUE)
  expect_match(topic_detail_html, "mdsLayer", fixed = TRUE)
  expect_match(topic_detail_html, "waterfallLayer", fixed = TRUE)
  expect_match(topic_detail_html, "pathLayer", fixed = TRUE)
  expect_match(topic_detail_html, "mdsTooltip", fixed = TRUE)
  expect_match(topic_detail_html, "body.embed .top", fixed = TRUE)
  expect_match(topic_detail_html, "wfTooltip", fixed = TRUE)
  expect_match(topic_detail_html, "pathTooltip", fixed = TRUE)
  expect_match(topic_detail_html, "mdsStats", fixed = TRUE)
  expect_match(topic_detail_html, "waterfallStats", fixed = TRUE)
  expect_match(topic_detail_html, "pathStats", fixed = TRUE)
  expect_match(topic_detail_html, "Mean document-to-topic probability", fixed = TRUE)
  expect_match(topic_detail_html, "marker", fixed = TRUE)
  expect_match(condition_detail_html, "Condition/Comparison MDS", fixed = TRUE)
  expect_match(condition_detail_html, "Topic Waterfall", fixed = TRUE)
  expect_match(condition_detail_html, "image/svg+xml", fixed = TRUE)
  expect_match(condition_detail_html, "XMLSerializer", fixed = TRUE)
  expect_match(condition_detail_html, "exportSvg", fixed = TRUE)
  expect_no_match(condition_detail_html, "mdsImage", fixed = TRUE)
  expect_no_match(condition_detail_html, "document.documentElement.outerHTML", fixed = TRUE)
  expect_no_match(condition_detail_html, "type:'text/html'", fixed = TRUE)
  expect_match(condition_detail_html, "mdsSvg", fixed = TRUE)
  expect_match(condition_detail_html, "mdsPointHit", fixed = TRUE)
  expect_match(condition_detail_html, "body.embed .top", fixed = TRUE)
  expect_match(condition_detail_html, "function drawMds()", fixed = TRUE)
  expect_no_match(condition_detail_html, "rgba(255,255,255,0.84)", fixed = TRUE)
  expect_no_match(condition_detail_html, "mdsHotspotLayer", fixed = TRUE)
  expect_no_match(condition_detail_html, "function drawMdsHotspots()", fixed = TRUE)
  expect_match(condition_detail_html, "groupColor", fixed = TRUE)
  expect_match(condition_detail_html, "selectGroup", fixed = TRUE)
  expect_match(condition_detail_html, "pathLabelTopicSpecific", fixed = TRUE)
  expect_match(condition_detail_html, "pathLabelGroupSpecific", fixed = TRUE)
  expect_match(condition_detail_html, "pathLabelBothSpecific", fixed = TRUE)
  expect_match(condition_detail_html, "Pathway name colors", fixed = TRUE)
  expect_match(condition_detail_html, "GROUP_MDS", fixed = TRUE)
  expect_match(condition_detail_html, "selectedGroupRows", fixed = TRUE)
  expect_match(condition_detail_html, "Document-to-topic probability", fixed = TRUE)
  expect_match(condition_detail_html, "waterfallLayer", fixed = TRUE)
  expect_match(condition_detail_html, "waterfallStats", fixed = TRUE)
  expect_match(theta_html, "theta_group_mds_k2.png", fixed = TRUE)

  mds_points <- data.table::fread(file.path(csv_dir, "theta_group_mds_points.csv"))
  expect_true("panel_label" %in% names(mds_points))
  expect_true(any(nzchar(mds_points$panel_label)))

  condition_svg <- list.files(
    file.path(root, "review", "condition_topic_reports", "pages", "assets"),
    pattern = "[.]svg$",
    full.names = TRUE
  )
  expect_true(length(condition_svg) >= 1L)

  score_mat <- data.table::fread(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv"))
  expect_equal(score_mat$method_setup, "cond fp aggr weight | LDA")
  expect_true("K2" %in% names(score_mat))
  expect_true("K3" %in% names(score_mat))
  expect_true(is.finite(score_mat$K2[[1L]]))
  expect_true(is.finite(score_mat$K3[[1L]]))

  pass_counts <- data.table::fread(file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  expect_equal(unique(pass_counts$selected_k), 2L)
})

test_that("Module 3 review colors nutrient-stress labels by stress type", {
  labels <- c(
    "HPAFII 0 BCAA Ctrl",
    "HPAFII 25 BCAA TGFb",
    "HPAFII 0.05 Glc Ctrl",
    "HPAFII 0 Gln.Arg Ctrl",
    "HPAFII 12.5uM Met.Cys TGFb",
    "HPAFII 0 Lys Ctrl",
    "HPAFII 0 Trp TGFb",
    "HPAFII 10 Arg Ctrl",
    "HPAFII 5 Gln Ctrl",
    "HPAFII 10 FBS Ctrl"
  )

  expect_equal(
    .m3tb_color_family(labels),
    c("BCAA", "BCAA", "Glc", "Gln.Arg", "Met.Cys", "Lys", "Trp", "Arg", "Gln", "Ctrl")
  )
  expect_gt(length(unique(.m3tb_group_color(labels))), 4L)
})

test_that("Module 3 benchmark infers nutrient-stress metric groups when repeated", {
  comparisons <- data.table::data.table(
    comparison_id = c(
      "HPAFII_0_BCAA_Ctrl_vs_10_FBS_Ctrl",
      "HPAFII_25_BCAA_Ctrl_vs_10_FBS_Ctrl",
      "HPAFII_0_Glc_Ctrl_vs_10_FBS_Ctrl",
      "HPAFII_0_FBS_Ctrl_vs_10_FBS_Ctrl",
      "HPAFII_0p4_FBS_Ctrl_vs_10_FBS_Ctrl"
    ),
    cond1_label = c(
      "HPAFII_0_BCAA_Ctrl",
      "HPAFII_25_BCAA_Ctrl",
      "HPAFII_0_Glc_Ctrl",
      "HPAFII_0_FBS_Ctrl",
      "HPAFII_0p4_FBS_Ctrl"
    ),
    cond2_label = rep("HPAFII_10_FBS_Ctrl", 5),
    comparison_label = c(
      "0_BCAA Ctrl vs 10_FBS Ctrl",
      "25_BCAA Ctrl vs 10_FBS Ctrl",
      "0_Glc Ctrl vs 10_FBS Ctrl",
      "0_FBS Ctrl vs 10_FBS Ctrl",
      "0.4_FBS Ctrl vs 10_FBS Ctrl"
    )
  )

  design <- .m3tb_design_table(comparisons)

  expect_true(any(design$context_type == "condition" & design$metric_group == "BCAA"))
  expect_true(any(design$context_type == "condition" & design$metric_group == "Ctrl"))
  expect_false(any(design$metric_group == "Full"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "BCAA::Target-Up"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "BCAA::Target-Down"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "FBS::Target-Up"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "FBS::Target-Down"))
  expect_false(any(design$context_type == "comparison" & grepl("^Ctrl::", design$metric_group)))
})

test_that("Module 3 topic reports use overlap denominator for pathway universe size", {
  topic_pathways <- data.table::data.table(
    topic = 2L,
    topic_num = 2L,
    pathway = "GO:CC: Intracellular Organelle Lumen",
    pathway_key = "GO:CC: Intracellular Organelle Lumen",
    padj = 0.01,
    overlap = "1/856",
    gene_in = 1L,
    gene_total = 856L,
    genes = "F5"
  )
  universe <- data.table::data.table(
    pathway = "GO:CC: Intracellular Organelle Lumen",
    pathway_key = "GO:CC: Intracellular Organelle Lumen",
    padj = 0.02,
    overlap = "1/856",
    overlap_hits = 1L,
    genes = "F5"
  )

  out <- craftgrn:::.m3tb_apply_universe_pathway_counts(topic_pathways, universe)

  expect_equal(out$gene_in, 1L)
  expect_equal(out$gene_total, 856L)
  expect_equal(out$gene_total_universe, 856L)
  expect_equal(out$gene_out, 855L)
})

test_that("Module 3 theta separation score does not mark singleton groups as perfect", {
  root <- tempfile("module3-theta-singleton-")
  dir.create(root, recursive = TRUE)
  theta <- matrix(
    c(0.90, 0.10, 0.60, 0.40, 0.20, 0.80),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("CondA::TF1", "CondB::TF1", "CondC::TF1"), c("Topic1", "Topic2"))
  )
  row <- data.table::data.table(
    method_order = 1L,
    context_type = "condition",
    setup = "condition_aggr_weight_lda",
    setup_label = "cond fp aggr weight",
    model_label = "LDA",
    method_setup = "cond fp aggr weight | LDA",
    selected_k = 2L,
    fp_mode = "aggregate_weight"
  )
  design <- data.table::data.table(
    context_type = "condition",
    comparison_label = c("CondA", "CondB", "CondC"),
    display_label = c("Cond A", "Cond B", "Cond C"),
    metric_group = c("CondA", "CondB", "CondC")
  )

  scored <- craftgrn:::.m3tb_score_theta_one(theta, row, design, root)

  expect_true(all(is.na(scored$per_label$theta_condition_label_score)))
  expect_equal(scored$score$n_scored_labels, 0L)
  expect_true(is.nan(scored$score$theta_condition_separation_score))
})

test_that("Module 3 shared-topic summary fills missing K values when cache is partial", {
  root <- tempfile("module3-shared-cache-")
  review <- file.path(root, "review")
  dir.create(file.path(review, "tables"), recursive = TRUE)

  make_topic_dir <- function(k) {
    topic_dir <- file.path(root, "topic_extraction", paste0("K", k))
    dir.create(topic_dir, recursive = TRUE)
    rows <- data.table::data.table(
      doc_id = c("CompA::Target-Up::TF1", "CompA::Target-Up::TF1"),
      topic_num = c(1L, k),
      tf = c("TF1", "TF1"),
      peak_id = c("P1", "P2"),
      gene_key = c("G1", "G2"),
      link_pass = TRUE,
      peak_pass = TRUE,
      gene_pass = TRUE
    )
    data.table::fwrite(rows, file.path(topic_dir, "topic_links_pass.csv"))
    data.table::fwrite(
      data.table::data.table(
        link_method = "gammafit",
        output_mode = "pass",
        n_scored_rows = 4L,
        n_pass_rows = 2L,
        pass_file = "topic_links_pass.csv"
      ),
      file.path(topic_dir, "topic_link_summary.csv")
    )
    file.path(topic_dir, "topic_links_pass.csv")
  }

  k2_path <- make_topic_dir(2L)
  k3_path <- make_topic_dir(3L)
  data.table::fwrite(
    data.table::data.table(
      n_topics = 1L,
      n_items = 2L,
      method_order = 1L,
      method_setup = "diff fp aggr | LDA",
      setup = "std_tf_diff_fp_aggr",
      model_label = "LDA",
      selected_k = 2L,
      unit = "Links"
    ),
    file.path(review, "tables", "topic_setup_shared_topic_counts.csv")
  )
  rows <- data.table::data.table(
    method_order = 1L,
    method_setup = "diff fp aggr | LDA",
    setup = "std_tf_diff_fp_aggr",
    model_label = "LDA",
    method = "comparison_aggr_lda",
    selected_k = c(2L, 3L),
    topic_extraction_dir = file.path(root, "topic_extraction"),
    topic_links_path = c(k2_path, k3_path)
  )

  old <- getOption("craftgrn.topic_review.fast_summary")
  options(craftgrn.topic_review.fast_summary = TRUE)
  on.exit(options(craftgrn.topic_review.fast_summary = old), add = TRUE)
  out <- craftgrn:::.m3tb_summarize_topic_links(root, rows, review_dir = review)

  expect_true(any(out$shared$unit == "Links" & out$shared$selected_k == 2L))
  expect_true(any(out$shared$unit == "Links" & out$shared$selected_k == 3L))
  expect_true(any(out$shared$unit == "Genes" & out$shared$selected_k == 3L))
  expect_true(any(out$shared$unit == "TFs" & out$shared$selected_k == 3L))
  expect_true(any(out$pass$count_basis == "Topic-link rows" & out$pass$selected_k == 3L))
})

test_that("Module 3 combined MDS plot separates comparison methods and directions", {
  mds <- data.table::data.table(
    k = 10L,
    method_order = rep(c(1L, 2L), each = 4L),
    method_setup = rep(c("diff fp aggr | LDA", "diff fp uniq | MultiVI"), each = 4L),
    panel_label = rep(c("diff fp aggr\nLDA", "diff fp uniq\nMultiVI"), each = 4L),
    doc_design = "comparison",
    group_label = paste0("G", seq_len(8L)),
    comparison_label = paste0("C", seq_len(8L), rep(c("::Target-Up", "::Target-Down"), 4L)),
    display_label = paste0("Condition ", seq_len(8L), rep(c(" Up", " Down"), 4L)),
    MDS1 = c(-1, -0.8, 0.8, 1, -0.9, -0.7, 0.7, 0.9),
    MDS2 = c(-0.2, 0.2, -0.2, 0.2, -0.3, 0.3, -0.3, 0.3),
    color = rep(c("#4E79A7", "#E15759"), 4L),
    shape_value = rep(c(16L, 25L), 4L)
  )
  p <- craftgrn:::.m3tb_condition_mds_plot(mds, title = "K10 condition/comparison MDS from theta")
  built <- ggplot2::ggplot_build(p)
  panels <- unique(as.character(built$layout$layout$split_panel_label))

  expect_setequal(
    panels,
    c(
      "diff fp aggr\nLDA\nUp",
      "diff fp aggr\nLDA\nDown",
      "diff fp uniq\nMultiVI\nUp",
      "diff fp uniq\nMultiVI\nDown"
    )
  )
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
  expect_equal(res$method_plan$topic_models_dir, file.path(root, "topic_models"))
  expect_equal(res$method_plan$topic_extraction_dir, file.path(root, "topic_extraction"))
  expect_equal(res$review_dir, file.path(root, "review"))
  expect_false(grepl("std_tf_", res$method_plan$topic_models_dir, fixed = TRUE))
  expect_false(grepl("lda$", res$method_plan$topic_models_dir))
})

test_that("Module 3 comparison design prefers concise comparison labels for display", {
  design <- craftgrn:::.m3tb_design_table(
    data.table::data.table(
      comparison_id = "Fib_BIRT_vs_BI",
      comparison_label = "Fib BIRT vs BI",
      comparison_display = "Fibroblast Dox72h BATF IRF4 RUNX3 Tbet vs Fibroblast Dox72h BATF IRF4",
      cond1_label = "Fib_BIRT",
      cond2_label = "Fib_BI"
    )
  )

  comparison_design <- design[context_type == "comparison"]
  expect_equal(
    comparison_design[comparison_label == "Fib_BIRT_vs_BI::Target-Up", display_label],
    "Fib BIRT vs BI Target-Up"
  )
  expect_false(any(grepl("Fibroblast Dox72h", comparison_design$display_label, fixed = TRUE)))
})

test_that("Module 3 standard layout finds root topic-link files", {
  root <- tempfile("module3-standard-topic-links-")
  dir.create(file.path(root, "topic_extraction"), recursive = TRUE)
  topic_links <- file.path(root, "topic_extraction", "topic_links_pass.csv")
  data.table::fwrite(
    data.table::data.table(doc_id = "Doc1", topic_num = 1L, tf = "TF1", peak_id = "P1", gene_key = "G1"),
    topic_links
  )
  row <- .m3tb_apply_output_layout(
    .module3_topic_method_plan(methods = "comparison_aggr_multivi", k_grid = 10L),
    root,
    "standard"
  )
  row[, selected_k := 10L]
  expect_equal(craftgrn:::.m3tb_find_topic_links(root, row[1]), topic_links)
})

test_that("Module 3 review HTML keeps method-specific reports in subfolders", {
  root <- tempfile("module3-topic-review-subfolders-")
  dir.create(root, recursive = TRUE)
  plan <- .m3tb_apply_output_layout(
    .module3_topic_method_plan(
      methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
      k_grid = 2L
    ),
    root,
    "benchmark"
  )

  make_fixture <- function(row) {
    vae_dir <- file.path(row$topic_models_dir[[1L]], "fixture_model", "vae_models")
    dir.create(vae_dir, recursive = TRUE, showWarnings = FALSE)
    theta <- data.table::data.table(
      doc_id = c("CondA::TF1", "CondB::TF1"),
      Topic1 = c(0.8, 0.2),
      Topic2 = c(0.2, 0.8)
    )
    phi <- data.table::data.table(
      term_id = c("GENE:G1", "GENE:G2", "PEAK:P1"),
      Topic1 = c(0.7, 0.2, 0.1),
      Topic2 = c(0.1, 0.7, 0.2)
    )
    data.table::fwrite(theta, file.path(vae_dir, "theta_K2.csv"))
    data.table::fwrite(phi, file.path(vae_dir, "phi_K2.csv"))
  }
  for (i in seq_len(nrow(plan))) make_fixture(plan[i])

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
    run_reports = TRUE,
    verbose = FALSE
  )

  review_dir <- res$review_dir
  expect_true(file.exists(file.path(review_dir, "topic_method_k_topic_mds_report.html")))
  expect_true(file.exists(file.path(review_dir, "topic_method_k_condition_mds_report.html")))
  expect_false(file.exists(file.path(review_dir, "topic_report_K2.html")))
  expect_false(file.exists(file.path(review_dir, "condition_topic_report_K2.html")))

  topic_reports <- list.files(file.path(review_dir, "topic_reports"), pattern = "_topic_report[.]html$", full.names = TRUE)
  condition_reports <- list.files(file.path(review_dir, "condition_topic_reports"), pattern = "_condition_topic_report[.]html$", full.names = TRUE)
  expect_equal(length(topic_reports), 2L)
  expect_equal(length(condition_reports), 2L)
  expect_equal(length(unique(basename(topic_reports))), 2L)
  expect_equal(length(unique(basename(condition_reports))), 2L)

  topic_index <- paste(readLines(file.path(review_dir, "topic_method_k_topic_mds_report.html"), warn = FALSE), collapse = "\n")
  condition_index <- paste(readLines(file.path(review_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  expect_match(topic_index, "topic_reports/", fixed = TRUE)
  expect_match(condition_index, "condition_topic_reports/", fixed = TRUE)
  expect_match(topic_index, "Method <select", fixed = TRUE)
  expect_match(condition_index, "Method <select", fixed = TRUE)
  expect_match(topic_index, "K <select", fixed = TRUE)
  expect_match(condition_index, "K <select", fixed = TRUE)
})

test_that("Module 3 review HTML groups multiple K values inside one method report", {
  root <- tempfile("module3-topic-review-k-dropdown-")
  dir.create(root, recursive = TRUE)
  plan <- .m3tb_apply_output_layout(
    .module3_topic_method_plan(
      methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
      k_grid = c(2L, 3L)
    ),
    root,
    "benchmark"
  )

  make_fixture <- function(row) {
    vae_dir <- file.path(row$topic_models_dir[[1L]], "fixture_model", "vae_models")
    dir.create(vae_dir, recursive = TRUE, showWarnings = FALSE)
    for (k in c(2L, 3L)) {
      theta <- data.table::data.table(
        doc_id = c("CondA::TF1", "CondB::TF1"),
        Topic1 = c(0.75, 0.25),
        Topic2 = c(0.25, 0.75)
      )
      phi <- data.table::data.table(
        term_id = c("GENE:G1", "GENE:G2", "PEAK:P1"),
        Topic1 = c(0.7, 0.2, 0.1),
        Topic2 = c(0.1, 0.7, 0.2)
      )
      if (k == 3L) {
        theta[, Topic3 := c(0.05, 0.05)]
        theta[, `:=`(Topic1 = Topic1 - 0.025, Topic2 = Topic2 - 0.025)]
        phi[, Topic3 := c(0.2, 0.2, 0.6)]
      }
      data.table::fwrite(theta, file.path(vae_dir, sprintf("theta_K%d.csv", k)))
      data.table::fwrite(phi, file.path(vae_dir, sprintf("phi_K%d.csv", k)))
    }
  }
  for (i in seq_len(nrow(plan))) make_fixture(plan[i])

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      condition_label = c("CondA", "CondB"),
      condition_group = c("CondA", "CondB")
    ),
    methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    k_grid = c(2L, 3L),
    output_layout = "benchmark",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = TRUE,
    verbose = FALSE
  )

  review_dir <- res$review_dir
  topic_reports <- list.files(file.path(review_dir, "topic_reports"), pattern = "_topic_report[.]html$", full.names = TRUE)
  condition_reports <- list.files(file.path(review_dir, "condition_topic_reports"), pattern = "_condition_topic_report[.]html$", full.names = TRUE)
  expect_equal(length(topic_reports), 2L)
  expect_equal(length(condition_reports), 2L)
  expect_equal(length(list.files(file.path(review_dir, "topic_reports", "pages"), pattern = "_topic_report[.]html$")), 4L)
  expect_equal(length(list.files(file.path(review_dir, "condition_topic_reports", "pages"), pattern = "_condition_topic_report[.]html$")), 4L)

  topic_html <- paste(readLines(topic_reports[[1L]], warn = FALSE), collapse = "\n")
  condition_html <- paste(readLines(condition_reports[[1L]], warn = FALSE), collapse = "\n")
  expect_match(topic_html, "K <select", fixed = TRUE)
  expect_match(topic_html, "K2", fixed = TRUE)
  expect_match(topic_html, "K3", fixed = TRUE)
  expect_match(topic_html, "pages/", fixed = TRUE)
  expect_match(condition_html, "K <select", fixed = TRUE)
  expect_match(condition_html, "K2", fixed = TRUE)
  expect_match(condition_html, "K3", fixed = TRUE)
  expect_match(condition_html, "pages/", fixed = TRUE)
})

test_that("Module 3 standard layout flattens one selected K extraction", {
  root <- tempfile("module3-topic-standard-single-k-")
  dir.create(root, recursive = TRUE)

  res <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      condition_label = c("CondA", "CondB"),
      condition_group = c("CondA", "CondB")
    ),
    methods = "condition_aggr_weight_lda",
    k_grid = 10L,
    output_layout = "standard",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = FALSE,
    verbose = FALSE
  )

  extract_roots <- .m3tb_extraction_output_dirs(res$method_plan[1], 10L, res$output_layout)
  expect_equal(extract_roots[[1L]], file.path(root, "topic_extraction"))

  res_multi_k <- run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      condition_label = c("CondA", "CondB"),
      condition_group = c("CondA", "CondB")
    ),
    methods = "condition_aggr_weight_lda",
    k_grid = c(10L, 12L),
    output_layout = "standard",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = FALSE,
    verbose = FALSE
  )
  extract_roots <- .m3tb_extraction_output_dirs(res_multi_k$method_plan[1], c(10L, 12L), res_multi_k$output_layout)
  expect_equal(extract_roots, file.path(root, "topic_extraction", c("K10", "K12")))
})

test_that("Module 3 topic benchmark exposes project-agnostic sample subsetting", {
  expect_false("celllines" %in% names(formals(run_module3_topic_benchmark)))
  expect_true("sample_subset" %in% names(formals(run_module3_topic_benchmark)))
  expect_true("analysis_label" %in% names(formals(run_module3_topic_benchmark)))
  expect_false("celllines" %in% names(formals(train_topic_models)))
  expect_true("sample_subset" %in% names(formals(train_topic_models)))
  expect_true("analysis_label" %in% names(formals(train_topic_models)))
  expect_true("flat_output" %in% names(formals(module3_train_topic_models)))
  expect_true(isTRUE(formals(module3_train_topic_models)$flat_output))
  expect_true("flatten_single_output" %in% names(formals(module3_extract_topics)))
  expect_true(isTRUE(formals(module3_extract_topics)$flatten_single_output))
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

test_that("Module 3 topic benchmark applies readable comparison labels", {
  theta <- matrix(
    c(
      0.8, 0.2,
      0.1, 0.9
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(theta) <- c(
    "Fib_BIRT_vs_BI::TF1::Target-Up",
    "Fib_BIRT_vs_BI::TF2::Target-Down"
  )
  colnames(theta) <- c("Topic1", "Topic2")
  design <- .m3tb_design_table(data.table::data.table(
    comparison_id = "Fib_BIRT_vs_BI",
    comparison_display = "Fibroblast BIRT vs BI",
    cond1_label = "Dox72h_BIRT",
    cond2_label = "Dox72h_BI"
  ))
  row <- data.table::data.table(
    method_order = 1L,
    context_type = "comparison",
    setup = "comparison_aggr_multivi",
    setup_label = "diff fp aggr",
    model_label = "MultiVI",
    method_setup = "diff fp aggr | MultiVI",
    selected_k = 2L,
    fp_mode = "aggregate"
  )
  csv_dir <- tempfile("theta-labels-")
  dir.create(csv_dir, recursive = TRUE)
  out <- .m3tb_score_theta_one(theta, row, design, csv_dir = csv_dir)

  expect_equal(
    out$per_label$comparison_label,
    c("Fib_BIRT_vs_BI::Target-Down", "Fib_BIRT_vs_BI::Target-Up")
  )
  expect_equal(
    out$per_label$display_label,
    c("Fibroblast BIRT vs BI Target-Down", "Fibroblast BIRT vs BI Target-Up")
  )
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

  stale_summary <- data.table::fread(file.path(out_dir, "topic_input_summary.csv"))
  stale_summary[, `:=`(doc_design = "condition", fp_term_mode = "aggregate_weight")]
  data.table::fwrite(stale_summary, file.path(out_dir, "topic_input_summary.csv"))
  rebuilt <- module3_construct_docs(
    filtered_dir = filtered_dir,
    output_dir = out_dir,
    tf_cluster_map = c(TF1 = "K01"),
    doc_mode = "tf",
    doc_design = "comparison",
    fp_term_mode = "aggregate",
    min_df = 1,
    threshold_gene_expr = 0,
    threshold_fp_score = 0,
    threshold_tf_expr = 0,
    direction_consistency = "none",
    require_tf_expr_either = FALSE,
    require_gene_expr_either = FALSE,
    require_fp_bound_either = FALSE,
    overwrite = FALSE,
    verbose = FALSE
  )
  expect_false(rebuilt$reused)
  expect_equal(rebuilt$summary$doc_design[[1L]], "comparison")
  expect_equal(rebuilt$summary$fp_term_mode[[1L]], "aggregate")
})

test_that("Module 3 production wrapper exposes compact defaults and QC report", {
  expect_true("run_topic_modeling" %in% getNamespaceExports("craftgrn"))
  expect_true("module3_construct_docs" %in% getNamespaceExports("craftgrn"))
  expect_true("build_module3_qc_report" %in% getNamespaceExports("craftgrn"))
  expect_true("project_config" %in% names(formals(run_topic_modeling)))
  expect_true("warplda_iterations" %in% names(formals(run_topic_modeling)))
  expect_true("vae_device" %in% names(formals(run_topic_modeling)))
  expect_true("vae_device" %in% names(formals(train_topic_models)))
  expect_true("vae_batch_size" %in% names(formals(train_topic_models)))
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

test_that("Module 3 topic wrapper resolves standard run settings from project config", {
  cfg <- list(
    topic_method = "comparison_aggr_weight_lda",
    topic_k_grid = c(8L, 10L),
    warplda_iterations = 25L,
    topic_link_output = "none",
    topic_count_method = "bin",
    topic_count_input = "pseudo_count_bin",
    topic_vae_device = "cuda",
    topic_vae_batch_size = 512L,
    topic_benchmark_enabled = TRUE,
    topic_benchmark_methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    topic_benchmark_k_grid = c(5L, 6L)
  )
  resolved <- .module3_resolve_topic_run_config(project_config = cfg)
  expect_equal(resolved$method, "comparison_aggr_weight_lda")
  expect_equal(resolved$k_grid, c(8L, 10L))
  expect_equal(resolved$warplda_iterations, 25L)
  expect_equal(resolved$topic_link_output, "none")
  expect_equal(resolved$count_method, "bin")
  expect_equal(resolved$count_input, "pseudo_count_bin")
  expect_equal(resolved$vae_device, "cuda")
  expect_equal(resolved$vae_batch_size, 512L)
  expect_true(resolved$benchmark$enabled)
  expect_equal(resolved$benchmark$methods, cfg$topic_benchmark_methods)
  expect_equal(resolved$benchmark$k_grid, c(5L, 6L))

  overridden <- .module3_resolve_topic_run_config(
    project_config = cfg,
    method = "condition_aggr_lda",
    k_grid = 12L,
    warplda_iterations = 3L,
    topic_link_output = "pass",
    count_method = "log",
    count_input = "pseudo_count_log",
    vae_device = "cpu",
    vae_batch_size = 128L
  )
  expect_equal(overridden$method, "condition_aggr_lda")
  expect_equal(overridden$k_grid, 12L)
  expect_equal(overridden$warplda_iterations, 3L)
  expect_equal(overridden$topic_link_output, "pass")
  expect_equal(overridden$count_method, "log")
  expect_equal(overridden$count_input, "pseudo_count_log")
  expect_equal(overridden$vae_device, "cpu")
  expect_equal(overridden$vae_batch_size, 128L)
})

test_that("Module 3 topic link defaults do not apply gene_prob max filtering", {
  resolved_default <- .module3_resolve_topic_run_config(project_config = list())
  expect_equal(resolved_default$extraction_args$link_topic_method, "gammafit")
  expect_equal(resolved_default$extraction_args$link_topic_prob_cutoff, 0.3)

  resolved_explicit <- .module3_resolve_topic_run_config(
    project_config = list(
      topic_link_method = "gene_prob",
      topic_link_prob_cutoff = "max"
    )
  )
  expect_equal(resolved_explicit$extraction_args$link_topic_method, "gene_prob")
  expect_equal(resolved_explicit$extraction_args$link_topic_prob_cutoff, "max")
})

test_that("Module 3 VAE device supports explicit auto mode", {
  resolved_default <- .module3_resolve_topic_run_config(project_config = list())
  expect_equal(resolved_default$vae_device, "auto")

  resolved_auto <- .module3_resolve_topic_run_config(
    project_config = list(topic_vae_device = "cuda"),
    vae_device = "auto"
  )
  expect_equal(resolved_auto$vae_device, "auto")

  py <- paste(readLines(system.file("python", "logistic_normal_vae_topics.py", package = "craftgrn"), warn = FALSE), collapse = "\n")
  expect_match(py, 'choices=["cpu", "cuda", "auto"]', fixed = TRUE)
  expect_match(py, "def _resolve_device", fixed = TRUE)
  expect_match(py, "requested_device", fixed = TRUE)
  expect_match(py, "resolved_device", fixed = TRUE)
})

test_that("Module 3 high-throughput topic input skips repeated-value diagnostics by default", {
  expect_false(formals(train_topic_models)$check_repeated_values)
  expect_false(formals(module3_construct_docs)$check_repeated_values)
  expect_true(formals(build_doc_term_joint)$check_repeated_values)
  expect_true(formals(build_doc_term_condition_union)$check_repeated_values)
})

test_that("Module 3 standard extraction keeps per-comparison pathway reports optional", {
  model_defaults <- paste(
    deparse(body(extract_regulatory_topics)),
    collapse = "\n"
  )
  expect_match(model_defaults, "pathway_per_comparison = FALSE", fixed = TRUE)
  expect_match(model_defaults, "run_pathway_enrichment = FALSE", fixed = TRUE)
})

test_that("Module 3 theta review PNG writer is headless safe", {
  skip_if_not(capabilities("cairo"))
  old_bitmap <- getOption("bitmapType")
  withr::defer(options(bitmapType = old_bitmap))
  options(bitmapType = "Xlib")
  root <- tempfile("module3-review-png-")
  score_result <- list(
    mds_points = data.table::data.table(
      k = 2L,
      method_setup = "diff fp aggr | MultiVI",
      display_label = c("Cond A", "Cond B"),
      MDS1 = c(-0.1, 0.1),
      MDS2 = c(0.2, -0.2),
      color = c("#4E79A7", "#E15759")
    ),
    scores = data.table::data.table(
      k = 2L,
      method_setup = "diff fp aggr | MultiVI",
      theta_condition_separation_score = 0.5
    )
  )
  out <- craftgrn:::.m3tb_write_review_pngs(score_result, root)
  expect_equal(nrow(out), 1L)
  expect_true(file.exists(file.path(root, "theta_phi_topic_distance_correlation_k2.png")))
  expect_true(file.exists(file.path(root, "theta_group_mds_k2.png")))
})
