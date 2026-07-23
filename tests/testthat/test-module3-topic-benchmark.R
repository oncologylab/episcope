test_that("Module 3 topic benchmark scores existing models and writes review reports", {
  old_scan <- getOption("craftgrn.topic_review.scan_link_tables")
  options(craftgrn.topic_review.scan_link_tables = TRUE)
  on.exit(options(craftgrn.topic_review.scan_link_tables = old_scan), add = TRUE)

  root <- tempfile("module3-topic-benchmark-")
  dir.create(root, recursive = TRUE)
  stale_review <- file.path(root, "review")
  dir.create(file.path(stale_review, "topic_reports"), recursive = TRUE)
  dir.create(
    file.path(stale_review, "condition_topic_reports", "pages", "assets"),
    recursive = TRUE
  )
  writeLines("stale", file.path(stale_review, "theta_phi_and_group_mds.html"))
  writeLines("stale", file.path(stale_review, "topic_method_k_topic_mds_report.html"))
  writeLines("stale", file.path(stale_review, "tf_std_six_setups_pass_state_counts.pdf"))
  writeLines("stale", file.path(stale_review, "topic_reports", "stale.html"))
  writeLines(
    "stale",
    file.path(stale_review, "condition_topic_reports", "pages", "assets", "stale.svg")
  )

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
  data.table::fwrite(
    data.table::data.table(
      topic = c(1L, 1L, 2L, 2L),
      term_id = c("GENE:G1", "GENE:G2", "GENE:G3", "PEAK:P1"),
      score = c(0.9, 0.8, 0.01, 0.02),
      in_topic = c(TRUE, TRUE, FALSE, FALSE)
    ),
    file.path(extraction_dir, "topic_terms.csv")
  )
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
  expect_true(file.exists(file.path(csv_dir, "topic_setup_item_coverage_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_link_assignment_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "topic_setup_shared_topic_counts.csv")))
  expect_true(file.exists(file.path(csv_dir, "theta_group_mds_points.csv")))
  expect_false(file.exists(file.path(root, "review", "tf_std_six_setups_pass_state_counts.pdf")))
  expect_false(file.exists(file.path(root, "review", "tf_std_six_setups_shared_topic_counts.pdf")))
  expect_false(file.exists(file.path(html_dir, "theta_phi_and_group_mds.html")))
  expect_false(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report.html")))
  expect_false(dir.exists(file.path(html_dir, "condition_topic_reports", "pages", "assets")))
  expect_true(file.exists(file.path(html_dir, "topic_method_k_condition_mds_report.html")))
  expect_false(file.exists(file.path(html_dir, "topic_method_k_topic_mds_report_global_term_group.html")))
  expect_false(file.exists(file.path(html_dir, "topic_method_k_condition_mds_report_global_term_group.html")))

  pass_counts <- data.table::fread(file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  item_counts <- data.table::fread(file.path(csv_dir, "topic_setup_item_coverage_counts.csv"))
  link_counts <- data.table::fread(file.path(csv_dir, "topic_setup_link_assignment_counts.csv"))
  shared_counts <- data.table::fread(file.path(csv_dir, "topic_setup_shared_topic_counts.csv"))
  expect_setequal(pass_counts$status, c("Pass", "Fail"))
  expect_setequal(pass_counts$unit, c("Terms", "Genes", "TF-gene-doc links", "Links", "TFs"))
  expect_equal(pass_counts[unit == "Terms" & status == "Pass", count], 2)
  expect_equal(pass_counts[unit == "Terms" & status == "Fail", count], 2)
  expect_equal(pass_counts[unit == "Terms" & status == "Pass", percent], 50)
  expect_equal(pass_counts[unit == "Genes" & status == "Pass", count], 4)
  expect_equal(pass_counts[unit == "Genes" & status == "Fail", count], 0)
  expect_equal(pass_counts[unit == "Genes" & status == "Pass", percent], 100)
  expect_equal(pass_counts[unit == "Links" & status == "Pass", count], 4)
  expect_equal(pass_counts[unit == "Links" & status == "Fail", count], 0)
  expect_equal(pass_counts[unit == "Links" & status == "Pass", percent], 100)
  expect_equal(pass_counts[unit == "TF-gene-doc links" & status == "Pass", count], 4)
  expect_equal(pass_counts[unit == "TF-gene-doc links" & status == "Fail", count], 0)
  expect_equal(pass_counts[unit == "TF-gene-doc links" & status == "Pass", percent], 100)
  expect_equal(pass_counts[unit == "TFs" & status == "Pass", count], 2)
  expect_equal(pass_counts[unit == "TFs" & status == "Fail", count], 0)
  expect_equal(pass_counts[unit == "TFs" & status == "Pass", percent], 100)
  expect_equal(item_counts, pass_counts, ignore_attr = TRUE)
  expect_equal(sum(link_counts[status == "Pass", count]), 4)
  expect_equal(sum(link_counts[status == "Fail", count]), 0)
  expect_equal(unique(link_counts$count_basis), "Standard TF-FP-gene links")
  expect_true("Pathways" %in% shared_counts$unit)

  condition_html <- paste(readLines(file.path(html_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  expect_match(condition_html, "Condition 1 <select", fixed = TRUE)
  expect_match(condition_html, "Condition 2 <select", fixed = TRUE)
  expect_match(condition_html, "Pathway <select", fixed = TRUE)
  expect_no_match(condition_html, "Report <select", fixed = TRUE)
  expect_match(condition_html, "Method <select", fixed = TRUE)
  expect_match(condition_html, "K <select", fixed = TRUE)
  expect_match(condition_html, "assets/p/", fixed = TRUE)
  expect_match(condition_html, "frame.src=embed(hit.src)", fixed = TRUE)
  expect_match(condition_html, "embed=1", fixed = TRUE)
  expect_match(
    condition_html,
    "options.find(x=>Number(x.trained_k)===trained)",
    fixed = TRUE
  )
  expect_match(condition_html, 'id="cond1Color"', fixed = TRUE)
  expect_match(condition_html, 'id="cond2Color"', fixed = TRUE)
  expect_no_match(condition_html, "Current selection", fixed = TRUE)
  expect_no_match(condition_html, "srcdoc", fixed = TRUE)

  condition_page_files <- list.files(
    file.path(html_dir, "assets", "p"),
    pattern = "[.]html$",
    full.names = TRUE
  )
  expect_false(dir.exists(file.path(html_dir, "topic_reports")))
  expect_false(dir.exists(file.path(html_dir, "condition_topic_reports")))
  expect_true(length(condition_page_files) >= 1L)
  condition_detail_html <- paste(readLines(condition_page_files[[1L]], warn = FALSE), collapse = "\n")
  expect_match(condition_detail_html, "Condition/Comparison MDS", fixed = TRUE)
  expect_match(condition_detail_html, "Topic Bar Plots", fixed = TRUE)
  expect_match(condition_detail_html, "craftgrn-module3-report-schema\" content=\"2", fixed = TRUE)
  expect_no_match(condition_detail_html, "Topic Waterfall", fixed = TRUE)
  expect_match(condition_detail_html, "id=\"tfSelect\"", fixed = TRUE)
  expect_match(condition_detail_html, "XMLSerializer", fixed = TRUE)
  expect_match(condition_detail_html, "exportSvg", fixed = TRUE)
  expect_match(condition_detail_html, "mdsSvg", fixed = TRUE)
  expect_match(condition_detail_html, "body.embed .top", fixed = TRUE)
  expect_match(condition_detail_html, "function drawMds()", fixed = TRUE)
  expect_match(condition_detail_html, "pathLabelTopicSpecific", fixed = TRUE)
  expect_match(condition_detail_html, "subgrnPanel", fixed = TRUE)
  expect_match(condition_detail_html, "Back to pathways", fixed = TRUE)
  expect_match(condition_detail_html, "GROUP_MDS", fixed = TRUE)
  expect_match(condition_detail_html, "TF_TOPIC", fixed = TRUE)

  mds_points <- data.table::fread(file.path(csv_dir, "theta_group_mds_points.csv"))
  expect_true("panel_label" %in% names(mds_points))
  expect_true(any(nzchar(mds_points$panel_label)))

  condition_svg <- list.files(
    file.path(root, "review", "assets", "p"),
    pattern = "[.]svg$",
    full.names = TRUE
  )
  expect_length(condition_svg, 0L)

  score_mat <- data.table::fread(file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv"))
  expect_equal(score_mat$method_setup, "cond fp aggr weight | LDA")
  expect_true("K2" %in% names(score_mat))
  expect_true("K3" %in% names(score_mat))
  expect_true(is.finite(score_mat$K2[[1L]]))
  expect_true(is.finite(score_mat$K3[[1L]]))

  pass_counts <- data.table::fread(file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  expect_equal(unique(pass_counts$selected_k), 2L)
})

test_that("Module 3 review uses explicit comparison groups instead of label substring rules", {
  comparisons <- data.table::data.table(
    comparison_id = c(
      "TCell_Resting_vs_Naive",
      "TCell_Restimulated_vs_Naive",
      "Fib_BATF_IRF4_vs_NoTF"
    ),
    cond1_label = c(
      "WT_D6_Resting",
      "WT_D6_Restimulated",
      "Dox72h_BATF_IRF4"
    ),
    cond2_label = c("Naive_Flox", "Naive_Flox", "Dox72h_NoTFs"),
    comparison_group = c(
      "TCell_D6_resting_reference",
      "TCell_D6_restimulated_reference",
      "Fib_Combo_OE"
    ),
    comparison_label = c(
      "WT D6 Resting vs Naive Flox",
      "WT D6 Restimulated vs Naive Flox",
      "Dox72h BATF IRF4 vs Dox72h NoTFs"
    )
  )

  design <- .m3tb_design_table(comparisons)

  expect_true(any(design$context_type == "comparison" & design$metric_group == "TCell_D6_resting_reference::Target-Up"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "TCell_D6_resting_reference::Target-Down"))
  expect_true(any(design$context_type == "comparison" & design$metric_group == "Fib_Combo_OE::Target-Up"))
  expect_false(any(design$metric_group == "Arg"))
  expect_false(any(grepl("^Arg::", design$metric_group)))
})

test_that("Module 3 review colors are generic and do not classify Target as Arg", {
  groups <- c(
    "TCell_D6_resting_reference::Target-Up",
    "TCell_D6_resting_reference::Target-Down",
    "Fib_Combo_OE::Target-Up",
    "Fib_Combo_OE::Target-Down"
  )

  expect_equal(.m3tb_color_family(groups), groups)
  expect_equal(length(unique(.m3tb_group_color(groups))), length(groups))
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

test_that("Module 3 condition reports fall back to topic-level pathways", {
  root <- tempfile("module3-condition-pathway-fallback-")
  extraction_dir <- file.path(root, "topic_extraction", "K2")
  dir.create(extraction_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(
      topic = c(1L, 2L),
      pathway = c("Reactome: Signal A", "GO:BP: Signal B"),
      pathway_key = c("Reactome: Signal A", "GO:BP: Signal B"),
      padj = c(0.01, 0.02),
      overlap = c("6/30", "7/40"),
      overlap_hits = c(6L, 7L),
      genes = c("G1;G2;G3;G4;G5;G6", "G7;G8;G9;G10;G11;G12;G13")
    ),
    file.path(extraction_dir, "topic_pathway_enrichment_peak_and_gene_dotplot.csv")
  )

  pathways <- craftgrn:::.m3tb_read_condition_pathway_tables(
    extraction_dir,
    model_dir = NULL,
    compute_universe = FALSE
  )

  expect_equal(nrow(pathways), 2L)
  expect_true(all(is.na(pathways$comparison_label)))
  expect_setequal(pathways$pathway, c("Reactome: Signal A", "GO:BP: Signal B"))
})

test_that("Module 3 condition reports read full per-comparison pathway tables", {
  root <- tempfile("module3-condition-pathway-per-comparison-")
  extraction_dir <- file.path(root, "topic_extraction", "K2")
  dir.create(extraction_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(
      comparison_id = c("CmpA", "CmpA"),
      direction_group = c("Up", "Up"),
      topic = c(1L, 1L),
      pathway = c("Reactome: Signal A", "GO:BP: Signal B"),
      pathway_norm_key = c("reactome:signal a", "go bp:signal b"),
      padj = c(0.01, 1),
      overlap = c("6/30", "0/40"),
      overlap_hits = c(6L, 0L),
      genes = c("G1;G2;G3;G4;G5;G6", "")
    ),
    file.path(extraction_dir, "per_comparison_topic_pathway_enrichment.csv")
  )

  pathways <- craftgrn:::.m3tb_read_condition_pathway_tables(
    extraction_dir,
    model_dir = NULL,
    compute_universe = FALSE
  )

  expect_equal(nrow(pathways), 1L)
  expect_setequal(pathways$pathway, "Reactome: Signal A")
  expect_setequal(pathways$comparison_label, "CmpA::Target-Up")
  expect_false("GO:BP: Signal B" %in% pathways$pathway)

  out_file <- tempfile(fileext = ".html")
  craftgrn:::.m3tb_condition_report_html(
    title = "Condition pathway test",
    group_mds = data.table::data.table(
      comparison_label = "CmpA::Target-Up",
      display_label = "CmpA Up",
      group_label = "CmpA Up",
      MDS1 = 0,
      MDS2 = 0,
      n_docs = 2L
    ),
    group_topic = data.table::data.table(
      comparison_label = "CmpA::Target-Up",
      display_label = "CmpA Up",
      n_docs = 2L,
      topic = "Topic1",
      topic_num = 1L,
      theta_mean = 0.8
    ),
    pathways = pathways,
    out_html = out_file
  )
  html <- paste(readLines(out_file, warn = FALSE), collapse = "\n")

  expect_no_match(html, "GO:BP: Signal B", fixed = TRUE)
  expect_no_match(html, "Filter: N_gene >= 5, adjusted p-value < 0.05", fixed = TRUE)
  expect_match(html, "Full per-condition or per-comparison pathway table", fixed = TRUE)
})

test_that("Module 3 review builds TF-topic rows from theta document IDs", {
  theta_condition <- matrix(
    c(0.8, 0.2, 0.3, 0.7),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("CondA::TF1", "CondA::TF2"), c("Topic1", "Topic2"))
  )
  condition_rows <- craftgrn:::.m3tb_tf_topic_rows(theta_condition, "condition")
  expect_setequal(condition_rows$comparison_label, "CondA")
  expect_setequal(condition_rows$tf, c("TF1", "TF2"))
  expect_equal(condition_rows[tf == "TF1" & topic == "Topic1", theta], 0.8)

  theta_comparison <- matrix(
    c(0.1, 0.9),
    nrow = 1L,
    dimnames = list("CmpA::TFX::Target-Up", c("Topic1", "Topic2"))
  )
  comparison_rows <- craftgrn:::.m3tb_tf_topic_rows(theta_comparison, "comparison")
  expect_setequal(comparison_rows$comparison_label, "CmpA::Target-Up")
  expect_setequal(comparison_rows$tf, "TFX")
  expect_equal(comparison_rows[topic == "Topic2", theta], 0.9)
})

test_that("Module 3 pathway sub-GRN payloads support TF-gene and TF-peak-gene views", {
  pathways <- data.table::data.table(
    comparison_id = "CmpA",
    direction_group = "Up",
    topic = 2L,
    pathway = "Reactome: Signal A",
    pathway_norm_key = "reactome:signal a",
    padj = 0.01,
    overlap_genes = "GeneA;GeneB"
  )
  edges_docs <- data.table::data.table(
    comparison_id = "CmpA",
    direction = "Target-Up",
    tf = c("TF1", "TF1", "TF2", "TF3", "TF1"),
    gene_key = c("GeneA", "GeneA", "GeneA", "GeneB", "GeneC"),
    peak_id = c("PeakShared", "Peak2", "PeakShared", "PeakShared", "Peak3"),
    delta_fp = c(2, 1, 3, 0.5, 9),
    delta_gene = c(1, 1, 1, 1, 1),
    log2fc_gene = c(1.2, 0.8, 1.5, 0.4, 2)
  )
  tf_membership <- data.table::data.table(
    comparison_id = "CmpA",
    direction = "Target-Up",
    tf = c("TF1", "TF3"),
    topic_num = c(2L, 1L),
    membership_pass = c(TRUE, TRUE)
  )

  payload <- craftgrn:::.m3tb_build_pathway_subgrn_payload(
    pathways = pathways,
    edges_docs = edges_docs,
    tf_membership = tf_membership,
    max_tf_gene_edges_per_context = 20L,
    max_tf_peak_gene_triplets_per_context = 20L
  )

  expect_equal(nrow(payload$manifest), 1L)
  expect_equal(payload$manifest$direction[[1L]], "Target-Up")
  expect_equal(payload$manifest$n_overlap_genes[[1L]], 2L)
  expect_equal(payload$manifest$n_tf_gene_edges[[1L]], 3L)
  expect_equal(payload$manifest$n_tf_peak_gene_triplets[[1L]], 4L)

  tf_gene <- payload$tf_gene_edges
  expect_setequal(tf_gene$gene_key, c("GeneA", "GeneB"))
  expect_false("GeneC" %in% tf_gene$gene_key)
  expect_equal(tf_gene[tf == "TF1" & gene_key == "GeneA", n_supporting_peaks], 2L)
  expect_true(tf_gene[tf == "TF1" & gene_key == "GeneA", topic_tf])
  expect_false(tf_gene[tf == "TF2" & gene_key == "GeneA", topic_tf])

  triplets <- payload$tf_peak_gene_triplets
  expect_equal(triplets[peak_id == "PeakShared", uniqueN(tf)], 3L)
  expect_equal(triplets[peak_id == "PeakShared", uniqueN(gene_key)], 2L)
  expect_true(all(triplets$gene_key %in% c("GeneA", "GeneB")))
})

test_that("Module 3 pathway sub-GRN contexts cover plotted pathway rows", {
  all_pathways <- data.table::data.table(
    comparison_id = c("CmpA", "CmpA", "CmpB", "CmpB"),
    direction_group = c("Up", "Up", "Down", "Down"),
    topic = c(1L, 1L, 1L, 2L),
    pathway = c("Path A", "Path B", "Path B", "Path C"),
    pathway_key = c("path:a", "path:b", "path:b", "path:c"),
    padj = c(0.01, 0.02, 0.03, 0.04),
    overlap_genes = c("GeneA", "GeneB", "GeneC", "GeneD")
  )
  condition_pathways <- data.table::data.table(
    comparison_id = "CmpA",
    direction_group = "Target-Up",
    topic = 1L,
    topic_num = 1L,
    pathway = "Path A",
    pathway_key = "path:a"
  )
  topic_pathways <- data.table::data.table(
    topic = 1L,
    topic_num = 1L,
    pathway = "Path B",
    pathway_key = "path:b"
  )

  selected <- craftgrn:::.m3tb_select_pathway_subgrn_contexts(
    all_pathways = all_pathways,
    condition_pathways = condition_pathways,
    topic_pathways = topic_pathways
  )

  selected_keys <- paste(selected$comparison_id, selected$direction_group, selected$topic, selected$pathway_key, sep = "|")
  expect_setequal(
    selected_keys,
    c("CmpA|Up|1|path:a", "CmpA|Up|1|path:b", "CmpB|Down|1|path:b")
  )
  expect_false(any(selected$pathway_key == "path:c"))
})

test_that("Module 3 condition report HTML exposes pathway sub-GRN controls", {
  out_file <- tempfile(fileext = ".html")
  craftgrn:::.m3tb_condition_report_html(
    title = "Condition subgrn test",
    group_mds = data.table::data.table(
      comparison_label = "CmpA::Target-Up",
      display_label = "CmpA Up",
      group_label = "CmpA Up",
      MDS1 = 0,
      MDS2 = 0,
      n_docs = 2L
    ),
    group_topic = data.table::data.table(
      comparison_label = "CmpA::Target-Up",
      display_label = "CmpA Up",
      n_docs = 2L,
      topic = "Topic2",
      topic_num = 2L,
      theta_mean = 0.8
    ),
    pathways = data.table::data.table(
      pathway_key = "reactome:signal a",
      pathway = "Reactome: Signal A",
      topic = 2L,
      topic_num = 2L,
      comparison_label = "CmpA::Target-Up",
      display_label = "CmpA Up",
      padj = 0.01,
      overlap = "2/30",
      gene_in = 2L,
      gene_total = 30L,
      gene_out = 28L,
      gene_total_universe = 30L,
      genes = "GeneA;GeneB"
    ),
    out_html = out_file,
    subgrn_manifest = data.table::data.table(
      subgrn_context_id = "ctx1",
      comparison_label = "CmpA::Target-Up",
      topic_num = 2L,
      pathway_key = "reactome:signal a",
      pathway = "Reactome: Signal A",
      payload_file = "pathway_subgrn_payloads/chunk_001.js"
    ),
    subgrn_payload_base = "../pathway_subgrn_payloads"
  )
  html <- paste(readLines(out_file, warn = FALSE), collapse = "\n")

  expect_match(html, "Sub GRN", fixed = TRUE)
  expect_match(html, "tfScopeSelect", fixed = TRUE)
  expect_match(html, "subgrnTopicThetaCutoff", fixed = TRUE)
  expect_match(html, "subgrnPrimaryTopicOnly", fixed = TRUE)
  expect_match(html, "subgrnTopicTheta", fixed = TRUE)
  expect_match(html, "networkModeSelect", fixed = TRUE)
  expect_match(html, "subgrnSpacingRange", fixed = TRUE)
  expect_match(html, "subgrnPaletteSelect", fixed = TRUE)
  expect_match(html, "subgrnShowArrows", fixed = TRUE)
  expect_match(html, "subgrnResetButton", fixed = TRUE)
  expect_match(html, "subgrnViewLayer", fixed = TRUE)
  expect_match(html, "<option value=\"clustered\">Clustered</option>", fixed = TRUE)
  expect_match(html, "<option value=\"spiral\">Spiral</option>", fixed = TRUE)
  expect_match(html, "openSubgrn", fixed = TRUE)
  expect_match(html, "pathway_subgrn_payloads/chunk_001.js", fixed = TRUE)
})

test_that("Module 3 pathway sub-GRN payload keeps theta and primary topic metadata", {
  pathways <- data.table::data.table(
    comparison_id = "CmpA",
    direction_group = "Up",
    topic = 2L,
    pathway = "Path A",
    pathway_key = "path:a",
    overlap_genes = "GeneA;GeneB",
    padj = 0.01
  )
  edges_docs <- data.table::data.table(
    comparison_id = "CmpA",
    direction = "Target-Up",
    tf = c("TF1", "TF2", "TF3"),
    gene_key = c("GeneA", "GeneA", "GeneB"),
    peak_id = c("Peak1", "Peak2", "Peak3"),
    delta_fp = c(2, 3, 4)
  )
  tf_membership <- data.table::data.table(
    comparison_id = "CmpA",
    direction = "Target-Up",
    tf = c("TF1", "TF1", "TF2", "TF2", "TF3"),
    topic_num = c(2L, 3L, 2L, 4L, 2L),
    theta = c(0.72, 0.20, 0.46, 0.51, 0.29),
    membership_pass = c(TRUE, FALSE, TRUE, TRUE, FALSE),
    primary_topic_num = c(2L, 2L, 4L, 4L, 2L)
  )

  payload <- craftgrn:::.m3tb_build_pathway_subgrn_compact_payload(
    pathways = pathways,
    edges_docs = edges_docs,
    tf_membership = tf_membership
  )

  tf_gene <- payload$tf_gene_edges
  expect_true(all(c("tf_topic_nums", "tf_topic_scores", "tf_primary_topic_num") %in% names(tf_gene)))
  expect_equal(tf_gene[tf == "TF1", tf_topic_nums], "2")
  expect_match(tf_gene[tf == "TF1", tf_topic_scores], "2:0.72", fixed = TRUE)
  expect_match(tf_gene[tf == "TF1", tf_topic_scores], "3:0.2", fixed = TRUE)
  expect_equal(tf_gene[tf == "TF2", tf_topic_nums], "2;4")
  expect_equal(tf_gene[tf == "TF2", tf_primary_topic_num], 4L)
  expect_equal(payload$manifest$n_topic_tfs[[1L]], 2L)
})

test_that("condition pathway sub-GRN reader enforces bounds and reuses cache", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")
  root <- tempfile("condition-subgrn-cache-")
  link_dir <- file.path(root, "condition_links")
  cache_dir <- file.path(root, "cache")
  dir.create(link_dir, recursive = TRUE)
  links <- data.frame(
    tf = c("TF1", "TF2", "TF3", "TF4", "TF1"),
    condition_id = "CondA",
    condition_label = "Condition A",
    gene_key = c("GeneA", "GeneA", "GeneA", "GeneA", "GeneB"),
    peak_id = paste0("Peak", seq_len(5L)),
    fp_score_condition = c(5, 4, 3, 2, 1)
  )
  link_path <- file.path(link_dir, "CondA_condition_links.parquet")
  arrow::write_parquet(links, link_path)
  data.table::fwrite(
    data.table::data.table(condition_id = "CondA", path = link_path),
    file.path(link_dir, "condition_links_manifest.csv")
  )
  pathways <- data.table::data.table(
    comparison_id = "CondA",
    topic_num = 1L,
    overlap_genes = "GeneA;GeneB"
  )
  membership <- data.table::data.table(
    comparison_id = "CondA",
    tf = paste0("TF", 1:4),
    topic_num = 1L,
    membership_pass = TRUE
  )
  withr::local_envvar(CRAFTGRN_PATHWAY_SUBGRN_MAX_TFS_PER_GENE = "2")
  first <- craftgrn:::.m3tb_read_condition_links_for_subgrn(
    link_dir,
    pathways,
    tf_membership = membership,
    cache_dir = cache_dir
  )
  second <- craftgrn:::.m3tb_read_condition_links_for_subgrn(
    link_dir,
    pathways,
    tf_membership = membership,
    cache_dir = cache_dir
  )
  expect_lte(first[gene_key == "GeneA", data.table::uniqueN(tf)], 2L)
  expect_equal(
    first[order(gene_key, tf), .(tf, gene_key, peak_id)],
    second[order(gene_key, tf), .(tf, gene_key, peak_id)]
  )
  cache_files <- list.files(cache_dir, pattern = "[.]rds$", full.names = TRUE)
  expect_length(cache_files, 1L)
  expect_identical(readRDS(cache_files[[1L]])$cache_version, 2L)
})

test_that("pathway sub-GRN chunks use compressed columnar browser payloads", {
  obj <- list(
    manifest = data.table::data.table(subgrn_context_id = "ctx1"),
    tf_gene_edges = data.table::data.table(tf = "TF1", gene_key = "GeneA"),
    tf_peak_gene_triplets = data.table::data.table(
      tf = "TF1", peak_id = "Peak1", gene_key = "GeneA"
    )
  )
  packed <- lapply(
    obj,
    craftgrn:::.module2_report_browser_browser_payload_to_columnar
  )
  encoded <- craftgrn:::.module2_report_browser_encode_browser_json_deflate_base64(
    packed
  )
  expect_gt(nchar(encoded), 20L)
  js <- craftgrn:::.m3tb_subgrn_js()
  expect_true(any(grepl("subgrnDecodePayload", js, fixed = TRUE)))
  expect_true(any(grepl("compressed_columnar", js, fixed = TRUE)))
  expect_true(any(grepl(
    "pathKey=x=>String(x||'').trim().toLowerCase().replace(/[^a-z0-9]+/g,'')",
    js,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "pathKey(m.pathway_key||m.pathway||'')===pk",
    js,
    fixed = TRUE
  )))
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

test_that("Fibroblast paper benchmark loader keeps Gold and Silver fibroblast targets", {
  path <- tempfile(fileext = ".csv")
  data.table::fwrite(data.table::data.table(
    id = paste0("row", 1:8),
    row_type = c("module_member", "module_member", rep("regulatory_edge", 6)),
    regulator = c("Batf_Irf4_module", "Batf_Irf4_module", "Batf_Irf4_module", "Batf_Irf4_module", "Tcell_module", "Batf_Irf4_module", "Batf_Irf4_module", "Batf_Irf4_module"),
    target = c("Batf", "Irf4", "Ccr7", "Ifngr1", "Sell", "LooseGene", "DownGene", "PathwayA"),
    target_type = c("TF", "TF", "gene", "gene", "gene", "gene", "gene", "pathway"),
    effect = c("member", "member", "up", "up", "up", "up", "down", "up"),
    sign_num = c(0, 0, 1, 1, 1, 1, -1, 1),
    context = c(
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "P14_Tcell_Batf_cKO",
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "NIH3T3_fibroblast_TF_overexpression_72h",
      "NIH3T3_fibroblast_TF_overexpression_72h"
    ),
    observation = "",
    evidence = "",
    confidence = c("gold", "gold", "gold", "silver", "gold", "bronze", "gold", "gold"),
    source = "test"
  ), path)

  bench <- craftgrn:::.m3fb_load_tsao_fibroblast_benchmark(path, confidence = c("gold", "silver"))
  bench_all <- craftgrn:::.m3fb_load_tsao_fibroblast_benchmark(path)

  expect_setequal(bench$module_members$target, c("Batf", "Irf4"))
  expect_setequal(bench$targets$target, c("Ccr7", "Ifngr1"))
  expect_false("LooseGene" %in% bench$targets$target)
  expect_false("DownGene" %in% bench$targets$target)
  expect_true(all(bench$targets$confidence %in% c("gold", "silver")))
  expect_true("LooseGene" %in% bench_all$targets$target)
})

test_that("Fibroblast paper benchmark scorer favors TF-target co-topic recovery", {
  bench <- list(
    module_members = data.table::data.table(
      regulator = "Batf_Irf4_Runx3_Tbx21_module",
      target = c("Batf", "Irf4", "Runx3", "Tbx21")
    ),
    targets = data.table::data.table(
      regulator = "Batf_Irf4_Runx3_Tbx21_module",
      target = c("Ccr7", "Ifngr1", "Il17ra", "Tgfb1"),
      confidence = c("gold", "gold", "silver", "silver")
    )
  )
  theta_good <- matrix(
    c(
      0.82, 0.10, 0.08,
      0.78, 0.12, 0.10,
      0.69, 0.21, 0.10,
      0.72, 0.18, 0.10
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      paste0("Fib_BIRT_vs_NoTF::Target-Up::", c("Batf", "Irf4", "Runx3", "Tbx21")),
      paste0("Topic", 1:3)
    )
  )
  theta_bad <- theta_good
  theta_bad[, ] <- c(
    0.82, 0.10, 0.08,
    0.10, 0.78, 0.12,
    0.10, 0.18, 0.72,
    0.72, 0.18, 0.10
  )
  phi <- matrix(
    c(
      0.90, 0.82, 0.80, 0.72,
      0.08, 0.10, 0.12, 0.18,
      0.02, 0.08, 0.08, 0.10
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Topic", 1:3), c("GENE:Ccr7", "GENE:Ifngr1", "GENE:Il17ra", "PEAK:Tgfb1"))
  )
  rows <- data.table::data.table(
    method_order = 1:2,
    method = c("comparison_aggr_multivi", "comparison_aggr_lda"),
    method_setup = c("diff fp aggr | MultiVI", "diff fp aggr | LDA"),
    model_label = c("MultiVI", "LDA"),
    selected_k = c(3L, 3L)
  )

  good <- craftgrn:::.m3fb_score_model_row(
    theta = theta_good,
    phi = phi,
    row = rows[1],
    benchmark = bench,
    comparisons = "Fib_BIRT_vs_NoTF",
    theta_cutoff = 0.3,
    gene_score_cutoff = 0.7
  )
  bad <- craftgrn:::.m3fb_score_model_row(
    theta = theta_bad,
    phi = phi,
    row = rows[2],
    benchmark = bench,
    comparisons = "Fib_BIRT_vs_NoTF",
    theta_cutoff = 0.3,
    gene_score_cutoff = 0.7
  )

  expect_gt(good$comparison_scores$model_score[[1L]], bad$comparison_scores$model_score[[1L]])
  expect_equal(good$selected_topics$topic_num, 1L)
  expect_equal(good$comparison_scores$n_required_tfs_covered[[1L]], 4L)
  expect_equal(good$comparison_scores$n_targets_cotopic[[1L]], 4L)
})

test_that("Fibroblast paper benchmark writer emits review tables and index", {
  root <- tempfile("fib-paper-review-")
  model_dir <- file.path(root, "run_004_comparison_aggr_multivi", "topic_models")
  dir.create(file.path(model_dir, "vae_models"), recursive = TRUE)
  review_dir <- file.path(root, "review")
  bench_path <- file.path(root, "tsao.csv")
  data.table::fwrite(data.table::data.table(
    id = paste0("row", 1:8),
    row_type = c(rep("module_member", 4), rep("regulatory_edge", 4)),
    regulator = "Batf_Irf4_Runx3_Tbx21_module",
    target = c("Batf", "Irf4", "Runx3", "Tbx21", "Ccr7", "Ifngr1", "Il17ra", "Tgfb1"),
    target_type = c(rep("TF", 4), rep("gene", 4)),
    effect = c(rep("member", 4), rep("up", 4)),
    sign_num = c(rep(0, 4), rep(1, 4)),
    context = "NIH3T3_fibroblast_TF_overexpression_72h",
    observation = "",
    evidence = "",
    confidence = c(rep("gold", 6), rep("silver", 2)),
    source = "test"
  ), bench_path)
  theta <- data.table::data.table(
    doc_id = paste0("Fib_BIRT_vs_NoTF::Target-Up::", c("Batf", "Irf4", "Runx3", "Tbx21")),
    Topic1 = c(0.8, 0.75, 0.7, 0.72),
    Topic2 = c(0.1, 0.15, 0.2, 0.18),
    Topic3 = c(0.1, 0.1, 0.1, 0.1)
  )
  phi <- data.table::data.table(
    term_id = paste0("Topic", 1:3),
    `GENE:Ccr7` = c(0.90, 0.08, 0.02),
    `GENE:Ifngr1` = c(0.82, 0.10, 0.08),
    `GENE:Il17ra` = c(0.80, 0.12, 0.08),
    `PEAK:Tgfb1` = c(0.72, 0.18, 0.10)
  )
  data.table::fwrite(theta, file.path(model_dir, "vae_models", "theta_K3.csv"))
  data.table::fwrite(phi, file.path(model_dir, "vae_models", "phi_K3.csv"))
  rows <- data.table::data.table(
    method_order = 1L,
    method = "comparison_aggr_multivi",
    method_setup = "diff fp aggr | MultiVI",
    model_label = "MultiVI",
    selected_k = 3L,
    model_dir = model_dir
  )

  out <- craftgrn:::.m3fb_score_existing_models(
    output_dir = root,
    model_rows = rows,
    benchmark_csv = bench_path,
    review_dir = review_dir,
    comparisons = "Fib_BIRT_vs_NoTF",
    verbose = FALSE
  )

  expect_true(file.exists(file.path(out$output_dir, "index.html")))
  expect_true(file.exists(file.path(out$output_dir, "tables", "fibroblast_model_k_leaderboard.csv")))
  expect_true(file.exists(file.path(out$output_dir, "figures", "fibroblast_model_k_leaderboard.pdf")))
  expect_equal(out$leaderboard$rank[[1L]], 1L)
  expect_equal(out$comparison_scores$n_targets_cotopic[[1L]], 4L)
})

test_that("Fibroblast paper benchmark maps targets to the expected comparison direction", {
  bench <- list(
    module_members = data.table::data.table(),
    targets = data.table::data.table(
      regulator = c("Batf_Irf4_module", "Batf_Irf4_module", "Batf_Irf4_Runx3_Tbx21_module", "Batf_Irf4_Runx3_Tbx21_module"),
      target = c("Ccr7", "Ifngr1", "Nfatc1", "Tgfb1"),
      target_key = c("CCR7", "IFNGR1", "NFATC1", "TGFB1"),
      confidence = c("silver", "silver", "gold", "bronze"),
      effect = "positive",
      sign_num = 1
    )
  )

  bi_targets <- craftgrn:::.m3fb_targets_for_comparison(bench, "Fib_BI_vs_NoTF")
  birt_targets <- craftgrn:::.m3fb_targets_for_comparison(bench, "Fib_BIRT_vs_NoTF")
  addon_targets <- craftgrn:::.m3fb_targets_for_comparison(bench, "Fib_BIRT_vs_BI")

  expect_setequal(bi_targets$target, c("Ccr7", "Ifngr1"))
  expect_setequal(birt_targets$target, c("Nfatc1", "Tgfb1"))
  expect_setequal(addon_targets$target, c("Nfatc1", "Tgfb1"))
  expect_equal(unique(craftgrn:::.m3fb_default_comparison_map()[comparison_id == "Fib_BIRT_vs_BI", direction]), "Target-Up")
})

test_that("Fibroblast paper benchmark full membership keeps Bronze and nonpassing genes for plots", {
  bench <- list(
    module_members = data.table::data.table(),
    targets = data.table::data.table(
      regulator = "Batf_Irf4_Runx3_Tbx21_module",
      target = c("Ccr7", "Akt2"),
      target_key = c("CCR7", "AKT2"),
      confidence = c("gold", "bronze"),
      effect = "positive",
      sign_num = 1
    )
  )
  theta <- matrix(
    c(0.8, 0.2, 0.7, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Fib_BIRT_vs_NoTF::Batf::Target-Up", "Fib_BIRT_vs_NoTF::Irf4::Target-Up"),
      c("Topic1", "Topic2")
    )
  )
  phi <- matrix(
    c(0.9, 0.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:Ccr7"))
  )
  rows <- data.table::data.table(
    method_order = 1L,
    method = "comparison_aggr_lda",
    method_setup = "diff fp aggr | LDA",
    model_label = "LDA",
    selected_k = 2L
  )

  full <- craftgrn:::.m3fb_full_benchmark_membership(
    theta = theta,
    phi = phi,
    row = rows,
    benchmark = bench,
    comparisons = "Fib_BIRT_vs_NoTF"
  )

  expect_true("Akt2" %in% full$item)
  expect_true("bronze" %in% full$confidence)
  expect_false(full[item == "Akt2" & topic_num == 1L, pass][[1L]])
  expect_equal(unique(full[item == "Akt2", expected_direction]), "Target-Up")
})

test_that("Fibroblast paper benchmark prepares top model membership and clean heatmap labels", {
  bench <- list(
    module_members = data.table::data.table(),
    targets = data.table::data.table(
      regulator = "Batf_Irf4_Runx3_Tbx21_module",
      target = c("Ccr7", "Akt2"),
      target_key = c("CCR7", "AKT2"),
      confidence = c("gold", "bronze"),
      effect = "positive",
      sign_num = 1
    )
  )
  theta <- matrix(
    c(0.8, 0.2, 0.7, 0.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Fib_BIRT_vs_NoTF::Batf::Target-Up", "Fib_BIRT_vs_NoTF::Irf4::Target-Up"),
      c("Topic1", "Topic2")
    )
  )
  phi <- matrix(
    c(0.9, 0.1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:Ccr7"))
  )
  row1 <- data.table::data.table(
    method_order = 1L,
    method = "comparison_aggr_lda",
    method_setup = "diff fp aggr | LDA",
    model_label = "LDA",
    selected_k = 2L
  )
  row2 <- data.table::copy(row1)
  row2[, `:=`(method_order = 2L, method_setup = "diff fp aggr | MultiVI", model_label = "MultiVI")]
  leaderboard <- data.table::data.table(
    method_order = 1:2,
    method = c("comparison_aggr_lda", "comparison_aggr_multivi"),
    method_setup = c("diff fp aggr | LDA", "diff fp aggr | MultiVI"),
    model_label = c("LDA", "MultiVI"),
    k = c(2L, 2L),
    rank = c(1L, 2L)
  )
  model_rows <- data.table::data.table(
    method_order = 1:2,
    method = c("comparison_aggr_lda", "comparison_aggr_multivi"),
    method_setup = c("diff fp aggr | LDA", "diff fp aggr | MultiVI"),
    model_label = c("LDA", "MultiVI"),
    selected_k = c(2L, 2L),
    model_dir = tempfile(c("m1", "m2"))
  )
  for (md in model_rows$model_dir) {
    dir.create(file.path(md, "vae_models"), recursive = TRUE)
    data.table::fwrite(data.table::data.table(doc_id = rownames(theta), theta), file.path(md, "vae_models", "theta_K2.csv"))
    data.table::fwrite(data.table::data.table(term_id = rownames(phi), phi), file.path(md, "vae_models", "phi_K2.csv"))
  }

  top_membership <- craftgrn:::.m3fb_full_membership_for_top_models(
    model_rows = model_rows,
    leaderboard = leaderboard,
    benchmark = bench,
    comparison_map = craftgrn:::.m3fb_default_comparison_map()[comparison_id == "Fib_BIRT_vs_NoTF"],
    top_n_models = 2L
  )
  plot_data <- craftgrn:::.m3fb_prepare_story_heatmap_data(
    leaderboard = leaderboard,
    membership_long = top_membership,
    top_n_models = 2L
  )

  expect_setequal(unique(top_membership$model_rank), c(1L, 2L))
  expect_false(any(grepl("\\[", plot_data$item_label)))
  expect_setequal(
    unique(plot_data[item_type %in% c("Required TF", "Target gene"), .(item_type, item_color)]$item_color),
    c("#2563eb", "#991b1b")
  )
})

test_that("Fibroblast paper SubGRN index exposes model topic and separate score controls", {
  out_dir <- tempfile("paper-subgrn-")
  dir.create(out_dir, recursive = TRUE)
  manifest <- data.table::data.table(
    subgrn_context_id = "paper_001",
    model_rank = 1L,
    method_setup = "diff fp aggr | LDA",
    model_label = "LDA",
    k = 12L,
    comparison_id = "Fib_BIRT_vs_NoTF",
    comparison_label = "Fib_BIRT_vs_NoTF::Target-Up",
    direction = "Target-Up",
    topic_num = 3L,
    pathway = "Paper benchmark targets",
    pathway_key = "paper:rank1:Fib_BIRT_vs_NoTF:Target-Up:3",
    overlap_genes = "Ccr7;Ifngr1",
    n_overlap_genes = 2L,
    n_tf_gene_edges = 1L,
    n_tf_peak_gene_triplets = 1L,
    payload_file = "paper_topic_subgrn_payload.js",
    selected_story_topic = TRUE
  )
  payload <- list(
    manifest = manifest,
    tf_gene_edges = data.table::data.table(
      subgrn_context_id = "paper_001",
      tf = "Batf",
      tf_upper = "BATF",
      gene_key = "Ccr7",
      abs_edge_score = 1,
      n_supporting_peaks = 1L,
      best_peak_id = "peak1",
      tf_topic_score = 0.8,
      gene_topic_score = 0.9
    ),
    tf_peak_gene_triplets = data.table::data.table(
      subgrn_context_id = "paper_001",
      tf = "Batf",
      tf_upper = "BATF",
      peak_id = "peak1",
      gene_key = "Ccr7",
      abs_edge_score = 1,
      tf_topic_score = 0.8,
      gene_topic_score = 0.9
    )
  )
  writeLines(
    paste0(
      "window.CRAFTGRN_PAPER_PAYLOAD=",
      jsonlite::toJSON(payload, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null"),
      ";"
    ),
    file.path(out_dir, "paper_topic_subgrn_payload.js"),
    useBytes = TRUE
  )

  index <- craftgrn:::.m3fb_write_subgrn_index(
    out_dir = out_dir,
    manifest = manifest,
    payload_file = "paper_topic_subgrn_payload.js",
    top_models = manifest[, .(model_rank, method_setup, model_label, k)]
  )
  html <- paste(readLines(index, warn = FALSE), collapse = "\n")

  expect_match(html, "paperModelSelect", fixed = TRUE)
  expect_match(html, "paperTopicSelect", fixed = TRUE)
  expect_match(html, "paperTfCutoff", fixed = TRUE)
  expect_match(html, "paperGeneCutoff", fixed = TRUE)
  expect_match(html, "paperEvidenceMode", fixed = TRUE)
  expect_match(html, "TF-peak-gene", fixed = TRUE)
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
  old_scan <- getOption("craftgrn.topic_review.scan_link_tables")
  options(
    craftgrn.topic_review.fast_summary = TRUE,
    craftgrn.topic_review.scan_link_tables = TRUE
  )
  on.exit(options(
    craftgrn.topic_review.fast_summary = old,
    craftgrn.topic_review.scan_link_tables = old_scan
  ), add = TRUE)
  out <- craftgrn:::.m3tb_summarize_topic_links(root, rows, review_dir = review)

  expect_true(any(out$shared$unit == "Links" & out$shared$selected_k == 2L))
  expect_true(any(out$shared$unit == "Links" & out$shared$selected_k == 3L))
  expect_true(any(out$shared$unit == "Genes" & out$shared$selected_k == 3L))
  expect_true(any(out$shared$unit == "TFs" & out$shared$selected_k == 3L))
  expect_true(any(out$link_assignment$count_basis == "Topic-link rows" & out$link_assignment$selected_k == 3L))
})

test_that("Module 3 review blocks oversized raw topic-link scans", {
  root <- tempfile("module3-link-scan-guard-")
  extraction_dir <- file.path(root, "topic_extraction", "K2")
  dir.create(extraction_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(
      tf = "TF1",
      peak_id = "P1",
      gene_key = "G1",
      topic_num = 1L,
      link_pass = TRUE
    ),
    file.path(extraction_dir, "topic_links_pass.csv")
  )
  rows <- data.table::data.table(
    method = "comparison_aggr_lda",
    method_order = 1L,
    method_setup = "diff fp aggr | LDA",
    setup = "std_tf_diff_fp_aggr",
    model_label = "LDA",
    selected_k = 2L,
    topic_extraction_dir = file.path(root, "topic_extraction")
  )
  withr::local_options(list(
    craftgrn.topic_review.scan_link_tables = TRUE,
    craftgrn.topic_review.max_link_file_bytes = 1
  ))

  expect_error(
    craftgrn:::.m3tb_summarize_topic_links(root, rows),
    "exceeds the safe review scan limit"
  )
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

test_that("Module 3 review HTML keeps condition reports in method subfolders", {
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
  expect_false(file.exists(file.path(review_dir, "topic_method_k_topic_mds_report.html")))
  expect_true(file.exists(file.path(review_dir, "topic_method_k_condition_mds_report.html")))
  expect_false(file.exists(file.path(review_dir, "topic_report_K2.html")))
  expect_false(file.exists(file.path(review_dir, "condition_topic_report_K2.html")))

  condition_pages <- list.files(file.path(review_dir, "assets", "p"), pattern = "[.]html$", full.names = TRUE)
  expect_false(dir.exists(file.path(review_dir, "topic_reports")))
  expect_false(dir.exists(file.path(review_dir, "condition_topic_reports")))
  expect_equal(length(condition_pages), 2L)
  expect_equal(length(unique(basename(condition_pages))), 2L)

  condition_index <- paste(readLines(file.path(review_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  expect_match(condition_index, "assets/p/", fixed = TRUE)
  expect_match(condition_index, "Method <select", fixed = TRUE)
  expect_match(condition_index, "K <select", fixed = TRUE)
})

test_that("Module 3 review validation rejects outdated HTML", {
  current <- tempfile(fileext = ".html")
  writeLines(c(
    '<meta name="craftgrn-module3-report-schema" content="2"/>',
    '<h2>Topic Bar Plots</h2>',
    '<input id="tfSearchInput"/>',
    '<select id="tfSelect" multiple size="3"></select>',
    '<script>const SUBGRN_MANIFEST=[];</script>'
  ), current)
  expect_true(craftgrn:::.m3tb_validate_current_report_html(current))

  stale <- tempfile(fileext = ".html")
  writeLines("<h2>Topic Waterfall</h2>", stale)
  expect_error(
    craftgrn:::.m3tb_validate_current_report_html(stale),
    "does not match the current schema",
    fixed = TRUE
  )
})

test_that("Module 3 condition-pair report uses schema 10 pathway and GRN modes", {
  report_root <- tempfile("module3-compact-report-")
  out_file <- file.path(report_root, "review", "assets", "p", "report.html")
  craftgrn:::.m3tb_condition_report_html(
    title = "Condition pair",
    group_mds = data.table::data.table(
      comparison_label = c("A", "B"),
      display_label = c("A", "B"),
      MDS1 = c(-1, 1),
      MDS2 = c(0, 0),
      n_docs = 2L
    ),
    group_topic = data.table::data.table(
      comparison_label = rep(c("A", "B"), each = 2L),
      topic = rep(c("Topic1", "Topic2"), 2L),
      topic_num = rep(1:2, 2L),
      theta_mean = c(0.8, 0.2, 0.3, 0.7)
    ),
    tf_topic = data.table::data.table(
      comparison_label = rep(c("A", "B"), each = 2L),
      tf = "Batf",
      tf_upper = "BATF",
      topic_num = rep(1:2, 2L),
      theta = c(0.8, 0.2, 0.3, 0.7)
    ),
    pathways = data.table::data.table(
      comparison_label = c("A", "B"),
      topic_num = 1L,
      pathway = "Path A",
      pathway_key = "path:a",
      padj = c(0.01, 0.02),
      combined_score = c(12, 4),
      gene_in = c(3L, 2L),
      overlap_genes = c("G1;G2;G3", "G1;G2")
    ),
    out_html = out_file,
    condition_payload = list(
      payload_file = "pair.js",
      payload_base = "../cp",
      conditions = c("A", "B"),
      n_tf_gene = 2L,
      n_tf_peak_gene = 2L
    ),
    report_state = list(
      condition_colors = list(A = "#E15759", B = "#4E79A7"),
      condition_order = c("B", "A"),
      defaults = list(condition_1 = "B", condition_2 = "A", topic = 2L)
    )
  )
  html <- paste(readLines(out_file, warn = FALSE), collapse = "\n")
  expect_match(html, '"base":"../cr"', fixed = TRUE)
  expect_match(html, '"payload_base":"../cp"', fixed = TRUE)
  expect_length(
    list.files(file.path(report_root, "review", "assets", "cr"), pattern = "[.]js$"),
    1L
  )
  expect_match(html, 'craftgrn-module3-report-schema" content="10', fixed = TRUE)
  expect_match(html, "Topic Activity", fixed = TRUE)
  expect_match(html, 'id="conditionTopicMetric"', fixed = TRUE)
  expect_match(html, 'value="rna_delta">Differential RNA activity', fixed = TRUE)
  expect_match(html, "pairwiseRnaTopicProfiles", fixed = TRUE)
  expect_match(html, "positive pairwise RNA activity", fixed = TRUE)
  expect_match(html, 'id="thetaAggregation"', fixed = TRUE)
  expect_match(html, "Matched TFs", fixed = TRUE)
  expect_match(html, "matchedTfKeys", fixed = TRUE)
  expect_match(html, "r.gexprA>=cfg.exprMin", fixed = TRUE)
  expect_no_match(html, "RNA expression share", fixed = TRUE)
  expect_match(html, "TF Probability", fixed = TRUE)
  expect_match(html, "id=\"tfButterflySvg\"", fixed = TRUE)
  expect_match(html, "P(topic | condition::TF)", fixed = TRUE)
  expect_match(html, "TF Activity", fixed = TRUE)
  expect_match(html, 'id="shortConditionNames" type="checkbox" checked', fixed = TRUE)
  expect_match(html, "mdsShortLabelMap", fixed = TRUE)
  expect_match(html, "mdsLabelCandidates", fixed = TRUE)
  expect_match(html, "mdsFinalLabelPenalty", fixed = TRUE)
  expect_match(html, "mdsRectOverlap", fixed = TRUE)
  expect_match(html, "mdsClampLabelCandidate", fixed = TRUE)
  expect_match(html, "mdsSegmentsCross", fixed = TRUE)
  expect_match(html, "mdsLeaderConflictPenalty", fixed = TRUE)
  expect_match(html, "mdsPruneLeaderLines", fixed = TRUE)
  expect_match(html, "fill-opacity", fixed = TRUE)
  expect_match(html, "selected?30:24", fixed = TRUE)
  expect_match(html, "pointInset=38", fixed = TRUE)
  expect_match(html, "pointInset=22", fixed = TRUE)
  expect_match(html, "fullY=[Math.max(0,yLoRaw-yPad)", fixed = TRUE)
  expect_match(html, "value=yd[0]+frac*(yd[1]-yd[0])", fixed = TRUE)
  expect_match(html, "Math.log2((shareA+pseudo)/(shareB+pseudo))", fixed = TRUE)
  expect_match(html, "log2 relative activity", fixed = TRUE)
  expect_false(grepl(
    "log2 relative FP activity (Condition 1 / Condition 2)",
    html,
    fixed = TRUE
  ))
  expect_match(html, "Math.asinh(d.deltaShare/soft)", fixed = TRUE)
  expect_match(html, "Normalized FP difference (pp):", fixed = TRUE)
  expect_match(html, "activityViewDomain", fixed = TRUE)
  expect_false(grepl("pseudo-log", html, fixed = TRUE))
  expect_match(html, "id=\"conditionProbabilityTitle\"", fixed = TRUE)
  expect_match(html, "id=\"tfConditionProbabilityTitle\"", fixed = TRUE)
  expect_match(html, "updateProbabilityPanelTitles", fixed = TRUE)
  expect_match(html, "Failed to load condition report payload after retries", fixed = TRUE)
  expect_match(html, "PAYLOAD_PROMISES.delete(file)", fixed = TRUE)
  expect_match(html, "function indexConditionPart", fixed = TRUE)
  expect_match(html, "EDGE_BY_COND.set(condition", fixed = TRUE)
  expect_match(html, "function ensureSelectedConditionEdges(){return loadSelectedConditionParts('edge')}", fixed = TRUE)
  expect_match(html, "const drawMdsUncached=drawMds", fixed = TRUE)
  expect_match(html, "const drawActivityUncached=drawActivity", fixed = TRUE)
  expect_match(html, "if(n&&!isFile)url.searchParams.set", fixed = TRUE)
  expect_match(html, "location.protocol!=='file:'", fixed = TRUE)
  expect_match(html, "conditions.reduce", fixed = TRUE)
  expect_match(html, "label:conditionLabel(o.value)", fixed = TRUE)
  expect_match(html, "shortConditionNames').addEventListener('change',refresh", fixed = TRUE)
  expect_match(html, "const pos=rows.filter(d=>d.delta>=0)", fixed = TRUE)
  expect_match(html, "neg=rows.filter(d=>d.delta<0)", fixed = TRUE)
  expect_match(html, "if(i%2===1)g.appendChild", fixed = TRUE)
  expect_match(
    html,
    "if(['cond1','cond2','condition','tf'].includes(changed))refreshConditionData()",
    fixed = TRUE
  )
  expect_match(html, "ranked-expression enrichment z-score", fixed = TRUE)
  expect_match(html, "differentialConditionLabels", fixed = TRUE)
  expect_match(html, "rowH=Math.min(40", fixed = TRUE)
  expect_match(html, "rowH=Math.min(38", fixed = TRUE)
  expect_match(html, "rowH=Math.min(22", fixed = TRUE)
  expect_match(html, "x:ax.cx+gap", fixed = TRUE)
  expect_match(html, "x:c2?ax.cx-gap-w1:ax.cx", fixed = TRUE)
  expect_match(html, "pathLabelTopicSpecific", fixed = TRUE)
  expect_match(html, "pathLabelConditionSpecific", fixed = TRUE)
  expect_match(html, ".pathLabelBothSpecific{fill:#171717}", fixed = TRUE)
  expect_match(html, "labelContrast", fixed = TRUE)
  expect_match(html, "id=\"cond1Select\"", fixed = TRUE)
  expect_match(html, "id=\"cond1Color\"", fixed = TRUE)
  expect_match(html, "id=\"cond2Color\"", fixed = TRUE)
  expect_match(html, "conditionColor(1)", fixed = TRUE)
  expect_match(html, "conditionColor(2)", fixed = TRUE)
  expect_match(html, "fill:conditionColor(1)", fixed = TRUE)
  expect_match(html, "fill:conditionColor(2)", fixed = TRUE)
  expect_match(html, "t.getComputedTextLength()<=maxWidth", fixed = TRUE)
  expect_false(grepl("syncConditionPlotColors", html, fixed = TRUE))
  expect_match(html, "q.textContent=v.toFixed(1)", fixed = TRUE)
  expect_match(html, "x1:cx-gap-xw,y1:axisY", fixed = TRUE)
  expect_match(html, "if(!c2)g.appendChild(el('line'", fixed = TRUE)
  expect_match(html, "if(selected)return conditionTopicRows()", fixed = TRUE)
  expect_match(html, "topicMode=!!selectedTf", fixed = TRUE)
  expect_match(html, "tfButterflyRowSelected", fixed = TRUE)
  expect_match(html, "if(changed==='pathway'&&pathwaySelect.value)openNetwork()", fixed = TRUE)
  expect_match(html, "countValue=tfSelect.value?(paired?r.tfCountUnion:r.tfCountA):0", fixed = TRUE)
  expect_match(html, "class:'pathTfCount'", fixed = TRUE)
  expect_no_match(html, "tfCountA+'/'+r.tfCountB+'/'+r.tfCountUnion", fixed = TRUE)
  expect_match(html, "value=\"topic\" selected>In Topic", fixed = TRUE)
  expect_match(html, "selectTf(d.tfu)", fixed = TRUE)
  expect_match(html, "return loadOverviewPayload()", fixed = TRUE)
  expect_match(html, "loadNetworkPayload()", fixed = TRUE)
  expect_match(html, "Loading compressed network data", fixed = TRUE)
  expect_match(html, "condition_payload_files", fixed = TRUE)
  expect_match(html, "loadReportData().then", fixed = TRUE)
  expect_match(html, "indexTfTopics()", fixed = TRUE)
  expect_match(html, "k==='tf_topic'?packed[k]:columnarRows", fixed = TRUE)
  expect_match(html, "id=\"networkTfPalette\"", fixed = TRUE)
  expect_match(html, 'id="networkTopTf" type="number" min="1" value="100"', fixed = TRUE)
  expect_match(html, 'id="topicPageStatus"', fixed = TRUE)
  expect_match(html, 'id="tfPageStatus"', fixed = TRUE)
  expect_match(html, 'id="pathPageStatus"', fixed = TRUE)
  expect_match(html, "grid-template-columns:minmax(620px,1fr) minmax(620px,1fr)", fixed = TRUE)
  expect_match(html, "const W=760,H=760,L=82,R=132", fixed = TRUE)
  expect_match(html, "box.width-t.offsetWidth-8", fixed = TRUE)
  expect_match(html, "const PAGE_SIZE={topic:20,tf:20,path:35}", fixed = TRUE)
  expect_match(html, 'id="pathwayScoreMethod"', fixed = TRUE)
  expect_match(html, 'id="pathwayDeOnly" type="checkbox"', fixed = TRUE)
  expect_match(html, "const displayRows=x=>deOnly?x.filter(r=>r.deOverlap>0):x", fixed = TRUE)
  expect_match(html, "Dot size = expressed overlap genes", fixed = TRUE)
  expect_match(html, "Dot size = topic-overlap genes", fixed = TRUE)
  expect_match(html, "Overall topic enrichment combined score", fixed = TRUE)
  expect_match(html, "Overall topic pathways", fixed = TRUE)
  expect_match(html, "if(!c1)return PATHWAYS.filter", fixed = TRUE)
  expect_no_match(html, ".slice(0,PAGE_SIZE.path)", fixed = TRUE)
  expect_no_match(html, ".filter(g=>!deOnly||", fixed = TRUE)
  expect_match(html, 'id="networkCorrectDirectionOnly" type="checkbox"', fixed = TRUE)
  expect_match(html, "function dynamicPathwayScore", fixed = TRUE)
  expect_match(html, "function decorateLinkDirection", fixed = TRUE)
  expect_match(html, "correctGenes", fixed = TRUE)
  expect_match(html, 'id="activityZoomIn"', fixed = TRUE)
  expect_match(html, '>Zoom out</button>', fixed = TRUE)
  expect_match(html, '>Zoom in</button>', fixed = TRUE)
  expect_false(grepl('id="activityZoomReset"', html, fixed = TRUE))
  expect_match(html, "strength=.78+.18*clamp", fixed = TRUE)
  expect_match(html, "condition1First=!c2||groupTopic(c1,topic)>=groupTopic(c2,topic)", fixed = TRUE)
  expect_no_match(html, "Overall topic significance", fixed = TRUE)
  expect_match(html, "significant pathways", fixed = TRUE)
  expect_match(html, "pathwayStars", fixed = TRUE)
  expect_match(html, "r:sel?19:14", fixed = TRUE)
  expect_match(html, "selectMdsCondition", fixed = TRUE)
  expect_match(html, "ACTIVE_CONDITION_SIDE", fixed = TRUE)
  expect_match(html, '"condition_colors":{"A":"#E15759","B":"#4E79A7"}', fixed = TRUE)
  expect_match(html, "REPORT_STATE.condition_order", fixed = TRUE)
  expect_match(html, "configuredCondition(defaults,'condition_2'", fixed = TRUE)
  expect_match(html, "defaults[key+'_suffix']", fixed = TRUE)
  expect_false(grepl("rows[i-1].scoreA>=rows[i-1].scoreB", html, fixed = TRUE))
  expect_match(html, "class:'pathAxisTickLabel'", fixed = TRUE)
  expect_match(html, "Math.pow(2,value)-1", fixed = TRUE)
  expect_match(html, "xTitle.textContent='MDS1'", fixed = TRUE)
  expect_match(html, "yTitle.textContent='MDS2'", fixed = TRUE)
  expect_false(grepl("chartScroll", html, fixed = TRUE))
  expect_false(grepl("pathScroll", html, fixed = TRUE))
  expect_false(grepl("pathShowAll", html, fixed = TRUE))
  expect_match(html, 'id="networkTabFilter"', fixed = TRUE)
  expect_match(html, 'id="networkTabLayout"', fixed = TRUE)
  expect_match(html, 'id="networkTabAppearance"', fixed = TRUE)
  expect_match(html, 'id="showPathwaysMode"', fixed = TRUE)
  expect_match(html, 'id="showGrnMode"', fixed = TRUE)
  expect_false(grepl('id="selectedContext"', html, fixed = TRUE))
  expect_false(grepl("contextBadge", html, fixed = TRUE))
  expect_false(grepl("updateContext()", html, fixed = TRUE))
  expect_match(html, 'id="networkPathwayFocus"', fixed = TRUE)
  expect_match(html, '<option value="highlight" selected>Highlight</option>', fixed = TRUE)
  expect_match(html, "function syncReportMode", fixed = TRUE)
  expect_match(html, "else if(networkOpen){updateNetworkHeading();drawNetwork()}", fixed = TRUE)
  expect_match(html, "byId('showPathwaysMode').onclick=closeNetwork", fixed = TRUE)
  expect_match(html, "byId('showGrnMode').onclick", fixed = TRUE)
  expect_match(html, 'id="networkThetaPreset"', fixed = TRUE)
  expect_match(html, 'id="networkPhiPreset"', fixed = TRUE)
  expect_match(html, 'id="networkPrimaryOnly"', fixed = TRUE)
  expect_match(html, 'id="networkSpacing"', fixed = TRUE)
  expect_match(html, '>Reset layout</button>', fixed = TRUE)
  expect_match(html, '>Red-blue (Default)<', fixed = TRUE)
  expect_match(html, '<option value="single">Single color</option>', fixed = TRUE)
  expect_match(html, 'id="networkTfSingleColor"', fixed = TRUE)
  expect_match(html, 'id="networkGeneSingleColor"', fixed = TRUE)
  expect_match(html, 'id="networkEdgeSingleColor"', fixed = TRUE)
  expect_match(html, 'id="networkTfMin" type="number" min="6" max="40"', fixed = TRUE)
  expect_match(html, 'id="networkTfMax" type="number" min="6" max="40"', fixed = TRUE)
  expect_match(html, "shape.setAttribute('fill'", fixed = TRUE)
  expect_match(html, "label.setAttribute('fill',contrast.fill)", fixed = TRUE)
  expect_match(html, '<option value="tailwind">Tailwind CSS inspired</option>', fixed = TRUE)
  expect_match(html, 'role="dialog"', fixed = TRUE)
  expect_match(html, ".networkPanel{position:absolute;inset:0", fixed = TRUE)
  expect_false(grepl(".networkPanel{position:fixed", html, fixed = TRUE))
  expect_match(html, '<option value="force" selected>Force</option>', fixed = TRUE)
  expect_match(html, '<option value="auto">Auto</option>', fixed = TRUE)
  expect_match(html, "fitNetworkView", fixed = TRUE)
  expect_match(html, "setLoading(true,'Loading topic GRN...')", fixed = TRUE)
  expect_match(html, "refreshConditionData()", fixed = TRUE)
  expect_match(html, "TF_EXPR_INDEX", fixed = TRUE)
  expect_match(html, "tfExpression(c1,r.tfu", fixed = TRUE)
  expect_match(html, "function setExpressionIndexValue", fixed = TRUE)
  expect_match(html, "indexEdgeExpressionRows(rows)", fixed = TRUE)
  expect_match(html, "function networkTfRnaValue", fixed = TRUE)
  expect_match(html, "function networkGeneRnaValue", fixed = TRUE)
  expect_match(html, "function networkExpressionLimit", fixed = TRUE)
  expect_match(html, "tfv=networkTfRnaValue(r.tfu)", fixed = TRUE)
  expect_match(html, "gv=networkGeneRnaValue(r.gene)", fixed = TRUE)
  expect_match(
    html,
    "RNA log2FC ('+conditionLabel(cond1Select.value)+' / '+conditionLabel(cond2Select.value)+')",
    fixed = TRUE
  )
  expect_no_match(html, "tfv=nodeLogfc(r.exprA,r.exprB)", fixed = TRUE)
  expect_no_match(html, "gv=nodeLogfc(r.gexprA,r.gexprB)", fixed = TRUE)
  expect_match(html, "showItemPage", fixed = TRUE)
  expect_true(craftgrn:::.m3tb_validate_current_report_html(out_file))

  index_file <- tempfile(fileext = ".html")
  craftgrn:::.m3cr_write_combined_report_index(
    index_file,
    data.table::data.table(
      label = "Condition pair",
      k = 2L,
      method_setup = "condition | LDA",
      run_id = "run_001",
      path = out_file
    ),
    report_state = list(
      condition_colors = list(A = "#E15759", B = "#4E79A7"),
      condition_order = c("B", "A"),
      defaults = list(condition_1 = "B")
    )
  )
  index_html <- paste(readLines(index_file, warn = FALSE), collapse = "\n")
  expect_match(index_html, "INITIAL_CONDITION_COLORS", fixed = TRUE)
  expect_match(index_html, '"condition_order":["B","A"]', fixed = TRUE)
  expect_match(index_html, '<label>Topic view <select id="conditionTopicMetric">', fixed = TRUE)
  expect_match(index_html, "metric:conditionTopicMetric.value", fixed = TRUE)
  expect_match(index_html, "sendState('metric')", fixed = TRUE)
  expect_match(index_html, "label:'Overall pathways'", fixed = TRUE)
  expect_match(index_html, "id==='cond2'&&!state.cond1", fixed = TRUE)
  expect_match(
    index_html,
    "document.getElementById('cond2Select').disabled=!state.cond1",
    fixed = TRUE
  )
  expect_match(
    index_html,
    "document.getElementById('tfSelect').disabled=!state.cond1",
    fixed = TRUE
  )
  expect_match(
    index_html,
    "const enabled=!!document.getElementById('cond1Select').value",
    fixed = TRUE
  )
  expect_match(index_html, "optionsSignature", fixed = TRUE)
  expect_match(index_html, "paletteSignature", fixed = TRUE)
  expect_false(grepl("localStorage", index_html, fixed = TRUE))
})

test_that("experimental gene-expression LDA is explicit and excluded from all", {
  experimental <- craftgrn:::.module3_topic_method_plan(
    methods = "condition_gene_expression_lda",
    k_grid = c(20L, 30L)
  )
  expect_equal(experimental$fp_mode, "gene_expression")
  expect_equal(experimental$backend, "warplda")
  expect_equal(experimental$weight_label, "gene_expression")
  expect_true(
    experimental$weight_label %in%
      eval(formals(craftgrn:::extract_regulatory_topics)$weight_label)
  )
  expect_true(experimental$experimental)
  standard <- craftgrn:::.module3_topic_method_plan(methods = "all", k_grid = 10L)
  expect_false("condition_gene_expression_lda" %in% standard$method)
})

test_that("experimental gene-expression MultiVI mirrors the LDA documents", {
  lda <- craftgrn:::.module3_topic_method_plan(
    methods = "condition_gene_expression_lda",
    k_grid = 30L
  )
  multivi <- craftgrn:::.module3_topic_method_plan(
    methods = "condition_gene_expression_multivi",
    k_grid = 30L
  )

  expect_equal(multivi$fp_mode, "gene_expression")
  expect_equal(multivi$backend, "vae")
  expect_equal(multivi$vae_variant, "multivi_encoder")
  expect_equal(multivi$weight_label, lda$weight_label)
  expect_equal(multivi$setup, lda$setup)
  expect_true(multivi$experimental)
  standard <- craftgrn:::.module3_topic_method_plan(methods = "all", k_grid = 30L)
  expect_false("condition_gene_expression_multivi" %in% standard$method)
})

test_that("condition pathway activity reuses fixed overall pathway genes", {
  pathways <- data.table::data.table(
    topic_num = 1L,
    pathway_key = "path:a",
    pathway = "Path A",
    padj = 0.01,
    genes = "G1;G2"
  )
  expression <- data.table::data.table(
    condition_id = rep(c("A", "B"), each = 6L),
    gene_key = rep(c("G1", "g2", "G3", "G4", "G5", "G6"), 2L),
    gene_expr = c(16, 8, 1, 2, 3, 4, 1, 2, 16, 8, 7, 6)
  )
  gene_topics <- data.table::data.table(
    gene_key = c("G1", "G2", "G3"),
    topic_num = c(1L, 1L, 2L)
  )
  out <- craftgrn:::.m3cr_pathway_expression_scores(
    pathways,
    expression,
    gene_topics,
    min_genes = 2L
  )
  expect_equal(nrow(out$activity), 2L)
  expect_true(out$activity[condition_id == "A", rank_enrichment] > 0)
  expect_true(out$activity[condition_id == "B", rank_enrichment] < 0)
  expect_true(out$activity[condition_id == "A", mean_gene_zscore] > out$activity[condition_id == "B", mean_gene_zscore])
  expect_true(all(out$activity$n_expression_universe == 6L))
  expect_setequal(out$gene_expression$gene_key, c("G1", "g2"))
  expect_false("G3" %in% out$gene_expression$gene_key)
  expect_true(all(c(
    "expression_rank", "gene_zscore", "n_expression_universe"
  ) %in% names(out$gene_expression)))
})

test_that("condition pathway activity handles no expression overlap", {
  out <- craftgrn:::.m3cr_pathway_expression_scores(
    data.table::data.table(
      topic_num = 1L,
      pathway_key = "path:a",
      genes = "G1;G2"
    ),
    data.table::data.table(
      condition_id = "A",
      gene_key = "G3",
      gene_expr = 4
    ),
    data.table::data.table(gene_key = "G3", topic_num = 1L)
  )
  expect_equal(nrow(out$activity), 0L)
  expect_equal(nrow(out$gene_expression), 0L)
  expect_named(
    out$activity,
    c(
      "topic_num", "pathway_key", "condition_id", "rank_enrichment",
      "mean_gene_zscore", "mean_log2_expression", "n_expression_genes",
      "n_expression_universe"
    )
  )
  expect_named(
    out$gene_expression,
    c(
      "condition_id", "gene_key", "gene_expr", "log2_expression",
      "expression_rank", "gene_zscore", "n_expression_universe"
    )
  )
})

test_that("condition report state carries configured link-direction cutoffs", {
  state <- craftgrn:::.module3_report_state(list(
    gene_log2fc_cutoff = 1.25,
    fp_delta_cutoff = 0.3,
    fp_filter_mode = "delta",
    threshold_gene_expr = 2,
    topic_condition_gene_expression_min = 10,
    tf_opposition_log2fc_cutoff = 0.8
  ))
  expect_equal(state$link_direction$gene_log2fc_cutoff, 1.25)
  expect_equal(state$link_direction$fp_cutoff, 0.3)
  expect_equal(state$link_direction$expression_min, 10)
  expect_equal(state$link_direction$tf_opposition_log2fc_cutoff, 0.8)
  expect_equal(state$link_direction$fp_filter_mode, "delta")
})

test_that("condition topic expression share uses unique assigned target genes", {
  expression <- data.table::data.table(
    condition_id = rep(c("A", "B"), each = 3L),
    gene_key = rep(c("G1", "G2", "G3"), 2L),
    gene_expr = c(15, 3, 1, 1, 3, 15)
  )
  topics <- data.table::data.table(
    gene_key = c("G1", "G2", "G3", "G3"),
    topic_num = c(1L, 1L, 2L, 1L),
    topic_score = c(0.9, 0.8, 0.9, 0.1)
  )

  out <- craftgrn:::.m3cr_condition_topic_expression_share(expression, topics)

  expect_equal(out[, sum(expression_topic_share), by = condition_id]$V1, c(1, 1))
  expect_gt(
    out[condition_id == "A" & topic_num == 1L, expression_topic_share],
    out[condition_id == "B" & topic_num == 1L, expression_topic_share]
  )
  expect_equal(
    out[condition_id == "A" & topic_num == 2L, n_expression_genes],
    1L
  )
})

test_that("condition pathway expression normalizes Ensembl gene identifiers", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Mm.eg.db")
  expression <- data.table::data.table(
    condition_id = "A",
    gene_key = c("ENSMUSG00000000001", "Gnai3"),
    gene_expr = c(3, 5)
  )
  observed <- craftgrn:::.m3cr_normalize_expression_gene_keys(expression)
  expect_setequal(observed$gene_key, c("Gnai3"))
})

test_that("Module 3 condition payload filenames are Windows-safe and deterministic", {
  name <- "run_001_condition_aggr_lda_K20_condition_pair"
  stem <- craftgrn:::.m3cr_payload_stem(name)
  expect_match(stem, "^cp_[a-f0-9]{12}$")
  expect_identical(stem, craftgrn:::.m3cr_payload_stem(name))
  expect_false(identical(stem, craftgrn:::.m3cr_payload_stem(paste0(name, "_other"))))
  expect_lte(nchar(paste0(stem, "_c99_e.js")), 30L)

  report_file <- craftgrn:::.m3cr_report_data_file(
    "run_001_condition_aggr_lda_K20_condition_topic_report.html"
  )
  expect_match(report_file, "^cr_[a-f0-9]{12}_d[.]js$")
  expect_lte(nchar(report_file), 20L)

  subgrn_stem <- craftgrn:::.m3tb_subgrn_payload_stem(
    "run_002_condition_aggr_multivi_K20_pathway_subgrn"
  )
  expect_match(subgrn_stem, "^sg_[a-f0-9]{12}$")
  expect_lte(nchar(paste0(subgrn_stem, "_chunk_999.js")), 30L)
})

test_that("Module 3 stage methods preserve the full method layout", {
  root <- tempfile("module3-stage-methods-")
  dir.create(root, recursive = TRUE)
  res <- craftgrn:::run_module3_topic_benchmark(
    filtered_dir = root,
    comparisons = data.table::data.table(condition_label = "A", condition_group = "A"),
    output_dir = root,
    methods = c("condition_aggr_lda", "condition_aggr_multivi"),
    training_methods = "condition_aggr_multivi",
    extraction_methods = "condition_aggr_multivi",
    k_grid = 10L,
    output_layout = "benchmark",
    run_training = FALSE,
    run_extraction = FALSE,
    run_reports = FALSE,
    path_length_check = FALSE,
    verbose = FALSE
  )
  expect_equal(res$method_plan$run_id, c("run_001", "run_002"))
  expect_equal(res$method_plan$method, c("condition_aggr_lda", "condition_aggr_multivi"))
  expect_error(
    craftgrn:::run_module3_topic_benchmark(
      filtered_dir = root,
      comparisons = data.table::data.table(condition_label = "A", condition_group = "A"),
      output_dir = root,
      methods = "condition_aggr_lda",
      training_methods = "condition_aggr_multivi",
      run_training = FALSE,
      run_extraction = FALSE,
      run_reports = FALSE,
      path_length_check = FALSE,
      verbose = FALSE
    ),
    "training_methods must be a subset"
  )
})

test_that("Module 3 condition review groups multiple K values inside one method report", {
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
  condition_pages <- list.files(file.path(review_dir, "assets", "p"), pattern = "[.]html$", full.names = TRUE)
  expect_false(dir.exists(file.path(review_dir, "topic_reports")))
  expect_false(dir.exists(file.path(review_dir, "condition_topic_reports")))
  expect_equal(length(condition_pages), 4L)

  condition_html <- paste(readLines(file.path(review_dir, "topic_method_k_condition_mds_report.html"), warn = FALSE), collapse = "\n")
  expect_match(condition_html, "K <select", fixed = TRUE)
  expect_match(condition_html, "K2", fixed = TRUE)
  expect_match(condition_html, "K3", fixed = TRUE)
  expect_match(condition_html, "assets/p/", fixed = TRUE)
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
  expect_true("topic_score_method" %in% names(formals(module3_extract_topics)))
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

test_that("Module 3 topic benchmark flattens model and extraction roots outside legacy layout", {
  expect_true(.m3tb_flat_output_for_layout("standard"))
  expect_true(.m3tb_flat_output_for_layout("benchmark"))
  expect_false(.m3tb_flat_output_for_layout("legacy"))
})

test_that("Module 3 benchmark extraction writes directly into each K folder", {
  root <- tempfile("module3-topic-benchmark-flat-extract-")
  dir.create(root, recursive = TRUE)

  setup <- .module3_topic_method_plan("comparison_aggr_multivi", k_grid = 10L)
  plan <- .m3tb_apply_output_layout(setup, root, "benchmark")
  model_name <- paste0(
    "GSE192390_vae_joint_tf_docs_",
    plan$weight_label[[1L]],
    "_multivi_encoder_Kgrid"
  )
  model_dir <- file.path(plan$topic_models_dir[[1L]], model_name)
  dir.create(file.path(model_dir, "vae_models"), recursive = TRUE)
  dir.create(file.path(model_dir, "rds"), recursive = TRUE)

  theta <- data.table::data.table(
    doc_id = c("CompA::KLF5::Target-Up", "CompA::KLF5::Target-Down"),
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  phi <- data.table::data.table(
    topic = c("Topic1", "Topic2"),
    `GENE:G1` = c(0.8, 0.2),
    `GENE:G2` = c(0.2, 0.8),
    `PEAK:G1` = c(0.7, 0.3),
    `PEAK:G2` = c(0.3, 0.7)
  )
  data.table::fwrite(theta, file.path(model_dir, "vae_models", "theta_K10.csv"))
  data.table::fwrite(phi, file.path(model_dir, "vae_models", "phi_K10.csv"))
  dtm <- Matrix::Matrix(
    matrix(c(2, 1, 2, 1, 1, 2, 1, 2), nrow = 2, byrow = TRUE),
    sparse = TRUE
  )
  rownames(dtm) <- theta$doc_id
  colnames(dtm) <- c("GENE:G1", "GENE:G2", "PEAK:G1", "PEAK:G2")
  edges_docs <- data.table::data.table(
    doc_id = theta$doc_id,
    tf = "KLF5",
    gene_key = c("G1", "G2"),
    peak_id = c("P1", "P2")
  )
  saveRDS(dtm, file.path(model_dir, "rds", "dtm.rds"))
  saveRDS(edges_docs, file.path(model_dir, "rds", "edges_docs.rds"))

  run_module3_topic_benchmark(
    filtered_dir = tempfile("unused-filtered-"),
    output_dir = root,
    comparisons = data.table::data.table(
      comparison_id = "CompA",
      comparison_display = "CompA",
      cond1_label = "Cond1",
      cond2_label = "Cond2"
    ),
    methods = "comparison_aggr_multivi",
    k_grid = 10L,
    output_layout = "benchmark",
    run_training = FALSE,
    run_extraction = TRUE,
    run_reports = FALSE,
    extraction_topic_report_args = list(
      extraction_steps = c("topic_terms", "tf_topic_assignment"),
      run_link_topic_scores = FALSE,
      run_pathway_enrichment = FALSE
    ),
    verbose = FALSE
  )

  k_dir <- file.path(plan$topic_extraction_dir[[1L]], "K10")
  expect_true(file.exists(file.path(k_dir, "topic_terms.csv")))
  expect_true(file.exists(file.path(k_dir, "topic_gene_peak_assignment.csv")))
  expect_true(file.exists(file.path(k_dir, "topic_term_assignment_summary.csv")))
  expect_equal(
    readLines(file.path(k_dir, "topic_term_assignment_method.txt")),
    "gammafit_maxprob"
  )
  assignment_summary <- data.table::fread(file.path(k_dir, "topic_term_assignment_summary.csv"))
  expect_equal(assignment_summary$model_family, "multivi")
  expect_equal(assignment_summary$gammafit_thrP, 0.5)
  expect_false(file.exists(file.path(k_dir, "raw_theta_documents_K2.pdf")))
  expect_false(dir.exists(file.path(k_dir, "tf_topic_assignment")))
  expect_false(dir.exists(file.path(k_dir, model_name)))
})

test_that("Module 3 path-length sanity check reports Windows-alias paths", {
  root <- tempfile("module3-path-length-")
  long_dir <- file.path(root, paste(rep("deepfolder", 4), collapse = .Platform$file.sep))
  dir.create(long_dir, recursive = TRUE)
  long_file <- file.path(long_dir, paste0(paste(rep("longname", 6), collapse = "_"), ".csv"))
  writeLines("x", long_file)

  report_file <- file.path(root, "path_length_sanity_check.csv")
  expect_warning(
    summary <- .m3tb_check_path_lengths(
      root = root,
      max_chars = 60L,
      path_alias = "Z:/project",
      out_file = report_file,
      verbose = TRUE
    ),
    "path length"
  )

  expect_true(file.exists(report_file))
  report <- data.table::fread(report_file)
  expect_true(any(report$over_limit))
  expect_true(any(grepl("^Z:/project", report$display_path)))
  expect_gt(summary$n_over_limit, 0L)
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
    c("method_order", "method_setup", "setup", "model_label", "selected_k", "status", "count", "count_basis")
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
  expect_true("topic_score_method" %in% names(formals(run_topic_modeling)))
  expect_true("topic_term_assignment_method" %in% names(formals(run_topic_modeling)))
  expect_true("condition_gene_weighting" %in% names(formals(run_topic_modeling)))
  expect_true("condition_gene_weighting" %in% names(formals(train_topic_models)))
  expect_true("condition_peak_weighting" %in% names(formals(run_topic_modeling)))
  expect_true("condition_peak_weighting" %in% names(formals(train_topic_models)))
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
    topic_score_method = "rowmax_phi",
    topic_term_assignment_method = "gammafit",
    topic_count_method = "bin",
    topic_count_input = "pseudo_count_bin",
    topic_condition_gene_weighting = "specificity",
    topic_condition_peak_weighting = "tf_expression",
    topic_condition_gene_expression_file = "condition_expression.csv",
    topic_condition_specificity_temperature = 0.5,
    topic_condition_specificity_floor = 0.1,
    topic_condition_specificity_expression_min = 10,
    topic_vae_device = "cuda",
    topic_vae_batch_size = 512L,
    pathway_databases = c("Reactome_2022", "KEGG_2021_Human"),
    report = list(
      condition_colors = list(A = "#E15759", B = "#4E79A7"),
      condition_order = c("B", "A"),
      defaults = list(condition_1 = "B")
    ),
    topic_benchmark_enabled = TRUE,
    topic_benchmark_methods = c("condition_aggr_weight_lda", "comparison_aggr_weight_lda"),
    topic_benchmark_k_grid = c(5L, 6L)
  )
  resolved <- .module3_resolve_topic_run_config(project_config = cfg)
  expect_equal(resolved$method, "comparison_aggr_weight_lda")
  expect_equal(resolved$k_grid, c(8L, 10L))
  expect_equal(resolved$warplda_iterations, 25L)
  expect_equal(resolved$topic_link_output, "none")
  expect_equal(resolved$extraction_args$topic_score_method, "rowmax_phi")
  expect_equal(resolved$extraction_args$topic_term_assignment_method, "gammafit")
  expect_equal(resolved$count_method, "bin")
  expect_equal(resolved$count_input, "pseudo_count_bin")
  expect_equal(resolved$condition_gene_weighting, "specificity")
  expect_equal(resolved$condition_peak_weighting, "tf_expression")
  expect_equal(resolved$condition_gene_expression_file, "condition_expression.csv")
  expect_equal(resolved$condition_specificity_temperature, 0.5)
  expect_equal(resolved$condition_specificity_floor, 0.1)
  expect_equal(resolved$condition_specificity_expression_min, 10)
  expect_equal(resolved$vae_device, "cuda")
  expect_equal(resolved$vae_batch_size, 512L)
  expect_equal(resolved$extraction_args$pathway_databases, cfg$pathway_databases)
  expect_equal(resolved$report_state$condition_colors$A, "#E15759")
  expect_equal(resolved$report_state$condition_order, c("B", "A"))
  expect_equal(resolved$report_state$defaults$condition_1, "B")
  expect_true(resolved$benchmark$enabled)
  expect_equal(resolved$benchmark$methods, cfg$topic_benchmark_methods)
  expect_equal(resolved$benchmark$k_grid, c(5L, 6L))

  overridden <- .module3_resolve_topic_run_config(
    project_config = cfg,
    method = "condition_aggr_lda",
    k_grid = 12L,
    warplda_iterations = 3L,
    topic_link_output = "pass",
    topic_score_method = "normtop_specificity",
    topic_term_assignment_method = "max_phi",
    count_method = "log",
    count_input = "pseudo_count_log",
    vae_device = "cpu",
    vae_batch_size = 128L
  )
  expect_equal(overridden$method, "condition_aggr_lda")
  expect_equal(overridden$k_grid, 12L)
  expect_equal(overridden$warplda_iterations, 3L)
  expect_equal(overridden$topic_link_output, "pass")
  expect_equal(overridden$extraction_args$topic_score_method, "normtop_specificity")
  expect_equal(overridden$extraction_args$topic_term_assignment_method, "max_phi")
  expect_equal(overridden$count_method, "log")
  expect_equal(overridden$count_input, "pseudo_count_log")
  expect_equal(overridden$vae_device, "cpu")
  expect_equal(overridden$vae_batch_size, 128L)
})

test_that("Module 3 topic link defaults do not apply gene_prob max filtering", {
  resolved_default <- .module3_resolve_topic_run_config(project_config = list())
  expect_equal(resolved_default$extraction_args$link_topic_method, "gammafit")
  expect_equal(resolved_default$extraction_args$topic_score_method, "normtop_specificity")
  expect_equal(resolved_default$extraction_args$topic_term_assignment_method, "gammafit_maxprob")
  expect_equal(resolved_default$extraction_args$thrP, 0.5)
  expect_equal(resolved_default$extraction_args$link_topic_prob_cutoff, 0.3)
  expect_equal(resolved_default$extraction_args$topic_tf_membership_cutoff, 0.3)
  expect_equal(resolved_default$extraction_args$topic_tf_primary_margin_cutoff, 0.1)
  expect_false(resolved_default$extraction_args$run_raw_theta_document_heatmap)
  expect_false("topic_tf_assignment_browser" %in% names(resolved_default$extraction_args))

  resolved_raw_theta <- .module3_resolve_topic_run_config(
    project_config = list(topic_raw_theta_document_heatmap = TRUE)
  )
  expect_true(resolved_raw_theta$extraction_args$run_raw_theta_document_heatmap)

  resolved_lda <- .module3_resolve_topic_run_config(
    project_config = list(),
    method = "condition_aggr_lda"
  )
  expect_equal(resolved_lda$extraction_args$topic_term_assignment_method, "gammafit_maxprob")
  expect_equal(resolved_lda$extraction_args$thrP, 0.70)

  resolved_multivi <- .module3_resolve_topic_run_config(
    project_config = list(),
    method = "condition_aggr_multivi"
  )
  expect_equal(resolved_multivi$extraction_args$thrP, 0.50)

  resolved_weighted <- .module3_resolve_topic_run_config(
    project_config = list(),
    method = "condition_aggr_weight_lda"
  )
  expect_equal(resolved_weighted$extraction_args$topic_term_assignment_method, "max_phi")

  resolved_mapped <- .module3_resolve_topic_run_config(
    project_config = list(topic_gammafit_thrP = list(lda = 0.79, multivi = 0.87)),
    method = "condition_aggr_multivi"
  )
  expect_equal(resolved_mapped$extraction_args$thrP, 0.87)

  resolved_explicit <- .module3_resolve_topic_run_config(
    project_config = list(
      topic_link_method = "gene_prob",
      topic_link_prob_cutoff = "max",
      topic_tf_membership_cutoff = 0.25,
      topic_tf_primary_margin_cutoff = 0.05
    )
  )
  expect_equal(resolved_explicit$extraction_args$link_topic_method, "gene_prob")
  expect_equal(resolved_explicit$extraction_args$link_topic_prob_cutoff, "max")
  expect_equal(resolved_explicit$extraction_args$topic_tf_membership_cutoff, 0.25)
  expect_equal(resolved_explicit$extraction_args$topic_tf_primary_margin_cutoff, 0.05)
  expect_false("topic_tf_assignment_browser" %in% names(resolved_explicit$extraction_args))

  resolved_theta_terms <- .module3_resolve_topic_run_config(
    project_config = list(topic_link_method = "theta_and_terms")
  )
  expect_equal(resolved_theta_terms$extraction_args$link_topic_method, "theta_and_terms")
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

test_that("Module 3 topic run config resolves pathway species", {
  resolved_mouse <- .module3_resolve_topic_run_config(
    project_config = list(ref_genome = "mm10")
  )
  expect_equal(resolved_mouse$extraction_args$pathway_species, "mouse")

  resolved_human <- .module3_resolve_topic_run_config(
    project_config = list(ref_genome = "mm10", pathway_species = "human")
  )
  expect_equal(resolved_human$extraction_args$pathway_species, "human")
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

test_that("Module 3 extraction inherits FP term mode from topic input summary", {
  model_defaults <- paste(
    deparse(body(extract_regulatory_topics)),
    collapse = "\n"
  )
  expect_match(model_defaults, "topic_input_summary_path", fixed = TRUE)
  expect_match(model_defaults, "topic_report_args$fp_term_mode", fixed = TRUE)
  expect_match(model_defaults, "Using fp_term_mode=", fixed = TRUE)
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

test_that("Module 3 review worker count is adaptive and memory bounded", {
  old_workers <- Sys.getenv("CRAFTGRN_REVIEW_WORKERS", unset = NA_character_)
  old_threads <- Sys.getenv("CRAFTGRN_REVIEW_THREADS", unset = NA_character_)
  withr::defer({
    if (is.na(old_workers)) Sys.unsetenv("CRAFTGRN_REVIEW_WORKERS") else Sys.setenv(CRAFTGRN_REVIEW_WORKERS = old_workers)
    if (is.na(old_threads)) Sys.unsetenv("CRAFTGRN_REVIEW_THREADS") else Sys.setenv(CRAFTGRN_REVIEW_THREADS = old_threads)
  })
  withr::local_options(list(
    craftgrn.review.max_workers = 8L,
    craftgrn.review.memory_reserve_gb = 32,
    craftgrn.review.memory_per_worker_gb = 6,
    craftgrn.review.auto_workers_in_tests = TRUE
  ))

  Sys.setenv(CRAFTGRN_REVIEW_WORKERS = "36", CRAFTGRN_REVIEW_THREADS = "1")
  expect_equal(
    craftgrn:::.m3tb_review_worker_count(96L, 128 * 1024^3, 36L),
    8L
  )
  expect_equal(
    craftgrn:::.m3tb_review_worker_count(4L, 128 * 1024^3, 36L),
    4L
  )
  expect_equal(
    craftgrn:::.m3tb_review_worker_count(96L, 40 * 1024^3, 36L),
    1L
  )

  Sys.unsetenv("CRAFTGRN_REVIEW_WORKERS")
  Sys.setenv(CRAFTGRN_REVIEW_THREADS = "20")
  expect_equal(
    craftgrn:::.m3tb_review_worker_count(96L, 128 * 1024^3, 36L),
    8L
  )

  Sys.unsetenv("CRAFTGRN_REVIEW_THREADS")
  expect_equal(
    craftgrn:::.m3tb_review_worker_count(8L, 128 * 1024^3, 36L),
    8L
  )
})

test_that("Module 3 review mapping retries sequentially after worker failure", {
  testthat::local_mocked_bindings(
    .m3tb_review_par_lapply = function(cl, x, fun) {
      stop("simulated worker connection loss")
    },
    .package = "craftgrn"
  )

  expect_warning(
    result <- craftgrn:::.m3tb_review_lapply(
      as.list(1:3),
      function(x) x * 2,
      workers = 2L
    ),
    "retrying 3 task\\(s\\) sequentially"
  )
  expect_equal(unlist(result), c(2, 4, 6))
})

test_that("Module 3 review parallel helper runs socket workers", {
  skip_on_cran()
  out <- craftgrn:::.m3tb_review_lapply(1:4, function(i) i * i, workers = 2L)
  expect_equal(unlist(out), c(1L, 4L, 9L, 16L))
})

test_that("Module 3 review socket workers inherit parent library paths", {
  skip_on_cran()
  parent_library <- tempfile("craftgrn-worker-library-")
  dir.create(parent_library)
  old_paths <- .libPaths()
  withr::defer(.libPaths(old_paths))
  .libPaths(c(parent_library, old_paths))

  out <- craftgrn:::.m3tb_review_lapply(
    1:2,
    function(i) list(
      libraries = normalizePath(.libPaths(), winslash = "/", mustWork = FALSE),
      package = normalizePath(
        getNamespaceInfo(asNamespace("craftgrn"), "path"),
        winslash = "/",
        mustWork = FALSE
      )
    ),
    workers = 2L
  )
  expected <- normalizePath(parent_library, winslash = "/", mustWork = FALSE)
  expected_package <- normalizePath(find.package("craftgrn", lib.loc = old_paths), winslash = "/")
  expect_true(all(vapply(out, function(x) expected %in% x$libraries, logical(1L))))
  expect_true(all(vapply(out, function(x) identical(x$package, expected_package), logical(1L))))
})

test_that("condition report payload keeps bounded strongest peak support", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")
  path <- tempfile(fileext = ".parquet")
  links <- data.table::data.table(
    condition_id = rep("condition_a", 5L),
    tf = c("TF1", "TF1", "TF1", "TF2", "TF2"),
    gene_key = c("Gene1", "Gene1", "Gene1", "Gene1", "Gene2"),
    peak_id = c("peak_a", "peak_a", "peak_b", "peak_a", "peak_c"),
    fp_score_condition = c(0.4, 0.8, 0.6, 0.7, 0.5),
    gene_expr_condition = rep(2, 5L),
    tf_expr_condition = rep(3, 5L)
  )
  arrow::write_parquet(links, path)

  out <- craftgrn:::.m3cr_collect_condition_edges(
    path,
    topic_genes = c("Gene1", "Gene2"),
    max_peaks_per_tf_gene = 1L
  )

  expect_equal(nrow(out$tf_peak_gene), 3L)
  expect_equal(
    out$tf_peak_gene[tf_upper == "TF1" & gene_key == "Gene1", peak_id],
    "peak_a"
  )
  expect_equal(
    out$tf_peak_gene[tf_upper == "TF1" & gene_key == "Gene1", fp_score],
    0.8
  )
  expect_equal(
    out$tf_peak_gene[gene_key == "Gene1", uniqueN(tf_upper)],
    2L
  )
})

test_that("Module 3 topic-link summary forces path arguments before socket workers", {
  old_workers <- Sys.getenv("CRAFTGRN_REVIEW_WORKERS", unset = NA_character_)
  withr::defer({
    if (is.na(old_workers)) Sys.unsetenv("CRAFTGRN_REVIEW_WORKERS") else Sys.setenv(CRAFTGRN_REVIEW_WORKERS = old_workers)
  })
  Sys.setenv(CRAFTGRN_REVIEW_WORKERS = "2")
  root <- tempfile("topic-review-force-")
  dir.create(root, recursive = TRUE)
  rows <- data.table::data.table(
    method = c("m1", "m2"),
    selected_k = c(2L, 3L),
    method_order = c(1L, 2L),
    method_setup = c("setup 1", "setup 2"),
    setup = c("setup1", "setup2"),
    model_label = c("LDA", "VAE"),
    topic_extraction_dir = file.path(root, c("missing1", "missing2"))
  )
  call_summary <- function(benchmark_dir) {
    craftgrn:::.m3tb_summarize_topic_links(
      output_dir = benchmark_dir,
      model_rows = rows,
      review_dir = file.path(benchmark_dir, "review")
    )
  }
  out <- call_summary(root)
  expect_s3_class(out$pass, "data.table")
  expect_equal(nrow(out$pass), 0L)
})

test_that("Module 3 review exposes complete raw and combined topic spaces", {
  root <- tempfile("topic-spaces-")
  extraction_dir <- file.path(root, "topic_extraction", "K4")
  model_dir <- file.path(root, "topic_models")
  dir.create(file.path(extraction_dir, "rds"), recursive = TRUE)
  dir.create(file.path(model_dir, "vae_models"), recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(
      raw_topic = 1:4,
      optimized_topic = c(1L, 1L, 2L, 3L)
    ),
    file.path(extraction_dir, "topic_optimization_map.csv")
  )
  data.table::fwrite(
    data.table::data.table(
      topic = 1:3,
      term_id = paste0("GENE:G", 1:3),
      in_topic = TRUE
    ),
    file.path(extraction_dir, "topic_terms.csv")
  )
  optimized_theta <- matrix(
    c(0.7, 0.2, 0.1, 0.1, 0.3, 0.6),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("A::TF1", "B::TF1"), paste0("Topic", 1:3))
  )
  saveRDS(
    optimized_theta,
    file.path(extraction_dir, "rds", "topic_theta_optimized.rds")
  )
  data.table::fwrite(
    data.table::data.table(
      doc_id = c("A::TF1", "B::TF1"),
      Topic1 = c(0.5, 0.1),
      Topic2 = c(0.2, 0.2),
      Topic3 = c(0.2, 0.3),
      Topic4 = c(0.1, 0.4)
    ),
    file.path(model_dir, "vae_models", "theta_K4.csv")
  )
  rows <- data.table::data.table(
    selected_k = 4L,
    method_order = 1L,
    method_setup = "cond fp aggr | LDA",
    model_label = "LDA",
    model_dir = model_dir,
    topic_extraction_dir = file.path(root, "topic_extraction")
  )

  spaces <- craftgrn:::.m3tb_topic_space_rows(rows)

  expect_equal(spaces$topic_space, c("raw", "combined"))
  expect_equal(spaces$topic_count, c(4L, 3L))
  expect_equal(
    spaces$report_label,
    c("Raw K4", "Combined K3 (from K4)")
  )
  expect_equal(
    craftgrn:::.m3tb_topic_space_theta(spaces[topic_space == "combined"]),
    optimized_theta
  )
  expect_match(
    craftgrn:::.m3tb_topic_space_file(
      extraction_dir,
      "raw",
      "condition_pathways"
    ),
    "per_condition_topic_pathway_enrichment_raw[.]csv$"
  )
  expect_match(
    craftgrn:::.m3tb_topic_space_file(
      extraction_dir,
      "raw",
      "overall_pathway_dotplot"
    ),
    "topic_pathway_enrichment_dotplot_raw[.]csv$"
  )
  expect_equal(
    craftgrn:::.m3tb_review_report_slug(spaces[topic_space == "combined"], 4L),
    "lda_K4_combined_K3"
  )
})

test_that("raw topic-space reports fall back to canonical non-optimized files", {
  extraction_dir <- tempfile("non-optimized-topics-")
  dir.create(extraction_dir)
  canonical <- c(
    terms = "topic_terms.csv",
    overall_pathways = "topic_pathway_enrichment_topic_terms.csv",
    overall_pathway_dotplot = "topic_pathway_enrichment_dotplot.csv",
    condition_pathways = "per_condition_topic_pathway_enrichment.csv"
  )
  file.create(file.path(extraction_dir, canonical))

  for (kind in names(canonical)) {
    expect_identical(
      craftgrn:::.m3tb_topic_space_file(extraction_dir, "raw", kind),
      file.path(extraction_dir, canonical[[kind]])
    )
  }
})
