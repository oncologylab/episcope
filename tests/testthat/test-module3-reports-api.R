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

test_that("Module 3 Enrichr helpers default to fast request scheduling", {
  expect_equal(
    formals(craftgrn:::plot_topic_pathway_enrichment_heatmap)$enrichr_sleep_time,
    0
  )
  expect_equal(
    formals(craftgrn:::plot_topic_pathway_enrichment_from_link_scores)$enrichr_sleep_time,
    0
  )
  expect_null(
    formals(craftgrn:::plot_topic_pathway_enrichment_from_link_scores)$enrichr_cache_dir
  )
  expect_equal(
    formals(craftgrn:::plot_topic_pathway_enrichment_from_link_scores)$enrichr_n_cores,
    1L
  )
  expect_equal(
    formals(craftgrn:::run_tfdocs_report_from_topic_base)$pathway_enrichr_sleep_time,
    0
  )
  expect_null(
    formals(craftgrn:::run_tfdocs_report_from_topic_base)$pathway_enrichr_cache_dir
  )
  expect_null(
    formals(craftgrn:::run_tfdocs_report_from_topic_base)$pathway_enrichr_n_cores
  )
  expect_equal(
    formals(craftgrn:::run_diff_grn_pathway_enrichment)$enrichr_sleep_time,
    0
  )
  expect_null(
    formals(craftgrn:::run_diff_grn_pathway_enrichment)$enrichr_cache_dir
  )
  expect_equal(craftgrn:::.normalize_enrichr_sleep_time(-1), 0)
  expect_equal(craftgrn:::.normalize_enrichr_sleep_time(0.25), 0.25)
  expect_equal(craftgrn:::.normalize_enrichr_n_cores(1), 1)
  expect_equal(
    craftgrn:::.module3_default_enrichr_cache_dir(file.path("root", "topic_extraction", "K10", "model")),
    file.path("root", "topic_extraction", "cache", "enrichr")
  )
  expect_match(craftgrn:::.enrichr_cache_key(c("B", "A", "A"), c("DB2", "DB1")), "^[a-f0-9]+$")
})

test_that("Module 3 pathway backend supports local enrichly enrichment", {
  skip_if_not_installed("enrichly")
  old_opt <- options(craftgrn.pathway_backend = NULL, craftgrn.enrichly.db_cache = NULL)
  on.exit(options(old_opt), add = TRUE)

  expect_equal(craftgrn:::.pathway_backend(), "enrichly")
  expect_equal(craftgrn:::.pathway_backend("enrichr"), "enrichr")
  expect_error(craftgrn:::.pathway_backend("bad_backend"), "pathway_backend")
  expect_false(identical(
    craftgrn:::.enrichr_cache_key(c("B", "A"), "Toy_DB", backend = "enrichly"),
    craftgrn:::.enrichr_cache_key(c("B", "A"), "Toy_DB", backend = "enrichr")
  ))

  db_cache <- tempfile("enrichly-db-cache-")
  dir.create(db_cache, recursive = TRUE)
  writeLines(
    c(
      paste(c("Term one", "desc", "A", "B", "C"), collapse = "\t"),
      paste(c("Term two", "desc", "X", "Y", "Z"), collapse = "\t")
    ),
    file.path(db_cache, "Toy_DB.txt")
  )
  options(craftgrn.enrichly.db_cache = db_cache)

  res <- craftgrn:::.run_enrichr_cached(
    genes = c("A", "B"),
    dbs = "Toy_DB",
    cache_dir = tempfile("enrich-cache-"),
    backend = "enrichly"
  )

  expect_named(res, "Toy_DB")
  expect_true(all(c("Term", "P.value", "Adjusted.P.value", "Genes") %in% names(res$Toy_DB)))
  expect_equal(res$Toy_DB$Term[[1L]], "Term one")
  expect_equal(res$Toy_DB$Genes[[1L]], "A;B")

  cfg <- tempfile("module3-pathway-backend-", fileext = ".yaml")
  writeLines(c("pathway_backend: enrichr", "topic_k: 3"), cfg)
  resolved <- craftgrn:::.module3_resolve_topic_run_config(project_config = cfg)
  expect_equal(resolved$pathway_backend, "enrichr")
})

test_that("Module 3 QC report writes differential TF summaries", {
  root <- tempfile("module3-qc-topic-")
  diff_dir <- tempfile("module3-qc-diff-")
  dir.create(file.path(root, "review", "tables"), recursive = TRUE)
  dir.create(diff_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(method = "condition_aggr_weight_lda", method_setup = "cond fp aggr weight | LDA"),
    file.path(root, "review", "tables", "module3_topic_method_plan.csv")
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
  dir.create(file.path(topic_dir, "review"), recursive = TRUE)
  dir.create(diff_dir, recursive = TRUE)
  writeLines("<!doctype html><html><body>Topic MDS</body></html>", file.path(topic_dir, "review", "topic_method_k_topic_mds_report.html"))
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

test_that("Module 3 differential GRN browser writes interactive TF-gene network payloads", {
  diff_dir <- tempfile("module3-diff-network-")
  dir.create(diff_dir, recursive = TRUE)
  links_up <- data.table::data.table(
    comparison_id = "CondA_vs_CondB",
    comparison_group = "Group1",
    cond1_id = "CondA",
    cond2_id = "CondB",
    cond1_label = "Condition A",
    cond2_label = "Condition B",
    tf = c("TF1", "TF1", "TF1", "TF2", "TF2", "TF3"),
    gene_key = c("GeneA", "GeneA", "GeneB", "GeneC", "GeneD", "GeneE"),
    peak_id = paste0("peak", seq_len(6)),
    delta_link_score = c(2.0, 1.0, 1.5, 1.2, 0.7, 0.2),
    link_score_cond1 = c(2.2, 1.3, 1.8, 1.4, 0.9, 0.5),
    link_score_cond2 = c(0.2, 0.3, 0.3, 0.2, 0.2, 0.3),
    log2FC_tf_expr = c(1.2, 1.2, 1.2, 0.7, 0.7, 0.2),
    log2FC_gene_expr = c(1.1, 1.1, 0.8, 0.6, 0.3, 0.1),
    distance_to_tss = c(100, 120, 500, 800, 900, 1000),
    candidate_source = "tss_window",
    r_gene = c(0.6, 0.62, 0.5, 0.48, 0.41, 0.3),
    r_rna_gene = c(0.7, 0.69, 0.6, 0.5, 0.43, 0.32)
  )
  links_down <- data.table::copy(links_up[1:2])
  links_down[, `:=`(
    tf = c("TF4", "TF4"),
    gene_key = c("GeneF", "GeneG"),
    peak_id = c("peak7", "peak8"),
    delta_link_score = c(-2.1, -1.7),
    log2FC_tf_expr = c(-1.1, -1.1),
    log2FC_gene_expr = c(-0.9, -0.7)
  )]
  data.table::fwrite(links_up, file.path(diff_dir, "CondA_vs_CondB_filtered_links_up.csv"))
  data.table::fwrite(links_down, file.path(diff_dir, "CondA_vs_CondB_filtered_links_down.csv"))

  out <- visualize_differential_grns(diff_dir, top_tf_n = 2L, top_link_n = 3L, verbose = FALSE)

  expect_true(file.exists(out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "comparisonSelect", fixed = TRUE)
  expect_match(html, "directionSelect", fixed = TRUE)
  expect_match(html, "topTfN", fixed = TRUE)
  expect_match(html, "topLinkN", fixed = TRUE)
  expect_match(html, "id=\"topTfN\" type=\"number\" min=\"0\"", fixed = TRUE)
  expect_match(html, "id=\"topLinkN\" type=\"number\" min=\"0\"", fixed = TRUE)
  expect_match(html, "function limitCount", fixed = TRUE)
  expect_match(html, "return v===0?Infinity", fixed = TRUE)
  expect_match(html, "if(state.drag)", fixed = TRUE)
  expect_match(html, "redrawPositions", fixed = TRUE)
  expect_match(html, "nodeLayer.addEventListener('mousedown'", fixed = TRUE)
  expect_match(html, "DEFAULT_DIRECTION='up'", fixed = TRUE)
  expect_match(html, "networkSvg", fixed = TRUE)
  expect_match(html, "CondA_vs_CondB", fixed = TRUE)
  expect_match(html, "supporting footprints", fixed = TRUE)
  expect_true(file.exists(file.path(dirname(out), "differential_grn_browser_manifest.csv")))
  expect_true(file.exists(file.path(dirname(out), "differential_grn_networks", "differential_grn_edges.csv")))
  manifest_dt <- data.table::fread(file.path(dirname(out), "differential_grn_browser_manifest.csv"))
  expect_equal(anyDuplicated(names(manifest_dt)), 0L)
  edge_dt <- data.table::fread(file.path(dirname(out), "differential_grn_networks", "differential_grn_edges.csv"))
  expect_equal(edge_dt[from == "TF1" & to == "GeneA", n_supporting_links], 2L)
  expect_equal(edge_dt[from == "TF1" & to == "GeneA", best_peak_ID], "peak1")
  expect_true(all(c("comparison_id", "direction", "from", "to", "edge_score", "n_supporting_links") %in% names(edge_dt)))
})
