topic_report_script <- normalizePath(file.path(
  Sys.getenv(
    "CRAFTGRN_REPO_DIR",
    unset = normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
  ),
  "dev/benchmark/03_build_topic_method_k_html_reports.R"
), mustWork = FALSE)

source_benchmark_script <- function(path) {
  testthat::skip_if_not(file.exists(path), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(path))))
  on.exit(setwd(old_wd), add = TRUE)
  source(path, local = parent.frame())
}

test_that("topic report helpers compute waterfall means by comparison label", {
  source_benchmark_script(topic_report_script)

  theta <- matrix(
    c(
      0.8, 0.2,
      0.6, 0.4,
      0.1, 0.9
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(theta) <- c(
    "HPAFII_A_vs_B::TF1::Up",
    "HPAFII_A_vs_B::TF2::Up",
    "HPAFII_A_vs_B::TF1::Down"
  )
  colnames(theta) <- c("Topic1", "Topic2")

  out <- .topic_report_waterfall(theta, doc_design = "comparison")

  expect_equal(
    out[comparison_label == "HPAFII_A_vs_B::Up" & topic == "Topic1", theta_mean],
    0.7
  )
  expect_equal(
    out[comparison_label == "HPAFII_A_vs_B::Down" & topic == "Topic2", theta_mean],
    0.9
  )
})

test_that("topic report helper filters pathways by overlap and adjusted p-value", {
  source_benchmark_script(topic_report_script)

  pathways <- data.table::data.table(
    topic = c(1L, 1L, 1L, 2L),
    pathway = c("keep strongest", "small overlap", "not significant", "topic two"),
    pathway_key = c("DB:keep strongest", "DB:small overlap", "DB:not significant", "DB:topic two"),
    padj = c(0.001, 0.002, 0.2, 0.03),
    overlap = c("7/40", "4/30", "10/100", "5/50"),
    overlap_hits = c(7L, 4L, 10L, 5L),
    genes = c("A;B;C;D;E;F;G", "A;B;C;D", "A;B;C;D;E;F;G;H;I;J", "K;L;M;N;O")
  )

  out <- .topic_report_filter_pathways(pathways, min_genes = 5L, padj_cutoff = 0.05, top_n = 30L)

  expect_equal(out$pathway, c("keep strongest", "topic two"))
  expect_equal(out$gene_in, c(7L, 5L))
  expect_equal(out$gene_total, c(40L, 50L))
})

test_that("topic report aborts when extraction pathway CSVs are missing", {
  source_benchmark_script(topic_report_script)

  out_dir <- tempfile("topic-extraction-")
  dir.create(out_dir, recursive = TRUE)

  expect_error(
    .read_pathway_files(out_dir),
    "Missing topic pathway dotplot CSVs"
  )
})

test_that("topic report completeness requires pathway CSV source", {
  source_benchmark_script(topic_report_script)

  extraction_dir <- tempfile("topic-extraction-")
  sub_dir <- file.path(extraction_dir, "report")
  dir.create(sub_dir, recursive = TRUE)
  file.create(file.path(sub_dir, c("topic_links.csv", "topic_terms.csv")))

  expect_false(.complete_extraction(extraction_dir))

  pathway <- data.table::data.table(topic = 1L, pathway = "DB: Path", padj = 0.01, overlap = "5/20")
  data.table::fwrite(pathway, file.path(sub_dir, "topic_pathway_enrichment_gene_only_dotplot.csv"))

  expect_true(.complete_extraction(extraction_dir))
})

test_that("topic report extracts gene universe from doc_term table", {
  source_benchmark_script(topic_report_script)

  doc_term <- data.table::data.table(
    doc_id = c("Doc1", "Doc1", "Doc2", "Doc2"),
    term_id = c("GENE:A", "PEAK:p1", "GENE:B", "GENE:A"),
    weight = c(1, 2, 3, 4)
  )
  path <- tempfile(fileext = ".csv")
  data.table::fwrite(doc_term, path)

  out <- .topic_report_gene_universe_from_doc_term(path)

  expect_equal(out, c("A", "B"))
})

test_that("topic report extracts gene universe from doc_term rds table", {
  source_benchmark_script(topic_report_script)

  doc_term <- data.table::data.table(
    doc_id = c("Doc1", "Doc1", "Doc2", "Doc2"),
    term_id = c("GENE:C", "PEAK:p1", "GENE:B", "GENE:C"),
    weight = c(1, 2, 3, 4)
  )
  path <- tempfile(fileext = ".rds")
  saveRDS(doc_term, path)

  out <- .topic_report_gene_universe_from_doc_term(path)

  expect_equal(out, c("B", "C"))
})

test_that("topic report places rds doc_term universe cache beside Kgrid directory", {
  source_benchmark_script(topic_report_script)

  doc_term_file <- file.path(tempdir(), "model_Kgrid", "rds", "doc_term.rds")
  dir.create(dirname(doc_term_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(data.table::data.table(term_id = "GENE:A"), doc_term_file)

  out <- .topic_report_universe_cache_dir(doc_term_file, fallback_dir = tempdir())

  expect_equal(basename(out), "model_Kgrid")
})

test_that("topic report keeps non-significant universe pathway rows for totals", {
  source_benchmark_script(topic_report_script)

  topic_pathways <- data.table::data.table(
    topic = c(1L, 1L),
    topic_label = c("Topic1", "Topic1"),
    pathway_key = c("DB:Sig", "DB:NonSigUniverse"),
    pathway = c("DB: Sig", "DB: NonSigUniverse"),
    padj = c(0.01, 0.02),
    neglog10_padj = c(2, 1.7),
    overlap = c("6/100", "8/100"),
    gene_in = c(6L, 8L),
    gene_total = c(100L, 100L),
    genes = c("A;B;C;D;E;F", "A;B;C;D;E;F;G;H")
  )
  universe <- data.table::data.table(
    pathway = c("DB: Sig", "DB: NonSigUniverse"),
    padj = c(0.001, 0.9),
    overlap_hits = c(25L, 31L)
  )

  out <- .topic_report_apply_universe_pathway_counts(topic_pathways, universe)

  expect_equal(out[pathway == "DB: Sig", gene_total_universe], 25L)
  expect_equal(out[pathway == "DB: NonSigUniverse", gene_total_universe], 31L)
  expect_equal(out[pathway == "DB: NonSigUniverse", gene_out], 23L)
})

test_that("topic report formats universe pathways with extraction formatter", {
  source_benchmark_script(topic_report_script)

  enrichr_result <- list(
    GO_Cellular_Component_2023 = data.frame(
      Term = "Nucleus (GO:0005634)",
      Adjusted.P.value = 0.001,
      P.value = 0.0001,
      Overlap = "2228/4487",
      Genes = "A;B",
      check.names = FALSE
    )
  )
  universe <- .topic_report_format_universe_enrichr(enrichr_result)

  expect_equal(universe$pathway, "GO:CC: Nucleus")
  expect_equal(universe$overlap_hits, 2228L)
})

test_that("topic report joins universe totals by exact pathway label", {
  source_benchmark_script(topic_report_script)

  topic_pathways <- data.table::data.table(
    topic = 1L,
    topic_label = "Topic1",
    pathway_key = "GO:CC: Nucleus",
    pathway = "GO:CC: Nucleus",
    padj = 0.01,
    neglog10_padj = 2,
    overlap = "176/4487",
    gene_in = 176L,
    gene_total = 4487L,
    genes = "A;B"
  )
  universe <- data.table::data.table(
    pathway = "GO:CC: Nucleus",
    padj = 0.001,
    overlap = "2228/4487",
    overlap_hits = 2228L
  )

  out <- .topic_report_apply_universe_pathway_counts(topic_pathways, universe)

  expect_equal(out$gene_total_universe, 2228L)
  expect_equal(out$gene_out, 2052L)
})

test_that("topic report aborts when universe pathway labels do not match exactly", {
  source_benchmark_script(topic_report_script)

  topic_pathways <- data.table::data.table(
    topic = 1L,
    topic_label = "Topic1",
    pathway_key = "GO:CC: Nucleus",
    pathway = "GO:CC: Nucleus",
    padj = 0.01,
    neglog10_padj = 2,
    overlap = "176/4487",
    gene_in = 176L,
    gene_total = 4487L,
    genes = "A;B"
  )
  universe <- data.table::data.table(
    pathway = "Nucleus (GO:0005634)",
    padj = 0.001,
    overlap = "2228/4487",
    overlap_hits = 2228L
  )

  expect_error(
    .topic_report_apply_universe_pathway_counts(topic_pathways, universe),
    "Missing full gene-universe pathway totals"
  )
})

test_that("topic report MDS returns one coordinate row per document", {
  source_benchmark_script(topic_report_script)

  theta <- matrix(
    c(
      0.8, 0.2,
      0.6, 0.4,
      0.1, 0.9
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(theta) <- c("DocA::TF1::Up", "DocA::TF2::Up", "DocB::TF1::Down")
  colnames(theta) <- c("Topic1", "Topic2")

  out <- .topic_report_doc_mds(theta, doc_design = "comparison")

  expect_equal(nrow(out), 3L)
  expect_true(all(c("doc_id", "comparison_label", "MDS1", "MDS2") %in% names(out)))
  expect_true(all(is.finite(out$MDS1)))
  expect_true(all(is.finite(out$MDS2)))
})

test_that("topic report topic map returns one coordinate row per topic", {
  source_benchmark_script(topic_report_script)

  theta <- matrix(
    c(
      0.8, 0.1, 0.1,
      0.6, 0.3, 0.1,
      0.1, 0.2, 0.7,
      0.2, 0.2, 0.6
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(theta) <- paste0("Doc", seq_len(4))
  colnames(theta) <- c("Topic1", "Topic2", "Topic3")
  gene_counts <- data.table::data.table(topic_num = c(1L, 3L), unique_genes = c(10L, 5L))

  out <- .topic_report_topic_mds(theta, gene_counts)

  expect_equal(nrow(out), 3L)
  expect_equal(out$topic, c("Topic1", "Topic2", "Topic3"))
  expect_equal(out$unique_genes, c(10L, 0L, 5L))
  expect_true(all(is.finite(out$MDS1)))
  expect_true(all(is.finite(out$MDS2)))
})

test_that("condition report MDS aggregates theta by condition or comparison label", {
  source_benchmark_script(topic_report_script)

  theta <- matrix(
    c(
      0.8, 0.2,
      0.6, 0.4,
      0.1, 0.9
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(theta) <- c(
    "HPAFII_A_vs_B::TF1::Up",
    "HPAFII_A_vs_B::TF2::Up",
    "HPAFII_A_vs_B::TF1::Down"
  )
  colnames(theta) <- c("Topic1", "Topic2")

  out <- .topic_report_group_mds(theta, doc_design = "comparison")

  expect_equal(out$comparison_label, c("HPAFII_A_vs_B::Down", "HPAFII_A_vs_B::Up"))
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Up", Topic1], 0.7)
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Down", Topic2], 0.9)
  expect_true(all(is.finite(out$MDS1)))
  expect_true(all(is.finite(out$MDS2)))
})

test_that("condition report MDS uses authoritative theta group MDS coordinates", {
  source_benchmark_script(topic_report_script)

  theta <- matrix(
    c(
      0.8, 0.2,
      0.6, 0.4,
      0.1, 0.9
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(theta) <- c(
    "HPAFII_A_vs_B::TF1::Up",
    "HPAFII_A_vs_B::TF2::Up",
    "HPAFII_A_vs_B::TF1::Down"
  )
  colnames(theta) <- c("Topic1", "Topic2")

  csv_dir <- tempfile("topic-csv-")
  dir.create(csv_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(
      group_label = c("HPAFII_A_vs_B::Up", "HPAFII_A_vs_B::Down"),
      display_label = c("A vs B Up", "A vs B Down"),
      MDS1 = c(10, 20),
      MDS2 = c(30, 40),
      k = c(10L, 10L),
      setup = c("std_tf_diff_fp_uniq", "std_tf_diff_fp_uniq"),
      doc_design = c("comparison", "comparison"),
      model_label = c("MultiVI", "MultiVI"),
      color_family = c("A", "A"),
      color = c("#111111", "#222222")
    ),
    file.path(csv_dir, "theta_group_mds_points.csv")
  )
  row <- data.table::data.table(
    project = "std_tf_diff_fp_uniq",
    selected_k = 10L,
    combo_id = "doc_tf_fp_uniq_model_mve_K10"
  )

  out <- .topic_report_authoritative_group_mds(theta, row, csv_dir = csv_dir)

  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Up", MDS1], 10)
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Down", MDS2], 40)
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Up", Topic1], 0.7)
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Down", Topic2], 0.9)
  expect_equal(out$color, c("#222222", "#111111"))
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Up", mds_label], "A_vs_B_Up")
  expect_equal(out[comparison_label == "HPAFII_A_vs_B::Down", shape_value], 25L)
})

test_that("condition report waterfall ranks topics within selected condition", {
  source_benchmark_script(topic_report_script)

  group_mds <- data.table::data.table(
    comparison_label = "HPAFII_A_vs_B::Up",
    doc_design = "comparison",
    MDS1 = 0,
    MDS2 = 0,
    n_docs = 2L,
    Topic1 = 0.7,
    Topic2 = 0.3
  )

  out <- .topic_report_group_topic_waterfall(group_mds)

  expect_equal(out$topic, c("Topic1", "Topic2"))
  expect_equal(out$theta_mean, c(0.7, 0.3))
  expect_equal(out$comparison_label, c("HPAFII_A_vs_B::Up", "HPAFII_A_vs_B::Up"))
})

test_that("topic report grid can be restricted to selected K values", {
  source_benchmark_script(topic_report_script)

  grid <- data.table::data.table(
    selected_k = c(5L, 10L, 20L),
    extract_id = c("K5", "K10", "K20")
  )

  out <- .filter_topic_report_grid(grid, selected_k = 10L)

  expect_equal(out$extract_id, "K10")
})

test_that("topic report HTML renders pathway bars and theta tooltips", {
  source_benchmark_script(topic_report_script)

  row <- data.table::data.table(
    project = "std_tf_cond_fp_aggr",
    extract_id = "example_extract",
    selected_k = 10L,
    combo_id = "doc_tf_fp_aggr_model_mve_K10"
  )
  topic_mds <- data.table::data.table(topic = "Topic1", topic_num = 1L, MDS1 = 0, MDS2 = 0, unique_genes = 10L)
  waterfall <- data.table::data.table(topic = "Topic1", topic_num = 1L, comparison_label = "A::Up", theta_mean = 0.5, n_docs = 2L)
  pathways <- data.table::data.table(
    topic = 1L,
    topic_label = "Topic1",
    pathway_key = "DB:Path",
    pathway = "DB: Path",
    padj = 0.01,
    neglog10_padj = 2,
    overlap = "7/40",
    gene_in = 7L,
    gene_total = 40L,
    gene_total_universe = 40L,
    gene_out = 33L,
    genes = "A;B;C;D;E;F;G"
  )
  out_file <- tempfile(fileext = ".html")

  .write_topic_report_html(row, topic_mds, waterfall, pathways, topic_mds[, .(topic_num, unique_genes)], out_file)
  html <- paste(readLines(out_file, warn = FALSE), collapse = "\n")

  expect_match(html, "pathAxis", fixed = TRUE)
  expect_match(html, "pathBarOut", fixed = TRUE)
  expect_match(html, "pathBarIn", fixed = TRUE)
  expect_false(grepl("pathCount", html, fixed = TRUE))
  expect_false(grepl("7 / 40", html, fixed = TRUE))
  expect_false(grepl("pathItem", html, fixed = TRUE))
  expect_false(grepl("+'\nmean theta: '", html, fixed = TRUE))
  expect_false(grepl("+'\ndocs: '", html, fixed = TRUE))
  expect_match(html, "mean document-to-topic probability", fixed = TRUE)
  expect_match(html, "\\ndocs: ", fixed = TRUE)
})

test_that("combined topic report embeds individual report HTML when available", {
  source_benchmark_script(topic_report_script)

  html_dir <- tempfile("topic-html-")
  report_dir <- file.path(html_dir, "topic_reports")
  dir.create(report_dir, recursive = TRUE)
  report_file <- file.path(report_dir, "cond_fp_uniq_MultiVI_K10_topic_report.html")
  writeLines("<!doctype html><html><body>embedded topic report</body></html>", report_file, useBytes = TRUE)

  out <- .write_combined_topic_report_html(html_dir, report_dir = report_dir)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")

  expect_match(html, "topic_reports/cond_fp_uniq_MultiVI_K10_topic_report.html", fixed = TRUE)
  expect_match(html, "frame.srcdoc=hit.html", fixed = TRUE)
  expect_match(html, "embedded topic report", fixed = TRUE)
})

test_that("combined condition report links condition reports by relative path", {
  source_benchmark_script(topic_report_script)

  html_dir <- tempfile("condition-topic-html-")
  report_dir <- file.path(html_dir, "condition_topic_reports")
  dir.create(report_dir, recursive = TRUE)
  report_file <- file.path(report_dir, "cond_fp_uniq_MultiVI_K10_condition_topic_report.html")
  writeLines("<!doctype html><html><body>embedded condition report</body></html>", report_file, useBytes = TRUE)

  out <- .write_combined_condition_topic_report_html(html_dir, report_dir = report_dir)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")

  expect_match(html, "condition_topic_reports/cond_fp_uniq_MultiVI_K10_condition_topic_report.html", fixed = TRUE)
})

test_that("condition-centered report HTML uses group MDS and topic waterfall payloads", {
  source_benchmark_script(topic_report_script)

  row <- data.table::data.table(
    project = "std_tf_diff_fp_uniq",
    extract_id = "example_extract",
    selected_k = 10L,
    combo_id = "doc_tf_fp_uniq_model_mve_K10",
    link_topic_method = "gammafit"
  )
  group_mds <- data.table::data.table(
    comparison_label = c("HPAFII_A_vs_B::Up", "HPAFII_A_vs_B::Down"),
    display_label = c("A vs B Up", "A vs B Down"),
    direction = c("Up", "Down"),
    doc_design = c("comparison", "comparison"),
    MDS1 = c(0, 1),
    MDS2 = c(0, 1),
    n_docs = c(2L, 1L),
    Topic1 = c(0.7, 0.1),
    Topic2 = c(0.3, 0.9)
  )
  waterfall <- .topic_report_group_topic_waterfall(group_mds)
  pathways <- data.table::data.table(
    topic = 1L,
    topic_label = "Topic1",
    pathway_key = "DB:Path",
    pathway = "DB: Path",
    padj = 0.01,
    neglog10_padj = 2,
    overlap = "7/40",
    gene_in = 7L,
    gene_total = 40L,
    gene_total_universe = 40L,
    gene_out = 33L,
    genes = "A;B;C;D;E;F;G"
  )
  out_file <- tempfile(fileext = ".html")

  .write_condition_topic_report_html(row, group_mds, waterfall, pathways, out_file)
  html <- paste(readLines(out_file, warn = FALSE), collapse = "\n")

  expect_match(html, "Condition/Comparison MDS", fixed = TRUE)
  expect_match(html, "Topic Waterfall", fixed = TRUE)
  expect_match(html, "GROUP_MDS_B64", fixed = TRUE)
  expect_match(html, "GROUP_TOPIC_B64", fixed = TRUE)
  expect_match(html, "selectGroup", fixed = TRUE)
  expect_false(grepl("+'\ndocs: '", html, fixed = TRUE))
  expect_match(html, "+'\\ndocs: '", fixed = TRUE)
})
