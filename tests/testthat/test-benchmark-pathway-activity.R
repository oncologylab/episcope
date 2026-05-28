source_tf_regulated_cluster_browser_script <- function() {
  script <- normalizePath(test_path("../..", "dev/benchmark/02_build_tf_regulated_cluster_html_browsers.R"), mustWork = FALSE)
  testthat::skip_if_not(file.exists(script), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(script))))
  on.exit(setwd(old_wd), add = TRUE)
  source(script, local = parent.frame())
}

test_that("pathway activity scores use matched genes as expression signatures", {
  source_tf_regulated_cluster_browser_script()

  enrich_dt <- data.table::data.table(
    major_cluster = "G01",
    database = "Reactome_2022",
    term = "Example pathway",
    overlap_genes = "A;B;MISSING",
    adjusted_p_value = 0.01
  )
  rna_dt <- data.table::data.table(
    HGNC = c("A", "B", "C"),
    Cond1 = c(7, 9, 5),
    Cond2 = c(3, 1, 5)
  )

  activity <- .score_pathway_activity_from_expression(
    pathway_dt = enrich_dt,
    rna_dt = rna_dt,
    min_genes = 2L,
    pseudocount = 1,
    methods = "mean_gene_zscore"
  )

  expect_equal(nrow(activity), 2L)
  expect_equal(unique(activity$n_genes_matched), 2L)
  expect_true(all(activity$score_method == "mean_gene_zscore"))
  expect_true(all(activity$score_label == "Aggregated activity score (relative gene expression)"))
  expect_gt(activity[condition == "Cond1", activity_score], 0)
  expect_lt(activity[condition == "Cond2", activity_score], 0)
})

test_that("pathway activity scores support matched-control relative expression scores", {
  source_tf_regulated_cluster_browser_script()

  enrich_dt <- data.table::data.table(
    major_cluster = "G01",
    database = "Reactome_2022",
    term = "Example pathway",
    overlap_genes = "A;B",
    adjusted_p_value = 0.01
  )
  rna_dt <- data.table::data.table(
    HGNC = c("A", "B", "C", "D", "E", "F"),
    Cond1 = c(12, 10, 5, 4, 3, 2),
    Cond2 = c(2, 3, 5, 4, 3, 2)
  )

  activity <- .score_pathway_activity_from_expression(
    pathway_dt = enrich_dt,
    rna_dt = rna_dt,
    min_genes = 2L,
    pseudocount = 1,
    methods = c("mean_gene_zscore", "expression_matched_module_score"),
    module_score_nbin = 3L,
    module_score_ctrl = 2L
  )

  expect_setequal(
    unique(activity$score_method),
    c("mean_gene_zscore", "expression_matched_module_score")
  )
  expect_true("Aggregated activity score (relative gene expression with expression-matched background)" %in% activity$score_label)
  expect_true(all(activity[score_method == "expression_matched_module_score", n_control_genes] > 0))
  expect_gt(
    activity[score_method == "expression_matched_module_score" & condition == "Cond1", activity_score],
    activity[score_method == "expression_matched_module_score" & condition == "Cond2", activity_score]
  )
})

test_that("pathway activity differentials follow comparison direction", {
  source_tf_regulated_cluster_browser_script()

  activity_dt <- data.table::data.table(
    major_cluster = "G01",
    database = "Reactome_2022",
    term = "Example pathway",
    pathway_key = "Reactome_2022\tExample pathway",
    n_genes_matched = 2L,
    score_method = "mean_gene_zscore",
    condition = c("Cond1", "Cond2"),
    activity_score = c(1.25, -0.75)
  )
  comparison_dt <- data.table::data.table(
    cond1_label = "Cond1",
    cond2_label = "Cond2"
  )

  diff_dt <- .compute_pathway_activity_differential(activity_dt, comparison_dt)

  expect_equal(nrow(diff_dt), 1L)
  expect_equal(diff_dt$delta_score, 2)
  expect_equal(diff_dt$comparison, "Cond1_vs_Cond2")
})

test_that("tf regulated cluster result labels default to all genes", {
  source_tf_regulated_cluster_browser_script()

  expect_equal(.result_label(100), "allGenes")
  expect_equal(.result_label(100, gene_mode = "top_n"), "top100union")
})

test_that("pathway activity outputs write expression-matched delta dotplot pdf", {
  source_tf_regulated_cluster_browser_script()

  out_dir <- tempfile("pathway_activity")
  csv_dir <- file.path(out_dir, "csv")
  pdf_dir <- file.path(out_dir, "pdf")
  dir.create(csv_dir, recursive = TRUE)
  dir.create(pdf_dir, recursive = TRUE)

  pathway_dt <- data.table::data.table(
    major_cluster = "G01",
    database = "Reactome_2022",
    term = "Example pathway",
    term_clean = "Example pathway",
    pathway_key = "Reactome_2022\tExample pathway",
    overlap_genes = "A;B",
    adjusted_p_value = 0.01
  )
  rna_path <- file.path(out_dir, "rna.csv")
  comparison_path <- file.path(out_dir, "comparisons.csv")
  data.table::fwrite(
    data.table::data.table(
      HGNC = c("A", "B", "C", "D", "E", "F"),
      Cond1 = c(12, 10, 5, 4, 3, 2),
      Cond2 = c(2, 3, 5, 4, 3, 2)
    ),
    rna_path
  )
  data.table::fwrite(
    data.table::data.table(cond1_label = "Cond1", cond2_label = "Cond2"),
    comparison_path
  )

  .write_pathway_activity_outputs(
    pathway_dt = pathway_dt,
    tag = "test",
    result_label = "allGenes",
    file_suffix = "_K03",
    csv_dir = csv_dir,
    pdf_dir = pdf_dir,
    rna_path = rna_path,
    comparison_path = comparison_path,
    min_genes = 2L,
    methods = c("mean_gene_zscore", "expression_matched_module_score"),
    module_score_nbin = 3L,
    module_score_ctrl = 2L
  )

  expect_true(file.exists(file.path(csv_dir, "test_allGenes_K03_pathway_activity_differential.csv")))
  expect_true(file.exists(file.path(pdf_dir, "test_allGenes_K03_pathway_activity_heatmap.pdf")))
  expect_true(file.exists(file.path(pdf_dir, "test_allGenes_K03_pathway_activity_delta_dotplot.pdf")))
})

test_that("pathway activity outputs only score significant pathways", {
  source_tf_regulated_cluster_browser_script()

  pathway_dt <- data.table::data.table(
    major_cluster = c("G01", "G02"),
    database = "Reactome_2022",
    term = c("Significant pathway", "Non-significant pathway"),
    term_clean = c("Significant pathway", "Non-significant pathway"),
    pathway_key = c(
      "Reactome_2022\tSignificant pathway",
      "Reactome_2022\tNon-significant pathway"
    ),
    overlap_genes = c("A;B", "C;D"),
    adjusted_p_value = c(0.01, 0.25)
  )
  rna_dt <- data.table::data.table(
    HGNC = c("A", "B", "C", "D", "E", "F"),
    Cond1 = c(12, 10, 11, 9, 3, 2),
    Cond2 = c(2, 3, 3, 4, 3, 2)
  )

  activity <- .score_pathway_activity_from_expression(
    pathway_dt = pathway_dt,
    rna_dt = rna_dt,
    min_genes = 2L,
    methods = "mean_gene_zscore",
    significance_cutoff = 0.05
  )

  expect_setequal(unique(activity$term), "Significant pathway")
  expect_false("Non-significant pathway" %in% activity$term)
})

test_that("pathway activity subnetwork edges use expression-matched score and matched pathway genes", {
  source_tf_regulated_cluster_browser_script()

  diff_dt <- data.table::data.table(
    major_cluster = c("G01", "G02"),
    database = "Reactome_2022",
    term = c("Pathway A", "Pathway B"),
    term_clean = c("Pathway A", "Pathway B"),
    pathway_key = c("Reactome_2022\tPathway A", "Reactome_2022\tPathway B"),
    matched_genes = c("GENE1;GENE2", "GENE3;GENE4"),
    adjusted_p_value = c(0.01, 0.02),
    score_method = c("expression_matched_module_score", "mean_gene_zscore"),
    score_label = c(
      "Aggregated activity score (relative gene expression with expression-matched background)",
      "Aggregated activity score (relative gene expression)"
    ),
    comparison = "Cond1_vs_Cond2",
    cond1 = "Cond1",
    cond2 = "Cond2",
    score_cond1 = c(1, 1),
    score_cond2 = c(0, 0),
    delta_score = c(1, 1),
    abs_delta_score = c(1, 1)
  )
  tf_gene_dt <- data.table::data.table(
    tf = c("TF1", "TF2", "TF3"),
    gene_norm = c("GENE1", "GENE2", "GENE3"),
    gene_symbol = c("GENE1", "GENE2", "GENE3"),
    fp_target_rna_r = c(0.9, 0.8, 0.7),
    tf_expression_target_r = c(0.6, 0.5, 0.4),
    n_supporting_links = c(3L, 2L, 1L),
    n_supporting_conditions = c(2L, 2L, 1L),
    best_peak_ID = c("p1", "p2", "p3")
  )
  rna_dt <- data.table::data.table(
    HGNC = c("TF1", "TF2", "TF3", "GENE1", "GENE2", "GENE3"),
    Cond1 = c(8, 4, 2, 16, 4, 8),
    Cond2 = c(2, 4, 8, 4, 4, 2)
  )

  edge_dt <- .build_pathway_activity_tf_gene_subnetwork_edges(
    diff_dt = diff_dt,
    tf_gene_dt = tf_gene_dt,
    rna_dt = rna_dt,
    top_n = 5L,
    max_edges_per_pathway = 10L
  )

  expect_setequal(unique(edge_dt$score_method), "expression_matched_module_score")
  expect_setequal(edge_dt$gene_norm, c("GENE1", "GENE2"))
  expect_false("GENE3" %in% edge_dt$gene_norm)
  expect_true(all(is.finite(edge_dt$gene_delta_log2)))
  expect_true(all(is.finite(edge_dt$tf_delta_log2)))
})

test_that("pathway activity condition subnetworks use high condition activity and relative expression", {
  source_tf_regulated_cluster_browser_script()

  activity_dt <- data.table::data.table(
    major_cluster = c("G01", "G02", "G01"),
    database = "Reactome_2022",
    term = c("Pathway A", "Pathway B", "Pathway C"),
    term_clean = c("Pathway A", "Pathway B", "Pathway C"),
    pathway_key = c("Reactome_2022\tPathway A", "Reactome_2022\tPathway B", "Reactome_2022\tPathway C"),
    matched_genes = c("GENE1;GENE2", "GENE3;GENE4", "GENE5;GENE6"),
    adjusted_p_value = c(0.01, 0.02, 0.03),
    score_method = c(
      "expression_matched_module_score",
      "expression_matched_module_score",
      "mean_gene_zscore"
    ),
    score_label = c(
      "Aggregated activity score (relative gene expression with expression-matched background)",
      "Aggregated activity score (relative gene expression with expression-matched background)",
      "Aggregated activity score (relative gene expression)"
    ),
    condition = "Cond1",
    activity_score = c(2, 1, 9),
    n_genes_matched = c(2L, 2L, 2L)
  )
  tf_gene_dt <- data.table::data.table(
    tf = c("TF1", "TF2", "TF3"),
    gene_norm = c("GENE1", "GENE2", "GENE5"),
    gene_symbol = c("GENE1", "GENE2", "GENE5"),
    fp_target_rna_r = c(0.9, 0.8, 0.7),
    tf_expression_target_r = c(0.6, 0.5, 0.4),
    n_supporting_links = c(3L, 2L, 1L),
    n_supporting_conditions = c(2L, 2L, 1L),
    best_peak_ID = c("p1", "p2", "p3")
  )
  rna_dt <- data.table::data.table(
    HGNC = c("TF1", "TF2", "TF3", "GENE1", "GENE2", "GENE5"),
    Cond1 = c(8, 4, 2, 16, 4, 8),
    Cond2 = c(2, 4, 8, 4, 4, 2)
  )

  edge_dt <- .build_pathway_activity_condition_tf_gene_subnetwork_edges(
    activity_dt = activity_dt,
    tf_gene_dt = tf_gene_dt,
    rna_dt = rna_dt,
    top_n = 2L,
    max_edges_per_pathway = 10L
  )

  expect_setequal(unique(edge_dt$score_method), "expression_matched_module_score")
  expect_setequal(edge_dt$pathway_key, "Reactome_2022\tPathway A")
  expect_setequal(edge_dt$cell_line, "Cond1")
  expect_setequal(edge_dt$gene_norm, c("GENE1", "GENE2"))
  expect_true(all(is.finite(edge_dt$gene_relative_expression)))
  expect_true(all(is.finite(edge_dt$tf_relative_expression)))
})

test_that("module2 composite TF connectivity uses direct TF targets and shared targets", {
  source_tf_regulated_cluster_browser_script()

  gene_dt <- data.table::data.table(
    tf = c("TF1", "TF1", "TF1", "TF2", "TF2", "TF2"),
    gene_norm = c("TF2", "GENEA", "GENEB", "TF1", "GENEA", "GENEC"),
    gene_symbol = c("TF2", "GENEA", "GENEB", "TF1", "GENEA", "GENEC")
  )

  mat <- .build_tf_tf_composite_connectivity_matrix(gene_dt, min_score = 0)

  expect_equal(rownames(mat), c("TF1", "TF2"))
  expect_equal(colnames(mat), c("TF1", "TF2"))
  expect_equal(mat["TF1", "TF2"], log2(2 + 0.5 * 1 + 1))
  expect_equal(mat["TF2", "TF1"], log2(2 + 0.5 * 1 + 1))
})
