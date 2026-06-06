test_that("Module 3 report API exposes only publication-facing helpers", {
  exports <- getNamespaceExports("craftgrn")
  expect_true("build_module3_qc_report" %in% exports)
  expect_true("report_differential_pathways" %in% exports)
  expect_true("report_master_tfs" %in% exports)
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

test_that("report_differential_pathways writes status tables for low-signal inputs", {
  filtered_dir <- file.path(tempdir(), paste0("module3-pathway-links-", sample.int(1e7, 1L)))
  output_dir <- file.path(tempdir(), paste0("module3-pathway-reports-", sample.int(1e7, 1L)))
  dir.create(filtered_dir, recursive = TRUE)
  data.table::fwrite(
    data.table::data.table(tf = "TF1", gene_key = "Gene1", peak_id = "peak1"),
    file.path(filtered_dir, "CondA_vs_CondB_filtered_links_up.csv")
  )

  res <- report_differential_pathways(
    filtered_dir = filtered_dir,
    output_dir = output_dir,
    make_dotplot = FALSE,
    min_genes = 5L,
    verbose = FALSE
  )

  expect_true(dir.exists(file.path(output_dir, "pathway_tables")))
  expect_true(file.exists(file.path(output_dir, "pathway_tables", "pathway_enrichment_manifest.csv")))
  expect_true(file.exists(res$selected_terms_file))
  manifest <- data.table::fread(file.path(output_dir, "pathway_tables", "pathway_enrichment_manifest.csv"))
  expect_equal(manifest$status[[1L]], "too_few_genes")
})

test_that("report_master_tfs writes report outputs with supporting tables separated", {
  testthat::skip_if_not_installed("ggrepel")
  testthat::skip_if_not_installed("scales")
  filtered_dir <- file.path(tempdir(), paste0("module3-master-links-", sample.int(1e7, 1L)))
  output_dir <- file.path(tempdir(), paste0("module3-master-reports-", sample.int(1e7, 1L)))
  dir.create(filtered_dir, recursive = TRUE)
  links <- data.table::data.table(
    tf = c("TF1", "TF1", "TF2", "TF3"),
    gene_key = c("GeneA", "GeneB", "TF1", "GeneC"),
    peak_id = paste0("peak", 1:4),
    delta_fp_score = c(1.2, 0.8, -1.1, 1.5),
    tf_expr_max = c(20, 20, 15, 10),
    log2FC_tf_expr = c(1.1, 1.1, -0.7, 0.4)
  )
  data.table::fwrite(links, file.path(filtered_dir, "CondA_vs_CondB_filtered_links_up.csv"))

  manifest <- report_master_tfs(
    filtered_dir = filtered_dir,
    output_dir = output_dir,
    waterfall_min_abs_net = 1L,
    verbose = FALSE
  )

  expect_true(nrow(manifest) >= 1L)
  expect_true(dir.exists(file.path(output_dir, "master_tf_tables")))
  expect_true(file.exists(file.path(output_dir, "master_tf_tables", "CondA_vs_CondB_master_tf_summary.csv")))
  expect_true(file.exists(file.path(output_dir, "CondA_vs_CondB_master_tf_summary.pdf")))
})
