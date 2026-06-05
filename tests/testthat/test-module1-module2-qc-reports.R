test_that("Module 1 QC report writes an HTML summary", {
  fixture <- module1_tiny_fixture()
  out_dir <- tempfile("module1-qc-")
  result <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    out_dir = out_dir,
    db = "JASPAR2024",
    r_cutoff = 0.8,
    write_outputs = TRUE,
    write_stats = FALSE,
    write_bed = FALSE,
    write_qc_report = FALSE,
    output_format = "csv",
    verbose = FALSE
  )

  expect_true(file.exists(file.path(out_dir, "cache", "module1_canonical_bound_fps.csv.gz")))
  expect_true(file.exists(file.path(out_dir, "cache", "module1_canonical_support_summary.csv")))

  html <- build_module1_qc_report(
    result,
    output_dir = file.path(out_dir, "reports"),
    scan_predicted_tfbs = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(html))
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_true(grepl("Module 1 QC Report", page, fixed = TRUE))
  expect_true(grepl("Run Parameters", page, fixed = TRUE))
  expect_true(grepl("Input Gates", page, fixed = TRUE))
  expect_true(grepl("Condition QC", page, fixed = TRUE))
  expect_true(grepl("Footprint Alignment And Input QC", page, fixed = TRUE))
  expect_true(grepl("Legacy Summary Coverage", page, fixed = TRUE))
  expect_true(grepl("02_fp_merge_summary.pdf", page, fixed = TRUE))
  expect_true(grepl("Total bound footprints per condition", page, fixed = TRUE))
  expect_true(grepl("Total expressed genes per condition", page, fixed = TRUE))
  expect_true(grepl("Motif-Supported Correlations", page, fixed = TRUE))
  expect_true(grepl("Prediction Output Integrity", page, fixed = TRUE))
  expect_true(grepl("Correlation Diagnostics", page, fixed = TRUE))
  expect_true(grepl("Footprint Motif Complexity", page, fixed = TRUE))
  expect_false(grepl("Top Predicted FPs", page, fixed = TRUE))
  expect_false(grepl("Top FPs by predicted TFBS", page, fixed = TRUE))
  expect_false(grepl("Embedded QC Artifacts", page, fixed = TRUE))
  expect_false(grepl("object class=\"embedded-file\"", page, fixed = TRUE))
  expect_true(grepl("Warning Checks", page, fixed = TRUE))
  expect_true(grepl("Correctness Checks", page, fixed = TRUE))
  expect_true(grepl("Workflow Funnel", page, fixed = TRUE))
  expect_true(grepl("Key Findings", page, fixed = TRUE))
  expect_true(grepl("Module 1 Interpretation", page, fixed = TRUE))
  expect_true(grepl("Recommended Review", page, fixed = TRUE))
  expect_true(grepl("Canonical support", page, fixed = TRUE))
  expect_true(grepl("what the filter did", page, fixed = TRUE))
  expect_true(grepl("Canonical-supported predicted FP check: <span class=\"status-pass\">PASS</span>", page, fixed = TRUE))
  expect_true(grepl("qc-plot-funnel", page, fixed = TRUE))
  expect_true(grepl("qc-plot-density", page, fixed = TRUE))
  expect_true(grepl("qc-plot-scatter", page, fixed = TRUE))
  expect_true(grepl("qc-plot-heatmap", page, fixed = TRUE))
  expect_true(grepl("qc-plot-lollipop", page, fixed = TRUE))
})

test_that("Module 1 QC reads chunked prediction statistics", {
  out_dir <- tempfile("module1-qc-stats-")
  chunk_dir <- file.path(out_dir, "module1_prediction_stats_chunks")
  dir.create(chunk_dir, recursive = TRUE)
  chunk_path <- file.path(chunk_dir, "module1_prediction_stats_chunk_0001.csv")
  stats <- tibble::tibble(
    fp_id = c("fp1", "fp2"),
    tf = c("TF_A", "TF_B"),
    best_r = c(0.9, 0.8)
  )
  readr::write_csv(stats, chunk_path)
  readr::write_csv(
    tibble::tibble(
      chunk_id = 1L,
      path = chunk_path,
      format = "csv",
      n_rows = nrow(stats)
    ),
    file.path(chunk_dir, "module1_prediction_stats_manifest.csv")
  )

  loaded <- .module1_qc_read_prediction_stats(list(), module1_dir = out_dir)

  expect_equal(nrow(loaded), 2L)
  expect_equal(loaded$tf, c("TF_A", "TF_B"))
  expect_true(all(loaded$pass))
})

test_that("QC metric reader accepts one-row integrity tables", {
  integrity <- tibble::tibble(
    n_missing_candidate_source = 0,
    n_missing_distance_to_tss = 0,
    n_bad_coordinate = 0,
    n_duplicate_fp_tf_keys = 0
  )

  expect_equal(.qc_metric_value(integrity, "n_missing_candidate_source"), 0)
  expect_equal(.qc_metric_value(integrity, "n_missing_distance_to_tss"), 0)
  expect_equal(.qc_metric_value(integrity, "n_bad_coordinate"), 0)
  expect_equal(.qc_metric_value(integrity, "n_duplicate_fp_tf_keys"), 0)
})

test_that("Module 2 QC report writes an HTML summary", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1, 9), s2 = c(2, 8), s3 = c(8, 2), s4 = c(9, 1)),
    fp_bound_condition = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1L, 1L), s2 = c(1L, 1L), s3 = c(1L, 1L), s4 = c(1L, 1L)),
    fp_annotation = tibble::tibble(fp_peak = c("chr1:100-140", "chr1:500-540"), atac_peak = c("chr1:90-160", "chr1:490-560"), motifs = c("M_A", "M_B"), tfs = c("TF_A", "TF_B")),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = c(1, 9, 1, 9), s2 = c(2, 8, 2, 8), s3 = c(8, 2, 8, 2), s4 = c(9, 1, 9, 1)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  predicted <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:500-540"),
    chr = "chr1",
    start = c(100L, 500L),
    end = c(140L, 540L),
    atac_peak = c("chr1:90-160", "chr1:490-560"),
    tf = c("TF_A", "TF_B")
  )
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  out_dir <- tempfile("module2-qc-")
  result <- predict_tf_targets(
    compact,
    predicted,
    gene_tss,
    project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)),
    output_dir = out_dir,
    max_distance_bp = 1000,
    n_cores = 1,
    output_format = "csv",
    write_qc_report = FALSE,
    verbose = FALSE
  )

  html <- build_module2_qc_report(
    result,
    multiomic_data = compact,
    output_dir = file.path(out_dir, "reports"),
    scan_large_tables = TRUE,
    validate_integrity = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(html))
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_true(grepl("Module 2 QC Report", page, fixed = TRUE))
  expect_true(grepl("Run Parameters", page, fixed = TRUE))
  expect_true(grepl("Input Handoff", page, fixed = TRUE))
  expect_true(grepl("Condition Context", page, fixed = TRUE))
  expect_true(grepl("TF-Target Correlation QC", page, fixed = TRUE))
  expect_true(grepl("FP-Target Correlation QC", page, fixed = TRUE))
  expect_true(grepl("Candidate Source QC", page, fixed = TRUE))
  expect_true(grepl("Condition Activity QC", page, fixed = TRUE))
  expect_true(grepl("Link Activity Summary", page, fixed = TRUE))
  expect_true(grepl("Total active links per condition", page, fixed = TRUE))
  expect_true(grepl("Active link counts per TF per condition", page, fixed = TRUE))
  expect_true(grepl("Warning Checks", page, fixed = TRUE))
  expect_true(grepl("Integrity Checks", page, fixed = TRUE))
  expect_true(grepl("Top TFs In Final Links", page, fixed = TRUE))
  expect_true(grepl("Top Target Genes In Final Links", page, fixed = TRUE))
  expect_false(grepl("Top Targets And FPs In Final Links", page, fixed = TRUE))
  expect_false(grepl("Top FPs by final links", page, fixed = TRUE))
  expect_false(grepl("Top FP-target correlation metrics", page, fixed = TRUE))
  expect_true(grepl("Key Findings", page, fixed = TRUE))
  expect_true(grepl("Module 2 Interpretation", page, fixed = TRUE))
  expect_true(grepl("Recommended Review", page, fixed = TRUE))
  expect_true(grepl("TF-target prefilter", page, fixed = TRUE))
  expect_true(grepl("FP-target evidence", page, fixed = TRUE))
  expect_true(grepl("qc-plot-flow", page, fixed = TRUE))
  expect_true(grepl("qc-plot-density", page, fixed = TRUE))
  expect_true(grepl("qc-plot-cumulative", page, fixed = TRUE))
  expect_true(grepl("qc-plot-scatter", page, fixed = TRUE))
  expect_true(grepl("qc-plot-heatmap", page, fixed = TRUE))
  expect_true(grepl("qc-plot-lollipop", page, fixed = TRUE))
})

test_that("Module 2 streamed condition activity is reported as skipped", {
  qc_summary <- tibble::tibble(
    metric = c(
      "n_tf_target_pairs_pass",
      "n_fp_target_pairs_pass",
      "n_module2_links",
      "n_active_link_conditions"
    ),
    value = c(2, 2, 2, NA_real_)
  )
  warning_checks <- .module2_qc_warning_checks(
    qc_summary = qc_summary,
    manifest_checks = tibble::tibble(table = "module2_links", file_exists = "yes"),
    candidate_scan = list(integrity = tibble::tibble(
      n_missing_candidate_source = 0,
      n_missing_distance_to_tss = 0
    )),
    link_scan = list(validation = tibble::tibble(
      metric = c("n_links_with_missing_tf_target_pass", "n_links_with_missing_fp_target_pass"),
      value = c(0, 0)
    )),
    tf_corr_scan = list(summary = tibble::tibble()),
    fp_corr_scan = list(summary = tibble::tibble()),
    condition_scan = list(summary = tibble::tibble())
  )

  condition_row <- warning_checks[warning_checks$check == "Condition activity table is available", ]
  expect_equal(condition_row$status, "SKIPPED")
  expect_match(condition_row$detail, "streamed output path")
})
