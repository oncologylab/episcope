test_that("Module 1 QC report writes an HTML summary", {
  fixture <- module1_tiny_fixture()
  out_dir <- tempfile("module1-qc-")
  result <- predict_tfbs(
    omics_data = as_multiomic_object(fixture$omics_data, verbose = FALSE),
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
  expect_true(grepl("Motif-Supported Correlations", page, fixed = TRUE))
  expect_true(grepl("Prediction Output Integrity", page, fixed = TRUE))
  expect_true(grepl("Correlation Diagnostics", page, fixed = TRUE))
  expect_true(grepl("Top Predicted FPs", page, fixed = TRUE))
  expect_true(grepl("Warning Checks", page, fixed = TRUE))
  expect_true(grepl("Correctness Checks", page, fixed = TRUE))
  expect_true(grepl("Workflow Funnel", page, fixed = TRUE))
  expect_true(grepl("Canonical-supported predicted FP check: <span class=\"status-pass\">PASS</span>", page, fixed = TRUE))
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
  compact <- as_multiomic_object(omics, verbose = FALSE)
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
  result <- link_tf_targets(
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
  expect_true(grepl("TF-Target Correlation QC", page, fixed = TRUE))
  expect_true(grepl("FP-Target Correlation QC", page, fixed = TRUE))
  expect_true(grepl("Candidate Source QC", page, fixed = TRUE))
  expect_true(grepl("Condition Activity QC", page, fixed = TRUE))
  expect_true(grepl("Warning Checks", page, fixed = TRUE))
  expect_true(grepl("Integrity Checks", page, fixed = TRUE))
  expect_true(grepl("Top TFs In Final Links", page, fixed = TRUE))
})
