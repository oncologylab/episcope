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
    project_config = list(ref_genome = "hg38", threshold_fp_tf_corr_r = 0.7),
    project_date = "2026-07-10",
    verbose = FALSE
  )

  expect_true(file.exists(html))
  page <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_true(grepl("Module 1 QC Report (2026-07-10)", page, fixed = TRUE))
  expect_true(grepl("1. Run Summary", page, fixed = TRUE))
  expect_true(grepl("Input parameters", page, fixed = TRUE))
  expect_true(grepl("# Conditions", page, fixed = TRUE))
  expect_true(grepl("Raw ATAC peaks", page, fixed = TRUE))
  expect_true(grepl("Raw footprints", page, fixed = TRUE))
  expect_true(grepl("Filtered footprints", page, fixed = TRUE))
  expect_true(grepl("Filtered footprints are the canonical-bound footprint set", page, fixed = TRUE))
  expect_true(grepl("Predicted unique TFBS", page, fixed = TRUE))
  expect_true(grepl("TFs with predicted binding", page, fixed = TRUE))
  expect_true(grepl("command override", page, fixed = TRUE))
  expect_true(grepl("2. Per-Condition QC", page, fixed = TRUE))
  expect_true(grepl("Footprint score distribution per condition", page, fixed = TRUE))
  expect_true(grepl("ATAC master table not provided", page, fixed = TRUE))
  expect_true(grepl("Gene expression distribution per condition", page, fixed = TRUE))
  expect_true(grepl("RNA PCA", page, fixed = TRUE))
  expect_false(grepl("ATAC PCA", page, fixed = TRUE))
  expect_true(grepl("Total filtered footprints detected per condition", page, fixed = TRUE))
  expect_true(grepl("Total predicted TFBS per condition", page, fixed = TRUE))
  expect_true(grepl("3. Correlation Summary", page, fixed = TRUE))
  expect_true(grepl("Selection: best(Pearson R, Spearman R) &gt;= 0.8", page, fixed = TRUE))
  expect_true(grepl("Gray bars are all evaluated pairs", page, fixed = TRUE))
  expect_true(grepl("Best R", page, fixed = TRUE))
  expect_true(grepl("4. Predicted Binding Sites per TF", page, fixed = TRUE))
  expect_true(grepl("5. Co-binding Summary", page, fixed = TRUE))
  expect_true(grepl("6. Top Predicted TFs at Motif-Containing Footprints", page, fixed = TRUE))
  expect_true(grepl("Technical Appendix", page, fixed = TRUE))
  expect_false(grepl("Top Predicted FPs", page, fixed = TRUE))
  expect_false(grepl("Top FPs by predicted TFBS", page, fixed = TRUE))
  expect_false(grepl("Embedded QC Artifacts", page, fixed = TRUE))
  expect_false(grepl("object class=\"embedded-file\"", page, fixed = TRUE))
  expect_true(grepl("Integrity and warning checks", page, fixed = TRUE))
  expect_false(grepl("Module 1 processing funnel", page, fixed = TRUE))
  expect_false(grepl("Input FPs", page, fixed = TRUE))
  expect_true(grepl("Module 1 Interpretation", page, fixed = TRUE))
  expect_true(grepl("Review Before Downstream Use", page, fixed = TRUE))
  expect_true(grepl("Canonical support", page, fixed = TRUE))
  expect_true(grepl("what the filter did", page, fixed = TRUE))
  expect_true(grepl("Canonical Support Check Detail", page, fixed = TRUE))
  expect_true(grepl("report-row", page, fixed = TRUE))
  expect_true(grepl("violin-plot", page, fixed = TRUE))
  expect_true(grepl("pca-plot", page, fixed = TRUE))
  expect_true(grepl("correlation-grid", page, fixed = TRUE))
  expect_true(grepl("binding-overall", page, fixed = TRUE))
  expect_true(grepl("cobind-chart", page, fixed = TRUE))
  expect_true(grepl("motif-chart", page, fixed = TRUE))
  expect_true(grepl("vertical-violin-plot", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"per-condition\"", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"distribution\"", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"pca\"", page, fixed = TRUE))
  expect_true(grepl("width:min(98vw,2880px)", page, fixed = TRUE))
  expect_true(grepl("option value=\"best\" selected", page, fixed = TRUE))
  expect_true(grepl("Colored retained-pair overlays are shown only for Best R", page, fixed = TRUE))
  expect_true(grepl("option selected>10</option>", page, fixed = TRUE))
  expect_true(grepl("@media(min-width:1600px)", page, fixed = TRUE))
  expect_true(grepl("How to read this section", page, fixed = TRUE))
  expect_true(grepl("Canonical motif TF reference", page, fixed = TRUE))
  expect_true(grepl("cobind-apply", page, fixed = TRUE))
  expect_true(grepl("qc-plot-density", page, fixed = TRUE))
  expect_true(grepl("qc-plot-lollipop", page, fixed = TRUE))
  html_ids <- regmatches(page, gregexpr("id=\"[^\"]+\"", page, perl = TRUE))[[1L]]
  expect_equal(length(html_ids), length(unique(html_ids)))
  expect_false(file.exists(file.path(out_dir, "reports", "module1_tfbs_explorer.html")))
  expect_true(file.exists(file.path(out_dir, "cache", "module1_qc_analysis.rds")))
})

test_that("Module 1 violin stages place quantile normalization last", {
  stages <- craftgrn:::.module1_qc_violin_stage_order(c(
    "quantile_normalized", "raw", "canonical_filtered"
  ))

  expect_equal(stages, c("raw", "canonical_filtered", "quantile_normalized"))
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

test_that("Module 1 QC recovers legacy raw footprint counts from compact caches", {
  cache_dir <- tempfile("module1-legacy-fp-cache-")
  dir.create(cache_dir, recursive = TRUE)
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("p1", "p2"),
      source_fp_peaks = c("s1;s2", "s2;s3"),
      n_source_fp_peaks = c(2L, 2L)
    ),
    file.path(cache_dir, "fp_sites_TEST.csv")
  )
  compact <- craftgrn:::as_multiomic_object(module1_tiny_fixture()$omics_data, verbose = FALSE)
  recovered <- craftgrn:::.module1_qc_raw_footprint_recovery(
    compact,
    project_config = list(fp_cache_dir = cache_dir, db = "TEST")
  )
  expect_equal(recovered$count, 3)

  cards <- craftgrn:::.module1_qc_run_cards(
    tibble::tibble(metric = character(), value = numeric()),
    compact,
    list(tf_summary_all = tibble::tibble()),
    predicted_rows = 1,
    legacy_raw_unavailable = TRUE
  )
  expect_equal(cards$value[cards$label == "Raw footprints"], "Unavailable for legacy run")
})

test_that("Module 1 QC reads singular provenance from Parquet compact caches", {
  skip_if_not_installed("arrow")
  cache_dir <- tempfile("module1-singular-fp-cache-")
  dir.create(cache_dir, recursive = TRUE)
  arrow::write_parquet(
    tibble::tibble(
      peak_ID = c("p1", "p2", "p3"),
      source_fp_peak = c("s1", "s2", "s2")
    ),
    file.path(cache_dir, "fp_sites_TEST.parquet")
  )
  compact <- craftgrn:::as_multiomic_object(module1_tiny_fixture()$omics_data, verbose = FALSE)
  recovered <- craftgrn:::.module1_qc_raw_footprint_recovery(
    compact,
    project_config = list(fp_cache_dir = cache_dir, db = "TEST")
  )
  expect_equal(recovered$count, 2)
})

test_that("Module 1 QC follows source projects and reads legacy footprint maps", {
  project_root <- tempfile("module1-legacy-fp-project-")
  source_root <- tempfile("module1-legacy-fp-source-")
  dir.create(file.path(project_root, "predict_tf_binding_sites"), recursive = TRUE)
  dir.create(file.path(source_root, "cache"), recursive = TRUE)
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("p1", "p2", "p3", "p4"),
      fp_peak_bak = c("s1", "s2", "s2", "s3")
    ),
    file.path(source_root, "cache", "fp_id_map_TEST.csv")
  )
  config_path <- file.path(project_root, "project.yaml")
  yaml::write_yaml(list(base_dir = project_root, source_project = source_root, db = "TEST"), config_path)
  compact <- craftgrn:::as_multiomic_object(module1_tiny_fixture()$omics_data, verbose = FALSE)

  recovered <- craftgrn:::.module1_qc_raw_footprint_recovery(
    compact,
    module1_dir = file.path(project_root, "predict_tf_binding_sites"),
    project_config = config_path
  )

  expect_equal(recovered$count, 3)
  expect_match(recovered$path, "fp_id_map_TEST[.]csv$")
})

test_that("Module 1 QC recovers and caches legacy raw score distributions", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  compact$qc$fp_score_distributions <- tibble::tibble()
  project_root <- tempfile("module1-distribution-project-")
  source_root <- tempfile("module1-distribution-source-")
  module1_dir <- file.path(project_root, "predict_tf_binding_sites")
  source_module1 <- file.path(source_root, "predict_tf_binding_sites")
  dir.create(module1_dir, recursive = TRUE)
  dir.create(source_module1, recursive = TRUE)
  raw <- fixture$omics_data$fp_score_condition_qn
  raw[, -1] <- lapply(raw[, -1, drop = FALSE], function(x) x * 3 + 1)
  readr::write_csv(raw, file.path(source_module1, "02_fp_score_raw_TEST.csv"))
  config_path <- file.path(project_root, "project.yaml")
  yaml::write_yaml(list(base_dir = project_root, source_project = source_root, db = "TEST"), config_path)

  distributions <- craftgrn:::.module1_qc_distribution_table(
    compact,
    module1_dir = module1_dir,
    project_config = config_path,
    verbose = FALSE
  )

  expect_setequal(unique(distributions$stage), c("raw", "quantile_normalized"))
  expect_equal(sort(unique(distributions$condition)), sort(colnames(compact$matrices$fp_score)))
  expect_true(file.exists(file.path(module1_dir, "module1_fp_score_distributions.csv")))
  cached <- craftgrn:::.module1_qc_read_distribution_cache(compact, module1_dir)
  expect_equal(nrow(cached), 2L * 4L * 101L)
})

test_that("Module 1 QC rejects mismatched raw score footprints", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  compact$qc$fp_score_distributions <- tibble::tibble()
  project_root <- tempfile("module1-distribution-mismatch-")
  module1_dir <- file.path(project_root, "predict_tf_binding_sites")
  dir.create(module1_dir, recursive = TRUE)
  raw <- fixture$omics_data$fp_score_condition_qn
  raw$peak_ID[[1L]] <- "wrong-footprint"
  readr::write_csv(raw, file.path(module1_dir, "02_fp_score_raw_TEST.csv"))

  distributions <- craftgrn:::.module1_qc_distribution_table(
    compact,
    module1_dir = module1_dir,
    project_config = list(base_dir = project_root, db = "TEST"),
    verbose = FALSE
  )

  expect_equal(unique(distributions$stage), "quantile_normalized")
  expect_false(file.exists(file.path(module1_dir, "module1_fp_score_distributions.csv")))
})

test_that("Module 1 QC recovers validated legacy canonical footprints", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)
  prepared <- craftgrn:::.module1_prepare_predict_omics(compact, verbose = FALSE)
  module1_dir <- tempfile("module1-canonical-recovery-")
  dir.create(file.path(module1_dir, "cache"), recursive = TRUE)
  used <- c(1L, 3L, 5L)
  explorer <- list(
    fingerprint = list(
      n_sites = nrow(compact$matrices$fp_score),
      fp_ids_digest = digest::digest(rownames(compact$matrices$fp_score), algo = "xxhash64")
    ),
    tf_names = craftgrn:::.module1_expressed_tfs(prepared),
    used_site_bits = craftgrn:::.module1_pack_site_indices(used, nrow(compact$matrices$fp_score))
  )
  saveRDS(explorer, file.path(module1_dir, "cache", "module1_qc_analysis.rds"))

  recovered <- craftgrn:::.module1_qc_recover_canonical_fp_cache(
    list(),
    module1_dir,
    compact,
    project_config = list(filter_to_canonical_bound = TRUE),
    verbose = FALSE
  )

  expect_equal(recovered$fp_id, rownames(compact$matrices$fp_score)[used])
  expect_true(file.exists(file.path(module1_dir, "cache", "module1_canonical_bound_fps.csv.gz")))
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
  expect_true(grepl("1. Run Summary", page, fixed = TRUE))
  expect_true(grepl("2. Per-Condition QC", page, fixed = TRUE))
  expect_true(grepl("3. Correlation Evidence", page, fixed = TRUE))
  expect_true(grepl("4. Candidate Evidence", page, fixed = TRUE))
  expect_true(grepl("5. Final Regulatory Links", page, fixed = TRUE))
  expect_true(grepl("Technical Appendix", page, fixed = TRUE))
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
  expect_true(grepl("class=\"qc-nav\"", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"module2-condition\"", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"module2-correlation\"", page, fixed = TRUE))
  expect_true(grepl("data-qc-group=\"module2-final-links\"", page, fixed = TRUE))
  expect_true(grepl("role=\"tab\"", page, fixed = TRUE))
  expect_true(grepl("role=\"tabpanel\"", page, fixed = TRUE))
  expect_true(grepl("aria-controls=", page, fixed = TRUE))
  expect_true(grepl("tabindex=\"0\"", page, fixed = TRUE))
  expect_true(grepl("tabindex=\"-1\"", page, fixed = TRUE))
  expect_true(grepl("ArrowLeft", page, fixed = TRUE))
  expect_true(grepl("How to read this section", page, fixed = TRUE))
  expect_true(grepl("status-badge", page, fixed = TRUE))
  expect_true(grepl("distinct row universes", page, fixed = TRUE))
  section_tags <- regmatches(page, gregexpr("<section id=", page, fixed = TRUE))[[1L]]
  expect_length(section_tags, 6L)
  html_ids <- regmatches(page, gregexpr("id=\"[^\"]+\"", page, perl = TRUE))[[1L]]
  expect_equal(length(html_ids), length(unique(html_ids)))
})

test_that("QC metric heatmaps preserve readable metric labels", {
  values <- tibble::tibble(
    source = "tss_window",
    `Target genes` = 12,
    `Median |distance|` = 2500
  )

  svg <- craftgrn:::.qc_metric_heatmap_svg(
    values,
    row_col = "source",
    value_cols = c("Target genes", "Median |distance|"),
    title = "Candidate evidence"
  )

  expect_match(svg, "Target genes", fixed = TRUE)
  expect_match(svg, "Median |distance|", fixed = TRUE)
  expect_match(svg, ">12</text>", fixed = TRUE)
  expect_match(svg, ">2,500</text>", fixed = TRUE)
  expect_false(grepl("fill=\"\"", svg, fixed = TRUE))
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
