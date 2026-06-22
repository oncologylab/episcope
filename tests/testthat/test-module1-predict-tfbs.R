test_that("predict_tfbs returns the public Module 1 contract", {
  fixture <- module1_tiny_fixture()
  out_dir <- file.path(tempdir(), paste0("craftgrn-module1-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  result <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    out_dir = out_dir,
    db = "JASPAR2024",
    r_cutoff = 0.8,
    write_outputs = FALSE,
    write_stats = FALSE,
    write_bed = FALSE,
    verbose = FALSE
  )

  expect_named(
    result,
    c(
      "omics_data",
      "high_confidence_footprints",
      "motif_supported_correlations",
      "prediction_stats",
      "predicted_tfbs",
      "prediction_stats_manifest",
      "reports",
      "parameters"
    )
  )
  expect_s3_class(result$prediction_stats, "data.frame")
  expect_s3_class(result$predicted_tfbs, "data.frame")
  expect_false(any(c("best_r", "best_method") %in% names(result$predicted_tfbs)))
  expect_true(all(c("fp_id", "chr", "start", "end", "atac_peak", "tf", "best_r", "best_method", "condition_support") %in% names(result$prediction_stats)))
  expect_false(fixture$expected_excluded_fp %in% result$high_confidence_footprints$fp_id)
  expect_false(fixture$expected_excluded_fp %in% result$prediction_stats$fp_id)
  expect_equal(result$parameters$db, "JASPAR2024")
  expect_equal(result$parameters$r_cutoff, 0.8)
  expect_true(result$parameters$filter_to_canonical_bound)
  expect_null(result$parameters$p_cutoff)
  expect_null(result$parameters$fdr_cutoff)
  expect_named(result$parameters$qc_summary, c("n_fp_input", "n_fp_bound_accessible", "n_expressed_tfs", "n_motif_supported_pairs", "n_canonical_pairs_pass", "n_canonical_bound_fps", "n_prediction_fps", "n_prediction_pairs", "n_prediction_stats", "n_predicted_tfbs"))
  expect_equal(result$parameters$qc_summary$n_prediction_fps, nrow(result$high_confidence_footprints))
  expect_false("TF_E" %in% result$motif_supported_correlations$tf)
})

test_that("predict_tfbs defaults do not write to the working directory", {
  expect_null(formals(predict_tfbs)$out_dir)
  expect_false(formals(predict_tfbs)$write_outputs)
  expect_null(formals(module1_predict_full_tfbs)$out_dir)
  expect_false(formals(module1_predict_full_tfbs)$write_outputs)

  fixture <- module1_tiny_fixture()
  result <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    r_cutoff = 0.8,
    verbose = FALSE
  )

  expect_s3_class(result$predicted_tfbs, "data.frame")
  expect_equal(result$parameters$write_outputs, FALSE)
  expect_length(result$reports, 0L)
})

test_that("predict_tfbs can retain all bound FPs after canonical-bound labeling", {
  fixture <- module1_tiny_fixture()
  strict <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    r_cutoff = 0.8,
    write_outputs = FALSE,
    write_stats = FALSE,
    verbose = FALSE
  )
  retained <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    r_cutoff = 0.8,
    filter_to_canonical_bound = FALSE,
    write_outputs = FALSE,
    write_stats = FALSE,
    verbose = FALSE
  )

  expect_gt(retained$parameters$qc_summary$n_prediction_fps, strict$parameters$qc_summary$n_prediction_fps)
  expect_gt(retained$parameters$qc_summary$n_prediction_pairs, strict$parameters$qc_summary$n_prediction_pairs)
  expect_equal(nrow(retained$high_confidence_footprints), strict$parameters$qc_summary$n_canonical_bound_fps)
  expect_false(fixture$expected_excluded_fp %in% retained$prediction_stats$fp_id)
})

test_that("motif-supported pair construction filters TF subsets before expansion", {
  fixture <- module1_tiny_fixture()
  fixture$omics_data$fp_annotation$tfs[[1L]] <- "TF_A, TF_B"
  pairs <- .module1_motif_supported_pairs(
    omics_data = fixture$omics_data,
    tf_subset = "TF_A"
  )

  expect_gt(nrow(pairs), 0L)
  expect_setequal(unique(pairs$tf), "TF_A")
  expect_false("TF_B" %in% pairs$tf)

  all_pairs <- .module1_motif_supported_pairs(fixture$omics_data)
  first_fp_pairs <- all_pairs[all_pairs$fp_id == fixture$omics_data$fp_annotation$fp_peak[[1L]], ]
  expect_setequal(first_fp_pairs$tf, c("TF_A", "TF_B"))
})

test_that("expressed TF selection is constrained to the TF list", {
  fixture <- module1_tiny_fixture()
  fixture$omics_data$rna_condition <- dplyr::bind_rows(
    fixture$omics_data$rna_condition,
    tibble::tibble(
      ensembl_gene_id = "gene_non_tf",
      HGNC = "HIGH_EXPR_NON_TF",
      cond_a = 100,
      cond_b = 100,
      cond_c = 100,
      cond_d = 100
    )
  )
  fixture$omics_data$rna_expressed <- dplyr::bind_rows(
    fixture$omics_data$rna_expressed,
    tibble::tibble(
      ensembl_gene_id = "gene_non_tf",
      HGNC = "HIGH_EXPR_NON_TF",
      cond_a = 1L,
      cond_b = 1L,
      cond_c = 1L,
      cond_d = 1L
    )
  )

  tfs <- .module1_expressed_tfs(fixture$omics_data)

  expect_false("HIGH_EXPR_NON_TF" %in% tfs)
  expect_true(all(tfs %in% fixture$omics_data$tf_list))
})

test_that("pair correlations match cor.test on complete condition data", {
  fixture <- module1_tiny_fixture()
  pairs <- tibble::tibble(
    fp_id = "chr1:100-140",
    atac_peak = "chr1:80-160",
    motifs = "M1",
    tf = "TF_A"
  )

  stats <- .module1_compute_pair_correlations(
    omics_data = fixture$omics_data,
    pairs = pairs,
    min_non_na = 3L
  )
  cond_cols <- .module1_condition_columns(fixture$omics_data)
  fp_values <- as.numeric(unlist(fixture$omics_data$fp_score_condition_qn[
    fixture$omics_data$fp_score_condition_qn$peak_ID == "chr1:100-140",
    cond_cols
  ], use.names = FALSE))
  tf_values <- as.numeric(unlist(fixture$omics_data$rna_condition[
    fixture$omics_data$rna_condition$HGNC == "TF_A",
    cond_cols
  ], use.names = FALSE))
  pearson <- suppressWarnings(stats::cor.test(fp_values, tf_values, method = "pearson"))
  spearman <- suppressWarnings(stats::cor.test(fp_values, tf_values, method = "spearman", exact = FALSE))

  expect_equal(stats$pearson_r, as.numeric(unname(pearson$estimate)), tolerance = 1e-10)
  expect_equal(stats$pearson_p, as.numeric(pearson$p.value), tolerance = 1e-10)
  expect_equal(stats$spearman_r, as.numeric(unname(spearman$estimate)), tolerance = 1e-10)
  expect_equal(stats$spearman_p, as.numeric(spearman$p.value), tolerance = 1e-10)
})

test_that("streamed prediction links match pairwise prediction links", {
  fixture <- module1_tiny_fixture()
  motif_supported_pairs <- .module1_motif_supported_pairs(fixture$omics_data)
  motif_supported_stats_raw <- .module1_compute_pair_correlations(
    omics_data = fixture$omics_data,
    pairs = motif_supported_pairs,
    min_non_na = 3L
  )
  motif_supported_stats <- .module1_merge_correlation_stats(
    pearson_stats = motif_supported_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
    spearman_stats = motif_supported_stats_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
    r_cutoff = 0.8
  )
  high_conf <- .module1_select_high_confidence_footprints(motif_supported_stats, r_cutoff = 0.8)
  pairs <- .module1_all_prediction_pairs(
    omics_data = fixture$omics_data,
    fp_ids = high_conf$fp_id
  )
  pairwise_raw <- .module1_compute_pair_correlations(
    omics_data = fixture$omics_data,
    pairs = pairs,
    min_non_na = 3L
  )
  pairwise_stats <- .module1_merge_correlation_stats(
    pearson_stats = pairwise_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
    spearman_stats = pairwise_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
    r_cutoff = 0.8
  )
  pairwise_links <- .module1_build_prediction_stats(
    correlation_stats = pairwise_stats,
    high_confidence_footprints = high_conf,
    omics_data = fixture$omics_data
  )

  streamed <- .module1_predict_tfbs_streamed(
    omics_data = fixture$omics_data,
    high_confidence_footprints = high_conf,
    r_cutoff = 0.8,
    min_non_na = 3L
  )

  expect_equal(streamed$prediction_pair_count, nrow(pairs))
  expect_equal(
    streamed$prediction_stats[order(streamed$prediction_stats$fp_id, streamed$prediction_stats$tf), names(pairwise_links)],
    pairwise_links[order(pairwise_links$fp_id, pairwise_links$tf), ]
  )
})

test_that("sparse C++ pair correlations match the R fallback", {
  fixture <- module1_tiny_fixture()
  pairs <- .module1_motif_supported_pairs(fixture$omics_data)
  cpp_stats <- .module1_compute_pair_correlations(
    omics_data = fixture$omics_data,
    pairs = pairs,
    min_non_na = 3L,
    cores = 2L
  )
  r_stats <- with_mocked_bindings(
    .module1_compute_pair_correlations(
      omics_data = fixture$omics_data,
      pairs = pairs,
      min_non_na = 3L,
      cores = 1L
    ),
    .module1_compute_pair_correlations_cpp = function(...) NULL
  )

  expect_equal(cpp_stats$pearson_r, r_stats$pearson_r, tolerance = 1e-10)
  expect_equal(cpp_stats$spearman_r, r_stats$spearman_r, tolerance = 1e-10)
  expect_equal(cpp_stats$pearson_p, r_stats$pearson_p, tolerance = 1e-10)
  expect_equal(cpp_stats$spearman_p, r_stats$spearman_p, tolerance = 1e-10)
})

test_that("predict_tfbs streams large link outputs when return_prediction_stats is disabled", {
  fixture <- module1_tiny_fixture()
  out_dir <- file.path(tempdir(), paste0("craftgrn-module1-stream-", as.integer(stats::runif(1L, 1, 1e9))))
  result <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    out_dir = out_dir,
    r_cutoff = 0.8,
    write_outputs = TRUE,
    write_stats = FALSE,
    output_format = "csv",
    return_prediction_stats = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(result$prediction_stats), 0L)
  expect_s3_class(result$prediction_stats_manifest, "data.frame")
  expect_true(file.exists(result$reports$prediction_stats_manifest))
  expect_true(dir.exists(result$reports$prediction_stats_chunks))
  expect_true(all(file.exists(result$prediction_stats_manifest$path)))
  expect_true(file.exists(result$reports$canonical_correlation_stats))
  expect_true(file.exists(result$reports$qc_summary))
  expect_true(file.exists(result$reports$predicted_tfbs_manifest))
  pred_manifest <- readr::read_csv(result$reports$predicted_tfbs_manifest, show_col_types = FALSE)
  expect_true(all(file.exists(pred_manifest$path)))
  expect_equal(sum(pred_manifest$n_rows), sum(result$prediction_stats_manifest$n_rows))
  qc <- readr::read_csv(result$reports$qc_summary, show_col_types = FALSE)
  expect_true(all(c("metric", "value") %in% names(qc)))
  expect_equal(sum(result$prediction_stats_manifest$n_rows), result$parameters$qc_summary$n_prediction_stats)
  expect_equal(sum(pred_manifest$n_rows), result$parameters$qc_summary$n_predicted_tfbs)
})

test_that("CraftGRN multiomic object uses semantic names", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)

  expect_true(craftgrn:::is_multiomic_object(compact))
  expect_named(compact$matrices, c("fp_score", "fp_bound", "gene_expr", "gene_on", "atac_score", "atac_open"))
  expect_named(compact$features, c("fp", "fp_motif", "atac", "gene"))
  expect_false(any(c("fp_score_condition_qn", "fp_bound_condition", "rna_condition", "rna_expressed", "fp_annotation") %in% names(compact)))
  expect_true(is.logical(compact$matrices$gene_on))
  expect_true(is.logical(compact$matrices$fp_bound))
  expect_silent(craftgrn:::validate_multiomic_object(compact))
})

test_that("predict_tfbs requires CraftGRN multiomic objects", {
  fixture <- module1_tiny_fixture()
  compact <- craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE)

  expect_error(
    predict_tfbs(fixture$omics_data, r_cutoff = 0.8, write_outputs = FALSE, verbose = FALSE),
    "CraftGRN multiomic object"
  )
  compact_result <- predict_tfbs(compact, r_cutoff = 0.8, write_outputs = FALSE, verbose = FALSE)
  expect_s3_class(compact_result$prediction_stats, "data.frame")
  expect_gt(nrow(compact_result$high_confidence_footprints), 0L)
})

test_that("predict_tfbs writes optional BED outputs", {
  fixture <- module1_tiny_fixture()
  out_dir <- file.path(tempdir(), paste0("craftgrn-module1-bed-", as.integer(stats::runif(1L, 1, 1e9))))
  result <- predict_tfbs(
    omics_data = craftgrn:::as_multiomic_object(fixture$omics_data, verbose = FALSE),
    out_dir = out_dir,
    r_cutoff = 0.8,
    write_outputs = TRUE,
    write_bed = TRUE,
    return_prediction_stats = TRUE,
    verbose = FALSE
  )

  expect_true(file.exists(result$reports$high_confidence_footprints_bed))
  expect_true(file.exists(result$reports$prediction_stats_bed))
  expect_gt(length(readLines(result$reports$high_confidence_footprints_bed)), 0L)
  expect_gt(length(readLines(result$reports$prediction_stats_bed)), 0L)
})
