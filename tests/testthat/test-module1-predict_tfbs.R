test_that("predict_tfbs returns the public Module 1 contract", {
  fixture <- module1_tiny_fixture()
  out_dir <- file.path(tempdir(), paste0("craftgrn-module1-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  result <- predict_tfbs(
    omics_data = fixture$omics_data,
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
      "tfbs_links",
      "tfbs_stats",
      "reports",
      "parameters"
    )
  )
  expect_s3_class(result$tfbs_links, "data.frame")
  expect_true(all(c("fp_id", "chr", "start", "end", "atac_peak", "tf", "best_r", "best_method", "condition_support") %in% names(result$tfbs_links)))
  expect_false(fixture$expected_excluded_fp %in% result$high_confidence_footprints$fp_id)
  expect_false(fixture$expected_excluded_fp %in% result$tfbs_links$fp_id)
  expect_null(result$tfbs_stats)
  expect_equal(result$parameters$db, "JASPAR2024")
  expect_equal(result$parameters$r_cutoff, 0.8)
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
  motif_supported_stats <- .module1_merge_tfbs_stats(
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
  pairwise_stats <- .module1_merge_tfbs_stats(
    pearson_stats = pairwise_raw[, c("fp_id", "atac_peak", "tf", "motifs", "pearson_r", "pearson_p", "pearson_p_adj"), drop = FALSE],
    spearman_stats = pairwise_raw[, c("fp_id", "atac_peak", "tf", "motifs", "spearman_r", "spearman_p", "spearman_p_adj"), drop = FALSE],
    r_cutoff = 0.8
  )
  pairwise_links <- .module1_build_tfbs_links(
    tfbs_stats = pairwise_stats,
    high_confidence_footprints = high_conf,
    omics_data = fixture$omics_data
  )

  streamed <- .module1_predict_links_streamed(
    omics_data = fixture$omics_data,
    high_confidence_footprints = high_conf,
    r_cutoff = 0.8,
    min_non_na = 3L
  )

  expect_equal(streamed$prediction_pair_count, nrow(pairs))
  expect_equal(
    streamed$tfbs_links[order(streamed$tfbs_links$fp_id, streamed$tfbs_links$tf), names(pairwise_links)],
    pairwise_links[order(pairwise_links$fp_id, pairwise_links$tf), ]
  )
})
