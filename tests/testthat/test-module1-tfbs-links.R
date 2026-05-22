test_that("Module 1 helpers select high-confidence footprints from motif-supported correlations", {
  motif_supported_stats <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:200-240", "chr2:100-150"),
    atac_peak = c("chr1:90-160", "chr1:180-260", "chr2:80-170"),
    tf = c("TF_A", "TF_B", "TF_D"),
    pearson_r = c(0.92, 0.20, 0.40),
    spearman_r = c(0.88, 0.85, 0.30)
  )

  high_conf <- .module1_select_high_confidence_footprints(
    motif_supported_stats,
    r_cutoff = 0.8
  )

  expect_equal(high_conf$fp_id, c("chr1:100-140", "chr1:200-240"))
  expect_true(all(c("fp_id", "chr", "start", "end", "atac_peak") %in% names(high_conf)))
  expect_equal(high_conf$start, c(100L, 200L))
})

test_that("Module 1 helpers merge Pearson and Spearman stats with best method", {
  pearson <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:100-140", "chr1:200-240"),
    atac_peak = c("chr1:90-160", "chr1:90-160", "chr1:180-260"),
    tf = c("TF_A", "TF_B", "TF_C"),
    pearson_r = c(0.91, 0.82, 0.10),
    pearson_p = c(0.01, 0.02, 0.80),
    pearson_p_adj = c(0.03, 0.04, 0.90)
  )
  spearman <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:100-140", "chr1:200-240"),
    atac_peak = c("chr1:90-160", "chr1:90-160", "chr1:180-260"),
    tf = c("TF_A", "TF_B", "TF_C"),
    spearman_r = c(0.88, 0.90, 0.20),
    spearman_p = c(0.02, 0.01, 0.70),
    spearman_p_adj = c(0.04, 0.03, 0.85)
  )

  merged <- .module1_merge_tfbs_stats(pearson, spearman, r_cutoff = 0.8)

  expect_equal(merged$best_method, c("pearson", "spearman", "spearman"))
  expect_equal(merged$pass, c(TRUE, TRUE, FALSE))
  expect_equal(merged$best_r, c(0.91, 0.90, 0.20), tolerance = 1e-8)
})

test_that("Module 1 helpers build sparse tfbs_links with multiple TFs per footprint", {
  fixture <- module1_tiny_fixture()
  high_conf <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:200-240"),
    chr = c("chr1", "chr1"),
    start = c(100L, 200L),
    end = c(140L, 240L),
    atac_peak = c("chr1:90-160", "chr1:180-260")
  )
  stats <- tibble::tibble(
    fp_id = c("chr1:100-140", "chr1:100-140", "chr1:200-240", fixture$expected_excluded_fp),
    atac_peak = c("chr1:90-160", "chr1:90-160", "chr1:180-260", "chr2:200-280"),
    tf = c("TF_A", "TF_B", "TF_C", "TF_A"),
    pearson_r = c(0.91, 0.82, 0.20, 0.95),
    pearson_p = c(0.01, 0.02, 0.80, 0.01),
    pearson_p_adj = c(0.03, 0.04, 0.90, 0.03),
    spearman_r = c(0.88, 0.90, 0.10, 0.96),
    spearman_p = c(0.02, 0.01, 0.70, 0.01),
    spearman_p_adj = c(0.04, 0.03, 0.85, 0.03),
    best_r = c(0.91, 0.90, 0.20, 0.96),
    best_method = c("pearson", "spearman", "pearson", "spearman"),
    pass = c(TRUE, TRUE, FALSE, TRUE)
  )

  links <- .module1_build_tfbs_links(
    tfbs_stats = stats,
    high_confidence_footprints = high_conf,
    omics_data = fixture$omics_data
  )

  expect_equal(nrow(links), 2L)
  expect_equal(links$fp_id, c("chr1:100-140", "chr1:100-140"))
  expect_equal(links$tf, c("TF_A", "TF_B"))
  expect_true(all(links$condition_support > 0L))
  expect_false(fixture$expected_excluded_fp %in% links$fp_id)
})
