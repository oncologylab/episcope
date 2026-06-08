test_that("gene TSS tables are normalized from file overrides", {
  tmp <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(gene = c("A", "B"), chrom = c("1", "chr2"), tss = c(10L, 20L), strand = c("+", "-")), tmp)
  out <- load_gene_tss(tmp, ref_genome = "hg38", verbose = FALSE)
  expect_equal(names(out), c("target_gene", "target_chr", "target_tss", "target_strand"))
  expect_equal(out$target_chr, c("chr1", "chr2"))
  expect_equal(out$target_strand, c("+", "-"))
})

test_that("automatic gene TSS loading follows ref_genome", {
  testthat::skip_if_not_installed("EnsDb.Hsapiens.v86")
  testthat::skip_if_not_installed("EnsDb.Mmusculus.v79")
  hg <- load_gene_tss(ref_genome = "hg38", genes = "KRAS", verbose = FALSE)
  mm <- load_gene_tss(ref_genome = "mm10", genes = "Kras", verbose = FALSE)
  expect_true(any(hg$target_gene == "KRAS"))
  expect_true(any(mm$target_gene == "Kras"))
  expect_true(all(grepl("^chr", hg$target_chr)))
  expect_true(all(grepl("^chr", mm$target_chr)))
})

test_that("Module 2 can resolve gene TSS from project config", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = "chr1:100-140", s1 = 1, s2 = 2, s3 = 8, s4 = 9),
    fp_bound_condition = tibble::tibble(peak_ID = "chr1:100-140", s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    fp_annotation = tibble::tibble(fp_peak = "chr1:100-140", atac_peak = "chr1:90-160", motifs = "M_A", tfs = "TF_A"),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2"), HGNC = c("TF_A", "GENE_UP"), s1 = c(1, 1), s2 = c(2, 2), s3 = c(8, 8), s4 = c(9, 9)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2"), HGNC = c("TF_A", "GENE_UP"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = "TF_A"
  )
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = "chr1:100-140", chr = "chr1", start = 100L, end = 140L, atac_peak = "chr1:90-160", tf = "TF_A")
  tss_path <- tempfile(fileext = ".csv")
  readr::write_csv(tibble::tibble(target_gene = "GENE_UP", target_chr = "chr1", target_tss = 120L, target_strand = "+"), tss_path)
  res <- predict_tf_targets(compact, pred, project_config = list(ref_genome = "hg38", gene_tss = tss_path, module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), max_distance_bp = 1000, n_cores = 1, verbose = FALSE)
  expect_true(check_predicted_links(res))
  expect_true(nrow(res$links) > 0)
})
