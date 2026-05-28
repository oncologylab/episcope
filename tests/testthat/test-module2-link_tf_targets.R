test_that("predicted_tfbs is compact and exportable", {
  tfbs_links <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:100-140"), chr = "chr1", start = 100L, end = 140L, atac_peak = "chr1:90-160", tf = c("TF_A", "TF_B"), best_r = c(0.9, 0.8), best_method = c("pearson", "spearman"), condition_support = c(2L, 1L))
  pred <- build_predicted_tfbs(tfbs_links)
  expect_true(all(c("tfbs_id", "fp_id", "chr", "start", "end", "atac_peak", "tf") %in% names(pred)))
  expect_false(any(c("best_r", "best_method") %in% names(pred)))
  tmp <- tempfile(fileext = ".bed")
  export_predicted_tfbs_bed(pred, out_file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(length(readLines(tmp)), 2L)
})

test_that("Module 2 computes TF-target first and restricts FP-target candidates", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1, 9), s2 = c(2, 8), s3 = c(8, 2), s4 = c(9, 1)),
    fp_bound_condition = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1L, 1L), s2 = c(1L, 1L), s3 = c(1L, 1L), s4 = c(1L, 1L)),
    fp_annotation = tibble::tibble(fp_peak = c("chr1:100-140", "chr1:500-540"), atac_peak = c("chr1:90-160", "chr1:490-560"), motifs = c("M_A", "M_B"), tfs = c("TF_A", "TF_B")),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = c(1, 9, 1, 9), s2 = c(2, 8, 2, 8), s3 = c(8, 2, 8, 2), s4 = c(9, 1, 9, 1)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  res <- link_tf_targets(compact, pred, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), max_distance_bp = 1000, n_cores = 1, verbose = FALSE)
  expect_true(validate_module2_links(res))
  expect_true(all(c("tf", "target_gene", "pass") %in% names(res$tf_target_corr)))
  expect_true(any(res$tf_target_corr$tf == "TF_A" & res$tf_target_corr$target_gene == "GENE_UP" & res$tf_target_corr$pass))
  expect_false(any(res$candidates$fp_id == "chr1:100-140" & res$candidates$target_gene == "GENE_DOWN"))
  expect_equal(nrow(unique(res$fp_target_corr[, c("fp_id", "target_gene")])), nrow(res$fp_target_corr))
  expect_true(all(res$links$module2_link_pass))
})


test_that("Module 2 streams predicted TFBS manifests", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1, 9), s2 = c(2, 8), s3 = c(8, 2), s4 = c(9, 1)),
    fp_bound_condition = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1L, 1L), s2 = c(1L, 1L), s3 = c(1L, 1L), s4 = c(1L, 1L)),
    fp_annotation = tibble::tibble(fp_peak = c("chr1:100-140", "chr1:500-540"), atac_peak = c("chr1:90-160", "chr1:490-560"), motifs = c("M_A", "M_B"), tfs = c("TF_A", "TF_B")),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = c(1, 9, 1, 9), s2 = c(2, 8, 2, 8), s3 = c(8, 2, 8, 2), s4 = c(9, 1, 9, 1)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  pred_dir <- tempfile("predicted-tfbs-")
  dir.create(pred_dir)
  pred_path <- file.path(pred_dir, "module1_predicted_tfbs_chunk_0001.csv")
  readr::write_csv(pred, pred_path)
  pred_manifest_path <- file.path(pred_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(tibble::tibble(chunk_id = 1L, path = pred_path, format = "csv", n_rows = nrow(pred)), pred_manifest_path)
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  out_dir <- tempfile("module2-stream-")
  res <- link_tf_targets(compact, pred_manifest_path, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), output_dir = out_dir, max_distance_bp = 1000, n_cores = 1, output_format = "csv", verbose = FALSE)
  expect_true(validate_module2_links(res))
  expect_true(file.exists(res$reports$links_manifest))
  loaded <- load_module2_links(out_dir)
  expect_true(validate_module2_links(loaded))
  q <- query_module2_links(loaded, tf = "TF_A")
  expect_true(nrow(q) > 0)
  expect_true(all(q$tf == "TF_A"))
})

test_that("Module 2 top TF target reports write self-contained HTML", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1, 9), s2 = c(2, 8), s3 = c(8, 2), s4 = c(9, 1)),
    fp_bound_condition = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1L, 1L), s2 = c(1L, 1L), s3 = c(1L, 1L), s4 = c(1L, 1L)),
    fp_annotation = tibble::tibble(fp_peak = c("chr1:100-140", "chr1:500-540"), atac_peak = c("chr1:90-160", "chr1:490-560"), motifs = c("M_A", "M_B"), tfs = c("TF_A", "TF_B")),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3"), HGNC = c("TF_A", "TF_B", "GENE_UP"), s1 = c(1, 9, 1), s2 = c(2, 8, 2), s3 = c(8, 2, 8), s4 = c(9, 1, 9)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3"), HGNC = c("TF_A", "TF_B", "GENE_UP"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  gene_tss <- tibble::tibble(target_gene = "GENE_UP", target_chr = "chr1", target_tss = 120L, target_strand = "+")
  res <- link_tf_targets(compact, pred, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), max_distance_bp = 1000, n_cores = 1, verbose = FALSE)
  out_dir <- tempfile("module2-reports-")
  man <- export_top_tf_targets(res, output_dir = out_dir, tfs = "TF_A", top_n = 1L, verbose = FALSE)
  expect_equal(nrow(man), 1L)
  expect_true(file.exists(man$path[[1L]]))
  html <- readLines(man$path[[1L]], warn = FALSE)
  expect_true(any(grepl("Export SVG", html, fixed = TRUE)))
  expect_true(any(grepl("const FULL_NODES", html, fixed = TRUE)))
})

test_that("Module 2 report builder writes direct HTML without pdf folders", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = "chr1:100-140", s1 = 5, s2 = 6),
    fp_bound_condition = tibble::tibble(peak_ID = "chr1:100-140", s1 = 1L, s2 = 1L),
    fp_annotation = tibble::tibble(fp_peak = "chr1:100-140", atac_peak = "chr1:90-160", motifs = "M_A", tfs = "TF_A"),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2"), HGNC = c("TF_A", "TF_B"), s1 = c(3, 4), s2 = c(5, 6)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2"), HGNC = c("TF_A", "TF_B"), s1 = 1L, s2 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- as_multiomic_object(omics, verbose = FALSE)
  module2 <- list(
    tf_target_corr = tibble::tibble(tf = c("TF_A", "TF_B"), target_gene = c("TF_B", "TF_A"), best_r = c(0.9, 0.1), pass = c(TRUE, FALSE)),
    fp_target_corr = tibble::tibble(fp_id = "chr1:100-140", target_gene = "TF_B", best_r = 0.8, pass = TRUE),
    links = tibble::tibble(link_id = "l1", tf = "TF_A", fp_id = "chr1:100-140", target_gene = "TF_B", candidate_id = "c1", module2_link_pass = TRUE),
    manifest = tibble::tibble()
  )
  out_dir <- tempfile("module2-report-builder-")
  man <- build_module2_reports(module2, multiomic_data = NULL, output_dir = out_dir, reports = c("direct_tf_tf", "tf_tf_connectivity"), k_values = 2L, verbose = FALSE)
  expect_equal(nrow(man), 2L)
  expect_true(all(file.exists(man$path)))
  expect_true(all(grepl("/html/", man$path, fixed = TRUE)))
  expect_false(dir.exists(file.path(out_dir, "pdf")))
})
