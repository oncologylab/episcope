test_that("predicted_tfbs is compact and exportable", {
  prediction_stats <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:100-140"), chr = "chr1", start = 100L, end = 140L, atac_peak = "chr1:90-160", tf = c("TF_A", "TF_B"), best_r = c(0.9, 0.8), best_method = c("pearson", "spearman"), condition_support = c(2L, 1L))
  pred <- output_predicted_tfbs(prediction_stats)
  expect_true(all(c("tfbs_id", "fp_id", "chr", "start", "end", "atac_peak", "tf") %in% names(pred)))
  expect_false(any(c("best_r", "best_method") %in% names(pred)))
  tmp <- tempfile(fileext = ".bed")
  export_predicted_tfbs_bed(pred, out_file = tmp)
  expect_true(file.exists(tmp))
  expect_equal(length(readLines(tmp)), 2L)

  tf_tmp <- tempfile(fileext = ".bed")
  export_predicted_tfbs_bed(pred, out_file = tf_tmp, tf = "TF_A")
  expect_equal(length(readLines(tf_tmp)), 1L)
  expect_match(readLines(tf_tmp), "TF_A|chr1:100-140", fixed = TRUE)

  split_dir <- tempfile("predicted-tfbs-bed-")
  manifest <- export_predicted_tfbs_bed(pred, out_dir = split_dir, tf = "TF_A", split_by = "tf")
  expect_equal(manifest$tf, "TF_A")
  expect_equal(manifest$n_rows, 1L)
  expect_true(file.exists(file.path(split_dir, "TF_A.bed")))
  expect_false(file.exists(file.path(split_dir, "TF_B.bed")))
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
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  res <- predict_tf_targets(compact, pred, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), max_distance_bp = 1000, n_cores = 1, verbose = FALSE)
  expect_true(check_predicted_links(res))
  expect_true(all(c("tf", "target_gene", "pass") %in% names(res$tf_target_corr)))
  expect_true(any(res$tf_target_corr$tf == "TF_A" & res$tf_target_corr$target_gene == "GENE_UP" & res$tf_target_corr$pass))
  expect_false(any(res$candidates$fp_id == "chr1:100-140" & res$candidates$target_gene == "GENE_DOWN"))
  expect_equal(nrow(unique(res$fp_target_corr[, c("fp_id", "target_gene")])), nrow(res$fp_target_corr))
  expect_true(all(res$links$module2_link_pass))
})

test_that("Module 2 correlation cutoffs tolerate rows without a best method", {
  x <- tibble::tibble(
    tf = c("TF_A", "TF_B"),
    target_gene = c("GENE_A", "GENE_B"),
    pearson_r = c(0.9, NA_real_),
    pearson_p = c(0.01, NA_real_),
    pearson_fdr = c(0.02, NA_real_),
    spearman_r = c(0.8, NA_real_),
    spearman_p = c(0.02, NA_real_),
    spearman_fdr = c(0.03, NA_real_)
  )

  out <- .module2_apply_corr_cutoffs(x, list(r = 0.3, p = 0.05, fdr = NULL))

  expect_true(out$pass[[1L]])
  expect_false(out$pass[[2L]])
  expect_equal(out$best_method[[1L]], "pearson")
  expect_true(is.na(out$best_method[[2L]]))
  expect_true(is.na(out$best_p[[2L]]))
})

test_that("Module 2 keeps regulatory-prior-only candidates outside TSS windows", {
  predicted_tfbs <- tibble::tibble(
    fp_id = "fp1",
    tf = "TF_A",
    chr = "chr1",
    start = 100L,
    end = 120L,
    atac_peak = "chr1:90-160"
  )
  tf_target_pass <- tibble::tibble(tf = "TF_A", target_gene = "GENE_FAR")
  gene_tss <- tibble::tibble(
    target_gene = "GENE_FAR",
    target_chr = "chr1",
    target_tss = 1000000L,
    target_strand = "+"
  )
  regulatory_prior <- tibble::tibble(
    fp_id = "fp1",
    target_gene = "GENE_FAR",
    prior_id = "prior1",
    prior_source = "test_prior",
    prior_score = 0.9,
    prior_status = "supported"
  )

  candidates <- .module2_build_candidates(
    predicted_tfbs = predicted_tfbs,
    tf_target_pass = tf_target_pass,
    gene_tss = gene_tss,
    regulatory_prior = regulatory_prior,
    max_distance_bp = 100L
  )

  expect_equal(nrow(candidates), 1L)
  expect_equal(candidates$candidate_source, "regulatory_prior")
  expect_false(candidates$within_tss_window)
  expect_true(candidates$prior_supported)
  expect_equal(candidates$prior_id, "prior1")
  expect_equal(candidates$prior_source, "test_prior")
  expect_equal(candidates$distance_to_tss, -999890)
})

test_that("Module 2 vectorized candidate merge preserves combined evidence", {
  fp_meta <- tibble::tibble(
    fp_id = "fp1",
    chr = "chr1",
    start = 100L,
    end = 120L,
    atac_peak = "chr1:90-160"
  )
  tf_target_pass <- tibble::tibble(tf = "TF_A", target_gene = "GENE_NEAR")
  gene_tss <- tibble::tibble(
    target_gene = "GENE_NEAR",
    target_chr = "chr1",
    target_tss = 110L,
    target_strand = "+"
  )
  regulatory_prior <- tibble::tibble(
    fp_id = "fp1",
    target_gene = "GENE_NEAR",
    prior_id = "prior1",
    prior_source = "test_prior",
    prior_score = 0.9,
    prior_status = "supported"
  )

  candidates <- .module2_build_candidates_from_fp_meta(
    fp_meta = fp_meta,
    tf_target_pass = tf_target_pass,
    gene_tss = gene_tss,
    regulatory_prior = regulatory_prior,
    max_distance_bp = 100L
  )

  expect_equal(nrow(candidates), 1L)
  expect_equal(candidates$candidate_source, "both")
  expect_true(candidates$within_tss_window)
  expect_true(candidates$prior_supported)
  expect_equal(candidates$prior_id, "prior1")
})

test_that("Module 2 candidate preflight uses a genomic-density upper bound", {
  fp_meta <- tibble::tibble(
    fp_id = c("fp1", "fp2"),
    chr = "chr1"
  )
  gene_tss <- tibble::tibble(
    target_gene = c("G1", "G2", "G3"),
    target_chr = "chr1",
    target_tss = c(0, 100, 1000)
  )
  prior <- tibble::tibble(fp_id = "fp1", target_gene = "G3")

  bound <- .module2_candidate_row_upper_bound(
    fp_meta = fp_meta,
    target_genes = gene_tss$target_gene,
    gene_tss = gene_tss,
    max_distance_bp = 100,
    regulatory_prior = prior
  )

  expect_equal(bound, 5)
})

test_that("Module 2 empty FP-target candidate sets keep the standard schema", {
  predicted_tfbs <- tibble::tibble(
    fp_id = "fp1",
    tf = "TF_A",
    chr = "chr1",
    start = 100L,
    end = 120L,
    atac_peak = "chr1:90-160"
  )
  tf_target_pass <- tibble::tibble(tf = "TF_A", target_gene = "GENE_FAR")
  gene_tss <- tibble::tibble(
    target_gene = "GENE_FAR",
    target_chr = "chr1",
    target_tss = 1000000L,
    target_strand = "+"
  )

  candidates <- .module2_build_candidates(
    predicted_tfbs = predicted_tfbs,
    tf_target_pass = tf_target_pass,
    gene_tss = gene_tss,
    max_distance_bp = 100L
  )

  expect_equal(nrow(candidates), 0L)
  expect_true(all(c(
    "candidate_id", "fp_id", "target_gene", "chr", "start", "end",
    "atac_peak", "target_chr", "target_tss", "target_strand",
    "distance_to_tss", "candidate_source", "within_tss_window",
    "prior_supported", "prior_id", "prior_source", "prior_score",
    "prior_status"
  ) %in% names(candidates)))
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
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  pred_dir <- tempfile("predicted-tfbs-")
  dir.create(pred_dir)
  pred_path <- file.path(pred_dir, "module1_predicted_tfbs_chunk_0001.csv")
  readr::write_csv(pred, pred_path)
  pred_manifest_path <- file.path(pred_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(tibble::tibble(chunk_id = 1L, path = pred_path, format = "csv", n_rows = nrow(pred)), pred_manifest_path)
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  out_dir <- tempfile("module2-stream-")
  res <- predict_tf_targets(compact, pred_manifest_path, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), output_dir = out_dir, max_distance_bp = 1000, n_cores = 1, output_format = "csv", verbose = FALSE)
  expect_true(check_predicted_links(res))
  expect_true(file.exists(res$reports$links_manifest))
  expect_true(file.exists(file.path(out_dir, "module2_qc_report.html")))
  expect_false(dir.exists(file.path(out_dir, "module2_fp_target_candidates_chunks")))
  expect_false(dir.exists(file.path(out_dir, "module2_fp_target_corr_chunks")))
  expect_false(file.exists(file.path(out_dir, "module2_fp_target_candidates_manifest.csv")))
  expect_false(file.exists(file.path(out_dir, "module2_fp_target_corr_manifest.csv")))
  expect_true(file.exists(file.path(out_dir, "data", "candidates", "fp_target_candidates.csv.gz")))
  expect_true(file.exists(file.path(out_dir, "data", "correlations", "fp_target_corr.csv.gz")))
  expect_true(file.exists(file.path(out_dir, "data", "links", "module2_links_manifest.csv")))
  expect_true(dir.exists(file.path(out_dir, "reports")))
  loaded <- load_predicted_links(out_dir)
  expect_true(check_predicted_links(loaded))
  q <- query_predicted_links(loaded, tf = "TF_A")
  expect_true(nrow(q) > 0)
  expect_true(all(q$tf == "TF_A"))
})

test_that("Module 2 streams one deduplicated FP-target candidate universe", {
  omics <- list(
    fp_score_condition_qn = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1, 9), s2 = c(2, 8), s3 = c(8, 2), s4 = c(9, 1)),
    fp_bound_condition = tibble::tibble(peak_ID = c("chr1:100-140", "chr1:500-540"), s1 = c(1L, 1L), s2 = c(1L, 1L), s3 = c(1L, 1L), s4 = c(1L, 1L)),
    fp_annotation = tibble::tibble(fp_peak = c("chr1:100-140", "chr1:500-540"), atac_peak = c("chr1:90-160", "chr1:490-560"), motifs = c("M_A", "M_B"), tfs = c("TF_A", "TF_B")),
    rna_condition = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = c(1, 9, 1, 9), s2 = c(2, 8, 2, 8), s3 = c(8, 2, 8, 2), s4 = c(9, 1, 9, 1)),
    rna_expressed = tibble::tibble(ensembl_gene_id = c("g1", "g2", "g3", "g4"), HGNC = c("TF_A", "TF_B", "GENE_UP", "GENE_DOWN"), s1 = 1L, s2 = 1L, s3 = 1L, s4 = 1L),
    tf_list = c("TF_A", "TF_B")
  )
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred_dir <- tempfile("predicted-tfbs-dedup-")
  dir.create(pred_dir)
  pred1 <- tibble::tibble(fp_id = "chr1:100-140", chr = "chr1", start = 100L, end = 140L, atac_peak = "chr1:90-160", tf = "TF_A")
  pred2 <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_B", "TF_B"))
  pred_path1 <- file.path(pred_dir, "module1_predicted_tfbs_chunk_0001.csv")
  pred_path2 <- file.path(pred_dir, "module1_predicted_tfbs_chunk_0002.csv")
  readr::write_csv(pred1, pred_path1)
  readr::write_csv(pred2, pred_path2)
  pred_manifest_path <- file.path(pred_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(tibble::tibble(
    chunk_id = c(1L, 2L),
    path = c(pred_path1, pred_path2),
    format = "csv",
    n_rows = c(nrow(pred1), nrow(pred2))
  ), pred_manifest_path)
  gene_tss <- tibble::tibble(target_gene = c("GENE_UP", "GENE_DOWN"), target_chr = "chr1", target_tss = c(120L, 520L), target_strand = "+")
  out_dir <- tempfile("module2-stream-dedup-")

  res <- predict_tf_targets(compact, pred_manifest_path, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), output_dir = out_dir, max_distance_bp = 1000, n_cores = 1, output_format = "csv", verbose = FALSE)

  expect_true(check_predicted_links(res))
  manifest <- readr::read_csv(file.path(out_dir, "module2_manifest.csv"), show_col_types = FALSE)
  cand_row <- manifest[manifest$table == "module2_fp_target_candidates", , drop = FALSE]
  expect_equal(nrow(cand_row), 1L)
  expect_equal(cand_row$format, "csv")
  expect_match(cand_row$path, "data/candidates/fp_target_candidates.csv.gz", fixed = TRUE)
  candidates <- .module2_read_predicted_chunk(cand_row$path[[1L]], cand_row$format[[1L]])
  expect_equal(nrow(candidates), nrow(unique(candidates[, c("fp_id", "target_gene")])))
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
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(fp_id = c("chr1:100-140", "chr1:500-540"), chr = "chr1", start = c(100L, 500L), end = c(140L, 540L), atac_peak = c("chr1:90-160", "chr1:490-560"), tf = c("TF_A", "TF_B"))
  gene_tss <- tibble::tibble(target_gene = "GENE_UP", target_chr = "chr1", target_tss = 120L, target_strand = "+")
  res <- predict_tf_targets(compact, pred, gene_tss, project_config = list(module2 = list(threshold_tf_target_corr_r = 0.8, threshold_fp_target_corr_r = 0.8)), max_distance_bp = 1000, n_cores = 1, verbose = FALSE)
  out_dir <- tempfile("module2-reports-")
  man <- report_top_tf_targets(res, output_dir = out_dir, tfs = "TF_A", top_n = 1L, verbose = FALSE)
  expect_equal(nrow(man), 1L)
  expect_true(file.exists(man$path[[1L]]))
  expect_equal(
    normalizePath(dirname(man$path[[1L]]), winslash = "/", mustWork = FALSE),
    normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  )
  html <- readLines(man$path[[1L]], warn = FALSE)
  expect_true(any(grepl("Export SVG", html, fixed = TRUE)))
  expect_true(any(grepl("const FULL_NODES", html, fixed = TRUE)))
  expect_false(any(grepl('id="viewMode"', html, fixed = TRUE)))
  expect_true(any(grepl('id="topRange"', html, fixed = TRUE)))
})

test_that("Module 2 regulon reports support combined TFs and condition payloads", {
  conditions <- c("condition_a", "condition_b")
  fp_ids <- c("chr1:100-140", "chr1:200-240", "chr1:300-340")
  omics <- list(
    fp_score_condition_qn = tibble::as_tibble(data.frame(
      peak_ID = fp_ids,
      condition_a = c(1, 2, 3),
      condition_b = c(4, 2, 1),
      check.names = FALSE
    )),
    fp_bound_condition = tibble::as_tibble(data.frame(
      peak_ID = fp_ids,
      condition_a = c(1L, 1L, 1L),
      condition_b = c(1L, 1L, 0L),
      check.names = FALSE
    )),
    fp_annotation = tibble::tibble(
      fp_peak = fp_ids,
      atac_peak = c("chr1:90-150", "chr1:190-250", "chr1:290-350"),
      motifs = c("M_A", "M_B", "M_C"),
      tfs = c("TF_A", "TF_B", "TF_C")
    ),
    rna_condition = tibble::as_tibble(data.frame(
      ensembl_gene_id = paste0("g", 1:5),
      HGNC = c("TF_A", "TF_B", "TF_C", "TARGET_1", "TARGET_2"),
      condition_a = c(4, 3, 2, 8, 1),
      condition_b = c(8, 2, 3, 2, 9),
      check.names = FALSE
    )),
    rna_expressed = tibble::as_tibble(data.frame(
      ensembl_gene_id = paste0("g", 1:5),
      HGNC = c("TF_A", "TF_B", "TF_C", "TARGET_1", "TARGET_2"),
      matrix(1L, 5, 2, dimnames = list(NULL, conditions)),
      check.names = FALSE
    )),
    tf_list = c("TF_A", "TF_B", "TF_C")
  )
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  module2 <- list(
    links = tibble::tibble(
      link_id = paste0("L", 1:3),
      tf = c("TF_A", "TF_B", "TF_C"),
      fp_id = fp_ids,
      target_gene = c("TARGET_1", "TARGET_2", "TARGET_1"),
      candidate_id = c("C1", "C2", "C3"),
      module2_link_pass = TRUE
    ),
    candidates = tibble::tibble(
      candidate_id = c("C1", "C2", "C3"),
      fp_id = fp_ids,
      target_gene = c("TARGET_1", "TARGET_2", "TARGET_1")
    ),
    tf_target_corr = tibble::tibble(
      tf = c("TF_A", "TF_B", "TF_C"),
      target_gene = c("TARGET_1", "TARGET_2", "TARGET_1"),
      best_r = c(0.72, 0.68, 0.64),
      pass = TRUE
    ),
    fp_target_corr = tibble::tibble(
      fp_id = fp_ids,
      target_gene = c("TARGET_1", "TARGET_2", "TARGET_1"),
      best_r = c(0.70, 0.86, 0.99),
      pass = TRUE
    )
  )
  out_dir <- tempfile("module2-combined-regulon-")
  manifest <- report_top_tf_targets(
    module2,
    output_dir = out_dir,
    tfs = c("TF_A", "TF_B"),
    top_n = 2L,
    default_top_n = 1L,
    verbose = FALSE,
    multiomic_data = compact,
    combine_tfs = TRUE,
    default_condition = "condition_a"
  )

  expect_equal(nrow(manifest), 1L)
  html <- paste(readLines(manifest$path[[1L]], warn = FALSE), collapse = "\n")
  expect_match(html, "REGULON_PAYLOAD_DEFLATE_BASE64", fixed = TRUE)
  expect_match(html, "const DEFAULT_TOP_N=1", fixed = TRUE)
  expect_match(html, "more than 500 nodes", fixed = TRUE)
  expect_no_match(basename(manifest$path[[1L]]), "log2fc1", fixed = TRUE)
  expect_match(html, "condition_a", fixed = TRUE)
  expect_match(html, 'id="nodePicker"', fixed = TRUE)
  expect_match(html, 'id="nodeSearch" type="search"', fixed = TRUE)
  expect_match(html, 'id="nodePaste"', fixed = TRUE)
  expect_match(html, "Show pasted nodes only", fixed = TRUE)
  expect_match(html, "parsePastedNodes", fixed = TRUE)
  expect_match(html, "Check all", fixed = TRUE)
  expect_match(html, "Uncheck all", fixed = TRUE)
  expect_match(html, "const DEFAULT_FP_R_CUTOFF=0.5", fixed = TRUE)
  expect_match(html, "Min FP R &gt;", fixed = TRUE)
  expect_match(html, "requestedTop", fixed = TRUE)
  expect_match(html, "Math.min(requested,eligible)", fixed = TRUE)
  expect_match(html, "selected:''", fixed = TRUE)
  expect_match(html, "eligibleTargets(fpCut()).length", fixed = TRUE)
  expect_match(html, "Compound spring", fixed = TRUE)
  expect_match(html, "Cond1 minus Cond2", fixed = TRUE)
  expect_match(html, 'id="cond1Select"', fixed = TRUE)
  expect_match(html, 'id="cond2Select"', fixed = TRUE)
  expect_match(html, 'id="appearance"', fixed = TRUE)
  expect_no_match(html, '<details class="appearance"', fixed = TRUE)
  expect_match(html, 'id="geneMin"', fixed = TRUE)
  expect_match(html, 'id="geneMax"', fixed = TRUE)
  expect_match(html, 'id="tfMin"', fixed = TRUE)
  expect_match(html, 'id="tfMax"', fixed = TRUE)
  expect_match(html, 'id="edgeOpacity"', fixed = TRUE)
  expect_match(html, 'id="edgeOpacityValue"', fixed = TRUE)
  expect_match(html, 'id="nodeOpacity"', fixed = TRUE)
  expect_match(html, 'id="nodeOpacityValue"', fixed = TRUE)
  expect_match(html, "none.textContent=noneLabel", fixed = TRUE)
  expect_match(html, "'All predictions'", fixed = TRUE)
  expect_match(html, "'None'", fixed = TRUE)
  expect_match(html, "!l.is_seed_edge", fixed = TRUE)
  expect_match(html, "predictionAutoAll", fixed = TRUE)
  expect_match(html, "displayed seed-TF targets", fixed = TRUE)
  expect_match(html, "c=>c!==ctl.c1.value", fixed = TRUE)
  expect_match(html, "e.delta=e.score1-e.score2", fixed = TRUE)
  expect_match(html, "Math.log2((a+p)/(b+p))", fixed = TRUE)
  expect_match(html, "g.mode==='comparison'", fixed = TRUE)
  expect_match(html, "Node size: maximum expression", fixed = TRUE)
  expect_match(html, 'id="topRange"', fixed = TRUE)
  expect_match(html, "allowed.has(l.to)", fixed = TRUE)
  expect_match(
    html,
    'id="spacing" type="range" min="0.1" max="2" step="0.01" value="0.5"',
    fixed = TRUE
  )
  expect_match(html, 'id="spacingValue" type="number"', fixed = TRUE)
  expect_match(html, "function nodeBoundary", fixed = TRUE)
  expect_match(html, "function fitView", fixed = TRUE)
  expect_match(html, "needsFit:true", fixed = TRUE)
  expect_match(html, "craftgrnPrepareIllustratorSvg", fixed = TRUE)
  expect_match(html, "data-export-role", fixed = TRUE)
  expect_match(html, "label-halo", fixed = TRUE)
  expect_match(html, "label-fill", fixed = TRUE)
  expect_match(html, "setAttribute('version','1.1')", fixed = TRUE)
  expect_match(html, "nodeBoundary(b,rb,-ux,-uy)-1", fixed = TRUE)
  expect_no_match(html, "rb+8", fixed = TRUE)
  expect_no_match(html, 'id="viewMode"', fixed = TRUE)
  expect_no_match(html, "B minus A", fixed = TRUE)
  expect_false(grepl("cytoscape", html, ignore.case = TRUE))
  ranked_links <- data.table::fread(file.path(
    out_dir,
    "csv",
    "TF_A_TF_B_Module2_top2_links.csv"
  ))
  ranked_targets <- unique(ranked_links[, .(target_gene = to, target_rank)])
  expect_equal(ranked_targets[target_gene == "TARGET_2", target_rank], 1L)
  expect_equal(ranked_targets[target_gene == "TARGET_1", target_rank], 2L)
  payload <- craftgrn:::.module2_regulon_build(
    module2 = module2,
    tfs = c("TF_A", "TF_B"),
    top_n = 2L,
    default_top_n = 1L,
    multiomic_data = compact,
    default_condition = "condition_a"
  )
  expect_no_warning(
    no_supporting_payload <- craftgrn:::.module2_regulon_build(
      module2 = module2,
      tfs = c("TF_A", "TF_B"),
      top_n = 2L,
      multiomic_data = compact,
      max_linked_tfs = 0L,
      max_linked_tfs_per_target = 0L
    )
  )
  expect_false("TF_C" %in% no_supporting_payload$nodes$id)
  expect_true(all(payload$condition_links$active %in% TRUE))
  expect_lt(nrow(payload$condition_links), nrow(payload$links) * length(conditions))

  mixed_case_module2 <- module2
  mixed_case_module2$links$target_gene[2L] <- "C2orf15"
  mixed_case_module2$candidates$target_gene[2L] <- "C2orf15"
  mixed_case_module2$tf_target_corr$target_gene[2L] <- "C2orf15"
  mixed_case_module2$fp_target_corr$target_gene[2L] <- "C2orf15"
  mixed_case_payload <- craftgrn:::.module2_regulon_build(
    module2 = mixed_case_module2,
    tfs = c("TF_A", "TF_B"),
    top_n = Inf,
    target_genes = c("TARGET_1", "C2ORF15"),
    max_linked_tfs = 0L,
    max_linked_tfs_per_target = 0L
  )
  expect_setequal(
    mixed_case_payload$nodes[is_target_gene %in% TRUE, id],
    c("TARGET_1", "C2ORF15")
  )

  selected_dir <- tempfile("module2-selected-regulon-")
  selected <- report_top_tf_targets(
    module2,
    output_dir = selected_dir,
    tfs = c("TF_A", "TF_B"),
    top_n = Inf,
    verbose = FALSE,
    combine_tfs = TRUE,
    target_genes = "TARGET_2"
  )
  expect_equal(nrow(selected), 1L)
  expect_match(basename(selected$path[[1L]]), "selected1", fixed = TRUE)
  selected_nodes <- data.table::fread(file.path(
    selected_dir,
    "csv",
    "TF_A_TF_B_Module2_selected1_nodes.csv"
  ))
  expect_setequal(selected_nodes[is_target_gene == TRUE, id], "TARGET_2")
})

test_that("regulon payload serialization preserves correlation precision", {
  value <- 0.300038281590023
  encoded <- craftgrn:::.module2_report_browser_encode_browser_json_deflate_base64(
    list(value = value),
    digits = 15L
  )
  decoded <- jsonlite::fromJSON(rawToChar(memDecompress(
    jsonlite::base64_dec(encoded),
    type = "gzip"
  )))
  expect_equal(decoded$value, value, tolerance = 1e-15)
})

test_that("Module 2 report APIs consume report-first output directories", {
  conditions <- paste0("s", 1:6)
  fp_ids <- sprintf("chr1:%d-%d", c(100, 200, 300, 400), c(140, 240, 340, 440))
  atac_ids <- sprintf("chr1:%d-%d", c(90, 190, 290, 390), c(150, 250, 350, 450))
  profiles <- rbind(
    TF_A = c(1, 2, 3, 4, 5, 6),
    TF_B = c(2, 3, 4, 5, 6, 7),
    TF_C = c(1, 3, 4, 6, 7, 9),
    TF_D = c(3, 4, 6, 7, 9, 10)
  )
  colnames(profiles) <- conditions
  omics <- list(
    fp_score_condition_qn = tibble::as_tibble(data.frame(
      peak_ID = fp_ids,
      profiles + 0.25,
      check.names = FALSE
    )),
    fp_bound_condition = tibble::as_tibble(data.frame(
      peak_ID = fp_ids,
      matrix(1L, 4, 6, dimnames = list(NULL, conditions)),
      check.names = FALSE
    )),
    fp_annotation = tibble::tibble(
      fp_peak = fp_ids,
      atac_peak = atac_ids,
      motifs = paste0("M_", LETTERS[1:4]),
      tfs = paste0("TF_", LETTERS[1:4])
    ),
    rna_condition = tibble::as_tibble(data.frame(
      ensembl_gene_id = paste0("g", 1:4),
      HGNC = rownames(profiles),
      profiles,
      check.names = FALSE
    )),
    rna_expressed = tibble::as_tibble(data.frame(
      ensembl_gene_id = paste0("g", 1:4),
      HGNC = rownames(profiles),
      matrix(1L, 4, 6, dimnames = list(NULL, conditions)),
      check.names = FALSE
    )),
    tf_list = rownames(profiles)
  )
  compact <- craftgrn:::as_multiomic_object(omics, verbose = FALSE)
  pred <- tibble::tibble(
    fp_id = fp_ids,
    chr = "chr1",
    start = c(100L, 200L, 300L, 400L),
    end = c(140L, 240L, 340L, 440L),
    atac_peak = atac_ids,
    tf = rownames(profiles)
  )
  gene_tss <- tibble::tibble(
    target_gene = rownames(profiles),
    target_chr = "chr1",
    target_tss = c(120L, 220L, 320L, 420L),
    target_strand = "+"
  )
  module2_dir <- tempfile("module2-report-first-")
  predict_tf_targets(
    compact,
    pred,
    gene_tss,
    project_config = list(module2 = list(
      threshold_tf_target_corr_r = 0.8,
      threshold_fp_target_corr_r = 0.8
    )),
    output_dir = module2_dir,
    max_distance_bp = 1000,
    n_cores = 1,
    output_format = "csv",
    verbose = FALSE,
    write_qc_report = FALSE
  )

  top <- report_top_tf_targets(
    module2_dir,
    output_dir = file.path(module2_dir, "top"),
    tfs = "TF_A",
    top_n = 2L,
    verbose = FALSE
  )
  direct <- report_direct_tf_tf_regulations(
    module2_dir,
    output_dir = file.path(module2_dir, "direct"),
    multiomic_data = compact,
    k_values = 1L,
    verbose = FALSE
  )
  connectivity <- report_tf_tf_coregulations(
    module2_dir,
    output_dir = file.path(module2_dir, "connectivity"),
    multiomic_data = compact,
    k_values = 1L,
    verbose = FALSE
  )

  expect_equal(nrow(top), 1L)
  expect_equal(nrow(direct), 1L)
  expect_equal(nrow(connectivity), 1L)
  expect_true(all(file.exists(c(top$path, direct$path, connectivity$path))))
  expect_match(paste(readLines(direct$path[[1L]], warn = FALSE), collapse = "\n"), "Direct TF heatmap", fixed = TRUE)
  expect_match(paste(readLines(connectivity$path[[1L]], warn = FALSE), collapse = "\n"), "Composite TF-TF heatmap", fixed = TRUE)
})

test_that("Module 2 report builder writes distinct direct and connectivity HTML", {
  tfs <- c("TF_A", "TF_B", "TF_C", "TF_D")
  pairs <- expand.grid(tf = tfs, target_gene = tfs, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$tf != pairs$target_gene, , drop = FALSE]
  n_pairs <- nrow(pairs)
  module2 <- list(
    tf_target_corr = tibble::tibble(
      tf = pairs$tf,
      target_gene = pairs$target_gene,
      best_r = seq(0.95, 0.84, length.out = n_pairs),
      pass = TRUE
    ),
    fp_target_corr = tibble::tibble(
      fp_id = paste0("fp", seq_len(n_pairs)),
      target_gene = pairs$target_gene,
      best_r = seq(0.83, 0.72, length.out = n_pairs),
      pass = TRUE
    ),
    links = tibble::tibble(
      link_id = paste0("l", seq_len(n_pairs)),
      tf = pairs$tf,
      fp_id = paste0("fp", seq_len(n_pairs)),
      target_gene = pairs$target_gene,
      candidate_id = paste0("c", seq_len(n_pairs)),
      module2_link_pass = TRUE
    ),
    manifest = tibble::tibble()
  )
  out_dir <- tempfile("module2-report-builder-")
  man <- dplyr::bind_rows(
    report_direct_tf_tf_regulations(module2, multiomic_data = NULL, output_dir = out_dir, k_values = 1L, verbose = FALSE),
    report_tf_tf_coregulations(module2, multiomic_data = NULL, output_dir = out_dir, k_values = 1L, verbose = FALSE)
  )
  expect_equal(nrow(man), 2L)
  expect_true(all(file.exists(man$path)))
  expect_true(all(
    normalizePath(dirname(man$path), winslash = "/", mustWork = FALSE) ==
      normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  ))
  direct_html <- paste(readLines(man$path[man$report == "direct_tf_tf"][[1L]], warn = FALSE), collapse = "\n")
  conn_html <- paste(readLines(man$path[man$report == "tf_tf_connectivity"][[1L]], warn = FALSE), collapse = "\n")
  expect_true(grepl("Direct TF heatmap", direct_html, fixed = TRUE))
  expect_true(grepl("Composite TF-TF heatmap", conn_html, fixed = TRUE))
  expect_true(grepl(">Single<", conn_html, fixed = TRUE))
  expect_true(grepl(">Overlay<", conn_html, fixed = TRUE))
  expect_false(grepl("Direct TF heatmap", conn_html, fixed = TRUE))
  expect_true(grepl("data:image/png;base64", direct_html, fixed = TRUE))
  expect_true(grepl("data:image/png;base64", conn_html, fixed = TRUE))
  expect_true(grepl("const HEATMAP_IMAGES", conn_html, fixed = TRUE))
  expect_true(grepl("craftgrnPrepareIllustratorSvg", direct_html, fixed = TRUE))
  expect_true(grepl("craftgrnPrepareIllustratorSvg", conn_html, fixed = TRUE))
  expect_true(grepl("exportArrowForLine", direct_html, fixed = TRUE))
  expect_true(grepl("exportPaneArrowForLine", conn_html, fixed = TRUE))
  expect_true(grepl("label-halo", direct_html, fixed = TRUE))
  expect_true(grepl("label-halo", conn_html, fixed = TRUE))
  expect_true(grepl('id="rightMaxEdges" type="number" min="0" value="120"', conn_html, fixed = TRUE))
  expect_true(grepl("[display=none],[visibility=hidden]", conn_html, fixed = TRUE))
  expect_true(grepl("setAttribute('version','1.1')", direct_html, fixed = TRUE))
  expect_true(grepl("setAttribute('version','1.1')", conn_html, fixed = TRUE))
  expect_false(dir.exists(file.path(out_dir, "pdf")))
})
