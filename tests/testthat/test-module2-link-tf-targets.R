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
  expect_false(dir.exists(file.path(out_dir, "pdf")))
})
