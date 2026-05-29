test_that("compact Module 3 bridge writes topic-compatible filtered links", {
  skip_if_not_installed("data.table")

  tmp <- withr::local_tempdir()
  cond <- c("Case", "Ctrl")
  fp_score <- matrix(c(5, 1, 1, 5), nrow = 2, byrow = TRUE, dimnames = list(c("fp_up", "fp_down"), cond))
  fp_bound <- fp_score > 0
  gene_expr <- matrix(c(20, 20, 40, 5, 5, 40), nrow = 3, byrow = TRUE, dimnames = list(c("TF1", "GENE_UP", "GENE_DOWN"), cond))
  gene_on <- gene_expr > 0
  multiomic <- list(
    schema = "craftgrn_multiomic_v1",
    project = list(),
    samples = data.frame(condition_id = cond),
    features = list(
      fp = data.frame(fp_id = rownames(fp_score)),
      gene = data.frame(gene_id = rownames(gene_expr))
    ),
    matrices = list(
      fp_score = fp_score,
      fp_bound = fp_bound,
      gene_expr = gene_expr,
      gene_on = gene_on
    ),
    refs = list(),
    qc = list(),
    paths = list()
  )
  class(multiomic) <- c("craftgrn_multiomic", "list")

  link_dt <- data.table::data.table(
    link_id = c("link_1", "link_2"),
    tf = "TF1",
    fp_id = c("fp_up", "fp_down"),
    target_gene = c("GENE_UP", "GENE_DOWN"),
    candidate_id = c("cand_1", "cand_2"),
    module2_link_pass = TRUE
  )
  cand_dt <- data.table::data.table(
    candidate_id = c("cand_1", "cand_2"),
    fp_id = c("fp_up", "fp_down"),
    target_gene = c("GENE_UP", "GENE_DOWN"),
    distance_to_tss = c(10, -10)
  )
  fp_corr <- data.table::data.table(
    fp_id = c("fp_up", "fp_down"),
    target_gene = c("GENE_UP", "GENE_DOWN"),
    best_r = c(0.5, 0.6),
    best_p = c(0.01, 0.02),
    best_fdr = c(0.01, 0.02),
    pass = TRUE
  )
  tf_corr <- data.table::data.table(
    tf = "TF1",
    target_gene = c("GENE_UP", "GENE_DOWN"),
    best_r = c(0.7, 0.8),
    best_p = c(0.01, 0.02),
    best_fdr = c(0.01, 0.02),
    pass = TRUE
  )

  data.table::fwrite(link_dt, file.path(tmp, "links.csv"))
  data.table::fwrite(cand_dt, file.path(tmp, "candidates.csv"))
  data.table::fwrite(fp_corr, file.path(tmp, "fp_corr.csv"))
  data.table::fwrite(tf_corr, file.path(tmp, "tf_corr.csv"))
  data.table::fwrite(data.table::data.table(chunk_id = 1L, path = file.path(tmp, "links.csv"), format = "csv", n_rows = 2L), file.path(tmp, "module2_links_manifest.csv"))
  data.table::fwrite(data.table::data.table(chunk_id = 1L, path = file.path(tmp, "candidates.csv"), format = "csv", n_rows = 2L), file.path(tmp, "module2_fp_target_candidates_manifest.csv"))
  data.table::fwrite(data.table::data.table(chunk_id = 1L, path = file.path(tmp, "fp_corr.csv"), format = "csv", n_rows = 2L), file.path(tmp, "module2_fp_target_corr_manifest.csv"))
  manifest <- data.table::data.table(
    table = c("module2_links", "module2_fp_target_candidates", "module2_fp_target_corr", "module2_tf_target_corr"),
    path = c(file.path(tmp, "module2_links_manifest.csv"), file.path(tmp, "module2_fp_target_candidates_manifest.csv"), file.path(tmp, "module2_fp_target_corr_manifest.csv"), file.path(tmp, "tf_corr.csv")),
    format = c("manifest", "manifest", "manifest", "csv"),
    n_rows = c(2L, 2L, 2L, 2L)
  )
  data.table::fwrite(manifest, file.path(tmp, "module2_manifest.csv"))

  out_dir <- file.path(tmp, "diff_links_filtered")
  res <- module3_prepare_differential_links(
    module2 = tmp,
    multiomic_data = multiomic,
    compar = data.frame(cond1_label = "Case", cond2_label = "Ctrl"),
    project_config = list(fp_filter_mode = "log2fc", fp_log2fc_cutoff = 1, gene_log2fc_cutoff = 1, threshold_expr = 0, threshold_fp_score = 0),
    output_dir = out_dir,
    n_cores = 1,
    overwrite = TRUE,
    verbose = FALSE
  )
  expect_equal(nrow(res), 1L)
  up <- data.table::fread(res$up_path[[1]])
  down <- data.table::fread(res$down_path[[1]])
  expect_equal(up$gene_key, "GENE_UP")
  expect_equal(down$gene_key, "GENE_DOWN")
  expect_true(all(c("tf", "gene_key", "peak_id", "log2FC_fp_score", "log2FC_gene_expr", "r_gene", "r_rna_gene", "distance_to_tss") %in% names(up)))
  expect_true(standardize_delta_links_one(res$up_path[[1]], keep_original = FALSE)$comparison_id[[1]] == "Case_vs_Ctrl")
})

test_that("Module 3 link_padded FP signal mode pads inactive condition scores", {
  skip_if_not_installed("data.table")

  cond <- c("Case", "Ctrl")
  fp_score <- matrix(c(5, 10), nrow = 1, dimnames = list("fp1", cond))
  fp_bound <- matrix(c(0, 1), nrow = 1, dimnames = list("fp1", cond))
  gene_expr <- matrix(c(20, 20, 5, 40), nrow = 2, byrow = TRUE, dimnames = list(c("TF1", "GENE1"), cond))
  gene_on <- gene_expr > 0
  multiomic <- list(
    matrices = list(
      fp_score = fp_score,
      fp_bound = fp_bound,
      gene_expr = gene_expr,
      gene_on = gene_on
    )
  )
  link_dt <- data.table::data.table(
    tf = "TF1",
    fp_id = "fp1",
    target_gene = "GENE1",
    r_gene = 0.5,
    p_gene = 0.01,
    p_adj_gene = 0.01,
    r_rna_gene = 0.7,
    p_rna_gene = 0.01,
    p_rna_adj_gene = 0.01
  )

  actual <- .module3_build_chunk_delta_prepared(
    link_dt,
    multiomic,
    "Case",
    "Ctrl",
    pseudocount = 1,
    fp_signal_mode = "actual"
  )
  padded <- .module3_build_chunk_delta_prepared(
    link_dt,
    multiomic,
    "Case",
    "Ctrl",
    pseudocount = 1,
    fp_signal_mode = "link_padded"
  )

  expect_equal(actual$fp_score_Case, 5)
  expect_equal(padded$fp_score_Case, 0)
  expect_equal(padded$fp_score_Ctrl, 10)
  expect_equal(padded$delta_fp_score, -10)
  expect_equal(padded$fp_signal_mode, "link_padded")
  expect_true(padded$active_Ctrl)
  expect_false(padded$active_Case)
})
