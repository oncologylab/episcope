test_that("Module 1 preparation writes quantile-normalized FP score CSV", {
  out_dir <- tempfile("craftgrn-module1-prep-output-")
  grn_set <- list(
    fp_score_condition_qn = tibble::tibble(
      peak_ID = c("fp1", "fp2"),
      sample_a = c(1.5, 2.5),
      sample_b = c(3.5, 4.5)
    )
  )

  out_path <- .write_fp_score_qn_csv(
    grn_set = grn_set,
    out_dir = out_dir,
    db = "JASPAR2024",
    verbose = FALSE
  )

  expect_equal(basename(out_path), "01_fp_scores_qn_JASPAR2024.csv")
  expect_true(file.exists(out_path))
  written <- data.table::fread(out_path, showProgress = FALSE)
  expect_equal(as.data.frame(written), as.data.frame(grn_set$fp_score_condition_qn))
})

test_that("condition matrix helpers preserve data.table metadata columns", {
  metadata <- data.table::data.table(
    id = c("sample_a", "sample_b"),
    condition_id = c("condition_a", "condition_b")
  )
  rna <- tibble::tibble(
    ensembl_gene_id = c("gene_1", "gene_2"),
    HGNC = c("Gene1", "Gene2"),
    sample_a = c(2, 0),
    sample_b = c(0, 3)
  )
  fp <- tibble::tibble(
    peak_ID = c("peak_1", "peak_2"),
    sample_a = c(1, 2),
    sample_b = c(3, 4)
  )

  rna_condition <- make_rna_condition(rna, metadata, "condition_id")
  rna_expressed <- make_rna_expressed(
    rna,
    metadata,
    "condition_id",
    threshold_gene_expr = 1
  )
  fp_condition <- make_fp_score_condition(fp, metadata, "condition_id")
  fp_bound <- make_fp_bound_condition(
    fp_bound_tbl = dplyr::mutate(fp, dplyr::across(-"peak_ID", ~ 1L)),
    fp_score_tbl = fp,
    atac_overlap_tbl = NULL,
    fp_annotation_tbl = NULL,
    metadata = metadata,
    label_col = "condition_id",
    threshold_fp_score = 0
  )

  expected_names <- c("condition_a", "condition_b")
  expect_equal(names(rna_condition)[-c(1L, 2L)], expected_names)
  expect_equal(names(rna_expressed)[-c(1L, 2L)], expected_names)
  expect_equal(names(fp_condition)[-1L], expected_names)
  expect_equal(names(fp_bound)[-1L], expected_names)
})
