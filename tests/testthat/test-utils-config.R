test_that("load_config resolves relative paths against the config directory", {
  project_dir <- tempfile("craftgrn-project-")
  dir.create(file.path(project_dir, "data"), recursive = TRUE)
  config_path <- file.path(project_dir, "project.yaml")
  writeLines(
    c(
      "base_dir: \".\"",
      "ref_genome: hg38",
      "db: JASPAR2024",
      "threshold_expr: 10",
      "threshold_fp_score: 2",
      "threshold_fp_tf_corr_r: 0.3",
      "link_window_bp: 100000",
      "threshold_rna_gene_corr_r: 0.3",
      "threshold_rna_gene_corr_p: null",
      "threshold_rna_gene_corr_fdr: null",
      "threshold_fp_gene_corr_r: 0.3",
      "threshold_fp_gene_corr_p: null",
      "threshold_fp_gene_corr_fdr: null",
      "sample_metadata: data/sample_metadata_strict.csv",
      "atac_master: data/atac_master.csv",
      "rna_mapped: data/rna_strict.csv"
    ),
    config_path
  )

  env <- new.env(parent = emptyenv())
  load_config(config_path, env = env)

  expect_equal(get("base_dir", envir = env), normalizePath(project_dir, winslash = "/", mustWork = FALSE))
  expect_equal(
    get("sample_metadata", envir = env),
    normalizePath(file.path(project_dir, "data/sample_metadata_strict.csv"), winslash = "/", mustWork = FALSE)
  )
  expect_equal(
    get("atac_master", envir = env),
    normalizePath(file.path(project_dir, "data/atac_master.csv"), winslash = "/", mustWork = FALSE)
  )
  expect_equal(
    get("rna_mapped", envir = env),
    normalizePath(file.path(project_dir, "data/rna_strict.csv"), winslash = "/", mustWork = FALSE)
  )
})
