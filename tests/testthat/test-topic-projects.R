topic_benchmark_projects_script <- normalizePath(file.path(
  Sys.getenv(
    "CRAFTGRN_REPO_DIR",
    unset = normalizePath(file.path(getwd(), "../.."), mustWork = FALSE)
  ),
  "dev/benchmark/03_topic_benchmark_projects.R"
), mustWork = FALSE)

source_benchmark_script <- function(path) {
  testthat::skip_if_not(file.exists(path), message = "dev benchmark script is not included in built package")
  old_wd <- setwd(dirname(dirname(dirname(path))))
  on.exit(setwd(old_wd), add = TRUE)
  source(path, local = parent.frame())
}

test_that("gene-only topic extractions still require pathway outputs", {
  source_benchmark_script(topic_benchmark_projects_script)

  extraction_dir <- file.path(tempdir(), "model_K12_thrP0p9_link_gene_gammafit")
  sub_dir <- file.path(extraction_dir, "report")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(sub_dir, c(
    "topic_links.csv",
    "topic_terms.csv",
    "topic_marker_term_heatmap.pdf",
    "HPAFII_topic_by_comparison.pdf",
    "report_doc_design.txt"
  )))

  expect_false(.extract_complete(extraction_dir))
  status <- .extraction_step_status(extraction_dir)
  expect_true(status$pathway_required)
  expect_false(status$pathway_complete_or_not_required)

  file.create(file.path(sub_dir, "topic_pathway_enrichment_gene_only_dotplot.pdf"))
  expect_true(.extract_complete(extraction_dir))
  status <- .extraction_step_status(extraction_dir)
  expect_true(status$pathway_complete_or_not_required)
})
