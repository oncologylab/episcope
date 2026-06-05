source_file <- function(...) {
  candidates <- c(
    file.path(testthat::test_path(), "..", "..", ...),
    file.path(getwd(), ...)
  )
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates)) candidates[[1L]] else NA_character_
}

test_that("pkgdown home relies on DESCRIPTION BugReports link only", {
  path <- source_file("_pkgdown.yml")
  testthat::skip_if_not(file.exists(path))
  cfg <- yaml::read_yaml(path)
  home_links <- if (is.null(cfg$home$links)) list() else cfg$home$links
  link_text <- vapply(home_links, function(x) {
    if (is.null(x$text)) "" else as.character(x$text)
  }, character(1L))
  expect_false(any(link_text == "Report a bug"))
})

test_that("pkgdown workflow fixes generated favicon markup and cleans stale pages", {
  path <- source_file(".github", "workflows", "pkgdown.yaml")
  testthat::skip_if_not(file.exists(path))
  workflow <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(workflow, "Fix generated favicon MIME type", fixed = TRUE)
  expect_match(workflow, "Remove local agent instruction pages", fixed = TRUE)
  expect_match(workflow, "AGENTS.html", fixed = TRUE)
  expect_match(workflow, "search.json", fixed = TRUE)
  expect_match(workflow, "sitemap.xml", fixed = TRUE)
  expect_match(workflow, "image/svg+xml", fixed = TRUE)
  expect_match(workflow, "clean: true", fixed = TRUE)
})

test_that("pkgdown reference uses publication-facing Module 1 and Module 2 names", {
  path <- source_file("_pkgdown.yml")
  testthat::skip_if_not(file.exists(path))
  cfg <- yaml::read_yaml(path)
  sections <- cfg$reference
  titles <- vapply(sections, function(x) as.character(x$title), character(1L))

  module1 <- sections[[which(titles == "TF binding site prediction")]]
  module2 <- sections[[which(titles == "Predict TF targets")]]

  expect_equal(module1$contents, c(
    "predict_tfbs",
    "build_module1_qc_report",
    "module1_prepare_tfbs_inputs",
    "module1_correlate_TF_to_canonical_tfbs",
    "module1_filter_canonical_bound_tfbs",
    "module1_predict_full_tfbs",
    "output_predicted_tfbs",
    "load_predicted_tfbs",
    "export_predicted_tfbs_bed"
  ))
  expect_equal(module2$contents, c(
    "predict_tf_targets",
    "build_module2_qc_report",
    "module2_identify_candidate_links",
    "module2_correlate_tf_targets",
    "module2_link_fp_targets",
    "module2_correlate_fp_targets",
    "output_predicted_links",
    "load_predicted_links",
    "query_predicted_links",
    "check_module2_links",
    "export_tf_target_bedpe",
    "report_top_tf_targets",
    "report_direct_tf_tf_regulations",
    "report_tf_tf_coregulations"
  ))
})
