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

  module1 <- sections[[which(titles == "TF binding site and motif utilities")]]
  module1_aux <- sections[[which(titles == "TFBS auxiliary functions")]]
  module2 <- sections[[which(titles == "Predict TF targets")]]
  module2_aux <- sections[[which(titles == "TF target auxiliary functions")]]

  expect_equal(module1$contents, c(
    "predict_tfbs",
    "module1_prepare_tfbs_inputs",
    "module1_correlate_TF_to_canonical_tfbs",
    "module1_filter_canonical_bound_tfbs",
    "module1_predict_full_tfbs",
    "build_module1_qc_report"
  ))
  expect_equal(module1_aux$contents, c(
    "output_predicted_tfbs",
    "load_predicted_tfbs",
    "export_predicted_tfbs_bed"
  ))
  expect_equal(module2$contents, c(
    "predict_tf_targets",
    "module2_identify_candidate_links",
    "module2_correlate_tf_targets",
    "module2_link_fp_targets",
    "module2_correlate_fp_targets",
    "module2_output_predicted_links",
    "build_module2_qc_report"
  ))
  expect_equal(module2_aux$contents, c(
    "load_predicted_links",
    "query_predicted_links",
    "check_predicted_links",
    "export_tf_target_bedpe",
    "report_top_tf_targets",
    "report_direct_tf_tf_regulations",
    "report_tf_tf_coregulations"
  ))
})

test_that("pkgdown reference uses clean Module 3 public API names", {
  path <- source_file("_pkgdown.yml")
  testthat::skip_if_not(file.exists(path))
  cfg <- yaml::read_yaml(path)
  sections <- cfg$reference
  titles <- vapply(sections, function(x) as.character(x$title), character(1L))

  module3 <- sections[[which(titles == "Regulatory topic modeling")]]
  vis <- sections[[which(titles == "Module 3 visualization utilities")]]

  expect_equal(module3$contents, c(
    "run_topic_modeling",
    "module3_prepare_differential_links",
    "module3_construct_docs",
    "module3_train_topic_models",
    "module3_extract_topics",
    "build_module3_qc_report"
  ))
  expect_equal(vis$contents, c(
    "visualize_topic_modeling_results",
    "visualize_differential_grns"
  ))

  all_contents <- unlist(lapply(sections, function(x) x$contents), use.names = FALSE)
  expect_false("run_regulatory_topics" %in% all_contents)
  expect_false("module3_prepare_topic_inputs" %in% all_contents)
  expect_false("report_differential_pathways" %in% all_contents)
  expect_false("report_master_tfs" %in% all_contents)
  expect_false("load_differential_links" %in% all_contents)
  expect_false("query_differential_links" %in% all_contents)
  expect_false("run_module3_topic_benchmark" %in% all_contents)
  expect_false("train_topic_models" %in% all_contents)
  expect_false("extract_regulatory_topics" %in% all_contents)
  expect_false("load_delta_links_many" %in% all_contents)
  expect_false("build_tf_cluster_map_from_motif" %in% all_contents)
  expect_false("default_diff_grn_pathway_databases" %in% all_contents)
  expect_false("run_diff_grn_pathway_analysis" %in% all_contents)
  expect_false("run_diff_grn_master_tf_summary" %in% all_contents)
})
