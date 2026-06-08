read_project_file <- function(path, marker) {
  candidates <- c(
    file.path(testthat::test_path(), "..", "..", path),
    file.path(testthat::test_path(), "..", "..", "..", path),
    file.path(getwd(), path)
  )
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  candidates <- candidates[file.exists(candidates)]
  contents <- vapply(candidates, function(candidate) {
    paste(readLines(candidate, warn = FALSE), collapse = "\n")
  }, character(1L))
  hit <- contents[grepl(marker, contents, fixed = TRUE)]
  testthat::skip_if_not(length(hit) > 0L)
  hit[[1L]]
}

test_that("README version badge tracks package metadata dynamically", {
  readme <- read_project_file("README.md", "# CraftGRN")
  expect_false(grepl("version-0.1.2", readme, fixed = TRUE))
  expect_match(
    readme,
    "https://img.shields.io/github/r-package/v/oncologylab/craftgrn",
    fixed = TRUE
  )
})

test_that("README and NEWS do not advertise retired public API names", {
  readme <- read_project_file("README.md", "# CraftGRN")
  news <- read_project_file("NEWS.md", "# craftgrn")
  retired <- c(
    "link_tf_targets",
    "run_regulatory_topics",
    "module3_prepare_topic_inputs",
    "train_topic_models",
    "extract_regulatory_topics",
    "run_module3_topic_benchmark",
    "run_diff_grn_pathway_analysis",
    "run_diff_grn_master_tf_summary"
  )
  for (fn in retired) {
    expect_false(grepl(paste0("`", fn, "\\(\\)`"), readme))
    expect_false(grepl(paste0("`", fn, "\\(\\)`"), news))
  }
})

test_that("README CraftGRN function references are exported or explicitly external", {
  readme <- read_project_file("README.md", "# CraftGRN")
  refs <- unique(unlist(regmatches(readme, gregexpr("`[A-Za-z0-9_.:]+\\(\\)`", readme))))
  refs <- sub("^`", "", sub("\\(\\)`$", "", refs))
  refs <- sub("^craftgrn::", "", refs)
  external <- c(
    "remotes::install_github",
    "pak::pak",
    "install.packages",
    "BiocManager::install"
  )
  local_refs <- setdiff(refs, external)
  expect_true(all(local_refs %in% getNamespaceExports("craftgrn")))
})
