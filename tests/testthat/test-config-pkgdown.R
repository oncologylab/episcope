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
  expect_match(workflow, "image/svg+xml", fixed = TRUE)
  expect_match(workflow, "clean: true", fixed = TRUE)
})
