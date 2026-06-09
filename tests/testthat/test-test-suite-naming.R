test_that("testthat files use approved kebab-case names", {
  files <- list.files(testthat::test_path(), pattern = "^(test|helper)-.*[.]R$", full.names = FALSE)
  expect_false(any(grepl("_", files, fixed = TRUE)))
  expect_true(all(grepl("^(test|helper)-[a-z0-9]+(-[a-z0-9]+)*[.]R$", files)))

  test_files <- grep("^test-", files, value = TRUE)
  name <- sub("[.]R$", "", sub("^test-", "", test_files))
  area <- vapply(strsplit(name, "-", fixed = TRUE), function(x) {
    if (length(x) >= 2L && identical(paste(x[1:2], collapse = "-"), "test-suite")) {
      "test-suite"
    } else {
      x[[1L]]
    }
  }, character(1L))
  allowed <- c("app", "benchmark", "config", "module1", "module2", "module3", "test-suite", "topic")
  expect_true(all(area %in% allowed))
})

test_that("test case names describe package behavior", {
  files <- list.files(testthat::test_path(), pattern = "^test-.*[.]R$", full.names = TRUE)
  lines <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
  descriptions <- sub('.*test_that[(]\"([^\"]+)\".*', "\\1", grep('test_that[(]\"', lines, value = TRUE))
  expect_false(any(grepl("^Test ", descriptions)))
  expect_false(any(grepl("^multiplication works$|^.* works$", descriptions)))
  expect_false(any(grepl("obsolete|deprecated|legacy placeholder", descriptions, ignore.case = TRUE)))
})

test_that("package tests do not require text2vec", {
  files <- list.files(testthat::test_path(), pattern = "^[a-z0-9-]+[.]R$", full.names = TRUE)
  lines <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
  package_name <- "text2vec"
  forbidden <- c(
    paste0("skip_if_not_installed(\"", package_name, "\")"),
    paste0("asNamespace(\"", package_name, "\")"),
    paste0("requireNamespace(\"", package_name, "\""),
    paste0(package_name, "::"),
    paste0("library(", package_name, ")"),
    paste0("require(", package_name, ")")
  )
  hits <- unlist(lapply(forbidden, grep, x = lines, value = TRUE, fixed = TRUE), use.names = FALSE)
  expect_length(hits, 0L)
})

test_that("package DESCRIPTION does not depend on text2vec", {
  desc <- as.list(utils::packageDescription("craftgrn"))
  dependency_fields <- intersect(c("Depends", "Imports", "Suggests", "Enhances", "LinkingTo"), names(desc))
  dependency_text <- paste(unlist(desc[dependency_fields], use.names = FALSE), collapse = "\n")
  expect_false(grepl("text2vec", dependency_text, fixed = TRUE))
})

test_that("Shiny app packages stay optional dependencies", {
  desc <- as.list(utils::packageDescription("craftgrn"))
  imports <- desc$Imports %||% ""
  suggests <- desc$Suggests %||% ""

  expect_false(grepl("shiny", imports, fixed = TRUE))
  expect_false(grepl("golem", imports, fixed = TRUE))
  expect_true(grepl("shiny", suggests, fixed = TRUE))
  expect_true(grepl("golem", suggests, fixed = TRUE))

  namespace_file <- system.file("NAMESPACE", package = "craftgrn")
  expect_true(nzchar(namespace_file))
  namespace_lines <- readLines(namespace_file, warn = FALSE)
  expect_false(any(grepl("import(shiny)", namespace_lines, fixed = TRUE)))
  expect_false(any(grepl("importFrom(shiny", namespace_lines, fixed = TRUE)))
  expect_false(any(grepl("importFrom(golem", namespace_lines, fixed = TRUE)))
})

test_that("Module 3 report packages install with the package", {
  desc <- as.list(utils::packageDescription("craftgrn"))
  imports <- desc$Imports %||% ""
  suggests <- desc$Suggests %||% ""
  required_report_packages <- c("pheatmap", "enrichR", "LDAvis")

  for (pkg in required_report_packages) {
    expect_true(grepl(pkg, imports, fixed = TRUE))
    expect_false(grepl(pkg, suggests, fixed = TRUE))
  }
})
