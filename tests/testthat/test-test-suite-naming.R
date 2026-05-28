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
