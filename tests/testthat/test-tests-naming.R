test_that("testthat files use kebab-case names", {
  files <- list.files(testthat::test_path(), pattern = "^(test|helper)-.*[.]R$", full.names = FALSE)
  expect_false(any(grepl("_", files, fixed = TRUE)))
  expect_true(all(grepl("^(test|helper)-[a-z0-9]+(-[a-z0-9]+)*[.]R$", files)))
})
