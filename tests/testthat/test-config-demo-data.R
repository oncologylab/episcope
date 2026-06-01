test_that("demo data registry is empty when no public bundle is configured", {
  info <- craftgrn_demo_data_info()

  expect_s3_class(info, "data.frame")
  expect_equal(nrow(info), 0L)
  expect_named(
    info,
    c(
      "name", "title", "version", "release_tag", "file", "project_dir",
      "url", "md5", "size", "description"
    )
  )
})

test_that("download reports no configured demo data", {
  expect_error(
    download_craftgrn_demo_data(destdir = tempdir(), verbose = FALSE),
    "No external CraftGRN demo data bundle is currently configured"
  )
})
