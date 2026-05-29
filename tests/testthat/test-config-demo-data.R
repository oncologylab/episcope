test_that("demo data registry exposes reproducible download metadata", {
  info <- craftgrn_demo_data_info()

  expect_s3_class(info, "data.frame")
  expect_equal(info$name, "nutrient_stress_chr22")
  expect_equal(info$release_tag, "demo-data-v0.1.0")
  expect_equal(info$file, "nutrient_stress_strict_JASPAR2024_chr22_demo_inputs.tar.gz")
  expect_equal(info$project_dir, "nutrient_stress_strict_JASPAR2024_chr22_demo")
  expect_equal(info$md5, "fa754783186b0e5119387b0405652331")
  expect_match(info$url, "https://github.com/oncologylab/craftgrn/releases/download/demo-data-v0.1.0/", fixed = TRUE)
})
