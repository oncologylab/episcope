# tests/testthat/test-utils_tflinks.R

testthat::test_that("species aliases are normalised correctly", {
  norm <- getFromNamespace(".tflink_normalise_species", "episcope")

  testthat::expect_identical(norm("Homo_sapiens"), "Homo_sapiens")
  testthat::expect_identical(norm("hs"), "Homo_sapiens")
  testthat::expect_identical(norm("mm"), "Mus_musculus")

  testthat::expect_error(
    norm("unknown_species"),
    "Unknown species",
    fixed = FALSE
  )
})

testthat::test_that("TFLink base/data dir helpers behave as expected", {
  old_opt <- getOption("episcope.tflink_data_dir")
  on.exit(options(episcope.tflink_data_dir = old_opt), add = TRUE)

  # Use a temp dir as a fake local episcope-data/TFLink
  tmp_base <- file.path(tempdir(), "TFLink_test_base")
  dir.create(tmp_base, recursive = TRUE, showWarnings = FALSE)

  episcope::tflink_set_data_dir(tmp_base)
  testthat::expect_identical(episcope::tflink_data_dir(), normalizePath(tmp_base, winslash = "/", mustWork = TRUE))

  # No species yet
  testthat::expect_identical(episcope::tflink_list_species(), character())
})

testthat::test_that("tflink_load reads from local chunk files and supports filters", {
  old_opt <- getOption("episcope.tflink_data_dir")
  on.exit(options(episcope.tflink_data_dir = old_opt), add = TRUE)

  # Create fake TFLink chunks under a temp base dir
  tmp_base <- file.path(tempdir(), "TFLink_test_local")
  species  <- "Homo_sapiens"
  species_dir <- file.path(tmp_base, species)
  dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

  # Two small chunks
  chunk1 <- data.frame(
    Name.TF     = c("HNF1A", "HNF4A"),
    Name.Target = c("GENE1", "GENE2"),
    Other       = 1:2,
    stringsAsFactors = FALSE
  )
  chunk2 <- data.frame(
    Name.TF     = c("HNF1A", "FOXA2"),
    Name.Target = c("GENE3", "GENE4"),
    Other       = 3:4,
    stringsAsFactors = FALSE
  )

  saveRDS(chunk1, file = file.path(species_dir, "TFLink_Homo_sapiens_chunk_001.rds"))
  saveRDS(chunk2, file = file.path(species_dir, "TFLink_Homo_sapiens_chunk_002.rds"))

  # Point episcope at this local "TFLink" base
  episcope::tflink_set_data_dir(tmp_base)

  # 1) Full load, no filters
  tbl_all <- episcope::tflink_load("Homo_sapiens", use_cache = FALSE)
  testthat::expect_true(is.data.frame(tbl_all))
  testthat::expect_identical(nrow(tbl_all), 4L)
  testthat::expect_true(all(c("Name.TF", "Name.Target") %in% names(tbl_all)))

  # 2) Filter by TF
  tbl_hnf <- episcope::tflink_load(
    species     = "hs",  # alias for Homo_sapiens
    filter_tf   = "HNF1A",
    select_cols = c("Name.TF", "Name.Target"),
    use_cache   = FALSE
  )

  testthat::expect_identical(sort(unique(tbl_hnf$Name.TF)), "HNF1A")
  testthat::expect_true(all(names(tbl_hnf) %in% c("Name.TF", "Name.Target")))
  testthat::expect_true(all(tbl_hnf$Name.Target %in% c("GENE1", "GENE3")))

  # 3) Filter by target
  tbl_gene2 <- episcope::tflink_load(
    species       = "Homo_sapiens",
    filter_target = "GENE2",
    use_cache     = FALSE
  )
  testthat::expect_identical(unique(tbl_gene2$Name.Target), "GENE2")

  # 4) Bad select_cols triggers an error
  testthat::expect_error(
    episcope::tflink_load(
      species     = "Homo_sapiens",
      select_cols = "NON_EXISTENT_COLUMN",
      use_cache   = FALSE
    ),
    "None of the requested columns present",
    fixed = FALSE
  )
})

testthat::test_that("tflink_load caching and tflink_clear_cache work", {
  old_opt <- getOption("episcope.tflink_data_dir")
  on.exit(options(episcope.tflink_data_dir = old_opt), add = TRUE)

  tmp_base   <- file.path(tempdir(), "TFLink_test_cache")
  species    <- "Homo_sapiens"
  species_dir <- file.path(tmp_base, species)
  dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

  chunk <- data.frame(
    Name.TF     = c("HNF1A", "HNF4A"),
    Name.Target = c("GENE1", "GENE2"),
    stringsAsFactors = FALSE
  )
  chunk_path <- file.path(species_dir, "TFLink_Homo_sapiens_chunk_001.rds")
  saveRDS(chunk, file = chunk_path)

  episcope::tflink_set_data_dir(tmp_base)
  episcope::tflink_clear_cache()

  # First call populates cache
  tbl1 <- episcope::tflink_load("Homo_sapiens", use_cache = TRUE)
  testthat::expect_identical(nrow(tbl1), 2L)

  # Modify file on disk
  chunk2 <- data.frame(
    Name.TF     = c("HNF1A", "HNF4A", "FOXA2"),
    Name.Target = c("GENE1", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )
  saveRDS(chunk2, file = chunk_path)

  # Because of cache, second call should still see the old 2-row table
  tbl2 <- episcope::tflink_load("Homo_sapiens", use_cache = TRUE)
  testthat::expect_identical(nrow(tbl2), 2L)

  # After clearing cache, it should see the updated 3 rows
  episcope::tflink_clear_cache()
  tbl3 <- episcope::tflink_load("Homo_sapiens", use_cache = TRUE)
  testthat::expect_identical(nrow(tbl3), 3L)
})

testthat::test_that("tflink_load fails early for unknown species", {
  # Unknown alias should fail before any file / network access
  testthat::expect_error(
    episcope::tflink_load("xx"),
    "Unknown species",
    fixed = FALSE
  )
})
