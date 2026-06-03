test_that("Module 1 preprocessing defaults use available cores", {
  expect_gt(.available_cores(logical = TRUE), 0L)
  expect_match(paste(deparse(formals(load_footprints)$n_workers), collapse = " "), "\\.available_cores")
  expect_match(paste(deparse(formals(align_footprints)$threads), collapse = " "), "\\.available_cores")
  expect_true(eval(formals(align_footprints)$write_fp_sites))
})

test_that("aligned footprint cache loader skips id_map by default", {
  cache_dir <- file.path(tempdir(), paste0("craftgrn-fp-cache-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_tag <- "TEST"

  readr::write_csv(
    tibble::tibble(peak_ID = c("chr1:1-10", "chr1:1-10", "chr1:20-30"), sample_a = c(1, 1, 2)),
    file.path(cache_dir, sprintf("fp_scores_%s.csv", cache_tag))
  )
  readr::write_csv(
    tibble::tibble(peak_ID = c("chr1:1-10", "chr1:1-10", "chr1:20-30"), sample_a = c(1L, 1L, 0L)),
    file.path(cache_dir, sprintf("fp_bounds_%s.csv", cache_tag))
  )
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("chr1:1-10", "chr1:20-30"),
      atac_peak = c("chr1:1-50", "chr1:1-50"),
      motifs_all = c("M1", "M2")
    ),
    file.path(cache_dir, sprintf("fp_sites_%s.csv", cache_tag))
  )
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("chr1:1-10", "chr1:20-30"),
      source_fp_peak = c("chr1:1-10", "chr1:20-30"),
      atac_peak = c("chr1:1-50", "chr1:1-50"),
      group_size = c(1, 1)
    ),
    file.path(cache_dir, sprintf("fp_id_map_%s.csv", cache_tag))
  )

  out <- load_fp_aligned_from_cache(
    cache_dir = cache_dir,
    cache_tag = cache_tag,
    output_mode = "distinct",
    verbose = FALSE
  )

  expect_equal(nrow(out$fp_score), 2L)
  expect_equal(nrow(out$fp_bound), 2L)
  expect_equal(nrow(out$fp_annotation), 2L)
  expect_equal(nrow(out$id_map), 0L)
})

test_that("aligned footprint cache loader loads id_map only when requested", {
  cache_dir <- file.path(tempdir(), paste0("craftgrn-fp-cache-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_tag <- "TEST"

  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1), file.path(cache_dir, sprintf("fp_scores_%s.csv", cache_tag)))
  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1L), file.path(cache_dir, sprintf("fp_bounds_%s.csv", cache_tag)))
  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", atac_peak = "chr1:1-50", motifs_all = "M1"), file.path(cache_dir, sprintf("fp_sites_%s.csv", cache_tag)))
  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", source_fp_peak = "chr1:1-10", atac_peak = "chr1:1-50", group_size = 1), file.path(cache_dir, sprintf("fp_id_map_%s.csv", cache_tag)))

  out <- load_fp_aligned_from_cache(
    cache_dir = cache_dir,
    cache_tag = cache_tag,
    output_mode = "distinct",
    load_id_map = TRUE,
    verbose = FALSE
  )

  expect_equal(nrow(out$id_map), 1L)
  expect_equal(out$id_map$group_size, 1)
})

test_that("aligned footprint cache loader reads Parquet cache when available", {
  testthat::skip_if_not_installed("arrow")
  cache_dir <- file.path(tempdir(), paste0("craftgrn-fp-cache-parquet-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_tag <- "TEST"

  arrow::write_parquet(
    tibble::tibble(peak_ID = c("chr1:1-10", "chr1:20-30"), sample_a = c(1, 2)),
    file.path(cache_dir, sprintf("fp_scores_%s.parquet", cache_tag))
  )
  arrow::write_parquet(
    tibble::tibble(peak_ID = c("chr1:1-10", "chr1:20-30"), sample_a = c(1L, 0L)),
    file.path(cache_dir, sprintf("fp_bounds_%s.parquet", cache_tag))
  )
  arrow::write_parquet(
    tibble::tibble(
      peak_ID = c("chr1:1-10", "chr1:20-30"),
      atac_peak = c("chr1:1-50", "chr1:1-50"),
      motifs_all = c("M1", "M2")
    ),
    file.path(cache_dir, sprintf("fp_sites_%s.parquet", cache_tag))
  )

  out <- load_fp_aligned_from_cache(
    cache_dir = cache_dir,
    cache_tag = cache_tag,
    output_mode = "distinct",
    cache_format = "auto",
    verbose = FALSE
  )

  expect_equal(nrow(out$fp_score), 2L)
  expect_equal(out$fp_annotation$motifs, c("M1", "M2"))
  expect_equal(nrow(out$id_map), 0L)
})

test_that("aligned footprint cache loader defaults to CSV when both formats exist", {
  testthat::skip_if_not_installed("arrow")
  cache_dir <- file.path(tempdir(), paste0("craftgrn-fp-cache-default-csv-", as.integer(stats::runif(1L, 1, 1e9))))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_tag <- "TEST"

  readr::write_csv(
    tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1),
    file.path(cache_dir, sprintf("fp_scores_%s.csv", cache_tag))
  )
  readr::write_csv(
    tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1L),
    file.path(cache_dir, sprintf("fp_bounds_%s.csv", cache_tag))
  )
  readr::write_csv(
    tibble::tibble(peak_ID = "chr1:1-10", atac_peak = "chr1:1-50", motifs_all = "CSV"),
    file.path(cache_dir, sprintf("fp_sites_%s.csv", cache_tag))
  )

  arrow::write_parquet(
    tibble::tibble(peak_ID = "chr1:1-10", sample_a = 2),
    file.path(cache_dir, sprintf("fp_scores_%s.parquet", cache_tag))
  )
  arrow::write_parquet(
    tibble::tibble(peak_ID = "chr1:1-10", sample_a = 0L),
    file.path(cache_dir, sprintf("fp_bounds_%s.parquet", cache_tag))
  )
  arrow::write_parquet(
    tibble::tibble(peak_ID = "chr1:1-10", atac_peak = "chr1:1-50", motifs_all = "PARQUET"),
    file.path(cache_dir, sprintf("fp_sites_%s.parquet", cache_tag))
  )

  out <- load_fp_aligned_from_cache(
    cache_dir = cache_dir,
    cache_tag = cache_tag,
    output_mode = "distinct",
    verbose = FALSE
  )

  expect_equal(out$fp_score$sample_a, 1)
  expect_equal(out$fp_bound$sample_a, 1L)
  expect_equal(out$fp_annotation$motifs, "CSV")
})

test_that("footprint alignment can write Parquet cache", {
  testthat::skip_if_not_installed("arrow")
  bench_dir <- file.path(tempdir(), paste0("craftgrn-align-parquet-", as.integer(stats::runif(1L, 1, 1e9))))
  raw_dir <- file.path(bench_dir, "raw")
  cache_dir <- file.path(bench_dir, "cache")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  score_path <- file.path(raw_dir, "M1_score.csv")
  bound_path <- file.path(raw_dir, "M1_bound.csv")
  annot_path <- file.path(raw_dir, "M1_annotation.csv")
  readr::write_csv(tibble::tibble(peak_ID = "chr1:10-20", sample_a = 1.1), score_path)
  readr::write_csv(tibble::tibble(peak_ID = "chr1:10-20", sample_a = 1L), bound_path)
  readr::write_csv(tibble::tibble(fp_peak = "chr1:10-20", atac_peak = "chr1:1-100", motifs = "M1"), annot_path)
  manifest <- tibble::tibble(motif = "M1", n_peaks = 1L, score = score_path, bound = bound_path, annot = annot_path)

  aligned <- align_footprints(
    manifest,
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 1,
    threads = 1L,
    cache_dir = cache_dir,
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = TRUE,
    output_mode = "distinct",
    return_id_map = FALSE,
    write_id_map = FALSE,
    cache_format = "parquet",
    verbose = FALSE
  )
  from_cache <- load_fp_aligned_from_cache(
    cache_dir = cache_dir,
    cache_tag = "TEST",
    output_mode = "distinct",
    cache_format = "parquet",
    verbose = FALSE
  )

  expect_true(file.exists(file.path(cache_dir, "fp_scores_TEST.parquet")))
  expect_false(file.exists(file.path(cache_dir, "fp_scores_TEST.csv")))
  expect_equal(from_cache$fp_score, aligned$fp_score)
  expect_equal(from_cache$fp_annotation[, c("fp_peak", "atac_peak", "motifs")], aligned$fp_annotation[, c("fp_peak", "atac_peak", "motifs")])
})

test_that("raw footprint loader reuses a complete manifest before scanning raw folders", {
  root_dir <- file.path(tempdir(), paste0("craftgrn-missing-root-", as.integer(stats::runif(1L, 1, 1e9))))
  out_dir <- file.path(tempdir(), paste0("craftgrn-fp-raw-cache-", as.integer(stats::runif(1L, 1, 1e9))), "fp_TEST")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  score_path <- file.path(out_dir, "M1_score.csv")
  bound_path <- file.path(out_dir, "M1_bound.csv")
  annot_path <- file.path(out_dir, "M1_annotation.csv")
  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1), score_path)
  readr::write_csv(tibble::tibble(peak_ID = "chr1:1-10", sample_a = 1L), bound_path)
  readr::write_csv(tibble::tibble(fp_peak = "chr1:1-10", atac_peak = "chr1:1-50", motifs = "M1"), annot_path)
  readr::write_csv(
    tibble::tibble(motif = "M1", n_peaks = 1L, score = score_path, bound = bound_path, annot = annot_path),
    file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
  )

  manifest <- load_footprints(
    root_dir = root_dir,
    db_name = "TEST",
    out_dir = out_dir,
    skip_existing = TRUE,
    verbose = FALSE
  )

  expect_true(isTRUE(attr(manifest, "from_cache")))
  expect_equal(manifest$motif, "M1")
})

test_that("raw footprint loader writes minimal per-motif cache tables", {
  root_dir <- file.path(tempdir(), paste0("craftgrn-raw-root-", as.integer(stats::runif(1L, 1, 1e9))))
  out_dir <- file.path(tempdir(), paste0("craftgrn-raw-cache-", as.integer(stats::runif(1L, 1, 1e9))), "fp_TEST")
  for (sid in c("S1", "S2")) {
    dir.create(file.path(root_dir, sid, "TEST", "M1"), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      data.table::data.table(
        TFBS_chr = c("chr1", "chr1"),
        TFBS_start = c(10L, 20L),
        TFBS_end = c(15L, 25L),
        TFBS_name = c("M1", "M1"),
        TFBS_score = c(8, 9),
        TFBS_strand = c("+", "-"),
        peak_chr = c("chr1", "chr1"),
        peak_start = c(1L, 1L),
        peak_end = c(100L, 100L),
        score = c(1.2, 0.1),
        bound = c(1L, 0L)
      ) |>
        stats::setNames(c(
          "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
          "TFBS_score", "TFBS_strand", "peak_chr", "peak_start", "peak_end",
          paste0(sid, "_ATAC_score"), paste0(sid, "_ATAC_bound")
        )),
      file.path(root_dir, sid, "TEST", "M1", "M1_overview.txt"),
      sep = "\t"
    )
  }

  manifest <- load_footprints(
    root_dir = root_dir,
    db_name = "TEST",
    out_dir = out_dir,
    sample_ids = c("S1", "S2"),
    n_workers = 1L,
    verbose = FALSE
  )

  score <- data.table::fread(manifest$score[[1L]], showProgress = FALSE)
  bound <- data.table::fread(manifest$bound[[1L]], showProgress = FALSE)
  annot <- data.table::fread(manifest$annot[[1L]], showProgress = FALSE)

  expect_equal(names(score), c("peak_ID", "S1", "S2"))
  expect_equal(names(bound), c("peak_ID", "S1", "S2"))
  expect_equal(names(annot), c("fp_peak", "atac_peak", "motifs"))
  expect_equal(nrow(score), 1L)
  expect_equal(score$S1, 1.2)
  expect_equal(bound$S2, 1L)
  expect_equal(annot$motifs, "M1")
})

test_that("raw footprint loader uses clear cache status wording", {
  root_dir <- file.path(tempdir(), paste0("craftgrn-log-root-", as.integer(stats::runif(1L, 1, 1e9))))
  out_dir <- file.path(tempdir(), paste0("craftgrn-log-cache-", as.integer(stats::runif(1L, 1, 1e9))), "fp_TEST")
  log_file <- file.path(tempdir(), paste0("craftgrn-log-", as.integer(stats::runif(1L, 1, 1e9)), ".txt"))
  dir.create(file.path(root_dir, "S1", "TEST", "M1"), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      TFBS_chr = "chr1",
      TFBS_start = 10L,
      TFBS_end = 15L,
      TFBS_name = "M1",
      peak_chr = "chr1",
      peak_start = 1L,
      peak_end = 100L,
      S1_ATAC_score = 2.5,
      S1_ATAC_bound = 1L
    ),
    file.path(root_dir, "S1", "TEST", "M1", "M1_overview.txt"),
    sep = "\t"
  )

  load_footprints(
    root_dir = root_dir,
    db_name = "TEST",
    out_dir = out_dir,
    sample_ids = "S1",
    n_workers = 1L,
    skip_existing = TRUE,
    verbose = FALSE,
    log_file = log_file
  )

  log_text <- paste(readLines(log_file), collapse = "\n")
  expect_match(log_text, "preparing cache outputs for this run", fixed = TRUE)
  expect_false(grepl("incomplete", log_text, fixed = TRUE))
  expect_false(grepl("recomputing", log_text, fixed = TRUE))
})

test_that("raw footprint loader recomputes incomplete cached motif triplets", {
  root_dir <- file.path(tempdir(), paste0("craftgrn-resume-root-", as.integer(stats::runif(1L, 1, 1e9))))
  out_dir <- file.path(tempdir(), paste0("craftgrn-resume-cache-", as.integer(stats::runif(1L, 1, 1e9))), "fp_TEST")
  log_file <- file.path(tempdir(), paste0("craftgrn-resume-log-", as.integer(stats::runif(1L, 1, 1e9)), ".txt"))
  dir.create(file.path(root_dir, "S1", "TEST", "M1"), recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      TFBS_chr = "chr1",
      TFBS_start = 10L,
      TFBS_end = 15L,
      TFBS_name = "M1",
      peak_chr = "chr1",
      peak_start = 1L,
      peak_end = 100L,
      S1_ATAC_score = 2.5,
      S1_ATAC_bound = 1L
    ),
    file.path(root_dir, "S1", "TEST", "M1", "M1_overview.txt"),
    sep = "\t"
  )
  readr::write_csv(tibble::tibble(peak_ID = "chr1:999-1000", S1 = 9), file.path(out_dir, "M1_score.csv"))

  manifest <- load_footprints(
    root_dir = root_dir,
    db_name = "TEST",
    out_dir = out_dir,
    sample_ids = "S1",
    n_workers = 1L,
    skip_existing = TRUE,
    verbose = FALSE,
    log_file = log_file
  )
  score <- data.table::fread(manifest$score[[1L]], showProgress = FALSE)
  bound <- data.table::fread(manifest$bound[[1L]], showProgress = FALSE)
  annot <- data.table::fread(manifest$annot[[1L]], showProgress = FALSE)
  log_text <- paste(readLines(log_file), collapse = "\n")

  expect_equal(nrow(manifest), 1L)
  expect_equal(score$peak_ID, "chr1:10-15")
  expect_equal(score$S1, 2.5)
  expect_equal(bound$S1, 1L)
  expect_equal(annot$fp_peak, "chr1:10-15")
  expect_match(log_text, "refreshing partial cache", fixed = TRUE)
  expect_false(grepl("incomplete", log_text, fixed = TRUE))
  expect_false(grepl("recomputing", log_text, fixed = TRUE))
})

test_that("raw footprint loader recomputes cached motif triplets with mismatched rows", {
  root_dir <- file.path(tempdir(), paste0("craftgrn-resume-mismatch-root-", as.integer(stats::runif(1L, 1, 1e9))))
  out_dir <- file.path(tempdir(), paste0("craftgrn-resume-mismatch-cache-", as.integer(stats::runif(1L, 1, 1e9))), "fp_TEST")
  dir.create(file.path(root_dir, "S1", "TEST", "M1"), recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      TFBS_chr = c("chr1", "chr1"),
      TFBS_start = c(10L, 20L),
      TFBS_end = c(15L, 25L),
      TFBS_name = c("M1", "M1"),
      peak_chr = c("chr1", "chr1"),
      peak_start = c(1L, 1L),
      peak_end = c(100L, 100L),
      S1_ATAC_score = c(2.5, 3.5),
      S1_ATAC_bound = c(1L, 1L)
    ),
    file.path(root_dir, "S1", "TEST", "M1", "M1_overview.txt"),
    sep = "\t"
  )
  readr::write_csv(tibble::tibble(peak_ID = c("chr1:10-15", "chr1:20-25"), S1 = c(2.5, 3.5)), file.path(out_dir, "M1_score.csv"))
  readr::write_csv(tibble::tibble(peak_ID = "chr1:10-15", S1 = 1L), file.path(out_dir, "M1_bound.csv"))
  readr::write_csv(tibble::tibble(fp_peak = "chr1:10-15", atac_peak = "chr1:1-100", motifs = "M1"), file.path(out_dir, "M1_annotation.csv"))

  manifest <- load_footprints(
    root_dir = root_dir,
    db_name = "TEST",
    out_dir = out_dir,
    sample_ids = "S1",
    n_workers = 1L,
    skip_existing = TRUE,
    verbose = FALSE
  )
  score <- data.table::fread(manifest$score[[1L]], showProgress = FALSE)
  bound <- data.table::fread(manifest$bound[[1L]], showProgress = FALSE)
  annot <- data.table::fread(manifest$annot[[1L]], showProgress = FALSE)

  expect_equal(nrow(score), 2L)
  expect_equal(nrow(bound), 2L)
  expect_equal(nrow(annot), 2L)
  expect_equal(manifest$n_peaks, 2L)
})

test_that("footprint alignment can skip returning and writing id_map", {
  bench_dir <- file.path(tempdir(), paste0("craftgrn-align-cache-", as.integer(stats::runif(1L, 1, 1e9))))
  raw_dir <- file.path(bench_dir, "raw")
  cache_dir <- file.path(bench_dir, "cache")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  score_path <- file.path(raw_dir, "M1_score.csv")
  bound_path <- file.path(raw_dir, "M1_bound.csv")
  annot_path <- file.path(raw_dir, "M1_annotation.csv")
  readr::write_csv(
    tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1.1, 1.1)),
    score_path
  )
  readr::write_csv(
    tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1L, 1L)),
    bound_path
  )
  readr::write_csv(
    tibble::tibble(
      fp_peak = c("chr1:10-20", "chr1:12-22"),
      atac_peak = c("chr1:1-100", "chr1:1-100"),
      motifs = c("M1", "M1")
    ),
    annot_path
  )
  manifest <- tibble::tibble(motif = "M1", n_peaks = 2L, score = score_path, bound = bound_path, annot = annot_path)

  aligned <- align_footprints(
    manifest,
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 1,
    threads = 1L,
    cache_dir = cache_dir,
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = TRUE,
    output_mode = "distinct",
    return_id_map = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(aligned$id_map), 0L)
  expect_false(file.exists(file.path(cache_dir, "fp_id_map_TEST.csv")))
  expect_true(file.exists(file.path(cache_dir, "fp_scores_TEST.csv")))
  expect_equal(nrow(aligned$fp_annotation), 1L)
  expect_equal(aligned$fp_annotation$fp_peak, "chr1:10-22")
})

test_that("footprint alignment can write cache without returning large tables", {
  bench_dir <- file.path(tempdir(), paste0("craftgrn-align-paths-", as.integer(stats::runif(1L, 1, 1e9))))
  raw_dir <- file.path(bench_dir, "raw")
  cache_dir <- file.path(bench_dir, "cache")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  score_path <- file.path(raw_dir, "M1_score.csv")
  bound_path <- file.path(raw_dir, "M1_bound.csv")
  annot_path <- file.path(raw_dir, "M1_annotation.csv")
  readr::write_csv(tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1.1, 1.1)), score_path)
  readr::write_csv(tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1L, 1L)), bound_path)
  readr::write_csv(
    tibble::tibble(
      fp_peak = c("chr1:10-20", "chr1:12-22"),
      atac_peak = c("chr1:1-100", "chr1:1-100"),
      motifs = c("M1", "M1")
    ),
    annot_path
  )
  manifest <- tibble::tibble(motif = "M1", n_peaks = 2L, score = score_path, bound = bound_path, annot = annot_path)

  aligned <- align_footprints(
    manifest,
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 1,
    threads = 1L,
    cache_dir = cache_dir,
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = TRUE,
    output_mode = "distinct",
    return_data = FALSE,
    return_id_map = FALSE,
    write_id_map = FALSE,
    write_fp_sites = TRUE,
    verbose = FALSE
  )

  expect_equal(aligned$output_mode, "distinct")
  expect_equal(aligned$counts$fp_score, 1L)
  expect_equal(aligned$counts$fp_bound, 1L)
  expect_equal(aligned$counts$fp_annotation, 0L)
  expect_true(file.exists(aligned$paths$fp_score))
  expect_true(file.exists(aligned$paths$fp_bound))
  expect_false("fp_annotation" %in% names(aligned$paths))
  expect_true(file.exists(aligned$paths$fp_sites))
  expect_equal(aligned$counts$fp_sites, 1L)
  fp_sites <- data.table::fread(aligned$paths$fp_sites, showProgress = FALSE)
  expect_equal(fp_sites$peak_ID, "chr1:10-22")
  expect_equal(fp_sites$n_source_fp_peaks, 2L)
  expect_equal(fp_sites$source_fp_peaks, "chr1:10-20;chr1:12-22")
  expect_equal(fp_sites$motifs_all, "M1")
  expect_false("fp_score" %in% names(aligned))
  expect_false(file.exists(file.path(cache_dir, "fp_id_map_TEST.csv")))
})

test_that("footprint alignment cache summary handles compact cache without annotation file", {
  bench_dir <- file.path(tempdir(), paste0("craftgrn-align-cache-hit-", as.integer(stats::runif(1L, 1, 1e9))))
  raw_dir <- file.path(bench_dir, "raw")
  cache_dir <- file.path(bench_dir, "cache")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  score_path <- file.path(raw_dir, "M1_score.csv")
  bound_path <- file.path(raw_dir, "M1_bound.csv")
  annot_path <- file.path(raw_dir, "M1_annotation.csv")
  readr::write_csv(tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1.1, 1.1)), score_path)
  readr::write_csv(tibble::tibble(peak_ID = c("chr1:10-20", "chr1:12-22"), sample_a = c(1L, 1L)), bound_path)
  readr::write_csv(tibble::tibble(fp_peak = c("chr1:10-20", "chr1:12-22"), atac_peak = "chr1:1-100", motifs = "M1"), annot_path)
  manifest <- tibble::tibble(motif = "M1", n_peaks = 2L, score = score_path, bound = bound_path, annot = annot_path)

  invisible(align_footprints(
    manifest,
    mid_slop = 10L,
    score_match_pct = 1,
    threads = 1L,
    cache_dir = cache_dir,
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = TRUE,
    return_data = FALSE,
    return_id_map = FALSE,
    write_id_map = FALSE,
    write_fp_sites = TRUE,
    verbose = FALSE
  ))

  summary <- align_footprints(
    manifest,
    cache_dir = cache_dir,
    cache_tag = "TEST",
    use_cache = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_equal(summary$counts$fp_score, 1L)
  expect_equal(summary$counts$fp_bound, 1L)
  expect_equal(summary$counts$fp_annotation, 0L)
  expect_equal(summary$counts$fp_sites, 1L)
})

test_that("chromosome-parallel footprint alignment matches single-process alignment", {
  bench_dir <- file.path(tempdir(), paste0("craftgrn-align-parallel-", as.integer(stats::runif(1L, 1, 1e9))))
  raw_dir <- file.path(bench_dir, "raw")
  cache_dir <- file.path(bench_dir, "cache")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  score_path <- file.path(raw_dir, "M1_score.csv")
  bound_path <- file.path(raw_dir, "M1_bound.csv")
  annot_path <- file.path(raw_dir, "M1_annotation.csv")
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("chr1:10-20", "chr1:12-22", "chr2:40-50", "chr2:42-52"),
      sample_a = c(1.1, 1.1, 3.2, 3.2),
      sample_b = c(2.0, 2.0, 4.0, 4.0)
    ),
    score_path
  )
  readr::write_csv(
    tibble::tibble(
      peak_ID = c("chr1:10-20", "chr1:12-22", "chr2:40-50", "chr2:42-52"),
      sample_a = c(1L, 1L, 1L, 1L),
      sample_b = c(1L, 1L, 0L, 0L)
    ),
    bound_path
  )
  readr::write_csv(
    tibble::tibble(
      fp_peak = c("chr1:10-20", "chr1:12-22", "chr2:40-50", "chr2:42-52"),
      atac_peak = c("chr1:1-100", "chr1:1-100", "chr2:1-100", "chr2:1-100"),
      motifs = c("M1", "M1", "M1", "M1")
    ),
    annot_path
  )
  manifest <- tibble::tibble(motif = "M1", n_peaks = 4L, score = score_path, bound = bound_path, annot = annot_path)

  single <- align_footprints(
    manifest,
    mid_slop = 10L,
    score_match_pct = 1,
    threads = 1L,
    cache_dir = file.path(cache_dir, "single"),
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = FALSE,
    return_id_map = FALSE,
    write_id_map = FALSE,
    parallel_by = "none",
    verbose = FALSE
  )
  parallel <- align_footprints(
    manifest,
    mid_slop = 10L,
    score_match_pct = 1,
    threads = 2L,
    cache_dir = file.path(cache_dir, "parallel"),
    cache_tag = "TEST",
    use_cache = FALSE,
    write_cache = FALSE,
    return_id_map = FALSE,
    write_id_map = FALSE,
    parallel_by = "chromosome",
    verbose = FALSE
  )

  expect_equal(parallel$fp_score, single$fp_score)
  expect_equal(parallel$fp_bound, single$fp_bound)
  expect_equal(parallel$fp_annotation, single$fp_annotation)
})

test_that("footprint alignment component helper preserves duplicate signature components", {
  m <- matrix(
    c(
      1L, 2L, 3L,
      1L, 2L, 3L,
      1L, 2L, 4L,
      9L, 9L, 9L
    ),
    ncol = 3L,
    byrow = TRUE
  )

  ids <- craftgrn:::`.assign_fp_score_components`(m, k_req = 2L, compress_duplicates = TRUE)

  expect_equal(ids[1L], ids[2L])
  expect_equal(ids[1L], ids[3L])
  expect_false(ids[1L] == ids[4L])
})
