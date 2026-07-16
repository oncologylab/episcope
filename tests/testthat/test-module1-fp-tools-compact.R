test_that("fp-tools compact caches import into aligned footprints", {
  root <- tempfile("fp-tools-")
  samples <- c("S1", "S2")
  for (sample_id in samples) {
    cache <- file.path(root, sample_id, "match_motifs", "cache")
    dir.create(cache, recursive = TRUE)
    sites <- data.table::data.table(
      motif = c("Batf_M1", "Batf_M1", "Irf4_M2"),
      TFBS_chr = "chr1",
      TFBS_start = c(10L, 20L, 20L),
      TFBS_end = c(15L, 25L, 25L),
      peak_chr = "chr1",
      peak_start = 1L,
      peak_end = 100L,
      score = if (sample_id == "S1") c(1, 3, 3) else c(2, 4, 4)
    )
    connection <- gzfile(file.path(cache, "motif_sites.tsv.gz"), open = "wt")
    utils::write.table(sites, connection, sep = "\t", row.names = FALSE, quote = FALSE)
    close(connection)
    writeLines(
      jsonlite::toJSON(list(thresholds = stats::setNames(list(2.5), sample_id)), auto_unbox = TRUE),
      file.path(cache, "thresholds.json")
    )
    data.table::fwrite(
      data.table::data.table(
        output_prefix = c("Batf_M1", "Irf4_M2"),
        name = c("Batf", "Irf4"),
        motif_id = c("M1", "M2")
      ),
      file.path(root, sample_id, "match_motifs", "motif_matches_results.txt"),
      sep = "\t"
    )
  }

  out <- .load_fp_tools_compact_footprints(root, samples, cache_dir = NULL, verbose = FALSE)
  expect_equal(nrow(out$fp_score), 2L)
  expect_equal(out$fp_score$S1, c(1, 3))
  expect_equal(out$fp_score$S2, c(2, 4))
  expect_equal(out$fp_bound$S1, c(FALSE, TRUE))
  expect_equal(sort(out$motif_db$gene_symbol), c("Batf", "Irf4"))
  expect_equal(nrow(out$fp_annotation), 3L)
  expect_equal(out$fp_sites$source_fp_peaks, out$fp_sites$peak_ID)
  expect_equal(out$fp_sites$n_source_fp_peaks, rep(1L, 2L))
  expect_equal(.module1_raw_footprint_count(out), 2)
})
