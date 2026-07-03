test_that("gene symbol resolver falls back to deterministic case keys", {
  res <- resolve_gene_symbols(
    c("Atp2b4", "ATP2B4", "missing gene", NA),
    species = "mouse",
    use_annotation = FALSE
  )

  expect_equal(res$input_gene, c("Atp2b4", "ATP2B4", "missing gene", NA))
  expect_equal(res$match_key[1:2], c("ATP2B4", "ATP2B4"))
  expect_equal(res$canonical_symbol[1:2], c("Atp2b4", "Atp2b4"))
  expect_equal(res$match_type[1:2], c("heuristic_case", "heuristic_case"))
  expect_false(res$matched[3])
  expect_true(is.na(res$match_key[4]))
})

test_that("human gene symbol resolver uses formal official symbols and aliases", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  res <- resolve_gene_symbols(
    c("STAT1", "stat1", "P53", "ENSG00000141510"),
    species = "human"
  )
  by_input <- stats::setNames(seq_len(nrow(res)), res$input_gene)

  expect_equal(res$canonical_symbol[by_input[["STAT1"]]], "STAT1")
  expect_equal(res$match_type[by_input[["STAT1"]]], "official_exact")
  expect_equal(res$canonical_symbol[by_input[["stat1"]]], "STAT1")
  expect_true(res$match_type[by_input[["stat1"]]] %in% c("official_case", "alias"))
  expect_equal(res$canonical_symbol[by_input[["P53"]]], "TP53")
  expect_true(res$matched[by_input[["P53"]]])
  expect_equal(res$canonical_symbol[by_input[["ENSG00000141510"]]], "TP53")
})

test_that("mouse gene symbol resolver handles dataset-style titlecase genes", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Mm.eg.db")

  res <- resolve_gene_symbols(
    c("Atp2b4", "ATP2B4", "Selplg", "SELPLG"),
    species = "mouse"
  )

  expect_equal(res$canonical_symbol, c("Atp2b4", "Atp2b4", "Selplg", "Selplg"))
  expect_true(all(res$matched))
  expect_true(all(res$match_type %in% c("official_exact", "official_case", "alias")))
})

test_that("gene symbol keys preserve original labels and canonical matching keys", {
  resolved <- .gene_symbol_key_table(
    c("Atp2b4", "ATP2B4", "Selplg"),
    species = "mouse",
    use_annotation = FALSE
  )

  expect_equal(resolved$gene, c("Atp2b4", "ATP2B4", "Selplg"))
  expect_equal(resolved$gene_key__, c("ATP2B4", "ATP2B4", "SELPLG"))
  expect_equal(resolved$gene_canonical, c("Atp2b4", "Atp2b4", "Selplg"))
})

test_that("pathway query canonicalizer uses formal symbols before enrichment", {
  testthat::skip_if_not_installed("AnnotationDbi")
  testthat::skip_if_not_installed("org.Hs.eg.db")

  out <- .canonicalize_pathway_genes(c("P53", "stat1", "missing gene"), pathway_species = "human")

  expect_true("TP53" %in% out)
  expect_true("STAT1" %in% out)
  expect_false("P53" %in% out)
  expect_false("missing gene" %in% out)
})

test_that("gene symbol conversion audit summarizes formal and fallback matches", {
  out_dir <- withr::local_tempdir()
  x <- data.table::data.table(
    gene = c("Stat1", "STAT1", "bad gene"),
    gene_canonical = c("Stat1", "Stat1", NA_character_),
    gene_match_type = c("official_exact", "heuristic_case", "unmatched"),
    gene_matched = c(TRUE, TRUE, FALSE),
    gene_ambiguous = c(FALSE, FALSE, FALSE)
  )

  audit <- .write_gene_symbol_conversion_audit(out_dir, list(test_genes = x))

  expect_true(file.exists(file.path(out_dir, "gene_symbol_conversion_audit.csv")))
  expect_true(file.exists(file.path(out_dir, "gene_symbol_conversion_details.csv")))
  expect_equal(audit$source, "test_genes")
  expect_equal(audit$total, 3L)
  expect_equal(audit$matched, 2L)
  expect_equal(audit$formal_matched, 1L)
  expect_equal(audit$heuristic_matched, 1L)
  expect_equal(audit$unmatched, 1L)
})
