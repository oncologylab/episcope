test_that("condition documents allocate exact post-transform modality tokens", {
  values <- data.table::data.table(
    condition_id = rep(c("CondA", "CondB"), each = 5),
    term_id = rep(
      c("GENE:G1", "GENE:G2", "PEAK:chr1:1-10", "PEAK:chr1:20-30", "PEAK:chr1:40-50"),
      2
    ),
    modality = rep(c("gene", "gene", "peak", "peak", "peak"), 2),
    token_score = c(9, 1, 8, 1, 1, 2, 8, 1, 2, 7)
  )
  observed <- craftgrn:::.module3_allocate_condition_document_tokens(
    values,
    gene_fraction = 0.6,
    token_budget = 100L
  )
  audit <- observed[, .(
    tokens = sum(pseudo_count_log),
    terms = .N,
    minimum = min(pseudo_count_log)
  ), by = .(condition_id, modality)]

  expect_true(all(audit$minimum >= 1L))
  expect_true(all(audit[modality == "gene", tokens] == 60L))
  expect_true(all(audit[modality == "peak", tokens] == 40L))
  expect_true(all(
    observed[, .(tokens = sum(pseudo_count_log)), by = condition_id]$tokens ==
      100L
  ))
  expect_identical(
    observed$pseudo_count_log,
    craftgrn:::.module3_allocate_condition_document_tokens(
      values,
      gene_fraction = 0.6,
      token_budget = 100L
    )$pseudo_count_log
  )
})

test_that("module3_construct_docs supports one document per condition", {
  root <- withr::local_tempdir()
  input_dir <- file.path(root, "condition_links")
  output_dir <- file.path(root, "topic_documents")
  dir.create(input_dir, recursive = TRUE)
  rows <- list(
    CondA = data.table::data.table(
      condition_id = "CondA",
      gene_key = rep(c("G1", "G2"), each = 2),
      peak_id = rep(c("chr1:1-10", "chr1:20-30"), 2),
      gene_expr_condition = rep(c(20, 15), each = 2),
      fp_score_condition = c(8, 3, 8, 3)
    ),
    CondB = data.table::data.table(
      condition_id = "CondB",
      gene_key = rep(c("G1", "G2"), each = 2),
      peak_id = rep(c("chr1:1-10", "chr1:20-30"), 2),
      gene_expr_condition = rep(c(12, 25), each = 2),
      fp_score_condition = c(2, 9, 2, 9)
    )
  )
  paths <- vapply(names(rows), function(condition) {
    path <- file.path(input_dir, paste0(condition, ".csv"))
    data.table::fwrite(rows[[condition]], path)
    path
  }, character(1L))
  data.table::fwrite(
    data.table::data.table(
      condition_id = names(rows),
      path = paths,
      format = "csv",
      n_links = vapply(rows, nrow, integer(1L))
    ),
    file.path(input_dir, "condition_links_manifest.csv")
  )

  result <- module3_construct_docs(
    filtered_dir = input_dir,
    output_dir = output_dir,
    input_source = "condition_links",
    doc_mode = "tf",
    doc_design = "condition",
    condition_document_unit = "condition",
    condition_gene_token_fraction = 0.6,
    condition_document_token_budget = 100L,
    fp_term_mode = "unique",
    gene_term_mode = "unique",
    include_tf_terms = FALSE,
    count_method = "log",
    threshold_gene_expr = 10,
    threshold_fp_score = 0,
    overwrite = TRUE,
    verbose = FALSE
  )
  doc_term <- readRDS(file.path(output_dir, "rds", "doc_term.rds"))
  summary <- data.table::fread(
    file.path(output_dir, "topic_input_summary.csv")
  )

  expect_false(result$reused)
  expect_setequal(unique(doc_term$doc_id), c("CondA", "CondB"))
  expect_false(any(grepl("::", doc_term$doc_id, fixed = TRUE)))
  expect_equal(summary$condition_document_unit, "condition")
  expect_equal(summary$n_documents, 2)
  expect_true(all(
    doc_term[, .(tokens = sum(pseudo_count_log)), by = doc_id]$tokens == 100L
  ))
})

test_that("TF-topic evidence uses unique correlation-weighted peaks", {
  links <- data.table::data.table(
    condition_id = c("C1", "C1", "C1", "C1", "C2", "C2"),
    tf = "TF1",
    peak_id = c("P1", "P1", "P2", "P2", "P1", "P2"),
    gene_key = c("G1", "G2", "G2", "G3", "G1", "G3"),
    topic_num = c(1L, 1L, 2L, 2L, 1L, 2L),
    tf_expression = 20
  )
  peaks <- data.table::data.table(
    condition_id = rep(c("C1", "C2"), each = 2),
    peak_id = rep(c("P1", "P2"), 2),
    topic_num = rep(1:2, 2),
    token_share = c(0.8, 0.2, 0.7, 0.3)
  )
  genes <- data.table::data.table(
    condition_id = c("C1", "C1", "C1", "C2", "C2"),
    gene_key = c("G1", "G2", "G3", "G1", "G3"),
    topic_num = c(1L, 1L, 2L, 1L, 2L),
    token_share = c(0.5, 0.3, 0.2, 0.6, 0.4)
  )
  prediction <- data.table::data.table(
    tf = "TF1",
    fp_id = c("P1", "P2"),
    best_r = c(0.9, 0.4)
  )
  observed <- craftgrn:::.module3_condition_tf_topic_evidence(
    aligned_links = links,
    peak_values = peaks,
    gene_values = genes,
    prediction_stats = prediction,
    topic_count = 2L,
    tf_expression_min = 10
  )
  evidence <- observed$evidence

  expect_equal(observed$join_audit$matched_percent, 100)
  expect_equal(
    evidence[condition_id == "C1" & topic_num == 1L, n_unique_peaks],
    1
  )
  expect_gt(
    evidence[condition_id == "C1" & topic_num == 1L, combined_probability],
    evidence[condition_id == "C1" & topic_num == 2L, combined_probability]
  )
  expect_true(all(evidence$global_primary_topic == 1L))
  expect_true(all(evidence$primary_topic_num == 1L))
  expect_true(all(evidence$assignment_source == "binding_target_evidence"))
  expect_equal(
    evidence[, sum(combined_probability), by = .(condition_id, tf)]$V1,
    c(1, 1),
    tolerance = 1e-10
  )
})

test_that("condition-document report state removes TF-document theta semantics", {
  html <- craftgrn:::.m3cr_condition_report_js()
  expect_match(
    html,
    "CONDITION_DOCUMENT_MODE",
    fixed = TRUE
  )
  expect_match(
    html,
    "P(topic | condition)",
    fixed = TRUE
  )
  expect_match(
    html,
    "TF-topic evidence probability",
    fixed = TRUE
  )
  expect_match(
    html,
    "aggregationLabel.style.display=CONDITION_DOCUMENT_MODE?'none':''",
    fixed = TRUE
  )
  expect_match(
    html,
    "topicCutoff.value='0.1'",
    fixed = TRUE
  )
})
