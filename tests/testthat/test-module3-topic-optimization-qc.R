test_that("Hellinger topic similarity averages Gene and Peak modalities", {
  phi <- rbind(
    Topic1 = c(0.4, 0.1, 0.4, 0.1),
    Topic2 = c(0.1, 0.4, 0.1, 0.4)
  )
  colnames(phi) <- c("GENE:A", "GENE:B", "PEAK:A", "PEAK:B")

  observed <- .m3_opt_hellinger_similarity(phi, 1:2, 3:4)
  expected <- sum(sqrt(c(0.8, 0.2) * c(0.2, 0.8)))

  expect_equal(unname(diag(observed)), c(1, 1))
  expect_equal(observed[1, 2], expected)
})

test_that("small topics merge into the larger deterministic representative", {
  phi <- rbind(
    Topic1 = c(0.45, 0.05, 0.45, 0.05),
    Topic2 = c(0.40, 0.10, 0.40, 0.10),
    Topic3 = c(0.05, 0.45, 0.05, 0.45)
  )
  colnames(phi) <- c("GENE:A", "GENE:B", "PEAK:A", "PEAK:B")
  theta <- rbind(
    `A::TF1` = c(0.8, 0.1, 0.1),
    `A::TF2` = c(0.1, 0.8, 0.1)
  )
  colnames(theta) <- rownames(phi)
  dtm <- Matrix::Matrix(
    1,
    nrow = 2,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(rownames(theta), colnames(phi))
  )

  result <- .m3_opt_merge_map(
    phi = phi,
    theta = theta,
    dtm = dtm,
    raw_topic_ids = 1:3,
    raw_links = c(1, 10, 10),
    raw_genes = c(1, 10, 10),
    gene_ids = 1:2,
    peak_ids = 3:4,
    min_genes = 2,
    min_links = 2,
    similarity_threshold = 0.98
  )

  expect_equal(result$mapping[[1]], 2)
  expect_equal(result$mapping[[2]], 2)
  expect_equal(result$audit$source_topic[[1]], 1)
  expect_equal(result$audit$representative_topic[[1]], 2)
})

test_that("small distinct topics are not forced into unrelated topics", {
  phi <- rbind(
    Topic1 = c(0.49, 0.01, 0.49, 0.01),
    Topic2 = c(0.01, 0.49, 0.01, 0.49),
    Topic3 = c(0.25, 0.25, 0.25, 0.25)
  )
  colnames(phi) <- c("GENE:A", "GENE:B", "PEAK:A", "PEAK:B")
  theta <- rbind(
    `A::TF1` = c(0.8, 0.1, 0.1),
    `A::TF2` = c(0.1, 0.8, 0.1)
  )
  colnames(theta) <- rownames(phi)
  dtm <- Matrix::Matrix(
    1,
    nrow = 2,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(rownames(theta), colnames(phi))
  )

  result <- .m3_opt_merge_map(
    phi = phi,
    theta = theta,
    dtm = dtm,
    raw_topic_ids = 1:3,
    raw_links = c(1, 10, 10),
    raw_genes = c(1, 10, 10),
    gene_ids = 1:2,
    peak_ids = 3:4,
    min_genes = 2,
    min_links = 2,
    similarity_threshold = 0.90
  )

  expect_equal(result$mapping, 1:3)
  expect_equal(nrow(result$audit), 0L)
})

test_that("topic merging can recover Gene and Peak agreement", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  colnames(phi) <- c("GENE:A", "PEAK:A")
  pair <- data.table::data.table(
    target_gene = "A",
    gene_term_id = "GENE:A",
    peak_term_id = "PEAK:A",
    gene_gammafit_topics = "1",
    peak_gammafit_topics = "2",
    assigned_topic = NA_character_
  )

  result <- .m3_opt_target_assignment(
    phi = phi,
    pair_assignment = pair,
    raw_topic_ids = 1:2,
    raw_to_group = c(2L, 2L)
  )

  expect_equal(result$optimized_gene_topic, 2)
  expect_equal(result$optimized_peak_topic, 2)
  expect_equal(result$optimized_assigned_topic, 2)
  expect_true(result$recovered_after_merge)
  expect_true("raw_assignment_status" %in% names(result))
})

test_that("gene-only optimization records its assignment mode and settings", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  colnames(phi) <- c("GENE:A", "GENE:B")
  theta <- rbind(
    `C1::TF1` = c(0.8, 0.2),
    `C2::TF1` = c(0.2, 0.8)
  )
  colnames(theta) <- rownames(phi)
  dtm <- Matrix::Matrix(
    matrix(c(5, 1, 1, 5), nrow = 2, byrow = TRUE),
    sparse = TRUE,
    dimnames = list(rownames(theta), colnames(phi))
  )
  assignments <- data.table::data.table(
    target_gene = c("A", "B"),
    assigned = TRUE,
    assigned_topic = c(1L, 2L),
    assignment_status = "assigned_gammafit_maxprob",
    gene_term_id = c("GENE:A", "GENE:B"),
    peak_term_id = c("GENE:A", "GENE:B"),
    gene_gammafit_topics = c("1", "2"),
    peak_gammafit_topics = c("1", "2")
  )
  topic_terms <- data.table::CJ(topic_num = 1:2, term_id = colnames(phi))
  topic_terms[, `:=`(
    topic = topic_num,
    term_group = "GENE",
    score = phi[cbind(topic_num, match(term_id, colnames(phi)))],
    in_topic = (topic_num == 1L & term_id == "GENE:A") |
      (topic_num == 2L & term_id == "GENE:B")
  )]
  topic_terms[, gammafit_candidate := in_topic]

  observed <- .module3_optimize_condition_topics(
    theta = theta,
    phi = phi,
    dtm = dtm,
    topic_terms = topic_terms,
    pair_assignment = assignments,
    assignment_mode = "gene_only",
    min_genes = 1L,
    min_links = 1L,
    similarity_threshold = 1,
    tf_topic_cutoff = 0,
    umap_max_links_per_condition = 10L,
    seed = 7L
  )

  expect_identical(observed$assignment_mode, "gene_only")
  expect_identical(observed$qc_seed, 7L)
  expect_equal(observed$pair_assignment$optimized_gene_topic,
               observed$pair_assignment$optimized_peak_topic)
})

test_that("optimized assignment summaries retain pre-merge failure reasons", {
  assignments <- data.table::data.table(
    gene_term_id = c("GENE:A", "GENE:B", "GENE:C"),
    peak_term_id = c("PEAK:A", "PEAK:B", "PEAK:C"),
    assigned = c(TRUE, TRUE, FALSE),
    assignment_status = c(
      "assigned_after_topic_merge",
      "assigned_gammafit_maxprob_agreement",
      "unassigned_after_topic_merge"
    ),
    raw_assignment_status = c(
      "maxprob_topic_disagreement",
      "assigned_gammafit_maxprob_agreement",
      "peak_no_gammafit_candidate"
    ),
    recovered_after_merge = c(TRUE, FALSE, FALSE),
    raw_assigned_topic = c(NA, 2L, 3L),
    common_candidate_count = c(2L, 1L, 0L),
    min_maxprob_probability = c(0.6, 0.7, NA),
    min_maxprob_margin = c(0.2, 0.3, NA)
  )

  summary <- .topic_gene_peak_assignment_summary(
    assignments,
    thrP = 0.7,
    assignment_method = "gammafit_maxprob_optimized",
    model_family = "lda"
  )

  expect_equal(summary$maxprob_topic_disagreement_count, 1L)
  expect_equal(summary$no_peak_gammafit_candidate_count, 1L)
  expect_equal(summary$recovered_after_topic_merge_count, 1L)
  expect_equal(summary$lost_after_topic_merge_count, 1L)
})

test_that("document-supported expression excludes genes absent from condition documents", {
  optimization <- list(
    qc = list(
      target_levels = c("A", "B"),
      assignments = data.table::data.table(
        doc_index = c(1L, 1L, 2L),
        target_index = c(1L, 2L, 1L),
        condition_id = c("C1", "C1", "C2"),
        optimized_topic = c(1L, 1L, 1L)
      )
    ),
    condition_gene_expression = data.table::data.table(
      condition_id = rep(c("C1", "C2"), each = 2L),
      target_gene = rep(c("A", "B"), 2L),
      expression = c(3, 15, 7, 1023)
    )
  )

  observed <- .m3_qc_document_supported_topic_expression_matrix(
    optimization,
    conditions = c("C1", "C2"),
    topics = 1L
  )

  expect_equal(observed["C1", "Topic 1"], mean(log2(c(3, 15) + 1)))
  expect_equal(observed["C2", "Topic 1"], log2(7 + 1))
  expect_false(observed["C2", "Topic 1"] == mean(log2(c(7, 1023) + 1)))
})

test_that("document expression topic shares use analyzed targets and average documents", {
  optimization <- list(
    qc = list(
      target_levels = c("A", "B"),
      assignments = data.table::data.table(
        doc_index = c(1L, 1L, 2L),
        target_index = c(1L, 2L, 1L),
        condition_id = c("C1", "C1", "C1"),
        optimized_topic = c(1L, 2L, 1L)
      )
    ),
    condition_gene_expression = data.table::data.table(
      condition_id = "C1",
      target_gene = c("A", "B", "UNUSED"),
      expression = c(3, 15, 1023)
    )
  )

  observed <- .m3_qc_document_expression_topic_share_matrix(
    optimization,
    conditions = "C1",
    topics = 1:2
  )
  first_document_topic_1 <- log2(3 + 1) / sum(log2(c(3, 15) + 1))
  first_document_topic_2 <- log2(15 + 1) / sum(log2(c(3, 15) + 1))

  expect_equal(observed["C1", "Topic 1"], (first_document_topic_1 + 1) / 2)
  expect_equal(observed["C1", "Topic 2"], first_document_topic_2 / 2)
  expect_equal(sum(observed["C1", ]), 1)
})

test_that("topic optimization controls validate in nested project config", {
  expect_silent(.validate_project_config(list(
    module3 = list(
      optimize_topics = TRUE,
      merge_min_genes = 50L,
      merge_min_links = 200L,
      merge_similarity_threshold = 0.9,
      run_topic_assignment_qc = TRUE,
      qc_umap_links_per_condition = 1000L,
      qc_top_tfs = 150L,
      qc_reference_condition = "ConditionB",
      qc_upregulated_log2fc_min = 1,
      qc_upregulated_pseudocount = 1,
      qc_seed = 1L
    )
  )))
  expect_error(
    .validate_project_config(list(
      module3 = list(merge_similarity_threshold = 1.1)
    )),
    "at most 1"
  )
  expect_error(
    .validate_project_config(list(
      module3 = list(qc_upregulated_log2fc_min = -1)
    )),
    "non-negative"
  )
})

test_that("topic run config resolves reproducible WarpLDA settings", {
  cfg <- list(
    topic_method = "condition_aggr_lda",
    topic_k = 30L,
    warplda_iterations = 20L,
    topic_warplda_sampler = "warp_omp",
    topic_warplda_beta = NULL,
    topic_warplda_seed = 123L,
    topic_count_method = "log",
    topic_count_scale = 50,
    topic_count_input = "pseudo_count_log"
  )

  expect_silent(.validate_project_config(cfg))
  observed <- .module3_resolve_topic_run_config(project_config = cfg)

  expect_identical(observed$warplda_iterations, 20L)
  expect_identical(observed$warplda_sampler, "warp_omp")
  expect_null(observed$warplda_beta)
  expect_identical(observed$warplda_seed, 123L)
  expect_identical(observed$count_method, "log")
  expect_equal(observed$count_scale, 50)
  expect_identical(observed$count_input, "pseudo_count_log")
})

test_that("repeated condition IDs retain concise stable display labels", {
  observed <- .m3_qc_short_condition_labels(c(
    "HPAFII_0_FBS_Ctrl",
    "HPAFII_0_FBS_Ctrl",
    "HPAFII_10_FBS_Ctrl"
  ))

  expect_equal(observed, c("0 FBS", "0 FBS", "10 FBS"))

  tcell <- c(
    "BATF_FloxedOut",
    "IRF4_FloxedOut",
    "WT_FloxedOut",
    "BATF_FloxReady",
    "IRF4_FloxReady",
    "WT_FloxReady"
  )
  expect_equal(
    .m3_qc_short_condition_labels(tcell),
    gsub("_", " ", tcell, fixed = TRUE)
  )
})

test_that("topic assignment retention labels distinguish links and genes", {
  expect_identical(
    .m3_qc_retention_labels(),
    list(
      links = c(
        "Total unique input links",
        "Links pass gamma filter",
        "Links pass max filter",
        "Links pass TF filter"
      ),
      genes = c(
        "Total unique input genes",
        "Genes pass gamma filter",
        "Genes pass max filter",
        "Genes pass TF filter"
      )
    )
  )
})

test_that("topic retention counts links and genes uniquely across conditions", {
  links <- data.table::data.table(
    condition_id = c("A", "B", "A", "B", "A"),
    tf_index = c(1L, 1L, 2L, 2L, 1L),
    target_index = c(1L, 1L, 1L, 1L, 2L)
  )

  expect_identical(.m3_qc_unique_link_count(links), 3L)
  expect_identical(
    .m3_qc_unique_target_count(c("G1", "G1", "G2", NA_character_, "")),
    2L
  )
})

test_that("TF heatmap pools conditions and deduplicates targets per topic", {
  n_tfs <- 160L
  assignments <- data.table::data.table(
    condition_id = rep(c("Condition_A", "Condition_B"), each = n_tfs * 2L),
    tf_index = rep(rep(seq_len(n_tfs), each = 2L), 2L),
    target_index = rep(seq_len(n_tfs * 2L), 2L),
    optimized_topic = rep(rep(1:2, n_tfs), 2L),
    optimized_aligned = TRUE
  )
  optimization <- list(qc = list(
    assignments = assignments,
    optimized_topic_ids = 1:2,
    tf_levels = paste0("TF", seq_len(n_tfs))
  ))

  matrix <- .m3_qc_pooled_tf_topic_matrix(
    optimization,
    top_n_tfs = 150
  )
  pages <- .m3_qc_tf_topic_pages(
    optimization,
    top_n_tfs = 150
  )

  expect_equal(dim(matrix), c(75, 2))
  expect_equal(attr(matrix, "selection_per_topic"), 75L)
  expect_equal(unname(matrix["TF1", ]), c(1, 1))
  expect_length(pages, 1)
  expect_s3_class(pages[[1]], "gtable")
})

test_that("QC counts unique TF-target pairs and TFs", {
  docs <- data.table::data.table(tf = c("TF1", "TF2"))
  pairs <- data.table::data.table(target_gene = c("G1", "G2"))
  assignments <- data.table::data.table(
    doc_index = c(1L, 1L, 1L, 2L),
    target_index = c(1L, 1L, 2L, 1L),
    condition_id = c("C1", "C1", "C1", "C1"),
    raw_aligned = TRUE,
    optimized_aligned = TRUE,
    raw_target_topic = 1L,
    optimized_topic = 1L
  )

  observed <- .m3_opt_qc_tables(
    assignments = assignments,
    docs = docs,
    pairs = pairs,
    raw_similarity = matrix(1, 1, 1),
    optimized_similarity = matrix(1, 1, 1),
    raw_topic_ids = 1L,
    optimized_topic_ids = 1L
  )

  expect_equal(observed$raw_counts$links, 3L)
  expect_equal(observed$raw_counts$genes, 2L)
  expect_equal(observed$raw_counts$tfs, 2L)
  expect_equal(observed$condition_topic$links, 3L)
  expect_equal(observed$condition_topic$tfs, 2L)
})

test_that("condition topics rank by mean theta before Jaccard", {
  theta <- rbind(
    `A::TF1` = c(0.60, 0.15, 0.05, 0.20),
    `A::TF2` = c(0.40, 0.20, 0.10, 0.30),
    `B::TF1` = c(0.05, 0.10, 0.55, 0.30),
    `B::TF2` = c(0.10, 0.20, 0.40, 0.30)
  )
  colnames(theta) <- paste0("Topic", 1:4)
  assignments <- data.table::data.table(
    doc_index = rep(seq_len(4L), each = 2L),
    condition_id = rep(c("A", "A", "B", "B"), each = 2L)
  )
  optimization <- list(
    theta = theta,
    qc = list(assignments = assignments)
  )
  jaccard <- rbind(
    A = c(0.10, 0.95, 0.20),
    B = c(0.90, 0.20, 0.10)
  )
  colnames(jaccard) <- as.character(1:3)

  condition_theta <- .m3_qc_condition_topic_theta(optimization)
  top <- .m3_qc_top_condition_topics(
    condition_theta,
    jaccard,
    top_n = 3L
  )

  expect_equal(unname(condition_theta["A", ]), c(0.50, 0.175, 0.075, 0.25))
  expect_equal(unname(condition_theta["B", ]), c(0.075, 0.15, 0.475, 0.30))
  expect_equal(top[condition_id == "A", topic_num], c(1L, 2L, 3L))
  expect_equal(top[condition_id == "B", topic_num], c(3L, 2L, 1L))
  expect_false(4L %in% top$topic_num)
  expect_equal(top[condition_id == "A" & topic_num == 2L, mean_jaccard], 0.95)
})

test_that("condition-topic expression uses one fixed expressed-gene pool", {
  edges <- data.table::data.table(
    condition_id = c("A", "A", "A", "B", "B", "A"),
    gene_key = c("G1", "G1", "G2", "G1", "G2", "G3"),
    gene_expr_condition = c(3, 3, 1, 0, 7, 15)
  )
  expression <- .m3_qc_condition_gene_expression(
    edges,
    genes = c("G1", "G2", "G3", "G4")
  )
  pairs <- data.table::data.table(
    target_gene = c("G1", "G2", "G3", "G4"),
    optimized_assigned_topic = 1:4
  )
  pairs[2L, optimized_assigned_topic := 1L]
  matrix <- .m3_qc_condition_topic_expression_matrix(
    expression,
    pairs,
    conditions = c("A", "B"),
    topics = 1:4
  )

  expect_equal(nrow(expression), 5L)
  expect_equal(unname(matrix["A", "Topic 1"]), 1.5)
  expect_equal(unname(matrix["B", "Topic 1"]), 1.5)
  expect_equal(unname(matrix["A", "Topic 3"]), 4)
  expect_equal(unname(matrix["B", "Topic 3"]), 0)
  expect_equal(unname(matrix[, "Topic 4"]), c(0, 0))
})

test_that("condition expression QC remains optional for legacy edges", {
  observed <- .m3_qc_condition_gene_expression(
    data.table::data.table(condition_id = "A", gene_key = "G1"),
    genes = "G1"
  )

  expect_named(observed, c("condition_id", "target_gene", "expression"))
  expect_equal(nrow(observed), 0L)
})

test_that("condition upregulation is relative to one explicit reference", {
  expression <- data.table::data.table(
    condition_id = rep(c("A", "B"), each = 3L),
    target_gene = rep(c("G1", "G2", "G3"), 2L),
    expression = c(7, 3, 1, 1, 2, 5)
  )
  upregulated <- .m3_qc_upregulated_genes_vs_reference(
    expression,
    reference_condition = "B",
    log2fc_min = 1,
    pseudocount = 1
  )
  pairs <- data.table::data.table(
    target_gene = c("G1", "G2", "G3"),
    optimized_assigned_topic = c(1L, 1L, 2L)
  )
  matrix <- .m3_qc_condition_topic_upregulated_matrix(
    upregulated,
    pairs,
    conditions = c("A", "B"),
    topics = 1:2
  )

  expect_equal(upregulated, data.table::data.table(
    condition_id = "A",
    target_gene = "G1"
  ))
  expect_equal(unname(matrix), matrix(c(1, 0, 0, 0), nrow = 2))
})

test_that("signed comparison audit identifies upregulated condition genes", {
  audit <- data.table::data.table(
    condition_1 = c("A", "C", "A"),
    condition_2 = c("B", "B", "B"),
    gene_key = c("G1", "G2", "G3"),
    expression_1 = c(7, 9, 3),
    expression_2 = c(1, 1, 2)
  )

  observed <- .m3_qc_upregulated_genes_from_comparisons(
    audit,
    reference_condition = "B",
    log2fc_min = 1
  )

  expect_equal(observed, data.table::data.table(
    condition_id = c("A", "C"),
    target_gene = c("G1", "G2")
  ))
})

test_that("comparison audit supplies deduplicated condition expression", {
  audit <- data.table::data.table(
    condition_1 = c("A", "A"),
    condition_2 = c("B", "C"),
    gene_key = c("G1", "G1"),
    expression_1 = c(4, 4),
    expression_2 = c(2, 7)
  )

  observed <- .m3_qc_condition_expression_from_comparisons(audit)

  expect_equal(observed, data.table::data.table(
    condition_id = c("A", "B", "C"),
    target_gene = rep("G1", 3),
    expression = c(4, 2, 7)
  ))
})

test_that("canonical long table supplies condition expression", {
  canonical <- data.table::data.table(
    condition_id = c("A", "A", "B", "B"),
    gene_key = c("G1", "G2", "G1", "G2"),
    expression = c(10, 0, 4, 8),
    expressed = c(TRUE, FALSE, TRUE, TRUE)
  )

  observed <- .m3_qc_condition_expression_from_comparisons(
    canonical,
    genes = "G1"
  )

  expect_equal(observed, data.table::data.table(
    condition_id = c("A", "B"),
    target_gene = c("G1", "G1"),
    expression = c(10, 4)
  ))
})

test_that("QC panel alignment reserves one shared legend width", {
  values <- data.frame(x = 1:2, y = 1, value = 1:2)
  short <- ggplot2::ggplot(values, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_continuous(name = "Count")
  long <- ggplot2::ggplot(values, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_continuous(name = "Mean log2(expression + 1)")

  aligned <- .m3_qc_align_plot_dimensions(list(short, long))

  expect_equal(aligned[[1L]]$widths, aligned[[2L]]$widths)
  expect_equal(aligned[[1L]]$heights, aligned[[2L]]$heights)
})

test_that("QC heatmap styling adapts without dropping below nine points", {
  small <- .m3_qc_heatmap_layout_spec(5, 10)
  large <- .m3_qc_heatmap_layout_spec(20, 40)

  expect_gte(small$axis_size, large$axis_size)
  expect_gte(large$axis_size, 9)
  expect_equal(.m3_qc_heatmap_fill_max(c(0, 1:100)), 95.65, tolerance = 0.01)
})

test_that("condition-topic layouts preserve pairs for odd counts", {
  layout <- .m3_qc_pair_layout(3L)

  expect_equal(layout$layout_matrix[1, ], rep(1L, 4L))
  expect_equal(layout$layout_matrix[2, ], 2:5)
  expect_equal(layout$layout_matrix[3, 1:2], 6:7)
  expect_true(all(is.na(layout$layout_matrix[3, 3:4])))
  expect_equal(layout$heights, c(1.35, 1, 1))
})

test_that("topic structure and count heatmaps use square cells", {
  similarity <- diag(3)
  rownames(similarity) <- colnames(similarity) <- paste0("Topic ", 1:3)
  counts <- data.table::data.table(
    raw_topic = 1:3,
    links = c(100, 50, 25),
    genes = c(20, 10, 5),
    tfs = c(10, 7, 4)
  )

  structure_plot <- .m3_qc_topic_structure_plot(
    similarity = similarity,
    counts = counts,
    topic_column = "raw_topic",
    topic_order = 1:3,
    title = "Structure",
    subtitle = "Square cells"
  )
  count_plot <- .m3_qc_count_heatmap(
    similarity,
    title = "Counts"
  )

  expect_equal(structure_plot$coordinates$ratio, 1)
  expect_equal(count_plot$coordinates$ratio, 1)
})

test_that("UMAP display jitter is small and deterministic", {
  coordinates <- data.table::data.table(
    UMAP1 = seq(0, 10, length.out = 20),
    UMAP2 = seq(-5, 5, length.out = 20)
  )

  first <- .m3_qc_jitter_umap(coordinates, seed = 17L)
  second <- .m3_qc_jitter_umap(coordinates, seed = 17L)

  expect_equal(first, second)
  expect_true(max(abs(first$UMAP1 - coordinates$UMAP1)) <= 0.05)
  expect_true(max(abs(first$UMAP2 - coordinates$UMAP2)) <= 0.05)
})

test_that("condition correlation combines normalized link and gene profiles", {
  links <- rbind(
    A = c(10, 0, 0),
    B = c(20, 0, 0),
    C = c(0, 0, 10)
  )
  genes <- rbind(
    A = c(5, 0, 0),
    B = c(10, 0, 0),
    C = c(0, 0, 5)
  )

  observed <- .m3_qc_condition_correlation(links, genes)
  plot <- .m3_qc_correlation_plot(observed, "Condition correlation")

  expect_equal(observed["A", "B"], 1)
  expect_lt(observed["A", "C"], 0)
  expect_equal(unname(diag(observed)), rep(1, 3))
  expect_equal(observed, t(observed))
  expect_equal(plot$coordinates$ratio, 1)
})

test_that("dense heatmap counts use compact integer labels", {
  expect_equal(
    .m3_qc_heatmap_count(c(999, 1157, 29640, 1500000)),
    c("999", "1k", "30k", "2M")
  )
})
