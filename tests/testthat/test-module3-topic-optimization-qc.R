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

test_that("coordinate Peak link universe counts repeated Peak support as one Gene", {
  dtm <- Matrix::Matrix(
    c(
      2, 3, 4,
      5, 6, 7
    ),
    nrow = 2L,
    byrow = TRUE,
    sparse = TRUE,
    dimnames = list(
      c("A::TF1", "B::TF1"),
      c("GENE:G1", "PEAK:P1", "PEAK:P2")
    )
  )
  pairs <- data.table::data.table(
    target_gene = c("G1", "G1"),
    gene_term_id = "GENE:G1",
    peak_term_id = c("PEAK:P1", "PEAK:P2")
  )

  observed <- .m3_opt_link_universe(dtm, pairs)

  expect_equal(observed$target_levels, "G1")
  expect_equal(unique(observed$links$target_index), 1L)
  expect_setequal(unique(observed$links$pair_index), 1:2)
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

test_that("topic merging prefers TF-term and theta correspondence", {
  phi <- rbind(
    Topic1 = c(0.45, 0.05, 0.45, 0.05),
    Topic2 = c(0.40, 0.10, 0.40, 0.10),
    Topic3 = c(0.30, 0.20, 0.30, 0.20)
  )
  colnames(phi) <- c(
    "GENE:TFX", "GENE:OTHER", "PEAK:TFX", "PEAK:OTHER"
  )
  theta <- rbind(`C1::TFX` = c(0.2, 0.1, 0.7))
  colnames(theta) <- rownames(phi)
  dtm <- Matrix::Matrix(
    1,
    nrow = 1,
    ncol = 4,
    sparse = TRUE,
    dimnames = list(rownames(theta), colnames(phi))
  )
  pair_assignment <- data.table::data.table(
    target_gene = "TFX",
    gene_term_id = "GENE:TFX",
    peak_term_id = "PEAK:TFX",
    gene_gammafit_topics = "1;2;3",
    peak_gammafit_topics = "1;2;3",
    assigned = TRUE,
    assigned_topic = 1L,
    assignment_status = "assigned_gammafit_maxprob_agreement"
  )

  legacy <- .m3_opt_merge_map(
    phi = phi,
    theta = theta,
    dtm = dtm,
    raw_topic_ids = 1:3,
    raw_links = c(1, 10, 10),
    raw_genes = c(1, 10, 10),
    gene_ids = 1:2,
    peak_ids = 3:4,
    pair_assignment = pair_assignment,
    min_genes = 2,
    min_links = 2,
    similarity_threshold = 0.9,
    prefer_tf_theta_correspondence = FALSE
  )
  preferred <- .m3_opt_merge_map(
    phi = phi,
    theta = theta,
    dtm = dtm,
    raw_topic_ids = 1:3,
    raw_links = c(1, 10, 10),
    raw_genes = c(1, 10, 10),
    gene_ids = 1:2,
    peak_ids = 3:4,
    pair_assignment = pair_assignment,
    min_genes = 2,
    min_links = 2,
    similarity_threshold = 0.9,
    prefer_tf_theta_correspondence = TRUE,
    tf_topic_cutoff = 0.3
  )

  expect_equal(legacy$mapping[[1L]], 2L)
  expect_equal(preferred$mapping[[1L]], 3L)
  expect_equal(preferred$audit$tf_theta_empty_before, 1L)
  expect_equal(preferred$audit$tf_theta_empty_after, 0L)
  expect_gt(
    preferred$audit$tf_theta_mean_after,
    legacy$audit$tf_theta_mean_after
  )
})

test_that("TF-term correspondence normalizes Tbet to Tbx21", {
  phi <- rbind(Topic1 = c(1, 1))
  colnames(phi) <- c("GENE:Tbx21", "PEAK:Tbx21")
  theta <- rbind(`C1::Tbet` = 1)
  colnames(theta) <- "Topic1"
  pair_assignment <- data.table::data.table(
    target_gene = "Tbx21",
    gene_term_id = "GENE:Tbx21",
    peak_term_id = "PEAK:Tbx21",
    gene_gammafit_topics = "1",
    peak_gammafit_topics = "1",
    assigned = TRUE,
    assigned_topic = 1L,
    assignment_status = "assigned_gammafit_maxprob_agreement"
  )

  observed <- .m3_opt_tf_theta_correspondence(
    theta = theta,
    phi = phi,
    pair_assignment = pair_assignment,
    raw_topic_ids = 1L,
    raw_to_group = 1L,
    cutoff = 0.3
  )

  expect_true(observed$available)
  expect_equal(observed$tf_term_assignments, 1L)
  expect_equal(observed$supported_tf_terms, 1L)
  expect_equal(observed$empty_tf_terms, 0L)
})

test_that("optimized theta uses the maximum raw member-topic value", {
  theta <- rbind(
    `C1::TF1` = c(0.2, 0.3, 0.5),
    `C2::TF1` = c(0.7, 0.1, 0.2)
  )
  colnames(theta) <- paste0("Topic", 1:3)
  mapping <- c(`3` = 3L, `1` = 2L, `2` = 2L)

  observed <- .m3_opt_max_theta(
    theta = theta,
    raw_topic_ids = 1:3,
    raw_to_group = mapping
  )

  expect_equal(colnames(observed), c("Topic2", "Topic3"))
  expect_equal(
    unname(observed),
    rbind(c(0.3, 0.5), c(0.7, 0.2))
  )
  expect_false(isTRUE(all.equal(observed[, "Topic2"], c(0.5, 0.8))))
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

test_that("TF-target optimization retains term-aligned links without theta gate", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2, 0.8, 0.2),
    Topic2 = c(0.2, 0.8, 0.2, 0.8)
  )
  colnames(phi) <- c(
    "GENE:A", "GENE:B", "TF1::A", "TF1::B"
  )
  theta <- rbind(
    `C1::TF1` = c(0.01, 0.99),
    `C2::TF1` = c(0.99, 0.01)
  )
  colnames(theta) <- rownames(phi)
  dtm <- Matrix::Matrix(
    matrix(
      c(5, 0, 5, 0, 0, 5, 0, 5),
      nrow = 2,
      byrow = TRUE
    ),
    sparse = TRUE,
    dimnames = list(rownames(theta), colnames(phi))
  )
  compact <- .assign_tf_target_terms_compact(
    phi,
    score_mat = phi,
    thrP = 0.7,
    min_terms = 1L
  )
  observed <- .module3_optimize_condition_topics(
    theta = theta,
    phi = phi,
    dtm = dtm,
    topic_terms = compact$topic_terms,
    pair_assignment = compact$link_assignment,
    assignment_mode = "tf_target",
    correspondence_assignment =
      .tf_target_gene_correspondence_assignment(compact$gene_assignment),
    require_theta_gate = FALSE,
    min_genes = 1L,
    min_links = 1L,
    similarity_threshold = 1,
    tf_topic_cutoff = 0.3,
    umap_max_links_per_condition = 10L,
    seed = 7L
  )

  expect_identical(observed$assignment_mode, "tf_target")
  expect_false(observed$require_theta_gate)
  expect_true(all(observed$qc$assignments$optimized_aligned))
  expect_true(any(!observed$qc$assignments$optimized_theta_pass))
  expect_setequal(observed$qc$target_levels, c("A", "B"))
  expect_true("pair_index" %in% names(observed$qc$assignments))
  expect_equal(nrow(observed$gene_assignment), 2L)
  expect_equal(nrow(observed$raw_correspondence_assignment), 2L)
  expect_setequal(
    observed$topic_terms[
      term_group == "GENE" & in_topic == TRUE,
      term_id
    ],
    c("GENE:A", "GENE:B")
  )
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

test_that("Gene-term UMAP uses assigned Gene phi terms only", {
  optimization <- list(
    raw_phi = rbind(
      Topic1 = c(0.8, 0.2, 0.6, 0.7),
      Topic2 = c(0.2, 0.8, 0.4, 0.3)
    ),
    raw_pair_assignment = data.table::data.table(
      target_gene = c("A", "B", "C"),
      assigned_topic = c(1L, 2L, NA_integer_),
      assigned = c(TRUE, TRUE, FALSE)
    )
  )
  colnames(optimization$raw_phi) <- c(
    "GENE:A", "GENE:B", "GENE:C", "PEAK:P1"
  )
  canonical_assignment <- data.table::data.table(
    term_id = c("GENE:A", "GENE:B", "PEAK:P1"),
    term_group = c("GENE", "GENE", "PEAK"),
    assigned_topic = c(1L, 2L, 1L),
    assigned = TRUE
  )

  observed <- .m3_qc_gene_term_umap_data(
    optimization,
    gene_term_assignment = canonical_assignment,
    seed = 7L
  )

  expect_identical(observed$target_gene, c("A", "B", "C"))
  expect_setequal(
    stats::na.omit(as.character(observed$topic)),
    c("Topic 1", "Topic 2")
  )
  expect_equal(nrow(observed), 3L)
  expect_false(observed[target_gene == "C", assigned])
  expect_true(all(is.finite(observed$UMAP1)))
  expect_true(all(is.finite(observed$UMAP2)))
  expect_identical(unique(observed$embedding_source), "phi")
})

test_that("Gene-term UMAP can use named condition-profile features", {
  genes <- LETTERS[1:8]
  phi <- rbind(
    Topic1 = c(rep(1, 4), rep(0, 4)),
    Topic2 = c(rep(0, 4), rep(1, 4))
  )
  colnames(phi) <- paste0("GENE:", genes)
  features <- cbind(
    Condition1 = seq(0.1, 0.8, length.out = 8),
    Condition2 = rev(seq(0.1, 0.8, length.out = 8)),
    Condition3 = rep(c(0.2, 0.8), 4)
  )
  rownames(features) <- genes
  supplied_coordinates <- cbind(
    UMAP1 = seq_along(genes),
    UMAP2 = -seq_along(genes)
  )
  rownames(supplied_coordinates) <- genes
  assignment <- data.table::data.table(
    term_id = paste0("GENE:", genes),
    term_group = "GENE",
    assigned_topic = rep(1:2, each = 4),
    assigned = TRUE
  )
  optimization <- list(
    raw_phi = phi,
    raw_topic_terms = assignment,
    gene_umap_features = features,
    gene_umap_coordinates = supplied_coordinates,
    gene_umap_feature_label = "condition-expression graph features"
  )

  observed <- .m3_qc_gene_term_umap_data(optimization, seed = 17L)

  expect_identical(observed$target_gene, genes)
  expect_identical(
    unique(observed$embedding_source),
    "condition-expression graph features"
  )
  expect_true(all(is.finite(observed$UMAP1)))
  expect_true(all(is.finite(observed$UMAP2)))
  expect_equal(observed$UMAP1, unname(supplied_coordinates[, "UMAP1"]))
  expect_equal(observed$UMAP2, unname(supplied_coordinates[, "UMAP2"]))
})

test_that("Peak-term UMAP keeps the strongest maximum-phi probabilities", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2, 9, 0.2, 0.5, 1),
    Topic2 = c(0.2, 0.8, 1, 0.8, 0.5, 9)
  )
  colnames(phi) <- c(
    "GENE:A", "GENE:B",
    "PEAK:chr1:1-2", "PEAK:chr1:3-4",
    "PEAK:chr1:5-6", "PEAK:chr1:7-8"
  )

  observed <- .m3_qc_peak_term_umap_data(
    list(raw_phi = phi),
    peak_term_assignment = data.table::data.table(
      term_id = c(
        "PEAK:chr1:1-2", "PEAK:chr1:3-4",
        "PEAK:chr1:5-6", "PEAK:chr1:7-8"
      ),
      term_group = "PEAK",
      assigned_topic = c(1L, 2L, 1L, 2L),
      assigned = TRUE
    ),
    top_n = 4L,
    seed = 19L
  )

  expect_identical(
    observed$peak_id,
    c("chr1:1-2", "chr1:7-8", "chr1:3-4", "chr1:5-6")
  )
  expect_equal(observed$max_probability, c(0.9, 0.9, 0.8, 0.5))
  expect_identical(observed$topic_num, c(1L, 2L, 2L, 1L))
  expect_identical(observed$selection_rank, 1:4)
  expect_true(all(is.finite(observed$UMAP1)))
  expect_true(all(is.finite(observed$UMAP2)))
})

test_that("Peak-term UMAP colors only GammaFit max-probability assignments", {
  phi <- rbind(
    Topic1 = c(9, 1, 6, 2),
    Topic2 = c(1, 9, 4, 8)
  )
  colnames(phi) <- paste0("PEAK:P", 1:4)
  assignment <- data.table::data.table(
    term_id = paste0("PEAK:P", 1:3),
    term_group = "PEAK",
    assigned_topic = c(2L, 1L, 2L),
    assigned = c(TRUE, FALSE, TRUE)
  )

  observed <- .m3_qc_peak_term_umap_data(
    list(raw_phi = phi),
    peak_term_assignment = assignment,
    top_n = 4L,
    seed = 29L
  )

  expect_identical(observed$peak_id, c("P1", "P2", "P4", "P3"))
  expect_identical(observed$assigned, c(TRUE, FALSE, FALSE, TRUE))
  expect_identical(observed$topic_num, c(2L, NA_integer_, NA_integer_, 2L))
  expect_identical(
    as.character(observed$topic),
    c("Topic 2", NA_character_, NA_character_, "Topic 2")
  )
})

test_that("Peak-term UMAP retains all available Peaks below its cap", {
  phi <- rbind(
    Topic1 = c(0.7, 0.3, 0.6),
    Topic2 = c(0.3, 0.7, 0.4)
  )
  colnames(phi) <- c("GENE:A", "PEAK:P1", "PEAK:P2")

  observed <- .m3_qc_peak_term_umap_data(
    list(
      raw_phi = phi,
      raw_topic_terms = data.table::data.table(
        term_id = c("PEAK:P1", "PEAK:P2"),
        term_group = "PEAK",
        assigned_topic = c(2L, 1L),
        assigned = c(TRUE, FALSE)
      )
    ),
    top_n = 10000L,
    seed = 23L
  )

  expect_identical(observed$peak_id, c("P1", "P2"))
  expect_equal(nrow(observed), 2L)
  expect_identical(observed$assigned, c(TRUE, FALSE))
  expect_error(
    .m3_qc_peak_term_umap_data(
      list(raw_phi = phi[, "GENE:A", drop = FALSE]),
      top_n = 10000L
    ),
    "No Peak terms"
  )
})

test_that("Peak-term UMAP accepts an unlimited Peak universe", {
  phi <- rbind(
    Topic1 = c(0.9, 0.4),
    Topic2 = c(0.1, 0.6)
  )
  colnames(phi) <- paste0("PEAK:P", 1:2)
  assignment <- data.table::data.table(
    term_id = rep(colnames(phi), each = 2L),
    term_group = "PEAK",
    topic_num = rep(1:2, times = 2L),
    in_topic = c(TRUE, FALSE, FALSE, TRUE)
  )

  observed <- .m3_qc_peak_term_umap_data(
    list(raw_phi = phi, raw_topic_terms = assignment),
    top_n = Inf,
    seed = 31L
  )

  expect_equal(nrow(observed), 2L)
  expect_identical(observed$selection_rank, 1:2)
  expect_setequal(observed$peak_id, c("P1", "P2"))
  expect_error(
    .m3_qc_peak_term_umap_data(
      list(raw_phi = phi, raw_topic_terms = assignment),
      top_n = 1.5
    ),
    "positive integer or Inf"
  )
  expect_error(
    .m3_qc_peak_term_umap_data(
      list(raw_phi = phi, raw_topic_terms = assignment),
      top_n = .Machine$integer.max + 1
    ),
    "positive integer or Inf"
  )
})

test_that("cross-K topic comparison detects a supported multimodal split", {
  terms <- c(paste0("GENE:G", 1:5), paste0("PEAK:P", 1:5))
  reference_phi <- rbind(
    Topic1 = c(rep(9, 4), 1, rep(9, 4), 1),
    Topic2 = c(rep(1, 4), 9, rep(1, 4), 9)
  )
  candidate_phi <- rbind(
    Topic1 = c(9, 9, 1, 1, 1, 9, 9, 1, 1, 1),
    Topic2 = c(1, 1, 9, 9, 1, 1, 1, 9, 9, 1),
    Topic3 = c(rep(1, 4), 9, rep(1, 4), 9)
  )
  colnames(reference_phi) <- colnames(candidate_phi) <- terms
  reference_assignment <- data.table::data.table(
    term_id = terms,
    term_group = rep(c("GENE", "PEAK"), each = 5L),
    topic_num = c(rep(1L, 4L), 2L, rep(1L, 4L), 2L),
    in_topic = TRUE
  )
  candidate_assignment <- data.table::data.table(
    term_id = terms,
    term_group = rep(c("GENE", "PEAK"), each = 5L),
    topic_num = c(1L, 1L, 2L, 2L, 3L, 1L, 1L, 2L, 2L, 3L),
    in_topic = TRUE
  )

  observed <- .m3_qc_compare_cross_k_topic_split(
    reference_phi,
    candidate_phi,
    reference_assignment,
    candidate_assignment
  )

  expect_identical(observed$summary$reference_topic, 1L)
  expect_identical(observed$summary$classification, "supported_split")
  expect_setequal(observed$summary$supported_children[[1L]], c(1L, 2L))
  expect_equal(
    observed$children[candidate_topic %in% 1:2, gene_fraction],
    c(0.5, 0.5)
  )
  expect_equal(
    observed$children[candidate_topic %in% 1:2, peak_fraction],
    c(0.5, 0.5)
  )
  expect_equal(nrow(observed$affinity), 2L * 3L * 2L)
})

test_that("cross-K topic comparison distinguishes persistence and partial splits", {
  terms <- c(paste0("GENE:G", 1:10), paste0("PEAK:P", 1:10))
  reference_phi <- rbind(Topic1 = rep(9, 20), Topic2 = rep(1, 20))
  candidate_phi <- rbind(
    Topic1 = c(rep(9, 18), 1, 1),
    Topic2 = c(rep(1, 18), 9, 9)
  )
  colnames(reference_phi) <- colnames(candidate_phi) <- terms
  reference_assignment <- data.table::data.table(
    term_id = terms,
    term_group = rep(c("GENE", "PEAK"), each = 10L),
    topic_num = 1L,
    in_topic = TRUE
  )
  persistent_assignment <- data.table::copy(reference_assignment)
  partial_assignment <- data.table::copy(reference_assignment)
  partial_assignment[term_id %in% c("GENE:G10", "PEAK:P10"), topic_num := 2L]

  persistent <- .m3_qc_compare_cross_k_topic_split(
    reference_phi,
    candidate_phi,
    reference_assignment,
    persistent_assignment
  )
  partial <- .m3_qc_compare_cross_k_topic_split(
    reference_phi,
    candidate_phi,
    reference_assignment,
    partial_assignment
  )

  expect_identical(persistent$summary$classification, "persistent_topic")
  expect_identical(partial$summary$classification, "partial_split")
})

test_that("topic UMAP continuation pages use available panel columns", {
  gene_umap <- data.table::data.table(
    target_gene = paste0("G", 1:4),
    topic_num = c(31L, 31L, 32L, 32L),
    assigned = TRUE,
    topic = factor(
      c("Topic 31", "Topic 31", "Topic 32", "Topic 32"),
      levels = c("Topic 31", "Topic 32")
    ),
    UMAP1 = c(-1, -0.5, 0.5, 1),
    UMAP2 = c(-1, -0.5, 0.5, 1)
  )
  observed <- .m3_qc_gene_term_umap_plot(
    optimization = list(),
    topic_palette = c(`Topic 31` = "#3366AA", `Topic 32` = "#D95F02"),
    gene_umap = gene_umap,
    topic_ids = 31:32
  )

  expect_equal(observed$facet$params$ncol, 2L)
  expect_equal(observed$facet$params$nrow, 1L)
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

test_that("Gene-only retention exposes the genes used at every stage", {
  optimization <- list(
    assignment_mode = "gene_only",
    raw_pair_assignment = data.table::data.table(
      target_gene = c("A", "B", "C", "D"),
      gene_gammafit_topics = c("1", "2", "", "1"),
      assigned_topic = c(1L, 2L, NA_integer_, 1L)
    ),
    qc = list(
      raw_topic_ids = 1:2,
      assignments = data.table::data.table(
        tf_index = 1L,
        target_index = 1:4,
        pair_index = 1:4,
        raw_aligned = c(TRUE, FALSE, TRUE, FALSE)
      )
    )
  )

  observed <- .m3_qc_gene_only_retention_gene_sets(optimization)

  expect_identical(observed$total, c("A", "B", "C", "D"))
  expect_identical(observed$gamma, c("A", "B", "D"))
  expect_identical(observed$max, c("A", "B", "D"))
  expect_identical(observed$tf, c("A", "C"))
  expect_silent(.m3_qc_retention_page(
    optimization,
    standardized = TRUE,
    peak_counts = c(10, 8, 7, 5)
  ))
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

test_that("condition-document correlation uses raw theta directly", {
  theta <- rbind(
    C2 = c(0.7, 0.2, 0.1),
    C1 = c(0.1, 0.2, 0.7),
    C3 = c(0.2, 0.6, 0.2)
  )
  colnames(theta) <- paste0("Topic", 1:3)
  optimization <- list(
    document_design = "condition",
    theta = theta,
    qc = list(assignments = data.table::data.table(
      doc_index = c(1L, 2L),
      condition_id = c("wrong_1", "wrong_2")
    ))
  )

  observed_theta <- .m3_qc_condition_topic_theta(optimization)
  expected_theta <- theta[c("C1", "C2", "C3"), , drop = FALSE]
  colnames(expected_theta) <- paste0("Topic ", 1:3)
  expect_equal(observed_theta, expected_theta)

  observed <- .m3_qc_theta_condition_correlation(observed_theta)
  expected <- stats::cor(t(expected_theta))
  expect_equal(observed, expected)
  expect_identical(rownames(observed), c("C1", "C2", "C3"))
})

test_that("condition::TF correlation uses TF documents shared by all conditions", {
  theta <- rbind(
    `C1::TF1` = c(0.8, 0.2),
    `C1::TF2` = c(0.1, 0.9),
    `C2::TF1` = c(0.7, 0.3),
    `C2::TF3` = c(0.2, 0.8)
  )
  colnames(theta) <- paste0("Topic", 1:2)

  observed <- .m3_qc_condition_tf_topic_theta(theta, matched_tfs = TRUE)

  expect_equal(
    unname(observed[, , drop = FALSE]),
    structure(rbind(c(0.8, 0.2), c(0.7, 0.3)), dimnames = NULL),
    ignore_attr = TRUE
  )
  expect_identical(rownames(observed), c("C1", "C2"))
  expect_identical(colnames(observed), c("Topic 1", "Topic 2"))
  expect_true(isTRUE(attr(observed, "matched_tfs")))
  expect_identical(attr(observed, "tf_count"), 1L)
})

test_that("pathway Gene sets retain expressed genes unique to one pathway", {
  sets <- data.table::data.table(
    pathway = c("P1", "P1", "P2", "P2", "P3"),
    gene = c("A", "B", "B", "C", "D")
  )

  observed <- .m3_qc_prepare_pathway_gene_sets(
    sets,
    expressed_genes = c("A", "B", "C")
  )

  expect_equal(observed$pathway, c("P1", "P2"))
  expect_equal(observed$gene, c("A", "C"))
  expect_identical(attr(observed, "pathway_order"), c("P1", "P2"))
})

test_that("dense heatmap counts use compact integer labels", {
  expect_equal(
    .m3_qc_heatmap_count(c(999, 1157, 29640, 1500000)),
    c("999", "1k", "30k", "2M")
  )
})

test_that("condition-document QC counts only confident primary TFs", {
  optimization <- list(
    raw_to_optimized = stats::setNames(1:2, 1:2),
    qc = list(
      raw_counts = data.table::data.table(
        raw_topic = 1:2,
        links = c(20L, 10L),
        genes = c(8L, 5L),
        tfs = c(99L, 99L)
      ),
      optimized_counts = data.table::data.table(
        optimized_topic = 1:2,
        links = c(20L, 10L),
        genes = c(8L, 5L),
        tfs = c(99L, 99L)
      ),
      condition_topic = data.table::CJ(
        condition_id = c("A", "B"),
        optimized_topic = 1:2
      )[, `:=`(links = 1L, genes = 1L, tfs = 99L)]
    )
  )
  evidence <- data.table::data.table(
    condition_id = rep(c("A", "B"), each = 6L),
    tf = rep(rep(c("TF1", "TF2", "TF3"), each = 2L), 2L),
    topic_num = rep(1:2, 6L),
    global_primary_topic = rep(c(1L, 1L, 1L, 1L, 2L, 2L), 2L),
    global_primary_confident = rep(
      c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
      2L
    ),
    condition_primary_topic = c(
      1L, 1L, 1L, 1L, 2L, 2L,
      2L, 2L, 1L, 1L, 2L, 2L
    ),
    condition_primary_confident = c(
      TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
      TRUE, TRUE, FALSE, FALSE, TRUE, TRUE
    )
  )

  observed <- .m3_qc_apply_primary_tf_counts(optimization, evidence)

  expect_equal(observed$qc$raw_counts$tfs, c(1L, 1L))
  expect_equal(observed$qc$optimized_counts$tfs, c(1L, 1L))
  expect_equal(observed$qc$raw_counts$links, c(20L, 10L))
  expect_equal(
    observed$qc$condition_topic[
      order(condition_id, optimized_topic),
      tfs
    ],
    c(1L, 1L, 0L, 2L)
  )
})

test_that("TF-target Gene sets use supported positive Module 2 links", {
  correlations <- data.table::data.table(
    tf = c("TF1", "TF1", "TF1", "TF2", "TF3", "TF4"),
    target_gene = c("A", "B", "C", "A", "B", "D"),
    best_r = c(0.90, 0.80, -0.90, 0.70, 0.95, 0.60),
    best_method = "pearson",
    best_fdr = 0.01,
    pass = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  links <- data.table::data.table(
    tf = c("TF1", "TF1", "TF1", "TF2", "TF4", "TF4"),
    target_gene = c("A", "A", "B", "A", "D", "D"),
    fp_id = c("P1", "P2", "P3", "P4", "P5", "P5"),
    module2_link_pass = TRUE
  )

  observed <- .m3_qc_prepare_tf_target_gene_sets(
    correlations,
    links,
    gene_universe = LETTERS[1:4],
    min_r = 0.5,
    top_n_tfs = 3L,
    top_n_targets = 2L
  )

  expect_equal(unique(observed$tf), c("TF1", "TF2", "TF4"))
  expect_equal(observed[tf == "TF1", target_gene], c("A", "B"))
  expect_equal(
    observed[tf == "TF1" & target_gene == "A", supporting_peak_count],
    2L
  )
  expect_false(any(observed$target_gene == "C"))
  expect_equal(observed, .m3_qc_prepare_tf_target_gene_sets(
    correlations[sample(.N)],
    links[sample(.N)],
    gene_universe = LETTERS[1:4],
    min_r = 0.5,
    top_n_tfs = 3L,
    top_n_targets = 2L
  ))
})

test_that("TF-target QC selection covers topics with a curated TF panel", {
  expect_identical(
    formals(.m3_qc_select_topic_tf_target_gene_sets)$top_n_targets,
    200L
  )
  expect_identical(
    formals(.write_module3_topic_assignment_qc)$tf_target_top_n,
    200L
  )
  expect_identical(
    formals(.write_module3_topic_assignment_qc)$tf_target_top_topics,
    3L
  )
  eligible <- data.table::data.table(
    tf = rep(c("TF_A", "TF_B", "TF_C"), each = 4L),
    target_gene = rep(c("G1", "G2", "G3", "G4"), 3L),
    best_r = c(
      0.9, 0.8, 0.6, 0.5,
      0.7, 0.6, 0.9, 0.8,
      0.8, 0.7, 0.7, 0.6
    ),
    supporting_peak_count = 1L
  )
  assignment <- data.table::data.table(
    target_gene = c("G1", "G2", "G3", "G4"),
    assigned_topic = c(1L, 1L, 2L, 2L),
    assigned = TRUE
  )
  panel <- data.table::data.table(
    tf = c("TF_A", "TF_B", "TF_C"),
    tf_function = c("Program A", "Program B", "Program C")
  )

  observed <- .m3_qc_select_topic_tf_target_gene_sets(
    eligible,
    assignment,
    tf_panel = panel,
    top_n_tfs = 3L,
    top_n_targets = 1L
  )

  expect_equal(data.table::uniqueN(observed$tf), 3L)
  expect_equal(sort(unique(observed$selected_topic)), 1:2)
  expect_equal(observed[, .N, by = tf]$N, rep(1L, 3L))
  expect_setequal(unique(observed$tf_function), panel$tf_function)

  all_topics <- .m3_qc_select_topic_tf_target_gene_sets(
    eligible,
    assignment,
    tf_panel = panel,
    top_n_tfs = 3L,
    top_n_targets = 4L,
    include_all_target_topics = TRUE
  )
  expect_equal(all_topics[, .N, by = tf]$N, rep(4L, 3L))
  expect_equal(all_topics[, data.table::uniqueN(selected_topic), by = tf]$V1,
               rep(2L, 3L))
  expect_equal(all_topics[, data.table::uniqueN(panel_topic), by = tf]$V1,
               rep(1L, 3L))
})

test_that("TF-target QC uses one evidence-ranked display pool", {
  eligible <- data.table::data.table(
    tf = rep(c("TF_A", "TF_B"), each = 6L),
    target_gene = c(paste0("G", 1:6), paste0("H", 1:6)),
    best_r = c(seq(0.99, 0.94, by = -0.01), seq(0.93, 0.88, by = -0.01)),
    supporting_peak_count = rep(6:1, 2L)
  )
  assignment <- data.table::data.table(
    target_gene = c(paste0("G", 1:6), paste0("H", 1:6)),
    assigned_topic = rep(c(1L, 1L, 2L, 2L, 3L, 3L), 2L),
    assigned = TRUE
  )

  observed <- .m3_qc_select_topic_tf_target_gene_sets(
    eligible,
    assignment,
    top_n_tfs = 2L,
    top_n_targets = 3L,
    include_all_target_topics = TRUE
  )

  expect_equal(observed[, .N, by = tf]$N, rep(3L, 2L))
  expect_setequal(
    observed[tf == "TF_A", target_gene],
    paste0("G", 1:3)
  )
  expect_setequal(
    observed[tf == "TF_B", target_gene],
    paste0("H", 1:3)
  )
  expect_true(all(observed[
    , any(selected_topic == panel_topic),
    by = tf
  ]$V1))
})

test_that("TF-target QC can select the union of top TFs per topic", {
  eligible <- data.table::data.table(
    tf = c(
      rep("TF_A", 4L), rep("TF_B", 3L), rep("TF_C", 2L),
      rep("TF_D", 4L), rep("TF_E", 3L), rep("TF_F", 2L)
    ),
    target_gene = c(
      paste0("A", 1:4), paste0("B", 1:3), paste0("C", 1:2),
      paste0("D", 1:4), paste0("E", 1:3), paste0("F", 1:2)
    ),
    best_r = 0.8,
    supporting_peak_count = 1L
  )
  assignment <- data.table::data.table(
    target_gene = eligible$target_gene,
    assigned_topic = c(rep(1L, 9L), rep(2L, 9L)),
    assigned = TRUE
  )

  observed <- .m3_qc_select_topic_tf_target_gene_sets(
    eligible,
    assignment,
    top_n_tfs = 4L,
    top_n_targets = 10L,
    include_all_target_topics = TRUE,
    selection_mode = "top_per_topic",
    top_tfs_per_topic = 2L
  )

  expect_setequal(unique(observed$tf), c("TF_A", "TF_B", "TF_D", "TF_E"))
  expect_equal(data.table::uniqueN(observed$display_rank), 4L)
})

test_that("TF-target UMAP highlights the three largest target topics", {
  gene_umap <- data.table::data.table(
    target_gene = paste0("G", 1:8),
    UMAP1 = seq(-1, 1, length.out = 8),
    UMAP2 = rev(seq(-1, 1, length.out = 8)),
    topic_num = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L)
  )
  targets <- data.table::data.table(
    display_rank = 1L,
    tf = "TF_A",
    selected_topic = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L),
    target_gene = paste0("G", 1:8),
    best_r = seq(0.9, 0.55, by = -0.05)
  )

  plot <- .m3_qc_tf_target_gene_umap_plot(gene_umap, targets)
  built <- ggplot2::ggplot_build(plot)

  expect_equal(nrow(built$data[[1L]]), 8L)
  expect_equal(nrow(built$data[[2L]]), 7L)
  expect_identical(
    plot$facet$params$labeller(data.frame(tf = factor("TF_A")))$tf,
    "TF_A\nT1/T2/T3 (3/2/2)"
  )
})

test_that("compact topic QC writes Gene, Peak, pathway, TF-target, and structure sections", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2, 0.7, 0.3),
    Topic2 = c(0.2, 0.8, 0.3, 0.7)
  )
  colnames(phi) <- c("GENE:A", "GENE:B", "PEAK:A", "PEAK:B")
  optimization <- list(
    document_design = "condition_tf",
    assignment_mode = "gene_peak",
    raw_phi = phi,
    raw_to_optimized = stats::setNames(1:2, 1:2),
    qc = list(
      raw_topic_ids = 1:2,
      optimized_topic_ids = 1:2,
      raw_topic_similarity = .m3_opt_hellinger_similarity(
        phi,
        1:2,
        3:4
      ),
      raw_counts = data.table::data.table(
        raw_topic = 1:2,
        peaks = c(1L, 1L),
        genes = c(1L, 1L)
      )
    )
  )
  assignment <- data.table::data.table(
    target_gene = c("A", "B"),
    assigned_topic = 1:2,
    assigned = TRUE
  )
  peak_assignment <- data.table::data.table(
    term_id = c("PEAK:A", "PEAK:B"),
    term_group = "PEAK",
    assigned_topic = 1:2,
    assigned = TRUE
  )
  pathway_gene_sets <- data.table::data.table(
    pathway = c("Pathway A", "Pathway B"),
    gene = c("A", "B")
  )
  tf_target_gene_sets <- data.table::data.table(
    tf = c("TF1", "TF2"),
    target_gene = c("A", "B"),
    best_r = c(0.8, 0.7),
    supporting_peak_count = 1L
  )
  output <- tempfile(fileext = ".pdf")

  expect_silent(.write_module3_topic_assignment_qc(
    optimization,
    out_file = output,
    gene_term_assignment = assignment,
    peak_term_assignment = peak_assignment,
    pathway_gene_sets = pathway_gene_sets,
    tf_target_gene_sets = tf_target_gene_sets,
    sections = c(
      "gene_phi_umap",
      "peak_phi_umap",
      "pathway_gene_umap",
      "tf_target_gene_umap",
      "raw_structure"
    ),
    peak_umap_top_n = 2L,
    sidebar_mode = "terms",
    seed = 7L
  ))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 1000)
  qpdf <- Sys.which("qpdf")
  if (nzchar(qpdf)) {
    expect_identical(
      system2(qpdf, c("--show-npages", output), stdout = TRUE),
      "5"
    )
  }
  expect_error(
    .write_module3_topic_assignment_qc(
      optimization,
      out_file = output,
      sections = "later_pages"
    ),
    "must be 'standard', 'all', or a subset"
  )
})

test_that("standard topic QC is the default compact input-design report", {
  phi <- rbind(
    Topic1 = c(0.7, 0.6, 0.2, 0.1, 0.65, 0.55, 0.25, 0.15),
    Topic2 = c(0.3, 0.4, 0.8, 0.9, 0.35, 0.45, 0.75, 0.85)
  )
  colnames(phi) <- c(
    paste0("GENE:", LETTERS[1:4]),
    paste0("PEAK:", LETTERS[1:4])
  )
  theta <- rbind(
    `C1::TF1` = c(0.8, 0.2),
    `C1::TF2` = c(0.6, 0.4),
    `C2::TF1` = c(0.3, 0.7),
    `C2::TF2` = c(0.1, 0.9)
  )
  pairs <- data.table::data.table(
    target_gene = LETTERS[1:4],
    gene_term_id = paste0("GENE:", LETTERS[1:4]),
    peak_term_id = paste0("PEAK:", LETTERS[1:4]),
    gene_gammafit_topics = c("1", "1", "2", "2"),
    peak_gammafit_topics = c("1", "1", "2", "2"),
    assigned_topic = c(1L, 1L, 2L, 2L),
    assigned = TRUE
  )
  assignments <- data.table::data.table(
    doc_index = rep(1:4, each = 4L),
    pair_index = rep(1:4, times = 4L),
    target_index = rep(1:4, times = 4L),
    tf_index = rep(c(1L, 2L, 1L, 2L), each = 4L),
    condition_id = rep(c("C1", "C1", "C2", "C2"), each = 4L),
    raw_aligned = rep(c(TRUE, TRUE, FALSE, TRUE), times = 4L)
  )
  optimization <- list(
    document_design = "condition_tf",
    assignment_mode = "gene_peak",
    theta = theta,
    raw_phi = phi,
    raw_pair_assignment = pairs,
    raw_to_optimized = stats::setNames(1:2, 1:2),
    qc = list(
      raw_topic_ids = 1:2,
      optimized_topic_ids = 1:2,
      raw_topic_similarity = .m3_opt_hellinger_similarity(phi, 1:2, 5:8),
      raw_counts = data.table::data.table(
        raw_topic = 1:2,
        peaks = c(2L, 2L),
        genes = c(2L, 2L)
      ),
      assignments = assignments
    )
  )
  gene_assignment <- data.table::data.table(
    target_gene = LETTERS[1:4],
    assigned_topic = c(1L, 1L, 2L, 2L),
    assigned = TRUE
  )
  peak_assignment <- data.table::data.table(
    term_id = paste0("PEAK:", LETTERS[1:4]),
    term_group = "PEAK",
    assigned_topic = c(1L, 1L, 2L, 2L),
    assigned = TRUE
  )
  tf_target_gene_sets <- data.table::data.table(
    tf = c("TF1", "TF1", "TF2", "TF2"),
    target_gene = LETTERS[1:4],
    best_r = c(0.9, 0.8, 0.7, 0.6),
    supporting_peak_count = 1L
  )
  output_gene_peak <- tempfile(fileext = ".pdf")
  output_gene_only <- tempfile(fileext = ".pdf")

  expect_identical(formals(.write_module3_topic_assignment_qc)$sections,
                   "standard")
  expect_true(is.infinite(
    eval(formals(.write_module3_topic_assignment_qc)$peak_umap_top_n)
  ))
  expect_silent(.write_module3_topic_assignment_qc(
    optimization,
    out_file = output_gene_peak,
    gene_term_assignment = gene_assignment,
    peak_term_assignment = peak_assignment,
    tf_target_gene_sets = tf_target_gene_sets,
    sidebar_mode = "terms",
    seed = 13L
  ))

  gene_only <- optimization
  gene_only$assignment_mode <- "gene_only"
  gene_only$raw_phi <- phi[, 1:4, drop = FALSE]
  gene_only$raw_pair_assignment[, peak_gammafit_topics := NULL]
  gene_only$qc$raw_topic_similarity <- diag(2L)
  gene_only$qc$raw_counts <- data.table::data.table(
    raw_topic = 1:2,
    genes = c(2L, 2L)
  )
  expect_silent(.write_module3_topic_assignment_qc(
    gene_only,
    out_file = output_gene_only,
    gene_term_assignment = gene_assignment,
    tf_target_gene_sets = tf_target_gene_sets,
    sidebar_mode = "terms",
    seed = 17L
  ))

  expect_gt(file.info(output_gene_peak)$size, 1000)
  expect_gt(file.info(output_gene_only)$size, 1000)
  qpdf <- Sys.which("qpdf")
  if (nzchar(qpdf)) {
    expect_identical(
      system2(qpdf, c("--show-npages", output_gene_peak), stdout = TRUE),
      "6"
    )
    expect_identical(
      system2(qpdf, c("--show-npages", output_gene_only), stdout = TRUE),
      "5"
    )
  }
})

test_that("Gene-only topics project Peaks without target-count duplication", {
  theta <- rbind(
    `C1::TF1` = c(Topic1 = 0.8, Topic2 = 0.2),
    `C1::TF2` = c(Topic1 = 0.2, Topic2 = 0.8)
  )
  phi <- rbind(
    Topic1 = c(`GENE:A` = 0.9, `GENE:B` = 0.1),
    Topic2 = c(`GENE:A` = 0.1, `GENE:B` = 0.9)
  )
  edges <- data.table::data.table(
    condition_id = "C1",
    doc_id = c("C1::TF1", "C1::TF1", "C1::TF1", "C1::TF2"),
    tf = c("TF1", "TF1", "TF1", "TF2"),
    peak_id = "P1",
    target_gene = c("A", "B", "B", "A"),
    tf_target_r = c(1, 0.5, 0.5, 1),
    fp_target_r = 1,
    tf_peak_r = 1,
    fp_score = 1,
    tf_expression = 1
  )

  projected <- .m3_project_peaks_from_gene_topics(
    theta,
    phi,
    edges,
    chunk_size = 2L
  )
  expect_equal(projected$audit$input_rows, 4L)
  expect_equal(projected$audit$deduplicated_rows, 3L)
  expect_equal(projected$audit$tf_peak_groups, 2L)
  expect_equal(projected$tf_peak_assignments$n_targets, c(2L, 1L))
  probability <- as.matrix(projected$condition_peak_probabilities[, .(
    Topic1,
    Topic2
  )])
  expect_equal(rowSums(probability), 1, tolerance = 1e-12)

  gene_probability <- .m3_opt_row_normalize(t(phi))
  tf1_a <- .m3_opt_row_normalize(
    theta[1, , drop = FALSE] * gene_probability["GENE:A", , drop = FALSE]
  )
  tf1_b <- .m3_opt_row_normalize(
    theta[1, , drop = FALSE] * gene_probability["GENE:B", , drop = FALSE]
  )
  tf2_a <- .m3_opt_row_normalize(
    theta[2, , drop = FALSE] * gene_probability["GENE:A", , drop = FALSE]
  )
  tf1 <- (tf1_a + 0.5 * tf1_b) / 1.5
  expected <- (0.75 * tf1 + tf2_a) / 1.75
  expect_equal(as.numeric(probability[1, ]), as.numeric(expected), tolerance = 1e-12)
})

test_that("condition Peak probabilities aggregate by evidence", {
  condition_peaks <- data.table::data.table(
    condition_id = c("C1", "C2", "C1"),
    peak_id = c("P1", "P1", "P2"),
    projection_evidence = c(2, 1, 1),
    Topic1 = c(0.75, 0.25, 0.1),
    Topic2 = c(0.25, 0.75, 0.9)
  )
  aggregated <- .m3_aggregate_projected_peak_topics(condition_peaks)
  expect_equal(aggregated[peak_id == "P1", Topic1], 7 / 12)
  expect_equal(aggregated[peak_id == "P1", Topic2], 5 / 12)
  expect_equal(aggregated[peak_id == "P1", primary_topic], 1L)
  expect_equal(aggregated[peak_id == "P2", primary_topic], 2L)
  expect_equal(aggregated[, Topic1 + Topic2], c(1, 1))
})

test_that("compact condition::TF QC adds a matched-theta correlation page", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2, 0.7, 0.3),
    Topic2 = c(0.2, 0.8, 0.3, 0.7)
  )
  colnames(phi) <- c("GENE:A", "GENE:B", "PEAK:A", "PEAK:B")
  theta <- rbind(
    `C1::TF1` = c(0.8, 0.2),
    `C1::TF2` = c(0.4, 0.6),
    `C2::TF1` = c(0.2, 0.8),
    `C2::TF3` = c(0.6, 0.4)
  )
  optimization <- list(
    document_design = "condition_tf",
    assignment_mode = "gene_peak",
    theta = theta,
    raw_phi = phi,
    raw_to_optimized = stats::setNames(1:2, 1:2),
    qc = list(
      raw_topic_ids = 1:2,
      optimized_topic_ids = 1:2,
      raw_topic_similarity = .m3_opt_hellinger_similarity(phi, 1:2, 3:4),
      raw_counts = data.table::data.table(
        raw_topic = 1:2,
        peaks = c(1L, 1L),
        genes = c(1L, 1L)
      )
    )
  )
  assignment <- data.table::data.table(
    target_gene = c("A", "B"),
    assigned_topic = 1:2,
    assigned = TRUE
  )
  pathway_gene_sets <- data.table::data.table(
    pathway = c("Pathway A", "Pathway B"),
    gene = c("A", "B")
  )
  output <- tempfile(fileext = ".pdf")

  expect_silent(.write_module3_topic_assignment_qc(
    optimization,
    out_file = output,
    gene_term_assignment = assignment,
    pathway_gene_sets = pathway_gene_sets,
    sections = c(
      "gene_phi_umap",
      "pathway_gene_umap",
      "raw_structure",
      "condition_correlation"
    ),
    sidebar_mode = "terms",
    seed = 9L
  ))
  qpdf <- Sys.which("qpdf")
  if (nzchar(qpdf)) {
    expect_identical(
      system2(qpdf, c("--show-npages", output), stdout = TRUE),
      "4"
    )
  }
})

test_that("gene-only topic QC uses a target-gene-only structure sidebar", {
  phi <- rbind(
    Topic1 = c(0.8, 0.2),
    Topic2 = c(0.2, 0.8)
  )
  colnames(phi) <- c("GENE:A", "GENE:B")
  optimization <- list(
    document_design = "condition",
    assignment_mode = "gene",
    raw_phi = phi,
    raw_to_optimized = stats::setNames(1:2, 1:2),
    qc = list(
      raw_topic_ids = 1:2,
      optimized_topic_ids = 1:2,
      raw_topic_similarity = diag(2),
      raw_counts = data.table::data.table(
        raw_topic = 1:2,
        genes = c(1L, 1L)
      )
    )
  )
  assignment <- data.table::data.table(
    target_gene = c("A", "B"),
    assigned_topic = 1:2,
    assigned = TRUE
  )
  output <- tempfile(fileext = ".pdf")

  expect_silent(.write_module3_topic_assignment_qc(
    optimization,
    out_file = output,
    gene_term_assignment = assignment,
    sections = c("gene_phi_umap", "raw_structure"),
    sidebar_mode = "auto",
    seed = 11L
  ))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 1000)
})
