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
    similarity_threshold = 1
  )

  expect_equal(result$mapping[[1]], 2)
  expect_equal(result$mapping[[2]], 2)
  expect_equal(result$audit$source_topic[[1]], 1)
  expect_equal(result$audit$representative_topic[[1]], 2)
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
      qc_seed = 1L
    )
  )))
  expect_error(
    .validate_project_config(list(
      module3 = list(merge_similarity_threshold = 1.1)
    )),
    "at most 1"
  )
})

test_that("repeated condition IDs retain concise stable display labels", {
  observed <- .m3_qc_short_condition_labels(c(
    "HPAFII_0_FBS_Ctrl",
    "HPAFII_0_FBS_Ctrl",
    "HPAFII_10_FBS_Ctrl"
  ))

  expect_equal(observed, c("0 FBS", "0 FBS", "10 FBS"))
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

  expect_equal(dim(matrix), c(150, 2))
  expect_equal(unname(matrix["TF1", ]), c(1, 1))
  expect_length(pages, 1)
  expect_s3_class(pages[[1]], "gtable")
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

test_that("QC panel alignment reserves one shared legend width", {
  values <- data.frame(x = 1:2, y = 1, value = 1:2)
  short <- ggplot2::ggplot(values, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_continuous(name = "Count")
  long <- ggplot2::ggplot(values, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_continuous(name = "Mean log2(expression + 1)")

  aligned <- .m3_qc_align_plot_widths(list(short, long))

  expect_equal(aligned[[1L]]$widths, aligned[[2L]]$widths)
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
    genes = c(20, 10, 5)
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
