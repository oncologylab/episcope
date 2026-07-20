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
      qc_top_tfs = 100L,
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

test_that("TF heatmaps use one page per condition", {
  assignments <- data.table::data.table(
    condition_id = "Condition_A",
    tf_index = rep(1:4, each = 2),
    target_index = seq_len(8),
    optimized_topic = rep(1:2, 4),
    optimized_aligned = TRUE
  )
  optimization <- list(qc = list(
    assignments = assignments,
    optimized_topic_ids = 1:2,
    tf_levels = paste0("TF", 1:4)
  ))

  pages <- .m3_qc_tf_topic_pages(
    optimization,
    top_n_tfs = 4
  )

  expect_length(pages, 1)
  expect_false(anyNA(pages[[1]]$data$row_label))
  expect_equal(length(unique(pages[[1]]$data$row_label)), 4)
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
