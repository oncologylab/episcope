test_that("condition background weighting penalizes uniform terms", {
  values <- data.table::data.table(
    condition_id = rep(c("C1", "C2", "C3"), each = 3L),
    feature_id = rep(c("uniform", "variable", "rare"), 3L),
    modality = "gene",
    value = c(10, 30, 20, 10, 2, 0, 10, 1, 0)
  )[value > 0]

  weighted <- craftgrn:::.m3_condition_background_weights(
    values,
    method = "idf_entropy",
    temperature = 0.5,
    floor = 0.1
  )
  specificity <- unique(weighted[, .(
    feature_id,
    background_specificity,
    background_multiplier
  )])

  expect_equal(
    specificity[feature_id == "uniform", background_multiplier],
    0.1,
    tolerance = 1e-8
  )
  expect_gt(
    specificity[feature_id == "variable", background_multiplier],
    specificity[feature_id == "uniform", background_multiplier]
  )
  expect_equal(
    specificity[feature_id == "rare", background_multiplier],
    1,
    tolerance = 1e-8
  )
})

test_that("condition modality balancing preserves equal target mass", {
  values <- data.table::data.table(
    condition_id = rep(c("C1", "C2"), each = 4L),
    feature_id = rep(c("G1", "G2", "P1", "P2"), 2L),
    modality = rep(c("gene", "gene", "peak", "peak"), 2L),
    value = c(10, 4, 100, 20, 3, 8, 40, 80)
  )
  weighted <- craftgrn:::.m3_condition_background_weights(
    values,
    method = "idf_entropy",
    floor = 0.1
  )
  balanced <- craftgrn:::.m3_balance_condition_modalities(
    weighted,
    target_mass = 500
  )
  totals <- balanced[, .(mass = sum(weight)),
    by = .(condition_id, modality)]

  expect_equal(totals$mass, rep(500, nrow(totals)), tolerance = 1e-8)
})

test_that("condition profile matrices select modalities deterministically", {
  values <- data.table::data.table(
    condition_id = rep(c("C2", "C1"), each = 2L),
    feature_id = c("G1", "P1", "G1", "P1"),
    modality = rep(c("gene", "peak"), 2L),
    weight = c(2, 4, 1, 3)
  )

  gene <- craftgrn:::.m3_condition_profile_matrix(values, "gene")
  peak <- craftgrn:::.m3_condition_profile_matrix(values, "peak")

  expect_equal(dim(gene), c(1L, 2L))
  expect_equal(dim(peak), c(1L, 2L))
  expect_identical(colnames(gene), c("C1", "C2"))
  expect_identical(colnames(peak), c("C1", "C2"))
  expect_equal(unname(gene[1, ]), c(1, 2))
  expect_equal(unname(peak[1, ]), c(3, 4))
})

test_that("SNN edges retain mutual shared-neighbor structure", {
  profile <- rbind(
    A1 = c(1, 0, 0),
    A2 = c(0.95, 0.05, 0),
    A3 = c(0.9, 0.1, 0),
    B1 = c(0, 0, 1),
    B2 = c(0, 0.05, 0.95),
    B3 = c(0, 0.1, 0.9)
  )
  edges <- craftgrn:::.m3_snn_edges(profile, neighbors = 2L, prune = 0)

  expect_true(nrow(edges) > 0L)
  expect_true(all(edges$from < edges$to))
  expect_true(all(edges$shared_neighbors >= 2L))
  expect_true(all(edges$weight > 0 & edges$weight <= 1))
})

test_that("Leiden resolution paths reuse one screen across requested K", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_ring(12L)
  igraph::E(graph)$weight <- 1
  path <- craftgrn:::.m3_leiden_resolution_path(
    graph,
    target_topics = c(2L, 3L),
    seed = 9L,
    resolution_grid = c(0.1, 0.5, 1),
    screen_iterations = 1L,
    final_iterations = 1L
  )

  expect_identical(names(path), c("2", "3"))
  expect_length(path[["2"]]$membership, 12L)
  expect_true(path[["2"]]$communities >= 1L)
  expect_equal(nrow(path[["2"]]$resolution_summary), 3L)
})

test_that("Leiden consensus returns exact requested community counts", {
  skip_if_not_installed("igraph")
  graph <- igraph::make_ring(20L)
  igraph::V(graph)$name <- paste0("V", seq_len(20L))
  igraph::E(graph)$weight <- 1
  profile <- cbind(
    x = cos(seq_len(20L) * pi / 10),
    y = sin(seq_len(20L) * pi / 10)
  )
  rownames(profile) <- igraph::V(graph)$name
  path <- craftgrn:::.m3_leiden_consensus_path(
    graph,
    vertex_profile = profile,
    target_topics = c(2L, 3L),
    seed = 7L,
    base_resolution = 1,
    leiden_iterations = 2L
  )

  expect_equal(data.table::uniqueN(path[["2"]]$membership), 2L)
  expect_equal(data.table::uniqueN(path[["3"]]$membership), 3L)
  expect_identical(names(path[["2"]]$membership), igraph::V(graph)$name)
})

test_that("WarpLDA alpha sum controls per-topic alpha", {
  skip_if_not(exists(".craftgrn_warplda_fit_cpp", asNamespace("craftgrn")))
  dtm <- Matrix::Matrix(
    matrix(c(4, 1, 0, 1, 4, 0, 0, 1, 4), nrow = 3, byrow = TRUE),
    sparse = TRUE
  )
  dtm <- methods::as(dtm, "dgCMatrix")
  rownames(dtm) <- paste0("D", 1:3)
  colnames(dtm) <- paste0("T", 1:3)
  out <- craftgrn:::run_warplda_models(
    dtm,
    K_grid = 3L,
    iterations = 10L,
    alpha_sum = 6,
    beta = 0.1,
    seed = 7L,
    sampler = "gibbs_sync",
    workers = 1L,
    threads_per_model = 1L,
    verbose = FALSE
  )

  expect_equal(out$metrics$alpha, 2)
  expect_equal(out$metrics$alpha_sum, 6)
})

test_that("condition-document QC uses the standard assignment report", {
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("scales")
  skip_if_not_installed("uwot")
  set.seed(11)
  phi <- matrix(stats::rexp(4 * 80), nrow = 4)
  theta <- matrix(stats::rexp(6 * 4), nrow = 6)
  phi <- phi / rowSums(phi)
  theta <- theta / rowSums(theta)
  rownames(phi) <- paste0("Topic", seq_len(nrow(phi)))
  colnames(phi) <- c(
    paste0("GENE:G", seq_len(40)),
    paste0("PEAK:P", seq_len(40))
  )
  rownames(theta) <- paste0("Condition_", seq_len(nrow(theta)))
  colnames(theta) <- rownames(phi)
  assignment <- data.table::data.table(
    term_id = colnames(phi),
    term_group = rep(c("GENE", "PEAK"), each = 40),
    assigned_topic = max.col(t(phi), ties.method = "first"),
    phi = apply(phi, 2L, max)
  )
  link_universe <- data.table::CJ(
    condition_id = rownames(theta),
    tf = c("TF1", "TF2", "TF3"),
    pair = seq_len(40)
  )
  link_universe[, `:=`(
    peak_id = paste0("P", pair),
    gene_key = paste0("G", pair),
    fp_score = 1 + pair / 40
  )]
  condition_expression <- data.table::CJ(
    condition_id = rownames(theta),
    target_gene = paste0("G", seq_len(40))
  )
  condition_expression[, expression := 1 + seq_len(.N) / .N]
  pathway_gene_sets <- data.table::data.table(
    pathway = rep(c("Pathway A", "Pathway B"), each = 5L),
    gene = paste0("G", seq_len(10L))
  )
  gene_umap_features <- matrix(
    seq_len(40L * nrow(theta)),
    nrow = 40L,
    dimnames = list(
      paste0("G", seq_len(40L)),
      rownames(theta)
    )
  )
  link_audit <- craftgrn:::.m3_condition_document_link_audit(
    link_universe,
    assignment
  )
  optimization <- craftgrn:::.m3_condition_document_qc_optimization(
    topic_space = list(phi = phi, theta = theta),
    assignment = assignment,
    link_audit = link_audit,
    condition_gene_expression = condition_expression,
    seed = 17L
  )
  expect_identical(optimization$document_design, "condition")
  expect_equal(unname(diag(optimization$qc$raw_topic_similarity)), rep(1, 4))
  expect_true("peaks" %in% names(optimization$qc$raw_counts))
  expect_gt(sum(optimization$qc$raw_counts$peaks), 0L)
  expect_gt(sum(nzchar(
    optimization$raw_pair_assignment$gene_gammafit_topics
  )), 0L)
  output_dir <- tempfile("condition_document_qc_")
  dir.create(output_dir)
  output <- file.path(output_dir, "topic_assignment_qc_K4.pdf")
  craftgrn:::.m3_write_condition_document_topic_qc(
    topic_space = list(phi = phi, theta = theta),
    assignment = assignment,
    link_universe = link_universe,
    condition_gene_expression = condition_expression,
    out_file = output,
    method_label = "Test",
    pathway_gene_sets = pathway_gene_sets,
    gene_umap_features = gene_umap_features,
    gene_umap_coordinates = gene_umap_features[, seq_len(2L), drop = FALSE],
    gene_umap_feature_label = "condition-expression graph features",
    seed = 17L,
    overwrite = TRUE
  )
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 1000)
  qpdf <- Sys.which("qpdf")
  if (nzchar(qpdf)) {
    expect_identical(
      system2(qpdf, c("--show-npages", output), stdout = TRUE),
      "5"
    )
  }
  pdftotext <- Sys.which("pdftotext")
  if (nzchar(pdftotext)) {
    pdf_text <- system2(pdftotext, c(output, "-"), stdout = TRUE)
    expect_false(any(grepl("Condition-topic UMAP", pdf_text, fixed = TRUE)))
    expect_false(any(grepl("Pathway genes", pdf_text, fixed = TRUE)))
    expect_true(any(grepl("Test K4: Gene topics", pdf_text, fixed = TRUE)))
    expect_true(any(grepl("Test K4: Peak topics", pdf_text, fixed = TRUE)))
    expect_true(any(grepl(
      "Test K4: Topic assignment retention",
      pdf_text,
      fixed = TRUE
    )))
    expect_true(any(grepl(
      "Pearson correlation of raw condition-document theta profiles",
      pdf_text,
      fixed = TRUE
    )))
    expect_false(any(grepl(
      "Condition-topic assignment profiles",
      pdf_text,
      fixed = TRUE
    )))
  }
})

test_that("target-gene SNN uses PCA features and retains every Gene", {
  skip_if_not_installed("FNN")
  skip_if_not_installed("igraph")
  gene_profile <- matrix(
    seq_len(15L * 5L),
    nrow = 15L,
    dimnames = list(
      paste0("G", seq_len(15L)),
      paste0("C", seq_len(5L))
    )
  )
  gene_profile <- gene_profile + outer(
    seq_len(nrow(gene_profile)),
    seq_len(ncol(gene_profile)),
    function(i, j) ((i * j) %% 7) / 3
  )

  observed <- craftgrn:::.m3_gene_snn_graph(
    gene_profile,
    neighbors = 5L,
    prune = 0,
    n_pcs = 3L
  )

  expect_equal(igraph::vcount(observed$graph), 15L)
  expect_equal(dim(observed$gene_pcs), c(15L, 3L))
  expect_identical(
    rownames(observed$gene_pcs),
    paste0("GENE:G", seq_len(15L))
  )
  expect_true(all(is.finite(observed$gene_pcs)))
  expect_equal(length(observed$variance_explained), 3L)
})

test_that("target-gene cluster similarity uses expression profiles", {
  assignment <- data.table::data.table(
    term_id = paste0("GENE:", c("A", "B", "C", "D")),
    term_group = "GENE",
    assigned_topic = c(1L, 1L, 2L, 2L),
    assigned = TRUE
  )
  values <- data.table::data.table(
    condition_id = rep(c("C1", "C2"), each = 4L),
    feature_id = rep(c("A", "B", "C", "D"), 2L),
    modality = "gene",
    value = c(10, 10, 0, 0, 0, 0, 10, 10)
  )

  observed <- craftgrn:::.m3_cluster_activity_similarity(
    assignment,
    values,
    topic_ids = 1:2
  )

  expect_equal(unname(diag(observed)), c(1, 1))
  expect_equal(observed[1, 2], 0)
  expect_identical(rownames(observed), c("Topic1", "Topic2"))
})

test_that("multimodal cluster similarity combines Gene and Peak activity", {
  assignment <- data.table::data.table(
    term_id = c("GENE:A", "GENE:B", "PEAK:P1", "PEAK:P2"),
    term_group = rep(c("GENE", "PEAK"), each = 2L),
    assigned_topic = c(1L, 2L, 1L, 2L),
    assigned = TRUE
  )
  values <- data.table::data.table(
    condition_id = rep(c("C1", "C2"), each = 4L),
    feature_id = rep(c("A", "B", "P1", "P2"), 2L),
    modality = rep(c("gene", "gene", "peak", "peak"), 2L),
    value = c(10, 0, 5, 0, 0, 10, 0, 5)
  )

  observed <- craftgrn:::.m3_cluster_activity_similarity(
    assignment,
    values,
    topic_ids = 1:2
  )

  expect_equal(unname(diag(observed)), c(1, 1))
  expect_equal(observed[1, 2], 0)
})

test_that("gene-only assignments inherit links from target topics", {
  phi <- rbind(
    Topic1 = c(0.8, 0.7, 0.2, 0.1, 0.6),
    Topic2 = c(0.2, 0.3, 0.8, 0.9, 0.4)
  )
  colnames(phi) <- paste0("GENE:G", seq_len(5L))
  theta <- rbind(
    C1 = c(0.8, 0.2),
    C2 = c(0.2, 0.8)
  )
  colnames(theta) <- rownames(phi)
  assignment <- data.table::data.table(
    term_id = colnames(phi),
    term_group = "GENE",
    assigned_topic = c(1L, 1L, 2L, 2L, 1L)
  )
  links <- data.table::data.table(
    condition_id = rep(c("C1", "C2"), each = 4L),
    tf = rep(c("TF1", "TF2"), each = 4L),
    gene_key = rep(paste0("G", seq_len(4L)), 2L),
    peak_id = paste0("P", seq_len(8L)),
    fp_score = seq_len(8L)
  )
  expression <- data.table::CJ(
    condition_id = c("C1", "C2"),
    target_gene = paste0("G", seq_len(5L))
  )
  expression[, expression := seq_len(.N)]

  link_audit <- craftgrn:::.m3_condition_document_link_audit(
    links,
    assignment
  )
  observed <- craftgrn:::.m3_condition_document_qc_optimization(
    topic_space = list(phi = phi, theta = theta),
    assignment = assignment,
    link_audit = link_audit,
    condition_gene_expression = expression,
    seed = 17L
  )

  expect_true(all(link_audit$aligned))
  expect_identical(observed$assignment_mode, "gene")
  expect_equal(unname(diag(observed$qc$raw_topic_similarity)), c(1, 1))
  expect_equal(sum(observed$qc$raw_counts$peaks), 0L)
  expect_equal(sum(observed$qc$raw_counts$genes), 5L)
})
