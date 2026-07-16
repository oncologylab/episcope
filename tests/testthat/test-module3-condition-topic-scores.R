test_that("condition theta mean reproduces direct document aggregation", {
  theta <- matrix(
    c(0.8, 0.2, 0.4, 0.6, 0.1, 0.9),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A::TF1", "A::TF2", "B::TF1"), c("Topic1", "Topic2"))
  )
  observed <- .condition_topic_theta_mean(theta)
  expect_equal(observed["A", ], c(Topic1 = 0.6, Topic2 = 0.4))
  expect_equal(observed["B", ], c(Topic1 = 0.1, Topic2 = 0.9))
})

test_that("condition target mass weights theta documents", {
  theta <- matrix(
    c(0.9, 0.1, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A::TF1", "A::TF2"), c("Topic1", "Topic2"))
  )
  mass <- data.table::data.table(
    condition_id = "A",
    doc_id = c("A::TF1", "A::TF2"),
    target_mass = c(9, 1),
    n_target_genes = c(9L, 1L)
  )
  observed <- .condition_topic_theta_weighted(theta, mass)
  expect_equal(observed["A", ], c(Topic1 = 0.82, Topic2 = 0.18), tolerance = 1e-12)
})

test_that("condition membership counts aggregate gene and peak once", {
  score <- matrix(
    c(
      0.8, 0.2, 0.6, 0.4,
      0.2, 0.8, 0.4, 0.6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Topic1", "Topic2"),
      c("GENE:G1", "GENE:G2", "PEAK:G1", "PEAK:G2")
    )
  )
  terms <- data.table::data.table(
    topic = c(1L, 2L, 1L, 2L),
    term_id = c("GENE:G1", "GENE:G2", "PEAK:G1", "PEAK:G2"),
    in_topic = TRUE
  )
  consensus <- .condition_topic_membership(score, terms, mode = "consensus")
  expect_equal(rowSums(consensus), c(G1 = 1, G2 = 1))
  expect_gt(consensus["G1", "Topic1"], consensus["G1", "Topic2"])
  expect_gt(consensus["G2", "Topic2"], consensus["G2", "Topic1"])
  projected <- .condition_topic_project(
    matrix(c(1, 0), nrow = 1, dimnames = list("A", c("G1", "G2"))),
    consensus
  )
  expect_equal(sum(projected), 1)
})

test_that("fixed phi fold-in favors the matching topic", {
  phi <- matrix(
    c(0.9, 0.1, 0.1, 0.9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Topic1", "Topic2"), c("GENE:G1", "GENE:G2"))
  )
  activity <- matrix(
    c(20, 1, 1, 20, 0, 0),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "Empty"), c("G1", "G2"))
  )
  observed <- .condition_topic_foldin(activity, phi, alpha = 0.1)
  expect_gt(observed["A", "Topic1"], observed["A", "Topic2"])
  expect_gt(observed["B", "Topic2"], observed["B", "Topic1"])
  expect_equal(observed["Empty", ], c(Topic1 = 0, Topic2 = 0))
  expect_equal(rowSums(observed)[c("A", "B")], c(A = 1, B = 1), tolerance = 1e-10)
})

test_that("positive condition contrast uses greater than or equal activity safely", {
  activity <- matrix(
    c(2, 1, 2, 3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Control", "Case"), c("G1", "G2"))
  )
  control <- .condition_topic_positive_contrast(activity, reference = "Control")
  expect_equal(control["Control", ], c(G1 = 0, G2 = 0))
  expect_equal(control["Case", ], c(G1 = 0, G2 = 2))
  expect_error(
    .condition_topic_positive_contrast(activity, reference = "Missing"),
    "reference is absent"
  )
})

test_that("condition expression distance preserves matrix dimensions", {
  expression <- matrix(
    c(1, 0, 0, 1, 1, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("G1", "G2"))
  )
  observed <- .condition_topic_expression_distance(expression)
  expect_equal(dim(observed), c(3L, 3L))
  expect_equal(diag(observed), c(A = 0, B = 0, C = 0))
  expect_equal(observed, t(observed))
})

test_that("condition-topic group metrics reward repeated-group separation", {
  separated <- rbind(
    A1 = c(0.9, 0.1), A2 = c(0.8, 0.2),
    B1 = c(0.1, 0.9), B2 = c(0.2, 0.8)
  )
  mixed <- rbind(
    A1 = c(0.9, 0.1), A2 = c(0.1, 0.9),
    B1 = c(0.8, 0.2), B2 = c(0.2, 0.8)
  )
  groups <- data.frame(
    condition_id = rownames(separated),
    condition_group = c("A", "A", "B", "B")
  )
  observed <- .condition_topic_group_metrics(
    list(separated = separated, mixed = mixed),
    groups
  )
  expect_gt(
    observed[method == "separated", condition_group_silhouette_mean],
    observed[method == "mixed", condition_group_silhouette_mean]
  )
  expect_equal(observed[method == "separated", nearest_same_group_fraction], 1)
  expect_gt(
    observed[method == "separated", between_within_js_ratio],
    observed[method == "mixed", between_within_js_ratio]
  )
})

test_that("sparse condition-topic assignment blanks unsupported topics", {
  score <- rbind(
    A = c(Topic1 = 0.6, Topic2 = 0.3, Topic3 = 0.1),
    B = c(Topic1 = 0.2, Topic2 = 0.7, Topic3 = 0.1)
  )
  relative <- .condition_topic_sparse_assignment(
    score,
    mode = "relative_max",
    relative_cutoff = 0.8
  )
  expect_equal(rowSums(relative$mask), c(A = 1, B = 1))
  expect_true(is.na(relative$display["A", "Topic2"]))
  expect_equal(relative$display["A", "Topic1"], 0.6)

  primary <- .condition_topic_sparse_assignment(score, mode = "primary_only")
  expect_equal(rowSums(primary$mask), c(A = 1, B = 1))

  flat <- matrix(
    c(0.40, 0.35, 0.25),
    nrow = 1,
    dimnames = list("C", c("Topic1", "Topic2", "Topic3"))
  )
  capped <- .condition_topic_sparse_assignment(
    flat,
    mode = "relative_max",
    relative_cutoff = 0.6,
    max_topics = 2
  )
  expect_equal(sum(capped$mask), 2)
  expect_false(capped$mask["C", "Topic3"])
})

test_that("positive-baseline assignment penalizes a shared background topic", {
  score <- rbind(
    A = c(Topic1 = 0.5, Topic2 = 0.4, Topic3 = 0.1),
    B = c(Topic1 = 0.5, Topic2 = 0.1, Topic3 = 0.4)
  )
  observed <- .condition_topic_sparse_assignment(
    score,
    mode = "positive_baseline",
    relative_cutoff = 0.8
  )
  expect_false(observed$mask["A", "Topic1"])
  expect_true(observed$mask["A", "Topic2"])
  expect_false(observed$mask["B", "Topic1"])
  expect_true(observed$mask["B", "Topic3"])
})

test_that("pair metrics report distinct and shared topic behavior", {
  methods <- list(
    separated = rbind(
      A = c(Topic1 = 0.9, Topic2 = 0.1),
      B = c(Topic1 = 0.1, Topic2 = 0.9)
    ),
    shared = rbind(
      A = c(Topic1 = 0.8, Topic2 = 0.2),
      B = c(Topic1 = 0.7, Topic2 = 0.3)
    )
  )
  observed <- .condition_topic_pair_metrics(
    methods,
    relative_cutoff = 0.7,
    max_topics = 2L
  )
  expect_true(observed[method == "separated", distinct_primary_topics])
  expect_false(observed[method == "shared", distinct_primary_topics])
  expect_gt(
    observed[method == "separated", js_distance],
    observed[method == "shared", js_distance]
  )
  expect_equal(observed$relative_cutoff, c(0.7, 0.7))
  expect_equal(observed$max_topics, c(2L, 2L))
})

test_that("pair metrics handle zero reference profiles", {
  methods <- list(
    one_zero = rbind(control = c(0, 0), treated = c(0, 1)),
    both_zero = rbind(control = c(0, 0), treated = c(0, 0))
  )
  colnames(methods$one_zero) <- colnames(methods$both_zero) <- c("Topic1", "Topic2")
  observed <- .condition_topic_pair_metrics(methods)
  expect_equal(observed[method == "one_zero", cosine_distance], 1)
  expect_equal(observed[method == "both_zero", cosine_distance], 0)
})

test_that("condition diagnostics tolerate undefined correlations", {
  methods <- list(centered = rbind(control = c(0, 0), treated = c(0, 1)))
  colnames(methods$centered) <- c("Topic1", "Topic2")
  expression <- rbind(control = c(1, 1), treated = c(1, 2))
  expect_warning(
    observed <- .condition_topic_diagnostics(methods, expression),
    NA
  )
  expect_true(is.na(observed$condition_correlation_min))
})

test_that("condition class colors cover known and unknown labels", {
  observed <- .condition_topic_group_colors(c("Ctrl", "TGFb", "Other"))
  expect_named(observed, c("Ctrl", "TGFb", "Other"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", observed)))
  expect_false(anyNA(observed))
})

test_that("condition-topic review methods have concise display labels", {
  observed <- .condition_topic_method_label(c(
    "theta_target_mass",
    "phi_foldin_one_vs_rest",
    "normtop_gamma_one_vs_rest",
    "phi_foldin_hvg",
    "normtop_gamma_hvg"
  ))
  expect_equal(observed, c(
    "Target-weighted theta",
    "Differential phi",
    "Differential specific-gene score",
    "HVG phi",
    "HVG specific-gene score"
  ))
})

test_that("HVG activity retains the most variable genes", {
  expression <- rbind(
    A = c(stable = 2, variable = 1, second = 0),
    B = c(stable = 2, variable = 9, second = 3)
  )
  observed <- .condition_topic_hvg_activity(expression, n_hvg = 2L)
  expect_equal(unname(observed[, "stable"]), c(0, 0))
  expect_equal(unname(observed[, "variable"]), unname(expression[, "variable"]))
  expect_equal(unname(observed[, "second"]), unname(expression[, "second"]))
})
