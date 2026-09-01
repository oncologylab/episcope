.m3_condition_background_weights <- function(values,
                                             method = c(
                                               "current_softmax",
                                               "idf",
                                               "entropy",
                                               "idf_entropy"
                                             ),
                                             temperature = 0.5,
                                             floor = 0.1) {
  .assert_pkg("data.table")
  method <- match.arg(method)
  temperature <- suppressWarnings(as.numeric(temperature[[1L]]))
  floor <- suppressWarnings(as.numeric(floor[[1L]]))
  if (!is.finite(temperature) || temperature <= 0) {
    .log_abort("`temperature` must be one positive finite number.")
  }
  if (!is.finite(floor) || floor < 0 || floor > 1) {
    .log_abort("`floor` must be between 0 and 1.")
  }
  x <- data.table::as.data.table(data.table::copy(values))
  .assert_has_cols(
    x,
    c("condition_id", "feature_id", "modality", "value"),
    context = "condition-term background weighting"
  )
  x <- x[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(feature_id) & nzchar(feature_id) &
      !is.na(modality) & nzchar(modality) &
      is.finite(value) & value > 0
  ]
  if (!nrow(x)) .log_abort("No positive condition-term values are available.")
  x[, feature_key__ := paste(modality, feature_id, sep = "\r")]
  condition_count <- data.table::uniqueN(x$condition_id)
  if (condition_count < 2L) {
    .log_abort("Background weighting requires at least two conditions.")
  }

  x[, log_value__ := log2(value + 1)]
  x[, centered__ := log_value__ / temperature -
      max(log_value__ / temperature), by = feature_key__]
  x[, softmax_mass__ := exp(centered__)]
  x[!is.finite(softmax_mass__), softmax_mass__ := 0]
  x[, `:=`(
    feature_mass__ = sum(softmax_mass__),
    feature_df__ = .N
  ), by = feature_key__]
  x[, relative_activity := feature_df__ * softmax_mass__ /
      pmax(feature_mass__, .Machine$double.eps)]
  specificity <- x[, {
    probability <- softmax_mass__ / pmax(
      sum(softmax_mass__),
      .Machine$double.eps
    )
    entropy <- -sum(probability * log(pmax(probability, .Machine$double.eps))) /
      log(condition_count)
    idf <- (condition_count - .N) / (condition_count - 1)
    data.table::data.table(
      condition_frequency = .N,
      condition_frequency_fraction = .N / condition_count,
      idf_specificity = pmin(1, pmax(0, idf)),
      entropy_specificity = pmin(1, pmax(0, 1 - entropy))
    )
  }, by = feature_key__]
  specificity[, background_specificity := switch(
    method,
    current_softmax = 1,
    idf = idf_specificity,
    entropy = entropy_specificity,
    idf_entropy = 0.5 * idf_specificity + 0.5 * entropy_specificity
  )]
  x <- specificity[x, on = "feature_key__"]
  x[, background_multiplier := floor +
      (1 - floor) * background_specificity]
  if (identical(method, "current_softmax")) {
    x[, background_multiplier := 1]
  }
  x[, multiplier := relative_activity * background_multiplier]
  x[, original_modality_mass__ := sum(value),
    by = .(condition_id, modality)]
  x[, weighted__ := value * multiplier]
  x[, weighted_modality_mass__ := sum(weighted__),
    by = .(condition_id, modality)]
  x[, weight_raw := weighted__ * original_modality_mass__ /
      pmax(weighted_modality_mass__, .Machine$double.eps)]
  x[, c(
    "feature_key__", "log_value__", "centered__", "softmax_mass__",
    "feature_mass__", "feature_df__", "original_modality_mass__",
    "weighted__", "weighted_modality_mass__"
  ) := NULL]
  x[]
}

.m3_balance_condition_modalities <- function(values,
                                             target_mass = 5e5) {
  .assert_pkg("data.table")
  x <- data.table::as.data.table(data.table::copy(values))
  .assert_has_cols(
    x,
    c("condition_id", "modality", "weight_raw"),
    context = "condition modality balancing"
  )
  target_mass <- suppressWarnings(as.numeric(target_mass[[1L]]))
  if (!is.finite(target_mass) || target_mass <= 0) {
    .log_abort("`target_mass` must be one positive finite number.")
  }
  totals <- x[, .(source_mass = sum(weight_raw)),
    by = .(condition_id, modality)]
  totals[, modality_scale := target_mass /
      pmax(source_mass, .Machine$double.eps)]
  x <- totals[x, on = c("condition_id", "modality")]
  x[, weight := weight_raw * modality_scale]
  x[, c("source_mass", "modality_scale") := NULL]
  x[]
}

.m3_condition_profile_matrix <- function(values,
                                         modality,
                                         value_column = "weight") {
  .assert_pkg("data.table")
  x <- data.table::as.data.table(values)
  .assert_has_cols(
    x,
    c("condition_id", "feature_id", "modality", value_column),
    context = "condition profile matrix"
  )
  modality_name <- as.character(modality[[1L]])
  x <- x[modality == modality_name]
  if (!nrow(x)) {
    .log_abort("No {modality_name} condition profiles are available.")
  }
  condition_levels <- sort(unique(as.character(values$condition_id)))
  feature_levels <- sort(unique(as.character(x$feature_id)))
  row_index <- match(x$feature_id, feature_levels)
  column_index <- match(x$condition_id, condition_levels)
  matrix_value <- suppressWarnings(as.numeric(x[[value_column]]))
  matrix_value[!is.finite(matrix_value) | matrix_value < 0] <- 0
  out <- matrix(
    0,
    nrow = length(feature_levels),
    ncol = length(condition_levels),
    dimnames = list(feature_levels, condition_levels)
  )
  out[cbind(row_index, column_index)] <- matrix_value
  out
}

.m3_center_l2_profiles <- function(x) {
  x <- log1p(as.matrix(x))
  x <- x - rowMeans(x)
  norm <- sqrt(rowSums(x * x))
  informative <- is.finite(norm) & norm > sqrt(.Machine$double.eps)
  x[informative, , drop = FALSE] / norm[informative]
}

.m3_snn_edges <- function(profile,
                          neighbors = 30L,
                          prune = 1 / 15) {
  .assert_pkg("FNN")
  profile <- as.matrix(profile)
  neighbors <- min(
    suppressWarnings(as.integer(neighbors[[1L]])),
    nrow(profile) - 1L
  )
  if (!is.finite(neighbors) || neighbors < 2L) {
    .log_abort("At least three informative terms are required for SNN.")
  }
  nearest <- FNN::get.knn(profile, k = neighbors, algorithm = "kd_tree")
  edges <- data.table::as.data.table(
    .topic_snn_edges_cpp(nearest$nn.index)
  )
  edges <- edges[is.finite(weight) & weight > prune]
  edges[]
}

.m3_multimodal_snn_graph <- function(gene_profile,
                                     peak_profile,
                                     peak_gene_pairs,
                                     neighbors = 30L,
                                     prune = 1 / 15,
                                     cross_edge_fraction = 0.25) {
  .assert_pkg("igraph")
  .assert_pkg("data.table")
  gene_profile <- .m3_center_l2_profiles(gene_profile)
  peak_profile <- .m3_center_l2_profiles(peak_profile)
  gene_edges <- .m3_snn_edges(gene_profile, neighbors, prune)
  peak_edges <- .m3_snn_edges(peak_profile, neighbors, prune)
  gene_names <- paste0("GENE:", rownames(gene_profile))
  peak_names <- paste0("PEAK:", rownames(peak_profile))
  gene_edges[, `:=`(
    from = gene_names[from],
    to = gene_names[to],
    modality = "gene"
  )]
  peak_edges[, `:=`(
    from = peak_names[from],
    to = peak_names[to],
    modality = "peak"
  )]
  normalize_edge_mass <- function(edges, target) {
    total <- sum(edges$weight)
    if (nrow(edges) && is.finite(total) && total > 0) {
      edges[, weight := weight * target / total]
    }
    edges
  }
  gene_edges <- normalize_edge_mass(gene_edges, 0.5)
  peak_edges <- normalize_edge_mass(peak_edges, 0.5)

  pairs <- data.table::as.data.table(data.table::copy(peak_gene_pairs))
  .assert_has_cols(pairs, c("peak_id", "gene_key"), "Peak-gene graph map")
  pairs[, `:=`(
    from = paste0("PEAK:", peak_id),
    to = paste0("GENE:", gene_key)
  )]
  pairs <- unique(pairs[
    from %in% peak_names & to %in% gene_names,
    .(from, to, peak_id, gene_key)
  ])
  peak_position <- match(pairs$peak_id, rownames(peak_profile))
  gene_position <- match(pairs$gene_key, rownames(gene_profile))
  concordance <- rowSums(
    peak_profile[peak_position, , drop = FALSE] *
      gene_profile[gene_position, , drop = FALSE]
  )
  pairs[, `:=`(
    weight = 0.1 + 0.9 * pmax(0, pmin(1, concordance)),
    shared_neighbors = NA_integer_,
    modality = "cross"
  )]
  pairs <- normalize_edge_mass(pairs, cross_edge_fraction)
  edges <- data.table::rbindlist(
    list(
      gene_edges[, .(from, to, weight, shared_neighbors, modality)],
      peak_edges[, .(from, to, weight, shared_neighbors, modality)],
      pairs[, .(from, to, weight, shared_neighbors, modality)]
    ),
    use.names = TRUE,
    fill = TRUE
  )
  vertices <- data.table::data.table(
    name = c(gene_names, peak_names),
    modality = c(
      rep("gene", length(gene_names)),
      rep("peak", length(peak_names))
    )
  )
  graph <- igraph::graph_from_data_frame(
    edges[, .(from, to, weight, modality)],
    directed = FALSE,
    vertices = vertices
  )
  components <- igraph::components(graph)
  if (components$no > 1L) {
    component_size <- components$csize
    hub_component <- which.max(component_size)
    hub_vertex <- which(components$membership == hub_component)[[1L]]
    other_components <- setdiff(seq_len(components$no), hub_component)
    component_vertex <- vapply(other_components, function(component) {
      which(components$membership == component)[[1L]]
    }, integer(1L))
    positive_weights <- edges$weight[is.finite(edges$weight) & edges$weight > 0]
    connectivity_weight <- if (length(positive_weights)) {
      stats::quantile(positive_weights, 0.1, names = FALSE) * 0.05
    } else {
      1e-8
    }
    connectivity_edges <- data.table::data.table(
      from = igraph::V(graph)$name[component_vertex],
      to = igraph::V(graph)$name[hub_vertex],
      weight = connectivity_weight,
      shared_neighbors = NA_integer_,
      modality = "connectivity"
    )
    edges <- data.table::rbindlist(
      list(edges, connectivity_edges),
      use.names = TRUE,
      fill = TRUE
    )
    graph <- igraph::graph_from_data_frame(
      edges[, .(from, to, weight, modality)],
      directed = FALSE,
      vertices = vertices
    )
  }
  list(
    graph = graph,
    edges = edges,
    gene_profile = gene_profile,
    peak_profile = peak_profile,
    components_before_connectivity = components$no
  )
}

.m3_gene_snn_graph <- function(gene_profile,
                               neighbors = 30L,
                               prune = 1 / 15,
                               n_pcs = 10L) {
  .assert_pkg("igraph")
  gene_profile <- .m3_center_l2_profiles(gene_profile)
  n_pcs <- min(
    suppressWarnings(as.integer(n_pcs[[1L]])),
    ncol(gene_profile) - 1L,
    nrow(gene_profile) - 1L
  )
  if (!is.finite(n_pcs) || n_pcs < 2L) {
    .log_abort("Target-gene SNN requires at least two principal components.")
  }
  pca <- stats::prcomp(
    gene_profile,
    center = FALSE,
    scale. = FALSE,
    rank. = n_pcs
  )
  gene_pcs <- pca$x[, seq_len(n_pcs), drop = FALSE]
  edges <- .m3_snn_edges(gene_pcs, neighbors, prune)
  gene_names <- paste0("GENE:", rownames(gene_profile))
  rownames(gene_pcs) <- gene_names
  edges[, `:=`(
    from = gene_names[from],
    to = gene_names[to],
    modality = "gene"
  )]
  vertices <- data.table::data.table(
    name = gene_names,
    modality = "gene"
  )
  graph <- igraph::graph_from_data_frame(
    edges[, .(from, to, weight, modality)],
    directed = FALSE,
    vertices = vertices
  )
  components <- igraph::components(graph)
  if (components$no > 1L) {
    hub_component <- which.max(components$csize)
    hub_vertex <- which(components$membership == hub_component)[[1L]]
    other_components <- setdiff(seq_len(components$no), hub_component)
    component_vertex <- vapply(other_components, function(component) {
      which(components$membership == component)[[1L]]
    }, integer(1L))
    positive_weights <- edges$weight[is.finite(edges$weight) & edges$weight > 0]
    connectivity_weight <- if (length(positive_weights)) {
      stats::quantile(positive_weights, 0.1, names = FALSE) * 0.05
    } else {
      1e-8
    }
    connectivity_edges <- data.table::data.table(
      from = gene_names[component_vertex],
      to = gene_names[hub_vertex],
      weight = connectivity_weight,
      shared_neighbors = NA_integer_,
      modality = "connectivity"
    )
    edges <- data.table::rbindlist(
      list(edges, connectivity_edges),
      use.names = TRUE,
      fill = TRUE
    )
    graph <- igraph::graph_from_data_frame(
      edges[, .(from, to, weight, modality)],
      directed = FALSE,
      vertices = vertices
    )
  }
  variance <- pca$sdev^2
  variance <- variance / sum(variance)
  list(
    graph = graph,
    edges = edges,
    gene_profile = gene_profile,
    gene_pcs = gene_pcs,
    variance_explained = variance[seq_len(n_pcs)],
    components_before_connectivity = components$no
  )
}

.m3_leiden_resolution_path <- function(graph,
                                       target_topics,
                                       seed = 123L,
                                       resolution_grid = exp(seq(
                                         log(1e-9),
                                         log(5),
                                         length.out = 30L
                                       )),
                                       screen_iterations = 3L,
                                       final_iterations = 10L) {
  .assert_pkg("igraph")
  target_topics <- sort(unique(suppressWarnings(as.integer(target_topics))))
  target_topics <- target_topics[is.finite(target_topics) & target_topics >= 2L]
  if (!length(target_topics)) {
    .log_abort("`target_topics` must contain at least one integer >= 2.")
  }
  fits <- lapply(seq_along(resolution_grid), function(i) {
    set.seed(as.integer(seed) + i - 1L)
    fit <- igraph::cluster_leiden(
      graph,
      objective_function = "modularity",
      weights = igraph::E(graph)$weight,
      resolution = resolution_grid[[i]],
      beta = 0.01,
      n_iterations = as.integer(screen_iterations)
    )
    list(
      fit = fit,
      resolution = resolution_grid[[i]],
      communities = length(fit)
    )
  })
  summary <- data.table::rbindlist(lapply(fits, function(x) {
    data.table::data.table(
      resolution = x$resolution,
      communities = x$communities
    )
  }))
  selected <- lapply(seq_along(target_topics), function(i) {
    target <- target_topics[[i]]
    winner <- which.min(vapply(fits, function(x) {
      abs(x$communities - target) + 1e-8 * x$resolution
    }, numeric(1L)))
    resolution <- fits[[winner]]$resolution
    set.seed(as.integer(seed) + 10000L + i)
    final <- igraph::cluster_leiden(
      graph,
      objective_function = "modularity",
      weights = igraph::E(graph)$weight,
      resolution = resolution,
      beta = 0.01,
      n_iterations = as.integer(final_iterations)
    )
    list(
      membership = igraph::membership(final),
      resolution = resolution,
      communities = length(final),
      requested_topics = target,
      resolution_summary = data.table::copy(summary)[
        , distance := abs(communities - target)
      ]
    )
  })
  names(selected) <- as.character(target_topics)
  selected
}

.m3_leiden_target_membership <- function(graph,
                                         target_topics,
                                         seed = 123L,
                                         resolution_grid = exp(seq(
                                           log(1e-9),
                                           log(5),
                                           length.out = 30L
                                         )),
                                         screen_iterations = 3L,
                                         final_iterations = 10L) {
  .m3_leiden_resolution_path(
    graph = graph,
    target_topics = target_topics,
    seed = seed,
    resolution_grid = resolution_grid,
    screen_iterations = screen_iterations,
    final_iterations = final_iterations
  )[[1L]]
}

.m3_leiden_consensus_path <- function(graph,
                                      vertex_profile,
                                      target_topics,
                                      seed = 123L,
                                      base_resolution = 1,
                                      leiden_iterations = 10L) {
  .assert_pkg("igraph")
  vertex_profile <- as.matrix(vertex_profile)
  vertex_names <- igraph::V(graph)$name
  if (is.null(rownames(vertex_profile)) ||
      any(!vertex_names %in% rownames(vertex_profile))) {
    .log_abort("`vertex_profile` must contain one named row per graph vertex.")
  }
  vertex_profile <- vertex_profile[vertex_names, , drop = FALSE]
  target_topics <- sort(unique(suppressWarnings(as.integer(target_topics))))
  target_topics <- target_topics[is.finite(target_topics) & target_topics >= 2L]
  if (!length(target_topics)) {
    .log_abort("`target_topics` must contain at least one integer >= 2.")
  }
  set.seed(as.integer(seed))
  base_fit <- igraph::cluster_leiden(
    graph,
    objective_function = "modularity",
    weights = igraph::E(graph)$weight,
    resolution = as.numeric(base_resolution),
    beta = 0.01,
    n_iterations = as.integer(leiden_iterations)
  )
  base_membership <- as.integer(igraph::membership(base_fit))
  base_count <- max(base_membership)
  if (max(target_topics) > base_count) {
    .log_abort(
      "Base Leiden fit has {base_count} communities, fewer than requested K={max(target_topics)}."
    )
  }
  centroid <- rowsum(vertex_profile, base_membership, reorder = TRUE)
  community_size <- tabulate(base_membership, nbins = base_count)
  centroid <- centroid / pmax(community_size, 1L)
  merge_tree <- if (base_count > 1L) {
    stats::hclust(stats::dist(centroid), method = "ward.D2")
  } else {
    NULL
  }
  selected <- lapply(seq_along(target_topics), function(i) {
    target <- target_topics[[i]]
    community_membership <- if (target == base_count) {
      seq_len(base_count)
    } else {
      as.integer(stats::cutree(merge_tree, k = target))
    }
    membership <- community_membership[base_membership]
    names(membership) <- vertex_names
    list(
      membership = membership,
      resolution = as.numeric(base_resolution),
      communities = data.table::uniqueN(membership),
      requested_topics = target,
      base_communities = base_count,
      merge_method = "leiden_centroid_ward"
    )
  })
  names(selected) <- as.character(target_topics)
  selected
}

.m3_graph_topic_affinity <- function(graph, membership) {
  .assert_pkg("Matrix")
  .assert_pkg("igraph")
  membership <- as.integer(membership[igraph::V(graph)$name])
  topic_count <- max(membership)
  adjacency <- igraph::as_adjacency_matrix(
    graph,
    attr = "weight",
    sparse = TRUE
  )
  indicator <- Matrix::sparseMatrix(
    i = seq_along(membership),
    j = membership,
    x = 1,
    dims = c(length(membership), topic_count)
  )
  affinity <- as.matrix(adjacency %*% indicator)
  totals <- rowSums(affinity)
  totals[!is.finite(totals) | totals <= 0] <- 1
  affinity <- affinity / totals
  rownames(affinity) <- igraph::V(graph)$name
  colnames(affinity) <- paste0("Topic", seq_len(topic_count))
  affinity
}

.m3_graph_condition_activity <- function(values,
                                         membership,
                                         value_column = "weight") {
  .assert_pkg("data.table")
  x <- data.table::as.data.table(data.table::copy(values))
  .assert_has_cols(
    x,
    c("condition_id", "feature_id", "modality", value_column),
    context = "graph condition-topic activity"
  )
  membership <- data.table::data.table(
    term_id = names(membership),
    topic = as.integer(membership)
  )
  x[, term_id := data.table::fifelse(
    modality == "gene",
    paste0("GENE:", feature_id),
    paste0("PEAK:", feature_id)
  )]
  x <- membership[x, on = "term_id", nomatch = 0L]
  x[, activity := suppressWarnings(as.numeric(get(value_column)))]
  out <- x[, .(activity = sum(activity)),
    by = .(condition_id, modality, topic)]
  out[, modality_probability := activity /
      pmax(sum(activity), .Machine$double.eps),
    by = .(condition_id, modality)]
  out <- out[, .(probability = mean(modality_probability)),
    by = .(condition_id, topic)]
  out[, probability := probability /
      pmax(sum(probability), .Machine$double.eps),
    by = condition_id]
  out[]
}

.m3_cluster_activity_similarity <- function(assignment,
                                            condition_term_values,
                                            topic_ids,
                                            value_column = "value") {
  assignment <- data.table::as.data.table(data.table::copy(assignment))
  values <- data.table::as.data.table(data.table::copy(
    condition_term_values
  ))
  .assert_has_cols(
    assignment,
    c("term_id", "term_group", "assigned_topic", "assigned"),
    context = "target-gene cluster assignment"
  )
  .assert_has_cols(
    values,
    c("condition_id", "feature_id", "modality", value_column),
    context = "cluster condition activity"
  )
  membership <- unique(assignment[
    assigned %in% TRUE & is.finite(assigned_topic),
    .(
      term_id = as.character(term_id),
      topic = as.integer(assigned_topic)
    )
  ])
  values[, `:=`(
    condition_id = as.character(condition_id),
    feature_id = as.character(feature_id),
    modality = tolower(as.character(modality)),
    activity = pmax(suppressWarnings(as.numeric(get(value_column))), 0)
  )]
  values[, term_id := data.table::fifelse(
    modality == "gene",
    paste0("GENE:", feature_id),
    paste0("PEAK:", feature_id)
  )]
  profile <- membership[values, on = "term_id", nomatch = 0L][
    is.finite(activity),
    .(activity = sum(activity)),
    by = .(topic, condition_id, modality)
  ]
  profile[, modality_profile := activity /
      pmax(sum(activity), .Machine$double.eps),
    by = .(topic, modality)]
  profile <- profile[, .(
    activity = mean(modality_profile)
  ), by = .(topic, condition_id)]
  conditions <- sort(unique(values$condition_id))
  grid <- data.table::CJ(
    topic = as.integer(topic_ids),
    condition_id = conditions,
    unique = TRUE
  )
  profile <- profile[grid, on = c("topic", "condition_id")]
  profile[is.na(activity), activity := 0]
  wide <- data.table::dcast(
    profile,
    topic ~ condition_id,
    value.var = "activity",
    fill = 0
  )
  topics <- wide$topic
  wide[, topic := NULL]
  matrix <- as.matrix(wide)
  rownames(matrix) <- paste0("Topic", topics)
  matrix <- .m3_opt_row_normalize(matrix)
  similarity <- tcrossprod(sqrt(matrix))
  diag(similarity) <- 1
  similarity
}

.m3_cd_qc_heatmap <- function(x,
                              title,
                              legend_title,
                              palette = c("#F7FBFF", "#6BAED6", "#08306B"),
                              show_values = FALSE) {
  x <- as.matrix(x)
  rows <- rownames(x) %||% paste0("R", seq_len(nrow(x)))
  columns <- colnames(x) %||% paste0("C", seq_len(ncol(x)))
  long <- data.table::as.data.table(as.table(x))
  data.table::setnames(long, c("row", "column", "value"))
  long[, `:=`(
    row = factor(as.character(row), levels = rev(rows)),
    column = factor(as.character(column), levels = columns),
    value = as.numeric(value)
  )]
  upper <- stats::quantile(long$value[is.finite(long$value)], 0.98,
    na.rm = TRUE, names = FALSE)
  if (!is.finite(upper) || upper <= 0) upper <- 1
  plot <- ggplot2::ggplot(
    long,
    ggplot2::aes(column, row, fill = pmin(value, upper))
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      limits = c(0, upper),
      name = legend_title
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    .m3_qc_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "right"
    )
  if (isTRUE(show_values)) {
    plot <- plot + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", value)),
      size = 2.2,
      family = "Helvetica",
      fontface = "bold"
    )
  }
  plot
}

.m3_cd_qc_condition_page <- function(theta,
                                      condition_colors,
                                      topic_colors,
                                      seed) {
  theta <- .m3_opt_row_normalize(theta)
  coordinates <- .m3_qc_umap_coordinates(theta, seed = seed)
  condition_id <- rownames(theta)
  primary_topic <- max.col(theta, ties.method = "first")
  labels <- .m3_qc_short_condition_labels(condition_id)
  points <- data.table::data.table(
    condition_id = condition_id,
    condition_label = labels,
    topic = factor(
      paste0("T", primary_topic),
      levels = paste0("T", seq_len(ncol(theta)))
    ),
    UMAP1 = as.numeric(coordinates[, 1L]),
    UMAP2 = as.numeric(coordinates[, 2L])
  )
  condition_colors <- condition_colors[condition_id]
  missing <- is.na(condition_colors) | !nzchar(condition_colors)
  if (any(missing)) {
    condition_colors[missing] <- .module3_bright_topic_palette(
      condition_id[missing]
    )
  }
  condition_colors <- stats::setNames(condition_colors, condition_id)
  p_condition <- ggplot2::ggplot(
    points,
    ggplot2::aes(UMAP1, UMAP2, color = condition_id)
  ) +
    ggplot2::geom_point(size = 4.2, alpha = 0.82) +
    ggplot2::scale_color_manual(values = condition_colors, guide = "none") +
    ggplot2::labs(
      title = "Condition theta UMAP",
      subtitle = "Color: condition",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    .m3_qc_theme() +
    ggplot2::coord_equal()
  p_topic <- ggplot2::ggplot(
    points,
    ggplot2::aes(UMAP1, UMAP2, color = topic)
  ) +
    ggplot2::geom_point(size = 4.2, alpha = 0.82) +
    ggplot2::scale_color_manual(
      values = stats::setNames(topic_colors, paste0("T", seq_along(topic_colors))),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Condition theta UMAP",
      subtitle = "Color: maximum-probability topic",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    .m3_qc_theme() +
    ggplot2::coord_equal()
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    label_layer <- ggrepel::geom_text_repel(
      ggplot2::aes(label = condition_label),
      family = "Helvetica",
      fontface = "bold",
      size = 3,
      color = "#20272E",
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.15,
      segment.size = 0.25,
      show.legend = FALSE
    )
    p_condition <- p_condition + label_layer
    p_topic <- p_topic + label_layer
  }
  condition_similarity <- tcrossprod(sqrt(theta))
  rownames(condition_similarity) <- colnames(condition_similarity) <- labels
  correlation <- .m3_cd_qc_heatmap(
    condition_similarity,
    "Condition theta affinity",
    "Affinity",
    palette = c("#FFF7EC", "#7FCDBB", "#084081")
  )
  list(
    page = gridExtra::arrangeGrob(
      p_condition,
      p_topic,
      correlation,
      layout_matrix = rbind(c(1L, 2L), c(3L, 3L)),
      heights = c(0.52, 0.48)
    ),
    coordinates = points
  )
}

.m3_cd_qc_gene_umap_pages <- function(phi,
                                       assignment,
                                       topic_colors,
                                       seed,
                                       topics_per_page = 20L) {
  optimization <- list(raw_phi = phi, raw_topic_terms = assignment)
  coordinates <- .m3_qc_gene_term_umap_data(
    optimization,
    seed = seed
  )
  topic_ids <- seq_len(nrow(phi))
  chunks <- split(
    topic_ids,
    ceiling(seq_along(topic_ids) / as.integer(topics_per_page))
  )
  pages <- lapply(seq_along(chunks), function(page_index) {
    plots <- lapply(chunks[[page_index]], function(topic_id) {
      highlighted <- coordinates[topic_num == topic_id]
      ggplot2::ggplot(coordinates, ggplot2::aes(UMAP1, UMAP2)) +
        ggplot2::geom_point(
          color = "#D9DEE3",
          size = 0.22,
          alpha = 0.38
        ) +
        ggplot2::geom_point(
          data = highlighted,
          color = topic_colors[[topic_id]],
          size = 0.34,
          alpha = 0.86
        ) +
        ggplot2::coord_equal() +
        ggplot2::labs(title = paste0("T", topic_id), x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 9, base_family = "Helvetica") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 9,
            face = "bold",
            hjust = 0.5
          ),
          aspect.ratio = 1
        )
    })
    gridExtra::arrangeGrob(
      grobs = plots,
      ncol = 4L,
      top = grid::textGrob(
        sprintf(
          "Assigned Gene phi UMAP by topic (%d/%d)",
          page_index,
          length(chunks)
        ),
        gp = grid::gpar(
          fontfamily = "Helvetica",
          fontface = "bold",
          fontsize = 11
        )
      )
    )
  })
  list(pages = pages, coordinates = coordinates)
}

.m3_cd_qc_profile_pages <- function(theta,
                                     assignment,
                                     aligned_pairs,
                                     topics_per_page = 25L) {
  theta <- .m3_opt_row_normalize(theta)
  condition_labels <- .m3_qc_short_condition_labels(rownames(theta))
  rownames(theta) <- condition_labels
  topic_ids <- seq_len(ncol(theta))
  topic_labels <- paste0("T", topic_ids)
  colnames(theta) <- topic_labels
  counts <- assignment[, .(
    count = data.table::uniqueN(term_id)
  ), by = .(topic = as.integer(assigned_topic), term_group)]
  aligned <- data.table::as.data.table(aligned_pairs)
  if (nrow(aligned) && all(c("aligned", "topic", "peak_id", "gene_key") %in%
      names(aligned))) {
    aligned <- aligned[aligned %in% TRUE & is.finite(topic), .(
      count = data.table::uniqueN(paste(peak_id, gene_key, sep = "\r"))
    ), by = .(topic = as.integer(topic))]
    aligned[, term_group := "ALIGNED"]
    counts <- data.table::rbindlist(list(counts, aligned), fill = TRUE)
  }
  count_grid <- data.table::CJ(
    topic = topic_ids,
    term_group = c("GENE", "PEAK", "ALIGNED")
  )
  counts <- counts[count_grid, on = c("topic", "term_group")]
  counts[is.na(count), count := 0]
  chunks <- split(
    topic_ids,
    ceiling(seq_along(topic_ids) / as.integer(topics_per_page))
  )
  lapply(seq_along(chunks), function(page_index) {
    selected <- chunks[[page_index]]
    selected_labels <- topic_labels[selected]
    theta_plot <- .m3_cd_qc_heatmap(
      theta[, selected, drop = FALSE],
      "Condition-topic probability",
      "Theta",
      palette = c("#FFF7EC", "#FEC44F", "#B10026")
    )
    count_data <- counts[topic %in% selected]
    count_data[, topic_label := factor(
      paste0("T", topic),
      levels = selected_labels
    )]
    count_plot <- ggplot2::ggplot(
      count_data,
      ggplot2::aes(topic_label, count, fill = term_group)
    ) +
      ggplot2::geom_col(position = "dodge", width = 0.78) +
      ggplot2::scale_fill_manual(values = c(
        GENE = "#D97824",
        PEAK = "#007C78",
        ALIGNED = "#4E79A7"
      )) +
      ggplot2::scale_y_continuous(labels = scales::label_comma()) +
      ggplot2::labs(
        title = "Assigned terms and aligned Peak-Gene pairs",
        x = "Topic",
        y = "Unique count",
        fill = NULL
      ) +
      .m3_qc_theme() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        legend.position = "top"
      )
    gridExtra::arrangeGrob(
      theta_plot,
      count_plot,
      ncol = 1L,
      heights = c(0.68, 0.32),
      top = grid::textGrob(
        sprintf(
          "Condition-topic assignment profiles (%d/%d)",
          page_index,
          length(chunks)
        ),
        gp = grid::gpar(
          fontfamily = "Helvetica",
          fontface = "bold",
          fontsize = 11
        )
      )
    )
  })
}

.m3_cd_qc_retention_page <- function(phi, assignment, aligned_pairs) {
  all_terms <- data.table::data.table(
    term_id = colnames(phi),
    term_group = .term_group(colnames(phi))
  )
  retained <- assignment[, .(
    retained = data.table::uniqueN(term_id)
  ), by = term_group]
  totals <- all_terms[, .(total = data.table::uniqueN(term_id)), by = term_group]
  retention <- retained[totals, on = "term_group"]
  retention[is.na(retained), retained := 0]
  retention <- data.table::melt(
    retention,
    id.vars = "term_group",
    measure.vars = c("total", "retained"),
    variable.name = "stage",
    value.name = "count"
  )
  aligned <- data.table::as.data.table(aligned_pairs)
  pair_rows <- if (nrow(aligned) && all(c("aligned", "peak_id", "gene_key") %in%
      names(aligned))) {
    data.table::data.table(
      term_group = "ALIGNED PAIRS",
      stage = c("total", "retained"),
      count = c(
        data.table::uniqueN(paste(aligned$peak_id, aligned$gene_key, sep = "\r")),
        aligned[aligned %in% TRUE,
          data.table::uniqueN(paste(peak_id, gene_key, sep = "\r"))]
      )
    )
  } else {
    data.table::data.table()
  }
  retention <- data.table::rbindlist(list(retention, pair_rows), fill = TRUE)
  retention[, stage := factor(
    stage,
    levels = c("total", "retained"),
    labels = c("Input", "Assigned/aligned")
  )]
  ggplot2::ggplot(
    retention,
    ggplot2::aes(stage, count, fill = stage)
  ) +
    ggplot2::geom_col(width = 0.66) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(count)),
      vjust = -0.25,
      family = "Helvetica",
      fontface = "bold",
      size = 3
    ) +
    ggplot2::facet_wrap(~term_group, scales = "free_y", ncol = 1L) +
    ggplot2::scale_fill_manual(values = c(
      Input = "#AAB4BE",
      `Assigned/aligned` = "#007C78"
    ), guide = "none") +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.14))
    ) +
    ggplot2::labs(
      title = "Topic assignment retention",
      x = NULL,
      y = "Unique count"
    ) +
    .m3_qc_theme()
}

.m3_condition_document_link_audit <- function(link_universe, assignment) {
  .assert_pkg("data.table")
  links <- data.table::as.data.table(data.table::copy(link_universe))
  .assert_has_cols(
    links,
    c("condition_id", "tf", "gene_key", "peak_id", "fp_score"),
    context = "condition-document TF-target link universe"
  )
  terms <- data.table::as.data.table(data.table::copy(assignment))
  .assert_has_cols(
    terms,
    c("term_id", "term_group", "assigned_topic"),
    context = "condition-document term assignment"
  )
  terms[, assigned_topic := suppressWarnings(as.integer(assigned_topic))]
  gene_rows <- terms[term_group == "GENE" & is.finite(assigned_topic)]
  peak_rows <- terms[term_group == "PEAK" & is.finite(assigned_topic)]
  gene_topic <- stats::setNames(
    gene_rows$assigned_topic,
    sub("^GENE:", "", gene_rows$term_id)
  )
  gene_only <- !nrow(peak_rows)
  links[, gene_topic__ := as.integer(gene_topic[gene_key])]
  if (gene_only) {
    links[, peak_topic__ := gene_topic__]
  } else {
    peak_topic <- stats::setNames(
      peak_rows$assigned_topic,
      sub("^PEAK:", "", peak_rows$term_id)
    )
    links[, peak_topic__ := as.integer(peak_topic[peak_id])]
  }
  links[, `:=`(
    gamma_pass = is.finite(gene_topic__) & is.finite(peak_topic__),
    max_pass = is.finite(gene_topic__) & is.finite(peak_topic__),
    aligned = is.finite(gene_topic__) & is.finite(peak_topic__) &
      gene_topic__ == peak_topic__,
    topic = data.table::fifelse(
      is.finite(gene_topic__) & is.finite(peak_topic__) &
        gene_topic__ == peak_topic__,
      gene_topic__,
      NA_integer_
    )
  )]
  links[, c("gene_topic__", "peak_topic__") := NULL]
  links[]
}

.m3_condition_document_qc_similarity <- function(aligned_links,
                                                  conditions,
                                                  topics) {
  out <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, as.character(topics))
  )
  if (!nrow(aligned_links)) return(out)
  condition_counts <- aligned_links[, .(
    condition_links = data.table::uniqueN(paste(tf, gene_key, sep = "\r")),
    condition_genes = data.table::uniqueN(gene_key)
  ), by = .(condition_id, topic)]
  topic_counts <- aligned_links[, .(
    topic_links = data.table::uniqueN(paste(tf, gene_key, sep = "\r")),
    topic_genes = data.table::uniqueN(gene_key)
  ), by = topic]
  values <- topic_counts[condition_counts, on = "topic"]
  values[, mean_jaccard := 0.5 * (
    condition_links / pmax(topic_links, 1L) +
      condition_genes / pmax(topic_genes, 1L)
  )]
  rows <- match(values$condition_id, conditions)
  columns <- match(values$topic, topics)
  keep <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[keep], columns[keep])] <- values$mean_jaccard[keep]
  out
}

.m3_condition_document_qc_optimization <- function(topic_space,
                                                    assignment,
                                                    link_audit,
                                                    condition_gene_expression,
                                                    seed = 20260804L,
                                                    sample_per_condition = 120L) {
  .assert_pkg("data.table")
  phi <- .m3_opt_row_normalize(as.matrix(topic_space$phi))
  theta <- .m3_opt_row_normalize(as.matrix(topic_space$theta))
  if (is.null(rownames(theta)) || is.null(colnames(phi))) {
    .log_abort("Condition-document theta rows and phi terms must be named.")
  }
  assigned <- data.table::as.data.table(data.table::copy(assignment))
  .assert_has_cols(
    assigned,
    c("term_id", "term_group", "assigned_topic"),
    context = "condition-document term assignment"
  )
  assigned[, assigned_topic := suppressWarnings(as.integer(assigned_topic))]
  keep_assigned <- assigned$term_id %in% colnames(phi) &
    is.finite(assigned$assigned_topic)
  assigned <- unique(assigned[
    keep_assigned,
    .(term_id, term_group, assigned_topic)
  ])
  all_terms <- data.table::data.table(
    term_id = colnames(phi),
    term_group = .term_group(colnames(phi))
  )
  term_rows <- assigned[all_terms, on = c("term_id", "term_group")]
  term_rows[, `:=`(
    assigned = is.finite(assigned_topic),
    candidate_count = as.integer(is.finite(assigned_topic))
  )]
  has_peak_terms <- any(term_rows$term_group == "PEAK")

  links <- data.table::as.data.table(data.table::copy(link_audit))
  .assert_has_cols(
    links,
    c(
      "condition_id", "tf", "gene_key", "peak_id", "fp_score",
      "gamma_pass", "max_pass", "aligned", "topic"
    ),
    context = "condition-document topic QC link audit"
  )
  conditions <- rownames(theta)
  topics <- seq_len(ncol(theta))
  tf_levels <- sort(unique(as.character(links$tf)))
  target_levels <- sub(
    "^GENE:",
    "",
    term_rows[term_group == "GENE", term_id]
  )
  tf_index <- stats::setNames(seq_along(tf_levels), tf_levels)
  target_index <- stats::setNames(seq_along(target_levels), target_levels)
  condition_index <- stats::setNames(seq_along(conditions), conditions)

  data.table::setorder(links, tf, gene_key, -aligned, -fp_score, condition_id)
  global_links <- links[, {
    aligned_rows <- which(aligned %in% TRUE)
    selected <- if (length(aligned_rows)) aligned_rows[[1L]] else 1L
    list(
      condition_id = as.character(condition_id[[selected]]),
      peak_id = as.character(peak_id[[selected]]),
      fp_score = as.numeric(fp_score[[selected]]),
      gamma_pass = any(gamma_pass %in% TRUE),
      max_pass = any(max_pass %in% TRUE),
      aligned = any(aligned %in% TRUE),
      topic = if (length(aligned_rows)) as.integer(topic[[selected]]) else NA_integer_
    )
  }, by = .(tf, gene_key)]
  aligned_condition <- links[aligned %in% TRUE]
  data.table::setorder(
    aligned_condition,
    condition_id,
    tf,
    gene_key,
    -fp_score,
    peak_id
  )
  aligned_condition <- aligned_condition[, .SD[1L], by = .(
    condition_id,
    tf,
    gene_key
  )]
  condition_gene_support <- aligned_condition[, .SD[1L], by = .(
    condition_id,
    gene_key
  )]
  sampled_condition <- data.table::rbindlist(lapply(
    conditions,
    function(condition) {
      rows <- aligned_condition[condition_id == condition]
      if (nrow(rows) > sample_per_condition) {
        set.seed(as.integer(seed) + match(condition, conditions))
        rows <- rows[sort(sample(seq_len(nrow(rows)), sample_per_condition))]
      }
      rows
    }
  ), use.names = TRUE, fill = TRUE)
  global_links[, row_source := "global"]
  condition_gene_support[, row_source := "condition_gene"]
  sampled_condition[, row_source := "umap_sample"]
  qc_links <- data.table::rbindlist(
    list(global_links, condition_gene_support, sampled_condition),
    use.names = TRUE,
    fill = TRUE
  )
  qc_links[, `:=`(
    pair_index = seq_len(.N),
    doc_index = unname(condition_index[condition_id]),
    tf_index = unname(tf_index[tf]),
    target_index = unname(target_index[gene_key]),
    raw_target_topic = as.integer(topic),
    raw_posterior_topic = as.integer(topic),
    raw_posterior_probability = data.table::fifelse(aligned, 1, 0),
    raw_posterior_margin = data.table::fifelse(aligned, 1, 0),
    raw_posterior_agrees = aligned,
    raw_theta_pass = aligned,
    raw_aligned = aligned,
    optimized_topic = as.integer(topic),
    optimized_posterior_topic = as.integer(topic),
    optimized_posterior_probability = data.table::fifelse(aligned, 1, 0),
    optimized_posterior_margin = data.table::fifelse(aligned, 1, 0),
    optimized_posterior_agrees = aligned,
    optimized_theta_pass = aligned,
    optimized_aligned = aligned,
    recovered_after_merge = FALSE,
    gene_token_count = 1,
    peak_token_count = 1
  )]
  qc_assignments <- qc_links[, .(
    pair_index,
    doc_index,
    tf_index,
    target_index,
    condition_id,
    gene_token_count,
    peak_token_count,
    raw_target_topic,
    raw_posterior_topic,
    raw_posterior_probability,
    raw_posterior_margin,
    raw_posterior_agrees,
    raw_theta_pass,
    raw_aligned,
    optimized_topic,
    optimized_posterior_topic,
    optimized_posterior_probability,
    optimized_posterior_margin,
    optimized_posterior_agrees,
    optimized_theta_pass,
    optimized_aligned,
    recovered_after_merge
  )]

  condition_topic <- aligned_condition[, .(
    links = data.table::uniqueN(paste(tf, gene_key, sep = "\r")),
    genes = data.table::uniqueN(gene_key),
    tfs = data.table::uniqueN(tf)
  ), by = .(condition_id, optimized_topic = as.integer(topic))]
  link_counts <- global_links[aligned %in% TRUE, .(
    links = data.table::uniqueN(paste(tf, gene_key, sep = "\r")),
    tfs = data.table::uniqueN(tf)
  ), by = .(topic = as.integer(topic))]
  term_counts <- term_rows[
    term_group %in% c("GENE", "PEAK") & assigned %in% TRUE &
      is.finite(assigned_topic),
    .(terms = data.table::uniqueN(term_id)),
    by = .(
      topic = as.integer(assigned_topic),
      term_group
    )
  ]
  term_counts <- data.table::dcast(
    term_counts,
    topic ~ term_group,
    value.var = "terms",
    fill = 0L
  )
  if (!"GENE" %in% names(term_counts)) term_counts[, GENE := 0L]
  if (!"PEAK" %in% names(term_counts)) term_counts[, PEAK := 0L]
  data.table::setnames(term_counts, c("GENE", "PEAK"), c("genes", "peaks"))
  counts <- merge(link_counts, term_counts, by = "topic", all = TRUE)
  raw_counts <- counts[data.table::data.table(raw_topic = topics),
    on = c("topic" = "raw_topic")]
  data.table::setnames(raw_counts, "topic", "raw_topic")
  optimized_counts <- counts[data.table::data.table(optimized_topic = topics),
    on = c("topic" = "optimized_topic")]
  data.table::setnames(optimized_counts, "topic", "optimized_topic")
  for (name in c("links", "genes", "tfs", "peaks")) {
    raw_counts[is.na(get(name)), (name) := 0L]
    optimized_counts[is.na(get(name)), (name) := 0L]
  }

  sample_rows <- which(qc_links$row_source == "umap_sample")
  sampled <- qc_links[sample_rows]
  gene_columns <- match(paste0("GENE:", sampled$gene_key), colnames(phi))
  peak_columns <- match(paste0("PEAK:", sampled$peak_id), colnames(phi))
  valid_sample <- is.finite(gene_columns) & is.finite(sampled$doc_index)
  if (has_peak_terms) valid_sample <- valid_sample & is.finite(peak_columns)
  sample_rows <- sample_rows[valid_sample]
  sampled <- sampled[valid_sample]
  gene_columns <- gene_columns[valid_sample]
  peak_columns <- peak_columns[valid_sample]
  gene_profile <- t(phi[, gene_columns, drop = FALSE])
  theta_profile <- theta[sampled$doc_index, , drop = FALSE]
  sample_probability <- if (has_peak_terms) {
    peak_profile <- t(phi[, peak_columns, drop = FALSE])
    .m3_opt_row_normalize(
      sqrt(pmax(gene_profile, 0) * pmax(peak_profile, 0)) *
        sqrt(pmax(theta_profile, 0))
    )
  } else {
    .m3_opt_row_normalize(
      pmax(gene_profile, 0) * sqrt(pmax(theta_profile, 0))
    )
  }

  assigned_gene_ids <- term_rows[
    term_group == "GENE" & candidate_count > 0L,
    sub("^GENE:", "", term_id)
  ]
  raw_pair_assignment <- data.table::data.table(
    target_gene = qc_links$gene_key,
    gene_gammafit_topics = data.table::fifelse(
      qc_links$gene_key %in% assigned_gene_ids,
      "1",
      ""
    ),
    peak_gammafit_topics = data.table::fifelse(
      qc_links$gamma_pass & has_peak_terms,
      "1",
      ""
    ),
    assigned_topic = data.table::fifelse(qc_links$max_pass, 1L, NA_integer_),
    assigned = qc_links$max_pass
  )
  pair_assignment <- term_rows[term_group == "GENE", .(
    target_gene = sub("^GENE:", "", term_id),
    optimized_assigned_topic = data.table::fifelse(
      assigned,
      assigned_topic,
      NA_integer_
    )
  )]
  expression <- data.table::as.data.table(data.table::copy(
    condition_gene_expression
  ))
  .assert_has_cols(
    expression,
    c("condition_id", "target_gene", "expression"),
    context = "condition-document Gene expression"
  )
  gene_term_columns <- match(
    term_rows[term_group == "GENE", term_id],
    colnames(phi)
  )
  peak_term_columns <- match(
    term_rows[term_group == "PEAK", term_id],
    colnames(phi)
  )
  similarity <- if (has_peak_terms) {
    .m3_opt_hellinger_similarity(
      phi,
      gene_term_columns,
      peak_term_columns
    )
  } else {
    gene_phi <- .m3_opt_row_normalize(
      phi[, gene_term_columns, drop = FALSE]
    )
    tcrossprod(sqrt(gene_phi))
  }
  diag(similarity) <- 1
  cross <- .m3_condition_document_qc_similarity(
    aligned_condition,
    conditions,
    topics
  )
  identity_map <- stats::setNames(topics, topics)
  list(
    document_design = "condition",
    assignment_mode = if (has_peak_terms) "gene_peak" else "gene",
    raw_to_optimized = identity_map,
    raw_theta = theta,
    theta = theta,
    raw_phi = phi,
    phi = phi,
    raw_score = phi,
    score = phi,
    raw_topic_terms = term_rows,
    sample_rows = sample_rows,
    raw_sample_probability = sample_probability,
    optimized_sample_probability = sample_probability,
    raw_sample_coordinates = NULL,
    optimized_sample_coordinates = NULL,
    raw_pair_assignment = raw_pair_assignment,
    pair_assignment = pair_assignment,
    condition_gene_expression = expression,
    qc = list(
      assignments = qc_assignments,
      tf_levels = tf_levels,
      target_levels = target_levels,
      raw_topic_ids = topics,
      optimized_topic_ids = topics,
      raw_counts = raw_counts,
      optimized_counts = optimized_counts,
      raw_topic_similarity = similarity,
      optimized_topic_similarity = similarity,
      condition_topic = condition_topic,
      condition_topic_similarity = list(mean = cross)
    )
  )
}

.m3_write_condition_document_topic_qc <- function(topic_space,
                                                   assignment,
                                                   link_universe,
                                                   condition_gene_expression,
                                                   out_file,
                                                   method_label,
                                                   condition_colors = NULL,
                                                   pathway_gene_sets = NULL,
                                                   pathway_colors = NULL,
                                                   tf_target_gene_sets = NULL,
                                                   tf_target_panel = NULL,
                                                   gene_umap_genes = NULL,
                                                   gene_umap_features = NULL,
                                                   gene_umap_coordinates = NULL,
                                                   gene_umap_feature_label = NULL,
                                                   condition_term_values = NULL,
                                                   qc_mode = c(
                                                     "topic_model",
                                                     "gene_clustering",
                                                     "multimodal_clustering"
                                                   ),
                                                   seed = 20260804L,
                                                   overwrite = FALSE) {
  .assert_pkg("data.table")
  qc_mode <- match.arg(qc_mode)
  if (!isTRUE(overwrite) && file.exists(out_file) && file.info(out_file)$size > 0) {
    return(invisible(out_file))
  }
  if (!is.list(topic_space) || is.null(topic_space$phi) ||
      is.null(topic_space$theta)) {
    .log_abort("`topic_space` must contain phi and theta matrices.")
  }
  assignment <- data.table::as.data.table(data.table::copy(assignment))
  .assert_has_cols(
    assignment,
    c("term_id", "term_group", "assigned_topic"),
    context = "condition-document topic QC assignment"
  )
  link_audit <- .m3_condition_document_link_audit(
    link_universe,
    assignment
  )
  optimization <- .m3_condition_document_qc_optimization(
    topic_space = topic_space,
    assignment = assignment,
    link_audit = link_audit,
    condition_gene_expression = condition_gene_expression,
    seed = seed
  )
  optimization$qc_mode <- qc_mode
  if (!identical(qc_mode, "topic_model")) {
    if (is.null(condition_term_values)) {
      condition_term_values <- data.table::data.table(
        condition_id = condition_gene_expression$condition_id,
        feature_id = condition_gene_expression$target_gene,
        modality = "gene",
        value = condition_gene_expression$expression
      )
    }
    similarity <- .m3_cluster_activity_similarity(
      assignment = assignment,
      condition_term_values = condition_term_values,
      topic_ids = optimization$qc$raw_topic_ids
    )
    optimization$qc$raw_topic_similarity <- similarity
    optimization$qc$optimized_topic_similarity <- similarity
  }
  optimization$gene_umap_features <- gene_umap_features
  optimization$gene_umap_coordinates <- gene_umap_coordinates
  optimization$gene_umap_feature_label <- gene_umap_feature_label
  .write_module3_topic_assignment_qc(
    optimization = optimization,
    out_file = out_file,
    title_prefix = sprintf("%s K%d", method_label, ncol(optimization$theta)),
    condition_colors = condition_colors,
    gene_term_assignment = assignment,
    gene_umap_genes = gene_umap_genes,
    pathway_gene_sets = pathway_gene_sets,
    pathway_colors = pathway_colors,
    tf_target_gene_sets = tf_target_gene_sets,
    tf_target_panel = tf_target_panel,
    top_n_tfs = 150L,
    seed = seed,
    sections = "standard",
    report_scope = "condition_correlation"
  )
  invisible(out_file)
}

utils::globalVariables(c(
  "aligned",
  "candidate_count",
  "column",
  "condition_genes",
  "condition_links",
  "gamma_pass",
  "gene_topic__",
  "max_pass",
  "mean_jaccard",
  "peak_topic__",
  "row_source",
  "topic_genes",
  "topic_links"
))
