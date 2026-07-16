# Condition-to-topic score comparisons for condition-document models.

.condition_topic_row_normalize <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  keep <- is.finite(totals) & totals > 0
  if (any(keep)) x[keep, ] <- x[keep, , drop = FALSE] / totals[keep]
  x[!keep, ] <- 0
  x
}

.condition_topic_js_distance <- function(x) {
  x <- .condition_topic_row_normalize(x)
  n <- nrow(x)
  out <- matrix(0, nrow = n, ncol = n, dimnames = list(rownames(x), rownames(x)))
  if (n < 2L) return(out)
  kl <- function(a, b) {
    keep <- a > 0
    sum(a[keep] * log(a[keep] / pmax(b[keep], 1e-15)))
  }
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      mid <- (x[i, ] + x[j, ]) / 2
      out[i, j] <- out[j, i] <- sqrt((kl(x[i, ], mid) + kl(x[j, ], mid)) / 2)
    }
  }
  out
}

.condition_topic_read_matrix_csv <- function(path, row_id_col = 1L) {
  .assert_pkg("data.table")
  if (!file.exists(path)) .log_abort("Condition-topic matrix file does not exist: {path}")
  dt <- data.table::fread(path, showProgress = FALSE)
  if (ncol(dt) < 2L) .log_abort("Condition-topic matrix file has no value columns: {path}")
  row_ids <- as.character(dt[[row_id_col]])
  value_cols <- setdiff(seq_along(dt), row_id_col)
  out <- as.matrix(dt[, ..value_cols])
  storage.mode(out) <- "double"
  rownames(out) <- row_ids
  out
}

.condition_topic_manifest <- function(condition_links_dir, conditions = NULL) {
  .assert_pkg("data.table")
  manifest_path <- .module3_condition_link_manifest_path(condition_links_dir)
  manifest <- data.table::fread(manifest_path, showProgress = FALSE)
  req <- c("condition_id", "path", "format")
  if (!all(req %in% names(manifest))) {
    .log_abort("Condition-link manifest must include condition_id, path, and format columns.")
  }
  manifest[, `:=`(
    condition_id = as.character(condition_id),
    path = as.character(path),
    format = as.character(format)
  )]
  relative <- !grepl("^/", manifest$path)
  manifest[relative, path := file.path(dirname(manifest_path), path)]
  if (!is.null(conditions)) manifest <- manifest[condition_id %in% conditions]
  if (!nrow(manifest)) .log_abort("No condition-link files remain for condition-topic scoring.")
  missing <- manifest[!file.exists(path), path]
  if (length(missing)) .log_abort("Condition-link file does not exist: {missing[[1L]]}")
  data.table::setorder(manifest, condition_id, path)
  manifest[]
}

.condition_topic_collect_activity <- function(condition_links_dir,
                                              conditions = NULL,
                                              verbose = TRUE) {
  .assert_pkg("data.table")
  manifest <- .condition_topic_manifest(condition_links_dir, conditions = conditions)
  gene_parts <- vector("list", nrow(manifest))
  doc_parts <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    if (isTRUE(verbose)) {
      .log_inform(
        "Condition-topic activity: reading {i}/{nrow(manifest)} ({manifest$condition_id[[i]]})."
      )
    }
    x <- .module3_read_table(
      manifest$path[[i]],
      manifest$format[[i]],
      columns = c(
        "condition_id", "doc_id", "tf", "gene_key", "target_gene",
        "gene_expr_condition"
      ),
      allow_missing_columns = TRUE
    )
    if (!"gene_key" %in% names(x) && "target_gene" %in% names(x)) {
      x[, gene_key := as.character(target_gene)]
    }
    if (!"condition_id" %in% names(x)) x[, condition_id := manifest$condition_id[[i]]]
    if (!"doc_id" %in% names(x) && "tf" %in% names(x)) {
      x[, doc_id := paste(condition_id, tf, sep = "::")]
    }
    req <- c("condition_id", "doc_id", "tf", "gene_key", "gene_expr_condition")
    if (!all(req %in% names(x))) {
      .log_abort(
        "Condition-link input is missing columns needed for condition-topic scoring: {paste(setdiff(req, names(x)), collapse = ', ')}"
      )
    }
    x[, `:=`(
      condition_id = as.character(condition_id),
      doc_id = as.character(doc_id),
      tf = as.character(tf),
      gene_key = as.character(gene_key),
      gene_expr_condition = suppressWarnings(as.numeric(gene_expr_condition))
    )]
    x <- x[
      !is.na(condition_id) & nzchar(condition_id) &
        !is.na(doc_id) & nzchar(doc_id) &
        !is.na(tf) & nzchar(tf) &
        !is.na(gene_key) & nzchar(gene_key) &
        is.finite(gene_expr_condition) & gene_expr_condition > 0
    ]
    if (!nrow(x)) next
    gene_expr <- x[, .(expression = max(gene_expr_condition)), by = .(condition_id, gene_key)]
    gene_tf <- unique(x[, .(condition_id, gene_key, tf)])[
      , .(n_regulating_tfs = .N), by = .(condition_id, gene_key)
    ]
    gene_parts[[i]] <- merge(
      gene_expr,
      gene_tf,
      by = c("condition_id", "gene_key"),
      all = TRUE,
      sort = FALSE
    )
    doc_gene <- unique(x[, .(condition_id, doc_id, gene_key, gene_expr_condition)])
    doc_parts[[i]] <- doc_gene[, .(
      target_mass = sum(log1p(gene_expr_condition)),
      n_target_genes = .N
    ), by = .(condition_id, doc_id)]
    rm(x, gene_expr, gene_tf, doc_gene)
    invisible(gc(FALSE))
  }
  gene_activity <- data.table::rbindlist(gene_parts, use.names = TRUE, fill = TRUE)
  doc_mass <- data.table::rbindlist(doc_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(gene_activity)) .log_abort("No positive condition-gene activity was found.")
  gene_activity <- gene_activity[, .(
    expression = max(expression, na.rm = TRUE),
    n_regulating_tfs = max(n_regulating_tfs, na.rm = TRUE)
  ), by = .(condition_id, gene_key)]
  doc_mass <- doc_mass[, .(
    target_mass = max(target_mass, na.rm = TRUE),
    n_target_genes = max(n_target_genes, na.rm = TRUE)
  ), by = .(condition_id, doc_id)]
  list(gene_activity = gene_activity[], doc_mass = doc_mass[], manifest = manifest)
}

.condition_topic_activity_matrices <- function(gene_activity, conditions, genes) {
  expression <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(genes),
    dimnames = list(conditions, genes)
  )
  tf_support <- expression
  rows <- match(gene_activity$condition_id, conditions)
  cols <- match(gene_activity$gene_key, genes)
  keep <- !is.na(rows) & !is.na(cols)
  if (any(keep)) {
    expression[cbind(rows[keep], cols[keep])] <- gene_activity$expression[keep]
    tf_support[cbind(rows[keep], cols[keep])] <- gene_activity$n_regulating_tfs[keep]
  }
  list(expression = log1p(expression), tf_support = tf_support)
}

.condition_topic_theta_mean <- function(theta) {
  conditions <- sub("::[^:]+$", "", rownames(theta))
  out <- do.call(rbind, lapply(unique(conditions), function(condition) {
    colMeans(theta[conditions == condition, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(out) <- unique(conditions)
  colnames(out) <- colnames(theta)
  .condition_topic_row_normalize(out)
}

.condition_topic_method_label <- function(method) {
  method <- as.character(method)
  labels <- c(
    theta_mean = "Original mean theta",
    theta_target_mass = "Target-weighted theta",
    phi_foldin_one_vs_rest = "Differential phi",
    normtop_gamma_one_vs_rest = "Differential specific-gene score",
    phi_foldin_hvg = "HVG phi",
    normtop_gamma_hvg = "HVG specific-gene score"
  )
  out <- unname(labels[method])
  missing <- is.na(out)
  out[missing] <- gsub("_", " ", method[missing], fixed = TRUE)
  out
}

.condition_topic_hvg_activity <- function(expression, n_hvg = 500L) {
  expression <- as.matrix(expression)
  n_hvg <- as.integer(n_hvg)
  if (!is.finite(n_hvg) || n_hvg < 1L) {
    .log_abort("n_hvg must be a positive integer.")
  }
  variability <- apply(expression, 2L, stats::var, na.rm = TRUE)
  variability[!is.finite(variability)] <- 0
  gene_names <- colnames(expression) %||% sprintf("gene_%d", seq_len(ncol(expression)))
  selected <- order(-variability, gene_names)[seq_len(min(n_hvg, ncol(expression)))]
  out <- matrix(0, nrow(expression), ncol(expression), dimnames = dimnames(expression))
  out[, selected] <- expression[, selected, drop = FALSE]
  out
}

.condition_topic_theta_weighted <- function(theta, doc_mass) {
  doc_ids <- rownames(theta)
  conditions <- sub("::[^:]+$", "", doc_ids)
  weights <- doc_mass$target_mass[match(doc_ids, doc_mass$doc_id)]
  weights[!is.finite(weights) | weights <= 0] <- 0
  out <- do.call(rbind, lapply(unique(conditions), function(condition) {
    idx <- which(conditions == condition)
    w <- weights[idx]
    if (!sum(w) > 0) return(colMeans(theta[idx, , drop = FALSE], na.rm = TRUE))
    colSums(theta[idx, , drop = FALSE] * w, na.rm = TRUE) / sum(w)
  }))
  rownames(out) <- unique(conditions)
  colnames(out) <- colnames(theta)
  .condition_topic_row_normalize(out)
}

.condition_topic_membership <- function(score_mat,
                                        topic_terms = NULL,
                                        mode = c("all", "gamma", "consensus")) {
  mode <- match.arg(mode)
  score_mat <- as.matrix(score_mat)
  topics <- rownames(score_mat)
  gene_terms <- colnames(score_mat)[startsWith(colnames(score_mat), "GENE:")]
  genes <- sub("^GENE:", "", gene_terms)
  gene_score <- t(score_mat[, gene_terms, drop = FALSE])
  rownames(gene_score) <- genes
  if (identical(mode, "all")) {
    return(.condition_topic_row_normalize(gene_score))
  }
  if (is.null(topic_terms) || !is.data.frame(topic_terms)) {
    .log_abort("Gamma and consensus condition-topic scores require topic_terms.")
  }
  terms <- data.table::as.data.table(topic_terms)
  if (!all(c("topic", "term_id", "in_topic") %in% names(terms))) {
    .log_abort("topic_terms must include topic, term_id, and in_topic columns.")
  }
  terms[, in_topic := .as_logical_flag(in_topic)]
  gene_mask <- matrix(0, nrow = length(genes), ncol = length(topics), dimnames = list(genes, topics))
  gene_pass <- terms[in_topic & startsWith(term_id, "GENE:")]
  gene_rows <- match(sub("^GENE:", "", gene_pass$term_id), genes)
  gene_cols <- match(paste0("Topic", as.integer(gene_pass$topic)), topics)
  keep <- !is.na(gene_rows) & !is.na(gene_cols)
  gene_mask[cbind(gene_rows[keep], gene_cols[keep])] <- 1
  gene_membership <- .condition_topic_row_normalize(gene_score * gene_mask)
  if (identical(mode, "gamma")) return(gene_membership)

  peak_terms <- paste0("PEAK:", genes)
  peak_idx <- match(peak_terms, colnames(score_mat))
  peak_score <- matrix(0, nrow = length(genes), ncol = length(topics), dimnames = list(genes, topics))
  has_peak <- !is.na(peak_idx)
  if (any(has_peak)) peak_score[has_peak, ] <- t(score_mat[, peak_idx[has_peak], drop = FALSE])
  peak_mask <- matrix(0, nrow = length(genes), ncol = length(topics), dimnames = list(genes, topics))
  peak_pass <- terms[in_topic & startsWith(term_id, "PEAK:")]
  peak_rows <- match(sub("^PEAK:", "", peak_pass$term_id), genes)
  peak_cols <- match(paste0("Topic", as.integer(peak_pass$topic)), topics)
  keep <- !is.na(peak_rows) & !is.na(peak_cols)
  peak_mask[cbind(peak_rows[keep], peak_cols[keep])] <- 1
  peak_membership <- .condition_topic_row_normalize(peak_score * peak_mask)
  available <- (rowSums(gene_membership) > 0) + (rowSums(peak_membership) > 0)
  consensus <- gene_membership + peak_membership
  keep <- available > 0
  consensus[keep, ] <- consensus[keep, , drop = FALSE] / available[keep]
  .condition_topic_row_normalize(consensus)
}

.condition_topic_project <- function(activity, membership) {
  genes <- intersect(colnames(activity), rownames(membership))
  if (!length(genes)) {
    return(matrix(0, nrow(activity), ncol(membership), dimnames = list(rownames(activity), colnames(membership))))
  }
  .condition_topic_row_normalize(activity[, genes, drop = FALSE] %*% membership[genes, , drop = FALSE])
}

.condition_topic_positive_contrast <- function(activity, reference = NULL) {
  activity <- as.matrix(activity)
  if (is.null(reference)) {
    out <- sweep(activity, 2L, colMeans(activity, na.rm = TRUE), "-")
  } else {
    if (!reference %in% rownames(activity)) .log_abort("Condition-topic reference is absent: {reference}")
    out <- sweep(activity, 2L, activity[reference, ], "-")
  }
  out[!is.finite(out) | out < 0] <- 0
  out
}

.condition_topic_foldin <- function(activity,
                                    phi,
                                    alpha = NULL,
                                    max_iter = 200L,
                                    tolerance = 1e-8) {
  activity <- as.matrix(activity)
  phi <- as.matrix(phi)
  genes <- intersect(colnames(activity), sub("^GENE:", "", colnames(phi)[startsWith(colnames(phi), "GENE:")]))
  topics <- rownames(phi)
  if (!length(genes)) {
    return(matrix(0, nrow(activity), nrow(phi), dimnames = list(rownames(activity), topics)))
  }
  phi <- phi[, paste0("GENE:", genes), drop = FALSE]
  phi[!is.finite(phi) | phi < 0] <- 0
  phi <- .condition_topic_row_normalize(phi)
  activity <- activity[, genes, drop = FALSE]
  k <- nrow(phi)
  alpha <- suppressWarnings(as.numeric(alpha))
  if (!length(alpha) || !all(is.finite(alpha)) || any(alpha <= 0)) alpha <- 50 / k
  if (length(alpha) == 1L) alpha <- rep(alpha, k)
  if (length(alpha) != k) .log_abort("Fold-in alpha must be scalar or have one value per topic.")
  out <- matrix(0, nrow(activity), k, dimnames = list(rownames(activity), topics))
  for (i in seq_len(nrow(activity))) {
    counts <- activity[i, ]
    keep <- counts > 0 & colSums(phi) > 0
    if (!any(keep)) next
    words <- phi[, keep, drop = FALSE]
    counts <- counts[keep]
    theta <- alpha / sum(alpha)
    for (iteration in seq_len(as.integer(max_iter))) {
      denom <- as.numeric(crossprod(theta, words))
      responsibility <- sweep(words, 2L, pmax(denom, 1e-300), "/") * theta
      topic_count <- rowSums(sweep(responsibility, 2L, counts, "*"))
      updated <- (topic_count + alpha) / (sum(topic_count) + sum(alpha))
      if (max(abs(updated - theta)) < tolerance) {
        theta <- updated
        break
      }
      theta <- updated
    }
    out[i, ] <- theta
  }
  .condition_topic_row_normalize(out)
}

.condition_topic_expression_distance <- function(expression) {
  x <- as.matrix(expression)
  norms <- sqrt(rowSums(x * x))
  norms[!is.finite(norms) | norms <= 0] <- 1
  similarity <- tcrossprod(x / norms)
  similarity[!is.finite(similarity)] <- 0
  similarity[similarity < -1] <- -1
  similarity[similarity > 1] <- 1
  distance <- 1 - similarity
  diag(distance) <- 0
  distance
}

.condition_topic_diagnostics <- function(methods, expression, contributing = NULL) {
  .assert_pkg("data.table")
  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x)) fun(x) else NA_real_
  }
  expression_distance <- .condition_topic_expression_distance(expression)
  expression_values <- expression_distance[upper.tri(expression_distance)]
  out <- lapply(names(methods), function(method) {
    score <- .condition_topic_row_normalize(methods[[method]])
    correlation <- suppressWarnings(stats::cor(t(score), use = "pairwise.complete.obs"))
    correlation_values <- correlation[upper.tri(correlation)]
    js <- .condition_topic_js_distance(score)
    js_values <- js[upper.tri(js)]
    entropy <- -rowSums(score * log(pmax(score, 1e-15))) / log(ncol(score))
    concordance <- suppressWarnings(stats::cor(js_values, expression_values, method = "spearman", use = "complete.obs"))
    genes <- if (!is.null(contributing) && method %in% names(contributing)) contributing[[method]] else rep(NA_integer_, nrow(score))
    data.table::data.table(
      method = method,
      n_conditions = nrow(score),
      n_topics = ncol(score),
      condition_correlation_median = finite_summary(correlation_values, stats::median),
      condition_correlation_min = finite_summary(correlation_values, min),
      js_distance_median = finite_summary(js_values, stats::median),
      js_distance_max = finite_summary(js_values, max),
      normalized_entropy_median = finite_summary(entropy, stats::median),
      primary_score_mean = mean(apply(score, 1L, max), na.rm = TRUE),
      primary_score_min = min(apply(score, 1L, max), na.rm = TRUE),
      primary_score_max = max(apply(score, 1L, max), na.rm = TRUE),
      expression_distance_spearman = concordance,
      contributing_genes_mean = mean(genes, na.rm = TRUE),
      contributing_genes_min = suppressWarnings(min(genes, na.rm = TRUE)),
      contributing_genes_max = suppressWarnings(max(genes, na.rm = TRUE))
    )
  })
  ans <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  numeric_cols <- names(ans)[vapply(ans, is.numeric, logical(1))]
  for (col in numeric_cols) ans[!is.finite(get(col)), (col) := NA_real_]
  ans[]
}

.condition_topic_group_metrics <- function(methods, condition_groups) {
  .assert_pkg("data.table")
  groups <- data.table::as.data.table(condition_groups)
  if (!all(c("condition_id", "condition_group") %in% names(groups))) {
    .log_abort("condition_groups must include condition_id and condition_group columns.")
  }
  groups <- unique(groups[, .(
    condition_id = as.character(condition_id),
    condition_group = as.character(condition_group)
  )])
  if (groups[, anyDuplicated(condition_id)]) {
    .log_abort("condition_groups must contain one row per condition_id.")
  }
  data.table::rbindlist(lapply(names(methods), function(method) {
    score <- .condition_topic_row_normalize(methods[[method]])
    group <- groups$condition_group[match(rownames(score), groups$condition_id)]
    keep <- !is.na(group) & nzchar(group)
    score <- score[keep, , drop = FALSE]
    group <- group[keep]
    distance <- .condition_topic_js_distance(score)
    pair_index <- which(upper.tri(distance), arr.ind = TRUE)
    pair_distance <- distance[pair_index]
    same_group <- group[pair_index[, 1L]] == group[pair_index[, 2L]]
    within <- pair_distance[same_group]
    between <- pair_distance[!same_group]

    group_size <- table(group)
    eligible <- names(group_size)[group_size >= 2L]
    silhouette <- rep(NA_real_, length(group))
    nearest_same <- rep(NA, length(group))
    for (i in seq_along(group)) {
      if (!group[[i]] %in% eligible) next
      same <- which(group == group[[i]] & seq_along(group) != i)
      other_groups <- setdiff(unique(group), group[[i]])
      a <- mean(distance[i, same], na.rm = TRUE)
      b <- min(vapply(other_groups, function(other_group) {
        mean(distance[i, group == other_group], na.rm = TRUE)
      }, numeric(1)), na.rm = TRUE)
      denominator <- max(a, b)
      silhouette[[i]] <- if (is.finite(denominator) && denominator > 0) {
        (b - a) / denominator
      } else {
        0
      }
      nearest <- order(distance[i, ], decreasing = FALSE)
      nearest <- nearest[nearest != i][1L]
      nearest_same[[i]] <- length(nearest) == 1L && group[[nearest]] == group[[i]]
    }

    within_mean <- mean(within, na.rm = TRUE)
    between_mean <- mean(between, na.rm = TRUE)
    data.table::data.table(
      method = method,
      n_grouped_conditions = nrow(score),
      n_condition_groups = data.table::uniqueN(group),
      n_repeated_condition_groups = length(eligible),
      n_within_group_pairs = length(within),
      n_between_group_pairs = length(between),
      within_group_js_mean = within_mean,
      within_group_js_median = stats::median(within, na.rm = TRUE),
      between_group_js_mean = between_mean,
      between_group_js_median = stats::median(between, na.rm = TRUE),
      between_within_js_ratio = between_mean / within_mean,
      between_within_js_margin = between_mean - within_mean,
      condition_group_silhouette_mean = mean(silhouette, na.rm = TRUE),
      nearest_same_group_fraction = mean(nearest_same, na.rm = TRUE)
    )
  }), use.names = TRUE, fill = TRUE)
}

.condition_topic_sparse_assignment <- function(score,
                                               mode = c("relative_max", "primary_only", "positive_baseline"),
                                               relative_cutoff = 0.8,
                                               max_topics = Inf) {
  mode <- match.arg(mode)
  score <- .condition_topic_row_normalize(score)
  relative_cutoff <- suppressWarnings(as.numeric(relative_cutoff)[1L])
  if (!is.finite(relative_cutoff) || relative_cutoff < 0 || relative_cutoff > 1) {
    .log_abort("relative_cutoff must be between 0 and 1.")
  }
  max_topics <- suppressWarnings(as.integer(max_topics)[1L])
  if (!is.finite(max_topics) || max_topics < 1L) max_topics <- ncol(score)
  evidence <- score
  if (identical(mode, "positive_baseline")) {
    evidence <- sweep(score, 2L, colMeans(score, na.rm = TRUE), "-")
    evidence[!is.finite(evidence) | evidence < 0] <- 0
  }
  mask <- matrix(FALSE, nrow(score), ncol(score), dimnames = dimnames(score))
  for (i in seq_len(nrow(score))) {
    primary <- which.max(score[i, ])
    if (identical(mode, "primary_only")) {
      selected <- primary
    } else {
      row_max <- max(evidence[i, ], na.rm = TRUE)
      selected <- if (is.finite(row_max) && row_max > 0) {
        which(evidence[i, ] >= relative_cutoff * row_max)
      } else {
        primary
      }
    }
    if (length(selected) > max_topics) {
      selected <- selected[order(evidence[i, selected], decreasing = TRUE)][seq_len(max_topics)]
    }
    if (!length(selected)) selected <- primary
    mask[i, selected] <- TRUE
  }
  assigned <- score
  assigned[!mask] <- 0
  display <- score
  display[!mask] <- NA_real_
  list(
    assigned = .condition_topic_row_normalize(assigned),
    display = display,
    mask = mask,
    evidence = evidence,
    retained_mass = rowSums(score * mask)
  )
}

.condition_topic_sparse_metrics <- function(assignments, condition_groups) {
  .assert_pkg("data.table")
  groups <- data.table::as.data.table(condition_groups)
  data.table::rbindlist(lapply(names(assignments), function(name) {
    assignment <- assignments[[name]]
    mask <- assignment$mask
    group <- groups$condition_group[match(rownames(mask), groups$condition_id)]
    keep <- !is.na(group) & nzchar(group)
    mask <- mask[keep, , drop = FALSE]
    group <- group[keep]
    pair_index <- which(upper.tri(matrix(FALSE, nrow(mask), nrow(mask))), arr.ind = TRUE)
    jaccard <- apply(pair_index, 1L, function(pair) {
      left <- mask[pair[[1L]], ]
      right <- mask[pair[[2L]], ]
      union_size <- sum(left | right)
      if (union_size > 0) sum(left & right) / union_size else 0
    })
    same_group <- group[pair_index[, 1L]] == group[pair_index[, 2L]]
    repeated <- names(table(group))[table(group) >= 2L]
    evaluable_pair <- same_group & group[pair_index[, 1L]] %in% repeated
    between <- !same_group
    within_jaccard <- mean(jaccard[evaluable_pair], na.rm = TRUE)
    between_jaccard <- mean(jaccard[between], na.rm = TRUE)
    data.table::data.table(
      sparse_method = name,
      topics_per_condition_mean = mean(rowSums(mask)),
      topics_per_condition_median = stats::median(rowSums(mask)),
      topics_per_condition_max = max(rowSums(mask)),
      retained_probability_mass_mean = mean(assignment$retained_mass[keep]),
      selected_topics_total = sum(colSums(mask) > 0),
      within_family_topic_jaccard = within_jaccard,
      between_family_topic_jaccard = between_jaccard,
      within_between_jaccard_margin = within_jaccard - between_jaccard
    )
  }), use.names = TRUE, fill = TRUE)
}

.condition_topic_pair_metrics <- function(methods,
                                          relative_cutoff = 0.7,
                                          max_topics = 2L) {
  .assert_pkg("data.table")
  data.table::rbindlist(lapply(names(methods), function(method) {
    score <- .condition_topic_row_normalize(methods[[method]])
    if (nrow(score) != 2L) {
      .log_abort("Pairwise condition-topic metrics require exactly two conditions.")
    }
    sparse <- .condition_topic_sparse_assignment(
      score,
      mode = "relative_max",
      relative_cutoff = relative_cutoff,
      max_topics = max_topics
    )
    primary <- colnames(score)[max.col(score, ties.method = "first")]
    selected_left <- sparse$mask[1L, ]
    selected_right <- sparse$mask[2L, ]
    selected_union <- sum(selected_left | selected_right)
    sparse_jaccard <- if (selected_union > 0L) {
      sum(selected_left & selected_right) / selected_union
    } else {
      0
    }
    left_norm <- sqrt(sum(score[1L, ]^2))
    right_norm <- sqrt(sum(score[2L, ]^2))
    cosine_distance <- if (left_norm > 0 && right_norm > 0) {
      1 - sum(score[1L, ] * score[2L, ]) / (left_norm * right_norm)
    } else if (left_norm > 0 || right_norm > 0) {
      1
    } else {
      0
    }
    entropy <- apply(score, 1L, function(x) {
      positive <- x[x > 0]
      if (!length(positive)) return(0)
      -sum(positive * log(positive)) / log(ncol(score))
    })
    data.table::data.table(
      method = method,
      condition_1 = rownames(score)[[1L]],
      condition_2 = rownames(score)[[2L]],
      js_distance = .condition_topic_js_distance(score)[1L, 2L],
      cosine_distance = cosine_distance,
      primary_topic_1 = primary[[1L]],
      primary_topic_2 = primary[[2L]],
      distinct_primary_topics = primary[[1L]] != primary[[2L]],
      max_topic_difference = max(abs(score[1L, ] - score[2L, ])),
      mean_topic_difference = mean(abs(score[1L, ] - score[2L, ])),
      normalized_entropy_mean = mean(entropy),
      sparse_topic_jaccard = sparse_jaccard,
      sparse_topics_condition_1 = sum(selected_left),
      sparse_topics_condition_2 = sum(selected_right),
      sparse_retained_mass_mean = mean(sparse$retained_mass),
      relative_cutoff = as.numeric(relative_cutoff),
      max_topics = as.integer(max_topics)
    )
  }), use.names = TRUE, fill = TRUE)
}

.condition_topic_group_colors <- function(groups) {
  groups <- unique(as.character(groups))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  fixed <- c(
    FBS = "#4477AA", BCAA = "#EE6677", `Gln.Arg` = "#228833",
    Gln = "#CCBB44", `Met.Cys` = "#66CCEE", Glc = "#AA3377",
    Lys = "#BBBBBB", Trp = "#EE7733", Arg = "#333333",
    Ctrl = "#111111", TGFb = "#CC6677", BATF_FloxedOut = "#AA3377",
    WT_FloxedOut = "#4477AA"
  )
  out <- fixed[groups]
  missing <- is.na(out)
  if (any(missing)) {
    palette <- grDevices::hcl.colors(sum(missing), palette = "Dark 3")
    out[missing] <- palette
  }
  stats::setNames(unname(out), groups)
}

.plot_sparse_condition_topic_assignments <- function(assignments,
                                                     out_file,
                                                     condition_groups,
                                                     title_prefix = NULL,
                                                     drop_empty_topics = FALSE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_abort("Package pheatmap is required for sparse condition-topic heatmaps.")
  }
  groups <- data.table::as.data.table(condition_groups)
  finite_values <- unlist(lapply(assignments, function(x) x$display[is.finite(x$display)]))
  color_max <- max(finite_values, na.rm = TRUE)
  if (!is.finite(color_max) || color_max <= 0) color_max <- 1
  colors <- grDevices::colorRampPalette(c("#FFFFFF", "#E8F2F8", "#BFDDEA", "#FFF7BC", "#FCAE61", "#D95F52"))(101)
  breaks <- seq(0, color_max, length.out = 102L)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(out_file, width = 12, height = 9.5, family = "Helvetica", onefile = TRUE, bg = "white")
  } else {
    grDevices::pdf(
      out_file, width = 12, height = 9.5, family = "Helvetica",
      onefile = TRUE, useDingbats = FALSE, bg = "white"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  for (name in names(assignments)) {
    assignment <- assignments[[name]]
    display <- t(assignment$display)
    if (isTRUE(drop_empty_topics)) {
      keep_topics <- colSums(assignment$mask) > 0
      display <- display[keep_topics, , drop = FALSE]
    }
    plot_matrix <- display
    plot_matrix[!is.finite(plot_matrix)] <- 0
    assigned <- assignment$assigned
    condition_tree <- if (nrow(assigned) > 1L) {
      stats::hclust(stats::as.dist(.condition_topic_js_distance(assigned)))
    } else {
      FALSE
    }
    topic_tree <- if (nrow(display) > 1L) {
      stats::hclust(stats::dist(t(assigned[, rownames(display), drop = FALSE])))
    } else {
      FALSE
    }
    annotation_col <- data.frame(
      `Condition class` = groups$condition_group[match(colnames(display), groups$condition_id)],
      check.names = FALSE
    )
    rownames(annotation_col) <- colnames(display)
    present <- unique(annotation_col[[1L]])
    annotation_colors <- list(`Condition class` = .condition_topic_group_colors(present))
    labels <- matrix("", nrow(display), ncol(display), dimnames = dimnames(display))
    labels[is.finite(display)] <- sprintf("%.2f", display[is.finite(display)])
    ph <- pheatmap::pheatmap(
      plot_matrix,
      color = colors,
      breaks = breaks,
      cluster_rows = topic_tree,
      cluster_cols = condition_tree,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 9,
      fontsize_row = 9,
      fontsize_col = 9,
      display_numbers = labels,
      number_color = "#111111",
      fontsize_number = 9,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      angle_col = 90,
      main = paste(c(title_prefix, .condition_topic_method_label(name)), collapse = " | "),
      border_color = "#FFFFFF",
      silent = TRUE
    )
    ph$gtable <- .condition_topic_bold_grob(ph$gtable)
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.draw(ph$gtable)
  }
  invisible(out_file)
}

.condition_topic_long_scores <- function(methods, method_info = NULL, contributing = NULL) {
  .assert_pkg("data.table")
  out <- data.table::rbindlist(lapply(names(methods), function(method) {
    score <- .condition_topic_row_normalize(methods[[method]])
    dt <- data.table::as.data.table(as.table(score))
    data.table::setnames(dt, c("condition_id", "topic", "score"))
    dt[, `:=`(
      method = method,
      condition_id = as.character(condition_id),
      topic = as.character(topic),
      topic_num = as.integer(sub("^Topic", "", as.character(topic))),
      score = as.numeric(score)
    )]
    dt[, topic_rank := data.table::frank(-score, ties.method = "first"), by = condition_id]
    dt[, `:=`(
      primary_topic = topic[which.max(score)],
      primary_score = max(score)
    ), by = condition_id]
    if (!is.null(contributing) && method %in% names(contributing)) {
      gene_count <- data.table::data.table(
        condition_id = rownames(score),
        contributing_genes = as.integer(contributing[[method]])
      )
      dt <- merge(dt, gene_count, by = "condition_id", all.x = TRUE, sort = FALSE)
    } else {
      dt[, contributing_genes := NA_integer_]
    }
    dt
  }), use.names = TRUE, fill = TRUE)
  if (!is.null(method_info) && nrow(method_info)) {
    out <- merge(out, method_info, by = "method", all.x = TRUE, sort = FALSE)
  }
  data.table::setorder(out, method, condition_id, topic_rank, topic_num)
  out[]
}

.condition_topic_build_methods <- function(theta,
                                           doc_mass,
                                           expression,
                                           tf_support,
                                           phi,
                                           alpha,
                                           memberships,
                                           reference_condition = NULL,
                                           n_hvg = 500L,
                                           foldin_max_iter = 200L,
                                           foldin_tolerance = 1e-8) {
  one_vs_rest <- .condition_topic_positive_contrast(expression)
  hvg_activity <- .condition_topic_hvg_activity(expression, n_hvg = n_hvg)
  control <- if (!is.null(reference_condition)) {
    .condition_topic_positive_contrast(expression, reference = reference_condition)
  } else {
    NULL
  }
  methods <- list(
    theta_mean = .condition_topic_theta_mean(theta),
    theta_target_mass = .condition_topic_theta_weighted(theta, doc_mass),
    phi_foldin_absolute = .condition_topic_foldin(
      expression, phi, alpha = alpha, max_iter = foldin_max_iter, tolerance = foldin_tolerance
    ),
    phi_foldin_one_vs_rest = .condition_topic_foldin(
      one_vs_rest, phi, alpha = alpha, max_iter = foldin_max_iter, tolerance = foldin_tolerance
    ),
    phi_foldin_hvg = .condition_topic_foldin(
      hvg_activity, phi, alpha = alpha, max_iter = foldin_max_iter, tolerance = foldin_tolerance
    ),
    normtop_all_absolute = .condition_topic_project(expression, memberships$all),
    normtop_all_one_vs_rest = .condition_topic_project(one_vs_rest, memberships$all),
    normtop_gamma_absolute = .condition_topic_project(expression, memberships$gamma),
    normtop_gamma_one_vs_rest = .condition_topic_project(one_vs_rest, memberships$gamma),
    normtop_gamma_hvg = .condition_topic_project(hvg_activity, memberships$gamma),
    normtop_consensus_absolute = .condition_topic_project(expression, memberships$consensus),
    normtop_consensus_one_vs_rest = .condition_topic_project(one_vs_rest, memberships$consensus),
    normtop_tf_support_one_vs_rest = .condition_topic_project(
      one_vs_rest * log1p(tf_support),
      memberships$gamma
    )
  )
  if (!is.null(control)) {
    methods$normtop_all_control <- .condition_topic_project(control, memberships$all)
    methods$normtop_gamma_control <- .condition_topic_project(control, memberships$gamma)
    methods$normtop_consensus_control <- .condition_topic_project(control, memberships$consensus)
  }
  methods
}

.condition_topic_bootstrap_stability <- function(theta,
                                                 doc_mass,
                                                 expression,
                                                 tf_support,
                                                 phi,
                                                 alpha,
                                                 memberships,
                                                 observed,
                                                 reference_condition = NULL,
                                                 n_hvg = 500L,
                                                 n_bootstrap = 10L,
                                                 seed = 1L) {
  n_bootstrap <- as.integer(n_bootstrap)
  if (!is.finite(n_bootstrap) || n_bootstrap < 1L) {
    return(data.table::data.table(
      method = names(observed),
      bootstrap_primary_stability_mean = NA_real_,
      bootstrap_primary_stability_min = NA_real_
    ))
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  observed_primary <- lapply(observed, function(x) {
    score <- .condition_topic_row_normalize(x)
    stats::setNames(colnames(score)[max.col(score, ties.method = "first")], rownames(score))
  })
  hits <- lapply(observed, function(x) {
    matrix(0, nrow = nrow(x), ncol = n_bootstrap, dimnames = list(rownames(x), NULL))
  })
  doc_condition <- sub("::[^:]+$", "", rownames(theta))
  doc_groups <- split(seq_len(nrow(theta)), doc_condition)
  for (bootstrap in seq_len(n_bootstrap)) {
    gene_weight <- tabulate(
      sample.int(ncol(expression), ncol(expression), replace = TRUE),
      nbins = ncol(expression)
    )
    doc_index <- unlist(lapply(doc_groups, function(idx) sample(idx, length(idx), replace = TRUE)), use.names = FALSE)
    boot <- .condition_topic_build_methods(
      theta = theta[doc_index, , drop = FALSE],
      doc_mass = doc_mass,
      expression = sweep(expression, 2L, gene_weight, "*"),
      tf_support = tf_support,
      phi = phi,
      alpha = alpha,
      memberships = memberships,
      reference_condition = reference_condition,
      n_hvg = n_hvg,
      foldin_max_iter = 75L,
      foldin_tolerance = 1e-6
    )
    for (method in intersect(names(hits), names(boot))) {
      score <- .condition_topic_row_normalize(boot[[method]])
      primary <- stats::setNames(colnames(score)[max.col(score, ties.method = "first")], rownames(score))
      condition_order <- rownames(hits[[method]])
      hits[[method]][, bootstrap] <- primary[condition_order] == observed_primary[[method]][condition_order]
    }
  }
  data.table::rbindlist(lapply(names(hits), function(method) {
    per_condition <- rowMeans(hits[[method]], na.rm = TRUE)
    data.table::data.table(
      method = method,
      bootstrap_primary_stability_mean = mean(per_condition, na.rm = TRUE),
      bootstrap_primary_stability_min = min(per_condition, na.rm = TRUE)
    )
  }))
}

.condition_topic_bold_grob <- function(grob) {
  if (inherits(grob, "text")) {
    gp <- as.list(grob$gp)
    gp$font <- NULL
    grob$gp <- do.call(grid::gpar, utils::modifyList(gp, list(
      font = 2,
      fontfamily = "Helvetica",
      fontsize = max(9, gp$fontsize %||% 9)
    )))
  }
  if (!is.null(grob$grobs)) grob$grobs <- lapply(grob$grobs, .condition_topic_bold_grob)
  if (!is.null(grob$children)) {
    grob$children <- do.call(grid::gList, lapply(grob$children, .condition_topic_bold_grob))
  }
  grob
}

.plot_condition_topic_method_comparison <- function(methods,
                                                    out_file,
                                                    title_prefix = NULL,
                                                    condition_groups = NULL) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_abort("Package pheatmap is required for condition-topic comparison heatmaps.")
  }
  values <- unlist(lapply(methods, as.numeric), use.names = FALSE)
  color_max <- max(values[is.finite(values)], na.rm = TRUE)
  if (!is.finite(color_max) || color_max <= 0) color_max <- 1
  colors <- grDevices::colorRampPalette(c("#E8F2F8", "#BFDDEA", "#FFF7BC", "#FCAE61", "#D95F52"))(101)
  breaks <- seq(0, color_max, length.out = 102L)
  annotation_col <- NULL
  annotation_colors <- NULL
  if (!is.null(condition_groups)) {
    groups <- data.table::as.data.table(condition_groups)
    annotation_col <- data.frame(
      `Condition class` = groups$condition_group[match(
        rownames(methods[[1L]]),
        groups$condition_id
      )],
      check.names = FALSE
    )
    rownames(annotation_col) <- rownames(methods[[1L]])
    present <- unique(annotation_col[[1L]])
    annotation_colors <- list(`Condition class` = .condition_topic_group_colors(present))
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      out_file,
      width = 12,
      height = 9.5,
      family = "Helvetica",
      onefile = TRUE,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      out_file,
      width = 12,
      height = 9.5,
      family = "Helvetica",
      onefile = TRUE,
      useDingbats = FALSE,
      bg = "white"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  for (method in names(methods)) {
    score <- .condition_topic_row_normalize(methods[[method]])
    mat <- t(score)
    title <- paste(c(title_prefix, .condition_topic_method_label(method)), collapse = " | ")
    ph <- pheatmap::pheatmap(
      mat,
      color = colors,
      breaks = breaks,
      cluster_rows = nrow(mat) > 1L,
      cluster_cols = ncol(mat) > 1L,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 9,
      fontsize_row = 9,
      fontsize_col = 9,
      display_numbers = matrix(sprintf("%.2f", mat), nrow = nrow(mat), dimnames = dimnames(mat)),
      number_color = "#111111",
      fontsize_number = 9,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      angle_col = 90,
      main = title,
      border_color = "#FFFFFF",
      silent = TRUE
    )
    ph$gtable <- .condition_topic_bold_grob(ph$gtable)
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    grid::grid.draw(ph$gtable)
  }
  invisible(out_file)
}

.plot_condition_topic_group_metrics <- function(diagnostics, out_file) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  diagnostics <- data.table::as.data.table(diagnostics)
  metric_labels <- c(
    condition_group_silhouette_mean = "Stress-family silhouette",
    nearest_same_group_fraction = "Nearest condition in same family",
    between_within_js_ratio = "Between / within JS distance",
    bootstrap_primary_stability_mean = "Primary-topic bootstrap stability"
  )
  metric_columns <- intersect(names(metric_labels), names(diagnostics))
  if (!length(metric_columns)) {
    .log_abort("No condition-group metrics are available to plot.")
  }
  plot_data <- data.table::melt(
    diagnostics,
    id.vars = c("method", "model_backend", "biological_separation_rank"),
    measure.vars = metric_columns,
    variable.name = "metric",
    value.name = "value"
  )
  plot_data[, `:=`(
    metric = factor(metric_labels[as.character(metric)], levels = unname(metric_labels[metric_columns])),
    method_label = .condition_topic_method_label(method),
    method_family = data.table::fcase(
      startsWith(method, "theta_"), "Theta",
      startsWith(method, "phi_foldin"), "Fixed phi",
      default = "NormTop"
    ),
    value_label = ifelse(is.finite(value), sprintf("%.2f", value), "NA")
  )]
  backends <- unique(plot_data$model_backend)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(
      out_file,
      width = 16,
      height = 10.5,
      family = "Helvetica",
      onefile = TRUE,
      bg = "white"
    )
  } else {
    grDevices::pdf(
      out_file,
      width = 16,
      height = 10.5,
      family = "Helvetica",
      onefile = TRUE,
      useDingbats = FALSE,
      bg = "white"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  for (backend in backends) {
    page <- plot_data[model_backend == backend]
    method_order <- unique(
      diagnostics[model_backend == backend][order(biological_separation_rank)]$method
    )
    page[, method_label := factor(
      method_label,
      levels = rev(.condition_topic_method_label(method_order))
    )]
    plot <- ggplot2::ggplot(page, ggplot2::aes(x = value, y = method_label, fill = method_family)) +
      ggplot2::geom_col(width = 0.72, color = "#303030", linewidth = 0.25) +
      ggplot2::geom_text(
        ggplot2::aes(label = value_label, hjust = ifelse(value >= 0, -0.08, 1.08)),
        family = "Helvetica",
        fontface = "bold",
        size = 3.2
      ) +
      ggplot2::facet_wrap(~metric, scales = "free_x", nrow = 1) +
      ggplot2::scale_fill_manual(values = c(
        Theta = "#4477AA", `Fixed phi` = "#EE6677", NormTop = "#228833"
      )) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
      ggplot2::labs(
        title = sprintf("HPAFII nutrient-family condition-topic scoring | %s K10", toupper(backend)),
        subtitle = "Higher is better for every panel; method rows are ordered by the combined rank",
        x = NULL,
        y = NULL,
        fill = "Method family"
      ) +
      ggplot2::theme_bw(base_size = 9, base_family = "Helvetica") +
      ggplot2::theme(
        text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 9, face = "bold"),
        axis.text = ggplot2::element_text(size = 9, face = "bold", color = "#111111"),
        strip.text = ggplot2::element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
    grid::grid.newpage()
    grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
    print(plot, newpage = FALSE)
  }
  invisible(out_file)
}

.run_condition_topic_scoring_pilot <- function(theta_file,
                                               phi_file,
                                               score_file,
                                               topic_terms_file,
                                               condition_links_dir,
                                               output_dir,
                                               model_fit_file = NULL,
                                               reference_condition = NULL,
                                               title_prefix = NULL,
                                               model_backend = "lda",
                                               condition_groups = NULL,
                                               activity = NULL,
                                               n_bootstrap = 10L,
                                               n_hvg = 500L,
                                               review_methods = c(
                                                 "theta_mean",
                                                 "phi_foldin_hvg"
                                               ),
                                               verbose = TRUE) {
  .assert_pkg("data.table")
  theta <- .condition_topic_read_matrix_csv(theta_file)
  phi <- .condition_topic_read_matrix_csv(phi_file)
  score_mat <- readRDS(score_file)
  topic_terms <- data.table::fread(topic_terms_file, showProgress = FALSE)
  conditions <- unique(sub("::[^:]+$", "", rownames(theta)))
  if (is.null(activity)) {
    activity <- .condition_topic_collect_activity(
      condition_links_dir,
      conditions = conditions,
      verbose = verbose
    )
  }
  if (!all(c("gene_activity", "doc_mass") %in% names(activity))) {
    .log_abort("activity must contain gene_activity and doc_mass tables.")
  }
  gene_terms <- colnames(score_mat)[startsWith(colnames(score_mat), "GENE:")]
  genes <- sub("^GENE:", "", gene_terms)
  matrices <- .condition_topic_activity_matrices(activity$gene_activity, conditions, genes)
  expression <- matrices$expression
  memberships <- list(
    all = .condition_topic_membership(score_mat, topic_terms, mode = "all"),
    gamma = .condition_topic_membership(score_mat, topic_terms, mode = "gamma"),
    consensus = .condition_topic_membership(score_mat, topic_terms, mode = "consensus")
  )
  alpha <- NULL
  if (!is.null(model_fit_file) && file.exists(model_fit_file)) {
    fit <- readRDS(model_fit_file)
    alpha <- fit$alpha %||% NULL
    rm(fit)
  }
  methods <- .condition_topic_build_methods(
    theta = theta,
    doc_mass = activity$doc_mass,
    expression = expression,
    tf_support = matrices$tf_support,
    phi = phi,
    alpha = alpha,
    memberships = memberships,
    reference_condition = reference_condition,
    n_hvg = n_hvg
  )
  review_methods <- unique(as.character(review_methods))
  missing_review_methods <- setdiff(review_methods, names(methods))
  if (length(missing_review_methods)) {
    .log_abort(
      "Unknown condition-topic review method(s): {paste(missing_review_methods, collapse = ', ')}"
    )
  }
  methods <- methods[review_methods]
  method_info <- data.table::data.table(
    method = names(methods),
    method_label = .condition_topic_method_label(names(methods)),
    model_backend = as.character(model_backend),
    score_type = data.table::fcase(
      grepl("^theta_", names(methods)), "model_probability",
      grepl("^phi_foldin", names(methods)) & identical(tolower(model_backend), "lda"), "model_probability",
      grepl("^phi_foldin", names(methods)), "fixed_phi_projection",
      default = "normalized_score"
    ),
    reference_mode = data.table::fcase(
      grepl("one_vs_rest", names(methods)), "one_vs_rest",
      grepl("_hvg$", names(methods)), "hvg",
      grepl("_control$", names(methods)), "reference_condition",
      default = "absolute"
    ),
    reference_condition = ifelse(grepl("_control$", names(methods)), reference_condition %||% NA_character_, NA_character_)
  )
  method_info[, n_hvg := ifelse(grepl("_hvg$", method), as.integer(n_hvg), NA_integer_)]
  contributing <- list()
  membership_map <- c(
    normtop_all_absolute = "all", normtop_all_one_vs_rest = "all", normtop_all_control = "all",
    normtop_gamma_absolute = "gamma", normtop_gamma_one_vs_rest = "gamma", normtop_gamma_control = "gamma",
    normtop_consensus_absolute = "consensus", normtop_consensus_one_vs_rest = "consensus",
    normtop_consensus_control = "consensus", normtop_tf_support_one_vs_rest = "gamma",
    normtop_gamma_hvg = "gamma"
  )
  for (method in names(methods)) {
    source <- if (grepl("one_vs_rest", method)) {
      .condition_topic_positive_contrast(expression)
    } else if (grepl("_hvg$", method)) {
      .condition_topic_hvg_activity(expression)
    } else if (grepl("_control$", method)) {
      .condition_topic_positive_contrast(expression, reference = reference_condition)
    } else {
      expression
    }
    membership_name <- unname(membership_map[method])
    if (!is.na(membership_name)) {
      membership <- memberships[[membership_name]]
      contributing[[method]] <- rowSums(source[, rownames(membership), drop = FALSE] > 0 &
        rep(rowSums(membership) > 0, each = nrow(source)))
    }
  }
  diagnostics <- .condition_topic_diagnostics(methods, expression, contributing = contributing)
  if (!is.null(condition_groups)) {
    group_metrics <- .condition_topic_group_metrics(methods, condition_groups)
    diagnostics <- merge(diagnostics, group_metrics, by = "method", all.x = TRUE, sort = FALSE)
  }
  if (isTRUE(verbose)) .log_inform("Condition-topic scoring: running {n_bootstrap} bootstrap stability iterations.")
  stability <- .condition_topic_bootstrap_stability(
    theta = theta,
    doc_mass = activity$doc_mass,
    expression = expression,
    tf_support = matrices$tf_support,
    phi = phi,
    alpha = alpha,
    memberships = memberships,
    observed = methods,
    reference_condition = reference_condition,
    n_hvg = n_hvg,
    n_bootstrap = n_bootstrap,
    seed = 1L
  )
  diagnostics <- merge(diagnostics, stability, by = "method", all.x = TRUE, sort = FALSE)
  diagnostics <- merge(diagnostics, method_info, by = "method", all.x = TRUE, sort = FALSE)
  scores <- .condition_topic_long_scores(methods, method_info = method_info, contributing = contributing)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, sprintf("condition_topic_score_method_comparison_K%d.pdf", ncol(theta)))
  score_file_out <- file.path(output_dir, sprintf("condition_topic_scores_K%d.csv", ncol(theta)))
  diagnostic_file <- file.path(output_dir, sprintf("condition_topic_score_diagnostics_K%d.csv", ncol(theta)))
  .plot_condition_topic_method_comparison(
    methods,
    pdf_file,
    title_prefix = title_prefix,
    condition_groups = condition_groups
  )
  data.table::fwrite(scores, score_file_out)
  data.table::fwrite(diagnostics, diagnostic_file)
  invisible(list(
    pdf = pdf_file,
    scores = score_file_out,
    diagnostics = diagnostic_file,
    methods = methods
  ))
}
