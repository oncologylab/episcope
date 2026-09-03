# Condition-topic optimization and assignment QC.

.m3_opt_row_normalize <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  totals[!is.finite(totals) | totals <= 0] <- 1
  x / totals
}

.m3_opt_topic_ids <- function(x, margin = "row") {
  labels <- if (identical(margin, "row")) rownames(x) else colnames(x)
  ids <- suppressWarnings(as.integer(sub("^Topic", "", labels)))
  if (is.null(labels) || anyNA(ids) || anyDuplicated(ids)) {
    ids <- seq_len(if (identical(margin, "row")) nrow(x) else ncol(x))
  }
  as.integer(ids)
}

.m3_opt_hellinger_similarity <- function(phi, gene_ids, peak_ids) {
  modality_similarity <- function(indices) {
    indices <- unique(as.integer(indices[is.finite(indices)]))
    if (!length(indices)) {
      .log_abort("Topic optimization requires matched terms for both Gene and Peak modalities.")
    }
    x <- phi[, indices, drop = FALSE]
    x <- .m3_opt_row_normalize(x)
    sim <- tcrossprod(sqrt(x))
    sim[!is.finite(sim)] <- 0
    sim[sim < 0] <- 0
    sim[sim > 1] <- 1
    sim
  }
  (modality_similarity(gene_ids) + modality_similarity(peak_ids)) / 2
}

.m3_opt_doc_metadata <- function(doc_ids) {
  parts <- strsplit(as.character(doc_ids), "::", fixed = TRUE)
  ok <- lengths(parts) == 2L
  if (!all(ok)) {
    .log_abort(
      "Condition-topic optimization requires condition::TF document IDs."
    )
  }
  data.table::data.table(
    doc_id = as.character(doc_ids),
    condition_id = vapply(parts, `[[`, character(1L), 1L),
    tf = vapply(parts, `[[`, character(1L), 2L)
  )
}

.m3_opt_parse_topics <- function(x, valid_topics) {
  values <- strsplit(ifelse(is.na(x), "", as.character(x)), ";", fixed = TRUE)
  lapply(values, function(value) {
    out <- suppressWarnings(as.integer(value))
    sort(unique(out[is.finite(out) & out %in% valid_topics]))
  })
}

.m3_opt_link_universe <- function(dtm, pair_assignment) {
  .assert_pkg("Matrix")
  terms <- colnames(dtm)
  if (is.null(terms) || is.null(rownames(dtm))) {
    .log_abort("DTM must have document and term names for topic optimization.")
  }
  pairs <- data.table::as.data.table(pair_assignment)
  pairs <- pairs[
    !is.na(gene_term_id) & !is.na(peak_term_id) &
      gene_term_id %in% terms & peak_term_id %in% terms
  ]
  if (!nrow(pairs)) {
    .log_abort("No matched Gene/Peak terms are available for topic optimization.")
  }
  gene_columns <- match(pairs$gene_term_id, terms)
  peak_columns <- match(pairs$peak_term_id, terms)
  pairs[, pair_index := seq_len(.N)]
  target_levels <- sort(unique(as.character(pairs$target_gene)))
  pair_target_index <- match(pairs$target_gene, target_levels)
  gene_matrix <- dtm[, gene_columns, drop = FALSE]
  entries <- Matrix::summary(gene_matrix)
  if (!nrow(entries)) {
    .log_abort("No condition-TF-target links are represented in the DTM.")
  }
  peak_global <- peak_columns[entries$j]
  peak_count <- as.numeric(dtm[cbind(entries$i, peak_global)])
  keep <- is.finite(entries$x) & entries$x > 0 &
    is.finite(peak_count) & peak_count > 0
  entries <- entries[keep, , drop = FALSE]
  peak_count <- peak_count[keep]
  doc_meta <- .m3_opt_doc_metadata(rownames(dtm))
  list(
    links = data.table::data.table(
      link_index = seq_len(nrow(entries)),
      doc_index = as.integer(entries$i),
      pair_index = as.integer(entries$j),
      target_index = as.integer(pair_target_index[entries$j]),
      gene_token_count = as.numeric(entries$x),
      peak_token_count = peak_count,
      condition_id = doc_meta$condition_id[entries$i]
    ),
    docs = doc_meta,
    pairs = pairs,
    target_levels = target_levels,
    gene_columns = gene_columns,
    peak_columns = peak_columns
  )
}

.m3_opt_link_universe_tf_target <- function(dtm, pair_assignment) {
  .assert_pkg("Matrix")
  terms <- colnames(dtm)
  if (is.null(terms) || is.null(rownames(dtm))) {
    .log_abort("DTM must have document and term names for topic optimization.")
  }
  pairs <- data.table::copy(data.table::as.data.table(pair_assignment))
  required <- c(
    "tf", "target_gene", "gene_term_id", "peak_term_id"
  )
  missing <- setdiff(required, names(pairs))
  if (length(missing)) {
    .log_abort(
      "TF-target optimization assignments are missing: {paste(missing, collapse = ', ')}."
    )
  }
  pairs <- pairs[
    !is.na(gene_term_id) & !is.na(peak_term_id) &
      gene_term_id %in% terms & peak_term_id %in% terms
  ]
  if (!nrow(pairs)) {
    .log_abort("No matched Gene and TF-target terms are available for optimization.")
  }
  pairs[, pair_index := seq_len(.N)]
  gene_columns <- match(pairs$gene_term_id, terms)
  link_columns <- match(pairs$peak_term_id, terms)
  entries <- Matrix::summary(dtm[, link_columns, drop = FALSE])
  if (!nrow(entries)) {
    .log_abort("No condition-specific TF-target terms are represented in the DTM.")
  }
  gene_count <- as.numeric(dtm[cbind(
    entries$i,
    gene_columns[entries$j]
  )])
  keep <- is.finite(entries$x) & entries$x > 0 &
    is.finite(gene_count) & gene_count > 0
  entries <- entries[keep, , drop = FALSE]
  gene_count <- gene_count[keep]
  doc_meta <- .m3_opt_doc_metadata(rownames(dtm))
  target_levels <- sort(unique(as.character(pairs$target_gene)))
  pair_target_index <- match(pairs$target_gene, target_levels)
  list(
    links = data.table::data.table(
      link_index = seq_len(nrow(entries)),
      doc_index = as.integer(entries$i),
      pair_index = as.integer(entries$j),
      target_index = as.integer(pair_target_index[entries$j]),
      gene_token_count = gene_count,
      peak_token_count = as.numeric(entries$x),
      condition_id = doc_meta$condition_id[entries$i]
    ),
    docs = doc_meta,
    pairs = pairs,
    target_levels = target_levels,
    gene_columns = gene_columns,
    peak_columns = link_columns
  )
}

.m3_opt_probability_metrics <- function(probability) {
  primary <- max.col(probability, ties.method = "first")
  primary_probability <- probability[cbind(seq_len(nrow(probability)), primary)]
  probability[cbind(seq_len(nrow(probability)), primary)] <- -Inf
  secondary <- apply(probability, 1L, max)
  secondary[!is.finite(secondary)] <- 0
  list(
    topic_index = as.integer(primary),
    probability = as.numeric(primary_probability),
    margin = as.numeric(primary_probability - secondary)
  )
}

.m3_opt_link_posteriors <- function(theta,
                                    phi,
                                    universe,
                                    group_index = seq_len(ncol(theta)),
                                    sample_rows = integer(),
                                    chunk_size = 50000L) {
  links <- universe$links
  n <- nrow(links)
  group_index <- as.integer(group_index)
  groups <- sort(unique(group_index))
  group_position <- match(group_index, groups)
  primary_group <- integer(n)
  primary_probability <- numeric(n)
  primary_margin <- numeric(n)
  sample_rows <- sort(unique(as.integer(sample_rows)))
  sample_probability <- matrix(
    0,
    nrow = length(sample_rows),
    ncol = length(groups)
  )
  sample_position <- integer(n)
  if (length(sample_rows)) sample_position[sample_rows] <- seq_along(sample_rows)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / as.integer(chunk_size)))
  for (rows in chunks) {
    link_rows <- links[rows]
    pair_index <- if ("pair_index" %in% names(link_rows)) {
      link_rows$pair_index
    } else {
      link_rows$target_index
    }
    theta_chunk <- theta[link_rows$doc_index, , drop = FALSE]
    gene_phi <- t(phi[, universe$gene_columns[pair_index], drop = FALSE])
    peak_phi <- t(phi[, universe$peak_columns[pair_index], drop = FALSE])
    gene_posterior <- .m3_opt_row_normalize(theta_chunk * gene_phi)
    peak_posterior <- .m3_opt_row_normalize(theta_chunk * peak_phi)
    probability <- .m3_opt_row_normalize(
      gene_posterior * link_rows$gene_token_count +
        peak_posterior * link_rows$peak_token_count
    )
    if (length(groups) != ncol(probability)) {
      collapsed <- matrix(0, nrow(probability), length(groups))
      for (j in seq_along(group_position)) {
        collapsed[, group_position[[j]]] <- collapsed[, group_position[[j]]] +
          probability[, j]
      }
      probability <- collapsed
    }
    metrics <- .m3_opt_probability_metrics(probability)
    primary_group[rows] <- groups[metrics$topic_index]
    primary_probability[rows] <- metrics$probability
    primary_margin[rows] <- metrics$margin
    sampled <- sample_position[rows] > 0L
    if (any(sampled)) {
      sample_probability[sample_position[rows[sampled]], ] <- probability[sampled, , drop = FALSE]
    }
  }
  colnames(sample_probability) <- paste0("Topic", groups)
  list(
    topic = primary_group,
    probability = primary_probability,
    margin = primary_margin,
    sample_probability = sample_probability
  )
}

.m3_project_peaks_from_gene_topics <- function(theta,
                                               phi,
                                               edges,
                                               chunk_size = 50000L) {
  .assert_pkg("data.table")
  theta <- as.matrix(theta)
  phi <- as.matrix(phi)
  if (is.null(rownames(theta)) || is.null(colnames(theta)) ||
      is.null(rownames(phi)) || is.null(colnames(phi))) {
    .log_abort(
      "Gene-topic Peak projection requires named theta and phi matrices."
    )
  }
  if (!setequal(colnames(theta), rownames(phi))) {
    .log_abort("Theta topics and phi topics do not match for Peak projection.")
  }
  phi <- phi[colnames(theta), , drop = FALSE]
  gene_terms <- startsWith(colnames(phi), "GENE:")
  if (!any(gene_terms)) {
    .log_abort("Gene-topic Peak projection requires GENE terms in phi.")
  }
  gene_ids <- sub("^GENE:", "", colnames(phi)[gene_terms])
  if (anyDuplicated(gene_ids)) {
    .log_abort("Gene-topic Peak projection found duplicated Gene terms in phi.")
  }
  chunk_size <- as.integer(chunk_size)
  if (!is.finite(chunk_size) || chunk_size < 1L) {
    .log_abort("Peak projection chunk_size must be a positive integer.")
  }

  required <- c(
    "condition_id", "doc_id", "tf", "peak_id", "target_gene",
    "tf_target_r", "fp_target_r", "tf_peak_r", "fp_score",
    "tf_expression"
  )
  x <- data.table::copy(data.table::as.data.table(edges))
  .assert_has_cols(x, required, context = "Gene-topic Peak projection edges")
  x <- x[, ..required]
  input_rows <- nrow(x)
  for (column in c("condition_id", "doc_id", "tf", "peak_id", "target_gene")) {
    data.table::set(x, j = column, value = as.character(x[[column]]))
  }
  text_ok <- Reduce(`&`, lapply(
    c("condition_id", "doc_id", "tf", "peak_id", "target_gene"),
    function(column) !is.na(x[[column]]) & nzchar(x[[column]])
  ))
  numeric_columns <- c(
    "tf_target_r", "fp_target_r", "tf_peak_r", "fp_score",
    "tf_expression"
  )
  for (column in numeric_columns) {
    data.table::set(
      x,
      j = column,
      value = suppressWarnings(as.numeric(x[[column]]))
    )
  }
  evidence_ok <- Reduce(`&`, lapply(
    numeric_columns,
    function(column) is.finite(x[[column]]) & x[[column]] > 0
  ))
  valid_input_rows <- sum(text_ok & evidence_ok)
  x <- x[text_ok & evidence_ok]
  if (!nrow(x)) {
    .log_abort("No valid positive-evidence edges remain for Peak projection.")
  }

  doc_index <- stats::setNames(seq_len(nrow(theta)), rownames(theta))
  gene_index <- stats::setNames(seq_along(gene_ids), gene_ids)
  x[, `:=`(
    doc_index = unname(doc_index[doc_id]),
    gene_index = unname(gene_index[target_gene]),
    target_evidence = tf_target_r * fp_target_r
  )]
  matched_rows <- x[
    is.finite(doc_index) & is.finite(gene_index) &
      is.finite(target_evidence) & target_evidence > 0
  ]
  if (!nrow(matched_rows)) {
    .log_abort("No Peak-projection edges match both theta documents and phi Genes.")
  }
  x <- matched_rows
  data.table::setorderv(
    x,
    c(
      "condition_id", "doc_id", "tf", "peak_id", "target_gene",
      "target_evidence", "tf_peak_r", "fp_score", "tf_expression"
    ),
    c(rep(1L, 5L), rep(-1L, 4L)),
    na.last = TRUE
  )
  x <- x[!duplicated(
    x,
    by = c("condition_id", "doc_id", "tf", "peak_id", "target_gene")
  )]
  deduplicated_rows <- nrow(x)
  x[, group_index := .GRP, by = .(condition_id, doc_id, tf, peak_id)]

  gene_probability <- .m3_opt_row_normalize(t(phi[, gene_terms, drop = FALSE]))
  group_count <- max(x$group_index)
  topic_count <- ncol(theta)
  group_score <- matrix(0, nrow = group_count, ncol = topic_count)
  starts <- seq.int(1L, nrow(x), by = chunk_size)
  for (start in starts) {
    rows <- start:min(start + chunk_size - 1L, nrow(x))
    probability <- .m3_opt_row_normalize(
      theta[x$doc_index[rows], , drop = FALSE] *
        gene_probability[x$gene_index[rows], , drop = FALSE]
    )
    contribution <- probability * x$target_evidence[rows]
    collapsed <- rowsum(
      contribution,
      x$group_index[rows],
      reorder = FALSE
    )
    group_rows <- as.integer(rownames(collapsed))
    group_score[group_rows, ] <-
      group_score[group_rows, , drop = FALSE] + collapsed
  }
  group_target_total <- x[, .(
    target_evidence_total = sum(target_evidence),
    target_evidence_mean = mean(target_evidence),
    n_targets = .N,
    tf_peak_r = max(tf_peak_r),
    fp_score = max(fp_score),
    tf_expression = max(tf_expression)
  ), by = .(group_index, condition_id, doc_id, tf, peak_id)]
  data.table::setorder(group_target_total, group_index)
  group_probability <- group_score / group_target_total$target_evidence_total
  group_probability <- .m3_opt_row_normalize(group_probability)
  group_evidence <- with(
    group_target_total,
    tf_peak_r * fp_score * tf_expression * target_evidence_mean
  )
  group_metrics <- .m3_opt_probability_metrics(group_probability)
  topic_ids <- .m3_opt_topic_ids(theta, "column")
  group_assignment <- group_target_total[, .(
    condition_id,
    doc_id,
    tf,
    peak_id,
    n_targets,
    projection_evidence = group_evidence,
    primary_topic = topic_ids[group_metrics$topic_index],
    primary_probability = group_metrics$probability,
    primary_margin = group_metrics$margin
  )]

  condition_peak <- unique(group_assignment[, .(condition_id, peak_id)])
  condition_peak[, condition_peak_index := seq_len(.N)]
  group_assignment[
    condition_peak,
    condition_peak_index := i.condition_peak_index,
    on = c("condition_id", "peak_id")
  ]
  weighted_group_probability <- group_probability *
    group_assignment$projection_evidence
  condition_peak_score <- rowsum(
    weighted_group_probability,
    group_assignment$condition_peak_index,
    reorder = FALSE
  )
  condition_peak_probability <- .m3_opt_row_normalize(condition_peak_score)
  condition_summary <- group_assignment[, .(
    n_tf_groups = .N,
    n_target_links = sum(n_targets),
    projection_evidence = sum(projection_evidence)
  ), by = .(condition_peak_index, condition_id, peak_id)]
  data.table::setorder(condition_summary, condition_peak_index)
  probability_names <- paste0("Topic", topic_ids)
  colnames(condition_peak_probability) <- probability_names
  condition_peak_probabilities <- cbind(
    condition_summary[, !"condition_peak_index"],
    data.table::as.data.table(condition_peak_probability)
  )
  condition_metrics <- .m3_opt_probability_metrics(condition_peak_probability)
  condition_peak_probabilities[, `:=`(
    primary_topic = topic_ids[condition_metrics$topic_index],
    primary_probability = condition_metrics$probability,
    primary_margin = condition_metrics$margin
  )]
  data.table::setcolorder(
    condition_peak_probabilities,
    c(
      "condition_id", "peak_id", "n_tf_groups", "n_target_links",
      "projection_evidence", "primary_topic", "primary_probability",
      "primary_margin", probability_names
    )
  )

  audit <- data.table::data.table(
    input_rows = input_rows,
    valid_positive_evidence_rows = valid_input_rows,
    matched_rows = nrow(matched_rows),
    deduplicated_rows = deduplicated_rows,
    tf_peak_groups = nrow(group_assignment),
    condition_peaks = nrow(condition_peak_probabilities),
    unique_peaks = data.table::uniqueN(condition_peak_probabilities$peak_id)
  )
  list(
    condition_peak_probabilities = condition_peak_probabilities[],
    tf_peak_assignments = group_assignment[],
    audit = audit[]
  )
}

.m3_aggregate_projected_peak_topics <- function(condition_peak_probabilities) {
  .assert_pkg("data.table")
  x <- data.table::copy(data.table::as.data.table(
    condition_peak_probabilities
  ))
  required <- c("condition_id", "peak_id", "projection_evidence")
  .assert_has_cols(
    x,
    required,
    context = "condition-level projected Peak probabilities"
  )
  topic_columns <- grep("^Topic[0-9]+$", names(x), value = TRUE)
  if (!length(topic_columns)) {
    .log_abort("Projected Peak probabilities do not contain Topic columns.")
  }
  topic_ids <- suppressWarnings(as.integer(sub("^Topic", "", topic_columns)))
  if (anyNA(topic_ids) || anyDuplicated(topic_ids)) {
    .log_abort("Projected Peak probability Topic columns are not unique.")
  }
  data.table::setorderv(x, c("peak_id", "condition_id"))
  valid <- !is.na(x$peak_id) & nzchar(as.character(x$peak_id)) &
    is.finite(x$projection_evidence) & x$projection_evidence > 0
  x <- x[valid]
  if (!nrow(x)) {
    .log_abort("No positive-evidence condition Peaks remain for aggregation.")
  }
  probability <- .m3_opt_row_normalize(as.matrix(x[, ..topic_columns]))
  x[, peak_index := .GRP, by = peak_id]
  peak_score <- rowsum(
    probability * x$projection_evidence,
    x$peak_index,
    reorder = FALSE
  )
  peak_probability <- .m3_opt_row_normalize(peak_score)
  peak_summary <- x[, .(
    n_conditions = data.table::uniqueN(condition_id),
    projection_evidence = sum(projection_evidence)
  ), by = .(peak_index, peak_id)]
  data.table::setorder(peak_summary, peak_index)
  colnames(peak_probability) <- topic_columns
  out <- cbind(
    peak_summary[, !"peak_index"],
    data.table::as.data.table(peak_probability)
  )
  metrics <- .m3_opt_probability_metrics(peak_probability)
  out[, `:=`(
    primary_topic = topic_ids[metrics$topic_index],
    primary_probability = metrics$probability,
    primary_margin = metrics$margin
  )]
  data.table::setcolorder(
    out,
    c(
      "peak_id", "n_conditions", "projection_evidence", "primary_topic",
      "primary_probability", "primary_margin", topic_columns
    )
  )
  out[]
}

.m3_opt_balanced_sample <- function(condition_id,
                                    max_per_condition = 10000L,
                                    seed = 20260716L) {
  max_per_condition <- as.integer(max_per_condition)
  seed <- as.integer(seed)
  parts <- split(seq_along(condition_id), condition_id)
  unlist(lapply(seq_along(parts), function(i) {
    rows <- parts[[i]]
    if (length(rows) <= max_per_condition) return(rows)
    set.seed(seed + i)
    sort(sample(rows, max_per_condition))
  }), use.names = FALSE)
}

.m3_opt_target_assignment <- function(phi,
                                      pair_assignment,
                                      raw_topic_ids,
                                      raw_to_group) {
  pairs <- data.table::copy(data.table::as.data.table(pair_assignment))
  raw_to_group <- stats::setNames(as.integer(raw_to_group), raw_topic_ids)
  candidate_gene <- .m3_opt_parse_topics(pairs$gene_gammafit_topics, raw_topic_ids)
  candidate_peak <- .m3_opt_parse_topics(pairs$peak_gammafit_topics, raw_topic_ids)
  group_ids <- sort(unique(raw_to_group))
  group_position <- stats::setNames(seq_along(group_ids), group_ids)
  raw_position <- stats::setNames(seq_along(raw_topic_ids), raw_topic_ids)
  phi_term <- .m3_opt_row_normalize(t(phi))

  select_group <- function(term_ids, candidates) {
    term_position <- match(term_ids, colnames(phi))
    selected <- rep(NA_integer_, length(term_position))
    selected_probability <- rep(NA_real_, length(term_position))
    for (i in seq_along(term_position)) {
      if (!is.finite(term_position[[i]]) || !length(candidates[[i]])) next
      raw_candidates <- candidates[[i]]
      candidate_groups <- unique(raw_to_group[as.character(raw_candidates)])
      values <- numeric(length(candidate_groups))
      for (j in seq_along(candidate_groups)) {
        members <- raw_topic_ids[raw_to_group == candidate_groups[[j]]]
        values[[j]] <- sum(phi_term[term_position[[i]], raw_position[as.character(members)]])
      }
      winner <- which.max(values)
      selected[[i]] <- candidate_groups[[winner]]
      selected_probability[[i]] <- values[[winner]] / sum(values)
    }
    list(topic = selected, probability = selected_probability)
  }

  gene <- select_group(pairs$gene_term_id, candidate_gene)
  peak <- select_group(pairs$peak_term_id, candidate_peak)
  raw_status <- if ("assignment_status" %in% names(pairs)) {
    as.character(pairs$assignment_status)
  } else {
    rep(NA_character_, nrow(pairs))
  }
  pairs[, raw_assignment_status := raw_status]
  pairs[, `:=`(
    raw_assigned_topic = suppressWarnings(as.integer(assigned_topic)),
    optimized_gene_topic = gene$topic,
    optimized_peak_topic = peak$topic,
    optimized_gene_probability = gene$probability,
    optimized_peak_probability = peak$probability
  )]
  pairs[, optimized_assigned_topic := data.table::fifelse(
    is.finite(optimized_gene_topic) &
      optimized_gene_topic == optimized_peak_topic,
    optimized_gene_topic,
    NA_integer_
  )]
  pairs[, recovered_after_merge := is.na(raw_assigned_topic) &
    is.finite(optimized_assigned_topic)]
  pairs[, assigned := is.finite(optimized_assigned_topic)]
  pairs[, assigned_topic := data.table::fifelse(
    assigned,
    as.character(optimized_assigned_topic),
    NA_character_
  )]
  pairs[, assignment_status := data.table::fcase(
    recovered_after_merge, "assigned_after_topic_merge",
    assigned, "assigned_gammafit_maxprob_agreement",
    default = "unassigned_after_topic_merge"
  )]
  pairs[]
}

.m3_opt_tf_identity_key <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[out == "TBET"] <- "TBX21"
  out
}

.m3_opt_max_theta <- function(theta, raw_topic_ids, raw_to_group) {
  if (ncol(theta) != length(raw_topic_ids) ||
      length(raw_to_group) != length(raw_topic_ids)) {
    .log_abort(
      "Raw theta columns, topic IDs, and optimization mapping must align."
    )
  }
  if (!is.null(names(raw_to_group))) {
    mapped_positions <- match(as.character(raw_topic_ids), names(raw_to_group))
    if (anyNA(mapped_positions)) {
      .log_abort(
        "The named optimization mapping does not cover every raw topic."
      )
    }
    raw_to_group <- raw_to_group[mapped_positions]
  }
  groups <- sort(unique(raw_to_group))
  out <- matrix(0, nrow(theta), length(groups))
  for (i in seq_along(groups)) {
    out[, i] <- apply(
      theta[, raw_to_group == groups[[i]], drop = FALSE],
      1L,
      max
    )
  }
  rownames(out) <- rownames(theta)
  colnames(out) <- paste0("Topic", groups)
  out
}

.m3_opt_tf_theta_correspondence <- function(theta,
                                            phi,
                                            pair_assignment,
                                            raw_topic_ids,
                                            raw_to_group,
                                            cutoff = 0.3) {
  empty_result <- list(
    available = FALSE,
    target_assignments = NA_integer_,
    tf_term_assignments = NA_integer_,
    supported_tf_terms = NA_integer_,
    empty_tf_terms = NA_integer_,
    support_rate = NA_real_,
    mean_theta = NA_real_
  )
  if (is.null(pair_assignment) ||
      !"target_gene" %in% names(pair_assignment) ||
      is.null(rownames(theta))) {
    return(empty_result)
  }
  docs <- .m3_opt_doc_metadata(rownames(theta))
  docs[, tf_key := .m3_opt_tf_identity_key(tf)]
  assignments <- .m3_opt_target_assignment(
    phi = phi,
    pair_assignment = pair_assignment,
    raw_topic_ids = raw_topic_ids,
    raw_to_group = raw_to_group
  )
  assignments[, tf_key := .m3_opt_tf_identity_key(target_gene)]
  tf_assignments <- unique(assignments[
    .as_logical_flag(assigned) & tf_key %in% docs$tf_key,
    .(
      tf_key,
      optimized_topic = suppressWarnings(as.integer(optimized_assigned_topic))
    )
  ])
  if (!nrow(tf_assignments)) {
    empty_result$target_assignments <- sum(.as_logical_flag(assignments$assigned))
    return(empty_result)
  }
  optimized_theta <- .m3_opt_max_theta(theta, raw_topic_ids, raw_to_group)
  optimized_topic_ids <- .m3_opt_topic_ids(optimized_theta, "column")
  theta_max <- vapply(seq_len(nrow(tf_assignments)), function(i) {
    topic_position <- match(
      tf_assignments$optimized_topic[[i]],
      optimized_topic_ids
    )
    doc_rows <- which(docs$tf_key == tf_assignments$tf_key[[i]])
    if (!length(doc_rows) || !is.finite(topic_position)) return(NA_real_)
    max(optimized_theta[doc_rows, topic_position], na.rm = TRUE)
  }, numeric(1L))
  supported <- is.finite(theta_max) & theta_max >= cutoff
  list(
    available = TRUE,
    target_assignments = sum(.as_logical_flag(assignments$assigned)),
    tf_term_assignments = length(theta_max),
    supported_tf_terms = sum(supported),
    empty_tf_terms = sum(!supported),
    support_rate = mean(supported),
    mean_theta = mean(theta_max[is.finite(theta_max)])
  )
}

.m3_opt_group_phi <- function(phi, theta, dtm, raw_topic_ids, raw_to_group) {
  doc_tokens <- as.numeric(Matrix::rowSums(dtm))
  topic_mass <- colSums(theta * doc_tokens)
  topic_mass[!is.finite(topic_mass) | topic_mass <= 0] <- 1
  groups <- sort(unique(raw_to_group))
  out <- matrix(0, nrow = length(groups), ncol = ncol(phi))
  for (i in seq_along(groups)) {
    members <- which(raw_to_group == groups[[i]])
    weights <- topic_mass[members]
    out[i, ] <- colSums(phi[members, , drop = FALSE] * weights) / sum(weights)
  }
  rownames(out) <- paste0("Topic", groups)
  colnames(out) <- colnames(phi)
  .m3_opt_row_normalize(out)
}

.m3_opt_merge_map <- function(phi,
                              theta,
                              dtm,
                              raw_topic_ids,
                              raw_links,
                              raw_genes,
                              gene_ids,
                              peak_ids,
                              pair_assignment = NULL,
                              min_genes = 150L,
                              min_links = 200L,
                              similarity_threshold = 0.65,
                              prefer_tf_theta_correspondence = TRUE,
                              tf_topic_cutoff = 0.3,
                              min_topics = 2L) {
  mapping <- stats::setNames(raw_topic_ids, raw_topic_ids)
  audit <- list()
  step <- 0L
  use_correspondence <- isTRUE(prefer_tf_theta_correspondence) &&
    !is.null(pair_assignment)
  repeat {
    group_phi <- .m3_opt_group_phi(phi, theta, dtm, raw_topic_ids, mapping)
    groups <- .m3_opt_topic_ids(group_phi, "row")
    links <- vapply(groups, function(g) sum(raw_links[mapping == g]), numeric(1L))
    genes <- vapply(groups, function(g) sum(raw_genes[mapping == g]), numeric(1L))
    similarity <- .m3_opt_hellinger_similarity(group_phi, gene_ids, peak_ids)
    diag(similarity) <- -Inf
    small <- which(links < min_links | genes < min_genes)
    if (length(groups) <= min_topics) break
    candidate_positions <- if (length(small)) {
      which(
        row(similarity) %in% small &
          similarity >= similarity_threshold,
        arr.ind = TRUE
      )
    } else {
      which(
        upper.tri(similarity) &
          similarity >= similarity_threshold,
        arr.ind = TRUE
      )
    }
    if (!nrow(candidate_positions)) break
    candidates <- lapply(seq_len(nrow(candidate_positions)), function(i) {
      source_pos <- candidate_positions[[i, 1L]]
      target_pos <- candidate_positions[[i, 2L]]
      source_group <- groups[[source_pos]]
      target_group <- groups[[target_pos]]
      source_size <- c(links[[source_pos]], genes[[source_pos]])
      target_size <- c(links[[target_pos]], genes[[target_pos]])
      if (source_size[[1L]] > target_size[[1L]] ||
          (source_size[[1L]] == target_size[[1L]] &&
            source_size[[2L]] > target_size[[2L]]) ||
          (identical(source_size, target_size) && source_group < target_group)) {
        tmp <- source_group
        source_group <- target_group
        target_group <- tmp
        tmp <- source_pos
        source_pos <- target_pos
        target_pos <- tmp
      }
      data.table::data.table(
        source_pos = source_pos,
        target_pos = target_pos,
        source_topic = source_group,
        representative_topic = target_group,
        reason = if (length(small)) {
          data.table::fcase(
            links[[source_pos]] < min_links & genes[[source_pos]] < min_genes,
            "small_links_and_genes",
            links[[source_pos]] < min_links,
            "small_links",
            default = "small_genes"
          )
        } else {
          "high_similarity"
        },
        similarity = similarity[source_pos, target_pos],
        source_links = links[[source_pos]],
        source_genes = genes[[source_pos]],
        representative_links_before = links[[target_pos]],
        representative_genes_before = genes[[target_pos]]
      )
    })
    candidates <- unique(data.table::rbindlist(candidates), by = c(
      "source_topic", "representative_topic"
    ))
    correspondence_before <- .m3_opt_tf_theta_correspondence(
      theta = theta,
      phi = phi,
      pair_assignment = pair_assignment,
      raw_topic_ids = raw_topic_ids,
      raw_to_group = mapping,
      cutoff = tf_topic_cutoff
    )
    candidates[, `:=`(
      target_assignments_before = correspondence_before$target_assignments,
      tf_term_assignments_before = correspondence_before$tf_term_assignments,
      tf_theta_supported_before = correspondence_before$supported_tf_terms,
      tf_theta_empty_before = correspondence_before$empty_tf_terms,
      tf_theta_support_rate_before = correspondence_before$support_rate,
      tf_theta_mean_before = correspondence_before$mean_theta,
      target_assignments_after = NA_integer_,
      tf_term_assignments_after = NA_integer_,
      tf_theta_supported_after = NA_integer_,
      tf_theta_empty_after = NA_integer_,
      tf_theta_support_rate_after = NA_real_,
      tf_theta_mean_after = NA_real_,
      correspondence_eligible = TRUE
    )]
    if (isTRUE(correspondence_before$available)) {
      for (i in seq_len(nrow(candidates))) {
        candidate_mapping <- mapping
        candidate_mapping[
          candidate_mapping == candidates$source_topic[[i]]
        ] <- candidates$representative_topic[[i]]
        after <- .m3_opt_tf_theta_correspondence(
          theta = theta,
          phi = phi,
          pair_assignment = pair_assignment,
          raw_topic_ids = raw_topic_ids,
          raw_to_group = candidate_mapping,
          cutoff = tf_topic_cutoff
        )
        candidates[i, `:=`(
          target_assignments_after = after$target_assignments,
          tf_term_assignments_after = after$tf_term_assignments,
          tf_theta_supported_after = after$supported_tf_terms,
          tf_theta_empty_after = after$empty_tf_terms,
          tf_theta_support_rate_after = after$support_rate,
          tf_theta_mean_after = after$mean_theta,
          correspondence_eligible = isTRUE(after$available) &&
            after$target_assignments >=
              correspondence_before$target_assignments &&
            after$supported_tf_terms >=
              correspondence_before$supported_tf_terms &&
            after$empty_tf_terms <= correspondence_before$empty_tf_terms
        )]
      }
    }
    if (use_correspondence && isTRUE(correspondence_before$available)) {
      candidates <- candidates[
        candidates[["correspondence_eligible"]] %in% TRUE
      ]
      if (!nrow(candidates)) break
      data.table::setorderv(
        candidates,
        c(
          "tf_theta_empty_after",
          "tf_theta_supported_after",
          "target_assignments_after",
          "tf_theta_support_rate_after",
          "tf_theta_mean_after",
          "similarity",
          "source_topic",
          "representative_topic"
        ),
        c(1L, -1L, -1L, -1L, -1L, -1L, 1L, 1L)
      )
    } else if (length(small)) {
      data.table::setorderv(
        candidates,
        c(
          "source_links",
          "source_genes",
          "source_topic",
          "similarity",
          "representative_topic"
        ),
        c(1L, 1L, 1L, -1L, 1L)
      )
    } else {
      data.table::setorderv(
        candidates,
        c("similarity", "source_topic", "representative_topic"),
        c(-1L, 1L, 1L)
      )
    }
    winner <- candidates[1L]
    step <- step + 1L
    audit[[step]] <- cbind(
      data.table::data.table(merge_step = step),
      winner[, -c("source_pos", "target_pos", "correspondence_eligible")]
    )
    mapping[mapping == winner[["source_topic"]]] <-
      winner[["representative_topic"]]
  }
  list(
    mapping = as.integer(mapping),
    audit = if (length(audit)) {
      data.table::rbindlist(audit, use.names = TRUE, fill = TRUE)
    } else {
      data.table::data.table(
        merge_step = integer(),
        source_topic = integer(),
        representative_topic = integer(),
        reason = character(),
        similarity = numeric(),
        source_links = numeric(),
        source_genes = numeric(),
        representative_links_before = numeric(),
        representative_genes_before = numeric(),
        target_assignments_before = integer(),
        tf_term_assignments_before = integer(),
        tf_theta_supported_before = integer(),
        tf_theta_empty_before = integer(),
        tf_theta_support_rate_before = numeric(),
        tf_theta_mean_before = numeric(),
        target_assignments_after = integer(),
        tf_term_assignments_after = integer(),
        tf_theta_supported_after = integer(),
        tf_theta_empty_after = integer(),
        tf_theta_support_rate_after = numeric(),
        tf_theta_mean_after = numeric()
      )
    }
  )
}

.m3_opt_theta <- function(theta, raw_topic_ids, raw_to_group) {
  groups <- sort(unique(raw_to_group))
  out <- matrix(0, nrow(theta), length(groups))
  for (i in seq_along(groups)) {
    out[, i] <- rowSums(theta[, raw_to_group == groups[[i]], drop = FALSE])
  }
  rownames(out) <- rownames(theta)
  colnames(out) <- paste0("Topic", groups)
  .m3_opt_row_normalize(out)
}

.m3_opt_topic_terms <- function(raw_topic_terms,
                                optimized_phi,
                                optimized_score,
                                raw_topic_ids,
                                raw_to_group,
                                optimized_pairs,
                                optimized_gene_pairs = NULL) {
  terms <- data.table::as.data.table(raw_topic_terms)
  if (!"topic_num" %in% names(terms)) {
    terms[, topic_num := suppressWarnings(as.integer(sub("^Topic", "", topic)))]
  }
  terms[, raw_topic_num := as.integer(topic_num)]
  map <- stats::setNames(raw_to_group, raw_topic_ids)
  terms[, optimized_topic_num := unname(map[as.character(raw_topic_num)])]
  terms <- terms[is.finite(optimized_topic_num)]
  if (!"term_group" %in% names(terms)) terms[, term_group := .term_group(term_id)]
  if (!"gammafit_candidate" %in% names(terms)) {
    terms[, gammafit_candidate := .as_logical_flag(in_topic)]
  }
  collapsed <- terms[, .(
    score = suppressWarnings(max(as.numeric(score), na.rm = TRUE)),
    gammafit_candidate = any(.as_logical_flag(gammafit_candidate)),
    raw_topic_members = paste(sort(unique(raw_topic_num)), collapse = ";")
  ), by = .(
    topic_num = optimized_topic_num,
    term_id,
    term_group
  )]
  collapsed[!is.finite(score), score := 0]
  topic_position <- match(collapsed$topic_num, .m3_opt_topic_ids(optimized_phi, "row"))
  term_position <- match(collapsed$term_id, colnames(optimized_phi))
  valid <- is.finite(topic_position) & is.finite(term_position)
  collapsed[, `:=`(phi = 0, optimized_score = 0)]
  collapsed$phi[valid] <- optimized_phi[cbind(topic_position[valid], term_position[valid])]
  collapsed$optimized_score[valid] <- optimized_score[
    cbind(topic_position[valid], term_position[valid])
  ]
  pair_terms <- data.table::as.data.table(optimized_pairs)
  gene_pair_terms <- data.table::as.data.table(
    optimized_gene_pairs %||% optimized_pairs
  )
  assignment_map <- data.table::rbindlist(list(
    unique(gene_pair_terms[
      is.finite(optimized_gene_topic),
      .(
        term_id = as.character(gene_term_id),
        term_group = "GENE",
        assigned_topic__ = as.integer(optimized_gene_topic)
      )
    ]),
    unique(pair_terms[
      is.finite(optimized_assigned_topic),
      .(
        term_id = as.character(peak_term_id),
        term_group = data.table::fifelse(
          grepl("^[^:]+::[^:]+$", peak_term_id),
          "TF_TARGET",
          "PEAK"
        ),
        assigned_topic__ = as.integer(optimized_assigned_topic)
      )
    ])
  ), use.names = TRUE, fill = TRUE)
  conflicts <- assignment_map[, data.table::uniqueN(assigned_topic__),
    by = term_id
  ][V1 > 1L]
  if (nrow(conflicts)) {
    .log_abort("Optimized term assignment maps a term to multiple topics.")
  }
  assignment_map <- unique(assignment_map, by = "term_id")
  assigned_topic_map <- stats::setNames(
    assignment_map$assigned_topic__,
    assignment_map$term_id
  )
  collapsed[, in_topic := data.table::fcase(
    term_group %in% c("GENE", "PEAK", "TF_TARGET"),
    is.finite(assigned_topic_map[term_id]) &
      topic_num == assigned_topic_map[term_id],
    default = gammafit_candidate
  )]
  missing_assigned <- assignment_map[
    !paste(term_id, assigned_topic__, sep = "\r") %in%
      paste(collapsed$term_id, collapsed$topic_num, sep = "\r")
  ]
  if (nrow(missing_assigned)) {
    topic_position <- match(
      missing_assigned$assigned_topic__,
      .m3_opt_topic_ids(optimized_phi, "row")
    )
    term_position <- match(missing_assigned$term_id, colnames(optimized_phi))
    valid <- is.finite(topic_position) & is.finite(term_position)
    missing_assigned <- missing_assigned[valid]
    topic_position <- topic_position[valid]
    term_position <- term_position[valid]
    raw_members <- vapply(
      missing_assigned$assigned_topic__,
      function(topic) {
        paste(raw_topic_ids[raw_to_group == topic], collapse = ";")
      },
      character(1L)
    )
    recovered <- data.table::data.table(
      topic = paste0("Topic", missing_assigned$assigned_topic__),
      topic_num = missing_assigned$assigned_topic__,
      term_id = missing_assigned$term_id,
      term_group = missing_assigned$term_group,
      score = optimized_score[cbind(topic_position, term_position)],
      phi = optimized_phi[cbind(topic_position, term_position)],
      in_topic = TRUE,
      assignment_method = "gammafit_maxprob_optimized",
      gammafit_candidate = TRUE,
      raw_topic_members = raw_members,
      assignment_status = "assigned_after_topic_merge"
    )
    collapsed <- data.table::rbindlist(
      list(collapsed, recovered),
      use.names = TRUE,
      fill = TRUE
    )
  }
  collapsed[, `:=`(
    topic = paste0("Topic", topic_num),
    assignment_method = data.table::fifelse(
      term_group %in% c("GENE", "PEAK", "TF_TARGET"),
      "gammafit_maxprob_optimized",
      "gammafit"
    ),
    score = optimized_score,
    assignment_status = data.table::fifelse(
      in_topic,
      "assigned_optimized_topic",
      "not_assigned_optimized_topic"
    )
  )]
  collapsed[, optimized_score := NULL]
  data.table::setcolorder(collapsed, c(
    "topic", "topic_num", "term_id", "term_group", "score", "phi",
    "in_topic", "assignment_method", "gammafit_candidate",
    "raw_topic_members", "assignment_status"
  ))
  data.table::setorder(collapsed, topic_num, term_group, -in_topic, -score, term_id)
  collapsed[]
}

.m3_opt_jaccard <- function(group, item) {
  group_levels <- unique(as.character(group))
  group <- match(as.character(group), group_levels)
  if (!is.integer(item) || anyNA(item) || any(item < 1L)) {
    item <- match(item, unique(item))
  }
  incidence <- Matrix::sparseMatrix(
    i = group,
    j = item,
    x = 1,
    dims = c(max(group), max(item))
  )
  incidence@x[] <- 1
  intersection <- as.matrix(Matrix::tcrossprod(incidence))
  size <- Matrix::rowSums(incidence)
  union <- outer(size, size, "+") - intersection
  out <- intersection / pmax(union, 1)
  diag(out) <- 1
  rownames(out) <- group_levels
  colnames(out) <- group_levels
  out
}

.m3_opt_cross_jaccard <- function(group_a, group_b, item) {
  levels_a <- unique(group_a)
  levels_b <- unique(group_b)
  item_levels <- unique(item)
  a <- Matrix::sparseMatrix(
    i = match(group_a, levels_a),
    j = match(item, item_levels),
    x = 1,
    dims = c(length(levels_a), length(item_levels))
  )
  b <- Matrix::sparseMatrix(
    i = match(group_b, levels_b),
    j = match(item, item_levels),
    x = 1,
    dims = c(length(levels_b), length(item_levels))
  )
  a@x[] <- 1
  b@x[] <- 1
  intersection <- as.matrix(Matrix::tcrossprod(a, b))
  union <- outer(Matrix::rowSums(a), Matrix::rowSums(b), "+") - intersection
  out <- intersection / pmax(union, 1)
  rownames(out) <- levels_a
  colnames(out) <- levels_b
  out
}

.m3_opt_qc_tables <- function(assignments,
                              docs,
                              pairs,
                              target_levels = NULL,
                              raw_similarity,
                              optimized_similarity,
                              raw_topic_ids,
                              optimized_topic_ids) {
  x <- data.table::copy(assignments)
  tf_levels <- unique(docs$tf)
  tf_index <- match(docs$tf, tf_levels)
  x[, `:=`(
    tf_index = tf_index[doc_index]
  )]
  raw <- x[raw_aligned == TRUE]
  optimized <- x[optimized_aligned == TRUE]
  raw_counts <- unique(raw[, .(
    raw_topic = raw_target_topic,
    tf_index,
    target_index
  )])[, .(
    links = .N,
    genes = data.table::uniqueN(target_index),
    tfs = data.table::uniqueN(tf_index)
  ), by = raw_topic]
  optimized_counts <- unique(optimized[, .(
    optimized_topic,
    tf_index,
    target_index
  )])[, .(
    links = .N,
    genes = data.table::uniqueN(target_index),
    tfs = data.table::uniqueN(tf_index)
  ), by = optimized_topic]
  condition_topic <- unique(optimized[, .(
    condition_id,
    optimized_topic,
    tf_index,
    target_index
  )])[, .(
    links = .N,
    genes = data.table::uniqueN(target_index),
    tfs = data.table::uniqueN(tf_index)
  ), by = .(condition_id, optimized_topic)]

  raw_unique <- unique(raw[, .(
    condition_id,
    raw_topic = raw_target_topic,
    tf_index,
    target_index
  )])
  condition_similarity <- if (nrow(raw_unique)) {
    link_item <- as.integer(
      (as.double(raw_unique$tf_index) - 1) * nrow(pairs) +
        raw_unique$target_index
    )
    link_sim <- .m3_opt_jaccard(raw_unique$condition_id, link_item)
    gene_unique <- unique(raw_unique[, .(condition_id, target_index)])
    gene_sim <- .m3_opt_jaccard(gene_unique$condition_id, gene_unique$target_index)
    (link_sim + gene_sim) / 2
  } else {
    matrix(numeric(), 0, 0)
  }

  optimized_unique <- unique(optimized[, .(
    condition_id,
    optimized_topic,
    tf_index,
    target_index
  )])
  condition_topic_similarity <- if (nrow(optimized_unique)) {
    link_item <- as.integer(
      (as.double(optimized_unique$tf_index) - 1) * nrow(pairs) +
        optimized_unique$target_index
    )
    link_cross <- .m3_opt_cross_jaccard(
      optimized_unique$condition_id,
      optimized_unique$optimized_topic,
      link_item
    )
    gene_unique <- unique(optimized_unique[, .(
      condition_id,
      optimized_topic,
      target_index
    )])
    gene_cross <- .m3_opt_cross_jaccard(
      gene_unique$condition_id,
      gene_unique$optimized_topic,
      gene_unique$target_index
    )
    list(link = link_cross, gene = gene_cross, mean = (link_cross + gene_cross) / 2)
  } else {
    list(link = matrix(numeric(), 0, 0), gene = matrix(numeric(), 0, 0), mean = matrix(numeric(), 0, 0))
  }
  if (is.null(target_levels)) {
    target_levels <- as.character(pairs$target_gene)
  }
  list(
    assignments = x,
    raw_counts = raw_counts,
    optimized_counts = optimized_counts,
    condition_topic = condition_topic,
    raw_topic_similarity = raw_similarity,
    optimized_topic_similarity = optimized_similarity,
    condition_similarity = condition_similarity,
    condition_topic_similarity = condition_topic_similarity,
    raw_topic_ids = raw_topic_ids,
    optimized_topic_ids = optimized_topic_ids
    ,
    tf_levels = tf_levels,
    target_levels = target_levels
  )
}

.module3_optimize_condition_topics <- function(theta,
                                               phi,
                                               dtm,
                                               topic_terms,
                                               pair_assignment,
                                               assignment_mode = c(
                                                 "gene_peak",
                                                 "gene_only",
                                                 "tf_target"
                                               ),
                                               correspondence_assignment = NULL,
                                               require_theta_gate = TRUE,
                                               condition_gene_expression = NULL,
                                               condition_upregulated_genes = NULL,
                                               upregulation_reference_condition = NULL,
                                               upregulated_log2fc_min = 1,
                                               min_genes = 150L,
                                               min_links = 200L,
                                               similarity_threshold = 0.65,
                                               prefer_tf_theta_correspondence = TRUE,
                                               tf_topic_cutoff = 0.3,
                                               umap_max_links_per_condition = 10000L,
                                               seed = 20260716L,
                                               chunk_size = 50000L) {
  assignment_mode <- match.arg(assignment_mode)
  min_genes <- suppressWarnings(as.integer(min_genes)[[1L]])
  min_links <- suppressWarnings(as.integer(min_links)[[1L]])
  similarity_threshold <- suppressWarnings(as.numeric(similarity_threshold)[[1L]])
  tf_topic_cutoff <- suppressWarnings(as.numeric(tf_topic_cutoff)[[1L]])
  umap_max_links_per_condition <- suppressWarnings(
    as.integer(umap_max_links_per_condition)[[1L]]
  )
  chunk_size <- suppressWarnings(as.integer(chunk_size)[[1L]])
  if (!is.finite(min_genes) || min_genes < 1L) {
    .log_abort("min_genes must be a positive integer.")
  }
  if (!is.finite(min_links) || min_links < 1L) {
    .log_abort("min_links must be a positive integer.")
  }
  if (!is.finite(similarity_threshold) ||
      similarity_threshold <= 0 || similarity_threshold > 1) {
    .log_abort("similarity_threshold must be greater than 0 and at most 1.")
  }
  if (!is.finite(tf_topic_cutoff) || tf_topic_cutoff < 0 || tf_topic_cutoff > 1) {
    .log_abort("tf_topic_cutoff must be between 0 and 1.")
  }
  if (!is.finite(umap_max_links_per_condition) ||
      umap_max_links_per_condition < 1L) {
    .log_abort("umap_max_links_per_condition must be a positive integer.")
  }
  if (!is.finite(chunk_size) || chunk_size < 1L) {
    .log_abort("chunk_size must be a positive integer.")
  }
  theta <- .m3_opt_row_normalize(theta)
  phi <- .m3_opt_row_normalize(phi)
  raw_topic_ids <- .m3_opt_topic_ids(phi, "row")
  if (!identical(raw_topic_ids, .m3_opt_topic_ids(theta, "column"))) {
    .log_abort("Theta and phi topic IDs do not match for topic optimization.")
  }
  universe <- if (identical(assignment_mode, "tf_target")) {
    .m3_opt_link_universe_tf_target(dtm, pair_assignment)
  } else {
    .m3_opt_link_universe(dtm, pair_assignment)
  }
  raw_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    universe,
    chunk_size = chunk_size
  )
  raw_target_topic <- suppressWarnings(as.integer(universe$pairs$assigned_topic))
  pair_index <- if ("pair_index" %in% names(universe$links)) {
    universe$links$pair_index
  } else {
    universe$links$target_index
  }
  target_topic_for_link <- raw_target_topic[pair_index]
  raw_posterior_agrees <- is.finite(target_topic_for_link) &
    raw_posterior$topic == target_topic_for_link
  raw_target_position <- match(target_topic_for_link, raw_topic_ids)
  raw_theta_pass <- rep(FALSE, nrow(universe$links))
  valid_raw_theta <- is.finite(raw_target_position)
  raw_theta_pass[valid_raw_theta] <- theta[cbind(
    universe$links$doc_index[valid_raw_theta],
    raw_target_position[valid_raw_theta]
  )] >= tf_topic_cutoff
  raw_aligned <- if (isTRUE(require_theta_gate)) {
    raw_posterior_agrees & raw_theta_pass
  } else {
    is.finite(target_topic_for_link)
  }
  raw_links <- vapply(raw_topic_ids, function(topic) {
    data.table::uniqueN(pair_index[
      raw_aligned & target_topic_for_link == topic
    ])
  }, integer(1L))
  raw_genes <- vapply(raw_topic_ids, function(topic) {
    data.table::uniqueN(universe$links$target_index[
      raw_aligned & target_topic_for_link == topic
    ])
  }, integer(1L))
  gene_ids <- match(universe$pairs$gene_term_id, colnames(phi))
  peak_ids <- match(universe$pairs$peak_term_id, colnames(phi))
  raw_similarity <- .m3_opt_hellinger_similarity(phi, gene_ids, peak_ids)
  merge <- .m3_opt_merge_map(
    phi = phi,
    theta = theta,
    dtm = dtm,
    raw_topic_ids = raw_topic_ids,
    raw_links = raw_links,
    raw_genes = raw_genes,
    gene_ids = gene_ids,
    peak_ids = peak_ids,
    pair_assignment = correspondence_assignment %||% pair_assignment,
    min_genes = min_genes,
    min_links = min_links,
    similarity_threshold = similarity_threshold,
    prefer_tf_theta_correspondence = prefer_tf_theta_correspondence,
    tf_topic_cutoff = tf_topic_cutoff
  )
  optimized_phi <- .m3_opt_group_phi(phi, theta, dtm, raw_topic_ids, merge$mapping)
  optimized_theta <- .m3_opt_theta(theta, raw_topic_ids, merge$mapping)
  optimized_topic_ids <- .m3_opt_topic_ids(optimized_phi, "row")
  optimized_pairs <- .m3_opt_target_assignment(
    phi = phi,
    pair_assignment = universe$pairs,
    raw_topic_ids = raw_topic_ids,
    raw_to_group = merge$mapping
  )
  optimized_gene_pairs <- if (is.null(correspondence_assignment)) {
    NULL
  } else {
    .m3_opt_target_assignment(
      phi = phi,
      pair_assignment = correspondence_assignment,
      raw_topic_ids = raw_topic_ids,
      raw_to_group = merge$mapping
    )
  }
  optimized_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    universe,
    group_index = merge$mapping,
    chunk_size = chunk_size
  )
  optimized_target_topic <- optimized_pairs$optimized_assigned_topic[pair_index]
  optimized_posterior_agrees <- is.finite(optimized_target_topic) &
    optimized_posterior$topic == optimized_target_topic
  optimized_target_position <- match(optimized_target_topic, optimized_topic_ids)
  optimized_theta_pass <- rep(FALSE, nrow(universe$links))
  valid_optimized_theta <- is.finite(optimized_target_position)
  optimized_theta_pass[valid_optimized_theta] <- optimized_theta[cbind(
    universe$links$doc_index[valid_optimized_theta],
    optimized_target_position[valid_optimized_theta]
  )] >= tf_topic_cutoff
  optimized_aligned <- if (isTRUE(require_theta_gate)) {
    optimized_posterior_agrees & optimized_theta_pass
  } else {
    is.finite(optimized_target_topic)
  }
  eligible_sample_rows <- which(optimized_aligned)
  if (!length(eligible_sample_rows)) {
    .log_abort("No aligned links remain after topic optimization.")
  }
  sample_rows <- .m3_opt_balanced_sample(
    universe$links$condition_id,
    max_per_condition = umap_max_links_per_condition,
    seed = seed
  )
  sample_universe <- universe
  sample_universe$links <- data.table::copy(universe$links[sample_rows])
  sample_universe$links[, link_index := seq_len(.N)]
  raw_sample_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    sample_universe,
    sample_rows = seq_len(nrow(sample_universe$links)),
    chunk_size = chunk_size
  )
  optimized_sample_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    sample_universe,
    group_index = merge$mapping,
    sample_rows = seq_len(nrow(sample_universe$links)),
    chunk_size = chunk_size
  )
  assignments <- data.table::copy(universe$links)
  assignments[, `:=`(
    raw_target_topic = target_topic_for_link,
    raw_posterior_topic = raw_posterior$topic,
    raw_posterior_probability = raw_posterior$probability,
    raw_posterior_margin = raw_posterior$margin,
    raw_posterior_agrees = raw_posterior_agrees,
    raw_theta_pass = raw_theta_pass,
    raw_aligned = raw_aligned,
    optimized_topic = optimized_target_topic,
    optimized_posterior_topic = optimized_posterior$topic,
    optimized_posterior_probability = optimized_posterior$probability,
    optimized_posterior_margin = optimized_posterior$margin,
    optimized_posterior_agrees = optimized_posterior_agrees,
    optimized_theta_pass = optimized_theta_pass,
    optimized_aligned = optimized_aligned,
    recovered_after_merge = optimized_aligned & !raw_aligned
  )]
  optimized_score <- score_terms_normtop(
    optimized_phi,
    method = "normtop_specificity"
  )
  optimized_terms <- .m3_opt_topic_terms(
    raw_topic_terms = topic_terms,
    optimized_phi = optimized_phi,
    optimized_score = optimized_score,
    raw_topic_ids = raw_topic_ids,
    raw_to_group = merge$mapping,
    optimized_pairs = optimized_pairs,
    optimized_gene_pairs = optimized_gene_pairs
  )
  optimized_similarity <- .m3_opt_hellinger_similarity(
    optimized_phi,
    gene_ids,
    peak_ids
  )
  raw_correspondence <- .m3_opt_tf_theta_correspondence(
    theta = theta,
    phi = phi,
    pair_assignment = correspondence_assignment %||% pair_assignment,
    raw_topic_ids = raw_topic_ids,
    raw_to_group = raw_topic_ids,
    cutoff = tf_topic_cutoff
  )
  optimized_correspondence <- .m3_opt_tf_theta_correspondence(
    theta = theta,
    phi = phi,
    pair_assignment = correspondence_assignment %||% pair_assignment,
    raw_topic_ids = raw_topic_ids,
    raw_to_group = merge$mapping,
    cutoff = tf_topic_cutoff
  )
  qc <- .m3_opt_qc_tables(
    assignments = assignments,
    docs = universe$docs,
    pairs = optimized_pairs,
    target_levels = universe$target_levels %||% NULL,
    raw_similarity = raw_similarity,
    optimized_similarity = optimized_similarity,
    raw_topic_ids = raw_topic_ids,
    optimized_topic_ids = optimized_topic_ids
  )
  list(
    document_design = "condition_tf",
    theta = optimized_theta,
    phi = optimized_phi,
    score = optimized_score,
    topic_terms = optimized_terms,
    pair_assignment = optimized_pairs,
    gene_assignment = optimized_gene_pairs,
    raw_theta = theta,
    raw_phi = phi,
    raw_topic_terms = topic_terms,
    raw_pair_assignment = pair_assignment,
    raw_correspondence_assignment = correspondence_assignment,
    assignment_mode = assignment_mode,
    require_theta_gate = isTRUE(require_theta_gate),
    condition_gene_expression = condition_gene_expression,
    condition_upregulated_genes = condition_upregulated_genes,
    upregulation_reference_condition = upregulation_reference_condition,
    upregulated_log2fc_min = upregulated_log2fc_min,
    raw_to_optimized = stats::setNames(merge$mapping, raw_topic_ids),
    merge_audit = merge$audit,
    raw_tf_theta_correspondence = raw_correspondence,
    optimized_tf_theta_correspondence = optimized_correspondence,
    sample_rows = sample_rows,
    raw_sample_probability = raw_sample_posterior$sample_probability,
    optimized_sample_probability = optimized_sample_posterior$sample_probability,
    tf_topic_cutoff = tf_topic_cutoff,
    min_genes = min_genes,
    min_links = min_links,
    similarity_threshold = similarity_threshold,
    prefer_tf_theta_correspondence = isTRUE(
      prefer_tf_theta_correspondence
    ),
    umap_max_links_per_condition = umap_max_links_per_condition,
    qc_seed = as.integer(seed),
    qc = qc
  )
}

.m3_qc_theme <- function(base_size = 9) {
  ggplot2::theme_minimal(base_family = "Helvetica", base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        colour = "#20272E"
      ),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, face = "bold", colour = "#4D5963"),
      plot.caption = ggplot2::element_text(size = 9, face = "bold", colour = "#4D5963"),
      axis.title = ggplot2::element_text(size = 9, face = "bold"),
      axis.text = ggplot2::element_text(size = 9, face = "bold"),
      legend.title = ggplot2::element_text(size = 9, face = "bold"),
      legend.text = ggplot2::element_text(size = 9, face = "bold"),
      panel.grid = ggplot2::element_line(colour = "#D5DDE2", linewidth = 0.25),
      plot.margin = ggplot2::margin(5, 7, 5, 7)
    )
}

.m3_qc_short_condition_labels <- function(x) {
  original <- as.character(x)
  unique_original <- unique(original)
  prefixes <- sub("_.*$", "", unique_original)
  strip_prefix <- length(unique_original) > 1L &&
    length(unique(prefixes)) == 1L
  out <- if (strip_prefix) sub("^[^_]+_", "", original) else original
  out <- sub("_Ctrl$", "", out)
  out <- sub("_TGFb$", " TGFb", out)
  out <- gsub("_", " ", out, fixed = TRUE)
  mapping <- unique(data.table::data.table(original, short = out))
  duplicated_label <- duplicated(mapping$short) |
    duplicated(mapping$short, fromLast = TRUE)
  mapping[duplicated_label, short := gsub("_", " ", original, fixed = TRUE)]
  mapping$short[match(original, mapping$original)]
}

.m3_qc_compact_count <- function(x) {
  x <- as.numeric(x)
  out <- format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  thousands <- x >= 1000 & x < 1e6
  millions <- x >= 1e6
  out[thousands] <- sprintf("%.1fk", x[thousands] / 1000)
  out[millions] <- sprintf("%.1fM", x[millions] / 1e6)
  sub("[.]0([kM])$", "\\1", out)
}

.m3_qc_heatmap_count <- function(x) {
  x <- as.numeric(x)
  out <- format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  thousands <- x >= 1000 & x < 1e6
  millions <- x >= 1e6
  out[thousands] <- paste0(round(x[thousands] / 1000), "k")
  out[millions] <- paste0(round(x[millions] / 1e6), "M")
  out
}

.m3_qc_cluster_order <- function(x, margin = c("row", "column")) {
  margin <- match.arg(margin)
  x <- as.matrix(x)
  if (identical(margin, "column")) x <- t(x)
  if (nrow(x) < 2L) return(seq_len(nrow(x)))
  work <- sqrt(pmax(x, 0))
  means <- rowMeans(work)
  sds <- apply(work, 1L, stats::sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  work <- sweep(sweep(work, 1L, means, "-"), 1L, sds, "/")
  work[!is.finite(work)] <- 0
  stats::hclust(stats::dist(work), method = "ward.D2")$order
}

.m3_qc_matrix_long <- function(x, row_order = NULL, column_order = NULL) {
  x <- as.matrix(x)
  if (is.null(rownames(x))) rownames(x) <- seq_len(nrow(x))
  if (is.null(colnames(x))) colnames(x) <- seq_len(ncol(x))
  if (is.null(row_order)) row_order <- seq_len(nrow(x))
  if (is.null(column_order)) column_order <- seq_len(ncol(x))
  out <- data.table::as.data.table(as.table(x))
  data.table::setnames(out, c("row_label", "column_label", "value"))
  out[, `:=`(
    row_label = factor(
      as.character(row_label),
      levels = rev(rownames(x)[row_order])
    ),
    column_label = factor(
      as.character(column_label),
      levels = colnames(x)[column_order]
    ),
    value = as.numeric(value)
  )]
  out
}

.m3_qc_heatmap_layout_spec <- function(n_conditions, n_topics) {
  n_conditions <- max(1L, as.integer(n_conditions))
  n_topics <- max(1L, as.integer(n_topics))
  axis_size <- max(
    9,
    min(11, 165 / n_conditions, 360 / n_topics)
  )
  list(
    axis_size = axis_size,
    title_size = max(9, min(11, axis_size + 1)),
    label_size = max(3.2, min(3.8, axis_size / 2.845)),
    tile_linewidth = if (n_conditions * n_topics >= 600L) 0.14 else 0.20,
    plot_margin = ggplot2::margin(2, 4, 2, 4)
  )
}

.m3_qc_heatmap_fill_max <- function(x, probability = 0.95) {
  values <- as.numeric(x)
  values <- values[is.finite(values) & values > 0]
  if (!length(values)) return(1)
  out <- as.numeric(stats::quantile(
    values,
    probs = probability,
    names = FALSE,
    type = 8
  ))
  if (!is.finite(out) || out <= 0) out <- max(values)
  out
}

.m3_qc_similarity_plot <- function(x,
                                   title,
                                   subtitle = NULL,
                                   row_order = NULL,
                                   column_order = NULL,
                                   autoscale = FALSE) {
  if (is.null(row_order)) row_order <- .m3_qc_cluster_order(x, "row")
  if (is.null(column_order)) {
    column_order <- .m3_qc_cluster_order(x, "column")
  }
  layout_spec <- .m3_qc_heatmap_layout_spec(nrow(x), ncol(x))
  long <- .m3_qc_matrix_long(x, row_order, column_order)
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(colour = "#D6DCE1", linewidth = 0.2)
  if (isTRUE(autoscale)) {
    values <- long$value[is.finite(long$value)]
    limits <- range(values)
    if (length(values) == 0L || diff(limits) <= 0) {
      limits <- c(0, max(c(values, 1), na.rm = TRUE))
    }
    p <- p + ggplot2::scale_fill_viridis_c(
      option = "C",
      limits = limits,
      oob = scales::squish,
      name = "Mean Jaccard"
    )
  } else {
    p <- p + ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Similarity"
    )
  }
  p +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    .m3_qc_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = layout_spec$axis_size),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = layout_spec$plot_margin
    )
}

.m3_qc_condition_correlation <- function(link_matrix, gene_matrix) {
  link_matrix <- as.matrix(link_matrix)
  gene_matrix <- as.matrix(gene_matrix)
  if (!identical(dim(link_matrix), dim(gene_matrix)) ||
      !identical(rownames(link_matrix), rownames(gene_matrix))) {
    .log_abort("Condition link and gene topic matrices must have matching rows.")
  }
  profiles <- cbind(
    .m3_opt_row_normalize(link_matrix),
    .m3_opt_row_normalize(gene_matrix)
  )
  if (nrow(profiles) == 1L) {
    out <- matrix(
      1,
      nrow = 1L,
      dimnames = list(rownames(profiles), rownames(profiles))
    )
    return(out)
  }
  out <- suppressWarnings(stats::cor(
    t(profiles),
    use = "pairwise.complete.obs"
  ))
  out[!is.finite(out)] <- 0
  diag(out) <- 1
  out
}

.m3_qc_theta_condition_correlation <- function(theta) {
  theta <- as.matrix(theta)
  if (!nrow(theta) || !ncol(theta)) {
    .log_abort("Condition theta must contain at least one row and one topic.")
  }
  if (is.null(rownames(theta)) || anyNA(rownames(theta)) ||
      any(!nzchar(rownames(theta))) || anyDuplicated(rownames(theta))) {
    .log_abort("Condition theta must have unique, non-empty condition row names.")
  }
  theta <- .m3_opt_row_normalize(theta)
  if (nrow(theta) == 1L) {
    return(matrix(
      1,
      nrow = 1L,
      dimnames = list(rownames(theta), rownames(theta))
    ))
  }
  out <- suppressWarnings(stats::cor(
    t(theta),
    use = "pairwise.complete.obs"
  ))
  out[!is.finite(out)] <- 0
  diag(out) <- 1
  out
}

.m3_qc_correlation_plot <- function(x,
                                    title,
                                    subtitle = NULL,
                                    order = NULL) {
  x <- as.matrix(x)
  if (is.null(order)) order <- .m3_qc_cluster_order(x, "row")
  layout_spec <- .m3_qc_heatmap_layout_spec(nrow(x), ncol(x))
  long <- .m3_qc_matrix_long(x, order, order)
  ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(colour = "#D6DCE1", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    .m3_qc_theme() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = layout_spec$axis_size),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = layout_spec$plot_margin
    )
}

.m3_qc_count_heatmap <- function(x,
                                 title,
                                 row_order = NULL,
                                 column_order = NULL,
                                 label_min = 1L,
                                 fill_max = NULL,
                                 square_cells = TRUE,
                                 show_legend = TRUE,
                                 show_labels = TRUE,
                                 label_size = NULL,
                                 show_x_axis = TRUE) {
  layout_spec <- .m3_qc_heatmap_layout_spec(nrow(x), ncol(x))
  if (is.null(label_size)) label_size <- layout_spec$label_size
  if (is.null(row_order)) row_order <- .m3_qc_cluster_order(x, "row")
  if (is.null(column_order)) column_order <- .m3_qc_cluster_order(x, "column")
  long <- .m3_qc_matrix_long(x, row_order, column_order)
  if (is.null(fill_max)) fill_max <- .m3_qc_heatmap_fill_max(long$value)
  if (!is.finite(fill_max) || fill_max <= 0) fill_max <- 1
  long[, label := ifelse(value >= label_min, .m3_qc_heatmap_count(value), "")]
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(
      colour = "#D6DCE1",
      linewidth = layout_spec$tile_linewidth
    )
  if (isTRUE(show_labels)) {
    p <- p + ggplot2::geom_text(
        ggplot2::aes(label = label, colour = value >= fill_max * 0.55),
        family = "Helvetica",
        fontface = "bold",
        size = label_size,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = c(`TRUE` = "white", `FALSE` = "#20272E")
      )
  }
  p <- p +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      trans = "sqrt",
      limits = c(0, fill_max),
      oob = scales::squish,
      labels = scales::label_comma(),
      name = "Count"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    .m3_qc_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = layout_spec$title_size),
      axis.text = ggplot2::element_text(size = layout_spec$axis_size),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = if (isTRUE(show_legend)) "right" else "none",
      plot.margin = layout_spec$plot_margin
    )
  if (!isTRUE(show_x_axis)) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  if (isTRUE(square_cells)) p <- p + ggplot2::coord_fixed(ratio = 1)
  p
}

.m3_qc_value_heatmap <- function(x,
                                 title,
                                 legend_title,
                                 row_order = NULL,
                                 column_order = NULL,
                                 limits = NULL,
                                 color_transform = "sqrt",
                                 show_x_axis = TRUE) {
  color_transform <- match.arg(color_transform, c("identity", "sqrt"))
  layout_spec <- .m3_qc_heatmap_layout_spec(nrow(x), ncol(x))
  if (is.null(row_order)) row_order <- .m3_qc_cluster_order(x, "row")
  if (is.null(column_order)) column_order <- .m3_qc_cluster_order(x, "column")
  long <- .m3_qc_matrix_long(x, row_order, column_order)
  if (is.null(limits)) {
    limits <- c(0, .m3_qc_heatmap_fill_max(long$value))
  }
  if (length(limits) != 2L || !all(is.finite(limits)) ||
      limits[[2L]] <= limits[[1L]]) {
    limits <- c(0, 1)
  }
  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(
      colour = "#D6DCE1",
      linewidth = layout_spec$tile_linewidth
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      trans = color_transform,
      limits = limits,
      oob = scales::squish,
      labels = scales::label_number(accuracy = 0.01),
      name = legend_title
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    .m3_qc_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = layout_spec$title_size),
      axis.text = ggplot2::element_text(size = layout_spec$axis_size),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = layout_spec$plot_margin
    )
  if (!isTRUE(show_x_axis)) {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  p
}

.m3_qc_topic_count_bar <- function(counts,
                                   topic_column,
                                   value_column,
                                   title,
                                   topic_order,
                                   fill) {
  x <- data.table::copy(counts)
  data.table::setnames(x, topic_column, "topic")
  data.table::setnames(x, value_column, "count")
  x[, topic := factor(
    paste0("Topic ", topic),
    levels = rev(topic_order)
  )]
  ggplot2::ggplot(x, ggplot2::aes(count, topic)) +
    ggplot2::geom_col(
      width = 0.82,
      fill = fill,
      colour = "#30383F",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .m3_qc_compact_count(count)),
      hjust = 1.06,
      colour = "white",
      family = "Helvetica",
      fontface = "bold",
      size = 3.2
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(title = title, x = "Full-universe count", y = NULL) +
    .m3_qc_theme() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(colour = NA),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(subtitle = " ")
}

.m3_qc_topic_structure_plot <- function(similarity,
                                        counts,
                                        topic_column,
                                        topic_order,
                                        title,
                                        subtitle,
                                        tf_title = "TFs",
                                        sidebar_columns = c(
                                          "links",
                                          "genes",
                                          "tfs"
                                        ),
                                        sidebar_labels = c(
                                          "Links",
                                          "Target genes",
                                          tf_title
                                        ),
                                        sidebar_colors = c(
                                          "#007C78",
                                          "#D97824",
                                          "#6A5ACD"
                                        )) {
  similarity <- as.matrix(similarity)
  topic_ids <- suppressWarnings(as.integer(sub(
    "^Topic ",
    "",
    rownames(similarity)
  )))
  ordered_ids <- topic_ids[topic_order]
  source_labels <- paste0("Topic ", ordered_ids)
  ordered_labels <- paste0("T", ordered_ids)
  k <- length(ordered_ids)
  reordered <- similarity[topic_order, topic_order, drop = FALSE]
  heatmap <- data.table::as.data.table(as.table(reordered))
  data.table::setnames(heatmap, c("row_label", "column_label", "similarity"))
  heatmap[, `:=`(
    x = match(as.character(column_label), source_labels),
    y = k - match(as.character(row_label), source_labels) + 1L,
    similarity = as.numeric(similarity)
  )]
  heatmap[, similarity := pmin(1, pmax(0, similarity))]

  count_data <- data.table::copy(counts)
  data.table::setnames(count_data, topic_column, "topic")
  count_data[, y := k - match(as.integer(topic), ordered_ids) + 1L]
  count_data <- count_data[is.finite(y)]
  if (!length(sidebar_columns) ||
      length(sidebar_columns) != length(sidebar_labels) ||
      length(sidebar_columns) != length(sidebar_colors)) {
    .log_abort("Topic QC sidebar columns, labels, and colors must align.")
  }
  missing_sidebar <- setdiff(sidebar_columns, names(count_data))
  for (column in missing_sidebar) count_data[, (column) := 0]
  bar_width <- max(3.5, min(6, k * 0.20))
  gap <- max(1, k * 0.04)
  sidebar_starts <- k + gap +
    (seq_along(sidebar_columns) - 1L) * (bar_width + gap)
  sidebar_data <- lapply(seq_along(sidebar_columns), function(i) {
    values <- suppressWarnings(as.numeric(count_data[[sidebar_columns[[i]]]]))
    values[!is.finite(values) | values < 0] <- 0
    value_max <- max(values, na.rm = TRUE)
    if (!is.finite(value_max) || value_max <= 0) value_max <- 1
    data.table::data.table(
      y = count_data$y,
      value = values,
      start = sidebar_starts[[i]],
      xmax = sidebar_starts[[i]] + values / value_max * bar_width,
      label = .m3_qc_compact_count(values),
      inside = values / value_max >= 0.28,
      color = sidebar_colors[[i]],
      panel = i
    )
  })
  right_edge <- utils::tail(sidebar_starts, 1L) + bar_width + 1.1
  title_y <- k + 1.15

  plot <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = heatmap,
      ggplot2::aes(x, y, fill = similarity),
      colour = "#D6DCE1",
      linewidth = 0.2,
      width = 1,
      height = 1
    )
  for (i in seq_along(sidebar_data)) {
    panel_data <- sidebar_data[[i]]
    plot <- plot +
      ggplot2::geom_rect(
        data = panel_data,
        ggplot2::aes(
          xmin = start,
          xmax = xmax,
          ymin = y - 0.40,
          ymax = y + 0.40
        ),
        fill = sidebar_colors[[i]],
        colour = "#30383F",
        linewidth = 0.2
      ) +
      ggplot2::geom_text(
        data = panel_data[inside == TRUE],
        ggplot2::aes(xmax - 0.12, y, label = label),
        hjust = 1,
        colour = "white",
        family = "Helvetica",
        fontface = "bold",
        size = 3.0
      ) +
      ggplot2::geom_text(
        data = panel_data[inside == FALSE],
        ggplot2::aes(xmax + 0.12, y, label = label),
        hjust = 0,
        colour = "#20272E",
        family = "Helvetica",
        fontface = "bold",
        size = 3.0
      ) +
      ggplot2::annotate(
        "text",
        x = sidebar_starts[[i]] + bar_width / 2,
        y = title_y,
        label = sidebar_labels[[i]],
        family = "Helvetica",
        fontface = "bold",
        size = 3.5
      )
  }
  plot +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Similarity"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(k),
      labels = ordered_labels,
      limits = c(0.5, right_edge + 0.25),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(k),
      labels = rev(ordered_labels),
      limits = c(0.5, title_y + 0.35),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL
    ) +
    .m3_qc_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.0, "in")
    )
}

.m3_qc_umap_coordinates <- function(probability, seed = 20260716L) {
  probability <- .m3_opt_row_normalize(probability)
  n <- nrow(probability)
  if (n < 3L) return(cbind(UMAP1 = seq_len(n), UMAP2 = 0))
  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(physical) || physical < 1L) physical <- 1L
  uwot::umap(
    sqrt(probability),
    metric = "euclidean",
    n_neighbors = min(30L, n - 1L),
    min_dist = 0.15,
    n_components = 2L,
    n_threads = min(20L, as.integer(physical)),
    n_sgd_threads = 1L,
    seed = as.integer(seed),
    verbose = FALSE,
    ret_model = FALSE
  )
}

.m3_qc_feature_umap_coordinates <- function(features, seed = 20260716L) {
  features <- as.matrix(features)
  storage.mode(features) <- "double"
  features[!is.finite(features)] <- 0
  n <- nrow(features)
  if (n < 3L) return(cbind(UMAP1 = seq_len(n), UMAP2 = 0))
  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (!is.finite(physical) || physical < 1L) physical <- 1L
  uwot::umap(
    features,
    metric = "euclidean",
    n_neighbors = min(30L, n - 1L),
    min_dist = 0.15,
    n_components = 2L,
    n_threads = min(20L, as.integer(physical)),
    n_sgd_threads = 1L,
    seed = as.integer(seed),
    verbose = FALSE,
    ret_model = FALSE
  )
}

.m3_qc_jitter_umap <- function(data,
                               seed = 20260718L,
                               fraction = 0.005) {
  out <- data.table::copy(data)
  if (!nrow(out) || !all(c("UMAP1", "UMAP2") %in% names(out))) return(out)
  x_span <- diff(range(out$UMAP1, finite = TRUE))
  y_span <- diff(range(out$UMAP2, finite = TRUE))
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  set.seed(as.integer(seed))
  out[, `:=`(
    UMAP1 = UMAP1 + stats::runif(.N, -x_span * fraction, x_span * fraction),
    UMAP2 = UMAP2 + stats::runif(.N, -y_span * fraction, y_span * fraction)
  )]
  out[]
}

.m3_qc_gene_term_umap_data <- function(optimization,
                                       gene_term_assignment = NULL,
                                       gene_universe = NULL,
                                       seed = 20260716L) {
  phi <- optimization$raw_phi
  if (is.null(phi) || !is.matrix(phi) || !nrow(phi) || !ncol(phi)) {
    .log_abort("Raw phi is required for the Gene-term assignment UMAP.")
  }
  assignment <- gene_term_assignment
  if (is.null(assignment) || !nrow(assignment)) {
    assignment <- optimization$raw_topic_terms
  }
  if (is.null(assignment) || !nrow(assignment)) {
    assignment <- optimization$raw_correspondence_assignment
  }
  if (is.null(assignment) || !nrow(assignment)) {
    assignment <- optimization$raw_pair_assignment
  }
  assignment <- data.table::as.data.table(assignment)
  if ("term_group" %in% names(assignment)) {
    assignment <- assignment[term_group == "GENE"]
  }
  if (!"target_gene" %in% names(assignment) &&
      "term_id" %in% names(assignment)) {
    assignment[, target_gene := sub("^GENE:", "", as.character(term_id))]
  }
  if (!"assigned_topic" %in% names(assignment)) {
    topic_column <- intersect(c("topic_num", "topic"), names(assignment))
    if (length(topic_column)) {
      assignment[, assigned_topic := get(topic_column[[1L]])]
    }
  }
  if (!"assigned" %in% names(assignment)) {
    if ("in_topic" %in% names(assignment)) {
      assignment[, assigned := .as_logical_flag(in_topic)]
    } else {
      assignment[, assigned := TRUE]
    }
  }
  .assert_has_cols(
    assignment,
    c("target_gene", "assigned_topic", "assigned"),
    context = "Gene-term assignment UMAP"
  )
  assignment[, assigned_topic := suppressWarnings(
    as.integer(as.character(assigned_topic))
  )]
  assignment <- unique(assignment[
    .as_logical_flag(assigned) & is.finite(assigned_topic),
    .(
      target_gene = as.character(target_gene),
      topic_num = as.integer(assigned_topic)
    )
  ])
  conflicting <- assignment[, .(n_topics = data.table::uniqueN(topic_num)),
                            by = target_gene][n_topics > 1L]
  if (nrow(conflicting)) {
    .log_abort(
      "{nrow(conflicting)} assigned Gene terms have multiple primary topics."
    )
  }
  assignment <- unique(assignment, by = "target_gene")
  raw_gene_term_ids <- colnames(phi)[startsWith(colnames(phi), "GENE:")]
  raw_gene_names <- sub("^GENE:", "", raw_gene_term_ids)
  if (!length(raw_gene_names) || anyDuplicated(raw_gene_names)) {
    .log_abort("Raw phi must contain unique Gene terms for the UMAP.")
  }
  if (is.null(gene_universe)) {
    gene_names <- raw_gene_names
  } else {
    gene_names <- unique(sub(
      "^GENE:",
      "",
      trimws(as.character(gene_universe))
    ))
    gene_names <- gene_names[!is.na(gene_names) & nzchar(gene_names)]
    missing_universe <- setdiff(gene_names, raw_gene_names)
    if (length(missing_universe)) {
      .log_warn(
        "Ignoring {length(missing_universe)} requested Gene terms absent from raw phi."
      )
    }
    gene_names <- gene_names[gene_names %in% raw_gene_names]
  }
  if (!length(gene_names)) {
    .log_abort("No Gene terms remain in the requested UMAP universe.")
  }
  gene_indices <- match(gene_names, raw_gene_names)
  phi_gene_indices <- match(raw_gene_term_ids[gene_indices], colnames(phi))
  missing_terms <- !assignment$target_gene %in% raw_gene_names
  if (any(missing_terms)) {
    .log_warn(
      "Ignoring {sum(missing_terms)} assigned Gene terms absent from raw phi."
    )
    assignment <- assignment[!missing_terms]
  }
  gene_data <- data.table::data.table(target_gene = gene_names)
  gene_data <- assignment[gene_data, on = "target_gene"]
  gene_data[, assigned := is.finite(topic_num)]
  embedding_source <- "phi"
  gene_features <- optimization$gene_umap_features
  if (!is.null(gene_features)) {
    gene_features <- as.matrix(gene_features)
    if (!nrow(gene_features) || !ncol(gene_features) ||
        is.null(rownames(gene_features)) ||
        anyNA(rownames(gene_features)) ||
        any(!nzchar(rownames(gene_features))) ||
        anyDuplicated(rownames(gene_features))) {
      .log_abort(
        "Gene UMAP features must have unique, non-empty Gene row names."
      )
    }
    feature_names <- sub("^GENE:", "", rownames(gene_features))
    if (anyDuplicated(feature_names)) {
      .log_abort("Gene UMAP feature names must remain unique after normalization.")
    }
    feature_indices <- match(gene_data$target_gene, feature_names)
    missing_features <- is.na(feature_indices)
    if (any(missing_features)) {
      .log_warn(
        "Ignoring {sum(missing_features)} Gene terms absent from Gene UMAP features."
      )
      gene_data <- gene_data[!missing_features]
      phi_gene_indices <- phi_gene_indices[!missing_features]
      feature_indices <- feature_indices[!missing_features]
    }
    if (!nrow(gene_data)) {
      .log_abort("No Gene terms match the Gene UMAP features.")
    }
    supplied_coordinates <- optimization$gene_umap_coordinates
    if (!is.null(supplied_coordinates)) {
      supplied_coordinates <- as.matrix(supplied_coordinates)
      if (ncol(supplied_coordinates) < 2L ||
          is.null(rownames(supplied_coordinates)) ||
          anyNA(rownames(supplied_coordinates)) ||
          any(!nzchar(rownames(supplied_coordinates))) ||
          anyDuplicated(rownames(supplied_coordinates))) {
        .log_abort(
          "Gene UMAP coordinates must have two columns and unique Gene row names."
        )
      }
      coordinate_names <- sub("^GENE:", "", rownames(supplied_coordinates))
      if (anyDuplicated(coordinate_names)) {
        .log_abort(
          "Gene UMAP coordinate names must remain unique after normalization."
        )
      }
      coordinate_indices <- match(gene_data$target_gene, coordinate_names)
      if (anyNA(coordinate_indices)) {
        .log_abort(
          "Gene UMAP coordinates are missing requested Gene terms."
        )
      }
      coordinates <- supplied_coordinates[
        coordinate_indices,
        seq_len(2L),
        drop = FALSE
      ]
      storage.mode(coordinates) <- "double"
      if (any(!is.finite(coordinates))) {
        .log_abort("Gene UMAP coordinates must be finite.")
      }
    } else {
      coordinates <- .m3_qc_feature_umap_coordinates(
        gene_features[feature_indices, , drop = FALSE],
        seed = seed
      )
    }
    embedding_source <- as.character(
      optimization$gene_umap_feature_label %||% "condition profiles"
    )[[1L]]
  } else {
    probability <- t(phi[, phi_gene_indices, drop = FALSE])
    probability <- .m3_opt_row_normalize(probability)
    coordinates <- .m3_qc_umap_coordinates(probability, seed = seed)
  }
  topic_levels <- paste0(
    "Topic ",
    as.integer(sub("^Topic", "", rownames(phi)))
  )
  data.table::data.table(
    target_gene = gene_data$target_gene,
    topic_num = gene_data$topic_num,
    assigned = gene_data$assigned,
    topic = factor(
      data.table::fifelse(
        gene_data$assigned,
        paste0("Topic ", gene_data$topic_num),
        NA_character_
      ),
      levels = topic_levels
    ),
    embedding_source = embedding_source,
    UMAP1 = as.numeric(coordinates[, 1L]),
    UMAP2 = as.numeric(coordinates[, 2L])
  )
}

.m3_qc_peak_term_umap_data <- function(optimization,
                                       peak_term_assignment = NULL,
                                       top_n = 10000L,
                                       seed = 20260716L) {
  phi <- optimization$raw_phi
  if (is.null(phi) || !is.matrix(phi) || !nrow(phi) || !ncol(phi)) {
    .log_abort("Raw phi is required for the Peak-term assignment UMAP.")
  }
  top_n <- suppressWarnings(as.numeric(top_n[[1L]]))
  unlimited <- is.infinite(top_n) && top_n > 0
  finite_integer <- is.finite(top_n) &&
    top_n >= 1 &&
    top_n <= .Machine$integer.max &&
    top_n == floor(top_n)
  if (!unlimited && !finite_integer) {
    .log_abort("Peak-term UMAP top_n must be a positive integer or Inf.")
  }
  top_n <- if (unlimited) .Machine$integer.max else as.integer(top_n)
  peak_term_ids <- colnames(phi)[startsWith(colnames(phi), "PEAK:")]
  peak_ids <- sub("^PEAK:", "", peak_term_ids)
  if (!length(peak_ids)) {
    .log_abort("No Peak terms are available in raw phi for the UMAP.")
  }
  if (anyNA(peak_ids) || any(!nzchar(peak_ids)) || anyDuplicated(peak_ids)) {
    .log_abort("Raw phi must contain unique, non-empty Peak terms for the UMAP.")
  }
  assignment <- peak_term_assignment
  if (is.null(assignment) || !nrow(assignment)) {
    assignment <- optimization$raw_topic_terms
  }
  if (is.null(assignment) || !nrow(assignment)) {
    assignment <- optimization$raw_pair_assignment
  }
  assignment <- data.table::as.data.table(assignment)
  if ("term_group" %in% names(assignment)) {
    assignment <- assignment[term_group == "PEAK"]
  }
  if (!"peak_id" %in% names(assignment)) {
    peak_column <- intersect(c("peak_term_id", "term_id"), names(assignment))
    if (length(peak_column)) {
      assignment[, peak_id := sub(
        "^PEAK:",
        "",
        as.character(get(peak_column[[1L]]))
      )]
    }
  }
  if (!"assigned_topic" %in% names(assignment)) {
    topic_column <- intersect(c("topic_num", "topic"), names(assignment))
    if (length(topic_column)) {
      assignment[, assigned_topic := get(topic_column[[1L]])]
    }
  }
  if (!"assigned" %in% names(assignment)) {
    if ("in_topic" %in% names(assignment)) {
      assignment[, assigned := .as_logical_flag(in_topic)]
    } else {
      assignment[, assigned := TRUE]
    }
  }
  .assert_has_cols(
    assignment,
    c("peak_id", "assigned_topic", "assigned"),
    context = "Peak-term assignment UMAP"
  )
  assignment[, assigned_topic := suppressWarnings(
    as.integer(as.character(assigned_topic))
  )]
  assignment <- unique(assignment[
    .as_logical_flag(assigned) & is.finite(assigned_topic),
    .(
      peak_id = as.character(peak_id),
      topic_num = as.integer(assigned_topic)
    )
  ])
  conflicting <- assignment[, .(n_topics = data.table::uniqueN(topic_num)),
                            by = peak_id][n_topics > 1L]
  if (nrow(conflicting)) {
    .log_abort(
      "{nrow(conflicting)} assigned Peak terms have multiple primary topics."
    )
  }
  assignment <- unique(assignment, by = "peak_id")
  peak_columns <- match(peak_term_ids, colnames(phi))
  probability <- phi[, peak_columns, drop = FALSE]
  probability[!is.finite(probability) | probability < 0] <- 0
  totals <- colSums(probability)
  totals[!is.finite(totals) | totals <= 0] <- 1
  probability <- sweep(probability, 2L, totals, "/")
  probability[!is.finite(probability) | probability < 0] <- 0
  topic_ids <- .m3_opt_topic_ids(phi, "row")
  topic_position <- max.col(t(probability), ties.method = "first")
  peak_data <- data.table::data.table(
    peak_id = peak_ids,
    peak_term_id = peak_term_ids,
    max_phi_topic = as.integer(topic_ids[topic_position]),
    max_probability = as.numeric(probability[cbind(
      topic_position,
      seq_along(topic_position)
    )])
  )
  data.table::setorder(peak_data, -max_probability, peak_id)
  peak_data <- peak_data[seq_len(min(top_n, .N))]
  peak_data[, selection_rank := seq_len(.N)]
  selected_columns <- match(peak_data$peak_term_id, peak_term_ids)
  coordinates <- .m3_qc_umap_coordinates(
    t(probability[, selected_columns, drop = FALSE]),
    seed = seed
  )
  peak_data <- assignment[peak_data, on = "peak_id"]
  topic_levels <- paste0("Topic ", topic_ids)
  peak_data[, `:=`(
    assigned = is.finite(topic_num),
    topic = factor(
      data.table::fifelse(
        is.finite(topic_num),
        paste0("Topic ", topic_num),
        NA_character_
      ),
      levels = topic_levels
    ),
    UMAP1 = as.numeric(coordinates[, 1L]),
    UMAP2 = as.numeric(coordinates[, 2L])
  )]
  peak_data[]
}

.m3_qc_primary_term_assignments <- function(assignment,
                                             context = "topic assignment") {
  assignment <- data.table::as.data.table(data.table::copy(assignment))
  if (!nrow(assignment)) {
    .log_abort("{context} must contain assigned Gene and Peak terms.")
  }
  if (!"term_id" %in% names(assignment)) {
    .log_abort("{context} must contain term_id.")
  }
  if (!"term_group" %in% names(assignment)) {
    assignment[, term_group := data.table::fcase(
      startsWith(term_id, "GENE:"), "GENE",
      startsWith(term_id, "PEAK:"), "PEAK",
      default = NA_character_
    )]
  }
  topic_column <- intersect(c("assigned_topic", "topic_num", "topic"), names(assignment))
  if (!length(topic_column)) {
    .log_abort("{context} must contain an assigned topic column.")
  }
  assignment[, primary_topic := suppressWarnings(as.integer(
    as.character(get(topic_column[[1L]]))
  ))]
  assigned_column <- intersect(c("assigned", "in_topic"), names(assignment))
  if (length(assigned_column)) {
    assignment[, primary_assigned := .as_logical_flag(get(assigned_column[[1L]]))]
  } else {
    assignment[, primary_assigned := TRUE]
  }
  output <- unique(assignment[
    primary_assigned &
      term_group %in% c("GENE", "PEAK") &
      !is.na(term_id) & nzchar(term_id) &
      is.finite(primary_topic),
    .(
      term_id = as.character(term_id),
      term_group = as.character(term_group),
      primary_topic = as.integer(primary_topic)
    )
  ])
  conflicting <- output[, .(topics = data.table::uniqueN(primary_topic)),
                        by = .(term_group, term_id)][topics > 1L]
  if (nrow(conflicting)) {
    .log_abort(
      "{context} contains {nrow(conflicting)} terms with multiple primary topics."
    )
  }
  unique(output, by = c("term_group", "term_id"))
}

.m3_qc_compare_cross_k_topic_split <- function(reference_phi,
                                                candidate_phi,
                                                reference_assignment,
                                                candidate_assignment,
                                                supported_fraction = 0.20,
                                                dominant_fraction = 0.70,
                                                partial_fraction = 0.10) {
  reference_phi <- as.matrix(reference_phi)
  candidate_phi <- as.matrix(candidate_phi)
  if (!nrow(reference_phi) || !nrow(candidate_phi) ||
      !ncol(reference_phi) || !ncol(candidate_phi) ||
      is.null(colnames(reference_phi)) || is.null(colnames(candidate_phi)) ||
      !identical(colnames(reference_phi), colnames(candidate_phi))) {
    .log_abort("Cross-K phi matrices must contain the same ordered terms.")
  }
  cutoffs <- c(supported_fraction, dominant_fraction, partial_fraction)
  if (any(!is.finite(cutoffs)) ||
      supported_fraction <= 0 || supported_fraction >= dominant_fraction ||
      dominant_fraction >= 1 ||
      partial_fraction <= 0 || partial_fraction > supported_fraction) {
    .log_abort("Cross-K split fractions are invalid or inconsistent.")
  }
  reference_topics <- .m3_opt_topic_ids(reference_phi, "row")
  candidate_topics <- .m3_opt_topic_ids(candidate_phi, "row")
  reference <- .m3_qc_primary_term_assignments(
    reference_assignment,
    context = "reference topic assignment"
  )
  candidate <- .m3_qc_primary_term_assignments(
    candidate_assignment,
    context = "candidate topic assignment"
  )
  reference_counts <- reference[, .(terms = .N),
                                by = .(term_group, primary_topic)]
  reference_counts[, modality_total := sum(terms), by = term_group]
  reference_counts[, modality_fraction := terms / modality_total]
  dominant_score <- reference_counts[, .(
    dominance_score = sum(modality_fraction),
    assigned_terms = sum(terms)
  ), by = primary_topic]
  data.table::setorder(
    dominant_score,
    -dominance_score,
    -assigned_terms,
    primary_topic
  )
  reference_topic <- dominant_score$primary_topic[[1L]]

  affinity <- data.table::rbindlist(lapply(c("GENE", "PEAK"), function(group) {
    term_index <- startsWith(colnames(reference_phi), paste0(group, ":"))
    if (!any(term_index)) return(NULL)
    reference_probability <- reference_phi[, term_index, drop = FALSE]
    candidate_probability <- candidate_phi[, term_index, drop = FALSE]
    reference_probability[!is.finite(reference_probability) |
                            reference_probability < 0] <- 0
    candidate_probability[!is.finite(candidate_probability) |
                            candidate_probability < 0] <- 0
    reference_probability <- reference_probability / pmax(
      rowSums(reference_probability),
      .Machine$double.eps
    )
    candidate_probability <- candidate_probability / pmax(
      rowSums(candidate_probability),
      .Machine$double.eps
    )
    similarity <- sqrt(reference_probability) %*%
      t(sqrt(candidate_probability))
    output <- data.table::as.data.table(as.table(similarity))
    data.table::setnames(
      output,
      c("reference_position", "candidate_position", "hellinger_affinity")
    )
    output[, `:=`(
      term_group = group,
      reference_topic = reference_topics[match(
        as.character(reference_position),
        rownames(reference_phi)
      )],
      candidate_topic = candidate_topics[match(
        as.character(candidate_position),
        rownames(candidate_phi)
      )]
    )]
    output[, .(
      term_group,
      reference_topic,
      candidate_topic,
      hellinger_affinity = as.numeric(hellinger_affinity)
    )]
  }), use.names = TRUE)

  source_terms <- reference[primary_topic == reference_topic, .(
    term_group,
    term_id,
    reference_topic = primary_topic
  )]
  candidate_terms <- candidate[, .(
    term_group,
    term_id,
    candidate_topic = primary_topic
  )]
  mapped <- merge(
    source_terms,
    candidate_terms,
    by = c("term_group", "term_id"),
    all.x = TRUE,
    sort = FALSE
  )
  redistribution <- mapped[, .(terms = .N),
                           by = .(term_group, reference_topic, candidate_topic)]
  source_totals <- source_terms[, .(source_terms = .N), by = term_group]
  redistribution <- source_totals[redistribution, on = "term_group"]
  redistribution[, fraction := terms / source_terms]
  data.table::setorder(redistribution, term_group, -fraction, candidate_topic)

  children <- data.table::dcast(
    redistribution[!is.na(candidate_topic)],
    candidate_topic ~ term_group,
    value.var = "fraction",
    fill = 0
  )
  for (column in c("GENE", "PEAK")) {
    if (!column %in% names(children)) children[, (column) := 0]
  }
  data.table::setnames(
    children,
    c("GENE", "PEAK"),
    c("gene_fraction", "peak_fraction")
  )
  children[, joint_fraction := pmin(gene_fraction, peak_fraction)]
  data.table::setorder(
    children,
    -joint_fraction,
    -gene_fraction,
    -peak_fraction,
    candidate_topic
  )
  supported_children <- children[
    gene_fraction >= supported_fraction & peak_fraction >= supported_fraction,
    candidate_topic
  ]
  partial_children <- children[
    gene_fraction >= partial_fraction & peak_fraction >= partial_fraction,
    candidate_topic
  ]
  dominant_retained <- any(
    children$gene_fraction > dominant_fraction &
      children$peak_fraction > dominant_fraction
  )
  gene_top <- children[order(-gene_fraction, candidate_topic), candidate_topic][1L]
  peak_top <- children[order(-peak_fraction, candidate_topic), candidate_topic][1L]
  modality_disagreement <- !identical(gene_top, peak_top)
  classification <- if (
    length(supported_children) >= 2L && !dominant_retained
  ) {
    "supported_split"
  } else if (length(partial_children) >= 2L || modality_disagreement) {
    "partial_split"
  } else {
    "persistent_topic"
  }
  summary <- data.table::data.table(
    reference_topic = as.integer(reference_topic),
    classification = classification,
    dominant_retained = dominant_retained,
    modality_disagreement = modality_disagreement,
    supported_children = list(as.integer(supported_children)),
    partial_children = list(as.integer(partial_children))
  )
  list(
    summary = summary,
    children = children[],
    redistribution = redistribution[],
    affinity = affinity[]
  )
}

.m3_qc_gene_term_umap_plot <- function(optimization,
                                       topic_palette,
                                       gene_term_assignment = NULL,
                                       gene_umap = NULL,
                                       topic_ids = NULL,
                                       seed = 20260716L) {
  if (is.null(gene_umap)) {
    gene_umap <- .m3_qc_gene_term_umap_data(
      optimization,
      gene_term_assignment = gene_term_assignment,
      seed = seed
    )
  }
  background <- gene_umap[, .(UMAP1, UMAP2)]
  display_umap <- gene_umap[
    .as_logical_flag(assigned) & is.finite(topic_num)
  ]
  if (!is.null(topic_ids)) {
    topic_ids <- sort(unique(suppressWarnings(as.integer(topic_ids))))
    topic_ids <- topic_ids[is.finite(topic_ids)]
    topic_levels <- paste0("Topic ", topic_ids)
    display_umap <- gene_umap[topic_num %in% topic_ids]
    display_umap[, topic := factor(as.character(topic), levels = topic_levels)]
  }
  point_layer <- function(data, mapping, ...) {
    args <- c(list(data = data, mapping = mapping), list(...))
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      args$raster.dpi <- 320
      return(do.call(ggrastr::geom_point_rast, args))
    }
    do.call(ggplot2::geom_point, args)
  }
  colors <- topic_palette[levels(display_umap$topic)]
  x_range <- range(gene_umap$UMAP1, finite = TRUE)
  y_range <- range(gene_umap$UMAP2, finite = TRUE)
  shared_span <- max(diff(x_range), diff(y_range)) * 1.04
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  x_limits <- x_mid + c(-0.5, 0.5) * shared_span
  y_limits <- y_mid + c(-0.5, 0.5) * shared_span
  topic_levels <- levels(display_umap$topic)
  panel_labels <- stats::setNames(
    sub("^Topic ", "T", topic_levels),
    topic_levels
  )
  panel_columns <- max(1L, min(5L, length(topic_levels)))
  panel_rows <- max(1L, min(
    6L,
    ceiling(length(topic_levels) / panel_columns)
  ))
  ggplot2::ggplot() +
    point_layer(
      background,
      ggplot2::aes(UMAP1, UMAP2),
      color = "#CBD2D7",
      alpha = 0.22,
      size = 0.08
    ) +
    point_layer(
      display_umap,
      ggplot2::aes(UMAP1, UMAP2, color = topic),
      alpha = 0.72,
      size = 0.16
    ) +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::facet_wrap(
      ggplot2::vars(topic),
      ncol = panel_columns,
      nrow = panel_rows,
      drop = FALSE,
      labeller = ggplot2::as_labeller(panel_labels)
    ) +
    ggplot2::coord_fixed(
      xlim = x_limits,
      ylim = y_limits,
      ratio = 1,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "UMAP1",
      y = "UMAP2"
    ) +
    .m3_qc_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "#B9C2C9",
        fill = NA,
        linewidth = 0.45
      ),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      aspect.ratio = 1,
      panel.spacing = grid::unit(0.08, "in")
    )
}

.m3_qc_prepare_pathway_gene_sets <- function(pathway_gene_sets,
                                             expressed_genes = NULL,
                                             drop_shared = TRUE) {
  x <- data.table::as.data.table(data.table::copy(pathway_gene_sets))
  .assert_has_cols(
    x,
    c("pathway", "gene"),
    context = "pathway Gene-term UMAP"
  )
  x[, `:=`(
    pathway = trimws(as.character(pathway)),
    gene = toupper(trimws(as.character(gene)))
  )]
  x <- unique(x[
    !is.na(pathway) & nzchar(pathway) & !is.na(gene) & nzchar(gene),
    .(pathway, gene)
  ])
  if (!nrow(x)) {
    .log_abort("No valid pathway genes are available for the Gene-term UMAP.")
  }
  pathway_order <- unique(x$pathway)
  if (!is.null(expressed_genes)) {
    expressed_genes <- unique(toupper(trimws(as.character(expressed_genes))))
    expressed_genes <- expressed_genes[
      !is.na(expressed_genes) & nzchar(expressed_genes)
    ]
    x <- x[gene %in% expressed_genes]
  }
  if (isTRUE(drop_shared) && nrow(x)) {
    memberships <- x[, .(pathways = data.table::uniqueN(pathway)), by = gene]
    x <- x[memberships[pathways == 1L], on = "gene", nomatch = 0L]
  }
  if (!nrow(x)) {
    .log_abort("No expressed pathway-unique genes remain for the Gene-term UMAP.")
  }
  pathway_order <- pathway_order[pathway_order %in% x$pathway]
  data.table::setattr(x, "pathway_order", pathway_order)
  x[]
}

.m3_qc_pathway_gene_sets_from_enrichly <- function(
    pathways,
    expressed_genes,
    database = "MSigDB_Hallmark_2020",
    pathway_species = "human",
    cache_dir = NULL) {
  .assert_pkg("enrichly")
  pathways <- unique(trimws(as.character(pathways)))
  pathways <- pathways[!is.na(pathways) & nzchar(pathways)]
  if (!length(pathways)) {
    .log_abort("At least one pathway name is required.")
  }
  manifest <- enrichly::enrichly_download(
    databases = database,
    species = pathway_species,
    cache_dir = .enrichly_db_cache_dir(cache_dir),
    overwrite = FALSE,
    verbose = FALSE
  )
  pathway_db <- enrichly::enrichly_load(
    manifest$path,
    databases = database
  )$sets
  pathway_db <- data.table::as.data.table(pathway_db)
  missing <- setdiff(pathways, unique(as.character(pathway_db$term)))
  if (length(missing)) {
    .log_abort(
      "Pathway database {database} is missing: {paste(missing, collapse = ', ')}."
    )
  }
  rows <- lapply(pathways, function(pathway) {
    genes <- .canonicalize_pathway_genes(
      pathway_db[term == pathway, gene],
      pathway_species = pathway_species
    )
    data.table::data.table(pathway = pathway, gene = genes)
  })
  .m3_qc_prepare_pathway_gene_sets(
    data.table::rbindlist(rows, use.names = TRUE),
    expressed_genes = expressed_genes,
    drop_shared = TRUE
  )
}

.m3_qc_peak_term_umap_plot <- function(optimization,
                                       topic_palette,
                                       peak_term_assignment = NULL,
                                       peak_umap = NULL,
                                       topic_ids = NULL,
                                       top_n = 10000L,
                                       seed = 20260716L) {
  if (is.null(peak_umap)) {
    peak_umap <- .m3_qc_peak_term_umap_data(
      optimization,
      peak_term_assignment = peak_term_assignment,
      top_n = top_n,
      seed = seed
    )
  }
  display_umap <- data.table::copy(peak_umap)
  display_umap[, target_gene := peak_id]
  .m3_qc_gene_term_umap_plot(
    optimization,
    topic_palette = topic_palette,
    gene_umap = display_umap,
    topic_ids = topic_ids,
    seed = seed
  )
}

.m3_qc_pathway_gene_umap_plot <- function(gene_umap,
                                           pathway_gene_sets,
                                           pathway_colors = NULL) {
  gene_umap <- data.table::as.data.table(data.table::copy(gene_umap))
  .assert_has_cols(
    gene_umap,
    c("target_gene", "UMAP1", "UMAP2"),
    context = "pathway Gene-term UMAP coordinates"
  )
  pathway_order <- attr(pathway_gene_sets, "pathway_order")
  pathway_gene_sets <- .m3_qc_prepare_pathway_gene_sets(
    pathway_gene_sets,
    drop_shared = TRUE
  )
  if (is.null(pathway_order)) {
    pathway_order <- attr(pathway_gene_sets, "pathway_order")
  }
  gene_umap[, gene := toupper(trimws(as.character(target_gene)))]
  highlighted <- merge(
    pathway_gene_sets,
    gene_umap,
    by = "gene",
    all = FALSE,
    sort = FALSE
  )
  available_pathways <- unique(as.character(highlighted$pathway))
  omitted <- setdiff(pathway_order, available_pathways)
  if (length(omitted)) {
    .log_warn(
      "Omitting pathway UMAP panels without matching Gene terms: {paste(omitted, collapse = ', ')}."
    )
  }
  pathway_order <- pathway_order[pathway_order %in% available_pathways]
  if (!length(pathway_order)) {
    .log_abort("No pathway genes match the assigned Gene-term UMAP universe.")
  }
  highlighted[, pathway := factor(pathway, levels = pathway_order)]
  counts <- highlighted[, .(genes = data.table::uniqueN(gene)), by = pathway]
  compact_pathway <- c(
    `Apical Junction` = "Apical junction",
    `Cholesterol Homeostasis` = "Cholesterol",
    `Unfolded Protein Response` = "UPR",
    `Oxidative Phosphorylation` = "OXPHOS",
    `Epithelial Mesenchymal Transition` = "EMT"
  )
  display_pathway <- as.character(counts$pathway)
  matched_labels <- compact_pathway[display_pathway]
  display_pathway[!is.na(matched_labels)] <- matched_labels[!is.na(matched_labels)]
  panel_labels <- stats::setNames(
    sprintf("%s (n=%s)", display_pathway, scales::comma(counts$genes)),
    as.character(counts$pathway)
  )
  default_colors <- c(
    "#0072B2", "#E69F00", "#CC79A7", "#009E73", "#D55E00"
  )
  if (is.null(pathway_colors)) {
    pathway_colors <- stats::setNames(
      rep(default_colors, length.out = length(pathway_order)),
      pathway_order
    )
  } else {
    pathway_colors <- as.character(pathway_colors)
    if (is.null(names(pathway_colors))) {
      names(pathway_colors) <- pathway_order[seq_along(pathway_colors)]
    }
    missing_colors <- setdiff(pathway_order, names(pathway_colors))
    if (length(missing_colors)) {
      pathway_colors <- c(
        pathway_colors,
        stats::setNames(
          rep(default_colors, length.out = length(missing_colors)),
          missing_colors
        )
      )
    }
  }
  point_layer <- function(data, mapping, ...) {
    args <- c(list(data = data, mapping = mapping), list(...))
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      args$raster.dpi <- 320
      return(do.call(ggrastr::geom_point_rast, args))
    }
    do.call(ggplot2::geom_point, args)
  }
  x_range <- range(gene_umap$UMAP1, finite = TRUE)
  y_range <- range(gene_umap$UMAP2, finite = TRUE)
  shared_span <- max(diff(x_range), diff(y_range)) * 1.04
  if (!is.finite(shared_span) || shared_span <= 0) shared_span <- 1
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  ggplot2::ggplot() +
    point_layer(
      gene_umap[, .(UMAP1, UMAP2)],
      ggplot2::aes(UMAP1, UMAP2),
      color = "#CBD2D7",
      alpha = 0.24,
      size = 0.10
    ) +
    point_layer(
      highlighted,
      ggplot2::aes(UMAP1, UMAP2, color = pathway),
      alpha = 0.88,
      size = 0.32
    ) +
    ggplot2::scale_color_manual(
      values = pathway_colors[pathway_order],
      guide = "none"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(pathway),
      ncol = 5L,
      drop = FALSE,
      labeller = ggplot2::as_labeller(panel_labels)
    ) +
    ggplot2::coord_fixed(
      xlim = x_mid + c(-0.5, 0.5) * shared_span,
      ylim = y_mid + c(-0.5, 0.5) * shared_span,
      ratio = 1,
      expand = FALSE
    ) +
    ggplot2::labs(
      x = "UMAP1",
      y = "UMAP2"
    ) +
    .m3_qc_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "#B9C2C9",
        fill = NA,
        linewidth = 0.45
      ),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      aspect.ratio = 1,
      panel.spacing = grid::unit(0.10, "in")
    )
}

.m3_qc_prepare_tf_target_gene_sets <- function(tf_target_corr,
                                                module2_links,
                                                gene_universe,
                                                min_r = 0.5,
                                                top_n_tfs = 30L,
                                                top_n_targets = 250L) {
  tf_target_corr <- data.table::as.data.table(data.table::copy(tf_target_corr))
  module2_links <- data.table::as.data.table(data.table::copy(module2_links))
  .assert_has_cols(
    tf_target_corr,
    c("tf", "target_gene", "best_r", "pass"),
    context = "Module 2 TF-target correlations"
  )
  .assert_has_cols(
    module2_links,
    c("tf", "target_gene", "fp_id", "module2_link_pass"),
    context = "Module 2 TF-Peak-target links"
  )
  min_r <- as.numeric(min_r)[1L]
  top_n_tfs <- as.integer(top_n_tfs)[1L]
  top_n_targets <- as.integer(top_n_targets)[1L]
  if (!is.finite(min_r) || min_r < -1 || min_r > 1) {
    .log_abort("`min_r` must be between -1 and 1.")
  }
  if (!is.finite(top_n_tfs) || top_n_tfs < 1L ||
      !is.finite(top_n_targets) || top_n_targets < 1L) {
    .log_abort("`top_n_tfs` and `top_n_targets` must be positive integers.")
  }
  genes <- unique(toupper(trimws(as.character(gene_universe))))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  if (!length(genes)) {
    .log_abort("`gene_universe` must contain at least one Gene.")
  }

  tf_target_corr[, `:=`(
    tf = trimws(as.character(tf)),
    tf_key = toupper(trimws(as.character(tf))),
    target_gene = toupper(trimws(as.character(target_gene))),
    best_r = as.numeric(best_r),
    pass = as.logical(pass)
  )]
  tf_target_corr <- tf_target_corr[
    pass %in% TRUE & is.finite(best_r) & best_r >= min_r &
      nzchar(tf_key) & target_gene %in% genes
  ]
  if (!nrow(tf_target_corr)) {
    .log_abort("No Module 2 TF-target correlations pass the requested filters.")
  }
  data.table::setorderv(
    tf_target_corr,
    c("tf_key", "target_gene", "best_r", "tf"),
    c(1L, 1L, -1L, 1L),
    na.last = TRUE
  )
  tf_target_corr <- unique(
    tf_target_corr,
    by = c("tf_key", "target_gene")
  )

  module2_links[, `:=`(
    tf_key = toupper(trimws(as.character(tf))),
    target_gene = toupper(trimws(as.character(target_gene))),
    fp_id = trimws(as.character(fp_id)),
    module2_link_pass = as.logical(module2_link_pass)
  )]
  eligible_keys <- tf_target_corr[, .(tf_key, target_gene)]
  module2_links <- module2_links[
    module2_link_pass %in% TRUE & nzchar(fp_id) & target_gene %in% genes
  ][eligible_keys, on = c("tf_key", "target_gene"), nomatch = 0L]
  if (!nrow(module2_links)) {
    .log_abort("No passing Module 2 TF-Peak-target links support the TF-target correlations.")
  }
  support <- unique(module2_links[, .(tf_key, target_gene, fp_id)])[
    , .(supporting_peak_count = data.table::uniqueN(fp_id)),
    by = .(tf_key, target_gene)
  ]
  eligible <- merge(
    tf_target_corr,
    support,
    by = c("tf_key", "target_gene"),
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(eligible)) {
    .log_abort("No supported TF-target pairs remain after Module 2 link matching.")
  }

  tf_summary <- eligible[, .(
    eligible_target_count = data.table::uniqueN(target_gene),
    median_best_r = stats::median(best_r, na.rm = TRUE),
    tf = sort(unique(tf))[1L]
  ), by = tf_key]
  data.table::setorderv(
    tf_summary,
    c("eligible_target_count", "median_best_r", "tf"),
    c(-1L, -1L, 1L),
    na.last = TRUE
  )
  tf_summary <- tf_summary[seq_len(min(top_n_tfs, .N))]
  tf_summary[, tf_rank := seq_len(.N)]
  eligible <- merge(
    eligible,
    tf_summary[, .(
      tf_key,
      tf,
      tf_rank,
      eligible_target_count,
      median_best_r
    )],
    by = "tf_key",
    all = FALSE,
    sort = FALSE,
    suffixes = c("_source", "")
  )
  data.table::setorderv(
    eligible,
    c("tf_rank", "best_r", "supporting_peak_count", "target_gene"),
    c(1L, -1L, -1L, 1L),
    na.last = TRUE
  )
  eligible[, target_rank := seq_len(.N), by = tf_rank]
  eligible <- eligible[target_rank <= top_n_targets]
  optional <- intersect(c("best_method", "best_fdr"), names(eligible))
  result <- eligible[, c(
    list(
      tf_rank = as.integer(tf_rank),
      tf = as.character(tf),
      target_rank = as.integer(target_rank),
      target_gene = as.character(target_gene),
      best_r = as.numeric(best_r),
      supporting_peak_count = as.integer(supporting_peak_count),
      eligible_target_count = as.integer(eligible_target_count),
      median_best_r = as.numeric(median_best_r)
    ),
    lapply(.SD, identity)
  ), .SDcols = optional]
  data.table::setorderv(result, c("tf_rank", "target_rank"))
  attr(result, "min_r") <- min_r
  attr(result, "top_n_tfs") <- top_n_tfs
  attr(result, "top_n_targets") <- top_n_targets
  result[]
}

.m3_qc_select_topic_tf_target_gene_sets <- function(tf_target_gene_sets,
                                                     gene_term_assignment,
                                                     tf_panel = NULL,
                                                     top_n_tfs = 30L,
                                                     top_n_targets = 500L) {
  tf_target_gene_sets <- data.table::as.data.table(
    data.table::copy(tf_target_gene_sets)
  )
  assignment <- data.table::as.data.table(
    data.table::copy(gene_term_assignment)
  )
  .assert_has_cols(
    tf_target_gene_sets,
    c("tf", "target_gene", "best_r", "supporting_peak_count"),
    context = "eligible high-confidence TF-target Gene sets"
  )
  topic_column <- intersect(
    c("assigned_topic", "topic_num", "topic"),
    names(assignment)
  )[1L]
  if (is.na(topic_column)) {
    .log_abort("Gene-term assignment must contain an assigned topic column.")
  }
  if (!"target_gene" %in% names(assignment)) {
    if (!"term_id" %in% names(assignment)) {
      .log_abort("Gene-term assignment must contain `target_gene` or `term_id`.")
    }
    assignment[, target_gene := sub("^GENE:", "", as.character(term_id))]
  }
  if ("term_group" %in% names(assignment)) {
    assignment <- assignment[term_group == "GENE"]
  }
  if ("assigned" %in% names(assignment)) {
    assignment <- assignment[as.logical(assigned) %in% TRUE]
  }
  assignment[, `:=`(
    target_gene = toupper(trimws(as.character(target_gene))),
    selected_topic = suppressWarnings(as.integer(get(topic_column)))
  )]
  assignment <- unique(assignment[
    nzchar(target_gene) & is.finite(selected_topic),
    .(target_gene, selected_topic)
  ])
  duplicated_genes <- assignment[, .N, by = target_gene][N > 1L, target_gene]
  if (length(duplicated_genes)) {
    .log_abort("Gene-term assignment maps a target Gene to multiple topics.")
  }
  tf_target_gene_sets[, target_gene := toupper(trimws(as.character(target_gene)))]
  eligible <- merge(
    tf_target_gene_sets,
    assignment,
    by = "target_gene",
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(eligible)) {
    .log_abort("No eligible TF targets match assigned Gene terms.")
  }
  top_n_tfs <- as.integer(top_n_tfs)[1L]
  top_n_targets <- as.integer(top_n_targets)[1L]
  if (!is.finite(top_n_tfs) || top_n_tfs < 1L ||
      !is.finite(top_n_targets) || top_n_targets < 1L) {
    .log_abort("`top_n_tfs` and `top_n_targets` must be positive integers.")
  }
  tf_topic <- eligible[, .(
    topic_target_count = data.table::uniqueN(target_gene),
    topic_median_r = stats::median(best_r, na.rm = TRUE)
  ), by = .(tf, selected_topic)]
  topic_size <- assignment[, .(
    assigned_topic_genes = data.table::uniqueN(target_gene)
  ), by = selected_topic]
  tf_topic <- merge(
    tf_topic,
    topic_size,
    by = "selected_topic",
    all.x = TRUE,
    sort = FALSE
  )
  tf_topic[, tf_target_total := sum(topic_target_count), by = tf]
  total_assigned_genes <- sum(topic_size$assigned_topic_genes)
  tf_topic[, `:=`(
    topic_target_share = topic_target_count / pmax(tf_target_total, 1),
    topic_prevalence = assigned_topic_genes / pmax(total_assigned_genes, 1)
  )]
  tf_topic[, topic_selection_score :=
    log1p(topic_target_count) * topic_median_r *
      topic_target_share / pmax(topic_prevalence, 1e-8)]
  if (!is.null(tf_panel)) {
    tf_panel <- data.table::as.data.table(data.table::copy(tf_panel))
    .assert_has_cols(tf_panel, "tf", context = "TF-target QC panel")
    if (!"tf_function" %in% names(tf_panel)) {
      tf_panel[, tf_function := NA_character_]
    }
    tf_panel[, `:=`(
      tf = toupper(trimws(as.character(tf))),
      tf_function = as.character(tf_function),
      panel_rank = seq_len(.N)
    )]
    tf_panel <- unique(tf_panel[nzchar(tf)], by = "tf")
    tf_topic[, tf := toupper(trimws(as.character(tf)))]
    tf_topic <- tf_topic[
      tf_panel[, .(tf, tf_function, panel_rank)],
      on = "tf",
      nomatch = 0L
    ]
    if (!nrow(tf_topic)) {
      .log_abort("None of the requested TF panel has eligible assigned targets.")
    }
    candidate_tfs <- unique(
      tf_topic[order(panel_rank), .(tf, panel_rank, tf_function)]
    )[seq_len(min(top_n_tfs, .N))]
  } else {
    tf_topic[, `:=`(
      tf_function = NA_character_,
      panel_rank = NA_integer_
    )]
    candidate_tfs <- unique(tf_topic[, .(
      tf,
      tf_target_total,
      mean_r = stats::median(topic_median_r, na.rm = TRUE)
    )])
    data.table::setorderv(
      candidate_tfs,
      c("tf_target_total", "mean_r", "tf"),
      c(-1L, -1L, 1L),
      na.last = TRUE
    )
    candidate_tfs <- candidate_tfs[seq_len(min(top_n_tfs, .N))]
    candidate_tfs[, `:=`(
      panel_rank = seq_len(.N),
      tf_function = NA_character_
    )]
  }
  tf_topic <- tf_topic[candidate_tfs, on = "tf", nomatch = 0L]
  topic_priority <- unique(tf_topic[, .(
    selected_topic,
    assigned_topic_genes
  )])
  data.table::setorderv(
    topic_priority,
    c("assigned_topic_genes", "selected_topic"),
    c(-1L, 1L),
    na.last = TRUE
  )
  topic_priority[, topic_priority := seq_len(.N)]
  coverage_topics <- topic_priority[
    seq_len(min(top_n_tfs, .N)),
    selected_topic
  ]
  remaining_tfs <- as.character(candidate_tfs$tf)
  remaining_topics <- as.integer(coverage_topics)
  coverage_rows <- list()
  while (length(remaining_tfs) && length(remaining_topics)) {
    candidates <- tf_topic[
      tf %in% remaining_tfs & selected_topic %in% remaining_topics
    ]
    if (!nrow(candidates)) break
    data.table::setorderv(
      candidates,
      c(
        "topic_selection_score", "topic_target_count", "topic_median_r",
        "panel_rank", "tf", "selected_topic"
      ),
      c(-1L, -1L, -1L, 1L, 1L, 1L),
      na.last = TRUE
    )
    chosen <- candidates[1L]
    coverage_rows[[length(coverage_rows) + 1L]] <- chosen
    remaining_tfs <- setdiff(remaining_tfs, chosen$tf)
    remaining_topics <- setdiff(remaining_topics, chosen$selected_topic)
  }
  selected_tfs <- data.table::rbindlist(coverage_rows, fill = TRUE)
  remaining_tfs <- setdiff(candidate_tfs$tf, selected_tfs$tf)
  if (length(remaining_tfs)) {
    additional <- tf_topic[tf %in% remaining_tfs]
    data.table::setorderv(
      additional,
      c(
        "tf", "topic_selection_score", "topic_target_count",
        "topic_median_r", "selected_topic"
      ),
      c(1L, -1L, -1L, -1L, 1L),
      na.last = TRUE
    )
    additional <- unique(additional, by = "tf")
    selected_tfs <- data.table::rbindlist(
      list(selected_tfs, additional),
      fill = TRUE
    )
  }
  selected_tfs <- merge(
    selected_tfs,
    topic_priority[, .(selected_topic, topic_priority)],
    by = "selected_topic",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setorderv(
    selected_tfs,
    c(
      "topic_priority", "topic_selection_score", "panel_rank", "tf"
    ),
    c(1L, -1L, 1L, 1L),
    na.last = TRUE
  )
  selected_tfs[, topic_tf_rank := seq_len(.N), by = selected_topic]
  selected_tfs[, display_rank := seq_len(.N)]
  selected <- merge(
    eligible,
    selected_tfs[, .(
      tf,
      selected_topic,
      display_rank,
      topic_tf_rank,
      topic_target_count,
      assigned_topic_genes,
      tf_function,
      topic_selection_score
    )],
    by = c("tf", "selected_topic"),
    all = FALSE,
    sort = FALSE
  )
  data.table::setorderv(
    selected,
    c("display_rank", "best_r", "supporting_peak_count", "target_gene"),
    c(1L, -1L, -1L, 1L),
    na.last = TRUE
  )
  selected[, target_rank := seq_len(.N), by = display_rank]
  selected <- selected[target_rank <= top_n_targets]
  data.table::setorderv(selected, c("display_rank", "target_rank"))
  attr(selected, "min_r") <- attr(tf_target_gene_sets, "min_r") %||% 0.5
  attr(selected, "top_n_tfs") <- top_n_tfs
  attr(selected, "top_n_targets") <- top_n_targets
  selected[]
}

.m3_qc_tf_target_gene_umap_plot <- function(gene_umap,
                                             tf_target_gene_sets,
                                             min_r = NULL) {
  gene_umap <- data.table::as.data.table(data.table::copy(gene_umap))
  tf_target_gene_sets <- data.table::as.data.table(
    data.table::copy(tf_target_gene_sets)
  )
  .assert_has_cols(
    gene_umap,
    c("target_gene", "UMAP1", "UMAP2"),
    context = "TF-target Gene-term UMAP coordinates"
  )
  .assert_has_cols(
    tf_target_gene_sets,
    c("display_rank", "tf", "selected_topic", "target_gene", "best_r"),
    context = "high-confidence TF-target Gene sets"
  )
  gene_umap[, gene := toupper(trimws(as.character(target_gene)))]
  tf_target_gene_sets[, target_gene := toupper(trimws(as.character(target_gene)))]
  highlighted <- merge(
    tf_target_gene_sets,
    gene_umap,
    by.x = "target_gene",
    by.y = "gene",
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(highlighted)) {
    .log_abort("No high-confidence TF targets match the Gene-term UMAP universe.")
  }
  tf_order <- unique(
    tf_target_gene_sets[order(display_rank), .(display_rank, tf)]
  )$tf
  tf_order <- tf_order[tf_order %in% unique(highlighted$tf)]
  highlighted[, tf := factor(tf, levels = tf_order)]
  counts <- highlighted[, .(targets = data.table::uniqueN(target_gene)), by = tf]
  topic_lookup <- unique(tf_target_gene_sets[, .(tf, selected_topic)])
  counts <- merge(counts, topic_lookup, by = "tf", all.x = TRUE, sort = FALSE)
  panel_labels <- stats::setNames(
    sprintf(
      "%s | T%d (n=%s)",
      as.character(counts$tf),
      as.integer(counts$selected_topic),
      scales::comma(counts$targets)
    ),
    as.character(counts$tf)
  )
  if (is.null(min_r)) min_r <- attr(tf_target_gene_sets, "min_r") %||% 0.5
  min_r <- as.numeric(min_r)[1L]
  if (!is.finite(min_r)) min_r <- 0.5
  max_r <- max(1, max(highlighted$best_r, na.rm = TRUE))
  point_layer <- function(data, mapping, ...) {
    args <- c(list(data = data, mapping = mapping), list(...))
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      args$raster.dpi <- 320
      return(do.call(ggrastr::geom_point_rast, args))
    }
    do.call(ggplot2::geom_point, args)
  }
  x_range <- range(gene_umap$UMAP1, finite = TRUE)
  y_range <- range(gene_umap$UMAP2, finite = TRUE)
  shared_span <- max(diff(x_range), diff(y_range)) * 1.04
  if (!is.finite(shared_span) || shared_span <= 0) shared_span <- 1
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  ggplot2::ggplot() +
    point_layer(
      gene_umap[, .(UMAP1, UMAP2)],
      ggplot2::aes(UMAP1, UMAP2),
      color = "#CDD3D8",
      alpha = 0.20,
      size = 0.08
    ) +
    point_layer(
      highlighted,
      ggplot2::aes(UMAP1, UMAP2, color = best_r),
      alpha = 0.90,
      size = 0.24
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("#F4D35E", "#E66B2E", "#7A1F5C", "#183153"),
      limits = c(min_r, max_r),
      oob = scales::squish,
      name = "Module 2 R"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(tf),
      ncol = 5L,
      nrow = 6L,
      drop = FALSE,
      labeller = ggplot2::as_labeller(panel_labels)
    ) +
    ggplot2::coord_fixed(
      xlim = x_mid + c(-0.5, 0.5) * shared_span,
      ylim = y_mid + c(-0.5, 0.5) * shared_span,
      ratio = 1,
      expand = FALSE
    ) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    .m3_qc_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "#B9C2C9",
        fill = NA,
        linewidth = 0.40
      ),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      aspect.ratio = 1,
      panel.spacing = grid::unit(0.04, "in"),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.65, "in")
    )
}

.m3_qc_umap_plot <- function(data,
                             colour_column,
                             title,
                             subtitle = NULL,
                             colors = NULL,
                             background = NULL,
                             label_column = NULL,
                             label_style = c("box", "text"),
                             seed = 20260716L,
                             compact = FALSE) {
  label_style <- match.arg(label_style)
  foreground <- data[!is.na(get(colour_column))]
  p <- ggplot2::ggplot()
  point_layer <- function(data, mapping, colour = NULL, alpha, size) {
    args <- list(
      data = data,
      mapping = mapping,
      alpha = alpha,
      size = size
    )
    if (!is.null(colour)) args$colour <- colour
    if (requireNamespace("ggrastr", quietly = TRUE)) {
      args$raster.dpi <- 320
      return(do.call(ggrastr::geom_point_rast, args))
    }
    do.call(ggplot2::geom_point, args)
  }
  if (!is.null(background) && nrow(background)) {
    p <- p + point_layer(
      background,
      ggplot2::aes(UMAP1, UMAP2),
      colour = "#B9C1C7",
      alpha = 0.08,
      size = 0.10
    )
  }
  p <- p + point_layer(
    foreground,
    ggplot2::aes(UMAP1, UMAP2, colour = .data[[colour_column]]),
    alpha = 0.25,
    size = 0.14
  )
  if (is.null(colors)) {
    values <- sort(unique(as.character(foreground[[colour_column]])))
    colors <- stats::setNames(.module3_bright_topic_palette(values), values)
  }
  if (!is.null(label_column) && nrow(foreground)) {
    .assert_pkg("ggrepel")
    centroids <- foreground[, .(
      UMAP1 = stats::median(UMAP1),
      UMAP2 = stats::median(UMAP2),
      label = as.character(get(label_column)[[1L]])
    ), by = .(colour_value__ = get(colour_column))]
    if (identical(label_style, "text")) {
      p <- p + ggrepel::geom_text_repel(
        data = centroids,
        ggplot2::aes(
          UMAP1,
          UMAP2,
          label = label,
          colour = colour_value__
        ),
        family = "Helvetica",
        fontface = "bold",
        size = 3.2,
        box.padding = 0.25,
        point.padding = 0.14,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = as.integer(seed),
        segment.colour = NA,
        show.legend = FALSE
      )
    } else {
      p <- p + ggrepel::geom_label_repel(
        data = centroids,
        ggplot2::aes(
          UMAP1,
          UMAP2,
          label = label,
          fill = colour_value__
        ),
        family = "Helvetica",
        fontface = "bold",
        size = 3.2,
        colour = "white",
        box.padding = 0.3,
        point.padding = 0.18,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = as.integer(seed),
        label.size = 0.22,
        segment.colour = NA,
        show.legend = FALSE
      )
    }
  }
  p <- p +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = "UMAP 1", y = "UMAP 2") +
    .m3_qc_theme() +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "#B9C2C9",
        fill = NA,
        linewidth = 0.45
      )
    )
  if (isTRUE(compact)) {
    p <- p + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 3, 2, 3)
    )
  }
  p
}

.m3_qc_arrange <- function(...,
                           ncol = 2L,
                           title = NULL,
                           layout_matrix = NULL,
                           widths = NULL,
                           heights = NULL) {
  .assert_pkg("gridExtra")
  plots <- list(...)
  if (!is.null(title)) {
    title <- paste(
      base::strwrap(gsub("_", " ", as.character(title)), width = 92L),
      collapse = "\n"
    )
  }
  top <- if (is.null(title)) NULL else grid::textGrob(
    title,
    gp = grid::gpar(
      fontfamily = "Helvetica",
      fontface = "bold",
      fontsize = 12,
      col = "#20272E"
    )
  )
  if (is.null(layout_matrix)) {
    gridExtra::arrangeGrob(
      grobs = plots,
      ncol = ncol,
      top = top,
      widths = widths,
      heights = heights
    )
  } else {
    gridExtra::arrangeGrob(
      grobs = plots,
      layout_matrix = layout_matrix,
      top = top,
      widths = widths,
      heights = heights
    )
  }
}

.m3_qc_align_plot_dimensions <- function(plots) {
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  shared_widths <- Reduce(
    grid::unit.pmax,
    lapply(grobs, function(grob) grob$widths)
  )
  shared_heights <- Reduce(
    grid::unit.pmax,
    lapply(grobs, function(grob) grob$heights)
  )
  lapply(grobs, function(grob) {
    grob$widths <- shared_widths
    grob$heights <- shared_heights
    grob
  })
}

.m3_qc_pair_layout <- function(pair_count) {
  pair_count <- suppressWarnings(as.integer(pair_count)[[1L]])
  if (!is.finite(pair_count) || pair_count < 1L) return(NULL)
  pair_rows <- ceiling(pair_count / 2L)
  layout <- matrix(
    NA_integer_,
    nrow = pair_rows + 1L,
    ncol = 4L
  )
  layout[1L, ] <- 1L
  for (i in seq_len(pair_count)) {
    row <- 1L + ceiling(i / 2L)
    column <- if (i %% 2L == 1L) 1L else 3L
    grob <- 2L + (i - 1L) * 2L
    layout[row, column:(column + 1L)] <- c(grob, grob + 1L)
  }
  list(
    layout_matrix = layout,
    heights = c(1.35, rep(1, pair_rows))
  )
}

.m3_qc_condition_topic_matrix <- function(counts, value, conditions, topics) {
  out <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, paste0("Topic ", topics))
  )
  if (!nrow(counts)) return(out)
  rows <- match(counts$condition_id, conditions)
  columns <- match(counts$optimized_topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- counts[[value]][valid]
  out
}

.m3_qc_condition_gene_expression <- function(edges_docs, genes = NULL) {
  x <- data.table::as.data.table(edges_docs)
  gene_column <- if ("gene_key" %in% names(x)) {
    "gene_key"
  } else if ("target_gene" %in% names(x)) {
    "target_gene"
  } else {
    NA_character_
  }
  required <- c("condition_id", "gene_expr_condition")
  if (is.na(gene_column) || !all(required %in% names(x))) {
    return(data.table::data.table(
      condition_id = character(),
      target_gene = character(),
      expression = numeric()
    ))
  }
  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    x <- x[get(gene_column) %in% genes]
  }
  x <- x[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(get(gene_column)) & nzchar(get(gene_column)) &
      is.finite(gene_expr_condition) & gene_expr_condition >= 0,
    .(
      expression = max(as.numeric(gene_expr_condition), na.rm = TRUE)
    ),
    by = .(
      condition_id = as.character(condition_id),
      target_gene = as.character(get(gene_column))
    )
  ]
  data.table::setorder(x, condition_id, target_gene)
  x[]
}

.m3_qc_condition_topic_expression_matrix <- function(condition_gene_expression,
                                                       pair_assignment,
                                                       conditions,
                                                       topics) {
  out <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, paste0("Topic ", topics))
  )
  expression <- data.table::as.data.table(condition_gene_expression)
  pairs <- data.table::as.data.table(pair_assignment)
  required_expression <- c("condition_id", "target_gene", "expression")
  required_pairs <- c("target_gene", "optimized_assigned_topic")
  if (!all(required_expression %in% names(expression)) ||
      !all(required_pairs %in% names(pairs))) {
    return(out)
  }
  pools <- unique(pairs[
    is.finite(optimized_assigned_topic),
    .(
      target_gene = as.character(target_gene),
      topic = as.integer(optimized_assigned_topic)
    )
  ])
  expressed <- unique(expression[
    is.finite(expression) & expression > 0,
    as.character(target_gene)
  ])
  pools <- pools[topic %in% topics & target_gene %in% expressed]
  if (!nrow(pools)) return(out)
  pool_sizes <- pools[, .(pool_size = data.table::uniqueN(target_gene)), by = topic]
  values <- merge(
    expression[condition_id %in% conditions],
    pools,
    by = "target_gene",
    all = FALSE,
    sort = FALSE
  )
  values[, expression_log2 := log2(pmax(as.numeric(expression), 0) + 1)]
  means <- values[, .(
    expression_sum = sum(expression_log2, na.rm = TRUE)
  ), by = .(condition_id, topic)]
  means <- merge(means, pool_sizes, by = "topic", all.x = TRUE, sort = FALSE)
  means[, mean_expression := expression_sum / pool_size]
  rows <- match(means$condition_id, conditions)
  columns <- match(means$topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- means$mean_expression[valid]
  out
}

.m3_qc_document_supported_topic_expression_matrix <- function(optimization,
                                                               conditions,
                                                               topics) {
  out <- matrix(
    NA_real_,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, paste0("Topic ", topics))
  )
  assignments <- data.table::as.data.table(optimization$qc$assignments)
  expression <- data.table::as.data.table(
    optimization$condition_gene_expression
  )
  target_levels <- as.character(optimization$qc$target_levels)
  required_assignments <- c(
    "doc_index", "target_index", "condition_id", "optimized_topic"
  )
  required_expression <- c("condition_id", "target_gene", "expression")
  if (!all(required_assignments %in% names(assignments)) ||
      !all(required_expression %in% names(expression)) ||
      !length(target_levels)) {
    return(out)
  }
  values <- assignments[
    condition_id %in% conditions &
      is.finite(target_index) & target_index >= 1L &
      target_index <= length(target_levels) &
      is.finite(optimized_topic) & optimized_topic %in% topics,
    .(
      doc_index = as.integer(doc_index),
      condition_id = as.character(condition_id),
      target_gene = target_levels[as.integer(target_index)],
      topic = as.integer(optimized_topic)
    )
  ]
  values <- unique(values)
  if (!nrow(values)) return(out)
  values <- merge(
    values,
    expression[condition_id %in% conditions],
    by = c("condition_id", "target_gene"),
    all = FALSE,
    sort = FALSE
  )
  values <- values[is.finite(expression) & expression >= 0]
  if (!nrow(values)) return(out)
  values[, expression_log2 := log2(as.numeric(expression) + 1)]
  means <- values[, .(
    mean_expression = mean(expression_log2),
    supporting_documents = data.table::uniqueN(doc_index),
    supporting_targets = data.table::uniqueN(target_gene)
  ), by = .(condition_id, topic)]
  rows <- match(means$condition_id, conditions)
  columns <- match(means$topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- means$mean_expression[valid]
  attr(out, "support") <- means
  out
}

.m3_qc_document_expression_topic_share_matrix <- function(optimization,
                                                           conditions,
                                                           topics) {
  out <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, paste0("Topic ", topics))
  )
  assignments <- data.table::as.data.table(optimization$qc$assignments)
  expression <- data.table::as.data.table(
    optimization$condition_gene_expression
  )
  target_levels <- as.character(optimization$qc$target_levels)
  required_assignments <- c(
    "doc_index", "target_index", "condition_id", "optimized_topic"
  )
  required_expression <- c("condition_id", "target_gene", "expression")
  if (!all(required_assignments %in% names(assignments)) ||
      !all(required_expression %in% names(expression)) ||
      !length(target_levels)) {
    return(out)
  }
  documents <- unique(assignments[
    condition_id %in% conditions,
    .(doc_index = as.integer(doc_index), condition_id = as.character(condition_id))
  ])
  document_counts <- documents[, .(document_count = .N), by = condition_id]
  values <- assignments[
    condition_id %in% conditions &
      is.finite(target_index) & target_index >= 1L &
      target_index <= length(target_levels) &
      is.finite(optimized_topic) & optimized_topic %in% topics,
    .(
      doc_index = as.integer(doc_index),
      condition_id = as.character(condition_id),
      target_gene = target_levels[as.integer(target_index)],
      topic = as.integer(optimized_topic)
    )
  ]
  values <- unique(values)
  if (!nrow(values)) return(out)
  values <- merge(
    values,
    expression[condition_id %in% conditions],
    by = c("condition_id", "target_gene"),
    all = FALSE,
    sort = FALSE
  )
  values <- values[is.finite(expression) & expression >= 0]
  if (!nrow(values)) return(out)
  values[, expression_mass := log2(as.numeric(expression) + 1)]
  document_topic <- values[, .(
    expression_mass = sum(expression_mass)
  ), by = .(doc_index, condition_id, topic)]
  document_topic[, topic_share := expression_mass / sum(expression_mass), by = doc_index]
  means <- document_topic[, .(
    topic_share_sum = sum(topic_share)
  ), by = .(condition_id, topic)]
  means <- merge(means, document_counts, by = "condition_id", all.x = TRUE)
  means[, mean_topic_share := topic_share_sum / document_count]
  rows <- match(means$condition_id, conditions)
  columns <- match(means$topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- means$mean_topic_share[valid]
  out
}

.m3_qc_upregulated_genes_from_comparisons <- function(comparison_genes,
                                                       genes = NULL,
                                                       reference_condition,
                                                       log2fc_min = 1,
                                                       pseudocount = 1) {
  expression <- .m3_qc_condition_expression_from_comparisons(
    comparison_genes,
    genes = genes
  )
  .m3_qc_upregulated_genes_vs_reference(
    expression,
    reference_condition = reference_condition,
    log2fc_min = log2fc_min,
    pseudocount = pseudocount
  )
}

.m3_qc_condition_expression_from_comparisons <- function(comparison_genes,
                                                           genes = NULL) {
  x <- data.table::as.data.table(comparison_genes)
  long_required <- c("condition_id", "gene_key", "expression")
  if (all(long_required %in% names(x))) {
    if (!is.null(genes)) {
      genes <- unique(as.character(genes))
      genes <- genes[!is.na(genes) & nzchar(genes)]
      x <- x[gene_key %in% genes]
    }
    out <- x[
      !is.na(condition_id) & nzchar(condition_id) &
        !is.na(gene_key) & nzchar(gene_key) &
        is.finite(expression) & expression >= 0,
      .(expression = max(as.numeric(expression))),
      by = .(
        condition_id = as.character(condition_id),
        target_gene = as.character(gene_key)
      )
    ]
    data.table::setorder(out, condition_id, target_gene)
    return(out[])
  }
  required <- c(
    "condition_1", "condition_2", "gene_key", "expression_1", "expression_2"
  )
  if (!all(required %in% names(x))) {
    return(data.table::data.table(
      condition_id = character(),
      target_gene = character(),
      expression = numeric()
    ))
  }
  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    x <- x[gene_key %in% genes]
  }
  out <- data.table::rbindlist(list(
    x[, .(
      condition_id = as.character(condition_1),
      target_gene = as.character(gene_key),
      expression = as.numeric(expression_1)
    )],
    x[, .(
      condition_id = as.character(condition_2),
      target_gene = as.character(gene_key),
      expression = as.numeric(expression_2)
    )]
  ))
  out <- out[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(target_gene) & nzchar(target_gene) &
      is.finite(expression) & expression >= 0,
    .(expression = max(expression)),
    by = .(condition_id, target_gene)
  ]
  data.table::setorder(out, condition_id, target_gene)
  out[]
}

.m3_qc_upregulated_genes_vs_reference <- function(condition_gene_expression,
                                                   reference_condition,
                                                   log2fc_min = 1,
                                                   pseudocount = 1) {
  expression <- data.table::copy(
    data.table::as.data.table(condition_gene_expression)
  )
  reference_condition <- as.character(reference_condition %||% "")[[1L]]
  log2fc_min <- suppressWarnings(as.numeric(log2fc_min)[[1L]])
  pseudocount <- suppressWarnings(as.numeric(pseudocount)[[1L]])
  empty <- data.table::data.table(
    condition_id = character(),
    target_gene = character()
  )
  required <- c("condition_id", "target_gene", "expression")
  if (!nzchar(reference_condition) ||
      !all(required %in% names(expression)) ||
      !is.finite(log2fc_min) || log2fc_min < 0 ||
      !is.finite(pseudocount) || pseudocount <= 0) {
    return(empty)
  }
  expression <- expression[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(target_gene) & nzchar(target_gene) &
      is.finite(expression) & expression >= 0,
    .(expression = max(as.numeric(expression))),
    by = .(
      condition_id = as.character(condition_id),
      target_gene = as.character(target_gene)
    )
  ]
  reference <- expression[
    condition_id == reference_condition,
    .(target_gene, reference_expression = expression)
  ]
  if (!nrow(reference)) return(empty)
  values <- merge(
    expression[condition_id != reference_condition],
    reference,
    by = "target_gene",
    all = FALSE,
    sort = FALSE
  )
  values[, log2fc_vs_reference := log2(
    (expression + pseudocount) / (reference_expression + pseudocount)
  )]
  out <- unique(values[
    is.finite(log2fc_vs_reference) & log2fc_vs_reference >= log2fc_min,
    .(condition_id, target_gene)
  ])
  data.table::setorder(out, condition_id, target_gene)
  out[]
}

.m3_qc_condition_topic_upregulated_matrix <- function(condition_upregulated_genes,
                                                        pair_assignment,
                                                        conditions,
                                                        topics) {
  out <- matrix(
    0,
    nrow = length(conditions),
    ncol = length(topics),
    dimnames = list(conditions, paste0("Topic ", topics))
  )
  upregulated <- data.table::as.data.table(condition_upregulated_genes)
  pairs <- data.table::as.data.table(pair_assignment)
  if (!all(c("condition_id", "target_gene") %in% names(upregulated)) ||
      !all(c("target_gene", "optimized_assigned_topic") %in% names(pairs))) {
    return(out)
  }
  assigned <- unique(pairs[
    is.finite(optimized_assigned_topic),
    .(
      target_gene = as.character(target_gene),
      topic = as.integer(optimized_assigned_topic)
    )
  ])
  values <- merge(
    unique(upregulated[condition_id %in% conditions]),
    assigned[topic %in% topics],
    by = "target_gene",
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(values)) return(out)
  counts <- values[, .(
    upregulated_genes = data.table::uniqueN(target_gene)
  ), by = .(condition_id, topic)]
  rows <- match(counts$condition_id, conditions)
  columns <- match(counts$topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- counts$upregulated_genes[valid]
  out
}

.m3_qc_pooled_tf_topic_matrix <- function(optimization,
                                           top_n_tfs = 150L) {
  top_n_tfs <- suppressWarnings(as.integer(top_n_tfs)[[1L]])
  if (!is.finite(top_n_tfs) || top_n_tfs < 1L) {
    .log_abort("top_n_tfs must be a positive integer.")
  }
  qc <- optimization$qc
  x <- qc$assignments[optimized_aligned == TRUE]
  topics <- as.integer(qc$optimized_topic_ids)
  if (!nrow(x) || !length(topics)) return(matrix(numeric(), 0L, 0L))
  per_topic_n <- max(1L, floor(top_n_tfs / length(topics)))
  topic_counts <- x[, .(
    targets = data.table::uniqueN(target_index)
  ), by = .(tf_index, optimized_topic)]
  topic_counts[, tf := as.character(qc$tf_levels[tf_index])]
  data.table::setorder(topic_counts, optimized_topic, -targets, tf)
  selected_ids <- unique(topic_counts[
    optimized_topic %in% topics,
    head(.SD, per_topic_n),
    by = optimized_topic
  ]$tf_index)
  totals <- x[, .(
    targets = data.table::uniqueN(target_index)
  ), by = tf_index]
  totals[, tf := as.character(qc$tf_levels[tf_index])]
  data.table::setorder(totals, -targets, tf)
  selected <- totals[tf_index %in% selected_ids]
  counts <- x[tf_index %in% selected$tf_index, .(
    targets = data.table::uniqueN(target_index)
  ), by = .(tf_index, optimized_topic)]
  out <- matrix(
    0,
    nrow = nrow(selected),
    ncol = length(topics),
    dimnames = list(selected$tf, paste0("Topic ", topics))
  )
  rows <- match(counts$tf_index, selected$tf_index)
  columns <- match(counts$optimized_topic, topics)
  valid <- is.finite(rows) & is.finite(columns)
  out[cbind(rows[valid], columns[valid])] <- counts$targets[valid]
  attr(out, "selection_per_topic") <- per_topic_n
  out
}

.m3_qc_condition_topic_theta <- function(optimization) {
  theta <- as.matrix(optimization$theta)
  if (identical(optimization$document_design %||% "condition_tf", "condition")) {
    if (!nrow(theta) || !ncol(theta)) {
      .log_abort("Condition-document theta must contain conditions and topics.")
    }
    if (is.null(rownames(theta)) || anyNA(rownames(theta)) ||
        any(!nzchar(rownames(theta))) || anyDuplicated(rownames(theta))) {
      .log_abort(
        "Condition-document theta must have one uniquely named row per condition."
      )
    }
    topic_ids <- .m3_opt_topic_ids(theta, "column")
    theta <- .m3_opt_row_normalize(theta)
    theta <- theta[order(rownames(theta)), , drop = FALSE]
    colnames(theta) <- paste0("Topic ", topic_ids)
    return(theta)
  }
  assignments <- optimization$qc$assignments
  if (!nrow(theta) || !nrow(assignments)) {
    return(matrix(numeric(), 0L, 0L))
  }
  doc_conditions <- unique(assignments[, .(doc_index, condition_id)])
  conflicts <- doc_conditions[, data.table::uniqueN(condition_id), by = doc_index][
    V1 > 1L
  ]
  if (nrow(conflicts)) {
    .log_abort("Each topic document must map to exactly one condition.")
  }
  doc_conditions <- doc_conditions[
    is.finite(doc_index) & doc_index >= 1L & doc_index <= nrow(theta)
  ]
  conditions <- sort(unique(as.character(doc_conditions$condition_id)))
  topic_ids <- .m3_opt_topic_ids(theta, "column")
  out <- vapply(conditions, function(condition) {
    rows <- doc_conditions[condition_id == condition, doc_index]
    colMeans(theta[rows, , drop = FALSE])
  }, numeric(ncol(theta)))
  out <- t(out)
  rownames(out) <- conditions
  colnames(out) <- paste0("Topic ", topic_ids)
  out
}

.m3_qc_condition_tf_topic_theta <- function(theta, matched_tfs = TRUE) {
  theta <- as.matrix(theta)
  if (!nrow(theta) || !ncol(theta)) {
    .log_abort("Condition::TF theta must contain documents and topics.")
  }
  document_ids <- rownames(theta)
  if (is.null(document_ids) || anyNA(document_ids) ||
      any(!nzchar(document_ids)) || anyDuplicated(document_ids) ||
      any(!grepl("::[^:]+$", document_ids))) {
    .log_abort(
      "Condition::TF theta must have unique '<condition>::<TF>' row names."
    )
  }
  document_index <- data.table::data.table(
    row_id = seq_len(nrow(theta)),
    condition_id = sub("::[^:]+$", "", document_ids),
    tf = sub("^.*::", "", document_ids)
  )
  if (document_index[, anyDuplicated(paste(condition_id, tf))] > 0L) {
    .log_abort("Condition::TF theta contains duplicated condition/TF documents.")
  }
  conditions <- sort(unique(document_index$condition_id))
  shared_tfs <- if (isTRUE(matched_tfs)) {
    Reduce(intersect, split(document_index$tf, document_index$condition_id))
  } else {
    character()
  }
  use_matched <- isTRUE(matched_tfs) && length(shared_tfs) > 0L
  if (use_matched) {
    document_index <- document_index[tf %in% shared_tfs]
  }
  theta <- .m3_opt_row_normalize(theta)
  out <- vapply(conditions, function(condition) {
    rows <- document_index[condition_id == condition, row_id]
    colMeans(theta[rows, , drop = FALSE])
  }, numeric(ncol(theta)))
  out <- t(out)
  rownames(out) <- conditions
  colnames(out) <- paste0("Topic ", .m3_opt_topic_ids(theta, "column"))
  attr(out, "matched_tfs") <- use_matched
  attr(out, "tf_count") <- if (use_matched) {
    length(shared_tfs)
  } else {
    length(unique(document_index$tf))
  }
  out
}

.m3_qc_top_condition_topics <- function(theta,
                                         jaccard,
                                         top_n = 3L) {
  top_n <- suppressWarnings(as.integer(top_n)[[1L]])
  if (!is.finite(top_n) || top_n < 1L) {
    .log_abort("top_n must be a positive integer.")
  }
  theta <- as.matrix(theta)
  jaccard <- as.matrix(jaccard)
  if (!nrow(theta) || !ncol(theta)) return(data.table::data.table())
  theta_long <- data.table::as.data.table(as.table(theta))
  data.table::setnames(
    theta_long,
    c("condition_id", "topic_label", "mean_theta")
  )
  theta_long[, `:=`(
    condition_id = as.character(condition_id),
    topic_label = as.character(topic_label),
    topic_num = suppressWarnings(as.integer(sub("^Topic ", "", topic_label))),
    mean_theta = as.numeric(get("mean_theta"))
  )]
  jaccard_long <- data.table::as.data.table(as.table(jaccard))
  data.table::setnames(
    jaccard_long,
    c("condition_id", "topic_num", "mean_jaccard")
  )
  jaccard_long[, `:=`(
    condition_id = as.character(condition_id),
    topic_num = suppressWarnings(as.integer(as.character(topic_num))),
    mean_jaccard = as.numeric(get("mean_jaccard"))
  )]
  out <- merge(
    theta_long,
    jaccard_long,
    by = c("condition_id", "topic_num"),
    all = FALSE,
    sort = FALSE
  )
  data.table::setorderv(
    out,
    c("condition_id", "mean_theta", "mean_jaccard", "topic_num"),
    c(1L, -1L, -1L, 1L)
  )
  out[, rank := seq_len(.N), by = condition_id]
  out[rank <= top_n]
}

.m3_qc_retention_plot <- function(labels,
                                  counts,
                                  title,
                                  fill = "#007C78") {
  x <- data.table::data.table(
    stage = factor(labels, levels = rev(labels)),
    count = as.numeric(counts)
  )
  count_max <- max(x$count, na.rm = TRUE)
  if (!is.finite(count_max) || count_max <= 0) count_max <- 1
  x[, label_inside := count >= 0.16 * count_max]
  ggplot2::ggplot(x, ggplot2::aes(count, stage)) +
    ggplot2::geom_col(
      width = 0.72,
      fill = fill,
      colour = "#30383F",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = scales::comma(count),
        hjust = data.table::fifelse(label_inside, 1.06, -0.08),
        colour = label_inside
      ),
      family = "Helvetica",
      fontface = "bold",
      size = 3.2
    ) +
    ggplot2::scale_colour_manual(
      values = c(`TRUE` = "white", `FALSE` = "#171717"),
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(title = title, x = "Full-universe count", y = NULL) +
    .m3_qc_theme()
}

.m3_qc_retention_labels <- function() {
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
}

.m3_qc_unique_link_count <- function(rows) {
  rows <- data.table::as.data.table(rows)
  if (!nrow(rows)) return(0L)
  data.table::uniqueN(rows, by = c("tf_index", "target_index"))
}

.m3_qc_unique_target_count <- function(values) {
  values <- values[!is.na(values)]
  if (is.character(values)) values <- values[nzchar(values)]
  data.table::uniqueN(values)
}

.m3_qc_retention_page <- function(optimization,
                                  title_prefix = NULL,
                                  standardized = FALSE) {
  qc <- optimization$qc
  assignments <- data.table::as.data.table(qc$assignments)
  pairs <- data.table::as.data.table(optimization$raw_pair_assignment)
  gene_only_assignment <- (
    optimization$assignment_mode %||% "gene_peak"
  ) %in% c("gene", "gene_only")
  .assert_has_cols(
    assignments,
    c("tf_index", "target_index", "raw_aligned"),
    context = "topic-assignment QC retention links"
  )
  .assert_has_cols(
    pairs,
    c(
      "target_gene", "gene_gammafit_topics", "assigned_topic"
    ),
    context = "topic-assignment QC retention targets"
  )
  if (!gene_only_assignment) {
    .assert_has_cols(
      pairs,
      "peak_gammafit_topics",
      context = "topic-assignment QC retention targets"
    )
  }
  pair_candidates <- .m3_opt_parse_topics(
    pairs$gene_gammafit_topics,
    qc$raw_topic_ids
  )
  peak_candidates <- if (gene_only_assignment) {
    pair_candidates
  } else {
    .m3_opt_parse_topics(
      pairs$peak_gammafit_topics,
      qc$raw_topic_ids
    )
  }
  target_gamma <- lengths(pair_candidates) > 0L &
    lengths(peak_candidates) > 0L
  target_max <- is.finite(suppressWarnings(as.integer(pairs$assigned_topic)))
  assignment_pair_index <- if ("pair_index" %in% names(assignments)) {
    assignments$pair_index
  } else {
    assignments$target_index
  }
  link_max <- target_max[assignment_pair_index]
  retention_labels <- .m3_qc_retention_labels()
  retention_links <- c(
    .m3_qc_unique_link_count(assignments),
    .m3_qc_unique_link_count(
      assignments[target_gamma[assignment_pair_index] %in% TRUE]
    ),
    .m3_qc_unique_link_count(assignments[link_max %in% TRUE]),
    .m3_qc_unique_link_count(assignments[raw_aligned == TRUE])
  )
  retention_genes <- if (
    identical(optimization$assignment_mode, "tf_target") &&
      !is.null(optimization$raw_correspondence_assignment)
  ) {
    gene_rows <- data.table::as.data.table(
      optimization$raw_correspondence_assignment
    )
    gene_gamma <- lengths(.m3_opt_parse_topics(
      gene_rows$gene_gammafit_topics,
      qc$raw_topic_ids
    )) > 0L
    gene_max <- .as_logical_flag(gene_rows$assigned)
    c(
      .m3_qc_unique_target_count(gene_rows$target_gene),
      .m3_qc_unique_target_count(gene_rows[gene_gamma, target_gene]),
      .m3_qc_unique_target_count(gene_rows[gene_max, target_gene]),
      .m3_qc_unique_target_count(
        assignments[raw_aligned == TRUE, target_index]
      )
    )
  } else {
    c(
      .m3_qc_unique_target_count(pairs$target_gene),
      .m3_qc_unique_target_count(pairs[target_gamma %in% TRUE, target_gene]),
      .m3_qc_unique_target_count(pairs[target_max %in% TRUE, target_gene]),
      .m3_qc_unique_target_count(
        assignments[raw_aligned == TRUE, target_index]
      )
    )
  }
  link_title <- if (isTRUE(standardized)) {
    "Unique TF-target pairs"
  } else {
    "TF-target links"
  }
  gene_title <- if (isTRUE(standardized)) {
    "Unique target genes"
  } else {
    "Target genes"
  }
  page_title <- if (isTRUE(standardized) && !is.null(title_prefix)) {
    paste0(title_prefix, ": Topic assignment retention")
  } else {
    "Topic assignment retention"
  }
  .m3_qc_arrange(
    .m3_qc_retention_plot(
      retention_labels$links,
      retention_links,
      link_title,
      fill = "#007C78"
    ),
    .m3_qc_retention_plot(
      retention_labels$genes,
      retention_genes,
      gene_title,
      fill = "#D97824"
    ),
    ncol = 1L,
    title = page_title
  )
}

.m3_qc_style_heatmap_grob <- function(grob) {
  if (!is.null(grob$gp)) {
    grob$gp$font <- NULL
    grob$gp$fontfamily <- "Helvetica"
    grob$gp$fontface <- "bold"
  }
  if (!is.null(grob$children) && length(grob$children)) {
    grob$children <- lapply(grob$children, .m3_qc_style_heatmap_grob)
  }
  if (!is.null(grob$grobs) && length(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, .m3_qc_style_heatmap_grob)
  }
  grob
}

.m3_qc_tf_topic_pages <- function(optimization,
                                  top_n_tfs = 150L) {
  mat <- .m3_qc_pooled_tf_topic_matrix(
    optimization,
    top_n_tfs = top_n_tfs
  )
  if (!nrow(mat) || !ncol(mat)) return(list())
  column_order <- .m3_qc_cluster_order(mat, "column")
  work <- sqrt(pmax(mat, 0))
  means <- rowMeans(work)
  sds <- apply(work, 1L, stats::sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  cluster_work <- sweep(sweep(work, 1L, means, "-"), 1L, sds, "/")
  cluster_work[!is.finite(cluster_work)] <- 0
  row_cluster <- if (nrow(cluster_work) >= 2L) {
    stats::hclust(stats::dist(cluster_work), method = "ward.D2")
  } else {
    FALSE
  }
  colors <- grDevices::colorRampPalette(
    c("#F5F7F7", "#2A788E", "#FDE725")
  )(100L)
  raw_max <- max(mat, na.rm = TRUE)
  legend_labels <- unique(c(
    0,
    pretty(c(0, raw_max), n = 4L),
    raw_max
  ))
  legend_labels <- legend_labels[
    is.finite(legend_labels) & legend_labels >= 0 & legend_labels <= raw_max
  ]
  legend_labels <- sort(unique(round(legend_labels)))
  legend_breaks <- sqrt(legend_labels)
  color_limits <- range(work, finite = TRUE)
  if (length(color_limits) != 2L || !all(is.finite(color_limits)) ||
      diff(color_limits) <= 0) {
    color_limits <- c(0, max(c(work, 1), na.rm = TRUE))
  }
  ph <- pheatmap::pheatmap(
    work[, column_order, drop = FALSE],
    cluster_rows = row_cluster,
    cluster_cols = FALSE,
    color = colors,
    breaks = seq(color_limits[[1L]], color_limits[[2L]], length.out = 101L),
    border_color = "#D6DCE1",
    treeheight_row = if (nrow(work) >= 2L) 42 else 0,
    treeheight_col = 0,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 9,
    fontsize_row = 5,
    fontsize_col = 7,
    angle_col = 90,
    main = "Unique assigned target genes",
    legend_breaks = legend_breaks,
    legend_labels = scales::comma(legend_labels),
    legend = TRUE,
    silent = TRUE
  )
  plot <- .m3_qc_style_heatmap_grob(ph$gtable)
  list(.m3_qc_arrange(
    plot,
    ncol = 1L,
    title = sprintf(
      "Top TFs balanced across topics (%d candidates/topic; %d unique TFs)",
      as.integer(attr(mat, "selection_per_topic") %||% 1L),
      nrow(mat)
    )
  ))
}

.m3_qc_apply_primary_tf_counts <- function(optimization,
                                           tf_topic_evidence) {
  evidence <- data.table::as.data.table(tf_topic_evidence)
  .assert_has_cols(
    evidence,
    c(
      "tf", "topic_num", "condition_id",
      "global_primary_topic", "condition_primary_topic",
      "global_primary_confident", "condition_primary_confident"
    ),
    context = "TF-topic evidence for assignment QC"
  )
  qc <- optimization$qc
  global_primary <- unique(evidence[global_primary_confident == TRUE, .(
    tf = as.character(tf),
    topic_num = as.integer(global_primary_topic)
  )])
  global_counts <- global_primary[, .(
    primary_tfs = data.table::uniqueN(tf)
  ), by = topic_num]
  update_primary_counts <- function(counts, topic_column) {
    x <- data.table::copy(counts)
    x <- global_counts[
      x,
      on = c("topic_num" = topic_column)
    ]
    x[, tfs := data.table::fcoalesce(as.integer(primary_tfs), 0L)]
    x[, primary_tfs := NULL]
    data.table::setnames(x, "topic_num", topic_column)
    x[]
  }
  qc$optimized_counts <- update_primary_counts(
    qc$optimized_counts,
    "optimized_topic"
  )
  identity_map <- identical(
    as.integer(optimization$raw_to_optimized),
    as.integer(names(optimization$raw_to_optimized))
  )
  if (identity_map) {
    qc$raw_counts <- update_primary_counts(qc$raw_counts, "raw_topic")
  }
  condition_primary <- unique(evidence[condition_primary_confident == TRUE, .(
    condition_id = as.character(condition_id),
    tf = as.character(tf),
    optimized_topic = as.integer(condition_primary_topic)
  )])[, .(
    primary_tfs = data.table::uniqueN(tf)
  ), by = .(condition_id, optimized_topic)]
  qc$condition_topic <- condition_primary[
    qc$condition_topic,
    on = c("condition_id", "optimized_topic")
  ]
  qc$condition_topic[, tfs := data.table::fcoalesce(
    as.integer(primary_tfs),
    0L
  )]
  qc$condition_topic[, primary_tfs := NULL]
  optimization$qc <- qc
  optimization
}

.write_module3_topic_assignment_qc <- function(optimization,
                                               out_file,
                                               title_prefix = NULL,
                                               condition_colors = NULL,
                                               tf_topic_evidence = NULL,
                                               gene_term_assignment = NULL,
                                               gene_umap_genes = NULL,
                                               peak_term_assignment = NULL,
                                               peak_umap_top_n = Inf,
                                               pathway_gene_sets = NULL,
                                               pathway_colors = NULL,
                                               tf_target_gene_sets = NULL,
                                               tf_target_panel = NULL,
                                               top_n_tfs = 150L,
                                               seed = 20260716L,
                                               sections = "standard",
                                               report_scope = c(
                                                 "auto",
                                                 "full",
                                                 "condition_correlation"
                                               ),
                                               sidebar_mode = c(
                                                 "auto",
                                                 "terms",
                                                 "links"
                                               )) {
  .assert_pkg("ggplot2")
  .assert_pkg("gridExtra")
  .assert_pkg("scales")
  gene_only_assignment <- (
    optimization$assignment_mode %||% "gene_peak"
  ) %in% c("gene", "gene_only")
  has_peak_terms <- !isTRUE(gene_only_assignment) &&
    any(startsWith(colnames(optimization$raw_phi), "PEAK:"))
  sections <- unique(as.character(sections))
  allowed_sections <- c(
    "gene_phi_umap",
    "peak_phi_umap",
    "pathway_gene_umap",
    "tf_target_gene_umap",
    "raw_structure",
    "condition_correlation"
  )
  standard_report <- identical(sections, "standard")
  full_report <- identical(sections, "all")
  if (!length(sections) ||
      (any(c("standard", "all") %in% sections) &&
       !standard_report && !full_report) ||
      any(!sections %in% c("standard", "all", allowed_sections))) {
    .log_abort(paste0(
      "Topic-assignment QC sections must be 'standard', 'all', or a subset of: ",
      paste(allowed_sections, collapse = ", "),
      "."
    ))
  }
  if (standard_report) {
    sections <- c(
      "gene_phi_umap",
      if (has_peak_terms) "peak_phi_umap",
      if (!is.null(tf_target_gene_sets)) "tf_target_gene_umap",
      "raw_structure",
      "condition_correlation"
    )
  }
  report_scope <- match.arg(report_scope)
  gene_clustering <- identical(
    optimization$qc_mode %||% "topic_model",
    "gene_clustering"
  )
  multimodal_clustering <- identical(
    optimization$qc_mode %||% "topic_model",
    "multimodal_clustering"
  )
  cluster_qc <- gene_clustering || multimodal_clustering
  if (identical(report_scope, "auto")) {
    report_scope <- if (identical(
      optimization$document_design %||% "condition_tf",
      "condition"
    )) {
      "condition_correlation"
    } else {
      "full"
    }
  }
  if ("pathway_gene_umap" %in% sections && is.null(pathway_gene_sets)) {
    .log_abort("pathway_gene_umap requires pathway_gene_sets.")
  }
  if ("tf_target_gene_umap" %in% sections && is.null(tf_target_gene_sets)) {
    .log_abort("tf_target_gene_umap requires tf_target_gene_sets.")
  }
  needs_gene_umap <- full_report || any(c(
    "gene_phi_umap",
    "pathway_gene_umap",
    "tf_target_gene_umap"
  ) %in% sections)
  needs_peak_umap <- "peak_phi_umap" %in% sections
  if (needs_gene_umap || needs_peak_umap) {
    .assert_pkg("uwot")
  }
  sidebar_mode <- match.arg(sidebar_mode)
  if (!is.null(tf_topic_evidence) && nrow(tf_topic_evidence)) {
    optimization <- .m3_qc_apply_primary_tf_counts(
      optimization,
      tf_topic_evidence
    )
  }
  qc <- optimization$qc
  title_prefix <- title_prefix %||% basename(dirname(out_file))
  topic_values <- sort(unique(c(
    as.integer(qc$raw_topic_ids),
    as.integer(qc$optimized_topic_ids %||% qc$raw_topic_ids)
  )))
  topic_labels <- paste0("Topic ", topic_values)
  topic_palette <- stats::setNames(
    .module3_bright_topic_palette(topic_labels),
    topic_labels
  )
  gene_umap_pages <- list()
  peak_umap_pages <- list()
  pathway_umap_pages <- list()
  tf_target_umap_pages <- list()
  gene_umap <- NULL
  if (needs_gene_umap) {
    gene_umap <- .m3_qc_gene_term_umap_data(
      optimization,
      gene_term_assignment = gene_term_assignment,
      gene_universe = gene_umap_genes,
      seed = seed + 1L
    )
  }
  if (full_report || "gene_phi_umap" %in% sections) {
    gene_umap_topic_chunks <- if (length(qc$raw_topic_ids) > 30L) {
      split(
        qc$raw_topic_ids,
        ceiling(seq_along(qc$raw_topic_ids) / 30L)
      )
    } else {
      list(qc$raw_topic_ids)
    }
    gene_umap_pages <- lapply(
      seq_along(gene_umap_topic_chunks),
      function(i) {
        page_title <- paste0(
          title_prefix,
          if (cluster_qc) ": Gene clusters" else ": Gene topics"
        )
        if (length(gene_umap_topic_chunks) > 1L) {
          page_title <- sprintf(
            "%s (%d/%d)",
            page_title,
            i,
            length(gene_umap_topic_chunks)
          )
        }
        .m3_qc_arrange(
          .m3_qc_gene_term_umap_plot(
            optimization,
            topic_palette = topic_palette,
            gene_term_assignment = gene_term_assignment,
            gene_umap = gene_umap,
            topic_ids = gene_umap_topic_chunks[[i]],
            seed = seed + 1L
          ),
          ncol = 1L,
          title = page_title
        )
      }
    )
  }
  if (needs_peak_umap) {
    peak_umap <- .m3_qc_peak_term_umap_data(
      optimization,
      peak_term_assignment = peak_term_assignment,
      top_n = peak_umap_top_n,
      seed = seed + 2L
    )
    peak_umap_topic_chunks <- if (length(qc$raw_topic_ids) > 30L) {
      split(
        qc$raw_topic_ids,
        ceiling(seq_along(qc$raw_topic_ids) / 30L)
      )
    } else {
      list(qc$raw_topic_ids)
    }
    peak_umap_pages <- lapply(
      seq_along(peak_umap_topic_chunks),
      function(i) {
        peak_scope <- if (is.infinite(peak_umap_top_n)) {
          NULL
        } else {
          sprintf(
            "top %s by maximum phi probability",
            scales::comma(nrow(peak_umap))
          )
        }
        page_title <- if (is.null(peak_scope)) {
          sprintf(
            "%s: Peak topics (GammaFit + max-phi colors)",
            title_prefix
          )
        } else {
          sprintf(
            "%s: Peak topics (%s; GammaFit + max-phi colors)",
            title_prefix,
            peak_scope
          )
        }
        if (length(peak_umap_topic_chunks) > 1L) {
          page_title <- sprintf(
            "%s (%d/%d)",
            page_title,
            i,
            length(peak_umap_topic_chunks)
          )
        }
        .m3_qc_arrange(
          .m3_qc_peak_term_umap_plot(
            optimization,
            topic_palette = topic_palette,
            peak_umap = peak_umap,
            topic_ids = peak_umap_topic_chunks[[i]],
            top_n = peak_umap_top_n,
            seed = seed + 2L
          ),
          ncol = 1L,
          title = page_title
        )
      }
    )
  }
  if (!is.null(pathway_gene_sets) &&
      (full_report || "pathway_gene_umap" %in% sections)) {
    pathway_umap_pages <- list(.m3_qc_arrange(
      .m3_qc_pathway_gene_umap_plot(
        gene_umap = gene_umap,
        pathway_gene_sets = pathway_gene_sets,
        pathway_colors = pathway_colors
      ),
      ncol = 1L,
      title = paste0(title_prefix, ": Pathway genes")
    ))
  }
  if (!is.null(tf_target_gene_sets) &&
      (full_report || "tf_target_gene_umap" %in% sections)) {
    displayed_tf_target_gene_sets <- if (all(c(
      "display_rank", "selected_topic"
    ) %in% names(tf_target_gene_sets))) {
      tf_target_gene_sets
    } else {
      .m3_qc_select_topic_tf_target_gene_sets(
        tf_target_gene_sets = tf_target_gene_sets,
        gene_term_assignment = gene_term_assignment,
        tf_panel = tf_target_panel,
        top_n_tfs = 30L,
        top_n_targets = 500L
      )
    }
    tf_target_umap_pages <- list(.m3_qc_arrange(
      .m3_qc_tf_target_gene_umap_plot(
        gene_umap = gene_umap,
        tf_target_gene_sets = displayed_tf_target_gene_sets
      ),
      ncol = 1L,
      title = paste0(title_prefix, ": High-confidence TF targets")
    ))
  }

  term_sidebar <- identical(sidebar_mode, "terms") ||
    (identical(sidebar_mode, "auto") && identical(
      optimization$document_design %||% "condition_tf",
      "condition"
    ))
  sidebar_columns <- if (gene_only_assignment) {
    "genes"
  } else if (term_sidebar) {
    c("peaks", "genes")
  } else {
    c("links", "genes", "tfs")
  }
  sidebar_labels <- if (gene_only_assignment) {
    "Target genes"
  } else if (term_sidebar) {
    c("Peaks", "Target genes")
  } else {
    c(
      "Links",
      "Target genes",
      if (!is.null(tf_topic_evidence) && nrow(tf_topic_evidence)) {
        "Confident primary TFs"
      } else {
        "TFs"
      }
    )
  }
  sidebar_colors <- if (gene_only_assignment) {
    "#D97824"
  } else if (term_sidebar) {
    c("#007C78", "#D97824")
  } else {
    c("#007C78", "#D97824", "#6A5ACD")
  }
  raw_topic_structure_plot <- NULL
  if (full_report || "raw_structure" %in% sections) {
    raw_topic_similarity <- qc$raw_topic_similarity
    rownames(raw_topic_similarity) <- colnames(raw_topic_similarity) <-
      paste0("Topic ", qc$raw_topic_ids)
    raw_topic_order <- .m3_qc_cluster_order(raw_topic_similarity, "row")
    structure_title <- if (gene_clustering) {
      "Gene-cluster expression structure"
    } else if (multimodal_clustering) {
      "Multimodal cluster activity structure"
    } else {
      "Raw topic assignment structure"
    }
    structure_subtitle <- if (gene_clustering) {
      paste0(
        "Hellinger similarity of cluster-level expression profiles across ",
        "conditions; bars show assigned target genes"
      )
    } else if (multimodal_clustering) {
      paste0(
        "Hellinger similarity of cluster-level Gene and Peak activity ",
        "profiles across conditions; bars show assigned terms"
      )
    } else {
      paste0(
        "Mean Hellinger similarity across separately normalized ",
        if (identical(optimization$assignment_mode, "tf_target")) {
          "Gene and TF-target phi"
        } else if (gene_only_assignment) {
          "Gene phi"
        } else {
          "Gene and Peak phi"
        },
        "; bars show full-universe assigned counts"
      )
    }
    raw_topic_structure_plot <- .m3_qc_topic_structure_plot(
      similarity = raw_topic_similarity,
      counts = qc$raw_counts,
      topic_column = "raw_topic",
      topic_order = raw_topic_order,
      title = structure_title,
      subtitle = structure_subtitle,
      sidebar_columns = sidebar_columns,
      sidebar_labels = sidebar_labels,
      sidebar_colors = sidebar_colors
    )
  }

  compact_condition_correlation_plot <- NULL
  compact_condition_correlation_title <- "Condition correlation by topic composition"
  if (!full_report && "condition_correlation" %in% sections) {
    condition_document <- identical(
      optimization$document_design %||% "condition_tf",
      "condition"
    )
    condition_theta <- if (condition_document) {
      .m3_qc_condition_topic_theta(optimization)
    } else {
      .m3_qc_condition_tf_topic_theta(
        optimization$theta,
        matched_tfs = TRUE
      )
    }
    condition_correlation <- .m3_qc_theta_condition_correlation(
      condition_theta
    )
    display_conditions <- .m3_qc_short_condition_labels(
      rownames(condition_correlation)
    )
    rownames(condition_correlation) <- colnames(condition_correlation) <-
      display_conditions
    condition_correlation_order <- .m3_qc_cluster_order(
      condition_correlation,
      "row"
    )
    condition_subtitle <- if (condition_document) {
      "Pearson correlation of raw condition-document theta profiles"
    } else if (isTRUE(attr(condition_theta, "matched_tfs"))) {
      sprintf(
        "Pearson correlation of mean raw theta across %d TF documents shared by all conditions",
        as.integer(attr(condition_theta, "tf_count"))
      )
    } else {
      "Pearson correlation of mean raw theta across available TF documents"
    }
    compact_condition_correlation_plot <- .m3_qc_correlation_plot(
      condition_correlation,
      compact_condition_correlation_title,
      condition_subtitle,
      order = condition_correlation_order
    )
  }

  if (!full_report) {
    pages <- c(
      gene_umap_pages,
      peak_umap_pages,
      pathway_umap_pages,
      tf_target_umap_pages
    )
    if ("raw_structure" %in% sections) {
      pages[[length(pages) + 1L]] <- .m3_qc_arrange(
        raw_topic_structure_plot,
        ncol = 1L,
        title = structure_title
      )
    }
    if ("condition_correlation" %in% sections) {
      pages[[length(pages) + 1L]] <- .m3_qc_arrange(
        compact_condition_correlation_plot,
        ncol = 1L,
        title = compact_condition_correlation_title
      )
    }
    if (standard_report) {
      pages[[length(pages) + 1L]] <- .m3_qc_retention_page(
        optimization,
        title_prefix = title_prefix,
        standardized = TRUE
      )
    }
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::cairo_pdf(
      out_file,
      width = 8.5,
      height = 11,
      family = "Helvetica",
      onefile = TRUE
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    for (page in pages) {
      grid::grid.newpage()
      grid::grid.draw(page)
    }
    return(invisible(out_file))
  }

  assignments <- qc$assignments
  identity_map <- identical(
    as.integer(optimization$raw_to_optimized),
    as.integer(names(optimization$raw_to_optimized))
  )
  optimized_aligned_sample <- NULL
  if (!identity_map) {
    sample_assignments <- assignments[optimization$sample_rows]
    optimized_coordinates <- optimization$optimized_sample_coordinates %||%
      .m3_qc_umap_coordinates(
        optimization$optimized_sample_probability,
        seed = seed
      )
    optimized_sample <- data.table::data.table(
      UMAP1 = optimized_coordinates[, 1L],
      UMAP2 = optimized_coordinates[, 2L],
      condition_id = sample_assignments$condition_id,
      tf = qc$tf_levels[sample_assignments$tf_index],
      target_gene = qc$target_levels[sample_assignments$target_index],
      topic = data.table::fifelse(
        sample_assignments$optimized_aligned,
        paste0("Topic ", sample_assignments$optimized_topic),
        NA_character_
      )
    )
    optimized_sample[, topic_short := sub("^Topic ", "T", topic)]
    optimized_sample <- .m3_qc_jitter_umap(
      optimized_sample,
      seed = seed + 21L
    )
    optimized_aligned_sample <- optimized_sample[!is.na(topic)]
  }
  condition_document <- identical(
    optimization$document_design %||% "condition_tf",
    "condition"
  )
  condition_theta_for_correlation <- if (condition_document) {
    .m3_qc_condition_topic_theta(optimization)
  } else {
    NULL
  }
  condition_values <- if (condition_document) {
    rownames(condition_theta_for_correlation)
  } else {
    sort(unique(assignments$condition_id))
  }

  optimized_topic_similarity <- qc$optimized_topic_similarity
  rownames(optimized_topic_similarity) <- colnames(optimized_topic_similarity) <-
    paste0("Topic ", qc$optimized_topic_ids)
  optimized_topic_order <- .m3_qc_cluster_order(
    optimized_topic_similarity,
    "row"
  )
  page3 <- if (!identity_map) {
    .m3_qc_arrange(
      .m3_qc_umap_plot(
        optimized_aligned_sample,
        "topic",
        "Optimized filtered aligned-link UMAP",
        "Sampled links retained after topic merge and alignment",
        colors = topic_palette,
        label_column = "topic_short",
        seed = seed + 4L
      ),
      .m3_qc_topic_structure_plot(
        similarity = optimized_topic_similarity,
        counts = qc$optimized_counts,
        topic_column = "optimized_topic",
        topic_order = optimized_topic_order,
        title = "Optimized topic assignment structure",
        subtitle = paste(
          "Mean Hellinger similarity after deterministic topic merging;",
          "bars show full-universe assigned counts"
        ),
        sidebar_columns = sidebar_columns,
        sidebar_labels = sidebar_labels,
        sidebar_colors = sidebar_colors
      ),
      ncol = 1L,
      heights = c(0.85, 1.15),
      title = "Optimized topic assignments"
    )
  } else {
    NULL
  }

  link_matrix <- .m3_qc_condition_topic_matrix(
    qc$condition_topic,
    "links",
    condition_values,
    qc$optimized_topic_ids
  )
  gene_matrix <- .m3_qc_condition_topic_matrix(
    qc$condition_topic,
    "genes",
    condition_values,
    qc$optimized_topic_ids
  )
  tf_matrix <- .m3_qc_condition_topic_matrix(
    qc$condition_topic,
    "tfs",
    condition_values,
    qc$optimized_topic_ids
  )
  display_conditions <- .m3_qc_short_condition_labels(rownames(link_matrix))
  rownames(link_matrix) <- display_conditions
  rownames(gene_matrix) <- display_conditions
  rownames(tf_matrix) <- display_conditions
  condition_correlation <- if (condition_document) {
    .m3_qc_theta_condition_correlation(condition_theta_for_correlation)
  } else {
    .m3_qc_condition_correlation(link_matrix, gene_matrix)
  }
  rownames(condition_correlation) <- colnames(condition_correlation) <-
    display_conditions
  condition_correlation_order <- .m3_qc_cluster_order(
    condition_correlation,
    "row"
  )
  condition_correlation_title <- if (gene_clustering) {
    "Condition correlation by gene-cluster activity"
  } else if (multimodal_clustering) {
    "Condition correlation by multimodal-cluster activity"
  } else {
    "Condition correlation by topic composition"
  }
  condition_correlation_subtitle <- if (gene_clustering) {
    "Pearson correlation of normalized expression mass across hard Gene clusters"
  } else if (multimodal_clustering) {
    paste0(
      "Pearson correlation of normalized Gene and Peak activity across ",
      "hard multimodal clusters"
    )
  } else if (condition_document) {
    "Pearson correlation of raw condition-document theta profiles"
  } else {
    paste(
      "Pearson correlation of separately normalized TF-target-link",
      "and target-gene topic profiles",
      sep = "\n"
    )
  }
  condition_correlation_plot <- .m3_qc_correlation_plot(
    condition_correlation,
    condition_correlation_title,
    condition_correlation_subtitle,
    order = condition_correlation_order
  )
  structure_pages <- if (
    max(length(condition_values), length(qc$raw_topic_ids)) > 24L
  ) {
    list(
      .m3_qc_arrange(
        raw_topic_structure_plot,
        ncol = 1L,
        title = structure_title
      ),
      .m3_qc_arrange(
        condition_correlation_plot,
        ncol = 1L,
        title = condition_correlation_title
      )
    )
  } else {
    list(.m3_qc_arrange(
      raw_topic_structure_plot,
      condition_correlation_plot,
      ncol = 1L,
      heights = c(1.1, 0.9),
      title = if (cluster_qc) {
        "Cluster and condition activity structure"
      } else {
        "Topic and condition assignment structure"
      }
    ))
  }
  if (identical(report_scope, "condition_correlation")) {
    pages <- c(
      gene_umap_pages,
      peak_umap_pages,
      pathway_umap_pages,
      tf_target_umap_pages,
      structure_pages
    )
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::cairo_pdf(
      out_file,
      width = 8.5,
      height = 11,
      family = "Helvetica",
      onefile = TRUE
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    for (page in pages) {
      grid::grid.newpage()
      grid::grid.draw(page)
    }
    return(invisible(out_file))
  }
  clustering_matrix <- cbind(
    .m3_opt_row_normalize(link_matrix),
    .m3_opt_row_normalize(gene_matrix),
    .m3_opt_row_normalize(tf_matrix)
  )
  row_order <- .m3_qc_cluster_order(clustering_matrix, "row")
  column_order <- .m3_qc_cluster_order(link_matrix + gene_matrix, "column")
  condition_theta <- .m3_qc_condition_topic_theta(optimization)
  condition_theta <- condition_theta[
    condition_values,
    paste0("Topic ", qc$optimized_topic_ids),
    drop = FALSE
  ]
  expression_matrix <- .m3_qc_document_expression_topic_share_matrix(
    optimization = optimization,
    conditions = condition_values,
    topics = qc$optimized_topic_ids
  )
  rownames(condition_theta) <- display_conditions
  rownames(expression_matrix) <- display_conditions
  expression_values <- as.numeric(expression_matrix)
  expression_values <- expression_values[is.finite(expression_values)]
  expression_limits <- if (length(unique(expression_values)) > 1L) {
    as.numeric(stats::quantile(
      expression_values,
      probs = c(0.05, 0.95),
      na.rm = TRUE,
      names = FALSE
    ))
  } else {
    NULL
  }
  if (!is.null(expression_limits) && diff(expression_limits) <= 0) {
    expression_limits <- NULL
  }
  profile_topic_chunks <- if (ncol(link_matrix) > 24L) {
    split(column_order, ceiling(seq_along(column_order) / 20L))
  } else {
    list(column_order)
  }
  profile_plot_groups <- lapply(profile_topic_chunks, function(selected) {
    selected_order <- seq_along(selected)
    .m3_qc_align_plot_dimensions(list(
      .m3_qc_count_heatmap(
        link_matrix[, selected, drop = FALSE],
        "Unique aligned TF-target links",
        row_order = row_order,
        column_order = selected_order,
        square_cells = FALSE,
        show_labels = FALSE
      ) + ggplot2::labs(x = "Topic"),
      .m3_qc_count_heatmap(
        tf_matrix[, selected, drop = FALSE],
        "Unique TFs",
        row_order = row_order,
        column_order = selected_order,
        square_cells = FALSE,
        show_labels = FALSE
      ) + ggplot2::labs(x = "Topic"),
      .m3_qc_value_heatmap(
        condition_theta[, selected, drop = FALSE],
        "Mean condition topic probability",
        "Mean theta",
        row_order = row_order,
        column_order = selected_order
      ) + ggplot2::labs(x = "Topic"),
      .m3_qc_value_heatmap(
        expression_matrix[, selected, drop = FALSE],
        "Mean expression of assigned target genes (topic share)",
        "Mean expression share",
        row_order = row_order,
        column_order = selected_order,
        limits = expression_limits,
        color_transform = "identity"
      ) + ggplot2::labs(x = "Topic")
    ))
  })
  plot_chunks <- if (length(condition_values) > 22L) {
    as.list(seq_len(4L))
  } else if (length(condition_values) > 12L) {
    list(1:2, 3:4)
  } else {
    list(seq_len(4L))
  }
  profile_page_specs <- unlist(lapply(
    seq_along(profile_plot_groups),
    function(group_index) {
      lapply(plot_chunks, function(indices) {
        list(group = group_index, indices = indices)
      })
    }
  ), recursive = FALSE)
  profile_pages <- lapply(seq_along(profile_page_specs), function(page_index) {
    spec <- profile_page_specs[[page_index]]
    page_title <- if (length(profile_page_specs) > 1L) {
      sprintf(
        "Condition-topic assignment profiles (%d/%d)",
        page_index,
        length(profile_page_specs)
      )
    } else {
      "Condition-topic assignment profiles"
    }
    do.call(.m3_qc_arrange, c(
      profile_plot_groups[[spec$group]][spec$indices],
      list(ncol = 1L, title = page_title)
    ))
  })

  cross <- qc$condition_topic_similarity$mean
  cross_display <- cross
  rownames(cross_display) <- .m3_qc_short_condition_labels(rownames(cross_display))
  cross_plot <- .m3_qc_similarity_plot(
    cross_display,
    "Condition-topic similarity",
    "Mean Jaccard of unique TF-target pairs and target genes",
    autoscale = TRUE
  )
  page5 <- .m3_qc_arrange(
    cross_plot,
    ncol = 1L,
    title = "Condition-topic similarity"
  )
  retention_page <- .m3_qc_retention_page(optimization)

  tf_pages <- .m3_qc_tf_topic_pages(
    optimization,
    top_n_tfs = top_n_tfs
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(
    out_file,
    width = 8.5,
    height = 11,
    family = "Helvetica",
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page in gene_umap_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  for (page in peak_umap_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  for (page in pathway_umap_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  for (page in tf_target_umap_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  for (page in structure_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  if (!identity_map) {
    grid::grid.newpage()
    grid::grid.draw(page3)
  }
  for (page in profile_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  grid::grid.newpage()
  grid::grid.draw(page5)
  for (page in tf_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
  grid::grid.newpage()
  grid::grid.draw(retention_page)
  invisible(out_file)
}

.write_module3_topic_optimization_outputs <- function(optimization,
                                                      out_dir,
                                                      raw_k,
                                                      title_prefix = NULL,
                                                      condition_colors = NULL,
                                                      top_n_tfs = 150L,
                                                      seed = 20260716L) {
  mapping <- data.table::data.table(
    raw_topic = as.integer(names(optimization$raw_to_optimized)),
    optimized_topic = as.integer(optimization$raw_to_optimized)
  )
  data.table::fwrite(
    mapping,
    file.path(out_dir, "topic_optimization_map.csv")
  )
  data.table::fwrite(
    optimization$merge_audit,
    file.path(out_dir, "topic_optimization_merge_audit.csv")
  )
  summary <- data.table::data.table(
    assignment_mode = optimization$assignment_mode %||% "gene_peak",
    require_theta_gate = isTRUE(optimization$require_theta_gate %||% TRUE),
    raw_topics = length(unique(mapping$raw_topic)),
    optimized_topics = length(unique(mapping$optimized_topic)),
    merged_topics = length(unique(mapping$raw_topic)) -
      length(unique(mapping$optimized_topic)),
    input_links = nrow(unique(optimization$qc$assignments[, .(
      condition_id, tf_index, target_index
    )])),
    raw_aligned_links = nrow(unique(
      optimization$qc$assignments[raw_aligned == TRUE, .(
        condition_id, tf_index, target_index
      )]
    )),
    optimized_aligned_links = nrow(unique(
      optimization$qc$assignments[optimized_aligned == TRUE, .(
        condition_id, tf_index, target_index
      )]
    )),
    raw_pair_assigned_genes = if (identical(
      optimization$assignment_mode,
      "tf_target"
    )) {
      data.table::uniqueN(optimization$raw_pair_assignment[
        .as_logical_flag(assigned),
        target_gene
      ])
    } else {
      sum(.as_logical_flag(optimization$raw_pair_assignment$assigned))
    },
    optimized_pair_assigned_genes = if (identical(
      optimization$assignment_mode,
      "tf_target"
    )) {
      data.table::uniqueN(optimization$pair_assignment[
        .as_logical_flag(assigned),
        target_gene
      ])
    } else {
      sum(.as_logical_flag(optimization$pair_assignment$assigned))
    },
    raw_assigned_tf_target_terms = if (identical(
      optimization$assignment_mode,
      "tf_target"
    )) {
      sum(.as_logical_flag(optimization$raw_pair_assignment$assigned))
    } else {
      NA_integer_
    },
    optimized_assigned_tf_target_terms = if (identical(
      optimization$assignment_mode,
      "tf_target"
    )) {
      sum(.as_logical_flag(optimization$pair_assignment$assigned))
    } else {
      NA_integer_
    },
    raw_assigned_genes = data.table::uniqueN(
      optimization$qc$assignments[raw_aligned == TRUE, target_index]
    ),
    optimized_assigned_genes = data.table::uniqueN(
      optimization$qc$assignments[optimized_aligned == TRUE, target_index]
    ),
    raw_tf_term_assignments =
      optimization$raw_tf_theta_correspondence$tf_term_assignments,
    raw_tf_theta_supported =
      optimization$raw_tf_theta_correspondence$supported_tf_terms,
    raw_tf_theta_empty =
      optimization$raw_tf_theta_correspondence$empty_tf_terms,
    optimized_tf_term_assignments =
      optimization$optimized_tf_theta_correspondence$tf_term_assignments,
    optimized_tf_theta_supported =
      optimization$optimized_tf_theta_correspondence$supported_tf_terms,
    optimized_tf_theta_empty =
      optimization$optimized_tf_theta_correspondence$empty_tf_terms,
    tf_topic_cutoff = optimization$tf_topic_cutoff,
    umap_sampled_links = length(optimization$sample_rows)
  )
  data.table::fwrite(
    summary,
    file.path(out_dir, "topic_optimization_summary.csv")
  )
  settings <- data.table::data.table(
    raw_k = as.integer(raw_k),
    assignment_mode = optimization$assignment_mode %||% "gene_peak",
    min_genes = as.integer(optimization$min_genes %||% NA_integer_),
    min_links = as.integer(optimization$min_links %||% NA_integer_),
    similarity_threshold = as.numeric(optimization$similarity_threshold %||% NA_real_),
    prefer_tf_theta_correspondence = isTRUE(
      optimization$prefer_tf_theta_correspondence
    ),
    tf_theta_merge_aggregation = "max",
    tf_topic_cutoff = as.numeric(optimization$tf_topic_cutoff),
    umap_max_links_per_condition = as.integer(
      optimization$umap_max_links_per_condition %||% length(optimization$sample_rows)
    ),
    qc_top_tfs = as.integer(top_n_tfs),
    upregulation_reference_condition = as.character(
      optimization$upregulation_reference_condition %||% NA_character_
    ),
    upregulated_log2fc_min = as.numeric(
      optimization$upregulated_log2fc_min %||% NA_real_
    ),
    qc_seed = as.integer(seed)
  )
  data.table::fwrite(
    settings,
    file.path(out_dir, "topic_assignment_qc_settings.csv")
  )
  data.table::fwrite(
    data.table::data.table(
      doc_index = seq_len(nrow(optimization$raw_theta)),
      doc_id = rownames(optimization$raw_theta)
    ),
    file.path(out_dir, "topic_assignment_document_lookup.csv")
  )
  data.table::fwrite(
    data.table::data.table(
      target_index = seq_along(optimization$qc$target_levels),
      target_gene = optimization$qc$target_levels
    ),
    file.path(out_dir, "topic_assignment_target_lookup.csv")
  )
  .save_all(out_dir, "topic_theta_optimized", optimization$theta)
  .save_all(out_dir, "topic_phi_optimized", optimization$phi)
  .save_all(out_dir, "topic_term_scores_optimized", optimization$score)
  assignments <- optimization$qc$assignments[, .(
    doc_index,
    pair_index = if ("pair_index" %in% names(optimization$qc$assignments)) {
      pair_index
    } else {
      NA_integer_
    },
    target_index,
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
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(
      assignments,
      file.path(out_dir, "topic_link_assignments.parquet"),
      compression = "zstd"
    )
  } else {
    saveRDS(
      assignments,
      file.path(out_dir, "topic_link_assignments.rds"),
      compress = "gzip"
    )
  }
  .write_module3_topic_assignment_qc(
    optimization,
    out_file = file.path(
      out_dir,
      sprintf("topic_assignment_qc_K%d.pdf", as.integer(raw_k))
    ),
    title_prefix = title_prefix,
    condition_colors = condition_colors,
    top_n_tfs = top_n_tfs,
    seed = seed
  )
  invisible(summary)
}
