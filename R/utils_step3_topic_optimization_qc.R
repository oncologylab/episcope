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
      target_index = as.integer(entries$j),
      gene_token_count = as.numeric(entries$x),
      peak_token_count = peak_count,
      condition_id = doc_meta$condition_id[entries$i]
    ),
    docs = doc_meta,
    pairs = pairs,
    gene_columns = gene_columns,
    peak_columns = peak_columns
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
    theta_chunk <- theta[link_rows$doc_index, , drop = FALSE]
    gene_phi <- t(phi[, universe$gene_columns[link_rows$target_index], drop = FALSE])
    peak_phi <- t(phi[, universe$peak_columns[link_rows$target_index], drop = FALSE])
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
                              min_genes = 50L,
                              min_links = 200L,
                              similarity_threshold = 0.90,
                              min_topics = 2L) {
  mapping <- stats::setNames(raw_topic_ids, raw_topic_ids)
  audit <- list()
  step <- 0L
  repeat {
    group_phi <- .m3_opt_group_phi(phi, theta, dtm, raw_topic_ids, mapping)
    groups <- .m3_opt_topic_ids(group_phi, "row")
    links <- vapply(groups, function(g) sum(raw_links[mapping == g]), numeric(1L))
    genes <- vapply(groups, function(g) sum(raw_genes[mapping == g]), numeric(1L))
    similarity <- .m3_opt_hellinger_similarity(group_phi, gene_ids, peak_ids)
    diag(similarity) <- -Inf
    small <- which(links < min_links | genes < min_genes)
    reason <- NULL
    if (length(groups) <= min_topics) break
    if (length(small)) {
      small <- small[order(links[small], genes[small], groups[small])]
      nearest <- vapply(small, function(source) {
        which.max(similarity[source, ])
      }, integer(1L))
      nearest_similarity <- similarity[cbind(small, nearest)]
      eligible <- which(
        is.finite(nearest_similarity) &
          nearest_similarity >= similarity_threshold
      )
      if (!length(eligible)) break
      source_pos <- small[eligible[[1L]]]
      target_pos <- nearest[eligible[[1L]]]
      reason <- data.table::fcase(
        links[[source_pos]] < min_links & genes[[source_pos]] < min_genes, "small_links_and_genes",
        links[[source_pos]] < min_links, "small_links",
        default = "small_genes"
      )
    } else {
      pair <- which(similarity == max(similarity), arr.ind = TRUE)
      pair <- pair[pair[, 1L] < pair[, 2L], , drop = FALSE]
      if (!nrow(pair) || similarity[pair[1L, 1L], pair[1L, 2L]] < similarity_threshold) break
      source_pos <- pair[1L, 1L]
      target_pos <- pair[1L, 2L]
      reason <- "high_similarity"
    }
    source_group <- groups[[source_pos]]
    target_group <- groups[[target_pos]]
    source_size <- c(links[[source_pos]], genes[[source_pos]])
    target_size <- c(links[[target_pos]], genes[[target_pos]])
    if (source_size[[1L]] > target_size[[1L]] ||
        (source_size[[1L]] == target_size[[1L]] && source_size[[2L]] > target_size[[2L]]) ||
        (identical(source_size, target_size) && source_group < target_group)) {
      tmp <- source_group
      source_group <- target_group
      target_group <- tmp
      tmp <- source_pos
      source_pos <- target_pos
      target_pos <- tmp
    }
    step <- step + 1L
    audit[[step]] <- data.table::data.table(
      merge_step = step,
      source_topic = source_group,
      representative_topic = target_group,
      reason = reason,
      similarity = similarity[source_pos, target_pos],
      source_links = links[source_pos],
      source_genes = genes[source_pos],
      representative_links_before = links[target_pos],
      representative_genes_before = genes[target_pos]
    )
    mapping[mapping == source_group] <- target_group
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
        representative_genes_before = numeric()
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
                                optimized_pairs) {
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
  assigned_target <- stats::setNames(
    optimized_pairs$optimized_assigned_topic,
    optimized_pairs$target_gene
  )
  collapsed[, target_key__ := data.table::fcase(
    term_group == "GENE", sub("^GENE:", "", term_id),
    term_group == "PEAK", sub("^PEAK:", "", term_id),
    default = NA_character_
  )]
  collapsed[, in_topic := data.table::fcase(
    term_group %in% c("GENE", "PEAK"),
    is.finite(assigned_target[target_key__]) &
      topic_num == assigned_target[target_key__],
    default = gammafit_candidate
  )]
  collapsed[, `:=`(
    topic = paste0("Topic", topic_num),
    assignment_method = data.table::fifelse(
      term_group %in% c("GENE", "PEAK"),
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
  collapsed[, c("target_key__", "optimized_score") := NULL]
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
    target_levels = pairs$target_gene
  )
}

.module3_optimize_condition_topics <- function(theta,
                                               phi,
                                               dtm,
                                               topic_terms,
                                               pair_assignment,
                                               assignment_mode = c("gene_peak", "gene_only"),
                                               condition_gene_expression = NULL,
                                               condition_upregulated_genes = NULL,
                                               upregulation_reference_condition = NULL,
                                               upregulated_log2fc_min = 1,
                                               min_genes = 50L,
                                               min_links = 200L,
                                               similarity_threshold = 0.90,
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
  universe <- .m3_opt_link_universe(dtm, pair_assignment)
  raw_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    universe,
    chunk_size = chunk_size
  )
  raw_target_topic <- suppressWarnings(as.integer(universe$pairs$assigned_topic))
  target_topic_for_link <- raw_target_topic[universe$links$target_index]
  raw_posterior_agrees <- is.finite(target_topic_for_link) &
    raw_posterior$topic == target_topic_for_link
  raw_target_position <- match(target_topic_for_link, raw_topic_ids)
  raw_theta_pass <- rep(FALSE, nrow(universe$links))
  valid_raw_theta <- is.finite(raw_target_position)
  raw_theta_pass[valid_raw_theta] <- theta[cbind(
    universe$links$doc_index[valid_raw_theta],
    raw_target_position[valid_raw_theta]
  )] >= tf_topic_cutoff
  raw_aligned <- raw_posterior_agrees & raw_theta_pass
  raw_links <- vapply(raw_topic_ids, function(topic) {
    sum(raw_aligned & target_topic_for_link == topic)
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
    min_genes = min_genes,
    min_links = min_links,
    similarity_threshold = similarity_threshold
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
  optimized_posterior <- .m3_opt_link_posteriors(
    theta,
    phi,
    universe,
    group_index = merge$mapping,
    chunk_size = chunk_size
  )
  optimized_target_topic <- optimized_pairs$optimized_assigned_topic[
    universe$links$target_index
  ]
  optimized_posterior_agrees <- is.finite(optimized_target_topic) &
    optimized_posterior$topic == optimized_target_topic
  optimized_target_position <- match(optimized_target_topic, optimized_topic_ids)
  optimized_theta_pass <- rep(FALSE, nrow(universe$links))
  valid_optimized_theta <- is.finite(optimized_target_position)
  optimized_theta_pass[valid_optimized_theta] <- optimized_theta[cbind(
    universe$links$doc_index[valid_optimized_theta],
    optimized_target_position[valid_optimized_theta]
  )] >= tf_topic_cutoff
  optimized_aligned <- optimized_posterior_agrees & optimized_theta_pass
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
    optimized_pairs = optimized_pairs
  )
  optimized_similarity <- .m3_opt_hellinger_similarity(
    optimized_phi,
    gene_ids,
    peak_ids
  )
  qc <- .m3_opt_qc_tables(
    assignments = assignments,
    docs = universe$docs,
    pairs = optimized_pairs,
    raw_similarity = raw_similarity,
    optimized_similarity = optimized_similarity,
    raw_topic_ids = raw_topic_ids,
    optimized_topic_ids = optimized_topic_ids
  )
  list(
    theta = optimized_theta,
    phi = optimized_phi,
    score = optimized_score,
    topic_terms = optimized_terms,
    pair_assignment = optimized_pairs,
    raw_theta = theta,
    raw_phi = phi,
    raw_topic_terms = topic_terms,
    raw_pair_assignment = pair_assignment,
    assignment_mode = assignment_mode,
    condition_gene_expression = condition_gene_expression,
    condition_upregulated_genes = condition_upregulated_genes,
    upregulation_reference_condition = upregulation_reference_condition,
    upregulated_log2fc_min = upregulated_log2fc_min,
    raw_to_optimized = stats::setNames(merge$mapping, raw_topic_ids),
    merge_audit = merge$audit,
    sample_rows = sample_rows,
    raw_sample_probability = raw_sample_posterior$sample_probability,
    optimized_sample_probability = optimized_sample_posterior$sample_probability,
    tf_topic_cutoff = tf_topic_cutoff,
    min_genes = min_genes,
    min_links = min_links,
    similarity_threshold = similarity_threshold,
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
                                        subtitle) {
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

  count_data <- data.table::copy(counts)
  data.table::setnames(count_data, topic_column, "topic")
  count_data[, y := k - match(as.integer(topic), ordered_ids) + 1L]
  count_data <- count_data[is.finite(y)]
  bar_width <- max(3.5, min(6, k * 0.20))
  gap <- max(1, k * 0.04)
  link_start <- k + gap
  gene_start <- link_start + bar_width + gap
  tf_start <- gene_start + bar_width + gap
  link_max <- max(count_data$links, na.rm = TRUE)
  gene_max <- max(count_data$genes, na.rm = TRUE)
  tf_max <- max(count_data$tfs, na.rm = TRUE)
  if (!is.finite(link_max) || link_max <= 0) link_max <- 1
  if (!is.finite(gene_max) || gene_max <= 0) gene_max <- 1
  if (!is.finite(tf_max) || tf_max <= 0) tf_max <- 1
  count_data[, `:=`(
    link_xmax = link_start + links / link_max * bar_width,
    gene_xmax = gene_start + genes / gene_max * bar_width,
    tf_xmax = tf_start + tfs / tf_max * bar_width,
    link_label = .m3_qc_compact_count(links),
    gene_label = .m3_qc_compact_count(genes),
    tf_label = .m3_qc_compact_count(tfs),
    link_inside = links / link_max >= 0.28,
    gene_inside = genes / gene_max >= 0.28,
    tf_inside = tfs / tf_max >= 0.28
  )]
  right_edge <- tf_start + bar_width + 1.1
  title_y <- k + 1.15

  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = heatmap,
      ggplot2::aes(x, y, fill = similarity),
      colour = "#D6DCE1",
      linewidth = 0.2,
      width = 1,
      height = 1
    ) +
    ggplot2::geom_rect(
      data = count_data,
      ggplot2::aes(
        xmin = link_start,
        xmax = link_xmax,
        ymin = y - 0.40,
        ymax = y + 0.40
      ),
      fill = "#007C78",
      colour = "#30383F",
      linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = count_data,
      ggplot2::aes(
        xmin = gene_start,
        xmax = gene_xmax,
        ymin = y - 0.40,
        ymax = y + 0.40
      ),
      fill = "#D97824",
      colour = "#30383F",
      linewidth = 0.2
    ) +
    ggplot2::geom_rect(
      data = count_data,
      ggplot2::aes(
        xmin = tf_start,
        xmax = tf_xmax,
        ymin = y - 0.40,
        ymax = y + 0.40
      ),
      fill = "#6A5ACD",
      colour = "#30383F",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      data = count_data[link_inside == TRUE],
      ggplot2::aes(link_xmax - 0.12, y, label = link_label),
      hjust = 1,
      colour = "white",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::geom_text(
      data = count_data[link_inside == FALSE],
      ggplot2::aes(link_xmax + 0.12, y, label = link_label),
      hjust = 0,
      colour = "#20272E",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::geom_text(
      data = count_data[gene_inside == TRUE],
      ggplot2::aes(gene_xmax - 0.12, y, label = gene_label),
      hjust = 1,
      colour = "white",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::geom_text(
      data = count_data[gene_inside == FALSE],
      ggplot2::aes(gene_xmax + 0.12, y, label = gene_label),
      hjust = 0,
      colour = "#20272E",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::geom_text(
      data = count_data[tf_inside == TRUE],
      ggplot2::aes(tf_xmax - 0.12, y, label = tf_label),
      hjust = 1,
      colour = "white",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::geom_text(
      data = count_data[tf_inside == FALSE],
      ggplot2::aes(tf_xmax + 0.12, y, label = tf_label),
      hjust = 0,
      colour = "#20272E",
      family = "Helvetica",
      fontface = "bold",
      size = 3.0
    ) +
    ggplot2::annotate(
      "text",
      x = link_start + bar_width / 2,
      y = title_y,
      label = "Links",
      family = "Helvetica",
      fontface = "bold",
      size = 3.5
    ) +
    ggplot2::annotate(
      "text",
      x = gene_start + bar_width / 2,
      y = title_y,
      label = "Target genes",
      family = "Helvetica",
      fontface = "bold",
      size = 3.5
    ) +
    ggplot2::annotate(
      "text",
      x = tf_start + bar_width / 2,
      y = title_y,
      label = "TFs",
      family = "Helvetica",
      fontface = "bold",
      size = 3.5
    ) +
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

.write_module3_topic_assignment_qc <- function(optimization,
                                               out_file,
                                               title_prefix = NULL,
                                               condition_colors = NULL,
                                               top_n_tfs = 150L,
                                               seed = 20260716L) {
  .assert_pkg("ggplot2")
  .assert_pkg("gridExtra")
  .assert_pkg("scales")
  .assert_pkg("uwot")
  qc <- optimization$qc
  assignments <- qc$assignments
  sample_rows <- optimization$sample_rows
  sample_assignments <- assignments[sample_rows]
  raw_coordinates <- optimization$raw_sample_coordinates %||%
    .m3_qc_umap_coordinates(
      optimization$raw_sample_probability,
      seed = seed
    )
  identity_map <- identical(
    as.integer(optimization$raw_to_optimized),
    as.integer(names(optimization$raw_to_optimized))
  )
  optimized_coordinates <- if (identity_map) {
    raw_coordinates
  } else if (!is.null(optimization$optimized_sample_coordinates)) {
    optimization$optimized_sample_coordinates
  } else {
    .m3_qc_umap_coordinates(
      optimization$optimized_sample_probability,
      seed = seed
    )
  }
  raw_sample <- data.table::data.table(
    UMAP1 = raw_coordinates[, 1L],
    UMAP2 = raw_coordinates[, 2L],
    condition_id = sample_assignments$condition_id,
    tf = qc$tf_levels[sample_assignments$tf_index],
    target_gene = qc$target_levels[sample_assignments$target_index],
    topic = data.table::fifelse(
      sample_assignments$raw_aligned,
      paste0("Topic ", sample_assignments$raw_target_topic),
      NA_character_
    )
  )
  raw_sample[, topic_short := sub("^Topic ", "T", topic)]
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
  raw_sample <- .m3_qc_jitter_umap(raw_sample, seed = seed + 20L)
  optimized_sample <- .m3_qc_jitter_umap(
    optimized_sample,
    seed = seed + 21L
  )
  condition_values <- sort(unique(assignments$condition_id))
  if (is.null(condition_colors) || !length(condition_colors)) {
    condition_colors <- .module3_bright_topic_palette(condition_values)
  } else {
    missing <- setdiff(condition_values, names(condition_colors))
    condition_colors <- c(
      condition_colors[intersect(names(condition_colors), condition_values)],
      .module3_bright_topic_palette(missing)
    )
  }
  title_prefix <- title_prefix %||% basename(dirname(out_file))
  raw_aligned_sample <- raw_sample[!is.na(topic)]
  optimized_aligned_sample <- optimized_sample[!is.na(topic)]

  raw_aligned_sample[, condition_label := .m3_qc_short_condition_labels(
    condition_id
  )]
  condition_labels <- .m3_qc_short_condition_labels(condition_values)
  condition_palette <- stats::setNames(
    condition_colors[condition_values],
    condition_labels
  )
  topic_values <- sort(unique(c(
    as.integer(qc$raw_topic_ids),
    as.integer(qc$optimized_topic_ids)
  )))
  topic_labels <- paste0("Topic ", topic_values)
  topic_palette <- stats::setNames(
    .module3_bright_topic_palette(topic_labels),
    topic_labels
  )
  page1 <- .m3_qc_arrange(
    .m3_qc_umap_plot(
      raw_aligned_sample,
      "topic",
      "Topic assignment",
      paste0(
        scales::comma(nrow(raw_aligned_sample)),
        " sampled links with posterior primary topic = assigned target topic"
      ),
      colors = topic_palette,
      label_column = "topic_short",
      seed = seed
    ),
    .m3_qc_umap_plot(
      raw_aligned_sample,
      "condition_label",
      "Condition",
      "The same aligned links and UMAP coordinates",
      colors = condition_palette,
      label_column = "condition_label",
      seed = seed + 1L
    ),
    ncol = 1L,
    heights = if (length(condition_values) > 18L) c(0.72, 1.28) else c(1, 1),
    title = paste0(title_prefix, ": filtered aligned-link UMAP")
  )

  raw_topic_similarity <- qc$raw_topic_similarity
  rownames(raw_topic_similarity) <- colnames(raw_topic_similarity) <-
    paste0("Topic ", qc$raw_topic_ids)
  raw_topic_order <- .m3_qc_cluster_order(raw_topic_similarity, "row")
  raw_topic_structure_plot <- .m3_qc_topic_structure_plot(
    similarity = raw_topic_similarity,
    counts = qc$raw_counts,
    topic_column = "raw_topic",
    topic_order = raw_topic_order,
    title = "Raw topic assignment structure",
    subtitle = paste(
      "Mean Hellinger similarity across separately normalized Gene and Peak phi;",
      "bars show full-universe assigned counts"
    )
  )

  optimized_topic_similarity <- qc$optimized_topic_similarity
  rownames(optimized_topic_similarity) <- colnames(optimized_topic_similarity) <-
    paste0("Topic ", qc$optimized_topic_ids)
  optimized_topic_order <- .m3_qc_cluster_order(
    optimized_topic_similarity,
    "row"
  )
  page3 <- .m3_qc_arrange(
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
      )
    ),
    ncol = 1L,
    heights = c(0.85, 1.15),
    title = "Optimized topic assignments"
  )

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
  condition_correlation <- .m3_qc_condition_correlation(
    link_matrix,
    gene_matrix
  )
  condition_correlation_order <- .m3_qc_cluster_order(
    condition_correlation,
    "row"
  )
  condition_correlation_plot <- .m3_qc_correlation_plot(
    condition_correlation,
    "Condition correlation by topic composition",
    paste(
      "Pearson correlation of separately normalized TF-target-link",
      "and target-gene topic profiles",
      sep = "\n"
    ),
    order = condition_correlation_order
  )
  structure_pages <- if (
    max(length(condition_values), length(qc$raw_topic_ids)) > 24L
  ) {
    list(
      .m3_qc_arrange(
        raw_topic_structure_plot,
        ncol = 1L,
        title = "Raw topic assignment structure"
      ),
      .m3_qc_arrange(
        condition_correlation_plot,
        ncol = 1L,
        title = "Condition correlation by topic composition"
      )
    )
  } else {
    list(.m3_qc_arrange(
      raw_topic_structure_plot,
      condition_correlation_plot,
      ncol = 1L,
      heights = c(1.1, 0.9),
      title = "Topic and condition assignment structure"
    ))
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
  condition_theta_ids <- condition_theta
  expression_matrix <- .m3_qc_document_expression_topic_share_matrix(
    optimization = optimization,
    conditions = condition_values,
    topics = qc$optimized_topic_ids
  )
  rownames(condition_theta) <- display_conditions
  rownames(expression_matrix) <- display_conditions
  page_link_matrix <- link_matrix
  page_tf_matrix <- tf_matrix
  page_condition_theta <- condition_theta
  page_expression_matrix <- expression_matrix
  page_column_order <- column_order
  expression_values <- as.numeric(page_expression_matrix)
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
  page4_plots <- .m3_qc_align_plot_dimensions(list(
    .m3_qc_count_heatmap(
      page_link_matrix,
      "Unique aligned TF-target links",
      row_order = row_order,
      column_order = page_column_order,
      square_cells = FALSE,
      show_labels = FALSE
    ) + ggplot2::labs(x = "Topic"),
    .m3_qc_count_heatmap(
      page_tf_matrix,
      "Unique TFs",
      row_order = row_order,
      column_order = page_column_order,
      square_cells = FALSE,
      show_labels = FALSE
    ) + ggplot2::labs(x = "Topic"),
    .m3_qc_value_heatmap(
      page_condition_theta,
      "Mean condition topic probability",
      "Mean theta",
      row_order = row_order,
      column_order = page_column_order
    ) + ggplot2::labs(x = "Topic"),
    .m3_qc_value_heatmap(
      page_expression_matrix,
      "Mean expression of assigned target genes (topic share)",
      "Mean expression share",
      row_order = row_order,
      column_order = page_column_order,
      limits = expression_limits,
      color_transform = "identity"
    ) + ggplot2::labs(x = "Topic")
  ))
  profile_chunks <- if (length(condition_values) > 22L) {
    as.list(seq_along(page4_plots))
  } else if (length(condition_values) > 12L) {
    list(1:2, 3:4)
  } else {
    list(seq_along(page4_plots))
  }
  profile_pages <- lapply(seq_along(profile_chunks), function(page_index) {
    indices <- profile_chunks[[page_index]]
    page_title <- if (length(profile_chunks) > 1L) {
      sprintf(
        "Condition-topic assignment profiles (%d/%d)",
        page_index,
        length(profile_chunks)
      )
    } else {
      "Condition-topic assignment profiles"
    }
    do.call(.m3_qc_arrange, c(
      page4_plots[indices],
      list(ncol = 1L, title = page_title)
    ))
  })

  cross <- qc$condition_topic_similarity$mean
  condition_top_topics <- .m3_qc_top_condition_topics(
    condition_theta_ids,
    cross,
    top_n = 3L
  )
  ordered_conditions <- condition_values[order(
    .m3_qc_short_condition_labels(condition_values),
    condition_values
  )]
  pair_plot_groups <- lapply(seq_along(ordered_conditions), function(i) {
    selected_condition_id <- ordered_conditions[[i]]
    condition_label <- condition_labels[
      match(selected_condition_id, condition_values)
    ]
    selected_topics <- condition_top_topics[
      condition_id == selected_condition_id
    ]
    selected_labels <- selected_topics$topic_label
    selected_short <- sub("^Topic ", "T", selected_labels)
    selected_annotations <- stats::setNames(selected_short, selected_labels)
    condition_highlight <- optimized_sample[
      condition_id == selected_condition_id
    ]
    condition_highlight[, selected := condition_label]
    topic_highlight <- optimized_sample[topic %in% selected_labels]
    topic_highlight[, "topic_annotation" := selected_annotations[topic]]
    list(
      .m3_qc_umap_plot(
        condition_highlight,
        "selected",
        paste0(condition_label, "\n(condition)"),
        sprintf(
          "Top mean-theta topic %s",
          if (length(selected_short)) selected_short[[1L]] else "not available"
        ),
        colors = stats::setNames(
          condition_colors[[selected_condition_id]],
          condition_label
        ),
        background = optimized_sample,
        compact = TRUE
      ),
      .m3_qc_umap_plot(
        topic_highlight,
        "topic",
        "Top theta topics\n(topic)",
        "Top mean-theta topics",
        colors = topic_palette[selected_labels],
        background = optimized_sample,
        label_column = "topic_annotation",
        label_style = "text",
        seed = seed + 30L + i,
        compact = TRUE
      )
    )
  })
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
  pairs_per_page <- if (length(ordered_conditions) > 18L) 4L else 6L
  pair_chunks <- split(
    seq_along(pair_plot_groups),
    ceiling(seq_along(pair_plot_groups) / pairs_per_page)
  )
  pair_pages <- lapply(seq_along(pair_chunks), function(page_index) {
    indices <- pair_chunks[[page_index]]
    plots <- unlist(pair_plot_groups[indices], recursive = FALSE)
    do.call(.m3_qc_arrange, c(
      plots,
      list(
        ncol = 4L,
        title = sprintf(
          "Condition-topic UMAP pairs (%d/%d)",
          page_index,
          length(pair_chunks)
        )
      )
    ))
  })

  pair_candidates <- .m3_opt_parse_topics(
    optimization$raw_pair_assignment$gene_gammafit_topics,
    qc$raw_topic_ids
  )
  peak_candidates <- .m3_opt_parse_topics(
    optimization$raw_pair_assignment$peak_gammafit_topics,
    qc$raw_topic_ids
  )
  target_gamma <- lengths(pair_candidates) > 0L & lengths(peak_candidates) > 0L
  target_max <- is.finite(
    suppressWarnings(as.integer(optimization$raw_pair_assignment$assigned_topic))
  )
  link_max <- target_max[assignments$target_index]
  retention_labels <- .m3_qc_retention_labels()
  retention_link_labels <- retention_labels$links
  retention_gene_labels <- retention_labels$genes
  retention_links <- c(
    .m3_qc_unique_link_count(assignments),
    .m3_qc_unique_link_count(
      assignments[target_gamma[target_index] %in% TRUE]
    ),
    .m3_qc_unique_link_count(assignments[link_max %in% TRUE]),
    .m3_qc_unique_link_count(assignments[raw_aligned == TRUE])
  )
  retention_genes <- c(
    .m3_qc_unique_target_count(
      optimization$raw_pair_assignment$target_gene
    ),
    .m3_qc_unique_target_count(
      optimization$raw_pair_assignment[target_gamma %in% TRUE, target_gene]
    ),
    .m3_qc_unique_target_count(
      optimization$raw_pair_assignment[target_max %in% TRUE, target_gene]
    ),
    .m3_qc_unique_target_count(
      assignments[raw_aligned == TRUE, target_index]
    )
  )
  retention_page <- .m3_qc_arrange(
    .m3_qc_retention_plot(
      retention_link_labels,
      retention_links,
      "TF-target links",
      fill = "#007C78"
    ),
    .m3_qc_retention_plot(
      retention_gene_labels,
      retention_genes,
      "Target genes",
      fill = "#D97824"
    ),
    ncol = 1L,
    title = "Topic assignment retention"
  )

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
  grid::grid.newpage()
  grid::grid.draw(page1)
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
  for (page in pair_pages) {
    grid::grid.newpage()
    grid::grid.draw(page)
  }
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
    raw_pair_assigned_genes = sum(.as_logical_flag(
      optimization$raw_pair_assignment$assigned
    )),
    optimized_pair_assigned_genes = sum(.as_logical_flag(
      optimization$pair_assignment$assigned
    )),
    raw_assigned_genes = data.table::uniqueN(
      optimization$qc$assignments[raw_aligned == TRUE, target_index]
    ),
    optimized_assigned_genes = data.table::uniqueN(
      optimization$qc$assignments[optimized_aligned == TRUE, target_index]
    ),
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
