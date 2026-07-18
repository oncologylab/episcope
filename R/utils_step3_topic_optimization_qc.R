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
      source_pos <- small[order(links[small], genes[small], groups[small])][[1L]]
      target_pos <- which.max(similarity[source_pos, ])
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
  raw_counts <- raw[, .(
    links = .N,
    genes = data.table::uniqueN(target_index)
  ), by = .(raw_topic = raw_target_topic)]
  optimized_counts <- optimized[, .(
    links = .N,
    genes = data.table::uniqueN(target_index)
  ), by = optimized_topic]
  condition_topic <- optimized[, .(
    links = .N,
    genes = data.table::uniqueN(target_index)
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
                                               min_genes = 50L,
                                               min_links = 200L,
                                               similarity_threshold = 0.90,
                                               tf_topic_cutoff = 0.3,
                                               umap_max_links_per_condition = 10000L,
                                               seed = 20260716L,
                                               chunk_size = 50000L) {
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
    raw_to_optimized = stats::setNames(merge$mapping, raw_topic_ids),
    merge_audit = merge$audit,
    sample_rows = sample_rows,
    raw_sample_probability = raw_sample_posterior$sample_probability,
    optimized_sample_probability = optimized_sample_posterior$sample_probability,
    tf_topic_cutoff = tf_topic_cutoff,
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
  out <- sub("^[^_]+_", "", original)
  out <- sub("_Ctrl$", "", out)
  out <- sub("_TGFb$", " TGFb", out)
  out <- gsub("_", " ", out, fixed = TRUE)
  mapping <- unique(data.table::data.table(original, short = out))
  duplicated_label <- duplicated(mapping$short) |
    duplicated(mapping$short, fromLast = TRUE)
  mapping[duplicated_label, short := original]
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

.m3_qc_similarity_plot <- function(x,
                                   title,
                                   subtitle = NULL,
                                   row_order = NULL,
                                   column_order = NULL) {
  if (is.null(row_order)) row_order <- .m3_qc_cluster_order(x, "row")
  if (is.null(column_order)) {
    column_order <- .m3_qc_cluster_order(x, "column")
  }
  long <- .m3_qc_matrix_long(x, row_order, column_order)
  ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(colour = "#D6DCE1", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Similarity"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    .m3_qc_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

.m3_qc_count_heatmap <- function(x,
                                 title,
                                 row_order = NULL,
                                 column_order = NULL,
                                 label_min = 1L) {
  if (is.null(row_order)) row_order <- .m3_qc_cluster_order(x, "row")
  if (is.null(column_order)) column_order <- .m3_qc_cluster_order(x, "column")
  long <- .m3_qc_matrix_long(x, row_order, column_order)
  fill_max <- max(long$value, na.rm = TRUE)
  if (!is.finite(fill_max) || fill_max <= 0) fill_max <- 1
  long[, label := ifelse(value >= label_min, .m3_qc_compact_count(value), "")]
  ggplot2::ggplot(
    long,
    ggplot2::aes(column_label, row_label, fill = value)
  ) +
    ggplot2::geom_tile(colour = "#D6DCE1", linewidth = 0.2) +
    ggplot2::geom_text(
      ggplot2::aes(label = label, colour = value >= fill_max * 0.55),
      family = "Helvetica",
      fontface = "bold",
      size = 3.2,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#20272E")) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      trans = "sqrt",
      limits = c(0, fill_max),
      labels = scales::label_comma(),
      name = "Count"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::coord_fixed(ratio = 1) +
    .m3_qc_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
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
  ordered_labels <- paste0("Topic ", ordered_ids)
  k <- length(ordered_ids)
  reordered <- similarity[topic_order, topic_order, drop = FALSE]
  heatmap <- data.table::as.data.table(as.table(reordered))
  data.table::setnames(heatmap, c("row_label", "column_label", "similarity"))
  heatmap[, `:=`(
    x = match(as.character(column_label), ordered_labels),
    y = k - match(as.character(row_label), ordered_labels) + 1L,
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
  link_max <- max(count_data$links, na.rm = TRUE)
  gene_max <- max(count_data$genes, na.rm = TRUE)
  if (!is.finite(link_max) || link_max <= 0) link_max <- 1
  if (!is.finite(gene_max) || gene_max <= 0) gene_max <- 1
  count_data[, `:=`(
    link_xmax = link_start + links / link_max * bar_width,
    gene_xmax = gene_start + genes / gene_max * bar_width,
    link_label = .m3_qc_compact_count(links),
    gene_label = .m3_qc_compact_count(genes),
    link_inside = links / link_max >= 0.28,
    gene_inside = genes / gene_max >= 0.28
  )]
  right_edge <- gene_start + bar_width + 1.1
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
                             seed = 20260716L) {
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
  p +
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

.m3_qc_funnel_plot <- function(labels,
                               counts,
                               title,
                               fill = "#007C78") {
  x <- data.table::data.table(
    stage = factor(labels, levels = rev(labels)),
    count = as.numeric(counts)
  )
  ggplot2::ggplot(x, ggplot2::aes(count, stage)) +
    ggplot2::geom_col(
      width = 0.72,
      fill = fill,
      colour = "#30383F",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(count)),
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
    .m3_qc_theme()
}

.m3_qc_tf_topic_pages <- function(optimization,
                                  top_n_tfs = 100L,
                                  rows_per_page = 34L) {
  x <- optimization$qc$assignments[optimized_aligned == TRUE]
  if (!nrow(x)) return(list())
  x <- x[, .(
    targets = data.table::uniqueN(target_index)
  ), by = .(condition_id, tf_index, optimized_topic)]
  totals <- x[, .(targets = sum(targets)), by = .(condition_id, tf_index)]
  data.table::setorder(totals, condition_id, -targets, tf_index)
  totals <- totals[, head(.SD, as.integer(top_n_tfs)), by = condition_id]
  x <- x[totals, on = .(condition_id, tf_index), nomatch = 0L]
  topics <- optimization$qc$optimized_topic_ids
  pages <- list()
  page_index <- 0L
  for (condition in unique(totals$condition_id)) {
    selected <- totals[condition_id == condition, tf_index]
    selected <- selected[!duplicated(selected)]
    mat <- matrix(
      0,
      nrow = length(selected),
      ncol = length(topics),
      dimnames = list(
        optimization$qc$tf_levels[selected],
        paste0("Topic ", topics)
      )
    )
    values <- x[condition_id == condition]
    mat[cbind(
      match(values$tf_index, selected),
      match(values$optimized_topic, topics)
    )] <- values$targets
    row_order <- .m3_qc_cluster_order(mat, "row")
    column_order <- .m3_qc_cluster_order(mat, "column")
    chunks <- split(row_order, ceiling(seq_along(row_order) / rows_per_page))
    for (i in seq_along(chunks)) {
      page_index <- page_index + 1L
      page_matrix <- mat[
        chunks[[i]],
        column_order,
        drop = FALSE
      ]
      pages[[page_index]] <- .m3_qc_count_heatmap(
        page_matrix,
        title = sprintf(
          "%s: top TFs by unique targets (%d/%d)",
          .m3_qc_short_condition_labels(condition),
          i,
          length(chunks)
        ),
        row_order = seq_len(nrow(page_matrix)),
        column_order = seq_len(ncol(page_matrix)),
        label_min = 1L
      ) +
        ggplot2::labs(
          subtitle = "Cells show unique target genes for each TF and optimized topic",
          fill = "Targets"
        )
    }
  }
  pages
}

.write_module3_topic_assignment_qc <- function(optimization,
                                               out_file,
                                               title_prefix = NULL,
                                               condition_colors = NULL,
                                               top_n_tfs = 100L,
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
  condition_palette <- stats::setNames(
    condition_colors[raw_aligned_sample[, unique(condition_id)]],
    .m3_qc_short_condition_labels(raw_aligned_sample[, unique(condition_id)])
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
    title = paste0(title_prefix, ": filtered aligned-link UMAP")
  )

  raw_topic_similarity <- qc$raw_topic_similarity
  rownames(raw_topic_similarity) <- colnames(raw_topic_similarity) <-
    paste0("Topic ", qc$raw_topic_ids)
  raw_topic_order <- .m3_qc_cluster_order(raw_topic_similarity, "row")
  page2 <- .m3_qc_arrange(
    .m3_qc_topic_structure_plot(
      similarity = raw_topic_similarity,
      counts = qc$raw_counts,
      topic_column = "raw_topic",
      topic_order = raw_topic_order,
      title = "Raw topic assignment structure",
      subtitle = paste(
        "Mean Hellinger similarity across separately normalized Gene and Peak phi;",
        "bars show full-universe assigned counts"
      )
    ),
    ncol = 1L
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
      "Optimized topics",
      "Sampled links retained after topic merge and alignment",
      label_column = "topic_short",
      seed = seed + 4L
    ),
    ncol = 1L,
    title = "Optimized filtered aligned-link UMAP"
  )
  page3_structure <- .m3_qc_arrange(
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
    ncol = 1L
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
  display_conditions <- .m3_qc_short_condition_labels(rownames(link_matrix))
  rownames(link_matrix) <- display_conditions
  rownames(gene_matrix) <- display_conditions
  clustering_matrix <- cbind(
    .m3_opt_row_normalize(link_matrix),
    .m3_opt_row_normalize(gene_matrix)
  )
  row_order <- .m3_qc_cluster_order(clustering_matrix, "row")
  column_order <- .m3_qc_cluster_order(link_matrix + gene_matrix, "column")
  page4 <- .m3_qc_arrange(
    .m3_qc_count_heatmap(
      link_matrix,
      "Aligned TF-target links",
      row_order = row_order,
      column_order = column_order
    ),
    .m3_qc_count_heatmap(
      gene_matrix,
      "Assigned target genes",
      row_order = row_order,
      column_order = column_order
    ),
    ncol = 1L,
    title = "Condition-topic assignment counts"
  )

  cross <- qc$condition_topic_similarity$mean
  pair_table <- data.table::as.data.table(as.table(cross))
  data.table::setnames(pair_table, c("condition_id", "topic", "similarity"))
  pair_table[, `:=`(
    condition_id = as.character(condition_id),
    topic = as.character(topic),
    topic_num = suppressWarnings(as.integer(as.character(topic)))
  )]
  data.table::setorder(pair_table, -similarity, condition_id, topic_num)
  chosen <- pair_table[0]
  used_conditions <- character()
  used_topics <- integer()
  for (i in seq_len(nrow(pair_table))) {
    candidate <- pair_table[i]
    if (candidate$condition_id %in% used_conditions ||
        candidate$topic_num %in% used_topics) next
    chosen <- data.table::rbindlist(list(chosen, candidate))
    used_conditions <- c(used_conditions, candidate$condition_id)
    used_topics <- c(used_topics, candidate$topic_num)
    if (nrow(chosen) >= 4L) break
  }
  pair_plots <- lapply(seq_len(nrow(chosen)), function(i) {
    candidate <- chosen[i]
    highlighted <- optimized_sample[
      condition_id == candidate$condition_id &
        topic == paste0("Topic ", candidate$topic_num)
    ]
    highlighted[, selected := candidate$condition_id]
    .m3_qc_umap_plot(
      highlighted,
      "selected",
      sprintf(
        "%s / Topic %d",
        .m3_qc_short_condition_labels(candidate$condition_id),
        candidate$topic_num
      ),
      sprintf("Mean Jaccard %.3f", candidate$similarity),
      colors = stats::setNames(
        condition_colors[[candidate$condition_id]],
        candidate$condition_id
      ),
      background = optimized_sample
    )
  })
  cross_display <- cross
  rownames(cross_display) <- .m3_qc_short_condition_labels(rownames(cross_display))
  cross_plot <- .m3_qc_similarity_plot(
    cross_display,
    "Condition-topic similarity",
    "Mean Jaccard of unique TF-target pairs and target genes"
  )
  page5_args <- c(
    list(cross_plot),
    pair_plots,
    list(title = "Top distinct condition-topic matches")
  )
  if (length(pair_plots) == 4L) {
    page5_args$layout_matrix <- rbind(
      c(1L, 1L),
      c(2L, 3L),
      c(4L, 5L)
    )
    page5_args$heights <- c(1.15, 1, 1)
  } else {
    page5_args$ncol <- 2L
  }
  page5 <- do.call(.m3_qc_arrange, page5_args)

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
  link_gamma <- target_gamma[assignments$target_index]
  link_max <- target_max[assignments$target_index]
  funnel_link_labels <- c(
    "Input universe",
    "Gene and Peak pass GammaFit",
    "Gene and Peak MaxProb agree",
    "Link posterior agrees",
    sprintf("TF theta >= %.2f", optimization$tf_topic_cutoff),
    "Retained/recovered after optimization"
  )
  funnel_gene_labels <- c(
    "Input target genes",
    "Gene and Peak pass GammaFit",
    "Gene and Peak MaxProb agree",
    "Genes with posterior-agreeing link",
    sprintf(
      "Genes with TF theta >= %.2f link",
      optimization$tf_topic_cutoff
    ),
    "Assigned/recovered after optimization"
  )
  funnel_links <- c(
    nrow(assignments),
    sum(link_gamma),
    sum(link_max),
    sum(assignments$raw_posterior_agrees),
    sum(assignments$raw_aligned),
    sum(assignments$optimized_aligned)
  )
  funnel_genes <- c(
    nrow(optimization$raw_pair_assignment),
    sum(target_gamma),
    sum(target_max),
    data.table::uniqueN(
      assignments[raw_posterior_agrees == TRUE, target_index]
    ),
    data.table::uniqueN(assignments[raw_aligned == TRUE, target_index]),
    data.table::uniqueN(assignments[optimized_aligned == TRUE, target_index])
  )
  funnel_page <- .m3_qc_arrange(
    .m3_qc_funnel_plot(
      funnel_gene_labels,
      funnel_genes,
      "Target-gene assignment funnel",
      fill = "#D97824"
    ),
    .m3_qc_funnel_plot(
      funnel_link_labels,
      funnel_links,
      "TF-target link assignment funnel",
      fill = "#007C78"
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
  grid::grid.newpage()
  grid::grid.draw(page2)
  if (!identity_map) {
    grid::grid.newpage()
    grid::grid.draw(page3)
    grid::grid.newpage()
    grid::grid.draw(page3_structure)
  }
  grid::grid.newpage()
  grid::grid.draw(page4)
  grid::grid.newpage()
  grid::grid.draw(page5)
  for (page in tf_pages) {
    grid::grid.newpage()
    print(page, newpage = FALSE)
  }
  grid::grid.newpage()
  grid::grid.draw(funnel_page)
  invisible(out_file)
}

.write_module3_topic_optimization_outputs <- function(optimization,
                                                      out_dir,
                                                      raw_k,
                                                      title_prefix = NULL,
                                                      condition_colors = NULL,
                                                      top_n_tfs = 100L,
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
    raw_topics = length(unique(mapping$raw_topic)),
    optimized_topics = length(unique(mapping$optimized_topic)),
    merged_topics = length(unique(mapping$raw_topic)) -
      length(unique(mapping$optimized_topic)),
    input_links = nrow(optimization$qc$assignments),
    raw_aligned_links = sum(optimization$qc$assignments$raw_aligned),
    optimized_aligned_links = sum(optimization$qc$assignments$optimized_aligned),
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
