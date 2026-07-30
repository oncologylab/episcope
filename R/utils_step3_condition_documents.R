# Condition-document topic inputs and post-model TF-topic evidence.

.module3_validate_condition_document_tokens <- function(gene_fraction,
                                                        token_budget) {
  gene_fraction <- suppressWarnings(as.numeric(gene_fraction[[1L]]))
  token_budget <- suppressWarnings(as.integer(token_budget[[1L]]))
  if (!is.finite(gene_fraction) || gene_fraction <= 0 || gene_fraction >= 1) {
    .log_abort("`condition_gene_token_fraction` must be strictly between 0 and 1.")
  }
  if (!is.finite(token_budget) || token_budget < 2L) {
    .log_abort("`condition_document_token_budget` must be an integer >= 2.")
  }
  list(gene_fraction = gene_fraction, token_budget = token_budget)
}

.module3_largest_remainder_counts <- function(score,
                                               total,
                                               tie_key = seq_along(score),
                                               minimum = 1L) {
  score <- as.numeric(score)
  score[!is.finite(score) | score < 0] <- 0
  total <- as.integer(total[[1L]])
  minimum <- as.integer(minimum[[1L]])
  n <- length(score)
  if (!n) return(integer())
  if (!is.finite(total) || total < n * minimum) {
    .log_abort(
      "Token allocation of {total} cannot retain {n} terms with minimum {minimum}."
    )
  }
  counts <- rep.int(minimum, n)
  remaining <- total - sum(counts)
  if (!remaining) return(counts)
  if (!any(score > 0)) score[] <- 1
  quota <- score / sum(score) * remaining
  extra <- floor(quota)
  counts <- counts + as.integer(extra)
  left <- total - sum(counts)
  if (left > 0L) {
    order_index <- order(
      -(quota - extra),
      as.character(tie_key),
      na.last = TRUE
    )
    counts[order_index[seq_len(left)]] <-
      counts[order_index[seq_len(left)]] + 1L
  }
  counts
}

.module3_allocate_condition_document_tokens <- function(values,
                                                        gene_fraction = 0.6,
                                                        token_budget = 1000000L) {
  .assert_pkg("data.table")
  x <- data.table::as.data.table(values)
  .assert_has_cols(
    x,
    c("condition_id", "term_id", "modality", "token_score"),
    context = "condition-document token allocation"
  )
  settings <- .module3_validate_condition_document_tokens(
    gene_fraction,
    token_budget
  )
  x[, modality := tolower(as.character(modality))]
  if (any(!x$modality %in% c("gene", "peak"))) {
    .log_abort("Condition-document modalities must be `gene` or `peak`.")
  }
  if (x[, anyDuplicated(paste(condition_id, term_id, sep = "\r"))] > 0L) {
    .log_abort("Condition-document values contain duplicate condition-term rows.")
  }
  conditions <- sort(unique(as.character(x$condition_id)))
  out <- vector("list", length(conditions))
  target_gene <- as.integer(round(
    settings$token_budget * settings$gene_fraction
  ))
  target_peak <- settings$token_budget - target_gene
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]]
    part <- data.table::copy(x[condition_id == condition])
    gene_index <- which(part$modality == "gene")
    peak_index <- which(part$modality == "peak")
    if (!length(gene_index) || !length(peak_index)) {
      .log_abort(
        "Condition {condition} must contain both Gene and Peak terms."
      )
    }
    if (length(gene_index) > target_gene || length(peak_index) > target_peak) {
      .log_abort(
        paste0(
          "Condition ", condition, " has more retained terms than its token ",
          "allocation. Increase `condition_document_token_budget` or adjust ",
          "`condition_gene_token_fraction`."
        )
      )
    }
    part[, pseudo_count_log := 0L]
    part$pseudo_count_log[gene_index] <-
      .module3_largest_remainder_counts(
        part$token_score[gene_index],
        target_gene,
        part$term_id[gene_index]
      )
    part$pseudo_count_log[peak_index] <-
      .module3_largest_remainder_counts(
        part$token_score[peak_index],
        target_peak,
        part$term_id[peak_index]
      )
    out[[i]] <- part
  }
  result <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  audit <- result[, .(
    n_terms = data.table::uniqueN(term_id),
    model_tokens = sum(as.double(pseudo_count_log))
  ), by = .(condition_id, modality)]
  audit[, expected_tokens := data.table::fifelse(
    modality == "gene",
    target_gene,
    target_peak
  )]
  if (audit[, any(model_tokens != expected_tokens)]) {
    .log_abort("Exact condition-document token allocation failed.")
  }
  totals <- audit[, .(model_tokens = sum(model_tokens)), by = condition_id]
  if (totals[, any(model_tokens != settings$token_budget)]) {
    .log_abort("One or more condition documents did not receive the exact token budget.")
  }
  data.table::setattr(result, "token_audit", audit)
  result[]
}

.module3_condition_specificity_values <- function(values,
                                                  temperature = 0.5,
                                                  uniform_floor = 0.1) {
  x <- data.table::as.data.table(values)
  temperature <- suppressWarnings(as.numeric(temperature[[1L]]))
  uniform_floor <- suppressWarnings(as.numeric(uniform_floor[[1L]]))
  if (!is.finite(temperature) || temperature <= 0) {
    .log_abort("`condition_specificity_temperature` must be positive.")
  }
  if (!is.finite(uniform_floor) || uniform_floor < 0 || uniform_floor > 1) {
    .log_abort("`condition_specificity_floor` must be between 0 and 1.")
  }
  x[, scaled__ := log2(value + 1) / temperature]
  x[, centered__ := scaled__ - max(scaled__, na.rm = TRUE), by = feature_id]
  x[, mass__ := exp(centered__)]
  x[!is.finite(mass__), mass__ := 0]
  x[, `:=`(
    mass_total__ = sum(mass__),
    observed_conditions__ = .N
  ), by = feature_id]
  x[, multiplier := uniform_floor +
      (1 - uniform_floor) *
      mass__ / pmax(mass_total__, .Machine$double.eps) *
      observed_conditions__]
  x[, original_mass__ := sum(value), by = .(condition_id, modality)]
  x[, weighted__ := value * multiplier]
  x[, weighted_mass__ := sum(weighted__), by = .(condition_id, modality)]
  x[, weight := weighted__ * original_mass__ /
      pmax(weighted_mass__, .Machine$double.eps)]
  x[, token_score := log1p(pmax(weight, 0))]
  x[, c(
    "scaled__", "centered__", "mass__", "mass_total__",
    "observed_conditions__", "original_mass__", "weighted__",
    "weighted_mass__"
  ) := NULL]
  x[]
}

.module3_collect_condition_document_features <- function(path,
                                                         format,
                                                         condition_id,
                                                         threshold_gene_expr,
                                                         threshold_fp_score) {
  format <- tolower(as.character(format[[1L]]))
  path <- as.character(path[[1L]])
  if (identical(format, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
    dataset <- arrow::open_dataset(path)
    schema_names <- names(dataset$schema)
    needed <- c(
      "gene_key", "gene_expr_condition",
      "peak_id", "fp_score_condition"
    )
    missing <- setdiff(needed, schema_names)
    if (length(missing)) {
      .log_abort(
        "Condition-link file is missing required columns: {paste(missing, collapse = ', ')}"
      )
    }
    gene <- dplyr::collect(
      dplyr::summarise(
        dplyr::group_by(
          dplyr::select(dataset, gene_key, gene_expr_condition),
          gene_key
        ),
        value = max(gene_expr_condition, na.rm = TRUE),
        .groups = "drop"
      )
    )
    peak <- dplyr::collect(
      dplyr::summarise(
        dplyr::group_by(
          dplyr::select(dataset, peak_id, fp_score_condition),
          peak_id
        ),
        value = max(fp_score_condition, na.rm = TRUE),
        .groups = "drop"
      )
    )
  } else {
    links <- .module3_read_table(
      path,
      format,
      columns = c(
        "gene_key", "gene_expr_condition",
        "peak_id", "fp_score_condition"
      ),
      allow_missing_columns = FALSE
    )
    gene <- links[, .(
      value = max(gene_expr_condition, na.rm = TRUE)
    ), by = gene_key]
    peak <- links[, .(
      value = max(fp_score_condition, na.rm = TRUE)
    ), by = peak_id]
  }
  gene <- data.table::as.data.table(gene)
  peak <- data.table::as.data.table(peak)
  data.table::setnames(gene, "gene_key", "feature_id")
  data.table::setnames(peak, "peak_id", "feature_id")
  gene[, `:=`(condition_id = condition_id, modality = "gene")]
  peak[, `:=`(condition_id = condition_id, modality = "peak")]
  gene <- gene[
    !is.na(feature_id) & nzchar(feature_id) &
      is.finite(value) & value >= threshold_gene_expr
  ]
  peak <- peak[
    !is.na(feature_id) &
      grepl("^chr[^:]+:[0-9]+-[0-9]+$", feature_id) &
      is.finite(value) & value > threshold_fp_score
  ]
  data.table::rbindlist(
    list(
      gene[, .(condition_id, feature_id, modality, value)],
      peak[, .(condition_id, feature_id, modality, value)]
    ),
    use.names = TRUE
  )
}

.module3_prepare_condition_document_inputs <- function(filtered_dir,
                                                       output_dir,
                                                       sample_subset = NULL,
                                                       analysis_label = NULL,
                                                       threshold_gene_expr = 10,
                                                       threshold_fp_score = 0,
                                                       condition_specificity_temperature = 0.5,
                                                       condition_specificity_floor = 0.1,
                                                       condition_gene_token_fraction = 0.6,
                                                       condition_document_token_budget = 1000000L,
                                                       input_signature = NA_character_,
                                                       save_full_doc_term_csv = FALSE,
                                                       verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("Matrix")
  manifest_path <- .module3_condition_link_manifest_path(filtered_dir)
  manifest <- data.table::fread(manifest_path, showProgress = FALSE)
  if (!all(c("condition_id", "path", "format") %in% names(manifest))) {
    .log_abort(
      "condition_links_manifest.csv must include condition_id, path, and format."
    )
  }
  subset <- if (is.null(sample_subset)) character() else {
    unique(as.character(sample_subset))
  }
  subset <- subset[!is.na(subset) & nzchar(subset)]
  if (length(subset)) manifest <- manifest[condition_id %in% subset]
  if (!nrow(manifest)) {
    .log_abort("No condition-link manifest rows remain for condition documents.")
  }
  manifest[, path := .module3_resolve_condition_manifest_paths(
    path,
    manifest_path
  )]
  manifest <- manifest[order(condition_id)]
  values <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    condition <- as.character(manifest$condition_id[[i]])
    if (isTRUE(verbose)) {
      .log_inform(
        "Preparing condition document {i}/{nrow(manifest)}: {condition}."
      )
    }
    values[[i]] <- .module3_collect_condition_document_features(
      manifest$path[[i]],
      manifest$format[[i]],
      condition,
      threshold_gene_expr,
      threshold_fp_score
    )
  }
  values <- data.table::rbindlist(values, use.names = TRUE, fill = TRUE)
  if (!nrow(values)) .log_abort("No condition-document terms passed input filters.")
  values <- .module3_condition_specificity_values(
    values,
    temperature = condition_specificity_temperature,
    uniform_floor = condition_specificity_floor
  )
  values[, term_id := data.table::fifelse(
    modality == "gene",
    paste0("GENE:", feature_id),
    paste0("PEAK:", feature_id)
  )]
  values <- .module3_allocate_condition_document_tokens(
    values,
    gene_fraction = condition_gene_token_fraction,
    token_budget = condition_document_token_budget
  )
  token_audit <- attr(values, "token_audit")
  values[, `:=`(
    doc_id = condition_id,
    pseudo_count = pseudo_count_log,
    pseudo_count_bin = 1L
  )]
  doc_term <- values[, .(
    condition_id,
    doc_id,
    term_id,
    modality,
    value,
    multiplier,
    weight,
    token_score,
    pseudo_count_log,
    pseudo_count,
    pseudo_count_bin
  )]
  data.table::setorder(doc_term, doc_id, term_id)
  edges_docs <- manifest[, .(
    condition_id = as.character(condition_id),
    doc_id = as.character(condition_id),
    source_path = as.character(path)
  )]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(
      values,
      file.path(output_dir, "condition_term_values.parquet"),
      compression = "zstd"
    )
  }
  write_doc_term_cache(
    doc_term,
    out_dir = output_dir,
    save_full_doc_term_csv = isTRUE(save_full_doc_term_csv)
  )
  .save_all(output_dir, "edges_docs", edges_docs)
  .save_all(output_dir, "doc_term", doc_term)
  dtm_obj <- build_sparse_dtm(doc_term, count_col = "pseudo_count_log")
  .save_all(output_dir, "dtm", dtm_obj$dtm)
  .save_all(
    output_dir,
    "dtm_index",
    list(
      doc_index = dtm_obj$doc_index,
      term_index = dtm_obj$term_index
    )
  )
  data.table::fwrite(
    token_audit,
    file.path(output_dir, "condition_document_token_audit.csv")
  )
  .write_module3_document_term_qc(
    doc_term = doc_term,
    output_dir = output_dir,
    count_column = "pseudo_count_log",
    title = paste0(
      "Module 3 document-term QC - ",
      as.character(analysis_label %||% "condition documents")[[1L]]
    ),
    document_unit = "condition",
    verbose = verbose
  )
  summary <- data.table::data.table(
    analysis_label = as.character(
      analysis_label %||% "condition_documents"
    )[[1L]],
    input_signature = as.character(input_signature[[1L]]),
    input_source = "condition_links",
    doc_design = "condition",
    doc_mode = "condition",
    document_unit = "condition",
    condition_document_unit = "condition",
    fp_term_mode = "unique",
    condition_gene_weighting = "specificity",
    condition_peak_weighting = "specificity",
    condition_specificity_temperature =
      as.numeric(condition_specificity_temperature),
    condition_specificity_floor = as.numeric(condition_specificity_floor),
    condition_gene_token_fraction =
      as.numeric(condition_gene_token_fraction),
    condition_peak_token_fraction =
      1 - as.numeric(condition_gene_token_fraction),
    condition_document_token_budget =
      as.integer(condition_document_token_budget),
    count_method = "log",
    count_scale = NA_real_,
    count_input_requested = "pseudo_count_log",
    count_input_effective = "pseudo_count_log",
    n_link_rows_loaded = sum(
      suppressWarnings(as.numeric(manifest$n_links)),
      na.rm = TRUE
    ),
    n_link_rows_after_subset = sum(
      suppressWarnings(as.numeric(manifest$n_links)),
      na.rm = TRUE
    ),
    n_link_rows_after_filter = sum(
      suppressWarnings(as.numeric(manifest$n_links)),
      na.rm = TRUE
    ),
    n_document_edge_rows = nrow(edges_docs),
    n_doc_term_rows = nrow(doc_term),
    n_documents = data.table::uniqueN(doc_term$doc_id),
    n_gene_terms = data.table::uniqueN(
      doc_term[modality == "gene", term_id]
    ),
    n_peak_terms = data.table::uniqueN(
      doc_term[modality == "peak", term_id]
    ),
    n_terms = data.table::uniqueN(doc_term$term_id),
    n_nonzero = nrow(doc_term),
    n_model_tokens = sum(as.double(doc_term$pseudo_count_log))
  )
  data.table::fwrite(
    summary,
    file.path(output_dir, "topic_input_summary.csv")
  )
  if (isTRUE(verbose)) {
    .log_inform(
      paste0(
        "Prepared ", summary$n_documents, " condition document(s), ",
        summary$n_gene_terms, " Gene terms, and ",
        summary$n_peak_terms, " coordinate Peak terms."
      )
    )
  }
  invisible(list(output_dir = output_dir, summary = summary, reused = FALSE))
}

.module3_row_normalize_nonnegative <- function(x) {
  x <- as.matrix(x)
  x[!is.finite(x) | x < 0] <- 0
  totals <- rowSums(x)
  totals[!is.finite(totals) | totals <= 0] <- 1
  x / totals
}

.module3_hellinger_similarity <- function(x) {
  x <- .module3_row_normalize_nonnegative(x)
  sqrt(x) %*% t(sqrt(x))
}

.module3_align_topic_fits <- function(reference,
                                      candidate,
                                      gene_weight = 0.4,
                                      peak_weight = 0.4,
                                      theta_weight = 0.2) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    .log_abort("Package {.pkg clue} is required for one-to-one topic alignment.")
  }
  ref_phi <- as.matrix(reference$phi)
  cand_phi <- as.matrix(candidate$phi)
  if (!identical(dim(ref_phi), dim(cand_phi)) ||
      !identical(colnames(ref_phi), colnames(cand_phi))) {
    .log_abort("Topic fits must have matching phi dimensions and term order.")
  }
  if (!identical(
    rownames(as.matrix(reference$theta)),
    rownames(as.matrix(candidate$theta))
  )) {
    .log_abort("Topic fits must have matching theta document order.")
  }
  gene_columns <- startsWith(colnames(ref_phi), "GENE:")
  peak_columns <- startsWith(colnames(ref_phi), "PEAK:")
  cross_similarity <- function(a, b) {
    a <- .module3_row_normalize_nonnegative(a)
    b <- .module3_row_normalize_nonnegative(b)
    sqrt(a) %*% t(sqrt(b))
  }
  gene_similarity <- cross_similarity(
    ref_phi[, gene_columns, drop = FALSE],
    cand_phi[, gene_columns, drop = FALSE]
  )
  peak_similarity <- cross_similarity(
    ref_phi[, peak_columns, drop = FALSE],
    cand_phi[, peak_columns, drop = FALSE]
  )
  ref_theta <- t(.module3_row_normalize_nonnegative(reference$theta))
  cand_theta <- t(.module3_row_normalize_nonnegative(candidate$theta))
  theta_similarity <- cross_similarity(ref_theta, cand_theta)
  similarity <- gene_weight * gene_similarity +
    peak_weight * peak_similarity +
    theta_weight * theta_similarity
  assignment <- as.integer(clue::solve_LSAP(similarity, maximum = TRUE))
  data.table::data.table(
    reference_topic = seq_len(nrow(similarity)),
    candidate_topic = assignment,
    similarity = similarity[cbind(seq_len(nrow(similarity)), assignment)],
    gene_similarity =
      gene_similarity[cbind(seq_len(nrow(similarity)), assignment)],
    peak_similarity =
      peak_similarity[cbind(seq_len(nrow(similarity)), assignment)],
    theta_similarity =
      theta_similarity[cbind(seq_len(nrow(similarity)), assignment)]
  )
}

.module3_select_medoid_topic_fit <- function(fits,
                                             seeds = names(fits)) {
  if (!is.list(fits) || length(fits) < 2L) {
    .log_abort("At least two topic fits are required for medoid selection.")
  }
  if (is.null(seeds) || length(seeds) != length(fits)) {
    seeds <- as.character(seq_along(fits))
  }
  pair_rows <- list()
  position <- 0L
  for (i in seq_along(fits)) {
    for (j in seq_along(fits)) {
      if (i >= j) next
      alignment <- .module3_align_topic_fits(fits[[i]], fits[[j]])
      position <- position + 1L
      pair_rows[[position]] <- data.table::data.table(
        seed_1 = as.character(seeds[[i]]),
        seed_2 = as.character(seeds[[j]]),
        median_similarity = stats::median(alignment$similarity),
        mean_similarity = mean(alignment$similarity),
        min_similarity = min(alignment$similarity)
      )
    }
  }
  pairs <- data.table::rbindlist(pair_rows)
  scores <- data.table::rbindlist(lapply(as.character(seeds), function(seed) {
    values <- pairs[
      seed_1 == seed | seed_2 == seed,
      mean_similarity
    ]
    data.table::data.table(
      seed = seed,
      medoid_score = mean(values),
      weakest_pair = min(values)
    )
  }))
  scores[, seed_order := match(seed, as.character(seeds))]
  data.table::setorder(scores, -medoid_score, -weakest_pair, seed_order)
  list(
    seed = scores$seed[[1L]],
    fit = fits[[match(scores$seed[[1L]], as.character(seeds))]],
    scores = scores[, seed_order := NULL],
    pairs = pairs
  )
}

.module3_tf_topic_component_normalize <- function(x,
                                                  value_column,
                                                  background_column,
                                                  output_column) {
  value_column <- as.character(value_column[[1L]])
  background_column <- as.character(background_column[[1L]])
  output_column <- as.character(output_column[[1L]])
  x[, adjusted__ := pmax(get(value_column), 0) /
      sqrt(pmax(get(background_column), .Machine$double.eps))]
  x[, (output_column) := {
    total <- sum(adjusted__, na.rm = TRUE)
    if (is.finite(total) && total > 0) adjusted__ / total else 0
  }, by = .(condition_id, tf)]
  x[, adjusted__ := NULL]
  x
}

.module3_condition_tf_topic_evidence <- function(aligned_links,
                                                 peak_values,
                                                 gene_values,
                                                 prediction_stats,
                                                 topic_count,
                                                 tf_expression_min = 10,
                                                 binding_weight = 0.7,
                                                 target_weight = 0.3,
                                                 membership_cutoff = 0.3,
                                                 primary_margin_cutoff = 0.1) {
  links <- data.table::as.data.table(aligned_links)
  peaks <- data.table::as.data.table(peak_values)
  genes <- data.table::as.data.table(gene_values)
  stats <- data.table::as.data.table(prediction_stats)
  .assert_has_cols(
    links,
    c(
      "condition_id", "tf", "peak_id", "gene_key",
      "topic_num", "tf_expression"
    ),
    context = "aligned condition links"
  )
  .assert_has_cols(
    peaks,
    c("condition_id", "peak_id", "topic_num", "token_share"),
    context = "condition Peak values"
  )
  .assert_has_cols(
    genes,
    c("condition_id", "gene_key", "topic_num", "token_share"),
    context = "condition Gene values"
  )
  .assert_has_cols(
    stats,
    c("tf", "fp_id", "best_r"),
    context = "Module 1 prediction statistics"
  )
  binding_weight <- as.numeric(binding_weight[[1L]])
  target_weight <- as.numeric(target_weight[[1L]])
  if (!is.finite(binding_weight) || !is.finite(target_weight) ||
      binding_weight < 0 || target_weight < 0 ||
      binding_weight + target_weight <= 0) {
    .log_abort("TF-topic component weights must be non-negative with a positive sum.")
  }
  total_weight <- binding_weight + target_weight
  binding_weight <- binding_weight / total_weight
  target_weight <- target_weight / total_weight
  stats <- stats[
    is.finite(best_r) & best_r >= 0.3,
    .(
      best_r = max(best_r, na.rm = TRUE)
    ),
    by = .(tf, fp_id)
  ]
  data.table::setnames(stats, "fp_id", "peak_id")
  links <- stats[
    links,
    on = c("tf", "peak_id")
  ]
  join_audit <- links[, .(
    total_aligned_rows = .N,
    matched_correlation_rows = sum(is.finite(best_r)),
    matched_percent = 100 * mean(is.finite(best_r))
  )]
  links <- links[
    is.finite(best_r) &
      is.finite(tf_expression) &
      tf_expression >= tf_expression_min
  ]
  if (!nrow(links)) {
    .log_abort("No TF-Peak evidence remains after correlation and expression filters.")
  }
  peak_topic <- unique(peaks[, .(
    condition_id,
    peak_id,
    peak_topic = as.integer(topic_num),
    peak_share = as.numeric(token_share)
  )])
  gene_topic <- unique(genes[, .(
    condition_id,
    gene_key,
    gene_topic = as.integer(topic_num),
    gene_share = as.numeric(token_share)
  )])
  links <- peak_topic[
    links,
    on = c("condition_id", "peak_id" = "peak_id"),
    nomatch = 0L
  ]
  links <- gene_topic[
    links,
    on = c("condition_id", "gene_key" = "gene_key"),
    nomatch = 0L
  ]
  links <- links[
    peak_topic == gene_topic &
      peak_topic == as.integer(topic_num)
  ]
  if (!nrow(links)) {
    .log_abort("No TF-Peak-Gene rows agree on a retained topic.")
  }
  links[, reliability := pmax(best_r, 0)^2]
  unique_peak <- unique(links[, .(
    condition_id,
    tf,
    topic_num = peak_topic,
    peak_id,
    peak_share,
    reliability
  )])
  binding <- unique_peak[, .(
    binding_raw = sum(peak_share * reliability),
    n_unique_peaks = data.table::uniqueN(peak_id)
  ), by = .(condition_id, tf, topic_num)]
  peak_background <- unique(peak_topic[, .(
    condition_id,
    topic_num = peak_topic,
    peak_id,
    peak_share
  )])[, .(
    peak_background = sum(peak_share)
  ), by = .(condition_id, topic_num)]
  target_support <- links[, .(
    strongest_peak_support = max(peak_share * reliability),
    gene_share = max(gene_share)
  ), by = .(condition_id, tf, topic_num = gene_topic, gene_key)]
  target <- target_support[, .(
    target_raw = sum(gene_share * strongest_peak_support),
    n_unique_targets = data.table::uniqueN(gene_key)
  ), by = .(condition_id, tf, topic_num)]
  gene_background <- unique(gene_topic[, .(
    condition_id,
    topic_num = gene_topic,
    gene_key,
    gene_share
  )])[, .(
    gene_background = sum(gene_share)
  ), by = .(condition_id, topic_num)]
  tf_condition <- unique(links[, .(condition_id, tf)])
  grid <- tf_condition[, .(
    topic_num = seq_len(as.integer(topic_count))
  ), by = .(condition_id, tf)]
  evidence <- binding[grid, on = c("condition_id", "tf", "topic_num")]
  evidence <- target[
    evidence,
    on = c("condition_id", "tf", "topic_num")
  ]
  evidence <- peak_background[
    evidence,
    on = c("condition_id", "topic_num")
  ]
  evidence <- gene_background[
    evidence,
    on = c("condition_id", "topic_num")
  ]
  for (column in c(
    "binding_raw", "target_raw", "n_unique_peaks", "n_unique_targets",
    "peak_background", "gene_background"
  )) {
    data.table::set(
      evidence,
      which(!is.finite(evidence[[column]])),
      column,
      0
    )
  }
  evidence <- .module3_tf_topic_component_normalize(
    evidence,
    "binding_raw",
    "peak_background",
    "binding_probability"
  )
  evidence <- .module3_tf_topic_component_normalize(
    evidence,
    "target_raw",
    "gene_background",
    "target_probability"
  )
  evidence[, component_available__ := target_raw > 0]
  evidence[, combined_probability := if (any(component_available__)) {
    binding_weight * binding_probability +
      target_weight * target_probability
  } else {
    binding_probability
  },
  by = .(condition_id, tf)]
  evidence[, combined_probability := {
    total <- sum(combined_probability)
    if (is.finite(total) && total > 0) combined_probability / total else 0
  }, by = .(condition_id, tf)]
  evidence[, component_available__ := NULL]
  data.table::setorder(
    evidence,
    condition_id,
    tf,
    -combined_probability,
    topic_num
  )
  evidence[, `:=`(
    condition_primary_topic = topic_num[[1L]],
    condition_primary_probability = combined_probability[[1L]],
    condition_second_probability =
      if (.N > 1L) combined_probability[[2L]] else 0
  ), by = .(condition_id, tf)]
  evidence[, membership_pass :=
    combined_probability >= membership_cutoff]
  global <- evidence[, .(
    global_probability = mean(combined_probability)
  ), by = .(tf, topic_num)]
  data.table::setorder(global, tf, -global_probability, topic_num)
  global[, `:=`(
    global_primary_topic = topic_num[[1L]],
    global_primary_probability = global_probability[[1L]],
    global_second_probability =
      if (.N > 1L) global_probability[[2L]] else 0
  ), by = tf]
  global[, global_primary_margin :=
    global_primary_probability - global_second_probability]
  evidence <- global[
    evidence,
    on = c("tf", "topic_num")
  ]
  evidence[, primary_confident :=
    topic_num == global_primary_topic &
      global_primary_probability >= membership_cutoff &
      global_primary_margin >= primary_margin_cutoff]
  evidence[, `:=`(
    topic = paste0("Topic", topic_num),
    theta = combined_probability,
    comparison_id = condition_id,
    direction = "All",
    primary_topic_num = as.integer(global_primary_topic),
    primary_topic = paste0("Topic", global_primary_topic),
    assignment_source = "binding_target_evidence"
  )]
  data.table::setorder(evidence, condition_id, tf, topic_num)
  list(evidence = evidence[], join_audit = join_audit[])
}
