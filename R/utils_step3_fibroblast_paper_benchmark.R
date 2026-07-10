# File: utils_step3_fibroblast_paper_benchmark.R
# Author: Yaoxiang Li
# Purpose: Fibroblast-only Tsao 2022 topic benchmark helpers.

.m3fb_default_comparison_map <- function() {
  data.table::data.table(
    comparison_id = c(
      "Fib_BI_vs_NoTF",
      "Fib_BIRT_Dox",
      "Fib_BIRT_vs_NoTF",
      "Fib_BIRT_vs_BI",
      "Fib_BIT_vs_BI",
      "Fib_BIR_vs_BI"
    ),
    direction = "Target-Up",
    required_tfs = c(
      "Batf;Irf4",
      "Batf;Irf4;Runx3;Tbx21",
      "Batf;Irf4;Runx3;Tbx21",
      "Runx3;Tbx21",
      "Tbx21",
      "Runx3"
    ),
    target_regulators = c(
      "Batf_Irf4_module",
      "Batf_Irf4_Runx3_Tbx21_module",
      "Batf_Irf4_Runx3_Tbx21_module",
      "Batf_Irf4_Runx3_Tbx21_module",
      "Batf_Irf4_Runx3_Tbx21_module",
      "Batf_Irf4_Runx3_Tbx21_module"
    ),
    comparison_weight = c(1.0, 1.5, 1.5, 1.0, 1.0, 1.0)
  )
}

.m3fb_split_semicolon <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  parts <- strsplit(x, ";", fixed = TRUE)
  lapply(parts, function(v) {
    v <- trimws(v)
    unique(v[nzchar(v)])
  })
}

.m3fb_gene_key <- function(x) {
  x <- as.character(x)
  x <- sub("^[A-Za-z]+:", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("_.*$", "", x)
  x <- sub(";.*$", "", x)
  x
}

.m3fb_load_tsao_fibroblast_benchmark <- function(path,
                                                  confidence = c("gold", "silver", "bronze"),
                                                  fibroblast_context = "NIH3T3_fibroblast_TF_overexpression_72h") {
  if (!file.exists(path)) .log_abort("Tsao 2022 benchmark CSV not found: {path}")
  dt <- data.table::fread(path, showProgress = FALSE)
  required <- c("row_type", "regulator", "target", "target_type", "effect", "sign_num", "context", "confidence")
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    .log_abort("Tsao 2022 benchmark CSV is missing required columns: {paste(missing, collapse = ', ')}")
  }
  dt[, `:=`(
    row_type = as.character(row_type),
    regulator = as.character(regulator),
    target = as.character(target),
    target_type = as.character(target_type),
    effect = as.character(effect),
    sign_num = suppressWarnings(as.numeric(sign_num)),
    context = as.character(context),
    confidence = tolower(as.character(confidence))
  )]
  confidence_keep <- tolower(as.character(confidence))
  members <- dt[
    row_type == "module_member" &
      context == fibroblast_context &
      target_type == "TF" &
      !is.na(target) & nzchar(target),
    .(regulator, target, target_key = toupper(target), confidence)
  ]
  members <- unique(members, by = c("regulator", "target_key"))

  targets <- dt[
    row_type == "regulatory_edge" &
      context == fibroblast_context &
      target_type == "gene" &
      confidence %in% confidence_keep &
      (effect %in% c("up", "induced", "positive") | sign_num > 0) &
      !is.na(target) & nzchar(target),
    .(regulator, target, target_key = toupper(target), confidence, effect, sign_num)
  ]
  targets <- unique(targets, by = c("regulator", "target_key", "confidence"))
  list(module_members = members, targets = targets, source = normalizePath(path, winslash = "/", mustWork = FALSE))
}

.m3fb_topic_nums_from_cols <- function(x) {
  out <- suppressWarnings(as.integer(sub("^Topic", "", as.character(x))))
  bad <- !is.finite(out)
  out[bad] <- seq_along(out)[bad]
  out
}

.m3fb_prepare_phi_gene_scores <- function(phi) {
  phi <- as.matrix(phi)
  storage.mode(phi) <- "numeric"
  if (is.null(rownames(phi))) rownames(phi) <- paste0("Topic", seq_len(nrow(phi)))
  if (is.null(colnames(phi))) colnames(phi) <- paste0("term_", seq_len(ncol(phi)))
  row_topic <- grepl("^Topic[0-9]+$", rownames(phi))
  col_topic <- grepl("^Topic[0-9]+$", colnames(phi))
  if (sum(col_topic) > sum(row_topic)) phi <- t(phi)
  phi[!is.finite(phi) | phi < 0] <- 0
  col_max <- apply(phi, 2L, max, na.rm = TRUE)
  col_max[!is.finite(col_max) | col_max <= 0] <- 1
  scaled <- sweep(phi, 2L, col_max, "/")
  topics <- rownames(scaled)
  topic_nums <- .m3fb_topic_nums_from_cols(topics)
  terms <- colnames(scaled)
  genes <- .m3fb_gene_key(terms)
  long <- data.table::as.data.table(as.table(scaled))
  data.table::setnames(long, c("topic", "term_id", "gene_score"))
  long[, `:=`(
    topic = as.character(topic),
    term_id = as.character(term_id),
    topic_num = topic_nums[match(as.character(topic), topics)],
    gene = genes[match(as.character(term_id), terms)],
    gene_key = toupper(genes[match(as.character(term_id), terms)])
  )]
  long <- long[nzchar(gene_key)]
  long[, .(
    gene_score = max(as.numeric(gene_score), na.rm = TRUE),
    term_ids = paste(unique(term_id[as.numeric(gene_score) == max(as.numeric(gene_score), na.rm = TRUE)]), collapse = ";")
  ), by = .(topic, topic_num, gene, gene_key)]
}

.m3fb_parse_theta_docs <- function(theta) {
  theta <- as.matrix(theta)
  storage.mode(theta) <- "numeric"
  if (is.null(colnames(theta))) colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  if (is.null(rownames(theta))) rownames(theta) <- paste0("doc_", seq_len(nrow(theta)))
  docs <- data.table::data.table(doc_id = rownames(theta))
  parts <- strsplit(docs$doc_id, "::", fixed = TRUE)
  is_direction <- function(x) x %in% c("Target-Up", "Target-Down", "Up", "Down")
  docs[, `:=`(
    comparison_id = vapply(parts, function(x) if (length(x) >= 1L) x[[1L]] else NA_character_, character(1L)),
    direction = vapply(parts, function(x) {
      if (length(x) < 2L) return(NA_character_)
      if (is_direction(x[[2L]])) return(x[[2L]])
      if (length(x) >= 3L && is_direction(x[[3L]])) return(x[[3L]])
      NA_character_
    }, character(1L)),
    tf = vapply(parts, function(x) {
      if (length(x) < 3L) return(NA_character_)
      if (is_direction(x[[2L]])) return(x[[3L]])
      x[[2L]]
    }, character(1L))
  )]
  topics <- colnames(theta)
  topic_nums <- .m3fb_topic_nums_from_cols(topics)
  long <- data.table::as.data.table(as.table(theta))
  data.table::setnames(long, c("doc_id", "topic", "theta"))
  long[, `:=`(
    doc_id = as.character(doc_id),
    topic = as.character(topic),
    theta = as.numeric(theta),
    topic_num = topic_nums[match(as.character(topic), topics)]
  )]
  merge(long, docs, by = "doc_id", all.x = TRUE, sort = FALSE)
}

.m3fb_required_tfs_for_comparison <- function(comparison_id, comparison_map = .m3fb_default_comparison_map()) {
  comparison_key <- as.character(comparison_id)[[1L]]
  hit <- comparison_map[comparison_id == comparison_key]
  if (!nrow(hit)) return(character())
  .m3fb_split_semicolon(hit$required_tfs[[1L]])[[1L]]
}

.m3fb_targets_for_comparison <- function(benchmark,
                                         comparison_id,
                                         comparison_map = .m3fb_default_comparison_map()) {
  targets <- data.table::copy(data.table::as.data.table(benchmark$targets))
  if (!nrow(targets)) return(targets)
  comparison_key <- as.character(comparison_id)[[1L]]
  hit <- comparison_map[comparison_id == comparison_key]
  if (!nrow(hit) || !"target_regulators" %in% names(hit)) return(targets)
  regulators <- .m3fb_split_semicolon(hit$target_regulators[[1L]])[[1L]]
  if (!length(regulators)) return(targets)
  targets[regulator %in% regulators]
}

.m3fb_score_model_row <- function(theta,
                                  phi,
                                  row,
                                  benchmark,
                                  comparisons = NULL,
                                  comparison_map = .m3fb_default_comparison_map(),
                                  theta_cutoff = 0.3,
                                  gene_score_cutoff = 0.7,
                                  max_story_topics = 3L) {
  row <- data.table::as.data.table(row)
  if (!nrow(row)) .log_abort("`row` must contain one model metadata row.")
  if (is.null(comparisons)) comparisons <- comparison_map$comparison_id
  comparison_map <- comparison_map[comparison_id %in% comparisons]
  if (!nrow(comparison_map)) {
    empty <- data.table::data.table()
    return(list(comparison_scores = empty, selected_topics = empty, membership_long = empty))
  }
  theta_long <- .m3fb_parse_theta_docs(theta)
  phi_genes <- .m3fb_prepare_phi_gene_scores(phi)
  targets <- data.table::copy(data.table::as.data.table(benchmark$targets))
  targets[, target_key := toupper(as.character(target))]
  targets[, confidence_weight := 1.0]
  all_scores <- list()
  all_topics <- list()
  all_membership <- list()

  for (i in seq_len(nrow(comparison_map))) {
    cmp <- comparison_map$comparison_id[[i]]
    direction <- comparison_map$direction[[i]]
    required_tfs <- .m3fb_split_semicolon(comparison_map$required_tfs[[i]])[[1L]]
    required_keys <- toupper(required_tfs)
    cmp_targets <- .m3fb_targets_for_comparison(list(targets = targets), cmp, comparison_map = comparison_map)
    cmp_theta <- theta_long[comparison_id == cmp & direction == direction & toupper(tf) %in% required_keys]
    tf_topic <- cmp_theta[, .(
      tf_theta = max(theta, na.rm = TRUE),
      tf_present = any(theta >= theta_cutoff, na.rm = TRUE)
    ), by = .(topic, topic_num, tf, tf_key = toupper(tf))]
    tf_topic <- tf_topic[tf_present == TRUE]
    gene_topic <- merge(
      phi_genes,
      cmp_targets[, .(target, target_key, confidence, confidence_weight)],
      by.x = "gene_key",
      by.y = "target_key",
      all = FALSE,
      sort = FALSE,
      allow.cartesian = TRUE
    )
    gene_topic[, gene_present := gene_score >= gene_score_cutoff]
    gene_topic <- gene_topic[gene_present == TRUE]

    topic_summary <- merge(
      tf_topic[, .(
        required_tfs_covered = paste(sort(unique(tf)), collapse = ";"),
        n_required_tfs_covered = data.table::uniqueN(tf_key),
        mean_tf_theta = mean(tf_theta, na.rm = TRUE)
      ), by = .(topic, topic_num)],
      gene_topic[, .(
        targets_covered = paste(sort(unique(target)), collapse = ";"),
        n_targets_covered = data.table::uniqueN(gene_key),
        weighted_targets_covered = sum(unique(confidence_weight), na.rm = TRUE),
        mean_gene_score = mean(gene_score, na.rm = TRUE)
      ), by = .(topic, topic_num)],
      by = c("topic", "topic_num"),
      all = FALSE,
      sort = FALSE
    )
    if (!nrow(topic_summary)) {
      topic_summary <- data.table::data.table(
        topic = character(),
        topic_num = integer(),
        required_tfs_covered = character(),
        n_required_tfs_covered = integer(),
        mean_tf_theta = numeric(),
        targets_covered = character(),
        n_targets_covered = integer(),
        weighted_targets_covered = numeric(),
        mean_gene_score = numeric()
      )
    }
    topic_summary[, `:=`(
      comparison_id = cmp,
      direction = direction,
      n_required_tfs = length(required_keys),
      n_targets_total = data.table::uniqueN(cmp_targets$target_key),
      weighted_targets_total = sum(unique(cmp_targets[, .(target_key, confidence_weight)])$confidence_weight, na.rm = TRUE),
      comparison_weight = as.numeric(comparison_map$comparison_weight[[i]])
    )]
    topic_summary[, topic_rank_score := n_targets_covered * 1000 + n_required_tfs_covered * 100 + mean_tf_theta + mean_gene_score]
    data.table::setorder(topic_summary, -topic_rank_score, topic_num)
    selected <- topic_summary[seq_len(min(max_story_topics, .N))]
    if (nrow(selected)) {
      selected_target_keys <- unique(unlist(.m3fb_split_semicolon(selected$targets_covered), use.names = FALSE))
      selected_tf_keys <- unique(unlist(.m3fb_split_semicolon(selected$required_tfs_covered), use.names = FALSE))
      selected_target_keys <- toupper(selected_target_keys)
      selected_tf_keys <- toupper(selected_tf_keys)
      n_targets <- sum(unique(cmp_targets[, .(target_key)])$target_key %in% selected_target_keys)
      weighted_targets <- sum(unique(cmp_targets[target_key %in% selected_target_keys, .(target_key, confidence_weight)])$confidence_weight, na.rm = TRUE)
      n_tfs <- sum(required_keys %in% selected_tf_keys)
    } else {
      n_targets <- 0L
      weighted_targets <- 0
      n_tfs <- 0L
    }
    tf_fraction <- if (length(required_keys)) n_tfs / length(required_keys) else NA_real_
    target_fraction <- if (data.table::uniqueN(cmp_targets$target_key)) n_targets / data.table::uniqueN(cmp_targets$target_key) else NA_real_
    weighted_target_fraction <- if (sum(unique(cmp_targets[, .(target_key, confidence_weight)])$confidence_weight, na.rm = TRUE) > 0) {
      weighted_targets / sum(unique(cmp_targets[, .(target_key, confidence_weight)])$confidence_weight, na.rm = TRUE)
    } else {
      NA_real_
    }
    model_score <- 0.60 * weighted_target_fraction + 0.25 * target_fraction + 0.15 * tf_fraction
    comp_score <- data.table::data.table(
      comparison_id = cmp,
      direction = direction,
      method_order = row$method_order[[1L]],
      method = if ("method" %in% names(row)) row$method[[1L]] else NA_character_,
      method_setup = if ("method_setup" %in% names(row)) row$method_setup[[1L]] else NA_character_,
      model_label = if ("model_label" %in% names(row)) row$model_label[[1L]] else NA_character_,
      k = as.integer(row$selected_k[[1L]]),
      n_selected_topics = nrow(selected),
      selected_topics = paste(selected$topic_num, collapse = ";"),
      n_required_tfs = length(required_keys),
      n_required_tfs_covered = n_tfs,
      required_tf_fraction = tf_fraction,
      n_targets_total = data.table::uniqueN(cmp_targets$target_key),
      n_targets_cotopic = n_targets,
      target_fraction = target_fraction,
      weighted_target_fraction = weighted_target_fraction,
      comparison_weight = as.numeric(comparison_map$comparison_weight[[i]]),
      model_score = model_score
    )
    if (nrow(selected)) selected[, `:=`(method_setup = comp_score$method_setup[[1L]], model_label = comp_score$model_label[[1L]], k = comp_score$k[[1L]])]
    all_scores[[length(all_scores) + 1L]] <- comp_score
    all_topics[[length(all_topics) + 1L]] <- selected
    membership <- data.table::rbindlist(list(
      tf_topic[, .(comparison_id = cmp, direction = direction, topic, topic_num, item_type = "TF", item = tf, score = tf_theta, pass = TRUE)],
      gene_topic[, .(comparison_id = cmp, direction = direction, topic, topic_num, item_type = "Target gene", item = target, score = gene_score, pass = TRUE)]
    ), use.names = TRUE, fill = TRUE)
    if (nrow(membership)) {
      membership[, `:=`(
        method_order = row$method_order[[1L]],
        method = if ("method" %in% names(row)) row$method[[1L]] else NA_character_,
        method_setup = comp_score$method_setup[[1L]],
        model_label = comp_score$model_label[[1L]],
        k = comp_score$k[[1L]]
      )]
    }
    all_membership[[length(all_membership) + 1L]] <- membership
  }
  list(
    comparison_scores = data.table::rbindlist(all_scores, use.names = TRUE, fill = TRUE),
    selected_topics = data.table::rbindlist(all_topics, use.names = TRUE, fill = TRUE),
    membership_long = data.table::rbindlist(all_membership, use.names = TRUE, fill = TRUE)
  )
}

.m3fb_score_existing_models <- function(output_dir,
                                        model_rows,
                                        benchmark_csv,
                                        review_dir,
                                        confidence = c("gold", "silver", "bronze"),
                                        theta_cutoff = 0.3,
                                        gene_score_cutoff = 0.7,
                                        max_story_topics = 3L,
                                        top_n_models = 5L,
                                        comparisons = NULL,
                                        verbose = TRUE) {
  model_rows <- data.table::as.data.table(model_rows)
  if (!nrow(model_rows)) return(NULL)
  benchmark <- .m3fb_load_tsao_fibroblast_benchmark(benchmark_csv, confidence = confidence)
  benchmark_all <- .m3fb_load_tsao_fibroblast_benchmark(benchmark_csv, confidence = c("gold", "silver", "bronze"))
  comparison_map <- .m3fb_default_comparison_map()
  if (!is.null(comparisons)) comparison_map <- comparison_map[comparison_id %in% comparisons]
  out_dir <- file.path(review_dir, "fibroblast_paper_topic_recovery")
  tables_dir <- file.path(out_dir, "tables")
  figures_dir <- file.path(out_dir, "figures")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(comparison_map, file.path(tables_dir, "fibroblast_expected_comparison_direction_map.csv"))
  data.table::fwrite(benchmark$module_members, file.path(tables_dir, "fibroblast_benchmark_module_members.csv"))
  data.table::fwrite(benchmark$targets, file.path(tables_dir, "fibroblast_benchmark_targets.csv"))
  data.table::fwrite(benchmark_all$targets, file.path(tables_dir, "fibroblast_benchmark_targets_all_confidence.csv"))

  results <- vector("list", nrow(model_rows))
  for (i in seq_len(nrow(model_rows))) {
    row <- model_rows[i]
    k <- as.integer(row$selected_k[[1L]])
    theta_file <- file.path(row$model_dir[[1L]], "vae_models", sprintf("theta_K%d.csv", k))
    phi_file <- file.path(row$model_dir[[1L]], "vae_models", sprintf("phi_K%d.csv", k))
    if (!file.exists(theta_file) || !file.exists(phi_file)) next
    if (isTRUE(verbose) && i == 1L) {
      .log_inform("Scoring fibroblast Tsao 2022 topic recovery from existing theta/phi files.")
    }
    theta <- .m3tb_read_probability_csv(theta_file, "doc_id")
    phi <- .m3tb_read_probability_csv(phi_file, "term_id")
    results[[i]] <- .m3fb_score_model_row(
      theta = theta,
      phi = phi,
      row = row,
      benchmark = benchmark,
      comparisons = comparison_map$comparison_id,
      comparison_map = comparison_map,
      theta_cutoff = theta_cutoff,
      gene_score_cutoff = gene_score_cutoff,
      max_story_topics = max_story_topics
    )
  }
  results <- results[!vapply(results, is.null, logical(1L))]
  if (!length(results)) return(NULL)
  comparison_scores <- data.table::rbindlist(lapply(results, `[[`, "comparison_scores"), use.names = TRUE, fill = TRUE)
  selected_topics <- data.table::rbindlist(lapply(results, `[[`, "selected_topics"), use.names = TRUE, fill = TRUE)
  membership_long <- data.table::rbindlist(lapply(results, `[[`, "membership_long"), use.names = TRUE, fill = TRUE)
  if (nrow(comparison_scores)) {
    leaderboard <- comparison_scores[, .(
      global_score = stats::weighted.mean(model_score, w = comparison_weight, na.rm = TRUE),
      mean_target_fraction = stats::weighted.mean(target_fraction, w = comparison_weight, na.rm = TRUE),
      mean_weighted_target_fraction = stats::weighted.mean(weighted_target_fraction, w = comparison_weight, na.rm = TRUE),
      mean_required_tf_fraction = stats::weighted.mean(required_tf_fraction, w = comparison_weight, na.rm = TRUE),
      total_selected_topics = sum(n_selected_topics, na.rm = TRUE),
      n_comparisons_scored = .N
    ), by = .(method_order, method, method_setup, model_label, k)]
    leaderboard[, multivi_preference := grepl("MultiVI", model_label, fixed = TRUE)]
    data.table::setorder(
      leaderboard,
      -global_score,
      -mean_weighted_target_fraction,
      total_selected_topics,
      -multivi_preference,
      k,
      method_order
    )
    leaderboard[, rank := seq_len(.N)]
  } else {
    leaderboard <- data.table::data.table()
  }
  full_membership <- .m3fb_full_membership_for_top_models(
    model_rows = model_rows,
    leaderboard = leaderboard,
    benchmark = benchmark_all,
    comparison_map = comparison_map,
    theta_cutoff = theta_cutoff,
    gene_score_cutoff = gene_score_cutoff,
    top_n_models = top_n_models
  )
  missing_genes <- .m3fb_missing_benchmark_genes(model_rows, benchmark_all)
  data.table::fwrite(leaderboard, file.path(tables_dir, "fibroblast_model_k_leaderboard.csv"))
  data.table::fwrite(comparison_scores, file.path(tables_dir, "fibroblast_comparison_topic_scores.csv"))
  data.table::fwrite(selected_topics, file.path(tables_dir, "fibroblast_selected_story_topics.csv"))
  data.table::fwrite(membership_long, file.path(tables_dir, "fibroblast_topic_gene_tf_membership_long.csv"))
  data.table::fwrite(full_membership, file.path(tables_dir, "fibroblast_expected_direction_topic_membership.csv"))
  data.table::fwrite(full_membership, file.path(tables_dir, "fibroblast_expected_direction_topic_membership_top_models.csv"))
  data.table::fwrite(missing_genes, file.path(tables_dir, "fibroblast_missing_benchmark_genes.csv"))
  .m3fb_plot_leaderboard(leaderboard, file.path(figures_dir, "fibroblast_model_k_leaderboard.pdf"))
  .m3fb_plot_story_heatmaps(
    leaderboard,
    full_membership,
    file.path(figures_dir, "fibroblast_story_topic_heatmaps.pdf"),
    top_n_models = top_n_models
  )
  .m3fb_write_index_html(
    out_dir = out_dir,
    leaderboard = leaderboard,
    comparison_scores = comparison_scores,
    selected_topics = selected_topics,
    full_membership = full_membership,
    model_rows = model_rows,
    theta_cutoff = theta_cutoff,
    gene_score_cutoff = gene_score_cutoff,
    benchmark_source = benchmark$source
  )
  subgrn_dir <- .m3fb_write_paper_subgrn(
    out_dir = file.path(out_dir, "paper_topic_subgrn"),
    leaderboard = leaderboard,
    selected_topics = selected_topics,
    model_rows = model_rows,
    full_membership = full_membership,
    theta_cutoff = theta_cutoff,
    gene_score_cutoff = gene_score_cutoff,
    top_n_models = top_n_models
  )
  list(
    leaderboard = leaderboard,
    comparison_scores = comparison_scores,
    selected_topics = selected_topics,
    membership_long = membership_long,
    full_membership = full_membership,
    missing_genes = missing_genes,
    subgrn_dir = subgrn_dir,
    output_dir = out_dir
  )
}

.m3fb_missing_benchmark_genes <- function(model_rows, benchmark) {
  model_rows <- data.table::as.data.table(model_rows)
  targets <- unique(data.table::as.data.table(benchmark$targets)[, .(target, target_key = toupper(target))])
  if (!nrow(model_rows) || !nrow(targets)) return(data.table::data.table())
  row <- model_rows[1]
  k <- as.integer(row$selected_k[[1L]])
  phi_file <- file.path(row$model_dir[[1L]], "vae_models", sprintf("phi_K%d.csv", k))
  if (!file.exists(phi_file)) return(data.table::data.table())
  header <- names(data.table::fread(phi_file, nrows = 0L, showProgress = FALSE))
  term_cols <- setdiff(header, header[[1L]])
  term_keys <- unique(toupper(.m3fb_gene_key(term_cols)))
  targets[, present_in_model_terms := target_key %in% term_keys]
  targets[]
}

.m3fb_full_membership_for_best_model <- function(model_rows,
                                                 leaderboard,
                                                 benchmark,
                                                 comparison_map,
                                                 theta_cutoff = 0.3,
                                                 gene_score_cutoff = 0.7) {
  .m3fb_full_membership_for_top_models(
    model_rows = model_rows,
    leaderboard = leaderboard,
    benchmark = benchmark,
    comparison_map = comparison_map,
    theta_cutoff = theta_cutoff,
    gene_score_cutoff = gene_score_cutoff,
    top_n_models = 1L
  )
}

.m3fb_full_membership_for_top_models <- function(model_rows,
                                                 leaderboard,
                                                 benchmark,
                                                 comparison_map,
                                                 theta_cutoff = 0.3,
                                                 gene_score_cutoff = 0.7,
                                                 top_n_models = 5L) {
  if (!nrow(leaderboard)) return(data.table::data.table())
  top_n_models <- max(1L, as.integer(top_n_models[[1L]]))
  top <- data.table::copy(leaderboard[rank <= top_n_models])
  if (!nrow(top)) return(data.table::data.table())
  model_rows <- data.table::as.data.table(model_rows)
  out <- vector("list", nrow(top))
  for (i in seq_len(nrow(top))) {
    hit <- top[i]
    row <- model_rows[
      method_setup == hit$method_setup[[1L]] &
        as.integer(selected_k) == as.integer(hit$k[[1L]])
    ]
    if (!nrow(row)) next
    row <- row[1]
    k <- as.integer(row$selected_k[[1L]])
    theta_file <- file.path(row$model_dir[[1L]], "vae_models", sprintf("theta_K%d.csv", k))
    phi_file <- file.path(row$model_dir[[1L]], "vae_models", sprintf("phi_K%d.csv", k))
    if (!file.exists(theta_file) || !file.exists(phi_file)) next
    one <- .m3fb_full_benchmark_membership(
      theta = .m3tb_read_probability_csv(theta_file, "doc_id"),
      phi = .m3tb_read_probability_csv(phi_file, "term_id"),
      row = row,
      benchmark = benchmark,
      comparisons = comparison_map$comparison_id,
      comparison_map = comparison_map,
      theta_cutoff = theta_cutoff,
      gene_score_cutoff = gene_score_cutoff
    )
    if (!nrow(one)) next
    one[, `:=`(
      model_rank = as.integer(hit$rank[[1L]]),
      global_score = if ("global_score" %in% names(hit)) as.numeric(hit$global_score[[1L]]) else NA_real_
    )]
    out[[i]] <- one
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

.m3fb_full_benchmark_membership <- function(theta,
                                            phi,
                                            row,
                                            benchmark,
                                            comparisons = NULL,
                                            comparison_map = .m3fb_default_comparison_map(),
                                            theta_cutoff = 0.3,
                                            gene_score_cutoff = 0.7) {
  row <- data.table::as.data.table(row)
  if (is.null(comparisons)) comparisons <- comparison_map$comparison_id
  comparison_map <- comparison_map[comparison_id %in% comparisons]
  if (!nrow(comparison_map)) return(data.table::data.table())
  theta_long <- .m3fb_parse_theta_docs(theta)
  phi_genes <- .m3fb_prepare_phi_gene_scores(phi)
  topics <- sort(unique(c(theta_long$topic_num, phi_genes$topic_num)))
  topics <- topics[is.finite(topics)]
  topic_table <- data.table::data.table(topic_num = as.integer(topics), topic = paste0("Topic", topics))
  rows <- list()
  for (i in seq_len(nrow(comparison_map))) {
    cmp <- comparison_map$comparison_id[[i]]
    direction <- comparison_map$direction[[i]]
    required_tfs <- .m3fb_split_semicolon(comparison_map$required_tfs[[i]])[[1L]]
    if (length(required_tfs)) {
      tf_grid <- data.table::CJ(
        comparison_id = cmp,
        direction = direction,
        item = required_tfs,
        topic_num = topic_table$topic_num,
        unique = TRUE
      )
      tf_grid[, `:=`(
        topic = paste0("Topic", topic_num),
        item_key = toupper(item),
        item_type = "Required TF",
        confidence = "gold",
        regulator = NA_character_,
        expected_direction = direction
      )]
      tf_scores <- theta_long[
        comparison_id == cmp & direction == direction & toupper(tf) %in% toupper(required_tfs),
        .(score = max(theta, na.rm = TRUE)),
        by = .(item_key = toupper(tf), topic_num)
      ]
      tf_grid <- merge(tf_grid, tf_scores, by = c("item_key", "topic_num"), all.x = TRUE, sort = FALSE)
      tf_grid[, pass := is.finite(score) & score >= theta_cutoff]
      rows[[length(rows) + 1L]] <- tf_grid
    }
    targets <- .m3fb_targets_for_comparison(benchmark, cmp, comparison_map = comparison_map)
    if (nrow(targets)) {
      target_grid <- data.table::CJ(
        comparison_id = cmp,
        direction = direction,
        item_key = unique(targets$target_key),
        topic_num = topic_table$topic_num,
        unique = TRUE
      )
      target_grid <- merge(
        target_grid,
        unique(targets[, .(item_key = target_key, item = target, item_type = "Target gene", confidence, regulator, expected_direction = direction)]),
        by = "item_key",
        all.x = TRUE,
        sort = FALSE
      )
      target_grid[, topic := paste0("Topic", topic_num)]
      gene_scores <- phi_genes[, .(score = max(gene_score, na.rm = TRUE)), by = .(item_key = gene_key, topic_num)]
      target_grid <- merge(target_grid, gene_scores, by = c("item_key", "topic_num"), all.x = TRUE, sort = FALSE)
      target_grid[, pass := is.finite(score) & score >= gene_score_cutoff]
      rows[[length(rows) + 1L]] <- target_grid
    }
  }
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(out)
  out[, `:=`(
    method_order = row$method_order[[1L]],
    method = if ("method" %in% names(row)) row$method[[1L]] else NA_character_,
    method_setup = if ("method_setup" %in% names(row)) row$method_setup[[1L]] else NA_character_,
    model_label = if ("model_label" %in% names(row)) row$model_label[[1L]] else NA_character_,
    k = as.integer(row$selected_k[[1L]])
  )]
  conf_order <- c(gold = 1L, silver = 2L, bronze = 3L)
  out[, confidence_order := data.table::fcoalesce(conf_order[tolower(confidence)], 9L)]
  data.table::setorder(out, comparison_id, item_type, confidence_order, item, topic_num)
  out[]
}

.m3fb_plot_leaderboard <- function(leaderboard, out_file) {
  if (!nrow(leaderboard)) return(invisible(NULL))
  top <- data.table::copy(leaderboard[rank <= 30L])
  top[, label := paste0(method_setup, " K", k)]
  top[, label := factor(label, levels = rev(label))]
  p <- ggplot2::ggplot(top, ggplot2::aes(x = label, y = global_score, fill = model_label)) +
    ggplot2::geom_col(width = 0.78, color = "black", linewidth = 0.2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Fibroblast paper topic recovery",
      x = "Method and K",
      y = "Global recovery score",
      fill = "Model"
    ) +
    ggplot2::theme_bw(base_family = "Helvetica", base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(out_file, p, width = 9.5, height = max(4.5, 0.22 * nrow(top) + 1.8), units = "in")
  invisible(out_file)
}

.m3fb_prepare_story_heatmap_data <- function(leaderboard, membership_long, top_n_models = 5L) {
  if (!nrow(leaderboard) || !nrow(membership_long)) return(data.table::data.table())
  top_n_models <- max(1L, as.integer(top_n_models[[1L]]))
  top <- data.table::copy(leaderboard[rank <= top_n_models])
  if (!nrow(top)) return(data.table::data.table())
  dt <- data.table::copy(membership_long)
  if (!"model_rank" %in% names(dt)) {
    rank_map <- top[, .(method_setup, k = as.integer(k), model_rank = as.integer(rank))]
    dt <- merge(dt, rank_map, by = c("method_setup", "k"), all.x = TRUE, sort = FALSE)
  }
  dt <- dt[model_rank %in% top$rank]
  if (!nrow(dt)) return(data.table::data.table())
  if (!"global_score" %in% names(dt)) {
    if (!"global_score" %in% names(top)) top[, global_score := NA_real_]
    score_map <- top[, .(method_setup, k = as.integer(k), global_score)]
    dt <- merge(dt, score_map, by = c("method_setup", "k"), all.x = TRUE, sort = FALSE)
  }
  if (!"confidence" %in% names(dt)) dt[, confidence := NA_character_]
  dt[, item_label := paste(item_type, item, sep = " ")]
  dt[, item_color := data.table::fcase(
    item_type == "Required TF", "#2563eb",
    item_type == "Target gene", "#991b1b",
    default = "#111111"
  )]
  dt[, topic_label := paste0("Topic ", topic_num)]
  dt[, comparison_label := paste(comparison_id, direction, sep = "::")]
  dt[, score_plot := data.table::fcoalesce(as.numeric(score), 0)]
  dt[, model_label_full := paste0("Rank ", model_rank, ": ", method_setup, " K", k)]
  data.table::setorder(dt, model_rank, comparison_id, direction, item_type, confidence_order, item, topic_num)
  dt[]
}

.m3fb_plot_story_heatmaps <- function(leaderboard, membership_long, out_file, top_n_models = 5L) {
  dt <- .m3fb_prepare_story_heatmap_data(leaderboard, membership_long, top_n_models = top_n_models)
  if (!nrow(dt)) return(invisible(NULL))
  grDevices::pdf(out_file, width = 13, height = 10.5, family = "Helvetica", onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  page_keys <- unique(dt[, .(model_rank, comparison_label)])
  data.table::setorder(page_keys, model_rank, comparison_label)
  for (i in seq_len(nrow(page_keys))) {
    key <- page_keys[i]
    page <- dt[model_rank == key$model_rank[[1L]] & comparison_label == key$comparison_label[[1L]]]
    item_levels <- unique(page[order(item_type, confidence_order, item)]$item_label)
    item_colors <- unique(page[, .(item_label, item_color)])
    item_colors <- item_colors[match(item_levels, item_label)]$item_color
    page[, item_label := factor(item_label, levels = rev(item_levels))]
    page[, topic_label := factor(topic_label, levels = paste0("Topic ", sort(unique(topic_num))))]
    page_model <- page[1]
    axis_colors <- rev(item_colors)
    p <- ggplot2::ggplot(page, ggplot2::aes(x = topic_label, y = item_label, fill = score_plot)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(pass, "*", "")),
        size = 2.5,
        fontface = "bold",
        color = "black"
      ) +
      ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#cc2f3c", limits = c(0, 1), oob = scales::squish) +
      ggplot2::labs(
        title = paste0(page_model$comparison_label[[1L]], " benchmark TF and target topic scores"),
        subtitle = paste0(
          page_model$model_label_full[[1L]],
          if (is.finite(page_model$global_score[[1L]])) sprintf(", score %.3f", page_model$global_score[[1L]]) else "",
          "; * marks score passing the assignment cutoff"
        ),
        x = "Topic",
        y = "Expected paper TF or target gene",
        fill = "Score"
      ) +
      ggplot2::theme_bw(base_family = "Helvetica", base_size = 9) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(face = "bold", size = 9),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.text.y = ggplot2::element_text(face = "bold", color = axis_colors),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(face = "bold")
      )
    print(p)
  }
  invisible(out_file)
}

.m3fb_write_index_html <- function(out_dir,
                                   leaderboard,
                                   comparison_scores,
                                   selected_topics,
                                   full_membership,
                                   model_rows,
                                   theta_cutoff,
                                   gene_score_cutoff,
                                   benchmark_source) {
  best <- if (nrow(leaderboard)) leaderboard[rank == 1L] else data.table::data.table()
  best_text <- if (nrow(best)) {
    sprintf(
      "%s, K%d, score %.3f",
      best$method_setup[[1L]],
      as.integer(best$k[[1L]]),
      as.numeric(best$global_score[[1L]])
    )
  } else {
    "No scored model"
  }
  top_rows <- if (nrow(leaderboard)) {
    rows <- leaderboard[rank <= 15L]
    paste0(
      "<tr><td>", rows$rank, "</td><td>", .m3tb_html_escape(rows$method_setup),
      "</td><td>K", rows$k, "</td><td>", sprintf("%.3f", rows$global_score),
      "</td><td>", sprintf("%.1f%%", 100 * rows$mean_weighted_target_fraction),
      "</td><td>", sprintf("%.1f%%", 100 * rows$mean_required_tf_fraction),
      "</td></tr>",
      collapse = ""
    )
  } else {
    "<tr><td colspan=\"6\">No rows</td></tr>"
  }
  story_rows <- if (nrow(selected_topics)) {
    rows <- selected_topics
    if (nrow(best)) rows <- rows[method_setup == best$method_setup[[1L]] & as.integer(k) == as.integer(best$k[[1L]])]
    rows <- rows[order(comparison_id, topic_num)]
    paste0(
      "<tr><td>", .m3tb_html_escape(rows$comparison_id), "::", .m3tb_html_escape(rows$direction),
      "</td><td>Topic ", rows$topic_num,
      "</td><td>", .m3tb_html_escape(rows$required_tfs_covered),
      "</td><td>", .m3tb_html_escape(rows$targets_covered),
      "</td></tr>",
      collapse = ""
    )
  } else {
    "<tr><td colspan=\"4\">No selected topics</td></tr>"
  }
  links <- .m3fb_condition_page_links(best, model_rows, out_dir)
  link_rows <- if (nrow(links)) {
    paste0(
      "<li><a href=\"", .m3tb_html_escape(links$rel_path), "\">",
      .m3tb_html_escape(links$label), "</a></li>",
      collapse = ""
    )
  } else {
    "<li>No existing condition-topic page was found for the selected model.</li>"
  }
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Fibroblast paper topic recovery</title>",
    "<style>body{font-family:Arial,Helvetica,sans-serif;font-weight:700;margin:24px;background:#fbfbf8;color:#111}h1{font-size:24px}h2{font-size:18px;margin-top:24px}table{border-collapse:collapse;width:100%;background:#fff}th,td{border:1px solid #d6d6d0;padding:6px 8px;font-size:12px;vertical-align:top}th{background:#eef3f1;text-align:left}.meta{color:#475569;font-size:12px;line-height:1.45}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px}.card{border:1px solid #d6d6d0;background:#fff;padding:10px}.big{font-size:18px}</style>",
    "</head><body>",
    "<h1>Fibroblast paper topic recovery</h1>",
    paste0("<div class=\"cards\"><div class=\"card\"><div>Selected model</div><div class=\"big\">", .m3tb_html_escape(best_text), "</div></div>",
           "<div class=\"card\"><div>Theta cutoff</div><div class=\"big\">", theta_cutoff, "</div></div>",
           "<div class=\"card\"><div>Gene topic score cutoff</div><div class=\"big\">", gene_score_cutoff, "</div></div></div>"),
    paste0("<p class=\"meta\">Benchmark source: ", .m3tb_html_escape(benchmark_source), "</p>"),
    "<h2>Leaderboard</h2>",
    "<table><thead><tr><th>Rank</th><th>Method</th><th>K</th><th>Score</th><th>Target recovery</th><th>TF recovery</th></tr></thead><tbody>",
    top_rows,
    "</tbody></table>",
    "<h2>Selected Story Topics</h2>",
    "<table><thead><tr><th>Comparison</th><th>Topic</th><th>TFs</th><th>Target genes</th></tr></thead><tbody>",
    story_rows,
    "</tbody></table>",
    "<h2>Existing Interactive Pages</h2><ul>",
    link_rows,
    "</ul>",
    "<h2>Files</h2><ul>",
    "<li><a href=\"tables/fibroblast_model_k_leaderboard.csv\">fibroblast_model_k_leaderboard.csv</a></li>",
    "<li><a href=\"tables/fibroblast_comparison_topic_scores.csv\">fibroblast_comparison_topic_scores.csv</a></li>",
    "<li><a href=\"tables/fibroblast_expected_comparison_direction_map.csv\">fibroblast_expected_comparison_direction_map.csv</a></li>",
    "<li><a href=\"tables/fibroblast_selected_story_topics.csv\">fibroblast_selected_story_topics.csv</a></li>",
    "<li><a href=\"tables/fibroblast_expected_direction_topic_membership.csv\">fibroblast_expected_direction_topic_membership.csv</a></li>",
    "<li><a href=\"figures/fibroblast_model_k_leaderboard.pdf\">fibroblast_model_k_leaderboard.pdf</a></li>",
    "<li><a href=\"figures/fibroblast_story_topic_heatmaps.pdf\">fibroblast_story_topic_heatmaps.pdf</a></li>",
    "</ul>",
    "</body></html>"
  )
  out_file <- file.path(out_dir, "index.html")
  writeLines(html, out_file, useBytes = TRUE)
  invisible(out_file)
}

.m3fb_condition_page_links <- function(best, model_rows, out_dir) {
  if (!nrow(best)) return(data.table::data.table())
  rows <- data.table::as.data.table(model_rows)[
    method_setup == best$method_setup[[1L]] &
      as.integer(selected_k) == as.integer(best$k[[1L]])
  ]
  if (!nrow(rows)) return(data.table::data.table())
  row <- rows[1]
  slug <- .m3tb_review_report_slug(row, as.integer(row$selected_k[[1L]]))
  review_dir <- dirname(out_dir)
  candidates <- file.path(review_dir, "condition_topic_reports", "pages", sprintf("%s_condition_topic_report.html", slug))
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(data.table::data.table())
  data.table::data.table(
    label = sprintf("%s K%d condition-topic report", best$method_setup[[1L]], as.integer(best$k[[1L]])),
    path = hits,
    rel_path = .m3tb_relative_path(hits, out_dir)
  )
}

.m3fb_parse_doc_parts <- function(doc_id) {
  doc_id <- as.character(doc_id)
  parts <- strsplit(doc_id, "::", fixed = TRUE)
  is_direction <- function(x) x %in% c("Target-Up", "Target-Down", "Up", "Down")
  data.table::data.table(
    doc_id = doc_id,
    comparison_id = vapply(parts, function(x) if (length(x) >= 1L) x[[1L]] else NA_character_, character(1L)),
    direction = vapply(parts, function(x) {
      if (length(x) < 2L) return(NA_character_)
      if (is_direction(x[[2L]])) return(.m3tb_subgrn_direction_label(x[[2L]]))
      if (length(x) >= 3L && is_direction(x[[3L]])) return(.m3tb_subgrn_direction_label(x[[3L]]))
      NA_character_
    }, character(1L)),
    tf_doc = vapply(parts, function(x) {
      if (length(x) < 3L) return(NA_character_)
      if (is_direction(x[[2L]])) return(x[[3L]])
      x[[2L]]
    }, character(1L))
  )
}

.m3fb_write_paper_subgrn <- function(out_dir,
                                     leaderboard,
                                     selected_topics,
                                     model_rows,
                                     full_membership = NULL,
                                     theta_cutoff = 0.3,
                                     gene_score_cutoff = 0.7,
                                     top_n_models = 5L,
                                     max_links_per_context = 300L) {
  if (!nrow(leaderboard) || is.null(full_membership) || !nrow(full_membership)) return(NA_character_)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  top_n_models <- max(1L, as.integer(top_n_models[[1L]]))
  top_models <- data.table::copy(leaderboard[rank <= top_n_models])
  model_rows <- data.table::as.data.table(model_rows)
  membership <- data.table::copy(data.table::as.data.table(full_membership))
  membership <- membership[model_rank %in% top_models$rank]
  if (!nrow(membership)) return(NA_character_)
  membership[, gene_key_upper := toupper(as.character(item_key))]
  target_context <- membership[item_type == "Target gene", .(
    overlap_genes = paste(unique(item[order(confidence_order, item)]), collapse = ";"),
    n_overlap_genes = data.table::uniqueN(item_key)
  ), by = .(model_rank, method_setup, model_label, k, comparison_id, direction, topic_num)]
  if (!nrow(target_context)) return(NA_character_)
  target_context[, `:=`(
    subgrn_context_id = sprintf("paper_%05d", seq_len(.N)),
    comparison_label = paste(comparison_id, direction, sep = "::"),
    pathway = "Paper benchmark target genes",
    pathway_key = paste("paper", model_rank, comparison_id, direction, topic_num, sep = ":")
  )]
  selected_flag <- data.table::copy(data.table::as.data.table(selected_topics))
  if (nrow(selected_flag)) {
    selected_flag <- merge(
      selected_flag,
      unique(top_models[, .(method_setup, k = as.integer(k), model_rank = as.integer(rank))]),
      by = c("method_setup", "k"),
      all.x = TRUE,
      sort = FALSE
    )
    selected_flag <- unique(selected_flag[, .(model_rank, comparison_id, direction, topic_num = as.integer(topic_num), selected_story_topic = TRUE)])
    target_context <- merge(target_context, selected_flag, by = c("model_rank", "comparison_id", "direction", "topic_num"), all.x = TRUE, sort = FALSE)
  } else {
    target_context[, selected_story_topic := FALSE]
  }
  target_context[is.na(selected_story_topic), selected_story_topic := FALSE]

  target_keys <- unique(membership[item_type == "Target gene", .(
    model_rank, comparison_id, direction, topic_num, gene_key_upper,
    target_gene = item,
    gene_topic_score = data.table::fcoalesce(as.numeric(score), 0),
    gene_topic_pass = as.logical(pass)
  )])
  target_keys <- merge(
    target_keys,
    target_context[, .(subgrn_context_id, model_rank, comparison_id, direction, topic_num)],
    by = c("model_rank", "comparison_id", "direction", "topic_num"),
    all.x = TRUE,
    sort = FALSE
  )
  all_tf_gene <- list()
  all_triplets <- list()
  for (i in seq_len(nrow(top_models))) {
    hit <- top_models[i]
    row <- model_rows[
      method_setup == hit$method_setup[[1L]] &
        as.integer(selected_k) == as.integer(hit$k[[1L]])
    ]
    if (!nrow(row)) next
    row <- row[1]
    extraction_dir <- .m3tb_find_extraction_subdir(row)
    if (is.na(extraction_dir) || !dir.exists(extraction_dir)) next
    link_path <- file.path(extraction_dir, "topic_links_pass.csv")
    if (!file.exists(link_path)) link_path <- file.path(extraction_dir, "topic_links.csv")
    if (!file.exists(link_path)) next
    header <- names(data.table::fread(link_path, nrows = 0L, showProgress = FALSE))
    cols <- intersect(c("doc_id", "tf", "peak_id", "gene_key", "topic_num", "peak_score", "gene_score", "link_pass", "peak_pass", "gene_pass"), header)
    links <- data.table::fread(link_path, select = cols, showProgress = FALSE)
    if (!nrow(links) || !"doc_id" %in% names(links)) next
    for (col in c("tf", "peak_id", "gene_key", "topic_num", "peak_score", "gene_score")) {
      if (!col %in% names(links)) links[, (col) := NA]
    }
    doc_meta <- .m3fb_parse_doc_parts(unique(links$doc_id))
    links <- merge(links, doc_meta, by = "doc_id", all.x = TRUE, sort = FALSE)
    if (!"tf" %in% names(links)) links[, tf := tf_doc]
    links[is.na(tf) | !nzchar(tf), tf := tf_doc]
    links[, `:=`(
      model_rank = as.integer(hit$rank[[1L]]),
      comparison_id = as.character(comparison_id),
      direction = as.character(direction),
      tf = as.character(tf),
      tf_upper = toupper(as.character(tf)),
      gene_key = as.character(gene_key),
      gene_key_upper = toupper(as.character(gene_key)),
      peak_id = as.character(peak_id),
      link_topic_num = suppressWarnings(as.integer(topic_num)),
      edge_score_row = abs(data.table::fcoalesce(suppressWarnings(as.numeric(peak_score)), 0)) +
        abs(data.table::fcoalesce(suppressWarnings(as.numeric(gene_score)), 0))
    )]
    links[, topic_num := NULL]
    links <- links[!is.na(comparison_id) & nzchar(comparison_id) &
      !is.na(direction) & nzchar(direction) &
      !is.na(tf) & nzchar(tf) &
      !is.na(gene_key_upper) & nzchar(gene_key_upper)]
    contexts <- target_keys[model_rank == hit$rank[[1L]]]
    if (!nrow(contexts)) next
    data.table::setkey(links, model_rank, comparison_id, direction, gene_key_upper)
    data.table::setkey(contexts, model_rank, comparison_id, direction, gene_key_upper)
    joined <- links[contexts, nomatch = 0L, allow.cartesian = TRUE]
    if (!nrow(joined)) next
    tf_path <- file.path(extraction_dir, "tf_topic_assignments.csv")
    tf_scores <- if (file.exists(tf_path)) data.table::fread(tf_path, showProgress = FALSE) else data.table::data.table()
    if (nrow(tf_scores) && all(c("comparison_id", "direction", "tf", "topic_num", "theta") %in% names(tf_scores))) {
      if (!"primary_topic_num" %in% names(tf_scores)) {
        if ("primary_topic" %in% names(tf_scores)) {
          tf_scores[, primary_topic_num := suppressWarnings(as.integer(sub("^Topic", "", as.character(primary_topic))))]
        } else {
          tf_scores[, primary_topic_num := NA_integer_]
        }
      }
      tf_scores[, `:=`(
        model_rank = as.integer(hit$rank[[1L]]),
        comparison_id = as.character(comparison_id),
        direction = as.character(direction),
        tf_upper = toupper(as.character(tf)),
        topic_num = as.integer(topic_num),
        tf_topic_score = suppressWarnings(as.numeric(theta)),
        tf_primary_topic_num = suppressWarnings(as.integer(primary_topic_num))
      )]
      tf_scores <- tf_scores[, .(model_rank, comparison_id, direction, tf_upper, topic_num, tf_topic_score, tf_primary_topic_num)]
      joined <- merge(joined, tf_scores, by = c("model_rank", "comparison_id", "direction", "tf_upper", "topic_num"), all.x = TRUE, sort = FALSE)
    } else {
      joined[, `:=`(tf_topic_score = NA_real_, tf_primary_topic_num = NA_integer_)]
    }
    joined[, `:=`(
      tf_topic_score = data.table::fcoalesce(suppressWarnings(as.numeric(tf_topic_score)), 0),
      tf_primary_topic_num = suppressWarnings(as.integer(tf_primary_topic_num)),
      gene_topic_score = data.table::fcoalesce(suppressWarnings(as.numeric(gene_topic_score)), 0),
      gene_topic_pass = as.logical(gene_topic_pass)
    )]
    data.table::setorder(joined, subgrn_context_id, -edge_score_row, tf_upper, gene_key)
    joined <- joined[, head(.SD, max_links_per_context), by = subgrn_context_id]
    tf_gene <- joined[, {
      best_i <- which.max(edge_score_row)
      if (!length(best_i) || !is.finite(edge_score_row[best_i])) best_i <- 1L
      .(
        model_rank = model_rank[[1L]],
        comparison_id = comparison_id[[1L]],
        direction = direction[[1L]],
        topic_num = topic_num[[1L]],
        tf = tf[[1L]],
        tf_upper = tf_upper[[1L]],
        gene_key = gene_key[[1L]],
        target_gene = target_gene[[1L]],
        edge_score = sum(edge_score_row, na.rm = TRUE),
        abs_edge_score = max(edge_score_row, na.rm = TRUE),
        n_supporting_peaks = data.table::uniqueN(peak_id[!is.na(peak_id) & nzchar(peak_id)]),
        best_peak_id = peak_id[[best_i]],
        tf_topic_score = max(tf_topic_score, na.rm = TRUE),
        tf_primary_topic_num = tf_primary_topic_num[[1L]],
        gene_topic_score = max(gene_topic_score, na.rm = TRUE),
        gene_topic_pass = any(gene_topic_pass %in% TRUE)
      )
    }, by = .(subgrn_context_id, tf_upper, gene_key)]
    triplets <- joined[, .(
      model_rank = model_rank[[1L]],
      comparison_id = comparison_id[[1L]],
      direction = direction[[1L]],
      topic_num = topic_num[[1L]],
      tf = tf[[1L]],
      tf_upper = tf_upper[[1L]],
      peak_id = peak_id[[1L]],
      gene_key = gene_key[[1L]],
      target_gene = target_gene[[1L]],
      edge_score = max(edge_score_row, na.rm = TRUE),
      abs_edge_score = max(edge_score_row, na.rm = TRUE),
      tf_topic_score = max(tf_topic_score, na.rm = TRUE),
      tf_primary_topic_num = tf_primary_topic_num[[1L]],
      gene_topic_score = max(gene_topic_score, na.rm = TRUE),
      gene_topic_pass = any(gene_topic_pass %in% TRUE)
    ), by = .(subgrn_context_id, tf_upper, peak_id, gene_key)]
    all_tf_gene[[length(all_tf_gene) + 1L]] <- tf_gene
    all_triplets[[length(all_triplets) + 1L]] <- triplets
  }
  tf_gene_edges <- data.table::rbindlist(all_tf_gene, use.names = TRUE, fill = TRUE)
  tf_peak_gene_triplets <- data.table::rbindlist(all_triplets, use.names = TRUE, fill = TRUE)
  if (!length(names(tf_gene_edges))) {
    tf_gene_edges <- data.table::data.table(
      subgrn_context_id = character(),
      model_rank = integer(),
      comparison_id = character(),
      direction = character(),
      topic_num = integer(),
      tf = character(),
      tf_upper = character(),
      gene_key = character(),
      target_gene = character(),
      edge_score = numeric(),
      abs_edge_score = numeric(),
      n_supporting_peaks = integer(),
      best_peak_id = character(),
      tf_topic_score = numeric(),
      tf_primary_topic_num = integer(),
      gene_topic_score = numeric(),
      gene_topic_pass = logical()
    )
  }
  if (!length(names(tf_peak_gene_triplets))) {
    tf_peak_gene_triplets <- data.table::data.table(
      subgrn_context_id = character(),
      model_rank = integer(),
      comparison_id = character(),
      direction = character(),
      topic_num = integer(),
      tf = character(),
      tf_upper = character(),
      peak_id = character(),
      gene_key = character(),
      target_gene = character(),
      edge_score = numeric(),
      abs_edge_score = numeric(),
      tf_topic_score = numeric(),
      tf_primary_topic_num = integer(),
      gene_topic_score = numeric(),
      gene_topic_pass = logical()
    )
  }
  manifest <- data.table::copy(target_context)
  counts1 <- if (nrow(tf_gene_edges)) tf_gene_edges[, .(n_tf_gene_edges = .N), by = subgrn_context_id] else data.table::data.table(subgrn_context_id = character(), n_tf_gene_edges = integer())
  counts2 <- if (nrow(tf_peak_gene_triplets)) tf_peak_gene_triplets[, .(n_tf_peak_gene_triplets = .N), by = subgrn_context_id] else data.table::data.table(subgrn_context_id = character(), n_tf_peak_gene_triplets = integer())
  manifest <- merge(manifest, counts1, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  manifest <- merge(manifest, counts2, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  manifest[is.na(n_tf_gene_edges), n_tf_gene_edges := 0L]
  manifest[is.na(n_tf_peak_gene_triplets), n_tf_peak_gene_triplets := 0L]
  payload_file <- "paper_topic_subgrn_payload.js"
  manifest[, payload_file := payload_file]
  obj <- list(
    manifest = manifest,
    tf_gene_edges = tf_gene_edges,
    tf_peak_gene_triplets = tf_peak_gene_triplets
  )
  json <- jsonlite::toJSON(obj, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null")
  js <- paste0("window.CRAFTGRN_PAPER_PAYLOAD=", json, ";\n")
  writeLines(js, file.path(out_dir, payload_file), useBytes = TRUE)
  data.table::fwrite(manifest, file.path(out_dir, "paper_topic_subgrn_manifest.csv"))
  data.table::fwrite(tf_gene_edges, file.path(out_dir, "paper_topic_subgrn_tf_gene_edges.csv"))
  data.table::fwrite(tf_peak_gene_triplets, file.path(out_dir, "paper_topic_subgrn_tf_peak_gene_triplets.csv"))
  .m3fb_write_subgrn_index(
    out_dir,
    manifest,
    payload_file,
    top_models = unique(manifest[, .(model_rank, method_setup, model_label, k)]),
    theta_cutoff = theta_cutoff,
    gene_score_cutoff = gene_score_cutoff
  )
  out_dir
}

.m3fb_write_subgrn_index <- function(out_dir,
                                     manifest,
                                     payload_file,
                                     top_models = NULL,
                                     theta_cutoff = 0.3,
                                     gene_score_cutoff = 0.7) {
  manifest <- data.table::as.data.table(manifest)
  if (is.null(top_models)) top_models <- unique(manifest[, .(model_rank, method_setup, model_label, k)])
  top_models <- data.table::as.data.table(top_models)
  manifest_json <- .m3tb_json_for_html(manifest)
  models_json <- .m3tb_json_for_html(top_models)
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Paper topic SubGRN</title>",
    "<style>",
    "html,body{width:100%;height:100%;margin:0;background:#fbfbf8;color:#111;font-family:Arial,Helvetica,sans-serif;font-weight:700}body{overflow:hidden}.page{height:100%;display:grid;grid-template-rows:auto 1fr}.top{border-bottom:1px solid #d6d6d0;background:#fff;padding:10px 12px}.title{font-size:18px;font-weight:700;margin-bottom:8px}.controls{display:flex;flex-wrap:wrap;gap:7px 10px;align-items:center;font-size:11px}.controls label{display:flex;align-items:center;gap:4px}.controls select,.controls input,.controls textarea{font:700 11px Arial,Helvetica,sans-serif;border:1px solid #9ca3af;border-radius:3px;background:#fff;color:#111;padding:4px 5px}.controls input[type=number]{width:58px}.controls input[type=range]{width:92px}.controls textarea{width:155px;height:28px;resize:vertical}.controls button{font:700 11px Arial,Helvetica,sans-serif;border:1px solid #111;background:#111;color:#fff;border-radius:3px;padding:5px 8px;cursor:pointer}.body{position:relative;min-height:0}.stats{position:absolute;left:10px;right:10px;top:8px;z-index:2;background:rgba(255,255,255,.9);border:1px solid #d6d6d0;padding:6px 8px;font-size:11px;color:#334155}.canvas{position:absolute;inset:0;padding-top:40px}.node{cursor:grab}.node.selected{stroke:#111;stroke-width:4px}.edge{stroke:#94a3b8;stroke-opacity:.55;stroke-linecap:round}.label{font-family:Arial,Helvetica,sans-serif;font-weight:700;font-size:10px;paint-order:stroke;stroke:#fff;stroke-width:3px;stroke-linejoin:round;fill:#111;pointer-events:none}.tfLabel{fill:#fff;stroke:#475569}.tooltip{position:absolute;display:none;background:rgba(17,17,17,.92);color:#fff;font:700 11px Arial,Helvetica,sans-serif;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:360px;line-height:1.35;z-index:5}",
    "</style></head><body>",
    "<div class=\"page\"><div class=\"top\"><div class=\"title\">Paper-validated fibroblast topic SubGRN browser</div>",
    "<div class=\"controls\">",
    "<label>Method+K <select id=\"paperModelSelect\"></select></label>",
    "<label>Comparison <select id=\"paperComparisonSelect\"></select></label>",
    "<label>Topic <select id=\"paperTopicSelect\"></select></label>",
    "<label>Evidence <select id=\"paperEvidenceMode\"><option value=\"paper_all\" selected>Paper genes</option><option value=\"tf_or_gene\">TF or target passes</option><option value=\"tf_and_gene\">TF and target pass</option><option value=\"tf_only\">TF passes</option><option value=\"gene_only\">Target passes</option></select></label>",
    paste0("<label>TF theta <input id=\"paperTfCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"", theta_cutoff, "\"/></label>"),
    paste0("<label>Target score <input id=\"paperGeneCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"", gene_score_cutoff, "\"/></label>"),
    "<label>Primary TF topic <input id=\"paperPrimaryOnly\" type=\"checkbox\"/></label>",
    "<label>Network <select id=\"paperNetworkMode\"><option value=\"tf_gene\" selected>TF-gene</option><option value=\"tf_peak_gene\">TF-peak-gene</option></select></label>",
    "<label>Layout <select id=\"paperLayoutSelect\"><option value=\"force\" selected>Force</option><option value=\"radial\">Radial</option><option value=\"columns\">Columns</option><option value=\"bipartite\">Bipartite</option><option value=\"hierarchy\">Hierarchy</option><option value=\"concentric\">Concentric</option><option value=\"circle\">Circle</option><option value=\"grid\">Grid</option><option value=\"spiral\">Spiral</option><option value=\"clustered\">Clustered</option></select></label>",
    "<label>Top links <input id=\"paperTopLinks\" type=\"number\" min=\"0\" value=\"300\"/></label>",
    "<label>Spacing <input id=\"paperSpacing\" type=\"range\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/></label>",
    "<label>Labels <input id=\"paperShowLabels\" type=\"checkbox\" checked/></label>",
    "<label>TF list <textarea id=\"paperTfList\" placeholder=\"optional\"></textarea></label>",
    "<label>Target list <textarea id=\"paperGeneList\" placeholder=\"optional\"></textarea></label>",
    "<button id=\"paperResetButton\" type=\"button\">Reset layout</button>",
    "<button id=\"paperExportButton\" type=\"button\">Export SVG</button>",
    "</div></div><div class=\"body\"><div id=\"paperStats\" class=\"stats\"></div><div class=\"canvas\"><div id=\"paperTooltip\" class=\"tooltip\"></div>",
    "<svg id=\"paperSvg\" viewBox=\"0 0 1600 900\" width=\"100%\" height=\"100%\"><defs><marker id=\"paperArrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"6\" markerHeight=\"6\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#64748b\"/></marker></defs><g id=\"paperView\"><g id=\"paperEdges\"></g><g id=\"paperNodes\"></g><g id=\"paperLabels\"></g></g></svg>",
    "</div></div></div>",
    "<script src=\"", .m3tb_html_escape(payload_file), "\"></script>",
    "<script>",
    paste0("const PAPER_MANIFEST=", manifest_json, ";"),
    paste0("const PAPER_MODELS=", models_json, ";"),
    .m3fb_paper_subgrn_js(),
    "</script>",
    "</body></html>"
  )
  writeLines(html, file.path(out_dir, "index.html"), useBytes = TRUE)
  invisible(file.path(out_dir, "index.html"))
}

.m3fb_paper_subgrn_js <- function() {
  c(
    "const PAPER_W=1600,PAPER_H=900;const paperPayload=window.CRAFTGRN_PAPER_PAYLOAD||{manifest:[],tf_gene_edges:[],tf_peak_gene_triplets:[]};let paperState={nodes:[],edges:[],drag:null,selected:'',view:{x:0,y:0,k:1}};",
    "function pEl(id){return document.getElementById(id)}function pNum(x,d){const v=Number(x);return Number.isFinite(v)?v:d}function pEsc(x){return String(x==null?'':x).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]))}function pUniq(a){return Array.from(new Set(a))}function pOpt(sel,val,txt){const o=document.createElement('option');o.value=val;o.textContent=txt;sel.appendChild(o)}function pTokens(x){return new Set(String(x||'').split(/[;,\\n\\t ]+/).map(s=>s.trim().toUpperCase()).filter(Boolean))}function pLimit(){const v=Math.floor(pNum(pEl('paperTopLinks').value,300));return v===0?Infinity:Math.max(1,v)}function pSpacing(){return Math.max(.5,Math.min(2,pNum(pEl('paperSpacing').value,1)))}",
    "function initPaperControls(){const m=pEl('paperModelSelect'),c=pEl('paperComparisonSelect'),t=pEl('paperTopicSelect');m.replaceChildren();PAPER_MODELS.slice().sort((a,b)=>pNum(a.model_rank,0)-pNum(b.model_rank,0)).forEach(x=>pOpt(m,String(x.model_rank),'Rank '+x.model_rank+': '+x.method_setup+' K'+x.k));function fillComparisons(){c.replaceChildren();const r=pNum(m.value,1),vals=pUniq(PAPER_MANIFEST.filter(x=>pNum(x.model_rank,0)===r).map(x=>x.comparison_label)).sort();vals.forEach(x=>pOpt(c,x,x));fillTopics()}function fillTopics(){t.replaceChildren();const r=pNum(m.value,1),cmp=c.value,rows=PAPER_MANIFEST.filter(x=>pNum(x.model_rank,0)===r&&String(x.comparison_label)===cmp).sort((a,b)=>pNum(a.topic_num,0)-pNum(b.topic_num,0));rows.forEach(x=>pOpt(t,x.subgrn_context_id,'Topic '+x.topic_num+(x.selected_story_topic?' selected':'')));renderPaper()}m.addEventListener('change',fillComparisons);c.addEventListener('change',fillTopics);['paperTopicSelect','paperEvidenceMode','paperTfCutoff','paperGeneCutoff','paperPrimaryOnly','paperNetworkMode','paperLayoutSelect','paperTopLinks','paperSpacing','paperShowLabels','paperTfList','paperGeneList'].forEach(id=>{const e=pEl(id);if(e){e.addEventListener('change',renderPaper);e.addEventListener('input',renderPaper)}});pEl('paperResetButton').addEventListener('click',()=>{paperState.view={x:0,y:0,k:1};renderPaper()});pEl('paperExportButton').addEventListener('click',exportPaperSvg);fillComparisons()}",
    "function currentContext(){const id=pEl('paperTopicSelect').value;return PAPER_MANIFEST.find(x=>x.subgrn_context_id===id)||PAPER_MANIFEST[0]||null}function passEvidence(e,ctx){const tf=pNum(e.tf_topic_score,0),gene=pNum(e.gene_topic_score,0),tfCut=pNum(pEl('paperTfCutoff').value,.3),geneCut=pNum(pEl('paperGeneCutoff').value,.7),mode=pEl('paperEvidenceMode').value,primary=pEl('paperPrimaryOnly').checked;if(primary&&pNum(e.tf_primary_topic_num,-1)!==pNum(ctx.topic_num,-2))return false;if(mode==='tf_and_gene')return tf>=tfCut&&gene>=geneCut;if(mode==='tf_or_gene')return tf>=tfCut||gene>=geneCut;if(mode==='tf_only')return tf>=tfCut;if(mode==='gene_only')return gene>=geneCut;return true}",
    "function rowsForContext(ctx){const mode=pEl('paperNetworkMode').value,src=mode==='tf_peak_gene'?paperPayload.tf_peak_gene_triplets:paperPayload.tf_gene_edges;const tfKeep=pTokens(pEl('paperTfList').value),geneKeep=pTokens(pEl('paperGeneList').value);let rows=(src||[]).filter(e=>String(e.subgrn_context_id)===String(ctx.subgrn_context_id)&&passEvidence(e,ctx));if(tfKeep.size)rows=rows.filter(e=>tfKeep.has(String(e.tf||e.tf_upper).toUpperCase())||tfKeep.has(String(e.tf_upper||'').toUpperCase()));if(geneKeep.size)rows=rows.filter(e=>geneKeep.has(String(e.gene_key||'').toUpperCase())||geneKeep.has(String(e.target_gene||'').toUpperCase()));rows=rows.sort((a,b)=>pNum(b.abs_edge_score,0)-pNum(a.abs_edge_score,0)||String(a.tf||'').localeCompare(String(b.tf||''))||String(a.gene_key||'').localeCompare(String(b.gene_key||'')));return Number.isFinite(pLimit())?rows.slice(0,pLimit()):rows}",
    "function pNode(map,id,type,label){if(!map.has(id))map.set(id,{id:id,type:type,label:label||id,count:0,score:0,x:PAPER_W/2,y:PAPER_H/2});return map.get(id)}function buildGraph(rows){const mode=pEl('paperNetworkMode').value,nodes=new Map(),edges=[];rows.forEach(e=>{const tf='TF:'+String(e.tf||e.tf_upper),gene='GENE:'+String(e.gene_key||''),score=pNum(e.abs_edge_score,0);pNode(nodes,tf,'TF',String(e.tf||e.tf_upper));pNode(nodes,gene,'Gene',String(e.target_gene||e.gene_key));if(mode==='tf_peak_gene'){const peak='PEAK:'+String(e.peak_id||'peak');pNode(nodes,peak,'Peak',String(e.peak_id||'peak'));edges.push({from:tf,to:peak,score:score,title:(e.tf||'TF')+' -> '+(e.peak_id||'peak')+'\\nTF theta: '+pNum(e.tf_topic_score,0).toFixed(3)+'\\nTarget score: '+pNum(e.gene_topic_score,0).toFixed(3)});edges.push({from:peak,to:gene,score:score,title:(e.peak_id||'peak')+' -> '+(e.target_gene||e.gene_key)+'\\nedge score: '+score.toFixed(3)})}else{edges.push({from:tf,to:gene,score:score,title:(e.tf||'TF')+' -> '+(e.target_gene||e.gene_key)+'\\nTF theta: '+pNum(e.tf_topic_score,0).toFixed(3)+'\\nTarget score: '+pNum(e.gene_topic_score,0).toFixed(3)+'\\npeaks: '+(e.n_supporting_peaks||0)+'\\nbest peak: '+(e.best_peak_id||'NA')})}});edges.forEach(e=>{const a=nodes.get(e.from),b=nodes.get(e.to);if(a){a.count++;a.score+=e.score}if(b){b.count++;b.score+=e.score}});return{nodes:Array.from(nodes.values()),edges:edges}}",
    "function byType(ns,type){return ns.filter(n=>n.type===type).sort((a,b)=>b.count-a.count||String(a.label).localeCompare(String(b.label)))}function setCol(items,x,top,bottom){const gap=Math.min(80*pSpacing(),(bottom-top)/Math.max(items.length-1,1));items.forEach((n,i)=>{n.x=x;n.y=(top+bottom)/2-gap*(items.length-1)/2+i*gap})}function ring(items,r,a0){items.forEach((n,i)=>{const a=a0+2*Math.PI*i/Math.max(1,items.length);n.x=PAPER_W/2+Math.cos(a)*r*pSpacing();n.y=PAPER_H/2+Math.sin(a)*r*pSpacing()})}function layoutGraph(g){const ns=g.nodes,es=g.edges,m=pEl('paperLayoutSelect').value;if(m==='columns'||m==='bipartite'||m==='hierarchy'){if(pEl('paperNetworkMode').value==='tf_peak_gene'){setCol(byType(ns,'TF'),260,80,830);setCol(byType(ns,'Peak'),800,80,830);setCol(byType(ns,'Gene'),1340,80,830)}else{setCol(byType(ns,'TF'),420,80,830);setCol(byType(ns,'Gene'),1180,80,830)}}else if(m==='radial'||m==='concentric'||m==='clustered'){ring(byType(ns,'TF'),220,-Math.PI/2);ring(byType(ns,'Peak'),320,0);ring(byType(ns,'Gene'),410,-Math.PI/2)}else if(m==='circle'){ring(ns.slice().sort((a,b)=>String(a.label).localeCompare(String(b.label))),360,-Math.PI/2)}else if(m==='grid'){const cols=Math.ceil(Math.sqrt(Math.max(1,ns.length))),gx=1180/Math.max(1,cols-1),rows=Math.ceil(ns.length/cols),gy=680/Math.max(1,rows-1);ns.slice().sort((a,b)=>a.type.localeCompare(b.type)||String(a.label).localeCompare(String(b.label))).forEach((n,i)=>{n.x=210+(i%cols)*gx;n.y=120+Math.floor(i/cols)*gy})}else if(m==='spiral'){ns.slice().sort((a,b)=>b.count-a.count).forEach((n,i)=>{const a=i*.72,r=(28+i*8)*pSpacing();n.x=PAPER_W/2+Math.cos(a)*r;n.y=PAPER_H/2+Math.sin(a)*r})}else{ring(byType(ns,'TF'),240,-Math.PI/2);ring(byType(ns,'Peak'),330,0);ring(byType(ns,'Gene'),420,-Math.PI/2);forceGraph(ns,es)}}function forceGraph(ns,es){const by=new Map(ns.map(n=>[n.id,n]));for(let iter=0;iter<120;iter++){ns.forEach(n=>{n.vx=n.vx||0;n.vy=n.vy||0});for(let i=0;i<ns.length;i++)for(let j=i+1;j<ns.length;j++){let a=ns[i],b=ns[j],dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(dx*dx+dy*dy,64),d=Math.sqrt(d2),f=Math.min(3,2400*pSpacing()/d2);dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f}es.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;let dx=b.x-a.x,dy=b.y-a.y,d=Math.max(Math.sqrt(dx*dx+dy*dy),1),f=(d-130*pSpacing())*.011;dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f});ns.forEach(n=>{n.vx+=(PAPER_W/2-n.x)*.0015;n.vy+=(PAPER_H/2-n.y)*.0015;n.x=Math.max(35,Math.min(PAPER_W-35,n.x+n.vx));n.y=Math.max(55,Math.min(PAPER_H-35,n.y+n.vy));n.vx*=.74;n.vy*=.74})}ns.forEach(n=>{delete n.vx;delete n.vy})}",
    "function el(name,attrs){const x=document.createElementNS('http://www.w3.org/2000/svg',name);Object.entries(attrs||{}).forEach(([k,v])=>x.setAttribute(k,v));return x}function nodeR(n){if(n.type==='TF')return Math.max(10,Math.min(30,10+Math.sqrt(n.count)*5));if(n.type==='Peak')return 5;return Math.max(6,Math.min(14,6+Math.sqrt(n.count)*2))}function lineEnd(a,b){const dx=b.x-a.x,dy=b.y-a.y,len=Math.max(Math.sqrt(dx*dx+dy*dy),1),r=b.r+3;return{x:b.x-dx/len*r,y:b.y-dy/len*r}}function applyView(){pEl('paperView').setAttribute('transform','translate('+paperState.view.x+' '+paperState.view.y+') scale('+paperState.view.k+')')}",
    "function renderPaper(){const ctx=currentContext(),stats=pEl('paperStats');if(!ctx){stats.textContent='No paper SubGRN context is available.';return}const rows=rowsForContext(ctx),g=buildGraph(rows);paperState.nodes=g.nodes;paperState.edges=g.edges;g.nodes.forEach(n=>{n.r=nodeR(n)});layoutGraph(g);drawPaper(g,ctx,rows)}function drawPaper(g,ctx,rows){const edgeLayer=pEl('paperEdges'),nodeLayer=pEl('paperNodes'),labelLayer=pEl('paperLabels'),stats=pEl('paperStats'),showLabels=pEl('paperShowLabels').checked,by=new Map(g.nodes.map(n=>[n.id,n]));edgeLayer.replaceChildren();nodeLayer.replaceChildren();labelLayer.replaceChildren();g.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const end=lineEnd(a,b),w=1+Math.min(3,Math.sqrt(pNum(e.score,0)));const line=el('line',{class:'edge',x1:a.x,y1:a.y,x2:end.x,y2:end.y,'stroke-width':w,'marker-end':'url(#paperArrow)'});line.dataset.title=e.title;edgeLayer.appendChild(line)});g.nodes.forEach(n=>{let shape;if(n.type==='TF')shape=el('rect',{class:'node',x:n.x-n.r*1.4,y:n.y-n.r*.75,width:n.r*2.8,height:n.r*1.5,rx:3,fill:'#2563eb'});else shape=el('circle',{class:'node',cx:n.x,cy:n.y,r:n.r,fill:n.type==='Peak'?'#f59e0b':'#991b1b'});shape.dataset.id=n.id;shape.dataset.title=n.label+'\\n'+n.type+'\\nedges: '+n.count;nodeLayer.appendChild(shape);if(showLabels){const tx=n.type==='TF'?n.x:n.x+n.r+6,ty=n.type==='TF'?n.y+4:n.y+3,anch=n.type==='TF'?'middle':'start';const text=el('text',{class:'label '+(n.type==='TF'?'tfLabel':''),x:tx,y:ty,'text-anchor':anch});text.textContent=n.label;labelLayer.appendChild(text)}});stats.textContent='Rank '+ctx.model_rank+' '+ctx.method_setup+' K'+ctx.k+' | '+ctx.comparison_label+' | Topic '+ctx.topic_num+' | nodes '+g.nodes.length+' | edges '+g.edges.length+' | rows '+rows.length+' | targets '+(ctx.overlap_genes||'');applyView()}",
    "function redrawPaper(){const by=new Map(paperState.nodes.map(n=>[n.id,n]));pEl('paperEdges').querySelectorAll('line').forEach((line,i)=>{const e=paperState.edges[i],a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const end=lineEnd(a,b);line.setAttribute('x1',a.x);line.setAttribute('y1',a.y);line.setAttribute('x2',end.x);line.setAttribute('y2',end.y)});pEl('paperNodes').querySelectorAll('.node').forEach(shape=>{const n=by.get(shape.dataset.id);if(!n)return;if(shape.tagName.toLowerCase()==='rect'){shape.setAttribute('x',n.x-n.r*1.4);shape.setAttribute('y',n.y-n.r*.75)}else{shape.setAttribute('cx',n.x);shape.setAttribute('cy',n.y)}});let j=0;pEl('paperLabels').querySelectorAll('text').forEach(text=>{const n=paperState.nodes[j++];if(!n)return;text.setAttribute('x',n.type==='TF'?n.x:n.x+n.r+6);text.setAttribute('y',n.type==='TF'?n.y+4:n.y+3)})}",
    "function exportPaperSvg(){const svg=pEl('paperSvg'),clone=svg.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('width',PAPER_W);clone.setAttribute('height',PAPER_H);const blob=new Blob([new XMLSerializer().serializeToString(clone)],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='paper_topic_subgrn.svg';document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(a.href)}",
    "{const svg=pEl('paperSvg'),nodes=pEl('paperNodes'),tip=pEl('paperTooltip');function svgPoint(ev){const p=svg.createSVGPoint();p.x=ev.clientX;p.y=ev.clientY;return p.matrixTransform(svg.getScreenCTM().inverse())}function world(ev){const p=svgPoint(ev);return{x:(p.x-paperState.view.x)/paperState.view.k,y:(p.y-paperState.view.y)/paperState.view.k}}svg.addEventListener('wheel',ev=>{ev.preventDefault();const p=svgPoint(ev),old=paperState.view.k,next=Math.max(.25,Math.min(6,old*(ev.deltaY<0?1.15:1/1.15)));paperState.view.x=p.x-(p.x-paperState.view.x)*next/old;paperState.view.y=p.y-(p.y-paperState.view.y)*next/old;paperState.view.k=next;applyView()},{passive:false});nodes.addEventListener('mousedown',ev=>{if(!(ev.target.classList&&ev.target.classList.contains('node')))return;ev.preventDefault();paperState.drag=paperState.nodes.find(n=>n.id===ev.target.dataset.id)||null});svg.addEventListener('mousemove',ev=>{if(paperState.drag){const p=world(ev);paperState.drag.x=Math.max(30,Math.min(PAPER_W-30,p.x));paperState.drag.y=Math.max(40,Math.min(PAPER_H-30,p.y));redrawPaper()}if(tip){tip.style.left=(ev.offsetX+12)+'px';tip.style.top=(ev.offsetY+12)+'px'}});svg.addEventListener('mouseup',()=>{paperState.drag=null});svg.addEventListener('mouseleave',()=>{paperState.drag=null;if(tip)tip.style.display='none'});svg.addEventListener('mouseover',ev=>{const title=ev.target.dataset?ev.target.dataset.title:'';if(title&&tip){tip.innerHTML=pEsc(title).replace(/\\n/g,'<br/>');tip.style.display='block'}});svg.addEventListener('mouseout',ev=>{if(ev.target.dataset&&ev.target.dataset.title&&tip)tip.style.display='none'})}",
    "initPaperControls();"
  )
}
