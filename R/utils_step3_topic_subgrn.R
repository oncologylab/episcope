.m3tb_subgrn_direction_label <- function(x) {
  x <- as.character(x)
  out <- x
  out[tolower(x) %in% c("up", "target-up", "target_up")] <- "Target-Up"
  out[tolower(x) %in% c("down", "target-down", "target_down")] <- "Target-Down"
  out
}

.m3tb_subgrn_direction_group <- function(x) {
  x <- .m3tb_subgrn_direction_label(x)
  out <- x
  out[out == "Target-Up"] <- "Up"
  out[out == "Target-Down"] <- "Down"
  out
}

.m3tb_subgrn_split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  parts <- strsplit(x, ";", fixed = TRUE)
  lapply(parts, function(v) {
    v <- trimws(v)
    unique(v[nzchar(v)])
  })
}

.m3tb_subgrn_context_table <- function(pathways, require_comparison = FALSE) {
  dt <- data.table::copy(data.table::as.data.table(pathways))
  if (!nrow(dt)) return(dt)
  if (!"pathway_key" %in% names(dt) && "pathway_norm_key" %in% names(dt)) {
    data.table::setnames(dt, "pathway_norm_key", "pathway_key")
  }
  if (!"topic" %in% names(dt) && "topic_num" %in% names(dt)) data.table::setnames(dt, "topic_num", "topic")
  if (!"topic_num" %in% names(dt) && "topic" %in% names(dt)) dt[, topic_num := topic]
  if (!"comparison_id" %in% names(dt)) dt[, comparison_id := NA_character_]
  if (!"direction_group" %in% names(dt)) dt[, direction_group := NA_character_]
  if ("comparison_label" %in% names(dt)) {
    missing_cmp <- is.na(dt$comparison_id) | !nzchar(as.character(dt$comparison_id))
    if (any(missing_cmp)) {
      label_chr <- as.character(dt$comparison_label)
      dt[missing_cmp, comparison_id := sub("::.*$", "", label_chr[missing_cmp])]
      missing_dir <- is.na(dt$direction_group) | !nzchar(as.character(dt$direction_group))
      dt[missing_dir, direction_group := ifelse(
        grepl("::", label_chr[missing_dir], fixed = TRUE),
        sub("^.*::", "", label_chr[missing_dir]),
        "All"
      )]
    }
  }
  dt[, `:=`(
    comparison_id = as.character(comparison_id),
    direction_group = .m3tb_subgrn_direction_group(direction_group),
    topic = as.integer(topic),
    topic_num = as.integer(topic_num),
    pathway = as.character(pathway),
    pathway_key = as.character(pathway_key)
  )]
  keep <- !is.na(dt$topic) & is.finite(dt$topic) & !is.na(dt$pathway_key) & nzchar(dt$pathway_key)
  if (isTRUE(require_comparison)) {
    keep <- keep &
      !is.na(dt$comparison_id) & nzchar(dt$comparison_id) &
      !is.na(dt$direction_group) & nzchar(dt$direction_group)
  }
  dt[keep]
}

.m3tb_select_pathway_subgrn_contexts <- function(all_pathways,
                                                 condition_pathways = NULL,
                                                 topic_pathways = NULL) {
  all_dt <- .m3tb_subgrn_context_table(all_pathways, require_comparison = TRUE)
  if (!nrow(all_dt)) return(all_dt)
  if ("overlap_genes" %in% names(all_dt)) {
    all_dt <- all_dt[!is.na(overlap_genes) & nzchar(overlap_genes)]
  }
  if (!nrow(all_dt)) return(all_dt)
  if (!"padj" %in% names(all_dt)) all_dt[, padj := NA_real_]
  if (!"overlap_hits" %in% names(all_dt)) {
    all_dt[, overlap_hits := lengths(.m3tb_subgrn_split_genes(if ("overlap_genes" %in% names(all_dt)) overlap_genes else ""))]
  }
  if (!"combined_score" %in% names(all_dt)) all_dt[, combined_score := NA_real_]
  all_dt[, `:=`(
    padj = suppressWarnings(as.numeric(padj)),
    overlap_hits = suppressWarnings(as.integer(overlap_hits)),
    combined_score = suppressWarnings(as.numeric(combined_score))
  )]
  data.table::setorder(all_dt, padj, -overlap_hits, -combined_score, comparison_id, direction_group, topic, pathway)
  all_dt <- unique(all_dt, by = c("comparison_id", "direction_group", "topic", "pathway_key"))

  selected <- list()
  if (!is.null(condition_pathways)) {
    cond <- .m3tb_subgrn_context_table(condition_pathways, require_comparison = TRUE)
    if (nrow(cond)) {
      cond <- unique(cond[, .(comparison_id, direction_group, topic, pathway_key)])
      selected[[length(selected) + 1L]] <- merge(
        all_dt,
        cond,
        by = c("comparison_id", "direction_group", "topic", "pathway_key"),
        all = FALSE,
        sort = FALSE
      )
    }
  }
  if (!is.null(topic_pathways)) {
    topic <- .m3tb_subgrn_context_table(topic_pathways, require_comparison = FALSE)
    if (nrow(topic)) {
      topic <- unique(topic[, .(topic, pathway_key)])
      selected[[length(selected) + 1L]] <- merge(
        all_dt,
        topic,
        by = c("topic", "pathway_key"),
        all = FALSE,
        sort = FALSE,
        allow.cartesian = TRUE
      )
    }
  }
  selected <- selected[!vapply(selected, is.null, logical(1L))]
  if (!length(selected)) return(all_dt)
  out <- data.table::rbindlist(selected, use.names = TRUE, fill = TRUE)
  out <- unique(out, by = c("comparison_id", "direction_group", "topic", "pathway_key"))
  data.table::setorder(out, padj, -overlap_hits, -combined_score, comparison_id, direction_group, topic, pathway)
  out[]
}

.m3tb_build_pathway_subgrn_payload <- function(pathways,
                                               edges_docs,
                                               tf_membership = NULL,
                                               max_tf_gene_edges_per_context = 300L,
                                               max_tf_peak_gene_triplets_per_context = 600L) {
  empty <- list(
    manifest = data.table::data.table(),
    tf_gene_edges = data.table::data.table(),
    tf_peak_gene_triplets = data.table::data.table()
  )
  pathways <- data.table::as.data.table(pathways)
  edges_docs <- data.table::as.data.table(edges_docs)
  if (!nrow(pathways) || !nrow(edges_docs)) return(empty)
  if (!"comparison_id" %in% names(edges_docs) && "condition_id" %in% names(edges_docs)) {
    edges_docs[, comparison_id := as.character(condition_id)]
  }
  if (!"direction" %in% names(edges_docs)) edges_docs[, direction := "All"]
  if (!"direction_group" %in% names(edges_docs)) edges_docs[, direction_group := "All"]
  if (!"pathway_key" %in% names(pathways) && "pathway_norm_key" %in% names(pathways)) {
    data.table::setnames(pathways, "pathway_norm_key", "pathway_key")
  }
  required_pathway <- c("comparison_id", "direction_group", "topic", "pathway", "pathway_key")
  if (!all(required_pathway %in% names(pathways))) return(empty)
  if (!"overlap_genes" %in% names(pathways)) pathways[, overlap_genes := ""]
  if (!"padj" %in% names(pathways)) pathways[, padj := NA_real_]
  pathways[, `:=`(
    comparison_id = as.character(comparison_id),
    direction = .m3tb_subgrn_direction_label(direction_group),
    topic_num = as.integer(topic),
    pathway = as.character(pathway),
    pathway_key = as.character(pathway_key),
    padj = suppressWarnings(as.numeric(padj))
  )]
  pathways <- pathways[!is.na(comparison_id) & nzchar(comparison_id) &
    !is.na(direction) & nzchar(direction) &
    is.finite(topic_num) &
    !is.na(pathway_key) & nzchar(pathway_key)]
  if (!nrow(pathways)) return(empty)
  data.table::setorder(pathways, comparison_id, direction, topic_num, pathway_key, padj)
  pathways[, subgrn_context_id := sprintf("subgrn_%06d", seq_len(.N))]
  pathways[, comparison_label := ifelse(direction == "All", comparison_id, paste(comparison_id, direction, sep = "::"))]
  gene_lists <- .m3tb_subgrn_split_genes(pathways$overlap_genes)
  context_genes <- data.table::rbindlist(lapply(seq_len(nrow(pathways)), function(i) {
    genes <- gene_lists[[i]]
    if (!length(genes)) return(NULL)
    data.table::data.table(
      subgrn_context_id = pathways$subgrn_context_id[[i]],
      comparison_id = pathways$comparison_id[[i]],
      direction = pathways$direction[[i]],
      topic_num = pathways$topic_num[[i]],
      pathway_key = pathways$pathway_key[[i]],
      gene_key = genes
    )
  }), use.names = TRUE, fill = TRUE)
  manifest <- pathways[, .(
    subgrn_context_id,
    comparison_id,
    direction,
    comparison_label,
    topic_num,
    pathway_key,
    pathway,
    padj,
    n_overlap_genes = lengths(gene_lists)
  )]
  if (!nrow(context_genes)) {
    manifest[, `:=`(n_tf_gene_edges = 0L, n_tf_peak_gene_triplets = 0L, n_topic_tfs = 0L)]
    return(list(manifest = manifest, tf_gene_edges = empty$tf_gene_edges, tf_peak_gene_triplets = empty$tf_peak_gene_triplets))
  }
  needed_edge_cols <- c("comparison_id", "direction", "direction_group", "tf", "gene_key", "peak_id")
  if (!all(c("comparison_id", "tf", "gene_key") %in% names(edges_docs))) return(empty)
  for (col in setdiff(needed_edge_cols, names(edges_docs))) edges_docs[, (col) := NA_character_]
  if (all(is.na(edges_docs$direction)) || !any(nzchar(as.character(edges_docs$direction)))) {
    edges_docs[, direction := .m3tb_subgrn_direction_label(direction_group)]
  } else {
    edges_docs[, direction := .m3tb_subgrn_direction_label(direction)]
  }
  score_cols <- c("delta_link_score", "delta_fp_score", "delta_gene_expr", "log2FC_gene_expr",
                  "delta_fp", "delta_gene", "log2fc_gene", "log2fc_fp",
                  "fc_mag_gene", "fp_score_condition")
  for (col in setdiff(score_cols, names(edges_docs))) edges_docs[, (col) := NA_real_]
  for (col in score_cols) edges_docs[, (col) := suppressWarnings(as.numeric(get(col)))]
  keep_edge_cols <- unique(c(
    "comparison_id", "direction", "direction_group", "tf", "gene_key", "peak_id",
    score_cols, "distance_to_tss", "log2fc_gene", "log2fc_tf", "delta_fp", "delta_gene"
  ))
  edges_docs <- edges_docs[, intersect(keep_edge_cols, names(edges_docs)), with = FALSE]
  edges_docs[, `:=`(
    comparison_id = as.character(comparison_id),
    tf = as.character(tf),
    tf_upper = toupper(as.character(tf)),
    gene_key = as.character(gene_key),
    peak_id = as.character(peak_id),
    edge_score_row = abs(data.table::fcoalesce(
      delta_link_score, delta_fp_score, delta_gene_expr, log2FC_gene_expr,
      delta_fp, delta_gene, log2fc_gene, log2fc_fp, fc_mag_gene,
      fp_score_condition, 0
    ))
  )]
  edges_docs <- edges_docs[!is.na(comparison_id) & nzchar(comparison_id) &
    !is.na(direction) & nzchar(direction) &
    !is.na(tf) & nzchar(tf) &
    !is.na(gene_key) & nzchar(gene_key)]
  edges_docs <- edges_docs[
    comparison_id %chin% unique(context_genes$comparison_id) &
      direction %chin% unique(context_genes$direction) &
      gene_key %chin% unique(context_genes$gene_key)
  ]
  if (!nrow(edges_docs)) return(list(manifest = manifest, tf_gene_edges = empty$tf_gene_edges, tf_peak_gene_triplets = empty$tf_peak_gene_triplets))
  data.table::setkey(edges_docs, comparison_id, direction, gene_key)
  data.table::setkey(context_genes, comparison_id, direction, gene_key)
  joined <- edges_docs[context_genes, nomatch = 0L, allow.cartesian = TRUE]
  if (!nrow(joined)) {
    manifest[, `:=`(n_tf_gene_edges = 0L, n_tf_peak_gene_triplets = 0L, n_topic_tfs = 0L)]
    return(list(manifest = manifest, tf_gene_edges = empty$tf_gene_edges, tf_peak_gene_triplets = empty$tf_peak_gene_triplets))
  }
  topic_tf <- data.table::data.table()
  if (!is.null(tf_membership)) {
    topic_tf <- data.table::as.data.table(tf_membership)
    if (!"comparison_id" %in% names(topic_tf) && "condition_id" %in% names(topic_tf)) {
      topic_tf[, comparison_id := as.character(condition_id)]
    }
    if (!"direction" %in% names(topic_tf)) topic_tf[, direction := "All"]
    if (nrow(topic_tf) && all(c("comparison_id", "direction", "tf", "topic_num") %in% names(topic_tf))) {
      if (!"membership_pass" %in% names(topic_tf)) topic_tf[, membership_pass := TRUE]
      topic_tf[, `:=`(
        comparison_id = as.character(comparison_id),
        direction = .m3tb_subgrn_direction_label(direction),
        tf_upper = toupper(as.character(tf)),
        topic_num = as.integer(topic_num),
        membership_pass = as.logical(membership_pass)
      )]
      topic_tf[is.na(direction) | !nzchar(direction), direction := "All"]
      topic_tf <- unique(topic_tf[isTRUE(membership_pass) | membership_pass == TRUE,
        .(comparison_id, direction, tf_upper, topic_num, topic_tf = TRUE)
      ])
    } else {
      topic_tf <- data.table::data.table()
    }
  }
  if (nrow(topic_tf)) {
    joined <- merge(
      joined,
      topic_tf,
      by = c("comparison_id", "direction", "tf_upper", "topic_num"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    joined[, topic_tf := NA]
  }
  joined[, topic_tf := isTRUE(topic_tf) | topic_tf == TRUE]
  max_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  min_abs_finite <- function(x) {
    x <- suppressWarnings(abs(as.numeric(x)))
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  for (col in c("distance_to_tss", "log2fc_gene", "log2fc_tf", "delta_fp", "delta_gene")) {
    if (!col %in% names(joined)) joined[, (col) := NA_real_]
    joined[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  tf_gene <- joined[, {
    best_i <- which.max(edge_score_row)
    if (!length(best_i) || !is.finite(edge_score_row[best_i])) best_i <- 1L
    data.table::data.table(
      tf = tf[[1L]],
      tf_upper = tf_upper[[1L]],
      gene_key = gene_key[[1L]],
      edge_score = sum(edge_score_row, na.rm = TRUE),
      abs_edge_score = max_finite(edge_score_row),
      n_supporting_peaks = data.table::uniqueN(peak_id[!is.na(peak_id) & nzchar(peak_id)]),
      best_peak_id = peak_id[[best_i]],
      best_distance_to_tss = min_abs_finite(distance_to_tss),
      topic_tf = any(topic_tf, na.rm = TRUE)
    )
  }, by = .(subgrn_context_id, tf_upper, gene_key)]
  tf_gene[is.na(n_supporting_peaks), n_supporting_peaks := 0L]
  tf_gene[, `:=`(
    tf_target_count = data.table::uniqueN(gene_key),
    tf_edge_score_sum = sum(abs_edge_score, na.rm = TRUE)
  ), by = .(subgrn_context_id, tf_upper)]
  data.table::setorder(tf_gene, subgrn_context_id, -tf_target_count, -tf_edge_score_sum, -abs_edge_score, tf_upper, gene_key)
  tf_gene[, tf_rank := match(tf_upper, unique(tf_upper)), by = subgrn_context_id]
  tf_gene[, edge_rank := seq_len(.N), by = subgrn_context_id]
  tf_gene_full_counts <- tf_gene[, .(n_tf_gene_edges = .N, n_topic_tfs = data.table::uniqueN(tf_upper[topic_tf])), by = subgrn_context_id]
  tf_gene <- tf_gene[edge_rank <= max(1L, as.integer(max_tf_gene_edges_per_context)[[1L]])]
  triplets <- joined[, .(
    tf = tf[[1L]],
    tf_upper = tf_upper[[1L]],
    peak_id = peak_id[[1L]],
    gene_key = gene_key[[1L]],
    edge_score = max_finite(edge_score_row),
    abs_edge_score = max_finite(edge_score_row),
    topic_tf = any(topic_tf, na.rm = TRUE),
    distance_to_tss = min_abs_finite(distance_to_tss)
  ), by = .(subgrn_context_id, tf_upper, peak_id, gene_key)]
  data.table::setorder(triplets, subgrn_context_id, -abs_edge_score, tf_upper, peak_id, gene_key)
  triplet_full_counts <- triplets[, .(n_tf_peak_gene_triplets = .N), by = subgrn_context_id]
  triplets[, triplet_rank := seq_len(.N), by = subgrn_context_id]
  triplets <- triplets[triplet_rank <= max(1L, as.integer(max_tf_peak_gene_triplets_per_context)[[1L]])]
  manifest <- merge(manifest, tf_gene_full_counts, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  manifest <- merge(manifest, triplet_full_counts, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  for (col in c("n_tf_gene_edges", "n_topic_tfs", "n_tf_peak_gene_triplets")) {
    manifest[is.na(get(col)), (col) := 0L]
    manifest[, (col) := as.integer(get(col))]
  }
  list(
    manifest = manifest,
    tf_gene_edges = tf_gene,
    tf_peak_gene_triplets = triplets
  )
}

.m3tb_build_pathway_subgrn_compact_payload <- function(pathways,
                                                       edges_docs,
                                                       tf_membership = NULL) {
  empty <- list(
    manifest = data.table::data.table(),
    context_genes = data.table::data.table(),
    tf_gene_edges = data.table::data.table(),
    tf_peak_gene_triplets = data.table::data.table()
  )
  pathways <- data.table::as.data.table(pathways)
  edges_docs <- data.table::as.data.table(edges_docs)
  if (!nrow(pathways) || !nrow(edges_docs)) return(empty)
  if (!"comparison_id" %in% names(edges_docs) && "condition_id" %in% names(edges_docs)) {
    edges_docs[, comparison_id := as.character(condition_id)]
  }
  if (!"direction" %in% names(edges_docs)) edges_docs[, direction := "All"]
  if (!"direction_group" %in% names(edges_docs)) edges_docs[, direction_group := "All"]
  if (!"pathway_key" %in% names(pathways) && "pathway_norm_key" %in% names(pathways)) {
    data.table::setnames(pathways, "pathway_norm_key", "pathway_key")
  }
  required_pathway <- c("comparison_id", "direction_group", "topic", "pathway", "pathway_key")
  if (!all(required_pathway %in% names(pathways))) return(empty)
  if (!"overlap_genes" %in% names(pathways)) pathways[, overlap_genes := ""]
  if (!"padj" %in% names(pathways)) pathways[, padj := NA_real_]
  pathways[, `:=`(
    comparison_id = as.character(comparison_id),
    direction = .m3tb_subgrn_direction_label(direction_group),
    topic_num = as.integer(topic),
    pathway = as.character(pathway),
    pathway_key = as.character(pathway_key),
    overlap_genes = as.character(overlap_genes),
    padj = suppressWarnings(as.numeric(padj))
  )]
  pathways <- pathways[!is.na(comparison_id) & nzchar(comparison_id) &
    !is.na(direction) & nzchar(direction) &
    is.finite(topic_num) &
    !is.na(pathway_key) & nzchar(pathway_key)]
  if (!nrow(pathways)) return(empty)
  data.table::setorder(pathways, comparison_id, direction, topic_num, pathway_key, padj)
  pathways[, subgrn_context_id := sprintf("subgrn_%06d", seq_len(.N))]
  pathways[, comparison_label := ifelse(direction == "All", comparison_id, paste(comparison_id, direction, sep = "::"))]
  gene_lists <- .m3tb_subgrn_split_genes(pathways$overlap_genes)
  context_genes <- data.table::rbindlist(lapply(seq_len(nrow(pathways)), function(i) {
    genes <- gene_lists[[i]]
    if (!length(genes)) return(NULL)
    data.table::data.table(
      subgrn_context_id = pathways$subgrn_context_id[[i]],
      comparison_id = pathways$comparison_id[[i]],
      direction = pathways$direction[[i]],
      topic_num = pathways$topic_num[[i]],
      gene_key = genes
    )
  }), use.names = TRUE, fill = TRUE)
  manifest <- pathways[, .(
    subgrn_context_id,
    comparison_id,
    direction,
    comparison_label,
    topic_num,
    pathway_key,
    pathway,
    padj,
    overlap_genes,
    n_overlap_genes = lengths(gene_lists)
  )]
  if (!nrow(context_genes)) {
    manifest[, `:=`(n_tf_gene_edges = 0L, n_tf_peak_gene_triplets = 0L, n_topic_tfs = 0L)]
    return(list(manifest = manifest, context_genes = context_genes, tf_gene_edges = empty$tf_gene_edges, tf_peak_gene_triplets = empty$tf_peak_gene_triplets))
  }

  needed_edge_cols <- c("comparison_id", "direction", "direction_group", "tf", "gene_key", "peak_id")
  if (!all(c("comparison_id", "tf", "gene_key") %in% names(edges_docs))) return(empty)
  for (col in setdiff(needed_edge_cols, names(edges_docs))) edges_docs[, (col) := NA_character_]
  if (all(is.na(edges_docs$direction)) || !any(nzchar(as.character(edges_docs$direction)))) {
    edges_docs[, direction := .m3tb_subgrn_direction_label(direction_group)]
  } else {
    edges_docs[, direction := .m3tb_subgrn_direction_label(direction)]
  }
  score_cols <- c("delta_link_score", "delta_fp_score", "delta_gene_expr", "log2FC_gene_expr",
                  "delta_fp", "delta_gene", "log2fc_gene", "log2fc_fp",
                  "fc_mag_gene", "fp_score_condition")
  for (col in setdiff(score_cols, names(edges_docs))) edges_docs[, (col) := NA_real_]
  for (col in score_cols) edges_docs[, (col) := suppressWarnings(as.numeric(get(col)))]
  keep_edge_cols <- unique(c(
    "comparison_id", "direction", "direction_group", "tf", "gene_key", "peak_id",
    score_cols, "distance_to_tss", "log2fc_gene", "log2fc_tf", "delta_fp", "delta_gene"
  ))
  edges_docs <- edges_docs[, intersect(keep_edge_cols, names(edges_docs)), with = FALSE]
  edges_docs[, `:=`(
    comparison_id = as.character(comparison_id),
    tf = as.character(tf),
    tf_upper = toupper(as.character(tf)),
    gene_key = as.character(gene_key),
    peak_id = as.character(peak_id),
    edge_score_row = abs(data.table::fcoalesce(
      delta_link_score, delta_fp_score, delta_gene_expr, log2FC_gene_expr,
      delta_fp, delta_gene, log2fc_gene, log2fc_fp, fc_mag_gene,
      fp_score_condition, 0
    ))
  )]
  edges_docs <- edges_docs[!is.na(comparison_id) & nzchar(comparison_id) &
    !is.na(direction) & nzchar(direction) &
    !is.na(tf) & nzchar(tf) &
    !is.na(gene_key) & nzchar(gene_key)]
  needed_genes <- unique(context_genes[, .(comparison_id, direction, gene_key)])
  data.table::setkey(edges_docs, comparison_id, direction, gene_key)
  data.table::setkey(needed_genes, comparison_id, direction, gene_key)
  joined <- edges_docs[needed_genes, nomatch = 0L, allow.cartesian = TRUE]
  if (!nrow(joined)) {
    manifest[, `:=`(n_tf_gene_edges = 0L, n_tf_peak_gene_triplets = 0L, n_topic_tfs = 0L)]
    return(list(manifest = manifest, context_genes = context_genes, tf_gene_edges = empty$tf_gene_edges, tf_peak_gene_triplets = empty$tf_peak_gene_triplets))
  }

  topic_tf <- data.table::data.table()
  if (!is.null(tf_membership)) {
    topic_tf <- data.table::as.data.table(tf_membership)
    if (!"comparison_id" %in% names(topic_tf) && "condition_id" %in% names(topic_tf)) {
      topic_tf[, comparison_id := as.character(condition_id)]
    }
    if (!"direction" %in% names(topic_tf)) topic_tf[, direction := "All"]
    if (nrow(topic_tf) && all(c("comparison_id", "direction", "tf", "topic_num") %in% names(topic_tf))) {
      if (!"membership_pass" %in% names(topic_tf)) topic_tf[, membership_pass := TRUE]
      if (!"theta" %in% names(topic_tf)) topic_tf[, theta := NA_real_]
      if (!"primary_topic_num" %in% names(topic_tf)) {
        if ("primary_topic" %in% names(topic_tf)) {
          topic_tf[, primary_topic_num := suppressWarnings(as.integer(sub("^Topic", "", as.character(primary_topic))))]
        } else {
          topic_tf[, primary_topic_num := NA_integer_]
        }
      }
      topic_tf[, `:=`(
        comparison_id = as.character(comparison_id),
        direction = .m3tb_subgrn_direction_label(direction),
        tf_upper = toupper(as.character(tf)),
        topic_num = as.integer(topic_num),
        theta = suppressWarnings(as.numeric(theta)),
        primary_topic_num = suppressWarnings(as.integer(primary_topic_num)),
        membership_pass = as.logical(membership_pass)
      )]
      topic_tf[is.na(direction) | !nzchar(direction), direction := "All"]
      topic_tf <- topic_tf[is.finite(topic_num)]
      if (nrow(topic_tf)) {
        topic_tf[, topic_score_token := ifelse(
          is.finite(theta),
          paste0(topic_num, ":", as.character(signif(theta, 6L))),
          paste0(topic_num, ":")
        )]
      }
      topic_tf <- topic_tf[, {
        pass_topics <- sort(unique(topic_num[isTRUE(membership_pass) | membership_pass == TRUE]))
        primary_vals <- primary_topic_num[is.finite(primary_topic_num)]
        if (length(primary_vals)) {
          primary_num <- primary_vals[[1L]]
        } else {
          theta_vals <- theta
          finite_theta <- is.finite(theta_vals)
          primary_num <- if (any(finite_theta)) topic_num[which.max(ifelse(finite_theta, theta_vals, -Inf))] else NA_integer_
        }
        .(
          tf_topic_nums = paste(pass_topics, collapse = ";"),
          tf_topic_scores = paste(topic_score_token[order(topic_num)], collapse = ";"),
          tf_primary_topic_num = as.integer(primary_num)
        )
      }, by = .(comparison_id, direction, tf_upper)]
    } else {
      topic_tf <- data.table::data.table()
    }
  }
  if (nrow(topic_tf)) {
    joined <- merge(joined, topic_tf, by = c("comparison_id", "direction", "tf_upper"), all.x = TRUE, sort = FALSE)
  } else {
    joined[, tf_topic_nums := ""]
    joined[, tf_topic_scores := ""]
    joined[, tf_primary_topic_num := NA_integer_]
  }
  joined[is.na(tf_topic_nums), tf_topic_nums := ""]
  joined[is.na(tf_topic_scores), tf_topic_scores := ""]
  joined[, tf_primary_topic_num := suppressWarnings(as.integer(tf_primary_topic_num))]

  max_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  min_abs_finite <- function(x) {
    x <- suppressWarnings(abs(as.numeric(x)))
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  for (col in c("distance_to_tss", "log2fc_gene", "log2fc_tf", "delta_fp", "delta_gene")) {
    if (!col %in% names(joined)) joined[, (col) := NA_real_]
    joined[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  tf_gene <- joined[, {
    best_i <- which.max(edge_score_row)
    if (!length(best_i) || !is.finite(edge_score_row[best_i])) best_i <- 1L
    data.table::data.table(
      tf = tf[[1L]],
      tf_upper = tf_upper[[1L]],
      gene_key = gene_key[[1L]],
      edge_score = sum(edge_score_row, na.rm = TRUE),
      abs_edge_score = max_finite(edge_score_row),
      n_supporting_peaks = data.table::uniqueN(peak_id[!is.na(peak_id) & nzchar(peak_id)]),
      best_peak_id = peak_id[[best_i]],
      best_distance_to_tss = min_abs_finite(distance_to_tss),
      tf_topic_nums = tf_topic_nums[[1L]],
      tf_topic_scores = tf_topic_scores[[1L]],
      tf_primary_topic_num = tf_primary_topic_num[[1L]]
    )
  }, by = .(comparison_id, direction, tf_upper, gene_key)]
  tf_gene[is.na(n_supporting_peaks), n_supporting_peaks := 0L]
  data.table::setorder(tf_gene, comparison_id, direction, -abs_edge_score, tf_upper, gene_key)

  triplets <- joined[, .(
    tf = tf[[1L]],
    tf_upper = tf_upper[[1L]],
    peak_id = peak_id[[1L]],
    gene_key = gene_key[[1L]],
    edge_score = max_finite(edge_score_row),
    abs_edge_score = max_finite(edge_score_row),
    tf_topic_nums = tf_topic_nums[[1L]],
    tf_topic_scores = tf_topic_scores[[1L]],
    tf_primary_topic_num = tf_primary_topic_num[[1L]],
    distance_to_tss = min_abs_finite(distance_to_tss)
  ), by = .(comparison_id, direction, tf_upper, peak_id, gene_key)]
  data.table::setorder(triplets, comparison_id, direction, -abs_edge_score, tf_upper, peak_id, gene_key)

  topic_match <- function(topic_num, topic_nums) {
    parts <- strsplit(as.character(topic_nums), ";", fixed = TRUE)
    vapply(seq_along(parts), function(i) as.character(topic_num[[i]]) %in% parts[[i]], logical(1L))
  }
  tf_gene_count <- tf_gene[, .(comparison_id, direction, gene_key, tf_upper, tf_topic_nums)]
  data.table::setkey(tf_gene_count, comparison_id, direction, gene_key)
  data.table::setkey(context_genes, comparison_id, direction, gene_key)
  count_join <- tf_gene_count[context_genes, nomatch = 0L, allow.cartesian = TRUE]
  if (nrow(count_join)) {
    count_join[, context_topic_tf := topic_match(topic_num, tf_topic_nums)]
    tf_gene_counts <- count_join[, .(
      n_tf_gene_edges = .N,
      n_topic_tfs = data.table::uniqueN(tf_upper[context_topic_tf])
    ), by = subgrn_context_id]
  } else {
    tf_gene_counts <- data.table::data.table(subgrn_context_id = character(), n_tf_gene_edges = integer(), n_topic_tfs = integer())
  }
  triplet_count <- triplets[, .(comparison_id, direction, gene_key, tf_upper)]
  data.table::setkey(triplet_count, comparison_id, direction, gene_key)
  triplet_join <- triplet_count[context_genes, nomatch = 0L, allow.cartesian = TRUE]
  triplet_counts <- if (nrow(triplet_join)) {
    triplet_join[, .(n_tf_peak_gene_triplets = .N), by = subgrn_context_id]
  } else {
    data.table::data.table(subgrn_context_id = character(), n_tf_peak_gene_triplets = integer())
  }
  manifest <- merge(manifest, tf_gene_counts, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  manifest <- merge(manifest, triplet_counts, by = "subgrn_context_id", all.x = TRUE, sort = FALSE)
  for (col in c("n_tf_gene_edges", "n_topic_tfs", "n_tf_peak_gene_triplets")) {
    manifest[is.na(get(col)), (col) := 0L]
    manifest[, (col) := as.integer(get(col))]
  }
  list(
    manifest = manifest,
    context_genes = context_genes,
    tf_gene_edges = tf_gene,
    tf_peak_gene_triplets = triplets
  )
}

.m3tb_subgrn_payload_stem <- function(payload_name) {
  key <- digest::digest(
    as.character(payload_name)[[1L]],
    algo = "xxhash64",
    serialize = FALSE
  )
  paste0("sg_", substr(key, 1L, 12L))
}

.m3tb_write_pathway_subgrn_payload <- function(extraction_dir,
                                               model_dir,
                                               payload_dir,
                                               payload_name,
                                               max_tf_gene_edges_per_context = 300L,
                                               max_tf_peak_gene_triplets_per_context = 600L,
                                               payload_contexts_per_file = 100L,
                                               max_payload_contexts = Inf,
                                               condition_pathways = NULL,
                                               topic_pathways = NULL) {
  dir.create(payload_dir, recursive = TRUE, showWarnings = FALSE)
  legacy_payload_name <- as.character(payload_name)[[1L]]
  payload_name <- .m3tb_subgrn_payload_stem(legacy_payload_name)
  unlink(
    file.path(payload_dir, paste0(legacy_payload_name, "_chunk_*.js")),
    force = TRUE
  )
  unlink(
    file.path(payload_dir, paste0(payload_name, "_chunk_*.js")),
    force = TRUE
  )
  pathways_path <- file.path(extraction_dir, "per_comparison_topic_pathway_enrichment.csv")
  if (!file.exists(pathways_path)) {
    pathways_path <- file.path(extraction_dir, "per_condition_topic_pathway_enrichment.csv")
  }
  edges_path <- file.path(model_dir, "rds", "edges_docs.rds")
  condition_links_dir <- .m3tb_subgrn_condition_links_dir(model_dir = model_dir, extraction_dir = extraction_dir)
  tf_assignment_path <- file.path(extraction_dir, "tf_topic_assignments.csv")
  tf_pass_path <- file.path(extraction_dir, "tf_topic_membership_pass.csv")
  tf_path <- if (file.exists(tf_assignment_path)) tf_assignment_path else tf_pass_path
  if (!file.exists(pathways_path) || (!file.exists(edges_path) && is.na(condition_links_dir))) {
    return(data.table::data.table())
  }
  pathways <- data.table::fread(pathways_path, showProgress = FALSE)
  if (!"comparison_id" %in% names(pathways) && "condition_id" %in% names(pathways)) {
    pathways[, comparison_id := condition_id]
  }
  if (is.null(condition_pathways)) {
    condition_pathways <- .m3tb_read_pathway_tables(
      extraction_dir,
      per_group = TRUE,
      model_dir = NULL,
      compute_universe = FALSE
    )
  }
  if (is.null(topic_pathways)) {
    topic_pathways <- .m3tb_read_pathway_tables(
      extraction_dir,
      per_group = FALSE,
      model_dir = NULL,
      compute_universe = FALSE
    )
  }
  pathways <- .m3tb_select_pathway_subgrn_contexts(
    all_pathways = pathways,
    condition_pathways = condition_pathways,
    topic_pathways = topic_pathways
  )
  pathways <- pathways[!is.na(overlap_genes) & nzchar(overlap_genes)]
  if (!nrow(pathways)) return(data.table::data.table())
  max_context_env <- Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_MAX_CONTEXTS", unset = "")
  max_payload_contexts <- suppressWarnings(as.integer(max_payload_contexts)[[1L]])
  if (nzchar(max_context_env)) {
    max_payload_contexts <- suppressWarnings(as.integer(max_context_env[[1L]]))
  }
  if (is.finite(max_payload_contexts) && max_payload_contexts <= 0L) {
    return(data.table::data.table())
  }
  top_per_group_env <- Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_TOP_PER_GROUP", unset = "")
  top_per_group <- suppressWarnings(as.integer(top_per_group_env[[1L]]))
  if (!nzchar(top_per_group_env) && !is.na(condition_links_dir)) {
    top_per_group <- 30L
  }
  if (is.finite(top_per_group) && top_per_group > 0L) {
    if (!"overlap_hits" %in% names(pathways)) {
      pathways[, overlap_hits := lengths(.m3tb_subgrn_split_genes(overlap_genes))]
    }
    if (!"combined_score" %in% names(pathways)) pathways[, combined_score := NA_real_]
    pathways[, `:=`(
      padj = suppressWarnings(as.numeric(padj)),
      overlap_hits = suppressWarnings(as.integer(overlap_hits)),
      combined_score = suppressWarnings(as.numeric(combined_score))
    )]
    data.table::setorder(pathways, comparison_id, direction_group, topic, padj, -overlap_hits, -combined_score, pathway)
    pathways[, subgrn_group_rank := seq_len(.N), by = .(comparison_id, direction_group, topic)]
    pathways <- pathways[subgrn_group_rank <= top_per_group][, subgrn_group_rank := NULL][]
  }
  if (is.finite(max_payload_contexts) && max_payload_contexts > 0L && nrow(pathways) > max_payload_contexts) {
    if (!"overlap_hits" %in% names(pathways)) {
      pathways[, overlap_hits := lengths(.m3tb_subgrn_split_genes(overlap_genes))]
    }
    if (!"combined_score" %in% names(pathways)) pathways[, combined_score := NA_real_]
    pathways[, `:=`(
      padj = suppressWarnings(as.numeric(padj)),
      overlap_hits = suppressWarnings(as.integer(overlap_hits)),
      combined_score = suppressWarnings(as.numeric(combined_score))
    )]
    data.table::setorder(pathways, padj, -overlap_hits, -combined_score, comparison_id, direction_group, topic, pathway)
    pathways <- pathways[seq_len(max_payload_contexts)]
  }
  tf_membership <- if (file.exists(tf_path)) data.table::fread(tf_path, showProgress = FALSE) else NULL
  edges_docs <- .m3tb_read_condition_links_for_subgrn(
    condition_links_dir,
    pathways,
    tf_membership = tf_membership,
    cache_dir = if (!is.na(condition_links_dir)) file.path(payload_dir, "cache") else NULL
  )
  if (!nrow(edges_docs)) {
    edges_docs <- readRDS(edges_path)
  }
  payload <- .m3tb_build_pathway_subgrn_compact_payload(
    pathways = pathways,
    edges_docs = edges_docs,
    tf_membership = tf_membership
  )
  manifest <- payload$manifest
  if (!nrow(manifest)) return(manifest)
  payload_contexts_per_file <- max(1L, as.integer(payload_contexts_per_file)[[1L]])
  if (!is.na(condition_links_dir)) {
    manifest[, payload_group := paste(comparison_id, direction, topic_num, sep = "__")]
    group_ids <- sort(unique(manifest$payload_group))
    group_map <- data.table::data.table(
      payload_group = group_ids,
      payload_chunk = seq_along(group_ids)
    )
    manifest <- merge(manifest, group_map, by = "payload_group", all.x = TRUE, sort = FALSE)
  } else {
    manifest[, payload_chunk := ceiling(seq_len(.N) / payload_contexts_per_file)]
  }
  manifest[, payload_file := sprintf("%s_chunk_%03d.js", payload_name, as.integer(payload_chunk))]
  for (pf in unique(manifest$payload_file)) {
    context_ids <- manifest[payload_file == pf, subgrn_context_id]
    chunk_genes <- unique(payload$context_genes[subgrn_context_id %chin% context_ids, .(comparison_id, direction, gene_key)])
    tf_gene_edges <- payload$tf_gene_edges
    tf_peak_gene_triplets <- payload$tf_peak_gene_triplets
    if (nrow(chunk_genes)) {
      data.table::setkey(tf_gene_edges, comparison_id, direction, gene_key)
      data.table::setkey(tf_peak_gene_triplets, comparison_id, direction, gene_key)
      data.table::setkey(chunk_genes, comparison_id, direction, gene_key)
      tf_gene_edges <- tf_gene_edges[chunk_genes, nomatch = 0L, allow.cartesian = TRUE]
      tf_peak_gene_triplets <- tf_peak_gene_triplets[chunk_genes, nomatch = 0L, allow.cartesian = TRUE]
    } else {
      tf_gene_edges <- tf_gene_edges[0]
      tf_peak_gene_triplets <- tf_peak_gene_triplets[0]
    }
    obj <- list(
      manifest = manifest[subgrn_context_id %chin% context_ids],
      tf_gene_edges = tf_gene_edges,
      tf_peak_gene_triplets = tf_peak_gene_triplets
    )
    obj$manifest[, payload_chunk := NULL]
    if ("payload_group" %in% names(obj$manifest)) obj$manifest[, payload_group := NULL]
    packed <- lapply(obj, .module2_report_browser_browser_payload_to_columnar)
    payload_b64 <- .module2_report_browser_encode_browser_json_deflate_base64(packed)
    js <- paste0(
      "window.CRAFTGRN_SUBGRN_PAYLOADS=window.CRAFTGRN_SUBGRN_PAYLOADS||{};\n",
      "window.CRAFTGRN_SUBGRN_PAYLOADS[", jsonlite::toJSON(pf, auto_unbox = TRUE), "]=",
      "{compressed_columnar:'", payload_b64, "'}",
      ";\n"
    )
    writeLines(js, file.path(payload_dir, pf), useBytes = TRUE)
  }
  manifest[, payload_chunk := NULL]
  if ("payload_group" %in% names(manifest)) manifest[, payload_group := NULL]
  manifest[]
}

.m3tb_subgrn_empty_manifest <- function() {
  data.table::data.table(
    subgrn_context_id = character(),
    comparison_id = character(),
    direction = character(),
    comparison_label = character(),
    topic_num = integer(),
    pathway_key = character(),
    pathway = character(),
    padj = numeric(),
    n_overlap_genes = integer(),
    n_tf_gene_edges = integer(),
    n_tf_peak_gene_triplets = integer(),
    n_topic_tfs = integer(),
    payload_file = character()
  )
}

.m3tb_subgrn_condition_links_dir <- function(model_dir, extraction_dir = NULL) {
  ancestor_candidates <- function(path, max_depth = 8L) {
    if (is.null(path) || !length(path) ||
        !nzchar(as.character(path)[[1L]])) {
      return(character())
    }
    current <- dirname(normalizePath(
      as.character(path)[[1L]],
      winslash = "/",
      mustWork = FALSE
    ))
    out <- character()
    for (depth in seq_len(max_depth)) {
      out <- c(out, file.path(current, "condition_links"))
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    out
  }
  candidates <- c(
    ancestor_candidates(model_dir),
    ancestor_candidates(extraction_dir)
  )
  candidates <- unique(candidates)
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) hit[[1L]] else NA_character_
}

.m3tb_read_condition_links_for_subgrn <- function(condition_links_dir,
                                                   pathways,
                                                   tf_membership = NULL,
                                                   cache_dir = NULL) {
  empty <- data.table::data.table()
  if (is.na(condition_links_dir) || !dir.exists(condition_links_dir)) return(empty)
  if (!requireNamespace("arrow", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) return(empty)
  pathways <- data.table::as.data.table(pathways)
  if (!nrow(pathways) || !"comparison_id" %in% names(pathways)) return(empty)
  if (!"overlap_genes" %in% names(pathways)) pathways[, overlap_genes := ""]
  if (!"topic_num" %in% names(pathways) && "topic" %in% names(pathways)) {
    pathways[, topic_num := suppressWarnings(as.integer(sub("^Topic", "", as.character(topic))))]
  }
  gene_lists <- .m3tb_subgrn_split_genes(pathways$overlap_genes)
  needed <- data.table::rbindlist(lapply(seq_len(nrow(pathways)), function(i) {
    genes <- gene_lists[[i]]
    if (!length(genes)) return(NULL)
    data.table::data.table(
      condition_id = as.character(pathways$comparison_id[[i]]),
      topic_num = suppressWarnings(as.integer(pathways$topic_num[[i]])),
      gene_key = genes
    )
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(needed)) return(empty)
  needed <- unique(needed[
    !is.na(condition_id) & nzchar(condition_id) &
      is.finite(topic_num) &
      !is.na(gene_key) & nzchar(gene_key)
  ])
  manifest_path <- file.path(condition_links_dir, "condition_links_manifest.csv")
  if (file.exists(manifest_path)) {
    manifest <- data.table::fread(manifest_path, showProgress = FALSE)
    if (!all(c("condition_id", "path") %in% names(manifest))) return(empty)
    manifest[, path := .module3_resolve_condition_manifest_paths(path, manifest_path)]
  } else {
    manifest <- data.table::data.table(
      condition_id = sub("_condition_links[.]parquet$", "", basename(list.files(condition_links_dir, "_condition_links[.]parquet$", full.names = TRUE))),
      path = list.files(condition_links_dir, "_condition_links[.]parquet$", full.names = TRUE)
    )
  }
  manifest <- manifest[condition_id %in% unique(needed$condition_id) & file.exists(path)]
  if (!nrow(manifest)) return(empty)
  max_tfs_per_gene <- suppressWarnings(as.integer(Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_MAX_TFS_PER_GENE", unset = "10")))
  max_peaks_per_tf_gene <- suppressWarnings(as.integer(Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_MAX_PEAKS_PER_TF_GENE", unset = "1")))
  if (!is.finite(max_tfs_per_gene) || max_tfs_per_gene < 1L) max_tfs_per_gene <- 10L
  if (!is.finite(max_peaks_per_tf_gene) || max_peaks_per_tf_gene < 1L) max_peaks_per_tf_gene <- 1L
  if (!is.null(cache_dir)) {
    cache_dir <- normalizePath(as.character(cache_dir)[[1L]], winslash = "/", mustWork = FALSE)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  priority_tfs <- data.table::data.table()
  if (!is.null(tf_membership)) {
    tm <- data.table::as.data.table(tf_membership)
    if (!"comparison_id" %in% names(tm) && "condition_id" %in% names(tm)) {
      tm[, comparison_id := condition_id]
    }
    if (!"topic_num" %in% names(tm) && "topic" %in% names(tm)) {
      tm[, topic_num := suppressWarnings(as.integer(sub("^Topic", "", as.character(topic))))]
    }
    if (!"membership_pass" %in% names(tm)) tm[, membership_pass := TRUE]
    if (nrow(tm) && all(c("comparison_id", "tf", "topic_num", "membership_pass") %in% names(tm))) {
      tm[, `:=`(
        condition_id = as.character(comparison_id),
        topic_num = suppressWarnings(as.integer(topic_num)),
        tf_upper_for_rank = toupper(as.character(tf)),
        membership_pass = as.logical(membership_pass)
      )]
      tm[is.na(membership_pass), membership_pass := TRUE]
      tm <- unique(tm[
        !is.na(condition_id) & nzchar(condition_id) &
          is.finite(topic_num) &
          !is.na(tf_upper_for_rank) & nzchar(tf_upper_for_rank) &
          membership_pass == TRUE,
        .(condition_id, topic_num, tf_upper_for_rank)
      ])
      if (nrow(tm)) {
        priority_tfs <- merge(
          unique(needed[, .(condition_id, topic_num, gene_key)]),
          tm,
          by = c("condition_id", "topic_num"),
          all = FALSE,
          allow.cartesian = TRUE,
          sort = FALSE
        )
        priority_tfs <- unique(priority_tfs[, .(condition_id, gene_key, tf_upper_for_rank)])
      }
    }
  }
  cols <- c(
    "link_id", "tf", "fp_id", "target_gene", "condition_id", "condition_label",
    "gene_key", "peak_id", "fp_score_condition", "gene_expr_condition",
    "tf_expr_condition", "active_link", "doc_id"
  )
  collapse_best_peaks <- function(parts) {
    parts <- parts[!vapply(parts, is.null, logical(1L))]
    if (!length(parts)) return(data.table::data.table())
    out <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
    if (!nrow(out)) return(out)
    out[, edge_score_for_rank := suppressWarnings(as.numeric(fp_score_condition))]
    out[!is.finite(edge_score_for_rank), edge_score_for_rank := 0]
    out[, tf_upper_for_rank := toupper(as.character(tf))]
    data.table::setorder(
      out,
      condition_id, gene_key, tf_upper_for_rank,
      -edge_score_for_rank, peak_id
    )
    out[, peak_rank_for_subgrn := seq_len(.N), by = .(
      condition_id, gene_key, tf_upper_for_rank
    )]
    out <- out[peak_rank_for_subgrn <= max_peaks_per_tf_gene]
    out[, peak_rank_for_subgrn := NULL]
    out
  }
  read_best_peaks <- function(path, genes) {
    ds <- arrow::open_dataset(path)
    query <- dplyr::select(
      dplyr::filter(ds, gene_key %in% genes),
      dplyr::any_of(cols)
    )
    if (max_peaks_per_tf_gene == 1L) {
      max_scores <- tryCatch({
        grouped <- dplyr::group_by(query, condition_id, gene_key, tf)
        summarized <- dplyr::summarise(
          grouped,
          max_score = max(fp_score_condition, na.rm = TRUE),
          .groups = "drop"
        )
        dplyr::collect(summarized)
      }, error = function(e) data.frame())
      if (!nrow(max_scores)) return(data.table::data.table())
      best_query <- dplyr::inner_join(
        query,
        arrow::arrow_table(max_scores),
        by = c("condition_id", "gene_key", "tf")
      )
      out <- tryCatch(
        dplyr::collect(dplyr::filter(best_query, fp_score_condition == max_score)),
        error = function(e) data.frame()
      )
      out <- data.table::as.data.table(out)
      if ("max_score" %in% names(out)) out[, max_score := NULL]
      return(collapse_best_peaks(list(out)))
    }
    reader <- tryCatch(
      arrow::as_record_batch_reader(query),
      error = function(e) NULL
    )
    if (is.null(reader)) return(data.table::data.table())
    parts <- list()
    repeat {
      batch <- reader$read_next_batch()
      if (is.null(batch)) break
      if (nrow(batch)) {
        part <- data.table::as.data.table(as.data.frame(batch))
        parts[[length(parts) + 1L]] <- collapse_best_peaks(list(part))
      }
      if (length(parts) >= 32L) {
        parts <- list(collapse_best_peaks(parts))
      }
    }
    collapse_best_peaks(parts)
  }
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    cond <- manifest$condition_id[[i]]
    genes <- needed[condition_id == cond, unique(gene_key)]
    if (!length(genes)) return(NULL)
    source_path <- normalizePath(manifest$path[[i]], winslash = "/", mustWork = TRUE)
    source_info <- file.info(source_path)
    cache_path <- if (is.null(cache_dir)) {
      NA_character_
    } else {
      file.path(
        cache_dir,
        sprintf("%s_peaks%d.rds", .safe_filename(cond), max_peaks_per_tf_gene)
      )
    }
    cached <- NULL
    if (!is.na(cache_path) && file.exists(cache_path)) {
      cached <- tryCatch(readRDS(cache_path), error = function(e) NULL)
      valid_cache <- is.list(cached) &&
        identical(as.integer(cached$cache_version), 2L) &&
        identical(as.character(cached$source_path), as.character(source_path)) &&
        isTRUE(all.equal(as.numeric(cached$source_size), as.numeric(source_info$size))) &&
        isTRUE(all.equal(as.numeric(cached$source_mtime), as.numeric(source_info$mtime))) &&
        identical(as.integer(cached$max_peaks_per_tf_gene), as.integer(max_peaks_per_tf_gene)) &&
        data.table::is.data.table(cached$rows)
      if (!isTRUE(valid_cache)) cached <- NULL
    }
    cached_rows <- if (is.null(cached)) data.table::data.table() else cached$rows
    cached_genes <- if (nrow(cached_rows)) unique(as.character(cached_rows$gene_key)) else character()
    missing_genes <- setdiff(genes, cached_genes)
    fresh_rows <- if (length(missing_genes)) read_best_peaks(source_path, missing_genes) else data.table::data.table()
    cache_cols <- c(
      "tf", "condition_id", "condition_label", "gene_key", "peak_id",
      "fp_score_condition"
    )
    if (nrow(fresh_rows)) {
      fresh_rows <- fresh_rows[, intersect(cache_cols, names(fresh_rows)), with = FALSE]
    }
    out <- collapse_best_peaks(list(cached_rows, fresh_rows))
    if (!is.na(cache_path) && length(missing_genes) && nrow(out)) {
      cache_obj <- list(
        cache_version = 2L,
        source_path = source_path,
        source_size = as.numeric(source_info$size),
        source_mtime = as.numeric(source_info$mtime),
        max_peaks_per_tf_gene = as.integer(max_peaks_per_tf_gene),
        rows = out
      )
      tmp_path <- tempfile(pattern = paste0(basename(cache_path), "."), tmpdir = cache_dir)
      saveRDS(cache_obj, tmp_path, compress = "gzip")
      if (!file.rename(tmp_path, cache_path)) unlink(tmp_path)
    }
    out <- out[gene_key %in% genes]
    if (!nrow(out)) return(NULL)
    if (!"comparison_id" %in% names(out)) out[, comparison_id := as.character(condition_id)]
    out[, `:=`(
      direction = "All",
      direction_group = "All"
    )]
    tf_rank <- out[, .(
      tf_best_score = max(edge_score_for_rank, na.rm = TRUE)
    ), by = .(condition_id, gene_key, tf_upper_for_rank)]
    data.table::setorder(tf_rank, condition_id, gene_key, -tf_best_score, tf_upper_for_rank)
    tf_rank[, tf_rank_for_subgrn := seq_len(.N), by = .(condition_id, gene_key)]
    if (nrow(priority_tfs)) {
      priority_cond <- priority_tfs[condition_id == cond]
      if (nrow(priority_cond)) {
        priority_cond[, priority_tf := TRUE]
        tf_rank <- merge(
          tf_rank,
          priority_cond,
          by = c("condition_id", "gene_key", "tf_upper_for_rank"),
          all.x = TRUE,
          sort = FALSE
        )
      } else {
        tf_rank[, priority_tf := FALSE]
      }
    } else {
      tf_rank[, priority_tf := FALSE]
    }
    tf_rank[is.na(priority_tf), priority_tf := FALSE]
    data.table::setorder(
      tf_rank,
      condition_id, gene_key, -priority_tf, -tf_best_score, tf_upper_for_rank
    )
    tf_rank[, tf_rank_for_subgrn := seq_len(.N), by = .(condition_id, gene_key)]
    tf_rank <- tf_rank[
      tf_rank_for_subgrn <= max_tfs_per_gene,
      .(condition_id, gene_key, tf_upper_for_rank)
    ]
    out <- merge(
      out,
      tf_rank,
      by = c("condition_id", "gene_key", "tf_upper_for_rank"),
      all = FALSE,
      sort = FALSE
    )
    out[, c("edge_score_for_rank", "tf_upper_for_rank") := NULL]
    out
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.m3tb_subgrn_css <- function() {
  paste0(
    ".subgrnPanel{position:absolute;inset:0;background:#fff;display:none;flex-direction:column;z-index:4}",
    ".subgrnTop{border-bottom:1px solid #e5e5df;padding:7px 9px;background:#fbfbfa}",
    ".subgrnTitle{font-size:14px;font-weight:700;line-height:1.2;margin-bottom:6px;color:#111}",
    ".subgrnControls{display:flex;gap:6px 8px;align-items:center;flex-wrap:wrap;font-size:11px}",
    ".subgrnControls select,.subgrnControls input{font:700 11px Arial,Helvetica,sans-serif;border:1px solid #aaa;background:#fff;color:#111;border-radius:3px;padding:4px 6px}",
    ".subgrnControls input{width:58px}.subgrnControls input[type=range]{width:82px}.subgrnControls button{font:700 11px Arial,Helvetica,sans-serif;border:1px solid #777;background:#222;color:#fff;border-radius:3px;padding:5px 8px;cursor:pointer}",
    ".subgrnStats{font-size:10px;line-height:1.25;color:#555;margin-top:5px}.subgrnCanvas{position:relative;flex:1;min-height:0}",
    ".subgrnCanvas.panning{cursor:grabbing}.subgrnNode{cursor:grab}.subgrnNode.selected{stroke:#111;stroke-width:4px}",
    ".subgrnNodeTf{stroke:none}.subgrnNodeGene{fill:#9ca3af}.subgrnNodePeak{fill:#f59e0b}.subgrnEdge{stroke:#9ca3af;stroke-width:1.1;stroke-opacity:.52;fill:none;stroke-linecap:round}.subgrnLabel{font-family:Arial,Helvetica,sans-serif;font-weight:700;font-size:12px;paint-order:stroke;stroke:#fff;stroke-width:3px;stroke-linejoin:round;fill:#111;pointer-events:none}.subgrnTfLabel{fill:#fff;stroke:#555}.subgrnTooltip{position:absolute;display:none;background:rgba(17,17,17,0.92);color:#fff;font:700 11px Arial,Helvetica,sans-serif;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:340px;line-height:1.35;z-index:6}"
  )
}

.m3tb_subgrn_html_panel <- function() {
  paste0(
    "<div id=\"subgrnPanel\" class=\"subgrnPanel\"><div class=\"subgrnTop\">",
    "<div id=\"subgrnTitle\" class=\"subgrnTitle\">Per Condition / Topic / Pathway Sub GRN</div>",
    "<div class=\"subgrnControls\">",
    "<label>Comparison <select id=\"subgrnComparisonSelect\"></select></label>",
    "<label>TF scope <select id=\"tfScopeSelect\"><option value=\"all\" selected>All TFs</option><option value=\"topic\">Topic TFs</option></select></label>",
    "<label>Topic theta >= <select id=\"subgrnTopicThetaPreset\"><option value=\"0.3\" selected>0.3</option><option value=\"0.5\">0.5</option><option value=\"0.7\">0.7</option><option value=\"custom\">Custom</option></select><input id=\"subgrnTopicThetaCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"0.3\"/></label>",
    "<label>Topic Phi >= <select id=\"subgrnTopicPhiPreset\"><option value=\"0\" selected>0</option><option value=\"0.1\">0.1</option><option value=\"0.2\">0.2</option><option value=\"0.3\">0.3</option><option value=\"0.5\">0.5</option><option value=\"0.7\">0.7</option><option value=\"custom\">Custom</option></select><input id=\"subgrnTopicPhiCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"0\"/></label>",
    "<label>Primary topic only <input id=\"subgrnPrimaryTopicOnly\" type=\"checkbox\"/></label>",
    "<label>Network <select id=\"networkModeSelect\"><option value=\"tf_gene\" selected>TF-gene</option><option value=\"tf_peak_gene\">TF-peak-gene</option></select></label>",
    "<label>Top TFs <input id=\"subgrnTopTfN\" type=\"number\" min=\"0\" value=\"100\"/></label>",
    "<label>Top links <input id=\"subgrnTopLinkN\" type=\"number\" min=\"0\" value=\"300\"/></label>",
    "<label>Layout <select id=\"subgrnLayoutSelect\"><option value=\"force\" selected>Force</option><option value=\"radial\">Radial</option><option value=\"columns\">Columns</option><option value=\"bipartite\">Bipartite</option><option value=\"hierarchy\">Hierarchy</option><option value=\"concentric\">Concentric</option><option value=\"circle\">Circle</option><option value=\"grid\">Grid</option><option value=\"spiral\">Spiral</option><option value=\"clustered\">Clustered</option></select></label>",
    "<label>Spacing <input id=\"subgrnSpacingRange\" type=\"range\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/><input id=\"subgrnSpacingValue\" type=\"number\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/></label>",
    "<label>TF box <input id=\"subgrnTfBoxMin\" type=\"number\" min=\"6\" max=\"80\" step=\"1\" value=\"14\"/><input id=\"subgrnTfBoxMax\" type=\"number\" min=\"8\" max=\"110\" step=\"1\" value=\"20\"/></label>",
    "<label>Labels <input id=\"subgrnShowLabels\" type=\"checkbox\" checked/></label>",
    "<label>Arrows <input id=\"subgrnShowArrows\" type=\"checkbox\" checked/></label>",
    "<label>Palette <select id=\"subgrnPaletteSelect\"><option value=\"default\" selected>Default</option><option value=\"npg\">NPG inspired</option><option value=\"aaas\">AAAS inspired</option><option value=\"nejm\">NEJM inspired</option><option value=\"lancet\">Lancet inspired</option><option value=\"jama\">JAMA inspired</option><option value=\"bmj\">BMJ inspired</option><option value=\"jco\">JCO inspired</option><option value=\"ucscgb\">UCSCGB inspired</option><option value=\"d3\">D3 inspired</option><option value=\"gephi\">Gephi inspired</option><option value=\"observable\">Observable inspired</option><option value=\"primer\">Primer inspired</option><option value=\"atlassian\">Atlassian inspired</option><option value=\"iterm\">iTerm inspired</option><option value=\"locuszoom\">LocusZoom inspired</option><option value=\"igv\">IGV inspired</option><option value=\"cosmic\">COSMIC inspired</option><option value=\"uchicago\">UChicago inspired</option><option value=\"startrek\">Star Trek inspired</option><option value=\"tron\">Tron inspired</option><option value=\"futurama\">Futurama inspired</option><option value=\"rickandmorty\">Rick and Morty inspired</option><option value=\"simpsons\">Simpsons inspired</option><option value=\"flatui\">Flat UI inspired</option><option value=\"frontiers\">Frontiers inspired</option><option value=\"gsea\">GSEA inspired</option><option value=\"bootstrap5\">Bootstrap 5 inspired</option><option value=\"material\">Material Design inspired</option><option value=\"tailwind\">Tailwind CSS inspired</option></select></label>",
    "<label>Select node <select id=\"subgrnNodeSelect\"></select></label>",
    "<button id=\"subgrnResetButton\" type=\"button\">Reset layout</button>",
    "<button id=\"subgrnBackButton\" type=\"button\">Back to pathways</button>",
    "<button id=\"subgrnExportButton\" type=\"button\">Export SVG</button>",
    "</div><div id=\"subgrnStats\" class=\"subgrnStats\"></div></div>",
    "<div id=\"subgrnCanvas\" class=\"subgrnCanvas\"><div id=\"subgrnTooltip\" class=\"subgrnTooltip\"></div><svg id=\"subgrnSvg\" viewBox=\"0 0 1600 900\"><defs><marker id=\"subgrnArrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"6\" markerHeight=\"6\" orient=\"auto-start-reverse\"><path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#6b7280\"/></marker></defs><g id=\"subgrnViewLayer\"><g id=\"subgrnEdgeLayer\"></g><g id=\"subgrnNodeLayer\"></g><g id=\"subgrnLabelLayer\"></g></g></svg></div></div>"
  )
}

.m3tb_subgrn_js <- function() {
  c(
    "window.CRAFTGRN_SUBGRN_PAYLOADS=window.CRAFTGRN_SUBGRN_PAYLOADS||{};let subgrnCurrentRows=[],subgrnCurrent=null;const SUBGRN_W=1600,SUBGRN_H=900;let subgrnState={nodes:[],edges:[],drag:null,pan:null,selected:'',view:{x:0,y:0,k:1}};",
    "const SUBGRN_COLOR_PRESETS={default:['#2563eb','#dc2626','#9ca3af','#f59e0b'],npg:['#E64B35','#4DBBD5','#7E6148','#3C5488'],aaas:['#EE0000','#3B4992','#008B45','#631879'],nejm:['#BC3C29','#0072B5','#E18727','#20854E'],lancet:['#ED0000','#00468B','#42B540','#0099B4'],jama:['#374E55','#DF8F44','#B24745','#6A6599'],bmj:['#2B8CBE','#E34A33','#969696','#636363'],jco:['#CD534C','#0073C2','#868686','#EFC000'],ucscgb:['#FF0000','#0000FF','#999999','#666666'],d3:['#FF7F0E','#1F77B4','#7F7F7F','#2CA02C'],gephi:['#FF7F00','#377EB8','#999999','#4DAF4A'],observable:['#4269D0','#EF553B','#AAAAAA','#6CC5B0'],primer:['#CF222E','#0969DA','#8C959F','#57606A'],atlassian:['#FF5630','#0052CC','#97A0AF','#42526E'],iterm:['#CC241D','#458588','#A89984','#665C54'],locuszoom:['#D7191C','#2C7BB6','#ABD9E9','#7570B3'],igv:['#E41A1C','#377EB8','#999999','#4DAF4A'],cosmic:['#D73027','#4575B4','#969696','#542788'],uchicago:['#800000','#155F83','#8A9045','#767676'],startrek:['#CC0C00','#5C88DA','#9C9C9C','#FFCC00'],tron:['#FF410D','#0085CA','#6EE2FF','#F7C530'],futurama:['#FF6F00','#1B9E77','#A6A6A6','#7570B3'],rickandmorty:['#FAFD7C','#00B0C8','#B2DF8A','#808080'],simpsons:['#FED439','#709AE1','#8A9197','#D2AF81'],flatui:['#E74C3C','#3498DB','#95A5A6','#2C3E50'],frontiers:['#E64B35','#4DBBD5','#B09C85','#3C5488'],gsea:['#E31A1C','#1F78B4','#BDBDBD','#33A02C'],bootstrap5:['#DC3545','#0D6EFD','#6C757D','#198754'],material:['#F44336','#2196F3','#9E9E9E','#4CAF50'],tailwind:['#EF4444','#3B82F6','#9CA3AF','#10B981']};",
    "function subgrnNum(x,d){const v=Number(x);return Number.isFinite(v)?v:d}function subgrnLimit(el,d){const v=Math.floor(subgrnNum(el.value,d));return v===0?Infinity:Math.max(1,v)}function subgrnEl(id){return document.getElementById(id)}function subgrnSpacing(){return Math.max(.5,Math.min(2,subgrnNum(subgrnEl('subgrnSpacingValue').value,1)))}function subgrnColors(){const p=subgrnEl('subgrnPaletteSelect').value;return SUBGRN_COLOR_PRESETS[p]||SUBGRN_COLOR_PRESETS.default}function subgrnHexRgb(hex){const x=String(hex||'#777777').replace('#','');const v=parseInt(x,16);return{r:(v>>16)&255,g:(v>>8)&255,b:v&255}}function subgrnMixRgb(a,b,t){return 'rgb('+Math.round(a.r+(b.r-a.r)*t)+','+Math.round(a.g+(b.g-a.g)*t)+','+Math.round(a.b+(b.b-a.b)*t)+')'}function subgrnEdgeColor(v,lo,hi){const t=hi>lo?(v-lo)/(hi-lo):.8;return subgrnMixRgb({r:226,g:226,b:226},subgrnHexRgb(subgrnColors()[3]),Math.max(0,Math.min(1,t)))}",
    "function subgrnContextRows(d){if(!SUBGRN_MANIFEST.length)return[];const tn=topicNum(d.topic_num||d.topic),pathKey=x=>String(x||'').trim().toLowerCase().replace(/[^a-z0-9]+/g,''),pk=pathKey(d.pathway_key||d.pathway_norm_key||d.pathway||'');let rows=SUBGRN_MANIFEST.filter(m=>+m.topic_num===tn&&pathKey(m.pathway_key||m.pathway||'')===pk);if(typeof groupSelect!=='undefined'&&groupSelect&&groupSelect.value&&String(groupSelect.value)!=='__overall__'){const gv=String(groupSelect.value),norm=s=>String(s||'').replace(/\\s+/g,'_');rows=rows.filter(m=>String(m.comparison_label||'')===gv||String(m.comparison_id||'')===gv||norm(m.comparison_label)===norm(gv)||norm(m.comparison_id)===norm(gv))}rows=rows.slice().sort((a,b)=>subgrnNum(b.n_overlap_genes,0)-subgrnNum(a.n_overlap_genes,0)||subgrnNum(a.padj,1)-subgrnNum(b.padj,1)||String(a.comparison_label||'').localeCompare(String(b.comparison_label||'')));return rows}",
    "function subgrnColumnarRows(x){const cols=x&&x.columns?x.columns:[],data=x&&x.data?x.data:[],n=data.length?data[0].length:0,rows=new Array(n);for(let i=0;i<n;i++){const row={};for(let j=0;j<cols.length;j++)row[cols[j]]=data[j][i];rows[i]=row}return rows}async function subgrnDecodePayload(x){if(!x||!x.compressed_columnar)return x;if(!('DecompressionStream' in window))throw new Error('This browser cannot decompress Sub-GRN data. Use a current Chrome, Edge, Firefox, or Safari.');const bin=atob(x.compressed_columnar),bytes=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)bytes[i]=bin.charCodeAt(i);const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate')),packed=JSON.parse(await new Response(stream).text()),out={};Object.keys(packed).forEach(k=>out[k]=subgrnColumnarRows(packed[k]));return out}",
    "function loadSubgrnPayload(file){return new Promise((resolve,reject)=>{if(!file)return resolve(null);const finish=()=>subgrnDecodePayload(window.CRAFTGRN_SUBGRN_PAYLOADS[file]||null).then(x=>{if(x)window.CRAFTGRN_SUBGRN_PAYLOADS[file]=x;resolve(x)}).catch(reject);if(window.CRAFTGRN_SUBGRN_PAYLOADS[file])return finish();const s=document.createElement('script');s.onload=finish;s.onerror=()=>reject(new Error('Failed to load '+file));s.src=(SUBGRN_PAYLOAD_BASE?SUBGRN_PAYLOAD_BASE.replace(/\\/$/,'')+'/':'')+file;document.head.appendChild(s)})}",
    "function openSubgrn(d){subgrnCurrentRows=subgrnContextRows(d);const panel=subgrnEl('subgrnPanel'),sel=subgrnEl('subgrnComparisonSelect'),title=subgrnEl('subgrnTitle'),stats=subgrnEl('subgrnStats');if(!panel||!sel)return;if(!subgrnCurrentRows.length){panel.style.display='flex';title.textContent='Per Condition / Topic / Pathway Sub GRN';stats.textContent='No sub-GRN payload is available for this pathway context.';subgrnClear();return}const support=(typeof tfCountKey==='function'&&typeof tfPathwayCounts!=='undefined')?tfPathwayCounts.get(tfCountKey(d)):null,preferred=support&&support.context_id?String(support.context_id):'';sel.replaceChildren();subgrnCurrentRows.forEach(r=>{const o=document.createElement('option');o.value=r.subgrn_context_id;let label=(r.comparison_label||r.comparison_id||'comparison')+' | '+(r.pathway||'pathway')+' | genes '+(r.n_overlap_genes||0);if(support&&support.context_id&&String(r.subgrn_context_id)===String(support.context_id)&&support.byTf&&support.byTf.length){label+=' | selected TF hits '+support.byTf.map(x=>(x.label||x.tf)+':'+(x.n||0)).join(', ')}o.textContent=label;sel.appendChild(o)});sel.value=(preferred&&subgrnCurrentRows.some(r=>String(r.subgrn_context_id)===preferred))?preferred:subgrnCurrentRows[0].subgrn_context_id;panel.style.display='flex';subgrnState.view={x:0,y:0,k:1};renderSubgrn()}",
    "function selectedSubgrnContext(){const id=subgrnEl('subgrnComparisonSelect').value;return subgrnCurrentRows.find(r=>r.subgrn_context_id===id)||subgrnCurrentRows[0]||null}function subgrnClear(){['subgrnEdgeLayer','subgrnNodeLayer','subgrnLabelLayer'].forEach(id=>{const x=subgrnEl(id);if(x)x.replaceChildren()})}",
    "function subgrnNode(id,type){return{id:id,type:type,score:0,count:0,x:800,y:450}}function subgrnAddNode(map,id,type){if(!map.has(id))map.set(id,subgrnNode(id,type));return map.get(id)}",
    "function subgrnGeneSet(ctx){const s=String(ctx.overlap_genes||ctx.genes||'');return new Set(s.split(';').map(x=>x.trim()).filter(Boolean))}",
    "function subgrnTopicTheta(e,ctx){const target=String(ctx.topic_num),tokens=String(e.tf_topic_scores||'').split(';').map(x=>x.trim()).filter(Boolean);for(const tok of tokens){const p=tok.split(':'),tn=p[0],tv=subgrnNum(p[1],NaN);if(tn===target&&Number.isFinite(tv))return tv}const parts=String(e.tf_topic_nums||'').split(';').map(x=>x.trim()).filter(Boolean);if(parts.includes(target)||e.topic_tf===true)return 1;return NaN}function subgrnTopicCutoff(){return Math.max(0,Math.min(1,subgrnNum(subgrnEl('subgrnTopicThetaCutoff').value,.3)))}function subgrnTopicPhiCutoff(){const x=subgrnEl('subgrnTopicPhiCutoff');return x?Math.max(0,Math.min(1,subgrnNum(x.value,0))):0}function subgrnTopicTf(e,ctx){const theta=subgrnTopicTheta(e,ctx);if(!Number.isFinite(theta)||theta<subgrnTopicCutoff())return false;if(subgrnEl('subgrnPrimaryTopicOnly').checked){const p=Number(e.tf_primary_topic_num);return Number.isFinite(p)&&p===Number(ctx.topic_num)}return true}function subgrnTopicGene(e){const cutoff=subgrnTopicPhiCutoff();if(cutoff<=0)return true;const v=Number(e.gene_topic_score);return Number.isFinite(v)&&v>=cutoff}",
    "function subgrnRowsForContext(payload,ctx){const mode=subgrnEl('networkModeSelect').value,scope=subgrnEl('tfScopeSelect').value,topTf=subgrnLimit(subgrnEl('subgrnTopTfN'),100),topLinks=subgrnLimit(subgrnEl('subgrnTopLinkN'),300),genes=subgrnGeneSet(ctx),selected=(typeof selectedTfSet==='function')?selectedTfSet():new Set();let rows=(mode==='tf_peak_gene'?payload.tf_peak_gene_triplets:payload.tf_gene_edges)||[];rows=rows.filter(e=>{const contextMatch=e.subgrn_context_id?e.subgrn_context_id===ctx.subgrn_context_id:(String(e.comparison_id||'')===String(ctx.comparison_id||'')&&String(e.direction||'')===String(ctx.direction||'')&&genes.has(String(e.gene_key||'')));const paperTf=e.paper_tf===true||String(e.paper_tf).toLowerCase()==='true';const tfPass=scope==='topic'?subgrnTopicTf(e,ctx):(scope==='paper'?paperTf&&subgrnTopicTf(e,ctx):true);return contextMatch&&subgrnTopicGene(e)&&tfPass});const tfScore=new Map();rows.forEach(e=>{const k=String(e.tf_upper||e.tf||'');const paperTf=e.paper_tf===true||String(e.paper_tf).toLowerCase()==='true';const v=tfScore.get(k)||{score:0,genes:new Set(),paper:false};v.score+=subgrnNum(e.abs_edge_score,0);v.genes.add(String(e.gene_key||''));v.paper=v.paper||paperTf;tfScore.set(k,v)});const tfRanks=Array.from(tfScore.entries()).sort((a,b)=>(b[1].paper?1:0)-(a[1].paper?1:0)||b[1].genes.size-a[1].genes.size||b[1].score-a[1].score||a[0].localeCompare(b[0])).map(x=>x[0]);const keepTf=new Set(Number.isFinite(topTf)?tfRanks.slice(0,topTf):tfRanks);selected.forEach(k=>keepTf.add(k));rows=rows.filter(e=>keepTf.has(String(e.tf_upper||e.tf||'').toUpperCase())).sort((a,b)=>((b.paper_tf===true||String(b.paper_tf).toLowerCase()==='true')?1:0)-((a.paper_tf===true||String(a.paper_tf).toLowerCase()==='true')?1:0)||selected.has(String(b.tf_upper||b.tf||'').toUpperCase())-selected.has(String(a.tf_upper||a.tf||'').toUpperCase())||subgrnNum(b.abs_edge_score,0)-subgrnNum(a.abs_edge_score,0)||String(a.tf||'').localeCompare(String(b.tf||''))||String(a.gene_key||'').localeCompare(String(b.gene_key||'')));if(Number.isFinite(topLinks)){if(selected.size){const key=e=>[String(e.tf_upper||e.tf||''),String(e.gene_key||''),String(e.peak_id||''),String(e.subgrn_context_id||'')].join('\\t'),selRows=rows.filter(e=>selected.has(String(e.tf_upper||e.tf||'').toUpperCase())),seen=new Set(selRows.map(key)),base=rows.filter(e=>!seen.has(key(e))).slice(0,Math.max(0,topLinks-selRows.length));return selRows.concat(base).slice(0,Math.max(topLinks,selRows.length))}return rows.slice(0,topLinks)}return rows}",
    "function buildSubgrnGraph(rows,ctx){const mode=subgrnEl('networkModeSelect').value,nodes=new Map(),edges=[],overlapGenes=subgrnGeneSet(ctx);overlapGenes.forEach(g=>{const n=subgrnAddNode(nodes,'GENE:'+String(g),'Gene');n.label=String(g);n.in_overlap=true});rows.forEach(e=>{const tf='TF:'+String(e.tf||e.tf_upper),gene='GENE:'+String(e.gene_key||'');const paperTf=e.paper_tf===true||String(e.paper_tf).toLowerCase()==='true';const tfNode=subgrnAddNode(nodes,tf,'TF');tfNode.label=String(e.tf||e.tf_upper);tfNode.paper_tf=tfNode.paper_tf||paperTf;const geneNode=subgrnAddNode(nodes,gene,'Gene');geneNode.label=String(e.gene_key||'');geneNode.in_overlap=true;if(mode==='tf_peak_gene'){const peak='PEAK:'+String(e.peak_id||'peak');subgrnAddNode(nodes,peak,'Peak').label=String(e.peak_id||'peak');edges.push({from:tf,to:peak,score:subgrnNum(e.abs_edge_score,0),paper_tf:paperTf,title:(e.tf||'TF')+' -> '+(e.peak_id||'peak')+'\\nedge score: '+subgrnNum(e.abs_edge_score,0).toFixed(3)});edges.push({from:peak,to:gene,score:subgrnNum(e.abs_edge_score,0),paper_tf:paperTf,title:(e.peak_id||'peak')+' -> '+(e.gene_key||'gene')+'\\nedge score: '+subgrnNum(e.abs_edge_score,0).toFixed(3)})}else{edges.push({from:tf,to:gene,score:subgrnNum(e.abs_edge_score,0),paper_tf:paperTf,title:(e.tf||'TF')+' -> '+(e.gene_key||'gene')+'\\nsupporting peaks: '+(e.n_supporting_peaks||0)+'\\nbest peak: '+(e.best_peak_id||'NA')+'\\nedge score: '+subgrnNum(e.abs_edge_score,0).toFixed(3)})}});edges.forEach(e=>{const a=nodes.get(e.from),b=nodes.get(e.to);if(a){a.count++;a.score+=e.score}if(b){b.count++;b.score+=e.score}});return{nodes:Array.from(nodes.values()),edges:edges}}",
    "function subgrnByType(ns,type){return ns.filter(n=>n.type===type).sort((a,b)=>b.count-a.count||String(a.label||a.id).localeCompare(String(b.label||b.id)))}function subgrnSetColumn(items,x,top,bottom){const desired=(items.some(n=>n.type==='TF')?86:items.some(n=>n.type==='Peak')?42:64)*subgrnSpacing(),gap=Math.min(desired,(bottom-top)/Math.max(items.length-1,1));items.forEach((n,i)=>{n.x=x;n.y=(top+bottom)/2-gap*(items.length-1)/2+i*gap})}function subgrnPlaceRing(items,ring,a0){items.forEach((n,i)=>{const a=a0+2*Math.PI*i/Math.max(items.length,1);n.x=SUBGRN_W/2+Math.cos(a)*ring*subgrnSpacing();n.y=SUBGRN_H/2+Math.sin(a)*ring*subgrnSpacing()})}function subgrnColumns(ns){const mode=subgrnEl('networkModeSelect').value;if(mode==='tf_peak_gene'){subgrnSetColumn(subgrnByType(ns,'TF'),260,72,828);subgrnSetColumn(subgrnByType(ns,'Peak'),800,72,828);subgrnSetColumn(subgrnByType(ns,'Gene'),1340,72,828)}else{subgrnSetColumn(subgrnByType(ns,'TF'),420,72,828);subgrnSetColumn(subgrnByType(ns,'Gene'),1180,72,828)}}function subgrnRadial(ns){subgrnPlaceRing(subgrnByType(ns,'TF'),210,-Math.PI/2);subgrnPlaceRing(subgrnByType(ns,'Peak'),300,0);subgrnPlaceRing(subgrnByType(ns,'Gene'),385,-Math.PI/2)}function subgrnCircle(ns){subgrnPlaceRing(ns.slice().sort((a,b)=>String(a.label||a.id).localeCompare(String(b.label||b.id))),340,-Math.PI/2)}function subgrnGrid(ns){const cols=Math.ceil(Math.sqrt(Math.max(1,ns.length))),gapX=1180/Math.max(1,cols-1),rows=Math.ceil(ns.length/cols),gapY=660/Math.max(1,rows-1);ns.slice().sort((a,b)=>a.type.localeCompare(b.type)||String(a.label||a.id).localeCompare(String(b.label||b.id))).forEach((n,i)=>{n.x=210+(i%cols)*gapX;n.y=120+Math.floor(i/cols)*gapY})}function subgrnSpiral(ns){const arr=ns.slice().sort((a,b)=>b.count-a.count||String(a.label||a.id).localeCompare(String(b.label||b.id)));arr.forEach((n,i)=>{const a=i*.72,r=(24+i*8)*subgrnSpacing();n.x=SUBGRN_W/2+Math.cos(a)*r;n.y=SUBGRN_H/2+Math.sin(a)*r})}function subgrnClustered(ns){const centers={TF:[420,450],Peak:[800,450],Gene:[1180,450]};['TF','Peak','Gene'].forEach(type=>{const arr=subgrnByType(ns,type),c=centers[type];arr.forEach((n,i)=>{const a=i*2.399963,r=(40+Math.sqrt(i)*42)*subgrnSpacing();n.x=c[0]+Math.cos(a)*r;n.y=c[1]+Math.sin(a)*r})})}",
    "function subgrnForce(ns,es){subgrnRadial(ns);const by=new Map(ns.map(n=>[n.id,n]));for(let iter=0;iter<180;iter++){ns.forEach(n=>{n.vx=n.vx||0;n.vy=n.vy||0});for(let i=0;i<ns.length;i++)for(let j=i+1;j<ns.length;j++){const a=ns[i],b=ns[j];let dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(dx*dx+dy*dy,64),d=Math.sqrt(d2),f=Math.min(3,2500*subgrnSpacing()*subgrnSpacing()/d2);dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f}es.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;let dx=b.x-a.x,dy=b.y-a.y,d=Math.max(Math.sqrt(dx*dx+dy*dy),1),target=120*subgrnSpacing(),f=(d-target)*.012;dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f});ns.forEach(n=>{n.vx+=(SUBGRN_W/2-n.x)*.002;n.vy+=(SUBGRN_H/2-n.y)*.002;n.x=Math.max(35,Math.min(SUBGRN_W-35,n.x+n.vx));n.y=Math.max(35,Math.min(SUBGRN_H-35,n.y+n.vy));n.vx*=.76;n.vy*=.76})}ns.forEach(n=>{delete n.vx;delete n.vy})}function layoutSubgrn(ns,es){const m=subgrnEl('subgrnLayoutSelect').value;if(m==='columns'||m==='bipartite'||m==='hierarchy')subgrnColumns(ns);else if(m==='radial'||m==='concentric')subgrnRadial(ns);else if(m==='circle')subgrnCircle(ns);else if(m==='grid')subgrnGrid(ns);else if(m==='spiral')subgrnSpiral(ns);else if(m==='clustered')subgrnClustered(ns);else subgrnForce(ns,es)}",
    "function subgrnRadiusScale(ns){const vals=ns.filter(n=>n.type==='TF').map(n=>subgrnNum(n.count,1)),lo=Math.min(...vals,1),hi=Math.max(...vals,1),b0=subgrnNum(subgrnEl('subgrnTfBoxMin').value,14),b1=subgrnNum(subgrnEl('subgrnTfBoxMax').value,20);return n=>{if(n.type==='Gene')return 6;if(n.type==='Peak')return 5.5;if(hi<=lo)return (b0+b1)/2;const t=(subgrnNum(n.count,1)-lo)/(hi-lo);return b0+Math.pow(Math.max(0,Math.min(1,t)),.85)*(b1-b0)}}function subgrnLineEnd(a,b){const dx=b.x-a.x,dy=b.y-a.y,len=Math.max(Math.sqrt(dx*dx+dy*dy),1),r=b.type==='TF'?b.r*1.2:b.r+2;return{x:b.x-dx/len*r,y:b.y-dy/len*r}}function subgrnLabelPos(n){if(n.type==='Gene')return{x:n.x+n.r+7,y:n.y,anchor:'start'};if(n.type==='Peak')return{x:n.x+8,y:n.y+4,anchor:'start'};return{x:n.x,y:n.y+4,anchor:'middle'}}function subgrnApplyView(){subgrnEl('subgrnViewLayer').setAttribute('transform','translate('+subgrnState.view.x+' '+subgrnState.view.y+') scale('+subgrnState.view.k+')')}function subgrnScreenPoint(ev){const svg=subgrnEl('subgrnSvg'),p=svg.createSVGPoint();p.x=ev.clientX;p.y=ev.clientY;return p.matrixTransform(svg.getScreenCTM().inverse())}function subgrnScreenToWorld(ev){const p=subgrnScreenPoint(ev);return{x:(p.x-subgrnState.view.x)/subgrnState.view.k,y:(p.y-subgrnState.view.y)/subgrnState.view.k}}",
    "function renderSubgrn(){const ctx=selectedSubgrnContext(),stats=subgrnEl('subgrnStats'),title=subgrnEl('subgrnTitle');if(!ctx){stats.textContent='No sub-GRN context selected.';subgrnClear();return}loadSubgrnPayload(ctx.payload_file).then(payload=>{if(!payload){stats.textContent='Sub-GRN payload could not be loaded.';return}title.textContent='Per Condition / Topic / Pathway Sub GRN: '+(ctx.comparison_label||ctx.comparison_id)+' / Topic '+ctx.topic_num+' / '+ctx.pathway;const rows=subgrnRowsForContext(payload,ctx),g=buildSubgrnGraph(rows,ctx);drawSubgrn(g,ctx,rows)}).catch(e=>{stats.textContent=e.message})}",
    "function drawSubgrn(g,ctx,rows){subgrnState.nodes=g.nodes;subgrnState.edges=g.edges;const rFor=subgrnRadiusScale(g.nodes);g.nodes.forEach(n=>{n.r=rFor(n)});layoutSubgrn(g.nodes,g.edges);const edgeLayer=subgrnEl('subgrnEdgeLayer'),nodeLayer=subgrnEl('subgrnNodeLayer'),labelLayer=subgrnEl('subgrnLabelLayer'),nodeSelect=subgrnEl('subgrnNodeSelect'),showLabels=subgrnEl('subgrnShowLabels').checked,stats=subgrnEl('subgrnStats'),by=new Map(g.nodes.map(n=>[n.id,n])),scores=g.edges.map(e=>subgrnNum(e.score,0)),lo=Math.min(...scores,0),hi=Math.max(...scores,1);edgeLayer.replaceChildren();nodeLayer.replaceChildren();labelLayer.replaceChildren();nodeSelect.replaceChildren();g.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const end=subgrnLineEnd(a,b),line=el('line',{class:'subgrnEdge',x1:a.x,y1:a.y,x2:end.x,y2:end.y,stroke:subgrnEdgeColor(subgrnNum(e.score,0),lo,hi),'data-from':e.from,'data-to':e.to});if(subgrnEl('subgrnShowArrows').checked)line.setAttribute('marker-end','url(#subgrnArrow)');line.dataset.title=e.title;edgeLayer.appendChild(line)});g.nodes.forEach(n=>{const colors=subgrnColors(),paperColor=(subgrnEl('paperTfColor')&&subgrnEl('paperTfColor').value)||colors[1],otherColor=(subgrnEl('otherTfColor')&&subgrnEl('otherTfColor').value)||colors[0];let shape;if(n.type==='TF')shape=el('rect',{class:'subgrnNode subgrnNodeTf',x:n.x-n.r*1.35,y:n.y-n.r*.78,width:n.r*2.7,height:n.r*1.56,rx:3,fill:n.paper_tf?paperColor:otherColor});else shape=el('circle',{class:'subgrnNode '+(n.type==='Peak'?'subgrnNodePeak':'subgrnNodeGene'),cx:n.x,cy:n.y,r:n.r,fill:n.type==='Peak'?colors[3]:colors[2]});shape.dataset.id=n.id;shape.dataset.title=(n.label||n.id)+'\\n'+(n.paper_tf?'Paper TF':n.type)+'\\nedges: '+n.count;shape.classList.toggle('selected',n.id===subgrnState.selected);nodeLayer.appendChild(shape);if(showLabels){const p=subgrnLabelPos(n),t=el('text',{class:'subgrnLabel '+(n.type==='TF'?'subgrnTfLabel':''),x:p.x,y:p.y,'text-anchor':p.anchor,'font-size':n.type==='TF'?Math.max(9,Math.min(19,n.r*.48)):9});t.dataset.id=n.id;t.textContent=n.label||n.id;labelLayer.appendChild(t)}});g.nodes.slice().sort((a,b)=>String(a.label||a.id).localeCompare(String(b.label||b.id),undefined,{sensitivity:'base'})||String(a.id).localeCompare(String(b.id))).forEach(n=>{const o=document.createElement('option');o.value=n.id;o.textContent=n.label||n.id;nodeSelect.appendChild(o)});if(subgrnState.selected&&Array.from(nodeSelect.options).some(o=>o.value===subgrnState.selected))nodeSelect.value=subgrnState.selected;const scope=subgrnEl('tfScopeSelect').value,primary=subgrnEl('subgrnPrimaryTopicOnly').checked?' primary-only':'';stats.textContent='Showing '+g.nodes.length+' nodes and '+g.edges.length+' drawn edges from '+rows.length+' selected payload row(s). Overlap genes: '+(ctx.n_overlap_genes||0)+'. TF scope: '+scope+(scope==='topic'?' theta >= '+subgrnTopicCutoff().toFixed(2)+primary:'')+'. Target phi >= '+subgrnTopicPhiCutoff().toFixed(2)+'. Mode: '+subgrnEl('networkModeSelect').value+'. Layout: '+subgrnEl('subgrnLayoutSelect').value+'.';subgrnApplyView()}",
    "function redrawSubgrnPositions(){const by=new Map(subgrnState.nodes.map(n=>[n.id,n]));subgrnEl('subgrnEdgeLayer').querySelectorAll('.subgrnEdge').forEach(line=>{const a=by.get(line.dataset.from),b=by.get(line.dataset.to);if(!a||!b)return;const end=subgrnLineEnd(a,b);line.setAttribute('x1',a.x);line.setAttribute('y1',a.y);line.setAttribute('x2',end.x);line.setAttribute('y2',end.y)});subgrnEl('subgrnNodeLayer').querySelectorAll('.subgrnNode').forEach(shape=>{const n=by.get(shape.dataset.id);if(!n)return;if(shape.tagName.toLowerCase()==='rect'){shape.setAttribute('x',n.x-n.r*1.35);shape.setAttribute('y',n.y-n.r*.78)}else{shape.setAttribute('cx',n.x);shape.setAttribute('cy',n.y)}});subgrnEl('subgrnLabelLayer').querySelectorAll('text').forEach(t=>{const n=by.get(t.dataset.id);if(!n)return;const p=subgrnLabelPos(n);t.setAttribute('x',p.x);t.setAttribute('y',p.y);t.setAttribute('text-anchor',p.anchor)})}",
    "function selectSubgrnNode(id){subgrnState.selected=id;subgrnEl('subgrnNodeLayer').querySelectorAll('.subgrnNode').forEach(c=>c.classList.toggle('selected',c.dataset.id===id))}function resetSubgrnLayout(){subgrnState.view={x:0,y:0,k:1};renderSubgrn()}function exportSubgrnSvg(){const svg=subgrnEl('subgrnSvg'),clone=svg.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('width',SUBGRN_W);clone.setAttribute('height',SUBGRN_H);const text=new XMLSerializer().serializeToString(clone),blob=new Blob([text],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='pathway_subgrn.svg';document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(a.href)}",
    "['subgrnComparisonSelect','tfScopeSelect','subgrnTopicThetaCutoff','subgrnTopicPhiCutoff','subgrnPrimaryTopicOnly','networkModeSelect','subgrnTopTfN','subgrnTopLinkN','subgrnLayoutSelect','subgrnShowLabels','subgrnShowArrows','subgrnPaletteSelect','subgrnTfBoxMin','subgrnTfBoxMax','paperTfColor','otherTfColor'].forEach(id=>{const x=subgrnEl(id);if(x){x.addEventListener('change',renderSubgrn);x.addEventListener('input',renderSubgrn)}});{const p=subgrnEl('subgrnTopicThetaPreset'),c=subgrnEl('subgrnTopicThetaCutoff');if(p&&c){p.addEventListener('change',()=>{if(p.value!=='custom')c.value=p.value;renderSubgrn()});c.addEventListener('input',()=>{if(!['0.3','0.5','0.7'].includes(String(c.value)))p.value='custom'})}}{const p=subgrnEl('subgrnTopicPhiPreset'),c=subgrnEl('subgrnTopicPhiCutoff');if(p&&c){p.addEventListener('change',()=>{if(p.value!=='custom')c.value=p.value;renderSubgrn()});c.addEventListener('input',()=>{if(!['0','0.1','0.2','0.3','0.5','0.7'].includes(String(c.value)))p.value='custom'})}}{const r=subgrnEl('subgrnSpacingRange'),v=subgrnEl('subgrnSpacingValue');if(r&&v){r.addEventListener('input',()=>{v.value=r.value;renderSubgrn()});v.addEventListener('change',()=>{r.value=v.value;renderSubgrn()})}}{const x=subgrnEl('subgrnNodeSelect');if(x)x.addEventListener('change',()=>selectSubgrnNode(x.value))}{const b=subgrnEl('subgrnResetButton');if(b)b.addEventListener('click',resetSubgrnLayout)}{const b=subgrnEl('subgrnBackButton');if(b)b.addEventListener('click',()=>{subgrnEl('subgrnPanel').style.display='none'})}{const b=subgrnEl('subgrnExportButton');if(b)b.addEventListener('click',exportSubgrnSvg)}",
    "{const svg=subgrnEl('subgrnSvg'),canvas=subgrnEl('subgrnCanvas'),tooltip=subgrnEl('subgrnTooltip'),nodeLayer=subgrnEl('subgrnNodeLayer');if(svg&&canvas){svg.addEventListener('wheel',ev=>{ev.preventDefault();const p=subgrnScreenPoint(ev),old=subgrnState.view.k,next=Math.max(.25,Math.min(6,old*(ev.deltaY<0?1.15:1/1.15)));subgrnState.view.x=p.x-(p.x-subgrnState.view.x)*next/old;subgrnState.view.y=p.y-(p.y-subgrnState.view.y)*next/old;subgrnState.view.k=next;subgrnApplyView()},{passive:false});nodeLayer.addEventListener('mousedown',ev=>{if(!(ev.target.classList&&ev.target.classList.contains('subgrnNode')))return;ev.preventDefault();ev.stopPropagation();const n=subgrnState.nodes.find(x=>x.id===ev.target.dataset.id);if(n){subgrnState.drag=n;subgrnEl('subgrnNodeSelect').value=n.id;selectSubgrnNode(n.id)}});svg.addEventListener('mousedown',ev=>{if(ev.target.classList&&ev.target.classList.contains('subgrnNode'))return;const p=subgrnScreenPoint(ev);subgrnState.pan={startX:p.x,startY:p.y,x:subgrnState.view.x,y:subgrnState.view.y};canvas.classList.add('panning')});svg.addEventListener('mousemove',ev=>{if(subgrnState.drag){const p=subgrnScreenToWorld(ev);subgrnState.drag.x=Math.max(30,Math.min(SUBGRN_W-30,p.x));subgrnState.drag.y=Math.max(30,Math.min(SUBGRN_H-30,p.y));redrawSubgrnPositions()}else if(subgrnState.pan){const p=subgrnScreenPoint(ev);subgrnState.view.x=subgrnState.pan.x+(p.x-subgrnState.pan.startX);subgrnState.view.y=subgrnState.pan.y+(p.y-subgrnState.pan.startY);subgrnApplyView()}if(tooltip){tooltip.style.left=(ev.offsetX+12)+'px';tooltip.style.top=(ev.offsetY+12)+'px'}});svg.addEventListener('mouseup',()=>{subgrnState.drag=null;subgrnState.pan=null;canvas.classList.remove('panning')});svg.addEventListener('mouseleave',()=>{subgrnState.drag=null;subgrnState.pan=null;canvas.classList.remove('panning');if(tooltip)tooltip.style.display='none'});svg.addEventListener('mouseover',ev=>{const title=ev.target.dataset?ev.target.dataset.title:'';if(!title||!tooltip)return;tooltip.innerHTML=esc(title).replace(/\\n/g,'<br/>');tooltip.style.display='block'});svg.addEventListener('mouseout',ev=>{if(ev.target.dataset&&ev.target.dataset.title&&tooltip)tooltip.style.display='none'})}}"
  )
}
