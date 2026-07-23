# Condition-pair report data and browser helpers for Module 3.

.m3cr_empty_payload_spec <- function() {
  list(
    payload_file = "",
    network_payload_file = "",
    condition_payload_files = list(),
    payload_base = "",
    conditions = character(),
    n_tf_gene = 0L,
    n_tf_peak_gene = 0L
  )
}

.m3cr_write_compressed_payload_file <- function(obj, payload_dir, payload_file) {
  packed <- lapply(obj, .module2_report_browser_browser_payload_to_columnar)
  payload_b64 <- .module2_report_browser_encode_browser_json_deflate_base64(packed)
  js <- paste0(
    "window.CRAFTGRN_CONDITION_PAYLOADS=window.CRAFTGRN_CONDITION_PAYLOADS||{};\n",
    "window.CRAFTGRN_CONDITION_PAYLOADS[",
    jsonlite::toJSON(payload_file, auto_unbox = TRUE),
    "]={compressed_columnar:'", payload_b64, "'};\n"
  )
  writeLines(js, file.path(payload_dir, payload_file), useBytes = TRUE)
  invisible(payload_file)
}

.m3cr_payload_stem <- function(payload_name) {
  key <- digest::digest(
    as.character(payload_name)[[1L]],
    algo = "xxhash64",
    serialize = FALSE
  )
  paste0("cp_", substr(key, 1L, 12L))
}

.m3cr_report_data_file <- function(out_html) {
  key <- digest::digest(
    basename(as.character(out_html)[[1L]]),
    algo = "xxhash64",
    serialize = FALSE
  )
  paste0("cr_", substr(key, 1L, 12L), "_d.js")
}

.m3cr_condition_link_manifest <- function(condition_links_dir) {
  if (is.na(condition_links_dir) || !dir.exists(condition_links_dir)) {
    return(data.table::data.table())
  }
  path <- file.path(condition_links_dir, "condition_links_manifest.csv")
  if (file.exists(path)) {
    out <- data.table::fread(path, showProgress = FALSE)
    if (!all(c("condition_id", "path") %in% names(out))) {
      return(data.table::data.table())
    }
    out[, path := as.character(path)]
    out[!grepl("^/", path), path := file.path(condition_links_dir, path)]
    return(unique(out[file.exists(path), .(
      condition_id = as.character(condition_id),
      path
    )]))
  }
  files <- sort(list.files(
    condition_links_dir,
    "_condition_links[.]parquet$",
    full.names = TRUE
  ))
  data.table::data.table(
    condition_id = sub("_condition_links[.]parquet$", "", basename(files)),
    path = files
  )
}

.m3cr_gene_topic_rows <- function(extraction_dir, topic_space = "combined") {
  path <- .m3tb_topic_space_file(
    extraction_dir,
    topic_space = topic_space,
    kind = "terms"
  )
  if (!file.exists(path)) return(data.table::data.table())
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  cols <- intersect(
    c("topic", "topic_num", "term_id", "score", "in_topic"),
    header
  )
  if (!all(c("term_id", "score") %in% cols) ||
      !any(c("topic", "topic_num") %in% cols)) {
    return(data.table::data.table())
  }
  dt <- data.table::fread(path, select = cols, showProgress = FALSE)
  if ("in_topic" %in% names(dt)) dt <- dt[as.logical(in_topic) %in% TRUE]
  topic_value <- if ("topic_num" %in% names(dt)) {
    suppressWarnings(as.integer(dt$topic_num))
  } else {
    suppressWarnings(as.integer(sub("^Topic\\s*", "", as.character(dt$topic))))
  }
  keep <- grepl("^GENE:", dt$term_id)
  dt <- data.table::data.table(
    topic_num = topic_value[keep],
    gene_key = sub("^GENE:", "", as.character(dt$term_id[keep])),
    topic_score = suppressWarnings(as.numeric(dt$score[keep]))
  )
  unique(dt[
    is.finite(topic_num) & !is.na(gene_key) & nzchar(gene_key) &
      is.finite(topic_score)
  ])
}

.m3cr_collect_condition_edges <- function(path,
                                          topic_genes = character(),
                                          max_peaks_per_tf_gene = 1L) {
  if (!requireNamespace("arrow", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    .log_abort("Packages `arrow` and `dplyr` are required for condition-pair report payloads.")
  }
  ds <- arrow::open_dataset(path)
  required <- c(
    "condition_id", "tf", "gene_key", "peak_id", "fp_score_condition",
    "gene_expr_condition", "tf_expr_condition"
  )
  available <- names(ds$schema)
  if (!all(required %in% available)) {
    .log_abort(
      "Condition-link file is missing report columns ({paste(setdiff(required, available), collapse = ', ')}): {path}"
    )
  }
  query <- dplyr::select(ds, dplyr::all_of(required))
  agg <- tryCatch({
    query |>
      dplyr::group_by(condition_id, tf, gene_key) |>
      dplyr::summarise(
        fp_sum = sum(fp_score_condition, na.rm = TRUE),
        fp_max = max(fp_score_condition, na.rm = TRUE),
        n_peaks = dplyr::n_distinct(peak_id),
        gene_expr = max(gene_expr_condition, na.rm = TRUE),
        tf_expr = max(tf_expr_condition, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::collect() |>
      data.table::as.data.table()
  }, error = function(e) NULL)
  if (is.null(agg)) {
    raw <- data.table::as.data.table(dplyr::collect(query))
    agg <- raw[, .(
      fp_sum = sum(fp_score_condition, na.rm = TRUE),
      fp_max = max(fp_score_condition, na.rm = TRUE),
      n_peaks = data.table::uniqueN(peak_id),
      gene_expr = max(gene_expr_condition, na.rm = TRUE),
      tf_expr = max(tf_expr_condition, na.rm = TRUE)
    ), by = .(condition_id, tf, gene_key)]
  }
  agg[, `:=`(
    condition_id = as.character(condition_id),
    tf = as.character(tf),
    tf_upper = toupper(as.character(tf)),
    gene_key = as.character(gene_key),
    fp_sum = suppressWarnings(as.numeric(fp_sum)),
    fp_max = suppressWarnings(as.numeric(fp_max)),
    n_peaks = suppressWarnings(as.integer(n_peaks)),
    gene_expr = suppressWarnings(as.numeric(gene_expr)),
    tf_expr = suppressWarnings(as.numeric(tf_expr))
  )]
  peak <- data.table::data.table()
  if (length(topic_genes)) {
    peak <- tryCatch({
      query |>
        dplyr::filter(gene_key %in% topic_genes) |>
        dplyr::collect() |>
        data.table::as.data.table()
    }, error = function(e) data.table::data.table())
    if (nrow(peak)) {
      peak[, `:=`(
        condition_id = as.character(condition_id),
        tf = as.character(tf),
        tf_upper = toupper(as.character(tf)),
        gene_key = as.character(gene_key),
        peak_id = as.character(peak_id),
        fp_score = suppressWarnings(as.numeric(fp_score_condition)),
        gene_expr = suppressWarnings(as.numeric(gene_expr_condition)),
        tf_expr = suppressWarnings(as.numeric(tf_expr_condition))
      )]
      peak <- peak[, .(
        fp_score = max(fp_score, na.rm = TRUE),
        gene_expr = max(gene_expr, na.rm = TRUE),
        tf_expr = max(tf_expr, na.rm = TRUE)
      ), by = .(condition_id, tf, tf_upper, gene_key, peak_id)]
      peak[!is.finite(fp_score), fp_score := 0]
      peak[!is.finite(gene_expr), gene_expr := 0]
      peak[!is.finite(tf_expr), tf_expr := 0]
      data.table::setorder(
        peak,
        condition_id, tf_upper, gene_key, -fp_score, peak_id
      )
      peak[, peak_rank := seq_len(.N), by = .(
        condition_id, tf_upper, gene_key
      )]
      peak <- peak[
        peak_rank <= max(1L, as.integer(max_peaks_per_tf_gene)[[1L]])
      ][, peak_rank := NULL]
    }
  }
  list(tf_gene = agg, tf_peak_gene = peak)
}

.m3cr_expression_audit_rows <- function(condition_links_dir) {
  paths <- c(
    file.path(condition_links_dir, "condition_gene_expression.csv"),
    file.path(condition_links_dir, "condition_gene_log2fc.csv"),
    file.path(condition_links_dir, "condition_comparison_gene_log2fc.csv")
  )
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return(data.table::data.table())
  if (basename(paths[[1L]]) == "condition_gene_expression.csv") {
    paths <- paths[[1L]]
  }
  rows <- lapply(paths, function(path) {
    dt <- data.table::fread(path, showProgress = FALSE)
    if (all(c("condition_id", "gene_key", "expression") %in% names(dt))) {
      return(dt[, .(
        condition_id = as.character(condition_id),
        gene_key = as.character(gene_key),
        gene_expr = suppressWarnings(as.numeric(expression))
      )])
    }
    if (!all(c("gene_key", "condition_1", "condition_2") %in% names(dt))) {
      return(NULL)
    }
    e1 <- intersect(c("expression_condition_1", "expression_1"), names(dt))[1L]
    e2 <- intersect(c("expression_condition_2", "expression_2"), names(dt))[1L]
    if (is.na(e1) || is.na(e2)) return(NULL)
    data.table::rbindlist(list(
      dt[, .(
        condition_id = as.character(condition_1),
        gene_key = as.character(gene_key),
        gene_expr = suppressWarnings(as.numeric(get(e1)))
      )],
      dt[, .(
        condition_id = as.character(condition_2),
        gene_key = as.character(gene_key),
        gene_expr = suppressWarnings(as.numeric(get(e2)))
      )]
    ))
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(out)
  out <- out[is.finite(gene_expr), .(gene_expr = max(gene_expr)), by = .(
    condition_id, gene_key
  )]
  out
}

.m3cr_normalize_expression_gene_keys <- function(expression) {
  expression <- data.table::copy(data.table::as.data.table(expression))
  if (!nrow(expression) || !"gene_key" %in% names(expression)) return(expression)
  keys <- unique(as.character(expression$gene_key))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  ensembl <- grepl("^ENS(MUS)?G[0-9]+([.][0-9]+)?$", keys, ignore.case = TRUE)
  if (!any(ensembl)) return(expression)

  ids <- keys[ensembl]
  bare_ids <- sub("[.][0-9]+$", "", ids)
  species <- ifelse(grepl("^ENSMUSG", bare_ids, ignore.case = TRUE), "mouse", "human")
  mapping <- data.table::rbindlist(lapply(unique(species), function(sp) {
    idx <- species == sp
    resolved <- resolve_gene_symbols(
      bare_ids[idx],
      species = sp,
      use_annotation = TRUE,
      allow_heuristic = FALSE
    )
    data.table::data.table(
      gene_key = ids[idx],
      canonical = ifelse(
        resolved$matched & !is.na(resolved$canonical_symbol),
        resolved$canonical_symbol,
        ids[idx]
      )
    )
  }), use.names = TRUE)
  idx <- match(expression$gene_key, mapping$gene_key)
  replace <- !is.na(idx)
  expression$gene_key[replace] <- mapping$canonical[idx[replace]]
  expression
}

.m3cr_condition_topic_expression_share <- function(gene_expr, gene_topics) {
  expression <- data.table::copy(data.table::as.data.table(gene_expr))
  topics <- data.table::copy(data.table::as.data.table(gene_topics))
  empty <- data.table::data.table(
    condition_id = character(), topic_num = integer(),
    expression_topic_share = numeric(), expression_mass = numeric(),
    n_expression_genes = integer()
  )
  if (!all(c("condition_id", "gene_key", "gene_expr") %in% names(expression)) ||
      !all(c("gene_key", "topic_num") %in% names(topics))) {
    return(empty)
  }
  expression <- .m3cr_normalize_expression_gene_keys(expression)
  expression[, `:=`(
    condition_id = as.character(condition_id),
    gene_match = toupper(as.character(gene_key)),
    gene_expr = suppressWarnings(as.numeric(gene_expr))
  )]
  expression <- expression[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(gene_match) & nzchar(gene_match) & is.finite(gene_expr),
    .(gene_expr = max(gene_expr, na.rm = TRUE)),
    by = .(condition_id, gene_match)
  ]
  if (!nrow(expression)) return(empty)

  if (!"topic_score" %in% names(topics)) topics[, topic_score := 1]
  topics[, `:=`(
    gene_match = toupper(as.character(gene_key)),
    topic_num = suppressWarnings(as.integer(topic_num)),
    topic_score = suppressWarnings(as.numeric(topic_score))
  )]
  topics[!is.finite(topic_score), topic_score := 0]
  topics <- topics[
    !is.na(gene_match) & nzchar(gene_match) & is.finite(topic_num)
  ]
  if (!nrow(topics)) return(empty)
  data.table::setorder(topics, gene_match, -topic_score, topic_num)
  topics <- topics[, .SD[1L], by = gene_match]

  values <- topics[
    expression,
    on = "gene_match",
    nomatch = 0L,
    allow.cartesian = FALSE
  ]
  if (!nrow(values)) return(empty)
  values[, expression_value := log2(pmax(gene_expr, 0) + 1)]
  out <- values[, .(
    expression_mass = sum(expression_value, na.rm = TRUE),
    n_expression_genes = data.table::uniqueN(gene_match)
  ), by = .(condition_id, topic_num)]
  out[, expression_topic_share := {
    total <- sum(expression_mass, na.rm = TRUE)
    if (is.finite(total) && total > 0) expression_mass / total else 0
  }, by = condition_id]
  data.table::setcolorder(out, names(empty))
  data.table::setorder(out, condition_id, topic_num)
  out[]
}

.m3cr_pathway_expression_scores <- function(topic_pathways,
                                             gene_expr,
                                             gene_topics = NULL,
                                             min_genes = 3L) {
  pathways <- data.table::copy(data.table::as.data.table(topic_pathways))
  expression <- data.table::copy(data.table::as.data.table(gene_expr))
  empty_activity <- data.table::data.table(
    topic_num = integer(), pathway_key = character(), condition_id = character(),
    rank_enrichment = numeric(), mean_gene_zscore = numeric(),
    mean_log2_expression = numeric(), n_expression_genes = integer(),
    n_expression_universe = integer()
  )
  empty_expression <- data.table::data.table(
    condition_id = character(), gene_key = character(), gene_expr = numeric(),
    log2_expression = numeric(), expression_rank = numeric(),
    gene_zscore = numeric(), n_expression_universe = integer()
  )
  if (!nrow(pathways) || !nrow(expression)) {
    return(list(activity = empty_activity, gene_expression = empty_expression))
  }
  if (!"topic_num" %in% names(pathways) && "topic" %in% names(pathways)) {
    pathways[, topic_num := suppressWarnings(as.integer(topic))]
  }
  if (!"pathway_key" %in% names(pathways)) pathways[, pathway_key := pathway]
  if (!"genes" %in% names(pathways) && "overlap_genes" %in% names(pathways)) {
    pathways[, genes := overlap_genes]
  }
  .assert_has_cols(
    pathways,
    c("topic_num", "pathway_key", "genes"),
    context = "condition pathway expression scoring"
  )
  .assert_has_cols(
    expression,
    c("condition_id", "gene_key", "gene_expr"),
    context = "condition pathway expression scoring"
  )
  pathways <- unique(pathways[
    is.finite(suppressWarnings(as.integer(topic_num))) &
      !is.na(pathway_key) & nzchar(as.character(pathway_key)) &
      !is.na(genes) & nzchar(as.character(genes)),
    .(
      topic_num = suppressWarnings(as.integer(topic_num)),
      pathway_key = as.character(pathway_key),
      genes = as.character(genes)
    )
  ])
  if (!nrow(pathways)) {
    return(list(activity = empty_activity, gene_expression = empty_expression))
  }
  gene_sets <- .m3tb_subgrn_split_genes(pathways$genes)
  membership <- data.table::rbindlist(lapply(seq_len(nrow(pathways)), function(i) {
    data.table::data.table(
      topic_num = pathways$topic_num[[i]],
      pathway_key = pathways$pathway_key[[i]],
      gene_key = gene_sets[[i]]
    )
  }), use.names = TRUE, fill = TRUE)
  membership <- unique(membership[
    !is.na(gene_key) & nzchar(as.character(gene_key))
  ])
  if (!nrow(membership)) {
    return(list(activity = empty_activity, gene_expression = empty_expression))
  }
  membership[, gene_match := toupper(as.character(gene_key))]
  expression <- .m3cr_normalize_expression_gene_keys(expression)
  expression[, `:=`(
    condition_id = as.character(condition_id),
    gene_key = as.character(gene_key),
    gene_match = toupper(as.character(gene_key)),
    gene_expr = suppressWarnings(as.numeric(gene_expr))
  )]
  expression <- expression[
    is.finite(gene_expr),
    .(gene_expr = max(gene_expr, na.rm = TRUE)),
    by = .(condition_id, gene_key, gene_match)
  ]
  if (!nrow(expression)) {
    return(list(activity = empty_activity, gene_expression = empty_expression))
  }
  expression[, log2_expression := log2(pmax(gene_expr, 0) + 1)]
  expression[, c("expression_rank", "n_expression_universe") := {
    n <- .N
    list(
      data.table::frank(log2_expression, ties.method = "average"),
      rep(as.integer(n), n)
    )
  }, by = condition_id]
  expression[, gene_zscore := {
    s <- stats::sd(log2_expression)
    if (!is.finite(s) || s <= 0) rep(0, .N) else
      (log2_expression - mean(log2_expression)) / s
  }, by = gene_key]
  scored <- expression[membership, on = "gene_match", nomatch = 0L, allow.cartesian = TRUE]
  if (!nrow(scored)) {
    return(list(activity = empty_activity, gene_expression = empty_expression))
  }
  min_genes <- max(1L, suppressWarnings(as.integer(min_genes)[[1L]]))
  activity <- scored[, .(
    rank_enrichment = {
      m <- data.table::uniqueN(gene_key)
      n <- as.integer(n_expression_universe[[1L]])
      expected <- m * (n + 1) / 2
      variance <- m * (n - m) * (n + 1) / 12
      if (is.finite(variance) && variance > 0) {
        (sum(expression_rank, na.rm = TRUE) - expected) / sqrt(variance)
      } else {
        0
      }
    },
    mean_gene_zscore = mean(gene_zscore, na.rm = TRUE),
    mean_log2_expression = mean(log2_expression, na.rm = TRUE),
    n_expression_genes = data.table::uniqueN(gene_key),
    n_expression_universe = as.integer(n_expression_universe[[1L]])
  ), by = .(topic_num, pathway_key, condition_id)]
  activity[n_expression_genes < min_genes, `:=`(
    rank_enrichment = NA_real_,
    mean_gene_zscore = NA_real_,
    mean_log2_expression = NA_real_
  )]
  pathway_genes <- unique(membership$gene_match)
  gene_expression <- expression[
    gene_match %in% pathway_genes,
    .(
      condition_id, gene_key, gene_expr, log2_expression,
      expression_rank, gene_zscore, n_expression_universe
    )
  ]
  data.table::setorder(activity, topic_num, pathway_key, condition_id)
  data.table::setorder(gene_expression, condition_id, gene_key)
  list(activity = activity, gene_expression = gene_expression)
}

.m3cr_write_condition_payload <- function(extraction_dir,
                                           model_dir,
                                           payload_dir,
                                           payload_name,
                                           payload_base,
                                           topic_space = "combined",
                                           topic_pathways = NULL,
                                           condition_pathways = NULL) {
  condition_links_dir <- .m3tb_subgrn_condition_links_dir(
    model_dir = model_dir,
    extraction_dir = extraction_dir
  )
  manifest <- .m3cr_condition_link_manifest(condition_links_dir)
  if (!nrow(manifest)) return(.m3cr_empty_payload_spec())
  gene_topics <- .m3cr_gene_topic_rows(extraction_dir, topic_space = topic_space)
  topic_genes <- unique(gene_topics$gene_key)
  max_peaks_per_tf_gene <- suppressWarnings(as.integer(Sys.getenv(
    "CRAFTGRN_CONDITION_REPORT_MAX_PEAKS_PER_TF_GENE",
    unset = "1"
  )))
  if (!is.finite(max_peaks_per_tf_gene) || max_peaks_per_tf_gene < 1L) {
    max_peaks_per_tf_gene <- 1L
  }
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    .m3cr_collect_condition_edges(
      manifest$path[[i]],
      topic_genes,
      max_peaks_per_tf_gene = max_peaks_per_tf_gene
    )
  })
  edge_gene_expr <- data.table::rbindlist(lapply(parts, function(part) {
    part$tf_gene[, .(
      gene_expr = max(gene_expr, na.rm = TRUE)
    ), by = .(condition_id, gene_key)]
  }), use.names = TRUE, fill = TRUE)
  gene_expr <- data.table::rbindlist(
    list(.m3cr_expression_audit_rows(condition_links_dir), edge_gene_expr),
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(gene_expr)) {
    gene_expr <- gene_expr[
      is.finite(gene_expr),
      .(gene_expr = max(gene_expr)),
      by = .(condition_id, gene_key)
    ]
  }
  pathway_scores <- .m3cr_pathway_expression_scores(
    topic_pathways = topic_pathways,
    gene_expr = gene_expr,
    gene_topics = gene_topics,
    min_genes = 3L
  )
  condition_pathways <- data.table::as.data.table(condition_pathways)
  if (nrow(condition_pathways)) {
    if (!"condition_id" %in% names(condition_pathways) &&
        "comparison_id" %in% names(condition_pathways)) {
      condition_pathways[, condition_id := as.character(comparison_id)]
    }
    if (!"genes" %in% names(condition_pathways)) {
      condition_pathways[, genes := NA_character_]
    }
    if ("overlap_genes" %in% names(condition_pathways)) {
      condition_pathways[
        is.na(genes) | !nzchar(genes),
        genes := as.character(overlap_genes)
      ]
    }
    condition_pathways <- condition_pathways[
      !is.na(condition_id) & nzchar(condition_id) &
        is.finite(topic_num) & is.finite(padj) & padj <= 0.05,
      .(
        condition_id = as.character(condition_id),
        topic_num = as.integer(topic_num),
        pathway = as.character(pathway),
        pathway_key = as.character(pathway_key),
        padj = as.numeric(padj),
        combined_score = as.numeric(combined_score),
        gene_in = as.integer(gene_in),
        gene_total_universe = as.integer(gene_total_universe),
        genes = as.character(genes)
      )
    ]
  } else {
    condition_pathways <- data.table::data.table(
      condition_id = character(),
      topic_num = integer(),
      pathway = character(),
      pathway_key = character(),
      padj = numeric(),
      combined_score = numeric(),
      gene_in = integer(),
      gene_total_universe = integer(),
      genes = character()
    )
  }
  condition_topic_expression <- .m3cr_condition_topic_expression_share(
    gene_expr = gene_expr,
    gene_topics = gene_topics
  )
  score_suffix <- if (identical(topic_space, "raw")) "_raw" else "_combined"
  score_file <- file.path(
    extraction_dir,
    paste0("topic_condition_pathway_expression_scores", score_suffix, ".csv")
  )
  data.table::fwrite(pathway_scores$activity, score_file)
  if (identical(topic_space, "combined")) {
    data.table::fwrite(
      pathway_scores$activity,
      file.path(extraction_dir, "topic_condition_pathway_expression_scores.csv")
    )
  }
  dir.create(payload_dir, recursive = TRUE, showWarnings = FALSE)
  payload_stem <- .m3cr_payload_stem(payload_name)
  unlink(file.path(payload_dir, paste0(payload_name, "*.js")), force = TRUE)
  unlink(file.path(payload_dir, paste0(payload_stem, "*.js")), force = TRUE)
  payload_file <- paste0(payload_stem, "_o.js")
  activity <- data.table::rbindlist(lapply(parts, function(part) {
    part$tf_gene[, .(
      fp_sum = sum(fp_sum, na.rm = TRUE),
      n_targets = data.table::uniqueN(gene_key),
      tf_expr = max(tf_expr, na.rm = TRUE)
    ), by = .(condition_id, tf, tf_upper)]
  }), use.names = TRUE, fill = TRUE)
  .m3cr_write_compressed_payload_file(
    list(
      tf_activity = activity,
      gene_topics = gene_topics,
      condition_topic_expression = condition_topic_expression,
      pathway_activity = pathway_scores$activity,
      pathway_gene_expression = pathway_scores$gene_expression
    ),
    payload_dir,
    payload_file
  )
  condition_files <- vector("list", nrow(manifest))
  names(condition_files) <- as.character(manifest$condition_id)
  for (i in seq_len(nrow(manifest))) {
    condition_key <- as.character(manifest$condition_id[[i]])
    prefix <- sprintf("%s_c%02d", payload_stem, i)
    edge_file <- paste0(prefix, "_e.js")
    peak_file <- paste0(prefix, "_p.js")
    pathway_file <- paste0(prefix, "_w.js")
    .m3cr_write_compressed_payload_file(
      list(tf_gene = parts[[i]]$tf_gene),
      payload_dir,
      edge_file
    )
    .m3cr_write_compressed_payload_file(
      list(tf_peak_gene = parts[[i]]$tf_peak_gene),
      payload_dir,
      peak_file
    )
    .m3cr_write_compressed_payload_file(
      list(
        condition_pathways = condition_pathways[
          condition_id == condition_key
        ]
      ),
      payload_dir,
      pathway_file
    )
    condition_files[[condition_key]] <- list(
      edge_file = edge_file,
      peak_file = peak_file,
      pathway_file = pathway_file
    )
  }
  n_tf_gene <- sum(vapply(parts, function(x) nrow(x$tf_gene), integer(1L)))
  n_tf_peak_gene <- sum(vapply(parts, function(x) nrow(x$tf_peak_gene), integer(1L)))
  list(
    payload_file = payload_file,
    network_payload_file = "",
    condition_payload_files = condition_files,
    payload_base = payload_base,
    conditions = sort(unique(as.character(manifest$condition_id))),
    n_tf_gene = n_tf_gene,
    n_tf_peak_gene = n_tf_peak_gene
  )
}

.m3cr_payload_spec_json <- function(spec) {
  jsonlite::toJSON(spec, auto_unbox = TRUE, null = "null", na = "null")
}

.m3cr_condition_report_css <- function() {
  paste0(
    "html,body{width:100%;height:100%;overflow:hidden}*{box-sizing:border-box}body{margin:0;background:#f3f4f2;color:#171717;font-family:Arial,Helvetica,sans-serif;font-weight:700;display:flex;flex-direction:column}",
    ".top{min-height:52px;flex:0 0 auto;background:#fff;border-bottom:1px solid #cfd2cc;display:flex;flex-wrap:wrap;gap:8px 12px;align-items:center;padding:6px 12px;white-space:nowrap;overflow:hidden}.top input[type=color]{width:30px;height:28px;padding:2px;border:1px solid #9ca3a0;border-radius:4px;background:#fff;cursor:pointer}",
    ".top label{display:flex;gap:5px;align-items:center;font-size:12px}.top label+label{padding-left:10px;border-left:1px solid #e2e4df}.top select{max-width:210px}",
    "select,input,button{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #9ca3a0;border-radius:4px;background:#fff;color:#171717;padding:6px 8px}button{background:#202322;color:#fff;border-color:#202322;cursor:pointer}button:hover{background:#000}button.secondary{background:#fff;color:#171717}button.secondary:hover{background:#eceeeb}button:focus-visible,select:focus-visible,input:focus-visible,summary:focus-visible{outline:3px solid #6ca0dc;outline-offset:1px}select.conditionTarget{outline:3px solid #202322;outline-offset:1px}",
    ".grid{flex:1;min-height:0;display:grid;grid-template-columns:minmax(620px,1fr) minmax(620px,1fr);gap:8px;padding:8px;overflow:hidden}",
    ".left{display:grid;grid-template-rows:minmax(0,1fr) minmax(0,1fr);gap:8px;min-height:0}.topPair,.bottomPair{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);min-height:0;height:100%}",
    ".mini{display:flex;flex-direction:column;min-width:0;min-height:0}.mini:first-child{border-right:1px solid #dfe2dc}",
    ".pane{background:#fff;border:1px solid #cfd2cc;display:flex;flex-direction:column;min-height:0;overflow:hidden}.paneHead{min-height:42px;padding:8px 11px;border-bottom:1px solid #dfe2dc;display:flex;justify-content:space-between;gap:10px;align-items:center}",
    ".pane h2,.mini h2{font-size:14px;line-height:1.1;margin:0;white-space:nowrap}.body{position:relative;flex:1;min-height:0;overflow:hidden}.meta{font-size:11px;color:#52605a}.bottomPair .meta{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap}.bottomPair .paneHead{padding:6px 8px;gap:4px}.contextWrap,.headTools{display:flex;align-items:center;justify-content:flex-end;gap:5px;min-width:0}.conditionTopicMetricControl{width:142px;min-width:0;font-size:10px;padding:4px 5px}body.embed .conditionTopicMetricControl,body.embed .conditionTopicMetricControl+.help{display:none}.contextBadge{background:#202322;color:#fff;border-radius:4px;padding:5px 8px;font-size:11px;max-width:520px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.iconButton{width:25px;height:25px;padding:0!important;display:inline-flex;align-items:center;justify-content:center;font-size:15px!important;line-height:1}.zoomButton{height:25px;padding:3px 8px!important;font-size:11px!important;line-height:1;white-space:nowrap}.pathScoreControl{font-size:11px;gap:4px!important}.pathScoreControl select{font-size:11px;padding:4px 6px}",
    ".help{display:inline-flex;width:19px;height:19px;border:1px solid #9ca3a0;border-radius:50%;align-items:center;justify-content:center;font-size:12px;color:#38413d;background:#fff;cursor:help}.tooltip{position:absolute;display:none;background:rgba(17,17,17,.96);color:#fff;font:700 12px Arial,Helvetica,sans-serif;padding:8px 9px;border-radius:4px;pointer-events:none;max-width:460px;line-height:1.35;z-index:30}",
    "svg{width:100%;height:100%;display:block}.fixedChart{position:absolute;inset:0;overflow:hidden;background:#fff}#mdsLayer text,#activityLayer text{font-size:24px;font-weight:700}#butterflyLayer text,#tfButterflyLayer text{font-size:21px;font-weight:700}.pager{display:inline-flex;align-items:center;gap:4px;white-space:nowrap}.pager button{width:25px;height:25px;padding:0;display:inline-flex;align-items:center;justify-content:center;font-size:15px;line-height:1;background:#fff;color:#171717;border-color:#9ca3a0}.pager button:hover:not(:disabled){background:#eceeeb}.pager button:disabled{opacity:.35;cursor:default}.pageStatus{min-width:42px;text-align:center;font-size:11px;color:#52605a}.pathLabel{font:700 11px Arial,Helvetica,sans-serif;fill:#171717}.pathLabelTopicSpecific{fill:#c42f3c}.pathLabelConditionSpecific{fill:#1f5fb8}.pathLabelBothSpecific{fill:#7626a8}.pathLegend{font:700 12px Arial,Helvetica,sans-serif}.pathRowSelected{fill:#f0f2ef}",
    ".networkPanel{position:absolute;inset:0;background:#fff;border:0;box-shadow:none;display:none;flex-direction:column;z-index:100}.networkTop{border-bottom:1px solid #cfd2cc;background:#f8f9f7}.networkHeading{min-height:54px;padding:7px 10px;display:flex;align-items:flex-start;justify-content:space-between;gap:12px;border-bottom:1px solid #dfe2dc}.networkHeadingText{min-width:0;flex:1}.networkTitle{font-size:15px;line-height:1.25;white-space:normal;overflow-wrap:anywhere}.networkContext{font-size:10px;line-height:1.25;color:#52605a;margin-top:3px;white-space:normal}.networkActions{display:flex;gap:6px;flex:0 0 auto}.networkActions button{padding:5px 8px}.networkTabs{height:31px;padding:3px 10px 0;display:flex;gap:3px;border-bottom:1px solid #cfd2cc}.networkTab{height:28px;padding:4px 13px;border:1px solid transparent;border-bottom:0;background:transparent;color:#26342e;border-radius:4px 4px 0 0}.networkTab:hover{background:#eceeeb}.networkTab.active{background:#fff;border-color:#cfd2cc;color:#111}.networkTabPanel{display:none;min-height:43px;padding:6px 10px;gap:6px 10px;align-items:center;flex-wrap:wrap;background:#fff}.networkTabPanel.active{display:flex}.networkControls label{display:flex;gap:4px;align-items:center;font-size:11px;white-space:nowrap}.networkControls select,.networkControls input{font-size:11px;padding:4px 6px}.networkControls input[type=number]{width:57px}.networkControls input[type=range]{width:92px;padding:0}.networkControls .wideSelect{min-width:132px;max-width:205px}.networkStats{font-size:10px;line-height:1.2;color:#46524d;padding:0 10px 5px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}.networkCanvas{position:relative;flex:1;min-height:0;background:#fff}.node{cursor:grab}.node.selected{stroke:#111;stroke-width:5}.nodeLabel{font:700 14px Arial,Helvetica,sans-serif;paint-order:stroke;stroke-width:4px;stroke-linejoin:round;pointer-events:none}.edge{fill:none;stroke-linecap:round;opacity:.58}",
    "body.network-open{overflow:hidden}body.embed .top{display:none}body.embed .grid{flex:1}",
    ".pathLabelTopicSpecific,.pathLabelConditionSpecific,.pathLabelBothSpecific{fill:#171717}.reportModeToggle{display:inline-flex;gap:2px}.reportModeToggle button{padding:4px 8px;background:#fff;color:#171717;border-color:#9ca3a0}.reportModeToggle button.active{background:#202322;color:#fff}.loadingOverlay{position:fixed;inset:0;z-index:1000;background:rgba(255,255,255,.86);display:flex;align-items:center;justify-content:center;gap:10px;font-size:13px}.loadingOverlay[hidden]{display:none}.loadingSpinner{width:26px;height:26px;border:4px solid #cfd5d0;border-top-color:#202322;border-radius:50%;animation:m3cr-spin .75s linear infinite}@keyframes m3cr-spin{to{transform:rotate(360deg)}}.networkControls{padding:5px 7px;gap:5px 8px}.networkHeading{min-height:46px;padding:5px 8px}.networkHeading button{padding:4px 7px}.singleColorInput{width:30px;height:27px;padding:2px}.singleColorInput[hidden]{display:none}",
    "@media(max-width:1450px){.contextBadge{display:none}.top{gap:7px}.top label+label{padding-left:7px}.bottomPair .help{display:none}.bottomPair .paneHead h2{font-size:13px}.paneHead:has(#shortConditionNames) .help,#mdsStats{display:none}.paneHead label:has(#shortConditionNames){font-size:10px!important}}",
    "@media(max-width:1240px){.grid{grid-template-columns:minmax(0,1fr) minmax(0,1fr)}.pathLabel{font-size:10.5px}}"
  )
}

.m3cr_network_panel_html <- function() {
  palette_options <- "<option value=\"rdbu\" selected>Red-blue</option><option value=\"prgn\">Purple-green</option><option value=\"brbg\">Brown-teal</option><option value=\"default\">Default</option><option value=\"npg\">NPG inspired</option><option value=\"aaas\">AAAS inspired</option><option value=\"nejm\">NEJM inspired</option><option value=\"lancet\">Lancet inspired</option><option value=\"jama\">JAMA inspired</option><option value=\"bmj\">BMJ inspired</option><option value=\"jco\">JCO inspired</option><option value=\"ucscgb\">UCSCGB inspired</option><option value=\"d3\">D3 inspired</option><option value=\"gephi\">Gephi inspired</option><option value=\"observable\">Observable inspired</option><option value=\"primer\">Primer inspired</option><option value=\"atlassian\">Atlassian inspired</option><option value=\"iterm\">iTerm inspired</option><option value=\"locuszoom\">LocusZoom inspired</option><option value=\"igv\">IGV inspired</option><option value=\"cosmic\">COSMIC inspired</option><option value=\"uchicago\">UChicago inspired</option><option value=\"startrek\">Star Trek inspired</option><option value=\"tron\">Tron inspired</option><option value=\"futurama\">Futurama inspired</option><option value=\"rickandmorty\">Rick and Morty inspired</option><option value=\"simpsons\">Simpsons inspired</option><option value=\"flatui\">Flat UI inspired</option><option value=\"frontiers\">Frontiers inspired</option><option value=\"gsea\">GSEA inspired</option><option value=\"bootstrap5\">Bootstrap 5 inspired</option><option value=\"material\">Material Design inspired</option><option value=\"tailwind\">Tailwind CSS inspired</option>"
  palette_options <- sub(
    ">Red-blue<",
    ">Red-blue (Default)<",
    palette_options,
    fixed = TRUE
  )
  palette_options <- sub(
    "<option value=\"default\">Default</option>",
    "",
    palette_options,
    fixed = TRUE
  )
  palette_options <- paste0(
    palette_options,
    "<option value=\"single\">Single color</option>"
  )
  paste0(
    "<div id=\"networkPanel\" class=\"networkPanel\" role=\"dialog\" aria-modal=\"true\" aria-labelledby=\"networkTitle\"><div class=\"networkTop\">",
    "<div class=\"networkHeading\"><div class=\"networkHeadingText\"><div id=\"networkTitle\" class=\"networkTitle\">Condition / Topic GRN</div><div id=\"networkContext\" class=\"networkContext\"></div></div><div class=\"networkActions\"><button id=\"networkFit\" class=\"secondary\" type=\"button\" title=\"Fit all displayed nodes into the network panel.\">Fit view</button><button id=\"networkExport\" class=\"secondary\" type=\"button\" title=\"Export the current GRN as SVG.\">Export SVG</button><button id=\"networkBack\" type=\"button\" title=\"Return to the pathway plot.\">Pathways</button></div></div>",
    "<div class=\"networkTabs\" role=\"tablist\" aria-label=\"Network controls\"><button class=\"networkTab active\" type=\"button\" data-network-tab=\"filter\">Filter</button><button class=\"networkTab\" type=\"button\" data-network-tab=\"layout\">Layout</button><button class=\"networkTab\" type=\"button\" data-network-tab=\"appearance\">Appearance</button></div>",
    "<div id=\"networkTabFilter\" class=\"networkControls networkTabPanel active\" data-network-panel=\"filter\">",
    "<label>TF scope <select id=\"networkTfScope\"><option value=\"topic\" selected>In Topic</option><option value=\"all\">All TFs</option></select></label>",
    "<label>Topic theta &gt;= <select id=\"networkThetaPreset\"><option value=\"0.3\" selected>0.3</option><option value=\"0.5\">0.5</option><option value=\"0.7\">0.7</option><option value=\"custom\">Custom</option></select><input id=\"networkThetaCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"0.3\"/></label>",
    "<label>Topic Phi &gt;= <select id=\"networkPhiPreset\"><option value=\"0\" selected>0</option><option value=\"0.1\">0.1</option><option value=\"0.2\">0.2</option><option value=\"0.3\">0.3</option><option value=\"0.5\">0.5</option><option value=\"0.7\">0.7</option><option value=\"custom\">Custom</option></select><input id=\"networkPhiCutoff\" type=\"number\" min=\"0\" max=\"1\" step=\"0.01\" value=\"0\"/></label>",
    "<label>Primary topic only <input id=\"networkPrimaryOnly\" type=\"checkbox\"/></label>",
    "<label title=\"Keep links whose target RNA change and footprint change pass project cutoffs in the same direction; exclude significant opposing TF RNA changes.\">Correct direction only <input id=\"networkCorrectDirectionOnly\" type=\"checkbox\"/></label>",
    "<label title=\"Highlight keeps the full topic GRN and dims non-pathway elements. Filter restricts the GRN to the selected pathway.\">Pathway focus <select id=\"networkPathwayFocus\"><option value=\"highlight\" selected>Highlight</option><option value=\"filter\">Filter</option></select></label>",
    "<label>Network <select id=\"networkMode\"><option value=\"tf_gene\" selected>TF-gene</option><option value=\"tf_peak_gene\">TF-peak-gene</option></select></label>",
    "<label>Top TFs <input id=\"networkTopTf\" type=\"number\" min=\"1\" value=\"100\"/></label><label>Top links <input id=\"networkTopLinks\" type=\"number\" min=\"1\" value=\"300\"/></label>",
    "<label>Select node <select id=\"networkNodeSelect\" class=\"wideSelect\"></select></label></div>",
    "<div id=\"networkTabLayout\" class=\"networkControls networkTabPanel\" data-network-panel=\"layout\">",
    "<label>Layout <select id=\"networkLayout\"><option value=\"force\" selected>Force</option><option value=\"auto\">Auto</option><option value=\"radial\">Radial</option><option value=\"columns\">Columns</option><option value=\"bipartite\">Bipartite</option><option value=\"hierarchy\">Hierarchy</option><option value=\"concentric\">Concentric</option><option value=\"circle\">Circle</option><option value=\"grid\">Grid</option><option value=\"spiral\">Spiral</option><option value=\"clustered\">Clustered</option></select></label>",
    "<label>Spacing <input id=\"networkSpacing\" type=\"range\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/><input id=\"networkSpacingValue\" type=\"number\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/></label><button id=\"networkReset\" class=\"secondary\" type=\"button\" title=\"Restore the selected layout and fit all nodes.\">Reset layout</button></div>",
    "<div id=\"networkTabAppearance\" class=\"networkControls networkTabPanel\" data-network-panel=\"appearance\">",
    "<label>TF palette <select id=\"networkTfPalette\">", palette_options, "</select><input id=\"networkTfSingleColor\" class=\"singleColorInput\" type=\"color\" value=\"#B2182B\" hidden title=\"Fixed TF color\"/></label><label>Gene palette <select id=\"networkGenePalette\">", palette_options, "</select><input id=\"networkGeneSingleColor\" class=\"singleColorInput\" type=\"color\" value=\"#2166AC\" hidden title=\"Fixed gene color\"/></label><label>Edge palette <select id=\"networkEdgePalette\">", palette_options, "</select><input id=\"networkEdgeSingleColor\" class=\"singleColorInput\" type=\"color\" value=\"#59615D\" hidden title=\"Fixed edge color\"/></label>",
    "<label title=\"Smallest TF node half-height, in SVG units.\">TF box min <input id=\"networkTfMin\" type=\"number\" min=\"6\" max=\"40\" value=\"14\"/></label><label title=\"Largest TF node half-height, in SVG units. Maximum 40.\">TF box max <input id=\"networkTfMax\" type=\"number\" min=\"6\" max=\"40\" value=\"20\"/></label>",
    "<label>Labels <input id=\"networkLabels\" type=\"checkbox\" checked/></label><label>Arrows <input id=\"networkArrows\" type=\"checkbox\" checked/></label></div>",
    "<div id=\"networkStats\" class=\"networkStats\"></div></div>",
    "<div id=\"networkCanvas\" class=\"networkCanvas\"><div id=\"networkTooltip\" class=\"tooltip\"></div>",
    "<svg id=\"networkSvg\" viewBox=\"0 0 1600 900\"><defs><marker id=\"networkArrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerUnits=\"userSpaceOnUse\" markerWidth=\"11\" markerHeight=\"11\" orient=\"auto\"><path d=\"M0,0 L10,5 L0,10 z\" fill=\"#59615d\"/></marker></defs><g id=\"networkView\"><g id=\"networkEdges\"></g><g id=\"networkNodes\"></g><g id=\"networkLabelsLayer\"></g></g></svg></div></div>"
  )
}

.m3cr_condition_report_js <- function() {
  r"---(
const NS='http://www.w3.org/2000/svg',PAL=['#4E79A7','#E15759','#59A14F','#F28E2B','#B07AA1','#76B7B2','#EDC948','#9C755F','#FF9DA7','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#1F78B4','#33A02C','#E31A1C','#6A3D9A','#B15928','#17BECF','#BCBD22','#7F7F7F'];
const NETWORK_COLOR_PRESETS={default:['#2563eb','#dc2626','#9ca3af','#f59e0b'],npg:['#E64B35','#4DBBD5','#7E6148','#3C5488'],aaas:['#EE0000','#3B4992','#008B45','#631879'],nejm:['#BC3C29','#0072B5','#E18727','#20854E'],lancet:['#ED0000','#00468B','#42B540','#0099B4'],jama:['#374E55','#DF8F44','#B24745','#6A6599'],bmj:['#2B8CBE','#E34A33','#969696','#636363'],jco:['#CD534C','#0073C2','#868686','#EFC000'],ucscgb:['#FF0000','#0000FF','#999999','#666666'],d3:['#FF7F0E','#1F77B4','#7F7F7F','#2CA02C'],gephi:['#FF7F00','#377EB8','#999999','#4DAF4A'],observable:['#4269D0','#EF553B','#AAAAAA','#6CC5B0'],primer:['#CF222E','#0969DA','#8C959F','#57606A'],atlassian:['#FF5630','#0052CC','#97A0AF','#42526E'],iterm:['#CC241D','#458588','#A89984','#665C54'],locuszoom:['#D7191C','#2C7BB6','#ABD9E9','#7570B3'],igv:['#E41A1C','#377EB8','#999999','#4DAF4A'],cosmic:['#D73027','#4575B4','#969696','#542788'],uchicago:['#800000','#155F83','#8A9045','#767676'],startrek:['#CC0C00','#5C88DA','#9C9C9C','#FFCC00'],tron:['#FF410D','#0085CA','#6EE2FF','#F7C530'],futurama:['#FF6F00','#1B9E77','#A6A6A6','#7570B3'],rickandmorty:['#FAFD7C','#00B0C8','#B2DF8A','#808080'],simpsons:['#FED439','#709AE1','#8A9197','#D2AF81'],flatui:['#E74C3C','#3498DB','#95A5A6','#2C3E50'],frontiers:['#E64B35','#4DBBD5','#B09C85','#3C5488'],gsea:['#E31A1C','#1F78B4','#BDBDBD','#33A02C'],bootstrap5:['#DC3545','#0D6EFD','#6C757D','#198754'],material:['#F44336','#2196F3','#9E9E9E','#4CAF50'],tailwind:['#EF4444','#3B82F6','#9CA3AF','#10B981']};
const byId=id=>document.getElementById(id),cond1Select=byId('cond1Select'),cond2Select=byId('cond2Select'),cond1Color=byId('cond1Color'),cond2Color=byId('cond2Color'),topicSelect=byId('topicSelect'),tfSelect=byId('tfSelect'),pathwaySelect=byId('pathwaySelect'),conditionTopicMetric=byId('conditionTopicMetric'),pathwayScoreMethod=byId('pathwayScoreMethod'),pathwayDeOnly=byId('pathwayDeOnly');
const PAGE_SIZE={topic:20,tf:20,path:35};
let PAYLOAD={tf_activity:[],tf_gene:[],tf_peak_gene:[],gene_topics:[],condition_topic_expression:[],gene_expr:[],pathway_activity:[],pathway_gene_expression:[]},PAYLOAD_PROMISES=new Map(),LOADED_EDGE_CONDITIONS=new Set(),LOADED_PEAK_CONDITIONS=new Set(),CONDITION_COLORS=Object.assign({},REPORT_STATE.condition_colors||{}),EDGE_BY_COND=new Map(),PEAK_BY_COND_PAIR=new Map(),GENE_EXPR_INDEX=new Map(),GENE_ACTIVITY_INDEX=new Map(),EXPR_TOPIC_INDEX=new Map(),PATH_ACTIVITY_INDEX=new Map(),TF_TOPIC_INDEX=new Map(),TF_BY_CONDITION=new Map(),TF_LABELS=new Map(),PATH_ROWS=[],PAGE_INDEX={topic:0,tf:0,path:0},DATA_VERSION=0,RNA_TOPIC_CACHE={key:'',left:new Map(),right:new Map(),genes:0},MDS_RENDER_KEY='',ACTIVITY_RENDER_KEY='',ACTIVE_CONDITION_SIDE=1,ACTIVITY_VIEW={k:1,ox:0,oy:0},networkOpen=false,lastNetworkTrigger=null,networkState={nodes:[],edges:[],view:{x:0,y:0,k:1},drag:null,pan:null,selected:''};
let LOADED_PATHWAY_CONDITIONS=new Set(),CONDITION_PATHWAYS_BY_COND=new Map(),TF_EXPR_INDEX=new Map(),LOADING_COUNT=1;
function el(n,a){const x=document.createElementNS(NS,n);Object.entries(a||{}).forEach(([k,v])=>x.setAttribute(k,v));return x}function esc(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')}function num(x,d=0){const v=Number(x);return Number.isFinite(v)?v:d}function clamp(x,a,b){return Math.max(a,Math.min(b,x))}function topicNum(x){return Number(String(x||'').replace(/^Topic/,''))}function tfKey(x){return String(x||'').trim().toUpperCase()}function pathKey(x){return String(x||'').trim().toLowerCase().replace(/[^a-z0-9]+/g,'')}function splitGenes(x){return String(x||'').split(/[;,]/).map(s=>s.trim()).filter(Boolean)}
function setLoading(active,message='Loading...'){const overlay=byId('loadingOverlay');if(!overlay)return;if(active)LOADING_COUNT++;else LOADING_COUNT=Math.max(0,LOADING_COUNT-1);if(active)byId('loadingMessage').textContent=message;const busy=LOADING_COUNT>0;overlay.hidden=!busy;document.body.setAttribute('aria-busy',busy?'true':'false')}
tfKey=x=>{const key=String(x||'').trim().toUpperCase();return key==='TBET'?'TBX21':key};
function tooltip(id,ev,msg){const t=byId(id),host=t.offsetParent||t.parentElement;t.innerHTML=esc(msg).replace(/\n/g,'<br/>');t.style.display='block';const box=host.getBoundingClientRect(),left=ev.clientX-box.left+12,top=ev.clientY-box.top+12;t.style.left=clamp(left,8,Math.max(8,box.width-t.offsetWidth-8))+'px';t.style.top=clamp(top,8,Math.max(8,box.height-t.offsetHeight-8))+'px'}function hideTip(id){byId(id).style.display='none'}
function mdsLabelRows(){return GROUP_MDS.map(d=>{const raw=String(d.mds_label||d.display_label||d.comparison_label||d.group_label||'').replace(/_/g,' ').trim(),full=String(d.display_label||d.comparison_label||d.group_label||raw).replace(/_/g,' ').trim();return{raw,full,tokens:full.split(/\s+/).filter(Boolean)}})}
function mdsFullLabelMap(){return new Map(mdsLabelRows().map(d=>[d.raw,d.full]))}
function mdsShortLabelMap(){const rows=mdsLabelRows(),usable=rows.filter(d=>d.tokens.length);if(usable.length<2)return new Map(rows.map(d=>[d.raw,d.full]));let prefix=0,suffix=0,minLen=Math.min(...usable.map(d=>d.tokens.length));while(prefix<minLen-1&&usable.every(d=>d.tokens[prefix].toLowerCase()===usable[0].tokens[prefix].toLowerCase()))prefix++;while(suffix<minLen-prefix-1&&usable.every(d=>d.tokens[d.tokens.length-1-suffix].toLowerCase()===usable[0].tokens[usable[0].tokens.length-1-suffix].toLowerCase()))suffix++;return new Map(rows.map(d=>{const end=d.tokens.length-suffix,label=d.tokens.slice(prefix,end).join(' ');return[d.raw,label||d.full]}))}
function mdsTextWidth(text){const canvas=mdsTextWidth.canvas||(mdsTextWidth.canvas=document.createElement('canvas')),ctx=canvas.getContext('2d');ctx.font='700 18px Arial, Helvetica, sans-serif';return Math.ceil(ctx.measureText(String(text)).width*1.4)+18}
function mdsRectOverlap(a,b){return Math.max(0,Math.min(a.x+a.w,b.x+b.w)-Math.max(a.x,b.x))*Math.max(0,Math.min(a.y+a.h,b.y+b.h)-Math.max(a.y,b.y))}
function mdsLabelCandidates(q){const out=[],w=q.w,h=q.h,x=q.pointX,y=q.pointY;[14,28,44,64,88,116].forEach(gap=>[[x+gap,y-h/2,'right'],[x-w-gap,y-h/2,'left'],[x+gap,y-h-gap,'top-right'],[x-w-gap,y-h-gap,'top-left'],[x+gap,y+gap,'bottom-right'],[x-w-gap,y+gap,'bottom-left'],[x-w/2,y-h-gap,'top'],[x-w/2,y+gap,'bottom']].forEach((d,i)=>out.push({x:d[0],y:d[1],w,h,location:d[2],rank:i+gap/14})));[18,48,84].forEach(gap=>[-140,-112,-84,-56,-28,0,28,56,84,112,140].forEach(dy=>{out.push({x:x+gap,y:y+dy-h/2,w,h,location:'right-track',rank:12+gap/14+Math.abs(dy)/28});out.push({x:x-w-gap,y:y+dy-h/2,w,h,location:'left-track',rank:12+gap/14+Math.abs(dy)/28})}));for(let gy=82;gy<=650;gy+=28){out.push({x:80,y:gy,w,h,location:'left-margin',rank:35});out.push({x:680-w,y:gy,w,h,location:'right-margin',rank:35})}return out}
function mdsLabelPenalty(candidate,q,placed,points,bounds){const area=candidate.w*candidate.h,insideW=Math.max(0,Math.min(candidate.x+candidate.w,bounds.x+bounds.w)-Math.max(candidate.x,bounds.x)),insideH=Math.max(0,Math.min(candidate.y+candidate.h,bounds.y+bounds.h)-Math.max(candidate.y,bounds.y)),outside=area-insideW*insideH,overlap=placed.reduce((s,d)=>s+mdsRectOverlap(candidate,d),0),pointOverlap=points.reduce((s,p)=>s+mdsRectOverlap(candidate,{x:p.x-p.r-3,y:p.y-p.r-3,w:2*p.r+6,h:2*p.r+6})*(p.id===q.id?0:1),0),cx=clamp(q.pointX,candidate.x,candidate.x+candidate.w),cy=clamp(q.pointY,candidate.y,candidate.y+candidate.h),leader=Math.hypot(cx-q.pointX,cy-q.pointY);return outside*150+overlap*90+pointOverlap*14+leader*.16+candidate.rank*.12}
function mdsFinalLabelPenalty(candidate,q,placed,points,bounds){const area=candidate.w*candidate.h,insideW=Math.max(0,Math.min(candidate.x+candidate.w,bounds.x+bounds.w)-Math.max(candidate.x,bounds.x)),insideH=Math.max(0,Math.min(candidate.y+candidate.h,bounds.y+bounds.h)-Math.max(candidate.y,bounds.y)),outside=area-insideW*insideH,overlap=placed.reduce((s,d)=>s+mdsRectOverlap(candidate,d),0);return overlap*1e7+outside*1e5+mdsLabelPenalty(candidate,q,placed,points,bounds)}
function mdsLayoutLabels(items,W,H,pad){const bounds={x:pad+6,y:pad+5,w:W-2*pad-12,h:H-2*pad-10},points=items.map(q=>({x:q.x,y:q.y,r:q.selected?30:24,id:q.id}));items.forEach(q=>{q.w=mdsTextWidth(q.text);q.h=36;q.nearest=Math.min(...points.filter(p=>p.id!==q.id).map(p=>Math.hypot(q.x-p.x,q.y-p.y)),Infinity)});const ordered=items.slice().sort((a,b)=>Number(b.selected)-Number(a.selected)||a.nearest-b.nearest||a.id.localeCompare(b.id)),placed=[];ordered.forEach(q=>{const candidates=mdsLabelCandidates(q),best=candidates.reduce((winner,d)=>mdsLabelPenalty(d,q,placed,points,bounds)<mdsLabelPenalty(winner,q,placed,points,bounds)?d:winner,candidates[0]);Object.assign(q,best);placed.push(q)});for(let pass=0;pass<10;pass++)ordered.forEach(q=>{const others=ordered.filter(d=>d!==q),candidates=mdsLabelCandidates(q),best=candidates.reduce((winner,d)=>mdsLabelPenalty(d,q,others,points,bounds)<mdsLabelPenalty(winner,q,others,points,bounds)?d:winner,candidates[0]);Object.assign(q,best)});const final=[];ordered.forEach(q=>{const candidates=mdsLabelCandidates(q),best=candidates.reduce((winner,d)=>mdsFinalLabelPenalty(d,q,final,points,bounds)<mdsFinalLabelPenalty(winner,q,final,points,bounds)?d:winner,candidates[0]);Object.assign(q,best);final.push(q)});return items}
function mdsLeaderEndpoints(q){const bx=clamp(q.pointX,q.x,q.x+q.w),by=clamp(q.pointY,q.y,q.y+q.h),dx=bx-q.pointX,dy=by-q.pointY,len=Math.max(1,Math.hypot(dx,dy)),r=q.selected?30:24;return{x1:q.pointX+dx/len*r,y1:q.pointY+dy/len*r,x2:bx,y2:by}}
function styleMdsPlot(){byId('mdsLayer').querySelectorAll('circle').forEach(p=>{const selected=p.getAttribute('stroke')==='#111';p.setAttribute('r',selected?30:24);p.setAttribute('fill-opacity',selected?'.76':'.48');p.setAttribute('stroke-width',selected?5:2.5)})}
function columnarRows(x){const cols=x&&x.columns?x.columns:[],data=x&&x.data?x.data:[],n=data.length?data[0].length:0,out=new Array(n);for(let i=0;i<n;i++){const r={};for(let j=0;j<cols.length;j++)r[cols[j]]=data[j][i];out[i]=r}return out}
async function decodePayload(x){if(!x||!x.compressed_columnar)return x||{};if(!('DecompressionStream' in window))throw new Error('A current browser with deflate support is required.');const bin=atob(x.compressed_columnar),bytes=new Uint8Array(bin.length);for(let i=0;i<bin.length;i++)bytes[i]=bin.charCodeAt(i);const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate')),packed=JSON.parse(await new Response(stream).text()),out={};Object.keys(packed).forEach(k=>out[k]=k==='tf_topic'?packed[k]:columnarRows(packed[k]));return out}
function loadPayloadFile(file){if(!file)return Promise.resolve({});if(PAYLOAD_PROMISES.has(file))return PAYLOAD_PROMISES.get(file);const promise=new Promise((resolve,reject)=>{window.CRAFTGRN_CONDITION_PAYLOADS=window.CRAFTGRN_CONDITION_PAYLOADS||{};const done=()=>{const packed=window.CRAFTGRN_CONDITION_PAYLOADS[file];delete window.CRAFTGRN_CONDITION_PAYLOADS[file];decodePayload(packed).then(resolve,reject)},attempt=n=>{if(window.CRAFTGRN_CONDITION_PAYLOADS[file])return done();const s=document.createElement('script'),base=CONDITION_PAYLOAD.payload_base.replace(/\/$/,'')+'/'+file,url=new URL(base,document.baseURI),isFile=url.protocol==='file:';s.async=true;if(n&&!isFile)url.searchParams.set('retry',Date.now());s.src=url.href;s.onload=()=>{s.remove();done()};s.onerror=()=>{s.remove();if(n<2)setTimeout(()=>attempt(n+1),300*(n+1));else reject(new Error('Failed to load condition report payload after retries: '+file+' ('+url.href+')'))};document.head.appendChild(s)};attempt(0)}).catch(e=>{PAYLOAD_PROMISES.delete(file);throw e});PAYLOAD_PROMISES.set(file,promise);return promise}
function loadReportData(){return new Promise((resolve,reject)=>{window.CRAFTGRN_CONDITION_REPORT_DATA=window.CRAFTGRN_CONDITION_REPORT_DATA||{};const file=REPORT_DATA_SPEC.file,done=()=>decodePayload(window.CRAFTGRN_CONDITION_REPORT_DATA[file]).then(resolve,reject);if(window.CRAFTGRN_CONDITION_REPORT_DATA[file])return done();const s=document.createElement('script');s.async=true;s.src=REPORT_DATA_SPEC.base.replace(/\/$/,'')+'/'+file;s.onload=done;s.onerror=()=>reject(new Error('Failed to load compressed report data: '+file));document.head.appendChild(s)})}
function mergePayload(x){Object.entries(x||{}).forEach(([k,v])=>{if(!Array.isArray(v))return;if(['tf_gene','tf_peak_gene','gene_expr'].includes(k))PAYLOAD[k]=(PAYLOAD[k]||[]).concat(v);else PAYLOAD[k]=v});indexPayload();DATA_VERSION++;return PAYLOAD}
function loadOverviewPayload(){return loadPayloadFile(CONDITION_PAYLOAD.payload_file).then(mergePayload)}
function loadNetworkPayload(){const file=CONDITION_PAYLOAD.network_payload_file;if(!file)return Promise.resolve(PAYLOAD);return loadPayloadFile(file).then(mergePayload)}
function conditionPayloadFiles(cond){return(CONDITION_PAYLOAD.condition_payload_files||{})[String(cond)]||{}}
function overallMode(){return!cond1Select.value}function selectedConditions(){return[cond1Select.value,cond2Select.value].filter(Boolean)}
function indexConditionPart(cond,kind,x){const condition=String(cond);if(kind==='peak'){const cm=new Map();(x.tf_peak_gene||[]).forEach(d=>{const k=tfKey(d.tf)+'\t'+d.gene_key,a=cm.get(k)||[];a.push(d);cm.set(k,a)});PEAK_BY_COND_PAIR.set(condition,cm)}else if(kind==='pathway'){CONDITION_PATHWAYS_BY_COND.set(condition,(x.condition_pathways||[]).filter(d=>String(d.condition_id)===condition))}else{EDGE_BY_COND.set(condition,(x.tf_gene||[]).filter(d=>String(d.condition_id)===condition))}DATA_VERSION++}
function loadConditionPart(cond,kind){const files=conditionPayloadFiles(cond),file=kind==='peak'?files.peak_file:kind==='pathway'?files.pathway_file:files.edge_file,loaded=kind==='peak'?LOADED_PEAK_CONDITIONS:kind==='pathway'?LOADED_PATHWAY_CONDITIONS:LOADED_EDGE_CONDITIONS;if(!file||loaded.has(String(cond)))return Promise.resolve(PAYLOAD);setLoading(true,kind==='pathway'?'Loading condition pathways...':kind==='peak'?'Loading peak support...':'Loading direct targets...');return loadPayloadFile(file).then(x=>{indexConditionPart(cond,kind,x);loaded.add(String(cond));PAYLOAD_PROMISES.delete(file);return PAYLOAD}).finally(()=>setLoading(false))}
function loadSelectedConditionParts(kind){const conditions=selectedConditions();if(location.protocol!=='file:')return Promise.all(conditions.map(c=>loadConditionPart(c,kind)));return conditions.reduce((promise,condition)=>promise.then(()=>loadConditionPart(condition,kind)),Promise.resolve())}
function ensureSelectedConditionEdges(){return loadSelectedConditionParts('edge')}
function ensureSelectedConditionPeaks(){return loadSelectedConditionParts('peak')}
function ensureSelectedConditionPathways(){return cond1Select.value&&!cond2Select.value?loadConditionPart(cond1Select.value,'pathway'):Promise.resolve(PAYLOAD)}
function indexPayload(){EDGE_BY_COND=new Map();PEAK_BY_COND_PAIR=new Map();GENE_EXPR_INDEX=new Map();GENE_ACTIVITY_INDEX=new Map();EXPR_TOPIC_INDEX=new Map();PATH_ACTIVITY_INDEX=new Map();TF_EXPR_INDEX=new Map();PAYLOAD.tf_gene.forEach(d=>{const c=String(d.condition_id),a=EDGE_BY_COND.get(c)||[];a.push(d);EDGE_BY_COND.set(c,a)});PAYLOAD.tf_peak_gene.forEach(d=>{const c=String(d.condition_id),cm=PEAK_BY_COND_PAIR.get(c)||new Map(),k=tfKey(d.tf)+'\t'+d.gene_key,a=cm.get(k)||[];a.push(d);cm.set(k,a);PEAK_BY_COND_PAIR.set(c,cm)});(PAYLOAD.tf_activity||[]).forEach(d=>TF_EXPR_INDEX.set(String(d.condition_id)+'\t'+tfKey(d.tf_upper||d.tf),num(d.tf_expr,NaN)));(PAYLOAD.condition_topic_expression||[]).forEach(d=>EXPR_TOPIC_INDEX.set(String(d.condition_id)+'\t'+Number(d.topic_num),num(d.expression_topic_share,0)));(PAYLOAD.pathway_gene_expression||[]).forEach(d=>{const key=String(d.condition_id)+'\t'+String(d.gene_key);GENE_EXPR_INDEX.set(key,num(d.gene_expr,NaN));GENE_ACTIVITY_INDEX.set(key,d)});(PAYLOAD.pathway_activity||[]).forEach(d=>PATH_ACTIVITY_INDEX.set(Number(d.topic_num)+'\t'+pathKey(d.pathway_key)+'\t'+String(d.condition_id),d))}
function indexTfTopics(){TF_TOPIC_INDEX=new Map();TF_BY_CONDITION=new Map();TF_LABELS=new Map();const cols=TF_TOPIC&&TF_TOPIC.columns||[],data=TF_TOPIC&&TF_TOPIC.data||[],at=name=>data[cols.indexOf(name)]||[],condition=at('comparison_label').length?at('comparison_label'):at('group_label'),tf=at('tf_upper').length?at('tf_upper'):at('tf'),display=at('tf_display').length?at('tf_display'):at('tf'),topic=at('topic_num').length?at('topic_num'):at('topic'),theta=at('theta'),n=theta.length;for(let i=0;i<n;i++){const c=String(condition[i]||''),key=tfKey(tf[i]),t=topicNum(topic[i]),v=num(theta[i],NaN);if(!c||!key||!Number.isFinite(t)||!Number.isFinite(v))continue;const idx=c+'\t'+key+'\t'+t,old=TF_TOPIC_INDEX.get(idx);if(old===undefined||v>old)TF_TOPIC_INDEX.set(idx,v);if(!TF_BY_CONDITION.has(c))TF_BY_CONDITION.set(c,new Set());TF_BY_CONDITION.get(c).add(key);if(!TF_LABELS.has(key))TF_LABELS.set(key,String(display[i]||tf[i]||key))}}
function setActiveConditionSide(side){ACTIVE_CONDITION_SIDE=Number(side)===2?2:1;cond1Select.classList.toggle('conditionTarget',ACTIVE_CONDITION_SIDE===1);cond2Select.classList.toggle('conditionTarget',ACTIVE_CONDITION_SIDE===2)}
function state(){return{cond1:cond1Select.value,cond2:cond2Select.value,activeConditionSide:ACTIVE_CONDITION_SIDE,color1:conditionColorForId(cond1Select.value),color2:conditionColorForId(cond2Select.value),colors:CONDITION_COLORS,topic:topicNum(topicSelect.value),tf:tfSelect.value,pathway:pathwaySelect.value,metric:conditionTopicMetric.value}}function sendState(ready=false){window.parent.postMessage({type:ready?'m3cr-ready':'m3cr-state',state:state(),options:controlOptions()},'*')}function controlOptions(){return{conditions:Array.from(cond1Select.options).map(o=>({value:o.value,label:conditionLabel(o.value)})).filter(x=>x.value),topics:Array.from(topicSelect.options).map(o=>({value:o.value,label:o.textContent})),tfs:Array.from(tfSelect.options).map(o=>({value:o.value,label:o.textContent})),pathways:Array.from(pathwaySelect.options).map(o=>({value:o.value,label:o.textContent}))}}
function fillSelect(sel,rows,current){sel.replaceChildren();rows.forEach(r=>{const o=document.createElement('option');o.value=String(r.value);o.textContent=String(r.label);sel.appendChild(o)});if(rows.some(r=>String(r.value)===String(current)))sel.value=String(current)}
function pageRows(rows,key){const size=PAGE_SIZE[key],pages=Math.max(1,Math.ceil(rows.length/size));PAGE_INDEX[key]=clamp(PAGE_INDEX[key],0,pages-1);const start=PAGE_INDEX[key]*size;updatePager(key,rows.length,pages);return rows.slice(start,start+size)}
function updatePager(key,total,pages){const prefix=key==='topic'?'topic':key==='tf'?'tf':'path',status=byId(prefix+'PageStatus'),prev=byId(prefix+'Prev'),next=byId(prefix+'Next');status.textContent=(PAGE_INDEX[key]+1)+' / '+pages;status.title=total+' items';prev.disabled=PAGE_INDEX[key]===0;next.disabled=PAGE_INDEX[key]>=pages-1}
function movePage(key,delta,draw){PAGE_INDEX[key]=Math.max(0,PAGE_INDEX[key]+delta);draw()}
function showItemPage(rows,key,predicate){const i=rows.findIndex(predicate);if(i>=0)PAGE_INDEX[key]=Math.floor(i/PAGE_SIZE[key])}
function resetPages(){PAGE_INDEX={topic:0,tf:0,path:0}}
function conditionLabels(){const out=new Map();GROUP_MDS.forEach(d=>out.set(String(d.comparison_label||d.group_label),String(d.display_label||d.comparison_label||d.group_label)));return out}
function configuredCondition(defaults,key,conditions,exclude=''){const exact=String(defaults[key]||'');if(exact!==exclude&&conditions.includes(exact))return exact;const suffix=String(defaults[key+'_suffix']||'');if(suffix){const match=conditions.find(x=>x!==exclude&&String(x).endsWith(suffix));if(match)return match}return''}
function initControls(){const labels=conditionLabels(),rank=new Map((REPORT_STATE.condition_order||[]).map((x,i)=>[String(x),i])),conditions=(CONDITION_PAYLOAD.conditions||Array.from(labels.keys())).slice().sort((a,b)=>(rank.get(String(a))??1e9)-(rank.get(String(b))??1e9)||(labels.get(a)||a).localeCompare(labels.get(b)||b,undefined,{sensitivity:'base'})),topics=Array.from(new Set(GROUP_TOPIC.map(d=>topicNum(d.topic_num||d.topic)))).filter(Number.isFinite).sort((a,b)=>a-b),tfs=Array.from(TF_LABELS.entries()).sort((a,b)=>a[1].localeCompare(b[1],undefined,{sensitivity:'base'})),defaults=REPORT_STATE.defaults||{},configured2=configuredCondition(defaults,'condition_2',conditions),default2=configured2||(conditions.length===2?conditions[1]:''),configured1=configuredCondition(defaults,'condition_1',conditions,default2),default1=configured1||conditions.find(x=>x!==default2)||conditions[0]||'';fillSelect(cond1Select,[{value:'',label:'Overall pathways'}].concat(conditions.map(x=>({value:x,label:labels.get(x)||x}))),default1);fillSelect(cond2Select,[{value:'',label:'None'}].concat(conditions.map(x=>({value:x,label:labels.get(x)||x}))),default2);fillSelect(topicSelect,topics.map(x=>({value:x,label:'Topic '+x})),topics.includes(Number(defaults.topic))?Number(defaults.topic):topics[0]||1);fillSelect(tfSelect,[{value:'',label:'All'}].concat(tfs.map(x=>({value:x[0],label:x[1]}))),defaults.tf||'');fillSelect(pathwaySelect,[{value:'',label:'None'}],'');byId('shortConditionNames').checked=defaults.short_condition_names!==false}
function matchedTfKeys(){const c1=String(cond1Select.value||''),c2=String(cond2Select.value||''),a=TF_BY_CONDITION.get(c1)||new Set(),b=TF_BY_CONDITION.get(c2)||new Set();return c2?Array.from(a).filter(tf=>b.has(tf)):Array.from(a)}
function groupTheta(cond,topic){const mode=byId('thetaAggregation')?byId('thetaAggregation').value:'all';if(mode==='matched'&&cond2Select.value){const keys=matchedTfKeys();if(keys.length)return keys.reduce((s,tf)=>s+tfTheta(cond,tf,topic),0)/keys.length}const r=GROUP_TOPIC.find(d=>String(d.comparison_label||d.group_label)===String(cond)&&topicNum(d.topic_num||d.topic)===Number(topic));return r?num(r.theta_mean):0}
function pairwiseRnaTopicProfiles(){const c1=String(cond1Select.value||''),c2=String(cond2Select.value||''),cfg=directionConfig(),key=[DATA_VERSION,c1,c2,cfg.exprMin,cfg.pseudocount].join('\t');if(RNA_TOPIC_CACHE.key===key)return RNA_TOPIC_CACHE;const assignments=new Map();(PAYLOAD.gene_topics||[]).forEach(d=>{const gene=String(d.gene_key||''),topic=Number(d.topic_num),score=num(d.topic_score,0),old=assignments.get(gene);if(gene&&Number.isFinite(topic)&&(!old||score>old.score||(score===old.score&&topic<old.topic)))assignments.set(gene,{topic,score})});const left=new Map(),right=new Map();let genes=0;assignments.forEach((assignment,gene)=>{const a=auditGeneExpr(c1,gene,NaN),b=c2?auditGeneExpr(c2,gene,NaN):NaN,expressedA=Number.isFinite(a)&&a>=cfg.exprMin,expressedB=Number.isFinite(b)&&b>=cfg.exprMin;if(!expressedA&&!expressedB)return;genes++;if(c2){const delta=Math.log2((Math.max(0,Number.isFinite(a)?a:0)+cfg.pseudocount)/(Math.max(0,Number.isFinite(b)?b:0)+cfg.pseudocount));if(delta>0)left.set(assignment.topic,(left.get(assignment.topic)||0)+delta);else if(delta<0)right.set(assignment.topic,(right.get(assignment.topic)||0)-delta)}else if(expressedA){left.set(assignment.topic,(left.get(assignment.topic)||0)+Math.log2(Math.max(0,a)+cfg.pseudocount))}});const normalize=values=>{const total=Array.from(values.values()).reduce((sum,value)=>sum+value,0);if(total>0)values.forEach((value,topic)=>values.set(topic,value/total));return values};RNA_TOPIC_CACHE={key,left:normalize(left),right:normalize(right),genes};return RNA_TOPIC_CACHE}
function groupTopic(cond,topic){if(conditionTopicMetric&&conditionTopicMetric.value==='rna_delta'){const profiles=pairwiseRnaTopicProfiles();return String(cond)===String(cond2Select.value)?num(profiles.right.get(Number(topic)),0):num(profiles.left.get(Number(topic)),0)}return groupTheta(cond,topic)}
function conditionTopicMetricLabel(){if(conditionTopicMetric&&conditionTopicMetric.value==='rna_delta')return cond2Select.value?'share of positive pairwise RNA activity':'assigned-gene RNA activity share';const matched=byId('thetaAggregation')&&byId('thetaAggregation').value==='matched'&&cond2Select.value;return matched?'mean theta across matched TF documents':'mean theta across available TF documents'}
function tfTheta(cond,tf,topic){return num(TF_TOPIC_INDEX.get(String(cond)+'\t'+tfKey(tf)+'\t'+Number(topic)),0)}
function tfPrimaryTopic(cond,tf){let best=NaN,bestV=-Infinity;Array.from(topicSelect.options).forEach(o=>{const t=Number(o.value),v=tfTheta(cond,tf,t);if(v>bestV||(v===bestV&&t<best)){best=t;bestV=v}});return best}
function selectTopTfTopic(tf){if(!tf)return;const conds=[cond1Select.value,cond2Select.value].filter(Boolean),topics=Array.from(topicSelect.options).map(o=>Number(o.value));let best=topics[0],bestV=-Infinity;topics.forEach(t=>{const v=conds.reduce((s,c)=>s+tfTheta(c,tf,t),0)/Math.max(1,conds.length);if(v>bestV){bestV=v;best=t}});topicSelect.value=String(best);pathwaySelect.value=''}
function colorScale(v,limit,palette='rdbu',single=false){const diverging={rdbu:['#2166AC','#F7F7F7','#B2182B'],prgn:['#762A83','#F7F7F7','#1B7837'],brbg:['#8C510A','#F5F5F5','#01665E']};let p=diverging[palette];if(!p){const q=NETWORK_COLOR_PRESETS[palette]||NETWORK_COLOR_PRESETS.default;p=[q[1],'#F5F5F5',q[0]]}if(single){const t=clamp(v/Math.max(limit,1e-9),0,1);return mix('#E5E7EB',p[2],t)}const t=clamp(v/Math.max(limit,1e-9),-1,1);return t<0?mix(p[1],p[0],-t):mix(p[1],p[2],t)}function mix(a,b,t){const rgb=x=>{const q=parseInt(x.slice(1),16);return[(q>>16)&255,(q>>8)&255,q&255]},x=rgb(a),y=rgb(b);return'rgb('+x.map((v,i)=>Math.round(v+(y[i]-v)*t)).join(',')+')'}
function conditionColorForId(id){return CONDITION_COLORS[String(id)]||'#7F7F7F'}function conditionColor(side){return conditionColorForId(side===2?cond2Select.value:cond1Select.value)}function syncSelectedColorInputs(){cond1Color.value=conditionColor(1);cond2Color.value=conditionColor(2)}
function ensureConditionColors(){Array.from(cond1Select.options).filter(o=>o.value).sort((a,b)=>a.textContent.localeCompare(b.textContent,undefined,{sensitivity:'base'})).forEach((o,i)=>{if(!/^#[0-9a-f]{6}$/i.test(String(CONDITION_COLORS[o.value]||'')))CONDITION_COLORS[o.value]=PAL[i%PAL.length]});syncSelectedColorInputs()}
function fillRgb(fill){const s=String(fill||'').trim();if(s.startsWith('#')){const q=parseInt(s.slice(1),16);return[(q>>16)&255,(q>>8)&255,q&255]}const m=s.match(/[\d.]+/g);return m&&m.length>=3?m.slice(0,3).map(Number):[127,127,127]}function labelContrast(fill){const rgb=fillRgb(fill).map(v=>{v/=255;return v<=.03928?v/12.92:Math.pow((v+.055)/1.055,2.4)}),lum=.2126*rgb[0]+.7152*rgb[1]+.0722*rgb[2],dark=lum>.42;return{fill:dark?'#111':'#fff',stroke:dark?'#fff':'#333'}}
function edgeRows(cond){return EDGE_BY_COND.get(String(cond))||[]}function auditGeneExpr(cond,gene,fallback){const v=GENE_EXPR_INDEX.get(String(cond)+'\t'+String(gene));return Number.isFinite(v)?v:fallback}function tfExpression(cond,tf,fallback=NaN){const v=TF_EXPR_INDEX.get(String(cond)+'\t'+tfKey(tf));return Number.isFinite(v)?v:fallback}function decorateLinkDirection(r){const cfg=directionConfig(),paired=!!cond2Select.value,geneLfc=paired?Math.log2((Math.max(0,Number.isFinite(r.gexprA)?r.gexprA:0)+cfg.pseudocount)/(Math.max(0,Number.isFinite(r.gexprB)?r.gexprB:0)+cfg.pseudocount)):NaN,tfLfc=paired?Math.log2((Math.max(0,Number.isFinite(r.exprA??r.tfExprA)?(r.exprA??r.tfExprA):0)+cfg.pseudocount)/(Math.max(0,Number.isFinite(r.exprB??r.tfExprB)?(r.exprB??r.tfExprB):0)+cfg.pseudocount)):NaN,fpChange=paired?(cfg.fpMode==='log2fc'?Math.log2((Math.max(0,r.a)+cfg.pseudocount)/(Math.max(0,r.b)+cfg.pseudocount)):r.a-r.b):r.a,targetPass=paired&&Math.abs(geneLfc)>=cfg.geneCut,fpPass=paired&&Math.abs(fpChange)>=cfg.fpCut,sameDirection=targetPass&&fpPass&&Math.sign(geneLfc)===Math.sign(fpChange),tfOpposing=paired&&Number.isFinite(tfLfc)&&Math.abs(tfLfc)>=cfg.tfCut&&Math.sign(tfLfc)===-Math.sign(geneLfc);return Object.assign(r,{delta:paired?r.a-r.b:r.a,geneLfc,tfLfc,fpChange,correctDirection:paired&&sameDirection&&!tfOpposing,tfOpposing})}function edgePair(){const c1=cond1Select.value,c2=cond2Select.value,m=new Map();edgeRows(c1).forEach(d=>m.set(tfKey(d.tf)+'\t'+d.gene_key,{tf:d.tf,tfu:tfKey(d.tf),gene:String(d.gene_key),a:num(d.fp_sum),b:0,exprA:tfExpression(c1,d.tf,num(d.tf_expr,NaN)),exprB:NaN,gexprA:num(d.gene_expr,NaN),gexprB:NaN,nA:num(d.n_peaks),nB:0}));if(c2)edgeRows(c2).forEach(d=>{const k=tfKey(d.tf)+'\t'+d.gene_key,r=m.get(k)||{tf:d.tf,tfu:tfKey(d.tf),gene:String(d.gene_key),a:0,b:0,exprA:tfExpression(c1,d.tf,NaN),exprB:NaN,gexprA:NaN,gexprB:NaN,nA:0,nB:0};r.b=num(d.fp_sum);r.exprB=tfExpression(c2,d.tf,num(d.tf_expr,NaN));r.gexprB=num(d.gene_expr,NaN);r.nB=num(d.n_peaks);m.set(k,r)});return Array.from(m.values()).map(r=>decorateLinkDirection(Object.assign(r,{exprA:tfExpression(c1,r.tfu,r.exprA),exprB:c2?tfExpression(c2,r.tfu,r.exprB):NaN,gexprA:auditGeneExpr(c1,r.gene,r.gexprA),gexprB:c2?auditGeneExpr(c2,r.gene,r.gexprB):NaN})))}
function selectMdsCondition(id){const target=ACTIVE_CONDITION_SIDE===2?cond2Select:cond1Select,other=ACTIVE_CONDITION_SIDE===2?cond1Select:cond2Select,previous=target.value;if(other.value===id){target.value=id;other.value=previous}else target.value=id;pathwaySelect.value='';resetPages();refresh();refreshConditionData()}
function drawMds(){const g=byId('mdsLayer');g.replaceChildren();const W=760,H=760,pad=72,pointInset=38,rows=GROUP_MDS,sx=scale(rows.map(d=>num(d.MDS1)),pad+pointInset,W-pad-pointInset),sy=scale(rows.map(d=>num(d.MDS2)),H-pad-pointInset,pad+pointInset),selected=new Map([[cond1Select.value,conditionColor(1)],[cond2Select.value,conditionColor(2)]]),shorten=byId('shortConditionNames')&&byId('shortConditionNames').checked,labelMap=shorten?mdsShortLabelMap():mdsFullLabelMap(),labels=[];g.appendChild(el('line',{x1:pad,y1:H-pad,x2:W-pad,y2:H-pad,stroke:'#111','stroke-width':1.5}));g.appendChild(el('line',{x1:pad,y1:pad,x2:pad,y2:H-pad,stroke:'#111','stroke-width':1.5}));const xTitle=el('text',{x:W/2,y:H-18,'text-anchor':'middle','font-size':15,'font-weight':700,fill:'#303834'});xTitle.textContent='MDS1';g.appendChild(xTitle);const yTitle=el('text',{x:20,y:H/2,'text-anchor':'middle','font-size':15,'font-weight':700,fill:'#303834',transform:'rotate(-90 20 '+(H/2)+')'});yTitle.textContent='MDS2';g.appendChild(yTitle);rows.forEach(d=>{const id=String(d.comparison_label||d.group_label),x=sx(num(d.MDS1)),y=sy(num(d.MDS2)),sel=selected.get(id),c=conditionColorForId(id),raw=String(d.mds_label||d.display_label||id).replace(/_/g,' ').trim(),text=labelMap.get(raw)||raw,p=el('circle',{cx:x,cy:y,r:sel?19:14,fill:c,stroke:sel?'#111':'#fff','stroke-width':sel?4:2.5,style:'cursor:pointer','data-condition-id':id});p.addEventListener('click',()=>selectMdsCondition(id));p.addEventListener('mousemove',ev=>tooltip('mdsTooltip',ev,(d.display_label||id)+'\ndocuments: '+num(d.n_docs)));p.addEventListener('mouseout',()=>hideTip('mdsTooltip'));g.appendChild(p);labels.push({d,id,pointX:x,pointY:y,x,y,c,p,text,selected:!!sel})});mdsLayoutLabels(labels,W,H,pad).forEach(q=>{const e=mdsLeaderEndpoints(q);g.appendChild(el('line',{x1:e.x1,y1:e.y1,x2:e.x2,y2:e.y2,stroke:'#59635e','stroke-width':1.25,opacity:.66}));const t=el('text',{x:q.x+5,y:q.y+25,'text-anchor':'start','font-size':18,'font-weight':700,fill:q.c,style:'cursor:pointer;paint-order:stroke;stroke:#fff;stroke-width:5px;stroke-linejoin:round'});t.textContent=q.text;t.addEventListener('click',()=>q.p.dispatchEvent(new MouseEvent('click')));t.addEventListener('mousemove',ev=>tooltip('mdsTooltip',ev,(q.d.display_label||q.id)+'\ndocuments: '+num(q.d.n_docs)));t.addEventListener('mouseout',()=>hideTip('mdsTooltip'));g.appendChild(t)});styleMdsPlot();byId('mdsStats').textContent=rows.length+' conditions'}function scale(vals,lo,hi){const x=vals.filter(Number.isFinite);if(!x.length)return()=>((lo+hi)/2);const a=Math.min(...x),b=Math.max(...x);return v=>a===b?(lo+hi)/2:lo+(v-a)/(b-a)*(hi-lo)}
function activityRows(){const c1=cond1Select.value,c2=cond2Select.value,detailed=LOADED_EDGE_CONDITIONS.has(String(c1))&&(!c2||LOADED_EDGE_CONDITIONS.has(String(c2)));let rows;if(!detailed){rows=(PAYLOAD.tf_activity||[]).filter(d=>[c1,c2].includes(String(d.condition_id))).reduce((m,d)=>{const key=tfKey(d.tf_upper||d.tf),q=m.get(key)||{tf:d.tf,tfu:key,a:0,b:0,y:0,exprA:0,exprB:0};if(String(d.condition_id)===String(c1)){q.a=num(d.fp_sum);q.y=Math.max(q.y,num(d.n_targets));q.exprA=num(d.tf_expr)}else{q.b=num(d.fp_sum);q.y=Math.max(q.y,num(d.n_targets));q.exprB=num(d.tf_expr)}m.set(key,q);return m},new Map());rows=Array.from(rows.values())}else{const by=new Map();edgePair().forEach(r=>{const q=by.get(r.tfu)||{tf:r.tf,tfu:r.tfu,a:0,b:0,targets:new Set(),exprA:[],exprB:[]};q.a+=r.a;q.b+=r.b;if(!c2||Math.abs(r.delta)>=1)q.targets.add(r.gene);if(Number.isFinite(r.exprA))q.exprA.push(r.exprA);if(Number.isFinite(r.exprB))q.exprB.push(r.exprB);by.set(r.tfu,q)});rows=Array.from(by.values()).map(q=>Object.assign(q,{y:q.targets.size,exprA:q.exprA.length?Math.max(...q.exprA):0,exprB:q.exprB.length?Math.max(...q.exprB):0}))}const totalA=rows.reduce((s,d)=>s+Math.max(0,d.a),0),totalB=rows.reduce((s,d)=>s+Math.max(0,d.b),0),positive=[];rows.forEach(d=>{if(d.a>0&&totalA>0)positive.push(d.a/totalA);if(c2&&d.b>0&&totalB>0)positive.push(d.b/totalB)});const pseudo=Math.max(1e-12,(positive.length?Math.min(...positive):1e-8)/2);return rows.map(d=>{const shareA=totalA>0?Math.max(0,d.a)/totalA:0,shareB=totalB>0?Math.max(0,d.b)/totalB:0;return Object.assign(d,{shareA,shareB,deltaShare:shareA-shareB,relativeLog2:Math.log2((shareA+pseudo)/(shareB+pseudo)),expr:Math.max(d.exprA||0,d.exprB||0),logfc:c2?Math.log2((d.exprA+1)/(d.exprB+1)):tfTheta(c1,d.tfu,topicSelect.value)})})}
function differentialConditionLabels(g,c1,c2,leftX,rightX,y=27){if(c2){const l=el('text',{x:leftX,y,'text-anchor':'middle','font-size':16,'font-weight':700,fill:conditionColor(1)});l.textContent=conditionLabel(c1);g.appendChild(l);const r=el('text',{x:rightX,y,'text-anchor':'middle','font-size':16,'font-weight':700,fill:conditionColor(2)});r.textContent=conditionLabel(c2);g.appendChild(r)}else{const t=el('text',{x:(leftX+rightX)/2,y,'text-anchor':'middle','font-size':16,'font-weight':700,fill:conditionColor(1)});t.textContent=conditionLabel(c1);g.appendChild(t)}}
function activityViewDomain(lo,hi,offset){const full=Math.max(hi-lo,1e-9),span=full/ACTIVITY_VIEW.k,minCenter=lo+span/2,maxCenter=hi-span/2,center=clamp((lo+hi)/2+offset*full,minCenter,maxCenter);return[center-span/2,center+span/2]}
function resetActivityZoom(){ACTIVITY_VIEW={k:1,ox:0,oy:0};ACTIVITY_RENDER_KEY=''}
function changeActivityZoom(factor){ACTIVITY_VIEW.k=clamp(ACTIVITY_VIEW.k*factor,1,16);if(ACTIVITY_VIEW.k===1){ACTIVITY_VIEW.ox=0;ACTIVITY_VIEW.oy=0}ACTIVITY_RENDER_KEY='';drawActivity()}
function panelMessage(g,message){const t=el('text',{x:380,y:380,'text-anchor':'middle','font-size':18,'font-weight':700,fill:'#59635e'});t.textContent=message;g.appendChild(t)}
function drawActivity(){const g=byId('activityLayer');g.replaceChildren();const c1=cond1Select.value,c2=cond2Select.value;if(!c1){panelMessage(g,'Select a condition to view TF activity.');byId('activityStats').textContent='Overall pathway mode';return}const W=760,H=760,L=82,R=132,T=58,B=76,pointInset=22,rows=activityRows(),absDelta=rows.map(d=>Math.abs(d.deltaShare)).filter(v=>v>0&&Number.isFinite(v)).sort((a,b)=>a-b),soft=Math.max(1e-8,absDelta.length?absDelta[Math.floor((absDelta.length-1)*.5)]:1e-4),displayX=d=>c2?-Math.asinh(d.deltaShare/soft):d.shareA,xvals=rows.map(displayX),xlim=c2?Math.max(...xvals.map(Math.abs),1):Math.max(...xvals,1e-6),fullX=c2?[-xlim,xlim]:[0,xlim],yvals=rows.map(d=>Math.log2(d.y+1)).filter(Number.isFinite),yLoRaw=yvals.length?Math.min(...yvals):0,yHiRaw=yvals.length?Math.max(...yvals):1,ySpan=Math.max(yHiRaw-yLoRaw,.5),yPad=Math.min(.4,.07*ySpan),fullY=[Math.max(0,yLoRaw-yPad),Math.max(Math.max(0,yLoRaw-yPad)+.1,yHiRaw+yPad)],xd=activityViewDomain(fullX[0],fullX[1],ACTIVITY_VIEW.ox),yd=activityViewDomain(fullY[0],fullY[1],ACTIVITY_VIEW.oy),sx=scale(xd,L+pointInset,W-R-pointInset),sy=scale(yd,H-B-pointInset,T+pointInset),lim=Math.max(...rows.map(d=>Math.abs(d.logfc)).filter(Number.isFinite),1),exprMax=Math.max(...rows.map(d=>Math.log2(d.expr+1)),1);differentialConditionLabels(g,c1,c2,L+(W-L-R)/4,L+3*(W-L-R)/4);g.appendChild(el('line',{x1:L,y1:H-B,x2:W-R,y2:H-B,stroke:'#111','stroke-width':1.5}));g.appendChild(el('line',{x1:L,y1:T,x2:L,y2:H-B,stroke:'#111','stroke-width':1.5}));[0,.25,.5,.75,1].forEach(frac=>{const value=yd[0]+frac*(yd[1]-yd[0]),yy=sy(value),count=Math.max(0,Math.round(Math.pow(2,value)-1));g.appendChild(el('line',{x1:L,y1:yy,x2:W-R,y2:yy,stroke:'#8a9690','stroke-width':.8,opacity:frac===0?.45:.2}));const tick=el('text',{x:L-10,y:yy+4,'text-anchor':'end','font-size':12,'font-weight':700,fill:'#46524d'});tick.textContent=String(count);g.appendChild(tick)});[0,.25,.5,.75,1].forEach(frac=>{const tv=xd[0]+frac*(xd[1]-xd[0]),xx=sx(tv),raw=c2?-100*soft*Math.sinh(tv):100*tv;g.appendChild(el('line',{x1:xx,y1:T,x2:xx,y2:H-B,stroke:Math.abs(tv)<1e-10?'#555':'#8a9690','stroke-width':Math.abs(tv)<1e-10?1.2:.8,'stroke-dasharray':Math.abs(tv)<1e-10?'4 4':'2 5',opacity:Math.abs(tv)<1e-10?1:.3}));const tick=el('text',{x:xx,y:H-B+22,'text-anchor':'middle','font-size':12,'font-weight':700,fill:'#46524d'});tick.textContent=(raw>0?'+':'')+raw.toFixed(Math.abs(raw)<.1?2:1);g.appendChild(tick)});rows.slice().sort((a,b)=>Number(tfKey(tfSelect.value)===a.tfu)-Number(tfKey(tfSelect.value)===b.tfu)).forEach(d=>{const dx=displayX(d),dy=Math.log2(d.y+1);if(dx<xd[0]||dx>xd[1]||dy<yd[0]||dy>yd[1])return;const px=sx(dx),py=sy(dy),r=2.2+3.3*Math.sqrt(Math.log2(d.expr+1)/exprMax),selected=tfKey(tfSelect.value)===d.tfu,base=c2?(d.deltaShare<0?conditionColor(2):conditionColor(1)):conditionColor(1),strength=.78+.18*clamp(Math.abs(d.logfc)/lim,0,1),p=el('circle',{cx:px,cy:py,r:selected?r+5:r,fill:selected?'#111':base,opacity:selected?1:strength,stroke:'#fff','stroke-width':selected?3:.8,style:'cursor:pointer'}),msg=d.tf+'\n'+conditionLabel(c1)+' normalized FP: '+(100*d.shareA).toFixed(3)+'%'+(c2?'\n'+conditionLabel(c2)+' normalized FP: '+(100*d.shareB).toFixed(3)+'%\n'+conditionLabel(c1)+' - '+conditionLabel(c2)+': '+(100*d.deltaShare).toFixed(3)+' percentage points\nlog2 relative activity: '+d.relativeLog2.toFixed(3):'')+'\ntarget genes: '+d.y+'\nTF expression: '+d.expr.toFixed(3)+(c2?'\nRNA log2FC: '+d.logfc.toFixed(3):'\nselected-topic theta: '+d.logfc.toFixed(3));p.addEventListener('click',()=>selectTf(d.tfu));p.addEventListener('mousemove',ev=>tooltip('activityTooltip',ev,msg));p.addEventListener('mouseout',()=>hideTip('activityTooltip'));g.appendChild(p);if(selected){const right=px<W-R-95,t=el('text',{x:px+(right?r+12:-r-12),y:py-9,'text-anchor':right?'start':'end','font-size':18,'font-weight':700,fill:'#111',style:'cursor:pointer;paint-order:stroke;stroke:#fff;stroke-width:5px;stroke-linejoin:round'});t.textContent=d.tf;t.addEventListener('click',()=>selectTf(d.tfu));g.appendChild(t)}});const xt=el('text',{x:(L+W-R)/2,y:H-14,'text-anchor':'middle','font-size':13,'font-weight':700});xt.textContent=c2?'Normalized FP difference (pp): '+conditionLabel(c1)+' - '+conditionLabel(c2):'Normalized FP: '+conditionLabel(c1);g.appendChild(xt);const yt=el('text',{x:24,y:(T+H-B)/2,'text-anchor':'middle','font-size':18,'font-weight':700,transform:'rotate(-90 24 '+((T+H-B)/2)+')'});yt.textContent=c2?'unique targets with abs(delta FP) >= 1':'unique target genes';g.appendChild(yt);byId('activityStats').textContent=rows.length+' TFs | '+ACTIVITY_VIEW.k.toFixed(1)+'x'}
function renderColorKey(){return Object.keys(CONDITION_COLORS).sort().map(k=>k+':'+CONDITION_COLORS[k]).join('|')}
const drawMdsUncached=drawMds;drawMds=()=>{const key=[cond1Select.value,cond2Select.value,byId('shortConditionNames')&&byId('shortConditionNames').checked,renderColorKey()].join('\t');if(key===MDS_RENDER_KEY)return;MDS_RENDER_KEY=key;drawMdsUncached()};
const drawActivityUncached=drawActivity;drawActivity=()=>{const key=[cond1Select.value,cond2Select.value,tfSelect.value,DATA_VERSION,conditionColor(1),conditionColor(2),ACTIVITY_VIEW.k,ACTIVITY_VIEW.ox,ACTIVITY_VIEW.oy].join('\t');if(key===ACTIVITY_RENDER_KEY)return;ACTIVITY_RENDER_KEY=key;drawActivityUncached()};
function singleAxisStart(rows,labelFn,minX,maxX){const longest=Math.max(0,...rows.map(d=>String(labelFn(d)||'').length));return clamp(34+longest*10.2,minX,maxX)}
function butterflyAxes(g,c1,c2,W,H,L,R,T,B,gap,singleStart=L){const cx=c2?W/2:singleStart,xw=c2?(W-L-R-2*gap)/2:W-singleStart-R,axisY=T-10;if(!c2)g.appendChild(el('line',{x1:cx,y1:T-12,x2:cx,y2:H-B,stroke:'#111','stroke-width':1.5}));if(c2){g.appendChild(el('line',{x1:cx-gap-xw,y1:axisY,x2:cx-gap,y2:axisY,stroke:conditionColor(1),'stroke-width':1.5}));g.appendChild(el('line',{x1:cx+gap,y1:axisY,x2:cx+gap+xw,y2:axisY,stroke:conditionColor(2),'stroke-width':1.5}))}else g.appendChild(el('line',{x1:cx,y1:axisY,x2:cx+xw,y2:axisY,stroke:conditionColor(1),'stroke-width':1.5}));[.25,.5,.75,1].forEach(v=>{if(c2){const xl=cx-gap-v*xw,xr=cx+gap+v*xw;g.appendChild(el('line',{x1:xl,y1:T-15,x2:xl,y2:H-B,stroke:'#8a9690','stroke-width':1,opacity:.25}));g.appendChild(el('line',{x1:xr,y1:T-15,x2:xr,y2:H-B,stroke:'#8a9690','stroke-width':1,opacity:.25}))}});const t1=el('text',{x:c2?cx-gap-xw/2:cx+xw/2,y:18,'text-anchor':'middle','font-size':16,'font-weight':700,fill:conditionColor(1)});t1.textContent=conditionLabel(c1);g.appendChild(t1);if(c2){const t2=el('text',{x:cx+gap+xw/2,y:18,'text-anchor':'middle','font-size':16,'font-weight':700,fill:conditionColor(2)});t2.textContent=conditionLabel(c2);g.appendChild(t2)}[0,.5,1].forEach(v=>{const addTick=(x,color)=>{g.appendChild(el('line',{x1:x,y1:axisY-4,x2:x,y2:axisY+4,stroke:color,'stroke-width':1.2}));const q=el('text',{x,y:axisY-7,'text-anchor':'middle','font-size':12,'font-weight':700,fill:color});q.textContent=v.toFixed(1);g.appendChild(q)};if(c2){addTick(cx-gap-v*xw,conditionColor(1));addTick(cx+gap+v*xw,conditionColor(2))}else addTick(cx+v*xw,conditionColor(1))});return{cx,xw}}
function conditionTopicRows(){const c1=cond1Select.value,c2=cond2Select.value,topics=Array.from(topicSelect.options).map(o=>Number(o.value));return topics.map(t=>({topic:t,a:groupTopic(c1,t),b:c2?groupTopic(c2,t):0})).sort((a,b)=>Math.max(b.a,b.b)-Math.max(a.a,a.b)||a.topic-b.topic)}
function drawButterfly(){const g=byId('butterflyLayer');g.replaceChildren();const c1=cond1Select.value,c2=cond2Select.value;if(!c1){panelMessage(g,'Select a condition to view topic activity.');byId('butterflyStats').textContent='Overall pathway mode';return}const W=760,H=760,L=24,R=24,T=58,B=30,gap=70,metric=conditionTopicMetricLabel(),all=conditionTopicRows(),vals=pageRows(all,'topic'),ax=butterflyAxes(g,c1,c2,W,H,L,R,T,B,gap,singleAxisStart(vals,d=>'Topic '+d.topic,165,245)),rowH=Math.min(40,(H-T-B)/Math.max(vals.length,1)),barH=Math.max(18,Math.min(28,rowH*.7));vals.forEach((d,i)=>{const y=T+i*rowH,w1=clamp(d.a,0,1)*ax.xw,w2=clamp(d.b,0,1)*ax.xw,sel=Number(topicSelect.value)===d.topic;if(i%2===1)g.appendChild(el('rect',{x:0,y:y-3,width:W,height:rowH,fill:'#f7f8f6'}));const label=el('text',{x:c2?ax.cx:ax.cx-10,y:y+barH*.79,'text-anchor':c2?'middle':'end','font-size':18,'font-weight':700,fill:sel?'#991B1B':'#171717',style:'cursor:pointer'});label.textContent='Topic '+d.topic;const choose=()=>{topicSelect.value=String(d.topic);pathwaySelect.value='';PAGE_INDEX.tf=0;PAGE_INDEX.path=0;refresh()};label.addEventListener('click',choose);g.appendChild(label);const a=el('rect',{x:c2?ax.cx-gap-w1:ax.cx,y,width:Math.max(2,w1),height:barH,fill:conditionColor(1),stroke:sel?'#111':'none','stroke-width':sel?3:0,style:'cursor:pointer'});a.addEventListener('click',choose);a.addEventListener('mousemove',ev=>tooltip('butterflyTooltip',ev,conditionLabel(c1)+' / Topic '+d.topic+'\n'+metric+': '+d.a.toFixed(4)));a.addEventListener('mouseout',()=>hideTip('butterflyTooltip'));g.appendChild(a);const av=el('text',{x:c2?Math.max(78,ax.cx-gap-w1-5):ax.cx+w1+5,y:y+barH*.73,'text-anchor':c2?'end':'start','font-size':12,'font-weight':700,fill:conditionColor(1),style:'pointer-events:none;paint-order:stroke;stroke:#fff;stroke-width:3px'});av.textContent=d.a.toFixed(3);g.appendChild(av);if(c2){const b=el('rect',{x:ax.cx+gap,y,width:Math.max(2,w2),height:barH,fill:conditionColor(2),stroke:sel?'#111':'none','stroke-width':sel?3:0,style:'cursor:pointer'});b.addEventListener('click',choose);b.addEventListener('mousemove',ev=>tooltip('butterflyTooltip',ev,conditionLabel(c2)+' / Topic '+d.topic+'\n'+metric+': '+d.b.toFixed(4)));b.addEventListener('mouseout',()=>hideTip('butterflyTooltip'));g.appendChild(b);const bv=el('text',{x:Math.min(W-78,ax.cx+gap+w2+5),y:y+barH*.73,'text-anchor':'start','font-size':12,'font-weight':700,fill:conditionColor(2),style:'pointer-events:none;paint-order:stroke;stroke:#fff;stroke-width:3px'});bv.textContent=d.b.toFixed(3);g.appendChild(bv)}});const start=all.length?PAGE_INDEX.topic*PAGE_SIZE.topic+1:0;byId('butterflyStats').textContent=start+'-'+(start+vals.length-1)+' of '+all.length}
function tfButterflyRows(){const c1=cond1Select.value,c2=cond2Select.value,topic=Number(topicSelect.value),selected=tfKey(tfSelect.value);if(selected)return conditionTopicRows().map(d=>({key:'topic:'+d.topic,label:'Topic '+d.topic,topic:d.topic,a:tfTheta(c1,selected,d.topic),b:c2?tfTheta(c2,selected,d.topic):0}));const rows=Array.from(TF_LABELS,([key,label])=>({key,label,a:tfTheta(c1,key,topic),b:c2?tfTheta(c2,key,topic):0,delta:tfTheta(c1,key,topic)-(c2?tfTheta(c2,key,topic):0)}));if(!c2)return rows.sort((a,b)=>b.a-a.a||a.label.localeCompare(b.label,undefined,{sensitivity:'base'})).slice(0,20);const pos=rows.filter(d=>d.delta>=0).sort((a,b)=>b.delta-a.delta||a.label.localeCompare(b.label,undefined,{sensitivity:'base'})).slice(0,10),neg=rows.filter(d=>d.delta<0).sort((a,b)=>a.delta-b.delta||a.label.localeCompare(b.label,undefined,{sensitivity:'base'})).slice(0,10);return pos.concat(neg)}
function tfButterflyRowSelected(d){return tfSelect.value?d.topic===Number(topicSelect.value):d.key===tfKey(tfSelect.value)}
function selectTfInTopic(tf){tfSelect.value=tfKey(tf);pathwaySelect.value='';refresh();refreshConditionData()}
function drawTfButterfly(){const g=byId('tfButterflyLayer');g.replaceChildren();const c1=cond1Select.value,c2=cond2Select.value;if(!c1){panelMessage(g,'Select a condition to view TF probability.');byId('tfButterflyStats').textContent='Overall pathway mode';return}const W=760,H=760,L=24,R=24,T=58,B=30,gap=74,selectedTf=tfKey(tfSelect.value),topicMode=!!selectedTf,all=tfButterflyRows(),rows=pageRows(all,'tf'),ax=butterflyAxes(g,c1,c2,W,H,L,R,T,B,gap,singleAxisStart(rows,d=>d.label,180,330)),rowH=Math.min(38,(H-T-B)/Math.max(rows.length,1)),barH=Math.max(18,Math.min(27,rowH*.7));rows.forEach((d,i)=>{const y=T+i*rowH,w1=clamp(d.a,0,1)*ax.xw,w2=clamp(d.b,0,1)*ax.xw,sel=tfButterflyRowSelected(d);if(i%2===1)g.appendChild(el('rect',{x:0,y:y-3,width:W,height:rowH,fill:'#f7f8f6'}));const label=el('text',{x:c2?ax.cx:ax.cx-10,y:y+barH*.79,'text-anchor':c2?'middle':'end','font-size':17,'font-weight':700,fill:sel?'#991B1B':'#171717',style:'cursor:pointer'});label.textContent=d.label;const choose=topicMode?()=>{topicSelect.value=String(d.topic);pathwaySelect.value='';PAGE_INDEX.topic=0;PAGE_INDEX.path=0;refresh()}:()=>selectTfInTopic(d.key);label.addEventListener('click',choose);g.appendChild(label);const detail=topicMode?tfLabel()+' / '+d.label:d.label+' / Topic '+topicSelect.value,a=el('rect',{x:c2?ax.cx-gap-w1:ax.cx,y,width:Math.max(2,w1),height:barH,fill:conditionColor(1),stroke:sel?'#111':'none','stroke-width':sel?3:0,style:'cursor:pointer'});a.addEventListener('click',choose);a.addEventListener('mousemove',ev=>tooltip('tfButterflyTooltip',ev,conditionLabel(c1)+' / '+detail+'\nP(topic | condition::TF): '+d.a.toFixed(4)));a.addEventListener('mouseout',()=>hideTip('tfButterflyTooltip'));g.appendChild(a);if(c2){const b=el('rect',{x:ax.cx+gap,y,width:Math.max(2,w2),height:barH,fill:conditionColor(2),stroke:sel?'#111':'none','stroke-width':sel?3:0,style:'cursor:pointer'});b.addEventListener('click',choose);b.addEventListener('mousemove',ev=>tooltip('tfButterflyTooltip',ev,conditionLabel(c2)+' / '+detail+'\nP(topic | condition::TF): '+d.b.toFixed(4)));b.addEventListener('mouseout',()=>hideTip('tfButterflyTooltip'));g.appendChild(b)}});const start=all.length?PAGE_INDEX.tf*PAGE_SIZE.tf+1:0;byId('tfButterflyStats').textContent=topicMode?tfLabel()+' | '+start+'-'+(start+rows.length-1)+' of '+all.length+' topics':start+'-'+(start+rows.length-1)+' of '+all.length+' TFs'}
function conditionLabel(id){if(!id)return'Overall topic pathways';const d=GROUP_MDS.find(x=>String(x.comparison_label||x.group_label)===String(id)),full=String(d&&(d.display_label||d.comparison_label)||id).replace(/_/g,' ');if(!d||!byId('shortConditionNames')||!byId('shortConditionNames').checked)return full;const raw=String(d.mds_label||d.display_label||d.comparison_label||id).replace(/_/g,' ').trim();return mdsShortLabelMap().get(raw)||full}
function topicGeneScore(gene,topic){let v=NaN;PAYLOAD.gene_topics.forEach(d=>{if(String(d.gene_key)===String(gene)&&Number(d.topic_num)===Number(topic))v=Math.max(Number.isFinite(v)?v:-Infinity,num(d.topic_score,NaN))});return v}function topicGenes(){const t=Number(topicSelect.value),cut=num(byId('networkPhiCutoff')&&byId('networkPhiCutoff').value,0);return new Set(PAYLOAD.gene_topics.filter(d=>Number(d.topic_num)===t&&num(d.topic_score)>=cut).map(d=>String(d.gene_key)))}
function selectedTfTargets(cond,tf){return new Set(edgeRows(cond).filter(d=>tfKey(d.tf)===tfKey(tf)).map(d=>String(d.gene_key)))}
function pathwayActivity(topic,key,condition){return PATH_ACTIVITY_INDEX.get(Number(topic)+'\t'+pathKey(key)+'\t'+String(condition))||null}
function pathwayScore(row){if(!row)return NaN;return num(row[pathwayScoreMethod.value],NaN)}
function directionConfig(){const x=REPORT_STATE.link_direction||{};return{geneCut:num(x.gene_log2fc_cutoff,1),tfCut:num(x.tf_opposition_log2fc_cutoff,num(x.gene_log2fc_cutoff,1)),fpCut:num(x.fp_cutoff,.5),fpMode:String(x.fp_filter_mode||'delta'),exprMin:num(x.expression_min,0),pseudocount:Math.max(1e-12,num(x.pseudocount,1))}}
function pairGeneLog2fc(gene){const c2=cond2Select.value;if(!c2)return NaN;const cfg=directionConfig(),a=auditGeneExpr(cond1Select.value,gene,NaN),b=auditGeneExpr(c2,gene,NaN);if(!Number.isFinite(a)&&!Number.isFinite(b))return NaN;return Math.log2((Math.max(0,Number.isFinite(a)?a:0)+cfg.pseudocount)/(Math.max(0,Number.isFinite(b)?b:0)+cfg.pseudocount))}
function pairGeneExpressed(gene){const cfg=directionConfig(),a=auditGeneExpr(cond1Select.value,gene,NaN),b=cond2Select.value?auditGeneExpr(cond2Select.value,gene,NaN):NaN;return(Number.isFinite(a)&&a>=cfg.exprMin)||(Number.isFinite(b)&&b>=cfg.exprMin)}
function pathwayGeneActivity(condition,gene){return GENE_ACTIVITY_INDEX.get(String(condition)+'\t'+String(gene))||null}
function dynamicPathwayScore(genes,condition){const rows=Array.from(new Set(genes)).map(g=>pathwayGeneActivity(condition,g)).filter(Boolean),m=rows.length;if(m<3)return{score:NaN,n:m};if(pathwayScoreMethod.value==='mean_gene_zscore')return{score:rows.reduce((s,d)=>s+num(d.gene_zscore,0),0)/m,n:m};const n=num(rows[0].n_expression_universe,0),expected=m*(n+1)/2,variance=m*(n-m)*(n+1)/12,score=variance>0?(rows.reduce((s,d)=>s+num(d.expression_rank,0),0)-expected)/Math.sqrt(variance):0;return{score,n:m}}
function pathwayStars(q){return q<=1e-4?'****':q<=1e-3?'***':q<=1e-2?'**':q<=.05?'*':''}
function orderedPathGenes(genes,condition){return genes.map(g=>({gene:g,expr:auditGeneExpr(condition,g,NaN)})).filter(d=>Number.isFinite(d.expr)).sort((a,b)=>b.expr-a.expr||a.gene.localeCompare(b.gene,undefined,{sensitivity:'base'}))}
function pathwayRows(){const topic=Number(topicSelect.value),c1=cond1Select.value,c2=cond2Select.value;if(!c1)return PATHWAYS.filter(d=>topicNum(d.topic_num||d.topic)===topic&&num(d.padj,1)<=.05).map(d=>{const genes=splitGenes(d.genes||d.overlap_genes),overlap=Math.max(num(d.gene_in,0),genes.length),score=Number.isFinite(Number(d.combined_score))?Math.max(0,Number(d.combined_score)):-Math.log10(Math.max(num(d.padj,1),1e-300));return{key:pathKey(d.pathway_key||d.pathway_norm_key||d.pathway),pathway:String(d.pathway),scoreA:score,scoreB:NaN,delta:score,padj:num(d.padj,1),overlap,deOverlap:0,deGenes:[],genesA:genes,genesB:genes,genes,overall:true}}).sort((a,b)=>a.padj-b.padj||b.overlap-a.overlap||a.pathway.localeCompare(b.pathway));const condition1First=!c2||groupTopic(c1,topic)>=groupTopic(c2,topic),deOnly=!!(c2&&pathwayDeOnly.checked),geneCut=directionConfig().geneCut;let rows=PATHWAYS.filter(d=>topicNum(d.topic_num||d.topic)===topic&&num(d.padj,1)<=.05).map(d=>{const key=pathKey(d.pathway_key||d.pathway_norm_key||d.pathway),genes=splitGenes(d.genes||d.overlap_genes).filter(pairGeneExpressed),deGenes=c2?genes.filter(g=>Math.abs(pairGeneLog2fc(g))>=geneCut):[],a=dynamicPathwayScore(genes,c1),b=c2?dynamicPathwayScore(genes,c2):{score:NaN,n:0};return{key,pathway:String(d.pathway),scoreA:a.score,scoreB:b.score,delta:c2?b.score-a.score:a.score,padj:num(d.padj,1),overlap:genes.length,deOverlap:deGenes.length,deGenes,genesA:genes,genesB:genes,genes,activityA:a,activityB:b}}).filter(r=>Number.isFinite(r.scoreA)&&(!c2||Number.isFinite(r.scoreB))&&r.overlap>=3);const displayRows=x=>deOnly?x.filter(r=>r.deOverlap>0):x,tf=tfSelect.value;if(tf){const ta=selectedTfTargets(c1,tf),tb=c2?selectedTfTargets(c2,tf):new Set();rows=rows.map(r=>{const tfA=r.genes.filter(g=>ta.has(g)),tfB=c2?r.genes.filter(g=>tb.has(g)):[],tfUnion=Array.from(new Set(tfA.concat(tfB)));return Object.assign(r,{tfA,tfB,tfUnion,tfCountA:tfA.length,tfCountB:tfB.length,tfCountUnion:tfUnion.length})}).filter(r=>r.tfCountUnion>0);if(!c2)return displayRows(rows.sort((a,b)=>b.tfCountA-a.tfCountA||a.padj-b.padj||a.pathway.localeCompare(b.pathway)));const c1Tf=rows.filter(d=>d.delta<=0).sort((a,b)=>b.tfCountUnion-a.tfCountUnion||a.delta-b.delta||a.padj-b.padj||a.pathway.localeCompare(b.pathway)),c2Tf=rows.filter(d=>d.delta>0).sort((a,b)=>b.tfCountUnion-a.tfCountUnion||b.delta-a.delta||a.padj-b.padj||a.pathway.localeCompare(b.pathway));return displayRows(condition1First?c1Tf.concat(c2Tf):c2Tf.concat(c1Tf))}if(!c2)return displayRows(rows.sort((a,b)=>b.scoreA-a.scoreA||a.padj-b.padj||a.pathway.localeCompare(b.pathway)));const c1Rows=rows.filter(d=>d.delta<=0).sort((a,b)=>a.delta-b.delta||a.padj-b.padj||a.pathway.localeCompare(b.pathway)),c2Rows=rows.filter(d=>d.delta>0).sort((a,b)=>b.delta-a.delta||a.padj-b.padj||a.pathway.localeCompare(b.pathway));return displayRows(condition1First?c1Rows.concat(c2Rows):c2Rows.concat(c1Rows))}
const expressionPathwayRows=pathwayRows;
pathwayRows=function(){const c1=cond1Select.value,c2=cond2Select.value;if(!c1||c2)return expressionPathwayRows();const topic=Number(topicSelect.value),source=CONDITION_PATHWAYS_BY_COND.get(String(c1))||[],tf=tfSelect.value,targets=tf?selectedTfTargets(c1,tf):new Set();let rows=source.filter(d=>topicNum(d.topic_num||d.topic)===topic&&num(d.padj,1)<=.05).map(d=>{const genes=splitGenes(d.genes||d.overlap_genes),overlap=Math.max(num(d.gene_in,0),genes.length),score=Number.isFinite(Number(d.combined_score))?Math.max(0,Number(d.combined_score)):-Math.log10(Math.max(num(d.padj,1),1e-300)),tfA=tf?genes.filter(g=>targets.has(g)):[];return{key:pathKey(d.pathway_key||d.pathway_norm_key||d.pathway),pathway:String(d.pathway),scoreA:score,scoreB:NaN,delta:score,padj:num(d.padj,1),overlap,deOverlap:0,deGenes:[],genesA:genes,genesB:[],genes,tfA,tfB:[],tfUnion:tfA,tfCountA:tfA.length,tfCountB:0,tfCountUnion:tfA.length,conditionEnrichment:true}});if(tf)rows=rows.filter(r=>r.tfCountA>0).sort((a,b)=>b.tfCountA-a.tfCountA||b.scoreA-a.scoreA||a.padj-b.padj||a.pathway.localeCompare(b.pathway));else rows.sort((a,b)=>b.scoreA-a.scoreA||a.padj-b.padj||b.overlap-a.overlap||a.pathway.localeCompare(b.pathway));return rows}
function updatePathwaySelect(rows){const current=pathwaySelect.value,opts=[{value:'',label:'None'}].concat(rows.slice().sort((a,b)=>a.pathway.localeCompare(b.pathway,undefined,{sensitivity:'base'})).map(r=>({value:r.key,label:r.pathway})));fillSelect(pathwaySelect,opts,opts.some(x=>x.value===current)?current:'')}
function pathwayLabelClass(pathway){const key=pathKey(pathway),topics=new Set(),conditions=new Set();PATHWAYS.forEach(d=>{if(pathKey(d.pathway_key||d.pathway_norm_key||d.pathway)!==key)return;topics.add(topicNum(d.topic_num||d.topic));const condition=String(d.comparison_label||d.condition_id||d.comparison_id||'');if(condition)conditions.add(condition)});const topicSpecific=topics.size===1,conditionSpecific=conditions.size===1;return topicSpecific&&conditionSpecific?'pathLabel pathLabelBothSpecific':topicSpecific?'pathLabel pathLabelTopicSpecific':conditionSpecific?'pathLabel pathLabelConditionSpecific':'pathLabel'}
function wrapText(g,text,x,y,maxWidth,cls='pathLabel'){const value=String(text).trim(),t=el('text',{x,y,'text-anchor':'end',class:cls,style:'cursor:pointer'});t.textContent=value;g.appendChild(t);if(!value||t.getComputedTextLength()<=maxWidth)return{node:t,lines:1};const words=value.split(/\s+/),lines=[];let line='';words.forEach(word=>{const candidate=line?line+' '+word:word;t.textContent=candidate;if(line&&t.getComputedTextLength()>maxWidth){lines.push(line);line=word}else line=candidate});if(line)lines.push(line);t.replaceChildren();lines.forEach((s,i)=>{const q=el('tspan',{x,dy:i===0?0:'.84em'});q.textContent=s;t.appendChild(q)});return{node:t,lines:lines.length}}
function addPathLegend(g,sizeLabel){const size=el('text',{x:24,y:20,class:'pathLegend',fill:'#46524d'});size.textContent='Size: '+sizeLabel.replace(/^Dot size\s*=\s*/i,'');g.appendChild(size)}
function drawPathways(){
  PATH_ROWS=pathwayRows();
  updatePathwaySelect(PATH_ROWS);
  const g=byId('pathLayer');
  g.replaceChildren();
  const W=1120,H=900,LABEL_X=630,PLOT_L=680,R=54,PLOT_R=W-R,T=112,B=68,
    AXIS_Y=82,BOTTOM_AXIS_Y=H-B+8,overall=overallMode(),paired=!!cond2Select.value,
    rows=pageRows(PATH_ROWS,'path'),rowH=Math.min(22,(H-T-B)/Math.max(rows.length,1)),
    scoreLabel=!overall&&!paired&&rows.some(r=>r.conditionEnrichment)?'condition-topic enrichment combined score':pathwayScoreMethod.value==='rank_enrichment'?'ranked-expression enrichment z-score':'mean standardized-expression score',
    axisScoreLabel=pathwayScoreMethod.value==='rank_enrichment'?'RNA rank-enrichment delta':'mean-expression delta';
  addPathLegend(g,overall?'Dot size = topic-overlap genes':'Dot size = expressed overlap genes');
  if(paired)differentialConditionLabels(g,cond1Select.value,cond2Select.value,PLOT_L+(PLOT_R-PLOT_L)/4,PLOT_L+3*(PLOT_R-PLOT_L)/4,43);
  else if(!overall)differentialConditionLabels(g,cond1Select.value,'',PLOT_L,PLOT_R,43);
  if(!rows.length){
    const t=el('text',{x:30,y:112,'font-size':17,'font-weight':700});
    t.textContent=overall?'No significant overall pathways are available for this topic.':tfSelect.value?'No significant overall topic pathways contain direct targets of the selected TF.':'No significant overall topic pathways have at least three expression-matched genes.';
    g.appendChild(t);
    byId('pathStats').textContent='0 pathways';
    return;
  }
  const scores=PATH_ROWS.map(r=>r.delta).filter(Number.isFinite),maxAbs=Math.max(...scores.map(Math.abs),.05),
    domain=paired?[-maxAbs,maxAbs]:[Math.min(0,...scores),Math.max(...scores,.05)],
    span=Math.max(domain[1]-domain[0],1e-9),sx=v=>PLOT_L+(v-domain[0])/span*(PLOT_R-PLOT_L),
    formatTick=v=>Math.abs(v)>=1?v.toFixed(1):v.toFixed(2),
    maxLogp=Math.max(...PATH_ROWS.map(r=>-Math.log10(Math.max(r.padj,1e-300))),1),
    axisText=overall?'Overall topic enrichment combined score':paired?conditionLabel(cond1Select.value)+' <- '+axisScoreLabel+' -> '+conditionLabel(cond2Select.value):scoreLabel+': '+conditionLabel(cond1Select.value);
  rows.forEach((r,i)=>{const y=T+i*rowH;if(i%2===1)g.appendChild(el('rect',{x:8,y:y-1,width:W-16,height:rowH,fill:'#f7f8f6'}))});
  [0,.25,.5,.75,1].forEach(frac=>{
    const value=domain[0]+frac*span,xx=sx(value),zero=Math.abs(value)<span/1000;
    g.appendChild(el('line',{x1:xx,y1:AXIS_Y,x2:xx,y2:BOTTOM_AXIS_Y,stroke:zero?'#4b5550':'#8a9690','stroke-width':zero?1.3:.9,'stroke-dasharray':zero?'5 4':'2 5',opacity:frac===0||frac===1||zero?1:.3}));
    [AXIS_Y+16,BOTTOM_AXIS_Y+17].forEach(ty=>{const tick=el('text',{x:xx,y:ty,'text-anchor':'middle','font-size':11,'font-weight':700,fill:'#46524d',class:'pathAxisTickLabel'});tick.textContent=formatTick(value);g.appendChild(tick)});
  });
  [AXIS_Y,BOTTOM_AXIS_Y].forEach(y=>g.appendChild(el('line',{x1:PLOT_L,y1:y,x2:PLOT_R,y2:y,stroke:'#111','stroke-width':1.6})));
  [64,H-7].forEach((y,i)=>{const axis=el('text',{x:(PLOT_L+PLOT_R)/2,y,'text-anchor':'middle','font-size':i?13:12,'font-weight':700,fill:'#303834'});axis.textContent=axisText;g.appendChild(axis)});
  rows.forEach((r,i)=>{
    const y=T+i*rowH,yy=y+rowH/2,selected=r.key===pathwaySelect.value;
    if(selected)g.appendChild(el('rect',{x:8,y:y-1,width:W-16,height:rowH,rx:3,class:'pathRowSelected'}));
    const lab=wrapText(g,r.pathway,LABEL_X,yy,LABEL_X-16,pathwayLabelClass(r.pathway)),stars=pathwayStars(r.padj),
      star=el('text',{x:LABEL_X+7,y:yy+4,'text-anchor':'start','font-size':12,'font-weight':700,fill:'#9A6700',style:'cursor:pointer'}),
      x=sx(r.delta),rad=Math.min(12,3.2+Math.sqrt(Math.max(1,r.overlap))*.85),
      opacity=.5+.5*clamp(-Math.log10(Math.max(r.padj,1e-300))/maxLogp,0,1),
      dot=el('circle',{cx:x,cy:yy,r:rad,fill:overall?'#4E79A7':paired&&r.delta>0?conditionColor(2):conditionColor(1),opacity,stroke:selected?'#111':'#fff','stroke-width':selected?4:2,style:'cursor:pointer'}),
      countValue=tfSelect.value?(paired?r.tfCountUnion:r.tfCountA):0,countText=tfSelect.value?String(countValue):'',
      orderedA=overall?[]:orderedPathGenes(r.genes,cond1Select.value),orderedB=paired?orderedPathGenes(r.genes,cond2Select.value):[],
      msg=overall?r.pathway+'\nOverall topic-pathway q: '+r.padj.toExponential(2)+' '+stars+'\nCombined score: '+r.delta.toFixed(3)+'\nTopic-overlap genes: '+r.overlap+'\nGenes: '+r.genes.join(', '):r.pathway+'\n'+(r.conditionEnrichment?'Condition-topic pathway q: ':'Overall topic-pathway q: ')+r.padj.toExponential(2)+' '+stars+'\n'+conditionLabel(cond1Select.value)+' '+scoreLabel+': '+r.scoreA.toFixed(3)+(paired?'\n'+conditionLabel(cond2Select.value)+' '+scoreLabel+': '+r.scoreB.toFixed(3)+'\n'+conditionLabel(cond2Select.value)+' - '+conditionLabel(cond1Select.value)+' difference: '+r.delta.toFixed(3):'')+'\nExpression-matched genes: '+r.overlap+'\n'+conditionLabel(cond1Select.value)+' genes by expression: '+orderedA.map(d=>d.gene+' ('+d.expr.toFixed(2)+')').join(', ')+(paired?'\n'+conditionLabel(cond2Select.value)+' genes by expression: '+orderedB.map(d=>d.gene+' ('+d.expr.toFixed(2)+')').join(', '):'')+(tfSelect.value?'\nUnique '+tfLabel()+' target genes: '+countText+'\n'+conditionLabel(cond1Select.value)+' targets: '+(r.tfA||[]).join(', ')+(paired?'\n'+conditionLabel(cond2Select.value)+' targets: '+(r.tfB||[]).join(', ')+'\nUnion targets: '+(r.tfUnion||[]).join(', '):''):'');
    star.textContent=stars;
    const lineStep=9.2;
    lab.node.setAttribute('y',yy-((lab.lines-1)*lineStep)/2+3.5);
    [dot,lab.node,star].forEach(n=>{
      n.addEventListener('click',()=>{lastNetworkTrigger=n;pathwaySelect.value=r.key;sendState(false);if(!overall)openNetwork()});
      n.addEventListener('mousemove',ev=>tooltip('pathTooltip',ev,msg));
      n.addEventListener('mouseout',()=>hideTip('pathTooltip'));
    });
    g.appendChild(star);
    g.appendChild(dot);
    if(countText){const right=x>PLOT_R-75,ct=el('text',{x:x+(right?-rad-5:rad+5),y:yy+4,'text-anchor':right?'end':'start','font-size':11,'font-weight':700,fill:'#111',class:'pathTfCount',style:'pointer-events:none;paint-order:stroke;stroke:#fff;stroke-width:3px'});ct.textContent=countText;g.appendChild(ct)}
  });
  const start=PATH_ROWS.length?PAGE_INDEX.path*PAGE_SIZE.path+1:0;
  byId('pathStats').textContent=PATH_ROWS.length+' significant pathways';
}
function tfLabel(){const o=tfSelect.options[tfSelect.selectedIndex];return o?o.textContent:tfSelect.value}
function selectTf(tf){tfSelect.value=tfKey(tf);selectTopTfTopic(tf);showItemPage(conditionTopicRows(),'topic',d=>d.topic===Number(topicSelect.value));showItemPage(tfButterflyRows(),'tf',tfButterflyRowSelected);PAGE_INDEX.path=0;refresh();refreshConditionData()}
function currentPath(){return PATH_ROWS.find(r=>r.key===pathwaySelect.value)||null}
function rowPairExpressed(r){const cfg=directionConfig();return(Number.isFinite(r.gexprA)&&r.gexprA>=cfg.exprMin)||(Number.isFinite(r.gexprB)&&r.gexprB>=cfg.exprMin)}
function networkPairRows(){const allowed=topicGenes(),p=currentPath(),filterPath=p&&byId('networkPathwayFocus').value==='filter';if(filterPath){const pg=new Set(p.genes);Array.from(allowed).forEach(g=>{if(!pg.has(g))allowed.delete(g)})}let rows=edgePair().filter(r=>allowed.has(r.gene)&&rowPairExpressed(r)),tf=tfSelect.value;if(tf){const seeds=new Set(rows.filter(r=>r.tfu===tfKey(tf)).map(r=>r.gene));rows=rows.filter(r=>seeds.has(r.gene))}if(byId('networkTfScope').value==='topic'){const cut=num(byId('networkThetaCutoff').value,.3),conds=[cond1Select.value,cond2Select.value].filter(Boolean),primary=byId('networkPrimaryOnly').checked,topic=Number(topicSelect.value);rows=rows.filter(r=>r.tfu===tfKey(tf)||conds.some(c=>tfTheta(c,r.tfu,topic)>=cut&&(!primary||tfPrimaryTopic(c,r.tfu)===topic)))}if(cond2Select.value&&byId('networkCorrectDirectionOnly').checked&&byId('networkMode').value==='tf_gene')rows=rows.filter(r=>r.correctDirection);return rows}
function peakPairRows(allowedPairs){const c1=cond1Select.value,c2=cond2Select.value,m=new Map(),pairKeys=Array.from(new Set(allowedPairs.map(r=>r.tfu+'\t'+r.gene)));function add(cond,side){const cm=PEAK_BY_COND_PAIR.get(String(cond))||new Map();pairKeys.forEach(pk=>(cm.get(pk)||[]).forEach(d=>{const k=tfKey(d.tf)+'\t'+d.peak_id+'\t'+d.gene_key,r=m.get(k)||{tf:d.tf,tfu:tfKey(d.tf),peak:String(d.peak_id),gene:String(d.gene_key),a:0,b:0,exprA:NaN,exprB:NaN,gexprA:NaN,gexprB:NaN};if(side==='a'){r.a=num(d.fp_score);r.exprA=num(d.tf_expr,NaN);r.gexprA=num(d.gene_expr,NaN)}else{r.b=num(d.fp_score);r.exprB=num(d.tf_expr,NaN);r.gexprB=num(d.gene_expr,NaN)}m.set(k,r)}))}add(c1,'a');if(c2)add(c2,'b');let rows=Array.from(m.values()).map(r=>decorateLinkDirection(Object.assign(r,{exprA:tfExpression(c1,r.tfu,r.exprA),exprB:c2?tfExpression(c2,r.tfu,r.exprB):NaN,gexprA:auditGeneExpr(c1,r.gene,r.gexprA),gexprB:c2?auditGeneExpr(c2,r.gene,r.gexprB):NaN})));if(c2&&byId('networkCorrectDirectionOnly').checked)rows=rows.filter(r=>r.correctDirection);return rows}
function nodeLogfc(a,b){if(!cond2Select.value)return Math.log2(Math.max(0,Number.isFinite(a)?a:0)+1);if(!Number.isFinite(a)&&!Number.isFinite(b))return NaN;return Math.log2((Number.isFinite(a)?a:0)+1)-Math.log2((Number.isFinite(b)?b:0)+1)}
function linkRank(a,b){return(b.tfu===tfKey(tfSelect.value))-(a.tfu===tfKey(tfSelect.value))||(b.correctDirection?1:0)-(a.correctDirection?1:0)||Math.abs(b.fpChange??b.delta)-Math.abs(a.fpChange??a.delta)||String(a.tf).localeCompare(String(b.tf),undefined,{sensitivity:'base'})||String(a.gene).localeCompare(String(b.gene),undefined,{sensitivity:'base'})}
function buildNetwork(){let rows=networkPairRows(),topTf=Math.max(1,num(byId('networkTopTf').value,100)),topLinks=Math.max(1,num(byId('networkTopLinks').value,300)),rank=new Map();rows.forEach(r=>{const q=rank.get(r.tfu)||{genes:new Set(),correctGenes:new Set(),score:0};q.genes.add(r.gene);if(r.correctDirection)q.correctGenes.add(r.gene);q.score+=Math.abs(r.fpChange??r.delta);rank.set(r.tfu,q)});const tfs=Array.from(rank.entries()).sort((a,b)=>(b[0]===tfKey(tfSelect.value))-(a[0]===tfKey(tfSelect.value))||b[1].correctGenes.size-a[1].correctGenes.size||b[1].genes.size-a[1].genes.size||b[1].score-a[1].score).slice(0,topTf).map(x=>x[0]),keep=new Set(tfs);rows=rows.filter(r=>keep.has(r.tfu)).sort(linkRank).slice(0,topLinks);const mode=byId('networkMode').value,source=mode==='tf_peak_gene'?peakPairRows(rows).sort(linkRank).slice(0,topLinks):rows,nodes=new Map(),edges=[];function add(id,type,label,value){if(!nodes.has(id))nodes.set(id,{id,type,label,value,count:0,x:800,y:450});return nodes.get(id)}source.forEach(r=>{const tf='TF:'+r.tfu,gene='GENE:'+r.gene,tfv=nodeLogfc(r.exprA,r.exprB),gv=nodeLogfc(r.gexprA,r.gexprB),directionText=cond2Select.value?'\ntarget RNA log2FC: '+r.geneLfc.toFixed(3)+'\nTF RNA log2FC: '+r.tfLfc.toFixed(3)+'\ncorrect direction: '+(r.correctDirection?'yes':'no'):'';add(tf,'TF',r.tf,tfv);add(gene,'Gene',r.gene,gv);if(mode==='tf_peak_gene'){const peak='PEAK:'+r.peak;add(peak,'Peak',r.peak,r.delta);edges.push({from:tf,to:peak,value:r.delta,title:r.tf+' -> '+r.peak+'\ndelta FP: '+r.delta.toFixed(3)+directionText});edges.push({from:peak,to:gene,value:r.delta,title:r.peak+' -> '+r.gene+'\ndelta FP: '+r.delta.toFixed(3)+directionText})}else edges.push({from:tf,to:gene,value:r.delta,title:r.tf+' -> '+r.gene+'\nFP '+conditionLabel(cond1Select.value)+': '+r.a.toFixed(3)+(cond2Select.value?'\nFP '+conditionLabel(cond2Select.value)+': '+r.b.toFixed(3)+'\ndelta FP: '+r.delta.toFixed(3)+directionText:'')})});edges.forEach(e=>{nodes.get(e.from).count++;nodes.get(e.to).count++});return{nodes:Array.from(nodes.values()),edges,rows}}
const buildNetworkBase=buildNetwork;
buildNetwork=function(){const g=buildNetworkBase(),p=currentPath(),highlight=!!(p&&byId('networkPathwayFocus').value==='highlight'),hit=new Set();if(highlight){const genes=new Set(p.genes);g.nodes.filter(n=>n.type==='Gene'&&genes.has(n.label)).forEach(n=>hit.add(n.id));for(let pass=0;pass<3;pass++)g.edges.slice().reverse().forEach(e=>{if(hit.has(e.to)){e.pathwayHit=true;hit.add(e.from)}})}g.nodes.forEach(n=>n.pathwayHit=!highlight||hit.has(n.id));g.edges.forEach(e=>e.pathwayHit=!highlight||e.pathwayHit===true);return g}
function placeNetwork(g){const mode=byId('networkMode').value,requested=byId('networkLayout').value,layout=requested==='auto'?(tfSelect.value?'focused':'columns'):requested,types=t=>g.nodes.filter(n=>n.type===t).sort((a,b)=>b.count-a.count||a.label.localeCompare(b.label,undefined,{sensitivity:'base'})),column=(arr,x,top=70,bottom=830)=>arr.forEach((n,i)=>{n.x=x;n.y=top+(i+1)*((bottom-top)/Math.max(1,arr.length+1))});if(layout==='focused'){const selected='TF:'+tfKey(tfSelect.value),chosen=g.nodes.find(n=>n.id===selected),otherTfs=types('TF').filter(n=>n.id!==selected);if(mode==='tf_peak_gene'){const allTfs=chosen?[chosen].concat(otherTfs):otherTfs;column(allTfs,180);if(chosen)chosen.y=450;column(types('Peak'),610);column(types('Gene'),1110)}else{if(chosen){chosen.x=190;chosen.y=450}column(types('Gene'),790);column(otherTfs,1380)}return}if(['columns','bipartite','hierarchy'].includes(layout)){column(types('TF'),mode==='tf_peak_gene'?210:300);if(mode==='tf_peak_gene')column(types('Peak'),790);column(types('Gene'),mode==='tf_peak_gene'?1370:1270);return}const arr=g.nodes.slice().sort((a,b)=>a.type.localeCompare(b.type)||a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));if(layout==='grid'){const cols=Math.ceil(Math.sqrt(arr.length));arr.forEach((n,i)=>{n.x=120+(i%cols)*1360/Math.max(1,cols-1);n.y=90+Math.floor(i/cols)*720/Math.max(1,Math.ceil(arr.length/cols)-1)});return}if(layout==='circle'){arr.forEach((n,i)=>{const a=2*Math.PI*i/Math.max(1,arr.length)-Math.PI/2;n.x=800+360*Math.cos(a);n.y=450+360*Math.sin(a)});return}if(layout==='spiral'){arr.forEach((n,i)=>{const a=i*.62,r=35+330*i/Math.max(1,arr.length-1);n.x=800+r*Math.cos(a);n.y=450+r*Math.sin(a)});return}if(layout==='clustered'){const tfs=types('TF'),byIdMap=new Map(g.nodes.map(n=>[n.id,n]));tfs.forEach((n,i)=>{const a=2*Math.PI*i/Math.max(1,tfs.length)-Math.PI/2;n.x=800+300*Math.cos(a);n.y=450+300*Math.sin(a)});const owners=new Map();g.edges.forEach(e=>{if(e.from.startsWith('TF:')&&!owners.has(e.to))owners.set(e.to,e.from)});arr.filter(n=>n.type!=='TF').forEach((n,i)=>{const o=byIdMap.get(owners.get(n.id)),a=i*2.4;n.x=(o?o.x:800)+78*Math.cos(a);n.y=(o?o.y:450)+78*Math.sin(a)});return}const rings=layout==='radial'?{TF:100,Peak:265,Gene:395}:{TF:225,Peak:315,Gene:395};['TF','Peak','Gene'].forEach(type=>types(type).forEach((n,i,a)=>{const ang=2*Math.PI*i/Math.max(1,a.length)-Math.PI/2;n.x=800+rings[type]*Math.cos(ang);n.y=450+rings[type]*Math.sin(ang)}));if(layout==='force')forceNetwork(g)}
function forceNetwork(g){const by=new Map(g.nodes.map(n=>[n.id,n]));for(let z=0;z<120;z++){g.nodes.forEach(n=>{n.vx=n.vx||0;n.vy=n.vy||0});for(let i=0;i<g.nodes.length;i++)for(let j=i+1;j<g.nodes.length;j++){const a=g.nodes[i],b=g.nodes[j],dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(dx*dx+dy*dy,100),d=Math.sqrt(d2),f=Math.min(3,2800/d2);a.vx+=dx/d*f;a.vy+=dy/d*f;b.vx-=dx/d*f;b.vy-=dy/d*f}g.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to),dx=b.x-a.x,dy=b.y-a.y,d=Math.max(1,Math.hypot(dx,dy)),f=(d-125)*.01;a.vx+=dx/d*f;a.vy+=dy/d*f;b.vx-=dx/d*f;b.vy-=dy/d*f});g.nodes.forEach(n=>{n.vx+=(800-n.x)*.002;n.vy+=(450-n.y)*.002;n.x=clamp(n.x+n.vx,35,1565);n.y=clamp(n.y+n.vy,35,865);n.vx*=.75;n.vy*=.75})}}
function fitNetworkView(){const nodes=networkState.nodes;if(!nodes.length){networkState.view={x:0,y:0,k:1};applyNetworkView();return}let x0=Infinity,x1=-Infinity,y0=Infinity,y1=-Infinity;nodes.forEach(n=>{const hw=n.type==='TF'?n.boxW/2:Math.max(n.r,24),hh=n.type==='TF'?n.r:Math.max(n.r,18);x0=Math.min(x0,n.x-hw-22);x1=Math.max(x1,n.x+hw+Math.min(150,String(n.label).length*8)+22);y0=Math.min(y0,n.y-hh-22);y1=Math.max(y1,n.y+hh+22)});const margin=72,bw=Math.max(120,x1-x0),bh=Math.max(100,y1-y0),k=clamp(Math.min((1600-2*margin)/bw,(900-2*margin)/bh),.3,1.8),cx=(x0+x1)/2,cy=(y0+y1)/2;networkState.view={x:800-cx*k,y:450-cy*k,k};applyNetworkView()}
function drawNetwork(){const g=buildNetwork();networkState.nodes=g.nodes;networkState.edges=g.edges;placeNetwork(g);const spacing=clamp(num(byId('networkSpacingValue').value,1),.5,2);g.nodes.forEach(n=>{n.x=800+(n.x-800)*spacing;n.y=450+(n.y-450)*spacing});const eg=byId('networkEdges'),ng=byId('networkNodes'),lg=byId('networkLabelsLayer'),sel=byId('networkNodeSelect');eg.replaceChildren();ng.replaceChildren();lg.replaceChildren();sel.replaceChildren();const by=new Map(g.nodes.map(n=>[n.id,n])),vals=g.edges.map(e=>Math.abs(e.value)),elim=Math.max(...vals,1),nlim=Math.max(...g.nodes.map(n=>Math.abs(num(n.value))).filter(Number.isFinite),1),single=!cond2Select.value,tfp=byId('networkTfPalette').value,gp=byId('networkGenePalette').value,ep=byId('networkEdgePalette').value,tf0=num(byId('networkTfMin').value,14),tf1=num(byId('networkTfMax').value,20),tfCounts=g.nodes.filter(n=>n.type==='TF').map(n=>n.count),cmin=Math.min(...tfCounts,1),cmax=Math.max(...tfCounts,1),selectedTfId=tfSelect.value?'TF:'+tfKey(tfSelect.value):'';g.nodes.forEach(n=>{const selectedTf=n.id===selectedTfId;n.r=n.type==='TF'?(cmax===cmin?(tf0+tf1)/2:tf0+(n.count-cmin)/(cmax-cmin)*(tf1-tf0)):n.type==='Peak'?7:9;if(selectedTf)n.r+=4;n.boxW=n.type==='TF'?Math.max(n.r*3.2,Math.min(190,22+String(n.label).length*8)):0});if(!g.nodes.length){const empty=el('text',{x:800,y:450,'text-anchor':'middle','font-size':22,'font-weight':700,fill:'#4b5550'});empty.textContent='No network edges pass the current filters.';lg.appendChild(empty)}g.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to),dx=b.x-a.x,dy=b.y-a.y,d=Math.max(1,Math.hypot(dx,dy)),endR=b.type==='TF'?Math.max(b.r+3,b.boxW/2+3):b.r+4,focus=!selectedTfId||e.from===selectedTfId,line=el('line',{class:'edge',x1:a.x,y1:a.y,x2:b.x-dx/d*endR,y2:b.y-dy/d*endR,stroke:colorScale(e.value,elim,ep,single),'stroke-width':1.2+5.5*Math.sqrt(Math.abs(e.value)/elim),opacity:focus?.82:.2,'data-title':e.title,'data-from':e.from,'data-to':e.to});if(byId('networkArrows').checked)line.setAttribute('marker-end','url(#networkArrow)');eg.appendChild(line)});g.nodes.forEach(n=>{const fill=n.type==='Peak'?'#E99B18':Number.isFinite(n.value)?colorScale(n.value,nlim,n.type==='TF'?tfp:gp,single):'#BDBDBD',selectedTf=n.id===selectedTfId,shape=n.type==='TF'?el('rect',{class:'node'+(selectedTf?' selected':''),x:n.x-n.boxW/2,y:n.y-n.r*.8,width:n.boxW,height:n.r*1.6,rx:4,fill,stroke:selectedTf?'#111':'#fff','stroke-width':selectedTf?5:1.5}):el('circle',{class:'node',cx:n.x,cy:n.y,r:n.r,fill,stroke:'#fff','stroke-width':1.5});shape.dataset.id=n.id;shape.dataset.title=n.label+'\n'+n.type+(selectedTf?' (selected TF)':'')+'\nRNA value: '+(Number.isFinite(n.value)?n.value.toFixed(3):'NA')+'\nedges: '+n.count;ng.appendChild(shape);if(byId('networkLabels').checked){const contrast=n.type==='TF'?labelContrast(fill):{fill:'#171717',stroke:'#fff'},t=el('text',{class:'nodeLabel','data-id':n.id,x:n.type==='TF'?n.x:n.x+n.r+7,y:n.y+5,'text-anchor':n.type==='TF'?'middle':'start',fill:contrast.fill,stroke:contrast.stroke});t.textContent=n.label;lg.appendChild(t)}});g.nodes.slice().sort((a,b)=>a.label.localeCompare(b.label,undefined,{sensitivity:'base'})||a.id.localeCompare(b.id)).forEach(n=>{const o=document.createElement('option');o.value=n.id;o.textContent=n.label;sel.appendChild(o)});networkState.selected=selectedTfId&&by.has(selectedTfId)?selectedTfId:networkState.selected;if(networkState.selected)sel.value=networkState.selected;bindNetworkPointer();markSelected();fitNetworkView();byId('networkStats').textContent='Showing '+g.nodes.length+' nodes and '+g.edges.length+' edges. Node color: '+(single?'absolute log2 expression':'RNA log2FC')+'. Edge color and width: '+(single?'FP score':'Condition 1 - Condition 2 delta FP')+'.'}
const drawNetworkBase=drawNetwork;
function syncSingleColorControls(){[['networkTfPalette','networkTfSingleColor'],['networkGenePalette','networkGeneSingleColor'],['networkEdgePalette','networkEdgeSingleColor']].forEach(([palette,color])=>{byId(color).hidden=byId(palette).value!=='single'})}
drawNetwork=function(){let lo=clamp(num(byId('networkTfMin').value,14),6,40),hi=clamp(num(byId('networkTfMax').value,20),6,40);if(lo>hi)hi=lo;byId('networkTfMin').value=lo;byId('networkTfMax').value=hi;syncSingleColorControls();drawNetworkBase();const selectedTfId=tfSelect.value?'TF:'+tfKey(tfSelect.value):'',highlight=!!(currentPath()&&byId('networkPathwayFocus').value==='highlight'),nodeById=new Map(networkState.nodes.map(n=>[n.id,n]));networkState.nodes.forEach(n=>{if(n.type==='TF'){n.r=Math.min(40,n.r);n.boxW=Math.min(190,n.boxW)}const shape=byId('networkNodes').querySelector('[data-id="'+CSS.escape(n.id)+'"]');if(!shape)return;if(n.type==='TF'){shape.setAttribute('width',n.boxW);shape.setAttribute('height',n.r*1.6);shape.setAttribute('x',n.x-n.boxW/2);shape.setAttribute('y',n.y-n.r*.8)}const palette=n.type==='TF'?byId('networkTfPalette').value:n.type==='Gene'?byId('networkGenePalette').value:'';if(palette==='single')shape.setAttribute('fill',n.type==='TF'?byId('networkTfSingleColor').value:byId('networkGeneSingleColor').value);shape.setAttribute('opacity',!highlight||n.pathwayHit||n.id===selectedTfId?1:.16)});byId('networkLabelsLayer').querySelectorAll('[data-id]').forEach(label=>{const n=nodeById.get(label.dataset.id);if(n){const palette=n.type==='TF'?byId('networkTfPalette').value:n.type==='Gene'?byId('networkGenePalette').value:'';if(palette==='single'){const fill=n.type==='TF'?byId('networkTfSingleColor').value:byId('networkGeneSingleColor').value,contrast=labelContrast(fill);label.setAttribute('fill',contrast.fill);label.setAttribute('stroke',contrast.stroke)}}label.setAttribute('opacity',!highlight||!n||n.pathwayHit||n.id===selectedTfId?1:.16)});byId('networkEdges').querySelectorAll('[data-from]').forEach((line,i)=>{const edge=networkState.edges[i];if(byId('networkEdgePalette').value==='single')line.setAttribute('stroke',byId('networkEdgeSingleColor').value);if(highlight&&!edge.pathwayHit)line.setAttribute('opacity','.07')});redrawNetwork()}
function bindNetworkPointer(){const svg=byId('networkSvg'),ng=byId('networkNodes'),tip=byId('networkTooltip');ng.onmousedown=ev=>{if(!ev.target.classList.contains('node'))return;ev.stopPropagation();networkState.drag=networkState.nodes.find(n=>n.id===ev.target.dataset.id);networkState.selected=ev.target.dataset.id;byId('networkNodeSelect').value=networkState.selected;markSelected()};svg.onmousedown=ev=>{if(ev.target.classList&&ev.target.classList.contains('node'))return;const p=screenPoint(ev);networkState.pan={x:p.x,y:p.y,vx:networkState.view.x,vy:networkState.view.y}};svg.onmousemove=ev=>{if(networkState.drag){const p=worldPoint(ev);networkState.drag.x=clamp(p.x,30,1570);networkState.drag.y=clamp(p.y,30,870);redrawNetwork()}else if(networkState.pan){const p=screenPoint(ev);networkState.view.x=networkState.pan.vx+p.x-networkState.pan.x;networkState.view.y=networkState.pan.vy+p.y-networkState.pan.y;applyNetworkView()}tip.style.left=(ev.offsetX+12)+'px';tip.style.top=(ev.offsetY+12)+'px'};svg.onmouseup=()=>{networkState.drag=null;networkState.pan=null};svg.onmouseleave=()=>{networkState.drag=null;networkState.pan=null;tip.style.display='none'};svg.onmouseover=ev=>{const m=ev.target.dataset&&ev.target.dataset.title;if(m){tip.innerHTML=esc(m).replace(/\n/g,'<br/>');tip.style.display='block'}};svg.onmouseout=ev=>{if(ev.target.dataset&&ev.target.dataset.title)tip.style.display='none'};svg.onwheel=ev=>{ev.preventDefault();const p=screenPoint(ev),old=networkState.view.k,next=clamp(old*(ev.deltaY<0?1.15:1/1.15),.25,6);networkState.view.x=p.x-(p.x-networkState.view.x)*next/old;networkState.view.y=p.y-(p.y-networkState.view.y)*next/old;networkState.view.k=next;applyNetworkView()}}
function screenPoint(ev){const p=byId('networkSvg').createSVGPoint();p.x=ev.clientX;p.y=ev.clientY;return p.matrixTransform(byId('networkSvg').getScreenCTM().inverse())}function worldPoint(ev){const p=screenPoint(ev),v=networkState.view;return{x:(p.x-v.x)/v.k,y:(p.y-v.y)/v.k}}function applyNetworkView(){const v=networkState.view;byId('networkView').setAttribute('transform','translate('+v.x+' '+v.y+') scale('+v.k+')')}function markSelected(){byId('networkNodes').querySelectorAll('.node').forEach(x=>x.classList.toggle('selected',x.dataset.id===networkState.selected))}function redrawNetwork(){const by=new Map(networkState.nodes.map(n=>[n.id,n]));byId('networkNodes').querySelectorAll('.node').forEach(s=>{const n=by.get(s.dataset.id);if(!n)return;if(s.tagName==='rect'){s.setAttribute('x',n.x-n.boxW/2);s.setAttribute('y',n.y-n.r*.8)}else{s.setAttribute('cx',n.x);s.setAttribute('cy',n.y)}});byId('networkLabelsLayer').querySelectorAll('[data-id]').forEach(t=>{const n=by.get(t.dataset.id);if(n){t.setAttribute('x',n.type==='TF'?n.x:n.x+n.r+6);t.setAttribute('y',n.y+4)}});byId('networkEdges').querySelectorAll('[data-from]').forEach(line=>{const a=by.get(line.dataset.from),b=by.get(line.dataset.to);if(!a||!b)return;const dx=b.x-a.x,dy=b.y-a.y,d=Math.max(1,Math.hypot(dx,dy)),endR=b.type==='TF'?Math.max(b.r+3,b.boxW/2+3):b.r+3;line.setAttribute('x1',a.x);line.setAttribute('y1',a.y);line.setAttribute('x2',b.x-dx/d*endR);line.setAttribute('y2',b.y-dy/d*endR)})}
function syncReportMode(){byId('showPathwaysMode').classList.toggle('active',!networkOpen);byId('showGrnMode').classList.toggle('active',networkOpen);byId('showGrnMode').disabled=overallMode();byId('networkPathwayFocus').disabled=!currentPath()}
function updateNetworkHeading(){const path=currentPath(),title=path?path.pathway:'Topic '+topicSelect.value+' GRN',context=conditionLabel(cond1Select.value)+(cond2Select.value?' vs '+conditionLabel(cond2Select.value):'')+' | Topic '+topicSelect.value+(tfSelect.value?' | TF '+tfLabel():' | All TFs')+(path?' | Pathway highlighted':'');byId('networkTitle').textContent=title;byId('networkTitle').title=title;byId('networkContext').textContent=context;syncReportMode()}
function openNetwork(){if(overallMode())return;networkOpen=true;document.body.classList.add('network-open');byId('networkPanel').style.display='flex';updateNetworkHeading();byId('networkStats').textContent='Loading compressed network data...';networkState.view={x:0,y:0,k:1};setLoading(true,'Loading topic GRN...');const loads=[ensureSelectedConditionEdges(),loadNetworkPayload()];if(byId('networkMode').value==='tf_peak_gene')loads.push(ensureSelectedConditionPeaks());Promise.all(loads).then(()=>{if(networkOpen){updateNetworkHeading();drawNetwork()}}).catch(e=>{byId('networkStats').textContent=e.message}).finally(()=>setLoading(false));byId('networkBack').focus()}function closeNetwork(){networkOpen=false;document.body.classList.remove('network-open');byId('networkPanel').style.display='none';syncReportMode();if(lastNetworkTrigger&&typeof lastNetworkTrigger.focus==='function')lastNetworkTrigger.focus()}
function exportSvg(svg,name){const clone=svg.cloneNode(true);clone.setAttribute('xmlns',NS);clone.setAttribute('width',1600);clone.setAttribute('height',900);const blob=new Blob([new XMLSerializer().serializeToString(clone)],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();URL.revokeObjectURL(a.href)}
function updateContext(){const x=byId('selectedContext'),scope=overallMode()?'Overall topic pathways':conditionLabel(cond1Select.value)+(cond2Select.value?' vs '+conditionLabel(cond2Select.value):''),full='Selected: '+scope+' | Topic '+topicSelect.value+(tfSelect.value?' | TF '+tfLabel():'')+(currentPath()?' | '+currentPath().pathway:'');x.textContent=(overallMode()?'Overall | ':'')+'Topic '+topicSelect.value+(tfSelect.value?' | '+tfLabel():' | All TFs');x.title=full}
function updateProbabilityPanelTitles(){const label=tfSelect.value?' - '+tfLabel():'',rna=conditionTopicMetric&&conditionTopicMetric.value==='rna_delta',matched=byId('thetaAggregation').value==='matched'&&cond2Select.value,n=matched?matchedTfKeys().length:0;byId('conditionProbabilityTitle').textContent=rna?'Differential RNA Activity':'Topic Activity';byId('conditionProbabilityTitle').title=rna?'The same expressed assigned genes are compared in both conditions; positive log2 expression differences are summed and normalized by topic.':matched?'Equal-weight mean theta over the same '+n+' TF documents in both selected conditions.':'Equal-weight mean theta over every available condition::TF document in each condition.';byId('tfConditionProbabilityTitle').textContent='TF Probability'+label;byId('tfConditionProbabilityTitle').title=tfSelect.value?'Condition::TF topic probabilities for '+tfLabel():'TF probabilities for the selected topic'}
function refresh(){['mdsTooltip','activityTooltip','butterflyTooltip','tfButterflyTooltip','pathTooltip','networkTooltip'].forEach(hideTip);ensureConditionColors();const overall=overallMode();if(overall){cond2Select.value='';tfSelect.value='';pathwayDeOnly.checked=false;if(networkOpen)closeNetwork()}else if(cond2Select.value===cond1Select.value)cond2Select.value='';const paired=!!cond2Select.value;cond2Select.disabled=overall;cond1Color.disabled=overall;cond2Color.disabled=overall||!paired;byId('thetaAggregation').disabled=overall;conditionTopicMetric.disabled=overall;tfSelect.disabled=overall;pathwayScoreMethod.disabled=overall;pathwayDeOnly.disabled=!paired;byId('networkCorrectDirectionOnly').disabled=!paired;if(!pathwayDeOnly.dataset.bound){pathwayDeOnly.dataset.bound='1';pathwayDeOnly.addEventListener('change',()=>{pathwaySelect.value='';PAGE_INDEX.path=0;drawPathways();updateContext();sendState(false);if(networkOpen){updateNetworkHeading();drawNetwork()}})}const directionOnly=byId('networkCorrectDirectionOnly');if(!directionOnly.dataset.bound){directionOnly.dataset.bound='1';directionOnly.addEventListener('change',()=>{if(networkOpen)drawNetwork()})}updateProbabilityPanelTitles();drawMds();drawActivity();drawButterfly();drawTfButterfly();drawPathways();updateContext();syncReportMode();if(networkOpen){updateNetworkHeading();drawNetwork()}sendState(false)}
function refreshConditionData(){if(overallMode())return Promise.resolve(PAYLOAD);const tasks=[];if(!cond2Select.value)tasks.push(ensureSelectedConditionPathways());if(tfSelect.value||networkOpen){const edgeConditions=[cond1Select.value,cond2Select.value].filter(Boolean),needsEdges=edgeConditions.some(x=>!LOADED_EDGE_CONDITIONS.has(String(x)));if(needsEdges)byId('activityStats').textContent='Loading direct target data...';tasks.push(ensureSelectedConditionEdges())}if(!tasks.length)return Promise.resolve(PAYLOAD);return Promise.all(tasks).then(()=>refresh()).catch(e=>{byId('activityStats').textContent=e.message})}
function bindActivityZoom(){const svg=byId('activitySvg');let drag=null;svg.addEventListener('wheel',ev=>{ev.preventDefault();changeActivityZoom(ev.deltaY<0?1.25:.8)},{passive:false});svg.addEventListener('mousedown',ev=>{if(ACTIVITY_VIEW.k<=1)return;drag={x:ev.clientX,y:ev.clientY,ox:ACTIVITY_VIEW.ox,oy:ACTIVITY_VIEW.oy};svg.style.cursor='grabbing'});window.addEventListener('mousemove',ev=>{if(!drag)return;ACTIVITY_VIEW.ox=clamp(drag.ox-(ev.clientX-drag.x)/600/ACTIVITY_VIEW.k,-.5,.5);ACTIVITY_VIEW.oy=clamp(drag.oy+(ev.clientY-drag.y)/600/ACTIVITY_VIEW.k,-.5,.5);ACTIVITY_RENDER_KEY='';drawActivity()});window.addEventListener('mouseup',()=>{drag=null;svg.style.cursor=''});svg.addEventListener('dblclick',()=>{resetActivityZoom();drawActivity()});byId('activityZoomIn').onclick=()=>changeActivityZoom(1.5);byId('activityZoomOut').onclick=()=>changeActivityZoom(1/1.5)}
function applyExternalState(s,changed){if(!s)return;if(changed==='cond1')setActiveConditionSide(1);if(changed==='cond2')setActiveConditionSide(2);if(s.activeConditionSide)setActiveConditionSide(s.activeConditionSide);if(s.cond1!==undefined&&Array.from(cond1Select.options).some(o=>o.value===String(s.cond1)))cond1Select.value=String(s.cond1);if(s.cond2!==undefined&&Array.from(cond2Select.options).some(o=>o.value===String(s.cond2)))cond2Select.value=String(s.cond2);if(s.topic&&Array.from(topicSelect.options).some(o=>o.value===String(s.topic)))topicSelect.value=String(s.topic);if(s.tf!==undefined&&Array.from(tfSelect.options).some(o=>o.value===String(s.tf)))tfSelect.value=String(s.tf);if(s.metric&&Array.from(conditionTopicMetric.options).some(o=>o.value===String(s.metric)))conditionTopicMetric.value=String(s.metric);if(['method','k','condition','cond1','cond2','topic','tf','metric'].includes(changed))pathwaySelect.value='';if(changed==='condition'||changed==='cond1'||changed==='cond2'){resetPages();resetActivityZoom()}if(changed==='metric'){PAGE_INDEX.topic=0;PAGE_INDEX.path=0}if(changed==='tf'&&tfSelect.value)selectTopTfTopic(tfSelect.value);if(changed==='topic'){PAGE_INDEX.tf=0;PAGE_INDEX.path=0}if(changed==='tf'){PAGE_INDEX.tf=0;PAGE_INDEX.path=0}showItemPage(conditionTopicRows(),'topic',d=>d.topic===Number(topicSelect.value));if(tfSelect.value)showItemPage(tfButterflyRows(),'tf',tfButterflyRowSelected);refresh();if(['cond1','cond2','condition','tf'].includes(changed))refreshConditionData();if(s.pathway!==undefined&&Array.from(pathwaySelect.options).some(o=>o.value===String(s.pathway))){pathwaySelect.value=String(s.pathway);showItemPage(PATH_ROWS,'path',d=>d.key===pathwaySelect.value);drawPathways();updateContext()}if(changed==='pathway'&&pathwaySelect.value)openNetwork();else if(networkOpen){updateNetworkHeading();refreshConditionData()}}
function init(){initControls();setActiveConditionSide(1);bindActivityZoom();[[cond1Select,1],[cond2Select,2]].forEach(([x,side])=>{x.addEventListener('focus',()=>setActiveConditionSide(side));x.addEventListener('change',()=>{setActiveConditionSide(side);pathwaySelect.value='';resetPages();resetActivityZoom();refresh()})});topicSelect.addEventListener('change',()=>{pathwaySelect.value='';PAGE_INDEX.tf=0;PAGE_INDEX.path=0;showItemPage(conditionTopicRows(),'topic',d=>d.topic===Number(topicSelect.value));refresh()});tfSelect.addEventListener('change',()=>{if(tfSelect.value)selectTopTfTopic(tfSelect.value);pathwaySelect.value='';PAGE_INDEX.tf=0;PAGE_INDEX.path=0;showItemPage(conditionTopicRows(),'topic',d=>d.topic===Number(topicSelect.value));if(tfSelect.value)showItemPage(tfButterflyRows(),'tf',tfButterflyRowSelected);refresh()});pathwayScoreMethod.addEventListener('change',()=>{pathwaySelect.value='';PAGE_INDEX.path=0;drawPathways();updateContext();sendState(false);if(networkOpen){updateNetworkHeading();drawNetwork()}});pathwaySelect.addEventListener('change',()=>{showItemPage(PATH_ROWS,'path',d=>d.key===pathwaySelect.value);drawPathways();updateContext();sendState(false);if(pathwaySelect.value){lastNetworkTrigger=pathwaySelect;openNetwork()}else if(networkOpen){updateNetworkHeading();drawNetwork()}});byId('topicPrev').onclick=()=>movePage('topic',-1,drawButterfly);byId('topicNext').onclick=()=>movePage('topic',1,drawButterfly);byId('tfPrev').onclick=()=>movePage('tf',-1,drawTfButterfly);byId('tfNext').onclick=()=>movePage('tf',1,drawTfButterfly);byId('pathPrev').onclick=()=>movePage('path',-1,drawPathways);byId('pathNext').onclick=()=>movePage('path',1,drawPathways);['networkTfScope','networkThetaCutoff','networkPhiCutoff','networkPrimaryOnly','networkMode','networkTopTf','networkTopLinks','networkLayout','networkTfPalette','networkGenePalette','networkEdgePalette','networkTfSingleColor','networkGeneSingleColor','networkEdgeSingleColor','networkTfMin','networkTfMax','networkLabels','networkArrows','networkPathwayFocus'].forEach(id=>{const x=byId(id);x.addEventListener('change',drawNetwork);x.addEventListener('input',drawNetwork)});[['networkThetaPreset','networkThetaCutoff'],['networkPhiPreset','networkPhiCutoff']].forEach(([presetId,inputId])=>{const preset=byId(presetId),input=byId(inputId);preset.addEventListener('change',()=>{if(preset.value!=='custom')input.value=preset.value;drawNetwork()});input.addEventListener('input',()=>{const match=Array.from(preset.options).find(o=>o.value!=='custom'&&Number(o.value)===Number(input.value));preset.value=match?match.value:'custom'})});const spacing=byId('networkSpacing'),spacingValue=byId('networkSpacingValue');spacing.addEventListener('input',()=>{spacingValue.value=spacing.value;drawNetwork()});spacingValue.addEventListener('input',()=>{spacing.value=clamp(num(spacingValue.value,1),.5,2);drawNetwork()});document.querySelectorAll('[data-network-tab]').forEach(tab=>tab.addEventListener('click',()=>{document.querySelectorAll('[data-network-tab]').forEach(x=>x.classList.toggle('active',x===tab));document.querySelectorAll('[data-network-panel]').forEach(x=>x.classList.toggle('active',x.dataset.networkPanel===tab.dataset.networkTab))}));byId('showPathwaysMode').onclick=closeNetwork;byId('showGrnMode').onclick=()=>{lastNetworkTrigger=byId('showGrnMode');openNetwork()};byId('networkBack').onclick=closeNetwork;byId('networkFit').onclick=fitNetworkView;byId('networkReset').onclick=()=>{networkState.view={x:0,y:0,k:1};drawNetwork()};byId('networkExport').onclick=()=>exportSvg(byId('networkSvg'),'condition_topic_grn.svg');byId('exportSvgButton').onclick=()=>exportSvg(networkOpen?byId('networkSvg'):byId('pathSvg'),'condition_topic_report.svg');byId('networkNodeSelect').onchange=()=>{networkState.selected=byId('networkNodeSelect').value;markSelected()};document.addEventListener('keydown',ev=>{if(ev.key==='Escape'&&networkOpen)closeNetwork()});let resizeTimer=null;window.addEventListener('resize',()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(()=>{drawMds();drawActivity();drawButterfly();drawTfButterfly();drawPathways();if(networkOpen)fitNetworkView()},100)});window.addEventListener('message',ev=>{if(ev.data&&ev.data.type==='m3cr-set-state')applyExternalState(ev.data.state,ev.data.changed);if(ev.data&&ev.data.type==='m3cr-active-condition')setActiveConditionSide(ev.data.side);if(ev.data&&ev.data.type==='m3cr-action'&&ev.data.action==='export')byId('exportSvgButton').click()});refresh();sendState(true);refreshConditionData()}
[cond1Select,cond2Select,tfSelect].forEach(x=>x.addEventListener('change',refreshConditionData));
conditionTopicMetric.addEventListener('change',()=>{pathwaySelect.value='';PAGE_INDEX.topic=0;PAGE_INDEX.path=0;refresh()});
byId('thetaAggregation').addEventListener('change',()=>{pathwaySelect.value='';PAGE_INDEX.topic=0;PAGE_INDEX.path=0;refresh()});
byId('shortConditionNames').addEventListener('change',refresh);new MutationObserver(styleMdsPlot).observe(byId('mdsLayer'),{childList:true});
byId('networkMode').addEventListener('change',()=>{if(byId('networkMode').value!=='tf_peak_gene')return;byId('networkStats').textContent='Loading peak data...';ensureSelectedConditionPeaks().then(()=>{if(networkOpen)drawNetwork()}).catch(e=>{byId('networkStats').textContent=e.message})});
cond1Color.addEventListener('input',()=>{if(cond1Select.value)CONDITION_COLORS[cond1Select.value]=cond1Color.value;refresh()});cond2Color.addEventListener('input',()=>{if(cond2Select.value)CONDITION_COLORS[cond2Select.value]=cond2Color.value;refresh()});
window.addEventListener('message',ev=>{if(!ev.data||ev.data.type!=='m3cr-set-state')return;const s=ev.data.state||{};if(s.colors&&typeof s.colors==='object')Object.entries(s.colors).forEach(([k,v])=>{if(/^#[0-9a-f]{6}$/i.test(String(v)))CONDITION_COLORS[k]=v});if(cond1Select.value&&/^#[0-9a-f]{6}$/i.test(String(s.color1||'')))CONDITION_COLORS[cond1Select.value]=s.color1;if(cond2Select.value&&/^#[0-9a-f]{6}$/i.test(String(s.color2||'')))CONDITION_COLORS[cond2Select.value]=s.color2;refresh()});
if(new URLSearchParams(location.search).get('embed')==='1')document.body.classList.add('embed');loadReportData().then(x=>{GROUP_MDS=x.group_mds||[];GROUP_TOPIC=x.group_topic||[];TF_TOPIC=x.tf_topic||{columns:[],data:[]};PATHWAYS=x.pathways||[];indexTfTopics();return loadOverviewPayload()}).then(()=>{init();setLoading(false)}).catch(e=>{byId('activityStats').textContent=e.message;indexPayload();indexTfTopics();init();setLoading(false)});
)---"
}

.m3cr_condition_report_html <- function(title,
                                         group_mds,
                                         group_topic,
                                         tf_topic,
                                         pathways,
                                         out_html,
                                         condition_payload,
                                         report_state = list()) {
  compact_layout <- identical(basename(dirname(out_html)), "p") &&
    identical(basename(dirname(dirname(out_html))), "assets")
  legacy_layout <- identical(basename(dirname(out_html)), "pages") &&
    identical(basename(dirname(dirname(out_html))), "condition_topic_reports")
  standard_layout <- compact_layout || legacy_layout
  review_dir <- if (standard_layout) dirname(dirname(dirname(out_html))) else dirname(out_html)
  report_data_dir <- if (compact_layout) {
    file.path(review_dir, "assets", "cr")
  } else {
    file.path(review_dir, "assets", "condition_report_data")
  }
  dir.create(report_data_dir, recursive = TRUE, showWarnings = FALSE)
  legacy_report_data_file <- paste0(
    sub("[.]html$", "", basename(out_html)),
    "_data.js"
  )
  report_data_file <- .m3cr_report_data_file(out_html)
  unlink(file.path(report_data_dir, legacy_report_data_file), force = TRUE)
  tf_keep <- intersect(c(
    "comparison_label", "group_label", "tf", "tf_display", "tf_upper",
    "topic", "topic_num", "theta"
  ), names(tf_topic))
  packed <- lapply(list(
    group_mds = data.table::as.data.table(group_mds),
    group_topic = data.table::as.data.table(group_topic),
    tf_topic = data.table::as.data.table(tf_topic)[, ..tf_keep],
    pathways = data.table::as.data.table(pathways)
  ), .module2_report_browser_browser_payload_to_columnar)
  report_data_b64 <- .module2_report_browser_encode_browser_json_deflate_base64(packed)
  writeLines(paste0(
    "window.CRAFTGRN_CONDITION_REPORT_DATA=window.CRAFTGRN_CONDITION_REPORT_DATA||{};\n",
    "window.CRAFTGRN_CONDITION_REPORT_DATA[",
    jsonlite::toJSON(report_data_file, auto_unbox = TRUE),
    "]={compressed_columnar:'", report_data_b64, "'};\n"
  ), file.path(report_data_dir, report_data_file), useBytes = TRUE)
  report_data_spec <- list(
    file = report_data_file,
    base = if (compact_layout) {
      "../cr"
    } else if (legacy_layout) {
      "../../assets/condition_report_data"
    } else {
      "assets/condition_report_data"
    }
  )
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    "<link rel=\"icon\" href=\"data:,\"/>",
    "<meta name=\"craftgrn-module3-report-schema\" content=\"10\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    paste0("<style>", .m3cr_condition_report_css(), "</style></head><body>"),
    "<div class=\"top\"><label>Condition 1 <select id=\"cond1Select\"></select><input id=\"cond1Color\" type=\"color\" value=\"#D95F59\" title=\"Condition 1 color\" aria-label=\"Condition 1 color\"/></label><label>Condition 2 <select id=\"cond2Select\"></select><input id=\"cond2Color\" type=\"color\" value=\"#4E79A7\" title=\"Condition 2 color\" aria-label=\"Condition 2 color\"/></label><label>Theta TF set <select id=\"thetaAggregation\"><option value=\"matched\" selected>Matched TFs</option><option value=\"all\">Available TFs</option></select></label><label>Topic <select id=\"topicSelect\"></select></label><label>TF <select id=\"tfSelect\"></select></label><label>Pathway <select id=\"pathwaySelect\"></select></label><button id=\"exportSvgButton\" type=\"button\">Export SVG</button></div>",
    "<div id=\"loadingOverlay\" class=\"loadingOverlay\" role=\"status\" aria-live=\"polite\"><span class=\"loadingSpinner\" aria-hidden=\"true\"></span><span id=\"loadingMessage\">Loading report data...</span></div>",
    "<main class=\"grid\"><div class=\"left\"><section class=\"pane\"><div class=\"topPair\">",
    "<div class=\"mini\"><div class=\"paneHead\"><h2>Condition MDS</h2><div class=\"headTools\"><label title=\"Remove the prefix and suffix shared by all condition labels. Full names remain in hover text.\" style=\"display:inline-flex;align-items:center;gap:3px;font-size:11px;white-space:nowrap\"><input id=\"shortConditionNames\" type=\"checkbox\" checked/> Short names</label><span class=\"help\" title=\"MDS of mean document theta. Click a point or label to select a condition.\">?</span><span id=\"mdsStats\" class=\"meta\"></span></div></div><div class=\"body\"><div id=\"mdsTooltip\" class=\"tooltip\"></div><svg id=\"mdsSvg\" viewBox=\"0 0 760 760\"><g id=\"mdsLayer\"></g></svg></div></div>",
    "<div class=\"mini\"><div class=\"paneHead\"><h2>TF Activity</h2><div class=\"headTools\"><button id=\"activityZoomOut\" class=\"secondary zoomButton\" type=\"button\" title=\"Zoom out\" aria-label=\"Zoom out\">Zoom out</button><button id=\"activityZoomIn\" class=\"secondary zoomButton\" type=\"button\" title=\"Zoom in\" aria-label=\"Zoom in\">Zoom in</button><span class=\"help\" title=\"TF FP activity is normalized by total retained FP activity in each condition. Use the buttons or mouse wheel to zoom, drag to pan, and double-click to reset.\">?</span><span id=\"activityStats\" class=\"meta\"></span></div></div><div class=\"body\"><div id=\"activityTooltip\" class=\"tooltip\"></div><svg id=\"activitySvg\" viewBox=\"0 0 760 760\"><g id=\"activityLayer\"></g></svg></div></div>",
    "</div></section>",
    "<section class=\"pane\"><div class=\"bottomPair\"><div class=\"mini\"><div class=\"paneHead\"><h2 id=\"conditionProbabilityTitle\">Topic Activity</h2><div class=\"headTools\"><select id=\"conditionTopicMetric\" class=\"conditionTopicMetricControl\" title=\"Choose fitted model theta or pairwise differential RNA activity\"><option value=\"theta\">Model theta</option><option value=\"rna_delta\">Differential RNA activity</option></select><span class=\"help\" title=\"Model theta is the fitted topic mixture. Differential RNA activity uses the same assigned genes in both selected conditions and sums positive log2 expression differences by topic.\">?</span><span id=\"butterflyStats\" class=\"meta\"></span><span class=\"pager\"><button id=\"topicPrev\" class=\"secondary\" type=\"button\" title=\"Previous topics\" aria-label=\"Previous topics\">&lt;</button><span id=\"topicPageStatus\" class=\"pageStatus\"></span><button id=\"topicNext\" class=\"secondary\" type=\"button\" title=\"Next topics\" aria-label=\"Next topics\">&gt;</button></span></div></div><div class=\"body\"><div id=\"butterflyTooltip\" class=\"tooltip\"></div><div class=\"fixedChart\"><svg id=\"butterflySvg\" viewBox=\"0 0 760 760\"><g id=\"butterflyLayer\"></g></svg></div></div></div><div class=\"mini\"><div class=\"paneHead\"><h2 id=\"tfConditionProbabilityTitle\">TF Probability</h2><div class=\"headTools\"><span class=\"help\" title=\"With TF set to All, this shows the ten largest positive and negative TF probability differences for the selected topic. With a TF selected, it preserves the condition-topic order.\">?</span><span id=\"tfButterflyStats\" class=\"meta\"></span><span class=\"pager\"><button id=\"tfPrev\" class=\"secondary\" type=\"button\" title=\"Previous rows\" aria-label=\"Previous rows\">&lt;</button><span id=\"tfPageStatus\" class=\"pageStatus\"></span><button id=\"tfNext\" class=\"secondary\" type=\"button\" title=\"Next rows\" aria-label=\"Next rows\">&gt;</button></span></div></div><div class=\"body\"><div id=\"tfButterflyTooltip\" class=\"tooltip\"></div><div class=\"fixedChart\"><svg id=\"tfButterflySvg\" viewBox=\"0 0 760 760\"><g id=\"tfButterflyLayer\"></g></svg></div></div></div></div></section></div>",
    "<section class=\"pane\"><div class=\"paneHead\"><h2>Pathways</h2><div class=\"contextWrap\"><span class=\"reportModeToggle\" role=\"group\" aria-label=\"Pathway or GRN view\"><button id=\"showPathwaysMode\" class=\"active\" type=\"button\" title=\"Show the pathway plot.\">Pathways</button><button id=\"showGrnMode\" type=\"button\" title=\"Show the topic GRN. The selected pathway can highlight or filter it.\">GRN</button></span><label class=\"pathScoreControl\">Score <select id=\"pathwayScoreMethod\"><option value=\"rank_enrichment\" selected>Ranked expression</option><option value=\"mean_gene_zscore\">Mean standardized expression</option></select></label><label class=\"pathScoreControl\" title=\"Show only pathways containing at least one expressed target with pairwise absolute RNA log2FC at least the configured cutoff. This display filter does not change pathway scores, dot sizes, overlap genes, or Sub-GRNs.\"><input id=\"pathwayDeOnly\" type=\"checkbox\"/> DE targets only</label><span id=\"pathStats\" class=\"meta\"></span><span class=\"pager\"><button id=\"pathPrev\" class=\"secondary\" type=\"button\" title=\"Previous pathways\" aria-label=\"Previous pathways\">&lt;</button><span id=\"pathPageStatus\" class=\"pageStatus\"></span><button id=\"pathNext\" class=\"secondary\" type=\"button\" title=\"Next pathways\" aria-label=\"Next pathways\">&gt;</button></span><div id=\"selectedContext\" class=\"contextBadge\"></div><span class=\"help\" title=\"Asterisks show adjusted pathway significance. Dot size is the number of overlapping genes. Single-condition mode uses the saved per-condition enrichment table; pairwise mode compares expression over the overall topic-pathway genes.\">?</span></div></div><div class=\"body\"><div id=\"pathTooltip\" class=\"tooltip\"></div><div class=\"fixedChart\"><svg id=\"pathSvg\" viewBox=\"0 0 1120 900\"><g id=\"pathLayer\"></g></svg></div>",
    .m3cr_network_panel_html(),
    "</div></section></main>",
    "<script>",
    "let GROUP_MDS=[],GROUP_TOPIC=[],TF_TOPIC=[],PATHWAYS=[];",
    paste0("const REPORT_STATE=", .m3cr_payload_spec_json(report_state), ";"),
    paste0("const REPORT_DATA_SPEC=", .m3tb_json_for_html(report_data_spec), ";"),
    paste0("const CONDITION_PAYLOAD=", .m3cr_payload_spec_json(condition_payload), ";"),
    .m3cr_condition_report_js(),
    "</script></body></html>"
  )
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_html, useBytes = TRUE)
  out_html
}

.m3cr_write_combined_report_index <- function(out_file, reports, report_state = list()) {
  reports <- data.table::as.data.table(reports)
  if (!nrow(reports)) return(NA_character_)
  reports[, src := vapply(
    path,
    .m3tb_relative_html_path,
    character(1L),
    base_dir = dirname(out_file)
  )]
  if (!"method_setup" %in% names(reports)) reports[, method_setup := label]
  if (!"run_id" %in% names(reports)) reports[, run_id := method_setup]
  if (!"report_key" %in% names(reports)) reports[, report_key := paste0("k", k)]
  if (!"report_label" %in% names(reports)) reports[, report_label := paste0("K", k)]
  if (!"trained_k" %in% names(reports)) reports[, trained_k := k]
  if (!"topic_count" %in% names(reports)) reports[, topic_count := k]
  if (!"topic_space" %in% names(reports)) reports[, topic_space := "raw"]
  reports[, k := suppressWarnings(as.integer(k))]
  reports[, space_order := data.table::fifelse(topic_space == "raw", 1L, 2L)]
  data.table::setorder(reports, method_setup, trained_k, space_order, run_id)
  json <- .m3tb_json_for_html(reports[, .(
    label, k, trained_k, topic_count, topic_space, report_key, report_label,
    method_setup, run_id, src
  )])
  js <- r"---(
const REPORTS=REPORT_DATA,REPORT_STATE=REPORT_STATE_DATA,methodSelect=document.getElementById('methodSelect'),kSelect=document.getElementById('kSelect'),conditionTopicMetric=document.getElementById('conditionTopicMetric'),cond1Color=document.getElementById('cond1Color'),cond2Color=document.getElementById('cond2Color'),frame=document.getElementById('frame');
const childControls=['cond1','cond2','topic','tf','pathway'],PALETTE=['#4E79A7','#E15759','#59A14F','#F28E2B','#B07AA1','#76B7B2','#EDC948','#9C755F','#FF9DA7','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#1F78B4','#33A02C','#E31A1C','#6A3D9A','#B15928','#17BECF','#BCBD22','#7F7F7F'],INITIAL_CONDITION_COLORS=Object.assign({},REPORT_STATE.condition_colors||{}),conditionColors=Object.assign({},INITIAL_CONDITION_COLORS);
let activeConditionSide=1,pendingPathway=null;
function uniq(x){return Array.from(new Set(x))}function fill(sel,rows,current){const signature=rows.map(r=>String(r.value)+'\t'+String(r.label)).join('\n');if(sel.dataset.optionsSignature!==signature){sel.replaceChildren();rows.forEach(r=>{const o=document.createElement('option');o.value=String(r.value);o.textContent=String(r.label);sel.appendChild(o)});sel.dataset.optionsSignature=signature}if(rows.some(r=>String(r.value)===String(current)))sel.value=String(current||'')}
function methods(){return uniq(REPORTS.map(r=>r.method_setup)).sort((a,b)=>String(a).localeCompare(String(b),undefined,{sensitivity:'base'}))}function methodRows(){return REPORTS.filter(r=>r.method_setup===methodSelect.value)}function selectedReport(){const rows=methodRows().filter(r=>String(r.report_key)===String(kSelect.value));return rows[0]||methodRows()[0]||REPORTS[0]}
function embed(src){return src+(String(src).includes('?')?'&':'?')+'embed=1'}function load(){const hit=selectedReport();if(!hit)return;frame.src=embed(hit.src);document.title=hit.label}
function reportOptions(){return methodRows().slice().sort((a,b)=>Number(a.trained_k)-Number(b.trained_k)||(a.topic_space==='raw'?-1:1)).map(r=>({value:r.report_key,label:r.report_label,topic_space:r.topic_space,trained_k:r.trained_k}))}function defaultReportOption(options){const d=REPORT_STATE.defaults||{},space=String(d.topic_space||'combined'),trained=Number(d.trained_k||30);return options.find(x=>x.topic_space===space&&Number(x.trained_k)===trained)||options.find(x=>x.topic_space==='combined'&&Number(x.trained_k)===30)||options.find(x=>Number(x.trained_k)===trained)||options.find(x=>Number(x.trained_k)===30)||options[0]}function initRunControls(){fill(methodSelect,methods().map(x=>({value:x,label:x})),methods()[0]);const options=reportOptions();fill(kSelect,options,defaultReportOption(options)?.value);load()}
function saveColors(){}function selectedColor(id){const v=document.getElementById(id+'Select').value;return conditionColors[v]||'#7F7F7F'}function syncColorInputs(){cond1Color.value=selectedColor('cond1');cond2Color.value=selectedColor('cond2')}
function buildPalette(rows){const box=document.getElementById('conditionPaletteList'),ordered=rows.slice().sort((a,b)=>a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));ordered.forEach((r,i)=>{if(!/^#[0-9a-f]{6}$/i.test(String(conditionColors[r.value]||'')))conditionColors[r.value]=PALETTE[i%PALETTE.length]});const signature=ordered.map(r=>r.value+'\t'+r.label+'\t'+conditionColors[r.value]).join('\n');if(box.dataset.paletteSignature!==signature){box.replaceChildren();ordered.forEach(r=>{const label=document.createElement('label'),input=document.createElement('input'),text=document.createElement('span');input.type='color';input.value=conditionColors[r.value];input.setAttribute('aria-label',r.label+' color');text.textContent=r.label;input.addEventListener('input',()=>{conditionColors[r.value]=input.value;box.dataset.paletteSignature='';saveColors();syncColorInputs();sendState('colors')});label.append(input,text);box.appendChild(label)});box.dataset.paletteSignature=signature}saveColors();syncColorInputs()}
function setActiveConditionSide(side,notify=true){activeConditionSide=Number(side)===2?2:1;document.getElementById('cond1Select').classList.toggle('conditionTarget',activeConditionSide===1);document.getElementById('cond2Select').classList.toggle('conditionTarget',activeConditionSide===2);if(notify)frame.contentWindow.postMessage({type:'m3cr-active-condition',side:activeConditionSide},'*')}
function sendState(changed){frame.contentWindow.postMessage({type:'m3cr-set-state',changed,state:{cond1:document.getElementById('cond1Select').value,cond2:document.getElementById('cond2Select').value,activeConditionSide,color1:cond1Color.value,color2:cond2Color.value,colors:conditionColors,topic:document.getElementById('topicSelect').value,tf:document.getElementById('tfSelect').value,pathway:document.getElementById('pathwaySelect').value,metric:conditionTopicMetric.value}},'*')}
function syncChild(msg){const opts=msg.options||{},state=msg.state||{};childControls.forEach(id=>{const sel=document.getElementById(id+'Select'),rows=opts[id==='cond1'||id==='cond2'?'conditions':id==='topic'?'topics':id==='tf'?'tfs':'pathways']||[];let use=rows;if(id==='cond1')use=[{value:'',label:'Overall pathways'}].concat(rows.filter(x=>x.value));if(id==='cond2')use=[{value:'',label:'None'}].concat(rows.filter(x=>x.value));if(id==='tf')use=[{value:'',label:'All'}].concat(rows.filter(x=>x.value));if(id==='pathway')use=[{value:'',label:'None'}].concat(rows.filter(x=>x.value));const requested=id==='pathway'&&pendingPathway!==null&&use.some(x=>String(x.value)===String(pendingPathway))?pendingPathway:state[id];fill(sel,use,requested);if(id==='pathway'&&pendingPathway!==null&&String(state[id]||'')===String(pendingPathway))pendingPathway=null;sel.disabled=!use.length||(id==='cond2'&&!state.cond1)||(id==='tf'&&!state.cond1)});if(state.activeConditionSide)setActiveConditionSide(state.activeConditionSide,false);const conditions=opts.conditions||[];Object.entries(state.colors||{}).forEach(([k,v])=>{if(!conditionColors[k]&&/^#[0-9a-f]{6}$/i.test(String(v)))conditionColors[k]=v});buildPalette(conditions);document.getElementById('cond2Select').disabled=!state.cond1;document.getElementById('tfSelect').disabled=!state.cond1;cond1Color.disabled=!state.cond1;cond2Color.disabled=!state.cond2;conditionTopicMetric.disabled=!state.cond1}
methodSelect.addEventListener('change',()=>{pendingPathway=null;const options=reportOptions();fill(kSelect,options,defaultReportOption(options)?.value);load()});kSelect.addEventListener('change',()=>{pendingPathway=null;load()});conditionTopicMetric.addEventListener('change',()=>{pendingPathway=null;sendState('metric')});childControls.forEach(id=>document.getElementById(id+'Select').addEventListener('change',()=>{if(id==='pathway')pendingPathway=document.getElementById('pathwaySelect').value;else pendingPathway=null;if(id==='cond1'){setActiveConditionSide(1,false);const enabled=!!document.getElementById('cond1Select').value;document.getElementById('cond2Select').disabled=!enabled;document.getElementById('tfSelect').disabled=!enabled;conditionTopicMetric.disabled=!enabled}if(id==='cond2')setActiveConditionSide(2,false);syncColorInputs();sendState(id)}));[['cond1Select',1],['cond2Select',2]].forEach(([id,side])=>document.getElementById(id).addEventListener('focus',()=>setActiveConditionSide(side)));cond1Color.addEventListener('input',()=>{const c=document.getElementById('cond1Select').value;if(c)conditionColors[c]=cond1Color.value;buildPalette(Array.from(document.getElementById('cond1Select').options).filter(o=>o.value).map(o=>({value:o.value,label:o.textContent})));sendState('colors')});cond2Color.addEventListener('input',()=>{const c=document.getElementById('cond2Select').value;if(c)conditionColors[c]=cond2Color.value;buildPalette(Array.from(document.getElementById('cond1Select').options).filter(o=>o.value).map(o=>({value:o.value,label:o.textContent})));sendState('colors')});document.getElementById('paletteButton').onclick=()=>document.getElementById('conditionPalette').classList.toggle('open');document.getElementById('paletteReset').onclick=()=>{Object.keys(conditionColors).forEach(k=>delete conditionColors[k]);Object.assign(conditionColors,INITIAL_CONDITION_COLORS);buildPalette(Array.from(document.getElementById('cond1Select').options).filter(o=>o.value).map(o=>({value:o.value,label:o.textContent})));sendState('colors')};document.getElementById('exportButton').onclick=()=>frame.contentWindow.postMessage({type:'m3cr-action',action:'export'},'*');window.addEventListener('message',ev=>{if(ev.data&&(ev.data.type==='m3cr-ready'||ev.data.type==='m3cr-state')){syncChild(ev.data);if(ev.data.type==='m3cr-ready')sendState('colors')}});setActiveConditionSide(1,false);initRunControls();
)---"
  js <- sub("REPORT_DATA", json, js, fixed = TRUE)
  js <- sub("REPORT_STATE_DATA", .m3cr_payload_spec_json(report_state), js, fixed = TRUE)
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    "<link rel=\"icon\" href=\"data:,\"/>",
    "<title>Condition Topic Report</title>",
    "<style>html,body{width:100%;height:100%;overflow:hidden}*{box-sizing:border-box}body{margin:0;background:#f3f4f2;color:#171717;font-family:Arial,Helvetica,sans-serif;font-weight:700}.top{height:52px;display:flex;gap:10px;align-items:center;padding:6px 12px;border-bottom:1px solid #cfd2cc;background:#fff;white-space:nowrap;overflow:hidden}.controlGroup{display:flex;gap:8px;align-items:center;padding-right:10px;border-right:1px solid #dfe2dc}.controlGroup:last-of-type{border-right:0}.top label{display:flex;gap:5px;align-items:center;font-size:12px}.top select,.top button{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #9ca3a0;border-radius:4px;background:#fff;color:#171717;padding:6px 8px}.top select.conditionTarget{outline:3px solid #202322;outline-offset:1px}.top input[type=color]{width:30px;height:28px;padding:2px;border:1px solid #9ca3a0;border-radius:4px;background:#fff;cursor:pointer}.top button{background:#202322;color:#fff;border-color:#202322;cursor:pointer}.top select{max-width:215px}.conditionPalette{position:fixed;z-index:50;top:50px;left:300px;right:300px;display:none;background:#fff;border:1px solid #9ca3a0;box-shadow:0 8px 26px rgba(0,0,0,.18);padding:10px}.conditionPalette.open{display:block}.paletteHead{display:flex;justify-content:space-between;align-items:center;margin-bottom:8px}.paletteHead strong{font-size:13px}.paletteHead button{font:700 11px Arial,Helvetica,sans-serif;padding:4px 7px}.conditionPaletteList{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:5px 12px}.conditionPaletteList label{display:flex;align-items:center;gap:7px;min-width:0;font-size:11px}.conditionPaletteList input{width:27px;height:23px;padding:1px;flex:0 0 auto}.conditionPaletteList span{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}iframe{display:block;width:100vw;height:calc(100vh - 52px);border:0;background:#fff}@media(max-width:2000px){.top{gap:6px;padding-left:8px;padding-right:8px}.controlGroup{gap:5px;padding-right:6px}.top label{gap:3px;font-size:11px}.top select{max-width:165px;padding:5px 6px}.top label:has(#pathwaySelect) select{max-width:180px}.top button{padding:6px 7px}}@media(max-width:1500px){.top{height:64px;gap:7px}.controlGroup{gap:6px;padding-right:7px}.top label{flex-direction:column;align-items:flex-start;gap:1px;font-size:10px}.top select{max-width:180px;padding:4px 6px}.top label:has(#pathwaySelect) select{max-width:210px}.top button{padding:6px 8px}.conditionPalette{top:62px;left:80px;right:80px}iframe{height:calc(100vh - 64px)}}</style></head><body>",
    "<div class=\"top\"><div class=\"controlGroup\"><label>Method <select id=\"methodSelect\"></select></label><label>K <select id=\"kSelect\"></select></label></div><div class=\"controlGroup\"><label>Condition 1 <select id=\"cond1Select\"></select><input id=\"cond1Color\" type=\"color\" value=\"#D95F59\" title=\"Condition 1 color\" aria-label=\"Condition 1 color\"/></label><label>Condition 2 <select id=\"cond2Select\"></select><input id=\"cond2Color\" type=\"color\" value=\"#4E79A7\" title=\"Condition 2 color\" aria-label=\"Condition 2 color\"/></label><button id=\"paletteButton\" type=\"button\">Condition colors</button></div><div class=\"controlGroup\"><label>Topic view <select id=\"conditionTopicMetric\"><option value=\"theta\">Model theta</option><option value=\"rna_delta\">Differential RNA activity</option></select></label><label>Topic <select id=\"topicSelect\"></select></label><label>TF <select id=\"tfSelect\"></select></label><label>Pathway <select id=\"pathwaySelect\"></select></label></div><button id=\"exportButton\" type=\"button\">Export SVG</button></div><div id=\"conditionPalette\" class=\"conditionPalette\"><div class=\"paletteHead\"><strong>Condition colors</strong><button id=\"paletteReset\" type=\"button\">Reset colors</button></div><div id=\"conditionPaletteList\" class=\"conditionPaletteList\"></div></div><iframe id=\"frame\" title=\"Condition topic analysis\"></iframe>",
    "<script>", js, "</script></body></html>"
  )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_file, useBytes = TRUE)
  out_file
}
