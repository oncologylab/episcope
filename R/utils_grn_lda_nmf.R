# Utilities for GRN topic modeling (pooled LDA / pooled NMF) and marker overlap scoring.
# Intended to be package-ready (functions live under `R/`), while scripts can `source()` it.

standardize_delta_links_one_tbl <- function(df, source_csv) {
  stopifnot(is.data.frame(df), is.character(source_csv), length(source_csv) == 1L)

  file <- basename(source_csv)

  # Expect: "<case>_vs_<ctrl>_delta_links.csv"
  core <- sub("_delta_links\\.csv$", "", file, perl = TRUE)
  parts <- strsplit(core, "_vs_", fixed = TRUE)[[1]]

  if (length(parts) < 2L) {
    cli::cli_abort("Cannot parse case/ctrl tags from filename: {file}")
  }
  # If there were multiple "_vs_" tokens (unexpected), collapse the remainder as ctrl
  case_tag <- parts[[1]]
  ctrl_tag <- paste(parts[-1], collapse = "_vs_")

  if (!is.character(case_tag) || length(case_tag) != 1L) case_tag <- as.character(case_tag[[1]])
  if (!is.character(ctrl_tag) || length(ctrl_tag) != 1L) ctrl_tag <- as.character(ctrl_tag[[1]])

  # Helper: fetch column "<prefix><tag>", else NA_real_/NA
  get_by <- function(prefix, tag, na_value = NA_real_) {
    tag <- as.character(tag)
    if (length(tag) != 1L) tag <- tag[[1]]
    nm <- paste0(prefix, tag)
    if (length(nm) != 1L) nm <- nm[[1]]
    if (nm %in% names(df)) df[[nm]] else na_value
  }

  need <- c("tf", "gene_key", "peak_id", "delta_link_score")
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    cli::cli_abort("Missing required column(s) in {file}: {paste(miss, collapse=', ')}")
  }

  out <- tibble::tibble(
    tf       = as.character(df[["tf"]]),
    gene_key = as.character(df[["gene_key"]]),
    peak_id  = as.character(df[["peak_id"]]),

    comparison_id = core,
    case_tag      = case_tag,
    ctrl_tag      = ctrl_tag,
    source_csv    = source_csv,

    link_score_case     = as.numeric(get_by("link_score_", case_tag, na_value = NA_real_)),
    link_score_ctrl     = as.numeric(get_by("link_score_", ctrl_tag, na_value = NA_real_)),
    link_sign_case      = as.character(get_by("link_sign_", case_tag, na_value = NA_character_)),
    link_sign_ctrl      = as.character(get_by("link_sign_", ctrl_tag, na_value = NA_character_)),
    fp_score_case       = {
      v <- get_by("fp_score_", case_tag, na_value = NA_real_)
      if (all(is.na(v))) v <- get_by("fp_bed_score_", case_tag, na_value = NA_real_)
      as.numeric(v)
    },
    fp_score_ctrl       = {
      v <- get_by("fp_score_", ctrl_tag, na_value = NA_real_)
      if (all(is.na(v))) v <- get_by("fp_bed_score_", ctrl_tag, na_value = NA_real_)
      as.numeric(v)
    },
    tf_expr_case        = as.numeric(get_by("tf_expr_", case_tag, na_value = NA_real_)),
    tf_expr_ctrl        = as.numeric(get_by("tf_expr_", ctrl_tag, na_value = NA_real_)),
    gene_expr_case      = as.numeric(get_by("gene_expr_", case_tag, na_value = NA_real_)),
    gene_expr_ctrl      = as.numeric(get_by("gene_expr_", ctrl_tag, na_value = NA_real_)),
    active_case         = as.logical(get_by("active_", case_tag, na_value = NA)),
    active_ctrl         = as.logical(get_by("active_", ctrl_tag, na_value = NA)),
    fp_bound_case       = as.integer(get_by("fp_bound_", case_tag, na_value = NA_integer_)),
    fp_bound_ctrl       = as.integer(get_by("fp_bound_", ctrl_tag, na_value = NA_integer_)),
    tf_expr_flag_case   = as.integer(get_by("tf_expr_flag_", case_tag, na_value = NA_integer_)),
    tf_expr_flag_ctrl   = as.integer(get_by("tf_expr_flag_", ctrl_tag, na_value = NA_integer_)),
    gene_expr_flag_case = as.integer(get_by("gene_expr_flag_", case_tag, na_value = NA_integer_)),
    gene_expr_flag_ctrl = as.integer(get_by("gene_expr_flag_", ctrl_tag, na_value = NA_integer_)),
    r_tf_case           = as.numeric(get_by("r_tf_", case_tag, na_value = NA_real_)),
    r_tf_ctrl           = as.numeric(get_by("r_tf_", ctrl_tag, na_value = NA_real_)),

    delta_link_score    = as.numeric(df[["delta_link_score"]]),
    delta_fp_score      = if ("delta_fp_score" %in% names(df)) as.numeric(df[["delta_fp_score"]]) else if ("delta_fp_bed_score" %in% names(df)) as.numeric(df[["delta_fp_bed_score"]]) else NA_real_,
    log2FC_fp_score     = if ("log2FC_fp_score" %in% names(df)) as.numeric(df[["log2FC_fp_score"]]) else if ("log2FC_fp_bed_score" %in% names(df)) as.numeric(df[["log2FC_fp_bed_score"]]) else NA_real_,
    delta_tf_expr       = if ("delta_tf_expr"      %in% names(df)) as.numeric(df[["delta_tf_expr"]])      else NA_real_,
    delta_gene_expr     = if ("delta_gene_expr"    %in% names(df)) as.numeric(df[["delta_gene_expr"]])    else NA_real_,

    log2FC_tf_expr      = if ("log2FC_tf_expr"     %in% names(df)) as.numeric(df[["log2FC_tf_expr"]])     else NA_real_,
    log2FC_gene_expr    = if ("log2FC_gene_expr"   %in% names(df)) as.numeric(df[["log2FC_gene_expr"]])   else NA_real_,

    r_gene              = if ("r_gene"             %in% names(df)) as.numeric(df[["r_gene"]])             else NA_real_,
    p_gene              = if ("p_gene"             %in% names(df)) as.numeric(df[["p_gene"]])             else NA_real_,
    p_adj_gene          = if ("p_adj_gene"         %in% names(df)) as.numeric(df[["p_adj_gene"]])         else NA_real_,

    r_tf                = if ("r_tf"               %in% names(df)) as.numeric(df[["r_tf"]])               else NA_real_,
    p_tf                = if ("p_tf"               %in% names(df)) as.numeric(df[["p_tf"]])               else NA_real_,
    p_adj_tf            = if ("p_adj_tf"           %in% names(df)) as.numeric(df[["p_adj_tf"]])           else NA_real_
  )

  out
}

load_delta_links_all_tidy <- function(delta_csvs, use_fread = TRUE) {
  if (!is.character(delta_csvs) || length(delta_csvs) == 0L) {
    cli::cli_abort("`delta_csvs` must be a non-empty character vector.")
  }
  delta_csvs <- delta_csvs[file.exists(delta_csvs)]
  if (length(delta_csvs) == 0L) cli::cli_abort("No existing delta CSVs.")

  read_one <- function(p) {
    if (isTRUE(use_fread) && requireNamespace("data.table", quietly = TRUE)) {
      df <- data.table::fread(p, showProgress = FALSE)
      df <- tibble::as_tibble(df)
    } else {
      df <- readr::read_csv(p, show_col_types = FALSE)
      df <- tibble::as_tibble(df)
    }
    standardize_delta_links_one_tbl(df, source_csv = p)
  }

  dplyr::bind_rows(lapply(delta_csvs, read_one))
}

read_union_edge_keys_from_filtered <- function(filtered_only, use_fread = TRUE) {
  filtered_only <- filtered_only[file.exists(filtered_only)]
  if (length(filtered_only) == 0L) cli::cli_abort("No existing *_filtered.csv files provided.")

  if (isTRUE(use_fread) && requireNamespace("data.table", quietly = TRUE)) {
    dt_list <- lapply(filtered_only, function(p) {
      x <- data.table::fread(p, select = c("tf", "gene_key", "peak_id"), showProgress = FALSE)
      core <- sub("_delta_links_filtered\\.csv$", "", basename(p), perl = TRUE)
      x[, comparison_id := core]
      x
    })
    u <- unique(data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE))
    return(u[, .(tf, gene_key, peak_id, comparison_id)])
  }

  # fallback (slower)
  dt_list <- lapply(filtered_only, function(p) {
    x <- readr::read_csv(p, show_col_types = FALSE)[, c("tf", "gene_key", "peak_id")]
    core <- sub("_delta_links_filtered\\.csv$", "", basename(p), perl = TRUE)
    x$comparison_id <- core
    x
  })
  dplyr::distinct(dplyr::bind_rows(dt_list), tf, gene_key, peak_id, comparison_id)
}

filter_edges_all_tidy_by_union <- function(edges_all_tidy, union_keys) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for fast filtering at this scale.")
  }

  dt_edges <- data.table::as.data.table(edges_all_tidy)
  dt_keys  <- data.table::as.data.table(union_keys)

  need_e <- c("tf", "gene_key", "peak_id", "comparison_id")
  miss_e <- setdiff(need_e, names(dt_edges))
  if (length(miss_e)) cli::cli_abort("edges_all_tidy missing: {paste(miss_e, collapse=', ')}")

  need_k <- c("tf", "gene_key", "peak_id", "comparison_id")
  miss_k <- setdiff(need_k, names(dt_keys))
  if (length(miss_k)) cli::cli_abort("union_keys missing: {paste(miss_k, collapse=', ')}")

  data.table::setkeyv(dt_edges, need_e)
  data.table::setkeyv(dt_keys, need_k)

  out <- dt_edges[dt_keys, nomatch = 0L] # keep rows present in union
  tibble::as_tibble(out)
}

aggregate_edges_by_tf_gene <- function(edges_tbl,
                                       sum_cols = NULL,
                                       first_cols = NULL) {
  stopifnot(is.data.frame(edges_tbl))
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    cli::cli_abort("Please install dplyr for aggregation.")
  }
  if (!nrow(edges_tbl)) return(edges_tbl)

  need <- c("comparison_id", "tf", "gene_key")
  miss <- setdiff(need, names(edges_tbl))
  if (length(miss)) {
    cli::cli_abort("edges_tbl missing: {paste(miss, collapse=', ')}")
  }

  if (is.null(sum_cols)) sum_cols <- character(0)
  if (is.null(first_cols)) first_cols <- character(0)

  sum_cols <- intersect(sum_cols, names(edges_tbl))
  first_cols <- setdiff(intersect(first_cols, names(edges_tbl)), sum_cols)

  # Grouping keys are already retained; exclude them from first_cols.
  group_cols <- c("comparison_id", "tf", "gene_key")
  first_cols <- setdiff(first_cols, group_cols)

  edges_tbl |>
    dplyr::group_by(.data$comparison_id, .data$tf, .data$gene_key) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(sum_cols), ~ sum(.x, na.rm = TRUE)),
      dplyr::across(dplyr::all_of(first_cols), ~ dplyr::first(.x)),
      .groups = "drop"
    )
}

compute_log2fc_from_case_ctrl <- function(edges_tbl, eps = 1e-6) {
  stopifnot(is.data.frame(edges_tbl))
  need <- c("tf_expr_case", "tf_expr_ctrl", "gene_expr_case", "gene_expr_ctrl",
            "fp_score_case", "fp_score_ctrl")
  miss <- setdiff(need, names(edges_tbl))
  if (length(miss)) {
    cli::cli_abort("edges_tbl missing: {paste(miss, collapse=', ')}")
  }
  eps <- as.numeric(eps)
  if (!is.finite(eps) || eps <= 0) eps <- 1e-6

  out <- edges_tbl
  out$log2fc_tf_expr <- log2((out$tf_expr_case + eps) / (out$tf_expr_ctrl + eps))
  out$log2fc_gene_expr <- log2((out$gene_expr_case + eps) / (out$gene_expr_ctrl + eps))
  out$log2fc_fp_score <- log2((out$fp_score_case + eps) / (out$fp_score_ctrl + eps))
  out
}

build_pooled_doc_term <- function(edges_all_tidy,
                                  which = c("abs", "gain", "loss"),
                                  min_abs_delta = NULL,
                                  top_terms_per_doc = 500L,
                                  min_df = 2L,
                                  fp_col = NULL,
                                  fp_transform = c("abs", "pos", "neg"),
                                  fp_weight = 0,
                                  tf_expr_col = NULL,
                                  tf_expr_transform = c("abs", "pos", "neg"),
                                  tf_expr_weight = 0,
                                  gene_expr_col = NULL,
                                  gene_expr_transform = c("abs", "pos", "neg"),
                                  gene_expr_weight = 0,
                                  scale_quantile = 0.9) {
  which <- match.arg(which)
  fp_transform <- match.arg(fp_transform)
  tf_expr_transform <- match.arg(tf_expr_transform)
  gene_expr_transform <- match.arg(gene_expr_transform)

  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for pooled LDA at this scale.")
  }

  dt <- data.table::as.data.table(edges_all_tidy)

  apply_transform <- function(x, how) {
    x <- as.numeric(x)
    if (how == "abs") return(abs(x))
    if (how == "pos") return(pmax(x, 0))
    pmax(-x, 0)
  }

  scale01 <- function(x, q = 0.9) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- 0
    s <- stats::quantile(x, probs = q, na.rm = TRUE, names = FALSE, type = 7)
    s <- if (is.finite(s) && s > 0) s else 1
    pmin(1, x / s)
  }

  # weight
  if (which == "abs") {
    dt[, w := abs(delta_link_score)]
  } else if (which == "gain") {
    dt[, w := pmax(delta_link_score, 0)]
  } else {
    dt[, w := pmax(-delta_link_score, 0)]
  }

  dt <- dt[is.finite(w) & w > 0]

  if (!is.null(min_abs_delta)) {
    dt <- dt[abs(delta_link_score) >= min_abs_delta]
  }

  # Optional: reweight by TF activity / TF expression / gene expression.
  # These are applied multiplicatively as: w := w * (1 + weight * scaled_component)
  if (!is.null(fp_col) && fp_weight != 0) {
    if (!(fp_col %in% names(dt))) cli::cli_abort("fp_col '{fp_col}' not found in edges table.")
    comp <- apply_transform(dt[[fp_col]], fp_transform)
    dt[, w := w * (1 + fp_weight * scale01(comp, q = scale_quantile))]
  }
  if (!is.null(tf_expr_col) && tf_expr_weight != 0) {
    if (!(tf_expr_col %in% names(dt))) cli::cli_abort("tf_expr_col '{tf_expr_col}' not found in edges table.")
    comp <- apply_transform(dt[[tf_expr_col]], tf_expr_transform)
    dt[, w := w * (1 + tf_expr_weight * scale01(comp, q = scale_quantile))]
  }
  if (!is.null(gene_expr_col) && gene_expr_weight != 0) {
    if (!(gene_expr_col %in% names(dt))) cli::cli_abort("gene_expr_col '{gene_expr_col}' not found in edges table.")
    comp <- apply_transform(dt[[gene_expr_col]], gene_expr_transform)
    dt[, w := w * (1 + gene_expr_weight * scale01(comp, q = scale_quantile))]
  }

  dt <- dt[is.finite(w) & w > 0]

  # doc_id = comparison + TF
  dt[, doc_id := paste(comparison_id, tf, sep = "::")]

  # aggregate duplicate tf-gene within doc
  dt <- dt[, .(w = sum(w)), by = .(doc_id, comparison_id, tf, gene_key)]

  # keep top terms per doc to control size
  if (!is.null(top_terms_per_doc) && is.finite(top_terms_per_doc) && top_terms_per_doc > 0) {
    data.table::setorder(dt, doc_id, -w)
    dt <- dt[, head(.SD, top_terms_per_doc), by = doc_id]
  }

  # min_df across docs
  df_tbl <- unique(dt[, .(doc_id, gene_key)])
  term_df <- df_tbl[, .N, by = gene_key]
  keep_terms <- term_df[N >= min_df, gene_key]
  dt <- dt[gene_key %in% keep_terms]

  dt[]
}

fit_pooled_lda <- function(doc_term_dt,
                           K = 20,
                           method = c("VEM", "Gibbs"),
                           seed = 1L,
                           gamma_cutoff = 0.2,
                           multi_sep = ";") {
  fit_pooled_lda_generic(
    doc_term_dt = doc_term_dt,
    term_col = "gene_key",
    meta_cols = c("comparison_id", "tf"),
    K = K,
    method = method,
    seed = seed,
    gamma_cutoff = gamma_cutoff,
    multi_sep = multi_sep
  )
}

fit_pooled_lda_generic <- function(doc_term_dt,
                                   term_col,
                                   meta_cols = character(0),
                                   K = 20,
                                   method = c("VEM", "Gibbs"),
                                   seed = 1L,
                                   gamma_cutoff = 0.2,
                                   multi_sep = ";") {
  method <- match.arg(method)

  if (!requireNamespace("slam", quietly = TRUE) ||
      !requireNamespace("topicmodels", quietly = TRUE)) {
    cli::cli_abort("Need packages: slam, topicmodels")
  }

  df <- as.data.frame(doc_term_dt)

  need <- c("doc_id", "w", term_col, meta_cols)
  miss <- setdiff(need, names(df))
  if (length(miss)) cli::cli_abort("doc_term_dt missing columns: {paste(miss, collapse=', ')}")

  # slam::simple_triplet_matrix does not allow duplicate (doc, term) pairs.
  # Aggregate any duplicates up-front.
  df_fit <- df[, c("doc_id", term_col, "w"), drop = FALSE]
  names(df_fit)[names(df_fit) == term_col] <- "term"
  df_fit$doc_id <- as.character(df_fit$doc_id)
  df_fit$term <- as.character(df_fit$term)
  df_fit$w <- as.numeric(df_fit$w)
  df_fit <- df_fit[is.finite(df_fit$w) & df_fit$w > 0 & !is.na(df_fit$doc_id) & df_fit$doc_id != "" & !is.na(df_fit$term) & df_fit$term != "", , drop = FALSE]
  if (!nrow(df_fit)) cli::cli_abort("After filtering, doc_term_dt has no positive finite weights.")

  df_fit <- stats::aggregate(w ~ doc_id + term, data = df_fit, FUN = sum)

  docs  <- sort(unique(df_fit$doc_id))
  terms <- sort(unique(df_fit$term))

  i <- match(df_fit$doc_id, docs)
  j <- match(df_fit$term, terms)

  s <- stats::quantile(df_fit$w, 0.9, na.rm = TRUE)
  s <- if (is.finite(s) && s > 0) s else 1
  v <- pmax(1L, as.integer(round(10 * df_fit$w / s)))

  stm <- slam::simple_triplet_matrix(
    i = i, j = j, v = v,
    nrow = length(docs), ncol = length(terms),
    dimnames = list(docs, terms)
  )

  ctrl <- if (method == "VEM") list(seed = seed) else list(seed = seed, burnin = 2000, iter = 2000, thin = 100)
  fit  <- topicmodels::LDA(stm, k = K, method = method, control = ctrl)
  post <- topicmodels::posterior(fit)

  gamma <- as.matrix(post$topics) # n_docs x K
  doc_id <- rownames(gamma)
  n_doc <- nrow(gamma)
  n_topic <- ncol(gamma)

  gamma_long <- tibble::tibble(
    doc_id = rep(doc_id, times = n_topic),
    topic  = rep(seq_len(n_topic), each = n_doc),
    gamma  = as.vector(gamma)
  )

  if (length(meta_cols)) {
    meta <- unique(df[, c("doc_id", meta_cols), drop = FALSE])
    meta <- tibble::as_tibble(meta)
    gamma_long <- dplyr::left_join(gamma_long, meta, by = "doc_id")
  } else {
    meta <- tibble::tibble(doc_id = unique(df$doc_id))
  }

  keep <- gamma_long[gamma_long$gamma >= gamma_cutoff, , drop = FALSE]
  if (nrow(keep)) {
    keep <- keep[order(keep$doc_id, -keep$gamma, keep$topic), , drop = FALSE]
    split_topics <- split(keep$topic, keep$doc_id)
    topics_str <- vapply(split_topics, function(v) paste(v, collapse = multi_sep), character(1))
    assign_tbl <- tibble::tibble(doc_id = names(topics_str), topics_str = unname(topics_str))
  } else {
    assign_tbl <- tibble::tibble(doc_id = character(0), topics_str = character(0))
  }

  doc_topics <- dplyr::left_join(meta, assign_tbl, by = "doc_id")

  list(
    fit = fit,
    doc_topic_long = gamma_long,
    doc_topics = doc_topics,
    term_col = term_col,
    meta_cols = meta_cols
  )
}

compute_topic_coherence_umass <- function(doc_term_dt,
                                          topic_term,
                                          top_n = 20L,
                                          smooth = 1L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for topic coherence.")
  }
  dt <- data.table::as.data.table(doc_term_dt)
  need <- c("doc_id", "gene_key", "w")
  miss <- setdiff(need, names(dt))
  if (length(miss)) cli::cli_abort("doc_term_dt missing columns: {paste(miss, collapse=', ')}")

  dt <- dt[is.finite(w) & w > 0]
  dt <- dt[!is.na(doc_id) & doc_id != "" & !is.na(gene_key) & gene_key != ""]
  if (!nrow(dt)) {
    return(tibble::tibble(topic = integer(0), topic_coherence = numeric(0), n_terms_used = integer(0)))
  }

  doc_ids <- unique(dt$doc_id)
  dt[, doc_idx := as.integer(factor(doc_id, levels = doc_ids))]

  term_docs <- dt[, .(docs = list(unique(doc_idx))), by = gene_key]
  docs_list <- setNames(term_docs$docs, term_docs$gene_key)

  term_names <- colnames(topic_term)
  if (is.null(term_names)) {
    cli::cli_abort("topic_term must have colnames matching gene_key.")
  }

  n_topics <- nrow(topic_term)
  out <- vector("list", n_topics)

  for (k in seq_len(n_topics)) {
    weights <- topic_term[k, ]
    ord <- order(weights, decreasing = TRUE, na.last = NA)
    terms <- term_names[ord]
    if (length(terms) > top_n) terms <- terms[seq_len(top_n)]
    terms <- terms[terms %in% names(docs_list)]
    n_terms_used <- length(terms)

    if (n_terms_used < 2L) {
      out[[k]] <- tibble::tibble(topic = k, topic_coherence = NA_real_, n_terms_used = n_terms_used)
      next
    }

    sum_coh <- 0
    n_pairs <- 0L
    for (i in 2:n_terms_used) {
      docs_i <- docs_list[[terms[i]]]
      for (j in 1:(i - 1)) {
        docs_j <- docs_list[[terms[j]]]
        d_j <- length(docs_j)
        if (d_j == 0L) next
        d_ij <- length(intersect(docs_i, docs_j))
        sum_coh <- sum_coh + log((d_ij + smooth) / d_j)
        n_pairs <- n_pairs + 1L
      }
    }

    coh <- if (n_pairs > 0L) sum_coh / n_pairs else NA_real_
    out[[k]] <- tibble::tibble(topic = k, topic_coherence = coh, n_terms_used = n_terms_used)
  }

  dplyr::bind_rows(out)
}

compute_lda_topic_coherence <- function(doc_term_dt,
                                        lda_fit,
                                        top_n = 20L,
                                        smooth = 1L) {
  if (!requireNamespace("topicmodels", quietly = TRUE)) {
    cli::cli_abort("Please install topicmodels for LDA coherence.")
  }
  topic_term <- topicmodels::posterior(lda_fit$fit)$terms
  compute_topic_coherence_umass(
    doc_term_dt = doc_term_dt,
    topic_term = topic_term,
    top_n = top_n,
    smooth = smooth
  )
}

build_gene_centric_doc_term <- function(edges_all_tidy,
                                        doc = c("gene", "comparison_gene"),
                                        term = c("tf", "peak_id"),
                                        which = c("abs", "gain", "loss"),
                                        min_abs_delta = NULL,
                                        top_terms_per_doc = 500L,
                                        min_df = 2L) {
  doc <- match.arg(doc)
  term <- match.arg(term)
  which <- match.arg(which)

  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for gene-centric LDA at this scale.")
  }

  dt <- data.table::as.data.table(edges_all_tidy)
  need <- c("comparison_id", "tf", "gene_key", "delta_link_score")
  miss <- setdiff(need, names(dt))
  if (length(miss)) cli::cli_abort("edges_all_tidy missing: {paste(miss, collapse=', ')}")

  if (which == "abs") {
    dt[, w := abs(delta_link_score)]
  } else if (which == "gain") {
    dt[, w := pmax(delta_link_score, 0)]
  } else {
    dt[, w := pmax(-delta_link_score, 0)]
  }

  dt <- dt[is.finite(w) & w > 0]
  if (!is.null(min_abs_delta)) dt <- dt[abs(delta_link_score) >= min_abs_delta]

  dt[, doc_id := if (doc == "gene") as.character(gene_key) else paste(comparison_id, gene_key, sep = "::")]
  dt[, term_id := if (term == "tf") as.character(tf) else as.character(peak_id)]

  dt <- dt[, .(w = sum(w)), by = .(doc_id, comparison_id, gene_key, term_id)]

  if (!is.null(top_terms_per_doc) && is.finite(top_terms_per_doc) && top_terms_per_doc > 0) {
    data.table::setorder(dt, doc_id, -w)
    dt <- dt[, head(.SD, top_terms_per_doc), by = doc_id]
  }

  df_tbl <- unique(dt[, .(doc_id, term_id)])
  term_df <- df_tbl[, .N, by = term_id]
  keep_terms <- term_df[N >= min_df, term_id]
  dt <- dt[term_id %in% keep_terms]

  # Return with the "term" column name matching the selected term type
  if (term == "tf") {
    data.table::setnames(dt, "term_id", "tf")
  } else {
    data.table::setnames(dt, "term_id", "peak_id")
  }

  dt[]
}

topic_marker_overlap_from_doc_topics <- function(doc_topics,
                                                entity_col,
                                                markers_epi,
                                                markers_mes,
                                                topics_sep = ";",
                                                universe_entities = NULL,
                                                top_n_overlap_entities = 30L,
                                                verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  if (!("topics_str" %in% names(doc_topics))) cli::cli_abort("doc_topics missing column: topics_str")
  if (!(entity_col %in% names(doc_topics))) cli::cli_abort("doc_topics missing entity_col: {entity_col}")

  markers_epi <- unique(as.character(markers_epi))
  markers_mes <- unique(as.character(markers_mes))

  dt_doc <- data.table::as.data.table(doc_topics)[, .(entity = as.character(get(entity_col)), topics_str = as.character(topics_str))]
  dt_doc <- dt_doc[!is.na(topics_str) & topics_str != "" & !is.na(entity) & entity != ""]
  if (!nrow(dt_doc)) cli::cli_abort("No doc_topics rows have topics_str.")

  dt_long <- dt_doc[, .(topic = unlist(strsplit(topics_str, split = topics_sep, fixed = TRUE))), by = entity]
  dt_long[, topic := as.integer(gsub("[^0-9]", "", as.character(topic)))]
  dt_long <- dt_long[is.finite(topic)]

  topic_entities <- unique(dt_long[, .(topic, entity)])

  if (is.null(universe_entities)) {
    universe_entities <- sort(unique(topic_entities$entity))
  } else {
    universe_entities <- unique(as.character(universe_entities))
  }
  N <- length(universe_entities)

  calc_stats <- function(topic_entity_vec, marker_vec) {
    topic_entity_vec <- unique(topic_entity_vec)
    marker_vec <- unique(marker_vec)

    topic_u  <- intersect(topic_entity_vec, universe_entities)
    marker_u <- intersect(marker_vec, universe_entities)

    a <- length(intersect(topic_u, marker_u))
    b <- length(setdiff(topic_u, marker_u))
    c <- length(setdiff(marker_u, topic_u))
    d <- N - a - b - c
    d <- max(d, 0)

    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

    list(
      overlap_n = a,
      marker_coverage = if (length(marker_u) > 0) a / length(marker_u) else NA_real_,
      fisher_p = unname(ft$p.value),
      fisher_or = unname(ft$estimate)
    )
  }

  topics <- sort(unique(topic_entities$topic))
  sizes <- topic_entities[, .(n_entities = data.table::uniqueN(entity)), by = topic]

  res_list <- lapply(topics, function(k) {
    ent_k <- topic_entities[topic == k, entity]
    epi <- calc_stats(ent_k, markers_epi)
    mes <- calc_stats(ent_k, markers_mes)

    overlap_epi <- intersect(unique(ent_k), markers_epi)
    overlap_mes <- intersect(unique(ent_k), markers_mes)

    tibble::tibble(
      topic = k,
      overlap_epi_n = epi$overlap_n,
      epi_marker_coverage = epi$marker_coverage,
      epi_fisher_p = epi$fisher_p,
      epi_fisher_or = as.numeric(epi$fisher_or),
      overlap_epi_entities = paste(utils::head(sort(overlap_epi), top_n_overlap_entities), collapse = ";"),
      overlap_mes_n = mes$overlap_n,
      mes_marker_coverage = mes$marker_coverage,
      mes_fisher_p = mes$fisher_p,
      mes_fisher_or = as.numeric(mes$fisher_or),
      overlap_mes_entities = paste(utils::head(sort(overlap_mes), top_n_overlap_entities), collapse = ";")
    )
  })

  stats_tbl <- dplyr::left_join(dplyr::bind_rows(res_list), tibble::as_tibble(sizes), by = "topic")

  best_epi <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_epi_n), epi_fisher_p, dplyr::desc(epi_marker_coverage)) |>
    dplyr::slice(1)

  best_mes <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_mes_n), mes_fisher_p, dplyr::desc(mes_marker_coverage)) |>
    dplyr::slice(1)

  if (isTRUE(verbose)) {
    message("[topic_marker_overlap_from_doc_topics] Universe N=", N, "; topics=", length(topics), "; entity_col=", entity_col)
    message("[topic_marker_overlap_from_doc_topics] Best EPI topic: ", best_epi$topic,
            " (overlap=", best_epi$overlap_epi_n,
            ", coverage=", signif(best_epi$epi_marker_coverage, 3),
            ", fisher_p=", signif(best_epi$epi_fisher_p, 3), ")")
    message("[topic_marker_overlap_from_doc_topics] Best MES topic: ", best_mes$topic,
            " (overlap=", best_mes$overlap_mes_n,
            ", coverage=", signif(best_mes$mes_marker_coverage, 3),
            ", fisher_p=", signif(best_mes$mes_fisher_p, 3), ")")
  }

  list(
    stats = stats_tbl,
    best_epi = best_epi,
    best_mes = best_mes,
    topic_entity = tibble::as_tibble(topic_entities)
  )
}

fit_pooled_nmf <- function(doc_term_dt,
                           K = 20,
                           method = c("brunet", "lee", "nsNMF"),
                           seed = 1L,
                           nrun = 10L,
                           gamma_cutoff = 0.2,
                           doc_topic_top_n = NULL,
                           multi_sep = ";",
                           allow_dense_fallback = TRUE,
                           dense_max_gb = 2,
                           nmf_options = "-p",
                           nmf_pbackend = NA,
                           nmf_stop = c("stationary", "connectivity", "none")) {
  method <- match.arg(method)
  nmf_stop <- match.arg(nmf_stop)

  if (!requireNamespace("NMF", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    cli::cli_abort("Need packages: NMF, Matrix")
  }

  stop_fun <- switch(
    nmf_stop,
    stationary = NMF::nmf.stop.stationary,
    connectivity = NMF::nmf.stop.connectivity,
    none = NULL
  )

  need <- c("doc_id", "comparison_id", "tf", "gene_key", "w")
  miss <- setdiff(need, names(doc_term_dt))
  if (length(miss)) cli::cli_abort("doc_term_dt missing columns: {paste(miss, collapse=', ')}")

  if (any(!is.finite(doc_term_dt$w)) || any(doc_term_dt$w < 0)) {
    cli::cli_abort("doc_term_dt$w must be finite and non-negative for NMF.")
  }

  docs  <- sort(unique(doc_term_dt$doc_id))
  terms <- sort(unique(doc_term_dt$gene_key))

  i <- match(doc_term_dt$doc_id, docs)
  j <- match(doc_term_dt$gene_key, terms)
  x <- as.numeric(doc_term_dt$w)

  V <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(docs), length(terms)),
    dimnames = list(docs, terms),
    giveCsparse = TRUE
  )

  nmf_args <- list(
    x = V,
    rank = as.integer(K),
    method = method,
    nrun = as.integer(nrun),
    seed = seed,
    .options = nmf_options,
    .pbackend = nmf_pbackend
  )
  if (!is.null(stop_fun)) nmf_args$stop <- stop_fun

  fit <- try(do.call(NMF::nmf, nmf_args), silent = TRUE)

  if (inherits(fit, "try-error")) {
    if (!isTRUE(allow_dense_fallback)) {
      cli::cli_abort("NMF::nmf() failed for sparse input; set allow_dense_fallback=TRUE (or install an NMF backend that supports sparse matrices).")
    }

    dense_gb <- (nrow(V) * ncol(V) * 8) / 1024^3
    if (!is.finite(dense_gb)) dense_gb <- Inf
    if (dense_gb > dense_max_gb) {
      cli::cli_abort(
        paste0(
          "NMF::nmf() failed for sparse input, and dense fallback would be ~",
          signif(dense_gb, 3),
          " GB (> dense_max_gb=",
          dense_max_gb,
          "). Reduce terms (e.g. lower top_terms_per_doc / raise min_df) or use a sparse-capable NMF implementation."
        )
      )
    }

    V_dense <- as.matrix(V)
    nmf_args_dense <- nmf_args
    nmf_args_dense$x <- V_dense
    fit <- do.call(NMF::nmf, nmf_args_dense)
  }

  W <- NMF::basis(fit) # docs x K
  doc_id <- rownames(W)

  rs <- rowSums(W)
  rs[!is.finite(rs) | rs <= 0] <- NA_real_
  gamma <- W / rs
  gamma[is.na(gamma)] <- 0

  n_doc <- nrow(gamma)
  n_topic <- ncol(gamma)

  gamma_long <- tibble::tibble(
    doc_id = rep(doc_id, times = n_topic),
    topic  = rep(seq_len(n_topic), each = n_doc),
    gamma  = as.vector(gamma)
  )

  meta <- unique(doc_term_dt[, c("doc_id", "comparison_id", "tf")])
  meta <- tibble::as_tibble(meta)
  gamma_long <- dplyr::left_join(gamma_long, meta, by = "doc_id")

  if (!is.null(doc_topic_top_n)) {
    doc_topic_top_n <- as.integer(doc_topic_top_n)
    if (!is.finite(doc_topic_top_n) || doc_topic_top_n <= 0) {
      cli::cli_abort("doc_topic_top_n must be a positive integer (or NULL).")
    }

    keep <- gamma_long |>
      dplyr::group_by(doc_id) |>
      dplyr::slice_max(order_by = gamma, n = doc_topic_top_n, with_ties = FALSE) |>
      dplyr::ungroup()
  } else {
    keep <- gamma_long[gamma_long$gamma >= gamma_cutoff, , drop = FALSE]
  }

  if (nrow(keep)) {
    keep <- keep[order(keep$doc_id, -keep$gamma, keep$topic), , drop = FALSE]
    split_topics <- split(keep$topic, keep$doc_id)
    topics_str <- vapply(split_topics, function(v) paste(v, collapse = multi_sep), character(1))
    assign_tbl <- tibble::tibble(doc_id = names(topics_str), topics_str = unname(topics_str))
  } else {
    assign_tbl <- tibble::tibble(doc_id = character(0), topics_str = character(0))
  }

  doc_topics <- dplyr::left_join(meta, assign_tbl, by = "doc_id")

  list(
    fit = fit,
    doc_topic_long = gamma_long,
    doc_topics = doc_topics
  )
}

topic_marker_overlap <- function(edges_tbl,
                                 doc_topics,
                                 markers_epi,
                                 markers_mes,
                                 abs_delta_min = NULL,
                                 gamma_cutoff_used = 0.2,
                                 topics_sep = ";",
                                 universe_genes = NULL,
                                 top_n_overlap_genes = 30L,
                                 verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }

  # ---- checks ----
  need_edges <- c("comparison_id", "tf", "gene_key")
  miss_edges <- setdiff(need_edges, names(edges_tbl))
  if (length(miss_edges)) {
    cli::cli_abort("edges_tbl missing required columns: {paste(miss_edges, collapse=', ')}")
  }

  need_docs <- c("comparison_id", "tf", "topics_str")
  miss_docs <- setdiff(need_docs, names(doc_topics))
  if (length(miss_docs)) {
    cli::cli_abort("doc_topics missing required columns: {paste(miss_docs, collapse=', ')}")
  }

  # ---- normalize markers ----
  markers_epi <- unique(as.character(markers_epi))
  markers_mes <- unique(as.character(markers_mes))

  # ---- doc -> topic mapping (explode topics_str safely) ----
  dt_doc <- data.table::as.data.table(doc_topics)[, .(comparison_id, tf, topics_str)]
  dt_doc <- dt_doc[!is.na(topics_str) & topics_str != ""]
  if (!nrow(dt_doc)) {
    cli::cli_abort("No doc_topics rows have topics_str (gamma_cutoff={gamma_cutoff_used} produced all NA?).")
  }

  dt_doc[, comparison_id := as.character(comparison_id)]
  dt_doc[, tf := as.character(tf)]
  dt_doc[, topics_str := as.character(topics_str)]

  # EXPAND rows: one row per (comparison_id, tf, topic)
  dt_doc_long <- dt_doc[
    ,
    .(topic = unlist(strsplit(topics_str, split = topics_sep, fixed = TRUE))),
    by = .(comparison_id, tf)
  ]
  dt_doc_long[, topic := as.integer(gsub("[^0-9]", "", as.character(topic)))]
  dt_doc_long <- dt_doc_long[is.finite(topic)]

  # ---- edges (optionally filter by abs(delta)) ----
  has_delta <- "delta_link_score" %in% names(edges_tbl)
  cols_e <- if (has_delta) c("comparison_id", "tf", "gene_key", "delta_link_score") else c("comparison_id", "tf", "gene_key")
  dt_e <- data.table::as.data.table(edges_tbl)[, ..cols_e]

  dt_e[, comparison_id := as.character(comparison_id)]
  dt_e[, tf := as.character(tf)]
  dt_e[, gene_key := as.character(gene_key)]

  if (!is.null(abs_delta_min)) {
    if (!has_delta) cli::cli_abort("abs_delta_min was provided but edges_tbl has no delta_link_score column.")
    dt_e <- dt_e[is.finite(delta_link_score) & abs(delta_link_score) >= abs_delta_min]
  }

  # ---- join doc-topics to edges to get genes per topic ----
  data.table::setkey(dt_doc_long, comparison_id, tf)
  data.table::setkey(dt_e, comparison_id, tf)

  # replicate edges across topics if multi-topic
  dt_join <- dt_doc_long[dt_e, nomatch = 0L, allow.cartesian = TRUE]
  if (!nrow(dt_join)) {
    cli::cli_abort("After join, no rows. Check that comparison_id/tf match between doc_topics and edges_tbl.")
  }

  # unique genes per topic (union across all docs assigned to that topic)
  dt_tgene <- unique(dt_join[, .(topic, gene_key)])

  # topic-level sizes
  topic_sizes <- dt_tgene[, .(n_genes = data.table::uniqueN(gene_key)), by = topic]
  topic_docs  <- dt_doc_long[, .(n_docs = .N, n_tfs = data.table::uniqueN(tf)), by = topic]

  # universe genes
  if (is.null(universe_genes)) {
    universe_genes <- unique(dt_e$gene_key)
  } else {
    universe_genes <- unique(as.character(universe_genes))
  }
  N <- length(universe_genes)

  calc_stats <- function(topic_gene_vec, marker_vec) {
    topic_gene_vec <- unique(topic_gene_vec)
    marker_vec <- unique(marker_vec)

    topic_u  <- intersect(topic_gene_vec, universe_genes)
    marker_u <- intersect(marker_vec, universe_genes)

    a <- length(intersect(topic_u, marker_u))
    b <- length(setdiff(topic_u, marker_u))
    c <- length(setdiff(marker_u, topic_u))
    d <- N - a - b - c
    d <- max(d, 0)

    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

    list(
      overlap_n = a,
      marker_coverage = if (length(marker_u) > 0) a / length(marker_u) else NA_real_,
      fisher_p = unname(ft$p.value),
      fisher_or = unname(ft$estimate)
    )
  }

  topics <- sort(unique(dt_tgene$topic))
  res_list <- lapply(topics, function(k) {
    genes_k <- dt_tgene[topic == k, gene_key]

    epi <- calc_stats(genes_k, markers_epi)
    mes <- calc_stats(genes_k, markers_mes)

    overlap_epi_genes <- intersect(unique(genes_k), markers_epi)
    overlap_mes_genes <- intersect(unique(genes_k), markers_mes)

    tibble::tibble(
      topic = k,

      overlap_epi_n = epi$overlap_n,
      epi_marker_coverage = epi$marker_coverage,
      epi_fisher_p = epi$fisher_p,
      epi_fisher_or = as.numeric(epi$fisher_or),
      overlap_epi_genes = paste(utils::head(sort(overlap_epi_genes), top_n_overlap_genes), collapse = ";"),

      overlap_mes_n = mes$overlap_n,
      mes_marker_coverage = mes$marker_coverage,
      mes_fisher_p = mes$fisher_p,
      mes_fisher_or = as.numeric(mes$fisher_or),
      overlap_mes_genes = paste(utils::head(sort(overlap_mes_genes), top_n_overlap_genes), collapse = ";")
    )
  })

  stats_tbl <- dplyr::left_join(dplyr::bind_rows(res_list), tibble::as_tibble(topic_sizes), by = "topic")
  stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_docs), by = "topic")

  best_epi <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_epi_n), epi_fisher_p, dplyr::desc(epi_marker_coverage)) |>
    dplyr::slice(1)

  best_mes <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_mes_n), mes_fisher_p, dplyr::desc(mes_marker_coverage)) |>
    dplyr::slice(1)

  if (isTRUE(verbose)) {
    message("[topic_marker_overlap] Universe genes N=", N,
            "; topics=", length(topics),
            "; abs_delta_min=", if (is.null(abs_delta_min)) "NULL" else abs_delta_min)
    message("[topic_marker_overlap] Best EPI topic: ", best_epi$topic,
            " (overlap=", best_epi$overlap_epi_n,
            ", coverage=", signif(best_epi$epi_marker_coverage, 3),
            ", fisher_p=", signif(best_epi$epi_fisher_p, 3), ")")
    message("[topic_marker_overlap] Best MES topic: ", best_mes$topic,
            " (overlap=", best_mes$overlap_mes_n,
            ", coverage=", signif(best_mes$mes_marker_coverage, 3),
            ", fisher_p=", signif(best_mes$mes_fisher_p, 3), ")")
  }

  list(
    stats = stats_tbl,
    best_epi = best_epi,
    best_mes = best_mes,
    topic_gene = tibble::as_tibble(dt_tgene)
  )
}

assign_edge_topics_lda <- function(edges_tbl,
                                   lda_fit,
                                   top_n = 1L,
                                   score_min = NULL,
                                   aggregate_peak = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("topicmodels", quietly = TRUE)) {
    cli::cli_abort("Need packages: data.table, topicmodels")
  }
  if (is.null(lda_fit$fit)) cli::cli_abort("lda_fit$fit is missing.")

  dt_edges <- data.table::as.data.table(edges_tbl)
  need <- c("comparison_id", "tf", "gene_key")
  miss <- setdiff(need, names(dt_edges))
  if (length(miss)) cli::cli_abort("edges_tbl missing: {paste(miss, collapse=', ')}")

  dt_edges[, comparison_id := as.character(comparison_id)]
  dt_edges[, tf := as.character(tf)]
  dt_edges[, gene_key := as.character(gene_key)]
  if (!("doc_id" %in% names(dt_edges))) {
    dt_edges[, doc_id := paste(comparison_id, tf, sep = "::")]
  }
  if ("peak_id" %in% names(dt_edges)) {
    dt_edges[, peak_id := as.character(peak_id)]
  } else {
    dt_edges[, peak_id := NA_character_]
  }
  dt_edges <- unique(dt_edges[, .(comparison_id, tf, gene_key, peak_id, doc_id)])
  dt_edges[, edge_id := paste(comparison_id, tf, gene_key, peak_id, sep = "|")]

  dt_gamma <- data.table::as.data.table(lda_fit$doc_topic_long)
  if (!all(c("doc_id", "topic", "gamma") %in% names(dt_gamma))) {
    cli::cli_abort("lda_fit$doc_topic_long must contain doc_id, topic, gamma.")
  }
  dt_gamma <- dt_gamma[is.finite(gamma) & gamma > 0]
  dt_gamma[, topic := as.integer(topic)]

  beta <- topicmodels::posterior(lda_fit$fit)$terms
  dt_beta <- data.table::as.data.table(as.table(beta))
  data.table::setnames(dt_beta, c("topic", "gene_key", "beta"))
  dt_beta[, topic := as.integer(as.character(topic))]
  dt_beta[, gene_key := as.character(gene_key)]
  dt_beta <- dt_beta[is.finite(beta) & beta > 0]

  dt_join <- dt_edges[dt_gamma, on = "doc_id", allow.cartesian = TRUE, nomatch = 0L]
  dt_join <- dt_join[dt_beta, on = .(gene_key, topic), nomatch = 0L]
  dt_join[, topic_score := gamma * beta]
  if (!is.null(score_min)) {
    dt_join <- dt_join[is.finite(topic_score) & topic_score >= score_min]
  }

  if (isTRUE(aggregate_peak)) {
    dt_join <- dt_join[, .(topic_score = sum(topic_score, na.rm = TRUE)),
                       by = .(comparison_id, tf, gene_key, topic)]
    dt_join[, peak_id := NA_character_]
    dt_join[, edge_id := paste(comparison_id, tf, gene_key, sep = "|")]
  } else {
    dt_join[, edge_id := paste(comparison_id, tf, gene_key, peak_id, sep = "|")]
  }

  data.table::setorder(dt_join, edge_id, -topic_score)
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
    dt_join <- dt_join[, head(.SD, top_n), by = edge_id]
  }

  tibble::as_tibble(dt_join)
}

assign_edge_topics_nmf <- function(edges_tbl,
                                   nmf_fit,
                                   top_n = 1L,
                                   score_min = NULL,
                                   aggregate_peak = FALSE) {
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("NMF", quietly = TRUE)) {
    cli::cli_abort("Need packages: data.table, NMF")
  }
  if (is.null(nmf_fit$fit)) cli::cli_abort("nmf_fit$fit is missing.")

  dt_edges <- data.table::as.data.table(edges_tbl)
  need <- c("comparison_id", "tf", "gene_key")
  miss <- setdiff(need, names(dt_edges))
  if (length(miss)) cli::cli_abort("edges_tbl missing: {paste(miss, collapse=', ')}")

  dt_edges[, comparison_id := as.character(comparison_id)]
  dt_edges[, tf := as.character(tf)]
  dt_edges[, gene_key := as.character(gene_key)]
  if (!("doc_id" %in% names(dt_edges))) {
    dt_edges[, doc_id := paste(comparison_id, tf, sep = "::")]
  }
  if ("peak_id" %in% names(dt_edges)) {
    dt_edges[, peak_id := as.character(peak_id)]
  } else {
    dt_edges[, peak_id := NA_character_]
  }
  dt_edges <- unique(dt_edges[, .(comparison_id, tf, gene_key, peak_id, doc_id)])
  dt_edges[, edge_id := paste(comparison_id, tf, gene_key, peak_id, sep = "|")]

  dt_gamma <- data.table::as.data.table(nmf_fit$doc_topic_long)
  if (!all(c("doc_id", "topic", "gamma") %in% names(dt_gamma))) {
    cli::cli_abort("nmf_fit$doc_topic_long must contain doc_id, topic, gamma.")
  }
  dt_gamma <- dt_gamma[is.finite(gamma) & gamma > 0]
  dt_gamma[, topic := as.integer(topic)]

  coef_mat <- NMF::coef(nmf_fit$fit)
  dt_beta <- data.table::as.data.table(as.table(coef_mat))
  data.table::setnames(dt_beta, c("topic", "gene_key", "beta"))
  dt_beta[, topic := as.integer(as.character(topic))]
  dt_beta[, gene_key := as.character(gene_key)]
  dt_beta <- dt_beta[is.finite(beta) & beta > 0]

  dt_join <- dt_edges[dt_gamma, on = "doc_id", allow.cartesian = TRUE, nomatch = 0L]
  dt_join <- dt_join[dt_beta, on = .(gene_key, topic), nomatch = 0L]
  dt_join[, topic_score := gamma * beta]
  if (!is.null(score_min)) {
    dt_join <- dt_join[is.finite(topic_score) & topic_score >= score_min]
  }

  if (isTRUE(aggregate_peak)) {
    dt_join <- dt_join[, .(topic_score = sum(topic_score, na.rm = TRUE)),
                       by = .(comparison_id, tf, gene_key, topic)]
    dt_join[, peak_id := NA_character_]
    dt_join[, edge_id := paste(comparison_id, tf, gene_key, sep = "|")]
  } else {
    dt_join[, edge_id := paste(comparison_id, tf, gene_key, peak_id, sep = "|")]
  }

  data.table::setorder(dt_join, edge_id, -topic_score)
  if (!is.null(top_n) && is.finite(top_n) && top_n > 0L) {
    dt_join <- dt_join[, head(.SD, top_n), by = edge_id]
  }

  tibble::as_tibble(dt_join)
}

assign_edge_topics_louvain <- function(edges_tbl,
                                       weight_col = "delta_link_score",
                                       min_abs_weight = NULL,
                                       seed = 1L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    cli::cli_abort("Please install igraph for Louvain baseline.")
  }
  need <- c("tf", "gene_key")
  miss <- setdiff(need, names(edges_tbl))
  if (length(miss)) cli::cli_abort("edges_tbl missing: {paste(miss, collapse=', ')}")
  if (!(weight_col %in% names(edges_tbl))) {
    cli::cli_abort("edges_tbl missing weight_col: {weight_col}")
  }

  dt <- data.table::as.data.table(edges_tbl)
  dt[, tf := as.character(tf)]
  dt[, gene_key := as.character(gene_key)]
  if ("comparison_id" %in% names(dt)) dt[, comparison_id := as.character(comparison_id)]
  if ("peak_id" %in% names(dt)) dt[, peak_id := as.character(peak_id)]

  dt[, weight := as.numeric(get(weight_col))]
  dt <- dt[is.finite(weight) & weight > 0]
  if (!is.null(min_abs_weight) && is.finite(min_abs_weight)) {
    dt <- dt[abs(weight) >= min_abs_weight]
  }
  if (!nrow(dt)) {
    return(tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      topic = integer(0),
      topic_score = numeric(0)
    ))
  }

  dt_graph <- dt[, .(weight = sum(abs(weight), na.rm = TRUE)), by = .(tf, gene_key)]
  g <- igraph::graph_from_data_frame(dt_graph, directed = FALSE)
  if (!igraph::ecount(g)) {
    return(tibble::tibble(
      comparison_id = character(0),
      tf = character(0),
      gene_key = character(0),
      peak_id = character(0),
      topic = integer(0),
      topic_score = numeric(0)
    ))
  }

  set.seed(seed)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  mem <- igraph::membership(cl)
  mem_dt <- data.table::data.table(node = names(mem), topic = as.integer(mem))

  dt_edges <- dt[, .(
    comparison_id = if ("comparison_id" %in% names(dt)) comparison_id else NA_character_,
    tf,
    gene_key,
    peak_id = if ("peak_id" %in% names(dt)) peak_id else NA_character_,
    topic_score = weight
  )]
  dt_edges <- unique(dt_edges)

  dt_edges[, topic := NA_integer_]
  mem_gene <- data.table::copy(mem_dt)
  data.table::setnames(mem_gene, "node", "gene_key")
  mem_tf <- data.table::copy(mem_dt)
  data.table::setnames(mem_tf, "node", "tf")
  dt_edges[mem_gene, topic := i.topic, on = "gene_key"]
  dt_edges[is.na(topic), mem_tf, topic := i.topic, on = "tf"]
  dt_edges <- dt_edges[!is.na(topic)]

  tibble::as_tibble(dt_edges)
}

topic_marker_overlap_from_edge_topics <- function(edge_topics,
                                                  markers_epi,
                                                  markers_mes,
                                                  universe_genes = NULL,
                                                  top_n_overlap_genes = 30L,
                                                  marker_mode = c("exclusive", "union"),
                                                  exclusive_method = c("mean_score", "sum_score", "edge_count"),
                                                  exclusive_top_n = 2L,
                                                  compute_gene_only = FALSE,
                                                  verbose = TRUE) {
  marker_mode <- match.arg(marker_mode)
  exclusive_method <- match.arg(exclusive_method)
  exclusive_top_n <- as.integer(exclusive_top_n)
  if (!is.finite(exclusive_top_n) || exclusive_top_n < 1L) {
    exclusive_top_n <- 1L
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  need <- c("topic", "gene_key")
  miss <- setdiff(need, names(edge_topics))
  if (length(miss)) cli::cli_abort("edge_topics missing: {paste(miss, collapse=', ')}")

  dt <- data.table::as.data.table(edge_topics)
  dt[, topic := as.integer(topic)]
  dt[, gene_key := as.character(gene_key)]
  if ("tf" %in% names(dt)) {
    dt[, tf := as.character(tf)]
  } else {
    dt[, tf := NA_character_]
  }
  if ("peak_id" %in% names(dt)) {
    dt[, peak_id := as.character(peak_id)]
  }
  if ("comparison_id" %in% names(dt)) {
    dt[, comparison_id := as.character(comparison_id)]
  }
  if (!("topic_score" %in% names(dt))) {
    dt[, topic_score := NA_real_]
  }
  dt <- dt[is.finite(topic) &
             ((is.character(gene_key) & !is.na(gene_key) & gene_key != "") |
                (is.character(tf) & !is.na(tf) & tf != ""))]

  topic_gene <- unique(dt[!is.na(gene_key) & gene_key != "", .(topic, gene_key)])
  topic_tf <- unique(dt[!is.na(tf) & tf != "", .(topic, tf)])
  topic_entities <- data.table::rbindlist(list(
    if (nrow(topic_gene)) data.table::data.table(topic = topic_gene$topic, entity = topic_gene$gene_key),
    if (nrow(topic_tf)) data.table::data.table(topic = topic_tf$topic, entity = topic_tf$tf)
  ), use.names = TRUE, fill = TRUE)
  topic_entities <- unique(topic_entities)

  topic_sizes <- topic_gene[, .(n_genes_unique = data.table::uniqueN(gene_key)), by = topic]
  topic_sizes_total <- dt[!is.na(gene_key) & gene_key != "", .(n_genes_total = .N), by = topic]
  topic_entities_sizes <- topic_entities[, .(n_entities_unique = data.table::uniqueN(entity)), by = topic]
  topic_edges <- if ("edge_id" %in% names(dt)) {
    dt[, .(n_edges = data.table::uniqueN(edge_id)), by = topic]
  } else {
    dt[, .(n_edges = .N), by = topic]
  }
  topic_edges_unique <- NULL
  if (all(c("tf", "gene_key") %in% names(dt))) {
    if ("peak_id" %in% names(dt)) {
      dt[, edge_key := paste(tf, gene_key, peak_id, sep = "|")]
    } else {
      dt[, edge_key := paste(tf, gene_key, sep = "|")]
    }
    topic_edges_unique <- dt[, .(n_edges_unique = data.table::uniqueN(edge_key)), by = topic]
  }
  topic_comparisons <- if ("comparison_id" %in% names(dt)) {
    dt[, .(n_comparisons = data.table::uniqueN(comparison_id)), by = topic]
  } else {
    data.table::data.table(topic = unique(dt$topic), n_comparisons = NA_integer_)
  }

  if (is.null(universe_genes)) {
    universe_genes_only <- unique(topic_gene$gene_key)
  } else {
    universe_genes_only <- unique(as.character(universe_genes))
  }
  universe_entities <- unique(c(universe_genes_only, topic_tf$tf))
  N <- length(universe_entities)

  markers_epi <- unique(as.character(markers_epi))
  markers_mes <- unique(as.character(markers_mes))
  markers_all <- unique(c(markers_epi, markers_mes))

  markers_epi_all_entities <- intersect(markers_epi, universe_entities)
  markers_mes_all_entities <- intersect(markers_mes, universe_entities)
  markers_epi_all_genes <- intersect(markers_epi, universe_genes_only)
  markers_mes_all_genes <- intersect(markers_mes, universe_genes_only)

  marker_assign <- NULL
  if (marker_mode == "exclusive") {
    dt_entities <- data.table::rbindlist(list(
      dt[!is.na(gene_key) & gene_key != "", .(topic, entity = gene_key, topic_score)],
      dt[!is.na(tf) & tf != "", .(topic, entity = tf, topic_score)]
    ), use.names = TRUE, fill = TRUE)
    dt_entities <- dt_entities[entity %in% markers_all]
    if (nrow(dt_entities)) {
      use_method <- exclusive_method
      if (use_method != "edge_count" && !any(is.finite(dt_entities$topic_score))) {
        use_method <- "edge_count"
      }
      if (use_method == "mean_score") {
        dt_entities <- dt_entities[, .(score = mean(topic_score, na.rm = TRUE)), by = .(topic, entity)]
      } else if (use_method == "sum_score") {
        dt_entities <- dt_entities[, .(score = sum(topic_score, na.rm = TRUE)), by = .(topic, entity)]
      } else {
        dt_entities <- dt_entities[, .(score = .N), by = .(topic, entity)]
      }
      data.table::setorder(dt_entities, entity, -score, topic)
      marker_assign <- dt_entities[is.finite(score), .SD[seq_len(min(.N, exclusive_top_n))], by = entity]
    } else {
      marker_assign <- data.table::data.table(entity = character(0), topic = integer(0), score = numeric(0))
    }
  }

  calc_stats <- function(topic_vec, marker_all, universe_vec, overlap_vec) {
    topic_u  <- intersect(unique(topic_vec), universe_vec)
    marker_u <- intersect(unique(marker_all), universe_vec)
    overlap_u <- intersect(unique(overlap_vec), topic_u)
    overlap_u <- intersect(overlap_u, marker_u)

    a <- length(overlap_u)
    b <- length(setdiff(topic_u, overlap_u))
    c <- length(setdiff(marker_u, overlap_u))
    d <- length(universe_vec) - a - b - c
    d <- max(d, 0)

    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

    list(
      overlap_n = a,
      marker_coverage = if (length(marker_u) > 0) a / length(marker_u) else NA_real_,
      fisher_p = unname(ft$p.value),
      fisher_or = unname(ft$estimate)
    )
  }

  topics <- sort(unique(topic_gene$topic))
  res_list <- lapply(topics, function(k) {
    genes_k <- topic_gene[topic == k, gene_key]
    entities_k <- topic_entities[topic == k, entity]
    if (marker_mode == "exclusive") {
      markers_epi_k <- marker_assign[topic == k & entity %in% markers_epi, entity]
      markers_mes_k <- marker_assign[topic == k & entity %in% markers_mes, entity]
      overlap_epi <- intersect(unique(entities_k), markers_epi_k)
      overlap_mes <- intersect(unique(entities_k), markers_mes_k)
      if (isTRUE(compute_gene_only)) {
        overlap_epi_genes_only <- intersect(unique(genes_k), markers_epi_k)
        overlap_mes_genes_only <- intersect(unique(genes_k), markers_mes_k)
      }
    } else {
      overlap_epi <- intersect(unique(entities_k), markers_epi_all_entities)
      overlap_mes <- intersect(unique(entities_k), markers_mes_all_entities)
      if (isTRUE(compute_gene_only)) {
        overlap_epi_genes_only <- intersect(unique(genes_k), markers_epi_all_genes)
        overlap_mes_genes_only <- intersect(unique(genes_k), markers_mes_all_genes)
      }
    }

    epi <- calc_stats(entities_k, markers_epi_all_entities, universe_entities, overlap_epi)
    mes <- calc_stats(entities_k, markers_mes_all_entities, universe_entities, overlap_mes)

    overlap_epi_genes <- overlap_epi
    overlap_mes_genes <- overlap_mes

    out_cols <- list(
      topic = k,
      overlap_epi_n = epi$overlap_n,
      epi_marker_coverage = epi$marker_coverage,
      epi_fisher_p = epi$fisher_p,
      epi_fisher_or = as.numeric(epi$fisher_or),
      overlap_epi_genes = paste(utils::head(sort(overlap_epi_genes), top_n_overlap_genes), collapse = ";"),
      overlap_mes_n = mes$overlap_n,
      mes_marker_coverage = mes$marker_coverage,
      mes_fisher_p = mes$fisher_p,
      mes_fisher_or = as.numeric(mes$fisher_or),
      overlap_mes_genes = paste(utils::head(sort(overlap_mes_genes), top_n_overlap_genes), collapse = ";")
    )
    if (isTRUE(compute_gene_only)) {
      epi_gene <- calc_stats(genes_k, markers_epi_all_genes, universe_genes_only, overlap_epi_genes_only)
      mes_gene <- calc_stats(genes_k, markers_mes_all_genes, universe_genes_only, overlap_mes_genes_only)
      out_cols <- c(out_cols, list(
        overlap_epi_n_genes_only = epi_gene$overlap_n,
        epi_marker_coverage_genes_only = epi_gene$marker_coverage,
        epi_fisher_p_genes_only = epi_gene$fisher_p,
        epi_fisher_or_genes_only = as.numeric(epi_gene$fisher_or),
        overlap_epi_genes_only = paste(utils::head(sort(overlap_epi_genes_only), top_n_overlap_genes), collapse = ";"),
        overlap_mes_n_genes_only = mes_gene$overlap_n,
        mes_marker_coverage_genes_only = mes_gene$marker_coverage,
        mes_fisher_p_genes_only = mes_gene$fisher_p,
        mes_fisher_or_genes_only = as.numeric(mes_gene$fisher_or),
        overlap_mes_genes_only = paste(utils::head(sort(overlap_mes_genes_only), top_n_overlap_genes), collapse = ";")
      ))
    }
    tibble::as_tibble(out_cols)
  })

  stats_tbl <- dplyr::left_join(dplyr::bind_rows(res_list), tibble::as_tibble(topic_sizes), by = "topic")
  stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_sizes_total), by = "topic")
  stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_entities_sizes), by = "topic")
  stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_edges), by = "topic")
  if (!is.null(topic_edges_unique)) {
    stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_edges_unique), by = "topic")
  } else {
    stats_tbl$n_edges_unique <- NA_integer_
  }
  stats_tbl <- dplyr::left_join(stats_tbl, tibble::as_tibble(topic_comparisons), by = "topic")
  stats_tbl$n_edges_avg <- ifelse(is.finite(stats_tbl$n_comparisons) & stats_tbl$n_comparisons > 0,
                                  stats_tbl$n_edges / stats_tbl$n_comparisons,
                                  NA_real_)
  stats_tbl <- stats_tbl |>
    dplyr::mutate(
      n_genes = n_genes_unique,
      marker_mode = marker_mode,
      exclusive_method = ifelse(marker_mode == "exclusive", exclusive_method, NA_character_),
      exclusive_top_n = ifelse(marker_mode == "exclusive", exclusive_top_n, NA_integer_)
    )

  best_epi <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_epi_n), epi_fisher_p, dplyr::desc(epi_marker_coverage)) |>
    dplyr::slice(1)
  best_mes <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_mes_n), mes_fisher_p, dplyr::desc(mes_marker_coverage)) |>
    dplyr::slice(1)

  if (isTRUE(verbose)) {
    message("[topic_marker_overlap_from_edge_topics] Universe entities N=", N,
            "; topics=", length(topics))
    message("[topic_marker_overlap_from_edge_topics] Best EPI topic: ", best_epi$topic,
            " (overlap=", best_epi$overlap_epi_n,
            ", coverage=", signif(best_epi$epi_marker_coverage, 3),
            ", fisher_p=", signif(best_epi$epi_fisher_p, 3), ")")
    message("[topic_marker_overlap_from_edge_topics] Best MES topic: ", best_mes$topic,
            " (overlap=", best_mes$overlap_mes_n,
            ", coverage=", signif(best_mes$mes_marker_coverage, 3),
            ", fisher_p=", signif(best_mes$mes_fisher_p, 3), ")")
  }

  list(
    stats = stats_tbl,
    best_epi = best_epi,
    best_mes = best_mes,
    topic_gene = tibble::as_tibble(topic_gene)
  )
}

topic_marker_overlap_from_edge_topics_by_group <- function(edge_topics,
                                                           group_col = "comparison_id",
                                                           markers_epi,
                                                           markers_mes,
                                                           top_n_overlap_genes = 30L,
                                                           marker_mode = c("exclusive", "union"),
                                                           exclusive_method = c("mean_score", "sum_score", "edge_count"),
                                                           exclusive_top_n = 2L,
                                                           compute_gene_only = FALSE,
                                                           verbose = TRUE) {
  marker_mode <- match.arg(marker_mode)
  exclusive_method <- match.arg(exclusive_method)
  exclusive_top_n <- as.integer(exclusive_top_n)
  if (!is.finite(exclusive_top_n) || exclusive_top_n < 1L) {
    exclusive_top_n <- 1L
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  dt <- data.table::as.data.table(edge_topics)
  if (!(group_col %in% names(dt))) {
    cli::cli_abort("edge_topics missing group_col: {group_col}")
  }

  groups <- unique(dt[[group_col]])
  res_stats <- vector("list", length(groups))
  res_best <- vector("list", length(groups))

  for (i in seq_along(groups)) {
    g <- groups[[i]]
    dt_g <- dt[dt[[group_col]] == g, , drop = FALSE]
    if (!nrow(dt_g)) next
    out <- topic_marker_overlap_from_edge_topics(
      edge_topics = tibble::as_tibble(dt_g),
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      universe_genes = NULL,
      top_n_overlap_genes = top_n_overlap_genes,
      marker_mode = marker_mode,
      exclusive_method = exclusive_method,
      exclusive_top_n = exclusive_top_n,
      compute_gene_only = compute_gene_only,
      verbose = FALSE
    )
    stats_g <- out$stats
    stats_g[[group_col]] <- g
    res_stats[[i]] <- stats_g

    best_epi <- out$best_epi
    best_mes <- out$best_mes
    res_best[[i]] <- tibble::tibble(
      !!group_col := g,
      best_epi_topic = if (nrow(best_epi)) best_epi$topic[[1]] else NA_integer_,
      best_epi_overlap = if (nrow(best_epi)) best_epi$overlap_epi_n[[1]] else NA_real_,
      best_epi_purity = if (nrow(best_epi) && is.finite(best_epi$n_entities_unique[[1]]) && best_epi$n_entities_unique[[1]] > 0) {
        best_epi$overlap_epi_n[[1]] / best_epi$n_entities_unique[[1]]
      } else {
        NA_real_
      },
      best_mes_topic = if (nrow(best_mes)) best_mes$topic[[1]] else NA_integer_,
      best_mes_overlap = if (nrow(best_mes)) best_mes$overlap_mes_n[[1]] else NA_real_,
      best_mes_purity = if (nrow(best_mes) && is.finite(best_mes$n_entities_unique[[1]]) && best_mes$n_entities_unique[[1]] > 0) {
        best_mes$overlap_mes_n[[1]] / best_mes$n_entities_unique[[1]]
      } else {
        NA_real_
      }
    )
  }

  stats_tbl <- dplyr::bind_rows(res_stats)
  best_tbl <- dplyr::bind_rows(res_best)

  if (isTRUE(verbose)) {
    message("[topic_marker_overlap_from_edge_topics_by_group] groups=", length(groups),
            "; group_col=", group_col)
  }

  list(stats = stats_tbl, best = best_tbl)
}

merge_topics_by_overlap <- function(topic_gene_tbl,
                                    overlap_thresh = 0.7,
                                    topic_col = "topic",
                                    gene_col = "gene_key") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  dt <- data.table::as.data.table(topic_gene_tbl)
  if (!(topic_col %in% names(dt)) || !(gene_col %in% names(dt))) {
    cli::cli_abort("topic_gene_tbl must contain {topic_col} and {gene_col}.")
  }
  dt <- dt[!is.na(get(topic_col)) & !is.na(get(gene_col))]
  dt[, topic_val := as.integer(get(topic_col))]
  dt[, gene_val := as.character(get(gene_col))]

  topics <- sort(unique(dt$topic_val))
  gene_sets <- lapply(topics, function(t) unique(dt[topic_val == t, gene_val]))
  names(gene_sets) <- as.character(topics)

  assigned <- setNames(rep(FALSE, length(topics)), as.character(topics))
  mapping <- data.table::data.table(topic = integer(0), topic_merged = integer(0))
  merged_gene_tbl <- list()

  for (i in seq_along(topics)) {
    t <- topics[[i]]
    t_key <- as.character(t)
    if (isTRUE(assigned[[t_key]])) next

    current_genes <- gene_sets[[t_key]]
    members <- c(t)
    assigned[[t_key]] <- TRUE

    for (j in seq_along(topics)) {
      t2 <- topics[[j]]
      t2_key <- as.character(t2)
      if (isTRUE(assigned[[t2_key]])) next
      genes2 <- gene_sets[[t2_key]]
      if (!length(genes2)) next
      overlap <- length(intersect(current_genes, genes2))
      denom1 <- length(current_genes)
      denom2 <- length(genes2)
      ratio1 <- if (denom1 > 0L) overlap / denom1 else 0
      ratio2 <- if (denom2 > 0L) overlap / denom2 else 0
      if (max(ratio1, ratio2) >= overlap_thresh) {
        members <- c(members, t2)
        assigned[[t2_key]] <- TRUE
        current_genes <- union(current_genes, genes2)
      }
    }

    mapping <- data.table::rbindlist(list(
      mapping,
      data.table::data.table(topic = members, topic_merged = t)
    ), use.names = TRUE, fill = TRUE)
    merged_gene_tbl[[t_key]] <- data.table::data.table(
      topic = t,
      gene_key = unique(current_genes)
    )
  }

  list(
    mapping = mapping,
    topic_genes = data.table::rbindlist(merged_gene_tbl, use.names = TRUE, fill = TRUE)
  )
}

merge_edge_topics_by_overlap <- function(edge_topics,
                                         overlap_thresh = 0.7,
                                         topic_col = "topic",
                                         gene_col = "gene_key") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  dt <- data.table::as.data.table(edge_topics)
  if (!(topic_col %in% names(dt)) || !(gene_col %in% names(dt))) {
    cli::cli_abort("edge_topics must contain {topic_col} and {gene_col}.")
  }
  dt[, topic_raw := get(topic_col)]
  dt[, topic_int := suppressWarnings(as.integer(get(topic_col)))]
  use_int <- !all(is.na(dt$topic_int))
  if (use_int) {
    dt[, (topic_col) := topic_int]
  } else {
    dt[, (topic_col) := as.character(topic_raw)]
  }
  topic_gene_tbl <- dt[, .(topic = get(topic_col), gene_key = get(gene_col))]
  topic_gene_tbl <- topic_gene_tbl[!is.na(topic) & !is.na(gene_key)]
  merge_res <- merge_topics_by_overlap(
    topic_gene_tbl = topic_gene_tbl,
    overlap_thresh = overlap_thresh,
    topic_col = "topic",
    gene_col = "gene_key"
  )
  map <- data.table::as.data.table(merge_res$mapping)
  if (use_int) {
    map[, topic := suppressWarnings(as.integer(topic))]
    map[, topic_merged := suppressWarnings(as.integer(topic_merged))]
  } else {
    map[, topic := as.character(topic)]
    map[, topic_merged := as.character(topic_merged)]
  }
  on_cols <- setNames("topic", topic_col)
  dt[map, topic_merged := i.topic_merged, on = on_cols]
  dt[is.na(topic_merged), topic_merged := topic_raw]
  dt[, topic := topic_merged]

  list(
    edge_topics = tibble::as_tibble(dt),
    mapping = tibble::as_tibble(map),
    topic_genes = tibble::as_tibble(merge_res$topic_genes)
  )
}

summarize_topic_benchmark_edge_topics <- function(topic_benchmark_dir,
                                                  fit_dir = file.path(topic_benchmark_dir, "fits"),
                                                  markers_epi,
                                                  markers_mes,
                                                  out_file = file.path(topic_benchmark_dir, "topic_benchmark_all_topics.csv"),
                                                  top_n_overlap_genes = 50L,
                                                  marker_mode = c("exclusive", "union"),
                                                  exclusive_method = c("mean_score", "sum_score", "edge_count"),
                                                  exclusive_top_n = 2L,
                                                  compute_gene_only = FALSE,
                                                  drop_gene_only = TRUE) {
  marker_mode <- match.arg(marker_mode)
  exclusive_method <- match.arg(exclusive_method)
  exclusive_top_n <- as.integer(exclusive_top_n)
  if (!is.finite(exclusive_top_n) || exclusive_top_n < 1L) {
    exclusive_top_n <- 1L
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table for this function (scale).")
  }
  edge_files <- list.files(fit_dir, pattern = "_edge_topics\\.csv$", full.names = TRUE)
  if (!length(edge_files)) {
    message("[summarize_topic_benchmark_edge_topics] No edge_topics files found in: ", fit_dir)
    return(tibble::tibble())
  }

  parse_meta <- function(fname) {
    base <- sub("_edge_topics\\.csv$", "", basename(fname))
    if (grepl("^lda_", base)) {
      m <- regexec("^lda_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", base)
      hit <- regmatches(base, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100, run_id = base))
      }
    }
    if (grepl("^nmf_", base)) {
      m <- regexec("^nmf_(.+)_K(\\d+)_([^_]+)_g(\\d+)$", base)
      hit <- regmatches(base, m)[[1]]
      if (length(hit)) {
        return(list(label = hit[2], method = hit[4], K = as.integer(hit[3]),
                    gamma_cutoff = as.numeric(hit[5]) / 100, run_id = base))
      }
    }
    if (grepl("^louvain_", base)) {
      label <- sub("^louvain_", "", base)
      return(list(label = label, method = "louvain", K = NA_integer_,
                  gamma_cutoff = NA_real_, run_id = base))
    }
    list(label = base, method = NA_character_, K = NA_integer_, gamma_cutoff = NA_real_, run_id = base)
  }

  stats_list <- lapply(edge_files, function(f) {
    meta <- parse_meta(f)
    edge_topics <- readr::read_csv(f, show_col_types = FALSE)
    if (!nrow(edge_topics)) return(NULL)

    out <- topic_marker_overlap_from_edge_topics(
      edge_topics = edge_topics,
      markers_epi = markers_epi,
      markers_mes = markers_mes,
      universe_genes = unique(edge_topics$gene_key),
      top_n_overlap_genes = top_n_overlap_genes,
      marker_mode = marker_mode,
      exclusive_method = exclusive_method,
      exclusive_top_n = exclusive_top_n,
      compute_gene_only = compute_gene_only,
      verbose = FALSE
    )

    stats_tbl <- out$stats |>
      dplyr::mutate(
        label = meta$label,
        method = meta$method,
        K = meta$K,
        gamma_cutoff = meta$gamma_cutoff,
        run_id = meta$run_id,
        purity_epi = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                            overlap_epi_n / n_entities_unique, NA_real_),
        purity_mes = ifelse(is.finite(n_entities_unique) & n_entities_unique > 0,
                            overlap_mes_n / n_entities_unique, NA_real_),
        purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE)
      )
    if (all(c("overlap_epi_n_genes_only", "overlap_mes_n_genes_only", "n_genes_unique") %in% names(stats_tbl))) {
      stats_tbl <- stats_tbl |>
        dplyr::mutate(
          purity_epi_genes_only = ifelse(is.finite(n_genes_unique) & n_genes_unique > 0,
                                         overlap_epi_n_genes_only / n_genes_unique, NA_real_),
          purity_mes_genes_only = ifelse(is.finite(n_genes_unique) & n_genes_unique > 0,
                                         overlap_mes_n_genes_only / n_genes_unique, NA_real_)
        )
    }
    stats_tbl
  })

  stats_tbl <- dplyr::bind_rows(Filter(Negate(is.null), stats_list))
  if (isTRUE(drop_gene_only) && nrow(stats_tbl)) {
    drop_cols <- c(
      "overlap_epi_n_genes_only",
      "epi_marker_coverage_genes_only",
      "epi_fisher_p_genes_only",
      "epi_fisher_or_genes_only",
      "overlap_epi_genes_only",
      "overlap_mes_n_genes_only",
      "mes_marker_coverage_genes_only",
      "mes_fisher_p_genes_only",
      "mes_fisher_or_genes_only",
      "overlap_mes_genes_only",
      "purity_epi_genes_only",
      "purity_mes_genes_only"
    )
    stats_tbl <- stats_tbl |> dplyr::select(-dplyr::any_of(drop_cols))
  }
  if (!is.null(out_file) && nrow(stats_tbl)) {
    readr::write_csv(stats_tbl, out_file)
  }
  stats_tbl
}

cluster_comparisons_by_tf_delta <- function(edges_all_tidy, min_abs_delta = 1) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Please install data.table.")
  }
  dt <- data.table::as.data.table(edges_all_tidy)
  dt <- dt[is.finite(delta_link_score) & abs(delta_link_score) >= min_abs_delta]
  dt[, w := abs(delta_link_score)]

  agg <- dt[, .(w = sum(w)), by = .(comparison_id, tf)]
  wide <- data.table::dcast(agg, comparison_id ~ tf, value.var = "w", fill = 0)

  rn <- wide$comparison_id
  X  <- as.matrix(wide[, -1, drop = FALSE])
  rownames(X) <- rn

  hc <- stats::hclust(stats::dist(X), method = "ward.D2")
  hc
}
