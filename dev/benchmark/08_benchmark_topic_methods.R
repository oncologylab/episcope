# Topic benchmarking (LDA variants) -------------------------------------------
# Standalone script to benchmark topic models on nutrient_stress edges.
library(episcope)

plot_dir <- file.path("/data/homes/yl814/episcope_test/nutrient_stress", "plot_lineage_plasticity_related_subnetworks")
markers_epi <- readr::read_tsv(file.path(plot_dir, "markers", "epithelial.markers_nutrient_stress_all_lines.txt")) |> dplyr::pull(HGNC)
markers_mes <- readr::read_tsv(file.path(plot_dir, "markers", "mesenchymal.markers_nutrient_stress_all_lines.txt")) |> dplyr::pull(HGNC)

source(file.path("R", "utils_grn_lda_nmf.R"))

base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"
lighting_folder <- file.path(
  base_dir,
  "lighting_fp_tf_corr_FDR_0.05_genehancer_jaspar2024_regulated_genes_1.5_delta_link_1_spearman"
)
rerun <- FALSE   # set TRUE to overwrite existing summaries/CSVs
run_nmf_demo <- FALSE
run_nmf_subset <- FALSE
use_parallel <- TRUE
run_lda_benchmark <- TRUE
run_nmf_benchmark <- FALSE
run_fast_methods <- TRUE
bench_cores <- as.integer(Sys.getenv("TOPIC_BENCH_CORES", unset = parallel::detectCores() - 1))
bench_cores <- max(1, bench_cores)

if (!exists("markers_epi") || !exists("markers_mes")) {
  stop("Please load markers_epi and markers_mes before running this script.")
}

# load and standardize edges
delta_csvs <- list.files(lighting_folder, "_delta_links\\.csv$", full.names = TRUE)
edges_all_tidy <- load_delta_links_all_tidy(delta_csvs)
filtered_csvs <- list.files(lighting_folder, "_delta_links_filtered\\.csv$", full.names = TRUE)
filtered_only <- filtered_csvs[grepl("_delta_links_filtered\\.csv$", filtered_csvs)]
union_keys <- read_union_edge_keys_from_filtered(filtered_only)
edges_filtered_tidy <- filter_edges_all_tidy_by_union(edges_all_tidy, union_keys)

# positive FC features (for alternative weighting)
edges_fc <- edges_filtered_tidy |>
  dplyr::mutate(
    fc_link_score = pmax(1e-6, link_score_case / pmax(link_score_ctrl, 1e-6)),
    fc_fp_bed     = pmax(1e-6, fp_bed_score_case / pmax(fp_bed_score_ctrl, 1e-6)),
    fc_gene_expr  = pmax(1e-6, 2^log2FC_gene_expr)
  )

# Combined (all comparisons), delta + optional expression weights
dt_docs_delta <- build_pooled_doc_term(
  edges_filtered_tidy,
  which = "abs",
  min_abs_delta = 1,
  top_terms_per_doc = 1000,
  min_df = 2,
  fp_col = "delta_fp_bed_score", fp_transform = "abs", fp_weight = 1,
  tf_expr_col = "log2FC_tf_expr", tf_expr_transform = "abs", tf_expr_weight = 1,
  gene_expr_col = "log2FC_gene_expr", gene_expr_transform = "pos", gene_expr_weight = 2,
  scale_quantile = 0.9
)
message("[delta docs] rows=", nrow(dt_docs_delta), " docs=", dplyr::n_distinct(dt_docs_delta$doc_id), " terms=", dplyr::n_distinct(dt_docs_delta$gene_key))

# Combined (all comparisons), FC-based weights
dt_docs_fc <- build_pooled_doc_term(
  edges_fc,
  which = "abs",
  min_abs_delta = 1,
  top_terms_per_doc = 1000,
  min_df = 2,
  fp_col = "fc_fp_bed",
  fp_transform = "abs",
  fp_weight = 1,
  gene_expr_col = "fc_gene_expr",
  gene_expr_transform = "abs",
  gene_expr_weight = 1,
  scale_quantile = 0.9
)
message("[fc docs] rows=", nrow(dt_docs_fc), " docs=", dplyr::n_distinct(dt_docs_fc$doc_id), " terms=", dplyr::n_distinct(dt_docs_fc$gene_key))

# Per-comparison builds (delta)
dt_docs_by_cmp <- lapply(split(edges_filtered_tidy, edges_filtered_tidy$comparison_id), function(df) {
  build_pooled_doc_term(
    df,
    which = "abs",
    min_abs_delta = 1,
    top_terms_per_doc = 1000,
    min_df = 2
  )
})
message("[per-comparison] built ", length(dt_docs_by_cmp), " doc-term tables")

# Per-condition (non-delta) edge scores and gene expression (wide -> long)
cond_cols_link <- grep("^link_score_", names(edges_all_tidy), value = TRUE)
cond_cols_expr <- grep("^gene_expr_", names(edges_all_tidy), value = TRUE)
edges_cond_long <- list()
if (length(cond_cols_link)) {
  edges_cond_long$link <- tidyr::pivot_longer(
    edges_all_tidy,
    cols = all_of(cond_cols_link),
    names_prefix = "link_score_",
    names_to = "condition",
    values_to = "w"
  ) |> dplyr::filter(is.finite(w) & w > 0)
  if (nrow(edges_cond_long$link)) {
    edges_cond_long$link$doc_id <- paste(edges_cond_long$link$condition, edges_cond_long$link$tf, sep = "::")
  }
}
if (length(cond_cols_expr)) {
  edges_cond_long$gene_expr <- tidyr::pivot_longer(
    edges_all_tidy,
    cols = all_of(cond_cols_expr),
    names_prefix = "gene_expr_",
    names_to = "condition",
    values_to = "w"
  ) |> dplyr::filter(is.finite(w) & w > 0)
  if (nrow(edges_cond_long$gene_expr)) {
    edges_cond_long$gene_expr$doc_id <- paste(edges_cond_long$gene_expr$condition, edges_cond_long$gene_expr$tf, sep = "::")
  }
}
message("[per-condition] link rows: ", ifelse(is.null(edges_cond_long$link), 0, nrow(edges_cond_long$link)),
        "; gene expr rows: ", ifelse(is.null(edges_cond_long$gene_expr), 0, nrow(edges_cond_long$gene_expr)))

# LDA/NMF benchmark output folder
benchmark_dir <- file.path(base_dir, "benchmark_topics_lda_nmf")
dir.create(benchmark_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Fast drop-in methods (T2/T3 + extras) ---------------------------------
topic_stats_from_edges <- function(topic_edge_tbl, markers_epi, markers_mes, top_n_overlap = 50L) {
  topic_gene <- unique(topic_edge_tbl[, c("topic", "gene_key")])
  topic_tf <- unique(topic_edge_tbl[, c("topic", "tf")])
  n_genes <- topic_gene |>
    dplyr::group_by(topic) |>
    dplyr::summarise(n_genes = dplyr::n_distinct(gene_key), .groups = "drop")
  n_tfs <- topic_tf |>
    dplyr::group_by(topic) |>
    dplyr::summarise(n_tfs = dplyr::n_distinct(tf), .groups = "drop")

  universe_genes <- unique(topic_gene$gene_key)
  N <- length(universe_genes)

  calc_stats <- function(topic_genes, marker_vec) {
    topic_genes <- unique(topic_genes)
    marker_vec <- unique(marker_vec)
    topic_u <- intersect(topic_genes, universe_genes)
    marker_u <- intersect(marker_vec, universe_genes)
    a <- length(intersect(topic_u, marker_u))
    b <- length(setdiff(topic_u, marker_u))
    c <- length(setdiff(marker_u, topic_u))
    d <- N - a - b - c
    d <- max(d, 0)
    ft <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    list(
      overlap_n = a,
      marker_coverage = if (length(marker_u)) a / length(marker_u) else NA_real_,
      fisher_p = unname(ft$p.value),
      fisher_or = unname(ft$estimate)
    )
  }

  topics <- sort(unique(topic_gene$topic))
  res_list <- lapply(topics, function(k) {
    genes_k <- topic_gene$gene_key[topic_gene$topic == k]
    epi <- calc_stats(genes_k, markers_epi)
    mes <- calc_stats(genes_k, markers_mes)
    overlap_epi_genes <- intersect(genes_k, markers_epi)
    overlap_mes_genes <- intersect(genes_k, markers_mes)
    tibble::tibble(
      topic = k,
      overlap_epi_n = epi$overlap_n,
      epi_marker_coverage = epi$marker_coverage,
      epi_fisher_p = epi$fisher_p,
      epi_fisher_or = as.numeric(epi$fisher_or),
      overlap_epi_genes = paste(utils::head(sort(overlap_epi_genes), top_n_overlap), collapse = ";"),
      overlap_mes_n = mes$overlap_n,
      mes_marker_coverage = mes$marker_coverage,
      mes_fisher_p = mes$fisher_p,
      mes_fisher_or = as.numeric(mes$fisher_or),
      overlap_mes_genes = paste(utils::head(sort(overlap_mes_genes), top_n_overlap), collapse = ";")
    )
  })

  stats_tbl <- dplyr::bind_rows(res_list) |>
    dplyr::left_join(n_genes, by = "topic") |>
    dplyr::left_join(n_tfs, by = "topic")

  stats_tbl$genes_full <- vapply(
    stats_tbl$topic,
    function(k) paste(sort(unique(topic_gene$gene_key[topic_gene$topic == k])), collapse = ";"),
    character(1)
  )
  stats_tbl$tf_full <- vapply(
    stats_tbl$topic,
    function(k) paste(sort(unique(topic_tf$tf[topic_tf$topic == k])), collapse = ";"),
    character(1)
  )

  stats_tbl$purity_epi <- ifelse(stats_tbl$n_genes > 0, stats_tbl$overlap_epi_n / stats_tbl$n_genes, NA_real_)
  stats_tbl$purity_mes <- ifelse(stats_tbl$n_genes > 0, stats_tbl$overlap_mes_n / stats_tbl$n_genes, NA_real_)
  stats_tbl$purity_majority <- pmax(stats_tbl$purity_epi, stats_tbl$purity_mes, na.rm = TRUE)

  best_epi <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_epi_n), epi_fisher_p, dplyr::desc(epi_marker_coverage)) |>
    dplyr::slice(1)
  best_mes <- stats_tbl |>
    dplyr::arrange(dplyr::desc(overlap_mes_n), mes_fisher_p, dplyr::desc(mes_marker_coverage)) |>
    dplyr::slice(1)

  list(stats = stats_tbl, best_epi = best_epi, best_mes = best_mes)
}

if (isTRUE(run_fast_methods)) {
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("igraph", quietly = TRUE)) {
    stop("Fast methods need data.table and igraph.")
  }

  dt_fast <- data.table::as.data.table(edges_filtered_tidy)
  dt_fast <- dt_fast[is.finite(delta_link_score)]
  dt_fast[, w := abs(delta_link_score)]
  dt_fast <- dt_fast[w > 0]
  dt_fast[, edge_id := paste(tf, gene_key, sep = "|")]
  agg <- dt_fast[, .(w = sum(w)), by = .(edge_id, tf, gene_key, comparison_id)]
  wide <- data.table::dcast(agg, edge_id + tf + gene_key ~ comparison_id, value.var = "w", fill = 0)

  edge_meta <- wide[, .(edge_id, tf, gene_key)]
  X <- as.matrix(wide[, setdiff(names(wide), c("edge_id", "tf", "gene_key")), with = FALSE])
  rownames(X) <- edge_meta$edge_id

  # T2: PCA + k-means (share_topic label per request)
  pca_k <- c(20, 30, 40)
  pc_n <- min(10, ncol(X))
  Xp <- stats::prcomp(X, center = TRUE, scale. = TRUE)$x[, seq_len(pc_n), drop = FALSE]
  for (k in pca_k) {
    km <- stats::kmeans(Xp, centers = k, nstart = 10)
    topic_edge <- tibble::tibble(
      topic = as.integer(km$cluster),
      tf = edge_meta$tf,
      gene_key = edge_meta$gene_key
    )
    res <- topic_stats_from_edges(topic_edge, markers_epi, markers_mes)
    out_stats <- res$stats
    out_stats$method <- "share_topic"
    out_stats$K <- k
    out_path <- file.path(benchmark_dir, sprintf("share_topic_K%d_summary.csv", k))
    readr::write_csv(out_stats, out_path)
  }

  # T2: raw k-means
  for (k in pca_k) {
    km <- stats::kmeans(X, centers = k, nstart = 10)
    topic_edge <- tibble::tibble(
      topic = as.integer(km$cluster),
      tf = edge_meta$tf,
      gene_key = edge_meta$gene_key
    )
    res <- topic_stats_from_edges(topic_edge, markers_epi, markers_mes)
    out_stats <- res$stats
    out_stats$method <- "share_topic_kmeans_raw"
    out_stats$K <- k
    out_path <- file.path(benchmark_dir, sprintf("share_topic_kmeans_raw_K%d_summary.csv", k))
    readr::write_csv(out_stats, out_path)
  }

  # T3: Louvain on TF-gene bipartite graph
  g <- igraph::graph_from_data_frame(
    dt_fast[, .(tf, gene_key, w)],
    directed = FALSE
  )
  igraph::E(g)$weight <- dt_fast$w
  comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  memb <- igraph::membership(comm)
  topic_edge <- tibble::tibble(
    topic = as.integer(memb[dt_fast$gene_key]),
    tf = dt_fast$tf,
    gene_key = dt_fast$gene_key
  )
  topic_edge <- topic_edge[!is.na(topic_edge$topic), , drop = FALSE]
  res <- topic_stats_from_edges(topic_edge, markers_epi, markers_mes)
  out_stats <- res$stats
  out_stats$method <- "louvain_bipartite"
  out_path <- file.path(benchmark_dir, "louvain_bipartite_summary.csv")
  readr::write_csv(out_stats, out_path)

  share_paths <- list.files(benchmark_dir, pattern = "^share_topic_.*_summary\\.csv$", full.names = TRUE)
  if (length(share_paths)) {
    share_all <- dplyr::bind_rows(lapply(share_paths, readr::read_csv, show_col_types = FALSE))
    share_all_path <- file.path(benchmark_dir, "share_topic_all_runs_summary.csv")
    readr::write_csv(share_all, share_all_path)
    message("[share_topic combined] ", share_all_path)
  }

  louvain_path <- file.path(benchmark_dir, "louvain_bipartite_summary.csv")
  if (file.exists(louvain_path)) {
    louvain_all_path <- file.path(benchmark_dir, "louvain_all_runs_summary.csv")
    readr::write_csv(readr::read_csv(louvain_path, show_col_types = FALSE), louvain_all_path)
    message("[louvain combined] ", louvain_all_path)
  }
}

# LDA benchmark: delta vs FC (combined)
variants <- list(
  list(label = "delta", docs = dt_docs_delta, edges = edges_filtered_tidy),
  list(label = "fc",    docs = dt_docs_fc,    edges = edges_fc)
)

lda_K <- c(20, 30, 40, 50, 60, 80, 100)
lda_methods <- c("VEM", "Gibbs")
gamma_cuts <- c(0.2, 0.5, 0.8)

augment_topic_stats <- function(stats_tbl, topic_entity_tbl, doc_topics, gene_col, tf_source_dt = NULL) {
  genes_by_topic <- topic_entity_tbl |>
    dplyr::group_by(topic) |>
    dplyr::summarise(
      genes_full = paste(sort(unique(.data[[gene_col]])), collapse = ";"),
      .groups = "drop"
    )

  tf_full_tbl <- NULL
  if (!is.null(doc_topics) && "topics_str" %in% names(doc_topics) && "tf" %in% names(doc_topics)) {
    tf_full_tbl <- doc_topics |>
      dplyr::filter(!is.na(topics_str) & topics_str != "") |>
      tidyr::separate_rows(topics_str, sep = ";") |>
      dplyr::mutate(topic = suppressWarnings(as.integer(gsub("[^0-9]", "", topics_str)))) |>
      dplyr::filter(is.finite(topic)) |>
      dplyr::group_by(topic) |>
      dplyr::summarise(tf_full = paste(sort(unique(tf)), collapse = ";"), .groups = "drop")
  }

  stats_aug <- stats_tbl
  stats_aug <- dplyr::left_join(stats_aug, genes_by_topic, by = "topic")
  if (!is.null(tf_full_tbl)) {
    stats_aug <- dplyr::left_join(stats_aug, tf_full_tbl, by = "topic")
  }
  stats_aug <- stats_aug |>
    dplyr::mutate(
      purity_epi = ifelse(is.finite(n_genes) & n_genes > 0, overlap_epi_n / n_genes, NA_real_),
      purity_mes = ifelse(is.finite(n_genes) & n_genes > 0, overlap_mes_n / n_genes, NA_real_),
      purity_majority = pmax(purity_epi, purity_mes, na.rm = TRUE)
    )
  stats_aug
}

write_best_topic_csv <- function(best_row, topic_entity_tbl, doc_topics, gene_col, tf_source_dt = NULL,
                                 label, out_csv) {
  topic_id <- as.integer(best_row$topic)
  topic_size <- as.integer(best_row$n_genes %||% nrow(topic_entity_tbl[topic_entity_tbl[[gene_col]] != "", , drop = FALSE]))
  if (!is.finite(topic_size) || topic_size <= 0) topic_size <- nrow(topic_entity_tbl)

  genes_topic <- unique(as.character(topic_entity_tbl[[gene_col]][topic_entity_tbl$topic == topic_id]))
  genes_topic <- genes_topic[!is.na(genes_topic) & genes_topic != ""]

  tf_list <- character(0)
  if (!is.null(tf_source_dt)) {
    if ("topics_str" %in% names(doc_topics) && "tf" %in% names(doc_topics)) {
      topics_tf <- doc_topics[!is.na(doc_topics$topics_str) & grepl(paste0("\\b", topic_id, "\\b"), doc_topics$topics_str), , drop = FALSE]
      tf_list <- unique(as.character(topics_tf$tf))
    } else if (all(c("doc_id", "tf") %in% names(tf_source_dt))) {
      tf_list <- unique(as.character(tf_source_dt$tf[tf_source_dt$doc_id %in% genes_topic]))
    }
  }

  purity_epi <- if (is.finite(topic_size) && topic_size > 0) as.numeric(best_row$overlap_epi_n) / topic_size else NA_real_
  purity_mes <- if (is.finite(topic_size) && topic_size > 0) as.numeric(best_row$overlap_mes_n) / topic_size else NA_real_
  purity_major <- max(purity_epi, purity_mes, na.rm = TRUE)

  df_out <- tibble::tibble(
    label = label,
    topic = topic_id,
    topic_size = topic_size,
    overlap_epi_n = best_row$overlap_epi_n,
    overlap_mes_n = best_row$overlap_mes_n,
    purity_epi = purity_epi,
    purity_mes = purity_mes,
    purity_majority = purity_major,
    epi_marker_coverage = best_row$epi_marker_coverage,
    mes_marker_coverage = best_row$mes_marker_coverage,
    epi_fisher_p = best_row$epi_fisher_p,
    mes_fisher_p = best_row$mes_fisher_p,
    overlap_epi_genes = best_row$overlap_epi_genes,
    overlap_mes_genes = best_row$overlap_mes_genes,
    genes_full = paste(sort(genes_topic), collapse = ";"),
    tf_full = paste(sort(tf_list), collapse = ";")
  )

  readr::write_csv(df_out, out_csv)
  message("[saved] ", out_csv)
  invisible(df_out)
}

grid <- expand.grid(
  gamma = gamma_cuts,
  v_idx = seq_along(variants),
  K = lda_K,
  method = lda_methods,
  stringsAsFactors = FALSE
)

run_one <- function(row) {
  g <- row$gamma
  gamma_tag <- sprintf("g%02d", as.integer(round(g * 100)))
  v <- variants[[row$v_idx]]
  k <- row$K
  m <- row$method

  summary_path <- file.path(benchmark_dir, sprintf("lda_%s_K%d_%s_%s_summary.csv", v$label, k, m, gamma_tag))
  best_epi_csv <- file.path(benchmark_dir, sprintf("lda_%s_K%d_%s_%s_best_epi.csv", v$label, k, m, gamma_tag))
  best_mes_csv <- file.path(benchmark_dir, sprintf("lda_%s_K%d_%s_%s_best_mes.csv", v$label, k, m, gamma_tag))
  if (!rerun && file.exists(summary_path)) {
    message("[skip existing] ", summary_path)
    return(NULL)
  }

  lda_fit <- fit_pooled_lda(
    v$docs,
    K = k,
    method = m,
    seed = 1,
    gamma_cutoff = g
  )

  out_curr <- topic_marker_overlap(
    edges_tbl   = v$edges,
    doc_topics  = lda_fit$doc_topics,
    markers_epi = markers_epi,
    markers_mes = markers_mes,
    abs_delta_min = 1,
    top_n_overlap_genes = 50
  )

  stats_aug <- augment_topic_stats(out_curr$stats, out_curr$topic_gene, lda_fit$doc_topics, "gene_key", v$docs)

  readr::write_csv(stats_aug, best_epi_csv)
  readr::write_csv(stats_aug, best_mes_csv)

  n_e <- as.numeric(out_curr$best_epi$n_genes %||% NA)
  n_m <- as.numeric(out_curr$best_mes$n_genes %||% NA)
  purity_e <- if (is.finite(n_e) && n_e > 0) out_curr$best_epi$overlap_epi_n / n_e else NA_real_
  purity_m <- if (is.finite(n_m) && n_m > 0) out_curr$best_mes$overlap_mes_n / n_m else NA_real_

  best_epi_full <- dplyr::filter(stats_aug, topic == out_curr$best_epi$topic)
  best_mes_full <- dplyr::filter(stats_aug, topic == out_curr$best_mes$topic)

  epi_overlap_epi_genes <- if (nrow(best_epi_full)) best_epi_full$overlap_epi_genes[[1]] else NA_character_
  epi_overlap_mes_genes <- if (nrow(best_epi_full)) best_epi_full$overlap_mes_genes[[1]] else NA_character_
  epi_genes_full <- if (nrow(best_epi_full)) best_epi_full$genes_full[[1]] else NA_character_
  epi_tf_full <- if (nrow(best_epi_full) && "tf_full" %in% names(best_epi_full)) best_epi_full$tf_full[[1]] else NA_character_
  epi_n_genes <- if (nrow(best_epi_full) && "n_genes" %in% names(best_epi_full)) best_epi_full$n_genes[[1]] else NA_real_
  epi_n_tfs <- if (nrow(best_epi_full) && "n_tfs" %in% names(best_epi_full)) best_epi_full$n_tfs[[1]] else NA_real_

  mes_overlap_epi_genes <- if (nrow(best_mes_full)) best_mes_full$overlap_epi_genes[[1]] else NA_character_
  mes_overlap_mes_genes <- if (nrow(best_mes_full)) best_mes_full$overlap_mes_genes[[1]] else NA_character_
  mes_genes_full <- if (nrow(best_mes_full)) best_mes_full$genes_full[[1]] else NA_character_
  mes_tf_full <- if (nrow(best_mes_full) && "tf_full" %in% names(best_mes_full)) best_mes_full$tf_full[[1]] else NA_character_
  mes_n_genes <- if (nrow(best_mes_full) && "n_genes" %in% names(best_mes_full)) best_mes_full$n_genes[[1]] else NA_real_
  mes_n_tfs <- if (nrow(best_mes_full) && "n_tfs" %in% names(best_mes_full)) best_mes_full$n_tfs[[1]] else NA_real_

  summary_row <- tibble::tibble(
    label = v$label,
    K = k,
    method = m,
    gamma_cutoff = g,
    gamma_tag = gamma_tag,
    best_epi_topic = out_curr$best_epi$topic,
    best_epi_overlap = out_curr$best_epi$overlap_epi_n,
    best_epi_purity = purity_e,
    best_epi_overlap_genes = epi_overlap_epi_genes,
    best_epi_overlap_mes_genes = epi_overlap_mes_genes,
    best_epi_genes_full = epi_genes_full,
    best_epi_tf_full = epi_tf_full,
    best_epi_n_genes = epi_n_genes,
    best_epi_n_tfs = epi_n_tfs,
    best_mes_topic = out_curr$best_mes$topic,
    best_mes_overlap = out_curr$best_mes$overlap_mes_n,
    best_mes_purity = purity_m,
    best_mes_overlap_genes = mes_overlap_epi_genes,
    best_mes_overlap_mes_genes = mes_overlap_mes_genes,
    best_mes_genes_full = mes_genes_full,
    best_mes_tf_full = mes_tf_full,
    best_mes_n_genes = mes_n_genes,
    best_mes_n_tfs = mes_n_tfs,
    summary_csv = summary_path,
    best_epi_csv = best_epi_csv,
    best_mes_csv = best_mes_csv
  )

  readr::write_csv(summary_row, summary_path)
  message("[saved] ", summary_path)
  summary_row
}

grid_list <- split(grid, seq_len(nrow(grid)))
summary_rows <- tibble::tibble()
if (isTRUE(run_lda_benchmark)) {
  res_list <- NULL
  if (use_parallel && .Platform$OS.type != "windows") {
    lda_cores <- max(1, min(length(grid_list), bench_cores))
    message("[LDA] using mc.cores=", lda_cores, " of ", parallel::detectCores())
    res_list <- parallel::mclapply(grid_list, run_one, mc.cores = lda_cores, mc.preschedule = FALSE)
  } else if (use_parallel && .Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("variants", "markers_epi", "markers_mes", "benchmark_dir", "run_one",
                  "gamma_cuts", "lda_K", "lda_methods", "rerun", "%||%",
                  "augment_topic_stats", "fit_pooled_lda", "topic_marker_overlap"),
      envir = environment()
    )
    res_list <- parallel::parLapply(cl, grid_list, run_one)
  } else {
    res_list <- lapply(grid_list, run_one)
  }
  summary_rows <- dplyr::bind_rows(Filter(Negate(is.null), res_list))
}

# Combine summaries
if (isTRUE(run_lda_benchmark)) {
  sum_paths <- list.files(benchmark_dir, pattern = "^lda_.*_summary\\.csv$", full.names = TRUE)
  if (length(sum_paths)) {
    all_sum <- dplyr::bind_rows(lapply(sum_paths, readr::read_csv, show_col_types = FALSE))
    all_sum_path <- file.path(benchmark_dir, "lda_all_runs_summary.csv")
    readr::write_csv(all_sum, all_sum_path)
    master_path <- file.path(benchmark_dir, "lda_all_runs_master.csv")
    if (file.exists(master_path)) {
      master_old <- readr::read_csv(master_path, show_col_types = FALSE)
      all_sum <- dplyr::distinct(dplyr::bind_rows(master_old, all_sum))
    }
    readr::write_csv(all_sum, master_path)
    message("[combined] ", all_sum_path)
    message("[master updated] ", master_path)
  }
}

# ---- NMF benchmark (delta, fc; combined docs) -------------------------------
nmf_variants <- list(
  list(label = "delta", docs = dt_docs_delta, edges = edges_filtered_tidy),
  list(label = "fc",    docs = dt_docs_fc,    edges = edges_fc)
)
nmf_K <- c(20, 30)
nmf_methods <- c("lee")
nmf_gamma <- c(0.3, 0.5)

nmf_grid <- expand.grid(
  gamma = nmf_gamma,
  v_idx = seq_along(nmf_variants),
  K = nmf_K,
  method = nmf_methods,
  stringsAsFactors = FALSE
)

run_one_nmf <- function(row) {
  g <- row$gamma
  gamma_tag <- sprintf("g%02d", as.integer(round(g * 100)))
  v <- nmf_variants[[row$v_idx]]
  k <- row$K
  m <- row$method

  summary_path <- file.path(benchmark_dir, sprintf("nmf_%s_K%d_%s_%s_summary.csv", v$label, k, m, gamma_tag))
  best_epi_csv <- file.path(benchmark_dir, sprintf("nmf_%s_K%d_%s_%s_best_epi.csv", v$label, k, m, gamma_tag))
  best_mes_csv <- file.path(benchmark_dir, sprintf("nmf_%s_K%d_%s_%s_best_mes.csv", v$label, k, m, gamma_tag))
  if (!rerun && file.exists(summary_path)) {
    message("[skip existing] ", summary_path)
    return(NULL)
  }

  # keep docs small for stability
  docs_run <- v$docs |>
    dplyr::group_by(doc_id) |>
    dplyr::slice_max(order_by = w, n = 500, with_ties = FALSE) |>
    dplyr::ungroup()

  nmf_fit <- fit_pooled_nmf(
    docs_run,
    K = k,
    method = m,
    seed = 1,
    nrun = 1,
    gamma_cutoff = g,
    nmf_stop = "none",
    nmf_options = "-p1"
  )

  out_curr <- topic_marker_overlap(
    edges_tbl   = v$edges,
    doc_topics  = nmf_fit$doc_topics,
    markers_epi = markers_epi,
    markers_mes = markers_mes,
    abs_delta_min = 1,
    top_n_overlap_genes = 50
  )

  stats_aug <- augment_topic_stats(out_curr$stats, out_curr$topic_gene, nmf_fit$doc_topics, "gene_key", docs_run)

  readr::write_csv(stats_aug, best_epi_csv)
  readr::write_csv(stats_aug, best_mes_csv)

  n_e <- as.numeric(out_curr$best_epi$n_genes %||% NA)
  n_m <- as.numeric(out_curr$best_mes$n_genes %||% NA)
  purity_e <- if (is.finite(n_e) && n_e > 0) out_curr$best_epi$overlap_epi_n / n_e else NA_real_
  purity_m <- if (is.finite(n_m) && n_m > 0) out_curr$best_mes$overlap_mes_n / n_m else NA_real_

  best_epi_full <- dplyr::filter(stats_aug, topic == out_curr$best_epi$topic)
  best_mes_full <- dplyr::filter(stats_aug, topic == out_curr$best_mes$topic)

  epi_overlap_epi_genes <- if (nrow(best_epi_full)) best_epi_full$overlap_epi_genes[[1]] else NA_character_
  epi_overlap_mes_genes <- if (nrow(best_epi_full)) best_epi_full$overlap_mes_genes[[1]] else NA_character_
  epi_genes_full <- if (nrow(best_epi_full)) best_epi_full$genes_full[[1]] else NA_character_
  epi_tf_full <- if (nrow(best_epi_full) && "tf_full" %in% names(best_epi_full)) best_epi_full$tf_full[[1]] else NA_character_

  mes_overlap_epi_genes <- if (nrow(best_mes_full)) best_mes_full$overlap_epi_genes[[1]] else NA_character_
  mes_overlap_mes_genes <- if (nrow(best_mes_full)) best_mes_full$overlap_mes_genes[[1]] else NA_character_
  mes_genes_full <- if (nrow(best_mes_full)) best_mes_full$genes_full[[1]] else NA_character_
  mes_tf_full <- if (nrow(best_mes_full) && "tf_full" %in% names(best_mes_full)) best_mes_full$tf_full[[1]] else NA_character_

  summary_row <- tibble::tibble(
    label = v$label,
    K = k,
    method = m,
    gamma_cutoff = g,
    gamma_tag = gamma_tag,
    best_epi_topic = out_curr$best_epi$topic,
    best_epi_overlap = out_curr$best_epi$overlap_epi_n,
    best_epi_purity = purity_e,
    best_epi_overlap_genes = epi_overlap_epi_genes,
    best_epi_overlap_mes_genes = epi_overlap_mes_genes,
    best_epi_genes_full = epi_genes_full,
    best_epi_tf_full = epi_tf_full,
    best_mes_topic = out_curr$best_mes$topic,
    best_mes_overlap = out_curr$best_mes$overlap_mes_n,
    best_mes_purity = purity_m,
    best_mes_overlap_genes = mes_overlap_epi_genes,
    best_mes_overlap_mes_genes = mes_overlap_mes_genes,
    best_mes_genes_full = mes_genes_full,
    best_mes_tf_full = mes_tf_full,
    summary_csv = summary_path,
    best_epi_csv = best_epi_csv,
    best_mes_csv = best_mes_csv
  )

  readr::write_csv(summary_row, summary_path)
  message("[saved] ", summary_path)
  summary_row
}

nmf_grid_list <- split(nmf_grid, seq_len(nrow(nmf_grid)))
nmf_summary_rows <- tibble::tibble()
if (isTRUE(run_nmf_benchmark)) {
  nmf_res_list <- NULL
  if (use_parallel && .Platform$OS.type != "windows") {
    nmf_cores <- max(1, bench_cores)
    message("[NMF] using mc.cores=", nmf_cores, " of ", parallel::detectCores(),
            " (jobs=", length(nmf_grid_list), ")")
    nmf_res_list <- parallel::mclapply(nmf_grid_list, run_one_nmf, mc.cores = nmf_cores, mc.preschedule = FALSE)
  } else if (use_parallel && .Platform$OS.type == "windows") {
    cl_nmf <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
    on.exit(parallel::stopCluster(cl_nmf), add = TRUE)
    parallel::clusterExport(
      cl_nmf,
      varlist = c("nmf_variants", "markers_epi", "markers_mes", "benchmark_dir", "run_one_nmf",
                  "rerun", "%||%", "augment_topic_stats", "fit_pooled_nmf", "topic_marker_overlap"),
      envir = environment()
    )
    nmf_res_list <- parallel::parLapply(cl_nmf, nmf_grid_list, run_one_nmf)
  } else {
    nmf_res_list <- lapply(nmf_grid_list, run_one_nmf)
  }

  nmf_summary_rows <- dplyr::bind_rows(Filter(Negate(is.null), nmf_res_list))
  nmf_sum_paths <- list.files(benchmark_dir, pattern = "^nmf_.*_summary\\.csv$", full.names = TRUE)
  if (length(nmf_sum_paths)) {
    nmf_all_sum <- dplyr::bind_rows(lapply(nmf_sum_paths, readr::read_csv, show_col_types = FALSE))
    nmf_all_sum_path <- file.path(benchmark_dir, "nmf_all_runs_summary.csv")
    readr::write_csv(nmf_all_sum, nmf_all_sum_path)
    nmf_master_path <- file.path(benchmark_dir, "nmf_all_runs_master.csv")
    if (file.exists(nmf_master_path)) {
      nmf_master_old <- readr::read_csv(nmf_master_path, show_col_types = FALSE)
      nmf_all_sum <- dplyr::distinct(dplyr::bind_rows(nmf_master_old, nmf_all_sum))
    }
    readr::write_csv(nmf_all_sum, nmf_master_path)
    message("[combined NMF] ", nmf_all_sum_path)
    message("[master NMF updated] ", nmf_master_path)
  }
}

# ---- Combined summary across all methods -----------------------------------
combine_all_methods_summary <- function(benchmark_dir) {
  paths <- list(
    lda = file.path(benchmark_dir, "lda_all_runs_summary.csv"),
    nmf = file.path(benchmark_dir, "nmf_all_runs_summary.csv"),
    share_topic = file.path(benchmark_dir, "share_topic_all_runs_summary.csv"),
    louvain = file.path(benchmark_dir, "louvain_all_runs_summary.csv")
  )
  tbls <- list()
  for (nm in names(paths)) {
    p <- paths[[nm]]
    if (file.exists(p)) {
      x <- readr::read_csv(p, show_col_types = FALSE)
      if (!"method_group" %in% names(x)) x$method_group <- nm
      tbls[[nm]] <- x
    }
  }
  if (!length(tbls)) return(invisible(NULL))
  common_cols <- Reduce(intersect, lapply(tbls, names))
  all_tbl <- dplyr::bind_rows(lapply(tbls, function(x) x[, common_cols]))
  out_path <- file.path(benchmark_dir, "all_methods_summary.csv")
  readr::write_csv(all_tbl, out_path)
  message("[combined all methods] ", out_path)
  invisible(all_tbl)
}

combine_all_methods_summary(benchmark_dir)

# ---- NMF demos --------------------------------------------------------------
if (isTRUE(run_nmf_demo)) {
  message("[nmf demo] running on tiny synthetic matrix to verify NMF workflow...")
  demo_edges <- tibble::tribble(
    ~comparison_id, ~tf, ~gene_key, ~delta_link_score,
    "C1", "TF1", "G1", 2,
    "C1", "TF1", "G2", 1.5,
    "C1", "TF2", "G1", 1.2,
    "C2", "TF1", "G3", 1.8,
    "C2", "TF2", "G2", 2.2,
    "C2", "TF2", "G3", 1.7
  )
  demo_docs <- build_pooled_doc_term(
    demo_edges,
    which = "abs",
    min_abs_delta = 0,
    top_terms_per_doc = 100,
    min_df = 1
  )
  demo_nmf <- fit_pooled_nmf(
    demo_docs,
    K = 2,
    method = "lee",
    seed = 1,
    nrun = 1,
    gamma_cutoff = 0.3,
    nmf_stop = "none",
    nmf_options = "-p1-t"
  )
  message("[nmf demo] topics:"); print(demo_nmf$doc_topics)
}

if (isTRUE(run_nmf_subset)) {
  message("[nmf subset] running NMF on pruned real docs for a quick sanity check...")
  dt_nmf <- dt_docs_delta |>
    dplyr::group_by(doc_id) |>
    dplyr::slice_max(order_by = w, n = 200, with_ties = FALSE) |>
    dplyr::ungroup()

  nmf_subset <- fit_pooled_nmf(
    dt_nmf,
    K = 20,
    method = "lee",
    seed = 1,
    nrun = 1,
    gamma_cutoff = 0.5,
    nmf_stop = "none",
    nmf_options = "-p1"
  )

  nmf_outdir <- file.path(benchmark_dir, "nmf_subset")
  dir.create(nmf_outdir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(nmf_subset$doc_topics, file.path(nmf_outdir, "nmf_subset_doc_topics.csv"))
  readr::write_csv(nmf_subset$doc_topic_long, file.path(nmf_outdir, "nmf_subset_doc_topic_long.csv"))
  message("[nmf subset] saved NMF outputs to ", nmf_outdir)
}
