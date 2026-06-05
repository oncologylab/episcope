# File: utils_step2_reports.R
# Purpose: Optional Module 2 HTML reports from relational outputs.


.module2_report_table_row <- function(module2, table_name) {
  man <- module2$manifest
  if (!is.data.frame(man) || !nrow(man)) return(NULL)
  hit <- man[as.character(man$table) == table_name, , drop = FALSE]
  if (!nrow(hit)) return(NULL)
  hit[1L, , drop = FALSE]
}

.module2_report_read_table <- function(module2, table_name, columns = NULL) {
  direct_map <- c(
    module1_predicted_tfbs = "predicted_tfbs",
    module2_fp_target_candidates = "candidates",
    module2_tf_target_corr = "tf_target_corr",
    module2_fp_target_corr = "fp_target_corr",
    module2_links = "links",
    module2_condition_activity = "condition_activity",
    module2_qc_summary = "qc_summary"
  )
  direct_name <- if (table_name %in% names(direct_map)) unname(direct_map[[table_name]]) else table_name
  x <- module2[[direct_name]]
  if (is.data.frame(x) && nrow(x) && !all(c("path", "format") %in% names(x))) {
    if (!is.null(columns)) x <- x[, intersect(columns, names(x)), drop = FALSE]
    return(tibble::as_tibble(x))
  }
  man_name <- paste0(table_name, "_manifest")
  chunk_man <- module2[[man_name]]
  if (is.data.frame(chunk_man) && all(c("path", "format") %in% names(chunk_man))) {
    rows <- lapply(seq_len(nrow(chunk_man)), function(i) {
      .module2_read_predicted_chunk(
        as.character(chunk_man$path[[i]]),
        as.character(chunk_man$format[[i]]),
        columns = columns
      )
    })
    return(tibble::as_tibble(data.table::rbindlist(rows, fill = TRUE)))
  }
  row <- .module2_report_table_row(module2, table_name)
  if (!is.null(row)) {
    path <- as.character(row$path[[1L]])
    fmt <- as.character(row$format[[1L]])
    if (identical(fmt, "manifest")) {
      chunk_man <- readr::read_csv(path, show_col_types = FALSE)
      rows <- lapply(seq_len(nrow(chunk_man)), function(i) {
        .module2_read_predicted_chunk(
          as.character(chunk_man$path[[i]]),
          as.character(chunk_man$format[[i]]),
          columns = columns
        )
      })
      return(tibble::as_tibble(data.table::rbindlist(rows, fill = TRUE)))
    }
    return(.module2_read_predicted_chunk(path, fmt, columns = columns))
  }
  .log_abort("Module 2 table not found: {table_name}")
}

.module2_report_read_links <- function(module2, columns = NULL, tf = NULL) {
  if (is.null(tf)) {
    return(.module2_report_read_table(module2, "module2_links", columns = columns))
  }
  query_module2_links(module2, tf = tf, pass_only = TRUE)
}

.module2_report_join_links <- function(module2, tf = NULL) {
  tf_expression_target_r <- fp_target_rna_r <- NULL
  links <- .module2_report_read_links(
    module2,
    columns = c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass"),
    tf = tf
  )
  links <- links[links$module2_link_pass %in% TRUE, , drop = FALSE]
  if (!nrow(links)) return(data.table::data.table())
  candidates <- .module2_report_read_table(
    module2,
    "module2_fp_target_candidates",
    columns = c("candidate_id", "fp_id", "target_gene", "chr", "start", "end", "atac_peak", "distance_to_tss")
  )
  tf_corr <- .module2_report_read_table(
    module2,
    "module2_tf_target_corr",
    columns = c("tf", "target_gene", "best_r", "pass")
  )
  fp_corr <- .module2_report_read_table(
    module2,
    "module2_fp_target_corr",
    columns = c("fp_id", "target_gene", "best_r", "pass")
  )
  out <- data.table::as.data.table(links)
  cand_dt <- data.table::as.data.table(candidates)
  tf_dt <- data.table::as.data.table(tf_corr)
  fp_dt <- data.table::as.data.table(fp_corr)
  out <- cand_dt[out, on = "candidate_id", nomatch = 0L]
  out <- tf_dt[out, on = c("tf", "target_gene"), nomatch = 0L]
  data.table::setnames(out, "best_r", "tf_expression_target_r", skip_absent = TRUE)
  out <- fp_dt[out, on = c("fp_id", "target_gene"), nomatch = 0L]
  data.table::setnames(out, "best_r", "fp_target_rna_r", skip_absent = TRUE)
  if (!"tf_expression_target_r" %in% names(out)) out[, tf_expression_target_r := NA_real_]
  if (!"fp_target_rna_r" %in% names(out)) out[, fp_target_rna_r := NA_real_]
  out[, `:=`(
    tf = toupper(as.character(tf)),
    gene_norm = toupper(as.character(target_gene)),
    gene_symbol = as.character(target_gene),
    peak_ID = as.character(fp_id),
    fp_target_rna_r = suppressWarnings(as.numeric(fp_target_rna_r)),
    tf_expression_target_r = suppressWarnings(as.numeric(tf_expression_target_r))
  )]
  out[]
}


.module2_report_direct_links <- function(module2) {
  tf_expression_target_r <- fp_target_rna_r <- NULL
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  tf_corr <- .module2_report_read_table(
    module2,
    "module2_tf_target_corr",
    columns = c("tf", "target_gene", "best_r", "pass")
  )
  tf_dt <- data.table::as.data.table(tf_corr)
  tf_dt[, `:=`(tf = toupper(as.character(tf)), target_gene = toupper(as.character(target_gene)))]
  tf_levels <- sort(unique(tf_dt$tf))
  tf_dt <- tf_dt[target_gene %in% tf_levels & pass %in% TRUE]
  if (!nrow(tf_dt)) return(data.table::data.table())
  data.table::setnames(tf_dt, "best_r", "tf_expression_target_r", skip_absent = TRUE)
  keep_pairs <- unique(tf_dt[, .(tf, target_gene)])
  link_man <- module2$module2_links_manifest
  if (!is.data.frame(link_man) || !nrow(link_man)) {
    link_tbl <- .module2_report_read_table(
      module2,
      "module2_links",
      columns = c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass")
    )
    link_rows <- list(data.table::as.data.table(link_tbl))
  } else {
    link_rows <- lapply(seq_len(nrow(link_man)), function(i) {
      one <- .module2_read_predicted_chunk(
        as.character(link_man$path[[i]]),
        as.character(link_man$format[[i]]),
        columns = c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass")
      )
      one <- data.table::as.data.table(one)
      one[, `:=`(tf = toupper(as.character(tf)), target_gene = toupper(as.character(target_gene)))]
      one <- one[module2_link_pass %in% TRUE & target_gene %in% tf_levels & tf != target_gene]
      if (!nrow(one)) return(data.table::data.table())
      one <- keep_pairs[one, on = c("tf", "target_gene"), nomatch = 0L]
      one[]
    })
  }
  links <- data.table::rbindlist(link_rows, use.names = TRUE, fill = TRUE)
  if (!nrow(links)) return(data.table::data.table())
  fp_corr <- .module2_report_read_table(
    module2,
    "module2_fp_target_corr",
    columns = c("fp_id", "target_gene", "best_r", "pass")
  )
  fp_dt <- data.table::as.data.table(fp_corr)
  fp_dt[, target_gene := toupper(as.character(target_gene))]
  fp_dt <- fp_dt[target_gene %in% tf_levels & pass %in% TRUE]
  data.table::setnames(fp_dt, "best_r", "fp_target_rna_r", skip_absent = TRUE)
  out <- tf_dt[links, on = c("tf", "target_gene"), nomatch = 0L]
  out <- fp_dt[out, on = c("fp_id", "target_gene"), nomatch = 0L]
  if (!"tf_expression_target_r" %in% names(out)) out[, tf_expression_target_r := NA_real_]
  if (!"fp_target_rna_r" %in% names(out)) out[, fp_target_rna_r := NA_real_]
  out[, `:=`(
    gene_norm = as.character(target_gene),
    gene_symbol = as.character(target_gene),
    peak_ID = as.character(fp_id),
    fp_target_rna_r = suppressWarnings(as.numeric(fp_target_rna_r)),
    tf_expression_target_r = suppressWarnings(as.numeric(tf_expression_target_r))
  )]
  out[]
}

.module2_report_condition_links <- function(link_dt, multiomic_data, conditions = NULL, fp_score_cutoff = 0) {
  condition_fp_score <- NULL
  if (!nrow(link_dt)) return(data.table::data.table())
  if (is.null(multiomic_data)) {
    out <- data.table::copy(link_dt)
    out[, condition := "all"]
    out[, condition_fp_score := NA_real_]
    return(out)
  }
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  validate_multiomic_object(multiomic_data)
  mats <- multiomic_data$matrices
  all_conditions <- colnames(mats$fp_score)
  if (is.null(conditions)) conditions <- all_conditions
  conditions <- intersect(as.character(conditions), all_conditions)
  if (!length(conditions)) .log_abort("No report conditions overlap multiomic_data condition names.")
  dt <- data.table::copy(link_dt)
  rows <- vector("list", length(conditions))
  fp_idx <- match(dt$fp_id, rownames(mats$fp_score))
  tf_idx <- match(dt$tf, rownames(mats$gene_expr))
  gene_idx <- match(dt$target_gene, rownames(mats$gene_expr))
  for (i in seq_along(conditions)) {
    cc <- conditions[[i]]
    fp_score <- as.numeric(mats$fp_score[fp_idx, cc])
    fp_bound <- as.logical(mats$fp_bound[fp_idx, cc])
    tf_expr <- as.numeric(mats$gene_expr[tf_idx, cc])
    gene_expr <- as.numeric(mats$gene_expr[gene_idx, cc])
    keep <- is.finite(fp_score) & fp_score >= fp_score_cutoff &
      fp_bound %in% TRUE &
      is.finite(tf_expr) & tf_expr > 0 &
      is.finite(gene_expr) & gene_expr > 0
    if (!any(keep)) next
    one <- dt[keep]
    one[, `:=`(
      condition = cc,
      condition_fp_score = fp_score[keep],
      tf_expr_flag = as.integer(tf_expr[keep] > 0),
      target_expr_flag = as.integer(gene_expr[keep] > 0),
      active_link = TRUE
    )]
    rows[[i]] <- one
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.module2_report_write_active_cache <- function(active_dt, out_dir, tag, result_label, filter_tag = "module2") {
  csv_dir <- file.path(out_dir, "csv")
  cache_dir <- file.path(out_dir, "cache", "condition_filtered_links", tag, result_label)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(active_dt)) .log_abort("No active Module 2 links available for report generation.")
  active_path <- file.path(csv_dir, sprintf("%s_significant_active_links_long.csv", tag))
  data.table::fwrite(active_dt, active_path)
  split <- split(active_dt, active_dt$condition)
  manifest <- lapply(names(split), function(condition) {
    condition_tag <- gsub("[^A-Za-z0-9]+", "_", condition)
    condition_tag <- gsub("^_|_$", "", condition_tag)
    path <- file.path(cache_dir, sprintf("%s_%s_%s_conditionFiltered_links_%s.csv", tag, result_label, condition_tag, filter_tag))
    data.table::fwrite(split[[condition]], path)
    data.table::data.table(condition = condition, filtered_links_csv = path)
  })
  manifest_dt <- data.table::rbindlist(manifest, use.names = TRUE)
  data.table::fwrite(
    manifest_dt,
    file.path(csv_dir, sprintf("%s_%s_conditionFiltered_tf_tf_composite_connectivity_manifest.csv", tag, result_label))
  )
  list(active_path = active_path, csv_dir = csv_dir, cache_dir = cache_dir, manifest = manifest_dt)
}

.module2_report_top_tables <- function(link_dt, seed_tf, top_n = 100) {
  edge_score <- abs_edge_score <- fp_target_rna_r <- tf_expression_target_r <- gene_norm <- node_id <- NULL
  seed_tf <- toupper(as.character(seed_tf)[[1L]])
  dt <- data.table::copy(link_dt)
  dt <- dt[tf == seed_tf]
  if (!nrow(dt)) return(NULL)
  dt[, edge_score := fp_target_rna_r * tf_expression_target_r]
  dt[, abs_edge_score := abs(edge_score)]
  target_dt <- dt[, .SD[which.max(abs_edge_score)], by = gene_norm]
  data.table::setorder(target_dt, -abs_edge_score, gene_norm)
  target_dt <- head(target_dt, top_n)
  target_genes <- unique(target_dt$gene_norm)
  extra <- link_dt[gene_norm %in% target_genes]
  extra[, edge_score := fp_target_rna_r * tf_expression_target_r]
  extra[, abs_edge_score := abs(edge_score)]
  edges <- extra[, .(
    from = tf,
    to = gene_norm,
    edge_class = ifelse(tf == seed_tf, "seed_to_target", "tf_to_target"),
    edge_score = max(edge_score, na.rm = TRUE),
    abs_edge_score = max(abs_edge_score, na.rm = TRUE),
    n_supporting_links = data.table::uniqueN(link_id),
    best_peak_ID = peak_ID[which.max(abs_edge_score)]
  ), by = .(tf, gene_norm)]
  nodes <- data.table::data.table(node_id = unique(c(edges$from, edges$to)))
  nodes[, `:=`(
    node_type = ifelse(node_id %in% unique(link_dt$tf), "TF", "Gene"),
    node_role = ifelse(node_id == seed_tf, "seed_tf", ifelse(node_id %in% target_genes, "selected_target", "supporting_tf")),
    is_seed_tf = node_id == seed_tf,
    is_target_gene = node_id %in% target_genes,
    perturbation_log2fc = NA_real_,
    strict_mean_expression = NA_real_
  )]
  list(nodes = nodes, edges = edges, targets = target_dt)
}


.module2_report_safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "item")
}

.module2_report_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

.module2_report_json <- function(x) {
  jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null", digits = 6)
}

.module2_report_write_html <- function(path, title, body, script) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  style <- paste(
    "body{margin:0;font-family:Arial,Helvetica,sans-serif;background:#f7f9fb;color:#18212f}",
    ".app{min-height:100vh;display:flex;flex-direction:column}",
    "header{background:#111827;color:white;padding:14px 18px;border-bottom:4px solid #2dd4bf}",
    "h1{font-size:20px;line-height:1.25;margin:0 0 6px 0;font-weight:700}.meta{font-size:12px;color:#cbd5e1}",
    ".toolbar{display:flex;flex-wrap:wrap;gap:10px;align-items:center;padding:10px 14px;background:white;border-bottom:1px solid #d8e0ea}",
    ".toolbar label{font-size:12px;font-weight:700;color:#344256;display:flex;gap:6px;align-items:center}",
    "select,input,button{font:inherit;font-size:12px;border:1px solid #b7c2d0;border-radius:6px;background:white;padding:5px 7px}",
    "button{cursor:pointer;background:#0f766e;color:white;border-color:#0f766e;font-weight:700}",
    ".main{display:grid;grid-template-columns:minmax(0,1fr) 320px;gap:0;min-height:0;flex:1}",
    ".canvasWrap{background:white;min-height:720px;position:relative;overflow:hidden}svg{width:100%;height:100%;min-height:720px;display:block}",
    ".side{border-left:1px solid #d8e0ea;background:#fbfdff;padding:12px;overflow:auto}.side h2{font-size:14px;margin:0 0 8px 0}",
    ".stat{font-size:12px;margin:5px 0;color:#344256}.node text{paint-order:stroke;stroke:white;stroke-width:4px;stroke-linejoin:round;font-weight:700}",
    ".edge{stroke:#64748b;stroke-opacity:.38}.edge.strong{stroke:#dc2626;stroke-opacity:.72}@media(max-width:900px){.main{grid-template-columns:1fr}.side{border-left:0;border-top:1px solid #d8e0ea}.canvasWrap,svg{min-height:560px}}",
    sep = "\n"
  )
  html <- paste0("<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>", .module2_report_html_escape(title), "</title><style>", style, "</style></head><body>", body, "<script>", script, "</script></body></html>")
  writeLines(html, path, useBytes = TRUE)
  path
}

.module2_report_write_network_html <- function(nodes, edges, out_html, title, subtitle = NULL) {
  is_seed_tf <- node_id <- node_type <- NULL
  nodes <- data.table::as.data.table(nodes)
  edges <- data.table::as.data.table(edges)
  seed_tf <- nodes[is_seed_tf %in% TRUE, node_id]
  if (!length(seed_tf)) seed_tf <- nodes[node_type == "TF", node_id]
  if (!length(seed_tf)) seed_tf <- title
  .module2_report_browser_write_network_html(
    nodes = nodes,
    edges = edges,
    case = list(seed_tf = seed_tf[[1L]]),
    out_html = out_html,
    title = title,
    tf_only = FALSE,
    max_edges_default = 0L
  )
}
.module2_report_build_direct_edges <- function(link_dt, min_supporting_peaks = 1L) {
  gene_norm <- from <- to <- fp_target_rna_r <- tf_expression_target_r <- NULL
  condition_fp_score <- n_supporting_peaks <- edge_score <- edge_rank <- NULL
  if (!nrow(link_dt) || !all(c("tf", "gene_norm", "peak_ID") %in% names(link_dt))) return(data.table::data.table())
  dt <- data.table::copy(link_dt)
  dt[, `:=`(from = toupper(as.character(tf)), to = toupper(as.character(gene_norm)), peak_ID = as.character(peak_ID))]
  dt <- dt[from != to]
  if (!nrow(dt)) return(data.table::data.table())
  dt[, link_score := suppressWarnings(as.numeric(fp_target_rna_r)) * suppressWarnings(as.numeric(tf_expression_target_r))]
  edge_dt <- dt[, .(n_supporting_peaks = data.table::uniqueN(peak_ID), edge_score = sum(link_score, na.rm = TRUE), mean_fp_target_rna_r = mean(fp_target_rna_r, na.rm = TRUE), mean_tf_expression_target_r = mean(tf_expression_target_r, na.rm = TRUE), mean_condition_fp_score = mean(condition_fp_score, na.rm = TRUE)), by = .(condition, from, to)]
  edge_dt <- edge_dt[n_supporting_peaks >= min_supporting_peaks & is.finite(edge_score) & edge_score > 0]
  data.table::setorder(edge_dt, condition, -edge_score, -n_supporting_peaks, from, to)
  edge_dt[, edge_rank := seq_len(.N), by = condition]
  edge_dt[]
}


.module2_report_condition_edges <- function(link_dt, multiomic_data, conditions = NULL, fp_score_cutoff = 0, verbose = TRUE) {
  gene_norm <- fp_target_rna_r <- tf_expression_target_r <- NULL
  if (!nrow(link_dt)) return(data.table::data.table())
  if (is.null(multiomic_data)) {
    dt <- data.table::copy(link_dt)
    dt[, `:=`(condition = "all", condition_fp_score = NA_real_)]
    return(.module2_report_build_direct_edges(dt, min_supporting_peaks = 1L))
  }
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  validate_multiomic_object(multiomic_data)
  mats <- multiomic_data$matrices
  all_conditions <- colnames(mats$fp_score)
  if (is.null(conditions)) conditions <- all_conditions
  conditions <- intersect(as.character(conditions), all_conditions)
  if (!length(conditions)) .log_abort("No report conditions overlap multiomic_data condition names.")
  dt <- data.table::copy(link_dt)
  fp_idx <- match(dt$fp_id, rownames(mats$fp_score))
  tf_idx <- match(dt$tf, rownames(mats$gene_expr))
  gene_idx <- match(dt$target_gene, rownames(mats$gene_expr))
  rows <- vector("list", length(conditions))
  for (i in seq_along(conditions)) {
    cc <- conditions[[i]]
    fp_score <- as.numeric(mats$fp_score[fp_idx, cc])
    fp_bound <- as.logical(mats$fp_bound[fp_idx, cc])
    tf_expr <- as.numeric(mats$gene_expr[tf_idx, cc])
    gene_expr <- as.numeric(mats$gene_expr[gene_idx, cc])
    keep <- is.finite(fp_score) & fp_score >= fp_score_cutoff &
      fp_bound %in% TRUE &
      is.finite(tf_expr) & tf_expr > 0 &
      is.finite(gene_expr) & gene_expr > 0
    if (!any(keep)) next
    one <- dt[keep, .(
      tf,
      gene_norm,
      peak_ID,
      fp_target_rna_r,
      tf_expression_target_r,
      condition = cc,
      condition_fp_score = fp_score[keep]
    )]
    rows[[i]] <- .module2_report_build_direct_edges(one, min_supporting_peaks = 1L)
    if (isTRUE(verbose)) .log_inform("Module 2 reports: aggregated direct TF-TF edges for condition {i}/{length(conditions)}: {cc}.")
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.module2_report_write_edge_cache <- function(edge_dt, output_dir, tag, result_label) {
  csv_dir <- file.path(output_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  edge_path <- file.path(csv_dir, sprintf("%s_%s_conditionFiltered_direct_tf_tf_edges.csv", tag, result_label))
  data.table::fwrite(edge_dt, edge_path)
  edge_path
}

.module2_report_norm_tf <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out[is.na(out)] <- ""
  out
}

.module2_report_cluster_rank <- function(cluster) {
  rank <- suppressWarnings(as.integer(sub("^[A-Za-z]+0*", "", as.character(cluster))))
  rank[!is.finite(rank)] <- .Machine$integer.max
  rank
}

.module2_report_order_tfs_for_heatmap <- function(tfs, mat, cluster_dt, axis = c("row", "col")) {
  axis <- match.arg(axis)
  cluster_rank <- NULL
  if (!length(tfs)) return(tfs)
  cl <- data.table::copy(cluster_dt)
  cl[, tf := .module2_report_norm_tf(tf)]
  cluster_map <- stats::setNames(cl$cluster, cl$tf)
  cluster <- unname(cluster_map[.module2_report_norm_tf(tfs)])
  cluster[is.na(cluster) | !nzchar(cluster)] <- "T00"
  strength <- if (identical(axis, "row")) rowSums(mat[tfs, , drop = FALSE], na.rm = TRUE) else colSums(mat[, tfs, drop = FALSE], na.rm = TRUE)
  ord <- data.table::data.table(tf = tfs, cluster = cluster, cluster_rank = .module2_report_cluster_rank(cluster), strength = as.numeric(strength))
  ord[order(cluster_rank, cluster, -strength, tf), tf]
}

.module2_report_draw_browser_heatmap_png <- function(mat, cluster_dt, out_png, condition, png_size = 900L, png_aspect = 0.60) {
  row_cluster <- col_cluster <- row_index <- col_index <- row_start <- row_end <- col_start <- col_end <- NULL
  join_key <- direct_tf_cluster <- heatmap_png <- h <- NULL
  if (!nrow(mat) || !ncol(mat)) return(data.table::data.table())
  if (!capabilities("png")) {
    .log_warn("PNG device is unavailable; skipping Module 2 browser heatmap.")
    return(data.table::data.table())
  }
  png_size <- suppressWarnings(as.integer(png_size[[1L]]))
  if (!is.finite(png_size) || is.na(png_size) || png_size < 300L) png_size <- 900L
  png_aspect <- suppressWarnings(as.numeric(png_aspect[[1L]]))
  if (!is.finite(png_aspect) || is.na(png_aspect) || png_aspect <= 0) png_aspect <- 0.60
  png_width <- max(300L, as.integer(round(png_size * png_aspect)))
  cl <- unique(data.table::copy(cluster_dt)[, .(tf = .module2_report_norm_tf(tf), cluster = as.character(cluster))])
  cluster_map <- stats::setNames(cl$cluster, cl$tf)
  row_cluster <- sub("^T", "R", unname(cluster_map[.module2_report_norm_tf(rownames(mat))]))
  col_cluster <- sub("^T", "C", unname(cluster_map[.module2_report_norm_tf(colnames(mat))]))
  row_cluster[is.na(row_cluster) | !nzchar(row_cluster)] <- "R00"
  col_cluster[is.na(col_cluster) | !nzchar(col_cluster)] <- "C00"
  finite_vals <- as.numeric(mat[is.finite(mat)])
  if (!length(finite_vals)) finite_vals <- 0
  fill_limits <- range(finite_vals, na.rm = TRUE)
  if (!is.finite(fill_limits[[1L]]) || !is.finite(fill_limits[[2L]]) || fill_limits[[1L]] == fill_limits[[2L]]) fill_limits <- c(0, max(1, fill_limits[[2L]]))
  pal <- grDevices::colorRampPalette(c("#2166AC", "#FEE08B", "#B2182B"))(256L)
  scaled <- (mat - fill_limits[[1L]]) / diff(fill_limits)
  scaled[!is.finite(scaled)] <- 0
  scaled <- pmin(1, pmax(0, scaled))
  color_mat <- matrix(pal[pmax(1L, pmin(256L, floor(scaled * 255) + 1L))], nrow = nrow(mat), dimnames = dimnames(mat))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  opened <- tryCatch({ grDevices::png(out_png, width = png_width, height = png_size, bg = "white", type = "cairo-png"); TRUE }, error = function(e) FALSE)
  if (!opened) grDevices::png(out_png, width = png_width, height = png_size, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::rasterImage(grDevices::as.raster(color_mat), 0, 0, 1, 1, interpolate = FALSE)
  row_bounds <- data.table::data.table(row_cluster = row_cluster, row_index = seq_along(row_cluster))[, .(row_start = min(row_index), row_end = max(row_index)), by = row_cluster]
  col_bounds <- data.table::data.table(col_cluster = col_cluster, col_index = seq_along(col_cluster))[, .(col_start = min(col_index), col_end = max(col_index)), by = col_cluster]
  row_bounds <- row_bounds[row_cluster != "R00"]
  col_bounds <- col_bounds[col_cluster != "C00"]
  if (!nrow(row_bounds) || !nrow(col_bounds)) return(data.table::data.table())
  col_bounds[, join_key := 1L]
  row_bounds[, join_key := 1L]
  out <- merge(col_bounds, row_bounds, by = "join_key", allow.cartesian = TRUE)
  out[, join_key := NULL]
  out[, `:=`(condition = as.character(condition), direct_tf_cluster = paste(col_cluster, row_cluster, sep = "-"), heatmap_png = out_png, x = (col_start - 1) / ncol(mat), y = (row_start - 1) / nrow(mat), w = (col_end - col_start + 1) / ncol(mat), h = (row_end - row_start + 1) / nrow(mat))]
  out[, .(condition, direct_tf_cluster, heatmap_png, x, y, w, h)]
}

.module2_report_heatmap_metadata <- function(edge_dt, k_label, report_name, png_size = 900L) {
  condition <- from <- to <- edge_score <- from_cluster <- to_cluster <- NULL
  if (!nrow(edge_dt)) return(data.table::data.table())
  dt <- data.table::copy(edge_dt)
  required <- c("condition", "from", "to", "edge_score", "from_cluster", "to_cluster")
  if (!all(required %in% names(dt))) return(data.table::data.table())
  tmp_dir <- tempfile("craftgrn-module2-heatmap-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(sort(unique(as.character(dt$condition))), function(cc) {
    one <- dt[condition == cc]
    if (!nrow(one)) return(NULL)
    tfs <- sort(unique(c(as.character(one$from), as.character(one$to))))
    mat <- matrix(0, nrow = length(tfs), ncol = length(tfs), dimnames = list(tfs, tfs))
    edge_sum <- one[, .(edge_score = sum(as.numeric(edge_score), na.rm = TRUE)), by = .(from, to)]
    mat[cbind(match(edge_sum$from, tfs), match(edge_sum$to, tfs))] <- edge_sum$edge_score
    if (!identical(report_name, "TF-TF connectivity")) mat <- log10(mat + 1)
    cl_from <- one[, .(tf = as.character(from), cluster = as.character(from_cluster))]
    cl_to <- one[, .(tf = as.character(to), cluster = as.character(to_cluster))]
    cl <- unique(data.table::rbindlist(list(cl_from, cl_to), use.names = TRUE, fill = TRUE))
    cl <- cl[nzchar(tf) & nzchar(cluster)]
    mat <- .module2_report_order_matrix_by_clusters(mat, cl)
    png_name <- sprintf("%s_%s_%s_heatmap.png", .module2_report_safe_filename(cc), .module2_report_safe_filename(report_name), k_label)
    .module2_report_draw_browser_heatmap_png(mat, cl, file.path(tmp_dir, png_name), cc, png_size = png_size)
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}


.module2_report_draw_connectivity_heatmap_png <- function(mat, cluster_dt, out_png, condition, png_size = 4800L, png_aspect = 0.60) {
  major_cluster <- row_index <- row_start <- row_end <- col_index <- col_start <- col_end <- h <- heatmap_png <- NULL
  if (!nrow(mat) || !ncol(mat)) return(data.table::data.table())
  if (!capabilities("png")) {
    .log_warn("PNG device is unavailable; skipping Module 2 connectivity browser heatmap.")
    return(data.table::data.table())
  }
  png_size <- suppressWarnings(as.integer(png_size[[1L]]))
  if (!is.finite(png_size) || is.na(png_size) || png_size < 300L) png_size <- 4800L
  png_aspect <- suppressWarnings(as.numeric(png_aspect[[1L]]))
  if (!is.finite(png_aspect) || is.na(png_aspect) || png_aspect <= 0) png_aspect <- 0.60
  png_width <- max(300L, as.integer(round(png_size * png_aspect)))
  cl <- unique(data.table::copy(cluster_dt)[, .(tf = .module2_report_norm_tf(tf), major_cluster = as.character(cluster))])
  cluster_map <- stats::setNames(cl$major_cluster, cl$tf)
  row_cluster <- unname(cluster_map[.module2_report_norm_tf(rownames(mat))])
  col_cluster <- unname(cluster_map[.module2_report_norm_tf(colnames(mat))])
  row_cluster[is.na(row_cluster) | !nzchar(row_cluster)] <- "T00"
  col_cluster[is.na(col_cluster) | !nzchar(col_cluster)] <- "T00"
  finite_vals <- as.numeric(mat[is.finite(mat)])
  if (!length(finite_vals)) finite_vals <- 0
  fill_limits <- range(finite_vals, na.rm = TRUE)
  if (!is.finite(fill_limits[[1L]]) || !is.finite(fill_limits[[2L]]) || fill_limits[[1L]] == fill_limits[[2L]]) fill_limits <- c(0, max(1, fill_limits[[2L]]))
  pal <- if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7L, name = "RdYlBu")))(100L)
  } else {
    grDevices::colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100L)
  }
  scaled <- (mat - fill_limits[[1L]]) / diff(fill_limits)
  scaled[!is.finite(scaled)] <- 0
  scaled <- pmin(1, pmax(0, scaled))
  color_mat <- matrix(pal[pmax(1L, pmin(length(pal), floor(scaled * (length(pal) - 1L)) + 1L))], nrow = nrow(mat), dimnames = dimnames(mat))
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  opened <- tryCatch({ grDevices::png(out_png, width = png_width, height = png_size, bg = "white", type = "cairo-png"); TRUE }, error = function(e) FALSE)
  if (!opened) grDevices::png(out_png, width = png_width, height = png_size, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::rasterImage(grDevices::as.raster(color_mat), 0, 0, 1, 1, interpolate = FALSE)
  row_bounds <- data.table::data.table(major_cluster = row_cluster, row_index = seq_along(row_cluster))[, .(row_start = min(row_index), row_end = max(row_index)), by = major_cluster]
  col_bounds <- data.table::data.table(major_cluster = col_cluster, col_index = seq_along(col_cluster))[, .(col_start = min(col_index), col_end = max(col_index)), by = major_cluster]
  out <- merge(col_bounds[major_cluster != "T00"], row_bounds[major_cluster != "T00"], by = "major_cluster", all = FALSE, sort = FALSE)
  if (!nrow(out)) return(data.table::data.table())
  out[, `:=`(condition = as.character(condition), heatmap_png = out_png, x = (col_start - 1) / ncol(mat), y = (row_start - 1) / nrow(mat), w = (col_end - col_start + 1) / ncol(mat), h = (row_end - row_start + 1) / nrow(mat))]
  out[, .(condition, major_cluster, heatmap_png, x, y, w, h)]
}

.module2_report_connectivity_heatmap_metadata <- function(edge_dt, k_label, png_size = 4800L) {
  condition <- edge_score <- from <- from_cluster <- major_cluster <- to <- to_cluster <- NULL
  if (!nrow(edge_dt)) return(data.table::data.table())
  dt <- data.table::copy(edge_dt)
  required <- c("condition", "from", "to", "edge_score", "from_cluster", "to_cluster")
  if (!all(required %in% names(dt))) return(data.table::data.table())
  if (!"major_cluster" %in% names(dt)) dt[, major_cluster := from_cluster]
  tmp_dir <- tempfile("craftgrn-module2-connectivity-heatmap-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(sort(unique(as.character(dt$condition))), function(cc) {
    one <- dt[condition == cc]
    if (!nrow(one)) return(NULL)
    tfs <- sort(unique(c(as.character(one$from), as.character(one$to))))
    mat <- matrix(0, nrow = length(tfs), ncol = length(tfs), dimnames = list(tfs, tfs))
    edge_sum <- one[, .(edge_score = sum(as.numeric(edge_score), na.rm = TRUE)), by = .(from, to)]
    mat[cbind(match(edge_sum$from, tfs), match(edge_sum$to, tfs))] <- edge_sum$edge_score
    mat[cbind(match(edge_sum$to, tfs), match(edge_sum$from, tfs))] <- edge_sum$edge_score
    cl_from <- one[, .(tf = as.character(from), cluster = as.character(from_cluster))]
    cl_to <- one[, .(tf = as.character(to), cluster = as.character(to_cluster))]
    cl <- unique(data.table::rbindlist(list(cl_from, cl_to), use.names = TRUE, fill = TRUE))
    cl <- cl[nzchar(tf) & nzchar(cluster)]
    mat <- .module2_report_order_matrix_by_clusters(mat, cl)
    png_name <- sprintf("%s_tf_tf_connectivity_heatmap_%s_browser.png", .module2_report_safe_filename(cc), k_label)
    .module2_report_draw_connectivity_heatmap_png(mat, cl, file.path(tmp_dir, png_name), cc, png_size = png_size)
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.module2_report_connectivity_payload <- function(edge_dt) {
  condition <- database <- edge_r <- edge_score <- from <- from_cluster <- major_cluster <- mean_condition_fp_score <- n_supporting_links <- network_id <- network_label <- pathway_adjusted_p_value <- pathway_key <- pathway_neglog10_padj <- pathway_term <- to <- NULL
  if (!nrow(edge_dt)) return(data.table::data.table())
  dt <- data.table::copy(edge_dt)
  if (!"major_cluster" %in% names(dt)) dt[, major_cluster := from_cluster]
  dt[is.na(major_cluster) | !nzchar(major_cluster), major_cluster := "T00"]
  dt[, network_id := paste(condition, major_cluster, "composite_tf_tf", sep = "__")]
  dt[, network_label := paste0(major_cluster, " composite TF-TF connectivity")]
  dt[, database := "Module2"]
  dt[, pathway_key := "composite_tf_tf_connectivity"]
  dt[, pathway_term := paste0(major_cluster, " composite TF-TF connectivity")]
  dt[, pathway_adjusted_p_value := 1]
  dt[, pathway_neglog10_padj := 0]
  dt[, edge_r := as.numeric(edge_score)]
  if (!"n_supporting_links" %in% names(dt)) dt[, n_supporting_links := 1L]
  if (!"mean_condition_fp_score" %in% names(dt)) dt[, mean_condition_fp_score := NA_real_]
  dt[, .(condition, network_id, network_label, major_cluster, database, pathway_key, pathway_term, pathway_adjusted_p_value, pathway_neglog10_padj, from, to, edge_score, edge_r, n_supporting_links, mean_condition_fp_score)]
}

.module2_report_write_tf_tf_connectivity_browser <- function(edge_dt, out_html, title) {
  base <- tools::file_path_sans_ext(basename(out_html))
  k_label <- sub("_.*$", "", base)
  out_suffix <- sub("^[^_]+_", "", base)
  payload_dt <- .module2_report_connectivity_payload(edge_dt)
  .module2_report_browser_write_condition_pathway_tf_gene_k_index(
    edge_dt = payload_dt,
    out_dir = dirname(out_html),
    title_prefix = title,
    k_label = k_label,
    max_edges_per_network = 120L,
    out_suffix = out_suffix,
    heatmap_dt = .module2_report_connectivity_heatmap_metadata(edge_dt, k_label = k_label)
  )
}

.module2_report_export_edge_browsers <- function(edge_dt, output_dir, html_subdir, report, out_suffix, title_template, report_name, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), verbose = TRUE) {
  html_dir <- output_dir
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(as.integer(k_values), function(k) {
    k_edge_dt <- .module2_report_assign_edge_clusters(edge_dt, k = k)
    if (!nrow(k_edge_dt)) return(NULL)
    html <- file.path(html_dir, sprintf("K%02d_%s.html", as.integer(k), out_suffix))
    .module2_report_write_tf_tf_browser(
      edge_dt = k_edge_dt,
      out_html = html,
      title = sprintf(title_template, tag, result_label, as.integer(k)),
      report_name = report_name
    )
    data.table::data.table(report = report, k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} {report} HTML file(s).")
  tibble::as_tibble(out)
}

.module2_report_direct_matrix_from_edges <- function(edge_dt) {
  edge_score <- from <- to <- NULL
  if (!nrow(edge_dt)) return(matrix(numeric(), 0L, 0L))
  tfs <- sort(unique(c(as.character(edge_dt$from), as.character(edge_dt$to))))
  if (length(tfs) < 2L) return(matrix(numeric(), 0L, 0L))
  mat <- matrix(0, nrow = length(tfs), ncol = length(tfs), dimnames = list(tfs, tfs))
  edge_sum <- edge_dt[, .(edge_score = sum(as.numeric(edge_score), na.rm = TRUE)), by = .(from, to)]
  mat[cbind(match(edge_sum$from, tfs), match(edge_sum$to, tfs))] <- edge_sum$edge_score
  log10(mat + 1)
}

.module2_report_cut_condition_clusters <- function(mat, k) {
  if (!nrow(mat)) return(data.table::data.table(tf = character(), cluster = character()))
  k <- as.integer(k[[1L]])
  if (!is.finite(k) || is.na(k) || k < 1L) k <- 1L
  k <- min(k, nrow(mat))
  labels <- if (nrow(mat) >= 2L && k > 1L) {
    hc <- tryCatch(stats::hclust(stats::dist(mat), method = "complete"), error = function(e) NULL)
    if (is.null(hc)) rep("T01", nrow(mat)) else paste0("T", sprintf("%02d", stats::cutree(hc, k = k)))
  } else {
    rep("T01", nrow(mat))
  }
  data.table::data.table(tf = rownames(mat), cluster = labels)
}

.module2_report_order_names_by_cluster <- function(names_vec, mat, cluster_dt, axis = c("row", "col")) {
  cluster_rank <- NULL
  axis <- match.arg(axis)
  names_vec <- as.character(names_vec)
  if (!length(names_vec) || !nrow(cluster_dt)) return(names_vec)
  cl <- unique(data.table::copy(cluster_dt)[, .(tf = .module2_report_norm_tf(tf), cluster = as.character(cluster))])
  cluster_map <- stats::setNames(cl$cluster, cl$tf)
  cluster_vec <- unname(cluster_map[.module2_report_norm_tf(names_vec)])
  cluster_vec[is.na(cluster_vec) | !nzchar(cluster_vec)] <- "T00"
  order_dt <- data.table::data.table(tf = names_vec, cluster = cluster_vec)
  order_dt[, cluster_rank := .module2_report_cluster_rank(cluster)]
  cluster_levels <- unique(order_dt[order(cluster_rank, cluster), cluster])
  unlist(lapply(cluster_levels, function(cluster_id) {
    members <- order_dt[cluster == cluster_id, tf]
    if (length(members) < 2L) return(members)
    idx <- match(members, if (identical(axis, "row")) rownames(mat) else colnames(mat))
    sub_mat <- if (identical(axis, "row")) mat[idx, , drop = FALSE] else t(mat[, idx, drop = FALSE])
    ord <- tryCatch(stats::hclust(stats::dist(sub_mat), method = "complete")$order, error = function(e) seq_along(members))
    members[ord]
  }), use.names = FALSE)
}

.module2_report_order_matrix_by_clusters <- function(mat, cluster_dt) {
  if (!nrow(mat) || !ncol(mat) || !nrow(cluster_dt)) return(mat)
  row_order <- .module2_report_order_names_by_cluster(rownames(mat), mat, cluster_dt, axis = "row")
  col_order <- .module2_report_order_names_by_cluster(colnames(mat), mat, cluster_dt, axis = "col")
  mat[row_order, col_order, drop = FALSE]
}

.module2_report_assign_edge_clusters <- function(edge_dt, k) {
  condition <- direct_tf_cluster <- to_cluster <- from_cluster <- NULL
  if (!nrow(edge_dt)) return(edge_dt)
  dt <- data.table::copy(edge_dt)
  if (!"condition" %in% names(dt)) dt[, condition := "all"]
  rows <- lapply(sort(unique(as.character(dt$condition))), function(cc) {
    one <- dt[condition == cc]
    mat <- .module2_report_direct_matrix_from_edges(one)
    cl <- .module2_report_cut_condition_clusters(mat, k = k)
    if (!nrow(cl)) return(data.table::data.table())
    out <- merge(one, cl, by.x = "from", by.y = "tf", all.x = TRUE, sort = FALSE)
    data.table::setnames(out, "cluster", "from_cluster")
    out <- merge(out, cl, by.x = "to", by.y = "tf", all.x = TRUE, sort = FALSE)
    data.table::setnames(out, "cluster", "to_cluster")
    out[is.na(from_cluster) | !nzchar(from_cluster), from_cluster := "T00"]
    out[is.na(to_cluster) | !nzchar(to_cluster), to_cluster := "T00"]
    out[, direct_tf_cluster := paste(sub("^T", "C", to_cluster), sub("^T", "R", from_cluster), sep = "-")]
    out[]
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.module2_report_composite_score_matrix_from_links <- function(link_dt) {
  gene_norm <- target_norm <- tf_norm <- NULL
  if (!nrow(link_dt) || !all(c("tf", "gene_norm") %in% names(link_dt))) return(matrix(numeric(), 0L, 0L))
  dt <- unique(data.table::copy(link_dt)[, .(tf_norm = .module2_report_norm_tf(tf), target_norm = .module2_report_norm_tf(gene_norm))])
  dt <- dt[nzchar(tf_norm) & nzchar(target_norm)]
  if (!nrow(dt)) return(matrix(numeric(), 0L, 0L))
  tf_levels <- sort(unique(dt$tf_norm))
  if (length(tf_levels) < 2L) return(matrix(numeric(), 0L, 0L))
  direct <- matrix(0, nrow = length(tf_levels), ncol = length(tf_levels), dimnames = list(tf_levels, tf_levels))
  tf2tf <- dt[target_norm %in% tf_levels]
  if (nrow(tf2tf)) direct[cbind(match(tf2tf$tf_norm, tf_levels), match(tf2tf$target_norm, tf_levels))] <- 1
  target_levels <- sort(unique(dt$target_norm))
  incidence <- matrix(0, nrow = length(tf_levels), ncol = length(target_levels), dimnames = list(tf_levels, target_levels))
  incidence[cbind(match(dt$tf_norm, tf_levels), match(dt$target_norm, target_levels))] <- 1
  comp <- direct + t(direct) + 0.5 * tcrossprod(incidence)
  diag(comp) <- 0
  comp
}

.module2_report_filter_composite_matrix <- function(score_mat, min_score = 2, max_tfs = Inf) {
  if (!nrow(score_mat) || !ncol(score_mat)) return(matrix(numeric(), 0L, 0L))
  mat <- score_mat
  mat[mat < min_score] <- 0
  keep <- pmax(rowSums(mat > 0, na.rm = TRUE), colSums(mat > 0, na.rm = TRUE)) > 0
  mat <- mat[keep, keep, drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(matrix(numeric(), 0L, 0L))
  if (is.finite(max_tfs) && !is.na(max_tfs) && max_tfs > 1L && nrow(mat) > max_tfs) {
    tf_strength <- rowSums(mat, na.rm = TRUE) + colSums(mat, na.rm = TRUE)
    keep_tf <- names(sort(tf_strength, decreasing = TRUE))[seq_len(max_tfs)]
    mat <- mat[keep_tf, keep_tf, drop = FALSE]
  }
  log2(mat + 1)
}

.module2_report_composite_edges_for_k <- function(link_dt, k, min_score = 2, max_tfs = Inf) {
  condition <- edge_rank <- edge_score <- from <- from_cluster <- major_cluster <- to <- to_cluster <- direct_tf_cluster <- NULL
  if (!nrow(link_dt)) return(data.table::data.table())
  dt <- data.table::copy(link_dt)
  if (!"condition" %in% names(dt)) dt[, condition := "all"]
  rows <- lapply(sort(unique(as.character(dt$condition))), function(cc) {
    score_mat <- .module2_report_composite_score_matrix_from_links(dt[condition == cc])
    mat <- .module2_report_filter_composite_matrix(score_mat, min_score = min_score, max_tfs = max_tfs)
    if (!nrow(mat)) return(data.table::data.table())
    cl <- .module2_report_cut_condition_clusters(mat, k = k)
    mat <- .module2_report_order_matrix_by_clusters(mat, cl)
    idx <- which(upper.tri(mat) & is.finite(mat) & mat > 0, arr.ind = TRUE)
    if (!nrow(idx)) return(data.table::data.table())
    edge_dt <- data.table::data.table(
      condition = cc,
      from = rownames(mat)[idx[, 1L]],
      to = colnames(mat)[idx[, 2L]],
      edge_score = as.numeric(mat[idx]),
      n_supporting_peaks = NA_integer_,
      mean_fp_target_rna_r = NA_real_,
      mean_tf_expression_target_r = NA_real_,
      mean_condition_fp_score = NA_real_
    )
    edge_dt <- merge(edge_dt, cl, by.x = "from", by.y = "tf", all.x = TRUE, sort = FALSE)
    data.table::setnames(edge_dt, "cluster", "from_cluster")
    edge_dt <- merge(edge_dt, cl, by.x = "to", by.y = "tf", all.x = TRUE, sort = FALSE)
    data.table::setnames(edge_dt, "cluster", "to_cluster")
    edge_dt[is.na(from_cluster) | !nzchar(from_cluster), from_cluster := "T00"]
    edge_dt[is.na(to_cluster) | !nzchar(to_cluster), to_cluster := "T00"]
    edge_dt[, direct_tf_cluster := paste(sub("^T", "C", to_cluster), sub("^T", "R", from_cluster), sep = "-")]
    edge_dt <- edge_dt[from_cluster == to_cluster]
    if (!nrow(edge_dt)) return(data.table::data.table())
    edge_dt[, major_cluster := from_cluster]
    data.table::setorder(edge_dt, condition, from_cluster, to_cluster, -edge_score, from, to)
    edge_dt[, edge_rank := seq_len(.N), by = condition]
    edge_dt[]
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.module2_report_write_tf_tf_browser <- function(edge_dt, out_html, title, report_name) {
  base <- tools::file_path_sans_ext(basename(out_html))
  k_label <- sub("_.*$", "", base)
  out_suffix <- sub("^[^_]+_", "", base)
  page_label <- if (identical(report_name, "TF-TF connectivity")) "TF-TF connectivity" else "direct TF-TF"
  edge_label <- if (identical(report_name, "TF-TF connectivity")) "condition-filtered TF-to-TF connectivity edges" else "condition-filtered direct TF-to-TF edges"
  .module2_report_browser_write_condition_direct_tf_tf_k_index(
    edge_dt = data.table::as.data.table(edge_dt),
    out_dir = dirname(out_html),
    title_prefix = title,
    k_label = k_label,
    max_edges_per_cluster = 250L,
    max_nodes_per_cluster = Inf,
    heatmap_dt = .module2_report_heatmap_metadata(edge_dt, k_label = k_label, report_name = report_name),
    network_label = report_name,
    page_label = page_label,
    edge_label = edge_label,
    out_suffix = out_suffix
  )
}
.module2_report_export_composite_browsers <- function(active_dt, fallback_edge_dt, output_dir, report, out_suffix, title_template, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), min_score = 2, verbose = TRUE) {
  html_dir <- output_dir
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(as.integer(k_values), function(k) {
    k_edge_dt <- .module2_report_composite_edges_for_k(active_dt, k = k, min_score = min_score)
    if (!nrow(k_edge_dt)) k_edge_dt <- .module2_report_assign_edge_clusters(fallback_edge_dt, k = k)
    if (!nrow(k_edge_dt)) return(NULL)
    html <- file.path(html_dir, sprintf("K%02d_%s.html", as.integer(k), out_suffix))
    .module2_report_write_tf_tf_connectivity_browser(
      edge_dt = k_edge_dt,
      out_html = html,
      title = sprintf(title_template, tag, result_label, as.integer(k))
    )
    data.table::data.table(report = report, k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} TF-TF connectivity HTML file(s).")
  tibble::as_tibble(out)
}

#' Build optional Module 2 HTML reports
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param multiomic_data Optional CraftGRN multiomic object for condition-filtered reports.
#' @param output_dir Report output directory.
#' @param reports Report families to build.
#' @param tfs TFs for top target reports.
#' @param conditions Optional condition subset.
#' @param k_values Cluster counts for TF-TF browsers.
#' @param top_n Number of top targets per TF.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @noRd
build_module2_reports <- function(module2, multiomic_data = NULL, output_dir = NULL, reports = c("top_tf_targets", "direct_tf_tf", "tf_tf_connectivity"), tfs = NULL, conditions = NULL, k_values = c(5L, 7L, 10L), top_n = 100L, verbose = TRUE) {
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  reports <- match.arg(reports, several.ok = TRUE)
  if (is.null(output_dir)) output_dir <- file.path(getwd(), "module2_reports")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out <- list()
  top_tfs <- tfs
  if ("top_tf_targets" %in% reports) {
    if (is.null(top_tfs)) {
      tf_corr <- .module2_report_read_table(module2, "module2_tf_target_corr", columns = c("tf"))
      top_tfs <- head(sort(unique(toupper(as.character(tf_corr$tf)))), 10L)
    }
    out[[length(out) + 1L]] <- export_top_tf_targets(
      module2 = module2,
      output_dir = output_dir,
      tfs = top_tfs,
      top_n = top_n,
      verbose = verbose
    )
  }
  if (any(c("direct_tf_tf", "tf_tf_connectivity") %in% reports)) {
    if (isTRUE(verbose)) .log_inform("Module 2 reports: building direct TF-TF report links.")
    direct_links <- .module2_report_direct_links(module2)
    if (!nrow(direct_links)) .log_abort("No direct TF-TF Module 2 links were found for report generation.")
    if (isTRUE(verbose)) .log_inform("Module 2 reports: aggregating condition-filtered TF-TF edges from {nrow(direct_links)} direct links.")
    tag <- "module2"
    result_label <- "top100"
    active_dt <- .module2_report_condition_links(direct_links, multiomic_data, conditions = conditions)
    if (!nrow(active_dt)) .log_abort("No condition-filtered Module 2 links were found for report generation.")
    edge_dt <- .module2_report_build_direct_edges(active_dt, min_supporting_peaks = 1L)
    if (!nrow(edge_dt)) .log_abort("No condition-filtered direct TF-TF edges were found for report generation.")
    .module2_report_write_edge_cache(edge_dt, output_dir, tag, result_label)
    if ("direct_tf_tf" %in% reports) {
      out[[length(out) + 1L]] <- .module2_report_export_edge_browsers(
        edge_dt = edge_dt,
        output_dir = output_dir,
        html_subdir = "direct_tf_tf_networks",
        report = "direct_tf_tf",
        out_suffix = "direct_tf_tf_network_browser",
        title_template = "%s %s K%02d direct TF-TF network",
        report_name = "Direct TF-TF",
        tag = tag,
        result_label = result_label,
        k_values = k_values,
        verbose = verbose
      )
    }
    if ("tf_tf_connectivity" %in% reports) {
      out[[length(out) + 1L]] <- .module2_report_export_composite_browsers(
        active_dt = active_dt,
        fallback_edge_dt = edge_dt,
        output_dir = output_dir,
        report = "tf_tf_connectivity",
        out_suffix = "tf_tf_connectivity_network_browser",
        title_template = "%s %s K%02d TF-TF connectivity",
        tag = tag,
        result_label = result_label,
        k_values = k_values,
        min_score = 2,
        verbose = verbose
      )
    }
  }
  manifest <- tibble::as_tibble(data.table::rbindlist(out, use.names = TRUE, fill = TRUE))
  readr::write_csv(manifest, file.path(output_dir, "module2_report_manifest.csv"))
  manifest
}

#' Export per-TF top target HTML reports
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param output_dir Output directory.
#' @param tfs TFs to report.
#' @param top_n Number of top targets per TF.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @noRd
export_top_tf_targets <- function(module2, output_dir, tfs, top_n = 100L, verbose = TRUE) {
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  dir.create(file.path(output_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
  html_dir <- output_dir
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  link_dt <- .module2_report_join_links(module2, tf = tfs)
  tfs <- toupper(as.character(tfs))
  rows <- lapply(tfs, function(tf) {
    tables <- .module2_report_top_tables(link_dt, tf, top_n = top_n)
    if (is.null(tables)) return(NULL)
    html <- file.path(html_dir, sprintf("%s_Module2_top%d_log2fc1_network.html", tf, as.integer(top_n)))
    .module2_report_write_network_html(
      nodes = tables$nodes,
      edges = tables$edges,
      out_html = html,
      title = sprintf("%s Module2 top-%d target network", tf, as.integer(top_n)),
      subtitle = "Top target genes and supporting TF-to-target edges"
    )
    data.table::fwrite(tables$nodes, file.path(output_dir, "csv", sprintf("%s_Module2_top%d_nodes.csv", tf, as.integer(top_n))))
    data.table::fwrite(tables$edges, file.path(output_dir, "csv", sprintf("%s_Module2_top%d_edges.csv", tf, as.integer(top_n))))
    data.table::data.table(report = "top_tf_targets", tf = tf, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} top TF-target HTML file(s).")
  tibble::as_tibble(out)
}

#' Export an interactive HTML browser of individual TF regulons
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param output_dir Output directory.
#' @param tfs TFs to report.
#' @param top_n Number of top targets per TF.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @export
report_top_tf_targets <- function(module2, output_dir, tfs, top_n = 100L, verbose = TRUE) {
  export_top_tf_targets(module2 = module2, output_dir = output_dir, tfs = tfs, top_n = top_n, verbose = verbose)
}

#' Export direct TF-TF browser reports
#'
#' @param cache Internal condition-filtered cache from build_module2_reports.
#' @param output_dir Output directory.
#' @param tag Report tag.
#' @param result_label Result label.
#' @param k_values Cluster counts.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @noRd
export_direct_tf_tf_browser <- function(cache, output_dir, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), verbose = TRUE) {
  html_dir <- output_dir
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  link_dt <- data.table::rbindlist(lapply(cache$manifest$filtered_links_csv, data.table::fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
  base_edges <- .module2_report_build_direct_edges(link_dt, min_supporting_peaks = 1L)
  rows <- lapply(as.integer(k_values), function(k) {
    edge_dt <- .module2_report_assign_edge_clusters(base_edges, k = k)
    if (!nrow(edge_dt)) return(NULL)
    html <- file.path(html_dir, sprintf("K%02d_direct_tf_tf_network_browser.html", as.integer(k)))
    .module2_report_write_tf_tf_browser(
      edge_dt = edge_dt,
      out_html = html,
      title = sprintf("%s %s K%02d direct TF-TF network", tag, result_label, as.integer(k)),
      report_name = "Direct TF-TF"
    )
    data.table::data.table(report = "direct_tf_tf", k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} direct TF-TF HTML file(s).")
  tibble::as_tibble(out)
}

#' Export an interactive HTML browser of direct TF-TF regulations
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param output_dir Output directory.
#' @param multiomic_data Optional CraftGRN multiomic object for condition-filtered reports.
#' @param k_values Cluster counts.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @export
report_direct_tf_tf_regulations <- function(module2, output_dir, multiomic_data = NULL, k_values = c(5L, 7L, 10L), verbose = TRUE) {
  build_module2_reports(
    module2 = module2,
    multiomic_data = multiomic_data,
    output_dir = output_dir,
    reports = "direct_tf_tf",
    k_values = k_values,
    verbose = verbose
  )
}

#' Export TF-TF connectivity browser reports
#'
#' @param cache Internal condition-filtered cache from build_module2_reports.
#' @param output_dir Output directory.
#' @param tag Report tag.
#' @param result_label Result label.
#' @param k_values Cluster counts.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @noRd
export_tf_tf_connectivity_browser <- function(cache, output_dir, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), verbose = TRUE) {
  html_dir <- output_dir
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  link_dt <- data.table::rbindlist(lapply(cache$manifest$filtered_links_csv, data.table::fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
  base_edges <- .module2_report_build_direct_edges(link_dt, min_supporting_peaks = 1L)
  rows <- lapply(as.integer(k_values), function(k) {
    edge_dt <- .module2_report_composite_edges_for_k(link_dt, k = k, min_score = 2)
    if (!nrow(edge_dt)) edge_dt <- .module2_report_assign_edge_clusters(base_edges, k = k)
    if (!nrow(edge_dt)) return(NULL)
    html <- file.path(html_dir, sprintf("K%02d_tf_tf_connectivity_network_browser.html", as.integer(k)))
    .module2_report_write_tf_tf_connectivity_browser(
      edge_dt = edge_dt,
      out_html = html,
      title = sprintf("%s %s K%02d TF-TF connectivity", tag, result_label, as.integer(k))
    )
    data.table::data.table(report = "tf_tf_connectivity", k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} TF-TF connectivity HTML file(s).")
  tibble::as_tibble(out)
}

#' Export an interactive HTML browser of TF-TF co-regulatory activities
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param output_dir Output directory.
#' @param multiomic_data Optional CraftGRN multiomic object for condition-filtered reports.
#' @param k_values Cluster counts.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @export
report_tf_tf_coregulations <- function(module2, output_dir, multiomic_data = NULL, k_values = c(5L, 7L, 10L), verbose = TRUE) {
  build_module2_reports(
    module2 = module2,
    multiomic_data = multiomic_data,
    output_dir = output_dir,
    reports = "tf_tf_connectivity",
    k_values = k_values,
    verbose = verbose
  )
}
