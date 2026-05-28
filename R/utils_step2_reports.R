# File: utils_step2_reports.R
# Purpose: Optional Module 2 HTML reports from compact relational outputs.

.module2_report_source_legacy <- function(script_name, env = parent.frame()) {
  roots <- unique(c(getwd(), dirname(getwd()), dirname(dirname(getwd()))))
  candidates <- file.path(roots, "dev", "benchmark", script_name)
  script_path <- candidates[file.exists(candidates)][1L]
  if (is.na(script_path)) {
    .log_abort("Legacy report writer not found. Expected dev/benchmark/{script_name}.")
  }
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(dirname(dirname(dirname(script_path))))
  sys.source(script_path, envir = env, keep.source = FALSE)
  invisible(script_path)
}

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

.module2_report_condition_links <- function(link_dt, multiomic_data, conditions = NULL, fp_score_cutoff = 0) {
  condition_fp_score <- NULL
  if (!nrow(link_dt)) return(data.table::data.table())
  if (is.null(multiomic_data)) {
    out <- data.table::copy(link_dt)
    out[, condition := "all"]
    out[, condition_fp_score := NA_real_]
    return(out)
  }
  if (!is_multiomic_object(multiomic_data)) multiomic_data <- as_multiomic_object(multiomic_data, verbose = FALSE)
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

#' Build optional Module 2 HTML reports
#'
#' @param module2 Module 2 result list, loaded output list, or output directory.
#' @param multiomic_data Optional compact multiomic object for condition-filtered reports.
#' @param output_dir Report output directory.
#' @param reports Report families to build.
#' @param tfs TFs for top target reports.
#' @param conditions Optional condition subset.
#' @param k_values Cluster counts for TF-TF browsers.
#' @param top_n Number of top targets per TF.
#' @param verbose Emit concise progress messages.
#' @return A tibble report manifest.
#' @export
build_module2_reports <- function(module2, multiomic_data = NULL, output_dir = NULL, reports = c("top_tf_targets", "direct_tf_tf", "tf_tf_connectivity"), tfs = NULL, conditions = NULL, k_values = c(5L, 7L, 10L), top_n = 100L, verbose = TRUE) {
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  reports <- match.arg(reports, several.ok = TRUE)
  if (is.null(output_dir)) output_dir <- file.path(getwd(), "module2_reports")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  joined <- .module2_report_join_links(module2, tf = tfs)
  if (!nrow(joined)) .log_abort("No passing Module 2 links were found for report generation.")
  active <- .module2_report_condition_links(joined, multiomic_data, conditions = conditions)
  tag <- "module2"
  result_label <- "top100"
  out <- list()
  if ("top_tf_targets" %in% reports) {
    out[[length(out) + 1L]] <- export_top_tf_targets(
      module2 = module2,
      output_dir = file.path(output_dir, "top_tf_targets"),
      tfs = if (is.null(tfs)) head(sort(unique(joined$tf)), 10L) else tfs,
      top_n = top_n,
      verbose = verbose
    )
  }
  if (any(c("direct_tf_tf", "tf_tf_connectivity") %in% reports)) {
    cache <- .module2_report_write_active_cache(active, output_dir, tag, result_label)
    if ("direct_tf_tf" %in% reports) {
      out[[length(out) + 1L]] <- export_direct_tf_tf_browser(
        cache = cache,
        output_dir = output_dir,
        tag = tag,
        result_label = result_label,
        k_values = k_values,
        verbose = verbose
      )
    }
    if ("tf_tf_connectivity" %in% reports) {
      out[[length(out) + 1L]] <- export_tf_tf_connectivity_browser(
        cache = cache,
        output_dir = output_dir,
        tag = tag,
        result_label = result_label,
        k_values = k_values,
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
#' @export
export_top_tf_targets <- function(module2, output_dir, tfs, top_n = 100L, verbose = TRUE) {
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  legacy_env <- new.env(parent = parent.frame())
  .module2_report_source_legacy("02_build_tf_target_perturbation_html_networks.R", env = legacy_env)
  dir.create(file.path(output_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
  link_dt <- .module2_report_join_links(module2, tf = tfs)
  tfs <- toupper(as.character(tfs))
  rows <- lapply(tfs, function(tf) {
    tables <- .module2_report_top_tables(link_dt, tf, top_n = top_n)
    if (is.null(tables)) return(NULL)
    case <- list(seed_tf = tf, cell_line = "Module2", output_prefix = paste0(tf, "_Module2_top", top_n))
    html <- file.path(output_dir, "pdf", sprintf("%s_Module2_top%d_log2fc1_network.html", tf, as.integer(top_n)))
    legacy_env$.write_network_html(
      nodes = tables$nodes,
      edges = tables$edges,
      case = case,
      out_html = html,
      title = sprintf("%s Module2 top-%d target network", tf, as.integer(top_n)),
      tf_only = FALSE,
      max_edges_default = 0L
    )
    data.table::fwrite(tables$nodes, file.path(output_dir, "csv", sprintf("%s_Module2_top%d_nodes.csv", tf, as.integer(top_n))))
    data.table::fwrite(tables$edges, file.path(output_dir, "csv", sprintf("%s_Module2_top%d_edges.csv", tf, as.integer(top_n))))
    data.table::data.table(report = "top_tf_targets", tf = tf, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} top TF-target HTML file(s).")
  tibble::as_tibble(out)
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
#' @export
export_direct_tf_tf_browser <- function(cache, output_dir, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), verbose = TRUE) {
  legacy_env <- new.env(parent = parent.frame())
  .module2_report_source_legacy("02_build_tf_regulated_cluster_html_browsers.R", env = legacy_env)
  pdf_dir <- file.path(output_dir, "pdf", "condition_filtered_direct_tf_tf")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  link_dt <- data.table::rbindlist(lapply(cache$manifest$filtered_links_csv, data.table::fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
  rows <- lapply(as.integer(k_values), function(k) {
    edge_dt <- legacy_env$.build_condition_direct_tf_tf_base_edges(link_dt, min_edge_score = 1L)
    if (!nrow(edge_dt)) return(NULL)
    edge_dt[, condition := if ("condition" %in% names(link_dt)) link_dt$condition[[1L]] else cache$manifest$condition[[1L]]]
    edge_dt[, `:=`(from_cluster = sprintf("T%02d", ((seq_len(.N) - 1L) %% k) + 1L), to_cluster = sprintf("T%02d", ((seq_len(.N) - 1L) %% k) + 1L))]
    html <- legacy_env$.write_condition_direct_tf_tf_k_index(
      edge_dt = edge_dt,
      out_dir = pdf_dir,
      title_prefix = tag,
      k_label = sprintf("K%02d", k),
      network_label = "Direct TF network",
      page_label = "direct TF-TF",
      edge_label = "direct TF-to-TF edges",
      out_suffix = "direct_tf_tf_network_browser"
    )
    data.table::data.table(report = "direct_tf_tf", k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} direct TF-TF HTML file(s).")
  tibble::as_tibble(out)
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
#' @export
export_tf_tf_connectivity_browser <- function(cache, output_dir, tag = "module2", result_label = "top100", k_values = c(5L, 7L, 10L), verbose = TRUE) {
  legacy_env <- new.env(parent = parent.frame())
  .module2_report_source_legacy("02_build_tf_regulated_cluster_html_browsers.R", env = legacy_env)
  pdf_dir <- file.path(output_dir, "pdf", "condition_filtered_tf_tf_connectivity")
  dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(as.integer(k_values), function(k) {
    link_dt <- data.table::rbindlist(lapply(cache$manifest$filtered_links_csv, data.table::fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
    edge_dt <- legacy_env$.build_condition_direct_tf_tf_base_edges(link_dt, min_edge_score = 1L)
    if (!nrow(edge_dt)) return(NULL)
    edge_dt[, condition := cache$manifest$condition[[1L]]]
    edge_dt[, `:=`(from_cluster = sprintf("T%02d", ((seq_len(.N) - 1L) %% k) + 1L), to_cluster = sprintf("T%02d", ((seq_len(.N) - 1L) %% k) + 1L))]
    html <- legacy_env$.write_condition_direct_tf_tf_k_index(
      edge_dt = edge_dt,
      out_dir = pdf_dir,
      title_prefix = tag,
      k_label = sprintf("K%02d", k),
      network_label = "TF-TF connectivity",
      page_label = "TF-TF connectivity",
      edge_label = "TF-TF connectivity edges",
      out_suffix = "tf_tf_connectivity_network_browser"
    )
    data.table::data.table(report = "tf_tf_connectivity", k = k, path = html)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (isTRUE(verbose)) .log_inform("Module 2 reports: wrote {nrow(out)} TF-TF connectivity HTML file(s).")
  tibble::as_tibble(out)
}
