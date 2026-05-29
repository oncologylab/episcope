# File: utils_module_qc_reports.R
# Author: Yaoxiang Li
# Purpose: Package-native HTML QC reports for Module 1 and Module 2.

.qc_html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

.qc_format_number <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), format(round(x, 3), big.mark = ",", scientific = FALSE, trim = TRUE), "NA")
}

.qc_percent <- function(num, den) {
  num <- suppressWarnings(as.numeric(num))
  den <- suppressWarnings(as.numeric(den))
  if (!is.finite(num) || !is.finite(den) || den <= 0) return("NA")
  paste0(format(round(100 * num / den, 2), trim = TRUE), "%")
}

.qc_safe_value <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(default)
  x <- suppressWarnings(as.numeric(x[[1L]]))
  if (is.finite(x)) x else default
}

.qc_metric_value <- function(qc_summary, metric, default = NA_real_) {
  if (!is.data.frame(qc_summary) || !all(c("metric", "value") %in% names(qc_summary))) return(default)
  hit <- qc_summary$value[as.character(qc_summary$metric) == metric]
  .qc_safe_value(hit, default = default)
}

.qc_table_html <- function(x, max_rows = 30L) {
  if (!is.data.frame(x) || !nrow(x)) return("<p class=\"empty\">No rows available.</p>")
  x <- as.data.frame(utils::head(x, max_rows), stringsAsFactors = FALSE)
  header <- paste0("<tr>", paste(sprintf("<th>%s</th>", .qc_html_escape(names(x))), collapse = ""), "</tr>")
  rows <- apply(x, 1L, function(row) {
    paste0("<tr>", paste(sprintf("<td>%s</td>", .qc_html_escape(row)), collapse = ""), "</tr>")
  })
  paste0("<table>", header, paste(rows, collapse = ""), "</table>")
}

.qc_cards_html <- function(x) {
  if (!is.data.frame(x) || !all(c("label", "value") %in% names(x))) return("")
  cards <- apply(as.data.frame(x, stringsAsFactors = FALSE), 1L, function(row) {
    paste0(
      "<div class=\"card\"><div class=\"card-label\">", .qc_html_escape(row[["label"]]),
      "</div><div class=\"card-value\">", .qc_html_escape(row[["value"]]),
      "</div></div>"
    )
  })
  paste0("<div class=\"cards\">", paste(cards, collapse = ""), "</div>")
}

.qc_bar_svg <- function(x, label_col, value_col, title = NULL, width = 860L, bar_height = 24L) {
  if (!is.data.frame(x) || !nrow(x) || !all(c(label_col, value_col) %in% names(x))) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[[value_col]] <- suppressWarnings(as.numeric(x[[value_col]]))
  x <- x[is.finite(x[[value_col]]) & x[[value_col]] >= 0, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  max_val <- max(x[[value_col]], na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  left <- 230L
  right <- 120L
  top <- if (is.null(title)) 12L else 44L
  row_gap <- 10L
  height <- top + nrow(x) * (bar_height + row_gap) + 20L
  plot_w <- width - left - right
  rows <- lapply(seq_len(nrow(x)), function(i) {
    y <- top + (i - 1L) * (bar_height + row_gap)
    val <- x[[value_col]][[i]]
    bw <- max(1, plot_w * val / max_val)
    label <- substr(as.character(x[[label_col]][[i]]), 1L, 34L)
    paste0(
      "<text x=\"0\" y=\"", y + 17L, "\" class=\"axis-label\">", .qc_html_escape(label), "</text>",
      "<rect x=\"", left, "\" y=\"", y, "\" width=\"", bw, "\" height=\"", bar_height, "\" rx=\"3\" class=\"bar\"/>",
      "<text x=\"", left + bw + 8L, "\" y=\"", y + 17L, "\" class=\"value-label\">", .qc_format_number(val), "</text>"
    )
  })
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    paste(rows, collapse = ""),
    "</svg>"
  )
}

.qc_hist_svg <- function(values, title = NULL, bins = 20L, width = 860L, height = 280L) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (!length(values)) return("<p class=\"empty\">No finite plot data available.</p>")
  h <- graphics::hist(values, breaks = bins, plot = FALSE)
  counts <- as.numeric(h$counts)
  max_count <- max(counts, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1
  left <- 56L
  right <- 18L
  top <- if (is.null(title)) 16L else 44L
  bottom <- 40L
  plot_w <- width - left - right
  plot_h <- height - top - bottom
  n <- length(counts)
  bar_w <- plot_w / max(1L, n)
  bars <- lapply(seq_len(n), function(i) {
    bh <- plot_h * counts[[i]] / max_count
    x <- left + (i - 1L) * bar_w
    y <- top + plot_h - bh
    paste0("<rect x=\"", x, "\" y=\"", y, "\" width=\"", max(1, bar_w - 1), "\" height=\"", bh, "\" class=\"bar\"/>")
  })
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    "<line x1=\"", left, "\" y1=\"", top + plot_h, "\" x2=\"", width - right, "\" y2=\"", top + plot_h, "\" class=\"axis\"/>",
    "<text x=\"", left, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(min(values))), "</text>",
    "<text x=\"", width - right - 90L, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(max(values))), "</text>",
    paste(bars, collapse = ""),
    "</svg>"
  )
}

.qc_relative_path <- function(path, from_dir) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  rp <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE), error = function(e) path)
  rd <- tryCatch(normalizePath(from_dir, winslash = "/", mustWork = FALSE), error = function(e) from_dir)
  split_path <- function(x) strsplit(x, "/", fixed = TRUE)[[1L]]
  rp_parts <- split_path(rp)
  rd_parts <- split_path(rd)
  common <- 0L
  max_common <- min(length(rp_parts), length(rd_parts))
  while (common < max_common && identical(rp_parts[[common + 1L]], rd_parts[[common + 1L]])) {
    common <- common + 1L
  }
  up <- rep("..", length(rd_parts) - common)
  down <- rp_parts[(common + 1L):length(rp_parts)]
  down <- down[!is.na(down) & nzchar(down)]
  rel <- paste(c(up, down), collapse = "/")
  if (!nzchar(rel)) rel <- basename(rp)
  utils::URLencode(rel, reserved = FALSE)
}

.qc_links_html <- function(paths, from_dir) {
  paths <- unique(as.character(paths))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return("<p class=\"empty\">No related report files found.</p>")
  items <- vapply(paths, function(path) {
    rel <- .qc_relative_path(path, from_dir)
    paste0("<li><a href=\"", .qc_html_escape(rel), "\">", .qc_html_escape(basename(path)), "</a></li>")
  }, character(1L))
  paste0("<ul class=\"links\">", paste(items, collapse = ""), "</ul>")
}

.qc_write_html <- function(path, title, sections) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  css <- paste(
    "body{margin:0;background:#f4f7fb;color:#172033;font-family:Arial,Helvetica,sans-serif}",
    ".wrap{max-width:1180px;margin:0 auto;padding:28px}",
    "header{background:#101827;color:white;border-bottom:5px solid #20b2aa}",
    "header .wrap{padding-top:24px;padding-bottom:24px}",
    "h1{font-size:28px;line-height:1.15;margin:0 0 8px 0}",
    "h2{font-size:19px;margin:0 0 12px 0;color:#101827}",
    "section{background:white;border:1px solid #d9e2ef;border-radius:8px;margin:18px 0;padding:18px;box-shadow:0 1px 2px rgba(17,24,39,.04)}",
    ".subtitle{color:#cbd5e1;font-size:13px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}",
    ".card{border:1px solid #d9e2ef;border-radius:7px;padding:12px;background:#fbfdff}.card-label{font-size:12px;color:#536173;font-weight:700}.card-value{font-size:24px;font-weight:800;margin-top:5px;color:#0f766e}",
    "table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid #e5edf6;padding:7px 8px;text-align:left}th{background:#edf4fb;color:#253246}",
    ".qc-plot{width:100%;height:auto;max-height:520px}.bar{fill:#168b87}.axis{stroke:#44546a;stroke-width:1}.axis-label,.value-label,.tick{font-size:12px;fill:#253246}.plot-title{font-size:15px;font-weight:700;fill:#101827}",
    ".empty{color:#69788c;font-style:italic}.status-pass{color:#0f766e;font-weight:800}.status-warn{color:#b45309;font-weight:800}.links{columns:2;line-height:1.8}a{color:#0f5f8f;text-decoration:none}a:hover{text-decoration:underline}",
    sep = "\n"
  )
  body <- paste(sections, collapse = "\n")
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>", .qc_html_escape(title), "</title><style>", css, "</style></head><body>",
    "<header><div class=\"wrap\"><h1>", .qc_html_escape(title), "</h1><div class=\"subtitle\">Generated by CraftGRN package QC reporting.</div></div></header>",
    "<main class=\"wrap\">", body, "</main></body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.qc_section <- function(title, body) {
  paste0("<section><h2>", .qc_html_escape(title), "</h2>", body, "</section>")
}

.qc_read_table_file <- function(path, format = NULL, columns = NULL) {
  if (is.null(path) || !file.exists(path)) return(tibble::tibble())
  fmt <- if (is.null(format) || is.na(format) || !nzchar(format)) {
    if (grepl("\\.parquet$", path, ignore.case = TRUE)) "parquet" else "csv"
  } else {
    as.character(format)
  }
  .module2_read_predicted_chunk(path, fmt, columns = columns)
}

.qc_read_manifest_chunks <- function(manifest_path) {
  if (is.null(manifest_path) || !file.exists(manifest_path)) return(tibble::tibble())
  tibble::as_tibble(readr::read_csv(manifest_path, show_col_types = FALSE))
}

.qc_find_module1_omics_rds <- function(module1_dir) {
  hits <- list.files(module1_dir, pattern = "^01_multiomic_data_object_.*\\.rds$", full.names = TRUE)
  if (length(hits)) hits[[1L]] else NA_character_
}

.module1_qc_scan_predicted_tfbs <- function(manifest, top_n = 20L, verbose = TRUE) {
  condition_support <- fp_id <- mean_condition_support <- n_fp <- n_predicted_tfbs <- tf <- N <- NULL
  if (!is.data.frame(manifest) || !nrow(manifest)) {
    return(list(tf_summary = tibble::tibble(), condition_support = tibble::tibble()))
  }
  top_n <- max(1L, as.integer(top_n[[1L]]))
  tf_rows <- vector("list", nrow(manifest))
  support_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    path <- as.character(manifest$path[[i]])
    fmt <- as.character(manifest$format[[i]])
    dt <- tryCatch(
      .qc_read_table_file(path, fmt, columns = c("tf", "fp_id", "condition_support")),
      error = function(e) .qc_read_table_file(path, fmt, columns = c("tf", "fp_id"))
    )
    dt <- data.table::as.data.table(dt)
    if (!nrow(dt)) next
    if (!"condition_support" %in% names(dt)) dt[, condition_support := NA_real_]
    tf_rows[[i]] <- dt[, .(
      n_predicted_tfbs = .N,
      n_fp = data.table::uniqueN(fp_id),
      mean_condition_support = mean(suppressWarnings(as.numeric(condition_support)), na.rm = TRUE)
    ), by = .(tf)]
    support_rows[[i]] <- dt[, .N, by = .(condition_support)]
    if (isTRUE(verbose)) .log_inform("Module 1 QC report: scanned predicted TFBS chunk {i}/{nrow(manifest)}.")
  }
  tf_summary <- data.table::rbindlist(tf_rows, use.names = TRUE, fill = TRUE)
  if (nrow(tf_summary)) {
    tf_summary <- tf_summary[, .(
      n_predicted_tfbs = sum(n_predicted_tfbs, na.rm = TRUE),
      n_fp = sum(n_fp, na.rm = TRUE),
      mean_condition_support = stats::weighted.mean(mean_condition_support, n_predicted_tfbs, na.rm = TRUE)
    ), by = tf]
    data.table::setorder(tf_summary, -n_predicted_tfbs, tf)
    tf_summary <- head(tf_summary, top_n)
  }
  support <- data.table::rbindlist(support_rows, use.names = TRUE, fill = TRUE)
  if (nrow(support)) {
    support <- support[, .(n_predicted_tfbs = sum(N, na.rm = TRUE)), by = condition_support]
    data.table::setorder(support, condition_support)
  }
  list(tf_summary = tibble::as_tibble(tf_summary), condition_support = tibble::as_tibble(support))
}

.module1_qc_input_metrics <- function(omics_data) {
  if (!is_multiomic_object(omics_data)) return(tibble::tibble(label = character(), value = character()))
  mats <- omics_data$matrices
  n_tf <- NA_integer_
  if (is.data.frame(omics_data$features$gene) && "is_tf" %in% names(omics_data$features$gene)) {
    n_tf <- sum(omics_data$features$gene$is_tf %in% TRUE, na.rm = TRUE)
  }
  tibble::tibble(
    label = c("Conditions", "Aligned FPs", "Bound FPs", "ATAC peaks", "Genes", "Expressed genes", "TF genes"),
    value = c(
      .qc_format_number(ncol(mats$fp_score)),
      .qc_format_number(nrow(mats$fp_score)),
      .qc_format_number(sum(rowSums(mats$fp_bound > 0, na.rm = TRUE) > 0)),
      .qc_format_number(nrow(mats$atac_score)),
      .qc_format_number(nrow(mats$gene_expr)),
      .qc_format_number(sum(rowSums(mats$gene_on > 0, na.rm = TRUE) > 0)),
      .qc_format_number(n_tf)
    )
  )
}

.module1_qc_motif_complexity <- function(omics_data) {
  fp_id <- motif <- n_tfs <- tf <- NULL
  if (!is_multiomic_object(omics_data) || !is.data.frame(omics_data$features$fp_motif)) {
    return(list(summary = tibble::tibble(), values = numeric()))
  }
  fm <- data.table::as.data.table(omics_data$features$fp_motif)
  if (!all(c("fp_id", "motif", "tf") %in% names(fm))) return(list(summary = tibble::tibble(), values = numeric()))
  counts <- fm[, .(n_motifs = data.table::uniqueN(motif), n_tfs = data.table::uniqueN(tf)), by = fp_id]
  summary <- counts[, .(
    n_fp = .N,
    median_tfs_per_fp = stats::median(n_tfs, na.rm = TRUE),
    p95_tfs_per_fp = as.numeric(stats::quantile(n_tfs, 0.95, na.rm = TRUE)),
    max_tfs_per_fp = max(n_tfs, na.rm = TRUE)
  )]
  list(summary = tibble::as_tibble(summary), values = counts$n_tfs)
}

#' Build a Module 1 QC HTML report
#'
#' Builds a compact HTML report for Module 1 input gates, canonical support,
#' prediction output scale, and related QC artifacts. The report can consume a
#' `predict_tfbs()` result, a step-by-step Module 1 result list, or a Module 1
#' output directory.
#'
#' @param module1 Module 1 result list or Module 1 output directory.
#' @param omics_data Optional compact multiomic object. Used when `module1` is
#'   an output directory or does not contain `omics_data`.
#' @param output_dir Directory where the HTML report should be written. If
#'   `NULL`, the report is written under `reports` inside the Module 1 output
#'   directory when available.
#' @param report_name HTML report filename.
#' @param scan_predicted_tfbs Logical; if `TRUE`, scan predicted TFBS chunks to
#'   summarize top TFs and condition support. This is comprehensive but can take
#'   extra time on full projects.
#' @param top_n Number of TFs to show in top-TF summaries.
#' @param verbose Emit concise progress messages.
#' @return Normalized path to the HTML report.
#' @export
build_module1_qc_report <- function(module1,
                                    omics_data = NULL,
                                    output_dir = NULL,
                                    report_name = "module1_qc_report.html",
                                    scan_predicted_tfbs = TRUE,
                                    top_n = 20L,
                                    verbose = TRUE) {
  module1_dir <- NULL
  if (is.character(module1) && length(module1) == 1L) {
    module1_dir <- module1
    module1 <- list()
  }
  if (is.null(omics_data) && is.list(module1) && is_multiomic_object(module1$omics_data)) {
    omics_data <- module1$omics_data
  }
  if (is.null(omics_data) && !is.null(module1_dir)) {
    rds <- .qc_find_module1_omics_rds(module1_dir)
    if (!is.na(rds) && file.exists(rds)) omics_data <- readRDS(rds)
  }
  if (is.null(output_dir)) {
    output_dir <- if (!is.null(module1_dir)) file.path(module1_dir, "reports") else file.path(getwd(), "module1_reports")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  qc_summary <- NULL
  if (is.list(module1) && is.list(module1$parameters$qc_summary)) {
    qc_summary <- tibble::tibble(metric = names(module1$parameters$qc_summary), value = as.numeric(unlist(module1$parameters$qc_summary, use.names = FALSE)))
  }
  if (is.null(qc_summary) && !is.null(module1_dir)) {
    qpath <- file.path(module1_dir, "module1_qc_summary.csv")
    if (file.exists(qpath)) qc_summary <- readr::read_csv(qpath, show_col_types = FALSE)
  }
  if (is.null(qc_summary)) qc_summary <- tibble::tibble(metric = character(), value = numeric())

  pred_manifest_path <- NULL
  if (is.list(module1)) pred_manifest_path <- module1$reports$predicted_tfbs_manifest %||% module1$predicted_tfbs_paths$manifest
  if ((is.null(pred_manifest_path) || !file.exists(pred_manifest_path)) && !is.null(module1_dir)) {
    pred_manifest_path <- file.path(module1_dir, "module1_predicted_tfbs_manifest.csv")
  }
  pred_manifest <- .qc_read_manifest_chunks(pred_manifest_path)
  support_check <- tibble::tibble()
  if (!is.null(module1_dir)) {
    support_path <- file.path(module1_dir, "module1_predicted_tfbs_canonical_support_check.csv")
    if (file.exists(support_path)) support_check <- readr::read_csv(support_path, show_col_types = FALSE)
  }

  predicted_scan <- if (isTRUE(scan_predicted_tfbs)) {
    .module1_qc_scan_predicted_tfbs(pred_manifest, top_n = top_n, verbose = verbose)
  } else {
    list(tf_summary = tibble::tibble(), condition_support = tibble::tibble())
  }
  motif_complexity <- .module1_qc_motif_complexity(omics_data)
  input_cards <- .module1_qc_input_metrics(omics_data)
  predicted_rows <- if (nrow(pred_manifest)) sum(suppressWarnings(as.numeric(pred_manifest$n_rows)), na.rm = TRUE) else .qc_metric_value(qc_summary, "n_predicted_tfbs")
  canonical_missing <- .qc_metric_value(support_check, "predicted_fp_without_motif_supported_predicted_tf", default = NA_real_)
  canonical_status <- if (is.finite(canonical_missing) && canonical_missing == 0) {
    "<span class=\"status-pass\">PASS</span>"
  } else if (is.finite(canonical_missing)) {
    "<span class=\"status-warn\">CHECK</span>"
  } else {
    "<span class=\"status-warn\">NOT AVAILABLE</span>"
  }

  cards <- rbind(
    input_cards,
    tibble::tibble(
      label = c("Predicted TFBS rows", "Predicted chunks", "Canonical support"),
      value = c(.qc_format_number(predicted_rows), .qc_format_number(nrow(pred_manifest)), gsub("<[^>]+>", "", canonical_status))
    )
  )
  funnel <- tibble::tibble(
    step = c("Input FPs", "Bound FPs", "Canonical-bound FPs", "Prediction FPs", "Predicted TFBS"),
    n = c(
      .qc_metric_value(qc_summary, "n_fp_input", default = if (is_multiomic_object(omics_data)) nrow(omics_data$matrices$fp_score) else NA_real_),
      .qc_metric_value(qc_summary, "n_fp_bound_accessible", default = if (is_multiomic_object(omics_data)) sum(rowSums(omics_data$matrices$fp_bound > 0, na.rm = TRUE) > 0) else NA_real_),
      .qc_metric_value(qc_summary, "n_canonical_bound_fps", default = .qc_metric_value(support_check, "predicted_unique_fp")),
      .qc_metric_value(qc_summary, "n_prediction_fps", default = .qc_metric_value(support_check, "predicted_unique_fp")),
      predicted_rows
    )
  )
  related <- character()
  if (!is.null(module1_dir)) {
    related <- c(
      list.files(module1_dir, pattern = "\\.(pdf|html)$", full.names = TRUE),
      pred_manifest_path
    )
  }
  sections <- list(
    .qc_section("Summary", .qc_cards_html(cards)),
    .qc_section("Correctness Checks", paste0(
      "<p>Canonical-supported predicted FP check: ", canonical_status, "</p>",
      .qc_table_html(support_check, max_rows = 20L)
    )),
    .qc_section("Workflow Funnel", .qc_bar_svg(funnel, "step", "n", title = "Module 1 row counts by processing stage")),
    .qc_section("Predicted TFBS Chunks", paste0(
      .qc_bar_svg(pred_manifest, "chunk_id", "n_rows", title = "Predicted TFBS rows per chunk"),
      .qc_table_html(pred_manifest, max_rows = 25L)
    )),
    .qc_section("Top Predicted TFs", paste0(
      .qc_bar_svg(predicted_scan$tf_summary, "tf", "n_predicted_tfbs", title = "Top TFs by predicted TFBS rows"),
      .qc_table_html(predicted_scan$tf_summary, max_rows = top_n)
    )),
    .qc_section("Motif Complexity", paste0(
      .qc_hist_svg(motif_complexity$values, title = "TFs per aligned FP from canonical motif support", bins = 30L),
      .qc_table_html(motif_complexity$summary, max_rows = 10L)
    )),
    .qc_section("Related Files", .qc_links_html(related, from_dir = output_dir))
  )
  out <- file.path(output_dir, report_name)
  path <- .qc_write_html(out, "Module 1 QC Report", sections)
  if (isTRUE(verbose)) .log_inform("Module 1 QC HTML report written: {path}")
  path
}

.module2_qc_table_from_module <- function(module2, table_name, columns = NULL) {
  tryCatch(.module2_report_read_table(module2, table_name, columns = columns), error = function(e) tibble::tibble())
}

.module2_qc_manifest_for <- function(module2, name) {
  direct <- module2[[paste0(name, "_manifest")]]
  if (is.data.frame(direct)) return(direct)
  man <- module2$manifest
  if (!is.data.frame(man) || !nrow(man)) return(tibble::tibble())
  hit <- man[as.character(man$table) == name, , drop = FALSE]
  if (!nrow(hit) || !identical(as.character(hit$format[[1L]]), "manifest")) return(tibble::tibble())
  .qc_read_manifest_chunks(as.character(hit$path[[1L]]))
}

.module2_qc_scan_candidates <- function(module2, bins = 40L, verbose = TRUE) {
  candidate_source <- distance_to_tss <- n_candidates <- NULL
  man <- .module2_qc_manifest_for(module2, "module2_fp_target_candidates")
  if (!nrow(man)) return(list(distance_values = numeric(), source_summary = tibble::tibble()))
  dist_sample <- numeric()
  src_rows <- vector("list", nrow(man))
  for (i in seq_len(nrow(man))) {
    dt <- data.table::as.data.table(.qc_read_table_file(as.character(man$path[[i]]), as.character(man$format[[i]]), columns = c("distance_to_tss", "candidate_source", "prior_supported")))
    if (!nrow(dt)) next
    src_rows[[i]] <- dt[, .(n_candidates = .N), by = .(candidate_source)]
    d <- suppressWarnings(as.numeric(dt$distance_to_tss))
    d <- d[is.finite(d)]
    if (length(d)) dist_sample <- c(dist_sample, if (length(d) > 50000L) utils::head(d, 50000L) else d)
    if (isTRUE(verbose)) .log_inform("Module 2 QC report: scanned candidate chunk {i}/{nrow(man)}.")
  }
  src <- data.table::rbindlist(src_rows, use.names = TRUE, fill = TRUE)
  if (nrow(src)) src <- src[, .(n_candidates = sum(n_candidates, na.rm = TRUE)), by = candidate_source]
  list(distance_values = dist_sample, source_summary = tibble::as_tibble(src))
}

.module2_qc_scan_links <- function(module2, top_n = 20L, validate_integrity = TRUE, verbose = TRUE) {
  fp_id <- module2_link_pass <- n_fp <- n_links <- n_target_genes <- pass <- target_gene <- tf <- NULL
  link_man <- .module2_qc_manifest_for(module2, "module2_links")
  if (!nrow(link_man) && is.data.frame(module2$links) && nrow(module2$links)) {
    dt <- data.table::as.data.table(module2$links)
    top <- dt[, .(n_links = .N, n_fp = data.table::uniqueN(fp_id), n_target_genes = data.table::uniqueN(target_gene)), by = tf]
    data.table::setorder(top, -n_links, tf)
    return(list(top_tf = tibble::as_tibble(head(top, top_n)), validation = tibble::tibble(metric = "n_links_scanned", value = nrow(dt))))
  }
  if (!nrow(link_man)) return(list(top_tf = tibble::tibble(), validation = tibble::tibble()))
  tf_keys <- character()
  fp_keys <- character()
  if (isTRUE(validate_integrity)) {
    tf_corr <- .module2_qc_table_from_module(module2, "module2_tf_target_corr", columns = c("tf", "target_gene", "pass"))
    fp_corr <- .module2_qc_table_from_module(module2, "module2_fp_target_corr", columns = c("fp_id", "target_gene", "pass"))
    tf_dt <- data.table::as.data.table(tf_corr)
    fp_dt <- data.table::as.data.table(fp_corr)
    if (nrow(tf_dt)) tf_keys <- paste(tf_dt[pass %in% TRUE]$tf, tf_dt[pass %in% TRUE]$target_gene, sep = "\r")
    if (nrow(fp_dt)) fp_keys <- paste(fp_dt[pass %in% TRUE]$fp_id, fp_dt[pass %in% TRUE]$target_gene, sep = "\r")
  }
  rows <- vector("list", nrow(link_man))
  n_scanned <- 0
  n_bad_tf <- 0
  n_bad_fp <- 0
  for (i in seq_len(nrow(link_man))) {
    dt <- data.table::as.data.table(.qc_read_table_file(as.character(link_man$path[[i]]), as.character(link_man$format[[i]]), columns = c("tf", "fp_id", "target_gene", "module2_link_pass")))
    if (!nrow(dt)) next
    dt <- dt[module2_link_pass %in% TRUE]
    n_scanned <- n_scanned + nrow(dt)
    rows[[i]] <- dt[, .(n_links = .N, n_fp = data.table::uniqueN(fp_id), n_target_genes = data.table::uniqueN(target_gene)), by = tf]
    if (isTRUE(validate_integrity) && length(tf_keys)) n_bad_tf <- n_bad_tf + sum(!paste(dt$tf, dt$target_gene, sep = "\r") %in% tf_keys)
    if (isTRUE(validate_integrity) && length(fp_keys)) n_bad_fp <- n_bad_fp + sum(!paste(dt$fp_id, dt$target_gene, sep = "\r") %in% fp_keys)
    if (isTRUE(verbose)) .log_inform("Module 2 QC report: scanned link chunk {i}/{nrow(link_man)}.")
  }
  top <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(top)) {
    top <- top[, .(n_links = sum(n_links, na.rm = TRUE), n_fp = sum(n_fp, na.rm = TRUE), n_target_genes = sum(n_target_genes, na.rm = TRUE)), by = tf]
    data.table::setorder(top, -n_links, tf)
    top <- head(top, top_n)
  }
  validation <- tibble::tibble(
    metric = c("n_links_scanned", "n_links_with_missing_tf_target_pass", "n_links_with_missing_fp_target_pass"),
    value = c(n_scanned, n_bad_tf, n_bad_fp)
  )
  list(top_tf = tibble::as_tibble(top), validation = validation)
}

.module2_qc_related_html <- function(module2_dir, output_dir) {
  if (is.null(module2_dir) || !dir.exists(module2_dir)) return(character())
  hits <- list.files(module2_dir, pattern = "\\.html$", full.names = TRUE, recursive = TRUE)
  hits[normalizePath(hits, winslash = "/", mustWork = FALSE) != normalizePath(file.path(output_dir, "module2_qc_report.html"), winslash = "/", mustWork = FALSE)]
}

#' Build a Module 2 QC HTML report
#'
#' Builds a compact HTML report for Module 2 correlation filters, candidate
#' construction, final TF-FP-target links, integrity checks, and related browser
#' reports.
#'
#' @param module2 Module 2 result list, loaded Module 2 list, or output
#'   directory.
#' @param multiomic_data Optional compact multiomic object used for context.
#' @param output_dir Directory where the HTML report should be written. If
#'   `NULL`, the report is written under `reports` inside the Module 2 output
#'   directory when available.
#' @param report_name HTML report filename.
#' @param scan_large_tables Logical; if `TRUE`, scan candidate and link chunks
#'   for top-TF, distance, and integrity summaries.
#' @param validate_integrity Logical; if `TRUE`, verify final links against
#'   passing TF-target and FP-target keys while scanning link chunks.
#' @param top_n Number of TFs to show in top-TF summaries.
#' @param verbose Emit concise progress messages.
#' @return Normalized path to the HTML report.
#' @export
build_module2_qc_report <- function(module2,
                                    multiomic_data = NULL,
                                    output_dir = NULL,
                                    report_name = "module2_qc_report.html",
                                    scan_large_tables = TRUE,
                                    validate_integrity = TRUE,
                                    top_n = 20L,
                                    verbose = TRUE) {
  module2_dir <- NULL
  if (is.character(module2) && length(module2) == 1L) {
    module2_dir <- if (dir.exists(module2)) module2 else dirname(module2)
    module2 <- load_module2_links(module2)
  } else if (is.list(module2) && is.character(module2$reports$manifest)) {
    module2_dir <- dirname(module2$reports$manifest)
  }
  if (is.null(output_dir)) {
    output_dir <- if (!is.null(module2_dir)) file.path(module2_dir, "reports") else file.path(getwd(), "module2_reports")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  qc_summary <- .module2_qc_table_from_module(module2, "module2_qc_summary")
  if (!nrow(qc_summary) && is.data.frame(module2$qc_summary)) qc_summary <- module2$qc_summary
  manifest <- if (is.data.frame(module2$manifest)) module2$manifest else tibble::tibble()
  manifest_checks <- if (nrow(manifest) && all(c("table", "path", "n_rows") %in% names(manifest))) {
    tibble::tibble(
      table = as.character(manifest$table),
      n_rows = .qc_format_number(manifest$n_rows),
      file_exists = ifelse(file.exists(as.character(manifest$path)), "yes", "no"),
      format = as.character(manifest$format)
    )
  } else {
    tibble::tibble()
  }
  candidate_scan <- if (isTRUE(scan_large_tables)) .module2_qc_scan_candidates(module2, verbose = verbose) else list(distance_values = numeric(), source_summary = tibble::tibble())
  link_scan <- if (isTRUE(scan_large_tables)) .module2_qc_scan_links(module2, top_n = top_n, validate_integrity = validate_integrity, verbose = verbose) else list(top_tf = tibble::tibble(), validation = tibble::tibble())
  cards <- tibble::tibble(
    label = c("Predicted TFBS", "TFs", "Target genes", "TF-target pass", "FP-target pass", "Final links"),
    value = c(
      .qc_format_number(.qc_metric_value(qc_summary, "n_predicted_tfbs")),
      .qc_format_number(.qc_metric_value(qc_summary, "n_tfs")),
      .qc_format_number(.qc_metric_value(qc_summary, "n_target_genes")),
      paste0(.qc_format_number(.qc_metric_value(qc_summary, "n_tf_target_pairs_pass")), " / ", .qc_format_number(.qc_metric_value(qc_summary, "n_tf_target_pairs_tested"))),
      paste0(.qc_format_number(.qc_metric_value(qc_summary, "n_fp_target_pairs_pass")), " / ", .qc_format_number(.qc_metric_value(qc_summary, "n_fp_target_pairs_tested"))),
      .qc_format_number(.qc_metric_value(qc_summary, "n_module2_links"))
    )
  )
  funnel <- tibble::tibble(
    step = c("Predicted TFBS", "TF-target tested", "TF-target pass", "FP-target candidates", "FP-target pass", "Final links"),
    n = c(
      .qc_metric_value(qc_summary, "n_predicted_tfbs"),
      .qc_metric_value(qc_summary, "n_tf_target_pairs_tested"),
      .qc_metric_value(qc_summary, "n_tf_target_pairs_pass"),
      .qc_metric_value(qc_summary, "n_fp_target_candidates"),
      .qc_metric_value(qc_summary, "n_fp_target_pairs_pass"),
      .qc_metric_value(qc_summary, "n_module2_links")
    )
  )
  integrity_status <- "NOT RUN"
  if (nrow(link_scan$validation)) {
    bad_tf <- .qc_metric_value(link_scan$validation, "n_links_with_missing_tf_target_pass", default = NA_real_)
    bad_fp <- .qc_metric_value(link_scan$validation, "n_links_with_missing_fp_target_pass", default = NA_real_)
    integrity_status <- if (is.finite(bad_tf) && is.finite(bad_fp) && bad_tf == 0 && bad_fp == 0) "PASS" else "CHECK"
  }
  integrity_class <- if (identical(integrity_status, "PASS")) "status-pass" else "status-warn"
  related <- .module2_qc_related_html(module2_dir, output_dir)
  sections <- list(
    .qc_section("Summary", .qc_cards_html(cards)),
    .qc_section("Workflow Funnel", .qc_bar_svg(funnel, "step", "n", title = "Module 2 row counts by processing stage")),
    .qc_section("Integrity Checks", paste0(
      "<p>Final-link relational integrity: <span class=\"", integrity_class, "\">", integrity_status, "</span></p>",
      .qc_table_html(link_scan$validation, max_rows = 20L),
      .qc_table_html(manifest_checks, max_rows = 20L)
    )),
    .qc_section("Candidate Distance To TSS", paste0(
      .qc_hist_svg(candidate_scan$distance_values, title = "Candidate FP distance to target TSS", bins = 60L),
      .qc_table_html(candidate_scan$source_summary, max_rows = 20L)
    )),
    .qc_section("Top TFs In Final Links", paste0(
      .qc_bar_svg(link_scan$top_tf, "tf", "n_links", title = "Top TFs by final Module 2 links"),
      .qc_table_html(link_scan$top_tf, max_rows = top_n)
    )),
    .qc_section("QC Summary Table", .qc_table_html(qc_summary, max_rows = 50L)),
    .qc_section("Related HTML Reports", .qc_links_html(related, from_dir = output_dir))
  )
  out <- file.path(output_dir, report_name)
  path <- .qc_write_html(out, "Module 2 QC Report", sections)
  if (isTRUE(verbose)) .log_inform("Module 2 QC HTML report written: {path}")
  path
}
