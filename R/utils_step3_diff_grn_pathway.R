# File: utils_step3_diff_grn_pathway.R
# Author: Yaoxiang Li
# Created: 2026-05-01
#
# Purpose:
# Provide reusable pathway enrichment and plotting utilities for Step3
# differential GRN outputs.
#
# Notes:
# - This file is project-agnostic. Project-specific condition ordering should be
#   supplied by callers through function arguments.
# - The standard input is a diff_links_filtered directory containing
#   *_filtered_links_up.csv and *_filtered_links_down.csv files.

.default_pathway_databases <- function() {
  c(
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "Reactome_2022",
    "WikiPathways_2024_Human",
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human"
  )
}

#' Report differential pathway enrichment
#'
#' @description
#' Runs pathway enrichment from Module 3 filtered differential links and writes
#' user-facing pathway report outputs. Supporting enrichment tables are written
#' under `output_dir/pathway_tables`.
#'
#' When `pathway_databases = NULL`, CraftGRN uses the package default Enrichr
#' databases: GO Biological Process 2023, GO Cellular Component 2023, GO
#' Molecular Function 2023, Reactome 2022, WikiPathways 2024 Human, MSigDB
#' Hallmark 2020, and KEGG 2021 Human.
#'
#' @param filtered_dir Directory containing files named
#'   \code{*_filtered_links_up.csv} and \code{*_filtered_links_down.csv}.
#' @param output_dir Directory where user-facing report outputs are written.
#' @param comparison_file Optional comparison design CSV.
#' @param pathway_databases Optional character vector of Enrichr database names.
#' @param tag Output filename and plot title tag.
#' @param overwrite If TRUE, recompute existing enrichment tables.
#' @param make_dotplot If TRUE, write the pathway term dotplot PDF.
#' @param gene_col Gene column in filtered link files.
#' @param padj_cut Adjusted p-value cutoff for significant enrichment outputs.
#' @param min_genes Minimum distinct genes required before running Enrichr.
#' @param verbose Emit concise progress messages.
#' @param top_n Number of top pathways selected per comparison and direction.
#' @param padj_display_cut Adjusted p-value cutoff used before top-n selection.
#' @param condition_order Optional project-specific condition order.
#' @param cell_line_order Optional project-specific cell line order.
#' @param stress_order Optional project-specific stress order.
#' @param direction_order Direction facet/order.
#' @param size_metric Dot size metric passed to the internal dotplot builder.
#' @param sig_color_break Color and border threshold passed to the internal
#'   dotplot builder.
#' @param panel_width Optional fixed dot-panel width in inches.
#' @param panel_height Optional fixed dot-panel height in inches.
#' @param font_family Optional PDF font family.
#'
#' @return Invisibly returns a list with output paths and selected terms.
#' @noRd
report_differential_pathways <- function(filtered_dir,
                                         output_dir,
                                         comparison_file = NULL,
                                         pathway_databases = NULL,
                                         tag = "module3",
                                         overwrite = FALSE,
                                         make_dotplot = TRUE,
                                         gene_col = "gene_key",
                                         padj_cut = 0.05,
                                         min_genes = 5L,
                                         verbose = TRUE,
                                         top_n = 5L,
                                         padj_display_cut = 0.05,
                                         condition_order = NULL,
                                         cell_line_order = NULL,
                                         stress_order = NULL,
                                         direction_order = c("up", "down"),
                                         size_metric = c("odds_ratio", "combined_score", "log10_combined_score", "overlap_percent"),
                                         sig_color_break = 1.3,
                                         panel_width = NULL,
                                         panel_height = NULL,
                                         font_family = NULL) {
  size_metric <- match.arg(size_metric)
  dbs <- if (is.null(pathway_databases)) .default_pathway_databases() else pathway_databases
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(output_dir, "pathway_tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  enrichment_files <- run_diff_grn_pathway_enrichment(
    filtered_dir = filtered_dir,
    out_dir = table_dir,
    dbs = dbs,
    gene_col = gene_col,
    padj_cut = padj_cut,
    min_genes = min_genes,
    overwrite = overwrite,
    verbose = verbose
  )
  enrich_dt <- read_diff_grn_pathway_enrichment(
    pathway_dir = table_dir,
    comparison_file = comparison_file
  )
  selected_terms <- if (nrow(enrich_dt)) {
    select_diff_grn_pathway_terms(
      enrich_dt = enrich_dt,
      top_n = top_n,
      padj_display_cut = padj_display_cut,
      condition_order = condition_order,
      cell_line_order = cell_line_order,
      stress_order = stress_order,
      direction_order = direction_order
    )
  } else {
    data.table::data.table(
      comparison_id = character(),
      direction = character(),
      database = character(),
      pathway = character(),
      adjusted_p = numeric(),
      overlap_hits = integer()
    )
  }
  term_path <- file.path(table_dir, sprintf("%s_pathway_plot_terms.csv", tag))
  data.table::fwrite(selected_terms, term_path)
  dotplot_path <- file.path(output_dir, sprintf("%s_pathway_term_dotplot.pdf", tag))
  if (isTRUE(make_dotplot) && nrow(enrich_dt)) {
    plot_diff_grn_pathway_dotplot(
      enrich_dt = enrich_dt,
      out_file = dotplot_path,
      tag = tag,
      top_n = top_n,
      padj_display_cut = padj_display_cut,
      condition_order = condition_order,
      cell_line_order = cell_line_order,
      stress_order = stress_order,
      direction_order = direction_order,
      size_metric = size_metric,
      sig_color_break = sig_color_break,
      panel_width = panel_width,
      panel_height = panel_height,
      font_family = font_family
    )
  }
  invisible(list(
    enrichment_files = enrichment_files,
    selected_terms_file = term_path,
    dotplot_file = if (isTRUE(make_dotplot)) dotplot_path else NA_character_,
    selected_terms = selected_terms
  ))
}

#' Run pathway enrichment for filtered differential GRN links
#'
#' @description
#' Reads filtered differential-link CSV files, extracts distinct target genes,
#' runs Enrichr with the default Enrichr background, and writes full and
#' significant enrichment tables for each comparison and direction.
#'
#' @param filtered_dir Directory containing files named
#'   \code{*_filtered_links_up.csv} and \code{*_filtered_links_down.csv}.
#' @param out_dir Directory where pathway enrichment CSVs are written.
#' @param dbs Character vector of Enrichr database names.
#' @param gene_col Gene column in filtered link files.
#' @param padj_cut Adjusted p-value cutoff for \code{*_sig.csv} outputs.
#' @param min_genes Minimum distinct genes required before running Enrichr.
#' @param overwrite If TRUE, recompute existing enrichment CSVs.
#' @param enrichr_sleep_time Seconds to pause between Enrichr requests.
#' @param enrichr_cache_dir Optional directory for cached Enrichr responses.
#' @param pathway_backend Pathway backend, either `"enrichly"` for local cached
#'   enrichment or `"enrichr"` for the Enrichr web API.
#' @param verbose Emit concise progress messages.
#'
#' @return Invisibly returns a character vector of full enrichment CSV paths.
#' @noRd
run_diff_grn_pathway_enrichment <- function(filtered_dir,
                                            out_dir,
                                            dbs = .default_pathway_databases(),
                                            gene_col = "gene_key",
                                            padj_cut = 0.05,
                                            min_genes = 5L,
                                            overwrite = FALSE,
                                            enrichr_sleep_time = 0,
                                            enrichr_cache_dir = NULL,
                                            pathway_backend = NULL,
                                            verbose = TRUE) {
  if (!dir.exists(filtered_dir)) .log_abort("`filtered_dir` not found: {filtered_dir}")
  pathway_backend <- .pathway_backend(pathway_backend)
  if (!.pathway_backend_available(pathway_backend)) {
    .log_warn("Skipping pathway enrichment because neither enrichly nor enrichR is installed.")
    return(invisible(character(0)))
  }
  enrichr_sleep_time <- .normalize_enrichr_sleep_time(enrichr_sleep_time)
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- file.path(out_dir, "cache", "enrichr")
  }
  if (identical(pathway_backend, "enrichr")) {
    .ensure_enrichr_ready(site = "Enrichr", verbose = verbose)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- sort(list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE))
  if (!length(files)) {
    .log_warn("No filtered link files found in: {filtered_dir}")
    return(invisible(character(0)))
  }
  manifest <- vector("list", length(files))
  out_files <- character(0)
  for (i in seq_along(files)) {
    path <- files[[i]]
    stem <- sub("\\.csv$", "", basename(path))
    out_path <- file.path(out_dir, paste0(stem, "_pathway_enrichment.csv"))
    sig_path <- sub("_pathway_enrichment\\.csv$", "_pathway_enrichment_sig.csv", out_path)
    if (!isTRUE(overwrite) && file.exists(out_path)) {
      out_files <- c(out_files, out_path)
      manifest[[i]] <- data.table::data.table(
        input_file = basename(path),
        output_file = basename(out_path),
        sig_output_file = basename(sig_path),
        status = "exists"
      )
      next
    }
    header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
    if (!(gene_col %in% header)) {
      manifest[[i]] <- data.table::data.table(
        input_file = basename(path),
        output_file = basename(out_path),
        status = "missing_gene_col"
      )
      next
    }
    genes <- unique(data.table::fread(path, select = gene_col, showProgress = FALSE)[[gene_col]])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) < as.integer(min_genes)) {
      manifest[[i]] <- data.table::data.table(
        input_file = basename(path),
        output_file = basename(out_path),
        status = "too_few_genes",
        n_genes = length(genes)
      )
      next
    }
    enr <- tryCatch(
      .run_enrichr_cached(
        genes = genes,
        dbs = dbs,
        sleep_time = enrichr_sleep_time,
        cache_dir = enrichr_cache_dir,
        backend = pathway_backend
      ),
      error = function(e) e
    )
    if (inherits(enr, "error")) {
      manifest[[i]] <- data.table::data.table(
        input_file = basename(path),
        output_file = basename(out_path),
        status = "enrichr_error",
        message = conditionMessage(enr)
      )
      next
    }
    rows <- lapply(names(enr), function(db_name) {
      tbl <- data.table::as.data.table(enr[[db_name]])
      if (!nrow(tbl) || !("Adjusted.P.value" %in% names(tbl))) return(NULL)
      data.table::data.table(
        source_file = basename(path),
        database = db_name,
        pathway = as.character(tbl$Term),
        adjusted_p = suppressWarnings(as.numeric(tbl$Adjusted.P.value)),
        p_value = if ("P.value" %in% names(tbl)) suppressWarnings(as.numeric(tbl$P.value)) else NA_real_,
        odds_ratio = if ("Odds.Ratio" %in% names(tbl)) suppressWarnings(as.numeric(tbl$Odds.Ratio)) else NA_real_,
        combined_score = if ("Combined.Score" %in% names(tbl)) suppressWarnings(as.numeric(tbl$Combined.Score)) else NA_real_,
        overlap = if ("Overlap" %in% names(tbl)) as.character(tbl$Overlap) else NA_character_,
        genes = if ("Genes" %in% names(tbl)) as.character(tbl$Genes) else NA_character_
      )
    })
    res <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    res <- res[is.finite(adjusted_p)]
    if (!nrow(res)) {
      manifest[[i]] <- data.table::data.table(
        input_file = basename(path),
        output_file = basename(out_path),
        status = "no_finite_padj"
      )
      next
    }
    data.table::setorder(res, adjusted_p, database, pathway)
    data.table::fwrite(res, out_path)
    data.table::fwrite(res[adjusted_p <= as.numeric(padj_cut)], sig_path)
    out_files <- c(out_files, out_path)
    manifest[[i]] <- data.table::data.table(
      input_file = basename(path),
      output_file = basename(out_path),
      sig_output_file = basename(sig_path),
      status = "written",
      n_genes = length(genes),
      n_rows = nrow(res),
      n_sig_rows = nrow(res[adjusted_p <= as.numeric(padj_cut)])
    )
    if (isTRUE(verbose)) .log_inform("Saved pathway enrichment: {out_path}")
  }
  manifest_dt <- data.table::rbindlist(manifest, use.names = TRUE, fill = TRUE)
  manifest_path <- file.path(out_dir, "pathway_enrichment_manifest.csv")
  data.table::fwrite(manifest_dt, manifest_path)
  if (isTRUE(verbose)) .log_inform("Saved pathway enrichment manifest: {manifest_path}")
  invisible(out_files)
}

#' Read differential GRN pathway enrichment tables
#'
#' @param pathway_dir Directory containing \code{*_pathway_enrichment.csv}.
#' @param comparison_file Optional comparison design CSV with \code{cond1_label}
#'   and \code{cond2_label}. If present, it is used to validate and annotate
#'   comparisons.
#' @param term_max_chars Maximum displayed term label length.
#'
#' @return A data.table with enrichment rows and comparison metadata.
#' @noRd
read_diff_grn_pathway_enrichment <- function(pathway_dir,
                                             comparison_file = NULL,
                                             term_max_chars = 58L) {
  if (!dir.exists(pathway_dir)) .log_abort("`pathway_dir` not found: {pathway_dir}")
  files <- sort(list.files(pathway_dir, "_pathway_enrichment\\.csv$", full.names = TRUE))
  files <- files[!grepl("_pathway_enrichment_sig\\.csv$", basename(files))]
  files <- files[basename(files) != "pathway_enrichment_manifest.csv"]
  if (!length(files)) return(data.table::data.table())
  compar_dt <- NULL
  if (!is.null(comparison_file)) {
    if (!file.exists(comparison_file)) .log_abort("`comparison_file` not found: {comparison_file}")
    compar_dt <- data.table::fread(comparison_file, showProgress = FALSE)
    if (!all(c("cond1_label", "cond2_label") %in% names(compar_dt))) {
      .log_abort("`comparison_file` must contain cond1_label and cond2_label.")
    }
    if (!("comparison_id" %in% names(compar_dt))) {
      compar_dt[, comparison_id := paste0(cond1_label, "_vs_", cond2_label)]
    }
    compar_dt <- unique(compar_dt[, .(comparison_id, cond1_label, cond2_label)])
  }
  out <- data.table::rbindlist(lapply(files, function(path) {
    dt <- data.table::fread(path, showProgress = FALSE)
    source_file <- sub("_pathway_enrichment\\.csv$", ".csv", basename(path))
    source_file <- sub("_delta_links_filtered\\.csv$", "_filtered_links.csv", source_file)
    direction <- if (grepl("_filtered_links_up\\.csv$", source_file)) "up" else "down"
    comparison_id <- sub("_filtered_links_(up|down)\\.csv$", "", source_file)
    parts <- .diff_grn_pathway_comparison_parts(comparison_id)
    dt[, source_file := source_file]
    dt[, comparison_id := comparison_id]
    dt[, cond1_label := parts$cond1_label]
    dt[, cond2_label := parts$cond2_label]
    dt[, cell_line := parts$cell_line]
    dt[, direction := direction]
    dt
  }), use.names = TRUE, fill = TRUE)
  if (!is.null(compar_dt)) {
    out[, c("cond1_label", "cond2_label") := NULL]
    out <- merge(out, compar_dt, by = "comparison_id", all.x = TRUE, sort = FALSE)
    out[is.na(cond1_label), cond1_label := sub("_vs_.*$", "", comparison_id)]
    out[is.na(cond2_label), cond2_label := sub("^.*_vs_", "", comparison_id)]
    out[, cell_line := sub("_.*$", "", cond1_label)]
  }
  if (!("adjusted_p" %in% names(out)) && "adjusted_p_value" %in% names(out)) {
    data.table::setnames(out, "adjusted_p_value", "adjusted_p")
  }
  for (col in c("p_value", "odds_ratio", "combined_score")) {
    if (!(col %in% names(out))) out[, (col) := NA_real_]
    out[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  for (col in c("overlap", "genes")) {
    if (!(col %in% names(out))) out[, (col) := NA_character_]
  }
  out[, overlap_hits := .diff_grn_pathway_parse_overlap(overlap)]
  out[, term_clean := .diff_grn_pathway_shorten_term(pathway, max_chars = term_max_chars)]
  out[, group_id := paste(comparison_id, direction, sep = " | ")]
  out[, stress_label := mapply(function(label, cell) {
    sub(paste0("^", cell, "_"), "", label)
  }, cond1_label, cell_line, USE.NAMES = FALSE)]
  out[, x_label := paste(cell_line, stress_label, sep = " ")]
  out[, pathway_key := paste(database, pathway, sep = " | ")]
  data.table::setorder(out, adjusted_p, -overlap_hits, database, pathway)
  out[]
}

#' Select pathway terms for a differential GRN dotplot
#'
#' @param enrich_dt Data table from \code{read_diff_grn_pathway_enrichment()}.
#' @param top_n Number of top pathways selected per comparison and direction.
#' @param padj_display_cut Adjusted p-value cutoff used before top-n selection.
#' @param condition_order Optional character vector or data.frame defining
#'   comparison or condition order. A data.frame can contain \code{comparison_id},
#'   \code{cond1_label}, and optionally \code{direction}.
#' @param cell_line_order Optional cell line order.
#' @param stress_order Optional stress label order.
#' @param direction_order Direction facet/order, usually \code{c("up", "down")}.
#'
#' @return A data.table containing selected pathway rows, all occurrences of
#'   selected pathways across conditions, and \code{is_cluster_top_pathway}.
#' @noRd
select_diff_grn_pathway_terms <- function(enrich_dt,
                                          top_n = 5L,
                                          padj_display_cut = 0.05,
                                          condition_order = NULL,
                                          cell_line_order = NULL,
                                          stress_order = NULL,
                                          direction_order = c("up", "down")) {
  top_n <- max(1L, as.integer(top_n))
  if (!("adjusted_p" %in% names(enrich_dt))) {
    return(data.table::data.table(
      comparison_id = character(),
      direction = character(),
      database = character(),
      pathway = character(),
      adjusted_p = numeric(),
      overlap_hits = integer()
    ))
  }
  dt <- data.table::as.data.table(enrich_dt)[is.finite(adjusted_p)]
  if (!nrow(dt)) return(dt)
  group_order <- .diff_grn_pathway_group_order(
    dt,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  group_levels <- group_order$group_id
  top_dt <- dt[, {
    ranked <- .SD[order(adjusted_p, -overlap_hits, database, pathway)]
    sig <- ranked[adjusted_p <= as.numeric(padj_display_cut)]
    if (nrow(sig)) head(sig, top_n) else head(ranked, top_n)
  }, by = group_id]
  if (!nrow(top_dt)) return(top_dt)
  selected_keys <- unique(top_dt$pathway_key)
  out <- dt[pathway_key %in% selected_keys]
  sig_keys <- unique(out[adjusted_p <= as.numeric(padj_display_cut), pathway_key])
  out <- out[pathway_key %in% sig_keys]
  if (!nrow(out)) return(out)
  top_key <- unique(top_dt[, .(group_id, pathway_key)])
  top_key[, is_cluster_top_pathway := TRUE]
  out <- merge(out, top_key, by = c("group_id", "pathway_key"), all.x = TRUE, sort = FALSE)
  out[is.na(is_cluster_top_pathway), is_cluster_top_pathway := FALSE]
  out[, neglog10_padj := -log10(pmax(adjusted_p, 1e-300))]
  visual_group_order <- data.table::copy(group_order)
  data.table::setorder(
    visual_group_order,
    direction_rank,
    condition_rank,
    cell_rank,
    stress_rank,
    cell_line,
    stress_label,
    direction
  )
  row_levels <- character(0)
  for (gid in visual_group_order$group_id) {
    owned <- out[group_id == gid & is_cluster_top_pathway == TRUE]
    if (!nrow(owned)) next
    data.table::setorder(owned, -neglog10_padj, -overlap_hits, database, pathway)
    for (key in owned$pathway_key) {
      if (!(key %in% row_levels)) row_levels <- c(row_levels, key)
    }
  }
  remaining <- setdiff(unique(out$pathway_key), row_levels)
  if (length(remaining)) {
    remain_rank <- out[pathway_key %in% remaining, .(
      max_logp = max(neglog10_padj, na.rm = TRUE),
      max_overlap = max(overlap_hits, na.rm = TRUE)
    ), by = pathway_key]
    data.table::setorder(remain_rank, -max_logp, -max_overlap, pathway_key)
    row_levels <- c(row_levels, remain_rank$pathway_key)
  }
  out[, pathway_key := factor(pathway_key, levels = rev(row_levels))]
  out <- merge(out, group_order[, .(group_id, group_rank = seq_len(.N))], by = "group_id", all.x = TRUE, sort = FALSE)
  data.table::setorder(out, group_rank, -is_cluster_top_pathway, adjusted_p, -overlap_hits, database, pathway)
  out[]
}

#' Plot selected differential GRN pathway terms
#'
#' @param enrich_dt Data table from \code{read_diff_grn_pathway_enrichment()}.
#' @param out_file PDF output path.
#' @param tag Optional plot title prefix.
#' @param top_n Number of top pathways selected per comparison and direction.
#' @param padj_display_cut Adjusted p-value cutoff used before top-n selection.
#' @param condition_order Optional project-specific condition order.
#' @param cell_line_order Optional project-specific cell line order.
#' @param stress_order Optional project-specific stress order.
#' @param direction_order Direction facet/order.
#' @param size_metric Dot size metric. Use \code{"odds_ratio"} for the
#'   standard Enrichr odds ratio, \code{"combined_score"} for raw Enrichr
#'   combined score, \code{"log10_combined_score"} for
#'   \code{log10(combined_score + 1)}, or \code{"overlap_percent"} for
#'   percent overlap within the input gene set.
#' @param sig_color_break Color and border threshold on
#'   \code{-log10(adjusted p-value)}. Dots below this threshold have no visible
#'   border; dots at or above this threshold use a border from the same color
#'   scale as the fill.
#' @param panel_width Optional fixed dot-panel width in inches.
#' @param panel_height Optional fixed dot-panel height in inches.
#' @param font_family Optional PDF font family.
#'
#' @return Invisibly returns \code{out_file}.
#' @noRd
plot_diff_grn_pathway_dotplot <- function(enrich_dt,
                                          out_file,
                                          tag = NULL,
                                          top_n = 5L,
                                          padj_display_cut = 0.05,
                                          condition_order = NULL,
                                          cell_line_order = NULL,
                                          stress_order = NULL,
                                          direction_order = c("up", "down"),
                                          size_metric = c("odds_ratio", "combined_score", "log10_combined_score", "overlap_percent"),
                                          sig_color_break = 1.3,
                                          panel_width = NULL,
                                          panel_height = NULL,
                                          font_family = NULL) {
  size_metric <- match.arg(size_metric)
  sig_color_break <- as.numeric(sig_color_break)[1L]
  if (!is.finite(sig_color_break) || sig_color_break < 0) {
    stop("sig_color_break must be a non-negative finite number.", call. = FALSE)
  }
  dt <- select_diff_grn_pathway_terms(
    enrich_dt,
    top_n = top_n,
    padj_display_cut = padj_display_cut,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  if (!nrow(dt)) return(invisible(FALSE))
  group_order <- .diff_grn_pathway_group_order(
    dt,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  group_levels <- group_order$group_id
  names(group_levels) <- group_order$x_label
  dt[, group_id := factor(group_id, levels = group_levels)]
  dt[, direction := factor(direction, levels = direction_order)]
  size_name <- switch(
    size_metric,
    odds_ratio = "Odds ratio",
    combined_score = "Combined score",
    log10_combined_score = "log10 combined score",
    overlap_percent = "Overlap percent"
  )
  if (!"cluster_size" %in% names(dt)) dt[, cluster_size := NA_real_]
  dt[, plot_size := switch(
    size_metric,
    odds_ratio = data.table::fifelse(is.finite(odds_ratio), odds_ratio, NA_real_),
    combined_score = data.table::fifelse(is.finite(combined_score), combined_score, NA_real_),
    log10_combined_score = data.table::fifelse(is.finite(combined_score), log10(combined_score + 1), NA_real_),
    overlap_percent = data.table::fifelse(is.finite(overlap_hits) & is.finite(cluster_size) & cluster_size > 0, 100 * overlap_hits / cluster_size, NA_real_)
  )]
  dt[!is.finite(plot_size) & is.finite(odds_ratio), plot_size := odds_ratio]
  dt[!is.finite(plot_size) & is.finite(combined_score), plot_size := log10(combined_score + 1)]
  dt[!is.finite(plot_size) & is.finite(overlap_hits) & is.finite(cluster_size) & cluster_size > 0, plot_size := 100 * overlap_hits / cluster_size]
  size_cap <- suppressWarnings(stats::quantile(dt$plot_size[is.finite(dt$plot_size)], probs = 0.95, na.rm = TRUE, names = FALSE))
  if (is.finite(size_cap)) dt[, plot_size := pmin(plot_size, size_cap)]
  dt[, is_sig := is.finite(neglog10_padj) & neglog10_padj >= sig_color_break]
  max_val <- max(dt$neglog10_padj[is.finite(dt$neglog10_padj)], sig_color_break, na.rm = TRUE)
  break_val <- sig_color_break
  family <- if (is.null(font_family)) .diff_grn_pathway_font_family() else font_family
  title <- if (is.null(tag) || !nzchar(tag)) "Differential GRN pathway terms" else sprintf("%s pathway terms", tag)
  p_terms <- ggplot2::ggplot(dt, ggplot2::aes(x = group_id, y = pathway_key, size = plot_size)) +
    ggplot2::geom_point(ggplot2::aes(fill = neglog10_padj), shape = 21, color = "transparent", alpha = 0.85) +
    ggplot2::geom_point(
      data = dt[is_sig == TRUE],
      ggplot2::aes(fill = neglog10_padj, color = neglog10_padj),
      shape = 21,
      stroke = 0.75,
      alpha = 0.95
    ) +
    ggplot2::facet_grid(. ~ direction, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(labels = function(x) group_order$x_label[match(x, group_order$group_id)]) +
    ggplot2::scale_y_discrete(labels = function(x) {
      idx <- match(x, as.character(dt$pathway_key))
      dt$term_clean[idx]
    }) +
    ggplot2::scale_size_continuous(name = size_name, range = c(0.8, 3.2)) +
    ggplot2::scale_fill_gradientn(
      colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
      values = .diff_grn_pathway_rescale(c(0, break_val, seq(break_val + 1e-6, max_val, length.out = 6))),
      limits = c(0, max_val),
      oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
      name = "-log10 adjusted p-value"
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
      values = .diff_grn_pathway_rescale(c(0, break_val, seq(break_val + 1e-6, max_val, length.out = 6))),
      limits = c(0, max_val),
      oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
      guide = "none"
    ) +
    ggplot2::labs(
      title = title,
      x = "Condition",
      y = "Pathway term",
      caption = sprintf("Union of top %d pathways per condition/direction; all occurrences shown.", top_n)
    ) +
    ggplot2::theme_bw(base_size = 9, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9, face = "bold", family = family),
      plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family),
      axis.title = ggplot2::element_text(size = 9, face = "bold", family = family),
      axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family),
      legend.title = ggplot2::element_text(size = 9, face = "bold", family = family),
      legend.text = ggplot2::element_text(size = 9, face = "bold", family = family),
      strip.text = ggplot2::element_text(size = 9, face = "bold", family = family)
    )
  n_term <- data.table::uniqueN(dt$pathway_key)
  n_group <- data.table::uniqueN(dt$group_id)
  if (is.null(panel_width)) panel_width <- max(2.0, 0.095 * n_group)
  if (is.null(panel_height)) panel_height <- max(1.1, 0.17 * n_term)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  .diff_grn_pathway_save_fixed_panel_pdf(
    p_terms,
    out_file,
    panel_width = panel_width,
    panel_height = panel_height,
    family = family
  )
  invisible(out_file)
}

#' Plot a differential GRN pathway network
#'
#' @description
#' Draws a bipartite network between comparison/direction groups and selected
#' pathway terms. The pathway selection matches
#' \code{select_diff_grn_pathway_terms()}, but edges are drawn only for rows
#' passing \code{padj_edge_cut}; non-significant pathway occurrences are kept in
#' the selection table but omitted from the network.
#'
#' @param enrich_dt Data table from \code{read_diff_grn_pathway_enrichment()}.
#' @param out_file PDF output path.
#' @param tag Optional plot title prefix.
#' @param top_n Number of top pathways selected per comparison and direction.
#' @param padj_display_cut Adjusted p-value cutoff used before top-n selection.
#' @param padj_edge_cut Adjusted p-value cutoff for drawing network edges.
#' @param condition_order Optional project-specific condition order.
#' @param cell_line_order Optional project-specific cell line order.
#' @param stress_order Optional project-specific stress order.
#' @param direction_order Direction order.
#' @param font_family Optional PDF font family.
#'
#' @return Invisibly returns \code{out_file}, or FALSE if no network is drawn.
#' @noRd
plot_diff_grn_pathway_network <- function(enrich_dt,
                                          out_file,
                                          tag = NULL,
                                          top_n = 5L,
                                          padj_display_cut = 0.05,
                                          padj_edge_cut = 0.05,
                                          condition_order = NULL,
                                          cell_line_order = NULL,
                                          stress_order = NULL,
                                          direction_order = c("up", "down"),
                                          font_family = NULL) {
  need_pkgs <- c("igraph", "ggraph")
  if (!all(vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
    .log_warn("Skipping pathway network because igraph and/or ggraph is not installed.")
    return(invisible(FALSE))
  }
  dt <- select_diff_grn_pathway_terms(
    enrich_dt,
    top_n = top_n,
    padj_display_cut = padj_display_cut,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  if (!nrow(dt)) return(invisible(FALSE))
  dt <- dt[is.finite(adjusted_p) & adjusted_p <= as.numeric(padj_edge_cut)]
  if (!nrow(dt)) {
    .log_warn("Skipping pathway network because no selected pathway rows pass padj_edge_cut.")
    return(invisible(FALSE))
  }
  group_order <- .diff_grn_pathway_group_order(
    dt,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  group_order[, group_node := paste0("G:", group_id)]
  dt[, term_short := .diff_grn_pathway_shorten_term(pathway, max_chars = 42L)]
  dt[, group_node := paste0("G:", group_id)]
  dt[, pathway_node := paste0("P:", term_short)]
  dt[, neglog10_padj := -log10(pmax(adjusted_p, 1e-300))]
  edge_dt <- unique(dt[, .(
    from = group_node,
    to = pathway_node,
    overlap_hits,
    neglog10_padj
  )])
  group_nodes <- unique(merge(
    dt[, .(group_node = unique(group_node)), by = group_id],
    group_order[, .(group_id, x_label, direction, group_node)],
    by = c("group_id", "group_node"),
    all.x = TRUE,
    sort = FALSE
  )[, .(
    name = group_node,
    label = paste(x_label, direction, sep = " "),
    node_type = "comparison"
  )])
  pathway_nodes <- unique(dt[, .(
    name = pathway_node,
    label = term_short,
    node_type = "pathway"
  )])
  node_dt <- data.table::rbindlist(list(group_nodes, pathway_nodes), use.names = TRUE, fill = TRUE)
  graph_obj <- igraph::graph_from_data_frame(edge_dt, vertices = node_dt, directed = FALSE)
  layout_dt <- ggraph::create_layout(graph_obj, layout = "stress")
  family <- if (is.null(font_family)) .diff_grn_pathway_font_family() else font_family
  title <- if (is.null(tag) || !nzchar(tag)) "Differential GRN pathway network" else sprintf("%s pathway network", tag)
  p_net <- ggraph::ggraph(layout_dt) +
    ggraph::geom_edge_link(ggplot2::aes(width = overlap_hits, color = neglog10_padj), alpha = 0.55) +
    ggraph::scale_edge_width(range = c(0.2, 1.8), name = "Overlap genes") +
    ggraph::scale_edge_color_gradientn(
      colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
      values = .diff_grn_pathway_rescale(c(0, 1.3, seq(1.3 + 1e-6, max(edge_dt$neglog10_padj, 1.3, na.rm = TRUE), length.out = 6))),
      limits = c(0, max(edge_dt$neglog10_padj, 1.3, na.rm = TRUE)),
      oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
      name = "-log10 adjusted p-value",
      guide = ggraph::guide_edge_colourbar()
    ) +
    ggraph::geom_node_point(ggplot2::aes(shape = node_type, fill = node_type), size = 3.2, color = "black") +
    ggplot2::scale_shape_manual(values = c(comparison = 22, pathway = 21), name = NULL) +
    ggplot2::scale_fill_manual(values = c(comparison = "#f0c36d", pathway = "#7fb3d5"), name = NULL) +
    ggraph::geom_node_text(ggplot2::aes(label = label), size = 3, family = family, fontface = "bold", repel = TRUE, max.overlaps = Inf) +
    ggplot2::labs(
      title = title,
      caption = sprintf("Edges shown only for adjusted p-value <= %s.", padj_edge_cut)
    ) +
    ggraph::theme_graph(base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9, face = "bold", family = family),
      plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family),
      legend.title = ggplot2::element_text(size = 9, face = "bold", family = family),
      legend.text = ggplot2::element_text(size = 9, face = "bold", family = family)
    )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p_net, width = 15, height = 12)
  invisible(out_file)
}

#' Run differential GRN pathway analysis from filtered links
#'
#' @description
#' Convenience wrapper for the standard downstream workflow: Enrichr enrichment,
#' pathway table loading, union-top pathway selection, and dotplot output.
#'
#' @param filtered_dir Directory containing filtered differential-link files.
#' @param comparison_file Optional comparison design CSV.
#' @param out_dir Output directory for pathway enrichment tables and plots.
#' @param tag Output filename and plot title tag.
#' @param overwrite If TRUE, overwrite existing enrichment tables.
#' @param make_dotplot If TRUE, write the pathway term dotplot PDF.
#' @param dbs Character vector of Enrichr database names.
#' @param gene_col Gene column in filtered link files.
#' @param padj_cut Adjusted p-value cutoff for significant enrichment outputs.
#' @param min_genes Minimum distinct genes required before running Enrichr.
#' @param verbose Emit concise progress messages.
#' @param top_n Number of top pathways selected per comparison and direction.
#' @param padj_display_cut Adjusted p-value cutoff used before top-n selection.
#' @param condition_order Optional project-specific condition order.
#' @param cell_line_order Optional project-specific cell line order.
#' @param stress_order Optional project-specific stress order.
#' @param direction_order Direction facet/order.
#' @param size_metric Dot size metric passed to
#'   \code{plot_diff_grn_pathway_dotplot()}.
#' @param sig_color_break Color and border threshold passed to
#'   \code{plot_diff_grn_pathway_dotplot()}.
#' @param panel_width Optional fixed dot-panel width in inches.
#' @param panel_height Optional fixed dot-panel height in inches.
#' @param font_family Optional PDF font family.
#' @param ... Reserved for future extensions.
#'
#' @return Invisibly returns a list with output paths and selected terms.
#' @noRd
run_diff_grn_pathway_analysis <- function(filtered_dir,
                                          comparison_file = NULL,
                                          out_dir,
                                          tag = "diff_grn",
                                          overwrite = FALSE,
                                          make_dotplot = TRUE,
                                          dbs = .default_pathway_databases(),
                                          gene_col = "gene_key",
                                          padj_cut = 0.05,
                                          min_genes = 5L,
                                          verbose = TRUE,
                                          enrichr_sleep_time = 0,
                                          enrichr_cache_dir = NULL,
                                          pathway_backend = NULL,
                                          top_n = 5L,
                                          padj_display_cut = 0.05,
                                          condition_order = NULL,
                                          cell_line_order = NULL,
                                          stress_order = NULL,
                                          direction_order = c("up", "down"),
                                          size_metric = c("odds_ratio", "combined_score", "log10_combined_score", "overlap_percent"),
                                          sig_color_break = 1.3,
                                          panel_width = NULL,
                                          panel_height = NULL,
                                          font_family = NULL,
                                          ...) {
  size_metric <- match.arg(size_metric)
  enrichment_files <- run_diff_grn_pathway_enrichment(
    filtered_dir = filtered_dir,
    out_dir = out_dir,
    dbs = dbs,
    gene_col = gene_col,
    padj_cut = padj_cut,
    min_genes = min_genes,
    overwrite = overwrite,
    enrichr_sleep_time = enrichr_sleep_time,
    enrichr_cache_dir = enrichr_cache_dir,
    pathway_backend = pathway_backend,
    verbose = verbose
  )
  enrich_dt <- read_diff_grn_pathway_enrichment(
    pathway_dir = out_dir,
    comparison_file = comparison_file
  )
  selected_terms <- select_diff_grn_pathway_terms(
    enrich_dt = enrich_dt,
    top_n = top_n,
    padj_display_cut = padj_display_cut,
    condition_order = condition_order,
    cell_line_order = cell_line_order,
    stress_order = stress_order,
    direction_order = direction_order
  )
  term_path <- file.path(out_dir, sprintf("%s_pathway_plot_terms.csv", tag))
  data.table::fwrite(selected_terms, term_path)
  dotplot_path <- file.path(out_dir, "pdf", sprintf("%s_pathway_term_dotplot.pdf", tag))
  if (isTRUE(make_dotplot)) {
    plot_diff_grn_pathway_dotplot(
      enrich_dt = enrich_dt,
      out_file = dotplot_path,
      tag = tag,
      top_n = top_n,
      padj_display_cut = padj_display_cut,
      condition_order = condition_order,
      cell_line_order = cell_line_order,
      stress_order = stress_order,
      direction_order = direction_order,
      size_metric = size_metric,
      sig_color_break = sig_color_break,
      panel_width = panel_width,
      panel_height = panel_height,
      font_family = font_family
    )
  }
  invisible(list(
    enrichment_files = enrichment_files,
    selected_terms_file = term_path,
    dotplot_file = if (isTRUE(make_dotplot)) dotplot_path else NA_character_,
    selected_terms = selected_terms
  ))
}

.diff_grn_pathway_comparison_parts <- function(comparison_id) {
  parts <- strsplit(comparison_id, "_vs_", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) {
    return(list(cond1_label = NA_character_, cond2_label = NA_character_, cell_line = NA_character_))
  }
  cell_line <- sub("_.*$", "", parts[[1L]])
  list(cond1_label = parts[[1L]], cond2_label = parts[[2L]], cell_line = cell_line)
}

.diff_grn_pathway_clean_term <- function(term) {
  term <- gsub("_", " ", as.character(term))
  term <- gsub("\\s*\\([^)]*\\)$", "", term)
  term <- gsub("^(HALLMARK|GOBP|GOMF|GOCC|KEGG|REACTOME|BIOCARTA|PID|WP) ", "", term)
  trimws(term)
}

.diff_grn_pathway_shorten_term <- function(term, max_chars = 54L) {
  term <- .diff_grn_pathway_clean_term(term)
  term <- gsub("\\s+", " ", term)
  too_long <- nchar(term) > max_chars
  term[too_long] <- paste0(substr(term[too_long], 1L, max_chars - 3L), "...")
  term
}

.diff_grn_pathway_parse_overlap <- function(overlap) {
  vals <- suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))
  vals[is.na(vals)] <- NA_integer_
  vals
}

.diff_grn_pathway_group_order <- function(enrich_dt,
                                          condition_order = NULL,
                                          cell_line_order = NULL,
                                          stress_order = NULL,
                                          direction_order = c("up", "down")) {
  ord <- unique(enrich_dt[, .(group_id, comparison_id, cond1_label, cell_line, stress_label, x_label, direction)])
  if (!is.null(condition_order)) {
    if (is.data.frame(condition_order)) {
      condition_order <- data.table::as.data.table(condition_order)
      if (!("condition_rank" %in% names(condition_order))) condition_order[, condition_rank := seq_len(.N)]
      by_cols <- intersect(c("group_id", "comparison_id", "cond1_label", "direction"), names(condition_order))
      if (length(by_cols)) {
        ord <- merge(ord, condition_order[, c(by_cols, "condition_rank"), with = FALSE], by = by_cols, all.x = TRUE, sort = FALSE)
      }
    } else {
      condition_order <- as.character(condition_order)
      ord[, condition_rank := match(comparison_id, condition_order)]
      ord[is.na(condition_rank), condition_rank := match(cond1_label, condition_order)]
      ord[is.na(condition_rank), condition_rank := match(group_id, condition_order)]
    }
  }
  if (!("condition_rank" %in% names(ord))) ord[, condition_rank := NA_integer_]
  if (!is.null(cell_line_order)) {
    ord[, cell_rank := match(cell_line, cell_line_order)]
  } else {
    ord[, cell_rank := data.table::frank(cell_line, ties.method = "dense")]
  }
  ord[is.na(cell_rank), cell_rank := max(cell_rank, na.rm = TRUE) + 1L]
  if (!is.null(stress_order)) {
    ord[, stress_rank := match(stress_label, stress_order)]
  } else {
    ord[, stress_rank := data.table::frank(stress_label, ties.method = "dense")]
  }
  ord[is.na(stress_rank), stress_rank := max(stress_rank, na.rm = TRUE) + 1L]
  ord[, direction_rank := match(direction, direction_order)]
  ord[is.na(direction_rank), direction_rank := length(direction_order) + 1L]
  data.table::setorder(ord, condition_rank, cell_rank, stress_rank, direction_rank, cell_line, stress_label, direction)
  ord
}

.diff_grn_pathway_font_family <- function() {
  pdf_preferred <- c("Helvetica", "ArialMT", "NimbusSan", "URWHelvetica", "sans")
  pdf_fonts <- names(grDevices::pdfFonts())
  pdf_hit <- pdf_preferred[pdf_preferred %in% pdf_fonts]
  if (length(pdf_hit)) return(pdf_hit[[1L]])
  "sans"
}

.diff_grn_pathway_rescale <- function(x) {
  rng <- range(x, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

.diff_grn_pathway_save_fixed_panel_pdf <- function(plot_obj, pdf_path, panel_width, panel_height, family = "sans") {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(pdf_path)) unlink(pdf_path)
  gt <- ggplot2::ggplotGrob(plot_obj)
  panel_idx <- grepl("^panel", gt$layout$name)
  panel_cols <- unique(gt$layout$l[panel_idx])
  panel_rows <- unique(gt$layout$t[panel_idx])
  gt$widths[panel_cols] <- grid::unit(panel_width, "in")
  gt$heights[panel_rows] <- grid::unit(panel_height, "in")
  width <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = family)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(gt)
  invisible(pdf_path)
}

utils::globalVariables(c(
  "adjusted_p",
  "cell_line",
  "cell_rank",
  "combined_score",
  "comparison_id",
  "cond1_label",
  "cond2_label",
  "condition_rank",
  "database",
  "direction",
  "direction_rank",
  "group_node",
  "group_id",
  "group_rank",
  "is_cluster_top_pathway",
  "label",
  "max_logp",
  "max_overlap",
  "name",
  "neglog10_padj",
  "node_type",
  "odds_ratio",
  "overlap",
  "overlap_hits",
  "pathway",
  "pathway_key",
  "pathway_node",
  "plot_size",
  "stress_label",
  "stress_rank",
  "term_short",
  "term_clean",
  "x_label"
))
