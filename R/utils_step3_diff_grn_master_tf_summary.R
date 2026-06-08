# File: utils_step3_diff_grn_master_tf_summary.R
# Author: Yaoxiang Li
# Created: 2026-05-01
#
# Purpose:
# Package-level downstream master TF summary utilities for Step3 differential
# GRN filtered links.

#' Report master TF activity
#'
#' @description
#' Builds standard Module 3 master-TF reports from filtered differential links.
#' User-facing plots are written directly under `output_dir`; supporting
#' per-comparison summary tables are written under `output_dir/master_tf_tables`.
#'
#' @param filtered_dir Directory containing filtered differential-link CSVs.
#' @param output_dir Directory where user-facing report outputs are written.
#' @param overwrite If TRUE, overwrite existing outputs.
#' @param connectivity_min_degree Deprecated. Min-distance heatmaps include
#'   source rows with at least one outgoing TF-to-TF link and target columns with
#'   at least one incoming TF-to-TF link.
#' @param waterfall_min_abs_net Minimum absolute net unique target count used to
#'   display TFs in the waterfall plot.
#' @param verbose Emit concise progress messages.
#'
#' @return Invisibly returns a data.table manifest of written/skipped outputs.
#' @noRd
report_master_tfs <- function(filtered_dir,
                              output_dir,
                              overwrite = TRUE,
                              connectivity_min_degree = 5L,
                              waterfall_min_abs_net = 20L,
                              verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(output_dir, "master_tf_tables")
  manifest <- run_diff_grn_master_tf_summary(
    filtered_dir = filtered_dir,
    out_dir = table_dir,
    overwrite = overwrite,
    connectivity_min_degree = connectivity_min_degree,
    waterfall_min_abs_net = waterfall_min_abs_net,
    verbose = verbose
  )
  if (nrow(manifest)) {
    path_cols <- intersect(
      c("summary_pdf", "waterfall_pdf", "mindist_pdf", "composite_pdf"),
      names(manifest)
    )
    for (col in path_cols) {
      paths <- manifest[[col]]
      paths <- paths[!is.na(paths) & nzchar(paths) & file.exists(paths)]
      for (src in paths) {
        dest <- file.path(output_dir, basename(src))
        if (!identical(normalizePath(src, winslash = "/", mustWork = FALSE), normalizePath(dest, winslash = "/", mustWork = FALSE))) {
          file.copy(src, dest, overwrite = TRUE)
        }
      }
    }
  }
  invisible(manifest)
}

#' Run master TF summaries from filtered differential GRN links
#'
#' @description
#' Reads paired \code{*_filtered_links_up.csv} and
#' \code{*_filtered_links_down.csv} files, writes per-comparison master TF
#' summary CSVs and plots. Connectivity outputs are written as two separate
#' PDFs: min-distance TF-to-TF connectivity and composite TF-to-TF connectivity.
#'
#' @param filtered_dir Directory containing filtered differential-link CSVs.
#' @param out_dir Output directory for master TF summary files.
#' @param overwrite If TRUE, overwrite existing outputs.
#' @param connectivity_min_degree Deprecated. Min-distance heatmaps include
#'   source rows with at least one outgoing TF-to-TF link and target columns with
#'   at least one incoming TF-to-TF link.
#' @param waterfall_min_abs_net Minimum absolute net unique target count used to
#'   display TFs in the waterfall plot.
#' @param verbose Emit concise progress messages.
#'
#' @return Invisibly returns a data.table manifest of written/skipped outputs.
#' @noRd
run_diff_grn_master_tf_summary <- function(filtered_dir,
                                           out_dir,
                                           overwrite = TRUE,
                                           connectivity_min_degree = 5L,
                                           waterfall_min_abs_net = 20L,
                                           verbose = TRUE) {
  if (!dir.exists(filtered_dir)) .log_abort("`filtered_dir` not found: {filtered_dir}")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  filtered_files <- sort(list.files(filtered_dir, "_filtered_links_(up|down)\\.csv$", full.names = TRUE))
  if (!length(filtered_files)) {
    .log_warn("No filtered differential-link files found in: {filtered_dir}")
    return(invisible(data.table::data.table()))
  }
  base_ids <- unique(sub("_(up|down)\\.csv$", "", basename(filtered_files)))
  manifest <- vector("list", length(base_ids))
  for (i in seq_along(base_ids)) {
    bid <- base_ids[[i]]
    pattern <- paste0("^", .diff_grn_master_escape_regex(bid), "_(up|down)\\.csv$")
    files <- filtered_files[grepl(pattern, basename(filtered_files))]
    contrast_id <- sub("_filtered_links.*$", "", bid)
    out_csv <- file.path(out_dir, paste0(contrast_id, "_master_tf_summary.csv"))
    out_pdf <- file.path(out_dir, paste0(contrast_id, "_master_tf_summary.pdf"))
    out_pdf_waterfall <- file.path(out_dir, paste0(contrast_id, "_tf_target_direction_waterfall.pdf"))
    out_pdf_mindist <- file.path(out_dir, paste0(contrast_id, "_tf_tf_connectivity_heatmap.pdf"))
    out_pdf_composite <- file.path(out_dir, paste0(contrast_id, "_tf_tf_composite_connectivity_heatmap.pdf"))
    expected <- c(out_csv, out_pdf, out_pdf_waterfall, out_pdf_mindist, out_pdf_composite)
    if (!isTRUE(overwrite) && all(file.exists(expected))) {
      if (isTRUE(verbose)) .log_inform("Skipping existing master TF summary: {contrast_id}")
      manifest[[i]] <- data.table::data.table(contrast_id = contrast_id, status = "exists")
      next
    }
    links_dt <- data.table::rbindlist(lapply(files, function(path) {
      dt <- data.table::fread(path, showProgress = FALSE)
      dt[, direction_group := if (grepl("_up\\.csv$", basename(path), ignore.case = TRUE)) "Up" else "Down"]
      dt
    }), use.names = TRUE, fill = TRUE)
    if (!nrow(links_dt)) {
      manifest[[i]] <- data.table::data.table(contrast_id = contrast_id, status = "empty_links")
      next
    }
    parts <- .diff_grn_master_contrast_parts(contrast_id)
    summary_dt <- summarize_diff_grn_master_tf_links(links_dt, cond1 = parts$cond1, cond2 = parts$cond2)
    if (!nrow(summary_dt)) {
      manifest[[i]] <- data.table::data.table(contrast_id = contrast_id, status = "empty_summary")
      next
    }
    data.table::fwrite(summary_dt, out_csv)
    plot_diff_grn_master_tf_summary(
      summary_dt = summary_dt,
      links_dt = links_dt,
      out_pdf = out_pdf,
      title_text = paste0("TF hubs (delta fp_score) - ", parts$label),
      cond1_label = parts$cond1,
      cond2_label = parts$cond2,
      connectivity_min_degree = connectivity_min_degree,
      waterfall_min_abs_net = waterfall_min_abs_net
    )
    if (isTRUE(verbose)) .log_inform("Saved master TF summary: {out_csv}")
    manifest[[i]] <- data.table::data.table(
      contrast_id = contrast_id,
      status = "written",
      n_links = nrow(links_dt),
      n_tfs = nrow(summary_dt),
      summary_csv = out_csv,
      summary_pdf = out_pdf,
      waterfall_pdf = out_pdf_waterfall,
      mindist_pdf = out_pdf_mindist,
      composite_pdf = out_pdf_composite
    )
  }
  out <- data.table::rbindlist(manifest, use.names = TRUE, fill = TRUE)
  invisible(out)
}

#' Summarize filtered links into master TF metrics
#'
#' @param links_dt Filtered differential-link table.
#' @param cond1 Optional condition-1 label used to find TF expression columns.
#' @param cond2 Optional condition-2 label used to find TF expression columns.
#'
#' @return A data.table with one row per TF.
#' @noRd
summarize_diff_grn_master_tf_links <- function(links_dt, cond1 = NULL, cond2 = NULL) {
  dt <- data.table::as.data.table(links_dt)
  tf_col <- if ("tf" %in% names(dt)) "tf" else if ("TF" %in% names(dt)) "TF" else NULL
  gene_col <- if ("gene_key" %in% names(dt)) "gene_key" else if ("gene" %in% names(dt)) "gene" else if ("target_gene" %in% names(dt)) "target_gene" else NULL
  delta_col <- if ("delta_fp_score" %in% names(dt)) "delta_fp_score" else if ("delta_fp_bed_score" %in% names(dt)) "delta_fp_bed_score" else if ("delta_fp" %in% names(dt)) "delta_fp" else NULL
  if (is.null(tf_col) || is.null(gene_col) || is.null(delta_col)) return(data.table::data.table())
  tf_expr_max <- NA_real_
  if (is.character(cond1) && nzchar(cond1) && is.character(cond2) && nzchar(cond2)) {
    tf_c1 <- paste0("tf_expr_", cond1)
    tf_c2 <- paste0("tf_expr_", cond2)
    if (all(c(tf_c1, tf_c2) %in% names(dt))) {
      tf_expr_max <- pmax(
        suppressWarnings(as.numeric(dt[[tf_c1]])),
        suppressWarnings(as.numeric(dt[[tf_c2]])),
        na.rm = TRUE
      )
    }
  }
  if (all(is.na(tf_expr_max)) && "tf_expr_max" %in% names(dt)) {
    tf_expr_max <- suppressWarnings(as.numeric(dt$tf_expr_max))
  }
  if (all(is.na(tf_expr_max)) && all(c("tf_expr_ctrl", "tf_expr_str") %in% names(dt))) {
    tf_expr_max <- pmax(
      suppressWarnings(as.numeric(dt$tf_expr_ctrl)),
      suppressWarnings(as.numeric(dt$tf_expr_str)),
      na.rm = TRUE
    )
  }
  tf_log2_fc <- if ("log2FC_tf_expr" %in% names(dt)) {
    suppressWarnings(as.numeric(dt$log2FC_tf_expr))
  } else if ("delta_tf_expr" %in% names(dt)) {
    suppressWarnings(as.numeric(dt$delta_tf_expr))
  } else if ("tf_log2_fc" %in% names(dt)) {
    suppressWarnings(as.numeric(dt$tf_log2_fc))
  } else {
    NA_real_
  }
  work <- data.table::data.table(
    TF = as.character(dt[[tf_col]]),
    gene = as.character(dt[[gene_col]]),
    delta = suppressWarnings(as.numeric(dt[[delta_col]])),
    tf_expr_max = tf_expr_max,
    tf_log2_fc = tf_log2_fc
  )
  work <- work[is.finite(delta) & !is.na(TF) & nzchar(TF)]
  if (!nrow(work)) return(data.table::data.table())
  work[, reg_sign := sign(delta)]
  work[, delta_abs := abs(delta)]
  out <- work[, .(
    tf_delta_sum = sum(delta, na.rm = TRUE),
    tf_delta_sum_activate = sum(ifelse(reg_sign > 0, delta, 0), na.rm = TRUE),
    tf_delta_sum_repress = sum(ifelse(reg_sign < 0, delta, 0), na.rm = TRUE),
    tf_delta_sum_abs = sum(delta_abs, na.rm = TRUE),
    tf_delta_sum_abs_activate = sum(ifelse(reg_sign > 0, delta_abs, 0), na.rm = TRUE),
    tf_delta_sum_abs_repress = sum(ifelse(reg_sign < 0, delta_abs, 0), na.rm = TRUE),
    tf_n_links = .N,
    tf_n_links_activate = sum(reg_sign > 0, na.rm = TRUE),
    tf_n_links_repress = sum(reg_sign < 0, na.rm = TRUE),
    tf_n_target_genes = data.table::uniqueN(gene),
    tf_n_target_genes_activate = data.table::uniqueN(gene[reg_sign > 0]),
    tf_n_target_genes_repress = data.table::uniqueN(gene[reg_sign < 0]),
    tf_expr_max = suppressWarnings(max(tf_expr_max, na.rm = TRUE)),
    tf_log2_fc = stats::median(tf_log2_fc, na.rm = TRUE),
    tf_hits_hub = NA_real_
  ), by = TF]
  out[!is.finite(tf_expr_max), tf_expr_max := NA_real_]
  out[!is.finite(tf_log2_fc), tf_log2_fc := NA_real_]
  out[, topic := 1L]
  out[, topic_mean_abs_delta := mean(tf_delta_sum_abs, na.rm = TRUE)]
  out[, topic_n_TFs := .N]
  out[, topic_n_links := sum(tf_n_links, na.rm = TRUE)]
  out[, topic_n_target_genes := data.table::uniqueN(work$gene)]
  out[, topic_rank := 1L]
  data.table::setorder(out, -tf_delta_sum_abs, TF)
  out[]
}

#' Plot master TF summary outputs
#'
#' @param summary_dt Master TF summary table.
#' @param links_dt Filtered links used for waterfall and connectivity plots.
#' @param out_pdf Master TF summary PDF output path.
#' @param title_text Plot title.
#' @param cond1_label Condition-1 label.
#' @param cond2_label Condition-2 label.
#' @param connectivity_min_degree Deprecated. Min-distance heatmaps include
#'   source rows with at least one outgoing TF-to-TF link and target columns with
#'   at least one incoming TF-to-TF link.
#' @param waterfall_min_abs_net Minimum absolute net target count for waterfall.
#' @param label_top_n Number of highest-priority TFs to annotate in the scatter
#'   plot.
#' @param label_top_fraction Minimum fraction of TFs to annotate in the scatter
#'   plot, after ranking by link count and absolute delta score.
#'
#' @return Invisibly returns \code{out_pdf}.
#' @noRd
plot_diff_grn_master_tf_summary <- function(summary_dt,
                                            links_dt = NULL,
                                            out_pdf,
                                            title_text = "TF hubs",
                                            cond1_label = NULL,
                                            cond2_label = NULL,
                                            connectivity_min_degree = 5L,
                                            waterfall_min_abs_net = 20L,
                                            label_top_n = 45L,
                                            label_top_fraction = 0.8) {
  for (pkg in c("ggrepel", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE)) .log_abort("Package `{pkg}` is required for master TF summary plots.")
  }
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  summary_dt <- data.table::as.data.table(summary_dt)
  need <- c("TF", "tf_expr_max", "tf_log2_fc", "tf_n_links", "tf_delta_sum", "tf_delta_sum_abs")
  if (!all(need %in% names(summary_dt))) .log_abort("Master TF summary missing required columns.")
  if (!("tf_n_target_genes" %in% names(summary_dt))) {
    summary_dt[, tf_n_target_genes := tf_n_links]
  }
  c1 <- if (is.null(cond1_label) || !nzchar(cond1_label)) "cond1" else cond1_label
  c2 <- if (is.null(cond2_label) || !nzchar(cond2_label)) "cond2" else cond2_label
  tf_dt <- summary_dt[, .(
    tf_links = sum(tf_n_links, na.rm = TRUE),
    tf_target_genes = sum(tf_n_target_genes, na.rm = TRUE),
    tf_sum_delta = sum(tf_delta_sum, na.rm = TRUE),
    tf_sum_abs_delta = sum(tf_delta_sum_abs, na.rm = TRUE),
    tf_expr_max = max(tf_expr_max, na.rm = TRUE),
    tf_log2_fc_med = stats::median(tf_log2_fc, na.rm = TRUE)
  ), by = TF]
  tf_dt[!is.finite(tf_expr_max), tf_expr_max := NA_real_]
  tf_dt[, tf_log2_expr_max := log2(pmax(tf_expr_max, 1e-9))]
  tf_dt[, x_sum := tf_sum_delta]
  tf_dt[, x_sum_signed_log10 := sign(x_sum) * log10(abs(x_sum) + 1)]
  tf_dt[, log2_fc_z_c := .diff_grn_master_clamp(.diff_grn_master_robust_z(tf_log2_fc_med), -2, 2)]
  label_top_n <- max(0L, as.integer(label_top_n)[1], ceiling(nrow(tf_dt) * as.numeric(label_top_fraction)[1]), na.rm = TRUE)
  if (label_top_n > 0L) {
    lab_dt <- tf_dt[order(-tf_target_genes, -tf_sum_abs_delta)][seq_len(min(label_top_n, .N))]
    lab_dt <- unique(data.table::rbindlist(list(
      lab_dt,
      tf_dt[order(-tf_sum_abs_delta, -tf_target_genes)][seq_len(min(max(5L, floor(label_top_n / 2L)), .N))]
    ), use.names = TRUE), by = "TF")
  } else {
    lab_dt <- tf_dt[0]
  }
  family <- .diff_grn_master_font_family()
  .scatter_plot <- function(x_col, x_label, y_log10 = FALSE, panel_title) {
    p <- ggplot2::ggplot(tf_dt, ggplot2::aes(x = .data[[x_col]], y = tf_target_genes)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(size = tf_log2_expr_max, color = log2_fc_z_c), shape = 16, alpha = 0.8) +
      ggrepel::geom_text_repel(data = lab_dt, ggplot2::aes(label = TF), size = 2.2, family = family, fontface = "bold", color = "black", segment.color = "black", max.overlaps = Inf, box.padding = 0.18, point.padding = 0.08, force = 1.2, force_pull = 0.35, min.segment.length = 0.03, max.time = 3) +
      ggplot2::scale_size_continuous(range = c(1, 4), name = "expr max (log2)") +
      ggplot2::scale_color_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "TF log2FC (z)") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.10))) +
      ggplot2::labs(
        title = panel_title,
        x = x_label,
        y = "number of unique target genes per TF",
        caption = paste0("delta fp_score = fp_score(", c1, ") - fp_score(", c2, ").")
      ) +
      ggplot2::theme_classic(base_size = 9, base_family = family) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        plot.title = ggplot2::element_text(size = 9, hjust = 0.5, face = "bold", family = family, color = "black"),
        axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
      )
    if (isTRUE(y_log10)) {
      p <- p +
        ggplot2::scale_y_log10(
          breaks = scales::breaks_log(n = 6),
          labels = scales::label_number(accuracy = 1),
          expand = ggplot2::expansion(mult = c(0.03, 0.08))
        )
    } else {
      p <- p +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.08)))
    }
    p
  }
  p_raw <- .scatter_plot(
    x_col = "x_sum",
    x_label = "sum delta fp_score per TF",
    y_log10 = FALSE,
    panel_title = paste0(title_text, " | raw axes")
  )
  p_log <- .scatter_plot(
    x_col = "x_sum_signed_log10",
    x_label = "signed log10(abs(sum delta fp_score) + 1) per TF",
    y_log10 = TRUE,
    panel_title = paste0(title_text, " | log axes")
  )
  grDevices::cairo_pdf(out_pdf, width = 20, height = 9.5, family = family)
  tryCatch({
    grid::grid.newpage()
    print(p_raw, vp = grid::viewport(x = 0.25, y = 0.5, width = 0.5, height = 1))
    print(p_log, vp = grid::viewport(x = 0.75, y = 0.5, width = 0.5, height = 1))
  }, finally = {
    grDevices::dev.off()
  })
  if (is.data.frame(links_dt) && nrow(links_dt)) {
    out_base <- sub("\\.pdf$", "", basename(out_pdf))
    out_base <- sub("_master_tf_summary$", "", out_base)
    plot_diff_grn_master_tf_waterfall(
      links_dt,
      out_pdf = file.path(dirname(out_pdf), paste0(out_base, "_tf_target_direction_waterfall.pdf")),
      title_text = title_text,
      waterfall_min_abs_net = waterfall_min_abs_net
    )
    plot_diff_grn_master_tf_connectivity(
      links_dt,
      out_pdf_mindist = file.path(dirname(out_pdf), paste0(out_base, "_tf_tf_connectivity_heatmap.pdf")),
      out_pdf_composite = file.path(dirname(out_pdf), paste0(out_base, "_tf_tf_composite_connectivity_heatmap.pdf")),
      title_text = title_text,
      connectivity_min_degree = connectivity_min_degree
    )
  }
  invisible(out_pdf)
}

.diff_grn_master_read_tf_cluster_annotation <- function(path = NULL) {
  path <- if (is.null(path) || !nzchar(path)) {
    .craftgrn_getenv("DIFF_GRN_TF_CLUSTER_ANNOTATION", unset = "")
  } else {
    as.character(path)[1]
  }
  if (!nzchar(path) || !file.exists(path)) {
    return(data.table::data.table(tf = character(), cluster = character()))
  }
  dt <- data.table::fread(path, showProgress = FALSE)
  tf_col <- intersect(c("tf", "TF", "Name"), names(dt))[1]
  cluster_col <- intersect(c("craftgrn_cluster", "jaspar_cluster", "CRAFT_GRN_cluster_DBD", "dbd_name", "cluster"), names(dt))[1]
  if (is.na(tf_col) || is.na(cluster_col)) {
    .log_warn("TF cluster annotation file lacks a TF or cluster column: {path}")
    return(data.table::data.table(tf = character(), cluster = character()))
  }
  out <- unique(dt[, .(
    tf = toupper(trimws(as.character(get(tf_col)))),
    cluster = trimws(as.character(get(cluster_col)))
  )])
  out[is.na(cluster) | !nzchar(cluster), cluster := "Unassigned"]
  out[!is.na(tf) & nzchar(tf)]
}

#' Plot TF target-direction waterfall
#'
#' @param links_dt Filtered links table.
#' @param out_pdf PDF output path.
#' @param title_text Plot title prefix.
#' @param waterfall_min_abs_net Minimum absolute net target count for display.
#'
#' @return Invisibly returns \code{out_pdf}.
#' @noRd
plot_diff_grn_master_tf_waterfall <- function(links_dt,
                                              out_pdf,
                                              title_text = "TF hubs",
                                              waterfall_min_abs_net = 20L) {
  link_dt <- .diff_grn_master_normalize_links(links_dt)
  if (!nrow(link_dt)) return(invisible(FALSE))
  gene_dt <- unique(link_dt[, .(TF, gene_key, direction, target_type)])
  tf_net <- gene_dt[, .(
    n_up_total = sum(direction == "Up", na.rm = TRUE),
    n_down_total = sum(direction == "Down", na.rm = TRUE)
  ), by = TF][, net_n := n_up_total - n_down_total]
  tf_order <- tf_net[abs(net_n) > as.numeric(waterfall_min_abs_net)][order(net_n, decreasing = TRUE)]$TF
  if (!length(tf_order)) {
    grDevices::pdf(out_pdf, width = 10, height = 6.5, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::plot.new()
    graphics::title(main = paste0(title_text, " | TF target-direction waterfall"))
    graphics::text(0.5, 0.5, paste0("No TFs pass cutoff: abs(unique(Up targets)-unique(Down targets)) > ", waterfall_min_abs_net))
    return(invisible(out_pdf))
  }
  bar_dt <- gene_dt[TF %in% tf_order, .(n = ifelse(.BY$direction == "Up", .N, -.N)), by = .(TF, target_type, direction)]
  bar_dt[, signed_log2_n := sign(n) * log2(abs(n) + 1)]
  bar_dt[, fill_group := paste0(direction, " | ", target_type)]
  bar_dt[, TF := factor(TF, levels = rev(tf_order))]
  family <- .diff_grn_master_font_family()
  panel_height <- max(1.6, 0.18 * length(tf_order))
  panel_width <- max(3.2, 0.55 * max(abs(bar_dt$signed_log2_n), na.rm = TRUE))
  p <- ggplot2::ggplot(bar_dt, ggplot2::aes(x = signed_log2_n, y = TF, fill = fill_group)) +
    ggplot2::geom_col(width = 0.62, color = "grey30", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c("Up | TF target" = "#B2182B", "Up | Gene target" = "#EF8A62", "Down | TF target" = "#2166AC", "Down | Gene target" = "#67A9CF")) +
    ggplot2::labs(title = paste0(sub("^TF hubs \\(delta fp_score\\) - ", "", title_text), "\nTF target-gene waterfall"), x = "signed log2(unique-gene count + 1)", y = "TF", fill = "Direction | Target", caption = paste0("x = sign(count) x log2(abs(count) + 1).\nTFs shown only if abs(unique(Up targets) - unique(Down targets)) > ", waterfall_min_abs_net, ".")) +
    ggplot2::theme_classic(base_size = 9, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      plot.title = ggplot2::element_text(size = 9, hjust = 0.5, face = "bold", family = family, color = "black"),
      axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
    )
  .diff_grn_master_save_fixed_panel_pdf(
    p,
    out_pdf,
    panel_width = panel_width,
    panel_height = panel_height,
    family = family
  )
  invisible(out_pdf)
}

#' Plot TF-to-TF connectivity heatmaps
#'
#' @param links_dt Filtered links table.
#' @param out_pdf_mindist Min-distance connectivity PDF path.
#' @param out_pdf_composite Composite connectivity PDF path.
#' @param title_text Plot title prefix.
#' @param connectivity_min_degree Deprecated. Min-distance heatmaps include
#'   source rows with at least one outgoing TF-to-TF link and target columns with
#'   at least one incoming TF-to-TF link.
#'
#' @return Invisibly returns a character vector of output PDF paths.
#' @noRd
plot_diff_grn_master_tf_connectivity <- function(links_dt,
                                                 out_pdf_mindist,
                                                 out_pdf_composite,
                                                 title_text = "TF hubs",
                                                 connectivity_min_degree = 5L) {
  link_dt <- .diff_grn_master_normalize_links(links_dt)
  tf2tf <- link_dt[target_type == "TF target", .(TF_src = TF, TF_tgt = gene_key)]
  if (!nrow(tf2tf)) {
    .diff_grn_master_empty_pdf(out_pdf_mindist, paste0(title_text, " | TF-to-TF connectivity heatmap"), "No TF->TF links found")
    .diff_grn_master_empty_pdf(out_pdf_composite, paste0(title_text, " | TF-to-TF composite connectivity"), "No TF->TF links found")
    return(invisible(c(out_pdf_mindist, out_pdf_composite)))
  }
  conn <- tf2tf[, .N, by = .(TF_src, TF_tgt)]
  tf_levels <- sort(unique(c(conn$TF_src, conn$TF_tgt)))
  mat <- matrix(0, nrow = length(tf_levels), ncol = length(tf_levels), dimnames = list(tf_levels, tf_levels))
  for (ii in seq_len(nrow(conn))) mat[conn$TF_src[ii], conn$TF_tgt[ii]] <- conn$N[ii]
  if (sum(rowSums(mat > 0, na.rm = TRUE) > 0) >= 1L && sum(colSums(mat > 0, na.rm = TRUE) > 0) >= 1L) {
    mindist <- .diff_grn_master_mindist_score(mat)
    plot_diff_grn_master_tf_connectivity_matrix(
      mindist,
      out_pdf_mindist,
      paste0(sub("^TF hubs \\(delta fp_score\\) - ", "", title_text), " | TF-to-TF min-distance connectivity"),
      "connectivity score",
      fill_limits = c(0, 1),
      fill_breaks = c(0, 0.5, 1),
      subtitle = .diff_grn_master_wrap_text(paste(
        "TF-to-TF min-distance connectivity (directed).",
        "Score(i,j)=1/# directed TF->TF steps; direct=1; two-step=0.5; unreachable=0."
      ), width = 96),
      caption = "TF inclusion cutoff: rows have >=1 outgoing TF-to-TF link; columns have >=1 incoming TF-to-TF link.",
      cell_size = 0.13,
      font_size = 8,
      row_annotation = .diff_grn_master_read_tf_cluster_annotation()
    )
  } else {
    .diff_grn_master_empty_pdf(out_pdf_mindist, paste0(title_text, " | TF-to-TF min-distance connectivity"), "No TFs pass directed TF-to-TF inclusion cutoff")
  }
  adj_dir_comp <- mat > 0
  diag(adj_dir_comp) <- FALSE
  comp_nodes <- intersect(
    rownames(mat)[rowSums(adj_dir_comp, na.rm = TRUE) > 0],
    colnames(mat)[colSums(adj_dir_comp, na.rm = TRUE) > 0]
  )
  nodes <- comp_nodes
  targets_by_tf <- lapply(nodes, function(tf_i) unique(link_dt[TF == tf_i, gene_key]))
  names(targets_by_tf) <- nodes
  comp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  for (ii in seq_along(nodes)) {
    for (jj in seq_along(nodes)) {
      if (ii == jj) next
      tf_i <- nodes[[ii]]
      tf_j <- nodes[[jj]]
      comp[ii, jj] <- as.numeric(mat[tf_i, tf_j] > 0) + as.numeric(mat[tf_j, tf_i] > 0) + 0.5 * length(intersect(targets_by_tf[[tf_i]], targets_by_tf[[tf_j]]))
    }
  }
  if (length(nodes)) {
    diag(comp) <- as.numeric(diag(mat[nodes, nodes, drop = FALSE]) > 0)
  }
  offdiag <- comp
  diag(offdiag) <- 0
  offdiag[offdiag < 5] <- 0
  keep_comp <- if (nrow(offdiag) && ncol(offdiag)) {
    pmax(apply(offdiag, 1, max, na.rm = TRUE), apply(offdiag, 2, max, na.rm = TRUE)) >= 5
  } else {
    logical()
  }
  if (sum(keep_comp) >= 2L) {
    comp_plot <- log2(pmax(offdiag[keep_comp, keep_comp, drop = FALSE], 0) + 1)
    plot_diff_grn_master_tf_connectivity_matrix(
      comp_plot,
      out_pdf_composite,
      paste0(sub("^TF hubs \\(delta fp_score\\) - ", "", title_text), " | TF-to-TF composite connectivity"),
      "log2(1+composite score)",
      subtitle = .diff_grn_master_wrap_text(
        "TF-to-TF composite connectivity. Score(i,j)=(i->j)+(j->i)+0.5*# shared targets; fill=log2(1+score).",
        width = 96
      ),
      caption = "Candidate TFs require >=1 outgoing and >=1 incoming TF-to-TF link.\nTF inclusion cutoff: max composite score >= 5. Display threshold: pairwise scores < 5 are shown as 0.",
      cell_size = 0.22,
      font_size = 8,
      row_annotation = .diff_grn_master_read_tf_cluster_annotation(),
      right_annotation = .diff_grn_master_direction_annotation(link_dt),
      x_label = NULL,
      y_label = NULL
    )
  } else {
    .diff_grn_master_empty_pdf(out_pdf_composite, paste0(title_text, " | TF-to-TF composite connectivity"), "No TFs pass composite cutoff")
  }
  invisible(c(out_pdf_mindist, out_pdf_composite))
}

.diff_grn_master_direction_annotation <- function(link_dt) {
  dt <- data.table::as.data.table(link_dt)
  need <- c("TF", "gene_key", "direction")
  if (!nrow(dt) || !all(need %chin% names(dt))) {
    return(data.table::data.table(tf = character(), direction_group = character()))
  }
  dt <- unique(dt[!is.na(TF) & nzchar(TF) & !is.na(gene_key) & nzchar(gene_key), .(
    tf = toupper(trimws(as.character(TF))),
    gene_key = as.character(gene_key),
    direction = trimws(as.character(direction))
  )])
  if (!nrow(dt)) return(data.table::data.table(tf = character(), direction_group = character()))
  out <- dt[, .(
    n_up = data.table::uniqueN(gene_key[direction == "Up"]),
    n_down = data.table::uniqueN(gene_key[direction == "Down"])
  ), by = tf]
  out[, direction_group := data.table::fifelse(
    n_up > 0 & n_down > 0,
    "Mixed",
    data.table::fifelse(n_up > 0, "Up", data.table::fifelse(n_down > 0, "Down", "None"))
  )]
  out[, .(tf, direction_group)]
}

#' Plot a TF connectivity matrix
#'
#' @param mat Numeric matrix.
#' @param out_pdf PDF output path.
#' @param title_text Plot title.
#' @param fill_title Fill legend title.
#' @param fill_limits Optional fill limits.
#' @param fill_breaks Optional fill breaks.
#' @param drop_empty Remove rows and columns with no nonzero off-diagonal
#'   connections before ordering and sizing the PDF.
#' @param subtitle Optional plot subtitle.
#' @param caption Optional bottom caption text.
#' @param cell_size Panel inches allocated to each heatmap row and column.
#' @param font_size Axis, legend, and title font size in points.
#' @param row_annotation Optional data frame with \code{tf} and \code{cluster}
#'   columns used to draw a row-side CRAFT-GRN cluster annotation.
#' @param col_annotation Optional data frame with \code{tf} and \code{cluster}
#'   columns used to draw a column-side cluster annotation.
#' @param right_annotation Optional data frame with \code{tf} and
#'   \code{direction_group} columns used to draw a right-side direction
#'   annotation.
#' @param row_annotation_title Title used for the row-side annotation legend.
#' @param col_annotation_title Title used for the column-side annotation.
#' @param right_annotation_title Title used for the right-side annotation
#'   legend.
#' @param row_cluster_gaps Add visible row gaps between adjacent row annotation
#'   clusters after clustering.
#' @param col_cluster_gaps Add visible column gaps between adjacent row
#'   annotation clusters after clustering. This is useful when rows and columns
#'   are the same TF universe.
#' @param cluster_gap_size Gap size, in tile units, inserted at each cluster
#'   boundary.
#' @param tile_border_color Color used for heatmap tile borders.
#' @param tile_border_linewidth Line width used for heatmap tile borders.
#' @param row_annotation_style Row annotation style, either legend based or
#'   inline cluster IDs.
#' @param row_annotation_side Side used for the row annotation strip.
#' @param row_annotation_width Width of the row annotation strip.
#' @param row_annotation_label_size Optional label size for inline row
#'   annotation IDs.
#' @param col_annotation_height Height of the column cluster annotation strip.
#' @param col_annotation_label_size Optional label size for inline column
#'   annotation IDs.
#' @param preserve_order Keep the input matrix row and column order instead of
#'   reclustering inside the plotting helper.
#' @param x_label Optional x-axis title.
#' @param y_label Optional y-axis title.
#' @param png_path Optional PNG output path written in addition to the PDF.
#'
#' @return Invisibly returns \code{out_pdf}.
#' @noRd
plot_diff_grn_master_tf_connectivity_matrix <- function(mat,
                                                        out_pdf,
                                                        title_text,
                                                        fill_title,
                                                        fill_limits = NULL,
                                                        fill_breaks = NULL,
                                                        drop_empty = TRUE,
                                                        subtitle = NULL,
                                                        caption = NULL,
                                                        cell_size = 0.18,
                                                        font_size = 9,
                                                        row_annotation = NULL,
                                                        col_annotation = NULL,
                                                        right_annotation = NULL,
                                                        row_annotation_title = "CRAFT-GRN cluster",
                                                        col_annotation_title = "Direct TF column cluster",
                                                        right_annotation_title = "Target direction",
                                                        row_cluster_gaps = FALSE,
                                                        col_cluster_gaps = FALSE,
                                                        cluster_gap_size = 0.35,
                                                        tile_border_color = "grey80",
                                                        tile_border_linewidth = 0.18,
                                                        row_annotation_style = c("legend", "inline_id"),
                                                        row_annotation_side = c("left", "right"),
                                                        row_annotation_width = 0.28,
                                                        row_annotation_label_size = NULL,
                                                        col_annotation_height = 0.9,
                                                        col_annotation_label_size = NULL,
                                                        preserve_order = FALSE,
                                                        x_label = "Targets",
                                                        y_label = "Regulators",
                                                        png_path = NULL) {
  if (!is.matrix(mat) || !length(mat)) return(invisible(FALSE))
  row_annotation_style <- match.arg(row_annotation_style)
  row_annotation_side <- match.arg(row_annotation_side)
  row_annotation_width <- max(0.05, as.numeric(row_annotation_width)[[1L]])
  if (!is.finite(row_annotation_width)) row_annotation_width <- 0.28
  if (is.null(row_annotation_label_size)) {
    row_annotation_label_size <- max(1.6, font_size / 4.2)
  }
  col_annotation_height <- max(0.05, as.numeric(col_annotation_height)[[1L]])
  if (!is.finite(col_annotation_height)) col_annotation_height <- 0.9
  if (is.null(col_annotation_label_size)) {
    col_annotation_label_size <- row_annotation_label_size
  }
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (isTRUE(drop_empty)) {
    nonzero <- is.finite(mat) & abs(mat) > 0
    if (nrow(mat) == ncol(mat) && identical(rownames(mat), colnames(mat))) {
      diag(nonzero) <- FALSE
    }
    row_keep <- rowSums(nonzero, na.rm = TRUE) > 0
    col_keep <- colSums(nonzero, na.rm = TRUE) > 0
    mat <- mat[row_keep, col_keep, drop = FALSE]
    if (!nrow(mat) || !ncol(mat)) {
      .diff_grn_master_empty_pdf(out_pdf, title_text, "No nonzero TF-to-TF connections after filtering")
      return(invisible(out_pdf))
    }
  }
  if (isTRUE(preserve_order)) {
    row_ord <- seq_len(nrow(mat))
  } else if (nrow(mat) >= 2L) {
    row_ord <- tryCatch(stats::hclust(stats::dist(mat), method = "complete")$order, error = function(e) seq_len(nrow(mat)))
  } else {
    row_ord <- seq_len(nrow(mat))
  }
  if (isTRUE(preserve_order)) {
    col_ord <- seq_len(ncol(mat))
  } else if (ncol(mat) >= 2L) {
    col_ord <- tryCatch(stats::hclust(stats::dist(t(mat)), method = "complete")$order, error = function(e) seq_len(ncol(mat)))
  } else {
    col_ord <- seq_len(ncol(mat))
  }
  mat <- mat[row_ord, col_ord, drop = FALSE]
  dt <- data.table::as.data.table(as.table(mat))
  names(dt) <- c("TF_row", "TF_col", "value")
  dt[, value := as.numeric(value)]
  row_levels <- rev(rownames(mat))
  col_levels <- colnames(mat)
  ann_dt <- data.table::data.table()
  if (!is.null(row_annotation) && nrow(data.table::as.data.table(row_annotation))) {
    ann_dt <- data.table::as.data.table(row_annotation)
    if (all(c("tf", "cluster") %chin% names(ann_dt))) {
      ann_dt[, tf := toupper(trimws(as.character(tf)))]
      ann_dt[, cluster := trimws(as.character(cluster))]
      ann_dt[is.na(cluster) | !nzchar(cluster), cluster := "Unassigned"]
      ann_dt <- unique(ann_dt[, .(tf, cluster)], by = "tf")
      ann_dt <- merge(
        data.table::data.table(TF_row = row_levels, row_index = seq_along(row_levels), tf = toupper(row_levels)),
        ann_dt,
        by = "tf",
        all.x = TRUE,
        sort = FALSE
      )
      ann_dt[is.na(cluster) | !nzchar(cluster), cluster := "Unassigned"]
    } else {
      ann_dt <- data.table::data.table()
    }
  }
  row_cluster_map <- if (nrow(ann_dt)) stats::setNames(ann_dt$cluster, ann_dt$TF_row) else character()
  make_gap_positions <- function(levels, cluster_map, do_gaps) {
    base_pos <- seq_along(levels)
    if (!isTRUE(do_gaps) || !length(cluster_map) || length(levels) < 2L) return(base_pos)
    clusters <- unname(cluster_map[levels])
    clusters[is.na(clusters) | !nzchar(clusters)] <- "Unassigned"
    boundary <- c(FALSE, clusters[-1L] != clusters[-length(clusters)])
    base_pos + cumsum(boundary) * cluster_gap_size
  }
  row_pos <- make_gap_positions(row_levels, row_cluster_map, row_cluster_gaps)
  col_pos <- make_gap_positions(col_levels, row_cluster_map, col_cluster_gaps)
  row_pos_map <- stats::setNames(row_pos, row_levels)
  col_pos_map <- stats::setNames(col_pos, col_levels)
  dt[, row_index := row_pos_map[as.character(TF_row)]]
  dt[, col_index := col_pos_map[as.character(TF_col)]]
  if (nrow(ann_dt)) {
    ann_dt[, row_index := row_pos_map[as.character(TF_row)]]
  }
  col_ann_dt <- data.table::data.table()
  if (!is.null(col_annotation) && nrow(data.table::as.data.table(col_annotation))) {
    col_ann_dt <- data.table::as.data.table(col_annotation)
    if (all(c("tf", "cluster") %chin% names(col_ann_dt))) {
      col_ann_dt[, tf := toupper(trimws(as.character(tf)))]
      col_ann_dt[, cluster := trimws(as.character(cluster))]
      col_ann_dt[is.na(cluster) | !nzchar(cluster), cluster := "Unassigned"]
      col_ann_dt <- unique(col_ann_dt[, .(tf, cluster)], by = "tf")
      col_ann_dt <- merge(
        data.table::data.table(TF_col = col_levels, col_index = seq_along(col_levels), tf = toupper(col_levels)),
        col_ann_dt,
        by = "tf",
        all.x = TRUE,
        sort = FALSE
      )
      col_ann_dt[is.na(cluster) | !nzchar(cluster), cluster := "Unassigned"]
      col_ann_dt[, col_index := col_pos_map[as.character(TF_col)]]
    } else {
      col_ann_dt <- data.table::data.table()
    }
  }
  if (is.null(fill_limits)) fill_limits <- c(0, max(dt$value, na.rm = TRUE))
  if (!all(is.finite(fill_limits)) || fill_limits[2] <= fill_limits[1]) fill_limits <- c(0, 1)
  if (is.null(fill_breaks)) fill_breaks <- pretty(fill_limits, n = 4)
  right_ann_dt <- data.table::data.table()
  if (!is.null(right_annotation) && nrow(data.table::as.data.table(right_annotation))) {
    right_ann_dt <- data.table::as.data.table(right_annotation)
    if (all(c("tf", "direction_group") %chin% names(right_ann_dt))) {
      right_ann_dt[, tf := toupper(trimws(as.character(tf)))]
      right_ann_dt[, direction_group := trimws(as.character(direction_group))]
      right_ann_dt[is.na(direction_group) | !nzchar(direction_group), direction_group := "None"]
      right_ann_dt <- unique(right_ann_dt[, .(tf, direction_group)], by = "tf")
      right_ann_dt <- merge(
        data.table::data.table(TF_row = row_levels, row_index = seq_along(row_levels), tf = toupper(row_levels)),
        right_ann_dt,
        by = "tf",
        all.x = TRUE,
        sort = FALSE
      )
      right_ann_dt[is.na(direction_group) | !nzchar(direction_group), direction_group := "None"]
    } else {
      right_ann_dt <- data.table::data.table()
    }
  }
  if (nrow(right_ann_dt)) {
    right_ann_dt[, row_index := row_pos_map[as.character(TF_row)]]
  }
  family <- .diff_grn_master_font_family()
  ann_cluster_levels <- sort(unique(c(
    if (nrow(ann_dt)) ann_dt$cluster else character(),
    if (nrow(col_ann_dt)) col_ann_dt$cluster else character()
  )))
  direction_levels <- c("Up", "Down", "Mixed", "None")
  direction_levels <- direction_levels[direction_levels %in% unique(right_ann_dt$direction_group)]
  has_left_annotation <- length(ann_cluster_levels) > 0L
  has_col_annotation <- nrow(col_ann_dt) > 0L
  has_right_annotation <- length(direction_levels) > 0L
  direction_label_chars <- if (has_right_annotation) {
    max(nchar(c(right_annotation_title, direction_levels)), na.rm = TRUE)
  } else {
    0
  }
  cluster_label_chars <- if (has_left_annotation) {
    max(nchar(c(row_annotation_title, ann_cluster_levels)), na.rm = TRUE)
  } else {
    0
  }
  col_max_pos <- max(col_pos, na.rm = TRUE)
  row_max_pos <- max(row_pos, na.rm = TRUE)
  annotation_gap <- 0.28
  right_ann_x <- col_max_pos + 0.55
  direction_legend_x <- col_max_pos + 1.35
  cluster_legend_x <- if (has_right_annotation) {
    direction_legend_x + 1.20 + 0.16 * direction_label_chars
  } else {
    col_max_pos + 1.45
  }
  right_expand <- if (has_left_annotation || has_right_annotation) {
    max(
      3.0,
      if (has_right_annotation) direction_legend_x - col_max_pos + 1.15 + 0.16 * direction_label_chars else 0,
      if (has_left_annotation) cluster_legend_x - col_max_pos + 1.15 + 0.16 * cluster_label_chars else 0
    )
  } else {
    0.05
  }
  left_annotation_active <- has_left_annotation && identical(row_annotation_side, "left")
  right_cluster_annotation_active <- has_left_annotation && identical(row_annotation_side, "right")
  row_annotation_x <- if (right_cluster_annotation_active) col_max_pos + 0.55 + annotation_gap + row_annotation_width / 2 else 0.25
  x_expand_add <- if (has_left_annotation || has_right_annotation) {
    c(if (left_annotation_active) row_annotation_width + annotation_gap + 0.37 else 0.05, right_expand + if (right_cluster_annotation_active) row_annotation_width + annotation_gap + 0.85 else 0)
  } else {
    c(0.05, 0.05)
  }
  top_expand <- if (has_col_annotation) col_annotation_height + 1.35 else if (has_left_annotation || has_right_annotation) 0.9 else 0.05
  y_expand_add <- c(0.05, top_expand)
  panel_x_span <- max(1, col_max_pos - min(col_pos, na.rm = TRUE) + sum(x_expand_add))
  panel_y_span <- max(1, row_max_pos - min(row_pos, na.rm = TRUE) + sum(y_expand_add))
  min_x_span <- if (has_left_annotation || has_right_annotation) 16.0 else 4.0
  min_y_span <- if (has_left_annotation || has_right_annotation) 5.2 else 4.0
  if (panel_x_span < min_x_span) {
    x_expand_add[2] <- x_expand_add[2] + min_x_span - panel_x_span
    panel_x_span <- min_x_span
  }
  if (panel_y_span < min_y_span) {
    y_expand_add[2] <- y_expand_add[2] + min_y_span - panel_y_span
    panel_y_span <- min_y_span
  }
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = col_index, y = row_index, fill = value)) +
    ggplot2::geom_tile(color = tile_border_color, linewidth = tile_border_linewidth) +
    ggplot2::scale_fill_gradientn(colors = c("#2166AC", "#FEE08B", "#B2182B"), limits = fill_limits, breaks = fill_breaks, oob = scales::squish, name = fill_title) +
    ggplot2::coord_fixed(ratio = 1, clip = if (has_left_annotation || has_right_annotation) "off" else "on") +
    ggplot2::scale_x_continuous(
      breaks = col_pos,
      labels = col_levels,
      expand = ggplot2::expansion(add = x_expand_add)
    ) +
    ggplot2::scale_y_continuous(
      breaks = row_pos,
      labels = row_levels,
      expand = ggplot2::expansion(add = y_expand_add)
    ) +
    ggplot2::labs(
      title = .diff_grn_master_wrap_text(title_text, width = 60),
      subtitle = .diff_grn_master_wrap_text(subtitle, width = 58),
      x = x_label,
      y = y_label,
      caption = .diff_grn_master_wrap_text(caption, width = 78)
    ) +
    ggplot2::theme_minimal(base_size = font_size, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black"),
      plot.title = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black", hjust = 0.5),
      axis.title.x = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black", margin = ggplot2::margin(t = 14)),
      axis.title.y = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = font_size, face = "bold", family = family, color = "black", margin = ggplot2::margin(t = 4)),
      axis.text.y = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black"),
      legend.title = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black"),
      legend.text = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black"),
      plot.caption = ggplot2::element_text(size = font_size, face = "bold", family = family, color = "black", hjust = 0, margin = ggplot2::margin(t = 12)),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, if (has_left_annotation || has_right_annotation) 100 else 8, if (has_left_annotation || has_right_annotation) 28 else 8, 8)
    )
  if (nrow(right_ann_dt)) {
    direction_cols <- c(Up = "#FC8D62", Down = "#66C2A5", Mixed = "#BDBDBD", None = "#F0F0F0")
    right_ann_dt[, direction_color := direction_cols[direction_group]]
    direction_legend_dt <- data.table::data.table(
      direction_group = direction_levels,
      direction_color = direction_cols[direction_levels],
      legend_y = seq_along(direction_levels)
    )
    p <- p +
      ggplot2::geom_tile(
        data = right_ann_dt,
        ggplot2::aes(x = right_ann_x, y = row_index),
        inherit.aes = FALSE,
        width = 0.28,
        height = 1,
        fill = right_ann_dt$direction_color,
        color = tile_border_color,
        linewidth = tile_border_linewidth
      ) +
      ggplot2::geom_tile(
        data = direction_legend_dt,
        ggplot2::aes(x = direction_legend_x, y = legend_y),
        inherit.aes = FALSE,
        width = 0.35,
        height = 0.8,
        fill = direction_legend_dt$direction_color,
        color = tile_border_color,
        linewidth = tile_border_linewidth
      ) +
      ggplot2::geom_text(
        data = direction_legend_dt,
        ggplot2::aes(x = direction_legend_x + 0.24, y = legend_y, label = direction_group),
        inherit.aes = FALSE,
        hjust = 0,
        size = max(2.2, font_size / 3),
        fontface = "bold",
        family = family
      ) +
      ggplot2::annotate(
        "text",
        x = direction_legend_x,
        y = length(direction_levels) + 0.7,
        label = right_annotation_title,
        hjust = 0,
        size = max(2.2, font_size / 3),
        fontface = "bold",
        family = family
      )
  }
  if (nrow(ann_dt)) {
    cluster_levels <- ann_cluster_levels
    cluster_cols <- .diff_grn_master_direct_tf_cluster_colors(cluster_levels)
    ann_dt[, cluster_color := cluster_cols[cluster]]
    p <- p +
      ggplot2::geom_tile(
        data = ann_dt,
        ggplot2::aes(x = row_annotation_x, y = row_index),
        inherit.aes = FALSE,
        width = row_annotation_width,
        height = 1,
        fill = ann_dt$cluster_color,
        color = tile_border_color,
        linewidth = tile_border_linewidth
      )
    if (identical(row_annotation_style, "inline_id")) {
      ann_label_dt <- data.table::copy(ann_dt)[order(row_index)]
      ann_label_dt[, run_id := data.table::rleid(cluster)]
      ann_label_dt <- ann_label_dt[
        ,
        .(
          row_index = mean(row_index, na.rm = TRUE),
          cluster_label = cluster[[1L]]
        ),
        by = run_id
      ]
      p <- p +
        ggplot2::geom_text(
          data = ann_label_dt,
          ggplot2::aes(x = row_annotation_x, y = row_index, label = cluster_label),
          inherit.aes = FALSE,
          size = row_annotation_label_size,
          fontface = "bold",
          family = family,
          color = "black"
        ) +
        ggplot2::annotate(
          "text",
          x = row_annotation_x,
          y = max(row_pos, na.rm = TRUE) + 0.55,
          label = row_annotation_title,
          hjust = 0.5,
          size = max(row_annotation_label_size * 0.72, font_size / 3.2),
          fontface = "bold",
          family = family
        )
    } else {
      legend_dt <- data.table::data.table(
        cluster = cluster_levels,
        cluster_color = cluster_cols[cluster_levels],
        legend_y = seq_along(cluster_levels)
      )
      legend_x <- cluster_legend_x
      p <- p +
        ggplot2::geom_tile(
          data = legend_dt,
          ggplot2::aes(x = legend_x, y = legend_y),
          inherit.aes = FALSE,
          width = 0.35,
          height = 0.8,
          fill = legend_dt$cluster_color,
          color = tile_border_color,
          linewidth = tile_border_linewidth
        ) +
        ggplot2::geom_text(
          data = legend_dt,
          ggplot2::aes(x = legend_x + 0.24, y = legend_y, label = cluster),
          inherit.aes = FALSE,
          hjust = 0,
          size = max(2.2, font_size / 3),
          fontface = "bold",
          family = family
        ) +
        ggplot2::annotate(
          "text",
          x = legend_x,
          y = length(cluster_levels) + 0.7,
          label = row_annotation_title,
          hjust = 0,
          size = max(2.2, font_size / 3),
          fontface = "bold",
          family = family
      )
    }
  }
  if (nrow(col_ann_dt)) {
    col_ann_dt[, cluster_color := "#FFFFFF"]
    col_annotation_y <- row_max_pos + 0.55 + annotation_gap + col_annotation_height / 2
    p <- p +
      ggplot2::geom_tile(
        data = col_ann_dt,
        ggplot2::aes(x = col_index, y = col_annotation_y),
        inherit.aes = FALSE,
        width = 1,
        height = col_annotation_height,
        fill = col_ann_dt$cluster_color,
        color = tile_border_color,
        linewidth = tile_border_linewidth
      )
    col_label_dt <- data.table::copy(col_ann_dt)[order(col_index)]
    col_label_dt[, run_id := data.table::rleid(cluster)]
    col_label_dt <- col_label_dt[
      ,
      .(
        col_index = mean(col_index, na.rm = TRUE),
        cluster_label = cluster[[1L]]
      ),
      by = run_id
    ]
    p <- p +
      ggplot2::geom_text(
        data = col_label_dt,
        ggplot2::aes(x = col_index, y = col_annotation_y, label = cluster_label),
        inherit.aes = FALSE,
        size = col_annotation_label_size,
        fontface = "bold",
        family = family,
        color = "black"
      ) +
      ggplot2::annotate(
        "text",
        x = min(col_pos, na.rm = TRUE),
        y = col_annotation_y + col_annotation_height * 0.75,
        label = col_annotation_title,
        hjust = 0,
        size = max(col_annotation_label_size * 0.72, font_size / 3.2),
        fontface = "bold",
        family = family
      )
  }
  .diff_grn_master_save_fixed_panel_pdf(
    p,
    out_pdf,
    panel_width = cell_size * panel_x_span,
    panel_height = cell_size * panel_y_span,
    family = family,
    png_path = png_path
  )
  invisible(out_pdf)
}

.diff_grn_master_normalize_links <- function(links_dt) {
  dt <- data.table::as.data.table(links_dt)
  tf_col <- if ("tf" %in% names(dt)) "tf" else if ("TF" %in% names(dt)) "TF" else NULL
  gene_col <- if ("gene_key" %in% names(dt)) "gene_key" else if ("gene" %in% names(dt)) "gene" else if ("target_gene" %in% names(dt)) "target_gene" else NULL
  peak_col <- if ("peak_id" %in% names(dt)) "peak_id" else NULL
  if (is.null(tf_col) || is.null(gene_col)) return(data.table::data.table())
  if ("direction_group" %in% names(dt)) {
    dt[, direction := ifelse(grepl("^up", tolower(trimws(as.character(direction_group)))), "Up", "Down")]
  } else if ("log2FC_gene_expr" %in% names(dt)) {
    dt[, direction := ifelse(suppressWarnings(as.numeric(log2FC_gene_expr)) >= 0, "Up", "Down")]
  } else {
    dt[, direction := "Up"]
  }
  cols <- c(tf_col, gene_col, peak_col, "direction")
  out <- unique(dt[, cols, with = FALSE])
  data.table::setnames(out, c(tf_col, gene_col), c("TF", "gene_key"))
  out <- out[!is.na(TF) & nzchar(TF) & !is.na(gene_key) & nzchar(gene_key)]
  tf_upper <- toupper(unique(out$TF))
  out[, target_type := ifelse(toupper(gene_key) %in% tf_upper, "TF target", "Gene target")]
  out[]
}

.diff_grn_master_mindist_score <- function(mat) {
  adj <- mat > 0
  diag(adj) <- FALSE
  n <- nrow(adj)
  dist_mat <- matrix(Inf, n, n, dimnames = dimnames(adj))
  for (src in seq_len(n)) {
    d <- rep(Inf, n)
    d[src] <- 0
    q <- src
    while (length(q)) {
      v <- q[[1]]
      q <- q[-1]
      nb <- which(adj[v, ] & is.infinite(d))
      if (length(nb)) {
        d[nb] <- d[v] + 1
        q <- c(q, nb)
      }
    }
    dist_mat[src, ] <- d
  }
  out <- matrix(0, n, n, dimnames = dimnames(mat))
  valid <- is.finite(dist_mat) & dist_mat > 0
  out[valid] <- 1 / dist_mat[valid]
  diag(out) <- as.numeric(diag(mat) > 0)
  out
}

.diff_grn_master_empty_pdf <- function(out_pdf, title_text, message) {
  grDevices::pdf(out_pdf, width = 10, height = 10, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot.new()
  graphics::title(main = title_text)
  graphics::text(0.5, 0.5, message)
  invisible(out_pdf)
}

.diff_grn_master_contrast_parts <- function(x) {
  stem <- sub("\\.(csv|pdf|html)$", "", basename(x), ignore.case = TRUE)
  stem <- sub("_master_tf_summary.*$", "", stem)
  stem <- sub("_filtered_links.*$", "", stem)
  stem <- sub("_delta_links.*$", "", stem)
  parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
  cond1 <- if (length(parts) >= 1L) parts[[1L]] else NA_character_
  cond2 <- if (length(parts) >= 2L) parts[[2L]] else NA_character_
  list(cond1 = cond1, cond2 = cond2, label = if (length(parts) == 2L) paste(cond1, "vs", cond2) else stem)
}

.diff_grn_master_escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

.diff_grn_master_robust_z <- function(x) {
  x <- as.numeric(x)
  m <- stats::median(x, na.rm = TRUE)
  madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(madv) || madv == 0) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
    return((x - m) / (sx + 1e-12))
  }
  (x - m) / (madv + 1e-12)
}

.diff_grn_master_clamp <- function(x, lo, hi) {
  pmax(pmin(x, hi), lo)
}

.diff_grn_master_wrap_text <- function(x, width = 90L) {
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

.diff_grn_master_font_family <- function() {
  pdf_preferred <- c("ArialMT", "Helvetica", "NimbusSan", "URWHelvetica", "sans")
  pdf_fonts <- names(grDevices::pdfFonts())
  hit <- pdf_preferred[pdf_preferred %in% pdf_fonts]
  if (length(hit)) return(hit[[1L]])
  "sans"
}

.diff_grn_master_direct_tf_cluster_colors <- function(cluster_levels) {
  cluster_levels <- as.character(cluster_levels)
  palette <- c(
    "#2B6CB0",
    "#fb5311",
    "#2F855A",
    "#6B46C1",
    "#B7791F",
    "#0F766E",
    "#9F1239",
    "#4A5568",
    "#1A365D",
    "#7C2D12"
  )
  idx <- suppressWarnings(as.integer(sub("^[A-Za-z]+0*", "", cluster_levels)))
  out <- grDevices::hcl.colors(length(cluster_levels), palette = "Set 2")
  use_fixed <- is.finite(idx) & idx >= 1L & idx <= length(palette) & grepl("^[RC][0-9]+$", cluster_levels)
  out[use_fixed] <- palette[idx[use_fixed]]
  stats::setNames(out, cluster_levels)
}

.diff_grn_master_save_fixed_panel_pdf <- function(plot_obj, pdf_path, panel_width, panel_height, family = "sans", png_path = NULL, png_dpi = 150) {
  gt <- ggplot2::ggplotGrob(plot_obj)
  panel_idx <- grepl("^panel", gt$layout$name)
  panel_cols <- unique(gt$layout$l[panel_idx])
  panel_rows <- unique(gt$layout$t[panel_idx])
  gt$widths[panel_cols] <- grid::unit(panel_width, "in")
  gt$heights[panel_rows] <- grid::unit(panel_height, "in")
  width <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
  height <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = family)
  grid::grid.newpage()
  grid::grid.draw(gt)
  grDevices::dev.off()
  if (!is.null(png_path) && nzchar(png_path)) {
    grDevices::png(png_path, width = width, height = height, units = "in", res = png_dpi, type = "cairo")
    grid::grid.newpage()
    grid::grid.draw(gt)
    grDevices::dev.off()
  }
  invisible(pdf_path)
}

utils::globalVariables(c(
  "TF",
  "TF_col",
  "TF_row",
  "TF_src",
  "TF_tgt",
  "cell_rank",
  "cluster",
  "cluster_color",
  ".data",
  "delta",
  "direction",
  "direction_color",
  "direction_group",
  "gene_key",
  "legend_y",
  "log2FC_gene_expr",
  "n",
  "n_down",
  "n_up",
  "net_n",
  "reg_sign",
  "signed_log2_n",
  "stress_rank",
  "target_type",
  "tf_delta_sum",
  "tf_delta_sum_abs",
  "tf_expr_max",
  "tf_log2_fc",
  "tf_log2_fc_med",
  "tf_log2_expr_max",
  "tf_links",
  "tf_n_links",
  "tf_n_target_genes",
  "tf_target_genes",
  "value",
  "x_sum",
  "x_sum_signed_log10"
))
