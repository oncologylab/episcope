#' Module 2 QC helpers
#'
#' Utilities for summarizing TF->TFBS->target link activity and producing
#' QC outputs for Module 2.
#'
#' @name module2_qc_helpers
#' @rdname module2_qc_helpers
#' @noRd
NULL

#' Build TF->target overview table with condition flags
#'
#' @param links Link table with TF/gene/peak keys and correlation statistics.
#' @param link_status Either a data.frame or a CSV path produced by
#'   [build_link_status_matrix()].
#' @param out_file Optional CSV path for the merged overview table.
#' @param verbose Emit status messages.
#'
#' @return Merged tibble (invisibly).
#' @export
build_tf_target_overview <- function(
  links,
  link_status,
  out_file = NULL,
  verbose = TRUE
) {
  if (!is.data.frame(links)) .log_abort("`links` must be a data.frame.")
  if (is.character(link_status)) {
    if (!file.exists(link_status)) .log_abort("link_status file not found: {link_status}")
    link_status <- readr::read_csv(link_status, show_col_types = FALSE)
  }
  if (!is.data.frame(link_status)) .log_abort("`link_status` must be a data.frame or CSV path.")

  key_cols <- c("TF", "gene_key", "peak_ID")
  if (!all(key_cols %in% names(links))) {
    .log_abort("`links` missing required columns: {paste(setdiff(key_cols, names(links)), collapse = ', ')}")
  }
  if (!all(key_cols %in% names(link_status))) {
    .log_abort("`link_status` missing required columns: {paste(setdiff(key_cols, names(link_status)), collapse = ', ')}")
  }

  out <- links |>
    dplyr::left_join(link_status, by = key_cols)

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out, out_file)
    if (isTRUE(verbose)) .log_inform("TF-target overview written: {out_file}")
  }

  invisible(out)
}

#' Summarize link activity per condition and per TF
#'
#' @param link_status Either a data.frame or a CSV path produced by
#'   [build_link_status_matrix()].
#' @param out_dir Output directory for summary CSVs.
#' @param db Database tag used in output filenames.
#' @param prefix Prefix for output filenames.
#' @param verbose Emit status messages.
#'
#' @return List with summary_total and summary_by_tf (invisible).
#' @export
summarize_link_activity <- function(
  link_status,
  out_dir,
  db,
  prefix = "step2",
  verbose = TRUE
) {
  if (is.character(link_status)) {
    if (!file.exists(link_status)) .log_abort("link_status file not found: {link_status}")
    link_status <- readr::read_csv(link_status, show_col_types = FALSE)
  }
  if (!is.data.frame(link_status)) .log_abort("`link_status` must be a data.frame or CSV path.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  key_cols <- c("TF", "gene_key", "peak_ID")
  if (!all(key_cols %in% names(link_status))) {
    .log_abort("`link_status` missing required columns: {paste(setdiff(key_cols, names(link_status)), collapse = ', ')}")
  }

  cond_cols <- setdiff(names(link_status), key_cols)
  if (!length(cond_cols)) .log_abort("No condition columns found in link_status.")

  status_long <- link_status |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(cond_cols),
      names_to = "condition",
      values_to = "active"
    )

  summary_total <- status_long |>
    dplyr::group_by(.data$condition) |>
    dplyr::summarise(n_links_active = sum(.data$active > 0, na.rm = TRUE), .groups = "drop")

  summary_by_tf_long <- status_long |>
    dplyr::filter(.data$active > 0) |>
    dplyr::group_by(.data$TF, .data$condition) |>
    dplyr::summarise(
      n_links_active = dplyr::n(),
      n_genes_active = dplyr::n_distinct(.data$gene_key),
      .groups = "drop"
    )

  summary_by_tf <- summary_by_tf_long |>
    dplyr::select("TF", "condition", "n_links_active") |>
    tidyr::pivot_wider(
      names_from = "condition",
      values_from = "n_links_active",
      values_fill = 0
    )

  total_row <- summary_total |>
    dplyr::mutate(TF = "__TOTAL__") |>
    tidyr::pivot_wider(
      names_from = "condition",
      values_from = "n_links_active",
      values_fill = 0
    ) |>
    dplyr::select(TF, dplyr::everything())

  summary_by_tf <- dplyr::bind_rows(total_row, summary_by_tf)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  by_tf_path <- file.path(out_dir, "link_activity_summary.csv")
  readr::write_csv(summary_by_tf, by_tf_path)
  total_path_legacy <- file.path(out_dir, sprintf("%s_link_activity_total_%s.csv", prefix, db))
  if (file.exists(total_path_legacy)) {
    file.remove(total_path_legacy)
  }
  by_tf_legacy <- file.path(out_dir, sprintf("%s_link_activity_by_tf_%s.csv", prefix, db))
  if (file.exists(by_tf_legacy)) {
    file.remove(by_tf_legacy)
  }

  if (isTRUE(verbose)) {
    .log_inform("Link activity by TF written: {by_tf_path}")
  }

  invisible(list(summary_total = summary_total, summary_by_tf = summary_by_tf))
}

#' Plot total active links per condition
#'
#' @param summary_total A data.frame with columns condition and n_links_active.
#' @param summary_by_tf Optional wide table with TF in first column and
#'   per-condition active-link counts in remaining columns.
#' @param out_dir Output directory for the PDF.
#' @param db Database tag used in the output filename.
#' @param prefix Prefix for output filename.
#' @param verbose Emit status messages.
#'
#' @return Path to the written PDF (invisible).
#' @export
plot_link_activity_qc <- function(
  summary_total,
  summary_by_tf = NULL,
  out_dir,
  db,
  prefix = "step2",
  verbose = TRUE
) {
  if (!is.data.frame(summary_total) || !all(c("condition", "n_links_active") %in% names(summary_total))) {
    .log_abort("`summary_total` must include condition and n_links_active columns.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .log_abort("Missing required package: ggplot2.")
  }
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  p <- ggplot2::ggplot(summary_total, ggplot2::aes(x = condition, y = n_links_active)) +
    ggplot2::geom_col(fill = "#3182bd") +
    ggplot2::labs(
      title = "Total active links per condition",
      x = "Condition",
      y = "Active links"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )

  p_tf <- NULL
  if (is.data.frame(summary_by_tf) && ncol(summary_by_tf) >= 2L) {
    tf_col <- names(summary_by_tf)[1]
    tf_df <- summary_by_tf
    if ("__TOTAL__" %in% tf_df[[tf_col]]) {
      tf_df <- tf_df[tf_df[[tf_col]] != "__TOTAL__", , drop = FALSE]
    }
    if (nrow(tf_df)) {
      cond_cols <- setdiff(names(tf_df), tf_col)
      tf_long <- tidyr::pivot_longer(
        tf_df,
        cols = dplyr::all_of(cond_cols),
        names_to = "condition",
        values_to = "n_links_active"
      )
      tf_tot <- tf_long |>
        dplyr::group_by(.data[[tf_col]]) |>
        dplyr::summarise(total_links = sum(.data$n_links_active, na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(.data$total_links))
      tf_levels <- tf_tot[[tf_col]]
      tf_long$TF_label <- factor(tf_long[[tf_col]], levels = rev(tf_levels))
      mat <- stats::xtabs(n_links_active ~ TF_label + condition, data = tf_long)
      if (nrow(mat) >= 2 && ncol(mat) >= 2) {
        hc_row <- stats::hclust(stats::dist(mat), method = "complete")
        hc_col <- stats::hclust(stats::dist(t(mat)), method = "complete")
        row_ord <- rownames(mat)[hc_row$order]
        col_ord <- colnames(mat)[hc_col$order]
      } else {
        row_ord <- rownames(mat)
        col_ord <- colnames(mat)
      }
      tf_long$TF_label <- factor(tf_long$TF_label, levels = row_ord)
      tf_long$condition <- factor(tf_long$condition, levels = col_ord)
      p_tf <- ggplot2::ggplot(tf_long, ggplot2::aes(x = .data$condition, y = .data$TF_label, fill = .data$n_links_active)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
        ggplot2::labs(
          title = "Active link counts per TF per condition (unsupervised clustering)",
          x = "Condition",
          y = "TF",
          fill = "Active links"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        )
    }
  }

  pdf_path <- file.path(out_dir, "link_activity_summary.pdf")
  n_cond <- max(1L, nrow(summary_total))
  bar_w <- max(10, min(18, 7 + 0.15 * n_cond))
  bar_h <- max(7, min(12, 5 + 0.08 * n_cond))
  grDevices::pdf(pdf_path, width = bar_w, height = bar_h)
  print(p)
  if (!is.null(p_tf)) {
    print(p_tf)
  }
  grDevices::dev.off()

  if (isTRUE(verbose)) .log_inform("Link activity summary saved: {pdf_path}")
  invisible(pdf_path)
}
