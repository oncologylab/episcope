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
  total_path <- file.path(out_dir, sprintf("%s_link_activity_total_%s.csv", prefix, db))
  by_tf_path <- file.path(out_dir, sprintf("%s_link_activity_by_tf_%s.csv", prefix, db))
  readr::write_csv(summary_total, total_path)
  readr::write_csv(summary_by_tf, by_tf_path)

  if (isTRUE(verbose)) {
    .log_inform("Link activity summary written: {total_path}")
    .log_inform("Link activity by TF written: {by_tf_path}")
  }

  invisible(list(summary_total = summary_total, summary_by_tf = summary_by_tf))
}

#' Plot total active links per condition
#'
#' @param summary_total A data.frame with columns condition and n_links_active.
#' @param out_dir Output directory for the PDF.
#' @param db Database tag used in the output filename.
#' @param prefix Prefix for output filename.
#' @param verbose Emit status messages.
#'
#' @return Path to the written PDF (invisible).
#' @export
plot_link_activity_qc <- function(
  summary_total,
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
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Total active links per condition",
      x = "Condition",
      y = "Active links"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  pdf_path <- file.path(out_dir, sprintf("%s_link_activity_qc_%s.pdf", prefix, db))
  grDevices::pdf(pdf_path, width = 8, height = 6)
  print(p)
  grDevices::dev.off()

  if (isTRUE(verbose)) .log_inform("Link activity QC saved: {pdf_path}")
  invisible(pdf_path)
}

#' Extract per-condition link tables from lighting output
#'
#' @param out_dir Directory containing per-condition link CSVs.
#' @param prefix Prefix used by [light_by_condition()].
#' @param read_tables If TRUE, read and return a combined table with
#'   a `condition` column.
#' @param verbose Emit status messages.
#'
#' @return Tibble with condition, path, and n_links (and optionally combined links).
#' @export
extract_link_info_by_condition <- function(
  out_dir,
  prefix = "step2",
  read_tables = FALSE,
  verbose = TRUE
) {
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!dir.exists(out_dir)) .log_abort("`out_dir` does not exist: {out_dir}")

  pattern <- sprintf("^%s_cond-.*_tf_gene_links\\.csv$", prefix)
  csvs <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  if (!length(csvs)) .log_abort("No per-condition link CSVs found in {out_dir} with prefix {prefix}.")

  parse_label <- function(path) {
    base <- basename(path)
    sub(sprintf("^%s_cond-", prefix), "", sub("_tf_gene_links\\.csv$", "", base))
  }

  manifest <- tibble::tibble(
    condition = vapply(csvs, parse_label, character(1)),
    path = csvs
  )

  manifest$n_links <- vapply(csvs, function(p) {
    tryCatch(nrow(readr::read_csv(p, show_col_types = FALSE)), error = function(e) NA_integer_)
  }, integer(1))

  if (isTRUE(read_tables)) {
    combined <- purrr::map2_dfr(csvs, manifest$condition, function(p, cond) {
      readr::read_csv(p, show_col_types = FALSE) |>
        dplyr::mutate(condition = cond)
    })
    if (isTRUE(verbose)) .log_inform("Loaded {nrow(combined)} per-condition links.")
    return(invisible(list(manifest = manifest, links = combined)))
  }

  if (isTRUE(verbose)) {
    .log_inform("Found {nrow(manifest)} per-condition link tables in {out_dir}.")
  }
  invisible(manifest)
}
