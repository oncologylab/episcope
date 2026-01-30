#' Extract TF-hub axes from LDA summary tables
#'
#' Given one or more `*_delta_links_filtered_lda_K{K}_summary.csv` files (the same
#' inputs used by `plot_from_summary_bulk()`), compute the per-TF values that drive:
#'   1) the hub-plot X axis (signed sums; or signed non-cancel),
#'   2) the hub-plot Y axis (#links, with a convenience log2p1 column).
#'
#' We return six tibbles:
#'   - x_activate   : signed x for activate-only (per TF × comparison)
#'   - x_repress    : signed x for repress-only  (per TF × comparison)
#'   - x_noncancel  : signed non-cancel x using OVERALL sign & abs(activate)+abs(repress)
#'   - y_activate   : #links for activate-only   (+ log2p1 convenience)
#'   - y_repress    : #links for repress-only    (+ log2p1 convenience)
#'   - y_noncancel  : #links for OVERALL         (+ log2p1 convenience)
#'
#' Notes:
#' - Signed x (activate/repress):  sign(sum_delta) * log2(abs(sum_delta)+1).
#' - Non-cancel x (overall): sign(overall_sum_delta) * log2(sum_abs_activate + sum_abs_repress + 1).
#' - Y uses raw counts; plotting applies log2 scale, so we also provide `*_log2p1`.
#'
#' @param summary_csvs Character vector of existing summary CSV paths.
#' @return A named list of six tibbles.
#' @examples
#' # res <- extract_tf_hub_axes(summary_files)
#' # res$x_activate |> dplyr::arrange(comparison, dplyr::desc(x_signed))
#' @export
extract_tf_hub_axes <- function(summary_csvs) {
  if (!is.character(summary_csvs) || length(summary_csvs) == 0L) {
    cli::cli_abort("`summary_csvs` must be a non-empty character vector of paths.")
  }
  summary_csvs <- summary_csvs[file.exists(summary_csvs)]
  if (!length(summary_csvs)) cli::cli_abort("No existing summary CSVs found.")

  .contrast_from_file <- function(f) {
    b <- basename(f)
    stem <- sub("_delta_links.*$", "", b)
    parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) == 2) paste(parts[1], "vs", parts[2]) else stem
  }

  # Per-file → per-TF aggregates that mirror make_tf_hubs_from_summary()
  one_file_summaries <- lapply(summary_csvs, function(f) {
    sz <- suppressWarnings(file.info(f)$size)
    if (!is.finite(sz) || is.na(sz) || sz == 0) return(NULL)

    S <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(S) || !NROW(S)) return(NULL)

    need <- c(
      "TF",
      "tf_n_links", "tf_n_links_activate", "tf_n_links_repress",
      "tf_delta_sum", "tf_delta_sum_activate", "tf_delta_sum_repress",
      "tf_delta_sum_abs", "tf_delta_sum_abs_activate", "tf_delta_sum_abs_repress"
    )
    if (!all(need %in% names(S))) {
      cli::cli_warn("Skipping {f}: missing required columns.")
      return(NULL)
    }

    # Aggregate across topics → per-TF for this comparison
    d <- dplyr::group_by(S, TF)
    df <- dplyr::summarise(
      d,
      tf_links_overall  = sum(tf_n_links,            na.rm = TRUE),
      tf_links_act      = sum(tf_n_links_activate,   na.rm = TRUE),
      tf_links_rep      = sum(tf_n_links_repress,    na.rm = TRUE),

      tf_sum_delta_overall = sum(tf_delta_sum,            na.rm = TRUE),
      tf_sum_delta_act     = sum(tf_delta_sum_activate,   na.rm = TRUE),
      tf_sum_delta_rep     = sum(tf_delta_sum_repress,    na.rm = TRUE),

      tf_sum_abs_overall   = sum(tf_delta_sum_abs,           na.rm = TRUE),
      tf_sum_abs_act       = sum(tf_delta_sum_abs_activate,  na.rm = TRUE),
      tf_sum_abs_rep       = sum(tf_delta_sum_abs_repress,   na.rm = TRUE),

      .groups = "drop"
    )

    # X (signed) for activate/repress
    x_signed_act <- sign(df$tf_sum_delta_act) * log2(abs(df$tf_sum_delta_act) + 1)
    x_signed_rep <- sign(df$tf_sum_delta_rep) * log2(abs(df$tf_sum_delta_rep) + 1)

    # X (signed non-cancel) overall: sign from OVERALL signed sum; magnitude from abs(act)+abs(rep)
    x_signed_nc <- {
      sgn_overall <- sign(df$tf_sum_delta_overall)
      mag_nc      <- df$tf_sum_abs_act + df$tf_sum_abs_rep
      sgn_overall * log2(mag_nc + 1)
    }

    # Y raw + convenience log2(1+y) (the plot uses log2 scale)
    y_overall <- df$tf_links_overall
    y_act     <- df$tf_links_act
    y_rep     <- df$tf_links_rep

    tibble::tibble(
      comparison = .contrast_from_file(f),
      TF         = df$TF,
      x_signed_activate   = x_signed_act,
      x_signed_repress    = x_signed_rep,
      x_signed_noncancel  = x_signed_nc,
      y_overall           = y_overall,
      y_overall_log2p1    = log2(pmax(y_overall, 0) + 1),
      y_activate          = y_act,
      y_activate_log2p1   = log2(pmax(y_act, 0) + 1),
      y_repress           = y_rep,
      y_repress_log2p1    = log2(pmax(y_rep, 0) + 1)
    )
  })

  one_file_summaries <- one_file_summaries[!vapply(one_file_summaries, is.null, logical(1))]
  if (!length(one_file_summaries)) {
    cli::cli_abort("All inputs were empty or malformed.")
  }
  all <- dplyr::bind_rows(one_file_summaries)

  # Split into the six requested tibbles (keep only the relevant cols)
  x_activate  <- all[, c("comparison","TF","x_signed_activate")]
  x_repress   <- all[, c("comparison","TF","x_signed_repress")]
  x_noncancel <- all[, c("comparison","TF","x_signed_noncancel")]

  y_activate  <- all[, c("comparison","TF","y_activate","y_activate_log2p1")]
  y_repress   <- all[, c("comparison","TF","y_repress","y_repress_log2p1")]
  y_noncancel <- all[, c("comparison","TF","y_overall","y_overall_log2p1")]

  # Nice names
  names(x_activate)[3]  <- "x_signed"
  names(x_repress)[3]   <- "x_signed"
  names(x_noncancel)[3] <- "x_signed"

  names(y_activate)[3:4]  <- c("y_links","y_links_log2p1")
  names(y_repress)[3:4]   <- c("y_links","y_links_log2p1")
  names(y_noncancel)[3:4] <- c("y_links","y_links_log2p1")

  list(
    x_activate   = x_activate,
    x_repress    = x_repress,
    x_noncancel  = x_noncancel,
    y_activate   = y_activate,
    y_repress    = y_repress,
    y_noncancel  = y_noncancel
  )
}

summary_files <- list.files(path = lighting_folder, pattern = "_lda_K20_summary.csv", full.names = TRUE)
summary_files_vs_culture <- grep("_vs_culture", summary_files, value = TRUE)
summary_files_vs_Ctrl    <- grep("_vs_Ctrl", summary_files, value = TRUE)

# --- Minimal 6-heatmap driver for one set of summary CSVs ---------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(matrixStats)
})
source("R/utils_ggheat.R")

plot_hub_heatmaps_simple <- function(
    summary_files,
    out_dir = "/data/homes/yl814/episcope_test/GSE87218_ATAC/heatmap_pca",
    tag     = c("vs_culture","vs_Ctrl")[1],
    x_thresh = 1,   # keep TFs if any |x_signed| > 1
    y_thresh = 1    # keep TFs if any y_links_log2p1 > 1
){
  stopifnot(length(summary_files) > 0)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_dir <- file.path(out_dir, tag); dir.create(out_dir, showWarnings = FALSE)

  res <- extract_tf_hub_axes(summary_files)

  pivot_mat <- function(df, value_col) {
    wide <- tidyr::pivot_wider(df[, c("comparison","TF", value_col)], id_cols = "TF",
                               names_from = "comparison", values_from = value_col, values_fill = 0)
    mat <- as.matrix(wide[, setdiff(names(wide), "TF"), drop = FALSE])
    rownames(mat) <- wide$TF
    storage.mode(mat) <- "numeric"
    mat
  }

  # Build matrices (columns = comparisons, rows = TFs)
  X_act <- pivot_mat(res$x_activate,  "x_signed")
  X_rep <- pivot_mat(res$x_repress,   "x_signed")
  X_nc  <- pivot_mat(res$x_noncancel, "x_signed")

  Y_act <- pivot_mat(res$y_activate,  "y_links_log2p1")
  Y_rep <- pivot_mat(res$y_repress,   "y_links_log2p1")
  Y_nc  <- pivot_mat(res$y_noncancel, "y_links_log2p1")

  # Simple row filters
  f_keep_x <- function(M) M[matrixStats::rowMaxs(abs(M), na.rm = TRUE) > x_thresh, , drop = FALSE]
  f_keep_y <- function(M) M[matrixStats::rowMaxs(M,       na.rm = TRUE) > y_thresh, , drop = FALSE]

  X_act <- f_keep_x(X_act); X_rep <- f_keep_x(X_rep); X_nc  <- f_keep_x(X_nc)
  Y_act <- f_keep_y(Y_act); Y_rep <- f_keep_y(Y_rep); Y_nc  <- f_keep_y(Y_nc)

  # Tiny plotting helper using ggheat() – no annotations, fixed palettes
  heat_div <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100)
  heat_seq <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(100)
  br_div   <- seq(-3, 3, length.out = length(heat_div) + 1)
  br_seq   <- seq(0, 6, length.out = length(heat_seq) + 1)

  plot_mat <- function(M, pdf_path, diverging = TRUE) {
    if (!nrow(M) || !ncol(M)) return(invisible(NULL))
    # order rows by variance for nicer look
    M <- M[order(matrixStats::rowVars(M, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
    ggheat(
      M,
      scale        = "none",
      color        = if (diverging) heat_div else heat_seq,
      breaks       = if (diverging) br_div   else br_seq,
      border_color = NA,
      show_rownames = FALSE,
      filename     = pdf_path,
      silent       = TRUE,
      width        = 10,
      height       = 16
    )
    message("Saved: ", pdf_path)
  }

  # Write 6 PDFs
  plot_mat(X_act, file.path(out_dir, "hub_X_activate_heatmap.pdf"),   TRUE)
  plot_mat(X_rep, file.path(out_dir, "hub_X_repress_heatmap.pdf"),    TRUE)
  plot_mat(X_nc,  file.path(out_dir, "hub_X_noncancel_heatmap.pdf"),  TRUE)

  plot_mat(Y_act, file.path(out_dir, "hub_Y_activate_heatmap.pdf"),   FALSE)
  plot_mat(Y_rep, file.path(out_dir, "hub_Y_repress_heatmap.pdf"),    FALSE)
  plot_mat(Y_nc,  file.path(out_dir, "hub_Y_noncancel_heatmap.pdf"),  FALSE)

  invisible(list(
    matrices = list(X_act=X_act, X_rep=X_rep, X_nc=X_nc, Y_act=Y_act, Y_rep=Y_rep, Y_nc=Y_nc),
    out_dir = out_dir
  ))
}

# -------------------------------------------------------------------
# Shared heatmap palette (your palette)
# -------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(matrixStats); library(stringr)
})


# Shared palettes (your palette for the heatmap)
heat_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100)
brks      <- seq(-3, 3, length.out = length(heat_cols) + 1)

# Annotation palettes
time_pal <- c(`1h`="#4daf4a", `4h`="#377eb8", `24h`="#984ea3", `6d`="#e41a1c")
cond_pal <- c(BG="#1b9e77", LPS="#d95f02", culture="#7570b3")

.trim_ws    <- function(x) gsub("^\\s+|\\s+$", "", x)
.left_of_vs <- function(x) sub("\\s+vs\\s+.*$", "", x)

derive_time_from_label <- function(label_chr) {
  g <- .trim_ws(label_chr)
  g <- stringr::str_replace_all(g, "_+", "_")
  out <- rep(NA_character_, length(g))

  out[grepl("(^|_)Ctrl($|_)", g, ignore.case = TRUE)] <- "0h"
  out[grepl("culture_5d", g, ignore.case = TRUE)]     <- "6d"

  idx_6d <- grepl("(^|_)(6d|d6)($|_)", g, ignore.case = TRUE)
  out[idx_6d & is.na(out)] <- "6d"

  has_h <- grepl("(^|_)([0-9]+)h($|_)", g, ignore.case = TRUE)
  if (any(has_h)) {
    hh <- sub("^.*(^|_)([0-9]+)h($|_).*?$", "\\2", g[has_h], perl = TRUE)
    fill <- has_h & is.na(out)
    out[fill] <- paste0(hh[match(which(fill), which(has_h))], "h")
  }

  has_d <- grepl("(^|_)d([0-9]+)($|_)", g, ignore.case = TRUE)
  if (any(has_d)) {
    dd <- as.integer(sub("^.*(^|_)d([0-9]+)($|_).*?$", "\\2", g[has_d], perl = TRUE))
    conv <- ifelse(dd == 1L, "24h", paste0(dd, "d"))
    fill <- has_d & is.na(out)
    out[fill] <- conv[match(which(fill), which(has_d))]
  }
  out
}

derive_condition_from_label <- function(label_chr) {
  g <- .trim_ws(label_chr)
  g <- stringr::str_replace_all(g, "_+", "_")
  out <- rep(NA_character_, length(g))

  out[is.na(out) & grepl("^BG(_|$)",    g, ignore.case = TRUE)] <- "BG"
  out[is.na(out) & grepl("^LPS(_|$)",   g, ignore.case = TRUE)] <- "LPS"
  out[is.na(out) & (grepl("^RPMI(_|$)", g, ignore.case = TRUE) | grepl("^culture(_|$)", g, ignore.case = TRUE))] <- "culture"
  out[is.na(out) & grepl("(^|_)Ctrl($|_)", g, ignore.case = TRUE)] <- "Ctrl"  # not used in Condition row, but safe
  out
}

make_col_ann_time_cond <- function(comparison_cols, tag) {
  cond1 <- .left_of_vs(comparison_cols)
  time_vec <- derive_time_from_label(cond1)
  cond_vec <- derive_condition_from_label(cond1)

  # For vs_culture set, hide "culture" (culture is the control)
  if (identical(tag, "vs_culture")) {
    cond_vec[cond_vec == "culture"] <- NA_character_
  }

  # For vs_Ctrl, show culture/LPS/BG as-is
  ann <- data.frame(
    row.names = comparison_cols,
    Time      = time_vec,
    Condition = cond_vec,
    stringsAsFactors = FALSE
  )
  ann
}

plot_hub_heatmaps_simple <- function(
    summary_files,
    out_dir  = "/data/homes/yl814/episcope_test/GSE87218_ATAC/heatmap_pca",
    tag      = c("vs_culture","vs_Ctrl")[1],
    x_thresh = 1,   # keep TFs if any |x_signed| > x_thresh
    y_thresh = 1    # keep TFs if any y_links_log2p1 > y_thresh
){
  stopifnot(length(summary_files) > 0)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Informative subfolder
  subdir  <- sprintf("%s__nComp%d__xgt%s_ygt%s", tag, length(summary_files), x_thresh, y_thresh)
  out_dir <- file.path(out_dir, subdir)
  dir.create(out_dir, showWarnings = FALSE)

  res <- extract_tf_hub_axes(summary_files)

  pivot_mat <- function(df, value_col) {
    wide <- tidyr::pivot_wider(
      df[, c("comparison","TF", value_col)],
      id_cols     = "TF",
      names_from  = "comparison",
      values_from = value_col,
      values_fill = 0
    )
    mat <- as.matrix(wide[, setdiff(names(wide), "TF"), drop = FALSE])
    rownames(mat) <- wide$TF
    storage.mode(mat) <- "numeric"
    mat
  }

  # Build TF × comparison matrices
  X_act <- pivot_mat(res$x_activate,  "x_signed")
  X_rep <- pivot_mat(res$x_repress,   "x_signed")
  X_nc  <- pivot_mat(res$x_noncancel, "x_signed")

  Y_act <- pivot_mat(res$y_activate,  "y_links_log2p1")
  Y_rep <- pivot_mat(res$y_repress,   "y_links_log2p1")
  Y_nc  <- pivot_mat(res$y_noncancel, "y_links_log2p1")

  # Row filters (per matrix)
  keep_x <- function(M) M[matrixStats::rowMaxs(abs(M), na.rm = TRUE) > x_thresh, , drop = FALSE]
  keep_y <- function(M) M[matrixStats::rowMaxs(M,       na.rm = TRUE) > y_thresh, , drop = FALSE]

  X_act <- keep_x(X_act); X_rep <- keep_x(X_rep); X_nc  <- keep_x(X_nc)
  Y_act <- keep_y(Y_act); Y_rep <- keep_y(Y_rep); Y_nc  <- keep_y(Y_nc)

  # Smaller plot size + larger fonts so it reads when zoomed out
  # height scales gently with #TF for readability but stays compact
  dyn_h <- function(nr) max(10, 0.15 * nr + 6)  # smaller than before
  base_width <- 6

  # Write a PDF + CSV with plotted order
  write_plot_and_csv <- function(M, ann, title_text, base_name) {
    if (!nrow(M) || !ncol(M)) {
      message("Skip (empty after filter): ", base_name)
      return(invisible(NULL))
    }
    # Order rows by variance for nicer visuals before clustering
    M <- M[order(matrixStats::rowVars(M, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]

    pdf_path <- file.path(out_dir, sprintf("%s__%s__TF%d__Comp%d.pdf", tag, base_name, nrow(M), ncol(M)))
    csv_path <- sub("\\.pdf$", ".csv", pdf_path)

    ph <- ggheat(
      M,
      annotation_col    = ann,
      annotation_colors = list(Time = time_pal, Condition = cond_pal),
      scale             = "none",
      color             = heat_cols,
      breaks            = brks,
      cluster_cols      = FALSE,
      clustering_method = "ward.D2",
      distfun           = function(x) stats::dist(x, method = "euclidean"),
      border_color      = NA,
      show_rownames     = TRUE,     # TF labels on rows
      show_colnames     = TRUE,     # comparison labels on columns
      main              = sprintf("%s | %s | TF=%d, Comp=%d | x>|%s| or y>%s",
                                  tag, title_text, nrow(M), ncol(M), x_thresh, y_thresh),
      filename          = pdf_path,
      silent            = TRUE,
      width             = base_width,
      height            = dyn_h(nrow(M)),
      fontsize_row      = 12,       # bigger row text
      fontsize_col      = 12,       # bigger column text
      fontsize          = 12        # overall elements (title, legend)
      # Note: true bold depends on ggheat/pheatmap support; fontsize ↑ for readability
    )
    message("Saved: ", pdf_path)

    # Extract plotted order for CSV
    row_ord <- tryCatch(ph$tree_row$order, error = function(...) NULL)
    col_ord <- tryCatch(ph$tree_col$order, error = function(...) NULL)
    if (is.null(row_ord)) row_ord <- seq_len(nrow(M))
    if (is.null(col_ord)) col_ord <- seq_len(ncol(M))

    ordered_rows <- rownames(M)[row_ord]
    ordered_cols <- colnames(M)[col_ord]

    out_df <- as.data.frame(M) |>
      tibble::rownames_to_column("TF") |>
      dplyr::slice(match(ordered_rows, TF)) |>
      dplyr::select(TF, dplyr::all_of(ordered_cols))

    write.csv(out_df, file = csv_path, row.names = FALSE)
    message("Saved: ", csv_path)
    invisible(list(pdf = pdf_path, csv = csv_path))
  }


  make_ann_and_plot <- function(M, title_text, base_name) {
    cols <- colnames(M)
    ann  <- make_col_ann_time_cond(cols, tag)
    # align columns with annotation
    ann <- ann[intersect(rownames(ann), colnames(M)), , drop = FALSE]
    # sort columns by Condition then Time
    cond_levels <- c("BG", "LPS", "culture")
    time_levels <- c("1h", "4h", "24h", "6d")

    ann$Condition <- factor(ann$Condition, levels = cond_levels, ordered = TRUE)
    ann$Time      <- factor(ann$Time,      levels = time_levels, ordered = TRUE)

    ord <- order(ann$Condition, ann$Time, na.last = TRUE)
    ann <- ann[ord, , drop = FALSE]

    # apply sorted columns to the matrix
    M <- M[, rownames(ann), drop = FALSE]
    # --------------------------------------------------------------------
    write_plot_and_csv(M, ann, title_text, base_name)
  }


  # Plot all six
  make_ann_and_plot(X_act, "hub_X_activate (signed)",                               "hub_X_activate")
  make_ann_and_plot(X_rep, "hub_X_repress (signed)",                                "hub_X_repress")
  make_ann_and_plot(X_nc,  "hub_X_noncancel (sign=overall, mag=|act|+|rep|)",      "hub_X_noncancel")

  make_ann_and_plot(Y_act, "hub_Y_activate (log2(links+1))",                        "hub_Y_activate")
  make_ann_and_plot(Y_rep, "hub_Y_repress (log2(links+1))",                         "hub_Y_repress")
  make_ann_and_plot(Y_nc,  "hub_Y_noncancel (overall log2(links+1))",               "hub_Y_noncancel")

  invisible(list(out_dir = out_dir, n_comp = length(summary_files)))
}


# vs culture
plot_hub_heatmaps_simple(
  summary_files = summary_files_vs_culture,
  out_dir = "/data/homes/yl814/episcope_test/GSE87218_ATAC/heatmap_pca",
  tag = "vs_culture",
  x_thresh = 1, y_thresh = 1
)

# vs Ctrl
plot_hub_heatmaps_simple(
  summary_files = summary_files_vs_Ctrl,
  out_dir = "/data/homes/yl814/episcope_test/GSE87218_ATAC/heatmap_pca",
  tag = "vs_Ctrl",
  x_thresh = 1, y_thresh = 1
)
