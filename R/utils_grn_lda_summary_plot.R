#' utils_grn_lda_summary_plot.R ?Plot from LDA summary tables (episcope)
#'
#' @author Yaoxiang Li
#' @family episcope-plots
#'
#' @description
#' Plot TFtopic bubble panels and TF hub scatterplots from
#' \code{*_filtered_lda_K{K}_summary.csv} files produced by
#' \code{summarize_lda_annotations_bulk()}.
#'
#' @section Outputs (per summary CSV):
#' \itemize{
#'   \item Bubbles (size = log2(expr max), color = signed log2(sum)):
#'         \code{*_bubble_byExpr_expr.pdf},
#'         \code{*_bubble_byExpr_expr_activate.pdf},
#'         \code{*_bubble_byExpr_expr_repress.pdf},
#'         \code{*_bubble_byExpr_expr_nocancel.pdf}
#'   \item Hubs (size = log2(expr max), fill = z(log2FC)):
#'         \code{<input>_summary.tf_hubs.pdf},
#'         \code{<input>_summary.tf_hubs_activate.pdf},
#'         \code{<input>_summary.tf_hubs_repress.pdf},
#'         \code{<input>_summary.tf_hubs_noncancel.pdf}
#' }
#'
#' @details
#' Bubble color is the signed log2 of the chosen delta-link sum (overall/activate/repress).
#' If \code{tf_expr_max} is unavailable, bubbles fall back to size by log2(#links).
#' Titles include the contrast parsed from the summary filename.
#'
#' @keywords internal
#' @noRd


# =============================
# Helpers & I/O (internal)
# =============================

.plog <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message("[utils_grn_lda_summary_plot] ", paste0(..., collapse = ""))
}

.load_summary <- function(x, verbose = TRUE) {
  if (is.character(x) && length(x) == 1L) {
    .plog("Reading CSV: ", x, verbose = verbose)
    readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .plog("Using in-memory tibble/data.frame (nrow=", nrow(x), ")", verbose = verbose)
    tibble::as_tibble(x)
  } else {
    cli::cli_abort("`x` must be a CSV path or a data.frame/tibble.")
  }
}

.path_swap_suffix <- function(path, old_suffix, new_suffix) {
  sub(paste0(old_suffix, "$"), new_suffix, path, perl = TRUE)
}

.path_add_postfix <- function(path, postfix, new_ext = ".pdf") {
  noext <- sub("\\.[^.]+$", "", path)
  paste0(noext, ".", postfix, new_ext)
}

.contrast_from_file <- function(f) {
  b <- basename(f)
  stem <- sub("_delta_links.*$", "", b)
  parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) == 2) paste(parts[1], "vs", parts[2]) else stem
}

.contrast_parts <- function(f) {
  b <- basename(f)
  stem <- sub("_delta_links.*$", "", b)
  parts <- strsplit(stem, "_vs_", fixed = TRUE)[[1]]
  cond1 <- if (length(parts) >= 1) parts[1] else NA_character_
  cond2 <- if (length(parts) >= 2) parts[2] else NA_character_
  list(cond1 = cond1, cond2 = cond2, label = if (length(parts) == 2) paste(cond1, "vs", cond2) else stem)
}

robust_z <- function(x) {
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

clamp <- function(x, lo, hi) pmax(pmin(x, hi), lo)


# ==========================================================
#' Expr-only bubble plots (overall / activate / repress)
#'
#' @param summary_df Tibble from \code{summarize_lda_annotations()}.
#' @param out_pdf Output PDF path.
#' @param variant One of \code{"overall"}, \code{"activate"}, \code{"repress"}.
#' @param top_tf_per_topic Top N TFs per topic (by expression) to show.
#' @param facet_ncol Number of columns in facet grid.
#' @param point_max_size Max bubble size (ggplot2::scale_size_area).
#' @param base_size Base text size.
#' @param color_cap Optional symmetric cap for color scale of signed log2 sum.
#' @param panel_height_base Unused (kept for API compatibility).
#' @param panel_height_per_tf Unused (kept for API compatibility).
#' @param main_title_prefix Title prefix.
#' @param color_sigma Clamp for z-scored color (TF log2FC).
#' @param no_cancel If TRUE, magnitude is sum of abs parts; sign from overall.
#' @param verbose Verbose logging.
#'
#' @return Invisibly returns \code{out_pdf}.
#' @export
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_x_continuous scale_size_area
#' @importFrom ggplot2 scale_fill_gradient2 facet_wrap scale_y_discrete labs guides guide_legend
#' @importFrom ggplot2 guide_colorbar theme_minimal theme element_blank element_rect element_text ggsave
#' @importFrom ggplot2 theme_classic element_text
#' @importFrom ggplot2 expansion
#' @importFrom dplyr mutate filter group_by slice_max ungroup distinct arrange desc pull select summarise
#' @importFrom tidyr separate_rows
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales log_trans rescale squish
make_expr_bubbles_from_summary <- function(
    summary_df, out_pdf,
    variant = c("overall","activate","repress"),
    top_tf_per_topic = 20,
    facet_ncol = 5,
    point_max_size = 7,
    base_size = 10,
    color_cap = NULL,
    panel_height_base = 0.6,
    panel_height_per_tf = 0.12,
    main_title_prefix = "Top TFs per topic (X = signed log2 sum delta_edge)",
    color_sigma = 2,
    no_cancel = FALSE,
    verbose = TRUE
){
  variant <- match.arg(variant)

  need <- c("TF","topic","topic_rank","tf_expr_max","tf_log2_fc",
            "tf_delta_sum","tf_delta_sum_activate","tf_delta_sum_repress",
            "tf_delta_sum_abs_activate","tf_delta_sum_abs_repress")
  if (!all(need %in% names(summary_df))) {
    cli::cli_abort("Summary is missing required columns. Have: {paste(names(summary_df), collapse=', ')}")
  }

  df0 <- dplyr::mutate(
    summary_df,
    topic = as.integer(topic),
    tf_log2_expr_max = log2(pmax(tf_expr_max, 1e-9)),
    tf_log2_fc       = suppressWarnings(as.numeric(tf_log2_fc))
  )

  sum_col <- switch(variant,
                    overall  = "tf_delta_sum",
                    activate = "tf_delta_sum_activate",
                    repress  = "tf_delta_sum_repress")
  df0$`._sum_raw` <- dplyr::coalesce(df0[[sum_col]], 0)
  df0$tf_signed_log_sum <- sign(df0$`._sum_raw`) * log2(abs(df0$`._sum_raw`) + 1)

  if (isTRUE(no_cancel)) {
    if (!all(c("tf_delta_sum_abs_activate","tf_delta_sum_abs_repress","tf_delta_sum") %in% names(df0))) {
      cli::cli_abort("no_cancel requires tf_delta_sum_abs_{activate,repress} and tf_delta_sum.")
    }
    if (variant == "overall") {
      mag  <- dplyr::coalesce(df0$tf_delta_sum_abs_activate, 0) + dplyr::coalesce(df0$tf_delta_sum_abs_repress, 0)
      sign_overall <- sign(dplyr::coalesce(df0$tf_delta_sum, 0))
      df0$tf_signed_log_sum <- sign_overall * log2(mag + 1)
    } else if (variant == "activate") {
      mag  <- dplyr::coalesce(df0$tf_delta_sum_abs_activate, 0)
      sign_overall <- sign(dplyr::coalesce(df0$tf_delta_sum, 0))
      df0$tf_signed_log_sum <- sign_overall * log2(mag + 1)
    } else {
      mag  <- dplyr::coalesce(df0$tf_delta_sum_abs_repress, 0)
      sign_overall <- sign(dplyr::coalesce(df0$tf_delta_sum, 0))
      df0$tf_signed_log_sum <- sign_overall * log2(mag + 1)
    }
  }

  z  <- robust_z(df0$tf_log2_fc)
  df0$tf_log2_fc_z_c <- clamp(z, -abs(color_sigma), abs(color_sigma))

  df_top <- df0 |>
    dplyr::filter(is.finite(tf_log2_expr_max)) |>
    dplyr::group_by(topic) |>
    dplyr::slice_max(order_by = tf_log2_expr_max, n = top_tf_per_topic, with_ties = FALSE) |>
    dplyr::ungroup()
  if (!nrow(df_top)) {
    warning("No rows to plot after top-per-topic selection: ", out_pdf)
    return(invisible(NULL))
  }

  rng <- range(df_top$tf_log2_expr_max, na.rm = TRUE, finite = TRUE)
  dlt <- diff(rng)
  if (!is.finite(dlt) || dlt == 0) {
    df_top$size_scaled <- 0.6
    size_breaks <- c(0.6)
    size_labs   <- sprintf("%.1f", rng[1])
  } else {
    df_top$size_scaled <- (df_top$tf_log2_expr_max - rng[1]) / dlt
    size_breaks <- c(0, 0.5, 1)
    size_labs   <- sprintf("%.1f", rng[1] + size_breaks * dlt)
  }

  topic_order <- df_top |>
    dplyr::distinct(topic, topic_rank) |>
    dplyr::arrange(topic_rank, topic) |>
    dplyr::pull(topic)

  y_order <- df_top |>
    dplyr::group_by(topic) |>
    dplyr::arrange(dplyr::desc(abs(tf_signed_log_sum)),
                   dplyr::desc(tf_log2_expr_max),
                   TF,
                   .by_group = TRUE) |>
    dplyr::mutate(TF_key = paste0(topic, "::", TF)) |>
    dplyr::pull(TF_key) |>
    unique()

  plot_df <- dplyr::mutate(
    df_top,
    topic_fac = factor(topic, levels = topic_order, ordered = TRUE),
    TF_fac    = factor(paste0(topic, "::", TF), levels = rev(y_order), ordered = TRUE)
  )

  cap <- if (is.null(color_cap)) max(abs(plot_df$tf_signed_log_sum), na.rm = TRUE) else color_cap
  if (!is.finite(cap) || cap <= 0) cap <- 1

  facet_labs_tbl <- plot_df |>
    dplyr::distinct(topic_fac, topic, topic_rank) |>
    dplyr::arrange(topic_rank, topic)
  facet_map <- setNames(
    paste0("Topic ", facet_labs_tbl$topic, " | Rank #", facet_labs_tbl$topic_rank),
    facet_labs_tbl$topic_fac
  )
  lab_fn <- function(x) {
    nm <- as.character(x)
    y  <- unname(facet_map[nm])
    y[is.na(y)] <- paste("Topic", nm[is.na(y)])
    y
  }

  x_lab <- if (isTRUE(no_cancel)) {
    "signed log2(sum abs(delta_edge)) per TF  (sign from OVERALL cond1-cond2)"
  } else {
    "signed log2(sum delta_edge) per TF"
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = tf_signed_log_sum, y = TF_fac)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    ggplot2::geom_point(
      ggplot2::aes(size = size_scaled, fill = tf_log2_fc_z_c),
      shape = 21, color = "grey10", stroke = 0.4, alpha = 0.95
    ) +
    ggplot2::scale_x_continuous(
      name   = x_lab,
      limits = c(-cap, cap),
      expand = ggplot2::expansion(mult = c(0.02, 0.06))
    ) +
    ggplot2::scale_size_area(
      max_size = point_max_size, limits = c(0, 1),
      breaks = size_breaks, labels = size_labs,
      name = "size = log2(expr max)"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-abs(color_sigma), abs(color_sigma)),
      oob = scales::squish, name = "TF log2FC (z)"
    ) +
    ggplot2::facet_wrap(
      ~ topic_fac, ncol = facet_ncol, scales = "free_y",
      labeller = ggplot2::labeller(topic_fac = lab_fn)
    ) +
    ggplot2::scale_y_discrete(labels = function(keys) sub("^\\d+::", "", keys)) +
    ggplot2::labs(y = "TF", title = main_title_prefix) +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, stroke = 0.4)),
      fill = ggplot2::guide_colorbar(order = 2)
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background   = ggplot2::element_rect(fill = "#eeeeee", color = NA),
      legend.position    = "bottom",
      plot.title         = ggplot2::element_text(hjust = 0.5),
      axis.title.y       = ggplot2::element_text(face = "bold"),
      axis.title.x       = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(out_pdf, p, width = 12, height = 9, dpi = 300)
  invisible(out_pdf)
}


# =============================
#' TF hub scatterplot (K-independent)
#'
#' @param summary_df Tibble from \code{summarize_lda_annotations()}.
#' @param out_pdf Output PDF path.
#' @param subset_mode One of \code{"overall"}, \code{"activate"}, \code{"repress"}.
#' @param x_sum_mode One of \code{"signed"}, \code{"abs"}, \code{"signed_noncancel"}.
#' @param width_in,height_in,dpi,base_size Plot sizing.
#' @param color_sigma Clamp for z-scored color (TF log2FC).
#' @param label_top_by_links,label_top_by_sum Integers for label seeding.
#' @param label_quantile_links,label_quantile_absx Quantiles for auto-labels.
#' @param label_cap Max number of labels.
#' @param min_x_for_label,min_y_for_label Filters for labeling.
#' @param label_size Label text size.
#' @param title_text Optional custom title.
#' @param cond1_label,cond2_label Optional condition labels (caption).
#' @param verbose Verbose logging.
#'
#' @return Invisibly returns \code{out_pdf}.
#' @export
#' @importFrom ggplot2 ggplot aes geom_vline geom_point labs scale_y_continuous scale_fill_gradient2
#' @importFrom ggplot2 theme_classic theme element_text guide_colorbar guide_legend guides ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales log_trans rescale squish
#' @importFrom dplyr group_by summarise mutate case_when filter pull arrange desc ungroup
make_tf_hubs_from_summary <- function(
    summary_df, out_pdf,
    subset_mode = c("overall","activate","repress"),
    x_sum_mode  = c("signed","abs","signed_noncancel"),
    width_in = 16, height_in = 9, dpi = 220, base_size = 18,
    color_sigma = 2,
    label_top_by_links = 150, label_top_by_sum = 150,
    label_quantile_links = 0.70, label_quantile_absx = 0.70,
    label_cap = 220, min_x_for_label = 0, min_y_for_label = 0,
    label_size = 2.0,
    title_text = NULL, cond1_label = NULL, cond2_label = NULL,
    verbose = TRUE
){
  subset_mode <- match.arg(subset_mode)
  x_sum_mode  <- match.arg(x_sum_mode)

  need <- c(
    "TF","tf_expr_max","tf_log2_fc",
    "tf_n_links","tf_n_links_activate","tf_n_links_repress",
    "tf_delta_sum","tf_delta_sum_activate","tf_delta_sum_repress",
    "tf_delta_sum_abs","tf_delta_sum_abs_activate","tf_delta_sum_abs_repress"
  )
  if (!all(need %in% names(summary_df))) {
    cli::cli_abort("Summary is missing required columns. Have: {paste(names(summary_df), collapse=', ')}")
  }

  c1 <- if (is.null(cond1_label) || !nzchar(cond1_label)) "cond1" else cond1_label
  c2 <- if (is.null(cond2_label) || !nzchar(cond2_label)) "cond2" else cond2_label

  col_links      <- switch(subset_mode,
                           overall  = "tf_n_links",
                           activate = "tf_n_links_activate",
                           repress  = "tf_n_links_repress")
  col_sum_abs    <- switch(subset_mode,
                           overall  = "tf_delta_sum_abs",
                           activate = "tf_delta_sum_abs_activate",
                           repress  = "tf_delta_sum_abs_repress")
  col_sum_signed <- switch(subset_mode,
                           overall  = "tf_delta_sum",
                           activate = "tf_delta_sum_activate",
                           repress  = "tf_delta_sum_repress")

  TF <- summary_df |>
    dplyr::group_by(TF) |>
    dplyr::summarise(
      tf_links           = sum(.data[[col_links]],      na.rm = TRUE),
      tf_sum_abs_delta   = sum(.data[[col_sum_abs]],    na.rm = TRUE),
      tf_sum_delta       = sum(.data[[col_sum_signed]], na.rm = TRUE),
      tf_sum_abs_act     = sum(dplyr::coalesce(tf_delta_sum_abs_activate, 0), na.rm = TRUE),
      tf_sum_abs_rep     = sum(dplyr::coalesce(tf_delta_sum_abs_repress,  0), na.rm = TRUE),
      tf_sum_abs_both    = sum(dplyr::coalesce(tf_delta_sum_abs_activate, 0) +
                                 dplyr::coalesce(tf_delta_sum_abs_repress,  0), na.rm = TRUE),
      tf_sum_delta_overall = sum(dplyr::coalesce(tf_delta_sum, 0), na.rm = TRUE),
      tf_expr_max        = max(tf_expr_max, na.rm = TRUE),
      tf_log2_fc_med     = stats::median(tf_log2_fc, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tf_links_plus1   = tf_links + 1L,
      tf_log2_expr_max = log2(pmax(tf_expr_max, 1e-9))
    )

  if (x_sum_mode == "abs") {
    TF$x_sum <- log2(TF$tf_sum_abs_delta + 1e-8)
  } else if (x_sum_mode == "signed") {
    TF$x_sum <- sign(TF$tf_sum_delta) * log2(abs(TF$tf_sum_delta) + 1)
  } else {
    sign_nc <- sign(TF$tf_sum_delta_overall)
    mag <- if (subset_mode == "activate") TF$tf_sum_abs_act else if (subset_mode == "repress") TF$tf_sum_abs_rep else TF$tf_sum_abs_both
    TF$x_sum <- sign_nc * log2(mag + 1)
  }
  TF$abs_x_for_label <- if (x_sum_mode == "abs") TF$x_sum else abs(TF$x_sum)

  TF$log2_fc_z_c <- clamp(robust_z(TF$tf_log2_fc_med), -abs(color_sigma), abs(color_sigma))

  by_links <- TF$TF[order(-TF$tf_links)][seq_len(min(label_top_by_links, nrow(TF)))]
  by_sum   <- TF$TF[order(-TF$tf_sum_abs_delta)][seq_len(min(label_top_by_sum, nrow(TF)))]
  q_links  <- stats::quantile(TF$tf_links,        probs = label_quantile_links, na.rm = TRUE)
  q_absx   <- stats::quantile(TF$abs_x_for_label, probs = label_quantile_absx,  na.rm = TRUE)
  by_q_links <- TF$TF[TF$tf_links >= q_links]
  by_q_absx  <- TF$TF[TF$abs_x_for_label >= q_absx]
  lab_ids <- unique(c(by_links, by_sum, by_q_links, by_q_absx))
  if (length(lab_ids) > label_cap) {
    TF$`._score` <- scales::rescale(TF$tf_links) + scales::rescale(TF$abs_x_for_label)
    lab_ids <- TF$TF[order(-TF$`._score`)][seq_len(label_cap)]
  }
  lab_df <- dplyr::filter(TF, TF %in% lab_ids,
                          abs_x_for_label >= min_x_for_label,
                          tf_links >= min_y_for_label)

  ymax <- max(TF$tf_links_plus1, na.rm = TRUE); if (!is.finite(ymax) || ymax <= 1) ymax <- 2
  kmax <- ceiling(log(ymax, base = 2))
  y_breaks <- 2^(0:kmax)
  y_labels <- pmax(y_breaks - 1, 0)

  have_expr <- any(is.finite(TF$tf_log2_expr_max))
  if (have_expr) {
    size_aes_name <- "tf_log2_expr_max"
    size_title    <- "expr max (log2)"
    TF$size_val   <- TF$tf_log2_expr_max
  } else {
    size_aes_name <- "size_fallback"
    size_title    <- "links (log2)"
    TF$size_val   <- log2(TF$tf_links + 1)
  }

  x_title <- switch(
    x_sum_mode,
    abs              = "log2(sum abs(delta link_score)) per TF",
    signed           = "signed sum delta link_score per TF (log2-scaled magnitude)",
    signed_noncancel = "signed log2(sum abs(delta link_score)) per TF  (sign from OVERALL cond1-cond2)"
  )
  if (is.null(title_text)) {
    title_text <- paste0("TF hubs  - ",
                         switch(subset_mode, overall="overall",
                                activate="activate-only", repress="repress-only"))
  }
  caption_text <- paste0(
    "delta link_score = link_score(", c1, ") - link_score(", c2, ").\n",
    "For 'signed_noncancel', sign = OVERALL (cond1-cond2), magnitude = subset sum|delta|.\n",
    "Y: log2(1 + #links). Size: log2(max TF RNA across ", c1, " & ", c2, "). Color: TF log2FC z-score."
  )

  p <- ggplot2::ggplot(TF, ggplot2::aes(x = x_sum, y = tf_links_plus1)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(size = size_val, fill = log2_fc_z_c),
                        shape = 21, color = "grey10", stroke = 0.4, alpha = 0.8) +
    ggrepel::geom_text_repel(
      data = lab_df, ggplot2::aes(label = TF),
      size = label_size, fontface = "bold",
      max.overlaps = Inf, box.padding = 0.18, point.padding = 0.12,
      min.segment.length = 0, segment.alpha = 0.5
    ) +
    ggplot2::scale_size_area(max_size = 8, name = size_title) +
    ggplot2::scale_fill_gradient2(
      low = "#4575b4", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-abs(color_sigma), abs(color_sigma)),
      oob = scales::squish, name = "TF log2FC (z)"
    ) +
    ggplot2::labs(
      title = title_text,
      x = x_title,
      y = "number of differential links per TF in [log2(1 + count)] scale",
      caption = caption_text
    ) +
    # BEGIN EDIT: use ggplot2::expansion, not scales::expansion
    ggplot2::scale_y_continuous(
      trans = scales::log_trans(base = 2),
      breaks = y_breaks, labels = y_labels,
      expand = ggplot2::expansion(mult = c(0.02, 0.06))
    ) +
    # END EDIT
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, override.aes = list(alpha = 1, stroke = 0.4)),
      fill = ggplot2::guide_colorbar(order = 2)
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = base_size - 4),
      legend.text  = ggplot2::element_text(size = base_size - 6),
      plot.title   = ggplot2::element_text(hjust = 0.5),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(out_pdf, p, width = width_in, height = height_in, dpi = dpi, limitsize = FALSE)
  invisible(out_pdf)
}


# =============================
#' Bulk driver over many summary CSVs
#'
#' @param summary_csvs Character vector of summary CSV paths.
#' @param top_tf_per_topic Integer; top TFs per topic for bubble plots.
#' @param pdf_width,pdf_height Bubble PDF size in inches.
#' @param hubs_width_in,hubs_height_in,hubs_dpi,hubs_base Hub plot sizing.
#' @param color_cap Optional symmetric cap for bubble color scale.
#' @param color_sigma Clamp for hub color (z-score).
#' @param expr_point_max_size Max bubble size.
#' @param panel_height_base,panel_height_per_tf Kept for API compatibility.
#' @param parallel Use future.apply.
#' @param plan Future plan string or function.
#' @param workers Integer; default \code{parallel::detectCores()-1}.
#' @param verbose Verbose logging.
#'
#' @return Invisibly returns a list with \code{summary_csvs} and unique \code{dirs}.
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom parallel detectCores
plot_from_summary_bulk <- function(
    summary_csvs,
    top_tf_per_topic = 20,
    pdf_width = 12, pdf_height = 9,
    hubs_width_in = 14, hubs_height_in = 9, hubs_dpi = 220, hubs_base = 18,
    color_cap = NULL, color_sigma = 2,
    expr_point_max_size = 7,
    panel_height_base = 0.6,
    panel_height_per_tf = 0.12,
    parallel = TRUE, plan = "multisession", workers = NULL,
    verbose = TRUE
){
  if (!is.character(summary_csvs) || length(summary_csvs) == 0L)
    cli::cli_abort("`summary_csvs` must be a non-empty character vector of paths.")
  summary_csvs <- summary_csvs[file.exists(summary_csvs)]
  if (!length(summary_csvs)) cli::cli_abort("No existing summary CSVs found.")

  if (isTRUE(parallel)) {
    strategy <- (function(x) {
      if (is.character(x)) {
        ns <- asNamespace("future")
        if (exists(x, envir = ns, mode = "function")) get(x, envir = ns) else future::multisession
      } else if (is.function(x)) x else future::multisession
    })(plan)
    if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    .plog("Launching future plan=", if (is.character(plan)) plan else "custom",
          ", workers=", workers, verbose = verbose)
    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    future::plan(strategy, workers = workers)
  }

  .title_with_contrast <- function(f, base) paste0(base, " - ", .contrast_from_file(f))

  .one <- function(f) {
    sz <- suppressWarnings(file.info(f)$size)
    if (!is.finite(sz) || is.na(sz) || sz == 0) {
      .plog("Skipping empty summary CSV (0 bytes): ", f, verbose = verbose)
      return(invisible(FALSE))
    }

    S <- tryCatch(
      .load_summary(f, verbose = verbose),
      error = function(e) {
        .plog("Failed to read summary CSV, skipping: ", f, " | err=", conditionMessage(e), verbose = verbose)
        NULL
      }
    )
    if (is.null(S) || NROW(S) == 0L) {
      .plog("Skipping summary with 0 rows: ", f, verbose = verbose)
      return(invisible(FALSE))
    }

    # ----- Bubble plots (overall / activate / repress) -----
    make_expr_bubbles_from_summary(
      S,
      out_pdf = .path_swap_suffix(f, "_summary.csv", "_bubble_byExpr_expr.pdf"),
      variant = "overall",
      top_tf_per_topic = top_tf_per_topic,
      facet_ncol = 5,
      point_max_size = expr_point_max_size,
      base_size = 10,
      color_cap = color_cap,
      panel_height_base = panel_height_base,
      panel_height_per_tf = panel_height_per_tf,
      main_title_prefix = .title_with_contrast(f, "Top TFs per topic (X = signed log2 sum delta_edge) - overall")
    )
    make_expr_bubbles_from_summary(
      S,
      out_pdf = .path_swap_suffix(f, "_summary.csv", "_bubble_byExpr_expr_activate.pdf"),
      variant = "activate",
      top_tf_per_topic = top_tf_per_topic,
      facet_ncol = 5,
      point_max_size = expr_point_max_size,
      base_size = 10,
      color_cap = color_cap,
      panel_height_base = panel_height_base,
      panel_height_per_tf = panel_height_per_tf,
      main_title_prefix = .title_with_contrast(f, "Top TFs per topic (X = signed log2 sum delta_edge) - activate")
    )
    make_expr_bubbles_from_summary(
      S,
      out_pdf = .path_swap_suffix(f, "_summary.csv", "_bubble_byExpr_expr_repress.pdf"),
      variant = "repress",
      top_tf_per_topic = top_tf_per_topic,
      facet_ncol = 5,
      point_max_size = expr_point_max_size,
      base_size = 10,
      color_cap = color_cap,
      panel_height_base = panel_height_base,
      panel_height_per_tf = panel_height_per_tf,
      main_title_prefix = .title_with_contrast(f, "Top TFs per topic (X = signed log2 sum delta_edge) - repress")
    )

    # ----- Bubble (no-cancel, overall) -----
    make_expr_bubbles_from_summary(
      S,
      out_pdf = .path_swap_suffix(f, "_summary.csv", "_bubble_byExpr_expr_nocancel.pdf"),
      variant = "overall",
      no_cancel = TRUE,
      top_tf_per_topic = top_tf_per_topic,
      facet_ncol = 5,
      point_max_size = expr_point_max_size,
      base_size = 10,
      color_cap = color_cap,
      panel_height_base = panel_height_base,
      panel_height_per_tf = panel_height_per_tf,
      main_title_prefix = .title_with_contrast(
        f, "Top TFs per topic (X = signed log2 sum abs(delta_edge)) - overall (no-cancel)"
      )
    )

    # ----- Hub plots (overall / activate / repress / noncancel) -----
    contrast <- .contrast_from_file(f)
    cp <- .contrast_parts(f)

    make_tf_hubs_from_summary(
      S,
      out_pdf   = .path_add_postfix(f, "tf_hubs"),
      subset_mode = "overall", x_sum_mode = "signed",
      width_in = hubs_width_in, height_in = hubs_height_in, dpi = hubs_dpi, base_size = hubs_base,
      color_sigma = color_sigma,
      label_top_by_links = 150, label_top_by_sum = 150,
      title_text = paste0("TF hubs - ", contrast, " - overall"),
      cond1_label = cp$cond1, cond2_label = cp$cond2,
      verbose = verbose
    )
    make_tf_hubs_from_summary(
      S,
      out_pdf   = .path_add_postfix(f, "tf_hubs_activate"),
      subset_mode = "activate", x_sum_mode = "signed",
      width_in = hubs_width_in, height_in = hubs_height_in, dpi = hubs_dpi, base_size = hubs_base,
      color_sigma = color_sigma,
      label_top_by_links = 150, label_top_by_sum = 150,
      title_text = paste0("TF hubs - ", contrast, " - activate-only"),
      cond1_label = cp$cond1, cond2_label = cp$cond2,
      verbose = verbose
    )
    make_tf_hubs_from_summary(
      S,
      out_pdf   = .path_add_postfix(f, "tf_hubs_repress"),
      subset_mode = "repress", x_sum_mode = "signed",
      width_in = hubs_width_in, height_in = hubs_height_in, dpi = hubs_dpi, base_size = hubs_base,
      color_sigma = color_sigma,
      label_top_by_links = 150, label_top_by_sum = 150,
      title_text = paste0("TF hubs - ", contrast, " - repress-only"),
      cond1_label = cp$cond1, cond2_label = cp$cond2,
      verbose = verbose
    )
    make_tf_hubs_from_summary(
      S,
      out_pdf   = .path_add_postfix(f, "tf_hubs_noncancel"),
      subset_mode = "overall", x_sum_mode = "signed_noncancel",
      width_in = hubs_width_in, height_in = hubs_height_in, dpi = hubs_dpi, base_size = hubs_base,
      color_sigma = color_sigma,
      label_top_by_links = 150, label_top_by_sum = 150,
      title_text = paste0("TF hubs - ", contrast, " - overall (no-cancel)"),
      cond1_label = cp$cond1, cond2_label = cp$cond2,
      verbose = verbose
    )

    invisible(TRUE)
  }

  if (isTRUE(parallel)) {
    future.apply::future_lapply(summary_csvs, .one, future.seed = TRUE)
  } else {
    lapply(summary_csvs, .one)
  }

  invisible(list(summary_csvs = summary_csvs, dirs = unique(dirname(summary_csvs))))
}
