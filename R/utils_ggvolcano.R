#' Publication-ready volcano plots
#'
#' Creates a volcano plot to visualize differential expression results.
#' This function is highly configurable to suit publication standards.
#'
#' @param data A data frame containing test statistics. Requires at least columns
#'   for variable names, log2 fold changes, and p-values.
#' @param labels Column name or vector for variable names (used for labeling).
#' @param logFC_col Column name for log2 fold changes.
#' @param pval_col Column name for nominal or adjusted p-values.
#' @param x_limits Limits of the x-axis (default auto-calculated).
#' @param y_limits Limits of the y-axis (default auto-calculated).
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param pval_cutoff P-value cutoff for significance.
#' @param logFC_cutoff Log2 fold-change cutoff for significance.
#' @param cutoff_line List of options for cutoff lines (`type`, `color`, `width`).
#' @param point_aes List of aesthetic options for points
#'   (`size`, `shape`, `color`, `alpha`, optionally `color_gradient`, etc.).
#' @param label_aes List of aesthetic options for labels
#'   (`size`, `color`, `face`, `parse`).
#' @param legend_aes List of aesthetic options for legend
#'   (`labels`, `position`, `label_size`, `icon_size`).
#' @param highlight_col Optional column name to use for custom highlighting
#'   instead of p-value/logFC (see `use_significance`).
#' @param highlight_colors Named or unnamed vector of colors used when
#'   `use_significance = FALSE`. If named, names must match the levels of
#'   `highlight_col`; otherwise they are matched in order.
#' @param highlight_labels Optional vector of labels for legend when
#'   `use_significance = FALSE`.
#' @param use_significance Logical; if TRUE (default), highlight based on
#'   logFC/p-value thresholds. If FALSE, highlight groups purely by
#'   `highlight_col` and associated colors/labels.
#' @param shade_options List of options for shading regions in the plot.
#' @param connector_aes List of aesthetic options for connectors
#'   (`line_width`, `arrow_type`, `arrow_ends`, `arrow_length`, `line_color`,
#'   `direction`, `draw_arrowheads`).
#' @param gridlines List with logical values indicating whether to draw
#'   gridlines (`major`, `minor`).
#' @param plot_border Add a border for plot axes (`"partial"` or `"full"`).
#' @param border_width Width of the border.
#' @param border_color Color of the border.
#' @param horizontal_line Numeric value(s) for drawing horizontal line(s).
#' @param horizontal_line_aes List of aesthetic options for the horizontal
#'   line(s) (`type`, `color`, `width`).
#'
#' @return A \code{ggplot2} object representing the volcano plot.
#'
#' @examples
#' # Basic usage with significance-based highlighting
#' # data <- read.csv(system.file("extdata", "example.csv", package = "ggvolcano"))
#' # ggvolcano(data,
#' #   logFC_col = "log2FoldChange",
#' #   pval_col = "pvalue"
#' # )
#'
#' @import ggplot2
#' @import ggrepel
#' @export
ggvolcano <- function(data,
                      labels = "",
                      logFC_col,
                      pval_col,
                      x_limits = c(
                        min(data[[logFC_col]], na.rm = TRUE) - 1.5,
                        max(data[[logFC_col]], na.rm = TRUE) + 1.5
                      ),
                      y_limits = c(
                        0,
                        max(-log10(data[[pval_col]]), na.rm = TRUE) + 5
                      ),
                      xlab = bquote( ~ Log[2] ~ "fold change"),
                      ylab = bquote( ~ -Log[10] ~ italic(P)),
                      title = "Volcano plot",
                      subtitle = "",
                      caption = paste0("total = ", nrow(data), " variables"),
                      pval_cutoff = 1e-6,
                      logFC_cutoff = 1.0,
                      cutoff_line = list(
                        type  = "longdash",
                        color = "black",
                        width = 0.4
                      ),
                      point_aes = list(
                        size  = 0.4,
                        shape = c(19, 19, 19, 19),
                        color = c("grey30", "#00CD6C", "#009ADE", "#FF1F5B"),
                        alpha = 0.7
                      ),
                      label_aes = list(
                        size  = 2.5,
                        color = "black",
                        face  = "plain",
                        parse = FALSE
                      ),
                      legend_aes = list(
                        labels = c(
                          "NS",
                          expression(Log[2] ~ FC),
                          "p-value",
                          expression(p - value ~ and ~ log[2] ~ FC)
                        ),
                        position   = "right",
                        label_size = 14,
                        icon_size  = 5
                      ),
                      highlight_col    = NULL,
                      highlight_colors = c("0" = "grey30", "1" = "red"),
                      highlight_labels = NULL,
                      use_significance = TRUE,
                      shade_options = NULL,
                      connector_aes = list(
                        line_width      = 0.5,
                        arrow_type      = "closed",
                        arrow_ends      = "first",
                        arrow_length    = grid::unit(0.01, "npc"),
                        line_color      = "grey10",
                        direction       = "both",
                        draw_arrowheads = TRUE
                      ),
                      gridlines = list(major = TRUE, minor = TRUE),
                      plot_border = "partial",
                      border_width = 0.8,
                      border_color = "black",
                      horizontal_line = NULL,
                      horizontal_line_aes = list(
                        type  = "longdash",
                        color = "black",
                        width = 0.4
                      ),
                      base_size = 20,
                      sample_frac = 1,
                      sample_keep_levels = NULL,
                      sample_seed = 1,
                      jitter = FALSE) {  # <--- NEW ARG

  if (!is.numeric(data[[logFC_col]])) {
    stop(paste(logFC_col, " is not numeric!", sep = ""))
  }
  if (!is.numeric(data[[pval_col]])) {
    stop(paste(pval_col, " is not numeric!", sep = ""))
  }

  data <- as.data.frame(data)

  ## --- define groups (significance vs custom highlight) -----------------------
  if (isTRUE(use_significance)) {
    data$significance <- "NS"
    data$significance[abs(data[[logFC_col]]) > logFC_cutoff] <- "FC"
    data$significance[data[[pval_col]] < pval_cutoff] <- "P"
    data$significance[
      (data[[pval_col]] < pval_cutoff) &
        (abs(data[[logFC_col]]) > logFC_cutoff)
    ] <- "FC_P"
    data$significance <- factor(
      data$significance,
      levels = c("NS", "FC", "P", "FC_P")
    )
    group_var <- data$significance
  } else {
    if (is.null(highlight_col)) {
      stop("When 'use_significance = FALSE', you must supply 'highlight_col'.")
    }
    if (!highlight_col %in% names(data)) {
      stop("Column '", highlight_col, "' not found in 'data'.")
    }
    group_var <- data[[highlight_col]]
    if (is.logical(group_var)) group_var <- ifelse(group_var, 1, 0)
    group_var <- as.factor(group_var)
  }
  data$group <- group_var

  ## --- optional down-sampling for speed --------------------------------------
  if (sample_frac < 1) {
    set.seed(sample_seed)
    if (is.null(sample_keep_levels)) {
      keep_idx <- sample(seq_len(nrow(data)),
                         size = ceiling(nrow(data) * sample_frac))
      data <- data[keep_idx, , drop = FALSE]
    } else {
      keep_flag  <- data$group %in% sample_keep_levels
      idx_keep   <- which(keep_flag)
      idx_sample <- which(!keep_flag)
      n_sample   <- ceiling(length(idx_sample) * sample_frac)
      if (n_sample > 0) {
        idx_sample <- sample(idx_sample, size = n_sample)
      }
      data <- data[c(idx_keep, idx_sample), , drop = FALSE]
    }
  }

  ## --- p == 0 handling -------------------------------------------------------
  if (min(data[[pval_col]], na.rm = TRUE) == 0) {
    warning(
      paste(
        "One or more p-values is 0.",
        "Converting to 10^-1 * lowest non-zero p-value..."
      ),
      call. = FALSE
    )
    nz <- data[[pval_col]] != 0
    min_nz <- min(data[nz, pval_col], na.rm = TRUE)
    data[!nz, pval_col] <- min_nz * 10^-1
  }

  data$labels <- labels
  data$xvals  <- data[[logFC_col]]
  data$yvals  <- data[[pval_col]]

  freq <- table(data$group)
  present_levels <- names(freq[freq > 0])

  ## --- theme ------------------------------------------------------------------
  base_theme <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(),
      plot.title        = ggplot2::element_text(
        angle = 0, size = base_size * 0.9, face = "bold", vjust = 1
      ),
      plot.subtitle     = ggplot2::element_text(
        angle = 0, size = base_size * 0.7, face = "plain", vjust = 1
      ),
      plot.caption      = ggplot2::element_text(
        angle = 0, size = base_size * 0.7, face = "plain", vjust = 1
      ),
      axis.text.x       = ggplot2::element_text(
        angle = 0, size = base_size * 0.8, vjust = 1
      ),
      axis.text.y       = ggplot2::element_text(
        angle = 0, size = base_size * 0.8, vjust = 0.5
      ),
      axis.title        = ggplot2::element_text(size = base_size * 0.9),
      legend.position   = legend_aes$position,
      legend.key        = ggplot2::element_blank(),
      legend.key.size   = grid::unit(0.5, "cm"),
      legend.text       = ggplot2::element_text(size = legend_aes$label_size),
      title             = ggplot2::element_text(size = legend_aes$label_size),
      legend.title      = ggplot2::element_blank()
    )

  ## --- colors / shapes / labels per group ------------------------------------
  if (isTRUE(use_significance)) {
    base_colors <- setNames(
      point_aes$color,
      c("NS", "FC", "P", "FC_P")
    )
    base_shapes <- setNames(
      point_aes$shape,
      c("NS", "FC", "P", "FC_P")
    )
    color_values <- base_colors[present_levels]
    shape_values <- base_shapes[present_levels]
    label_vec    <- legend_aes$labels
  } else {
    levs <- present_levels
    if (!is.null(names(highlight_colors)) &&
        all(levs %in% names(highlight_colors))) {
      color_values <- highlight_colors[levs]
    } else {
      color_values <- rep(highlight_colors, length.out = length(levs))
      names(color_values) <- levs
    }
    shape_values <- rep(point_aes$shape[1], length(levs))
    names(shape_values) <- levs

    if (!is.null(highlight_labels)) {
      label_vec <- highlight_labels
    } else if (!is.null(legend_aes$labels)) {
      label_vec <- legend_aes$labels
    } else {
      label_vec <- levs
    }
  }

  if (!is.null(label_vec)) {
    if (length(label_vec) < length(present_levels)) {
      label_vec <- c(
        label_vec,
        present_levels[(length(label_vec) + 1L):length(present_levels)]
      )
    }
    label_vec <- label_vec[seq_along(present_levels)]
  } else {
    label_vec <- present_levels
  }

  ## --- position (jitter vs no jitter) ----------------------------------------
  pos_points <- if (isTRUE(jitter)) {
    ggplot2::position_jitter(width = 0.01, height = 0.1)
  } else {
    ggplot2::position_identity()
  }

  ## --- core plot -------------------------------------------------------------
  if (is.null(point_aes$color_gradient)) {
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = xvals, y = -log10(yvals))
    ) +
      base_theme +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          order = 1,
          override.aes = list(
            shape = unname(shape_values[present_levels]),
            size  = legend_aes$icon_size
          )
        )
      ) +
      ggplot2::geom_point(
        ggplot2::aes(color = group, shape = group),
        alpha    = point_aes$alpha,
        size     = point_aes$size,
        na.rm    = TRUE,
        position = pos_points
      ) +
      ggplot2::scale_color_manual(
        values = color_values[present_levels],
        breaks = present_levels,
        labels = label_vec,
        drop   = FALSE
      ) +
      ggplot2::scale_shape_manual(
        values = shape_values[present_levels],
        breaks = present_levels,
        guide  = "none",
        drop   = FALSE
      )
  } else {
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = xvals, y = -log10(yvals))
    ) +
      base_theme +
      ggplot2::geom_point(
        ggplot2::aes(color = yvals, shape = group),
        alpha    = point_aes$alpha,
        size     = point_aes$size,
        na.rm    = TRUE,
        position = pos_points
      ) +
      ggplot2::scale_colour_gradient(
        low    = point_aes$color_gradient[1],
        high   = point_aes$color_gradient[2],
        limits = point_aes$color_gradient_limits,
        breaks = point_aes$color_gradient_breaks,
        labels = point_aes$color_gradient_labels
      ) +
      ggplot2::scale_shape_manual(
        values = shape_values[present_levels],
        breaks = present_levels,
        guide  = "none",
        drop   = FALSE
      )
  }

  ## --- axes, cutoffs, titles -------------------------------------------------
  plot <- plot +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::xlim(x_limits[1], x_limits[2]) +
    ggplot2::ylim(y_limits[1], y_limits[2]) +
    ggplot2::geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype   = cutoff_line$type,
      colour     = cutoff_line$color,
      linewidth  = cutoff_line$width
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pval_cutoff),
      linetype   = cutoff_line$type,
      colour     = cutoff_line$color,
      linewidth  = cutoff_line$width
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      caption  = caption
    )

  if (!is.null(horizontal_line)) {
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = -log10(horizontal_line),
        linetype   = horizontal_line_aes$type,
        colour     = horizontal_line_aes$color,
        linewidth  = horizontal_line_aes$width
      )
  }

  if (plot_border == "full") {
    plot <- plot +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          colour = border_color,
          fill   = NA,
          size   = border_width
        )
      )
  } else if (plot_border == "partial") {
    plot <- plot +
      ggplot2::theme(
        axis.line       = ggplot2::element_line(
          linewidth = border_width,
          colour    = border_color
        ),
        panel.border    = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank()
      )
  }

  if (isTRUE(gridlines$major)) {
    plot <- plot + ggplot2::theme(panel.grid.major = ggplot2::element_line())
  } else {
    plot <- plot + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }
  if (isTRUE(gridlines$minor)) {
    plot <- plot + ggplot2::theme(panel.grid.minor = ggplot2::element_line())
  } else {
    plot <- plot + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }

  plot <- plot + ggrepel::geom_text_repel(
    data = subset(
      data,
      data[[pval_col]] < pval_cutoff &
        abs(data[[logFC_col]]) > logFC_cutoff
    ),
    ggplot2::aes(
      label = subset(
        data,
        data[[pval_col]] < pval_cutoff &
          abs(data[[logFC_col]]) > logFC_cutoff
      )[["labels"]]
    ),
    size          = label_aes$size,
    segment.color = connector_aes$line_color,
    segment.size  = connector_aes$line_width,
    arrow = if (isTRUE(connector_aes$draw_arrowheads)) {
      grid::arrow(
        length = connector_aes$arrow_length,
        type   = connector_aes$arrow_type,
        ends   = connector_aes$arrow_ends
      )
    } else {
      NULL
    },
    colour            = label_aes$color,
    fontface          = label_aes$face,
    parse             = label_aes$parse,
    na.rm             = TRUE,
    direction         = connector_aes$direction,
    max.overlaps      = 15,
    min.segment.length = 0
  )

  if (!is.null(shade_options)) {
    plot <- plot +
      ggplot2::stat_density2d(
        data   = subset(data, rownames(data) %in% shade_options$variables),
        fill   = shade_options$fill,
        alpha  = shade_options$alpha,
        geom   = "polygon",
        contour = TRUE,
        size   = shade_options$size,
        bins   = shade_options$bins,
        show.legend = FALSE,
        na.rm  = TRUE
      )
  }

  plot <- plot + ggplot2::coord_cartesian(clip = "off")
  return(plot)
}

utils::globalVariables(c("xvals", "yvals", "significance", "group"))

