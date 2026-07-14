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
  if (!is.data.frame(qc_summary)) return(default)
  if (all(c("metric", "value") %in% names(qc_summary))) {
    hit <- qc_summary$value[as.character(qc_summary$metric) == metric]
    return(.qc_safe_value(hit, default = default))
  }
  if (metric %in% names(qc_summary) && nrow(qc_summary)) {
    return(.qc_safe_value(qc_summary[[metric]], default = default))
  }
  default
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

.qc_bullets_html <- function(items) {
  items <- as.character(items)
  items <- items[!is.na(items) & nzchar(items)]
  if (!length(items)) return("")
  paste0("<ul class=\"qc-bullets\">", paste0("<li>", .qc_html_escape(items), "</li>", collapse = ""), "</ul>")
}

.qc_callout_html <- function(title, items, class = "info") {
  body <- .qc_bullets_html(items)
  if (!nzchar(body)) return("")
  paste0(
    "<div class=\"callout callout-", .qc_html_escape(class), "\">",
    "<h3>", .qc_html_escape(title), "</h3>",
    body,
    "</div>"
  )
}

.qc_nice_percent <- function(num, den) {
  pct <- .qc_percent(num, den)
  if (identical(pct, "NA")) "not available" else pct
}

.qc_weighted_quantile <- function(x, w, probs = c(0.5, 0.9)) {
  x <- suppressWarnings(as.numeric(x))
  w <- suppressWarnings(as.numeric(w))
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) {
    idx <- which(cw >= p)[[1L]]
    x[[idx]]
  }, numeric(1L))
}

.qc_expand_weighted_values <- function(values, weights, max_total = 6000L) {
  values <- suppressWarnings(as.numeric(values))
  weights <- suppressWarnings(as.numeric(weights))
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  if (!length(values)) return(numeric())
  max_total <- suppressWarnings(as.integer(max_total[[1L]]))
  if (!is.finite(max_total) || max_total < 1L) max_total <- 6000L
  weight_total <- sum(weights)
  if (!is.finite(weight_total) || weight_total <= 0) return(numeric())
  times <- pmax(1L, round(weights / weight_total * max_total))
  times <- pmin(times, max_total)
  rep(values, times)
}

.qc_metric_sentence <- function(label, pass, tested) {
  pass_num <- suppressWarnings(as.numeric(pass))
  tested_num <- suppressWarnings(as.numeric(tested))
  if (!is.finite(pass_num) || !is.finite(tested_num) || tested_num <= 0) {
    return(paste0(label, ": not available."))
  }
  paste0(
    label, ": ", .qc_format_number(pass), " of ", .qc_format_number(tested),
    " passed (", .qc_nice_percent(pass, tested), ")."
  )
}

.qc_metric_sentence_or <- function(label, pass, tested, fallback) {
  pass_num <- suppressWarnings(as.numeric(pass))
  tested_num <- suppressWarnings(as.numeric(tested))
  if (!is.finite(pass_num) || !is.finite(tested_num) || tested_num <= 0) {
    return(fallback)
  }
  .qc_metric_sentence(label, pass, tested)
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

.qc_key_value_table <- function(x) {
  if (is.null(x) || !length(x)) return(tibble::tibble(item = character(), value = character()))
  vals <- lapply(x, function(v) {
    if (is.null(v) || !length(v)) return("NA")
    if (length(v) > 1L) return(paste(as.character(v), collapse = ", "))
    as.character(v)
  })
  tibble::tibble(item = names(vals), value = unname(unlist(vals, use.names = FALSE)))
}

.qc_status_table <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(tibble::tibble(check = character(), status = character(), detail = character()))
  need <- c("check", "status", "detail")
  for (nm in setdiff(need, names(x))) x[[nm]] <- ""
  x[, need, drop = FALSE]
}

.qc_status <- function(pass = NA, detail = "") {
  if (is.na(pass)) {
    list(status = "NOT RUN", detail = detail)
  } else if (isTRUE(pass)) {
    list(status = "PASS", detail = detail)
  } else {
    list(status = "CHECK", detail = detail)
  }
}

.qc_truncate_label <- function(x, max_chars = 34L) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  max_chars <- max(4L, as.integer(max_chars[[1L]]))
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 3L), "...")
  x
}

.qc_manifest_checks <- function(manifest) {
  if (!is.data.frame(manifest) || !nrow(manifest) || !"path" %in% names(manifest)) {
    return(tibble::tibble(table = character(), path = character(), n_rows = character(), file_exists = character(), format = character()))
  }
  n_rows <- if ("n_rows" %in% names(manifest)) .qc_format_number(manifest$n_rows) else rep("NA", nrow(manifest))
  fmt <- if ("format" %in% names(manifest)) as.character(manifest$format) else rep("NA", nrow(manifest))
  table <- if ("table" %in% names(manifest)) as.character(manifest$table) else if ("chunk_id" %in% names(manifest)) paste0("chunk_", manifest$chunk_id) else basename(as.character(manifest$path))
  tibble::tibble(
    table = table,
    path = basename(as.character(manifest$path)),
    n_rows = n_rows,
    file_exists = ifelse(file.exists(as.character(manifest$path)), "yes", "no"),
    format = fmt
  )
}

.qc_pass_rate_table <- function(rows) {
  if (!is.data.frame(rows) || !nrow(rows)) return(tibble::tibble())
  rows$pass_rate <- vapply(seq_len(nrow(rows)), function(i) {
    .qc_percent(rows$pass[[i]], rows$tested[[i]])
  }, character(1L))
  rows
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
    label <- .qc_truncate_label(x[[label_col]][[i]], 34L)
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

.qc_plot_grid <- function(...) {
  plots <- list(...)
  plots <- plots[lengths(plots) > 0L]
  if (!length(plots)) return("")
  cards <- vapply(plots, function(plot) paste0("<figure class=\"plot-card\">", plot, "</figure>"), character(1L))
  paste0("<div class=\"plot-grid\">", paste(cards, collapse = ""), "</div>")
}

.qc_plot_card_html <- function(plot, wide = FALSE) {
  plot <- as.character(plot)
  if (!length(plot) || !nzchar(plot[[1L]])) return("")
  cls <- if (isTRUE(wide)) "plot-card plot-card-wide" else "plot-card"
  paste0("<figure class=\"", cls, "\">", plot[[1L]], "</figure>")
}

.qc_finite_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.finite(x)]
}

.qc_rescale <- function(x, from = NULL, to = c(0, 1)) {
  x <- suppressWarnings(as.numeric(x))
  if (is.null(from)) {
    finite <- x[is.finite(x)]
    if (!length(finite)) return(rep(mean(to), length(x)))
    from <- range(finite, na.rm = TRUE)
  }
  if (!all(is.finite(from)) || from[[1L]] == from[[2L]]) return(rep(mean(to), length(x)))
  to[[1L]] + (x - from[[1L]]) / (from[[2L]] - from[[1L]]) * (to[[2L]] - to[[1L]])
}

.qc_color_gradient <- function(x) {
  x <- pmin(1, pmax(0, suppressWarnings(as.numeric(x))))
  low <- c(247, 247, 247)
  mid <- c(103, 169, 207)
  high <- c(5, 48, 97)
  vapply(x, function(v) {
    if (!is.finite(v)) v <- 0
    if (v <= 0.5) {
      z <- v / 0.5
      rgb <- low + (mid - low) * z
    } else {
      z <- (v - 0.5) / 0.5
      rgb <- mid + (high - mid) * z
    }
    grDevices::rgb(rgb[[1L]], rgb[[2L]], rgb[[3L]], maxColorValue = 255)
  }, character(1L))
}

.qc_funnel_svg <- function(x, label_col, value_col, title = NULL, width = 860L, height = 320L) {
  if (!is.data.frame(x) || !nrow(x) || !all(c(label_col, value_col) %in% names(x))) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[[value_col]] <- suppressWarnings(as.numeric(x[[value_col]]))
  x <- x[is.finite(x[[value_col]]) & x[[value_col]] >= 0, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  max_val <- max(x[[value_col]], na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  top <- if (is.null(title)) 16L else 44L
  bottom <- 18L
  left <- 20L
  plot_left <- 350L
  right <- 150L
  gap <- 10L
  stage_h <- max(24, (height - top - bottom - gap * (nrow(x) - 1L)) / nrow(x))
  plot_w <- width - plot_left - right
  rows <- lapply(seq_len(nrow(x)), function(i) {
    val <- x[[value_col]][[i]]
    y <- top + (i - 1L) * (stage_h + gap)
    bw <- max(38, plot_w * sqrt(val / max_val))
    bx <- plot_left + (plot_w - bw) / 2
    fill <- .qc_color_gradient((i - 1L) / max(1L, nrow(x) - 1L))
    label <- .qc_truncate_label(x[[label_col]][[i]], 42L)
    paste0(
      "<rect x=\"", bx, "\" y=\"", y, "\" width=\"", bw, "\" height=\"", stage_h, "\" rx=\"5\" fill=\"", fill, "\" class=\"funnel-bar\"/>",
      "<text x=\"", left, "\" y=\"", y + stage_h / 2 + 4L, "\" class=\"axis-label\">", .qc_html_escape(label), "</text>",
      "<text x=\"", width - right + 8L, "\" y=\"", y + stage_h / 2 + 4L, "\" class=\"value-label\">", .qc_html_escape(.qc_format_number(val)), "</text>"
    )
  })
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-funnel\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
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

.qc_density_svg <- function(values, title = NULL, width = 860L, height = 280L) {
  values <- .qc_finite_numeric(values)
  if (!length(values)) return("<p class=\"empty\">No finite plot data available.</p>")
  if (length(unique(values)) < 2L) {
    values <- c(values, values + 1e-6)
  }
  dens <- stats::density(values, n = 160L, na.rm = TRUE)
  left <- 56L
  right <- 18L
  top <- if (is.null(title)) 16L else 44L
  bottom <- 40L
  plot_w <- width - left - right
  plot_h <- height - top - bottom
  x <- .qc_rescale(dens$x, to = c(left, left + plot_w))
  y <- .qc_rescale(dens$y, to = c(top + plot_h, top))
  base_y <- top + plot_h
  points <- paste(paste0(round(x, 2), ",", round(y, 2)), collapse = " ")
  area <- paste0(left, ",", base_y, " ", points, " ", left + plot_w, ",", base_y)
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-density\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    "<line x1=\"", left, "\" y1=\"", base_y, "\" x2=\"", width - right, "\" y2=\"", base_y, "\" class=\"axis\"/>",
    "<polygon points=\"", area, "\" class=\"density-area\"/>",
    "<polyline points=\"", points, "\" class=\"density-line\"/>",
    "<text x=\"", left, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(min(values))), "</text>",
    "<text x=\"", width - right - 90L, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(max(values))), "</text>",
    "</svg>"
  )
}

.qc_cumulative_svg <- function(values, title = NULL, width = 860L, height = 280L, max_points = 300L) {
  values <- sort(.qc_finite_numeric(values))
  if (!length(values)) return("<p class=\"empty\">No finite plot data available.</p>")
  if (length(values) > max_points) {
    keep <- unique(round(seq(1L, length(values), length.out = max_points)))
    values <- values[keep]
    total_n <- max(keep)
    yvals <- keep / total_n
  } else {
    yvals <- seq_along(values) / length(values)
  }
  left <- 56L
  right <- 18L
  top <- if (is.null(title)) 16L else 44L
  bottom <- 40L
  plot_w <- width - left - right
  plot_h <- height - top - bottom
  x <- .qc_rescale(values, to = c(left, left + plot_w))
  y <- .qc_rescale(yvals, from = c(0, 1), to = c(top + plot_h, top))
  points <- paste(paste0(round(x, 2), ",", round(y, 2)), collapse = " ")
  base_y <- top + plot_h
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-cumulative\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    "<line x1=\"", left, "\" y1=\"", base_y, "\" x2=\"", width - right, "\" y2=\"", base_y, "\" class=\"axis\"/>",
    "<line x1=\"", left, "\" y1=\"", top, "\" x2=\"", left, "\" y2=\"", base_y, "\" class=\"axis\"/>",
    "<polyline points=\"", points, "\" class=\"line-strong\"/>",
    "<text x=\"", left, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(min(values))), "</text>",
    "<text x=\"", width - right - 90L, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(max(values))), "</text>",
    "<text x=\"", left + 6L, "\" y=\"", top + 14L, "\" class=\"tick\">1.0</text>",
    "</svg>"
  )
}

.qc_lollipop_svg <- function(x, label_col, value_col, title = NULL, width = 860L, height = NULL, max_rows = 20L) {
  if (!is.data.frame(x) || !nrow(x) || !all(c(label_col, value_col) %in% names(x))) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(utils::head(x, max_rows), stringsAsFactors = FALSE)
  x[[value_col]] <- suppressWarnings(as.numeric(x[[value_col]]))
  x <- x[is.finite(x[[value_col]]) & x[[value_col]] >= 0, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  if (is.null(height)) height <- max(180L, 48L + nrow(x) * 24L)
  max_val <- max(x[[value_col]], na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  left <- 350L
  right <- 130L
  top <- if (is.null(title)) 14L else 44L
  plot_w <- width - left - right
  row_h <- max(18, (height - top - 22L) / nrow(x))
  rows <- lapply(seq_len(nrow(x)), function(i) {
    y <- top + (i - 0.5) * row_h
    val <- x[[value_col]][[i]]
    px <- left + plot_w * val / max_val
    label <- .qc_truncate_label(x[[label_col]][[i]], 34L)
    paste0(
      "<text x=\"0\" y=\"", y + 4L, "\" class=\"axis-label\">", .qc_html_escape(label), "</text>",
      "<line x1=\"", left, "\" y1=\"", y, "\" x2=\"", px, "\" y2=\"", y, "\" class=\"stem\"/>",
      "<circle cx=\"", px, "\" cy=\"", y, "\" r=\"5\" class=\"point-accent\"/>",
      "<text x=\"", px + 9L, "\" y=\"", y + 4L, "\" class=\"value-label\">", .qc_html_escape(.qc_format_number(val)), "</text>"
    )
  })
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-lollipop\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    "<line x1=\"", left, "\" y1=\"", top, "\" x2=\"", left, "\" y2=\"", height - 16L, "\" class=\"axis-light\"/>",
    paste(rows, collapse = ""),
    "</svg>"
  )
}

.qc_scatter_svg <- function(x, x_col, y_col, label_col = NULL, size_col = NULL, title = NULL, width = 860L, height = 320L, max_points = 120L) {
  if (!is.data.frame(x) || !nrow(x) || !all(c(x_col, y_col) %in% names(x))) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[[x_col]] <- suppressWarnings(as.numeric(x[[x_col]]))
  x[[y_col]] <- suppressWarnings(as.numeric(x[[y_col]]))
  x <- x[is.finite(x[[x_col]]) & is.finite(x[[y_col]]), , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  if (nrow(x) > max_points) x <- x[seq_len(max_points), , drop = FALSE]
  left <- 58L
  right <- 24L
  top <- if (is.null(title)) 18L else 44L
  bottom <- 44L
  plot_w <- width - left - right
  plot_h <- height - top - bottom
  px <- .qc_rescale(x[[x_col]], to = c(left, left + plot_w))
  py <- .qc_rescale(x[[y_col]], to = c(top + plot_h, top))
  radius <- rep(4.5, nrow(x))
  if (!is.null(size_col) && size_col %in% names(x)) {
    radius <- .qc_rescale(suppressWarnings(as.numeric(x[[size_col]])), to = c(3.5, 9))
  }
  pts <- vapply(seq_len(nrow(x)), function(i) {
    label <- if (!is.null(label_col) && label_col %in% names(x)) as.character(x[[label_col]][[i]]) else ""
    paste0(
      "<circle cx=\"", round(px[[i]], 2), "\" cy=\"", round(py[[i]], 2), "\" r=\"", round(radius[[i]], 2), "\" class=\"point\"/>",
      if (nzchar(label) && i <= 8L) paste0("<text x=\"", round(px[[i]] + 7, 2), "\" y=\"", round(py[[i]] + 4, 2), "\" class=\"point-label\">", .qc_html_escape(.qc_truncate_label(label, 18L)), "</text>") else ""
    )
  }, character(1L))
  base_y <- top + plot_h
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-scatter\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    "<line x1=\"", left, "\" y1=\"", base_y, "\" x2=\"", width - right, "\" y2=\"", base_y, "\" class=\"axis\"/>",
    "<line x1=\"", left, "\" y1=\"", top, "\" x2=\"", left, "\" y2=\"", base_y, "\" class=\"axis\"/>",
    paste(pts, collapse = ""),
    "<text x=\"", left, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(x_col), "</text>",
    "<text x=\"", left + 6L, "\" y=\"", top + 14L, "\" class=\"tick\">", .qc_html_escape(y_col), "</text>",
    "</svg>"
  )
}

.qc_metric_heatmap_svg <- function(x, row_col, value_cols, title = NULL, width = 860L, cell_h = 24L, max_rows = 24L) {
  value_cols <- intersect(value_cols, names(x))
  if (!is.data.frame(x) || !nrow(x) || !(row_col %in% names(x)) || !length(value_cols)) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(utils::head(x, max_rows), stringsAsFactors = FALSE)
  mat <- as.data.frame(lapply(x[value_cols], function(v) suppressWarnings(as.numeric(v))), stringsAsFactors = FALSE)
  keep <- rowSums(is.finite(as.matrix(mat))) > 0
  x <- x[keep, , drop = FALSE]
  mat <- mat[keep, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  left <- 190L
  top <- if (is.null(title)) 34L else 58L
  right <- 20L
  col_w <- (width - left - right) / length(value_cols)
  height <- top + nrow(x) * cell_h + 24L
  scaled <- as.data.frame(lapply(mat, function(v) .qc_rescale(v, to = c(0, 1))), stringsAsFactors = FALSE)
  cells <- list()
  k <- 1L
  for (i in seq_len(nrow(x))) {
    y <- top + (i - 1L) * cell_h
    cells[[k]] <- paste0("<text x=\"0\" y=\"", y + 16L, "\" class=\"axis-label\">", .qc_html_escape(.qc_truncate_label(x[[row_col]][[i]], 28L)), "</text>")
    k <- k + 1L
    for (j in seq_along(value_cols)) {
      val <- mat[[value_cols[[j]]]][[i]]
      fill <- .qc_color_gradient(scaled[[value_cols[[j]]]][[i]])
      cells[[k]] <- paste0(
        "<rect x=\"", left + (j - 1L) * col_w, "\" y=\"", y, "\" width=\"", max(1, col_w - 2), "\" height=\"", cell_h - 2L, "\" rx=\"2\" fill=\"", fill, "\" class=\"heat-cell\"/>",
        "<text x=\"", left + (j - 1L) * col_w + 5L, "\" y=\"", y + 16L, "\" class=\"heat-label\">", .qc_html_escape(.qc_format_number(val)), "</text>"
      )
      k <- k + 1L
    }
  }
  headers <- vapply(seq_along(value_cols), function(j) {
    paste0("<text x=\"", left + (j - 1L) * col_w + 5L, "\" y=\"", top - 9L, "\" class=\"tick\">", .qc_html_escape(substr(value_cols[[j]], 1L, 20L)), "</text>")
  }, character(1L))
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-heatmap\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    paste(headers, collapse = ""),
    paste(cells, collapse = ""),
    "</svg>"
  )
}

.qc_matrix_heatmap_svg <- function(x, row_col, value_cols, title = NULL, width = 1040L, cell_h = 18L, max_rows = 28L) {
  value_cols <- intersect(value_cols, names(x))
  if (!is.data.frame(x) || !nrow(x) || !(row_col %in% names(x)) || !length(value_cols)) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(utils::head(x, max_rows), stringsAsFactors = FALSE)
  mat <- suppressWarnings(as.matrix(data.frame(lapply(x[value_cols], as.numeric), check.names = FALSE)))
  keep <- rowSums(is.finite(mat) & mat > 0, na.rm = TRUE) > 0
  x <- x[keep, , drop = FALSE]
  mat <- mat[keep, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  left <- 175L
  right <- 20L
  top <- if (is.null(title)) 86L else 112L
  col_w <- max(24L, (width - left - right) / ncol(mat))
  width <- max(width, left + right + col_w * ncol(mat))
  height <- top + nrow(mat) * cell_h + 18L
  scaled <- .qc_rescale(log1p(mat), to = c(0, 1))
  dim(scaled) <- dim(mat)
  cells <- vector("list", nrow(mat) * (ncol(mat) + 1L))
  k <- 1L
  for (i in seq_len(nrow(mat))) {
    y <- top + (i - 1L) * cell_h
    cells[[k]] <- paste0("<text x=\"0\" y=\"", y + 13L, "\" class=\"axis-label\">", .qc_html_escape(.qc_truncate_label(x[[row_col]][[i]], 26L)), "</text>")
    k <- k + 1L
    for (j in seq_len(ncol(mat))) {
      val <- mat[i, j]
      fill <- .qc_color_gradient(scaled[i, j])
      cells[[k]] <- paste0(
        "<rect x=\"", left + (j - 1L) * col_w, "\" y=\"", y, "\" width=\"", max(1, col_w - 1), "\" height=\"", cell_h - 1L, "\" rx=\"1\" fill=\"", fill, "\" class=\"heat-cell\"/>",
        if (ncol(mat) <= 12L && nrow(mat) <= 14L) {
          paste0("<text x=\"", left + (j - 1L) * col_w + 4L, "\" y=\"", y + 13L, "\" class=\"heat-label\">", .qc_html_escape(.qc_format_number(val)), "</text>")
        } else {
          ""
        }
      )
      k <- k + 1L
    }
  }
  headers <- vapply(seq_len(ncol(mat)), function(j) {
    x0 <- left + (j - 0.5) * col_w
    paste0(
      "<text transform=\"translate(", x0, ",", top - 8L, ") rotate(-45)\" class=\"tick\">",
      .qc_html_escape(.qc_truncate_label(value_cols[[j]], 18L)),
      "</text>"
    )
  }, character(1L))
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-heatmap qc-plot-matrix-heatmap\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    paste(headers, collapse = ""),
    paste(cells, collapse = ""),
    "</svg>"
  )
}

.qc_flow_svg <- function(x, label_col, value_col, title = NULL, width = 860L, height = 260L) {
  if (!is.data.frame(x) || !nrow(x) || !all(c(label_col, value_col) %in% names(x))) return("<p class=\"empty\">No plot data available.</p>")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[[value_col]] <- suppressWarnings(as.numeric(x[[value_col]]))
  x <- x[is.finite(x[[value_col]]) & x[[value_col]] >= 0, , drop = FALSE]
  if (!nrow(x)) return("<p class=\"empty\">No finite plot data available.</p>")
  top <- if (is.null(title)) 20L else 48L
  bottom <- 44L
  node_w <- 86L
  gap <- (width - 40L - node_w * nrow(x)) / max(1L, nrow(x) - 1L)
  max_val <- max(x[[value_col]], na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  center_y <- top + (height - top - bottom) / 2
  max_h <- height - top - bottom
  node_h <- pmax(16, max_h * sqrt(x[[value_col]] / max_val))
  nodes <- lapply(seq_len(nrow(x)), function(i) {
    nx <- 20L + (i - 1L) * (node_w + gap)
    ny <- center_y - node_h[[i]] / 2
    fill <- .qc_color_gradient((i - 1L) / max(1L, nrow(x) - 1L))
    label <- .qc_truncate_label(x[[label_col]][[i]], 16L)
    paste0(
      "<rect x=\"", nx, "\" y=\"", ny, "\" width=\"", node_w, "\" height=\"", node_h[[i]], "\" rx=\"6\" fill=\"", fill, "\" class=\"flow-node\"/>",
      "<text x=\"", nx, "\" y=\"", height - 24L, "\" class=\"tick\">", .qc_html_escape(label), "</text>",
      "<text x=\"", nx, "\" y=\"", height - 8L, "\" class=\"tick\">", .qc_html_escape(.qc_format_number(x[[value_col]][[i]])), "</text>"
    )
  })
  bands <- character()
  if (nrow(x) > 1L) {
    bands <- vapply(seq_len(nrow(x) - 1L), function(i) {
      x1 <- 20L + (i - 1L) * (node_w + gap) + node_w
      x2 <- 20L + i * (node_w + gap)
      bh <- max(8, min(node_h[[i]], node_h[[i + 1L]]) * 0.62)
      y1 <- center_y - bh / 2
      y2 <- center_y + bh / 2
      paste0(
        "<path d=\"M", x1, " ", y1, " C", x1 + gap * 0.45, " ", y1, " ", x2 - gap * 0.45, " ", y1, " ", x2, " ", y1,
        " L", x2, " ", y2, " C", x2 - gap * 0.45, " ", y2, " ", x1 + gap * 0.45, " ", y2, " ", x1, " ", y2, " Z\" class=\"flow-band\"/>"
      )
    }, character(1L))
  }
  title_svg <- if (is.null(title)) "" else paste0("<text x=\"0\" y=\"22\" class=\"plot-title\">", .qc_html_escape(title), "</text>")
  paste0(
    "<svg class=\"qc-plot qc-plot-flow\" viewBox=\"0 0 ", width, " ", height, "\" role=\"img\">",
    title_svg,
    paste(bands, collapse = ""),
    paste(nodes, collapse = ""),
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
    ":root{--ink:#17202a;--muted:#5d6673;--line:#d8dee8;--panel:#ffffff;--soft:#f7f9fb;--accent:#2166ac;--accent2:#1b9e77;--navy:#172a45;--paper-red:#b2182b;--paper-gold:#c99400;--warn:#9a5b00;--skip:#6b7280}",
    "*{box-sizing:border-box}body{margin:0;background:#ffffff;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:14px;line-height:1.5}",
    ".wrap{width:min(96vw,1920px);max-width:none;margin:0 auto;padding:24px 30px}",
    "header{background:#fff;border-bottom:1px solid #cfd6df}",
    "header .wrap{padding-top:24px;padding-bottom:18px}",
    "h1{font-size:28px;line-height:1.12;margin:0 0 5px 0;letter-spacing:0;color:#111827;font-weight:750}",
    "h2{font-size:17px;line-height:1.25;margin:0 0 12px 0;color:#111827;font-weight:750;border-bottom:1px solid var(--line);padding-bottom:7px}",
    "section{background:var(--panel);border:1px solid var(--line);border-radius:3px;margin:16px 0;padding:16px 18px 18px 18px;min-width:0}",
    "h3{font-size:14px;line-height:1.25;margin:0 0 8px 0;color:#111827;font-weight:750}",
    ".subtitle{color:var(--muted);font-size:13px}.metric-note{margin:7px 0 0;color:var(--muted);font-size:11px}.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:0;margin-top:4px;border:1px solid var(--line);border-bottom:0;border-right:0}",
    ".card{border-right:1px solid var(--line);border-bottom:1px solid var(--line);border-radius:0;padding:9px 11px;background:#fff}.card-label{font-size:11px;color:#4b5563;font-weight:750;text-transform:uppercase;letter-spacing:.02em}.card-value{font-size:21px;font-weight:760;margin-top:3px;color:var(--navy);line-height:1.1;font-variant-numeric:tabular-nums}",
    ".callout{border:1px solid #d7dde6;border-left:3px solid var(--accent);border-radius:2px;background:#fafbfc;padding:12px 14px;margin:10px 0}.callout-warn{border-left-color:var(--warn);background:#fffaf2}.qc-bullets{margin:0;padding-left:18px}.qc-bullets li{margin:4px 0}",
    ".report-row{padding:4px 0 2px}.report-row+.report-row{border-top:2px solid #aeb9c8;margin-top:22px;padding-top:22px}.run-setup{display:grid;grid-template-columns:minmax(0,2fr) minmax(240px,1fr);gap:18px;align-items:start}.run-setup>*{min-width:0}.condition-summary{border:1px solid var(--line);padding:14px;background:var(--soft)}.condition-count{font-size:30px;font-weight:800;color:var(--navy)}.condition-list{margin-top:8px;color:var(--muted);font-size:12px;overflow-wrap:anywhere}",
    ".plot-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:16px;align-items:start;margin-top:10px}.plot-grid-three{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:14px;align-items:start}.plot-card{margin:0;background:#fff;border:1px solid #e1e6ee;border-radius:4px;padding:14px 16px;overflow:hidden;max-width:620px;width:100%}.plot-grid .plot-card,.plot-card-wide{margin-top:12px;max-width:none;width:100%}.plot-grid-three .plot-card{max-width:none;padding:10px}.qc-tabset{min-width:0}.qc-view-tabs{display:flex;flex-wrap:wrap;gap:6px;margin:4px 0 12px;padding:4px;border-bottom:1px solid var(--line)}.qc-view-tabs button{appearance:none;border:1px solid transparent;border-radius:4px;background:transparent;color:#526071;font:inherit;font-size:12px;font-weight:750;padding:7px 11px;cursor:pointer}.qc-view-tabs button:hover{background:#f0f4f8}.qc-view-tabs button.active{background:#e8f0f8;border-color:#9fb4cc;color:#173d69}.qc-view{display:none}.qc-view.active{display:block}.qc-view-intro{margin:0 0 8px;color:var(--muted);font-size:12px}.qc-guide{margin:0 0 12px;border:0;padding:0}.qc-guide summary{display:inline-block;border:1px solid #c8d2df;border-radius:4px;background:#f7f9fb;padding:6px 9px;font-size:12px}.qc-guide .callout{margin-top:8px}",
    ".table-scroll{max-width:100%;overflow-x:auto;-webkit-overflow-scrolling:touch}table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px;border:1px solid #dfe5ed}th,td{border:1px solid #dfe5ed;padding:7px 9px;text-align:left;vertical-align:top}th{background:#f0f3f7;color:#253246;font-weight:750;white-space:nowrap}td{overflow-wrap:anywhere;word-break:break-word}tr:nth-child(even) td{background:#fbfcfd}",
    ".qc-plot{display:block;width:100%;height:auto;max-height:none}.bar{fill:var(--accent)}.axis{stroke:#44546a;stroke-width:1.1}.axis-light{stroke:#c8d0dc;stroke-width:1.1}.axis-label,.value-label,.tick{font-size:18px;fill:#253246}.value-label{font-weight:700}.plot-title{font-size:20px;font-weight:760;fill:#111827}",
    ".density-area{fill:#67a9cf;opacity:.2}.density-line,.line-strong{fill:none;stroke:var(--accent);stroke-width:2.6}.stem{stroke:var(--accent);stroke-width:2.1}.point{fill:var(--accent);stroke:white;stroke-width:1.2;opacity:.84}.point-accent{fill:var(--paper-red);stroke:white;stroke-width:1.2;opacity:.88}.point-label{font-size:14px;fill:#253246}.heat-label{font-size:12px;fill:#111827}.pca-grid,.violin-grid{stroke:#e5e9ef;stroke-width:1}.pca-point{stroke:#fff;stroke-width:1.5;opacity:.92}.pca-tick{font-size:11px;fill:#526071}.axis-title{font-size:13px;font-weight:700;fill:#111827}.condition-legend{display:flex;flex-wrap:wrap;gap:6px 12px;margin:10px 2px 2px}.condition-tag{font-size:11px;color:#435066;display:inline-flex;align-items:center;gap:5px}.condition-tag i{width:9px;height:9px;border-radius:50%;display:inline-block}.violin-shape{opacity:.82;stroke:#374151;stroke-width:.35}.violin-median{stroke:#111827;stroke-width:1.2}.violin-label{font-size:11px;fill:#253246}.qc-plot-matrix-heatmap .axis-label,.qc-plot-matrix-heatmap .tick{font-size:12px}.qc-plot-matrix-heatmap .plot-title{font-size:16px}.qc-plot-matrix-heatmap .heat-label{font-size:10px}.flow-band{fill:#67a9cf;opacity:.18}.flow-node{stroke:white;stroke-width:1}",
    ".empty{color:#69788c;font-style:italic}.status-pass{color:#1d6b46;font-weight:800}.status-warn{color:var(--warn);font-weight:800}.status-skip{color:var(--skip);font-weight:800}.links{columns:2;line-height:1.8}a{color:#244f82;text-decoration:none}a:hover{text-decoration:underline}",
    ".qc-nav{position:sticky;top:0;z-index:5;display:flex;flex-wrap:wrap;gap:6px;margin:-16px -18px 14px;padding:9px 18px;background:rgba(255,255,255,.96);border-bottom:1px solid var(--line)}.qc-nav a{padding:4px 8px;border-radius:3px;background:#f0f4f8;font-size:12px;font-weight:750}.report-link a{display:inline-block;padding:8px 11px;border:1px solid #9fb4cc;border-radius:4px;background:#eef4fa;font-weight:750}.correlation-controls{display:flex;gap:14px;flex-wrap:wrap;padding:10px 12px;background:var(--soft);border:1px solid var(--line);border-radius:3px}.correlation-controls label{font-weight:750}.correlation-controls select{margin-left:5px;padding:5px;border:1px solid #b8c4d2;border-radius:3px;background:#fff}.correlation-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin-top:14px}.correlation-card{margin:0;border:1px solid #d9dee6;border-radius:4px;padding:8px;background:#fff}.corr-title{font-size:15px;font-weight:800;fill:#111827}.corr-gridline{stroke:#ebedf0;stroke-width:1}.corr-axis{stroke:#111827;stroke-width:1.7}.corr-tick{font-size:10px;fill:#111827}.corr-axis-title{font-size:12px;font-weight:800;fill:#111827}details{margin-top:14px;border-top:1px solid var(--line);padding-top:10px}summary{cursor:pointer;font-weight:750;color:#253246;margin-bottom:8px}section{scroll-margin-top:14px}",
    "@media(min-width:1600px){.cards{grid-template-columns:repeat(8,minmax(0,1fr))}.correlation-grid{grid-template-columns:repeat(4,minmax(0,1fr))}section{padding:18px 22px 22px}}@media(max-width:900px){.plot-grid-three{grid-template-columns:1fr}.run-setup{grid-template-columns:1fr}.cards{grid-template-columns:repeat(2,minmax(0,1fr))}}@media(max-width:760px){.wrap{width:100%;padding:18px}.plot-grid,.correlation-grid{grid-template-columns:1fr}.cards{grid-template-columns:1fr}h1{font-size:25px}.links{columns:1}}",
    sep = "\n"
  )
  body <- paste(sections, collapse = "\n")
  body <- gsub("<table", "<div class=\"table-scroll\"><table", body, fixed = TRUE)
  body <- gsub("</table>", "</table></div>", body, fixed = TRUE)
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    "<title>", .qc_html_escape(title), "</title><style>", css, "</style></head><body>",
    "<header><div class=\"wrap\"><h1>", .qc_html_escape(title), "</h1><div class=\"subtitle\">CraftGRN workflow quality-control diagnostics.</div></div></header>",
    "<main class=\"wrap\">", body, "</main></body></html>"
  )
  writeLines(html, path, useBytes = TRUE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.qc_section <- function(title, body) {
  id <- tolower(gsub("[^A-Za-z0-9]+", "-", sub("^[0-9]+[.]?[[:space:]]*", "", title)))
  id <- gsub("(^-+|-+$)", "", id)
  paste0("<section id=\"", .qc_html_escape(id), "\"><h2>", .qc_html_escape(title), "</h2>", body, "</section>")
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

.qc_as_numeric_matrix <- function(x) {
  if (is.null(x)) return(matrix(numeric(), nrow = 0L, ncol = 0L))
  if (is.matrix(x)) {
    storage.mode(x) <- "double"
    return(x)
  }
  if (is.data.frame(x)) {
    rn <- rownames(x)
    out <- suppressWarnings(as.matrix(data.frame(lapply(x, as.numeric), check.names = FALSE)))
    rownames(out) <- rn
    return(out)
  }
  out <- suppressWarnings(as.matrix(x))
  storage.mode(out) <- "double"
  out
}

.qc_condition_columns <- function(omics_data) {
  if (!is_multiomic_object(omics_data)) return(character())
  mats <- omics_data$matrices
  primary <- colnames(mats$fp_score)
  primary <- primary[!is.na(primary) & nzchar(primary)]
  if (length(primary)) return(unique(primary))
  primary <- colnames(mats$gene_expr)
  primary <- primary[!is.na(primary) & nzchar(primary)]
  if (length(primary)) return(unique(primary))
  primary <- colnames(mats$fp_bound)
  primary <- primary[!is.na(primary) & nzchar(primary)]
  if (length(primary)) return(unique(primary))
  candidates <- list(mats$fp_score, mats$fp_bound, mats$gene_expr, mats$gene_on, mats$atac_score, mats$atac_open)
  cols <- unlist(lapply(candidates, colnames), use.names = FALSE)
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  cols
}

.qc_condition_metadata <- function(omics_data, conditions) {
  conditions <- as.character(conditions)
  out <- data.frame(condition = conditions, stringsAsFactors = FALSE)
  if (!is_multiomic_object(omics_data) || !is.data.frame(omics_data$samples) || !nrow(omics_data$samples)) {
    return(tibble::as_tibble(out))
  }
  samples <- as.data.frame(omics_data$samples, stringsAsFactors = FALSE)
  key_cols <- intersect(c("id", "condition_id", "condition", "sample", "name"), names(samples))
  key_col <- NA_character_
  for (nm in key_cols) {
    vals <- as.character(samples[[nm]])
    if (all(conditions %in% vals)) {
      key_col <- nm
      break
    }
  }
  if (is.na(key_col) && length(key_cols)) key_col <- key_cols[[1L]]
  if (is.na(key_col)) return(tibble::as_tibble(out))
  match_idx <- match(conditions, as.character(samples[[key_col]]))
  meta_cols <- intersect(c("id", "sample_id", "cell", "stress_type", "run_grn", "sample", "condition_id", "name"), names(samples))
  for (nm in meta_cols) {
    out[[nm]] <- samples[[nm]][match_idx]
  }
  tibble::as_tibble(out)
}

.qc_condition_matrix_indices <- function(mat, conditions, condition_meta = NULL) {
  mat <- .qc_as_numeric_matrix(mat)
  if (!nrow(mat) || !ncol(mat) || is.null(colnames(mat))) {
    return(rep(list(integer()), length(conditions)))
  }
  cn <- colnames(mat)
  direct <- match(conditions, cn)
  if (all(!is.na(direct))) return(as.list(as.integer(direct)))
  idx_list <- vector("list", length(conditions))
  for (i in seq_along(conditions)) {
    hits <- integer()
    direct_i <- match(conditions[[i]], cn)
    if (!is.na(direct_i)) hits <- c(hits, direct_i)
    if (is.data.frame(condition_meta) && nrow(condition_meta) >= i) {
      map_cols <- intersect(c("id", "sample_id", "sample", "name", "condition_id", "condition"), names(condition_meta))
      vals <- unique(as.character(unlist(condition_meta[i, map_cols, drop = FALSE], use.names = FALSE)))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      hits <- c(hits, match(vals, cn, nomatch = 0L))
    }
    idx_list[[i]] <- unique(as.integer(hits[hits > 0L]))
  }
  idx_list
}

.qc_col_count_positive <- function(mat, conditions, condition_meta = NULL) {
  mat <- .qc_as_numeric_matrix(mat)
  out <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
  if (!nrow(mat) || !ncol(mat) || is.null(colnames(mat))) return(out)
  idx_list <- .qc_condition_matrix_indices(mat, conditions, condition_meta = condition_meta)
  for (i in seq_along(conditions)) {
    idx <- idx_list[[i]]
    if (!length(idx)) next
    if (length(idx) == 1L) {
      out[[i]] <- sum(mat[, idx] > 0, na.rm = TRUE)
    } else {
      out[[i]] <- sum(rowSums(mat[, idx, drop = FALSE] > 0, na.rm = TRUE) > 0, na.rm = TRUE)
    }
  }
  out
}

.qc_col_median <- function(mat, conditions, condition_meta = NULL) {
  mat <- .qc_as_numeric_matrix(mat)
  out <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
  if (!nrow(mat) || !ncol(mat) || is.null(colnames(mat))) return(out)
  idx_list <- .qc_condition_matrix_indices(mat, conditions, condition_meta = condition_meta)
  for (i in seq_along(conditions)) {
    idx <- idx_list[[i]]
    if (!length(idx)) next
    vals <- as.numeric(mat[, idx, drop = FALSE])
    vals <- vals[is.finite(vals)]
    if (length(vals)) out[[i]] <- stats::median(vals)
  }
  out
}

.qc_col_p95 <- function(mat, conditions, condition_meta = NULL) {
  mat <- .qc_as_numeric_matrix(mat)
  out <- stats::setNames(rep(NA_real_, length(conditions)), conditions)
  if (!nrow(mat) || !ncol(mat) || is.null(colnames(mat))) return(out)
  idx_list <- .qc_condition_matrix_indices(mat, conditions, condition_meta = condition_meta)
  for (i in seq_along(conditions)) {
    idx <- idx_list[[i]]
    if (!length(idx)) next
    vals <- as.numeric(mat[, idx, drop = FALSE])
    vals <- vals[is.finite(vals)]
    if (length(vals)) out[[i]] <- as.numeric(stats::quantile(vals, 0.95, na.rm = TRUE))
  }
  out
}

.qc_multiomic_condition_qc <- function(omics_data) {
  empty <- list(
    cards = tibble::tibble(label = character(), value = character()),
    condition_summary = tibble::tibble(),
    design_summary = tibble::tibble(),
    note = "No CraftGRN multiomic object was available."
  )
  if (!is_multiomic_object(omics_data)) return(empty)
  conditions <- .qc_condition_columns(omics_data)
  if (!length(conditions)) return(empty)
  mats <- omics_data$matrices
  has_atac <- .module1_qc_has_atac(omics_data)
  condition_summary <- .qc_condition_metadata(omics_data, conditions)
  condition_summary$n_bound_fp <- as.numeric(.qc_col_count_positive(mats$fp_bound, conditions, condition_meta = condition_summary))
  if (has_atac) {
    condition_summary$n_open_atac <- as.numeric(.qc_col_count_positive(mats$atac_open, conditions, condition_meta = condition_summary))
  }
  condition_summary$n_expressed_gene <- as.numeric(.qc_col_count_positive(mats$gene_on, conditions, condition_meta = condition_summary))
  condition_summary <- condition_summary[order(condition_summary$condition), , drop = FALSE]
  run_vals <- if ("run_grn" %in% names(condition_summary)) {
    condition_summary$run_grn %in% TRUE | as.character(condition_summary$run_grn) %in% c("TRUE", "true", "1", "yes")
  } else {
    rep(NA, nrow(condition_summary))
  }
  cards <- tibble::tibble(
    label = c("Conditions", "Run-GRN conditions", "Median bound FPs", "Median expressed genes"),
    value = c(
      .qc_format_number(nrow(condition_summary)),
      .qc_format_number(sum(run_vals %in% TRUE, na.rm = TRUE)),
      .qc_format_number(stats::median(condition_summary$n_bound_fp, na.rm = TRUE)),
      .qc_format_number(stats::median(condition_summary$n_expressed_gene, na.rm = TRUE))
    )
  )
  if (has_atac) {
    cards <- dplyr::bind_rows(
      cards[1:3, , drop = FALSE],
      tibble::tibble(label = "Median open ATAC", value = .qc_format_number(stats::median(condition_summary$n_open_atac, na.rm = TRUE))),
      cards[4, , drop = FALSE]
    )
  }
  design_rows <- list()
  design_cols <- intersect(c("cell", "stress_type", "run_grn"), names(condition_summary))
  for (nm in design_cols) {
    dt <- data.table::as.data.table(condition_summary[, c(nm, "condition"), drop = FALSE])
    data.table::setnames(dt, nm, "value")
    dt[, value := as.character(value)]
    dt[is.na(value) | !nzchar(value), value := "missing"]
    summary_dt <- dt[, .(n_conditions = .N), by = "value"]
    summary_dt[["field"]] <- nm
    design_rows[[nm]] <- summary_dt[, c("field", "value", "n_conditions")]
  }
  design_summary <- data.table::rbindlist(design_rows, use.names = TRUE, fill = TRUE)
  if (nrow(design_summary)) {
    data.table::setorderv(design_summary, c("field", "n_conditions", "value"), order = c(1L, -1L, 1L))
  }
  list(
    cards = cards,
    condition_summary = tibble::as_tibble(condition_summary),
    design_summary = tibble::as_tibble(design_summary),
    note = "Condition summaries are computed from the CraftGRN multiomic matrices used by Module 1 and Module 2."
  )
}

.qc_condition_qc_html <- function(condition_qc, top_n = 30L) {
  if (!is.list(condition_qc) || !is.data.frame(condition_qc$condition_summary) || !nrow(condition_qc$condition_summary)) {
    return(paste0(
      .qc_callout_html("Condition QC Coverage", c(
        "Total bound footprints per condition replaces the bound-FP panel from 03_fp_norm_bound_summary.pdf.",
        "Total expressed genes per condition replaces 05_gene_expr_flag_summary.pdf."
      )),
      "<p class=\"empty\">No condition QC data available.</p>"
    ))
  }
  cs <- condition_qc$condition_summary
  cs_bound <- cs[order(-suppressWarnings(as.numeric(cs$n_bound_fp)), cs$condition), , drop = FALSE]
  cs_gene <- cs[order(-suppressWarnings(as.numeric(cs$n_expressed_gene)), cs$condition), , drop = FALSE]
  heat_cols <- intersect(c("n_bound_fp", "n_open_atac", "n_expressed_gene"), names(cs))
  heat_cols <- heat_cols[vapply(cs[heat_cols], function(v) any(is.finite(suppressWarnings(as.numeric(v)))), logical(1L))]
  paste0(
    .qc_callout_html("Condition QC Coverage", c(
      "This section carries the same information as the old bound-footprint and expressed-gene summary PDFs, but is computed directly from the CraftGRN multiomic object.",
      "Total bound footprints per condition replaces the bound-FP panel from 03_fp_norm_bound_summary.pdf.",
      "Total expressed genes per condition replaces 05_gene_expr_flag_summary.pdf.",
      condition_qc$note
    )),
    .qc_cards_html(condition_qc$cards),
    .qc_plot_grid(
      .qc_lollipop_svg(cs_bound, "condition", "n_bound_fp", title = "Total bound footprints per condition", max_rows = top_n),
      .qc_lollipop_svg(cs_gene, "condition", "n_expressed_gene", title = "Total expressed genes per condition", max_rows = top_n),
      .qc_metric_heatmap_svg(cs, "condition", heat_cols, title = "Condition input signal matrix", max_rows = top_n)
    ),
    .qc_table_html(condition_qc$design_summary, max_rows = 30L),
    .qc_table_html(cs, max_rows = top_n)
  )
}

.module1_qc_alignment_summary <- function(omics_data) {
  empty <- list(counts = tibble::tibble(), group_summary = tibble::tibble(), group_values = numeric())
  if (!is_multiomic_object(omics_data)) return(empty)
  fp <- omics_data$features$fp
  fp_motif <- omics_data$features$fp_motif
  atac <- omics_data$features$atac
  gene <- omics_data$features$gene
  has_atac <- .module1_qc_has_atac(omics_data)
  n_aligned <- if (is.data.frame(fp)) nrow(fp) else nrow(omics_data$matrices$fp_score)
  n_atac <- if (!has_atac) NA_real_ else if (is.data.frame(atac)) nrow(atac) else nrow(omics_data$matrices$atac_score)
  n_gene <- if (is.data.frame(gene)) nrow(gene) else nrow(omics_data$matrices$gene_expr)
  n_fm <- if (is.data.frame(fp_motif)) nrow(fp_motif) else NA_integer_
  n_fm_fp <- if (is.data.frame(fp_motif) && "fp_id" %in% names(fp_motif)) length(unique(fp_motif$fp_id)) else NA_integer_
  n_motif <- if (is.data.frame(fp_motif) && "motif" %in% names(fp_motif)) length(unique(fp_motif$motif)) else NA_integer_
  n_tf <- if (is.data.frame(fp_motif) && "tf" %in% names(fp_motif)) length(unique(fp_motif$tf)) else NA_integer_
  group_values <- numeric()
  if (is.data.frame(fp) && "group_size" %in% names(fp)) {
    group_values <- suppressWarnings(as.numeric(fp$group_size))
    group_values <- group_values[is.finite(group_values) & group_values > 0]
  }
  raw_count <- if (length(group_values)) sum(group_values) else NA_real_
  count_values <- c(raw_count, n_aligned, n_atac, n_gene, n_fm, n_fm_fp, n_motif, n_tf)
  counts <- tibble::tibble(
    metric = c(
      "raw footprints represented by aligned groups",
      "aligned footprints",
      "ATAC peaks",
      "genes",
      "FP-motif rows",
      "motif-supported aligned FPs",
      "unique motifs",
      "unique TFs"
    ),
    n = suppressWarnings(as.numeric(count_values)),
    value = .qc_format_number(count_values),
    note = c(
      if (length(group_values)) "computed from features$fp$group_size" else "not retained in this CraftGRN multiomic object",
      "compact aligned FP site count",
      if (has_atac) "compact ATAC feature count" else "ATAC master table not provided",
      "compact RNA feature count",
      "one row per FP-motif support row",
      "unique FP ids in the FP-motif table",
      "unique motif ids in the FP-motif table",
      "unique TF symbols in the FP-motif table"
    )
  )
  group_summary <- if (length(group_values)) {
    tibble::tibble(
      metric = c("median group size", "p95 group size", "max group size"),
      value = c(
        .qc_format_number(stats::median(group_values, na.rm = TRUE)),
        .qc_format_number(as.numeric(stats::quantile(group_values, 0.95, na.rm = TRUE))),
        .qc_format_number(max(group_values, na.rm = TRUE))
      ),
      detail = "Raw footprint count per aligned FP from features$fp$group_size."
    )
  } else {
    tibble::tibble(
      metric = "raw merge group size",
      value = "not available",
      detail = "The current CraftGRN multiomic object in this project does not retain raw source footprint counts per aligned FP."
    )
  }
  list(counts = counts, group_summary = group_summary, group_values = group_values)
}

.module1_qc_legacy_summary_coverage <- function(alignment_summary) {
  has_group <- length(alignment_summary$group_values) > 0L
  tibble::tibble(
    old_artifact = c("02_fp_merge_summary.pdf", "03_fp_norm_bound_summary.pdf", "05_gene_expr_flag_summary.pdf"),
    old_information = c(
      "Raw vs aligned vs ATAC counts and merge group-size distribution",
      "FP score distribution and bound footprints per condition",
      "Expressed genes per condition"
    ),
    current_html_section = c(
      "Footprint Alignment And Input QC",
      "Condition QC",
      "Condition QC"
    ),
    current_status = c(
      if (has_group) "available" else "partly available: raw group-size was not retained in this CraftGRN multiomic object",
      "available from fp_score and fp_bound matrices",
      "available from gene_on matrix"
    )
  )
}

.module1_qc_scan_predicted_tfbs <- function(manifest, top_n = 20L, verbose = TRUE, canonical_stats = NULL, canonical_fps = character()) {
  chr <- condition_support <- end <- file_exists <- fp_id <- has_canonical_support <- mean_condition_support <- n_bad_coord <- n_bad_coordinate <- n_duplicate_fp_tf_keys <- n_duplicate_key <- n_fp <- n_predicted_tfbs <- n_rows_scanned <- n_tf <- pass <- start <- tf <- N <- NULL
  if (!is.data.frame(manifest) || !nrow(manifest)) {
    return(list(
      tf_summary = tibble::tibble(),
      tf_summary_all = tibble::tibble(),
      condition_support = tibble::tibble(),
      integrity = tibble::tibble(),
      canonical_support_check = tibble::tibble()
    ))
  }
  canonical_fps <- unique(as.character(canonical_fps))
  canonical_fps <- canonical_fps[!is.na(canonical_fps) & nzchar(canonical_fps)]
  if (is.data.frame(canonical_stats) && nrow(canonical_stats) && all(c("fp_id", "tf") %in% names(canonical_stats))) {
    canon_dt <- data.table::as.data.table(canonical_stats)
    if (!"pass" %in% names(canon_dt)) canon_dt[, pass := TRUE]
    canon_dt <- unique(canon_dt[pass %in% TRUE & !is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf), .(fp_id = as.character(fp_id), tf = as.character(tf))])
    canonical_fps <- unique(c(canonical_fps, canon_dt$fp_id))
  }
  top_n <- max(1L, as.integer(top_n[[1L]]))
  tf_rows <- vector("list", nrow(manifest))
  support_rows <- vector("list", nrow(manifest))
  integrity_rows <- vector("list", nrow(manifest))
  canonical_support_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    path <- as.character(manifest$path[[i]])
    fmt <- as.character(manifest$format[[i]])
    dt <- tryCatch(
      .qc_read_table_file(path, fmt, columns = c("tf", "fp_id", "chr", "start", "end", "condition_support")),
      error = function(e) .qc_read_table_file(path, fmt, columns = c("tf", "fp_id"))
    )
    dt <- data.table::as.data.table(dt)
    if (!nrow(dt)) next
    if (!"condition_support" %in% names(dt)) dt[, condition_support := NA_real_]
    if (!"chr" %in% names(dt)) dt[, chr := NA_character_]
    if (!"start" %in% names(dt)) dt[, start := NA_real_]
    if (!"end" %in% names(dt)) dt[, end := NA_real_]
    dt[, `:=`(
      tf = as.character(tf),
      fp_id = as.character(fp_id),
      start = suppressWarnings(as.numeric(start)),
      end = suppressWarnings(as.numeric(end))
    )]
    tf_rows[[i]] <- dt[, .(
      n_predicted_tfbs = .N,
      n_fp = data.table::uniqueN(fp_id),
      mean_condition_support = mean(suppressWarnings(as.numeric(condition_support)), na.rm = TRUE)
    ), by = .(tf)]
    if (length(canonical_fps)) {
      pred_fps <- unique(dt[, .(fp_id)])
      pred_fps[, has_canonical_support := fp_id %in% canonical_fps]
      canonical_support_rows[[i]] <- pred_fps
    }
    support_rows[[i]] <- dt[, .N, by = .(condition_support)]
    bad_coord <- sum(is.na(chr) | !nzchar(chr) | !is.finite(start) | !is.finite(end) | end <= start)
    duplicate_key <- sum(duplicated(dt[, .(fp_id, tf)]))
    integrity_rows[[i]] <- data.table::data.table(
      chunk_id = if ("chunk_id" %in% names(manifest)) manifest$chunk_id[[i]] else i,
      n_rows_scanned = nrow(dt),
      n_bad_coordinate = bad_coord,
      n_duplicate_fp_tf_keys = duplicate_key,
      file_exists = file.exists(path)
    )
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
  }
  tf_summary_all <- tf_summary
  tf_summary <- head(tf_summary, top_n)
  support <- data.table::rbindlist(support_rows, use.names = TRUE, fill = TRUE)
  if (nrow(support)) {
    support <- support[, .(n_predicted_tfbs = sum(N, na.rm = TRUE)), by = condition_support]
    data.table::setorder(support, condition_support)
  }
  integrity <- data.table::rbindlist(integrity_rows, use.names = TRUE, fill = TRUE)
  if (nrow(integrity)) {
    integrity <- integrity[, .(
      n_chunks_scanned = .N,
      n_rows_scanned = sum(n_rows_scanned, na.rm = TRUE),
      n_bad_coordinate = sum(n_bad_coordinate, na.rm = TRUE),
      n_duplicate_fp_tf_keys = sum(n_duplicate_fp_tf_keys, na.rm = TRUE),
      n_missing_files = sum(!file_exists, na.rm = TRUE)
    )]
  }
  canonical_support_check <- tibble::tibble()
  canonical_support <- data.table::rbindlist(canonical_support_rows, use.names = TRUE, fill = TRUE)
  if (nrow(canonical_support)) {
    canonical_support <- canonical_support[, .(has_canonical_support = any(has_canonical_support, na.rm = TRUE)), by = fp_id]
    canonical_support_check <- tibble::tibble(
      metric = c(
        "predicted_unique_fp",
        "predicted_fp_with_motif_supported_predicted_tf",
        "predicted_fp_without_motif_supported_predicted_tf"
      ),
      value = c(
        nrow(canonical_support),
        sum(canonical_support$has_canonical_support %in% TRUE, na.rm = TRUE),
        sum(!(canonical_support$has_canonical_support %in% TRUE), na.rm = TRUE)
      )
    )
  }
  list(
    tf_summary = tibble::as_tibble(tf_summary),
    tf_summary_all = tibble::as_tibble(tf_summary_all),
    condition_support = tibble::as_tibble(support),
    integrity = tibble::as_tibble(integrity),
    canonical_support_check = canonical_support_check
  )
}

.module1_qc_input_metrics <- function(omics_data) {
  if (!is_multiomic_object(omics_data)) return(tibble::tibble(label = character(), value = character()))
  mats <- omics_data$matrices
  has_atac <- .module1_qc_has_atac(omics_data)
  n_tf <- NA_integer_
  if (is.data.frame(omics_data$features$gene) && "is_tf" %in% names(omics_data$features$gene)) {
    n_tf <- sum(omics_data$features$gene$is_tf %in% TRUE, na.rm = TRUE)
  }
  labels <- c("Conditions", "Aligned FPs", "Bound FPs", "Genes", "Expressed genes", "TF genes")
  values <- c(
    ncol(mats$fp_score), nrow(mats$fp_score),
    sum(rowSums(mats$fp_bound > 0, na.rm = TRUE) > 0),
    nrow(mats$gene_expr), sum(rowSums(mats$gene_on > 0, na.rm = TRUE) > 0), n_tf
  )
  if (has_atac) {
    labels <- append(labels, "ATAC peaks", after = 3L)
    values <- append(values, nrow(mats$atac_score), after = 3L)
  }
  tibble::tibble(label = labels, value = .qc_format_number(values))
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

.module1_qc_tf_prediction_overview <- function(tf_summary_all) {
  if (!is.data.frame(tf_summary_all) || !nrow(tf_summary_all)) {
    return(tibble::tibble(metric = character(), value = character()))
  }
  n_pred <- suppressWarnings(as.numeric(tf_summary_all$n_predicted_tfbs))
  n_fp <- suppressWarnings(as.numeric(tf_summary_all$n_fp))
  cond <- suppressWarnings(as.numeric(tf_summary_all$mean_condition_support))
  tibble::tibble(
    metric = c(
      "TFs with predicted TFBS",
      "median predicted TFBS per TF",
      "p95 predicted TFBS per TF",
      "median predicted FPs per TF",
      "median mean condition support"
    ),
    value = c(
      .qc_format_number(nrow(tf_summary_all)),
      .qc_format_number(stats::median(n_pred, na.rm = TRUE)),
      .qc_format_number(as.numeric(stats::quantile(n_pred, 0.95, na.rm = TRUE))),
      .qc_format_number(stats::median(n_fp, na.rm = TRUE)),
      .qc_format_number(stats::median(cond, na.rm = TRUE))
    )
  )
}

.module1_qc_condition_support_overview <- function(condition_support, total_conditions = NA_real_) {
  if (!is.data.frame(condition_support) || !nrow(condition_support) || !all(c("condition_support", "n_predicted_tfbs") %in% names(condition_support))) {
    return(tibble::tibble(metric = character(), value = character()))
  }
  support <- suppressWarnings(as.numeric(condition_support$condition_support))
  rows <- suppressWarnings(as.numeric(condition_support$n_predicted_tfbs))
  keep <- is.finite(support) & is.finite(rows) & rows >= 0
  support <- support[keep]
  rows <- rows[keep]
  total <- sum(rows, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(tibble::tibble(metric = character(), value = character()))
  qs <- .qc_weighted_quantile(support, rows, probs = c(0.5, 0.9))
  total_conditions <- suppressWarnings(as.numeric(total_conditions[[1L]]))
  half_cutoff <- if (is.finite(total_conditions) && total_conditions > 0) ceiling(total_conditions / 2) else NA_real_
  all_cutoff <- if (is.finite(total_conditions) && total_conditions > 0) total_conditions else NA_real_
  n_zero <- sum(rows[support <= 0], na.rm = TRUE)
  n_one <- sum(rows[support == 1], na.rm = TRUE)
  n_half <- if (is.finite(half_cutoff)) sum(rows[support >= half_cutoff], na.rm = TRUE) else NA_real_
  n_all <- if (is.finite(all_cutoff)) sum(rows[support >= all_cutoff], na.rm = TRUE) else NA_real_
  tibble::tibble(
    metric = c(
      "predicted TFBS rows scanned",
      "weighted mean condition support",
      "weighted median condition support",
      "weighted p90 condition support",
      "support in zero conditions",
      "support in one condition",
      "support in at least half of conditions",
      "support in all conditions"
    ),
    value = c(
      .qc_format_number(total),
      .qc_format_number(stats::weighted.mean(support, rows, na.rm = TRUE)),
      .qc_format_number(qs[[1L]]),
      .qc_format_number(qs[[2L]]),
      paste0(.qc_format_number(n_zero), " (", .qc_percent(n_zero, total), ")"),
      paste0(.qc_format_number(n_one), " (", .qc_percent(n_one, total), ")"),
      paste0(.qc_format_number(n_half), " (", .qc_percent(n_half, total), ")"),
      paste0(.qc_format_number(n_all), " (", .qc_percent(n_all, total), ")")
    )
  )
}

.module1_qc_parameters <- function(module1, scan_predicted_tfbs = FALSE) {
  params <- if (is.list(module1)) module1$parameters else NULL
  keep <- c(
    "db", "r_cutoff", "p_cutoff", "fdr_cutoff", "filter_to_canonical_bound",
    "min_non_na", "cores", "output_format", "write_stats", "write_bed"
  )
  out <- params[intersect(keep, names(params))]
  out$scan_predicted_tfbs <- isTRUE(scan_predicted_tfbs)
  .qc_key_value_table(out)
}

.module1_qc_flatten_values <- function(x, prefix = "") {
  if (is.null(x) || !length(x)) return(list())
  out <- list()
  for (nm in names(x)) {
    key <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    value <- x[[nm]]
    if (is.list(value) && !is.data.frame(value)) {
      out <- c(out, .module1_qc_flatten_values(value, key))
    } else if (length(value) <= 20L) {
      out[[key]] <- if (!length(value) || is.null(value)) "NA" else paste(as.character(value), collapse = ", ")
    }
  }
  out
}

.module1_qc_parameter_provenance <- function(module1, omics_data, project_config = NULL, scan_predicted_tfbs = FALSE) {
  yaml_values <- .module1_relevant_config(project_config)
  if (!length(yaml_values) && is.list(module1$parameters$project_config)) yaml_values <- module1$parameters$project_config
  if (!length(yaml_values) && is_multiomic_object(omics_data) && is.list(omics_data$project$config)) yaml_values <- omics_data$project$config
  yaml_values <- .module1_qc_flatten_values(yaml_values)
  yaml_map <- c(
    threshold_fp_tf_corr_r = "r_cutoff",
    threshold_fp_tf_corr_p = "p_cutoff",
    threshold_fp_tf_corr_fdr = "fdr_cutoff"
  )
  yaml_names <- names(yaml_values)
  mapped <- unname(yaml_map[yaml_names])
  yaml_names[!is.na(mapped)] <- mapped[!is.na(mapped)]
  names(yaml_values) <- yaml_names

  params <- if (is.list(module1$parameters)) module1$parameters else list()
  params$qc_summary <- NULL
  params$expressed_tfs <- if (length(params$expressed_tfs)) paste(sort(unique(as.character(params$expressed_tfs))), collapse = ", ") else NULL
  params$project_config <- NULL
  params$scan_predicted_tfbs <- isTRUE(scan_predicted_tfbs)
  params <- .module1_qc_flatten_values(params)
  keys <- unique(c(names(yaml_values), names(params)))
  if (!length(keys)) return(tibble::tibble())
  yaml_col <- vapply(keys, function(k) yaml_values[[k]] %||% "", character(1L))
  command_col <- vapply(keys, function(k) params[[k]] %||% "", character(1L))
  command_present <- nzchar(command_col)
  tibble::tibble(
    parameter = keys,
    yaml_value = yaml_col,
    command_value = command_col,
    effective_value = ifelse(command_present, command_col, yaml_col),
    source = ifelse(command_present & nzchar(yaml_col) & command_col != yaml_col, "command override", ifelse(command_present, "command/default", "YAML"))
  )
}

.module1_qc_report_date <- function(module1, omics_data, project_config = NULL, project_date = NULL) {
  config_values <- .module1_relevant_config(project_config)
  candidates <- list(
    project_date,
    config_values$project_date,
    if (is.list(module1$parameters)) module1$parameters$project_date else NULL,
    if (is_multiomic_object(omics_data)) omics_data$project$project_date else NULL,
    Sys.Date()
  )
  value <- candidates[[which(vapply(candidates, function(x) !is.null(x) && length(x) && !is.na(x[[1L]]) && nzchar(as.character(x[[1L]])), logical(1L)))[[1L]]]]
  parsed <- suppressWarnings(as.Date(as.character(value[[1L]])))
  if (is.na(parsed)) .log_abort("`project_date` must be coercible to YYYY-MM-DD.")
  format(parsed, "%Y-%m-%d")
}

.module1_qc_has_atac <- function(omics_data) {
  if (!is_multiomic_object(omics_data)) return(FALSE)
  presence <- omics_data$qc$input_presence$atac
  if (is.logical(presence) && length(presence) == 1L && !is.na(presence)) return(isTRUE(presence))
  is.matrix(omics_data$matrices$atac_score) && nrow(omics_data$matrices$atac_score) > 0L
}

.module1_qc_raw_footprint_recovery <- function(omics_data, module1_dir = NULL, project_config = NULL) {
  direct_paths <- character()
  if (is_multiomic_object(omics_data) && is.list(omics_data$paths)) {
    direct_paths <- unlist(omics_data$paths[intersect(c("fp_sites", "id_map"), names(omics_data$paths))], use.names = FALSE)
  }
  config <- .module1_relevant_config(project_config)
  cache_dirs <- c(config$fp_cache_dir)
  if (!is.null(config$base_dir) && nzchar(config$base_dir)) cache_dirs <- c(cache_dirs, file.path(config$base_dir, "cache"))
  if (!is.null(module1_dir) && nzchar(module1_dir)) cache_dirs <- c(cache_dirs, file.path(module1_dir, "cache"), file.path(dirname(module1_dir), "cache"))
  cache_dirs <- unique(as.character(cache_dirs))
  cache_dirs <- cache_dirs[!is.na(cache_dirs) & nzchar(cache_dirs)]
  tags <- unique(as.character(c(config$fp_cache_tag, config$db)))
  tags <- tags[!is.na(tags) & nzchar(tags)]
  cache_paths <- character()
  for (cache_dir in cache_dirs) {
    for (tag in tags) {
      cache_paths <- c(
        cache_paths,
        file.path(cache_dir, paste0("fp_sites_", tag, c(".parquet", ".csv"))),
        file.path(cache_dir, paste0("fp_id_map_", tag, c(".parquet", ".csv")))
      )
    }
  }
  paths <- unique(c(as.character(direct_paths), cache_paths))
  paths <- paths[!is.na(paths) & nzchar(paths) & file.exists(paths)]
  for (path in paths) {
    columns <- if (grepl("fp_id_map_", basename(path), fixed = TRUE)) "source_fp_peak" else c("source_fp_peaks", "n_source_fp_peaks")
    table <- tryCatch(.qc_read_table_file(path, columns = columns), error = function(e) tibble::tibble())
    aligned <- if ("source_fp_peak" %in% names(table)) list(id_map = table) else list(fp_sites = table)
    count <- .module1_raw_footprint_count(aligned)
    if (is.finite(count)) return(list(count = count, path = path))
  }
  list(count = NA_real_, path = NA_character_)
}

.module1_qc_run_cards <- function(qc_summary, omics_data, predicted_scan, predicted_rows, legacy_raw_unavailable = FALSE) {
  input_counts <- if (is_multiomic_object(omics_data) && is.list(omics_data$qc$input_counts)) omics_data$qc$input_counts else list()
  has_atac <- .module1_qc_has_atac(omics_data)
  atac_fallback <- if (has_atac && is.data.frame(omics_data$features$atac)) nrow(omics_data$features$atac) else NA_real_
  n_raw_atac <- .qc_metric_value(qc_summary, "n_raw_atac_peaks", default = input_counts$n_raw_atac_peaks %||% atac_fallback)
  n_raw_fp <- .qc_metric_value(qc_summary, "n_raw_footprints", default = input_counts$n_raw_footprints %||% NA_real_)
  n_aligned <- .qc_metric_value(qc_summary, "n_aligned_footprints", default = if (is_multiomic_object(omics_data)) nrow(omics_data$matrices$fp_score) else NA_real_)
  n_filtered <- .qc_metric_value(qc_summary, "n_canonical_bound_fps")
  expressed_fallback <- NA_real_
  if (is_multiomic_object(omics_data) && is.data.frame(omics_data$features$gene) && "is_tf" %in% names(omics_data$features$gene)) {
    expressed_fallback <- sum(
      omics_data$features$gene$is_tf %in% TRUE & rowSums(omics_data$matrices$gene_on > 0, na.rm = TRUE) > 0,
      na.rm = TRUE
    )
  }
  n_expressed <- .qc_metric_value(qc_summary, "n_expressed_tfs", default = expressed_fallback)
  scanned_tf_n <- if (is.data.frame(predicted_scan$tf_summary_all) && nrow(predicted_scan$tf_summary_all)) nrow(predicted_scan$tf_summary_all) else NA_real_
  n_pred_tfs <- .qc_metric_value(qc_summary, "n_tfs_with_predicted_binding", default = scanned_tf_n)
  format_value <- function(value, missing = "Not recorded") {
    if (is.finite(suppressWarnings(as.numeric(value)))) .qc_format_number(value) else missing
  }
  values <- c(
    if (has_atac) format_value(n_raw_atac) else "Not provided",
    format_value(n_raw_fp, if (isTRUE(legacy_raw_unavailable)) "Unavailable for legacy run" else "Not recorded"),
    format_value(n_aligned), format_value(n_filtered), format_value(predicted_rows),
    format_value(n_expressed), format_value(n_pred_tfs)
  )
  tibble::tibble(
    label = c(
      "Raw ATAC peaks", "Raw footprints", "Aligned footprints",
      "Filtered footprints", "Predicted unique TFBS",
      "Expressed TFs", "TFs with predicted binding"
    ),
    value = values
  )
}

.module1_qc_gate_table <- function(qc_summary, omics_data = NULL) {
  n_aligned <- if (is_multiomic_object(omics_data)) nrow(omics_data$matrices$fp_score) else .qc_metric_value(qc_summary, "n_fp_input")
  n_bound <- .qc_metric_value(qc_summary, "n_fp_bound_accessible")
  n_motif <- .qc_metric_value(qc_summary, "n_motif_supported_pairs")
  n_canon_pair <- .qc_metric_value(qc_summary, "n_canonical_pairs_pass")
  n_canon_fp <- .qc_metric_value(qc_summary, "n_canonical_bound_fps")
  n_pred_pair <- .qc_metric_value(qc_summary, "n_prediction_pairs")
  n_pred_tfbs <- .qc_metric_value(qc_summary, "n_predicted_tfbs")
  rows <- tibble::tibble(
    gate = c(
      "Aligned FP to bound/accessibility-supported FP",
      "Motif-supported FP-TF pairs to canonical-bound pairs",
      "Bound/accessibility-supported FP to canonical-bound FP",
      "Prediction FP-TF pairs to predicted TFBS"
    ),
    tested = c(n_aligned, n_motif, n_bound, n_pred_pair),
    pass = c(n_bound, n_canon_pair, n_canon_fp, n_pred_tfbs),
    removed = c(n_aligned - n_bound, n_motif - n_canon_pair, n_bound - n_canon_fp, n_pred_pair - n_pred_tfbs)
  )
  finite_rows <- is.finite(suppressWarnings(as.numeric(rows$tested))) &
    is.finite(suppressWarnings(as.numeric(rows$pass)))
  .qc_pass_rate_table(rows[finite_rows, , drop = FALSE])
}

.module1_qc_read_canonical_stats <- function(module1, module1_dir = NULL) {
  if (is.list(module1) && is.data.frame(module1$motif_supported_correlations) && nrow(module1$motif_supported_correlations)) {
    return(tibble::as_tibble(module1$motif_supported_correlations))
  }
  if (!is.null(module1_dir)) {
    path <- file.path(module1_dir, "module1_canonical_prediction_stats.csv.gz")
    if (file.exists(path)) return(tibble::as_tibble(readr::read_csv(path, show_col_types = FALSE)))
  }
  tibble::tibble()
}

.module1_qc_read_canonical_fp_cache <- function(module1, module1_dir = NULL) {
  if (is.list(module1) && is.data.frame(module1$high_confidence_footprints) && nrow(module1$high_confidence_footprints)) {
    keep <- intersect(c("fp_id", "n_canonical_bound_tfs", "has_canonical_bound"), names(module1$high_confidence_footprints))
    return(tibble::as_tibble(module1$high_confidence_footprints[, keep, drop = FALSE]))
  }
  if (!is.null(module1_dir)) {
    paths <- c(
      file.path(module1_dir, "cache", "module1_canonical_bound_fps.csv.gz"),
      file.path(module1_dir, "module1_high_confidence_footprints.csv")
    )
    for (path in paths) {
      if (file.exists(path)) {
        out <- tibble::as_tibble(readr::read_csv(path, show_col_types = FALSE))
        keep <- intersect(c("fp_id", "n_canonical_bound_tfs", "has_canonical_bound"), names(out))
        if (length(keep)) return(out[, keep, drop = FALSE])
      }
    }
  }
  tibble::tibble()
}

.module1_qc_normalize_prediction_stats <- function(x) {
  x <- tibble::as_tibble(x)
  if (nrow(x) && !"pass" %in% names(x)) x$pass <- TRUE
  x
}

.module1_qc_read_prediction_stats <- function(module1, module1_dir = NULL) {
  if (is.list(module1) && is.data.frame(module1$prediction_stats) && nrow(module1$prediction_stats)) {
    return(.module1_qc_normalize_prediction_stats(module1$prediction_stats))
  }
  if (!is.null(module1_dir)) {
    path <- file.path(module1_dir, "module1_prediction_stats.csv.gz")
    if (file.exists(path)) return(.module1_qc_normalize_prediction_stats(readr::read_csv(path, show_col_types = FALSE)))
    manifest_candidates <- c(
      file.path(module1_dir, "module1_prediction_stats_manifest.csv"),
      file.path(module1_dir, "module1_prediction_stats_chunks", "module1_prediction_stats_manifest.csv")
    )
    existing_manifest <- manifest_candidates[file.exists(manifest_candidates)]
    manifest_path <- if (length(existing_manifest)) existing_manifest[[1L]] else manifest_candidates[[1L]]
    manifest <- .qc_read_manifest_chunks(manifest_path)
    if (nrow(manifest)) {
      max_rows <- getOption("craftgrn.qc_prediction_stats_max_rows", 1000000L)
      max_rows <- suppressWarnings(as.integer(max_rows[[1L]]))
      if (!is.finite(max_rows) || max_rows < 1L) max_rows <- 1000000L
      keep_cols <- c(
        "fp_id", "tf", "pearson_r", "pearson_p", "pearson_fdr",
        "spearman_r", "spearman_p", "spearman_fdr", "best_r",
        "best_method", "best_p", "best_fdr", "pass"
      )
      rows <- vector("list", nrow(manifest))
      n_loaded <- 0L
      for (i in seq_len(nrow(manifest))) {
        if (n_loaded >= max_rows) break
        chunk <- tryCatch(
          .qc_read_table_file(
            as.character(manifest$path[[i]]),
            as.character(manifest$format[[i]]),
            columns = keep_cols
          ),
          error = function(e) .qc_read_table_file(
            as.character(manifest$path[[i]]),
            as.character(manifest$format[[i]]),
            columns = NULL
          )
        )
        if (nrow(chunk)) chunk <- chunk[, intersect(keep_cols, names(chunk)), drop = FALSE]
        if (!nrow(chunk)) next
        remaining <- max_rows - n_loaded
        if (nrow(chunk) > remaining) chunk <- utils::head(chunk, remaining)
        rows[[i]] <- chunk
        n_loaded <- n_loaded + nrow(chunk)
      }
      out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
      if (nrow(out)) return(.module1_qc_normalize_prediction_stats(out))
    }
  }
  tibble::tibble()
}

.module1_qc_corr_summary <- function(x, label = "correlation") {
  best_r <- fp_id <- n_pass <- pass <- tf <- NULL
  if (!is.data.frame(x) || !nrow(x)) {
    return(list(summary = tibble::tibble(), top_tf = tibble::tibble(), best_r = numeric()))
  }
  dt <- data.table::as.data.table(x)
  if (!"pass" %in% names(dt)) dt[, pass := FALSE]
  if (!"best_r" %in% names(dt)) dt[, best_r := NA_real_]
  dt[, best_r := suppressWarnings(as.numeric(best_r))]
  summary <- tibble::tibble(
    table = label,
    n_rows = nrow(dt),
    n_pass = sum(dt$pass %in% TRUE, na.rm = TRUE),
    pass_rate = .qc_percent(sum(dt$pass %in% TRUE, na.rm = TRUE), nrow(dt)),
    median_best_r = stats::median(dt$best_r, na.rm = TRUE),
    p95_best_r = as.numeric(stats::quantile(dt$best_r, 0.95, na.rm = TRUE)),
    max_best_r = suppressWarnings(max(dt$best_r, na.rm = TRUE))
  )
  top_tf <- tibble::tibble()
  if ("tf" %in% names(dt)) {
    top_tf_dt <- dt[pass %in% TRUE, .(
      n_pass = .N,
      n_fp = if ("fp_id" %in% names(dt)) data.table::uniqueN(fp_id) else NA_integer_,
      median_best_r = stats::median(best_r, na.rm = TRUE)
    ), by = tf]
    if (nrow(top_tf_dt)) {
      data.table::setorder(top_tf_dt, -n_pass, tf)
      top_tf <- tibble::as_tibble(head(top_tf_dt, 20L))
    }
  }
  list(summary = summary, top_tf = top_tf, best_r = dt$best_r[is.finite(dt$best_r)])
}

.module1_qc_warning_checks <- function(pred_manifest, support_check, predicted_scan, canonical_stats, prediction_stats) {
  manifest_check <- if (is.data.frame(pred_manifest) && nrow(pred_manifest)) {
    all(file.exists(as.character(pred_manifest$path)))
  } else {
    NA
  }
  canonical_missing <- .qc_metric_value(support_check, "predicted_fp_without_motif_supported_predicted_tf", default = NA_real_)
  scan_integrity <- predicted_scan$integrity
  bad_coord <- .qc_metric_value(scan_integrity, "n_bad_coordinate", default = NA_real_)
  bad_duplicate <- .qc_metric_value(scan_integrity, "n_duplicate_fp_tf_keys", default = NA_real_)
  status_rows <- list(
    list("Predicted TFBS manifest files exist", .qc_status(manifest_check, "Every listed chunk path should exist.")),
    list("Predicted FPs retain canonical support", .qc_status(if (is.finite(canonical_missing)) canonical_missing == 0 else NA, paste0("Missing canonical support count: ", .qc_format_number(canonical_missing)))),
    list("Predicted TFBS coordinates are valid", .qc_status(if (is.finite(bad_coord)) bad_coord == 0 else NA, paste0("Bad coordinate rows: ", .qc_format_number(bad_coord)))),
    list("Predicted FP-TF keys are not duplicated within chunks", .qc_status(if (is.finite(bad_duplicate)) bad_duplicate == 0 else NA, paste0("Duplicate keys detected within chunks: ", .qc_format_number(bad_duplicate)))),
    list("Motif-supported correlation table is available", .qc_status(is.data.frame(canonical_stats) && nrow(canonical_stats) > 0, paste0("Rows: ", .qc_format_number(nrow(canonical_stats))))),
    list("Prediction correlation statistics are available", .qc_status(is.data.frame(prediction_stats) && nrow(prediction_stats) > 0, paste0("Rows: ", .qc_format_number(nrow(prediction_stats)))))
  )
  tibble::tibble(
    check = vapply(status_rows, `[[`, character(1L), 1L),
    status = vapply(status_rows, function(x) x[[2L]]$status, character(1L)),
    detail = vapply(status_rows, function(x) x[[2L]]$detail, character(1L))
  )
}

.module1_qc_canonical_status <- function(support_check, canonical_stats, qc_summary) {
  canonical_missing <- .qc_metric_value(support_check, "predicted_fp_without_motif_supported_predicted_tf", default = NA_real_)
  canonical_unique <- .qc_metric_value(support_check, "predicted_unique_fp", default = NA_real_)
  n_canonical_pass <- .qc_metric_value(qc_summary, "n_canonical_pairs_pass", default = NA_real_)
  n_canonical_pairs <- .qc_metric_value(qc_summary, "n_motif_supported_pairs", default = NA_real_)
  if (is.finite(canonical_missing) && canonical_missing == 0) {
    return(list(
      card = "PASS",
      html = "<span class=\"status-pass\">PASS</span>",
      detail = paste0("All ", .qc_format_number(canonical_unique), " predicted FPs retained canonical support."),
      checked = TRUE
    ))
  }
  if (is.finite(canonical_missing)) {
    return(list(
      card = "CHECK",
      html = "<span class=\"status-warn\">CHECK</span>",
      detail = paste0(.qc_format_number(canonical_missing), " predicted FPs did not have canonical support in the scanned support table."),
      checked = TRUE
    ))
  }
  if (is.data.frame(canonical_stats) && nrow(canonical_stats)) {
    return(list(
      card = "Not scanned",
      html = "<span class=\"status-warn\">NOT SCANNED</span>",
      detail = "Canonical correlation statistics were available, but the predicted TFBS chunks were not scanned against them.",
      checked = FALSE
    ))
  }
  if (is.finite(n_canonical_pass) && is.finite(n_canonical_pairs)) {
    return(list(
      card = "Recorded",
      html = "<span class=\"status-pass\">RECORDED</span>",
      detail = paste0("Run summary recorded ", .qc_format_number(n_canonical_pass), " passing canonical pairs out of ", .qc_format_number(n_canonical_pairs), "."),
      checked = FALSE
    ))
  }
  list(
    card = "Not checked",
    html = "<span class=\"status-warn\">NOT CHECKED</span>",
    detail = "Canonical-support files were not found in this output directory, so the report cannot independently verify canonical support for this saved run.",
    checked = FALSE
  )
}

.module1_qc_content_sections <- function(qc_summary,
                                         predicted_rows,
                                         pred_manifest,
                                         support_check,
                                         predicted_scan,
                                         warning_checks,
                                         condition_overview = tibble::tibble(),
                                         omics_data = NULL,
                                         canonical_status = NULL) {
  n_input <- .qc_metric_value(
    qc_summary,
    "n_fp_input",
    default = if (is_multiomic_object(omics_data)) nrow(omics_data$matrices$fp_score) else NA_real_
  )
  n_bound <- .qc_metric_value(
    qc_summary,
    "n_fp_bound_accessible",
    default = if (is_multiomic_object(omics_data)) sum(rowSums(omics_data$matrices$fp_bound > 0, na.rm = TRUE) > 0) else NA_real_
  )
  n_canonical_fp <- .qc_metric_value(
    qc_summary,
    "n_canonical_bound_fps",
    default = .qc_metric_value(support_check, "predicted_unique_fp", default = NA_real_)
  )
  n_prediction_pairs <- .qc_metric_value(qc_summary, "n_prediction_pairs")
  n_canonical_pairs <- .qc_metric_value(qc_summary, "n_motif_supported_pairs")
  n_canonical_pass <- .qc_metric_value(qc_summary, "n_canonical_pairs_pass")
  canonical_missing <- .qc_metric_value(support_check, "predicted_fp_without_motif_supported_predicted_tf", default = NA_real_)
  canonical_unique <- .qc_metric_value(support_check, "predicted_unique_fp", default = n_canonical_fp)
  top_tf <- if (is.data.frame(predicted_scan$tf_summary) && nrow(predicted_scan$tf_summary)) {
    paste0(
      as.character(predicted_scan$tf_summary$tf[[1L]]), " has ",
      .qc_format_number(predicted_scan$tf_summary$n_predicted_tfbs[[1L]]),
      " predicted TFBS rows."
    )
  } else {
    "Top predicted TF summary was not scanned."
  }
  n_check <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "CHECK", na.rm = TRUE) else NA_integer_
  n_not_run <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "NOT RUN", na.rm = TRUE) else NA_integer_
  n_skipped <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "SKIPPED", na.rm = TRUE) else NA_integer_
  support_zero <- if (is.data.frame(condition_overview) && nrow(condition_overview)) {
    hit <- condition_overview$value[condition_overview$metric == "support in zero conditions"]
    if (length(hit)) as.character(hit[[1L]]) else "not available"
  } else {
    "not available"
  }
  canonical_pair_sentence <- .qc_metric_sentence_or(
    "Canonical support",
    n_canonical_pass,
    n_canonical_pairs,
    "Canonical support counts were not recorded in this saved output directory."
  )
  canonical_fp_sentence <- if (is.finite(suppressWarnings(as.numeric(n_canonical_fp)))) {
    paste0("After canonical support, ", .qc_format_number(n_canonical_fp), " FPs were eligible for all expressed TF prediction.")
  } else {
    "Canonical-bound FP eligibility was not recorded in this saved output directory."
  }
  prediction_sentence <- if (is.finite(suppressWarnings(as.numeric(n_prediction_pairs)))) {
    paste0("The prediction stage evaluated ", .qc_format_number(n_prediction_pairs), " FP-TF pairs and wrote ", .qc_format_number(predicted_rows), " predicted TFBS rows across ", .qc_format_number(nrow(pred_manifest)), " chunks.")
  } else {
    paste0("The saved output contains ", .qc_format_number(predicted_rows), " predicted TFBS rows across ", .qc_format_number(nrow(pred_manifest)), " chunks; the total evaluated FP-TF pair count was not recorded.")
  }
  canonical_check_sentence <- if (is.finite(canonical_missing) && is.finite(canonical_unique)) {
    paste0("The canonical-supported FP check found ", .qc_format_number(canonical_missing), " predicted FPs without canonical support out of ", .qc_format_number(canonical_unique), " predicted FPs.")
  } else if (is.list(canonical_status)) {
    canonical_status$detail
  } else {
    "Canonical-supported FP check details were not recorded in this saved output directory."
  }
  key <- .qc_callout_html("What Happened", c(
    paste0("Module 1 started from ", .qc_format_number(n_input), " aligned FPs and kept ", .qc_format_number(n_bound), " bound/accessibility-supported FPs."),
    canonical_pair_sentence,
    if (is.list(canonical_status)) paste0("Canonical support check: ", canonical_status$card, ". ", canonical_status$detail) else "",
    canonical_fp_sentence,
    prediction_sentence,
    top_tf
  ))
  interp <- .qc_callout_html("Module 1 Interpretation", c(
    "Input Gates shows what the filter did at each stage: bound/accessibility support, motif-supported canonical evidence, canonical-bound FP selection, and final predicted TFBS selection.",
    "Canonical support means a predicted FP has at least one motif-supported TF-FP pair that passed the configured correlation cutoffs before broader TF prediction.",
    canonical_check_sentence,
    paste0("Condition support summarizes how broadly a predicted TFBS is active across samples; rows with zero support are ", support_zero, "."),
    "Prediction correlation diagnostics should be read as the broad TFBS prediction step after canonical-bound FPs have already been selected."
  ))
  review_class <- if (is.finite(n_check) && n_check > 0) "warn" else "info"
  review <- .qc_callout_html("Review Before Downstream Use", c(
    paste0("Warning Checks contains ", .qc_format_number(n_check), " CHECK rows, ", .qc_format_number(n_not_run), " NOT RUN rows, and ", .qc_format_number(n_skipped), " SKIPPED optional rows."),
    "Review TFs with very high predicted counts to make sure they reflect biology rather than overly permissive cutoffs.",
    "Review Condition Support to confirm that predicted TFBS are not dominated by a very small number of conditions or by support-zero rows.",
    "Review Prediction Output Integrity before using predicted_tfbs as the Module 2 handoff."
  ), class = review_class)
  list(key = key, interpretation = interp, review = review)
}

.module1_qc_sample_annotations <- function(omics_data, sample_ids) {
  sample_ids <- unique(as.character(sample_ids))
  out <- tibble::tibble(sample_id = sample_ids, display_label = sample_ids, biological_group = sample_ids, group_field = "sample_id", matched = FALSE)
  if (!is_multiomic_object(omics_data) || !is.data.frame(omics_data$samples) || !nrow(omics_data$samples)) return(out)
  samples <- as.data.frame(omics_data$samples, stringsAsFactors = FALSE)
  key_priority <- c("id", "sample_id", "condition_id", "module1_condition_id", "name", "sample", "id_atac", "id_fp", "id_rna", "strict_match_rna", "broad_match_rna")
  label_priority <- c("condition_display", "name", "condition_id", "sample", "sample_id", "id")
  group_priority <- c("condition_base", "comparison_group", "cell_stress_type", "condition", "cell_context", "cell", "condition_id", "name")
  columns_for <- function(keys) {
    unlist(lapply(keys, function(key) names(samples)[tolower(names(samples)) == key]), use.names = FALSE)
  }
  key_columns <- unique(columns_for(key_priority))
  label_columns <- unique(columns_for(label_priority))
  group_columns <- unique(columns_for(group_priority))
  for (i in seq_along(sample_ids)) {
    row_index <- NA_integer_
    for (column in key_columns) {
      hit <- which(as.character(samples[[column]]) == sample_ids[[i]])
      if (length(hit)) {
        row_index <- hit[[1L]]
        break
      }
    }
    if (!is.finite(row_index)) next
    first_value <- function(columns, fallback) {
      for (column in columns) {
        value <- as.character(samples[[column]][[row_index]])
        if (!is.na(value) && nzchar(value)) return(list(value = value, column = column))
      }
      list(value = fallback, column = "sample_id")
    }
    label <- first_value(label_columns, sample_ids[[i]])
    group <- first_value(group_columns, label$value)
    out$display_label[[i]] <- label$value
    out$biological_group[[i]] <- group$value
    out$group_field[[i]] <- group$column
    out$matched[[i]] <- TRUE
  }
  out
}

.module1_qc_transform_matrix <- function(mat, assay = c("rna", "atac", "fp")) {
  assay <- match.arg(assay)
  mat <- .qc_as_numeric_matrix(mat)
  if (!length(mat)) return(list(matrix = mat, transformation = "identity"))
  row_index <- unique(round(seq(1, nrow(mat), length.out = min(nrow(mat), 10000L))))
  probe <- as.numeric(mat[row_index, , drop = FALSE])
  probe <- probe[is.finite(probe)]
  has_negative <- length(probe) && any(probe < 0)
  p99 <- if (length(probe)) as.numeric(stats::quantile(probe, 0.99, na.rm = TRUE, names = FALSE)) else NA_real_
  use_log <- !identical(assay, "fp") && !has_negative && is.finite(p99) && p99 > 50
  if (use_log) mat <- log2(pmax(mat, 0) + 1)
  list(matrix = mat, transformation = if (use_log) "log2(x + 1)" else "identity")
}

.module1_qc_pca_scores <- function(mat, label, annotations = NULL, assay = c("rna", "atac", "fp"), max_features = 2000L) {
  assay <- match.arg(assay)
  prepared <- .module1_qc_transform_matrix(mat, assay = assay)
  mat <- prepared$matrix
  if (nrow(mat) < 2L || ncol(mat) < 3L) return(tibble::tibble())
  variance <- apply(mat, 1L, stats::var, na.rm = TRUE)
  keep_n <- min(as.integer(max_features), sum(is.finite(variance) & variance > 0))
  keep <- head(order(variance, decreasing = TRUE, na.last = NA), keep_n)
  if (length(keep) < 2L) return(tibble::tibble())
  selected <- mat[keep, , drop = FALSE]
  for (i in seq_len(nrow(selected))) {
    bad <- !is.finite(selected[i, ])
    if (any(bad)) {
      replacement <- stats::median(selected[i, !bad], na.rm = TRUE)
      if (!is.finite(replacement)) replacement <- 0
      selected[i, bad] <- replacement
    }
  }
  fit <- tryCatch(stats::prcomp(t(selected), center = TRUE, scale. = TRUE), error = function(e) NULL)
  if (is.null(fit) || ncol(fit$x) < 2L) return(tibble::tibble())
  explained <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  sample_id <- rownames(fit$x)
  if (!is.data.frame(annotations) || !nrow(annotations)) {
    annotations <- tibble::tibble(sample_id = sample_id, display_label = sample_id, biological_group = sample_id, group_field = "sample_id", matched = FALSE)
  }
  annotations <- annotations[match(sample_id, annotations$sample_id), , drop = FALSE]
  tibble::tibble(
    sample_id = sample_id,
    display_label = ifelse(is.na(annotations$display_label), sample_id, annotations$display_label),
    biological_group = ifelse(is.na(annotations$biological_group), sample_id, annotations$biological_group),
    metadata_matched = annotations$matched %in% TRUE,
    PC1 = fit$x[, 1L],
    PC2 = fit$x[, 2L],
    panel = label,
    transformation = prepared$transformation,
    n_features = length(keep),
    x_label = paste0("PC1 (", round(explained[[1L]], 1), "%)"),
    y_label = paste0("PC2 (", round(explained[[2L]], 1), "%)")
  )
}

.qc_ggplot_svg <- function(plot, width = 8, height = 5, css_class = "qc-plot") {
  path <- tempfile(fileext = ".svg")
  device_open <- FALSE
  on.exit({
    if (device_open) try(grDevices::dev.off(), silent = TRUE)
    unlink(path)
  }, add = TRUE)
  rendered <- tryCatch({
    grDevices::svg(path, width = width, height = height, bg = "white", onefile = TRUE, pointsize = 11)
    device_open <- TRUE
    print(plot)
    grDevices::dev.off()
    device_open <- FALSE
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }, error = function(e) "")
  if (!nzchar(rendered)) return("<p class=\"empty\">Plot rendering failed.</p>")
  start <- regexpr("<svg", rendered, fixed = TRUE)[[1L]]
  if (start < 1L) return("<p class=\"empty\">Plot rendering failed.</p>")
  rendered <- substr(rendered, start, nchar(rendered))
  id_matches <- regmatches(rendered, gregexpr("id=\"[^\"]+\"", rendered, perl = TRUE))[[1L]]
  if (length(id_matches) && !identical(id_matches, character())) {
    ids <- unique(gsub("^id=\"|\"$", "", id_matches))
    ids <- ids[order(nchar(ids), decreasing = TRUE)]
    prefix <- paste0("qcsvg-", gsub("[^A-Za-z0-9]", "", basename(path)), "-")
    for (id in ids) {
      rendered <- gsub(paste0("id=\"", id, "\""), paste0("id=\"", prefix, id, "\""), rendered, fixed = TRUE)
      rendered <- gsub(paste0("#", id), paste0("#", prefix, id), rendered, fixed = TRUE)
    }
  }
  sub("<svg ", paste0("<svg class=\"", css_class, "\" "), rendered, fixed = TRUE)
}

.module1_qc_pca_svg <- function(x, title, color_map) {
  PC1 <- PC2 <- biological_group <- display_label <- NULL
  if (!is.data.frame(x) || !nrow(x)) return(paste0("<h3>", .qc_html_escape(title), "</h3><p class=\"empty\">PCA requires at least three samples and two variable features.</p>"))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(x, ggplot2::aes(x = PC1, y = PC2, color = biological_group, label = display_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "#d7dde6", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, color = "#d7dde6", linewidth = 0.35) +
    ggplot2::geom_point(size = 3.6, alpha = 0.92) +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::labs(
      subtitle = paste0(x$transformation[[1L]], "; ", .qc_format_number(x$n_features[[1L]]), " variable features"),
      x = x$x_label[[1L]], y = x$y_label[[1L]], color = "Biological group"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#526071"),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 34, 8, 8)
  )
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(color = "#334155", size = 2.65, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.25, min.segment.length = 0, seed = 1, show.legend = FALSE)
  } else {
    p <- p + ggplot2::geom_text(color = "#334155", size = 2.4, check_overlap = TRUE, vjust = -0.8, show.legend = FALSE)
  }
  paste0(
    "<div class=\"plot-heading\"><h3>", .qc_html_escape(title), "</h3></div>",
    .qc_ggplot_svg(p, width = 10.4, height = 6.2, css_class = "qc-plot pca-plot")
  )
}

.module1_qc_condition_legend <- function(color_map) {
  if (!length(color_map)) return("")
  items <- vapply(seq_along(color_map), function(i) {
    paste0("<span class=\"condition-tag\"><i style=\"background:", unname(color_map[[i]]), "\"></i>", .qc_html_escape(names(color_map)[[i]]), "</span>")
  }, character(1L))
  paste0("<div class=\"condition-legend\">", paste(items, collapse = ""), "</div>")
}

.module1_qc_violin_svg <- function(x, title, stage_labels = NULL, width = 900L) {
  value <- y <- group <- stage <- q25 <- q75 <- median <- NULL
  need <- c("condition", "stage", "probability", "value")
  if (!is.data.frame(x) || !nrow(x) || !all(need %in% names(x))) return("<p class=\"empty\">Distribution data are unavailable.</p>")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x <- x[is.finite(x$probability) & is.finite(x$value), , drop = FALSE]
  conditions <- unique(as.character(x$condition))
  stages <- unique(as.character(x$stage))
  if (!length(conditions) || !length(stages)) return("<p class=\"empty\">Distribution data are unavailable.</p>")
  labels <- c(raw = "Raw", quantile_normalized = "Quantile-normalized", canonical_filtered = "Filtered", gene_expression = "Gene expression")
  if (!is.null(stage_labels)) labels[names(stage_labels)] <- stage_labels
  colors <- c(raw = "#27c1c7", quantile_normalized = "#ed8582", canonical_filtered = "#315f97", gene_expression = "#31a354")
  density_parts <- list()
  summary_parts <- list()
  k <- 1L
  stage_offset <- stats::setNames((seq_along(stages) - (length(stages) + 1) / 2) * 0.22, stages)
  violin_half_height <- if (length(stages) == 1L) 0.27 else 0.085
  for (ci in seq_along(conditions)) {
    condition <- conditions[[ci]]
    base_y <- length(conditions) - ci + 1L
    for (si in seq_along(stages)) {
      stage <- stages[[si]]
      d <- x[x$condition == condition & x$stage == stage, , drop = FALSE]
      d <- d[order(d$probability), , drop = FALSE]
      if (nrow(d) < 5L) next
      probs <- seq(0.005, 0.995, length.out = 201L)
      sample_values <- stats::approx(d$probability, d$value, xout = probs, rule = 2, ties = "ordered")$y
      center <- base_y + stage_offset[[stage]]
      den <- tryCatch(stats::density(sample_values, n = 128L, na.rm = TRUE), error = function(e) NULL)
      if (is.null(den) || !length(den$y) || max(den$y) <= 0) next
      scaled_density <- violin_half_height * den$y / max(den$y)
      density_parts[[k]] <- tibble::tibble(
        condition = condition,
        stage = stage,
        group = paste(condition, stage, sep = "::"),
        value = c(den$x, rev(den$x)),
        y = c(center + scaled_density, rev(center - scaled_density))
      )
      quartiles <- stats::approx(d$probability, d$value, xout = c(0.25, 0.5, 0.75), rule = 2, ties = "ordered")$y
      summary_parts[[k]] <- tibble::tibble(condition = condition, stage = stage, y = center, q25 = quartiles[[1L]], median = quartiles[[2L]], q75 = quartiles[[3L]])
      k <- k + 1L
    }
  }
  density_data <- dplyr::bind_rows(density_parts)
  summary_data <- dplyr::bind_rows(summary_parts)
  if (!nrow(density_data)) return("<p class=\"empty\">Distribution data are unavailable.</p>")
  stage_labels <- stats::setNames(unname(labels[stages]), stages)
  stage_colors <- stats::setNames(vapply(stages, function(stage) colors[[stage]] %||% "#64748b", character(1L)), stages)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = density_data, ggplot2::aes(x = value, y = y, group = group, fill = stage), color = "#ffffff", linewidth = 0.25, alpha = 0.72) +
    ggplot2::geom_segment(data = summary_data, ggplot2::aes(x = q25, xend = q75, y = y, yend = y), color = "#172033", linewidth = 1.2) +
    ggplot2::geom_point(data = summary_data, ggplot2::aes(x = median, y = y, fill = stage), shape = 21, color = "#172033", stroke = 0.4, size = 2.1) +
    ggplot2::scale_y_continuous(breaks = rev(seq_along(conditions)), labels = conditions, expand = ggplot2::expansion(add = c(0.35, 0.5))) +
    ggplot2::scale_fill_manual(values = stage_colors, labels = stage_labels, drop = FALSE) +
    ggplot2::labs(
      x = if (identical(stages, "gene_expression")) "Expression value" else "Footprint score",
      y = NULL, fill = NULL, color = NULL,
      caption = "Violin shapes are deterministic reconstructions from saved quantile summaries; center dots and lines show the median and interquartile range."
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8.5, color = "#111827"),
      panel.grid.major.y = ggplot2::element_line(color = "#eef1f5", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.key.height = grid::unit(3.5, "mm"),
      legend.text = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 8, color = "#64748b", hjust = 0),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
  paste0(
    "<div class=\"plot-heading\"><h3>", .qc_html_escape(title), "</h3></div>",
    .qc_ggplot_svg(p, width = width / 100, height = max(4.4, 2.1 + length(conditions) * 0.25), css_class = "qc-plot violin-plot raincloud-plot distribution-plot")
  )
}

.module1_qc_distribution_table <- function(omics_data, canonical_fps = character()) {
  if (!is_multiomic_object(omics_data)) return(tibble::tibble())
  out <- if (is.data.frame(omics_data$qc$fp_score_distributions)) omics_data$qc$fp_score_distributions else tibble::tibble()
  probs <- seq(0, 1, length.out = 101L)
  if (!nrow(out)) {
    out <- dplyr::bind_rows(lapply(seq_len(ncol(omics_data$matrices$fp_score)), function(i) {
      values <- omics_data$matrices$fp_score[, i]
      tibble::tibble(
        condition = colnames(omics_data$matrices$fp_score)[[i]],
        stage = "quantile_normalized",
        probability = probs,
        value = as.numeric(stats::quantile(values, probs, na.rm = TRUE, names = FALSE)),
        n = length(values)
      )
    }))
  }
  keep <- match(canonical_fps, rownames(omics_data$matrices$fp_score), nomatch = 0L)
  keep <- unique(keep[keep > 0L])
  if (length(keep)) {
    filtered <- dplyr::bind_rows(lapply(seq_len(ncol(omics_data$matrices$fp_score)), function(i) {
      values <- omics_data$matrices$fp_score[keep, i]
      tibble::tibble(
        condition = colnames(omics_data$matrices$fp_score)[[i]],
        stage = "canonical_filtered",
        probability = probs,
        value = as.numeric(stats::quantile(values, probs, na.rm = TRUE, names = FALSE)),
        n = length(values)
      )
    }))
    out <- dplyr::bind_rows(out, filtered)
  }
  tibble::as_tibble(out)
}

.module1_qc_distribution_summary <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) return(tibble::tibble())
  data.table::as.data.table(x)[, .(
    median = value[which.min(abs(probability - 0.5))],
    p95 = value[which.min(abs(probability - 0.95))],
    n = max(n, na.rm = TRUE)
  ), by = .(condition, stage)] |>
    tibble::as_tibble()
}

.module1_qc_correlation_browser <- function(canonical_stats, prediction_stats, bins = 40L, hist_rows = NULL) {
  parts <- list()
  add_scope <- function(x, scope) {
    if (!is.data.frame(x) || !nrow(x) || !"tf" %in% names(x)) return()
    for (method in c("pearson_r", "spearman_r")) {
      if (!method %in% names(x)) next
      pass <- if ("pass" %in% names(x)) x$pass %in% TRUE else suppressWarnings(as.numeric(x[[method]])) >= 0.5
      dt <- data.table::data.table(tf = as.character(x$tf), value = suppressWarnings(as.numeric(x[[method]])), pass = pass)
      dt <- dt[is.finite(value) & value >= -1 & value <= 1]
      if (!nrow(dt)) next
      dt[, bin := pmin(bins - 1L, pmax(0L, floor((value + 1) / 2 * bins)))]
      per_tf <- dt[, .(n = .N, n_pass = sum(pass)), by = .(tf, bin)]
      overall <- dt[, .(n = .N, n_pass = sum(pass)), by = bin][, tf := "Overall"]
      out <- data.table::rbindlist(list(overall, per_tf), use.names = TRUE)
      out[, `:=`(scope = scope, method = method)]
      parts[[length(parts) + 1L]] <<- out[, .(scope, method, tf, bin, n, n_pass)]
    }
  }
  if (is.data.frame(hist_rows) && nrow(hist_rows)) {
    rows <- data.table::as.data.table(hist_rows)
  } else {
    add_scope(canonical_stats, "Canonical motif-supported")
    add_scope(prediction_stats, "Full prediction")
    if (!length(parts)) return("<p class=\"empty\">Correlation histograms are unavailable in this saved run.</p>")
    rows <- data.table::rbindlist(parts, use.names = TRUE)
  }
  if (!"n_pass" %in% names(rows)) rows[, n_pass := 0]
  rows[, n_pass := pmin(as.numeric(n), pmax(0, as.numeric(n_pass)))]
  json <- jsonlite::toJSON(rows, dataframe = "rows", auto_unbox = TRUE)
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  paste0(
    "<div class=\"correlation-browser\"><div class=\"correlation-controls\"><label>Stage <select id=\"corr-scope\"></select></label><label>TF 1 <select id=\"corr-tf-1\"></select></label><label>TF 2 <select id=\"corr-tf-2\"></select></label><label>TF 3 <select id=\"corr-tf-3\"></select></label></div><div id=\"corr-grid\" class=\"correlation-grid\"></div></div>",
    "<script>(()=>{const R=", json, ",scope=document.getElementById('corr-scope'),sels=[1,2,3].map(i=>document.getElementById('corr-tf-'+i)),grid=document.getElementById('corr-grid');const scopes=[...new Set(R.map(x=>x.scope))],alpha=[...new Set(R.filter(x=>x.tf!=='Overall').map(x=>x.tf))].sort(),rank=[...alpha].sort((a,b)=>R.filter(x=>x.tf===b).reduce((s,x)=>s+(+x.n_pass||0),0)-R.filter(x=>x.tf===a).reduce((s,x)=>s+(+x.n_pass||0),0)||a.localeCompare(b));scopes.forEach(x=>scope.add(new Option(x,x)));sels.forEach((s,i)=>{alpha.forEach(x=>s.add(new Option(x,x)));s.value=rank[Math.min(i,rank.length-1)]||alpha[i]||''});function hist(method,name){const rows=R.filter(x=>x.scope===scope.value&&x.method===method&&x.tf===name),w=360,h=285,l=58,b=52,t=42,pw=w-l-12,ph=h-b-t,m=Math.max(1,...rows.map(x=>+x.n)),color=method==='pearson_r'?'#67a9cf':'#ef8a62',title=name+' '+(method==='pearson_r'?'Pearson R':'Spearman R');let s='<figure class=\"correlation-card\"><svg class=\"qc-plot\" viewBox=\"0 0 '+w+' '+h+'\"><text x=\"'+(w/2)+'\" y=\"22\" text-anchor=\"middle\" class=\"corr-title\">'+title+'</text>';for(let i=0;i<=2;i++){const v=m*i/2,y=t+ph-ph*i/2;s+='<line x1=\"'+l+'\" y1=\"'+y+'\" x2=\"'+(w-12)+'\" y2=\"'+y+'\" class=\"corr-gridline\"/><text x=\"'+(l-7)+'\" y=\"'+(y+4)+'\" text-anchor=\"end\" class=\"corr-tick\">'+Math.round(v).toLocaleString()+'</text>'}rows.forEach(x=>{const bw=pw/40,px=l+(+x.bin)*bw,bh=ph*(+x.n)/m,passh=ph*(+x.n_pass||0)/m;s+='<rect x=\"'+px+'\" y=\"'+(t+ph-bh)+'\" width=\"'+Math.max(1,bw-.5)+'\" height=\"'+bh+'\" fill=\"#dddddd\"><title>Total: '+(+x.n).toLocaleString()+'</title></rect><rect x=\"'+px+'\" y=\"'+(t+ph-passh)+'\" width=\"'+Math.max(1,bw-.5)+'\" height=\"'+passh+'\" fill=\"'+color+'\"><title>Passing: '+(+x.n_pass||0).toLocaleString()+'</title></rect>'});[-1,-.5,0,.5,1].forEach(v=>{const x=l+(v+1)/2*pw;s+='<text x=\"'+x+'\" y=\"'+(h-29)+'\" text-anchor=\"middle\" class=\"corr-tick\">'+v+'</text>'});s+='<line x1=\"'+l+'\" y1=\"'+(t+ph)+'\" x2=\"'+(w-12)+'\" y2=\"'+(t+ph)+'\" class=\"corr-axis\"/><line x1=\"'+l+'\" y1=\"'+t+'\" x2=\"'+l+'\" y2=\"'+(t+ph)+'\" class=\"corr-axis\"/><text x=\"'+(l+pw/2)+'\" y=\"'+(h-5)+'\" text-anchor=\"middle\" class=\"corr-axis-title\">R value</text><text x=\"15\" y=\"'+(t+ph/2)+'\" transform=\"rotate(-90 15 '+(t+ph/2)+')\" text-anchor=\"middle\" class=\"corr-axis-title\">Count</text></svg></figure>';return s}function distinct(){const used=new Set();sels.forEach(s=>{if(used.has(s.value)){const next=alpha.find(x=>!used.has(x));if(next)s.value=next}used.add(s.value)})}function draw(){distinct();const names=['Overall',...sels.map(s=>s.value)];grid.innerHTML=names.map(n=>hist('pearson_r',n)+hist('spearman_r',n)).join('')}scope.addEventListener('change',draw);sels.forEach(s=>s.addEventListener('change',draw));draw()})()</script>"
  )
}

#' Build a Module 1 QC HTML report
#'
#' Builds a comprehensive HTML report for Module 1 run parameters, input gates,
#' motif-supported canonical support, prediction output integrity, correlation
#' diagnostics, condition-level CraftGRN multiomic input QC, metadata-aware PCA,
#' footprint alignment summaries, warning checks, and interactive Binding,
#' Co-binding, and Motif review. The report can consume a
#' `predict_tfbs()` result, a step-by-step Module 1 result list, or a Module 1
#' output directory.
#'
#' @param module1 Module 1 result list or Module 1 output directory.
#' @param omics_data Optional CraftGRN multiomic object. Used when `module1` is
#'   an output directory or does not contain `omics_data`.
#' @param output_dir Directory where the HTML report should be written. If
#'   `NULL`, the report is written under `reports` inside the Module 1 output
#'   directory when available.
#' @param report_name HTML report filename.
#' @param scan_predicted_tfbs Logical; if `TRUE`, scan predicted TFBS chunks to
#'   summarize top TFs and condition support. This is comprehensive but can take
#'   extra time on full projects.
#' @param top_n Number of TFs to show in top-TF summaries.
#' @param project_config Optional YAML path or list used for parameter
#'   provenance in the Run Summary.
#' @param project_date Optional project date displayed in the report title.
#' @param build_tfbs_explorer Logical; embed the interactive Binding,
#'   Co-binding, and Motif explorer sections in the QC report.
#' @param default_condition Optional initial condition for TFBS report previews
#'   and the linked explorer.
#' @param verbose Emit concise progress messages.
#' @return Normalized path to the HTML report.
#' @export
build_module1_qc_report <- function(module1,
                                    omics_data = NULL,
                                    output_dir = NULL,
                                    report_name = "module1_qc_report.html",
                                    scan_predicted_tfbs = TRUE,
                                    top_n = 20L,
                                    project_config = NULL,
                                    project_date = NULL,
                                    build_tfbs_explorer = TRUE,
                                    default_condition = NULL,
                                    verbose = TRUE) {
  module1_input <- module1
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

  canonical_stats <- if (isTRUE(scan_predicted_tfbs)) {
    .module1_qc_read_canonical_stats(module1, module1_dir = module1_dir)
  } else {
    tibble::tibble()
  }
  canonical_fp_cache <- .module1_qc_read_canonical_fp_cache(module1, module1_dir = module1_dir)
  canonical_fps <- if (is.data.frame(canonical_fp_cache) && nrow(canonical_fp_cache) && "fp_id" %in% names(canonical_fp_cache)) {
    as.character(canonical_fp_cache$fp_id)
  } else {
    character()
  }
  prediction_stats <- if (isTRUE(scan_predicted_tfbs)) {
    .module1_qc_read_prediction_stats(module1, module1_dir = module1_dir)
  } else {
    tibble::tibble()
  }
  predicted_scan <- if (isTRUE(scan_predicted_tfbs)) {
    .module1_qc_scan_predicted_tfbs(
      pred_manifest,
      top_n = top_n,
      verbose = verbose,
      canonical_stats = canonical_stats,
      canonical_fps = canonical_fps
    )
  } else {
    list(tf_summary = tibble::tibble(), tf_summary_all = tibble::tibble(), condition_support = tibble::tibble(), integrity = tibble::tibble(), canonical_support_check = tibble::tibble())
  }
  if (!nrow(support_check) && is.data.frame(predicted_scan$canonical_support_check)) {
    support_check <- predicted_scan$canonical_support_check
  }
  canonical_corr <- .module1_qc_corr_summary(canonical_stats, label = "Motif-supported FP-TF correlations")
  prediction_corr <- .module1_qc_corr_summary(prediction_stats, label = "Prediction FP-TF correlations")
  motif_complexity <- .module1_qc_motif_complexity(omics_data)
  alignment_summary <- .module1_qc_alignment_summary(omics_data)
  legacy_coverage <- .module1_qc_legacy_summary_coverage(alignment_summary)
  condition_qc <- .qc_multiomic_condition_qc(omics_data)
  input_cards <- .module1_qc_input_metrics(omics_data)
  tf_overview <- .module1_qc_tf_prediction_overview(predicted_scan$tf_summary_all)
  total_conditions <- if (is_multiomic_object(omics_data)) ncol(omics_data$matrices$fp_score) else NA_real_
  condition_overview <- .module1_qc_condition_support_overview(predicted_scan$condition_support, total_conditions = total_conditions)
  predicted_rows <- if (nrow(pred_manifest)) sum(suppressWarnings(as.numeric(pred_manifest$n_rows)), na.rm = TRUE) else .qc_metric_value(qc_summary, "n_predicted_tfbs")
  canonical_status <- .module1_qc_canonical_status(support_check, canonical_stats, qc_summary)

  input_counts <- if (is_multiomic_object(omics_data) && is.list(omics_data$qc$input_counts)) omics_data$qc$input_counts else list()
  raw_fp_known <- .qc_metric_value(qc_summary, "n_raw_footprints", default = input_counts$n_raw_footprints %||% NA_real_)
  raw_recovery <- list(count = NA_real_, path = NA_character_)
  if (!is.finite(raw_fp_known)) {
    raw_recovery <- .module1_qc_raw_footprint_recovery(omics_data, module1_dir = module1_dir, project_config = project_config)
    if (is.finite(raw_recovery$count)) {
      qc_summary <- dplyr::bind_rows(qc_summary, tibble::tibble(metric = "n_raw_footprints", value = raw_recovery$count))
      raw_fp_known <- raw_recovery$count
    }
  }
  legacy_raw_unavailable <- !is.finite(raw_fp_known) && !is.null(module1_dir)
  cards <- .module1_qc_run_cards(qc_summary, omics_data, predicted_scan, predicted_rows, legacy_raw_unavailable = legacy_raw_unavailable)
  parameter_table <- .module1_qc_parameter_provenance(module1, omics_data, project_config, scan_predicted_tfbs)
  report_date <- .module1_qc_report_date(module1, omics_data, project_config, project_date)
  report_title <- paste0("Module 1 QC Report (", report_date, ")")
  explorer_source <- if (!is.null(module1_dir)) module1_dir else module1_input
  explorer_cache <- NULL
  if (isTRUE(build_tfbs_explorer) && is_multiomic_object(omics_data)) {
    cache_path <- .module1_explorer_cache_path(explorer_source, file.path(output_dir, report_name))
    explorer_cache <- tryCatch(
      .module1_build_explorer_cache(explorer_source, omics_data, cache_path, rebuild = FALSE, verbose = verbose),
      error = function(e) {
        .log_warn("Embedded Module 1 TFBS Explorer was not built: {conditionMessage(e)}")
        NULL
      }
    )
  } else if (is_multiomic_object(omics_data)) {
    cache_path <- .module1_explorer_cache_path(explorer_source, file.path(output_dir, report_name))
    if (file.exists(cache_path)) explorer_cache <- tryCatch(readRDS(cache_path), error = function(e) NULL)
  }
  if (is.list(explorer_cache)) {
    if (!length(canonical_fps) && length(explorer_cache$used_site_bits) && is_multiomic_object(omics_data)) {
      used_idx <- .module1_unpack_site_indices(explorer_cache$used_site_bits, nrow(omics_data$matrices$fp_score))
      canonical_fps <- rownames(omics_data$matrices$fp_score)[used_idx]
    }
    predicted_rows <- sum(explorer_cache$tf_counts[, "Overall"], na.rm = TRUE)
    cards$value[cards$label == "Filtered footprints"] <- .qc_format_number(length(canonical_fps))
    cards$value[cards$label == "Predicted unique TFBS"] <- .qc_format_number(predicted_rows)
    cards$value[cards$label == "TFs with predicted binding"] <- .qc_format_number(length(explorer_cache$tf_names))
    if (!is.finite(.qc_metric_value(qc_summary, "n_canonical_bound_fps"))) {
      qc_summary <- dplyr::bind_rows(qc_summary, tibble::tibble(metric = "n_canonical_bound_fps", value = length(canonical_fps)))
    }
    if (!is.finite(.qc_metric_value(qc_summary, "n_tfs_with_predicted_binding"))) {
      qc_summary <- dplyr::bind_rows(qc_summary, tibble::tibble(metric = "n_tfs_with_predicted_binding", value = length(explorer_cache$tf_names)))
    }
  }
  default_condition_use <- if (is_multiomic_object(omics_data)) {
    .module1_explorer_default_condition(.qc_condition_columns(omics_data), default_condition, project_config, omics_data)
  } else {
    ""
  }
  explorer_components <- if (is.list(explorer_cache) && isTRUE(build_tfbs_explorer)) {
    .module1_tfbs_explorer_components(explorer_cache, default_condition_use, standalone = FALSE)
  } else {
    NULL
  }
  gate_table <- .module1_qc_gate_table(qc_summary, omics_data = omics_data)
  manifest_checks <- .qc_manifest_checks(pred_manifest)
  warning_checks <- .module1_qc_warning_checks(
    pred_manifest = pred_manifest,
    support_check = support_check,
    predicted_scan = predicted_scan,
    canonical_stats = canonical_stats,
    prediction_stats = prediction_stats
  )
  related <- character()
  if (!is.null(module1_dir)) {
    related <- c(
      list.files(module1_dir, pattern = "\\.html$", full.names = TRUE),
      pred_manifest_path
    )
  }
  condition_support_plots <- .qc_plot_grid(
    .qc_lollipop_svg(
      predicted_scan$condition_support,
      "condition_support",
      "n_predicted_tfbs",
      title = "Predicted TFBS rows by condition support"
    ),
    .qc_cumulative_svg(
      .qc_expand_weighted_values(
        predicted_scan$condition_support$condition_support,
        predicted_scan$condition_support$n_predicted_tfbs
      ),
      title = "Condition-support cumulative profile"
    )
  )
  alignment_counts_for_plot <- alignment_summary$counts
  if (is.data.frame(alignment_counts_for_plot) && nrow(alignment_counts_for_plot)) {
    alignment_counts_for_plot <- alignment_counts_for_plot[is.finite(suppressWarnings(as.numeric(alignment_counts_for_plot$n))), , drop = FALSE]
  }
  alignment_plots <- .qc_plot_grid(
    .qc_lollipop_svg(alignment_counts_for_plot, "metric", "n", title = "Footprint merge: peak counts"),
    .qc_hist_svg(alignment_summary$group_values, title = "Footprint merge: group size distribution"),
    .qc_cumulative_svg(alignment_summary$group_values, title = "Footprint merge: group size cumulative curve")
  )
  footprint_distribution_plots <- .qc_plot_grid(
    .qc_density_svg(motif_complexity$values, title = "Canonical TFs per aligned FP density"),
    .qc_cumulative_svg(motif_complexity$values, title = "Cumulative TF multiplicity per aligned FP")
  )
  content <- .module1_qc_content_sections(
    qc_summary = qc_summary,
    predicted_rows = predicted_rows,
    pred_manifest = pred_manifest,
    support_check = support_check,
    predicted_scan = predicted_scan,
    warning_checks = warning_checks,
    condition_overview = condition_overview,
    omics_data = omics_data,
    canonical_status = canonical_status
  )
  distribution_table <- .module1_qc_distribution_table(omics_data, canonical_fps = canonical_fps)
  sample_ids <- if (is_multiomic_object(omics_data)) {
    unique(unlist(lapply(
      list(omics_data$matrices$gene_expr, omics_data$matrices$atac_score, omics_data$matrices$fp_score),
      colnames
    ), use.names = FALSE))
  } else {
    character()
  }
  sample_annotations <- .module1_qc_sample_annotations(omics_data, sample_ids)
  display_map <- stats::setNames(sample_annotations$display_label, sample_annotations$sample_id)
  if (nrow(distribution_table)) {
    readable <- unname(display_map[as.character(distribution_table$condition)])
    distribution_table$condition <- ifelse(is.na(readable) | !nzchar(readable), as.character(distribution_table$condition), readable)
  }
  correlation_browser <- .module1_qc_correlation_browser(
    canonical_stats,
    prediction_stats,
    hist_rows = if (is.list(explorer_cache)) explorer_cache$correlation_hist else NULL
  )
  distribution_summary <- .module1_qc_distribution_summary(distribution_table)
  if (nrow(distribution_summary)) distribution_summary$condition_stage <- paste(distribution_summary$condition, distribution_summary$stage, sep = " | ")
  gene_distribution <- if (is_multiomic_object(omics_data)) {
    probs <- seq(0, 1, length.out = 101L)
    gene_prepared <- .module1_qc_transform_matrix(omics_data$matrices$gene_expr, assay = "rna")
    dplyr::bind_rows(lapply(seq_len(ncol(gene_prepared$matrix)), function(i) {
      values <- suppressWarnings(as.numeric(gene_prepared$matrix[, i]))
      values <- values[is.finite(values)]
      sample_id <- colnames(gene_prepared$matrix)[[i]]
      readable <- unname(display_map[sample_id])
      if (!length(readable) || is.na(readable) || !nzchar(readable)) readable <- sample_id
      tibble::tibble(
        condition = readable,
        stage = "gene_expression",
        probability = probs,
        value = as.numeric(stats::quantile(values, probs, names = FALSE)),
        n = length(values),
        transformation = gene_prepared$transformation
      )
    }))
  } else {
    tibble::tibble()
  }
  pca_rna <- if (is_multiomic_object(omics_data)) .module1_qc_pca_scores(omics_data$matrices$gene_expr, "RNA PCA", annotations = sample_annotations, assay = "rna") else tibble::tibble()
  pca_atac <- if (is_multiomic_object(omics_data)) .module1_qc_pca_scores(omics_data$matrices$atac_score, "ATAC PCA", annotations = sample_annotations, assay = "atac") else tibble::tibble()
  filtered_idx <- if (is_multiomic_object(omics_data)) match(canonical_fps, rownames(omics_data$matrices$fp_score), nomatch = 0L) else integer()
  filtered_idx <- unique(filtered_idx[filtered_idx > 0L])
  pca_fp <- if (is_multiomic_object(omics_data) && length(filtered_idx)) {
    .module1_qc_pca_scores(omics_data$matrices$fp_score[filtered_idx, , drop = FALSE], "Filtered FP PCA", annotations = sample_annotations, assay = "fp")
  } else {
    tibble::tibble()
  }
  has_atac <- .module1_qc_has_atac(omics_data)
  pca_groups <- sort(unique(unlist(lapply(list(pca_rna, pca_atac, pca_fp), function(x) {
    if (is.data.frame(x) && "biological_group" %in% names(x)) as.character(x$biological_group) else character()
  }), use.names = FALSE)))
  pca_groups <- pca_groups[!is.na(pca_groups) & nzchar(pca_groups)]
  pca_colors <- stats::setNames(grDevices::hcl.colors(max(3L, length(pca_groups)), palette = "Dark 3")[seq_along(pca_groups)], pca_groups)
  pca_plots <- paste0(
    "<div class=\"qc-tabset\"><div class=\"qc-view-tabs\" role=\"tablist\" aria-label=\"PCA assay\">",
    "<button type=\"button\" class=\"active\" data-qc-group=\"pca\" data-qc-view=\"rna\" aria-selected=\"true\">RNA</button>",
    if (has_atac) "<button type=\"button\" data-qc-group=\"pca\" data-qc-view=\"atac\" aria-selected=\"false\">ATAC</button>" else "",
    "<button type=\"button\" data-qc-group=\"pca\" data-qc-view=\"fp\" aria-selected=\"false\">Filtered footprints</button></div>",
    "<div class=\"qc-view active\" data-qc-panel-group=\"pca\" data-qc-panel=\"rna\"><figure class=\"plot-card plot-card-wide\">", .module1_qc_pca_svg(pca_rna, "RNA PCA", pca_colors), "</figure></div>",
    if (has_atac) paste0("<div class=\"qc-view\" data-qc-panel-group=\"pca\" data-qc-panel=\"atac\"><figure class=\"plot-card plot-card-wide\">", .module1_qc_pca_svg(pca_atac, "ATAC PCA", pca_colors), "</figure></div>") else "",
    "<div class=\"qc-view\" data-qc-panel-group=\"pca\" data-qc-panel=\"fp\"><figure class=\"plot-card plot-card-wide\">", .module1_qc_pca_svg(pca_fp, "Filtered FP PCA", pca_colors), "</figure></div></div>",
    .module1_qc_condition_legend(pca_colors)
  )
  filtered_by_condition <- if (is_multiomic_object(omics_data) && length(filtered_idx)) {
    condition_ids <- colnames(omics_data$matrices$fp_bound)
    readable <- unname(display_map[condition_ids])
    tibble::tibble(
      condition = ifelse(is.na(readable) | !nzchar(readable), condition_ids, readable),
      n_filtered_fp = as.numeric(colSums(omics_data$matrices$fp_bound[filtered_idx, , drop = FALSE] > 0, na.rm = TRUE))
    )
  } else {
    tibble::tibble()
  }
  predicted_by_condition <- if (is.list(explorer_cache) && !is.null(explorer_cache$tf_counts)) {
    readable <- unname(display_map[explorer_cache$conditions])
    tibble::tibble(
      condition = ifelse(is.na(readable) | !nzchar(readable), explorer_cache$conditions, readable),
      n_predicted_tfbs = as.numeric(colSums(explorer_cache$tf_counts[, explorer_cache$conditions, drop = FALSE]))
    )
  } else {
    tibble::tibble()
  }
  explorer_empty <- "<p class=\"empty\">Interactive TFBS data are unavailable in this saved run.</p>"
  explorer_style <- if (is.list(explorer_components)) paste0("<style>", explorer_components$css, "</style>") else ""
  n_conditions <- if (is_multiomic_object(omics_data)) ncol(omics_data$matrices$fp_score) else .qc_metric_value(qc_summary, "n_conditions")
  condition_names <- if (is_multiomic_object(omics_data)) colnames(omics_data$matrices$fp_score) else character()
  summary_cards <- dplyr::bind_rows(
    tibble::tibble(label = "# Conditions", value = .qc_format_number(n_conditions)),
    cards
  )
  setup_html <- paste0(
    "<details><summary>Input parameters and condition roster</summary><div class=\"run-setup\"><div><h3>Input parameters</h3><p class=\"subtitle\">Effective YAML and command-line values with their provenance.</p>",
    .qc_table_html(parameter_table, max_rows = 100L),
    "</div><aside class=\"condition-summary\"><div class=\"card-label\">Condition roster</div><div class=\"condition-count\">", .qc_format_number(n_conditions), "</div><div class=\"condition-list\">", .qc_html_escape(paste(condition_names, collapse = " | ")), "</div></aside></div></details>"
  )
  missing_fp_stages <- setdiff(c("raw", "quantile_normalized", "canonical_filtered"), unique(as.character(distribution_table$stage)))
  distribution_notice <- if (length(missing_fp_stages)) {
    paste0("<p class=\"empty\">", .qc_html_escape(paste0(
      "Not recorded in this saved multiomic object: ", paste(gsub("_", " ", missing_fp_stages), collapse = ", "),
      ". New preprocessing runs retain compact raw and normalized distribution summaries."
    )), "</p>")
  } else {
    ""
  }
  distribution_views <- paste0(
    "<div class=\"qc-tabset\"><div class=\"qc-view-tabs\" role=\"tablist\" aria-label=\"Distribution assay\">",
    "<button type=\"button\" class=\"active\" data-qc-group=\"distribution\" data-qc-view=\"footprints\" aria-selected=\"true\">Footprint scores</button>",
    "<button type=\"button\" data-qc-group=\"distribution\" data-qc-view=\"expression\" aria-selected=\"false\">Gene expression</button></div>",
    "<div class=\"qc-view active\" data-qc-panel-group=\"distribution\" data-qc-panel=\"footprints\">", distribution_notice,
    "<figure class=\"plot-card plot-card-wide\">", .module1_qc_violin_svg(distribution_table, if (has_atac) "ATAC-aligned footprint score distribution per condition" else "Footprint score distribution per condition", width = 1200L), "</figure></div>",
    "<div class=\"qc-view\" data-qc-panel-group=\"distribution\" data-qc-panel=\"expression\"><figure class=\"plot-card plot-card-wide\">",
    .module1_qc_violin_svg(gene_distribution, paste0("Gene expression distribution per condition (", if (nrow(gene_distribution)) gene_distribution$transformation[[1L]] else "unavailable", ")"), width = 1200L),
    "</figure></div></div>"
  )
  per_condition_script <- paste0(
    "<script>(()=>{document.querySelectorAll('[data-qc-view]').forEach(button=>button.addEventListener('click',()=>{",
    "const group=button.dataset.qcGroup,view=button.dataset.qcView;document.querySelectorAll('[data-qc-view]').forEach(x=>{if(x.dataset.qcGroup===group){const active=x.dataset.qcView===view;x.classList.toggle('active',active);x.setAttribute('aria-selected',active?'true':'false')}});",
    "document.querySelectorAll('[data-qc-panel]').forEach(x=>{if(x.dataset.qcPanelGroup===group)x.classList.toggle('active',x.dataset.qcPanel===view)})}))})()</script>"
  )
  sections <- list(
    .qc_section("1. Run Summary", paste0(
      explorer_style,
      "<nav class=\"qc-nav\"><a href=\"#run-summary\">Summary</a><a href=\"#per-condition-qc\">Per-condition QC</a><a href=\"#correlation-summary\">Correlations</a><a href=\"#predicted-binding-sites-per-tf\">Binding</a><a href=\"#co-binding-summary\">Co-binding</a><a href=\"#top-predicted-tfs-at-motif-containing-footprints\">Motifs</a></nav>",
      "<div class=\"report-row\">", .qc_cards_html(summary_cards),
      "<p class=\"metric-note\">Filtered footprints are the canonical-bound footprint set. Predicted unique TFBS counts distinct TF-footprint pairs across the union of all conditions; each pair is counted once.</p>",
      if (has_atac) "" else "<p class=\"metric-note\">ATAC master table not provided; footprint-bound calls and footprint scores were used without an ATAC-open gate.</p>",
      "</div>",
      "<div class=\"report-row\">", setup_html, "</div>"
    )),
    .qc_section("2. Per-Condition QC", paste0(
      "<p class=\"qc-view-intro\">Choose one review mode to keep dense condition-level diagnostics readable.</p><div class=\"qc-tabset\"><div class=\"qc-view-tabs\" role=\"tablist\" aria-label=\"Per-condition QC review mode\">",
      "<button type=\"button\" class=\"active\" data-qc-group=\"per-condition\" data-qc-view=\"distributions\" aria-selected=\"true\">Distributions</button>",
      "<button type=\"button\" data-qc-group=\"per-condition\" data-qc-view=\"pca\" aria-selected=\"false\">PCA</button>",
      "<button type=\"button\" data-qc-group=\"per-condition\" data-qc-view=\"counts\" aria-selected=\"false\">Counts</button></div>",
      "<div class=\"qc-view active\" data-qc-panel-group=\"per-condition\" data-qc-panel=\"distributions\">", distribution_views, "</div>",
      "<div class=\"qc-view\" data-qc-panel-group=\"per-condition\" data-qc-panel=\"pca\"><h3>PCA of samples based on the most variable ", if (has_atac) "RNA, ATAC, and filtered footprints" else "RNA and filtered footprints", "</h3>", pca_plots, "</div>",
      "<div class=\"qc-view\" data-qc-panel-group=\"per-condition\" data-qc-panel=\"counts\"><div class=\"plot-grid\"><figure class=\"plot-card\">", .qc_bar_svg(filtered_by_condition, "condition", "n_filtered_fp", title = "Total filtered footprints detected per condition"), "</figure><figure class=\"plot-card\">", .qc_bar_svg(predicted_by_condition, "condition", "n_predicted_tfbs", title = "Total predicted TFBS per condition"), "</figure></div></div></div>",
      "<details><summary>Per-condition values</summary>",
      .qc_table_html(distribution_summary, max_rows = 300L),
      .qc_table_html(condition_qc$condition_summary, max_rows = 100L), "</details>", per_condition_script
    )),
    .qc_section("3. Correlation Summary", paste0(
      correlation_browser,
      "<details><summary>Correlation values</summary>",
      .qc_table_html(dplyr::bind_rows(canonical_corr$summary, prediction_corr$summary), max_rows = 20L),
      .qc_table_html(canonical_corr$top_tf, max_rows = top_n), "</details>"
    )),
    .qc_section("4. Predicted Binding Sites per TF", paste0(
      "<details class=\"qc-guide\"><summary>How to read this section</summary>", .qc_callout_html("Overall vs condition", c(
        "Overall counts include unique predicted TF-footprint pairs.",
        paste0("The initial condition is ", default_condition_use, "; condition-specific counts require a bound footprint and expressed TF."),
        "Change condition or top N and click a TF to highlight it across both panels."
      )), "</details>",
      if (is.list(explorer_components)) explorer_components$binding else explorer_empty,
      "<details><summary>Prediction summary values</summary>", .qc_table_html(tf_overview, max_rows = 20L), "</details>"
    )),
    .qc_section("5. Co-binding Summary", paste0(
      "<details class=\"qc-guide\"><summary>How to read this section</summary>", .qc_callout_html("Exact co-binding", c(
        "The explorer supports overall or condition-specific co-binding.",
        "For multiple focal TFs, a site must be shared by every selected TF.",
        "Stacked bars separate co-bound sites from the remainder of each partner TF total."
      )), "</details>",
      if (is.list(explorer_components)) explorer_components$cobinding else explorer_empty
    )),
    .qc_section("6. Top Predicted TFs at Motif-Containing Footprints", paste0(
      "<details class=\"qc-guide\"><summary>How to read this section</summary>", .qc_callout_html("Motif-centered review", c(
        "Select a motif to rank predicted TFs found at footprint sites carrying that motif.",
        "The true top 20 are shown first; canonical motif TF references follow with their actual rank or not-predicted status.",
        "Canonical motif identity does not force a TF to rank first because broader prediction uses TF-footprint correlation evidence."
      )), "</details>",
      if (is.list(explorer_components)) explorer_components$motif else explorer_empty,
      if (is.list(explorer_components)) explorer_components$script else ""
    )),
    .qc_section("Technical Appendix", paste0(
      content$review,
      "<details><summary>Workflow gates and footprint alignment</summary>", content$interpretation,
      .qc_table_html(gate_table, max_rows = 20L), alignment_plots,
      .qc_table_html(alignment_summary$counts, max_rows = 20L),
      .qc_table_html(alignment_summary$group_summary, max_rows = 10L), "</details>",
      "<details><summary>Condition support and motif complexity</summary>", condition_support_plots,
      .qc_table_html(condition_overview, max_rows = 20L), footprint_distribution_plots,
      .qc_table_html(motif_complexity$summary, max_rows = 10L), "</details>",
      "<details><summary>Integrity and warning checks</summary>",
      .qc_callout_html("Canonical Support Check Detail", canonical_status$detail, class = if (isTRUE(canonical_status$checked)) "info" else "warn"),
      .qc_table_html(support_check, max_rows = 20L),
      .qc_table_html(predicted_scan$integrity, max_rows = 20L),
      .qc_table_html(manifest_checks, max_rows = 25L),
      .qc_table_html(.qc_status_table(warning_checks), max_rows = 30L), "</details>",
      "<details><summary>Predicted TFBS chunks and related files</summary>",
      .qc_lollipop_svg(pred_manifest, "chunk_id", "n_rows", title = "Predicted TFBS rows per chunk"),
      .qc_table_html(pred_manifest, max_rows = 25L), .qc_links_html(related, from_dir = output_dir), "</details>"
    ))
  )
  out <- file.path(output_dir, report_name)
  path <- .qc_write_html(out, report_title, sections)
  stale_explorer <- file.path(output_dir, "module1_tfbs_explorer.html")
  if (file.exists(stale_explorer)) unlink(stale_explorer)
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
  if (!nrow(hit)) return(tibble::tibble())
  if (identical(as.character(hit$format[[1L]]), "manifest")) {
    return(.qc_read_manifest_chunks(as.character(hit$path[[1L]])))
  }
  hit$chunk_id <- seq_len(nrow(hit))
  tibble::as_tibble(hit[, intersect(c("chunk_id", "path", "format", "n_rows"), names(hit)), drop = FALSE])
}

.module2_qc_scan_candidates <- function(module2, bins = 40L, verbose = TRUE) {
  abs_distance_to_tss <- candidate_source <- distance_to_tss <- fp_id <- median_abs_distance_to_tss <- n_candidates <- n_fp <- n_missing_candidate_source <- n_missing_distance_to_tss <- n_missing_fp_or_target <- n_prior_supported <- n_rows_scanned <- n_target_genes <- p95_abs_distance_to_tss <- prior_supported <- target_gene <- N <- NULL
  man <- .module2_qc_manifest_for(module2, "module2_fp_target_candidates")
  if (!nrow(man)) {
    return(list(
      distance_values = numeric(),
      source_summary = tibble::tibble(),
      top_target = tibble::tibble(),
      top_fp = tibble::tibble(),
      integrity = tibble::tibble()
    ))
  }
  dist_sample <- numeric()
  src_rows <- vector("list", nrow(man))
  target_rows <- vector("list", nrow(man))
  fp_rows <- vector("list", nrow(man))
  integrity_rows <- vector("list", nrow(man))
  for (i in seq_len(nrow(man))) {
    dt <- data.table::as.data.table(.qc_read_table_file(as.character(man$path[[i]]), as.character(man$format[[i]]), columns = c("fp_id", "target_gene", "distance_to_tss", "candidate_source", "prior_supported")))
    if (!nrow(dt)) next
    if (!"prior_supported" %in% names(dt)) dt[, prior_supported := FALSE]
    if (!"candidate_source" %in% names(dt)) dt[, candidate_source := NA_character_]
    d_all <- suppressWarnings(as.numeric(dt$distance_to_tss))
    dt[, abs_distance_to_tss := abs(d_all)]
    src_rows[[i]] <- dt[, .(
      n_candidates = .N,
      n_fp = data.table::uniqueN(fp_id),
      n_target_genes = data.table::uniqueN(target_gene),
      n_prior_supported = sum(prior_supported %in% TRUE, na.rm = TRUE),
      median_abs_distance_to_tss = stats::median(abs_distance_to_tss, na.rm = TRUE),
      p95_abs_distance_to_tss = as.numeric(stats::quantile(abs_distance_to_tss, 0.95, na.rm = TRUE))
    ), by = .(candidate_source)]
    target_rows[[i]] <- dt[, .(n_candidates = .N, n_fp = data.table::uniqueN(fp_id)), by = target_gene]
    fp_rows[[i]] <- dt[, .(n_candidates = .N, n_target_genes = data.table::uniqueN(target_gene)), by = fp_id]
    integrity_rows[[i]] <- data.table::data.table(
      chunk_id = if ("chunk_id" %in% names(man)) man$chunk_id[[i]] else i,
      n_rows_scanned = nrow(dt),
      n_missing_candidate_source = sum(is.na(dt$candidate_source) | !nzchar(dt$candidate_source)),
      n_missing_distance_to_tss = sum(!is.finite(d_all)),
      n_missing_fp_or_target = sum(is.na(dt$fp_id) | !nzchar(dt$fp_id) | is.na(dt$target_gene) | !nzchar(dt$target_gene))
    )
    d <- suppressWarnings(as.numeric(dt$distance_to_tss))
    d <- d[is.finite(d)]
    if (length(d)) dist_sample <- c(dist_sample, if (length(d) > 50000L) utils::head(d, 50000L) else d)
    if (isTRUE(verbose)) .log_inform("Module 2 QC report: scanned candidate chunk {i}/{nrow(man)}.")
  }
  src <- data.table::rbindlist(src_rows, use.names = TRUE, fill = TRUE)
  if (nrow(src)) {
    src <- src[, .(
      n_candidates = sum(n_candidates, na.rm = TRUE),
      n_fp = sum(n_fp, na.rm = TRUE),
      n_target_genes = sum(n_target_genes, na.rm = TRUE),
      n_prior_supported = sum(n_prior_supported, na.rm = TRUE),
      median_abs_distance_to_tss = stats::median(median_abs_distance_to_tss, na.rm = TRUE),
      p95_abs_distance_to_tss = max(p95_abs_distance_to_tss, na.rm = TRUE)
    ), by = candidate_source]
    data.table::setorder(src, -n_candidates, candidate_source)
  }
  top_target <- data.table::rbindlist(target_rows, use.names = TRUE, fill = TRUE)
  if (nrow(top_target)) {
    top_target <- top_target[, .(n_candidates = sum(n_candidates, na.rm = TRUE), n_fp = sum(n_fp, na.rm = TRUE)), by = target_gene]
    data.table::setorder(top_target, -n_candidates, target_gene)
    top_target <- head(top_target, 20L)
  }
  top_fp <- data.table::rbindlist(fp_rows, use.names = TRUE, fill = TRUE)
  if (nrow(top_fp)) {
    top_fp <- top_fp[, .(n_candidates = sum(n_candidates, na.rm = TRUE), n_target_genes = sum(n_target_genes, na.rm = TRUE)), by = fp_id]
    data.table::setorder(top_fp, -n_candidates, fp_id)
    top_fp <- head(top_fp, 20L)
  }
  integrity <- data.table::rbindlist(integrity_rows, use.names = TRUE, fill = TRUE)
  if (nrow(integrity)) {
    integrity <- integrity[, .(
      n_chunks_scanned = .N,
      n_rows_scanned = sum(n_rows_scanned, na.rm = TRUE),
      n_missing_candidate_source = sum(n_missing_candidate_source, na.rm = TRUE),
      n_missing_distance_to_tss = sum(n_missing_distance_to_tss, na.rm = TRUE),
      n_missing_fp_or_target = sum(n_missing_fp_or_target, na.rm = TRUE)
    )]
  }
  list(
    distance_values = dist_sample,
    source_summary = tibble::as_tibble(src),
    top_target = tibble::as_tibble(top_target),
    top_fp = tibble::as_tibble(top_fp),
    integrity = tibble::as_tibble(integrity)
  )
}

.module2_qc_scan_links <- function(module2, top_n = 20L, validate_integrity = TRUE, verbose = TRUE) {
  fp_id <- module2_link_pass <- n_fp <- n_links <- n_target_genes <- n_tf <- pass <- target_gene <- tf <- NULL
  link_man <- .module2_qc_manifest_for(module2, "module2_links")
  if (!nrow(link_man) && is.data.frame(module2$links) && nrow(module2$links)) {
    dt <- data.table::as.data.table(module2$links)
    top <- dt[, .(n_links = .N, n_fp = data.table::uniqueN(fp_id), n_target_genes = data.table::uniqueN(target_gene)), by = tf]
    data.table::setorder(top, -n_links, tf)
    top_target <- dt[, .(n_links = .N, n_tf = data.table::uniqueN(tf), n_fp = data.table::uniqueN(fp_id)), by = target_gene]
    data.table::setorder(top_target, -n_links, target_gene)
    top_fp <- dt[, .(n_links = .N, n_tf = data.table::uniqueN(tf), n_target_genes = data.table::uniqueN(target_gene)), by = fp_id]
    data.table::setorder(top_fp, -n_links, fp_id)
    return(list(
      top_tf = tibble::as_tibble(head(top, top_n)),
      top_target = tibble::as_tibble(head(top_target, top_n)),
      top_fp = tibble::as_tibble(head(top_fp, top_n)),
      validation = tibble::tibble(metric = "n_links_scanned", value = nrow(dt))
    ))
  }
  if (!nrow(link_man)) return(list(top_tf = tibble::tibble(), top_target = tibble::tibble(), top_fp = tibble::tibble(), validation = tibble::tibble()))
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
  target_rows <- vector("list", nrow(link_man))
  fp_rows <- vector("list", nrow(link_man))
  n_scanned <- 0
  n_bad_tf <- 0
  n_bad_fp <- 0
  for (i in seq_len(nrow(link_man))) {
    dt <- data.table::as.data.table(.qc_read_table_file(as.character(link_man$path[[i]]), as.character(link_man$format[[i]]), columns = c("tf", "fp_id", "target_gene", "module2_link_pass")))
    if (!nrow(dt)) next
    dt <- dt[module2_link_pass %in% TRUE]
    n_scanned <- n_scanned + nrow(dt)
    rows[[i]] <- dt[, .(n_links = .N, n_fp = data.table::uniqueN(fp_id), n_target_genes = data.table::uniqueN(target_gene)), by = tf]
    target_rows[[i]] <- dt[, .(n_links = .N, n_tf = data.table::uniqueN(tf), n_fp = data.table::uniqueN(fp_id)), by = target_gene]
    fp_rows[[i]] <- dt[, .(n_links = .N, n_tf = data.table::uniqueN(tf), n_target_genes = data.table::uniqueN(target_gene)), by = fp_id]
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
  top_target <- data.table::rbindlist(target_rows, use.names = TRUE, fill = TRUE)
  if (nrow(top_target)) {
    top_target <- top_target[, .(n_links = sum(n_links, na.rm = TRUE), n_tf = sum(n_tf, na.rm = TRUE), n_fp = sum(n_fp, na.rm = TRUE)), by = target_gene]
    data.table::setorder(top_target, -n_links, target_gene)
    top_target <- head(top_target, top_n)
  }
  top_fp <- data.table::rbindlist(fp_rows, use.names = TRUE, fill = TRUE)
  if (nrow(top_fp)) {
    top_fp <- top_fp[, .(n_links = sum(n_links, na.rm = TRUE), n_tf = sum(n_tf, na.rm = TRUE), n_target_genes = sum(n_target_genes, na.rm = TRUE)), by = fp_id]
    data.table::setorder(top_fp, -n_links, fp_id)
    top_fp <- head(top_fp, top_n)
  }
  validation <- tibble::tibble(
    metric = c("n_links_scanned", "n_links_with_missing_tf_target_pass", "n_links_with_missing_fp_target_pass"),
    value = c(n_scanned, n_bad_tf, n_bad_fp)
  )
  list(top_tf = tibble::as_tibble(top), top_target = tibble::as_tibble(top_target), top_fp = tibble::as_tibble(top_fp), validation = validation)
}

.module2_qc_related_html <- function(module2_dir, output_dir) {
  if (is.null(module2_dir) || !dir.exists(module2_dir)) return(character())
  hits <- list.files(module2_dir, pattern = "\\.html$", full.names = TRUE, recursive = TRUE)
  hits[normalizePath(hits, winslash = "/", mustWork = FALSE) != normalizePath(file.path(output_dir, "module2_qc_report.html"), winslash = "/", mustWork = FALSE)]
}

.module2_qc_parameters <- function(module2, scan_large_tables = FALSE, validate_integrity = FALSE) {
  params <- if (is.list(module2)) module2$parameters else NULL
  keep <- c("max_distance_bp", "n_cores", "output_format", "streamed", "write_qc_report", "qc_report_validate")
  out <- params[intersect(keep, names(params))]
  out$scan_large_tables <- isTRUE(scan_large_tables)
  out$validate_integrity <- isTRUE(validate_integrity)
  .qc_key_value_table(out)
}

.module2_qc_input_handoff <- function(module2, multiomic_data = NULL, qc_summary = NULL) {
  omics_cards <- tibble::tibble(label = character(), value = character())
  if (is_multiomic_object(multiomic_data)) {
    mats <- multiomic_data$matrices
    safe_nrow <- function(x) if (is.null(x)) NA_integer_ else nrow(x)
    safe_ncol <- function(x) if (is.null(x)) NA_integer_ else ncol(x)
    omics_cards <- tibble::tibble(
      label = c("Conditions", "FP score rows", "Gene expression rows", "ATAC score rows"),
      value = c(
        .qc_format_number(safe_ncol(mats$fp_score)),
        .qc_format_number(safe_nrow(mats$fp_score)),
        .qc_format_number(safe_nrow(mats$gene_expr)),
        .qc_format_number(safe_nrow(mats$atac_score))
      )
    )
  }
  module_cards <- tibble::tibble(
    label = c("Predicted TFBS rows", "TFs", "Target genes", "Final links"),
    value = c(
      .qc_format_number(.qc_metric_value(qc_summary, "n_predicted_tfbs")),
      .qc_format_number(.qc_metric_value(qc_summary, "n_tfs")),
      .qc_format_number(.qc_metric_value(qc_summary, "n_target_genes")),
      .qc_format_number(.qc_metric_value(qc_summary, "n_module2_links"))
    )
  )
  rbind(module_cards, omics_cards)
}

.module2_qc_corr_summary <- function(module2, table_name, label, group_cols = character(), top_n = 20L, scan = TRUE) {
  best_r <- n_pass <- pass <- NULL
  if (!isTRUE(scan)) {
    return(list(summary = tibble::tibble(), top = tibble::tibble(), best_r = numeric()))
  }
  x <- .module2_qc_table_from_module(
    module2,
    table_name,
    columns = unique(c(group_cols, "pearson_r", "pearson_p", "pearson_fdr", "spearman_r", "spearman_p", "spearman_fdr", "best_r", "best_p", "best_fdr", "best_method", "pass"))
  )
  if (!is.data.frame(x) || !nrow(x)) {
    return(list(summary = tibble::tibble(), top = tibble::tibble(), best_r = numeric()))
  }
  dt <- data.table::as.data.table(x)
  if (!"pass" %in% names(dt)) dt[, pass := FALSE]
  if (!"best_r" %in% names(dt)) dt[, best_r := NA_real_]
  dt[, best_r := suppressWarnings(as.numeric(best_r))]
  summary <- tibble::tibble(
    table = label,
    n_rows = nrow(dt),
    n_pass = sum(dt$pass %in% TRUE, na.rm = TRUE),
    pass_rate = .qc_percent(sum(dt$pass %in% TRUE, na.rm = TRUE), nrow(dt)),
    median_best_r = stats::median(dt$best_r, na.rm = TRUE),
    p95_best_r = as.numeric(stats::quantile(dt$best_r, 0.95, na.rm = TRUE)),
    max_best_r = suppressWarnings(max(dt$best_r, na.rm = TRUE))
  )
  top <- tibble::tibble()
  first_group <- group_cols[[1L]] %||% NA_character_
  second_group <- group_cols[[2L]] %||% NA_character_
  if (!is.na(first_group) && first_group %in% names(dt)) {
    by_cols <- first_group
    top_dt <- dt[pass %in% TRUE, .(
      n_pass = .N,
      n_unique_partner = if (!is.na(second_group) && second_group %in% names(dt)) data.table::uniqueN(get(second_group)) else NA_integer_,
      median_best_r = stats::median(best_r, na.rm = TRUE)
    ), by = by_cols]
    if (nrow(top_dt)) {
      data.table::setorder(top_dt, -n_pass)
      top <- tibble::as_tibble(head(top_dt, top_n))
    }
  }
  list(summary = summary, top = top, best_r = dt$best_r[is.finite(dt$best_r)])
}

.module2_qc_empty_condition_activity <- function(source = "not_run") {
  list(
    summary = tibble::tibble(),
    tf_summary = tibble::tibble(),
    tf_matrix = tibble::tibble(),
    signal_summary = tibble::tibble(),
    source = source
  )
}

.module2_qc_stream_condition_activity <- function(module2, multiomic_data = NULL, top_n = 20L, verbose = TRUE) {
  condition <- fp_id <- module2_link_pass <- n_active <- target_gene <- tf <- NULL
  if (!is_multiomic_object(multiomic_data)) return(.module2_qc_empty_condition_activity("missing_multiomic_data"))
  link_man <- .module2_qc_manifest_for(module2, "module2_links")
  if (!nrow(link_man)) return(.module2_qc_empty_condition_activity("missing_link_manifest"))
  mats <- multiomic_data$matrices
  conditions <- .qc_condition_columns(multiomic_data)
  gene_on <- if (!is.null(mats$gene_on)) .qc_as_numeric_matrix(mats$gene_on) else .qc_as_numeric_matrix(mats$gene_expr) > 0
  fp_bound <- .qc_as_numeric_matrix(mats$fp_bound)
  fp_score <- .qc_as_numeric_matrix(mats$fp_score)
  if (!nrow(gene_on) || !nrow(fp_bound) || !nrow(fp_score)) {
    return(.module2_qc_empty_condition_activity("missing_matrices"))
  }
  conditions <- conditions[
    conditions %in% colnames(gene_on) &
      conditions %in% colnames(fp_bound) &
      conditions %in% colnames(fp_score)
  ]
  if (!length(conditions)) return(.module2_qc_empty_condition_activity("missing_condition_columns"))
  gene_cols <- match(conditions, colnames(gene_on))
  fp_cols <- match(conditions, colnames(fp_bound))
  score_cols <- match(conditions, colnames(fp_score))
  active_counts <- stats::setNames(numeric(length(conditions)), conditions)
  n_scanned <- 0
  tf_rows <- vector("list", length(conditions) * nrow(link_man))
  k <- 1L
  for (i in seq_len(nrow(link_man))) {
    dt <- data.table::as.data.table(.qc_read_table_file(
      as.character(link_man$path[[i]]),
      as.character(link_man$format[[i]]),
      columns = c("tf", "fp_id", "target_gene", "module2_link_pass")
    ))
    if (!nrow(dt)) next
    if (!"module2_link_pass" %in% names(dt)) dt[, module2_link_pass := TRUE]
    dt <- dt[module2_link_pass %in% TRUE]
    if (!nrow(dt)) next
    n_scanned <- n_scanned + nrow(dt)
    tf_idx <- match(dt$tf, rownames(gene_on))
    target_idx <- match(dt$target_gene, rownames(gene_on))
    fp_idx <- match(dt$fp_id, rownames(fp_bound))
    valid <- !is.na(tf_idx) & !is.na(target_idx) & !is.na(fp_idx)
    for (j in seq_along(conditions)) {
      tf_on <- gene_on[tf_idx, gene_cols[[j]]] > 0
      target_on <- gene_on[target_idx, gene_cols[[j]]] > 0
      fp_is_bound <- fp_bound[fp_idx, fp_cols[[j]]] > 0
      fp_has_score <- fp_score[fp_idx, score_cols[[j]]] > 0
      active <- valid & tf_on & target_on & fp_is_bound & fp_has_score
      active[is.na(active)] <- FALSE
      n_active_j <- sum(active)
      active_counts[[j]] <- active_counts[[j]] + n_active_j
      if (n_active_j > 0L) {
        tf_rows[[k]] <- data.table::data.table(tf = dt$tf[active], condition = conditions[[j]])[
          ,
          .(n_active = .N),
          by = .(tf, condition)
        ]
        k <- k + 1L
      }
    }
    if (isTRUE(verbose)) .log_inform("Module 2 QC report: scanned condition activity chunk {i}/{nrow(link_man)}.")
  }
  summary <- tibble::tibble(
    condition = conditions,
    n_rows = n_scanned,
    n_active = as.numeric(active_counts),
    active_rate = vapply(active_counts, .qc_percent, character(1L), den = n_scanned)
  )
  summary <- summary[order(-summary$n_active, summary$condition), , drop = FALSE]
  tf_long <- data.table::rbindlist(tf_rows, use.names = TRUE, fill = TRUE)
  tf_summary <- tibble::tibble()
  tf_matrix <- tibble::tibble()
  if (nrow(tf_long)) {
    tf_summary_dt <- tf_long[, .(n_active = sum(n_active, na.rm = TRUE), n_active_conditions = data.table::uniqueN(condition[n_active > 0])), by = tf]
    data.table::setorder(tf_summary_dt, -n_active, tf)
    tf_summary <- tibble::as_tibble(head(tf_summary_dt, top_n))
    keep_tf <- tf_summary$tf
    tf_matrix_dt <- tf_long[tf %in% keep_tf, .(n_active = sum(n_active, na.rm = TRUE)), by = .(tf, condition)]
    tf_matrix_dt <- data.table::dcast(tf_matrix_dt, tf ~ condition, value.var = "n_active", fill = 0)
    for (cc in conditions) {
      if (!cc %in% names(tf_matrix_dt)) tf_matrix_dt[[cc]] <- 0
    }
    data.table::setorderv(tf_matrix_dt, "tf")
    ord <- match(tf_matrix_dt$tf, keep_tf)
    tf_matrix_dt <- tf_matrix_dt[order(ord)]
    tf_matrix <- tibble::as_tibble(tf_matrix_dt[, c("tf", conditions), with = FALSE])
  }
  signal_summary <- tibble::tibble(
    metric = c("links scanned", "conditions scanned", "active condition-link events", "source"),
    value = c(
      .qc_format_number(n_scanned),
      .qc_format_number(length(conditions)),
      .qc_format_number(sum(active_counts, na.rm = TRUE)),
      "streamed link scan"
    )
  )
  list(
    summary = tibble::as_tibble(summary),
    tf_summary = tf_summary,
    tf_matrix = tf_matrix,
    signal_summary = signal_summary,
    source = "streamed"
  )
}

.module2_qc_condition_activity <- function(module2, multiomic_data = NULL, scan = TRUE, top_n = 20L, verbose = TRUE) {
  active <- condition <- n_active <- NULL
  if (!isTRUE(scan)) return(.module2_qc_empty_condition_activity("not_run"))
  activity <- .module2_qc_table_from_module(
    module2,
    "module2_condition_activity",
    columns = c("condition", "active", "tf_expr", "target_expr", "fp_score", "fp_bound", "atac_score")
  )
  if (!is.data.frame(activity) || !nrow(activity)) {
    return(.module2_qc_stream_condition_activity(module2, multiomic_data = multiomic_data, top_n = top_n, verbose = verbose))
  }
  dt <- data.table::as.data.table(activity)
  if (!"active" %in% names(dt)) dt[, active := FALSE]
  dt[, active := active %in% TRUE]
  summary <- dt[, .(
    n_rows = .N,
    n_active = sum(active, na.rm = TRUE),
    active_rate = .qc_percent(sum(active, na.rm = TRUE), .N)
  ), by = condition]
  data.table::setorder(summary, -n_active, condition)
  signal_cols <- intersect(c("tf_expr", "target_expr", "fp_score", "atac_score"), names(dt))
  signal_summary <- tibble::tibble(metric = character(), median_active_value = numeric())
  if (length(signal_cols)) {
    signal_summary <- tibble::tibble(
      metric = signal_cols,
      median_active_value = vapply(signal_cols, function(nm) {
        stats::median(suppressWarnings(as.numeric(dt[active %in% TRUE][[nm]])), na.rm = TRUE)
      }, numeric(1L))
    )
  }
  list(
    summary = tibble::as_tibble(head(summary, top_n)),
    tf_summary = tibble::tibble(),
    tf_matrix = tibble::tibble(),
    signal_summary = signal_summary,
    source = "materialized"
  )
}

.module2_qc_warning_checks <- function(qc_summary, manifest_checks, candidate_scan, link_scan, tf_corr_scan, fp_corr_scan, condition_scan) {
  bad_tf <- .qc_metric_value(link_scan$validation, "n_links_with_missing_tf_target_pass", default = NA_real_)
  bad_fp <- .qc_metric_value(link_scan$validation, "n_links_with_missing_fp_target_pass", default = NA_real_)
  missing_manifest <- if (is.data.frame(manifest_checks) && nrow(manifest_checks)) {
    sum(as.character(manifest_checks$file_exists) != "yes", na.rm = TRUE)
  } else {
    NA_real_
  }
  missing_source <- .qc_metric_value(candidate_scan$integrity, "n_missing_candidate_source", default = NA_real_)
  missing_dist <- .qc_metric_value(candidate_scan$integrity, "n_missing_distance_to_tss", default = NA_real_)
  tf_pass <- .qc_metric_value(qc_summary, "n_tf_target_pairs_pass", default = NA_real_)
  fp_pass <- .qc_metric_value(qc_summary, "n_fp_target_pairs_pass", default = NA_real_)
  n_links <- .qc_metric_value(qc_summary, "n_module2_links", default = NA_real_)
  status_rows <- list(
    list("Output manifest files exist", .qc_status(if (is.finite(missing_manifest)) missing_manifest == 0 else NA, paste0("Missing files: ", .qc_format_number(missing_manifest)))),
    list("Final links have passing TF-target support", .qc_status(if (is.finite(bad_tf)) bad_tf == 0 else NA, paste0("Missing support rows: ", .qc_format_number(bad_tf)))),
    list("Final links have passing FP-target support", .qc_status(if (is.finite(bad_fp)) bad_fp == 0 else NA, paste0("Missing support rows: ", .qc_format_number(bad_fp)))),
    list("Candidate source is annotated", .qc_status(if (is.finite(missing_source)) missing_source == 0 else NA, paste0("Missing candidate_source rows: ", .qc_format_number(missing_source)))),
    list("Candidate distance to TSS is available", .qc_status(if (is.finite(missing_dist)) missing_dist == 0 else NA, paste0("Missing distance rows: ", .qc_format_number(missing_dist)))),
    list("TF-target correlation produced passing pairs", .qc_status(if (is.finite(tf_pass)) tf_pass > 0 else nrow(tf_corr_scan$summary) > 0, paste0("Passing pairs: ", .qc_format_number(tf_pass)))),
    list("FP-target correlation produced passing pairs", .qc_status(if (is.finite(fp_pass)) fp_pass > 0 else nrow(fp_corr_scan$summary) > 0, paste0("Passing pairs: ", .qc_format_number(fp_pass)))),
    list("Final Module 2 links are present", .qc_status(if (is.finite(n_links)) n_links > 0 else NA, paste0("Final links: ", .qc_format_number(n_links)))),
    list(
      "Condition activity table is available",
      if (is.data.frame(condition_scan$summary) && nrow(condition_scan$summary) > 0) {
        .qc_status(TRUE, paste0("Condition rows: ", .qc_format_number(nrow(condition_scan$summary)), "; source: ", condition_scan$source %||% "unknown"))
      } else if (is.finite(.qc_metric_value(qc_summary, "n_active_link_conditions", default = NA_real_))) {
        .qc_status(FALSE, paste0("Condition rows: ", .qc_format_number(nrow(condition_scan$summary))))
      } else {
        list(status = "SKIPPED", detail = "Condition activity could not be computed because the materialized table or streamed output path inputs were unavailable.")
      }
    )
  )
  tibble::tibble(
    check = vapply(status_rows, `[[`, character(1L), 1L),
    status = vapply(status_rows, function(x) x[[2L]]$status, character(1L)),
    detail = vapply(status_rows, function(x) x[[2L]]$detail, character(1L))
  )
}

.module2_qc_content_sections <- function(qc_summary,
                                         link_scan,
                                         warning_checks,
                                         condition_scan = NULL) {
  n_predicted <- .qc_metric_value(qc_summary, "n_predicted_tfbs")
  n_tf_test <- .qc_metric_value(qc_summary, "n_tf_target_pairs_tested")
  n_tf_pass <- .qc_metric_value(qc_summary, "n_tf_target_pairs_pass")
  n_fp_candidates <- .qc_metric_value(qc_summary, "n_fp_target_candidates")
  n_fp_test <- .qc_metric_value(qc_summary, "n_fp_target_pairs_tested")
  n_fp_pass <- .qc_metric_value(qc_summary, "n_fp_target_pairs_pass")
  n_links <- .qc_metric_value(qc_summary, "n_module2_links")
  bad_tf <- .qc_metric_value(link_scan$validation, "n_links_with_missing_tf_target_pass", default = NA_real_)
  bad_fp <- .qc_metric_value(link_scan$validation, "n_links_with_missing_fp_target_pass", default = NA_real_)
  n_check <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "CHECK", na.rm = TRUE) else NA_integer_
  n_not_run <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "NOT RUN", na.rm = TRUE) else NA_integer_
  n_skipped <- if (is.data.frame(warning_checks) && nrow(warning_checks)) sum(warning_checks$status == "SKIPPED", na.rm = TRUE) else NA_integer_
  top_tf <- if (is.data.frame(link_scan$top_tf) && nrow(link_scan$top_tf)) {
    paste0(
      as.character(link_scan$top_tf$tf[[1L]]), " has ",
      .qc_format_number(link_scan$top_tf$n_links[[1L]]),
      " final links in the scanned output."
    )
  } else {
    "Top final-link TF summary was not scanned."
  }
  top_condition <- if (is.list(condition_scan) && is.data.frame(condition_scan$summary) && nrow(condition_scan$summary)) {
    paste0(
      "Condition activity scan found ",
      .qc_format_number(sum(suppressWarnings(as.numeric(condition_scan$summary$n_active)), na.rm = TRUE)),
      " active condition-link events in the displayed conditions; ",
      as.character(condition_scan$summary$condition[[1L]]),
      " had the highest active-link count."
    )
  } else {
    "Condition activity was not scanned or could not be computed from the available output files."
  }
  key <- .qc_callout_html("What Happened", c(
    paste0("Module 2 consumed ", .qc_format_number(n_predicted), " predicted TFBS rows from Module 1."),
    .qc_metric_sentence("TF-target prefilter", n_tf_pass, n_tf_test),
    paste0("FP-target candidate construction produced ", .qc_format_number(n_fp_candidates), " candidate rows and tested ", .qc_format_number(n_fp_test), " FP-target correlations."),
    .qc_metric_sentence("FP-target evidence", n_fp_pass, n_fp_test),
    paste0("Final assembly wrote ", .qc_format_number(n_links), " TF-FP-target links after requiring predicted TFBS, passing TF-target support, and passing FP-target support."),
    top_tf,
    top_condition
  ))
  interp <- .qc_callout_html("Module 2 Interpretation", c(
    "Module 2 is intentionally ordered to reduce over-calculation: TF-target expression correlation is the first prefilter, then FP-target evidence is tested only for candidate FP-gene pairs.",
    "FP-target evidence combines the configured TSS window or regulatory prior with FP score to target expression correlation.",
    "Final links represent a relational join: a TF is predicted to bind an FP from Module 1, the TF-target expression pair passes, and the FP-target pair passes.",
    paste0("Relational integrity check found ", .qc_format_number(bad_tf), " links without TF-target support and ", .qc_format_number(bad_fp), " links without FP-target support.")
  ))
  review_class <- if (is.finite(n_check) && n_check > 0) "warn" else "info"
  review <- .qc_callout_html("Review Before Downstream Use", c(
    paste0("Warning Checks contains ", .qc_format_number(n_check), " CHECK rows, ", .qc_format_number(n_not_run), " NOT RUN rows, and ", .qc_format_number(n_skipped), " SKIPPED optional rows."),
    "Review Candidate Source QC to confirm whether links mainly come from the TSS window, regulatory prior input, or both.",
    "Review Candidate Distance To TSS to make sure the selected window and regulatory prior match the expected biology.",
    "Review Top TFs In Final Links and Condition Activity QC for overly broad regulators or condition imbalance before using Module 2 output in Module 3."
  ), class = review_class)
  list(key = key, interpretation = interp, review = review)
}

#' Build a Module 2 QC HTML report
#'
#' Builds a comprehensive HTML report for Module 2 run parameters, CraftGRN multiomic input
#' handoff, TF-target and FP-target correlation filters, candidate source and
#' distance-to-TSS evidence, final TF-FP-target links, condition activity,
#' CraftGRN multiomic condition context, warning checks, integrity checks, and
#' related browser reports.
#'
#' @param module2 Module 2 result list, loaded Module 2 list, or output
#'   directory.
#' @param multiomic_data Optional CraftGRN multiomic object used for context.
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
  candidate_scan <- if (isTRUE(scan_large_tables)) {
    .module2_qc_scan_candidates(module2, verbose = verbose)
  } else {
    list(distance_values = numeric(), source_summary = tibble::tibble(), top_target = tibble::tibble(), top_fp = tibble::tibble(), integrity = tibble::tibble())
  }
  link_scan <- if (isTRUE(scan_large_tables)) {
    .module2_qc_scan_links(module2, top_n = top_n, validate_integrity = validate_integrity, verbose = verbose)
  } else {
    list(top_tf = tibble::tibble(), top_target = tibble::tibble(), top_fp = tibble::tibble(), validation = tibble::tibble())
  }
  tf_corr_scan <- .module2_qc_corr_summary(
    module2,
    table_name = "module2_tf_target_corr",
    label = "TF expression to target expression",
    group_cols = c("tf", "target_gene"),
    top_n = top_n,
    scan = isTRUE(scan_large_tables)
  )
  fp_corr_scan <- .module2_qc_corr_summary(
    module2,
    table_name = "module2_fp_target_corr",
    label = "FP score to target expression",
    group_cols = c("target_gene", "fp_id"),
    top_n = top_n,
    scan = isTRUE(scan_large_tables)
  )
  condition_scan <- .module2_qc_condition_activity(
    module2,
    multiomic_data = multiomic_data,
    scan = isTRUE(scan_large_tables),
    top_n = top_n,
    verbose = verbose
  )
  parameter_table <- .module2_qc_parameters(module2, scan_large_tables = scan_large_tables, validate_integrity = validate_integrity)
  input_cards <- .module2_qc_input_handoff(module2, multiomic_data = multiomic_data, qc_summary = qc_summary)
  condition_qc <- .qc_multiomic_condition_qc(multiomic_data)
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
  warning_checks <- .module2_qc_warning_checks(
    qc_summary = qc_summary,
    manifest_checks = manifest_checks,
    candidate_scan = candidate_scan,
    link_scan = link_scan,
    tf_corr_scan = tf_corr_scan,
    fp_corr_scan = fp_corr_scan,
    condition_scan = condition_scan
  )
  content <- .module2_qc_content_sections(
    qc_summary = qc_summary,
    link_scan = link_scan,
    warning_checks = warning_checks,
    condition_scan = condition_scan
  )
  related <- .module2_qc_related_html(module2_dir, output_dir)
  tf_corr_plots <- .qc_plot_grid(
    .qc_density_svg(tf_corr_scan$best_r, title = "TF-target best r density"),
    .qc_scatter_svg(
      tf_corr_scan$top,
      x_col = "n_unique_partner",
      y_col = "n_pass",
      label_col = "tf",
      size_col = "median_best_r",
      title = "TF-target passing breadth"
    ),
    .qc_metric_heatmap_svg(
      tf_corr_scan$top,
      row_col = "tf",
      value_cols = c("n_pass", "n_unique_partner", "median_best_r"),
      title = "Top TF-target correlation metrics"
    )
  )
  fp_corr_plots <- .qc_plot_grid(
    .qc_density_svg(fp_corr_scan$best_r, title = "FP-target best r density"),
    .qc_scatter_svg(
      fp_corr_scan$top,
      x_col = "n_unique_partner",
      y_col = "n_pass",
      label_col = "target_gene",
      size_col = "median_best_r",
      title = "Target genes with passing FP evidence"
    ),
    .qc_metric_heatmap_svg(
      fp_corr_scan$top,
      row_col = "target_gene",
      value_cols = c("n_pass", "n_unique_partner", "median_best_r"),
      title = "Top FP-target target-gene metrics"
    )
  )
  candidate_plots <- .qc_plot_grid(
    .qc_density_svg(candidate_scan$distance_values, title = "Candidate signed distance density"),
    .qc_cumulative_svg(abs(candidate_scan$distance_values), title = "Cumulative absolute distance to TSS"),
    .qc_metric_heatmap_svg(
      candidate_scan$source_summary,
      row_col = "candidate_source",
      value_cols = c("n_candidates", "n_fp", "n_target_genes", "n_prior_supported", "median_abs_distance_to_tss", "p95_abs_distance_to_tss"),
      title = "Candidate source evidence matrix"
    )
  )
  final_link_plots <- .qc_plot_grid(
    .qc_scatter_svg(
      link_scan$top_tf,
      x_col = "n_target_genes",
      y_col = "n_links",
      label_col = "tf",
      size_col = "n_fp",
      title = "Final-link breadth by TF"
    ),
    .qc_metric_heatmap_svg(
      link_scan$top_tf,
      row_col = "tf",
      value_cols = c("n_links", "n_fp", "n_target_genes"),
      title = "Top final-link TF metrics"
    ),
    .qc_lollipop_svg(link_scan$top_tf, "tf", "n_links", title = "Top TFs by final links")
  )
  sections <- list(
    .qc_section("Run Parameters", .qc_table_html(parameter_table, max_rows = 30L)),
    .qc_section("Summary", .qc_cards_html(cards)),
    .qc_section("Key Findings", content$key),
    .qc_section("Input Handoff", .qc_cards_html(input_cards)),
    .qc_section("Condition Context", .qc_condition_qc_html(condition_qc, top_n = max(30L, as.integer(top_n[[1L]])))),
    .qc_section("Workflow Funnel", paste0(
      content$interpretation,
      .qc_flow_svg(funnel, "step", "n", title = "Module 2 relational flow")
    )),
    .qc_section("TF-Target Correlation QC", paste0(
      .qc_table_html(tf_corr_scan$summary, max_rows = 10L),
      tf_corr_plots,
      .qc_table_html(tf_corr_scan$top, max_rows = top_n)
    )),
    .qc_section("Candidate Source QC", paste0(
      candidate_plots,
      .qc_table_html(candidate_scan$source_summary, max_rows = 20L),
      .qc_table_html(candidate_scan$integrity, max_rows = 10L)
    )),
    .qc_section("FP-Target Correlation QC", paste0(
      .qc_table_html(fp_corr_scan$summary, max_rows = 10L),
      fp_corr_plots,
      .qc_table_html(fp_corr_scan$top, max_rows = top_n)
    )),
    .qc_section("Integrity Checks", paste0(
      "<p>Final-link relational integrity: <span class=\"", integrity_class, "\">", integrity_status, "</span></p>",
      .qc_table_html(link_scan$validation, max_rows = 20L),
      .qc_table_html(manifest_checks, max_rows = 20L)
    )),
    .qc_section("Candidate Distance To TSS", paste0(
      .qc_plot_grid(
        .qc_density_svg(candidate_scan$distance_values, title = "Signed FP to target TSS distance"),
        .qc_cumulative_svg(abs(candidate_scan$distance_values), title = "Absolute distance cumulative curve")
      ),
      .qc_table_html(candidate_scan$source_summary, max_rows = 20L)
    )),
    .qc_section("Top TFs In Final Links", paste0(
      final_link_plots,
      .qc_table_html(link_scan$top_tf, max_rows = top_n)
    )),
    .qc_section("Top Target Genes In Final Links", paste0(
      .qc_plot_grid(
        .qc_lollipop_svg(link_scan$top_target, "target_gene", "n_links", title = "Top target genes by final links")
      ),
      .qc_table_html(link_scan$top_target, max_rows = top_n)
    )),
    .qc_section("Recommended Review", content$review),
    .qc_section("Condition Activity QC", paste0(
      .qc_callout_html("Link Activity Summary", c(
        "This section replaces the old link_activity_summary.pdf with native HTML plots.",
        "An active condition-link event requires the final Module 2 link to pass, the TF and target gene to be expressed in that condition, and the FP to be bound with positive FP score.",
        "Active link counts per TF per condition are reported when condition activity can be computed from the CraftGRN multiomic object."
      )),
      .qc_plot_grid(
        .qc_lollipop_svg(condition_scan$summary, "condition", "n_active", title = "Total active links per condition")
      ),
      .qc_plot_card_html(
        .qc_matrix_heatmap_svg(
          condition_scan$tf_matrix,
          "tf",
          setdiff(names(condition_scan$tf_matrix), "tf"),
          title = "Active link counts per TF per condition"
        ),
        wide = TRUE
      ),
      .qc_table_html(condition_scan$summary, max_rows = top_n),
      .qc_table_html(condition_scan$tf_summary, max_rows = top_n),
      .qc_table_html(condition_scan$signal_summary, max_rows = 20L)
    )),
    .qc_section("Warning Checks", .qc_table_html(.qc_status_table(warning_checks), max_rows = 30L)),
    .qc_section("QC Summary Table", .qc_table_html(qc_summary, max_rows = 50L)),
    .qc_section("Related HTML Reports", .qc_links_html(related, from_dir = output_dir))
  )
  out <- file.path(output_dir, report_name)
  path <- .qc_write_html(out, "Module 2 QC Report", sections)
  if (isTRUE(verbose)) .log_inform("Module 2 QC HTML report written: {path}")
  path
}
