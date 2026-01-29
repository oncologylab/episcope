#' - New: preserves the unscaled input matrix and makes it available (in heatmap order) as `matrix_ordered`.
#' - Changed default clustering method to `"ward.D2"`
#' - Updated default `border_color` to `NA` (no cell borders).
#' - New: `annotation_row_side` ("left" [default] or "right") to place row annotations.


#' Pretty heatmaps
#'
#' Implementation of heatmaps that offers more control
#' over dimensions and appearance.
#'
#' @name ggheat-package
#'
#' @import grid
#' @import RColorBrewer
#' @import scales
#' @import gtable
#' @importFrom stats hclust dist cutree kmeans sd cor as.dist
#' @import grDevices
#' @import graphics
#'
#' @keywords internal
"_PACKAGE"

library(grid)
library(RColorBrewer)
library(scales)
library(gtable)
library(grDevices)
library(graphics)

lo <- function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, annotation_names_row, annotation_names_col, main, fontsize, fontsize_row, fontsize_col, angle_col, gaps_row, gaps_col, annotation_row_side = c("left","right"), ...) {

  annotation_row_side <- match.arg(annotation_row_side)

  # Get height of colnames and length of rownames
  if (!is.null(coln[1]) | (!is.na2(annotation_row) & annotation_names_row)) {
    if (!is.null(coln[1])) {
      t <- coln
    } else {
      t <- ""
    }
    tw <- strwidth(t, units = "in", cex = fontsize_col / fontsize)
    if (annotation_names_row) {
      t <- c(t, colnames(annotation_row))
      tw <- c(tw, strwidth(colnames(annotation_row), units = "in"))
    }
    longest_coln <- which.max(tw)
    gp <- list(fontsize = ifelse(longest_coln <= length(coln), fontsize_col, fontsize), ...)
    coln_height <- unit(1, "grobheight", textGrob(t[longest_coln], rot = angle_col, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  } else {
    coln_height <- unit(5, "bigpts")
  }

  if (!is.null(rown[1])) {
    t <- rown
    tw <- strwidth(t, units = "in", cex = fontsize_row / fontsize)
    if (annotation_names_col) {
      t <- c(t, colnames(annotation_col))
      tw <- c(tw, strwidth(colnames(annotation_col), units = "in"))
    }
    longest_rown <- which.max(tw)
    gp <- list(fontsize = ifelse(longest_rown <= length(rown), fontsize_row, fontsize), ...)
    rown_width <- unit(1, "grobwidth", textGrob(t[longest_rown], rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  } else {
    rown_width <- unit(5, "bigpts")
  }

  gp <- list(fontsize = fontsize, ...)
  # Legend position
  if (!is.na2(legend)) {
    longest_break <- which.max(nchar(names(legend)))
    longest_break <- unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
    title_length <- unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
    legend_width <- unit(12, "bigpts") + longest_break * 1.2
    legend_width <- max(title_length, legend_width)
  } else {
    legend_width <- unit(0, "bigpts")
  }

  # Set main title height
  if (is.na(main)) {
    main_height <- unit(0, "npc")
  } else {
    main_height <- unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }

  # Column annotations
  textheight <- unit(fontsize, "bigpts")

  if (!is.na2(annotation_col)) {
    # Column annotation height
    annot_col_height <- ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

    # Width of the correponding legend
    t <- c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width <- unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_col_legend_width <- unit(0, "npc")
    }
  } else {
    annot_col_height <- unit(0, "bigpts")
    annot_col_legend_width <- unit(0, "bigpts")
  }

  # Row annotations
  if (!is.na2(annotation_row)) {
    # Row annotation width
    annot_row_width <- ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")

    # Width of the correponding legend
    t <- c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width <- unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_row_legend_width <- unit(0, "npc")
    }
  } else {
    annot_row_width <- unit(0, "bigpts")
    annot_row_legend_width <- unit(0, "bigpts")
  }

  annot_legend_width <- max(annot_row_legend_width, annot_col_legend_width)

  # Tree height
  treeheight_col <- unit(treeheight_col, "bigpts") + unit(5, "bigpts")
  treeheight_row <- unit(treeheight_row, "bigpts") + unit(5, "bigpts")

  # Set cell sizes
  if (is.na(cellwidth)) {
    mat_width <- unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width
  } else {
    mat_width <- unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(4, "bigpts")
  }

  if (is.na(cellheight)) {
    mat_height <- unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
  } else {
    mat_height <- unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(4, "bigpts")
  }

  # Produce gtable
  if (annotation_row_side == "left") {
    # ORIGINAL order:
    # 1: row dendro | 2: row annotation | 3: matrix | 4: row names | 5: legend | 6: annot legend
    gt <- gtable(
      widths = unit.c(treeheight_row, annot_row_width, mat_width, rown_width, legend_width, annot_legend_width),
      heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height),
      vp = viewport(gp = do.call(gpar, gp))
    )
  } else {
    # RIGHT side order:
    # 1: row dendro | 2: matrix | 3: row annotation | 4: row names | 5: legend | 6: annot legend
    gt <- gtable(
      widths = unit.c(treeheight_row, mat_width, annot_row_width, rown_width, legend_width, annot_legend_width),
      heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height),
      vp = viewport(gp = do.call(gpar, gp))
    )
  }


  cw <- convertWidth(mat_width - (length(gaps_col) * unit(4, "bigpts")), "bigpts", valueOnly = T) / ncol
  ch <- convertHeight(mat_height - (length(gaps_row) * unit(4, "bigpts")), "bigpts", valueOnly = T) / nrow

  # Return minimal cell dimension in bigpts to decide if borders are drawn
  mindim <- min(cw, ch)

  res <- list(gt = gt, mindim = mindim)

  return(res)
}

find_coordinates <- function(n, gaps, m = 1:n) {
  if (length(gaps) == 0) {
    return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc")))
  }

  if (max(gaps) > n) {
    stop("Gaps do not match with matrix size")
  }

  size <- (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))

  gaps2 <- apply(sapply(gaps, function(gap, x) {
    x > gap
  }, m), 1, sum)
  coord <- m * size + (gaps2 * unit("4", "bigpts"))

  return(list(coord = coord, size = size))
}

draw_dendrogram <- function(hc, gaps, horizontal = T) {
  h <- hc$height / max(hc$height) / 1.05
  m <- hc$merge
  o <- hc$order
  n <- length(o)

  m[m > 0] <- n + m[m > 0]
  m[m < 0] <- abs(m[m < 0])

  dist <- matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y")))
  dist[1:n, 1] <- 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)

  for (i in 1:nrow(m)) {
    dist[n + i, 1] <- (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
    dist[n + i, 2] <- h[i]
  }

  draw_connection <- function(x1, x2, y1, y2, y) {
    res <- list(
      x = c(x1, x1, x2, x2),
      y = c(y1, y, y, y2)
    )

    return(res)
  }

  x <- rep(NA, nrow(m) * 4)
  y <- rep(NA, nrow(m) * 4)
  id <- rep(1:nrow(m), rep(4, nrow(m)))

  for (i in 1:nrow(m)) {
    c <- draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    k <- (i - 1) * 4 + 1
    x[k:(k + 3)] <- c$x
    y[k:(k + 3)] <- c$y
  }

  x <- find_coordinates(n, gaps, x * n)$coord
  y <- unit(y, "npc")

  if (!horizontal) {
    a <- x
    x <- unit(1, "npc") - y
    y <- unit(1, "npc") - a
  }
  res <- polylineGrob(x = x, y = y, id = id)

  return(res)
}

draw_matrix <- function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color) {
  n <- nrow(matrix)
  m <- ncol(matrix)

  coord_x <- find_coordinates(m, gaps_cols)
  coord_y <- find_coordinates(n, gaps_rows)

  x <- coord_x$coord - 0.5 * coord_x$size
  y <- unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)

  coord <- expand.grid(y = y, x = x)

  res <- gList()

  res[["rect"]] <- rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))

  if (attr(fmat, "draw")) {
    res[["text"]] <- textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
  }

  res <- gTree(children = res)

  return(res)
}

draw_colnames <- function(coln, gaps, vjust_col, hjust_col, angle_col, ...) {
  coord <- find_coordinates(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size

  res <- textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = vjust_col, hjust = hjust_col, rot = angle_col, gp = gpar(...))

  return(res)
}

draw_rownames <- function(rown, gaps, ...) {
  coord <- find_coordinates(length(rown), gaps)
  y <- unit(1, "npc") - (coord$coord - 0.5 * coord$size)

  res <- textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))

  return(res)
}

draw_legend <- function(color, breaks, legend, ...) {
  color <- color[!is.infinite(breaks)]
  breaks <- breaks[!is.infinite(breaks)]

  height <- min(unit(1, "npc"), unit(150, "bigpts"))

  legend_pos <- (legend - min(breaks)) / (max(breaks) - min(breaks))
  legend_pos <- height * legend_pos + (unit(1, "npc") - height)

  breaks <- (breaks - min(breaks)) / (max(breaks) - min(breaks))
  breaks <- height * breaks + (unit(1, "npc") - height)

  h <- breaks[-1] - breaks[-length(breaks)]

  rect <- rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
  text <- textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))

  res <- grobTree(rect, text)

  return(res)
}

convert_annotations <- function(annotation, annotation_colors) {
  new <- annotation
  for (i in 1:ncol(annotation)) {
    a <- annotation[, i]
    b <- annotation_colors[[colnames(annotation)[i]]]
    if (is.character(a) | is.factor(a)) {
      a <- as.character(a)

      if (length(setdiff(setdiff(a, NA), names(b))) > 0) {
        stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
      }
      new[, i] <- b[a]
    } else {
      a <- cut(a, breaks = 100)
      new[, i] <- colorRampPalette(b)(100)[a]
    }
  }
  return(as.matrix(new))
}

draw_annotations <- function(converted_annotations, border_color, gaps, fontsize, horizontal) {
  n <- ncol(converted_annotations)
  m <- nrow(converted_annotations)

  coord_x <- find_coordinates(m, gaps)

  x <- coord_x$coord - 0.5 * coord_x$size

  # y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
  y <- cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
  y <- unit(y, "bigpts")

  if (horizontal) {
    coord <- expand.grid(x = x, y = y)
    res <- rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = unit(fontsize, "bigpts"), gp = gpar(fill = converted_annotations, col = border_color))
  } else {
    a <- x
    x <- unit(1, "npc") - y
    y <- unit(1, "npc") - a

    coord <- expand.grid(y = y, x = x)
    res <- rectGrob(x = coord$x, y = coord$y, width = unit(fontsize, "bigpts"), height = coord_x$size, gp = gpar(fill = converted_annotations, col = border_color))
  }

  return(res)
}

draw_annotation_names <- function(annotations, fontsize, horizontal, hjust_col, vjust_col, angle_col) {
  n <- ncol(annotations)

  x <- unit(3, "bigpts")

  y <- cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1
  y <- unit(y, "bigpts")

  if (horizontal) {
    res <- textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
  } else {
    a <- x
    x <- unit(1, "npc") - y
    y <- unit(1, "npc") - a

    res <- textGrob(colnames(annotations), x = x, y = y, vjust = vjust_col, hjust = hjust_col, rot = angle_col, gp = gpar(fontsize = fontsize, fontface = 2))
  }

  return(res)
}

draw_annotation_legend <- function(annotation, annotation_colors, border_color, ...) {
  y <- unit(1, "npc")
  text_height <- unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))

  res <- gList()
  for (i in names(annotation)) {
    res[[i]] <- textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))

    y <- y - 1.5 * text_height
    if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
      n <- length(annotation_colors[[i]])
      yy <- y - (1:n - 1) * 2 * text_height

      res[[paste(i, "r")]] <- rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]]))
      res[[paste(i, "t")]] <- textGrob(names(annotation_colors[[i]]), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(...))

      y <- y - n * 2 * text_height
    } else {
      yy <- y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
      h <- 8 * text_height * 0.25

      res[[paste(i, "r")]] <- rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
      res[[paste(i, "r2")]] <- rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = NA))

      txt <- rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
      yy <- y - c(1, 7) * text_height
      res[[paste(i, "t")]] <- textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
      y <- y - 8 * text_height
    }
    y <- y - 1.5 * text_height
  }

  res <- gTree(children = res)

  return(res)
}

draw_main <- function(text, ...) {
  res <- textGrob(text, gp = gpar(fontface = "bold", ...))

  return(res)
}

vplayout <- function(x, y) {
  return(viewport(layout.pos.row = x, layout.pos.col = y))
}

heatmap_motor <- function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, annotation_names_row, annotation_names_col, main, fontsize, fontsize_row, fontsize_col, hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, gaps_col, gaps_row, labels_row, labels_col, annotation_row_side = c("left","right"),...) {

  annotation_row_side <- match.arg(annotation_row_side)

  # Set layout
  # lo <- lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, gaps_col = gaps_col, ...)

  # res <- lo$gt
  # mindim <- lo$mindim

  # Set layout
  lo_ret <- lo(coln = labels_col, rown = labels_row,
               nrow = nrow(matrix), ncol = ncol(matrix),
               cellwidth = cellwidth, cellheight = cellheight,
               treeheight_col = treeheight_col, treeheight_row = treeheight_row,
               legend = legend,
               annotation_col = annotation_col, annotation_row = annotation_row,
               annotation_colors = annotation_colors, annotation_legend = annotation_legend,
               annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col,
               main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
               angle_col = angle_col, gaps_row = gaps_row, gaps_col = gaps_col,
               annotation_row_side = annotation_row_side, ...)

  res <- lo_ret$gt
  mindim <- lo_ret$mindim

  # Column indices depending on side
  if (annotation_row_side == "left") {
    col_matrix <- 3
    col_row_annot <- 2
    col_row_names <- 4
    col_title <- 3
    col_col_tree <- 3
    col_colnames <- 3
    col_colannot <- 3
    col_colannot_names <- 4
  } else { # right
    col_matrix <- 2
    col_row_annot <- 3
    col_row_names <- 4
    col_title <- 2
    col_col_tree <- 2
    col_colnames <- 2
    col_colannot <- 2
    col_colannot_names <- 4
  }

  if (!is.na(filename)) {
    if (is.na(height)) {
      height <- convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if (is.na(width)) {
      width <- convertWidth(gtable_width(res), "inches", valueOnly = T)
    }

    # Get file type
    r <- regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) stop("Improper filename")
    ending <- substr(filename, r + 1, r + attr(r, "match.length"))

    f <- switch(ending,
      pdf = function(x, ...) pdf(x, ...),
      png = function(x, ...) png(x, units = "in", res = 300, ...),
      jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
      jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
      tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
      bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
      stop("File type should be: pdf, png, bmp, jpg, tiff")
    )

    # print(sprintf("height:%f width:%f", height, width))

    # gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)

    f(filename, height = height, width = width)
    # gt <- heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, fontsize_number = fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
    gt <- heatmap_motor(matrix,
                        cellwidth = cellwidth, cellheight = cellheight, border_color = border_color,
                        tree_col = tree_col, tree_row = tree_row,
                        treeheight_col = treeheight_col, treeheight_row = treeheight_row,
                        breaks = breaks, color = color, legend = legend,
                        annotation_row = annotation_row, annotation_col = annotation_col,
                        annotation_colors = annotation_colors, annotation_legend = annotation_legend,
                        annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col,
                        filename = NA, main = main,
                        fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                        hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col,
                        fmat = fmat, fontsize_number = fontsize_number, number_color = number_color,
                        labels_row = labels_row, labels_col = labels_col,
                        gaps_col = gaps_col, gaps_row = gaps_row,
                        annotation_row_side = annotation_row_side, ...)
    grid.draw(gt)
    dev.off()

    return(gt)
  }

  # Omit border color if cell size is too small
  if (mindim < 3) border_color <- NA

  # Draw title
  if (!is.na(main)) {
    elem <- draw_main(main, fontsize = 1.3 * fontsize, ...)
    res <- gtable_add_grob(res, elem, t = 1, l = col_title, name = "main", clip = "off")
  }

  # Draw tree for the columns
  if (!is.na2(tree_col) & treeheight_col != 0) {
    elem <- draw_dendrogram(tree_col, gaps_col, horizontal = TRUE)
    res <- gtable_add_grob(res, elem, t = 2, l = col_col_tree, name = "col_tree")
  }

  # Draw tree for the rows
  if (!is.na2(tree_row) & treeheight_row != 0) {
    elem <- draw_dendrogram(tree_row, gaps_row, horizontal = FALSE)
    res <- gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }

  # Draw matrix
  elem <- draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
  res <- gtable_add_grob(res, elem, t = 4, l = col_matrix, clip = "off", name = "matrix")

  # Draw colnames
  if (length(labels_col) != 0) {
    pars <- list(labels_col, gaps = gaps_col, fontsize = fontsize_col,
                 hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, ...)
    elem <- do.call(draw_colnames, pars)
    res <- gtable_add_grob(res, elem, t = 5, l = col_colnames, clip = "off", name = "col_names")
  }

  # Draw rownames
  if (length(labels_row) != 0) {
    pars <- list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
    elem <- do.call(draw_rownames, pars)
    res <- gtable_add_grob(res, elem, t = 4, l = col_row_names, clip = "off", name = "row_names")
  }

  # Draw annotation tracks on cols
  if (!is.na2(annotation_col)) {
    converted_annotation <- convert_annotations(annotation_col, annotation_colors)
    elem <- draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = TRUE)
    res <- gtable_add_grob(res, elem, t = 3, l = col_colannot, clip = "off", name = "col_annotation")

    if (annotation_names_col) {
      elem <- draw_annotation_names(annotation_col, fontsize, horizontal = TRUE)
      res <- gtable_add_grob(res, elem, t = 3, l = col_colannot_names, clip = "off",
                             name = "col_annotation_names")
    }
  }

  # Draw annotation tracks on rows
  if (!is.na2(annotation_row)) {
    converted_annotation <- convert_annotations(annotation_row, annotation_colors)
    elem <- draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = FALSE)
    res <- gtable_add_grob(res, elem, t = 4, l = col_row_annot, clip = "off", name = "row_annotation")

    if (annotation_names_row) {
      elem <- draw_annotation_names(annotation_row, fontsize, horizontal = FALSE,
                                    hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col)
      res <- gtable_add_grob(res, elem, t = 5, l = col_row_annot + 1, clip = "off",
                             name = "row_annotation_names")
    }
  }

  # Draw annotation legend
  annotation <- c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
  annotation <- annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]

  if (length(annotation) > 0 & annotation_legend) {
    elem <- draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
    t <- ifelse(is.null(labels_row), 4, 3)
    res <- gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
  }

  # Draw legend
  if (!is.na2(legend)) {
    elem <- draw_legend(color, breaks, legend, fontsize = fontsize, ...)
    t <- ifelse(is.null(labels_row), 4, 3)
    res <- gtable_add_grob(res, elem, t = t, l = 5, b = 5, clip = "off", name = "legend")
  }

  return(res)
}

generate_breaks <- function(x, n, center = F) {
  if (center) {
    m <- max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res <- seq(-m, m, length.out = n + 1)
  } else {
    res <- seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }

  return(res)
}

scale_vec_colours <- function(x, col = rainbow(10), breaks = NA, na_col) {
  res <- col[as.numeric(cut(x, breaks = breaks, include.lowest = T))]
  res[is.na(res)] <- na_col
  return(res)
}

scale_colours <- function(mat, col = rainbow(10), breaks = NA, na_col) {
  mat <- as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks, na_col = na_col), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cluster_mat <- function(mat, distance, method) {
  if (!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))) {
    stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if (!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & !(inherits(distance, "dist"))) {
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  if (distance[1] == "correlation") {
    d <- as.dist(1 - cor(t(mat)))
  } else {
    if (inherits(distance, "dist")) {
      d <- distance
    } else {
      d <- dist(mat, method = distance)
    }
  }

  return(hclust(d, method = method))
}

scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = T)
  s <- apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat <- function(mat, scale) {
  if (!(scale %in% c("none", "row", "column"))) {
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat <- switch(scale,
    none = mat,
    row = scale_rows(mat),
    column = t(scale_rows(t(mat)))
  )
  return(mat)
}

generate_annotation_colours <- function(annotation, annotation_colors, drop) {
  if (is.na2(annotation_colors)) {
    annotation_colors <- list()
  }
  count <- 0
  for (i in 1:length(annotation)) {
    annotation[[i]] <- annotation[[i]][!is.na(annotation[[i]])]
    if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
      if (is.factor(annotation[[i]]) & !drop) {
        count <- count + length(levels(annotation[[i]]))
      } else {
        count <- count + length(unique(annotation[[i]]))
      }
    }
  }

  factor_colors <- dscale(factor(1:count), hue_pal(l = 75))

  oldseed <- NULL
  if (exists(".Random.seed")) {
    oldseed <- get(".Random.seed", pos = .GlobalEnv)
  }

  set.seed(3453)

  cont_counter <- 2
  for (i in 1:length(annotation)) {
    if (!(names(annotation)[i] %in% names(annotation_colors))) {
      if (is.character(annotation[[i]]) | is.factor(annotation[[i]])) {
        n <- length(unique(annotation[[i]]))
        if (is.factor(annotation[[i]]) & !drop) {
          n <- length(levels(annotation[[i]]))
        }
        ind <- sample(1:length(factor_colors), n)
        annotation_colors[[names(annotation)[i]]] <- factor_colors[ind]
        l <- levels(as.factor(annotation[[i]]))
        l <- l[l %in% unique(annotation[[i]])]
        if (is.factor(annotation[[i]]) & !drop) {
          l <- levels(annotation[[i]])
        }

        names(annotation_colors[[names(annotation)[i]]]) <- l
        factor_colors <- factor_colors[-ind]
      } else {
        annotation_colors[[names(annotation)[i]]] <- brewer_pal("seq", cont_counter)(5)[1:4]
        cont_counter <- cont_counter + 1
      }
    }
  }

  if (!is.null(oldseed)) {
    assign(".Random.seed", oldseed, pos = .GlobalEnv)
  } else {
    remove(.Random.seed, pos = .GlobalEnv)
  }

  return(annotation_colors)
}

kmeans_ggheat <- function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...) {
  # Filter data
  if (!is.na(sd_limit)) {
    s <- apply(mat, 1, sd)
    mat <- mat[s > sd_limit, ]
  }

  # Cluster data
  set.seed(1245678)
  km <- kmeans(mat, k, iter.max = 100)
  mat2 <- km$centers

  # Compose rownames
  t <- table(km$cluster)
  rownames(mat2) <- sprintf("cl%s_size_%d", names(t), t)

  # Draw heatmap
  ggheat(mat2, ...)
}

find_gaps <- function(tree, cutree_n) {
  v <- cutree(tree, cutree_n)[tree$order]
  gaps <- which((v[-1] - v[-length(v)]) != 0)
}

is.na2 <- function(x) {
  if (is.list(x) | length(x) > 1) {
    return(FALSE)
  }
  if (length(x) == 0) {
    return(TRUE)
  }

  return(is.na(x))
}

identity2 <- function(x, ...) {
  return(x)
}

#' A function to draw clustered heatmaps.
#'
#' A function to draw clustered heatmaps where one has better control over some graphical
#' parameters such as cell size, etc.
#'
#' The function also allows to aggregate the rows using kmeans clustering. This is
#' advisable if number of rows is so big that R cannot handle their hierarchical
#' clustering anymore, roughly more than 1000. Instead of showing all the rows
#' separately one can cluster the rows in advance and show only the cluster centers.
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to aggregate the
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one
#' element longer than color vector. Used for mapping values to colors. Useful, if needed
#' to map certain values to certain colors, to certain values. If value is NA then the
#' breaks are calculated automatically. When breaks do not cover the range of values,
#' then any value larger than \code{max(breaks)} will have the largest color and any value
#' lower than\code{ min(breaks)} will get the lowest color.
#' @param border_color color of cell borders on heatmap; default = NA (no borders). Set to, e.g., `"grey80"` if you want colored borders.
#' @param cellwidth individual cell width in points. If left as NA, then the values
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA,
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in
#' either the row direction or the column direction, or none. Corresponding values are
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object,
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible
#' values are \code{"correlation"} for Pearson correlation and all the distances
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Defaults to `"ward.D2"`, and accepts the same values as \code{\link{hclust}}.
#' @param clustering_callback callback function to modify the clustering. Is
#' called with two parameters: original \code{hclust} object and the matrix
#' used for clustering. Must return a \code{hclust} object.
#' @param cutree_rows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the
#' argument is ignored
#' @param cutree_cols similar to \code{cutree_rows}, but for columns
#' @param treeheight_row the height of a tree for rows, if these are clustered.
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered.
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation_row data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The
#' rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotation_row_side character, "left" (default) or "right" ?where to draw row annotations.
#' @param annotation_col similar to annotation_row, but for columns.
#' @param annotation deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors list for specifying annotation_row and
#' annotation_col track colors manually. It is  possible to define the colors
#' for only some of the features. Check examples for  details.
#' @param annotation_legend boolean value showing if the legend for annotation
#' tracks should be drawn.
#' @param annotation_names_row boolean value showing if the names for row annotation
#' tracks should be drawn.
#' @param annotation_names_col boolean value showing if the names for column annotation
#' tracks should be drawn.
#' @param drop_levels logical to determine if unused levels are also shown in
#' the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot
#' @param fontsize_row fontsize for rownames (Default: fontsize)
#' @param fontsize_col fontsize for colnames (Default: fontsize)
#' @param angle_col angle of the column labels, right now one can choose only from few
#' predefined options (0, 45, 90, 270 and 315)
#' @param display_numbers logical determining if the numeric values are also printed to
#' the cells. If this is a matrix (with same dimensions as original matrix), the contents
#' of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in cells.
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential
#' notation (see more in \code{\link{sprintf}}).
#' @param number_color color of the text
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param gaps_row vector of row indices that show where to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutree_row}
#'  to see how to introduce gaps to clustered rows.
#' @param gaps_col similar to gaps_row, but for columns.
#' @param labels_row custom labels for rows that are used instead of rownames.
#' @param labels_col similar to labels_row, but for columns.
#' @param filename file path where to save the picture. Filetype is decided by
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param na_col specify the color of the NA cell in the matrix.
#' @param \dots graphical parameters for the text used in plot. Parameters passed to
#' \code{\link{grid.text}}, see \code{\link{gpar}}.
#'
#' @return
#' Invisibly a \code{ggheat} object that is a list with components
#' \itemize{
#'     \item \code{tree_row} the clustering of rows as \code{\link{hclust}} object
#'     \item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was
#' specified
#'     \item \code{gtable} a \code{\link[gtable]{gtable}} object containing the heatmap,
#'     can be used for combining the heatmap with other plots
#' }
#'
#' @author  Yaoxiang Li <liyaoxiang@@outlook.com>
#' @examples
#' # Create test matrix
#' test <- matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] <- test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] <- test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] <- test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) <- paste("Test", 1:10, sep = "")
#' rownames(test) <- paste("Gene", 1:20, sep = "")
#'
#' # Draw heatmaps
#' ggheat(test)
#' ggheat(test, kmeans_k = 2)
#' ggheat(test, scale = "row", clustering_distance_rows = "correlation")
#' ggheat(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#' ggheat(test, cluster_row = FALSE)
#' ggheat(test, legend = FALSE)
#'
#' # Show text within cells
#' ggheat(test, display_numbers = TRUE)
#' ggheat(test, display_numbers = TRUE, number_format = "%.1e")
#' ggheat(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
#' ggheat(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c(
#'   "0",
#'   "1e-4", "1e-3", "1e-2", "1e-1", "1"
#' ))
#'
#' # Fix cell sizes and save to file with correct size
#' ggheat(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#' # ggheat(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
#'
#' # Generate annotations for rows and columns
#' annotation_col <- data.frame(
#'   CellType = factor(rep(c("CT1", "CT2"), 5)),
#'   Time = 1:5
#' )
#' rownames(annotation_col) <- paste("Test", 1:10, sep = "")
#'
#' annotation_row <- data.frame(
#'   GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
#' )
#' rownames(annotation_row) <- paste("Gene", 1:20, sep = "")
#'
#' # Display row and color annotations
#' ggheat(test, annotation_col = annotation_col)
#' ggheat(test, annotation_col = annotation_col, annotation_legend = FALSE)
#' ggheat(test, annotation_col = annotation_col, annotation_row = annotation_row)
#'
#' # Change angle of text in the columns
#' ggheat(test, annotation_col = annotation_col, annotation_row = annotation_row, angle_col = "45")
#' ggheat(test, annotation_col = annotation_col, angle_col = "0")
#'
#' # Specify colors
#' ann_colors <- list(
#'   Time = c("white", "firebrick"),
#'   CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#'   GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
#' )
#'
#' ggheat(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
#' ggheat(test,
#'   annotation_col = annotation_col, annotation_row = annotation_row,
#'   annotation_colors = ann_colors
#' )
#' ggheat(test, annotation_col = annotation_col, annotation_colors = ann_colors[2])
#'
#' # Gaps in heatmaps
#' ggheat(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
#' ggheat(test,
#'   annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14),
#'   cutree_col = 2
#' )
#'
#' # Show custom strings as row/col names
#' labels_row <- c(
#'   "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
#'   "", "", "Il10", "Il15", "Il1b"
#' )
#'
#' ggheat(test, annotation_col = annotation_col, labels_row = labels_row)
#'
#' # Specifying clustering from distance matrix
#' drows <- dist(test, method = "minkowski")
#' dcols <- dist(t(test), method = "minkowski")
#' ggheat(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
#'
#' # Modify ordering of the clusters using clustering callback option
#' callback <- function(hc, mat) {
#'   sv <- svd(t(mat))$v[, 1]
#'   dend <- reorder(as.dendrogram(hc), wts = sv)
#'   as.hclust(dend)
#' }
#'
#' ggheat(test, clustering_callback = callback)
#'
#' \dontrun{
#' # Same using dendsort package
#' library(dendsort)
#'
#' callback <- function(hc, ...) {
#'   dendsort(hc)
#' }
#' ggheat(test, clustering_callback = callback)
#' }
#'
#' @export
ggheat <-
  function(mat,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           kmeans_k = NA,
           breaks = NA,
           border_color = NA,
           cellwidth = NA,
           cellheight = NA,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           clustering_callback = identity2,
           cutree_rows = NA,
           cutree_cols = NA,
           treeheight_row = ifelse((inherits(cluster_rows, "hclust")) ||
                                     cluster_rows, 50, 0),
           treeheight_col = ifelse((inherits(cluster_cols, "hclust")) ||
                                     cluster_cols, 50, 0),
           legend = TRUE,
           legend_breaks = NA,
           legend_labels = NA,
           annotation_row = NA,
           annotation_row_side = c("left","right"),
           annotation_col = NA,
           annotation = NA,
           annotation_colors = NA,
           annotation_legend = TRUE,
           annotation_names_row = TRUE,
           annotation_names_col = TRUE,
           drop_levels = TRUE,
           show_rownames = T,
           show_colnames = T,
           main = NA,
           fontsize = 10,
           fontsize_row = fontsize,
           fontsize_col = fontsize,
           angle_col = c("270", "0", "45", "90", "315"),
           display_numbers = F,
           number_format = "%.2f",
           number_color = "grey30",
           fontsize_number = 0.8 * fontsize,
           gaps_row = NULL,
           gaps_col = NULL,
           labels_row = NULL,
           labels_col = NULL,
           filename = NA,
           width = NA,
           height = NA,
           silent = FALSE,
           na_col = "#DDDDDD",
           ...) {


    annotation_row_side <- match.arg(annotation_row_side)

    # Set labels
  if (is.null(labels_row)) {
    labels_row <- rownames(mat)
  }
  if (is.null(labels_col)) {
    labels_col <- colnames(mat)
  }

  # Preprocess matrix
  mat <- as.matrix(mat)
  orig_mat <- mat    # Updated: capture the unscaled data


  if (scale != "none") {
    mat <- scale_mat(mat, scale)
    if (is.na2(breaks)) {
      breaks <- generate_breaks(mat, length(color), center = T)
    }
  }


  # Kmeans
  if (!is.na(kmeans_k)) {
    # Cluster data
    km <- kmeans(mat, kmeans_k, iter.max = 100)
    mat <- km$centers

    # Compose rownames
    t <- table(km$cluster)
    labels_row <- sprintf("Cluster: %s Size: %d", names(t), t)
  } else {
    km <- NA
  }

  # Format numbers to be displayed in cells
  if (is.matrix(display_numbers) | is.data.frame(display_numbers)) {
    if (nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)) {
      stop("If display_numbers provided as matrix, its dimensions have to match with mat")
    }

    display_numbers <- as.matrix(display_numbers)
    fmat <- matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
    fmat_draw <- TRUE
  } else {
    if (display_numbers) {
      fmat <- matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw <- TRUE
    } else {
      fmat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
      fmat_draw <- FALSE
    }
  }

  # Do clustering
  if ((inherits(cluster_rows, "hclust")) || cluster_rows) {
    if (inherits(cluster_rows, "hclust")) {
      tree_row <- cluster_rows
    } else {
      tree_row <- cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
      tree_row <- clustering_callback(tree_row, mat)
    }
    mat <- mat[tree_row$order, , drop = FALSE]
    fmat <- fmat[tree_row$order, , drop = FALSE]
    labels_row <- labels_row[tree_row$order]
    if (!is.na(cutree_rows)) {
      gaps_row <- find_gaps(tree_row, cutree_rows)

      # Updated: automatic annotate row clusters when `cutree_rows` provided
      ## get the actual cluster IDs
      row_clusters <- cutree(tree_row, cutree_rows)
      ## build a default annotation df of those cluster IDs
      default_row_annot <- data.frame(cluster = factor(row_clusters))
      rownames(default_row_annot) <- names(row_clusters)
      ## merge into any userprovided annotation_row
      if (is.na2(annotation_row)) {
        annotation_row <- default_row_annot
      } else {
        annotation_row <- cbind(
          annotation_row[rownames(default_row_annot), , drop = FALSE],
          default_row_annot
        )
      }
      # End of Updated

    } else {
      gaps_row <- NULL
    }
  } else {
    tree_row <- NA
    treeheight_row <- 0
  }

  if ((inherits(cluster_cols, "hclust")) || cluster_cols) {
    if (inherits(cluster_cols, "hclust")) {
      tree_col <- cluster_cols
    } else {
      tree_col <- cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
      tree_col <- clustering_callback(tree_col, t(mat))
    }
    mat <- mat[, tree_col$order, drop = FALSE]
    fmat <- fmat[, tree_col$order, drop = FALSE]
    labels_col <- labels_col[tree_col$order]
    if (!is.na(cutree_cols)) {
      gaps_col <- find_gaps(tree_col, cutree_cols)
    } else {
      gaps_col <- NULL
    }
  } else {
    tree_col <- NA
    treeheight_col <- 0
  }

  attr(fmat, "draw") <- fmat_draw

  # Colors and scales
  if (!is.na2(legend_breaks) & !is.na2(legend_labels)) {
    if (length(legend_breaks) != length(legend_labels)) {
      stop("Lengths of legend_breaks and legend_labels must be the same")
    }
  }


  if (is.na2(breaks)) {
    breaks <- generate_breaks(as.vector(mat), length(color))
  }

  if (!is.infinite(min(breaks))) {
    breaks <- c(-Inf, breaks)
    color <- c(color[1], color)
  }
  if (!is.infinite(max(breaks))) {
    breaks <- c(breaks, Inf)
    color <- c(color, color[length(color)])
  }
  non_inf_breaks <- breaks[!is.infinite(breaks)]

  if (legend & is.na2(legend_breaks)) {
    legend <- grid.pretty(range(as.vector(non_inf_breaks)))
    names(legend) <- legend
  } else if (legend & !is.na2(legend_breaks)) {
    legend <- legend_breaks[legend_breaks >= min(non_inf_breaks) & legend_breaks <= max(non_inf_breaks)]

    if (!is.na2(legend_labels)) {
      legend_labels <- legend_labels[legend_breaks >= min(non_inf_breaks) & legend_breaks <= max(non_inf_breaks)]
      names(legend) <- legend_labels
    } else {
      names(legend) <- legend
    }
  } else {
    legend <- NA
  }

  mat <- scale_colours(mat, col = color, breaks = breaks, na_col = na_col)

  # Preparing annotations
  if (is.na2(annotation_col) & !is.na2(annotation)) {
    annotation_col <- annotation
  }
  # Select only the ones present in the matrix
  if (!is.na2(annotation_col)) {
    annotation_col <- annotation_col[colnames(mat), , drop = F]
  }

  if (!is.na2(annotation_row)) {
    annotation_row <- annotation_row[rownames(mat), , drop = F]
  }

  annotation <- c(annotation_row, annotation_col)
  annotation <- annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
  if (length(annotation) != 0) {
    annotation_colors <- generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
  } else {
    annotation_colors <- NA
  }

  if (!show_rownames) {
    labels_row <- NULL
  }

  if (!show_colnames) {
    labels_col <- NULL
  }

  # Set colname rotating parameters
  angle_col <- as.character(angle_col)
  angle_col <- match.arg(angle_col)

  if (angle_col == "0") {
    angle_col <- 0
    hjust_col <- 0.5
    vjust_col <- 1
  }
  if (angle_col == "45") {
    angle_col <- 45
    hjust_col <- 1
    vjust_col <- 1
  }
  if (angle_col == "90") {
    angle_col <- 90
    hjust_col <- 1
    vjust_col <- 0.5
  }
  if (angle_col == "270") {
    angle_col <- 270
    hjust_col <- 0
    vjust_col <- 0.5
  }
  if (angle_col == "315") {
    angle_col <- 315
    hjust_col <- 0
    vjust_col <- 1
  }

  # Draw heatmap
  pdf(file = NULL)
  gt <-
    heatmap_motor(
      mat,
      border_color = border_color,
      cellwidth = cellwidth,
      cellheight = cellheight,
      treeheight_col = treeheight_col,
      treeheight_row = treeheight_row,
      tree_col = tree_col,
      tree_row = tree_row,
      filename = filename,
      width = width,
      height = height,
      breaks = breaks,
      color = color,
      legend = legend,
      annotation_row = annotation_row,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      annotation_legend = annotation_legend,
      annotation_names_row = annotation_names_row,
      annotation_names_col = annotation_names_col,
      main = main,
      fontsize = fontsize,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col,
      hjust_col = hjust_col,
      vjust_col = vjust_col,
      angle_col = angle_col,
      fmat = fmat,
      fontsize_number = fontsize_number,
      number_color = number_color,
      gaps_row = gaps_row,
      gaps_col = gaps_col,
      labels_row = labels_row,
      labels_col = labels_col,
      annotation_row_side = annotation_row_side,
      ...
    )

  dev.off()

  if (is.na(filename) & !silent) {
    grid.newpage()
    grid.draw(gt)
  }
  # Updated: build the ordered, un-scaled matrix
  matrix_ordered <- orig_mat[
    if ((inherits(cluster_rows,"hclust"))||cluster_rows) tree_row$order else seq_len(nrow(orig_mat)),
    if ((inherits(cluster_cols,"hclust"))||cluster_cols) tree_col$order else seq_len(ncol(orig_mat)),
    drop = FALSE
  ]
  # Updated: extract the row order and clusters
  if (!is.na(cutree_rows)) {
    ord <- if ((inherits(cluster_rows,"hclust"))||cluster_rows) tree_row$order else seq_len(nrow(orig_mat))
    clusts <- cutree(tree_row, cutree_rows)[rownames(orig_mat)[ord]]
    matrix_ordered <- as.data.frame(matrix_ordered)
    matrix_ordered$cluster_id <- factor(clusts)
  }
  invisible(structure(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt,
                           matrix_ordered = matrix_ordered # Updated: original values, in heatmap order
                           ), class = "ggheat"))
}


#' @method grid.draw ggheat
#' @export
grid.draw.ggheat <- function(x, recording = TRUE) {
  grid.draw(x$gtable)
}

#' @method print ggheat
#' @export
print.ggheat <- function(x, ...) {
  grid.draw(x)
}
