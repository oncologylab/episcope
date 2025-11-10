# BEGIN FILE: utils_grn_link_network_topic.R
# nocov start

#' Validate and normalize a topic-centric GRN links table
#'
#' Many \code{episcope} functions expect a "links" table in long format where each
#' row represents a TF→enhancer→gene candidate relationship (optionally scoped to a
#' topic and a condition) with one or more scores (e.g., link score, footprint score).
#'
#' Required columns:
#' \itemize{
#'   \item \strong{tf}: transcription factor symbol (character)
#'   \item \strong{gene}: gene symbol or ID (character)
#'   \item \strong{enhancer}: enhancer ID (character; e.g., "chr1:100-200")
#'   \item \strong{topic}: topic/group label (character or factor)
#' }
#' Optional columns:
#' \itemize{
#'   \item \strong{condition}: condition / group label
#'   \item numeric score columns (e.g., \strong{link_score}, \strong{fp_score})
#' }
#'
#' @param links data.frame-like
#' @param cols optional name mapping list(tf, gene, enhancer, topic, condition)
#' @param must_have_scores character() of numeric score columns required
#' @return standardized data.frame
#' @export
grn_topic_validate_links <- function(links,
                                     cols = list(),
                                     must_have_scores = character(0)) {
  if (!is.data.frame(links)) cli::cli_abort("`links` must be a data.frame-like object.")

  req_keys <- c("tf","gene","enhancer","topic")
  opt_keys <- c("condition")
  all_keys <- c(req_keys, opt_keys)

  name_map <- stats::setNames(all_keys, all_keys)
  if (length(cols)) {
    nm <- names(cols)
    bad <- setdiff(nm, all_keys)
    if (length(bad)) cli::cli_abort("Unknown column mapping key(s): {bad}")
    for (k in nm) name_map[[k]] <- cols[[k]]
  }

  for (k in req_keys) {
    if (!name_map[[k]] %in% colnames(links)) {
      cli::cli_abort("Missing required column: '{name_map[[k]]}' for key '{k}'.")
    }
  }

  out <- links
  for (k in all_keys) {
    src <- name_map[[k]]
    if (src %in% colnames(out) && src != k) {
      if (k %in% colnames(out)) cli::cli_abort("Cannot rename column to '{k}': name already exists.")
      colnames(out)[match(src, colnames(out))] <- k
    }
  }

  char_cols <- intersect(c("tf","gene","enhancer","topic","condition"), colnames(out))
  for (cc in char_cols) out[[cc]] <- as.character(out[[cc]])

  if (length(must_have_scores)) {
    missing_scores <- setdiff(must_have_scores, colnames(out))
    if (length(missing_scores)) {
      cli::cli_abort("Missing required score column(s): {missing_scores}.")
    }
    for (sc in must_have_scores) {
      out[[sc]] <- suppressWarnings(as.numeric(out[[sc]]))
      if (any(!is.finite(out[[sc]]), na.rm = TRUE)) {
        cli::cli_warn("Non-finite values detected in score column '{sc}'.")
      }
    }
  }
  out
}

#' Build a topic-centric TF→gene network by aggregating enhancer links
#' @param links validated links table
#' @param topic character or NULL
#' @param score_col character; numeric score column to aggregate
#' @param agg_fun one of "sum","mean","max","median"
#' @param min_score numeric filter on aggregated score
#' @param top_n_edges integer or NULL
#' @param min_tf_occ integer; min distinct genes per TF
#' @param min_gene_occ integer; min distinct TFs per gene
#' @param normalize_by_targets logical; divide by contributing enhancers
#' @return list(nodes, edges, meta)
#' @export
grn_topic_build_network <- function(links,
                                    topic = NULL,
                                    score_col = "link_score",
                                    agg_fun = c("sum","mean","max","median"),
                                    min_score = -Inf,
                                    top_n_edges = NULL,
                                    min_tf_occ = 0,
                                    min_gene_occ = 0,
                                    normalize_by_targets = FALSE) {

  if (!is.data.frame(links)) cli::cli_abort("`links` must be a data.frame.")
  need <- c("tf","gene","enhancer","topic")
  if (!all(need %in% colnames(links))) {
    cli::cli_abort("`links` missing columns: {setdiff(need, colnames(links))}. See ?grn_topic_validate_links.")
  }
  if (!score_col %in% colnames(links)) {
    cli::cli_abort("Score column '{score_col}' not found in `links`.")
  }

  agg_fun <- match.arg(agg_fun)

  dat <- links
  if (!is.null(topic)) {
    dat <- dat[dat$topic == topic, , drop = FALSE]
    if (!nrow(dat)) cli::cli_abort("No rows for topic '{topic}'.")
  }

  dat_sc <- dat[, c("tf","gene","enhancer", score_col)]
  sc <- suppressWarnings(as.numeric(dat_sc[[score_col]]))
  sc[!is.finite(sc)] <- NA_real_
  dat_sc[[score_col]] <- sc

  key <- paste0(dat_sc$tf, "||", dat_sc$gene)

  n_enh <- stats::aggregate(dat_sc$enhancer,
                            by = list(key = key),
                            FUN = function(x) length(unique(x)))
  colnames(n_enh)[2] <- "n_enhancers"

  agg_fun_obj <- switch(
    agg_fun,
    sum    = function(x) sum(x, na.rm = TRUE),
    mean   = function(x) mean(x, na.rm = TRUE),
    max    = function(x) max(x, na.rm = TRUE),
    median = function(x) stats::median(x, na.rm = TRUE)
  )
  agg_df <- stats::aggregate(dat_sc[[score_col]],
                             by = list(key = key),
                             FUN = agg_fun_obj)
  colnames(agg_df)[2] <- "score"

  out <- merge(agg_df, n_enh, by = "key", all.x = TRUE, all.y = FALSE)
  sp <- strsplit(out$key, "\\|\\|")
  out$tf   <- vapply(sp, function(s) s[[1]], character(1))
  out$gene <- vapply(sp, function(s) s[[2]], character(1))
  out$key <- NULL

  if (isTRUE(normalize_by_targets)) {
    out$score <- out$score / pmax(1L, out$n_enhancers)
  }

  out <- out[is.finite(out$score) & out$score >= min_score, , drop = FALSE]
  if (!nrow(out)) cli::cli_abort("No edges remain after min_score filtering.")

  if (!is.null(top_n_edges)) {
    ord <- order(out$score, decreasing = TRUE)
    out <- out[ord, , drop = FALSE]
    n_keep <- min(nrow(out), as.integer(top_n_edges))
    out <- out[seq_len(n_keep), , drop = FALSE]
  }

  tf_deg   <- stats::aggregate(out$gene, by = list(tf = out$tf),
                               FUN = function(x) length(unique(x)))
  colnames(tf_deg)[2] <- "deg"
  keep_tf <- tf_deg$tf[tf_deg$deg >= min_tf_occ]

  gene_deg <- stats::aggregate(out$tf, by = list(gene = out$gene),
                               FUN = function(x) length(unique(x)))
  colnames(gene_deg)[2] <- "deg"
  keep_gene <- gene_deg$gene[gene_deg$deg >= min_gene_occ]

  out <- out[out$tf %in% keep_tf & out$gene %in% keep_gene, , drop = FALSE]
  if (!nrow(out)) cli::cli_abort("No edges remain after degree pruning.")

  tf_nodes   <- data.frame(id = unique(out$tf),   label = unique(out$tf),   type = "TF",   stringsAsFactors = FALSE)
  gene_nodes <- data.frame(id = unique(out$gene), label = unique(out$gene), type = "gene", stringsAsFactors = FALSE)
  topic_nodes <- NULL
  if (!is.null(topic)) {
    topic_nodes <- data.frame(id = paste0("topic:", topic), label = topic, type = "topic", stringsAsFactors = FALSE)
  }

  nodes <- rbind(tf_nodes, gene_nodes)
  if (!is.null(topic_nodes)) nodes <- rbind(nodes, topic_nodes)

  edges <- data.frame(
    from = out$tf,
    to   = out$gene,
    score = out$score,
    n_enhancers = out$n_enhancers,
    stringsAsFactors = FALSE
  )

  meta <- list(
    topic = topic,
    score_col = score_col,
    agg_fun = agg_fun,
    min_score = min_score,
    top_n_edges = top_n_edges,
    min_tf_occ = min_tf_occ,
    min_gene_occ = min_gene_occ,
    normalize_by_targets = normalize_by_targets,
    n_edges = nrow(edges),
    n_nodes = nrow(nodes)
  )

  list(nodes = nodes, edges = edges, meta = meta)
}

#' Compute a delta topic network between two conditions
#'
#' For a given topic, aggregate edges within \code{baseline} and \code{contrast}
#' conditions, then compute \code{delta = contrast_score - baseline_score}.
#'
#' @param links links table
#' @param topic character topic
#' @param condition_col column with condition labels
#' @param baseline baseline label
#' @param contrast contrast label
#' @inheritParams grn_topic_build_network
#' @return list(nodes, edges, meta) with score_base, score_contrast, delta
#' @export
grn_topic_delta_network <- function(links,
                                    topic,
                                    condition_col = "condition",
                                    baseline,
                                    contrast,
                                    score_col = "link_score",
                                    agg_fun = c("sum","mean","max","median"),
                                    min_score = -Inf,
                                    normalize_by_targets = FALSE) {
  if (!condition_col %in% colnames(links)) {
    cli::cli_abort("Condition column '{condition_col}' not found in `links`.")
  }
  agg_fun <- match.arg(agg_fun)

  net_base <- grn_topic_build_network(
    links[links[[condition_col]] == baseline, , drop = FALSE],
    topic = topic, score_col = score_col, agg_fun = agg_fun,
    min_score = -Inf, top_n_edges = NULL,
    min_tf_occ = 0, min_gene_occ = 0,
    normalize_by_targets = normalize_by_targets
  )

  net_con  <- grn_topic_build_network(
    links[links[[condition_col]] == contrast, , drop = FALSE],
    topic = topic, score_col = score_col, agg_fun = agg_fun,
    min_score = -Inf, top_n_edges = NULL,
    min_tf_occ = 0, min_gene_occ = 0,
    normalize_by_targets = normalize_by_targets
  )

  key_b <- paste0(net_base$edges$from, "||", net_base$edges$to)
  key_c <- paste0(net_con$edges$from,  "||", net_con$edges$to)

  df_b <- data.frame(key = key_b,
                     from = net_base$edges$from,
                     to   = net_base$edges$to,
                     score_base = net_base$edges$score,
                     stringsAsFactors = FALSE)

  df_c <- data.frame(key = key_c,
                     from = net_con$edges$from,
                     to   = net_con$edges$to,
                     score_contrast = net_con$edges$score,
                     stringsAsFactors = FALSE)

  keys_all <- unique(c(df_b$key, df_c$key))
  res <- data.frame(
    key = keys_all,
    from = NA_character_,
    to   = NA_character_,
    score_base = NA_real_,
    score_contrast = NA_real_,
    stringsAsFactors = FALSE
  )
  if (nrow(df_b)) {
    idx <- match(df_b$key, res$key)
    res$from[idx] <- df_b$from
    res$to[idx] <- df_b$to
    res$score_base[idx] <- df_b$score_base
  }
  if (nrow(df_c)) {
    idx <- match(df_c$key, res$key)
    need_from <- is.na(res$from[idx]) & !is.na(df_c$from)
    res$from[idx][need_from] <- df_c$from[need_from]
    need_to <- is.na(res$to[idx]) & !is.na(df_c$to)
    res$to[idx][need_to] <- df_c$to[need_to]
    res$score_contrast[idx] <- df_c$score_contrast
  }

  res$score_base[!is.finite(res$score_base)] <- 0
  res$score_contrast[!is.finite(res$score_contrast)] <- 0
  res$delta <- res$score_contrast - res$score_base

  keep <- abs(res$delta) >= min_score | res$score_base >= min_score | res$score_contrast >= min_score
  res <- res[keep, , drop = FALSE]

  edges <- res[, c("from","to","score_base","score_contrast","delta")]
  nodes <- rbind(
    data.frame(id = unique(edges$from), label = unique(edges$from), type = "TF", stringsAsFactors = FALSE),
    data.frame(id = unique(edges$to),   label = unique(edges$to),   type = "gene", stringsAsFactors = FALSE),
    data.frame(id = paste0("topic:", topic), label = topic, type = "topic", stringsAsFactors = FALSE)
  )
  nodes <- nodes[!duplicated(nodes$id), , drop = FALSE]

  meta <- list(
    topic = topic,
    score_col = score_col,
    agg_fun = agg_fun,
    condition_col = condition_col,
    baseline = baseline,
    contrast = contrast,
    n_edges = nrow(edges),
    n_nodes = nrow(nodes)
  )

  list(nodes = nodes, edges = edges, meta = meta)
}

#' Render a topic network to interactive HTML (visNetwork)
#' @param net list from grn_topic_build_network or grn_topic_delta_network
#' @param out_html path or NULL
#' @param use_delta logical; if TRUE and delta in edges, use |delta| for width & color by sign
#' @param title optional title
#' @param physics logical
#' @param seed integer for layout reproducibility
#' @return out_html invisibly if written, else the widget
#' @export
grn_topic_render_html <- function(net,
                                  out_html = NULL,
                                  use_delta = FALSE,
                                  title = NULL,
                                  physics = TRUE,
                                  seed = 123) {
  if (!is.list(net) || is.null(net$nodes) || is.null(net$edges)) {
    cli::cli_abort("`net` must be a list with `nodes` and `edges` data.frames.")
  }

  nodes <- net$nodes
  edges <- net$edges

  nodes$shape <- ifelse(nodes$type == "TF", "box",
                        ifelse(nodes$type == "gene", "ellipse", "database"))
  nodes$color.background <- ifelse(nodes$type == "TF", "#e8f0fe",
                                   ifelse(nodes$type == "gene", "#e6f4ea", "#fff4e5"))
  nodes$color.border <- ifelse(nodes$type == "TF", "#1a73e8",
                               ifelse(nodes$type == "gene", "#137333", "#c26401"))
  nodes$font.size <- 16L
  nodes$title <- paste0(nodes$type, ": ", nodes$label)

  width_val <- edges$score
  color <- rep("#555555", nrow(edges))

  if (isTRUE(use_delta) && "delta" %in% colnames(edges)) {
    width_val <- abs(edges$delta)
    color <- ifelse(edges$delta >= 0, "#1a73e8", "#d93025")
  }

  if (length(width_val) && any(is.finite(width_val))) {
    v <- width_val
    v[!is.finite(v)] <- 0
    rng <- range(v, na.rm = TRUE)
    widths <- if (rng[1] == rng[2]) rep(3, length(v)) else 1 + 9 * (v - rng[1]) / (rng[2] - rng[1])
  } else {
    widths <- rep(3, nrow(edges))
  }

  edges_vis <- data.frame(
    from = edges$from,
    to   = edges$to,
    width = widths,
    color = color,
    title = if ("delta" %in% colnames(edges)) {
      paste0("Base=", round(edges$score_base,3),
             "<br>Contrast=", round(edges$score_contrast,3),
             "<br>Delta=", round(edges$delta,3))
    } else {
      paste0("Score=", round(edges$score,3),
             if ("n_enhancers" %in% colnames(edges))
               paste0("<br>n_enhancers=", edges$n_enhancers) else "")
    },
    stringsAsFactors = FALSE
  )

  vis <- visNetwork::visNetwork(nodes, edges_vis, main = title) |>
    visNetwork::visIgraphLayout(randomSeed = as.integer(seed)) |>
    visNetwork::visEdges(smooth = FALSE, arrows = "to") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                           nodesIdSelection = TRUE) |>
    visNetwork::visInteraction(navigationButtons = TRUE, dragNodes = TRUE) |>
    visNetwork::visPhysics(enabled = isTRUE(physics))

  if (is.null(out_html)) {
    return(vis)
  } else {
    htmlwidgets::saveWidget(vis, file = out_html, selfcontained = TRUE)
    return(invisible(out_html))
  }
}

#' Plot a static topic network with ggplot2
#' @param net network list
#' @param out_file optional PDF path
#' @param width,height inches for PDF
#' @param seed reproducible layout
#' @param title optional
#' @param use_delta logical
#' @return ggplot
#' @export
grn_topic_plot_static <- function(net,
                                  out_file = NULL,
                                  width = 7,
                                  height = 5,
                                  seed = 123,
                                  title = NULL,
                                  use_delta = FALSE) {
  if (!is.list(net) || is.null(net$nodes) || is.null(net$edges)) {
    cli::cli_abort("`net` must be a list with `nodes` and `edges` data.frames.")
  }
  nodes <- net$nodes
  edges <- net$edges

  g <- igraph::graph_from_data_frame(
    d = data.frame(from = edges$from, to = edges$to, stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = data.frame(name = nodes$id, type = nodes$type, stringsAsFactors = FALSE)
  )
  set.seed(as.integer(seed))
  lay <- igraph::layout_with_fr(g)
  pos <- as.data.frame(lay)
  colnames(pos) <- c("x","y")
  pos$id <- igraph::V(g)$name

  nodes2 <- merge(nodes, pos, by.x = "id", by.y = "id", all.x = TRUE, all.y = FALSE)
  e2 <- merge(edges, nodes2[, c("id","x","y")], by.x = "from", by.y = "id", all.x = TRUE, all.y = FALSE)
  colnames(e2)[match(c("x","y"), colnames(e2))] <- c("x_from","y_from")
  e2 <- merge(e2, nodes2[, c("id","x","y")], by.x = "to", by_y = "id", all.x = TRUE, all.y = FALSE)
  # fix possible wrong colname
  if (!"y.x" %in% colnames(e2) && !"x_to" %in% colnames(e2)) {
    colnames(e2)[(ncol(e2)-1):ncol(e2)] <- c("x_to","y_to")
  }

  width_val <- if (isTRUE(use_delta) && "delta" %in% colnames(e2)) abs(e2$delta) else e2$score
  width_val <- suppressWarnings(as.numeric(width_val))
  width_val[!is.finite(width_val)] <- 0
  rng <- range(width_val, na.rm = TRUE)
  edge_size <- if (rng[1] == rng[2]) rep(0.4, length(width_val)) else (0.2 + 1.8 * (width_val - rng[1]) / (rng[2] - rng[1]))
  edge_color <- if (isTRUE(use_delta) && "delta" %in% colnames(e2)) ifelse(e2$delta >= 0, "#1a73e8", "#d93025") else "#555555"

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = e2,
      ggplot2::aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      linewidth = edge_size,
      alpha = 0.8,
      lineend = "round",
      colour = edge_color
    ) +
    ggplot2::geom_point(
      data = nodes2,
      ggplot2::aes(x = x, y = y, shape = type),
      size = 4.5,
      stroke = 0.6,
      colour = "#222222",
      fill = scales::alpha("#FFFFFF", 0.9)
    ) +
    ggplot2::geom_text(
      data = nodes2,
      ggplot2::aes(x = x, y = y, label = label),
      size = 3.6,
      vjust = -1.1
    ) +
    ggplot2::scale_shape_manual(values = c(TF = 22, gene = 21, topic = 23)) +
    ggplot2::labs(title = if (is.null(title)) {
      if (!is.null(net$meta$topic)) paste0("Topic network — ", net$meta$topic) else "Topic network"
    } else title,
    x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(out_file)) {
    grDevices::pdf(out_file, width = width, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
    invisible(p)
  } else {
    p
  }
}

#' Convenience: build, render HTML, and (optionally) save a static PDF
#' @inheritParams grn_topic_build_network
#' @inheritParams grn_topic_render_html
#' @inheritParams grn_topic_plot_static
#' @param out_pdf optional PDF path to save a static snapshot
#' @return list(net, html_path, pdf_path)
#' @export
grn_topic_build_and_render <- function(links,
                                       topic = NULL,
                                       score_col = "link_score",
                                       agg_fun = c("sum","mean","max","median"),
                                       min_score = -Inf,
                                       top_n_edges = NULL,
                                       min_tf_occ = 0,
                                       min_gene_occ = 0,
                                       normalize_by_targets = FALSE,
                                       out_html = NULL,
                                       out_pdf = NULL,
                                       title = NULL,
                                       physics = TRUE,
                                       seed = 123) {
  net <- grn_topic_build_network(
    links = links,
    topic = topic,
    score_col = score_col,
    agg_fun = agg_fun,
    min_score = min_score,
    top_n_edges = top_n_edges,
    min_tf_occ = min_tf_occ,
    min_gene_occ = min_gene_occ,
    normalize_by_targets = normalize_by_targets
  )

  html_path <- NULL
  if (!is.null(out_html)) {
    html_path <- grn_topic_render_html(
      net, out_html = out_html, title = title, physics = physics, seed = seed
    )
  }

  pdf_path <- NULL
  if (!is.null(out_pdf)) {
    grn_topic_plot_static(net, out_file = out_pdf, width = 7, height = 5, seed = seed, title = title)
    pdf_path <- out_pdf
  }

  list(net = net, html_path = html_path, pdf_path = pdf_path)
}

# small infix for defaults without importing rlang
`%||%` <- function(a, b) if (!is.null(a)) a else b

# nocov end
# END FILE: utils_grn_link_network_topic.R
