# utils_plot_tf_network_delta.R — Simplified Δ-network plot (all links) w/ crowd controls
# Author: Yaoxiang Li (episcope) — refined by ChatGPT
# Updated: 2025-10-31

#' Plot an interactive Δ-network for all TF–peak–gene links
#' (single-panel Δ = stress − control) with crowd-control options.
#'
#' The plot draws TFs (boxes), genes (dots), and optional peaks (triangles).
#' Edge color encodes Δ magnitude/sign (red = gain/activation; blue = loss/repression);
#' **arrow head styles** encode sign: **arrow** for activation (+), **blunt bar** for repression (−).
#'
#' @param data Tibble/data.frame or path to CSV/TSV with per-row TF–peak–gene entry.
#'   Must contain at least: TF/tf, gene_key/gene, peak_id/peak, and two per-condition
#'   link scores (e.g., \code{link_score_*}). Signs/expressions are optional.
#' @param out_html Optional path to write HTML; if `NULL` returns the htmlwidget.
#' @param layout_algo One of: "fr","kk","lgl","dh","graphopt","circle","grid","sphere","star","random".
#' @param physics Logical; if `TRUE`, enables a repulsive physics preset to relieve crowding.
#' @param add_direct Logical; also draw faint direct TF→gene Δ edges in addition to TF→PEAK→gene.
#' @param edge_filter_min Keep edges with \code{|score_ctrl|} or \code{|score_str|} ≥ this value (magnitude-only).
#' @param min_delta_abs Keep edges with \code{|Δ|} ≥ this value (Δ = str − ctrl).
#' @param keep_top_edges_per_tf Optional integer; per TF keep only the strongest \code{|Δ|} edges.
#' @param peak_mode One of: "show_all" (default), "shared_only" (hide unique-peak triangles),
#'   "hide" (draw TF→gene only; collapses multiple peaks to strongest Δ per TF–gene).
#' @param show_peaks Convenience boolean to show peak triangles (\code{TRUE}, default) or hide them
#'   entirely (\code{FALSE}). If set to \code{FALSE}, this overrides \code{peak_mode} to "hide".
#' @param gene_fc_thresh Threshold for node border up/down coloring (reference set by `de_reference`).
#' @param de_reference "str_over_ctrl" (default) or "ctrl_over_str" for log2FC border coloring.
#' @param motif_db "jaspar2024" (default) or "hocomocov13" to recognize TF symbols.
#' @param width,height Dimensions passed to \code{visNetwork}.
#' @param seed Integer seed for reproducible igraph layouts (ignored when \code{physics=TRUE} and nodes are moved).
#' @param verbose Logical; print progress messages.
#' @param score_ctrl_col,score_str_col Optional explicit column names of per-condition link scores.
#' @param sign_ctrl_col,sign_str_col Optional explicit column names of per-condition link signs ("+"|"-"|0 or numeric).
#' @param tf_expr_ctrl_col,tf_expr_str_col,gene_expr_ctrl_col,gene_expr_str_col Optional
#'   expression columns (used for node tooltips, TF label sizing, and border coloring).
#'
#' @return A \code{visNetwork} htmlwidget (or, invisibly, \code{out_html} path if saved).
#'
#' @export
plot_tf_network_delta <- function(
    data,
    plot_title,
    out_html       = NULL,
    layout_algo    = c("fr","kk","lgl","dh","graphopt","circle","grid","sphere","star","random"),
    add_direct     = TRUE,
    edge_filter_min = 0,
    min_delta_abs   = 0,
    keep_top_edges_per_tf = NULL,
    peak_mode      = c("show_all","shared_only","hide"),
    show_peaks     = TRUE,     # NEW: simple on/off for peak triangles
    gene_fc_thresh = 1.5,
    de_reference   = c("str_over_ctrl","ctrl_over_str"),
    motif_db       = c("jaspar2024","hocomocov13"),
    physics        = FALSE,
    width          = "100%",
    height         = "800px",
    seed           = 1L,
    verbose        = TRUE,
    # explicit mappings
    score_ctrl_col  = NULL,
    score_str_col   = NULL,
    sign_ctrl_col   = NULL,
    sign_str_col    = NULL,
    tf_expr_ctrl_col   = NULL,
    tf_expr_str_col    = NULL,
    gene_expr_ctrl_col = NULL,
    gene_expr_str_col  = NULL
){
  layout_algo <- match.arg(layout_algo)
  de_reference <- match.arg(de_reference)
  motif_db     <- match.arg(motif_db)
  peak_mode    <- match.arg(peak_mode)
  # Boolean override: if show_peaks == FALSE, force peak_mode = "hide"
  if (identical(show_peaks, FALSE)) peak_mode <- "hide"
  set.seed(as.integer(seed) %% .Machine$integer.max)

  # ── Styling ────────────────────────────────────────────────────────────────
  col_edge_pos    <- "#df1d16"  # Δ > 0
  col_edge_neg    <- "#255dae"  # Δ < 0
  col_edge_zero   <- "#C8C8C8"
  col_gene_fill   <- "#6ad961"  #
  col_border_none <- "#BFBFBF"
  col_border_up   <- "#df1d16"
  col_border_down <- "#255dae"
  tf_font_range   <- c(22, 44)
  gene_size_range <- c(1, 5)

  `%||%` <- function(a, b) if (is.null(a)) b else a
  to_rgba <- function(hex, alpha){
    hex <- gsub("^#", "", hex)
    r <- strtoi(substr(hex,1,2),16L); g <- strtoi(substr(hex,3,4),16L); b <- strtoi(substr(hex,5,6),16L)
    sprintf("rgba(%d,%d,%d,%.3f)", r,g,b, max(0, min(1, alpha)))
  }
  robust_z <- function(x){
    x <- as.numeric(x)
    m <- stats::median(x, na.rm = TRUE)
    madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(madv) || madv == 0) {
      sx <- stats::sd(x, na.rm = TRUE); if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
      return((x - m)/(sx + 1e-12))
    }
    (x - m)/(madv + 1e-12)
  }
  clamp01 <- function(x) pmax(0, pmin(1, x))
  parse_sign <- function(x){
    if (is.null(x)) return(NA_integer_)
    if (is.numeric(x)) return(ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L)))
    x <- trimws(as.character(x))
    ifelse(x %in% c("+","plus","pos","positive"), 1L,
           ifelse(x %in% c("-","minus","neg","negative"), -1L,
                  ifelse(x %in% c("0","zero","0.0","."), 0L, NA_integer_)))
  }
  .make_arrows <- function(sign_vec) {
    # + or NA → normal arrow; − → T-bar (blunt)
    kinds <- ifelse(is.na(sign_vec) | sign_vec >= 0L, "arrow", "bar")
    I(lapply(kinds, function(tp) list(to = list(enabled = TRUE, type = tp))))
  }

  # Known TF universe (from installed package files)
  .tf_cache <- new.env(parent = emptyenv())
  get_tf_syms <- function() get0("KNOWN_TF_SYMBOLS", envir = .tf_cache, inherits = FALSE)
  set_tf_syms <- function(x, db) { assign("KNOWN_TF_SYMBOLS", x, envir = .tf_cache); assign("KNOWN_TF_DB", db, envir = .tf_cache); invisible(TRUE) }
  load_tf_universe <- function(db = motif_db){
    cur_db <- get0("KNOWN_TF_DB", envir = .tf_cache, inherits = FALSE)
    cur    <- get_tf_syms()
    if (!is.null(cur) && identical(cur_db, db)) return(invisible(TRUE))
    path <- if (db == "jaspar2024") {
      system.file("extdata","genome","JASPAR2024.txt", package = "episcope")
    } else {
      system.file("extdata","genome","HOCOMOCO_v13_motif_cluster_definition_with_sub_cluster.txt", package = "episcope")
    }
    if (!nzchar(path) || !file.exists(path)) cli::cli_abort("Motif cluster definition not found for motif_db='{db}'.")
    df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    syms <- unique(trimws(unlist(strsplit(df$HGNC, "::", fixed = TRUE), use.names = FALSE)))
    set_tf_syms(syms, db); invisible(TRUE)
  }
  is_known_tf <- function(x){
    if (is.null(get_tf_syms())) load_tf_universe(motif_db)
    x <- trimws(as.character(x))
    x %in% (get_tf_syms() %||% character(0))
  }

  # ── Load & map data ───────────────────────────────────────────────────────
  if (is.character(data) && length(data) == 1L && file.exists(data)) {
    if (verbose) message("[Δnet] Reading: ", data)
    ext <- tolower(tools::file_ext(data))
    DF <- if (ext %in% c("tsv","txt")) readr::read_tsv(data, show_col_types = FALSE) else readr::read_csv(data, show_col_types = FALSE)
  } else if (is.data.frame(data)) {
    DF <- tibble::as_tibble(data)
  } else cli::cli_abort("`data` must be a data.frame/tibble or an existing CSV/TSV path.")

  ns <- names(DF)

  # Auto-detect per-condition link scores if not provided
  if (is.null(score_ctrl_col) || is.null(score_str_col)) {
    score_cols <- ns[grepl("^link_score_", ns)]
    score_cols <- score_cols[vapply(score_cols, function(cn) is.numeric(DF[[cn]]) || is.integer(DF[[cn]]), logical(1))]
    if (length(score_cols) < 2L) {
      cli::cli_abort(c(
        "Could not auto-detect two per-condition score columns starting with 'link_score_'.",
        "x" = "Pass `score_ctrl_col` and `score_str_col` explicitly."
      ))
    }
    score_cols <- sort(score_cols)
    ctrl_idx <- grep("ctrl|control", score_cols, ignore.case = TRUE)
    if (length(ctrl_idx) == 1L) {
      score_ctrl_col <- score_cols[ctrl_idx]
      score_str_col  <- setdiff(score_cols, score_ctrl_col)[1]
    } else {
      score_ctrl_col <- score_cols[1]
      score_str_col  <- score_cols[2]
      if (verbose) message("[Δnet] Using first two score columns as (ctrl,str): ",
                           score_ctrl_col, " | ", score_str_col,
                           " (override via score_ctrl_col/score_str_col if needed)")
    }
  }

  suf_ctrl <- sub("^link_score_", "", score_ctrl_col)
  suf_str  <- sub("^link_score_", "", score_str_col)
  guess <- function(prefix, suf) {
    cand <- paste0(prefix, suf)
    if (cand %in% ns) cand else {
      rx <- paste0("^", gsub("([.*+?^${}()|[\\]\\/])", "\\\\\\1", prefix))
      hits <- ns[grepl(rx, ns)]
      norm <- function(x) gsub("[^A-Za-z0-9]+","_", x)
      sufN <- norm(suf)
      cand2 <- hits[norm(sub("^.*?_", "", hits)) == sufN]
      if (length(cand2)) cand2[1] else NA_character_
    }
  }

  sign_ctrl_col <- sign_ctrl_col %||% guess("link_sign_", suf_ctrl)
  sign_str_col  <- sign_str_col  %||% guess("link_sign_", suf_str)
  tf_expr_ctrl_col    <- tf_expr_ctrl_col    %||% guess("tf_expr_",    suf_ctrl)
  tf_expr_str_col     <- tf_expr_str_col     %||% guess("tf_expr_",    suf_str)
  gene_expr_ctrl_col  <- gene_expr_ctrl_col  %||% guess("gene_expr_",  suf_ctrl)
  gene_expr_str_col   <- gene_expr_str_col   %||% guess("gene_expr_",  suf_str)

  tf_col   <- if ("TF" %in% ns) "TF" else if ("tf" %in% ns) "tf" else "TF"
  gene_col <- if ("gene_key" %in% ns) "gene_key" else if ("gene" %in% ns) "gene" else "gene_key"
  peak_col <- if ("peak_id" %in% ns) "peak_id" else if ("peak" %in% ns) "peak" else "peak_id"
  if (!tf_col %in% ns || !gene_col %in% ns || !peak_col %in% ns) {
    cli::cli_abort("Missing required columns: need TF/tf, gene_key/gene, and peak_id/peak.")
  }
  if (!score_ctrl_col %in% ns || !score_str_col %in% ns) {
    cli::cli_abort("Missing score columns: {score_ctrl_col} / {score_str_col}. Provide correct values via score_ctrl_col / score_str_col.")
  }

  take <- unique(na.omit(c(
    tf_col, gene_col, peak_col,
    score_ctrl_col, score_str_col,
    sign_ctrl_col, sign_str_col,
    tf_expr_ctrl_col, tf_expr_str_col,
    gene_expr_ctrl_col, gene_expr_str_col
  )))
  ed <- DF[, take, drop = FALSE]
  names(ed)[match(tf_col,   names(ed))] <- "TF"
  names(ed)[match(gene_col, names(ed))] <- "gene_key"
  names(ed)[match(peak_col, names(ed))] <- "peak_id"
  names(ed)[match(score_ctrl_col, names(ed))] <- "score_ctrl"
  names(ed)[match(score_str_col,  names(ed))] <- "score_str"
  if (!is.na(sign_ctrl_col)   && sign_ctrl_col   %in% names(ed)) names(ed)[match(sign_ctrl_col,   names(ed))] <- "sign_ctrl"
  if (!is.na(sign_str_col)    && sign_str_col    %in% names(ed)) names(ed)[match(sign_str_col,    names(ed))] <- "sign_str"
  if (!is.na(tf_expr_ctrl_col)   && tf_expr_ctrl_col   %in% names(ed)) names(ed)[match(tf_expr_ctrl_col,   names(ed))] <- "tf_expr_ctrl"
  if (!is.na(tf_expr_str_col)    && tf_expr_str_col    %in% names(ed)) names(ed)[match(tf_expr_str_col,    names(ed))] <- "tf_expr_str"
  if (!is.na(gene_expr_ctrl_col) && gene_expr_ctrl_col %in% names(ed)) names(ed)[match(gene_expr_ctrl_col, names(ed))] <- "gene_expr_ctrl"
  if (!is.na(gene_expr_str_col)  && gene_expr_str_col  %in% names(ed)) names(ed)[match(gene_expr_str_col,  names(ed))] <- "gene_expr_str"

  ed$TF        <- trimws(as.character(ed$TF))
  ed$gene_key  <- trimws(as.character(ed$gene_key))
  ed$peak_id   <- trimws(as.character(ed$peak_id))
  ed$score_ctrl <- suppressWarnings(as.numeric(ed$score_ctrl))
  ed$score_str  <- suppressWarnings(as.numeric(ed$score_str))
  if (!"sign_ctrl" %in% names(ed)) ed$sign_ctrl <- NA
  if (!"sign_str"  %in% names(ed)) ed$sign_str  <- NA
  if (!"tf_expr_ctrl"   %in% names(ed)) ed$tf_expr_ctrl   <- NA_real_
  if (!"tf_expr_str"    %in% names(ed)) ed$tf_expr_str    <- NA_real_
  if (!"gene_expr_ctrl" %in% names(ed)) ed$gene_expr_ctrl <- NA_real_
  if (!"gene_expr_str"  %in% names(ed)) ed$gene_expr_str  <- NA_real_

  ed$sign_ctrl <- parse_sign(ed$sign_ctrl)
  ed$sign_str  <- parse_sign(ed$sign_str)

  ed <- ed[is.finite(ed$score_ctrl) | is.finite(ed$score_str), , drop = FALSE]
  ed <- ed[!is.na(ed$TF) & !is.na(ed$gene_key) & !is.na(ed$peak_id), , drop = FALSE]

  pass_ctrl_mag <- ifelse(is.finite(ed$score_ctrl) & abs(ed$score_ctrl) >= edge_filter_min, TRUE, FALSE)
  pass_str_mag  <- ifelse(is.finite(ed$score_str)  & abs(ed$score_str)  >= edge_filter_min, TRUE, FALSE)
  ed$any_pass   <- pass_ctrl_mag | pass_str_mag

  # Δ score and filters
  delta <- suppressWarnings(as.numeric(ed$score_str) - as.numeric(ed$score_ctrl))
  delta[!is.finite(delta)] <- 0
  ed$.delta <- delta

  if (min_delta_abs > 0) {
    ed <- ed[abs(ed$.delta) >= min_delta_abs, , drop = FALSE]
  }
  if (is.finite(keep_top_edges_per_tf) && keep_top_edges_per_tf > 0) {
    ed <- ed |>
      dplyr::group_by(TF) |>
      dplyr::slice_max(order_by = abs(.delta), n = keep_top_edges_per_tf, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  if (!nrow(ed)) cli::cli_abort("No edges left after filtering (try lowering thresholds).")

  # Optional peak handling
  if (peak_mode != "show_all") {
    tf_per_peak <- ed |>
      dplyr::group_by(peak_id) |>
      dplyr::summarise(n_tf = dplyr::n_distinct(TF), .groups = "drop")
    shared_peaks <- tf_per_peak$peak_id[tf_per_peak$n_tf > 1]
    if (peak_mode == "shared_only") {
      ed <- ed[ed$peak_id %in% shared_peaks, , drop = FALSE]
    } else { # hide
      ed <- ed |>
        dplyr::group_by(TF, gene_key) |>
        dplyr::slice_max(order_by = abs(.delta), n = 1, with_ties = FALSE) |>
        dplyr::ungroup()
    }
  }

  # cap for width/alpha scaling
  cap <- {
    x <- abs(ed$.delta[is.finite(ed$.delta)])
    q <- stats::quantile(x[x > 0], 0.995, na.rm = TRUE, names = FALSE)
    if (!is.finite(q) || q <= 0) q <- suppressWarnings(max(x, na.rm = TRUE))
    if (!is.finite(q) || q <= 0) 1 else q
  }

  # ── Nodes (TF/gene/peak) ──────────────────────────────────────────────────
  load_tf_universe(motif_db)
  node_ids <- unique(c(ed$TF, ed$gene_key))
  node_types <- ifelse(is_known_tf(node_ids), "TF", "gene")
  nodes <- tibble::tibble(id = node_ids, type = node_types)

  if (peak_mode != "hide") {
    peak_nodes <- tibble::tibble(id = paste0("PEAK:", unique(ed$peak_id)), type = "peak")
    nodes <- dplyr::bind_rows(nodes, peak_nodes)
  }

  # expression summaries for styling
  tf_ctrl <- dplyr::summarise(dplyr::group_by(ed, TF), tf_expr_ctrl = stats::median(tf_expr_ctrl, na.rm = TRUE)); names(tf_ctrl)[1] <- "id"
  tf_str  <- dplyr::summarise(dplyr::group_by(ed, TF), tf_expr_str  = stats::median(tf_expr_str,  na.rm = TRUE)); names(tf_str)[1] <- "id"
  g_ctrl  <- dplyr::summarise(dplyr::group_by(ed, gene_key), gene_expr_ctrl = stats::median(gene_expr_ctrl, na.rm = TRUE)); names(g_ctrl)[1] <- "id"
  g_str   <- dplyr::summarise(dplyr::group_by(ed, gene_key), gene_expr_str  = stats::median(gene_expr_str,  na.rm = TRUE)); names(g_str)[1] <- "id"

  nodes <- nodes |>
    dplyr::left_join(tf_ctrl, by = "id") |>
    dplyr::left_join(tf_str,  by = "id") |>
    dplyr::left_join(g_ctrl,  by = "id") |>
    dplyr::left_join(g_str,   by = "id") |>
    dplyr::mutate(
      node_raw_ctrl = dplyr::if_else(type == "TF", tf_expr_ctrl, gene_expr_ctrl),
      node_raw_str  = dplyr::if_else(type == "TF", tf_expr_str,  gene_expr_str),
      node_z_ctrl = robust_z(node_raw_ctrl),
      node_z_str  = robust_z(node_raw_str)
    )

  eps <- 1e-9
  if (de_reference == "ctrl_over_str") {
    l2fc_gene <- log2((nodes$gene_expr_ctrl + eps)/(nodes$gene_expr_str + eps))
    l2fc_tf   <- log2((nodes$tf_expr_ctrl  + eps)/(nodes$tf_expr_str  + eps))
  } else {
    l2fc_gene <- log2((nodes$gene_expr_str + eps)/(nodes$gene_expr_ctrl + eps))
    l2fc_tf   <- log2((nodes$tf_expr_str  + eps)/(nodes$tf_expr_ctrl  + eps))
  }
  nodes$l2fc <- ifelse(nodes$type == "TF", l2fc_tf, l2fc_gene)

  tf_idx_nodes <- nodes$type == "TF"
  max_abs_tf <- suppressWarnings(max(abs(nodes$l2fc[tf_idx_nodes]), na.rm = TRUE)); if (!is.finite(max_abs_tf) || max_abs_tf <= 0) max_abs_tf <- 1
  tf_signed_all <- rep(NA_real_, nrow(nodes))
  tf_signed_all[tf_idx_nodes] <- clamp01(0.5 + 0.5 * (nodes$l2fc[tf_idx_nodes] / max_abs_tf))
  ramp_div <- function(w){ w[!is.finite(w)] <- 0.5; rgb <- grDevices::colorRamp(c(col_edge_neg, "#FFFFFF", col_edge_pos))(w)/255; grDevices::rgb(rgb[,1], rgb[,2], rgb[,3]) }

  gene_z_joint <- pmax(abs(nodes$node_z_ctrl), abs(nodes$node_z_str), na.rm = TRUE)
  if (!any(is.finite(gene_z_joint))) gene_size <- rep(mean(gene_size_range), nrow(nodes)) else {
    s <- (gene_z_joint / max(gene_z_joint, na.rm = TRUE))^0.9
    gene_size <- gene_size_range[1] + s * (gene_size_range[2] - gene_size_range[1])
  }
  thr <- log2(gene_fc_thresh)
  border_col <- ifelse(!is.finite(nodes$l2fc) | abs(nodes$l2fc) < thr, col_border_none,
                       ifelse(nodes$l2fc > 0, col_border_up, col_border_down))
  border_width <- ifelse(nodes$type == "TF", 1,
                         ifelse(nodes$type == "peak", 0.5,
                                dplyr::case_when(
                                  is.finite(abs(nodes$l2fc)) & (2^abs(nodes$l2fc)) > 8 ~ 4,
                                  is.finite(abs(nodes$l2fc)) & (2^abs(nodes$l2fc)) > 4 ~ 3,
                                  is.finite(abs(nodes$l2fc)) & (2^abs(nodes$l2fc)) > 2 ~ 2,
                                  TRUE ~ 1
                                )))

  tf_median_expr <- pmax(nodes$tf_expr_ctrl %||% 0, nodes$tf_expr_str %||% 0, na.rm = TRUE)
  tf_median_expr <- tf_median_expr[tf_idx_nodes]
  tf_font <- {
    if (!length(tf_median_expr) || !any(is.finite(tf_median_expr)) || max(tf_median_expr, na.rm = TRUE) == 0) {
      rep(mean(tf_font_range), sum(tf_idx_nodes))
    } else {
      s <- (tf_median_expr / max(tf_median_expr, na.rm = TRUE))^1
      tf_font_range[1] + s * (tf_font_range[2] - tf_font_range[1])
    }
  }

  vn_nodes <- data.frame(
    id    = nodes$id,
    # show labels for TFs; keep genes unlabeled (or set to nodes$id if you want labels)
    label = ifelse(nodes$type %in% c("TF","gene"), nodes$id, ""),
    group = nodes$type,
    color.background = ifelse(nodes$type == "TF", ramp_div(tf_signed_all), col_gene_fill),
    color.border     = border_col,
    borderWidth      = border_width,
    value = ifelse(nodes$type == "gene", gene_size, NA),
    size  = ifelse(nodes$type == "peak", 6, NA),
    # make both TF and gene nodes boxes; peaks remain triangles
    shape = ifelse(nodes$type == "peak", "triangle", "box"),
    physics = FALSE,
    stringsAsFactors = FALSE
  )

  # (optional) slightly round box corners globally
  # widget <- widget |> visNetwork::visNodes(shapeProperties = list(borderRadius = 6))

  vn_nodes[is.na(vn_nodes$label), "label"] <- ""
  vn_nodes[is.na(vn_nodes$value), "value"] <- 1

  vn_nodes$title <- ifelse(
    nodes$type == "TF",
    sprintf("TF %s<br>log2FC(TF RNA): %s<br>expr_ctrl=%s<br>expr_str=%s",
            nodes$id,
            ifelse(is.finite(nodes$l2fc), sprintf("%.3f", nodes$l2fc), "NA"),
            ifelse(is.finite(nodes$tf_expr_ctrl), sprintf("%.3f", nodes$tf_expr_ctrl), "NA"),
            ifelse(is.finite(nodes$tf_expr_str),  sprintf("%.3f", nodes$tf_expr_str),  "NA")),
    ifelse(
      nodes$type == "peak",
      sprintf("Peak %s", gsub("^PEAK:","", nodes$id)),
      sprintf("Gene %s<br>expr_ctrl=%s<br>expr_str=%s",
              nodes$id,
              ifelse(is.finite(nodes$gene_expr_ctrl), sprintf("%.3f", nodes$gene_expr_ctrl), "NA"),
              ifelse(is.finite(nodes$gene_expr_str),  sprintf("%.3f", nodes$gene_expr_str),  "NA"))
    )
  )
  # ---- Adjust TF label color for dark backgrounds ----
  is_dark <- function(hex) {
    rgb <- grDevices::col2rgb(hex) / 255
    lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]  # perceptual brightness
    lum < 0.45  # threshold: lower = darker
  }

  tf_rows <- vn_nodes$group == "TF"
  if (any(tf_rows)) {
    bg_col <- vn_nodes$color.background[tf_rows]
    dark_mask <- is_dark(bg_col)
    vn_nodes$`font.color`[tf_rows] <- ifelse(dark_mask, "#FFFFFF", "#000000")
  }
  # ---- Make gene labels bold ----
  gene_rows <- vn_nodes$group == "gene"
  if (any(gene_rows)) {
    vn_nodes$`font.bold`[gene_rows] <- TRUE
  }


  # ── Edges ──────────────────────────────────────────────────────────────────
  mag <- pmin(abs(ed$.delta), cap) / max(cap, 1e-9)
  alpha <- 0.10 + mag * (0.78 - 0.10)  # lower max alpha to reduce overpaint
  width_bin <- dplyr::case_when(abs(ed$.delta) > 8 ~ 4, abs(ed$.delta) > 4 ~ 3,
                                abs(ed$.delta) > 2 ~ 2, abs(ed$.delta) > 0 ~ 1, TRUE ~ 1)
  base_hex <- ifelse(ed$.delta >= 0, col_edge_pos, col_edge_neg)
  col_rgba <- mapply(to_rgba, base_hex, alpha, USE.NAMES = FALSE)

  link_sign <- ed$sign_ctrl
  na_idx <- is.na(link_sign) | link_sign == 0L
  link_sign[na_idx] <- ed$sign_str[na_idx]

  e_tf_gene <- NULL
  if (isTRUE(add_direct)) {
    col_rgba_dir <- gsub(",\\s*[0-9.]+\\)$", ",0.28)", col_rgba)
    e_tf_gene <- tibble::tibble(
      from  = ed$TF,
      to    = ed$gene_key,
      width = width_bin,
      color = col_rgba_dir,
      dashes = FALSE,
      arrows = .make_arrows(link_sign),
      title  = sprintf("%s → %s (direct)<br>Δ (str − ctrl) = %.3f<br>ctrl=%.3f, str=%.3f",
                       ed$TF, ed$gene_key, ed$.delta, ed$score_ctrl, ed$score_str)
    )
  }

  if (peak_mode == "hide") {
    vn_edges <- e_tf_gene
  } else {
    # TF→PEAK connectors: solid, pale (no dotted lines)
    e_tf_pk <- tibble::tibble(
      from  = ed$TF,
      to    = paste0("PEAK:", ed$peak_id),
      width = 1,
      color = to_rgba(col_edge_zero, 0.35),
      dashes = FALSE,
      arrows = "",
      title  = sprintf("%s → peak %s", ed$TF, ed$peak_id)
    )
    # PEAK→gene carries Δ & the sign-based head (arrow for + ; blunt bar for −)
    e_pk_gene <- tibble::tibble(
      from  = paste0("PEAK:", ed$peak_id),
      to    = ed$gene_key,
      width = width_bin,
      color = col_rgba,
      dashes = FALSE,
      arrows = .make_arrows(link_sign),   # <-- HERE
      title  = sprintf("peak %s → %s<br>Δ (str − ctrl) = %.3f<br>ctrl=%.3f, str=%.3f",
                       ed$peak_id, ed$gene_key, ed$.delta, ed$score_ctrl, ed$score_str)
    )
    vn_edges <- dplyr::bind_rows(
      e_tf_pk,
      e_pk_gene,
      e_tf_gene %||% tibble::tibble(
        from=character(0), to=character(0), width=numeric(0), color=character(0),
        dashes=logical(0), arrows=character(0), title=character(0)
      )
    )
  }

  vn_edges$dashes <- as.logical(vn_edges$dashes)

  # ── Layout via igraph ─────────────────────────────────────────────────────
  g <- igraph::graph_from_data_frame(d = vn_edges[, c("from","to")], directed = TRUE,
                                     vertices = vn_nodes[, c("id","group"), drop = FALSE])
  layout_fun <- switch(
    layout_algo,
    fr       = igraph::layout_with_fr,
    kk       = igraph::layout_with_kk,
    lgl      = igraph::layout_with_lgl,
    dh       = igraph::layout_with_dh,
    graphopt = igraph::layout_with_graphopt,
    circle   = igraph::layout_in_circle,
    grid     = igraph::layout_on_grid,
    sphere   = igraph::layout_on_sphere,
    star     = igraph::layout_as_star,
    random   = igraph::layout_randomly
  )
  xy <- layout_fun(g); xy <- as.data.frame(xy); names(xy) <- c("x","y"); xy$id <- igraph::V(g)$name
  vn_nodes$x <- xy$x[match(vn_nodes$id, xy$id)]
  vn_nodes$y <- xy$y[match(vn_nodes$id, xy$id)]
  vn_nodes$physics <- isTRUE(physics)

  title_txt <- if (is.character(plot_title) && length(plot_title) == 1L) {
    sprintf("%s | Δ-network ", plot_title)
  } else "Δ-network"

  widget <- visNetwork::visNetwork(vn_nodes, vn_edges, width = width, height = height, main = title_txt) |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                           nodesIdSelection = FALSE) |>
    visNetwork::visEdges(smooth = FALSE) |>
    visNetwork::visNodes(borderWidth = 1, font = list(bold = TRUE),
                         margin = list(top=6,right=10,bottom=6,left=10))

  if (isTRUE(physics)) {
    widget <- widget |>
      visNetwork::visPhysics(
        enabled = TRUE,
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -50,
          centralGravity = 0.01,
          springLength = 120,
          springConstant = 0.06,
          avoidOverlap = 0.6
        ),
        stabilization = list(enabled = TRUE, iterations = 800)
      ) |>
      # Freeze physics right after stabilization so nodes stop moving
      visNetwork::visEvents(
        stabilizationIterationsDone = "function () { this.setOptions({ physics: false }); }"
      )
  } else {
    widget <- widget |>
      visNetwork::visPhysics(enabled = FALSE) |>
      visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
  }

  if (is.null(out_html)) {
    return(widget)
  } else {
    dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
    htmlwidgets::saveWidget(widget, file = out_html, selfcontained = FALSE)
    if (verbose) message("[Δnet] \u2713 wrote: ", normalizePath(out_html, mustWork = FALSE))
    invisible(out_html)
  }
}
