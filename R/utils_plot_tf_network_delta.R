# utils_plot_tf_network_delta.R ?Simplified -network plot (all links) w/ crowd controls
# Author: Yaoxiang Li
# Updated: 2026-02-17

#' Plot an interactive -network for all TFpeakgene links
#' (single-panel  = case ?control) with crowd-control options.
#'
#' The plot draws TFs (rounded boxes), genes (circles), and optional peaks (triangles).
#' Edge color encodes  magnitude/sign (red = gain/activation; blue = loss/repression);
#' arrow head styles encode sign: arrow for activation (+), blunt bar for repression (?.
#' Node fill color encodes expression direction (up / down / neutral).
#'
#' @param data Tibble/data.frame or path to CSV/TSV with per-row TFpeakgene entry.
#'   Must contain at least: TF/tf, gene_key/gene, peak_id/peak, and two per-condition
#'   link scores (e.g., \code{link_score_*}). Signs/expressions are optional.
#' @param plot_title Title string shown above the network.
#' @param out_html Optional path to write HTML; if `NULL` returns the htmlwidget.
#' @param layout_algo One of: "fr","kk","lgl","dh","graphopt","circle","grid","sphere","star","random".
#' @param physics Logical; if `TRUE`, enables a repulsive physics preset to relieve crowding.
#' @param add_direct Logical; also draw faint direct TFgene  edges in addition to TFPEAKgene.
#' @param edge_filter_min Keep edges with \code{|score_ctrl|} or \code{|score_str|} ?this value.
#' @param min_delta_abs Keep edges with \code{||} ?this value ( = case ?ctrl).
#' @param keep_top_edges_per_tf Optional integer; per TF keep only the strongest \code{||} edges.
#' @param peak_mode One of: "show_all" (default), "shared_only", "hide".
#' @param show_peaks Convenience boolean; if FALSE, forces \code{peak_mode = "hide"}.
#' @param gene_fc_thresh Threshold on |log2FC| used to call nodes up/down vs neutral.
#' @param col_node_up Fill color for up-regulated nodes (TFs/genes).
#' @param col_node_down Fill color for down-regulated nodes.
#' @param col_node_neutral Fill color for neutral nodes.
#' @param font_color Base label color (all nodes).
#' @param edge_width_min,edge_width_max Minimum and maximum edge thickness for ||.
#' @param de_reference "str_over_ctrl" (default) or "ctrl_over_str" for log2FC sign.
#' @param motif_db "jaspar2024" (default) or "hocomocov13" to recognize TF symbols.
#' @param width,height Dimensions passed to \code{visNetwork}.
#' @param seed Integer seed for reproducible igraph layouts.
#' @param verbose Logical; print progress messages.
#' @param score_ctrl_col,score_str_col Optional explicit column names of per-condition link scores.
#' @param sign_ctrl_col,sign_str_col Optional explicit column names of per-condition link signs.
#' @param tf_expr_ctrl_col,tf_expr_str_col,gene_expr_ctrl_col,gene_expr_str_col Optional expression columns.
#' @param log2fc_tf_col,log2fc_gene_col Precomputed log2FC columns used for node direction/color and TF tooltip.
#'   Defaults are `log2FC_tf_expr` and `log2FC_gene_expr`. This function will not recompute these values.
#'
#' @return A \code{visNetwork} htmlwidget (or, invisibly, \code{out_html} path if saved).
#'
#' @export
plot_tf_network_delta <- function(
    data,
    plot_title,
    out_html       = NULL,
    html_title     = NULL,
    layout_algo    = c("fr","kk","lgl","dh","graphopt","circle","grid","sphere","star","random"),
    add_direct     = TRUE,
    edge_filter_min = 0,
    min_delta_abs   = 0,
    keep_top_edges_per_tf = NULL,
    peak_mode      = c("show_all","shared_only","hide"),
    show_peaks     = TRUE,
    gene_fc_thresh = 1.5,
    col_node_up      = "#dfb734",
    col_node_down    = "#283673",
    col_node_neutral = "#555555",
    font_color       = "#FFFFFF",
    edge_width_min   = 2,
    edge_width_max   = 8,
    size_by        = c("log2fc", "expr_max"),
    label_inside   = TRUE,
    label_max_chars = 18,
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
    gene_expr_str_col  = NULL,
    log2fc_tf_col      = "log2FC_tf_expr",
    log2fc_gene_col    = "log2FC_gene_expr"
){
  layout_algo  <- match.arg(layout_algo)
  de_reference <- match.arg(de_reference)
  motif_db     <- match.arg(motif_db)
  peak_mode    <- match.arg(peak_mode)
  size_by      <- match.arg(size_by)
  if (identical(show_peaks, FALSE)) peak_mode <- "hide"
  set.seed(as.integer(seed) %% .Machine$integer.max)

  #  Styling
  col_edge_pos    <- "#df1d16"  #  > 0
  col_edge_neg    <- "#255dae"  #  < 0
  col_edge_zero   <- "#C8C8C8"
  tf_font_range   <- c(24, 56)
  tf_box_aspect   <- 1.6
  gene_size_range <- c(10, 50)
  tf_size_range   <- gene_size_range
  gene_font_range <- c(
    max(10, round(gene_size_range[1] * 0.4)),
    max(18, round(gene_size_range[2] * 0.6))
  )

  `%||%` <- function(a, b) if (is.null(a)) b else a
  to_rgba <- function(hex, alpha){
    hex <- gsub("^#", "", hex)
    r <- strtoi(substr(hex,1,2),16L)
    g <- strtoi(substr(hex,3,4),16L)
    b <- strtoi(substr(hex,5,6),16L)
    sprintf("rgba(%d,%d,%d,%.3f)", r,g,b, max(0, min(1, alpha)))
  }
  robust_z <- function(x){
    x <- as.numeric(x)
    m <- stats::median(x, na.rm = TRUE)
    madv <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(madv) || madv == 0) {
      sx <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
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
    kinds <- ifelse(is.na(sign_vec) | sign_vec >= 0L, "arrow", "bar")
    I(lapply(kinds, function(tp) list(to = list(enabled = TRUE, type = tp))))
  }

  # Known TF universe (from installed package files)
  .tf_cache <- new.env(parent = emptyenv())
  get_tf_syms <- function() get0("KNOWN_TF_SYMBOLS", envir = .tf_cache, inherits = FALSE)
  set_tf_syms <- function(x, db) {
    assign("KNOWN_TF_SYMBOLS", x, envir = .tf_cache)
    assign("KNOWN_TF_DB", x = db, envir = .tf_cache)
    invisible(TRUE)
  }
  load_tf_universe <- function(db = motif_db){
    cur_db <- get0("KNOWN_TF_DB", envir = .tf_cache, inherits = FALSE)
    cur    <- get_tf_syms()
    if (!is.null(cur) && identical(cur_db, db)) return(invisible(TRUE))
    ref_genome <- get0("ref_genome", envir = .GlobalEnv, inherits = FALSE)
    if (is.null(ref_genome) || !nzchar(ref_genome)) ref_genome <- "hg38"
    ref_genome <- tolower(as.character(ref_genome))
    path <- if (db == "jaspar2024") {
      if (ref_genome == "mm10") {
        system.file("extdata","genome","JASPAR2024_mm10.txt", package = "episcope")
      } else {
        system.file("extdata","genome","JASPAR2024_hg38.txt", package = "episcope")
      }
    } else {
      system.file("extdata","genome","HOCOMOCOv13.txt", package = "episcope")
    }
    if (!nzchar(path) || !file.exists(path)) {
      # Fallbacks for legacy filenames
      if (db == "jaspar2024") {
        path <- system.file("extdata","genome","JASPAR2024.txt", package = "episcope")
      }
    }
    if (!nzchar(path) || !file.exists(path)) {
      .log_abort("Motif cluster definition not found for motif_db='{db}'.")
    }
    df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    if ("HGNC" %in% names(df) && !("gene_symbol" %in% names(df))) {
      df <- dplyr::rename(df, gene_symbol = HGNC)
    }
    if (!"gene_symbol" %in% names(df)) {
      .log_abort("Motif DB file missing gene_symbol column: {path}")
    }
    syms <- unique(trimws(unlist(strsplit(df$gene_symbol, "::", fixed = TRUE), use.names = FALSE)))
    set_tf_syms(syms, db); invisible(TRUE)
  }
  is_known_tf <- function(x){
    if (is.null(get_tf_syms())) load_tf_universe(motif_db)
    x <- trimws(as.character(x))
    x %in% (get_tf_syms() %||% character(0))
  }

  #  Load & map data
  if (is.character(data) && length(data) == 1L && file.exists(data)) {
    if (verbose) message("[net] Reading: ", data)
    ext <- tolower(tools::file_ext(data))
    DF <- if (ext %in% c("tsv","txt")) {
      readr::read_tsv(data, show_col_types = FALSE)
    } else {
      readr::read_csv(data, show_col_types = FALSE)
    }
  } else if (is.data.frame(data)) {
    DF <- tibble::as_tibble(data)
  } else {
    .log_abort("`data` must be a data.frame/tibble or an existing CSV/TSV path.")
  }

  ns <- names(DF)

  # Auto-detect per-condition link scores if not provided
  if (is.null(score_ctrl_col) || is.null(score_str_col)) {
    score_cols <- ns[grepl("^link_score_", ns)]
    score_cols <- score_cols[vapply(
      score_cols,
      function(cn) is.numeric(DF[[cn]]) || is.integer(DF[[cn]]),
      logical(1)
    )]
    if (length(score_cols) < 2L) {
      .log_abort(c(
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
      if (verbose) {
        message("[net] Using first two score columns as (ctrl,case): ",
                score_ctrl_col, " | ", score_str_col,
                " (override via score_ctrl_col/score_str_col if needed)")
      }
    }
  }

  suf_ctrl <- sub("^link_score_", "", score_ctrl_col)
  suf_str  <- sub("^link_score_", "", score_str_col)
  guess <- function(prefix, suf) {
    cand <- paste0(prefix, suf)
    if (cand %in% ns) {
      cand
    } else {
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
    .log_abort("Missing required columns: need TF/tf, gene_key/gene, and peak_id/peak.")
  }
  if (!score_ctrl_col %in% ns || !score_str_col %in% ns) {
    .log_abort("Missing score columns: {score_ctrl_col} / {score_str_col}. Provide correct values via score_ctrl_col / score_str_col.")
  }
  if (!is.character(log2fc_tf_col) || length(log2fc_tf_col) != 1L || !nzchar(log2fc_tf_col) || !(log2fc_tf_col %in% ns)) {
    .log_abort("Missing required precomputed TF log2FC column: {log2fc_tf_col}. Provide correct `log2fc_tf_col`; this plot does not recompute TF log2FC.")
  }
  if (!is.character(log2fc_gene_col) || length(log2fc_gene_col) != 1L || !nzchar(log2fc_gene_col) || !(log2fc_gene_col %in% ns)) {
    .log_abort("Missing required precomputed gene log2FC column: {log2fc_gene_col}. Provide correct `log2fc_gene_col`; this plot does not recompute gene log2FC.")
  }

  take <- unique(na.omit(c(
    tf_col, gene_col, peak_col,
    score_ctrl_col, score_str_col,
    sign_ctrl_col, sign_str_col,
    tf_expr_ctrl_col, tf_expr_str_col,
    gene_expr_ctrl_col, gene_expr_str_col,
    log2fc_tf_col, log2fc_gene_col
  )))
  ed <- DF[, take, drop = FALSE]
  names(ed)[match(tf_col,   names(ed))] <- "TF"
  names(ed)[match(gene_col, names(ed))] <- "gene_key"
  names(ed)[match(peak_col, names(ed))] <- "peak_id"
  names(ed)[match(score_ctrl_col, names(ed))] <- "score_ctrl"
  names(ed)[match(score_str_col,  names(ed))] <- "score_str"
  if (!is.na(sign_ctrl_col)   && sign_ctrl_col   %in% names(ed)) {
    names(ed)[match(sign_ctrl_col,   names(ed))] <- "sign_ctrl"
  }
  if (!is.na(sign_str_col)    && sign_str_col    %in% names(ed)) {
    names(ed)[match(sign_str_col,    names(ed))] <- "sign_str"
  }
  if (!is.na(tf_expr_ctrl_col)   && tf_expr_ctrl_col   %in% names(ed)) {
    names(ed)[match(tf_expr_ctrl_col,   names(ed))] <- "tf_expr_ctrl"
  }
  if (!is.na(tf_expr_str_col)    && tf_expr_str_col    %in% names(ed)) {
    names(ed)[match(tf_expr_str_col,    names(ed))] <- "tf_expr_str"
  }
  if (!is.na(gene_expr_ctrl_col) && gene_expr_ctrl_col %in% names(ed)) {
    names(ed)[match(gene_expr_ctrl_col, names(ed))] <- "gene_expr_ctrl"
  }
  if (!is.na(gene_expr_str_col)  && gene_expr_str_col  %in% names(ed)) {
    names(ed)[match(gene_expr_str_col,  names(ed))] <- "gene_expr_str"
  }
  if (!is.na(log2fc_tf_col)   && log2fc_tf_col   %in% names(ed)) {
    names(ed)[match(log2fc_tf_col,   names(ed))] <- "log2fc_tf_expr"
  }
  if (!is.na(log2fc_gene_col) && log2fc_gene_col %in% names(ed)) {
    names(ed)[match(log2fc_gene_col, names(ed))] <- "log2fc_gene_expr"
  }

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
  if (!"log2fc_tf_expr" %in% names(ed)) .log_abort("Missing required precomputed TF log2FC values after mapping.")
  if (!"log2fc_gene_expr" %in% names(ed)) .log_abort("Missing required precomputed gene log2FC values after mapping.")
  ed$log2fc_tf_expr <- suppressWarnings(as.numeric(ed$log2fc_tf_expr))
  ed$log2fc_gene_expr <- suppressWarnings(as.numeric(ed$log2fc_gene_expr))

  ed$sign_ctrl <- parse_sign(ed$sign_ctrl)
  ed$sign_str  <- parse_sign(ed$sign_str)

  ed <- ed[is.finite(ed$score_ctrl) | is.finite(ed$score_str), , drop = FALSE]
  ed <- ed[!is.na(ed$TF) & !is.na(ed$gene_key) & !is.na(ed$peak_id), , drop = FALSE]

  pass_ctrl_mag <- ifelse(is.finite(ed$score_ctrl) & abs(ed$score_ctrl) >= edge_filter_min, TRUE, FALSE)
  pass_str_mag  <- ifelse(is.finite(ed$score_str)  & abs(ed$score_str)  >= edge_filter_min, TRUE, FALSE)
  ed$any_pass   <- pass_ctrl_mag | pass_str_mag

  #  score and filters
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

  if (!nrow(ed)) .log_abort("No edges left after filtering (try lowering thresholds).")

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

  #  Nodes (TF/gene/peak)
  load_tf_universe(motif_db)
  node_ids   <- unique(c(ed$TF, ed$gene_key))
  node_types <- ifelse(is_known_tf(node_ids), "TF", "gene")
  nodes <- tibble::tibble(id = node_ids, type = node_types)

  if (peak_mode != "hide") {
    peak_nodes <- tibble::tibble(id = paste0("PEAK:", unique(ed$peak_id)), type = "peak")
    nodes <- dplyr::bind_rows(nodes, peak_nodes)
  }

  .median_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    stats::median(x)
  }

  tf_ctrl <- dplyr::summarise(dplyr::group_by(ed, TF),
                              tf_expr_ctrl = .median_finite(tf_expr_ctrl))
  names(tf_ctrl)[1] <- "id"
  tf_str  <- dplyr::summarise(dplyr::group_by(ed, TF),
                              tf_expr_str  = .median_finite(tf_expr_str))
  names(tf_str)[1] <- "id"
  g_ctrl  <- dplyr::summarise(dplyr::group_by(ed, gene_key),
                              gene_expr_ctrl = .median_finite(gene_expr_ctrl))
  names(g_ctrl)[1] <- "id"
  g_str   <- dplyr::summarise(dplyr::group_by(ed, gene_key),
                              gene_expr_str  = .median_finite(gene_expr_str))
  names(g_str)[1] <- "id"
  tf_l2fc <- dplyr::summarise(
    dplyr::group_by(ed, TF),
    l2fc_tf = .median_finite(log2fc_tf_expr)
  )
  names(tf_l2fc)[1] <- "id"
  g_l2fc <- dplyr::summarise(
    dplyr::group_by(ed, gene_key),
    l2fc_gene = .median_finite(log2fc_gene_expr)
  )
  names(g_l2fc)[1] <- "id"

  nodes <- nodes |>
    dplyr::left_join(tf_ctrl, by = "id") |>
    dplyr::left_join(tf_str,  by = "id") |>
    dplyr::left_join(g_ctrl,  by = "id") |>
    dplyr::left_join(g_str,   by = "id") |>
    dplyr::left_join(tf_l2fc, by = "id") |>
    dplyr::left_join(g_l2fc,  by = "id") |>
    dplyr::mutate(
      # For nodes classified as TF by motif_db but never used as TF in `ed$TF`,
      # borrow gene_expr_* as TF expression so we don't end up with NA.
      tf_expr_ctrl = dplyr::if_else(
        type == "TF" & is.na(tf_expr_ctrl) & !is.na(gene_expr_ctrl),
        gene_expr_ctrl,
        tf_expr_ctrl
      ),
      tf_expr_str = dplyr::if_else(
        type == "TF" & is.na(tf_expr_str) & !is.na(gene_expr_str),
        gene_expr_str,
        tf_expr_str
      ),
      node_raw_ctrl = dplyr::if_else(type == "TF", tf_expr_ctrl, gene_expr_ctrl),
      node_raw_str  = dplyr::if_else(type == "TF", tf_expr_str,  gene_expr_str),
      node_z_ctrl   = robust_z(node_raw_ctrl),
      node_z_str    = robust_z(node_raw_str),
      l2fc = dplyr::if_else(type == "TF", l2fc_tf, l2fc_gene)
    )

  # TF font size based on |log2FC| ---------------------------------

  tf_idx_nodes <- nodes$type == "TF"
  tf_fc_mag    <- abs(nodes$l2fc[tf_idx_nodes])

  tf_font <- {
    if (!length(tf_fc_mag) || !any(is.finite(tf_fc_mag)) ||
        max(tf_fc_mag, na.rm = TRUE) == 0) {
      rep(mean(tf_font_range), sum(tf_idx_nodes))
    } else {
      s <- tf_fc_mag / max(tf_fc_mag, na.rm = TRUE)
      tf_font_range[1] + s * (tf_font_range[2] - tf_font_range[1])
    }
  }


  # node size based on log2FC or expression -------------------------

  expr_ctrl <- nodes$node_raw_ctrl
  expr_str  <- nodes$node_raw_str
  expr_max  <- pmax(expr_ctrl, expr_str)
  if (size_by == "log2fc") {
    size_signal <- abs(nodes$l2fc)
  } else {
    size_signal <- log1p(pmax(expr_max, 0))
  }

  scale_range <- function(signal, range) {
    if (!length(signal) || !any(is.finite(signal)) || max(signal, na.rm = TRUE) == 0) {
      return(rep(mean(range), length(signal)))
    }
    s <- signal / max(signal, na.rm = TRUE)
    range[1] + s * (range[2] - range[1])
  }

  gene_idx <- nodes$type == "gene"
  tf_idx   <- nodes$type == "TF"

  gene_size <- rep(NA_real_, nrow(nodes))
  tf_value  <- rep(NA_real_, nrow(nodes))
  gene_size[gene_idx] <- scale_range(size_signal[gene_idx], gene_size_range)
  tf_value[tf_idx]    <- scale_range(size_signal[tf_idx], tf_size_range)
  size_scaled <- rep(NA_real_, nrow(nodes))
  size_scaled[gene_idx] <- gene_size[gene_idx]
  size_scaled[tf_idx]   <- tf_value[tf_idx]


  thr <- log2(gene_fc_thresh)
  node_state <- ifelse(
    !is.finite(nodes$l2fc) | abs(nodes$l2fc) < thr, "neutral",
    ifelse(nodes$l2fc > 0, "up", "down")
  )
  fill_col <- ifelse(
    node_state == "up",   col_node_up,
    ifelse(node_state == "down", col_node_down, col_node_neutral)
  )

  border_width <- ifelse(nodes$type == "peak", 0.5, 0)

  label_text <- nodes$id
  if (isTRUE(label_inside) && is.finite(label_max_chars) && label_max_chars > 0) {
    label_text <- ifelse(
      nchar(label_text) > label_max_chars,
      paste0(substr(label_text, 1, label_max_chars - 3), "..."),
      label_text
    )
  }

  vn_nodes <- data.frame(
    id    = nodes$id,
    label = ifelse(nodes$type %in% c("TF","gene"), label_text, ""),
    group = nodes$type,
    color.background = fill_col,
    color.border     = fill_col,
    borderWidth      = border_width,
    value = ifelse(nodes$type == "gene", gene_size,
                   ifelse(nodes$type == "TF", tf_value, NA)),
    size  = ifelse(nodes$type == "peak", 6, NA),
    # Shapes: label-inside vs label-outside based on label_inside
    shape = ifelse(
      nodes$type == "peak", "triangle",
      ifelse(
        nodes$type == "TF",
        if (isTRUE(label_inside)) "box" else "square",
        if (isTRUE(label_inside)) "circle" else "dot"
      )
    ),
    physics = FALSE,
    stringsAsFactors = FALSE
  )

  vn_nodes[is.na(vn_nodes$label), "label"] <- ""
  vn_nodes[is.na(vn_nodes$value), "value"] <- 1

  vn_nodes$title <- ifelse(
    nodes$type == "TF",
    sprintf("TF %s<br>log2FC(TF RNA, case/ctrl): %s<br>expr_ctrl=%s<br>expr_case=%s<br>expr_max=%s<br>size_signal=%s<br>size_scaled=%s",
            nodes$id,
            ifelse(is.finite(nodes$l2fc), sprintf("%.3f", nodes$l2fc), "NA"),
            ifelse(is.finite(nodes$tf_expr_ctrl), sprintf("%.3f", nodes$tf_expr_ctrl), "NA"),
            ifelse(is.finite(nodes$tf_expr_str),  sprintf("%.3f", nodes$tf_expr_str),  "NA"),
            ifelse(is.finite(pmax(nodes$tf_expr_ctrl, nodes$tf_expr_str)),
                   sprintf("%.3f", pmax(nodes$tf_expr_ctrl, nodes$tf_expr_str)), "NA"),
            ifelse(is.finite(size_signal), sprintf("%.3f", size_signal), "NA"),
            ifelse(is.finite(size_scaled), sprintf("%.3f", size_scaled), "NA")),
    ifelse(
      nodes$type == "peak",
      sprintf("Peak %s", gsub("^PEAK:","", nodes$id)),
      sprintf("Gene %s<br>expr_ctrl=%s<br>expr_case=%s<br>expr_max=%s<br>size_signal=%s<br>size_scaled=%s",
              nodes$id,
              ifelse(is.finite(nodes$gene_expr_ctrl), sprintf("%.3f", nodes$gene_expr_ctrl), "NA"),
              ifelse(is.finite(nodes$gene_expr_str),  sprintf("%.3f", nodes$gene_expr_str),  "NA"),
              ifelse(is.finite(pmax(nodes$gene_expr_ctrl, nodes$gene_expr_str)),
                     sprintf("%.3f", pmax(nodes$gene_expr_ctrl, nodes$gene_expr_str)), "NA"),
              ifelse(is.finite(size_signal), sprintf("%.3f", size_signal), "NA"),
              ifelse(is.finite(size_scaled), sprintf("%.3f", size_scaled), "NA"))
    )
  )

  node_sizes <- data.frame(
    id = nodes$id,
    type = nodes$type,
    expr_ctrl = expr_ctrl,
    expr_case = expr_str,
    expr_max = expr_max,
    size_signal = size_signal,
    size_scaled = size_scaled,
    stringsAsFactors = FALSE
  )

  tf_rows   <- vn_nodes$group == "TF"
  gene_rows <- vn_nodes$group == "gene"

  vn_nodes$`font.size`        <- 20
  vn_nodes$`font.bold`        <- FALSE
  vn_nodes$`font.color`       <- font_color
  vn_nodes$`font.strokeWidth` <- 2
  vn_nodes$`font.strokeColor` <- vn_nodes$color.background

  if (any(tf_rows)) {
    vn_nodes$`font.bold`[tf_rows]  <- TRUE
    if (isTRUE(label_inside)) {
      tf_sizes <- size_scaled[tf_rows]
      tf_sizes[!is.finite(tf_sizes)] <- NA_real_
      tf_labels <- label_text[tf_rows]
      label_chars <- pmax(1, nchar(tf_labels))
      tf_font <- tf_font * 1.15
      max_font_by_label <- (tf_sizes * 2) / (0.5 * label_chars)
      tf_font <- pmin(tf_font, max_font_by_label)
      tf_font <- pmax(tf_font, 10)
    }
    vn_nodes$`font.size`[tf_rows]  <- tf_font
  }
  if (any(gene_rows)) {
    vn_nodes$`font.bold`[gene_rows] <- TRUE
    if (isTRUE(label_inside)) {
      gene_sizes <- size_scaled[gene_rows]
      gene_sizes[!is.finite(gene_sizes)] <- NA_real_
      if (any(is.finite(gene_sizes))) {
        min_s <- min(gene_sizes, na.rm = TRUE)
        max_s <- max(gene_sizes, na.rm = TRUE)
        if (max_s <= min_s) {
          gene_font <- rep(12, sum(gene_rows))
        } else {
          s <- (gene_sizes - min_s) / (max_s - min_s)
          gene_font <- gene_font_range[1] + s * (gene_font_range[2] - gene_font_range[1])
        }
      } else {
        gene_font <- rep(12, sum(gene_rows))
      }
      gene_labels <- label_text[gene_rows]
      label_chars <- pmax(1, nchar(gene_labels))
      max_font_by_label <- (gene_sizes * 2) / (0.5 * label_chars)
      gene_font <- pmin(gene_font, max_font_by_label)
      gene_font <- pmax(gene_font, 8)
      vn_nodes$`font.size`[gene_rows] <- gene_font
    } else {
      vn_nodes$`font.size`[gene_rows] <- 22
    }
  }

  if (isTRUE(label_inside)) {
    size_px <- ifelse(is.finite(size_scaled), round(size_scaled * 2), NA_integer_)
    vn_nodes$widthConstraint <- rep(NA_integer_, nrow(vn_nodes))
    vn_nodes$heightConstraint <- rep(NA_integer_, nrow(vn_nodes))
    gene_rows <- vn_nodes$group == "gene"
    tf_rows <- vn_nodes$group == "TF"
    if (any(gene_rows)) {
      vn_nodes$widthConstraint[gene_rows] <- size_px[gene_rows]
      vn_nodes$heightConstraint[gene_rows] <- size_px[gene_rows]
    }
    if (any(tf_rows)) {
      tf_height <- size_px[tf_rows]
      tf_width <- round(tf_height * tf_box_aspect)
      vn_nodes$widthConstraint[tf_rows] <- tf_width
      vn_nodes$heightConstraint[tf_rows] <- tf_height
    }
  }

  #  Edges
  mag <- pmin(abs(ed$.delta), cap) / max(cap, 1e-9)
  alpha <- 0.10 + mag * (0.78 - 0.10)

  edge_width <- edge_width_min + mag * (edge_width_max - edge_width_min)

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
      width = edge_width,
      color = col_rgba_dir,
      dashes = FALSE,
      arrows = .make_arrows(link_sign),
      title  = sprintf("%s ?%s (direct)<br> (case ?ctrl) = %.3f<br>ctrl=%.3f, case=%.3f",
                       ed$TF, ed$gene_key, ed$.delta, ed$score_ctrl, ed$score_str)
    )
  }

  if (peak_mode == "hide") {
    vn_edges <- e_tf_gene
  } else {
    e_tf_pk <- tibble::tibble(
      from  = ed$TF,
      to    = paste0("PEAK:", ed$peak_id),
      width = 1,
      color = to_rgba(col_edge_zero, 0.35),
      dashes = FALSE,
      arrows = "",
      title  = sprintf("%s ?peak %s", ed$TF, ed$peak_id)
    )
    e_pk_gene <- tibble::tibble(
      from  = paste0("PEAK:", ed$peak_id),
      to    = ed$gene_key,
      width = edge_width,
      color = col_rgba,
      dashes = FALSE,
      arrows = .make_arrows(link_sign),
      title  = sprintf("peak %s ?%s<br> (case ?ctrl) = %.3f<br>ctrl=%.3f, case=%.3f",
                       ed$peak_id, ed$gene_key, ed$.delta, ed$score_ctrl, ed$score_str)
    )
    vn_edges <- dplyr::bind_rows(
      e_tf_pk,
      e_pk_gene,
      e_tf_gene %||% tibble::tibble(
        from   = character(0), to = character(0),
        width  = numeric(0),  color = character(0),
        dashes = logical(0),  arrows = character(0),
        title  = character(0)
      )
    )
  }

  vn_edges$dashes <- as.logical(vn_edges$dashes)

  #  Layout via igraph
  g <- igraph::graph_from_data_frame(
    d = vn_edges[, c("from","to")],
    directed = TRUE,
    vertices = vn_nodes[, c("id","group"), drop = FALSE]
  )
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
  xy <- layout_fun(g)
  xy <- as.data.frame(xy)
  names(xy) <- c("x","y")
  xy$id <- igraph::V(g)$name
  vn_nodes$x <- xy$x[match(vn_nodes$id, xy$id)]
  vn_nodes$y <- xy$y[match(vn_nodes$id, xy$id)]
  vn_nodes$physics <- isTRUE(physics)

  title_txt <- if (is.character(plot_title) && length(plot_title) == 1L) {
    sprintf("%s | -network ", plot_title)
  } else {
    "-network"
  }

  tf_ids <- sort(nodes$id[nodes$type == "TF"])
  gene_ids <- sort(nodes$id[nodes$type == "gene"])
  select_ids <- c(tf_ids, gene_ids)

  widget <- visNetwork::visNetwork(
    vn_nodes, vn_edges,
    width = width, height = height, main = title_txt
  ) |>
    visNetwork::visLayout(randomSeed = as.integer(seed)) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = list(enabled = TRUE, values = select_ids),
      manipulation = TRUE
    ) |>
    visNetwork::visEdges(smooth = FALSE) |>
    visNetwork::visNodes(
      font = list(bold = TRUE),
      margin = list(top = 6, right = 10, bottom = 6, left = 10),
      scaling = list(
        min = min(gene_size_range[1], tf_size_range[1]),
        max = max(gene_size_range[2], tf_size_range[2])
      ),
      shapeProperties = list(borderRadius = 8)
    )

  if (isTRUE(physics)) {
    widget <- widget |>
      visNetwork::visPhysics(
        enabled = TRUE,
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -55,
          centralGravity = 0.01,
          springLength = 160,
          springConstant = 0.07,
          avoidOverlap = 0.8
        ),
        stabilization = list(enabled = TRUE, iterations = 1500)
      ) |>
      visNetwork::visEvents(
        stabilizationIterationsDone =
          "function () { this.setOptions({ physics: false }); }"
      )
  } else {
    widget <- widget |>
      visNetwork::visPhysics(enabled = FALSE) |>
      visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE)
  }

  attr(widget, "node_sizes") <- node_sizes
  attr(widget, "node_size_by") <- size_by

  if (is.null(html_title) || !nzchar(html_title)) {
    html_title <- plot_title
  }

  if (is.null(out_html)) {
    return(widget)
  } else {
    dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
    htmlwidgets::saveWidget(widget, file = out_html, selfcontained = FALSE)
    .set_html_title(out_html, html_title)
    if (verbose) message("[net] ?wrote: ", normalizePath(out_html, mustWork = FALSE))
    invisible(out_html)
  }
}

.set_html_title <- function(html_path, title_text) {
  if (is.null(title_text) || !nzchar(title_text)) return(invisible(FALSE))
  if (!file.exists(html_path)) return(invisible(FALSE))
  html_txt <- readLines(html_path, warn = FALSE)
  if (!any(grepl("<title>", html_txt, fixed = TRUE))) return(invisible(FALSE))
  html_txt <- sub("<title>.*</title>", paste0("<title>", title_text, "</title>"), html_txt)
  writeLines(html_txt, html_path)
  invisible(TRUE)
}
