# BEGIN FILE: utils_grn_link_delta_loader.R
# nocov start

# ------------------------------ Constants ------------------------------------
.col_tf_fill_low   <- "#4575b4"  # blue
.col_tf_fill_high  <- "#d73027"  # red
.col_edge_neg      <- "#6f9ed4"  # Δ < 0 (blue-ish)
.col_edge_pos      <- "#d86b5e"  # Δ > 0 (red-ish)

.gene_size_range   <- c(8, 24)   # visNetwork node value for genes
.gene_size_gamma   <- 0.7
.tf_font_range     <- c(18, 30)  # px
.tf_font_gamma     <- 0.8
.col_border_none   <- "#BDBDBD"
.col_border_up     <- "#c0392b"
.col_border_down   <- "#2980b9"
.peak_node_value   <- 6
.tf_label_halo_w   <- 3
.tf_label_halo_c   <- "#FFFFFF80"
.peak_ring_inner_offset_px <- 40

.tf_peak_dash_pattern <- c(2, 8)  # dotted

# ------------------------------ Utilities ------------------------------------
.llog <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    msg <- paste0(format(Sys.time(), "%H:%M:%S"), " ", paste(list(...), collapse = ""))
    utils::message(msg)
  }
}

clamp <- function(x, lo, hi) pmin(hi, pmax(lo, x))

robust_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || all(!is.finite(x))) return(rep(NA_real_, length(x)))
  m <- stats::median(x, na.rm = TRUE)
  mad <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) return(rep(0, length(x)))
  (x - m) / mad
}

.as_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(y %in% c("1","true","t","yes","y"))
  }
  rep(NA, length(x))
}

.as_sign <- function(x) {
  # Accept {-1, 0, +1}, or { -, + }, or {repressor, activator}
  if (is.numeric(x)) {
    y <- sign(x)
    y[!is.finite(y)] <- 0
    return(as.integer(y))
  }
  if (is.character(x)) {
    y <- tolower(trimws(x))
    out <- integer(length(y))
    out[y %in% c("-","neg","repressor","down")] <- -1L
    out[y %in% c("+","pos","activator","up")]   <-  1L
    out[!(y %in% c("-","neg","repressor","down","+","pos","activator","up"))] <- 0L
    return(out)
  }
  rep(0L, length(x))
}

.compute_pass_flags <- function(df, thr, score_col, active_col = NULL) {
  sc <- suppressWarnings(as.numeric(df[[score_col]]))
  pass <- is.finite(sc) & (sc >= thr)
  if (!is.null(active_col) && active_col %in% names(df)) {
    ac <- .as_bool(df[[active_col]])
    pass <- pass & !isFALSE(ac)  # if NA, keep result driven by score
  }
  pass
}

.compute_caps <- function(v_pos, v_neg, q = 0.995) {
  v_pos <- v_pos[is.finite(v_pos) & v_pos > 0]
  v_neg <- -v_neg[is.finite(v_neg) & v_neg < 0]
  cpos <- if (length(v_pos)) stats::quantile(v_pos, q, names = FALSE) else 1
  cneg <- if (length(v_neg)) stats::quantile(v_neg, q, names = FALSE) else 1
  if (!is.finite(cpos) || cpos <= 0) cpos <- 1
  if (!is.finite(cneg) || cneg <= 0) cneg <- 1
  list(pos = cpos, neg = cneg)
}

# ------------------------------ Filename & columns ---------------------------
.guess_pair_from_path <- function(path) {
  base <- basename(path)
  base <- sub("_delta_links.*$", "", base, perl = TRUE)
  bits <- strsplit(base, "_vs_", fixed = TRUE)[[1]]
  list(stress_tag = bits[[1]], ctrl_tag = bits[[2]])
}

.find_tagged_col <- function(df, prefix, tag) {
  exact <- paste0(prefix, tag)
  if (exact %in% names(df)) return(exact)
  esc_tag <- gsub("([][{}()+*^$.|?\\\\-])", "\\\\\\1", tag, perl = TRUE)
  esc_pre <- gsub("([][{}()+*^$.|?\\\\-])", "\\\\\\1", prefix, perl = TRUE)
  pat <- paste0("^", esc_pre, esc_tag, "$")
  hits <- grep(pat, names(df), value = TRUE, perl = TRUE)
  if (length(hits) >= 1) return(hits[1])
  stop("Cannot find column for prefix='", prefix, "' and tag='", tag, "'. Available: ",
       paste(names(df), collapse = ", "))
}

.resolve_pair_columns <- function(C, stress_tag, ctrl_tag) {
  list(
    sc_str   = .find_tagged_col(C, "link_score_", stress_tag),
    sc_ctrl  = .find_tagged_col(C, "link_score_", ctrl_tag),
    tf_str   = .find_tagged_col(C, "tf_expr_",    stress_tag),
    tf_ctrl  = .find_tagged_col(C, "tf_expr_",    ctrl_tag),
    ge_str   = .find_tagged_col(C, "gene_expr_",  stress_tag),
    ge_ctrl  = .find_tagged_col(C, "gene_expr_",  ctrl_tag),
    rtf_str  = if (any(grepl("^r_tf_", names(C)))) .find_tagged_col(C, "r_tf_", stress_tag) else NULL,
    rtf_ctrl = if (any(grepl("^r_tf_", names(C)))) .find_tagged_col(C, "r_tf_", ctrl_tag) else NULL
  )
}

# BEGIN EDIT: robust topic column detection
.find_topic_col <- function(x, topic_id = NULL) {
  nms <- names(x)
  if (length(nms) == 0) return(NA_character_)

  # Exact preferred names first
  preferred <- c("topic", "main_topic")
  hit <- intersect(preferred, nms)
  if (length(hit)) {
    return(hit[1])
  }

  # Any column whose name includes "topic" (case-insensitive)
  cand <- grep("topic", nms, value = TRUE, ignore.case = TRUE)
  # also allow plural forms
  cand <- unique(c(cand, grep("topics", nms, value = TRUE, ignore.case = TRUE)))

  # Also common LDA patterns
  cand <- unique(c(cand, grep("^lda_.*_topic(s)?$", nms, value = TRUE, ignore.case = TRUE)))

  if (!length(cand)) return(NA_character_)

  # If topic_id supplied, prefer a candidate column that plausibly contains it
  if (!is.null(topic_id) && length(topic_id) == 1 && is.finite(topic_id)) {
    tid <- as.integer(topic_id)
    for (cn in cand) {
      v <- x[[cn]]
      # Accept integer, numeric or character with delimited ints
      if (is.integer(v) || is.numeric(v)) {
        if (any(is.finite(v) & as.integer(v) == tid, na.rm = TRUE)) return(cn)
      } else if (is.character(v)) {
        # split on ; , | space
        toks <- unlist(strsplit(paste0(v, collapse = ";"), "[;,| ]+", perl = TRUE), use.names = FALSE)
        toks <- suppressWarnings(as.integer(toks))
        if (any(is.finite(toks) & toks == tid, na.rm = TRUE)) return(cn)
      }
    }
  }

  # Otherwise return the first reasonable candidate
  cand[1]
}
# END EDIT


# ------------------------------ Column detection for comp CSV ----------------
.detect_cols <- function(C, cond1_tag, cond2_tag) {
  # cond1_tag = control (ctrl), cond2_tag = stress (str)
  # Heuristics for TF/gene/peak columns
  tf_col <- if ("TF" %in% names(C)) "TF" else if ("tf" %in% names(C)) "tf" else stop("No TF column.")
  gene_col <- if ("gene_key" %in% names(C)) "gene_key" else if ("gene" %in% names(C)) "gene" else stop("No gene column.")
  peak_col <- if ("peak_id" %in% names(C)) "peak_id" else if ("peak" %in% names(C)) "peak" else stop("No peak column.")

  topic_col <- .find_topic_col(C)
  if (is.na(topic_col)) stop("No topic column found in comparison CSV.")

  # Scores/signs/active flags per side
  sc_ctrl <- if (any(grepl("^link_score_", names(C)))) .find_tagged_col(C, "link_score_", cond1_tag) else "score_ctrl"
  sc_str  <- if (any(grepl("^link_score_", names(C)))) .find_tagged_col(C, "link_score_", cond2_tag) else "score_str"

  sg_ctrl <- if ("link_sign_ctrl" %in% names(C)) "link_sign_ctrl" else if ("sign_ctrl" %in% names(C)) "sign_ctrl" else "sign_ctrl"
  sg_str  <- if ("link_sign_str"  %in% names(C)) "link_sign_str"  else if ("sign_str"  %in% names(C)) "sign_str"  else "sign_str"

  ac_ctrl <- if ("active_ctrl" %in% names(C)) "active_ctrl" else NULL
  ac_str  <- if ("active_str"  %in% names(C)) "active_str"  else NULL

  tf_ctrl <- if (any(grepl("^tf_expr_", names(C)))) .find_tagged_col(C, "tf_expr_", cond1_tag) else NULL
  tf_str  <- if (any(grepl("^tf_expr_", names(C)))) .find_tagged_col(C, "tf_expr_", cond2_tag) else NULL
  ge_ctrl <- if (any(grepl("^gene_expr_", names(C)))) .find_tagged_col(C, "gene_expr_", cond1_tag) else NULL
  ge_str  <- if (any(grepl("^gene_expr_", names(C)))) .find_tagged_col(C, "gene_expr_", cond2_tag) else NULL

  rtf_ctrl <- if (any(grepl("^r_tf_", names(C)))) .find_tagged_col(C, "r_tf_", cond1_tag) else NULL
  rtf_str  <- if (any(grepl("^r_tf_", names(C)))) .find_tagged_col(C, "r_tf_", cond2_tag) else NULL

  list(tf = tf_col, gene = gene_col, peak = peak_col, topic = topic_col,
       score_ctrl = sc_ctrl, score_str = sc_str,
       sign_ctrl = sg_ctrl, sign_str = sg_str,
       active_ctrl = ac_ctrl, active_str = ac_str,
       tf_expr_ctrl = tf_ctrl, tf_expr_str = tf_str,
       gene_expr_ctrl = ge_ctrl, gene_expr_str = ge_str,
       r_tf_ctrl = rtf_ctrl, r_tf_str = rtf_str)
}

# ------------------------------ Public Entry (Simple) ------------------------
#' Render Δ-topic single-panel network (simple; tolerant of suffixed columns)
#' @param comp_csv path to annotated comparison CSV
#' @param summary_csv path to matching topic summary CSV
#' @param topic_id integer topic
#' @param out_html optional html path; auto if NULL
#' @param top_n_tfs_per_topic integer whitelist from summary
#' @param edge_filter_min numeric; min |Δ| or per-side score to keep
#' @param edge_filter_on "either","stress","control","both"
#' @param gene_fc_thresh numeric (fold-change)
#' @param de_reference "str_over_ctrl" or "ctrl_over_str"
#' @param verbose logical
#' @return out_html invisibly
#' @export
render_link_network_delta_topic_simple <- function(
    comp_csv, summary_csv, topic_id,
    out_html = NULL,
    top_n_tfs_per_topic = 20,
    edge_filter_min = 2,
    edge_filter_on = c("either","stress","control","both"),
    gene_fc_thresh = 1.5,
    de_reference = c("str_over_ctrl","ctrl_over_str"),
    verbose = TRUE
){
  edge_filter_on <- match.arg(edge_filter_on)
  de_reference   <- match.arg(de_reference)
  stopifnot(file.exists(comp_csv), file.exists(summary_csv))

  .llog("[tfhubΔ] Reading: ", comp_csv,  verbose = verbose)
  C <- readr::read_csv(comp_csv, show_col_types = FALSE)
  .llog("[tfhubΔ] Reading: ", summary_csv, verbose = verbose)
  S <- readr::read_csv(summary_csv, show_col_types = FALSE)

  tags <- .guess_pair_from_path(comp_csv)
  stress_tag <- tags$stress_tag
  ctrl_tag   <- tags$ctrl_tag
  .llog(sprintf("[tfhubΔ] Δ-topic: stress='%s', control='%s' (Δ = stress − control).",
                stress_tag, ctrl_tag), verbose = verbose)

  comp_topic_col <- .find_topic_col(C)
  sum_topic_col  <- .find_topic_col(S)
  if (is.na(comp_topic_col)) stop("No topic column found in comparison CSV.")
  if (is.na(sum_topic_col))  stop("No topic column found in summary CSV.")

  if (comp_topic_col != "topic") C$topic <- suppressWarnings(as.integer(C[[comp_topic_col]]))
  if (sum_topic_col  != "topic") S$topic <- suppressWarnings(as.integer(S[[sum_topic_col]]))

  # TF whitelist for the topic
  S_sub <- S[S$topic == as.integer(topic_id), , drop = FALSE]
  tf_col_sum <- if ("TF" %in% names(S_sub)) "TF" else if ("tf" %in% names(S_sub)) "tf" else NA_character_
  if (is.na(tf_col_sum)) stop("Summary lacks TF column ('TF' or 'tf').")
  tf_keep <- unique(S_sub[[tf_col_sum]])
  if (length(tf_keep) > top_n_tfs_per_topic) tf_keep <- tf_keep[seq_len(top_n_tfs_per_topic)]
  .llog(sprintf("[tfhubΔ] Topic %s: keeping up to %d TF(s).", topic_id, length(tf_keep)), verbose = verbose)

  # Subset comparison to topic + TF whitelist
  Csub <- C[C$topic == as.integer(topic_id), , drop = FALSE]
  tf_col_comp <- if ("tf" %in% names(Csub)) "tf" else if ("TF" %in% names(Csub)) "TF" else NA_character_
  if (is.na(tf_col_comp)) stop("Comparison lacks TF column ('tf' or 'TF').")
  Csub <- Csub[Csub[[tf_col_comp]] %in% tf_keep, , drop = FALSE]

  # |Δ link| (if absent, compute from suffixed cols)
  if ("delta_link_score" %in% names(Csub)) {
    abs_d <- abs(suppressWarnings(as.numeric(Csub$delta_link_score)))
  } else {
    sc_str  <- .find_tagged_col(Csub, "link_score_", stress_tag)
    sc_ctrl <- .find_tagged_col(Csub, "link_score_", ctrl_tag)
    abs_d <- abs(
      suppressWarnings(as.numeric(Csub[[sc_str]])) -
        suppressWarnings(as.numeric(Csub[[sc_ctrl]]))
    )
  }
  keep <- is.finite(abs_d) & abs_d >= edge_filter_min
  Csub <- Csub[keep, , drop = FALSE]
  if (!nrow(Csub)) stop("No edges left after edge_filter_min=", edge_filter_min)

  tmp_dir <- file.path(tempdir(), "episcope_norm")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  comp_norm <- file.path(tmp_dir, basename(comp_csv))
  sum_norm  <- file.path(tmp_dir, basename(summary_csv))
  readr::write_csv(Csub, comp_norm)
  readr::write_csv(S,    sum_norm)

  out <- render_link_network_delta_topic(
    comp_csv = comp_norm,
    summary_csv = sum_norm,
    topic_id = topic_id,
    out_html = out_html,
    top_n_tfs_per_topic = top_n_tfs_per_topic,
    edge_filter_min = edge_filter_min,
    edge_filter_on  = edge_filter_on,
    gene_fc_thresh  = gene_fc_thresh,
    de_reference    = de_reference,
    verbose = verbose,
    topic_rank = if ("topic_rank" %in% names(S_sub))
      unique(suppressWarnings(as.integer(S_sub$topic_rank)))[1] else NA_integer_
  )
  invisible(out)
}

# ------------------------------ Core builder ---------------------------------
# BEGIN EDIT: rebuild .build_nodes_edges_from_comp with graceful topic handling
.build_nodes_edges_from_comp <- function(comp_df, summary_df, topic_id,
                                         top_n_tfs_per_topic = 20,
                                         edge_filter_min = 2,
                                         edge_filter_on = c("either","stress","control","both"),
                                         gene_fc_thresh = 1.5,
                                         tf_color_low = .col_tf_fill_low,
                                         tf_color_high = .col_tf_fill_high,
                                         tf_color_gamma = 1,
                                         cond1_tag = NULL,
                                         cond2_tag = NULL,
                                         de_reference = c("ctrl_over_str","str_over_ctrl"),
                                         verbose = TRUE) {
  edge_filter_on <- match.arg(edge_filter_on)
  de_reference   <- match.arg(de_reference)

  # detect columns (this will now call the new .find_topic_col)
  cols <- .detect_cols(comp_df, cond1_tag = cond1_tag, cond2_tag = cond2_tag)
  tf_col   <- cols$tf; gene_col <- cols$gene; peak_col <- cols$peak
  sc_ctrl  <- cols$score_ctrl; sc_str <- cols$score_str
  sg_ctrl  <- cols$sign_ctrl;  sg_str <- cols$sign_str
  ac_ctrl  <- cols$active_ctrl; ac_str <- cols$active_str
  tfc      <- cols$tf_expr_ctrl; tfs <- cols$tf_expr_str
  gxc      <- cols$gene_expr_ctrl; gxs <- cols$gene_expr_str

  # Find topic column in summary (required for TF whitelist)
  if (!("TF" %in% names(summary_df))) cli::cli_abort("summary CSV lacks TF column.")
  topic_col_in_S <- .find_topic_col(summary_df, topic_id = topic_id)
  if (is.na(topic_col_in_S)) cli::cli_abort("summary CSV lacks a usable topic column.")

  S_topic <- summary_df[summary_df[[topic_col_in_S]] == as.integer(topic_id), , drop = FALSE]
  if ("tf_delta_sum_abs" %in% names(S_topic)) {
    ord <- order(-suppressWarnings(as.numeric(S_topic$tf_delta_sum_abs)))
    S_topic <- S_topic[ord, , drop = FALSE]
  }
  if (nrow(S_topic) > top_n_tfs_per_topic) S_topic <- S_topic[seq_len(top_n_tfs_per_topic), , drop = FALSE]
  tf_whitelist <- unique(S_topic[["TF"]])
  if (!length(tf_whitelist)) cli::cli_abort("No TFs for topic {topic_id} found in summary CSV.")

  # --- Robust topic filtering on comp_df ---
  topic_col <- .find_topic_col(comp_df, topic_id = topic_id)

  if (!is.na(topic_col)) {
    # Handle both single-valued and delimited topic columns
    v <- comp_df[[topic_col]]
    use_rows <- rep(FALSE, nrow(comp_df))

    if (is.numeric(v) || is.integer(v)) {
      use_rows <- as.integer(v) == as.integer(topic_id)
      use_rows[is.na(use_rows)] <- FALSE
    } else {
      # character or factor: allow "7" or "6;7;12" etc.
      vv <- as.character(v)
      # Fast path: exact single-value match
      use_rows <- vv == as.character(topic_id)
      use_rows[is.na(use_rows)] <- FALSE
      # If nothing matched, try delimited
      if (!any(use_rows)) {
        delim_pat <- "(^|[;,| ]+)%s([;,| ]+|$)"
        pat <- sprintf(delim_pat, gsub("([.^$|()*+?{\\[\\]\\\\])", "\\\\\\1", as.character(topic_id)))
        use_rows <- grepl(pat, vv, perl = TRUE)
        use_rows[is.na(use_rows)] <- FALSE
      }
    }

    comp_exp <- comp_df[use_rows, , drop = FALSE]
    if (!nrow(comp_exp)) {
      cli::cli_warn("No rows matched topic {topic_id} in column '{topic_col}'. Using all rows as a fallback.")
      comp_exp <- comp_df
    }
  } else {
    cli::cli_warn("No topic column found in comparison CSV. Using all rows (no topic filtering).")
    comp_exp <- comp_df
  }
  # ------------------------------------------

  # build edge table
  TF        <- comp_exp[[tf_col]]
  gene_key  <- comp_exp[[gene_col]]
  peak_id   <- comp_exp[[peak_col]]
  score_ctrl <- suppressWarnings(as.numeric(comp_exp[[sc_ctrl]]))
  score_str  <- suppressWarnings(as.numeric(comp_exp[[sc_str]]))
  sign_ctrl  <- .as_sign(comp_exp[[sg_ctrl]]) ; sign_ctrl[!is.finite(sign_ctrl)] <- 0L
  sign_str   <- .as_sign(comp_exp[[sg_str]])  ; sign_str[!is.finite(sign_str)]  <- 0L
  active_ctrl <- if (!is.null(ac_ctrl)) .as_bool(comp_exp[[ac_ctrl]]) else rep(NA, nrow(comp_exp))
  active_str  <- if (!is.null(ac_str))  .as_bool(comp_exp[[ac_str]])  else rep(NA, nrow(comp_exp))

  r_tf_ctrl <- if (!is.null(cols$r_tf_ctrl)) suppressWarnings(as.numeric(comp_exp[[cols$r_tf_ctrl]])) else rep(NA_real_, nrow(comp_exp))
  r_tf_str  <- if (!is.null(cols$r_tf_str )) suppressWarnings(as.numeric(comp_exp[[cols$r_tf_str ]])) else rep(NA_real_, nrow(comp_exp))
  r_any_col <- grep("^r_tf($|_)", names(comp_exp), value = TRUE)
  r_any_col <- setdiff(r_any_col, c(cols$r_tf_ctrl, cols$r_tf_str))
  r_tf_any  <- if (length(r_any_col)) suppressWarnings(as.numeric(comp_exp[[r_any_col[1]]])) else rep(NA_real_, nrow(comp_exp))
  r_tf      <- ifelse(is.finite(r_tf_ctrl), r_tf_ctrl, ifelse(is.finite(r_tf_str), r_tf_str, r_tf_any))

  tf_expr_ctrl   <- if (!is.null(tfc)) suppressWarnings(as.numeric(comp_exp[[tfc]])) else rep(NA_real_, nrow(comp_exp))
  tf_expr_str    <- if (!is.null(tfs)) suppressWarnings(as.numeric(comp_exp[[tfs]])) else rep(NA_real_, nrow(comp_exp))
  gene_expr_ctrl <- if (!is.null(gxc)) suppressWarnings(as.numeric(comp_exp[[gxc]])) else rep(NA_real_, nrow(comp_exp))
  gene_expr_str  <- if (!is.null(gxs)) suppressWarnings(as.numeric(comp_exp[[gxs]])) else rep(NA_real_, nrow(comp_exp))

  df <- data.frame(TF = TF, gene_key = gene_key, peak_id = peak_id,
                   score_ctrl = score_ctrl, score_str = score_str,
                   sign_ctrl = sign_ctrl, sign_str = sign_str,
                   active_ctrl = active_ctrl, active_str = active_str,
                   r_tf_ctrl = r_tf_ctrl, r_tf_str = r_tf_str, r_tf_any = r_tf_any, r_tf = r_tf,
                   tf_expr_ctrl = tf_expr_ctrl, tf_expr_str = tf_expr_str,
                   gene_expr_ctrl = gene_expr_ctrl, gene_expr_str = gene_expr_str,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$TF) & !is.na(df$gene_key) & !is.na(df$peak_id), , drop = FALSE]
  df <- df[!duplicated(df[, c("TF","gene_key","peak_id")]), , drop = FALSE]

  # TF whitelist
  df <- df[df$TF %in% tf_whitelist, , drop = FALSE]
  if (!nrow(df)) cli::cli_abort("No (TF,gene,peak) edges remain after TF whitelist.")

  # edge-level filters
  pass_ctrl <- .compute_pass_flags(df, edge_filter_min, "score_ctrl", "active_ctrl")
  pass_str  <- .compute_pass_flags(df, edge_filter_min, "score_str",  "active_str")
  # <<< NEW: retain flags on df so downstream code (e.g., .aggregate_peak_gene) can use them
  df$pass_ctrl <- pass_ctrl
  df$pass_str  <- pass_str
  keep <- switch(edge_filter_on,
                 control = pass_ctrl,
                 stress  = pass_str,
                 both    = pass_ctrl & pass_str,
                 either  = pass_ctrl | pass_str)
  keep[is.na(keep)] <- FALSE
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) cli::cli_abort("All edges filtered by 'edge_filter_on' and threshold.")

  # Node table (TFs + genes) and expression summaries
  node_ids <- sort(unique(c(df$TF, df$gene_key)))
  type <- ifelse(node_ids %in% unique(df$TF), "TF", "gene")
  nodes <- data.frame(id = node_ids, type = type, stringsAsFactors = FALSE)

  tf_ctrl_tab <- unique(df[, c("TF","tf_expr_ctrl")]); colnames(tf_ctrl_tab) <- c("id","tf_expr_ctrl")
  tf_str_tab  <- unique(df[, c("TF","tf_expr_str")]);  colnames(tf_str_tab)  <- c("id","tf_expr_str")
  g_ctrl_tab  <- unique(df[, c("gene_key","gene_expr_ctrl")]); colnames(g_ctrl_tab) <- c("id","gene_expr_ctrl")
  g_str_tab   <- unique(df[, c("gene_key","gene_expr_str")]);  colnames(g_str_tab)  <- c("id","gene_expr_str")

  nodes <- merge(nodes, tf_ctrl_tab, by = "id", all.x = TRUE)
  nodes <- merge(nodes, tf_str_tab,  by = "id", all.x = TRUE)
  nodes <- merge(nodes, g_ctrl_tab,  by = "id", all.x = TRUE)
  nodes <- merge(nodes, g_str_tab,   by = "id", all.x = TRUE)

  deg_ctrl <- stats::aggregate(df$TF[pass_ctrl], by = list(id = df$TF[pass_ctrl]), FUN = length)
  colnames(deg_ctrl)[2] <- "n_active_ctrl"
  deg_str  <- stats::aggregate(df$TF[pass_str],  by = list(id = df$TF[pass_str]),  FUN = length)
  colnames(deg_str)[2] <- "n_active_str"
  nodes <- merge(nodes, deg_ctrl, by = "id", all.x = TRUE)
  nodes <- merge(nodes, deg_str,  by = "id", all.x = TRUE)
  nodes$n_active_ctrl[is.na(nodes$n_active_ctrl)] <- 0L
  nodes$n_active_str [is.na(nodes$n_active_str )] <- 0L

  nodes$node_raw_ctrl <- ifelse(nodes$type == "TF", nodes$tf_expr_ctrl, nodes$gene_expr_ctrl)
  nodes$node_raw_str  <- ifelse(nodes$type == "TF", nodes$tf_expr_str,  nodes$gene_expr_str)
  nodes$node_z_ctrl   <- robust_z(nodes$node_raw_ctrl)
  nodes$node_z_str    <- robust_z(nodes$node_raw_str)

  eps <- 1e-9; thr <- log2(gene_fc_thresh)
  l2fc_gene <- if (de_reference == "ctrl_over_str") {
    log2((nodes$gene_expr_ctrl + eps) / (nodes$gene_expr_str + eps))
  } else {
    log2((nodes$gene_expr_str + eps) / (nodes$gene_expr_ctrl + eps))
  }
  l2fc_tf <- if (de_reference == "ctrl_over_str") {
    log2((nodes$tf_expr_ctrl + eps) / (nodes$tf_expr_str + eps))
  } else {
    log2((nodes$tf_expr_str + eps) / (nodes$tf_expr_ctrl + eps))
  }
  nodes$l2fc <- ifelse(nodes$type == "TF", l2fc_tf, l2fc_gene)

  list(nodes = nodes, ed = df)
}
# END EDIT

# ------------------------------ Layout & Peak pseudonodes --------------------
.layout_topic_circle <- function(nodes, ed) {
  # Place TFs evenly on a circle; genes on a larger concentric ring near their TFs
  tfs   <- nodes$id[nodes$type == "TF"]
  genes <- nodes$id[nodes$type == "gene"]

  n_tf <- length(tfs); n_gene <- length(genes)
  r_tf <- 200; r_gene <- 320

  theta <- if (n_tf > 0) seq(0, 2*pi, length.out = n_tf + 1)[- (n_tf + 1)] else numeric(0)
  tf_xy <- if (n_tf > 0) data.frame(id = tfs, x = r_tf * cos(theta), y = r_tf * sin(theta), stringsAsFactors = FALSE) else data.frame(id=character(),x=numeric(),y=numeric())

  # genes pulled towards the mean angle of their TF neighbors
  gene_theta <- rep(0, n_gene)
  if (n_gene > 0 && n_tf > 0) {
    # compute nearest TF angle for each gene by frequency in edges
    pairs <- unique(ed[, c("TF","gene_key")])
    pairs <- pairs[pairs$TF %in% tfs & pairs$gene_key %in% genes, , drop = FALSE]
    # map TF -> angle
    tf_angle <- stats::setNames(theta, tfs)
    for (i in seq_len(n_gene)) {
      g <- genes[i]
      tset <- pairs$TF[pairs$gene_key == g]
      if (length(tset)) {
        gene_theta[i] <- mean(tf_angle[tset])
      } else {
        gene_theta[i] <- 2*pi*i/n_gene
      }
    }
  } else if (n_gene > 0) {
    gene_theta <- 2*pi*seq_len(n_gene)/n_gene
  }
  gene_xy <- if (n_gene > 0) data.frame(id = genes, x = r_gene * cos(gene_theta), y = r_gene * sin(gene_theta), stringsAsFactors = FALSE) else data.frame(id=character(),x=numeric(),y=numeric())

  coords <- rbind(tf_xy, gene_xy)
  coords
}

.add_peak_pseudonodes <- function(nodes, coords, ed) {
  # shared peaks: PEAK:<id> appear once if same peak connects >=2 TFs or >=2 genes
  pk_id <- unique(ed$peak_id)
  if (!length(pk_id)) return(list(nodes = nodes, coords = coords, peak_ids = character(), bridge_pairs = data.frame()))
  peak_ids <- paste0("PEAK:", pk_id)

  # place each peak midway between its TF and gene centroid
  add_nodes <- data.frame(
    id = peak_ids,
    type = "peak",
    stringsAsFactors = FALSE
  )
  # rough placement: mean of connected gene & TF coords
  tf_xy <- coords[match(unique(ed$TF), coords$id), c("x","y")]
  gene_xy <- coords[match(unique(ed$gene_key), coords$id), c("x","y")]
  cx <- mean(c(tf_xy$x, gene_xy$x), na.rm = TRUE)
  cy <- mean(c(tf_xy$y, gene_xy$y), na.rm = TRUE)
  add_coords <- data.frame(id = peak_ids, x = cx, y = cy, stringsAsFactors = FALSE)

  nodes2  <- rbind(nodes, add_nodes)
  coords2 <- rbind(coords, add_coords)
  list(nodes = nodes2, coords = coords2, peak_ids = peak_ids, bridge_pairs = data.frame())
}

.add_unique_peak_pseudonodes <- function(nodes, coords, ed, w_gene = 0.08, jitter_sd_px = 2) {
  # Treat all peaks as unique for simplicity (UPEAK:)
  pk_id <- unique(ed$peak_id)
  if (!length(pk_id)) return(list(nodes = nodes, coords = coords, upeak_ids = character(), bridge_pairs = data.frame()))
  upeak_ids <- paste0("UPEAK:", pk_id)

  # place near their target genes
  gid <- unique(ed$gene_key)
  gxy <- coords[match(gid, coords$id), c("x","y")]
  # naive assignment: cycle genes
  base_xy <- if (nrow(gxy)) gxy else data.frame(x = 0, y = 0)
  rep_xy <- base_xy[rep(seq_len(nrow(base_xy)), length.out = length(upeak_ids)), , drop = FALSE]
  set.seed(1)
  add_coords <- data.frame(
    id = upeak_ids,
    x = rep_xy$x * (1 - w_gene) + stats::rnorm(nrow(rep_xy), 0, jitter_sd_px),
    y = rep_xy$y * (1 - w_gene) + stats::rnorm(nrow(rep_xy), 0, jitter_sd_px),
    stringsAsFactors = FALSE
  )
  add_nodes <- data.frame(id = upeak_ids, type = "peak", stringsAsFactors = FALSE)
  nodes2  <- rbind(nodes, add_nodes)
  coords2 <- rbind(coords, add_coords)
  list(nodes = nodes2, coords = coords2, upeak_ids = upeak_ids, bridge_pairs = data.frame())
}

.bump_peak_nodes_radially <- function(coords, nodes, push_px = 22) {
  pk <- coords$id[coords$id %in% nodes$id[which(nodes$type == "peak")]]
  if (!length(pk)) return(coords)
  for (i in match(pk, coords$id)) {
    vx <- coords$x[i]; vy <- coords$y[i]
    r <- sqrt(vx*vx + vy*vy) + push_px
    ang <- atan2(vy, vx)
    coords$x[i] <- r * cos(ang); coords$y[i] <- r * sin(ang)
  }
  coords
}

# ------------------------------ Edge styling ---------------------------------
# Convert multi-edges into "fanned" edges by adding a small curvature (vis option: smooth)
.fan_multi_edges <- function(df) {
  if (!nrow(df)) return(df)
  key <- paste(df$from, df$to)
  tab <- table(key)
  multi <- names(tab[tab > 1])
  if (!length(multi)) return(df)
  df$physics <- FALSE
  df$smooth  <- TRUE
  df
}

.force_edge_cols <- function(df) {
  base <- data.frame(id = character(0), from = character(0), to = character(0),
                     arrows = "to", color = "#888888", width = 1,
                     dashes = FALSE, smooth = FALSE, title = "", stringsAsFactors = FALSE)
  miss <- setdiff(colnames(base), colnames(df))
  for (m in miss) df[[m]] <- base[[m]]
  df
}

# Topic Δ direct edges TF→gene
.style_edges_delta_topic <- function(ed, cap_abs = 1) {
  if (!nrow(ed)) return(data.frame())
  d <- suppressWarnings(as.numeric(ed$score_str) - as.numeric(ed$score_ctrl))
  absd <- abs(d)
  w <- ifelse(absd > 8, 4, ifelse(absd > 4, 3, ifelse(absd > 2, 2, ifelse(absd > 0, 1, 1))))
  col <- ifelse(d >= 0,
                grDevices::rgb(216,107,94, maxColorValue = 255),
                grDevices::rgb(111,158,212, maxColorValue = 255))
  data.frame(
    id = paste("dir", ed$TF, ed$gene_key, sep="|"),
    from = ed$TF, to = ed$gene_key,
    arrows = "to", color = col, width = w,
    dashes = FALSE, smooth = FALSE,
    title = paste0("Δ=", round(d, 3)),
    stringsAsFactors = FALSE
  )
}

# TF→PEAK gray dual (ctrl/str) collapsed to a single dotted gray hint
.style_tf_peak_edges_gray_dual <- function(ed) {
  if (!nrow(ed)) return(data.frame())
  data.frame(
    id = paste("tfpeak", ed$TF, ed$peak_id, sep="|"),
    from = ed$TF, to = paste0("PEAK:", ed$peak_id),
    arrows = "to", color = "#9E9E9E", width = 1,
    dashes = TRUE, smooth = FALSE,
    title = "TF→peak (dual)",
    stringsAsFactors = FALSE
  )
}

.style_tf_upeak_edges_gray_dual <- function(ed) {
  if (!nrow(ed)) return(data.frame())
  data.frame(
    id = paste("tfu", ed$TF, ed$peak_id, sep="|"),
    from = ed$TF, to = paste0("UPEAK:", ed$peak_id),
    arrows = "to", color = "#B0B0B0", width = 1,
    dashes = TRUE, smooth = FALSE,
    title = "TF→unique peak (dual)",
    stringsAsFactors = FALSE
  )
}

# Aggregate peak→gene per side, then return Δ edges
# BEGIN EDIT: robust peak→gene aggregator used by Δ-topic plot
.aggregate_peak_gene <- function(ed) {
  # Expect columns: peak_id, gene_key, score_ctrl, score_str, sign_ctrl, sign_str,
  # and (now) pass_ctrl, pass_str. Be tolerant if sign/flags are missing.
  has <- function(nm) nm %in% names(ed)

  # Safe vectors
  sc_ctrl <- suppressWarnings(as.numeric(ed$score_ctrl))
  sc_str  <- suppressWarnings(as.numeric(ed$score_str))
  sg_ctrl <- if (has("sign_ctrl")) .as_sign(ed$sign_ctrl) else rep(1L, nrow(ed))
  sg_str  <- if (has("sign_str"))  .as_sign(ed$sign_str)  else rep(1L, nrow(ed))
  pc      <- if (has("pass_ctrl")) ed$pass_ctrl else rep(TRUE, nrow(ed))
  ps      <- if (has("pass_str"))  ed$pass_str  else rep(TRUE, nrow(ed))

  # Only count edges that pass on each side
  w_ctrl  <- ifelse(isTRUE(pc), 1, 0)
  w_str   <- ifelse(isTRUE(ps),  1, 0)

  # Signed scores by side
  signed_ctrl <- (sc_ctrl * sg_ctrl) * w_ctrl
  signed_str  <- (sc_str  * sg_str)  * w_str

  df <- data.frame(
    peak_id = ed$peak_id,
    gene_key = ed$gene_key,
    score_ctrl = signed_ctrl,
    score_str  = signed_str,
    n_edges_ctrl = as.integer(w_ctrl),
    n_edges_str  = as.integer(w_str),
    stringsAsFactors = FALSE
  )

  # Aggregate per peak→gene (sum of signed, counts of passing edges)
  agg <- df |>
    dplyr::group_by(.data$peak_id, .data$gene_key) |>
    dplyr::summarise(
      score_ctrl = sum(.data$score_ctrl, na.rm = TRUE),
      score_str  = sum(.data$score_str,  na.rm = TRUE),
      n_edges_ctrl = sum(.data$n_edges_ctrl, na.rm = TRUE),
      n_edges_str  = sum(.data$n_edges_str,  na.rm = TRUE),
      .groups = "drop"
    )

  # Keep column names expected by styling helpers
  agg
}
# END EDIT

.style_edges_delta_pg <- function(pg_df, cap_abs = 1) {
  if (!nrow(pg_df)) return(data.frame())
  d <- suppressWarnings(as.numeric(pg_df$score_str) - as.numeric(pg_df$score_ctrl))
  absd <- abs(d)
  w <- ifelse(absd > 8, 4, ifelse(absd > 4, 3, ifelse(absd > 2, 2, ifelse(absd > 0, 1, 1))))
  col <- ifelse(d >= 0,
                grDevices::rgb(216,107,94, maxColorValue = 255),
                grDevices::rgb(111,158,212, maxColorValue = 255))
  data.frame(
    id    = paste("pg", pg_df$peak_id, pg_df$gene_key, sep="|"),
    from  = paste0("PEAK:", pg_df$peak_id),
    to    = pg_df$gene_key,
    arrows = "to", color = col, width = w,
    dashes = FALSE, smooth = FALSE,
    title = paste0("Δ=", round(d, 3)),
    stringsAsFactors = FALSE
  )
}

.style_edges_delta_upg <- function(pg_df, cap_abs = 1) {
  if (!nrow(pg_df)) return(data.frame())
  d <- suppressWarnings(as.numeric(pg_df$score_str) - as.numeric(pg_df$score_ctrl))
  absd <- abs(d)
  w <- ifelse(absd > 8, 4, ifelse(absd > 4, 3, ifelse(absd > 2, 2, ifelse(absd > 0, 1, 1))))
  col <- ifelse(d >= 0,
                grDevices::rgb(216,107,94, maxColorValue = 255),
                grDevices::rgb(111,158,212, maxColorValue = 255))
  data.frame(
    id    = paste("upg", pg_df$peak_id, pg_df$gene_key, sep="|"),
    from  = paste0("UPEAK:", pg_df$peak_id),
    to    = pg_df$gene_key,
    arrows = "to", color = col, width = w,
    dashes = FALSE, smooth = FALSE,
    title = paste0("Δ=", round(d, 3)),
    stringsAsFactors = FALSE
  )
}

.drop_direct_edges_covered_by_peaks <- function(dir_edges, edsub, peak_nodes) {
  if (!nrow(dir_edges)) return(dir_edges)
  has_peak <- paste(edsub$TF, edsub$gene_key)  # presence implies covered
  keep <- !(paste(dir_edges$from, dir_edges$to) %in% unique(has_peak))
  dir_edges[keep, , drop = FALSE]
}

.neutral_bridge_edges <- function(pairs) {
  # optional thin gray bridges from PEAK/UPEAK to gene
  if (!nrow(pairs)) return(data.frame())
  data.frame(
    id = paste("br", pairs$from, pairs$to, sep="|"),
    from = pairs$from, to = pairs$to,
    arrows = "to", color = "#BDBDBD", width = 1,
    dashes = FALSE, smooth = FALSE, title = "bridge",
    stringsAsFactors = FALSE
  )
}

# ------------------------------ vis helpers (no-ops/overlays) ----------------
.attach_vis_jumpers <- function(widget) widget
.layer_stack_nodes_and_edges <- function(widget) widget
.attach_draw_peak_ring <- function(widget, r_ring) widget
.attach_constrain_peak_ring <- function(widget, r_ring) widget

# ------------------------------ Main Renderer --------------------------------
#' Δ-topic (stress − control) single-panel network
#' See long description in the source code header.
#' @inheritParams render_link_network_delta_topic_simple
#' @param topic_rank optional integer rank to include in title
#' @return out_html invisibly
#' @export
render_link_network_delta_topic <- function(
    comp_csv, summary_csv, topic_id,
    out_html = NULL,
    top_n_tfs_per_topic = 20,
    edge_filter_min = 2,
    edge_filter_on = c("either","stress","control","both"),
    gene_fc_thresh = 1.5,
    de_reference = c("str_over_ctrl","ctrl_over_str"),
    verbose = TRUE,
    topic_rank = NULL
){
  edge_filter_on <- match.arg(edge_filter_on)
  de_reference   <- match.arg(de_reference)

  stopifnot(file.exists(comp_csv), file.exists(summary_csv))
  .llog("Reading comparison CSV: ", comp_csv, verbose = verbose)
  C <- readr::read_csv(comp_csv, show_col_types = FALSE)
  .llog("Reading summary CSV: ", summary_csv, verbose = verbose)
  S <- readr::read_csv(summary_csv, show_col_types = FALSE)

  base <- basename(comp_csv)
  hdr  <- sub("_delta_links.*$", "", base)
  bits <- strsplit(hdr, "_vs_", fixed = TRUE)[[1]]
  stress_tag <- bits[[1]] %||% "stress"
  ctrl_tag   <- bits[[2]] %||% "control"
  .llog("Δ-topic: STRESS_vs_CONTROL: stress='", stress_tag,
        "', control='", ctrl_tag, "'. Δ = stress − control.", verbose = verbose)

  core <- .build_nodes_edges_from_comp(
    comp_df = C, summary_df = S, topic_id = topic_id,
    top_n_tfs_per_topic = top_n_tfs_per_topic,
    edge_filter_min = edge_filter_min,
    edge_filter_on  = edge_filter_on,
    cond1_tag = ctrl_tag,   # CONTROL → 'ctrl'
    cond2_tag = stress_tag, # STRESS  → 'str'
    de_reference = de_reference,
    verbose = verbose
  )

  nodes <- core$nodes
  ed    <- core$ed

  coords <- .layout_topic_circle(nodes, ed = ed[, c("TF","gene_key","peak_id")])
  pk  <- .add_peak_pseudonodes(nodes, coords, ed)       # PEAK:
  nodes  <- pk$nodes;  coords <- pk$coords
  upk <- .add_unique_peak_pseudonodes(nodes, coords, ed, w_gene = 0.08, jitter_sd_px = 2)  # UPEAK:
  nodes  <- upk$nodes; coords <- upk$coords
  coords <- .bump_peak_nodes_radially(coords, nodes, push_px = 22)

  # Node styling → visNetwork nodes data.frame
  tf_idx  <- nodes$type == "TF"
  tf_l2fc <- nodes$l2fc[tf_idx]
  max_abs <- suppressWarnings(max(abs(tf_l2fc), na.rm = TRUE)); if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
  tf_signed <- rep(NA_real_, nrow(nodes)); tf_signed[tf_idx] <- tf_l2fc / max_abs
  ramp_div <- function(u_signed){
    u_signed[!is.finite(u_signed)] <- 0
    w <- 0.5 + 0.5 * clamp(u_signed, -1, 1)
    grDevices::rgb(grDevices::colorRamp(c(.col_edge_neg, "#FFFFFF", .col_edge_pos))(w)/255)
  }

  gene_z_joint <- pmax(abs(nodes$node_z_ctrl), abs(nodes$node_z_str), na.rm = TRUE)
  gene_size <- {
    x <- pmax(0, suppressWarnings(as.numeric(gene_z_joint)))
    if (!any(is.finite(x)) || max(x, na.rm = TRUE) == 0) rep(mean(.gene_size_range), nrow(nodes)) else {
      s <- (x / max(x, na.rm = TRUE))^.gene_size_gamma
      .gene_size_range[1] + s * (.gene_size_range[2] - .gene_size_range[1])
    }
  }

  thr <- log2(gene_fc_thresh)
  border_col <- ifelse(!is.finite(nodes$l2fc) | abs(nodes$l2fc) < thr, .col_border_none,
                       ifelse(nodes$l2fc > 0, .col_border_up, .col_border_down))
  fc_abs <- 2^abs(nodes$l2fc)
  bw_gene <- ifelse(!is.finite(fc_abs), 1,
                    ifelse(fc_abs > 8, 4, ifelse(fc_abs > 4, 3, ifelse(fc_abs > 2, 2, 1))))
  border_width <- ifelse(nodes$type == "TF", 1, bw_gene)
  border_width[nodes$type == "peak"] <- 0.5

  pos <- coords[match(nodes$id, coords$id), c("x","y")]
  pos$x[!is.finite(pos$x)] <- 0; pos$y[!is.finite(pos$y)] <- 0

  x_tf <- pmax(nodes$tf_expr_ctrl, nodes$tf_expr_str, na.rm = TRUE)
  x_tf[!is.finite(x_tf)] <- 0
  if (!any(is.finite(x_tf)) || max(x_tf, na.rm = TRUE) == 0) {
    tf_font <- rep(mean(.tf_font_range), nrow(nodes))
  } else {
    s <- (x_tf / max(x_tf, na.rm = TRUE))^.tf_font_gamma
    tf_font <- .tf_font_range[1] + s * (.tf_font_range[2] - .tf_font_range[1])
  }

  nd <- data.frame(
    id    = nodes$id,
    label = ifelse(nodes$type == "TF", nodes$id, ""),
    x = pos$x, y = pos$y,
    color.background = ifelse(nodes$type == "TF", ramp_div(tf_signed), "#D9EAD3"),
    color.border     = border_col,
    borderWidth      = border_width,
    value = ifelse(nodes$type == "TF", 1,
                   ifelse(nodes$type == "peak", .peak_node_value, gene_size)),
    title = ifelse(
      nodes$type == "TF",
      sprintf("TF %s%s%s<br>log2FC(TF RNA) (ref=%s): %.3f<br>#active edges: ctrl=%s, str=%s",
              nodes$id,
              ifelse(is.finite(nodes$tf_expr_ctrl), sprintf("<br>expr_ctrl=%.3f", nodes$tf_expr_ctrl), ""),
              ifelse(is.finite(nodes$tf_expr_str ), sprintf("<br>expr_str =%.3f", nodes$tf_expr_str ), ""),
              de_reference,
              ifelse(is.finite(nodes$l2fc), nodes$l2fc, NA_real_),
              nodes$n_active_ctrl, nodes$n_active_str),
      ifelse(
        nodes$type == "peak",
        sprintf("Peak %s", gsub("^(UPEAK|PEAK):","", nodes$id)),
        sprintf("Gene %s%s%s",
                nodes$id,
                ifelse(is.finite(nodes$gene_expr_ctrl), sprintf("<br>expr_ctrl=%.3f", nodes$gene_expr_ctrl), ""),
                ifelse(is.finite(nodes$gene_expr_str ), sprintf("<br>expr_str =%.3f", nodes$gene_expr_str ), ""))
      )
    ),
    shape = ifelse(nodes$type == "TF", "box",
                   ifelse(nodes$type == "peak", "triangle", "dot")),
    physics = FALSE, fixed = FALSE,
    "font.color"       = "#2b2b2b",
    "font.size"        = ifelse(nodes$type == "TF", tf_font, 0),
    "font.face"        = "arial",
    "font.strokeWidth" = ifelse(nodes$type == "TF", .tf_label_halo_w, 0),
    "font.strokeColor" = ifelse(nodes$type == "TF", .tf_label_halo_c, NA),
    stringsAsFactors = FALSE
  )

  delta_all <- abs(suppressWarnings(as.numeric(ed$score_str) - as.numeric(ed$score_ctrl)))
  cap_joint <- stats::quantile(delta_all[delta_all > 0], 0.995, na.rm = TRUE, names = FALSE)
  if (!is.finite(cap_joint) || cap_joint <= 0) cap_joint <- suppressWarnings(max(delta_all, na.rm = TRUE))
  if (!is.finite(cap_joint) || cap_joint <= 0) cap_joint <- 1

  shared_ids <- gsub("^PEAK:",  "", pk$peak_ids)
  unique_ids <- gsub("^UPEAK:", "", upk$upeak_ids)

  ed_shared <- if (length(shared_ids)) ed[ed$peak_id %in% shared_ids, , drop = FALSE] else ed[0,]
  ed_unique <- if (length(unique_ids)) ed[ed$peak_id %in% unique_ids, , drop = FALSE] else ed[0,]

  e_tfpk_shared <- .style_tf_peak_edges_gray_dual(ed_shared)
  e_tfu         <- .style_tf_upeak_edges_gray_dual(ed_unique)

  pg_all   <- .aggregate_peak_gene(ed)
  pg_share <- if (length(shared_ids)) pg_all[paste0("PEAK:", pg_all$peak_id)  %in% pk$peak_ids, , drop = FALSE] else pg_all[0,]
  pg_uniq  <- if (length(unique_ids)) pg_all[paste0("UPEAK:", pg_all$peak_id) %in% upk$upeak_ids, , drop = FALSE] else pg_all[0,]

  pg_delta_all <- abs(suppressWarnings(as.numeric(pg_all$score_str) - as.numeric(pg_all$score_ctrl)))
  cap_pg <- stats::quantile(pg_delta_all[pg_delta_all > 0], 0.995, na.rm = TRUE, names = FALSE)
  if (!is.finite(cap_pg) || cap_pg <= 0) cap_pg <- suppressWarnings(max(pg_delta_all, na.rm = TRUE))
  if (!is.finite(cap_pg) || cap_pg <= 0) cap_pg <- 1

  e_pg_delta  <- .style_edges_delta_pg (pg_share, cap_abs = cap_pg)
  e_upg_delta <- .style_edges_delta_upg(pg_uniq,  cap_abs = cap_pg)

  e_delta_direct <- .style_edges_delta_topic(ed, cap_abs = cap_joint)
  e_delta_direct <- .drop_direct_edges_covered_by_peaks(
    dir_edges = e_delta_direct,
    edsub     = ed,
    peak_nodes = c(pk$peak_ids, upk$upeak_ids)
  )
  if (nrow(e_delta_direct)) {
    e_delta_direct$id     <- paste("dir", e_delta_direct$id, sep="|")
    e_delta_direct$dashes <- FALSE
    e_delta_direct$color  <- gsub(",\\s*[0-9.]+\\)$", ",0.35)", e_delta_direct$color)
  }

  edv <- rbind(
    .force_edge_cols(e_tfpk_shared),
    .force_edge_cols(e_tfu),
    .force_edge_cols(e_pg_delta),
    .force_edge_cols(e_upg_delta),
    .force_edge_cols(e_delta_direct)
  )
  edv <- .fan_multi_edges(edv)

  main_hdr <- sprintf(
    "%s | Topic %d%s – Δ network (%s − %s; edges per TF–gene–peak)",
    hdr, as.integer(topic_id),
    if (!is.null(topic_rank) && is.finite(topic_rank)) sprintf(" [rank %d]", as.integer(topic_rank)) else "",
    stress_tag, ctrl_tag
  )

  tf_xy <- nd[nd$shape == "box", c("x","y")]
  r_tf_med <- stats::median(sqrt(tf_xy$x^2 + tf_xy$y^2), na.rm = TRUE)
  r_ring   <- max(10, r_tf_med - .peak_ring_inner_offset_px)

  allowed_ids <- setdiff(
    unique(nd$id),
    unique(nd$id[grepl("^(PEAK|UPEAK):", nd$id) | nd$shape == "triangle"])
  )

  widget <- visNetwork::visNetwork(nd, edv, width = "100%", height = "700px", main = main_hdr) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                           nodesIdSelection = list(enabled = TRUE, useLabels = TRUE, values = allowed_ids)) %>%
    visNetwork::visEdges(scaling = list(min = 1, max = 10), smooth = FALSE) %>%
    visNetwork::visNodes(borderWidth = 1, font = list(bold = TRUE),
                         margin = list(top=6,right=10,bottom=6,left=10)) %>%
    visNetwork::visPhysics(enabled = FALSE) %>%
    visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
    visNetwork::visExport(type = "pdf",
                          name = sprintf("topic-%d-delta", as.integer(topic_id)),
                          label = "Save as PDF") %>%
    .attach_vis_jumpers() %>%
    .layer_stack_nodes_and_edges() %>%
    .attach_draw_peak_ring(r_ring) %>%
    .attach_constrain_peak_ring(r_ring)

  # BEGIN EDIT: wrap legend/title into the widget and save it
  widget_out <- htmlwidgets::prependContent(
    widget,
    htmltools::h3(main_hdr),
    htmltools::div(
      class = "grn-legend",
      htmltools::HTML(paste0(
        "<span>Legend:</span>",
        "<span><span class='grad'></span>TF (box) fill by signed log2FC(TF RNA) (ref=", de_reference, "): blue (−) → white (0) → red (+)</span>",
        "<span><span class='dot'></span>Gene (dot) size ∝ max |Z(expr)| across sides</span>",
        "<span>Edge color: red if Δ&gt;0, blue if Δ&lt;0 (Δ = ", stress_tag, " − ", ctrl_tag, ")</span>",
        "<span>Edge width ∝ |Δ link_score| (bins: &gt;8⇒4, &gt;4⇒3, &gt;2⇒2, &gt;0⇒1)</span>",
        "<span>Triangle = peak pseudonode (shared PEAK: or unique UPEAK:)</span>",
        "<span>Arrow tip: pointed = activation; blunt bar = repression (if sign provided)</span>"
      ))
    ),
    htmltools::tags$style(htmltools::HTML("
    body { font-family: Helvetica, Arial, sans-serif; margin: 16px; }
    .grn-legend { position: relative !important; float: none !important; clear: both;
      width: min(1200px, 96vw); margin: 12px auto 6px; box-sizing: border-box;
      display: flex; gap: 18px; align-items: center; justify-content: center;
      flex-wrap: wrap; font-size: 12px; color: #555; white-space: normal !important; }
    .grn-legend .grad  { display:inline-block; width:40px; height:14px; background:linear-gradient(90deg,#4575b4,#FFFFFF,#d73027); border:1px solid #BFBFBF; margin-right:6px;}
    .grn-legend .dot   { display:inline-block; width:12px; height:12px; background:#D9D9D9; border-radius:50%; border:1px solid #BFBFBF; margin-right:6px;}
  "))
  )

  if (is.null(out_html)) {
    file_base <- tools::file_path_sans_ext(basename(comp_csv))
    rank_tag  <- if (!is.null(topic_rank) && is.finite(topic_rank))
      sprintf("_rank-%02d", as.integer(topic_rank)) else ""
    out_html <- file.path(dirname(comp_csv),
                          sprintf("%s_topic-%d%s_subnetwork_delta.html",
                                  file_base, as.integer(topic_id), rank_tag))
  }
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(nd)) cli::cli_abort("No nodes to plot.")
  if (!nrow(edv)) cli::cli_warn("No edges after filtering; plot will show isolated nodes.")

  ok <- try(htmlwidgets::saveWidget(widget_out, file = out_html, selfcontained = TRUE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    htmlwidgets::saveWidget(
      widget_out, file = out_html, selfcontained = FALSE,
      libdir = paste0(tools::file_path_sans_ext(basename(out_html)), "_files")
    )
  }
  .llog("✓ wrote: ", normalizePath(out_html, mustWork = FALSE), verbose = verbose)
  # END EDIT

  invisible(out_html)
}

# Safe infix (duplicate here for self-contained file)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# nocov end
# END FILE: utils_grn_link_delta_loader.R
