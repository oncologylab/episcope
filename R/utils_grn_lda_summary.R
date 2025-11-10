#' utils_grn_lda_summary.R — Summarize LDA-annotated delta-link tables (episcope)
#'
#' @author Yaoxiang Li
#' @family episcope-utils
#' @keywords internal
#'
#' @description
#' Summarize *_filtered_lda_K{K}.csv files that include an \code{LDA_K{K}_topic}
#' column (multi-topic like "6;11" allowed). Produces per-topic metrics and
#' per-(TF, topic) metrics, plus a topic rank.
#'
#' @section Topic-level metrics:
#' \itemize{
#'   \item \code{topic_mean_abs_delta}
#'   \item \code{topic_n_TFs}
#'   \item \code{topic_n_links}  (renamed from topic_n_edges; same calculation)
#'   \item \code{topic_rank}     (descending by \code{topic_mean_abs_delta})
#' }
#'
#' @section TF × topic metrics:
#' \itemize{
#'   \item \code{tf_delta_sum} (signed, using strict regulatory sign)
#'   \item \code{tf_delta_sum_activate}, \code{tf_delta_sum_repress}
#'   \item \code{tf_delta_sum_abs}, \code{tf_delta_sum_abs_activate}, \code{tf_delta_sum_abs_repress}
#'   \item \code{tf_n_links}, \code{tf_n_links_activate}, \code{tf_n_links_repress}
#'   \item \code{tf_expr_max} (max across any \code{tf_expr_*} columns, if present)
#'   \item \code{tf_log2_fc}  (from \code{log2FC_tf_expr}, if present)
#'   \item \code{tf_hits_hub} (per-topic HITS hub score scaled to [0,1])
#' }
#'
#' @details
#' \strong{STRICT sign policy.} Exactly two columns named \code{link_sign_*}
#' must be present (e.g., \code{link_sign_cond1}, \code{link_sign_cond2}).
#' After filtering by |delta|, every kept row must have valid signs and the two
#' signs must agree; otherwise an error is thrown.
#'
#' Signed aggregation uses
#' \deqn{\mathrm{delta\_oriented} = \mathrm{delta\_link\_score} \times \mathrm{regulatory\_sign},}
#' where \eqn{\mathrm{regulatory\_sign} \in \{+1,-1\}} derived from \code{link_sign_*}
#' (“+”→+1, “−”→−1). Topic-level metrics are based on |delta| and are unaffected by sign.
#'
#' @section Imports:
#' Uses \pkg{cli}, \pkg{readr}, \pkg{tibble}, \pkg{dplyr}, \pkg{tidyr},
#' \pkg{future}, \pkg{future.apply}, \pkg{parallel}, \pkg{igraph}, \pkg{stats}.
#'
#' @noRd


# =============================
# Internal utilities (helpers)
# =============================

.slog <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message("[utils_grn_lda_summary] ", paste0(..., collapse = ""))
}

.load_annotated <- function(x, verbose = TRUE) {
  if (is.character(x) && length(x) == 1L) {
    .slog("Reading CSV: ", x, verbose = verbose)
    tbl <- readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .slog("Using in-memory tibble/data.frame (nrow=", nrow(x), ")", verbose = verbose)
    tbl <- tibble::as_tibble(x)
  } else {
    cli::cli_abort("`x` must be a CSV path or a data.frame/tibble.")
  }
  tbl
}

.pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    cli::cli_abort("Missing required column for {label}: one of {paste(candidates, collapse = ', ')}.")
  }
  hit[[1]]
}

.detect_topic_col <- function(df) {
  ks <- grep("^LDA_K\\d+_topic$", names(df), value = TRUE)
  if (length(ks) == 0L) cli::cli_abort("No LDA topic column found (expected 'LDA_K{K}_topic').")
  if (length(ks) == 1L) return(ks[[1]])
  Knum <- as.integer(gsub("^LDA_K(\\d+)_topic$", "\\1", ks))
  ks[order(Knum, decreasing = TRUE)][[1]]
}

.default_summary_path <- function(in_csv) {
  paste0(sub("\\.csv$", "", in_csv, ignore.case = TRUE), "_summary.csv")
}

# Map "+" / "−" → +1 / −1; others → NA_real_
.str_sign_to_num <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  out[!is.na(x) & x == "+"] <-  1
  out[!is.na(x) & x == "-"] <- -1
  out
}

# STRICT: require EXACTLY two link_sign_* columns and perfect agreement per kept row.
# Returns numeric vector (+1/-1) of the agreed sign for each row of `df`.
.compute_reg_sign_strict <- function(df, tf_col, gene_col) {
  sign_cols <- grep("^link_sign_", names(df), value = TRUE)
  if (length(sign_cols) != 2L) {
    cli::cli_abort("STRICT sign check requires exactly two link_sign_* columns; found {length(sign_cols)}: {paste(sign_cols, collapse = ', ')}.")
  }
  sign_cols <- sort(sign_cols) # deterministic
  s1_chr <- df[[sign_cols[1]]]
  s2_chr <- df[[sign_cols[2]]]
  s1 <- .str_sign_to_num(s1_chr)
  s2 <- .str_sign_to_num(s2_chr)

  invalid_mask  <- is.na(s1) | is.na(s2)
  disagree_mask <- (!invalid_mask) & (s1 != s2)

  bad_n <- sum(invalid_mask | disagree_mask, na.rm = TRUE)
  if (bad_n > 0) {
    ex_idx <- which(invalid_mask | disagree_mask)[seq_len(min(10L, bad_n))]
    tf_v   <- if (!missing(tf_col)   && tf_col   %in% names(df)) df[[tf_col]][ex_idx]   else NA_character_
    gene_v <- if (!missing(gene_col) && gene_col %in% names(df)) df[[gene_col]][ex_idx] else NA_character_
    msg <- paste0(
      "Found ", bad_n, " row(s) with missing/invalid or disagreeing link_sign_*.\n",
      "First examples (up to 10):\n",
      paste0("  TF=", tf_v, " | gene_key=", gene_v,
             " | ", sign_cols[1], "=", s1_chr[ex_idx],
             " | ", sign_cols[2], "=", s2_chr[ex_idx], collapse = "\n")
    )
    stop(msg)
  }
  s1
}

# HITS helper (per-topic hub score per TF; scaled to [0,1])
.compute_tf_hits_by_topic <- function(expanded_df, tf_col, gene_col) {
  if (!nrow(expanded_df)) {
    out <- tibble::tibble(topic = integer(0))
    out[[tf_col]] <- character(0)
    out[["tf_hits_hub"]] <- numeric(0)
    return(out)
  }

  topics <- sort(unique(expanded_df$topic))
  out_list <- vector("list", length(topics))

  for (i in seq_along(topics)) {
    tt <- topics[i]
    sub <- expanded_df[expanded_df$topic == tt, , drop = FALSE]

    el <- unique(sub[, c(tf_col, gene_col), drop = FALSE])
    el <- el[stats::complete.cases(el), , drop = FALSE]
    if (!nrow(el)) {
      out <- tibble::tibble(topic = integer(0))
      out[[tf_col]] <- character(0)
      out[["tf_hits_hub"]] <- numeric(0)
      out_list[[i]] <- out
      next
    }

    g <- igraph::graph_from_data_frame(
      d = stats::setNames(el, c("from","to")),
      directed = TRUE
    )

    # Preferred: HITS with no weights
    hs <- tryCatch(igraph::hits(g, weights = NULL)$hub,
                   error = function(e) rep(0, igraph::gorder(g)))
    vnames <- igraph::V(g)$name
    names(hs) <- vnames

    tf_names <- unique(el[[tf_col]])
    hs_tf <- hs[match(tf_names, names(hs))]
    hs_tf[!is.finite(hs_tf)] <- 0

    # Fallback: normalized out-degree if all zeros
    if (all(hs_tf == 0, na.rm = TRUE)) {
      deg_out <- igraph::degree(g, mode = "out")
      names(deg_out) <- vnames
      hs_tf <- as.numeric(deg_out[match(tf_names, names(deg_out))])
      hs_tf[!is.finite(hs_tf)] <- 0
    }

    # Scale to [0,1] within topic
    maxv <- max(hs_tf, na.rm = TRUE)
    if (!is.finite(maxv) || maxv <= 0) hs_tf[] <- 0 else hs_tf <- hs_tf / maxv

    n <- length(tf_names)
    out <- tibble::tibble(
      topic = rep.int(tt, n),
      !!rlang::sym(tf_col) := tf_names,
      tf_hits_hub = as.numeric(hs_tf)
    )
    out_list[[i]] <- out
  }

  dplyr::bind_rows(out_list)
}

# =============================
# Public API: single file
# =============================

#' Summarize one LDA-annotated delta-links table (TF × topic)
#'
#' @param x Annotated tibble or CSV path (must include an \code{LDA_K{K}_topic} column).
#' @param topic_col Optional topic column (auto-detects the latest \code{^LDA_K\\d+_topic$}).
#' @param multi_sep Separator for multi-topic strings (default \code{";"}).
#' @param edge_filter_min Keep links with \code{|delta_link_score| >=} this (default \code{0.25}).
#' @param save_csv Write \code{"*_summary.csv"} if `x` is a path (default \code{TRUE}).
#' @param out_file Optional override for output CSV path.
#' @param verbose Verbose logging (default \code{TRUE}).
#' @param use_tf_r_weight Logical; when \code{TRUE}, multiply \code{delta_link_score}
#'   by a TF-level \code{r} column (after |delta| filtering, before sign orientation).
#' @param tf_r_col Optional explicit TF r column (e.g., \code{"r_tf_ctrl"}). If \code{NULL},
#'   auto-detects.
#'
#' @return Tibble with per-(TF, topic) rows and topic-level ranks.
#' @export
#' @importFrom readr read_csv write_csv
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr mutate filter group_by summarise n_distinct arrange desc row_number left_join
#' @importFrom dplyr select distinct bind_rows ungroup across rename
#' @importFrom tidyr separate_rows
#' @importFrom rlang sym
summarize_lda_annotations <- function(x,
                                      topic_col = NULL,
                                      multi_sep = ";",
                                      edge_filter_min = 0.25,
                                      save_csv = TRUE,
                                      out_file = NULL,
                                      verbose = TRUE,
                                      use_tf_r_weight = FALSE,
                                      tf_r_col = NULL) {
  df <- .load_annotated(x, verbose = verbose)

  tf_col    <- .pick_col(df, c("tf","TF"), "TF")
  gene_col  <- .pick_col(df, c("gene_key","gene","target"), "gene_key")
  delta_col <- .pick_col(df, c("delta_link_score","d_edge_score","delta"), "delta_link_score")
  if (is.null(topic_col)) topic_col <- .detect_topic_col(df)

  # Filter by |delta|
  delta_vec <- suppressWarnings(as.numeric(df[[delta_col]]))
  edges_f <- tibble::as_tibble(df)
  edges_f$delta <- delta_vec
  edges_f <- edges_f[is.finite(edges_f$delta) & abs(edges_f$delta) >= edge_filter_min, , drop = FALSE]

  if (!nrow(edges_f)) {
    warning("No links pass |delta| >= ", edge_filter_min, " — summary will be empty.")
    empty <- tibble::tibble(
      TF = character(0), topic = integer(0),
      topic_mean_abs_delta = numeric(0), topic_n_TFs = integer(0),
      topic_n_links = integer(0), topic_rank = integer(0),
      tf_delta_sum = numeric(0),
      tf_delta_sum_activate = numeric(0), tf_delta_sum_repress = numeric(0),
      tf_delta_sum_abs = numeric(0),
      tf_delta_sum_abs_activate = numeric(0), tf_delta_sum_abs_repress = numeric(0),
      tf_n_links = integer(0),
      tf_n_links_activate = integer(0), tf_n_links_repress = integer(0),
      tf_expr_max = numeric(0), tf_log2_fc = numeric(0),
      tf_hits_hub = numeric(0)
    )
    if (isTRUE(save_csv) && is.character(x) && length(x) == 1L) {
      if (is.null(out_file)) out_file <- .default_summary_path(x)
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(empty, out_file)
      .slog("✓ wrote empty summary: ", out_file, verbose = verbose)
    }
    return(empty)
  }

  # Optional: weight delta by a TF r column (after filtering, before orientation)
  if (isTRUE(use_tf_r_weight)) {
    picked_tf_r_col <- tf_r_col
    if (is.null(picked_tf_r_col)) {
      hits <- grep("^r_tf(_|$)", names(edges_f), value = TRUE)
      if (!length(hits)) hits <- grep("(^|_)tf_r($|_)", names(edges_f), value = TRUE)
      if (length(hits)) picked_tf_r_col <- hits[[1]]
    }
    if (is.null(picked_tf_r_col) || !(picked_tf_r_col %in% names(edges_f))) {
      cli::cli_abort("use_tf_r_weight=TRUE but no suitable 'r_tf_*' column found. Provide tf_r_col= or include one in the input.")
    }
    r_w <- suppressWarnings(as.numeric(edges_f[[picked_tf_r_col]]))
    r_w[!is.finite(r_w)] <- 0
    edges_f$delta <- edges_f$delta * r_w
    .slog("Applied TF r weighting using column: ", picked_tf_r_col, verbose = verbose)
  }

  # STRICT regulatory sign (+1 / −1)
  reg_sign <- .compute_reg_sign_strict(edges_f, tf_col = tf_col, gene_col = gene_col)

  # Oriented/abs deltas
  edges_f$reg_sign       <- reg_sign
  edges_f$delta_oriented <- edges_f$delta * edges_f$reg_sign
  edges_f$delta_abs      <- abs(edges_f$delta)

  # Expand multi-topic membership
  topic_chr <- trimws(as.character(edges_f[[topic_col]]))
  keep_topic <- !is.na(topic_chr) & nzchar(topic_chr)
  edges_keep <- edges_f[keep_topic, , drop = FALSE]
  if (!nrow(edges_keep)) {
    warning("No rows with valid topics after expanding '", topic_col, "'.")
    if (isTRUE(save_csv) && is.character(x) && length(x) == 1L) {
      if (is.null(out_file)) out_file <- .default_summary_path(x)
      readr::write_csv(tibble::tibble(), out_file)
      .slog("✓ wrote empty summary: ", out_file, verbose = verbose)
    }
    return(tibble::tibble())
  }

  expanded <- edges_keep |>
    dplyr::mutate(topic_str = topic_chr[keep_topic]) |>
    tidyr::separate_rows(topic_str, sep = multi_sep, convert = TRUE) |>
    dplyr::mutate(topic = as.integer(topic_str)) |>
    dplyr::filter(is.finite(.data$topic))

  if (!nrow(expanded)) {
    warning("No valid integer topics after split.")
    if (isTRUE(save_csv) && is.character(x) && length(x) == 1L) {
      if (is.null(out_file)) out_file <- .default_summary_path(x)
      readr::write_csv(tibble::tibble(), out_file)
      .slog("✓ wrote empty summary: ", out_file, verbose = verbose)
    }
    return(tibble::tibble())
  }

  # Topic-level summary (|delta|)
  topic_summary <- expanded |>
    dplyr::group_by(.data$topic) |>
    dplyr::summarise(
      topic_mean_abs_delta = mean(delta_abs, na.rm = TRUE),
      topic_n_TFs          = dplyr::n_distinct(.data[[tf_col]]),
      topic_n_links        = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$topic_mean_abs_delta)) |>
    dplyr::mutate(topic_rank = dplyr::row_number())

  # TF × topic summaries
  tf_topic <- expanded |>
    dplyr::group_by(.data$topic, .data[[tf_col]]) |>
    dplyr::summarise(
      tf_delta_sum                  = sum(delta_oriented, na.rm = TRUE),
      tf_delta_sum_activate         = sum(ifelse(reg_sign > 0,  delta_oriented, 0), na.rm = TRUE),
      tf_delta_sum_repress          = sum(ifelse(reg_sign < 0,  delta_oriented, 0), na.rm = TRUE),
      tf_delta_sum_abs              = sum(delta_abs, na.rm = TRUE),
      tf_delta_sum_abs_activate     = sum(ifelse(reg_sign > 0,  delta_abs, 0), na.rm = TRUE),
      tf_delta_sum_abs_repress      = sum(ifelse(reg_sign < 0,  delta_abs, 0), na.rm = TRUE),
      tf_n_links                    = dplyr::n(),
      tf_n_links_activate           = sum(reg_sign > 0, na.rm = TRUE),
      tf_n_links_repress            = sum(reg_sign < 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(TF = !!rlang::sym(tf_col))

  # TF hub (HITS) per topic
  hits_tbl <- .compute_tf_hits_by_topic(
    expanded_df = expanded[, c(tf_col, gene_col, "topic")],
    tf_col = tf_col,
    gene_col = gene_col
  ) |>
    dplyr::rename(TF = dplyr::all_of(tf_col))

  tf_topic <- dplyr::left_join(tf_topic, hits_tbl, by = c("topic", "TF"))
  if (!("tf_hits_hub" %in% names(tf_topic))) tf_topic$tf_hits_hub <- NA_real_
  tf_topic$tf_hits_hub[!is.finite(tf_topic$tf_hits_hub)] <- 0

  # Per-TF expression max across tf_expr_* (optional)
  tf_expr_cols <- grep("^tf_expr_", names(df), value = TRUE)
  tf_expr_tbl <- if (length(tf_expr_cols)) {
    # Take first non-NA numeric per TF and then rowwise max
    df_tmp <- tibble::as_tibble(df)
    df_tmp$TF <- df_tmp[[tf_col]]
    agg <- df_tmp |>
      dplyr::group_by(.data$TF) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(tf_expr_cols),
          ~ suppressWarnings(dplyr::first(na.omit(as.numeric(.x)))),
          .names = "{.col}"
        ),
        .groups = "drop"
      )
    if (!nrow(agg)) {
      tibble::tibble(TF = unique(df[[tf_col]]), tf_expr_max = NA_real_)
    } else {
      agg |>
        dplyr::rowwise() |>
        dplyr::mutate(tf_expr_max = max(c(dplyr::c_across(dplyr::all_of(tf_expr_cols))), na.rm = TRUE)) |>
        dplyr::ungroup() |>
        dplyr::select(.data$TF, .data$tf_expr_max)
    }
  } else {
    tibble::tibble(TF = unique(df[[tf_col]]), tf_expr_max = NA_real_)
  }

  # Per-TF log2FC of TF expression (optional)
  tf_fc_tbl <- if ("log2FC_tf_expr" %in% names(df)) {
    df_tmp <- tibble::as_tibble(df)
    df_tmp$TF <- df_tmp[[tf_col]]
    df_tmp |>
      dplyr::group_by(.data$TF) |>
      dplyr::summarise(
        tf_log2_fc = suppressWarnings(dplyr::first(na.omit(as.numeric(.data$log2FC_tf_expr)))),
        .groups = "drop"
      )
  } else {
    tibble::tibble(TF = unique(df[[tf_col]]), tf_log2_fc = NA_real_)
  }

  tf_topic_summary <- tf_topic |>
    dplyr::left_join(topic_summary, by = "topic") |>
    dplyr::left_join(tf_expr_tbl,   by = "TF") |>
    dplyr::left_join(tf_fc_tbl,     by = "TF") |>
    dplyr::arrange(.data$topic_rank, dplyr::desc(abs(.data$tf_delta_sum))) |>
    dplyr::select(
      .data$TF, .data$topic,
      .data$topic_mean_abs_delta, .data$topic_n_TFs, .data$topic_n_links, .data$topic_rank,
      .data$tf_delta_sum,
      .data$tf_delta_sum_activate, .data$tf_delta_sum_repress,
      .data$tf_delta_sum_abs,
      .data$tf_delta_sum_abs_activate, .data$tf_delta_sum_abs_repress,
      .data$tf_n_links, .data$tf_n_links_activate, .data$tf_n_links_repress,
      .data$tf_expr_max, .data$tf_log2_fc,
      .data$tf_hits_hub
    )

  if (isTRUE(save_csv) && is.character(x) && length(x) == 1L) {
    if (is.null(out_file)) out_file <- .default_summary_path(x)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(tf_topic_summary, out_file)
    .slog("✓ wrote: ", out_file, verbose = verbose)
  }

  tf_topic_summary
}

#' Bulk/parallel summarization over many annotated CSVs
#' (writes "*_summary.csv" next to each input)
#'
#' @param annotated_csvs Character vector of annotated CSV paths.
#' @param topic_col Optional topic column name (auto-detect by default).
#' @param multi_sep Topic separator (default ";").
#' @param edge_filter_min |delta| threshold (default 0.25).
#' @param parallel Use future.apply (default TRUE).
#' @param plan Future plan (string or function; default "multisession").
#' @param workers Integer; default: parallel::detectCores()-1.
#' @param write_manifest Write mapping CSV (default TRUE).
#' @param manifest_file Optional manifest path (default next to first input).
#' @param verbose Verbose logging (default TRUE).
#' @param use_tf_r_weight Pass-through to `summarize_lda_annotations()`. Default FALSE.
#' @param tf_r_col Pass-through to `summarize_lda_annotations()`. Default NULL.
#'
#' @return List: results, summary_paths, summary_dirs, manifest, manifest_path
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom parallel detectCores
#' @importFrom tibble tibble
#' @importFrom readr write_csv
summarize_lda_annotations_bulk <- function(annotated_csvs,
                                           topic_col = NULL,
                                           multi_sep = ";",
                                           edge_filter_min = 0.25,
                                           parallel = TRUE,
                                           plan = "multisession",
                                           workers = NULL,
                                           write_manifest = TRUE,
                                           manifest_file = NULL,
                                           verbose = TRUE,
                                           use_tf_r_weight = FALSE,
                                           tf_r_col = NULL) {
  if (!is.character(annotated_csvs) || length(annotated_csvs) == 0L) {
    cli::cli_abort("`annotated_csvs` must be a non-empty character vector of CSV paths.")
  }
  exist_mask <- file.exists(annotated_csvs)
  if (!all(exist_mask)) {
    missing <- annotated_csvs[!exist_mask]
    warning("These annotated CSVs do not exist and will be skipped:\n  - ",
            paste(missing, collapse = "\n  - "))
    annotated_csvs <- annotated_csvs[exist_mask]
  }
  if (length(annotated_csvs) == 0L) cli::cli_abort("No existing annotated CSVs to process.")

  summary_paths <- vapply(annotated_csvs, .default_summary_path, character(1))
  summary_dirs  <- unique(dirname(summary_paths))

  if (isTRUE(parallel)) {
    strategy <- (function(x) {
      if (is.character(x)) {
        ns <- asNamespace("future")
        if (exists(x, envir = ns, mode = "function")) get(x, envir = ns) else future::multisession
      } else if (is.function(x)) x else future::multisession
    })(plan)
    if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    .slog("Launching future plan=", if (is.character(plan)) plan else "custom",
          ", workers=", workers, verbose = verbose)
    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    future::plan(strategy, workers = workers)
  }

  run_one <- function(f) {
    summarize_lda_annotations(
      f,
      topic_col = topic_col,
      multi_sep = multi_sep,
      edge_filter_min = edge_filter_min,
      save_csv = TRUE,
      verbose = verbose,
      use_tf_r_weight = use_tf_r_weight,
      tf_r_col = tf_r_col
    )
  }

  results <- if (isTRUE(parallel)) {
    future.apply::future_lapply(annotated_csvs, run_one, future.seed = TRUE)
  } else {
    lapply(annotated_csvs, run_one)
  }
  names(results) <- basename(annotated_csvs)

  manifest <- tibble::tibble(
    annotated_csv = annotated_csvs,
    summary_csv   = summary_paths,
    exists        = file.exists(summary_paths)
  )

  manifest_path <- NULL
  if (isTRUE(write_manifest)) {
    manifest_path <- if (is.null(manifest_file)) {
      file.path(dirname(annotated_csvs[[1]]), "delta_links_lda_summary_manifest.csv")
    } else {
      manifest_file
    }
    dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(manifest, manifest_path)
    .slog("✓ wrote manifest: ", manifest_path, verbose = verbose)
  }

  list(
    results       = results,
    summary_paths = summary_paths,
    summary_dirs  = summary_dirs,
    manifest      = manifest,
    manifest_path = manifest_path
  )
}
