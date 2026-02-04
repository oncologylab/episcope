#' utils_grn_filter.R ?Filter TFgene link delta tables (episcope)
#'
#' @author Yaoxiang Li
#' @family episcope-utils
#' @keywords internal
#' @noRd


# =============================
# Internal utilities (helpers)
# =============================

# Local logger (kept separate from utils_grn_diff)
.flog <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) .log_inform(c(v = paste0("[utils_grn_filter] ", paste0(..., collapse = ""))))
}

# Load a delta-links table from CSV or in-memory tibble/data.frame
.load_delta_links <- function(x, verbose = TRUE) {
  if (is.character(x) && length(x) == 1L) {
    .flog("Reading CSV: ", x, verbose = verbose)
    tbl <- readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .flog("Using in-memory tibble/data.frame (nrow=", nrow(x), ")", verbose = verbose)
    tbl <- tibble::as_tibble(x)
  } else {
    cli::cli_abort("`x` must be a CSV file path or a data.frame/tibble.")
  }
  tbl
}

# Utility: get column names that start with a given prefix
.cols_with_prefix <- function(df, prefix) {
  grep(paste0("^", prefix), names(df), value = TRUE)
}

# Utility: row-wise "any > threshold" over a set of numeric columns
.any_gt <- function(df, cols, thr) {
  if (length(cols) == 0L) {
    return(rep(TRUE, nrow(df)))  # pass-through if none
  }
  m <- as.matrix(df[, cols, drop = FALSE])
  if (nrow(m) == 0L) {
    return(logical(0))
  }
  apply(m, 1L, function(v) any(v > thr, na.rm = TRUE))
}

# Split links into up/down by gene log2FC and footprint delta, with optional TF opposition filter
.split_links_by_direction <- function(df,
                                      gene_log2_col = "log2FC_gene_expr",
                                      fp_delta_col = "delta_fp_bed_score",
                                      tf_log2_col = "log2FC_tf_expr",
                                      gene_log2_abs_min = NULL,
                                      fp_delta_min = NULL,
                                      tf_log2_abs_min = NULL,
                                      verbose = TRUE) {
  if (!nrow(df)) return(list(up = df, down = df))
  nms <- names(df)
  if (!gene_log2_col %in% nms || !fp_delta_col %in% nms) {
    .log_warn("Missing {gene_log2_col} or {fp_delta_col}; skipping up/down split.")
    return(list(up = df[0, , drop = FALSE], down = df[0, , drop = FALSE]))
  }
  if (is.null(gene_log2_abs_min) || !is.finite(gene_log2_abs_min)) {
    .log_warn("gene_log2_abs_min not set; skipping up/down split.")
    return(list(up = df[0, , drop = FALSE], down = df[0, , drop = FALSE]))
  }
  if (is.null(fp_delta_min) || !is.finite(fp_delta_min)) {
    .log_warn("fp_delta_min not set; skipping up/down split.")
    return(list(up = df[0, , drop = FALSE], down = df[0, , drop = FALSE]))
  }

  gene_l2 <- suppressWarnings(as.numeric(df[[gene_log2_col]]))
  fp_del <- suppressWarnings(as.numeric(df[[fp_delta_col]]))
  tf_l2 <- if (tf_log2_col %in% nms) suppressWarnings(as.numeric(df[[tf_log2_col]])) else NA_real_

  gene_pass <- is.finite(gene_l2) & (abs(gene_l2) >= gene_log2_abs_min)

  if (is.finite(tf_log2_abs_min) && tf_log2_col %in% nms) {
    tf_pass <- is.finite(tf_l2) & (abs(tf_l2) >= tf_log2_abs_min)
    opp_dir <- tf_pass & gene_pass & (sign(tf_l2) == -sign(gene_l2))
    if (any(opp_dir, na.rm = TRUE)) {
      df <- df[!opp_dir, , drop = FALSE]
      gene_l2 <- gene_l2[!opp_dir]
      fp_del <- fp_del[!opp_dir]
      gene_pass <- gene_pass[!opp_dir]
    }
  }

  up <- gene_pass & (gene_l2 > 0) & is.finite(fp_del) & (fp_del >= fp_delta_min)
  down <- gene_pass & (gene_l2 < 0) & is.finite(fp_del) & (fp_del <= -fp_delta_min)

  list(
    up = df[which(up), , drop = FALSE],
    down = df[which(down), , drop = FALSE]
  )
}

# Build default output path from an input CSV
.default_filtered_path <- function(in_csv) {
  base <- sub("\\.csv$", "", in_csv, ignore.case = TRUE)
  paste0(base, "_filtered.csv")
}

# Derive activator/repressor from link_sign_* columns
.link_mode_from_signs <- function(df) {
  cols <- .cols_with_prefix(df, "link_sign_")
  if (length(cols) == 0L) return(rep(NA_character_, nrow(df)))
  m  <- as.data.frame(df[, cols, drop = FALSE], stringsAsFactors = FALSE)
  nn <- rowSums(!is.na(m))
  plus_all  <- (rowSums(m == "+", na.rm = TRUE) == nn) & (nn > 0)
  minus_all <- (rowSums(m == "-", na.rm = TRUE) == nn) & (nn > 0)
  ifelse(plus_all & !minus_all, "activator",
         ifelse(minus_all & !plus_all, "repressor", NA_character_))
}

# Common dir for a set of paths
.common_dir <- function(paths) {
  u <- unique(dirname(paths))
  if (length(u) == 1L) u else getwd()
}

# Compute default filtered paths
.default_filtered_paths <- function(csvs) {
  vapply(csvs, function(p) .default_filtered_path(p), character(1))
}


# =================================
# Public API: filtering
# =================================

#' Filter delta-links table by multi-criteria (expression, FP, link, delta)
#'
#' Filter *_delta_links.csv (or in-memory tables) using:
#'   1) any gene_expr_*    > gene_expr_min
#'   2) any tf_expr_*      > tf_expr_min
#'   3) any fp_bed_score_* > fp_min
#'   4) any link_score_*   > link_min
#'   5) |delta_link_score| > abs_delta_min
#'   6) Optional absolute |log2FC| thresholds:
#'      - log2FC_gene_expr (via apply_de_gene + de_gene_log2_abs_min)
#'      - log2FC_tf_expr   (via apply_de_tf   + de_tf_log2_abs_min)
#'
#' Columns are NOT renamed; output preserves original names. If a criterions
#' columns are missing, it is skipped with a warning.
#'
#' @param x Tibble/data.frame from \code{compare_links_two_conditions()}, or path to a CSV.
#' @param gene_expr_min numeric; default 4
#' @param tf_expr_min numeric; default 4
#' @param fp_min numeric; default 4
#' @param link_min numeric; default 2
#' @param abs_delta_min numeric; default 2
#' @param save_csv logical; write filtered CSV. Default FALSE
#' @param out_file optional character; output path. If NULL and `x` is a file,
#'   writes alongside input as "*_delta_links_filtered.csv".
#' @param verbose logical; default TRUE
#' @param apply_de_gene logical; require |log2FC_gene_expr| ?de_gene_log2_abs_min
#' @param apply_de_tf   logical; require |log2FC_tf_expr|   ?de_tf_log2_abs_min
#' @param de_gene_log2_abs_min numeric absolute log2FC threshold.
#' @param de_tf_log2_abs_min   numeric absolute log2FC threshold.
#' @param de_gene_fc_min numeric absolute fold-change threshold (non-log2) for genes.
#' @param de_tf_fc_min   numeric absolute fold-change threshold (non-log2) for TFs.
#' @param fp_delta_min numeric absolute delta threshold for footprint scores.
#' @param fp_delta_min numeric absolute delta threshold for footprint scores
#'   (column `delta_fp_bed_score`).
#' @param enforce_link_expr_sign logical; if TRUE, drop rows whose delta_link_score
#'   direction contradicts expected gene direction given link_sign_* (default TRUE).
#' @param expr_dir_col character; gene-change column to use for direction; default "log2FC_gene_expr"
#'
#' @return Tibble filtered to rows passing all criteria.
#' @export
#'
#' @examples
#' \dontrun{
#' # Single file
#' out <- episcope::filter_links_deltas(
#'   "path/to/HPAFII_0_FBS_vs_HPAFII_10_FBS_delta_links.csv",
#'   gene_expr_min = 4, tf_expr_min = 4, fp_min = 4, link_min = 2, abs_delta_min = 2,
#'   apply_de_gene = TRUE, de_gene_log2_abs_min = 0.585, save_csv = TRUE
#' )
#' }
filter_links_deltas <- function(x,
                                gene_expr_min = 4,
                                tf_expr_min   = 4,
                                fp_min        = 4,
                                link_min      = 2,
                                abs_delta_min = 2,
                                apply_de_gene = FALSE,
                                apply_de_tf   = FALSE,
                                de_gene_log2_abs_min = NULL,
                                de_tf_log2_abs_min   = NULL,
                                de_gene_fc_min = NULL,
                                de_tf_fc_min   = NULL,
                                fp_delta_min = NULL,
                                enforce_link_expr_sign = TRUE,
                                expr_dir_col = "log2FC_gene_expr",
                                save_csv = FALSE,
                                out_file = NULL,
                                verbose = TRUE) {

  df <- .load_delta_links(x, verbose = verbose)
  n0 <- nrow(df)
  nms <- names(df)

  if (is.null(de_gene_log2_abs_min) && !is.null(de_gene_fc_min)) {
    de_gene_log2_abs_min <- log2(as.numeric(de_gene_fc_min))
  }
  if (is.null(de_tf_log2_abs_min) && !is.null(de_tf_fc_min)) {
    de_tf_log2_abs_min <- log2(as.numeric(de_tf_fc_min))
  }

  # Identify column groups by prefix
  cols_gene <- .cols_with_prefix(df, "gene_expr_")
  cols_tf   <- .cols_with_prefix(df, "tf_expr_")
  cols_fp   <- .cols_with_prefix(df, "fp_bed_score_")
  cols_link <- .cols_with_prefix(df, "link_score_")
  has_delta <- ("delta_link_score" %in% nms)

  .flog("Columns detected ?gene_expr_*: ", length(cols_gene),
        ", tf_expr_*: ", length(cols_tf),
        ", fp_bed_score_*: ", length(cols_fp),
        ", link_score_*: ", length(cols_link),
        ", delta_link_score: ", has_delta, ".", verbose = verbose)

  .flog("Thresholds ?gene_expr_min=", gene_expr_min,
        ", tf_expr_min=", tf_expr_min,
        ", fp_min=", fp_min,
        ", link_min=", link_min,
        ", abs_delta_min=", abs_delta_min, ".", verbose = verbose)

  # Per-criterion passes
  pass_gene <- .any_gt(df, cols_gene, gene_expr_min)
  if (length(cols_gene) == 0L) warning("No gene_expr_* columns found; skipping gene expression filter.")
  .flog("Pass gene_expr: ", sum(pass_gene), "/", n0, " rows.", verbose = verbose)

  pass_tf <- .any_gt(df, cols_tf, tf_expr_min)
  if (length(cols_tf) == 0L) warning("No tf_expr_* columns found; skipping TF expression filter.")
  .flog("Pass tf_expr: ", sum(pass_tf), "/", n0, " rows.", verbose = verbose)

  pass_fp <- .any_gt(df, cols_fp, fp_min)
  if (length(cols_fp) == 0L) warning("No fp_bed_score_* columns found; skipping FP filter.")
  .flog("Pass fp_bed_score: ", sum(pass_fp), "/", n0, " rows.", verbose = verbose)

  pass_link <- .any_gt(df, cols_link, link_min)
  if (length(cols_link) == 0L) warning("No link_score_* columns found; skipping link_score filter.")
  .flog("Pass link_score: ", sum(pass_link), "/", n0, " rows.", verbose = verbose)

  if (has_delta) {
    dv <- df[["delta_link_score"]]
    pass_delta <- (abs(dv) > abs_delta_min) & is.finite(dv)
  } else {
    warning("delta_link_score column not found; skipping delta filter.")
    pass_delta <- rep(TRUE, n0)
  }
  .flog("Pass delta_link_score: ", sum(pass_delta), "/", n0, " rows.", verbose = verbose)

  # Optional consistency between delta_link_score and gene change
  if (isTRUE(enforce_link_expr_sign)) {
    mode <- .link_mode_from_signs(df)  # "activator", "repressor", or NA
    s_link <- suppressWarnings(sign(as.numeric(df[["delta_link_score"]])))
    s_link[!is.finite(s_link)] <- 0L

    gene_col <- if (expr_dir_col %in% nms) {
      expr_dir_col
    } else if ("delta_gene_expr" %in% nms) {
      "delta_gene_expr"
    } else {
      NULL
    }

    if (is.null(gene_col)) {
      warning("enforce_link_expr_sign=TRUE but no gene direction column found; skipping consistency filter.")
      pass_sign <- rep(TRUE, n0)
    } else {
      lg <- suppressWarnings(as.numeric(df[[gene_col]]))
      s_gene <- sign(lg); s_gene[!is.finite(s_gene)] <- 0L
      prod <- s_link * s_gene
      pass_sign <- rep(TRUE, n0)               # default: pass if unknown/neutral
      idx_act <- which(mode == "activator")
      idx_rep <- which(mode == "repressor")
      if (length(idx_act)) pass_sign[idx_act] <- (prod[idx_act] ==  1) | (s_link[idx_act] == 0) | (s_gene[idx_act] == 0)
      if (length(idx_rep)) pass_sign[idx_rep] <- (prod[idx_rep] == -1) | (s_link[idx_rep] == 0) | (s_gene[idx_rep] == 0)
      .flog("Pass link-gene sign consistency: ", sum(pass_sign), "/", n0, " rows.", verbose = verbose)
    }
  } else {
    pass_sign <- rep(TRUE, n0)
  }

  # Absolute log2FC filters
  if (isTRUE(apply_de_gene)) {
    if (!("log2FC_gene_expr" %in% nms)) {
      warning("Column 'log2FC_gene_expr' not found; skipping differential GENE filter.")
      pass_de_gene <- rep(TRUE, n0)
    } else if (is.null(de_gene_log2_abs_min)) {
      warning("apply_de_gene=TRUE but de_gene_log2_abs_min is NULL; skipping differential GENE filter.")
      pass_de_gene <- rep(TRUE, n0)
    } else {
      lg <- suppressWarnings(as.numeric(df[["log2FC_gene_expr"]]))
      thr_abs <- as.numeric(de_gene_log2_abs_min)
      pass_de_gene <- is.finite(lg) & (abs(lg) >= thr_abs)
      .flog("Pass DE gene (|log2FC| ?", thr_abs, "): ",
            sum(pass_de_gene), "/", n0, " rows.", verbose = verbose)
    }
  } else {
    pass_de_gene <- rep(TRUE, n0)
  }

  if (isTRUE(apply_de_tf)) {
    if (!("log2FC_tf_expr" %in% nms)) {
      warning("Column 'log2FC_tf_expr' not found; skipping differential TF filter.")
      pass_de_tf <- rep(TRUE, n0)
    } else if (is.null(de_tf_log2_abs_min)) {
      warning("apply_de_tf=TRUE but de_tf_log2_abs_min is NULL; skipping differential TF filter.")
      pass_de_tf <- rep(TRUE, n0)
    } else {
      lt <- suppressWarnings(as.numeric(df[["log2FC_tf_expr"]]))
      thr_abs <- as.numeric(de_tf_log2_abs_min)
      pass_de_tf <- is.finite(lt) & (abs(lt) >= thr_abs)
      .flog("Pass DE TF (|log2FC| ?", thr_abs, "): ",
            sum(pass_de_tf), "/", n0, " rows.", verbose = verbose)
    }
  } else {
    pass_de_tf <- rep(TRUE, n0)
  }

  # Combine (logical AND)
  keep <- pass_gene & pass_tf & pass_fp & pass_link & pass_delta & pass_de_gene & pass_de_tf & pass_sign

  out <- df[keep, , drop = FALSE]
  .flog("Kept ", nrow(out), " / ", n0, " rows after all filters.", verbose = verbose)

  # Write CSV if requested
  if (isTRUE(save_csv)) {
    if (is.null(out_file)) {
      if (is.character(x) && length(x) == 1L) {
        out_file <- .default_filtered_path(x)
      } else {
        cli::cli_abort("save_csv=TRUE but `out_file` is NULL and `x` is not a file path.")
      }
    }
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out, out_file)
    .flog("?wrote: ", out_file, verbose = verbose)
  }

  out
}


#' Bulk/parallel filter over many delta-link CSVs
#'
#' Applies \code{filter_links_deltas()} across many \code{*_delta_links.csv} files,
#' optionally in parallel, and writes a manifest mapping source -> filtered outputs.
#'
#' @param delta_csvs Character vector of paths to *_delta_links.csv files
#'                   (e.g., \code{specs$out_file} from \code{build_cellwise_contrasts_from_index()}).
#' @inheritParams filter_links_deltas
#' @param parallel Logical; run in parallel via \pkg{future.apply}. Default TRUE.
#' @param plan Future plan; string (e.g. "multisession") or function. Default "multisession".
#' @param workers Integer; default: availableCores()-1.
#' @param save_csv Logical; write *_filtered.csv next to each input (default TRUE).
#' @param write_manifest Logical; write a manifest CSV mapping source->filtered (default TRUE).
#' @param manifest_file Optional path for manifest; if NULL, written to the common dir as
#'        "delta_links_filtered_manifest.csv".
#' @param split_direction Logical; write up/down filtered links per delta file.
#' @param write_combined Logical; write combined filtered CSV (default TRUE).
#' @return A list with fields: results (list of tibbles), filtered_paths, filtered_dirs, manifest (tibble), manifest_path.
#' @export
#'
#' @examples
#' \dontrun{
#' delta_csvs <- list.files("inst/extdata/lighting", "_delta_links.csv", full.names = TRUE)
#' de_gene_log2_abs_min <- 0.585
#' bulk <- episcope::filter_links_deltas_bulk(
#'   delta_csvs,
#'   gene_expr_min = 4, tf_expr_min = 4, fp_min = 4, link_min = 2, abs_delta_min = 2,
#'   apply_de_gene = TRUE,
#'   de_gene_log2_abs_min = de_gene_log2_abs_min,
#'   enforce_link_expr_sign = TRUE,
#'   expr_dir_col = "log2FC_gene_expr",
#'   workers = 20
#' )
#' bulk$filtered_paths
#' bulk$manifest_path
#' }
filter_links_deltas_bulk <- function(delta_csvs,
                                     gene_expr_min = 4,
                                     tf_expr_min   = 4,
                                     fp_min        = 4,
                                     link_min      = 2,
                                     abs_delta_min = 2,
                                     apply_de_gene = FALSE,
                                     apply_de_tf   = FALSE,
                                     de_gene_log2_abs_min = NULL,
                                     de_tf_log2_abs_min   = NULL,
                                     de_gene_fc_min = NULL,
                                     de_tf_fc_min   = NULL,
                                     fp_delta_min = NULL,
                                     split_direction = FALSE,
                                     write_combined = TRUE,
                                     enforce_link_expr_sign = FALSE,
                                     expr_dir_col = "log2FC_gene_expr",
                                     parallel = TRUE,
                                     plan = "multisession",
                                     workers = NULL,
                                     save_csv = TRUE,
                                     write_manifest = TRUE,
                                     manifest_file = NULL,
                                     filtered_dir = NULL,
                                     verbose = TRUE) {
  # sanity
  if (!is.character(delta_csvs) || length(delta_csvs) == 0L) {
    cli::cli_abort("`delta_csvs` must be a non-empty character vector of CSV paths.")
  }
  exist_mask <- file.exists(delta_csvs)
  if (!all(exist_mask)) {
    missing <- delta_csvs[!exist_mask]
    warning(
      "These input CSVs do not exist and will be skipped:\n  - ",
      paste(missing, collapse = "\n  - ")
    )
    delta_csvs <- delta_csvs[exist_mask]
  }
  if (length(delta_csvs) == 0L) cli::cli_abort("No existing delta CSVs to process.")

  # expected outputs (for manifest)
  if (!is.null(filtered_dir) && nzchar(filtered_dir)) {
    filtered_paths <- file.path(filtered_dir, basename(.default_filtered_paths(delta_csvs)))
  } else {
    filtered_paths <- .default_filtered_paths(delta_csvs)
  }
  filtered_dirs  <- unique(dirname(filtered_paths))

  if (is.null(de_gene_log2_abs_min) && !is.null(de_gene_fc_min)) {
    de_gene_log2_abs_min <- log2(as.numeric(de_gene_fc_min))
  }
  if (is.null(de_tf_log2_abs_min) && !is.null(de_tf_fc_min)) {
    de_tf_log2_abs_min <- log2(as.numeric(de_tf_fc_min))
  }

  n_files <- length(delta_csvs)
  if (isTRUE(verbose)) {
    .log_inform("Filtering {n_files} delta-link file{?s}...", n_files = n_files)
    if (isTRUE(split_direction)) {
      .log_inform("Writing up/down filtered links per delta file.")
    }
  }

  # parallel plan
  if (isTRUE(parallel)) {
    strategy <- (function(x) {
      if (is.character(x)) {
        ns <- asNamespace("future")
        if (exists(x, envir = ns, mode = "function")) get(x, envir = ns) else future::multisession
      } else if (is.function(x)) x else future::multisession
    })(plan)
    if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    .flog("Launching future plan=", if (is.character(plan)) plan else "custom",
          ", workers=", workers, verbose = verbose)
    oplan <- future::plan()
    on.exit(future::plan(oplan), add = TRUE)
    future::plan(strategy, workers = workers)
  }

  worker_verbose <- if (isTRUE(parallel)) FALSE else verbose
  run_one <- function(f, out_path) {
    out_dir_use <- dirname(out_path)
    if (!is.null(filtered_dir) && nzchar(filtered_dir)) {
      out_dir_use <- filtered_dir
    }
    df <- filter_links_deltas(
      f,
      gene_expr_min = gene_expr_min,
      tf_expr_min   = tf_expr_min,
      fp_min        = fp_min,
      link_min      = link_min,
      abs_delta_min = abs_delta_min,
      apply_de_gene = apply_de_gene,
      apply_de_tf   = apply_de_tf,
      de_gene_log2_abs_min = de_gene_log2_abs_min,
      de_tf_log2_abs_min   = de_tf_log2_abs_min,
      de_gene_fc_min = de_gene_fc_min,
      de_tf_fc_min   = de_tf_fc_min,
      fp_delta_min = fp_delta_min,
      enforce_link_expr_sign = enforce_link_expr_sign,
      expr_dir_col = expr_dir_col,
      save_csv = isTRUE(save_csv) && isTRUE(write_combined),
      out_file = out_path,
      verbose = worker_verbose
    )
    if (isTRUE(split_direction)) {
      base <- sub("\\.csv$", "", basename(f), ignore.case = TRUE)
      up_path <- file.path(out_dir_use, paste0(base, "_filtered_up.csv"))
      down_path <- file.path(out_dir_use, paste0(base, "_filtered_down.csv"))
      split <- .split_links_by_direction(
        df,
        gene_log2_abs_min = de_gene_log2_abs_min,
        fp_delta_min = fp_delta_min,
        tf_log2_abs_min = de_tf_log2_abs_min,
        verbose = worker_verbose
      )
      readr::write_csv(split$up, up_path)
      readr::write_csv(split$down, down_path)
      attr(df, "filtered_up_path") <- up_path
      attr(df, "filtered_down_path") <- down_path
    }
    df
  }

  results <- if (isTRUE(parallel)) {
    if (isTRUE(verbose) && requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(along = seq_len(n_files))
        do_one_wrap <- function(i) {
          on.exit(p(), add = TRUE)
          run_one(delta_csvs[[i]], filtered_paths[[i]])
        }
        future.apply::future_lapply(seq_along(delta_csvs), do_one_wrap)
      })
    } else {
      future.apply::future_lapply(seq_along(delta_csvs), function(i) run_one(delta_csvs[[i]], filtered_paths[[i]]))
    }
  } else {
    lapply(seq_along(delta_csvs), function(i) run_one(delta_csvs[[i]], filtered_paths[[i]]))
  }
  names(results) <- basename(delta_csvs)

  # manifest
  manifest <- tibble::tibble(
    source_delta = delta_csvs,
    filtered_csv = filtered_paths,
    exists       = file.exists(filtered_paths)
  )
  if (isTRUE(split_direction)) {
    up_paths <- vapply(results, function(x) attr(x, "filtered_up_path"), character(1))
    down_paths <- vapply(results, function(x) attr(x, "filtered_down_path"), character(1))
    manifest$filtered_up <- up_paths
    manifest$filtered_down <- down_paths
    manifest$up_exists <- file.exists(up_paths)
    manifest$down_exists <- file.exists(down_paths)
  }

  manifest_path <- NULL
  if (isTRUE(write_manifest)) {
    if (is.null(manifest_file)) {
      if (!is.null(filtered_dir) && nzchar(filtered_dir)) {
        manifest_path <- file.path(filtered_dir, "delta_links_filtered_manifest.csv")
      } else {
        common_dir <- .common_dir(delta_csvs)
        manifest_path <- file.path(common_dir, "delta_links_filtered_manifest.csv")
      }
    } else {
      manifest_path <- manifest_file
    }
    dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(manifest, manifest_path)
    .flog("?wrote manifest: ", manifest_path, verbose = verbose)
  }

  list(
    results        = results,
    filtered_paths = filtered_paths,
    filtered_dirs  = filtered_dirs,
    manifest       = manifest,
    manifest_path  = manifest_path
  )
}
