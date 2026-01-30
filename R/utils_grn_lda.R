#' utils_grn_lda.R ?Annotate delta-link tables with LDA topic (episcope)
#'
#' @author Yaoxiang Li
#' @family episcope-utils
#' @keywords internal
#'
#' @details
#' Documents = TFs; terms = gene_key; weights = f(delta_link_score) with one of:
#' \itemize{
#'   \item \code{which="abs"}: \eqn{|delta_link_score|}
#'   \item \code{which="gain"}: \eqn{max(delta_link_score, 0)}
#'   \item \code{which="loss"}: \eqn{max(-delta_link_score, 0)}
#' }
#' Multi-topic assignment per TF keeps \emph{all} topics with posterior gamma
#' \eqn{\ge} \code{gamma_cutoff} (default 0.2), joined back as a
#' semicolon-separated string in a new column \code{LDA_K{K}_topic}.
#'
#' @section Imports:
#' Uses \pkg{readr}, \pkg{tibble}, \pkg{dplyr}, \pkg{tidyr},
#' \pkg{topicmodels}, \pkg{slam}, \pkg{future}, \pkg{future.apply},
#' \pkg{parallel}, \pkg{stats}, \pkg{cli}.
#'
#' @noRd


# =============================
# Internal utilities (helpers)
# =============================

# Local logger (kept separate from other utils)
.ldalog <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) message("[utils_grn_lda] ", paste0(..., collapse = ""))
}

# Load a filtered delta-links table from CSV or in-memory tibble/data.frame
.load_table <- function(x, verbose = TRUE) {
  if (is.character(x) && length(x) == 1L) {
    .ldalog("Reading CSV: ", x, verbose = verbose)
    tbl <- readr::read_csv(x, show_col_types = FALSE)
  } else if (inherits(x, "data.frame")) {
    .ldalog("Using in-memory tibble/data.frame (nrow=", nrow(x), ")", verbose = verbose)
    tbl <- tibble::as_tibble(x)
  } else {
    cli::cli_abort("`x` must be a CSV file path or a data.frame/tibble.")
  }
  tbl
}

# Pick the first existing column from a set of synonyms
.pick_col <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    cli::cli_abort("Required column missing for {label}: one of {paste(candidates, collapse=', ')}.")
  }
  hit[[1]]
}

# Default output path: append _lda_K{K}.csv
.default_lda_out_path <- function(in_csv, K) {
  base <- sub("\\.csv$", "", in_csv, ignore.case = TRUE)
  paste0(base, "_lda_K", K, ".csv")
}

# Topicmodels control (simple, seed only for VEM; basic defaults for Gibbs)
.tm_control <- function(method = c("VEM","Gibbs"), seed = 1L) {
  method <- match.arg(method)
  if (method == "VEM") {
    list(seed = seed)
  } else {
    list(seed = seed, burnin = 2000, iter = 2000, thin = 100)
  }
}

# Convert topic label to plain integer (works for "1" or "Topic1")
.safe_topic_int <- function(x) {
  as.integer(gsub("[^0-9]", "", as.character(x)))
}

# Core LDA fit that returns per-TF multi-topic assignments (gamma >= cutoff)
# Avoid .data pronoun by renaming to fixed column names up-front
.fit_lda_assign_topics <- function(df,
                                   which = c("abs","gain","loss"),
                                   min_df = 2L,
                                   K = 20,
                                   seed = 1L,
                                   method = c("VEM","Gibbs"),
                                   gamma_cutoff = 0.2,
                                   multi_sep = ";",
                                   verbose = TRUE) {
  which  <- match.arg(which)
  method <- match.arg(method)

  # Column resolution (be lenient to casing/synonyms)
  tf_col_in    <- .pick_col(df, c("tf","TF"), "TF")
  gene_col_in  <- .pick_col(df, c("gene_key","gene","target"), "gene_key")
  delta_col_in <- .pick_col(df, c("delta_link_score","d_edge_score","delta"), "delta_link_score")

  # Normalize column names to avoid rlang::.data inside dplyr verbs
  dat0 <- df[, c(tf_col_in, gene_col_in, delta_col_in)]
  names(dat0) <- c("TF", "gene_key", "dval")

  # Build weighted edges
  w <- switch(which,
              abs  = abs(dat0$dval),
              gain = pmax(dat0$dval, 0),
              loss = pmax(-dat0$dval, 0))
  dat <- tibble::tibble(TF = as.character(dat0$TF),
                        gene_key = as.character(dat0$gene_key),
                        w = as.numeric(w))
  dat <- dat[is.finite(dat$w) & dat$w > 0, , drop = FALSE]

  if (!nrow(dat)) {
    warning("No positive weights available for LDA (after '", which, "' transform). Returning NA assignments.")
    all_tfs <- unique(df[[tf_col_in]])
    return(tibble::tibble(!!tf_col_in := all_tfs,
                          topics_str = NA_character_,
                          gamma_max  = NA_real_))
  }

  # Prune ultra-rare terms by document frequency (across TFs)
  df_tbl <- dat[, c("TF","gene_key")]
  df_tbl <- unique(df_tbl)
  term_df <- as.data.frame(table(df_tbl$gene_key), stringsAsFactors = FALSE)
  names(term_df) <- c("gene_key","df")
  keep_terms <- term_df$gene_key[term_df$df >= min_df]
  dat <- dat[dat$gene_key %in% keep_terms, , drop = FALSE]

  if (!nrow(dat)) {
    warning("After min_df=", min_df, " filter, no terms remain. Returning NA assignments.")
    all_tfs <- unique(df[[tf_col_in]])
    return(tibble::tibble(!!tf_col_in := all_tfs,
                          topics_str = NA_character_,
                          gamma_max  = NA_real_))
  }

  # Dedup (TF, gene_key), sum weights
  dat <- stats::aggregate(w ~ TF + gene_key, data = dat, sum)

  # Build triplet matrix of integer counts (scaled from weights)
  docs  <- sort(unique(dat$TF))
  terms <- sort(unique(dat$gene_key))
  i <- match(dat$TF, docs)
  j <- match(dat$gene_key, terms)

  s <- stats::quantile(dat$w, 0.9, na.rm = TRUE)
  s <- if (is.finite(s) && s > 0) s else 1
  v <- pmax(1L, as.integer(round(10 * dat$w / s)))

  stm <- slam::simple_triplet_matrix(i, j, v,
                                     nrow = length(docs), ncol = length(terms),
                                     dimnames = list(docs, terms))

  # Fit LDA
  ctrl <- .tm_control(method = method, seed = seed)
  fit  <- topicmodels::LDA(stm, k = K, method = method, control = ctrl)

  post <- topicmodels::posterior(fit)  # list: topics (gamma), terms (phi)

  # Topics (per TF): coerce to long without rlang::.data
  gamma_mat <- as.matrix(post$topics)                # rows = docs (TF), cols = topics
  tf_vec    <- rownames(gamma_mat)
  topic_vec <- colnames(gamma_mat)
  topics_long <- tibble::tibble(
    !!tf_col_in := rep(tf_vec, times = ncol(gamma_mat)),
    topic       = rep(.safe_topic_int(topic_vec), each = nrow(gamma_mat)),
    gamma       = as.vector(gamma_mat)
  )

  # Per-TF max gamma (for diagnostics)
  # use base tapply to avoid dplyr NSE
  gamma_max_vals <- tapply(topics_long$gamma, topics_long[[tf_col_in]], max, na.rm = TRUE)
  gamma_max_tbl  <- tibble::tibble(!!tf_col_in := names(gamma_max_vals),
                                   gamma_max  = as.numeric(gamma_max_vals))

  # Multi-topic selection per TF using gamma >= cutoff
  keep_idx <- which(topics_long$gamma >= gamma_cutoff)
  if (length(keep_idx)) {
    t_keep <- topics_long[keep_idx, , drop = FALSE]
    # order by gamma desc within TF
    ord <- order(t_keep[[tf_col_in]], -t_keep$gamma, t_keep$topic)
    t_keep <- t_keep[ord, , drop = FALSE]

    # collapse topics per TF
    split_topics <- split(t_keep$topic, t_keep[[tf_col_in]])
    topics_str   <- vapply(split_topics, function(v) paste(v, collapse = multi_sep), character(1))
    multi_assign <- tibble::tibble(!!tf_col_in := names(topics_str),
                                   topics_str   = unname(topics_str))
  } else {
    multi_assign <- tibble::tibble(!!tf_col_in := character(0), topics_str = character(0))
  }

  # Ensure we include TFs with no topics >= cutoff (NA topics_str)
  all_tfs <- tibble::tibble(!!tf_col_in := docs)
  out <- dplyr::left_join(all_tfs, multi_assign, by = tf_col_in)
  out <- dplyr::left_join(out, gamma_max_tbl, by = tf_col_in)

  out
}


# =================================
# Public API: single table
# =================================

#' Annotate a delta-links table with LDA topic(s) per TF
#'
#' @param x Tibble/data.frame (filtered delta links) or CSV path.
#' @param K Integer, number of topics (default 20).
#' @param which One of 'abs' (default), 'gain', 'loss' to form weights from delta_link_score.
#' @param min_df Integer; minimum TF-document frequency for a gene term (default 2).
#' @param method 'VEM' (default) or 'Gibbs' (passed to topicmodels::LDA).
#' @param seed Integer RNG seed (default 1).
#' @param gamma_cutoff Numeric; keep all topics with gamma >= cutoff (default 0.2).
#' @param multi_sep String to separate multiple topic IDs (default ';').
#' @param save_csv Logical; if TRUE and `x` is a CSV path, write alongside input
#'        with suffix \code{_lda_K{K}.csv} (default TRUE).
#' @param out_file Optional override path when saving.
#' @param verbose Logical; default TRUE.
#'
#' @return Tibble identical to input plus one extra column \code{LDA_K{K}_topic}
#'         containing one or more topic IDs (semicolon-separated) or NA.
#' @export
#' @importFrom readr read_csv write_csv
#' @importFrom tibble as_tibble
#' @importFrom dplyr left_join rename
annotate_links_with_lda_topic <- function(x,
                                          K = 20,
                                          which = c("abs","gain","loss"),
                                          min_df = 2L,
                                          method = c("VEM","Gibbs"),
                                          seed = 1L,
                                          gamma_cutoff = 0.2,
                                          multi_sep = ";",
                                          save_csv = TRUE,
                                          out_file = NULL,
                                          verbose = TRUE) {
  df <- .load_table(x, verbose = verbose)

  # Resolve TF column (for joining)
  tf_col <- .pick_col(df, c("tf","TF"), "TF")

  # Fit and get per-TF multi-topic assignment
  A <- .fit_lda_assign_topics(df,
                              which = which, min_df = min_df,
                              K = K, seed = seed, method = method,
                              gamma_cutoff = gamma_cutoff, multi_sep = multi_sep,
                              verbose = verbose)

  topic_col <- paste0("LDA_K", K, "_topic")

  # Join back by TF; preserve original columns (avoid rlang::.data)
  A_small <- A[, c(tf_col, "topics_str"), drop = FALSE]
  names(A_small)[names(A_small) == "topics_str"] <- topic_col
  out <- dplyr::left_join(df, A_small, by = tf_col)

  # Save if requested and x was a path
  if (isTRUE(save_csv) && is.character(x) && length(x) == 1L) {
    if (is.null(out_file)) out_file <- .default_lda_out_path(x, K)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out, out_file)
    .ldalog("?wrote: ", out_file, verbose = verbose)
  }

  out
}


# =================================
# Public API: bulk / parallel
# =================================

#' Bulk/parallel LDA annotation over many filtered CSVs
#'
#' @param filtered_csvs Character vector of filtered CSV paths
#'        (e.g., \code{bulk$filtered_paths} from filter step).
#' @inheritParams annotate_links_with_lda_topic
#' @param parallel Logical; run via future.apply. Default TRUE.
#' @param plan Future plan; string or function. Default "multisession".
#' @param workers Integer; default parallel::detectCores()-1.
#' @param write_manifest Logical; write a manifest CSV mapping source->annotated (default TRUE).
#' @param manifest_file Optional path for manifest; if NULL, writes next to the first input as
#'        "delta_links_lda_manifest.csv".
#' @param verbose Logical; default TRUE.
#'
#' @return A list with:
#'   \itemize{
#'     \item results: named list of annotated tibbles
#'     \item assigned_paths: character vector of output CSVs
#'     \item assigned_dirs: unique directories
#'     \item manifest: tibble (source_filtered, annotated_csv, exists)
#'     \item manifest_path: path where manifest was written (or NULL)
#'   }
#' @export
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom parallel detectCores
#' @importFrom tibble tibble
#' @importFrom readr write_csv
annotate_links_with_lda_topic_bulk <- function(filtered_csvs,
                                               K = 20,
                                               which = c("abs","gain","loss"),
                                               min_df = 2L,
                                               method = c("VEM","Gibbs"),
                                               seed = 1L,
                                               gamma_cutoff = 0.2,
                                               multi_sep = ";",
                                               save_csv = TRUE,
                                               parallel = TRUE,
                                               plan = "multisession",
                                               workers = NULL,
                                               write_manifest = TRUE,
                                               manifest_file = NULL,
                                               verbose = TRUE) {
  if (!is.character(filtered_csvs) || length(filtered_csvs) == 0L) {
    cli::cli_abort("`filtered_csvs` must be a non-empty character vector of CSV paths.")
  }
  exist_mask <- file.exists(filtered_csvs)
  if (!all(exist_mask)) {
    missing <- filtered_csvs[!exist_mask]
    warning("These filtered CSVs do not exist and will be skipped:\n  - ",
            paste(missing, collapse = "\n  - "))
    filtered_csvs <- filtered_csvs[exist_mask]
  }
  if (length(filtered_csvs) == 0L) cli::cli_abort("No existing filtered CSVs to process.")

  assigned_paths <- vapply(filtered_csvs, function(p) .default_lda_out_path(p, K), character(1))
  assigned_dirs  <- unique(dirname(assigned_paths))

  # Parallel plan
  if (isTRUE(parallel)) {
    strategy <- (function(x) {
      if (is.character(x)) {
        ns <- asNamespace("future")
        if (exists(x, envir = ns, mode = "function")) {
          get(x, envir = ns)
        } else {
          future::multisession
        }
      } else if (is.function(x)) x else future::multisession
    })(plan)
    if (is.null(workers)) workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    .ldalog("Launching future plan=", if (is.character(plan)) plan else "custom",
            ", workers=", workers, verbose = verbose)
    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    future::plan(strategy, workers = workers)
  }

  run_one <- function(f) {
    annotate_links_with_lda_topic(
      f,
      K = K, which = which, min_df = min_df,
      method = method, seed = seed,
      gamma_cutoff = gamma_cutoff, multi_sep = multi_sep,
      save_csv = save_csv, verbose = verbose
    )
  }

  results <- if (isTRUE(parallel)) {
    future.apply::future_lapply(filtered_csvs, run_one, future.seed = TRUE)
  } else {
    lapply(filtered_csvs, run_one)
  }
  names(results) <- basename(filtered_csvs)

  manifest <- tibble::tibble(
    source_filtered = filtered_csvs,
    annotated_csv   = assigned_paths,
    exists          = file.exists(assigned_paths)
  )

  manifest_path <- NULL
  if (isTRUE(write_manifest)) {
    manifest_path <- if (is.null(manifest_file)) {
      file.path(dirname(filtered_csvs[[1]]), "delta_links_lda_manifest.csv")
    } else {
      manifest_file
    }
    dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(manifest, manifest_path)
    .ldalog("?wrote manifest: ", manifest_path, verbose = verbose)
  }

  list(
    results        = results,
    assigned_paths = assigned_paths,
    assigned_dirs  = assigned_dirs,
    manifest       = manifest,
    manifest_path  = manifest_path
  )
}
