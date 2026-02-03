#' JASPAR-like Position Frequency Matrix (PFM) / Position Probability Matrix (PPM) clustering utilities
#'
#' Alignment-based motif similarity using mean column-wise Pearson correlation
#' across overlaps (direct + reverse-complement). Includes pairwise metrics,
#' similarity matrix, and hclust + cut-to-K clustering.
#'
#' @details
#' Motifs are represented as 4 x L probability matrices (Position Probability Matrix; PPM)
#' with rows A,C,G,T. Position Frequency Matrices (PFMs; counts) are accepted and
#' converted to PPM internally.
#' Per-column correlation is Pearson over 4 values; if sd(v1)==0 OR sd(v2)==0,
#' correlation is defined as 0 for that column. Alignment offsets are tested
#' from -(w2-1) to (w1-1) inclusive. Best alignment chosen by:
#' (1) maximize cor, (2) maximize Ncor, (3) maximize overlap length w.
#'
#' Normalizations:
#'   Ncor1 = cor * (w / w1)
#'   Ncor2 = cor * (w / w2)
#'   Ncor  = cor * (w / max(w1, w2))
#'
#' @name utils_motif_clustering
NULL

# -----------------------------
# Parsing JASPAR PFM text
# -----------------------------

#' Parse JASPAR-like Position Frequency Matrix (PFM; counts) text into PPMs
#'
#' @param jaspar_text Character scalar: full text with one or more motifs in
#'   JASPAR PFM format. Also supports JASPAR-like blocks with four unlabeled
#'   numeric rows (assumed A,C,G,T order), as seen in HOCOMOCO exports.
#' @param strict Logical; if TRUE, stop on malformed motifs. If FALSE, skip
#'   malformed motifs with a warning.
#' @return Named list. Each element is a 4 x L numeric matrix of probabilities
#'   with rows A,C,G,T.
#' @export
parse_jaspar_pfm_text <- function(jaspar_text, strict = TRUE) {
  stopifnot(is.character(jaspar_text), length(jaspar_text) == 1L)

  lines <- unlist(strsplit(jaspar_text, "\n", fixed = TRUE), use.names = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  hdr_idx <- which(startsWith(lines, ">"))
  if (length(hdr_idx) == 0L) {
    stop("No motifs found: expected lines starting with '>'")
  }
  hdr_idx <- c(hdr_idx, length(lines) + 1L)

  motifs <- list()

  for (i in seq_len(length(hdr_idx) - 1L)) {
    start <- hdr_idx[i]
    end <- hdr_idx[i + 1L] - 1L

    header <- sub("^>", "", lines[start])
    # Typically: "MA0004.1\tArnt" or "MA0004.1 Arnt"
    header_parts <- strsplit(header, "[\t ]+", perl = TRUE)[[1]]
    motif_id <- header_parts[1]
    motif_name <- if (length(header_parts) >= 2L) header_parts[2] else motif_id
    motif_key <- paste0(motif_id, "|", motif_name)

    block <- lines[(start + 1L):end]
    if (length(block) < 4L) next

    bases <- c("A", "C", "G", "T")
    parse_nums <- function(x) {
      nums <- regmatches(
        x,
        gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", x)
      )[[1]]
      if (!length(nums)) return(numeric(0))
      as.numeric(nums)
    }

    labeled_present <- vapply(
      bases,
      function(b) any(grepl(paste0("^", b, "\\b"), block, ignore.case = TRUE)),
      logical(1)
    )

    bad <- FALSE
    mat_counts <- NULL

    if (all(labeled_present)) {
      mat_counts <- lapply(bases, function(b) {
        row_line <- block[grepl(paste0("^", b, "\\b"), block, ignore.case = TRUE)][1]
        nums <- parse_nums(row_line)
        if (!length(nums)) {
          if (isTRUE(strict)) {
            stop(sprintf("Missing numeric values for motif %s row %s", motif_key, b))
          }
          warning(sprintf("Skipping motif %s: missing numeric values for row %s", motif_key, b), call. = FALSE)
          bad <<- TRUE
          return(numeric(0))
        }
        nums
      })
    } else if (any(labeled_present)) {
      # Mixed labeled/unlabeled rows: treat as malformed.
      missing_bases <- bases[!labeled_present]
      if (isTRUE(strict)) {
        stop(sprintf("Missing base row(s) '%s' for motif %s", paste(missing_bases, collapse = ","), motif_key))
      }
      warning(
        sprintf("Skipping motif %s: missing base row(s) '%s'", motif_key, paste(missing_bases, collapse = ",")),
        call. = FALSE
      )
      bad <- TRUE
    } else {
      # Unlabeled 4-row format (e.g., HOCOMOCO exports).
      numeric_rows <- lapply(block, parse_nums)
      numeric_rows <- numeric_rows[vapply(numeric_rows, length, integer(1)) > 0L]
      if (length(numeric_rows) < 4L) {
        if (isTRUE(strict)) {
          stop(sprintf("Missing numeric rows for motif %s", motif_key))
        }
        warning(sprintf("Skipping motif %s: missing numeric rows", motif_key), call. = FALSE)
        bad <- TRUE
      } else {
        mat_counts <- numeric_rows[1:4]
      }
    }
    if (isTRUE(bad)) next

    lens <- vapply(mat_counts, length, integer(1))
    if (length(unique(lens)) != 1L) {
      if (isTRUE(strict)) {
        stop(sprintf("Inconsistent row lengths for motif %s: %s", motif_key, paste(lens, collapse = ",")))
      }
      warning(sprintf("Skipping motif %s: inconsistent row lengths (%s)", motif_key, paste(lens, collapse = ",")),
              call. = FALSE)
      next
    }

    counts <- do.call(rbind, mat_counts)
    rownames(counts) <- bases
    probs <- .pwm_counts_to_probs(counts)
    motifs[[motif_key]] <- probs
  }

  motifs
}

#' Read a JASPAR PFM file into PPMs
#'
#' @param file_path Path to a JASPAR PFM file.
#' @return Named list of 4 x L probability matrices (PPMs).
#' @export
read_jaspar_pfm <- function(file_path, strict = TRUE) {
  if (!grepl("^https?://", file_path) && !file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  jaspar_text <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
  parse_jaspar_pfm_text(jaspar_text, strict = strict)
}

#' Load PFMs for a known motif database (converted to PPMs)
#'
#' @param db Character; supported values: "JASPAR2024", "HOCOMOCOv13".
#' @return Named list of 4 x L probability matrices (PPMs).
#' @export
load_motif_pfms <- function(db) {
  db <- as.character(db)
  url <- switch(
    db,
    "JASPAR2024" = "https://raw.githubusercontent.com/oncologylab/episcope-data/refs/heads/main/TFBS_db/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
    "HOCOMOCOv13" = "https://raw.githubusercontent.com/oncologylab/episcope-data/refs/heads/main/TFBS_db/H13CORE_jaspar_format.txt",
    stop("Unsupported db: ", db)
  )
  strict <- !identical(db, "HOCOMOCOv13")
  read_jaspar_pfm(url, strict = strict)
}

#' @rdname load_motif_pfms
#' @export
load_motif_pwms <- function(db) {
  cli::cli_warn("`load_motif_pwms()` is deprecated; use `load_motif_pfms()` (PFM -> PPM).")
  load_motif_pfms(db)
}

# -----------------------------
# Position Probability Matrix (PPM) utilities
# -----------------------------

.pwm_counts_to_probs <- function(counts) {
  stopifnot(is.matrix(counts), nrow(counts) == 4L)
  col_sums <- colSums(counts)
  probs <- counts
  for (j in seq_len(ncol(counts))) {
    if (col_sums[j] <= 0) {
      probs[, j] <- 0.25
    } else {
      probs[, j] <- counts[, j] / col_sums[j]
    }
  }
  probs
}

#' Reverse-complement a Position Probability Matrix (PPM; rows A,C,G,T)
#'
#' @param pwm 4 x L numeric matrix (PPM)
#' @return 4 x L numeric matrix reversed and complemented
#' @export
pwm_revcomp <- function(pwm) {
  stopifnot(is.matrix(pwm), nrow(pwm) == 4L)
  bases <- rownames(pwm)
  if (is.null(bases)) stop("PWM must have rownames A,C,G,T")
  idx <- match(c("A", "C", "G", "T"), bases)
  if (anyNA(idx)) stop("PWM rownames must include A,C,G,T")
  pwm_std <- pwm[idx, , drop = FALSE]
  pwm_rev <- pwm_std[, ncol(pwm_std):1, drop = FALSE]
  pwm_rc <- pwm_rev[c(4, 3, 2, 1), , drop = FALSE]
  rownames(pwm_rc) <- c("A", "C", "G", "T")
  pwm_rc
}

.prepare_pwm <- function(pwm) {
  stopifnot(is.matrix(pwm), nrow(pwm) == 4L)
  bases <- rownames(pwm)
  if (is.null(bases) || anyNA(match(c("A", "C", "G", "T"), bases))) {
    rownames(pwm) <- c("A", "C", "G", "T")
  } else {
    pwm <- pwm[match(c("A", "C", "G", "T"), bases), , drop = FALSE]
  }
  pwm
}

.col_corr4_fast <- function(v1, v2) {
  s1 <- stats::sd(v1)
  s2 <- stats::sd(v2)
  if (s1 == 0 || s2 == 0) return(0)
  n <- 4
  m1 <- mean(v1)
  m2 <- mean(v2)
  cov12 <- sum((v1 - m1) * (v2 - m2)) / (n - 1)
  as.numeric(cov12 / (s1 * s2))
}

# -----------------------------
# Pairwise alignment + metrics
# -----------------------------

#' Compute best alignment similarity for a pair of motifs (PPMs)
#'
#' @param pwm1 4 x L1 probability matrix (PPM)
#' @param pwm2 4 x L2 probability matrix (PPM)
#' @param min_overlap Integer >=1. Minimum overlap length for alignments.
#' @return Named list with cor, Ncor, Ncor1, Ncor2, w1, w2, w, strand, offset.
#' @export
best_pwm_alignment <- function(pwm1, pwm2, min_overlap = 1L) {
  stopifnot(is.matrix(pwm1), is.matrix(pwm2), nrow(pwm1) == 4L, nrow(pwm2) == 4L)
  stopifnot(is.numeric(min_overlap), length(min_overlap) == 1L, min_overlap >= 1L)

  pwm1 <- .prepare_pwm(pwm1)
  pwm2 <- .prepare_pwm(pwm2)
  w1 <- ncol(pwm1)
  w2 <- ncol(pwm2)
  offsets <- seq.int(from = -(w2 - 1L), to = (w1 - 1L), by = 1L)

  eval_one <- function(p2, strand_label) {
    best <- list(
      cor = -Inf, Ncor = -Inf, Ncor1 = -Inf, Ncor2 = -Inf,
      w1 = w1, w2 = w2, w = 0L, strand = strand_label, offset = 0L
    )

    for (off in offsets) {
      start1 <- 0L
      end1 <- w1 - 1L
      start2 <- off
      end2 <- off + (w2 - 1L)

      ov_start <- max(start1, start2)
      ov_end <- min(end1, end2)
      w <- ov_end - ov_start + 1L
      if (w < min_overlap) next

      idx1 <- (ov_start - start1) + 1L
      idx2 <- (ov_start - start2) + 1L
      cols1 <- idx1:(idx1 + w - 1L)
      cols2 <- idx2:(idx2 + w - 1L)

      cors <- numeric(w)
      for (k in seq_len(w)) {
        cors[k] <- .col_corr4_fast(pwm1[, cols1[k]], p2[, cols2[k]])
      }
      cor_mean <- mean(cors)

      Ncor1 <- cor_mean * (w / w1)
      Ncor2 <- cor_mean * (w / w2)
      Ncor <- cor_mean * (w / max(w1, w2))

      is_better <- (cor_mean > best$cor) ||
        (cor_mean == best$cor && Ncor > best$Ncor) ||
        (cor_mean == best$cor && Ncor == best$Ncor && w > best$w)

      if (is_better) {
        best$cor <- cor_mean
        best$Ncor <- Ncor
        best$Ncor1 <- Ncor1
        best$Ncor2 <- Ncor2
        best$w <- as.integer(w)
        best$strand <- strand_label
        best$offset <- as.integer(off)
      }
    }

    best
  }

  best_d <- eval_one(pwm2, "D")
  best_r <- eval_one(pwm_revcomp(pwm2), "R")

  pick_r <- (best_r$cor > best_d$cor) ||
    (best_r$cor == best_d$cor && best_r$Ncor > best_d$Ncor) ||
    (best_r$cor == best_d$cor && best_r$Ncor == best_d$Ncor && best_r$w > best_d$w)

  if (pick_r) best_r else best_d
}

# -----------------------------
# All-pairs comparisons + clustering
# -----------------------------

#' Compute all pairwise comparisons (upper triangle) for a set of motifs (PPMs)
#'
#' @param motifs Named list of 4xL probability matrices.
#' @param min_overlap Minimum overlap length for alignments.
#' @param cores Integer number of workers (default: all available cores).
#' @param verbose Logical; emit concise progress messages (default TRUE).
#' @return data.frame with pairwise metrics similar to JASPAR pairwise_comparisons.tab.
#' @export
pairwise_comparisons <- function(motifs, min_overlap = 1L, cores = 1L, verbose = TRUE) {
  stopifnot(is.list(motifs), length(motifs) >= 2L)
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) cores <- 1L

  ids <- names(motifs)
  n <- length(ids)
  if (isTRUE(verbose)) .log_inform("Pairwise PWM comparisons: {n} motifs.")

  # Prepare motif ordering/rownames once
  motifs <- lapply(motifs, .prepare_pwm)

  pair_idx <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      pair_idx[[length(pair_idx) + 1L]] <- c(i, j)
    }
  }

  worker_fun <- function(idx_vec) {
    out <- vector("list", length(idx_vec))
    for (k in seq_along(idx_vec)) {
      ij <- idx_vec[[k]]
      i <- ij[1]
      j <- ij[2]
      m1 <- motifs[[i]]
      m2 <- motifs[[j]]

      best <- best_pwm_alignment(m1, m2, min_overlap = min_overlap)

      key1 <- ids[i]
      key2 <- ids[j]
      name1 <- if (grepl("\\|", key1)) strsplit(key1, "\\|", fixed = FALSE)[[1]][2] else key1
      name2 <- if (grepl("\\|", key2)) strsplit(key2, "\\|", fixed = FALSE)[[1]][2] else key2

      out[[k]] <- data.frame(
        id1 = key1,
        id2 = key2,
        name1 = name1,
        name2 = name2,
        cor = best$cor,
        Ncor = best$Ncor,
        Ncor1 = best$Ncor1,
        Ncor2 = best$Ncor2,
        w1 = best$w1,
        w2 = best$w2,
        w = best$w,
        strand = best$strand,
        offset = best$offset,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, out)
  }

  if (cores == 1L) {
    res <- worker_fun(pair_idx)
  } else {
    # Split pairs into chunks; avoid printing from workers.
    n_pairs <- length(pair_idx)
    n_chunks <- min(cores, n_pairs)
    chunk_id <- rep(seq_len(n_chunks), length.out = n_pairs)
    chunks <- split(pair_idx, chunk_id)

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      res_list <- parallel::parLapply(cl, chunks, worker_fun)
    } else {
      res_list <- parallel::mclapply(chunks, worker_fun, mc.cores = cores)
    }
    res <- do.call(rbind, res_list)
  }

  res
}

#' Build an NxN similarity matrix from pairwise comparisons (using Ncor)
#'
#' @param motifs Named list of motifs (defines ordering).
#' @param pw data.frame from pairwise_comparisons().
#' @return numeric matrix NxN with diag=1 and symmetric entries = Ncor.
#' @export
build_similarity_matrix <- function(motifs, pw) {
  ids <- names(motifs)
  n <- length(ids)
  sim <- matrix(0, nrow = n, ncol = n, dimnames = list(ids, ids))
  diag(sim) <- 1

  for (k in seq_len(nrow(pw))) {
    i <- match(pw$id1[k], ids)
    j <- match(pw$id2[k], ids)
    if (is.na(i) || is.na(j)) next
    sim[i, j] <- pw$Ncor[k]
    sim[j, i] <- pw$Ncor[k]
  }

  sim
}

#' Cut an hclust tree to target K clusters (by searching cut height)
#'
#' @param hc hclust object.
#' @param target_k Integer number of clusters desired.
#' @return Integer vector of cluster labels aligned to hc$order labels.
#' @export
cutree_to_k <- function(hc, target_k) {
  stopifnot(inherits(hc, "hclust"))
  stopifnot(is.numeric(target_k), length(target_k) == 1L, target_k >= 2L)

  heights <- sort(unique(hc$height))
  heights <- c(-Inf, heights, Inf)

  best_labels <- NULL
  best_diff <- Inf

  for (h in heights) {
    labs <- stats::cutree(hc, h = h)
    k <- length(unique(labs))
    diff <- abs(k - target_k)
    if (diff < best_diff) {
      best_diff <- diff
      best_labels <- labs
      if (diff == 0L) break
    }
  }

  best_labels
}

#' Run a JASPAR-like PFM/PPM clustering workflow
#'
#' @param jaspar_text Full JASPAR PFM text.
#' @param min_overlap Minimum overlap length for alignments.
#' @param linkage Linkage method passed to hclust().
#' @param target_clusters Target number of clusters (K).
#' @param cores Number of workers for pairwise comparisons.
#' @param verbose Logical; emit concise progress messages.
#' @return List with motifs, pairwise, similarity, hclust, clusters, clusters_table.
#' @export
run_jaspar_like_clustering <- function(
  jaspar_text = NULL,
  db = NULL,
  min_overlap = 1L,
  linkage = "average",
  target_clusters = 233L,
  cores = 1L,
  verbose = TRUE
) {
  if (is.null(jaspar_text) && is.null(db)) {
    stop("Provide either `jaspar_text` or `db`.")
  }
  motifs <- if (!is.null(jaspar_text)) {
    parse_jaspar_pfm_text(jaspar_text)
  } else {
    load_motif_pfms(db)
  }
  if (isTRUE(verbose)) .log_inform("Parsed {length(motifs)} motifs.")
  if (length(motifs) < 2L) stop("Need at least 2 motifs")

  pw <- pairwise_comparisons(motifs, min_overlap = min_overlap, cores = cores, verbose = verbose)
  sim <- build_similarity_matrix(motifs, pw)
  dist_mat <- stats::as.dist(1 - sim)

  hc <- stats::hclust(dist_mat, method = linkage)
  cl <- cutree_to_k(hc, target_clusters)

  ids <- names(motifs)
  cluster_ids <- split(ids, cl)
  cluster_df <- data.frame(
    cluster = sprintf("cluster_%03d", seq_along(cluster_ids)),
    id = vapply(cluster_ids, function(x) paste(x, collapse = ","), character(1)),
    name = vapply(cluster_ids, function(x) {
      nm <- vapply(x, function(k) {
        parts <- strsplit(k, "\\|", fixed = FALSE)[[1]]
        if (length(parts) >= 2L) parts[2] else parts[1]
      }, character(1))
      paste(unique(nm), collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) .log_inform("QC plots (global distributions): done.")
  list(
    motifs = motifs,
    pairwise = pw,
    similarity = sim,
    hclust = hc,
    clusters = cl,
    clusters_table = cluster_df
  )
}

.get_ref_similarity <- function(db, cache_dir, cores = 1L, verbose = TRUE) {
  if (is.null(db)) return(NULL)
  cache_dir <- if (is.null(cache_dir)) file.path(tempdir(), "motif_clustering_cache") else cache_dir
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(cache_dir, sprintf("ref_similarity_%s.rds", db))
  if (file.exists(cache_path)) {
    if (isTRUE(verbose)) .log_inform("Using cached reference similarity: {cache_path}")
    return(readRDS(cache_path))
  }
  if (isTRUE(verbose)) .log_inform("Computing reference similarity for {db}.")
  motifs <- load_motif_pfms(db)
  pw <- pairwise_comparisons(motifs, min_overlap = 1L, cores = cores, verbose = verbose)
  sim <- build_similarity_matrix(motifs, pw)
  saveRDS(sim, cache_path)
  sim
}

# -----------------------------
# Performance guidance / Rcpp skeleton (comments)
# -----------------------------

# For large motif sets (~1000), consider:
# - Parallelizing pairwise comparisons (already supported via cores)
# - Rcpp acceleration for best_pwm_alignment inner loops
# - Precomputing per-column mean/sd and using vectorized dot products
#
# Optional Rcpp skeleton (outline only):
# // [[Rcpp::export]]
# Rcpp::List best_pwm_alignment_cpp(Rcpp::NumericMatrix pwm1, Rcpp::NumericMatrix pwm2, int min_overlap);
#
# The C++ implementation should:
# - Precompute column means/sds (4xL).
# - Iterate offsets and compute correlation via fixed-n formula.
# - Apply tie-breaks: cor, then Ncor, then w.

# -----------------------------
# Footprint-driven motif co-occurrence clustering
# -----------------------------

#' Build motif-by-peak co-occurrence stats from aligned footprints
#'
#' @param fp_annotation Tibble/data.frame with `fp_peak` and `motifs` columns.
#' @param id_map Optional tibble/data.frame with `peak_ID` and `group_size`.
#' @param peak_col Character; aligned peak column name in fp_annotation.
#' @param motif_col Character; motif column name in fp_annotation.
#' @param min_peak_support Numeric; minimum weighted peak support per motif.
#' @param weight_by_group_size Logical; weight by id_map$group_size if present.
#' @param verbose Logical; emit concise progress messages.
#' @return List with motif_counts, pair_counts, motifs_keep, peaks_n, motifs_n.
#' @export
build_motif_peak_stats <- function(
  fp_annotation,
  id_map = NULL,
  peak_col = "fp_peak",
  motif_col = "motifs",
  min_peak_support = 100,
  weight_by_group_size = TRUE,
  verbose = TRUE
) {
  stopifnot(is.data.frame(fp_annotation))
  if (!all(c(peak_col, motif_col) %in% names(fp_annotation))) {
    stop("fp_annotation must contain columns: ", peak_col, ", ", motif_col)
  }

  dt <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::as.data.table(fp_annotation[, c(peak_col, motif_col), drop = FALSE])
  } else {
    fp_annotation[, c(peak_col, motif_col), drop = FALSE]
  }

  if (isTRUE(verbose)) {
    .log_inform("Motif co-occurrence: building motif-by-peak index.")
  }

  if (inherits(dt, "data.table")) {
    data.table::setnames(dt, c(peak_col, motif_col), c("fp_peak", "motif"))
    dt <- unique(dt, by = c("fp_peak", "motif"))
  } else {
    names(dt) <- c("fp_peak", "motif")
    dt <- unique(dt)
  }

  if (isTRUE(weight_by_group_size) && !is.null(id_map)) {
    if (all(c("peak_ID", "group_size") %in% names(id_map))) {
      wmap <- unique(id_map[, c("peak_ID", "group_size"), drop = FALSE])
      names(wmap) <- c("fp_peak", "group_size")
      if (inherits(dt, "data.table") && requireNamespace("data.table", quietly = TRUE)) {
        wmap <- data.table::as.data.table(wmap)
        dt <- wmap[dt, on = "fp_peak"]
      } else {
        dt <- merge(dt, wmap, by = "fp_peak", all.x = TRUE)
      }
    }
  }

  if (!("group_size" %in% names(dt))) {
    dt$group_size <- 1
  } else {
    dt$group_size[is.na(dt$group_size)] <- 1
  }

  if (inherits(dt, "data.table")) {
    motif_counts <- dt[, .(weight = sum(group_size)), by = motif]
  } else {
    motif_counts <- aggregate(group_size ~ motif, dt, sum)
    names(motif_counts) <- c("motif", "weight")
  }

  motifs_keep <- motif_counts$motif[motif_counts$weight >= min_peak_support]
  if (isTRUE(verbose)) {
    .log_inform("Motifs kept: {length(motifs_keep)} (min_peak_support={min_peak_support}).")
  }

  dt <- dt[dt$motif %in% motifs_keep, , drop = FALSE]

  if (inherits(dt, "data.table")) {
    peak_groups <- dt[, .(motifs = list(sort(unique(motif))), weight = max(group_size)), by = fp_peak]
  } else {
    peak_split <- split(dt, dt$fp_peak)
    peak_groups <- lapply(peak_split, function(x) {
      list(motifs = sort(unique(x$motif)), weight = max(x$group_size))
    })
    peak_groups <- data.frame(
      fp_peak = names(peak_groups),
      motifs = I(lapply(peak_groups, `[[`, "motifs")),
      weight = vapply(peak_groups, `[[`, numeric(1), "weight"),
      stringsAsFactors = FALSE
    )
  }

  make_pairs <- function(motifs, w) {
    if (length(motifs) < 2L) return(NULL)
    cmb <- utils::combn(motifs, 2L)
    data.frame(
      motif1 = cmb[1, ],
      motif2 = cmb[2, ],
      weight = rep(w, ncol(cmb)),
      stringsAsFactors = FALSE
    )
  }

  if (inherits(peak_groups, "data.table")) {
    pair_list <- mapply(make_pairs, peak_groups$motifs, peak_groups$weight, SIMPLIFY = FALSE)
    pair_dt <- data.table::rbindlist(pair_list, use.names = TRUE, fill = TRUE)
    if (nrow(pair_dt)) {
      pair_counts <- pair_dt[, .(weight = sum(weight)), by = .(motif1, motif2)]
    } else {
      pair_counts <- data.table::data.table(motif1 = character(0), motif2 = character(0), weight = numeric(0))
    }
  } else {
    pair_list <- mapply(make_pairs, peak_groups$motifs, peak_groups$weight, SIMPLIFY = FALSE)
    pair_dt <- do.call(rbind, pair_list)
    if (is.null(pair_dt)) {
      pair_counts <- data.frame(motif1 = character(0), motif2 = character(0), weight = numeric(0))
    } else {
      pair_counts <- aggregate(weight ~ motif1 + motif2, pair_dt, sum)
    }
  }

  list(
    motif_counts = motif_counts,
    pair_counts = pair_counts,
    motifs_keep = motifs_keep,
    peaks_n = length(unique(dt$fp_peak)),
    motifs_n = length(motifs_keep)
  )
}

#' Compute motif co-occurrence similarity matrix
#'
#' @param motif_stats Output of build_motif_peak_stats().
#' @param method Character; "jaccard" or "cosine".
#' @param cores Integer number of workers (default 1L).
#' @return Numeric similarity matrix with diag=1.
#' @export
compute_motif_similarity <- function(motif_stats, method = c("jaccard", "cosine"), cores = 1L) {
  method <- match.arg(method)
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  motifs <- motif_stats$motifs_keep
  n <- length(motifs)
  sim <- matrix(0, nrow = n, ncol = n, dimnames = list(motifs, motifs))
  diag(sim) <- 1

  counts <- motif_stats$motif_counts
  if (!all(c("motif", "weight") %in% names(counts))) {
    stop("motif_counts must have columns motif, weight.")
  }
  count_map <- stats::setNames(counts$weight, counts$motif)

  pairs <- motif_stats$pair_counts
  if (!nrow(pairs)) return(sim)

  compute_block <- function(idx) {
    m1 <- pairs$motif1[idx]
    m2 <- pairs$motif2[idx]
    w12 <- pairs$weight[idx]
    c1 <- count_map[m1]
    c2 <- count_map[m2]
    ok <- is.finite(c1) & is.finite(c2) & c1 > 0 & c2 > 0
    if (any(ok)) {
      w12 <- w12[ok]
      m1 <- m1[ok]
      m2 <- m2[ok]
      c1 <- c1[ok]
      c2 <- c2[ok]
      s <- if (method == "jaccard") {
        w12 / (c1 + c2 - w12)
      } else {
        w12 / sqrt(c1 * c2)
      }
      data.frame(m1 = m1, m2 = m2, s = s, stringsAsFactors = FALSE)
    } else {
      data.frame(m1 = character(0), m2 = character(0), s = numeric(0))
    }
  }

  if (cores == 1L) {
    df <- compute_block(seq_len(nrow(pairs)))
  } else {
    n_pairs <- nrow(pairs)
    idx <- split(seq_len(n_pairs), rep_len(seq_len(min(cores, n_pairs)), n_pairs))
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      df_list <- parallel::parLapply(cl, idx, compute_block)
    } else {
      df_list <- parallel::mclapply(idx, compute_block, mc.cores = cores)
    }
    df <- do.call(rbind, df_list)
  }

  if (nrow(df)) {
    sim[cbind(df$m1, df$m2)] <- df$s
    sim[cbind(df$m2, df$m1)] <- df$s
  }

  sim
}

.save_plot_safe <- function(path, plot_fun, width = 8, height = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- plot_fun(ggplot2 = TRUE)
    ggplot2::ggsave(path, p, width = width, height = height, limitsize = FALSE)
  } else {
    grDevices::pdf(path, width = width, height = height)
    plot_fun(ggplot2 = FALSE)
    grDevices::dev.off()
  }
}

.theme_pub <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, hjust = 1),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.3)
    )
}

.extract_jaspar_id <- function(x) {
  m <- regexpr("MA[0-9]+\\.[0-9]+", x)
  out <- ifelse(m == -1L, NA_character_, regmatches(x, m))
  out
}

.extract_motif_key <- function(x, db = NULL) {
  if (!is.character(x) || !length(x)) return(x)
  db <- if (!is.null(db) && nzchar(db)) toupper(as.character(db)) else NA_character_
  if (identical(db, "JASPAR2024")) {
    return(.extract_jaspar_id(x))
  }
  if (identical(db, "HOCOMOCOv13")) {
    return(x)
  }
  # Fallback: if JASPAR-style IDs are present, use them; otherwise keep raw IDs.
  jaspar_ids <- .extract_jaspar_id(x)
  if (any(!is.na(jaspar_ids))) return(jaspar_ids)
  x
}

.make_tf_map_from_ids <- function(ids) {
  has_name <- grepl("\\|", ids)
  if (!any(has_name)) return(NULL)
  parts <- strsplit(ids[has_name], "\\|", fixed = FALSE)
  id_only <- vapply(parts, function(x) x[1], character(1))
  name_only <- vapply(parts, function(x) x[2], character(1))
  stats::setNames(name_only, id_only)
}

.normalize_similarity_ids <- function(sim) {
  if (is.null(rownames(sim))) return(list(sim = sim, tf_map = NULL))
  tf_map <- .make_tf_map_from_ids(rownames(sim))
  ids <- sub("\\|.*$", "", rownames(sim))
  rownames(sim) <- ids
  colnames(sim) <- ids
  if (anyDuplicated(ids)) {
    cli::cli_warn("ref_similarity has duplicated IDs after normalization; keeping first occurrence.")
    keep <- !duplicated(ids)
    sim <- sim[keep, keep, drop = FALSE]
  }
  list(sim = sim, tf_map = tf_map)
}

.make_cache_key <- function(prefix, ...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(as.character(parts))]
  key <- paste(parts, collapse = "__")
  key <- gsub("[^A-Za-z0-9]+", "_", key)
  paste0(prefix, "_", key)
}

.choose_best_k <- function(hc, dist_obj, k_grid = NULL, verbose = TRUE) {
  n <- length(hc$labels)
  if (is.null(k_grid) || !length(k_grid)) {
    max_k <- min(600L, n - 1L)
    if (max_k < 20L) {
      k_grid <- seq.int(2L, max_k, by = 1L)
    } else {
      k_grid <- seq.int(20L, max_k, by = 5L)
    }
  } else {
    k_grid <- unique(as.integer(k_grid))
    k_grid <- k_grid[k_grid >= 2L & k_grid < n]
  }
  if (!length(k_grid)) return(list(k = NA_integer_, scores = NULL))
  if (!requireNamespace("cluster", quietly = TRUE)) {
    if (isTRUE(verbose)) {
      cli::cli_warn("cluster package not available; using median K from grid.")
    }
    return(list(k = k_grid[ceiling(length(k_grid) / 2)], scores = NULL))
  }
  sil_scores <- vapply(k_grid, function(k) {
    cl <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(cl, dist_obj)
    mean(sil[, "sil_width"])
  }, numeric(1))
  k_df <- data.frame(k = k_grid, silhouette = sil_scores)
  k_info <- .k_curve_with_second_derivative(k_df)
  best_k <- .select_best_k_from_scores(k_info$curve, k_info$d2)
  list(k = best_k, scores = k_info$curve, d2 = k_info$d2)
}

.cluster_peak_summary <- function(fp_annotation, cl_map, top_n = 5L) {
  if (!is.data.frame(fp_annotation)) return(NULL)
  if (!all(c("fp_peak", "motifs") %in% names(fp_annotation))) return(NULL)
  if (!all(c("motif", "cluster_id") %in% names(cl_map))) return(NULL)

  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(fp_annotation[, c("fp_peak", "motifs")])
    data.table::setnames(dt, c("fp_peak", "motifs"), c("fp_peak", "motif"))
    dt <- unique(dt, by = c("fp_peak", "motif"))
    cmap <- data.table::as.data.table(cl_map)
    dt <- cmap[dt, on = "motif", nomatch = 0]
    if (!nrow(dt)) return(NULL)
    dt_peak <- dt[, .(n_motifs = data.table::uniqueN(motif)), by = .(cluster_id, fp_peak)]
    summary <- dt_peak[, .(
      n_peaks = .N,
      shared_peaks = sum(n_motifs >= 2L),
      avg_motifs_per_peak = mean(n_motifs)
    ), by = cluster_id]
    summary[, pct_shared_peaks := ifelse(n_peaks > 0, shared_peaks / n_peaks, 0)]
    top_peaks <- dt_peak[order(-n_motifs)][, .SD[seq_len(min(.N, top_n))], by = cluster_id]
    top_peaks[, peak_label := paste0(fp_peak, "(", n_motifs, ")")]
    top_peaks <- top_peaks[, .(top_shared_peaks = paste(peak_label, collapse = ",")), by = cluster_id]
    summary <- merge(summary, top_peaks, by = "cluster_id", all.x = TRUE)
    return(as.data.frame(summary))
  }

  dt <- unique(fp_annotation[, c("fp_peak", "motifs")])
  names(dt) <- c("fp_peak", "motif")
  dt <- merge(dt, cl_map, by = "motif")
  if (!nrow(dt)) return(NULL)
  dt_peak <- aggregate(motif ~ cluster_id + fp_peak, dt, function(x) length(unique(x)))
  names(dt_peak)[names(dt_peak) == "motif"] <- "n_motifs"
  summary <- aggregate(n_motifs ~ cluster_id, dt_peak, function(x) {
    c(n_peaks = length(x), shared_peaks = sum(x >= 2L), avg_motifs_per_peak = mean(x))
  })
  summary <- data.frame(
    cluster_id = summary$cluster_id,
    n_peaks = summary$n_motifs[, "n_peaks"],
    shared_peaks = summary$n_motifs[, "shared_peaks"],
    avg_motifs_per_peak = summary$n_motifs[, "avg_motifs_per_peak"]
  )
  summary$pct_shared_peaks <- ifelse(summary$n_peaks > 0, summary$shared_peaks / summary$n_peaks, 0)
  dt_peak <- dt_peak[order(-dt_peak$n_motifs), ]
  top_list <- split(dt_peak, dt_peak$cluster_id)
  top_peaks <- lapply(top_list, function(x) {
    head_x <- head(x, top_n)
    paste0(head_x$fp_peak, "(", head_x$n_motifs, ")", collapse = ",")
  })
  summary$top_shared_peaks <- unname(top_peaks[as.character(summary$cluster_id)])
  summary
}

.make_cluster_name <- function(motif_ids, tf_names, mode = c("rich", "compact"), max_parts = 4L, max_chars = 45L) {
  mode <- match.arg(mode)
  tf_names <- tf_names[!is.na(tf_names) & nzchar(tf_names)]
  if (!length(tf_names)) return(NA_character_)
  tf_names <- unique(unlist(strsplit(tf_names, ",", fixed = TRUE)))
  tf_names <- trimws(tf_names)
  tf_names <- tf_names[nzchar(tf_names)]
  tf_names <- vapply(tf_names, function(x) strsplit(x, "::", fixed = TRUE)[[1]][1], character(1))
  tf_upper <- toupper(tf_names)

  if (identical(mode, "compact")) {
    max_parts <- min(max_parts, 3L)
    max_chars <- min(max_chars, 35L)
  } else {
    max_parts <- max(max_parts, 5L)
    max_chars <- max(max_chars, 55L)
  }

  homeobox_set <- c("HOX", "PAX", "NKX", "LHX", "DLX", "ALX", "EMX", "EN", "GSX",
                    "RAX", "VSX", "PRRX", "LBX", "LMX", "MNX", "SHOX", "TLX", "GBX")

  family_key <- function(x) {
    if (grepl("^ZNF", x)) return("ZNF")
    if (grepl("^KLF", x)) return("KLF")
    if (grepl("^E2F", x)) return("E2F")
    if (grepl("^FOX", x)) return("FOX")
    if (grepl("^NR", x)) return(x)
    if (grepl("^HOX", x) || grepl("^POU6F", x) || grepl("^POU", x)) return("HOX")
    if (grepl("^ALX", x)) return("ALX")
    if (grepl("^DLX", x)) return("DLX")
    if (grepl("^LHX", x)) return("LHX")
    if (grepl("^NKX", x)) return("NKX")
    if (grepl("^PAX", x)) return("PAX")
    if (grepl("^POU", x)) return("POU")
    if (grepl("^PRRX", x)) return("PRRX")
    if (grepl("^RAX", x)) return("RAX")
    if (grepl("^VSX", x)) return("VSX")
    if (grepl("^VAX", x)) return("VAX")
    if (grepl("^EMX", x)) return("EMX")
    if (grepl("^EN", x)) return("EN")
    if (grepl("^GBX", x)) return("GBX")
    if (grepl("^GSX", x)) return("GSX")
    if (grepl("^LBX", x)) return("LBX")
    if (grepl("^LMX", x)) return("LMX")
    if (grepl("^MEOX", x)) return("MEOX")
    if (grepl("^MNX", x)) return("MNX")
    if (grepl("^SHOX", x)) return("SHOX")
    if (grepl("^TLX", x)) return("TLX")
    if (grepl("^UNCX", x)) return("UNCX")
    x
  }

  families <- vapply(tf_upper, family_key, character(1))
  fam_tab <- sort(table(families), decreasing = TRUE)
  fam_keep <- names(fam_tab)[seq_len(min(max_parts, length(fam_tab)))]
  if (identical(mode, "rich")) {
    homeobox_hits <- intersect(names(fam_tab), homeobox_set)
    if (length(homeobox_hits) >= 3L) {
      fam_keep <- unique(c("HOMEObox", homeobox_hits, fam_keep))
    }
    fam_keep <- unique(fam_keep)[seq_len(min(max_parts, length(unique(fam_keep))))]
  }
  name_raw <- paste0("C_", paste(fam_keep, collapse = ";"))
  if (nchar(name_raw) > max_chars) {
    name_raw <- paste0(substr(name_raw, 1L, max_chars - 2L), "..")
  }
  name_raw
}

.apply_cluster_to_motif_db <- function(motif_db, cluster_df, db = NULL) {
  if (is.null(motif_db) || !is.data.frame(motif_db)) return(motif_db)
  if (is.null(cluster_df) || !nrow(cluster_df)) return(motif_db)
  if (!"motif" %in% names(motif_db)) return(motif_db)
  if (!all(c("cluster_id", "cluster_name", "id") %in% names(cluster_df))) return(motif_db)

  split_ids <- strsplit(cluster_df$id, ",", fixed = TRUE)
  motif_vec <- unlist(split_ids, use.names = FALSE)
  if (!length(motif_vec)) return(motif_db)
  cluster_id <- rep(cluster_df$cluster_id, lengths(split_ids))
  cluster_name <- rep(cluster_df$cluster_name, lengths(split_ids))
  map <- data.frame(
    motif_key = motif_vec,
    sub_cluster = as.integer(cluster_id),
    sub_cluster_name = as.character(cluster_name),
    stringsAsFactors = FALSE
  )
  map <- map[!duplicated(map$motif_key), , drop = FALSE]
  motif_key <- .extract_motif_key(motif_db$motif, db = db)
  motif_key[is.na(motif_key)] <- motif_db$motif[is.na(motif_key)]
  motif_db$motif_key <- motif_key
  out <- dplyr::left_join(motif_db, map, by = "motif_key")
  out$motif_key <- NULL
  out
}

.save_similarity_heatmap <- function(sim_mat, clusters, out_path, title, show_labels = FALSE) {
  if (is.null(sim_mat) || !is.matrix(sim_mat)) return(invisible(NULL))
  if (is.null(clusters)) return(invisible(NULL))
  cl_vec <- clusters[rownames(sim_mat)]
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann <- data.frame(cluster = factor(cl_vec))
    rownames(ann) <- rownames(sim_mat)
    pheatmap::pheatmap(
      sim_mat,
      annotation_row = ann,
      annotation_col = ann,
      show_rownames = isTRUE(show_labels),
      show_colnames = isTRUE(show_labels),
      clustering_method = "average",
      color = grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
      filename = out_path,
      main = title,
      fontsize = 7
    )
  } else {
    grDevices::pdf(out_path, width = 8, height = 8)
    stats::heatmap(
      sim_mat,
      Rowv = NA,
      Colv = NA,
      scale = "none",
      labRow = if (isTRUE(show_labels)) rownames(sim_mat) else NA,
      labCol = if (isTRUE(show_labels)) colnames(sim_mat) else NA,
      main = title,
      col = grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)
    )
    grDevices::dev.off()
  }
  invisible(TRUE)
}

.cluster_similarity_matrix <- function(sim_mat, clusters, cluster_labels = NULL) {
  if (is.null(sim_mat) || !is.matrix(sim_mat)) return(NULL)
  if (is.null(rownames(sim_mat)) || is.null(names(clusters))) return(NULL)
  common <- intersect(rownames(sim_mat), names(clusters))
  if (!length(common)) return(NULL)
  sim_mat <- sim_mat[common, common, drop = FALSE]
  cl <- clusters[common]
  cl_ids <- sort(unique(cl))
  n <- length(cl_ids)
  out <- matrix(NA_real_, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    idx_i <- common[cl == cl_ids[i]]
    for (j in i:n) {
      idx_j <- common[cl == cl_ids[j]]
      val <- mean(sim_mat[idx_i, idx_j, drop = FALSE], na.rm = TRUE)
      out[i, j] <- val
      out[j, i] <- val
    }
  }
  labels <- if (!is.null(cluster_labels)) {
    cluster_labels[as.character(cl_ids)]
  } else {
    paste0("cluster_", cl_ids)
  }
  missing <- is.na(labels) | labels == ""
  if (any(missing)) labels[missing] <- paste0("cluster_", cl_ids[missing])
  rownames(out) <- labels
  colnames(out) <- labels
  out
}

.save_silhouette_plot <- function(dist_obj, clusters, out_path, title) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(invisible(NULL))
  cl <- as.integer(clusters)
  sil <- cluster::silhouette(cl, dist_obj)
  sil_df <- as.data.frame(sil[, c("cluster", "sil_width")])
  names(sil_df) <- c("cluster", "sil_width")
  sil_df$cluster <- as.integer(sil_df$cluster)

  sil_df <- sil_df[order(sil_df$cluster, sil_df$sil_width), , drop = FALSE]
  sil_df$y <- seq_len(nrow(sil_df))
  avg_sil <- mean(sil_df$sil_width, na.rm = TRUE)

  .save_plot_safe(
    out_path,
    function(ggplot2 = TRUE) {
      if (ggplot2) {
        ggplot2::ggplot(sil_df, ggplot2::aes(y = y, x = 0, xend = sil_width, yend = y, color = factor(cluster))) +
          ggplot2::geom_segment(linewidth = 0.6, alpha = 0.85) +
          ggplot2::geom_vline(xintercept = avg_sil, linetype = "dashed", color = "#d95f02") +
          ggplot2::scale_color_viridis_d(option = "C", guide = "none") +
          ggplot2::labs(
            title = title,
            x = "Silhouette coefficient",
            y = "Clustered motifs (sorted within cluster)",
            caption = paste0("Average silhouette = ", sprintf("%.3f", avg_sil))
          ) +
          .theme_pub()
      } else {
        plot(sil_df$sil_width, sil_df$y, type = "h",
             main = title, xlab = "Silhouette coefficient",
             ylab = "Clustered motifs (sorted within cluster)")
        abline(v = avg_sil, col = "red", lty = 2)
      }
    },
    width = 10,
    height = 6
  )
  invisible(TRUE)
}

.save_embedding_plot <- function(dist_obj, clusters, out_path, title, top_clusters = 12L) {
  if (is.null(dist_obj) || is.null(clusters)) return(invisible(NULL))
  cl <- as.integer(clusters)
  coords <- stats::cmdscale(dist_obj, k = 2, eig = FALSE)
  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    cluster = as.factor(cl),
    stringsAsFactors = FALSE
  )
  cl_sizes <- sort(table(cl), decreasing = TRUE)
  top_ids <- names(cl_sizes)[seq_len(min(top_clusters, length(cl_sizes)))]
  df$cluster_top <- ifelse(df$cluster %in% top_ids, as.character(df$cluster), "Other")

  .save_plot_safe(
    out_path,
    function(ggplot2 = TRUE) {
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = cluster_top)) +
          ggplot2::geom_point(size = 0.5, alpha = 0.6) +
          ggplot2::scale_color_viridis_d(option = "C") +
          ggplot2::labs(
            title = title,
            x = "MDS dimension 1",
            y = "MDS dimension 2",
            caption = paste0("Top ", length(top_ids), " clusters colored; others in 'Other'.")
          ) +
          .theme_pub()
      } else {
        plot(df$x, df$y, pch = 16, cex = 0.3,
             main = title, xlab = "MDS dimension 1", ylab = "MDS dimension 2")
      }
    },
    width = 8,
    height = 6
  )
  invisible(TRUE)
}

.k_curve_with_second_derivative <- function(k_df) {
  k_df <- k_df[order(k_df$k), , drop = FALSE]
  if (nrow(k_df) < 3L) {
    return(list(curve = k_df, d2 = data.frame(k = numeric(0), d2 = numeric(0))))
  }
  k <- k_df$k
  s <- k_df$silhouette
  d1 <- diff(s) / diff(k)
  k_mid <- (k[-1] + k[-length(k)]) / 2
  d2 <- diff(d1) / diff(k_mid)
  k_d2 <- k[2:(length(k) - 1L)]
  list(
    curve = k_df,
    d2 = data.frame(k = k_d2, d2 = d2)
  )
}

.select_best_k_from_scores <- function(k_df, d2_df) {
  if (!is.data.frame(k_df) || !nrow(k_df)) return(NA_integer_)
  k_df <- k_df[order(k_df$k), , drop = FALSE]
  rank_sil <- rank(-k_df$silhouette, ties.method = "min", na.last = "keep")
  if (!is.data.frame(d2_df) || !nrow(d2_df)) {
    best <- k_df[order(rank_sil, k_df$k), , drop = FALSE][1, ]
    return(best$k)
  }
  d2_df <- d2_df[order(d2_df$k), , drop = FALSE]
  rank_d2_map <- rank(-d2_df$d2, ties.method = "min", na.last = "keep")
  d2_rank <- rep(max(rank_d2_map, na.rm = TRUE) + 1L, nrow(k_df))
  match_idx <- match(k_df$k, d2_df$k)
  has_d2 <- !is.na(match_idx)
  d2_rank[has_d2] <- rank_d2_map[match_idx[has_d2]]
  combo <- rank_sil + d2_rank
  best <- k_df[order(combo, k_df$k), , drop = FALSE][1, ]
  best$k
}

.plot_best_k_curve <- function(k_df, d2_df, out_path, title, line_color) {
  if (is.null(k_df) || !nrow(k_df)) return(invisible(FALSE))
  k_df <- k_df[order(k_df$k), , drop = FALSE]
  top_n <- min(5L, nrow(k_df))
  top_k <- k_df[order(-k_df$silhouette, k_df$k), ][seq_len(top_n), , drop = FALSE]
  top_k$rank <- seq_len(nrow(top_k))
  d2_df <- if (is.null(d2_df)) data.frame(k = numeric(0), d2 = numeric(0)) else d2_df
  top_d2 <- if (nrow(d2_df)) {
    d2_df[order(-d2_df$d2, d2_df$k), ][seq_len(min(5L, nrow(d2_df))), , drop = FALSE]
  } else {
    NULL
  }
  if (!is.null(top_d2) && nrow(top_d2)) top_d2$rank <- seq_len(nrow(top_d2))

  colors_top <- c("#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e")

  .save_plot_safe(
    out_path,
    function(ggplot2 = TRUE) {
      if (ggplot2 && requireNamespace("gridExtra", quietly = TRUE)) {
        p1 <- ggplot2::ggplot(k_df, ggplot2::aes(x = k, y = silhouette)) +
          ggplot2::geom_line(color = line_color, linewidth = 0.7) +
          ggplot2::geom_point(color = line_color, size = 1.4) +
          ggplot2::geom_point(
            data = top_k,
            ggplot2::aes(x = k, y = silhouette, color = factor(rank)),
            size = 2.6
          ) +
          ggplot2::geom_text(
            data = top_k,
            ggplot2::aes(x = k, y = silhouette, label = k, color = factor(rank)),
            vjust = -0.9,
            size = 3
          ) +
          ggplot2::scale_color_manual(values = colors_top[seq_len(nrow(top_k))]) +
          ggplot2::labs(
            title = "Mean silhouette",
            x = "Number of clusters (K)",
            y = "Mean silhouette"
          ) +
          .theme_pub()

        p2 <- ggplot2::ggplot(d2_df, ggplot2::aes(x = k, y = d2)) +
          ggplot2::geom_line(color = "#6a51a3", linewidth = 0.7) +
          ggplot2::geom_point(color = "#6a51a3", size = 1.4) +
          ggplot2::labs(
            title = "Second derivative",
            x = "Number of clusters (K)",
            y = "Second derivative"
          ) +
          .theme_pub()
        if (!is.null(top_d2) && nrow(top_d2)) {
          p2 <- p2 +
            ggplot2::geom_point(
              data = top_d2,
              ggplot2::aes(x = k, y = d2, color = factor(rank)),
              size = 2.6
            ) +
            ggplot2::geom_text(
              data = top_d2,
              ggplot2::aes(x = k, y = d2, label = k, color = factor(rank)),
              vjust = -0.9,
              size = 3
            ) +
            ggplot2::scale_color_manual(values = colors_top[seq_len(nrow(top_d2))])
        }
        gridExtra::grid.arrange(p1, p2, ncol = 2, top = title)
      } else {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        par(mfrow = c(1, 2))
        plot(k_df$k, k_df$silhouette, type = "b",
             main = "Mean silhouette",
             xlab = "Number of clusters (K)", ylab = "Mean silhouette")
        cols <- colors_top[seq_len(nrow(top_k))]
        points(top_k$k, top_k$silhouette, pch = 19, col = cols)
        text(top_k$k, top_k$silhouette, labels = top_k$k, pos = 3, col = cols)
        if (nrow(d2_df)) {
          plot(d2_df$k, d2_df$d2, type = "b",
               main = "Second derivative",
               xlab = "Number of clusters (K)", ylab = "Second derivative")
          if (!is.null(top_d2) && nrow(top_d2)) {
            cols2 <- colors_top[seq_len(nrow(top_d2))]
            points(top_d2$k, top_d2$d2, pch = 19, col = cols2)
            text(top_d2$k, top_d2$d2, labels = top_d2$k, pos = 3, col = cols2)
          }
        }
        title(main = title, outer = TRUE, line = -2)
      }
    },
    width = 12,
    height = 5
  )
  invisible(TRUE)
}
#' Run footprint-driven motif clustering with quality control (QC) plots
#'
#' @param fp_aligned Output from align_footprints().
#' @param out_dir Output directory to save results and plots.
#' @param base_dir Optional base directory used to build a default out_dir.
#' @param mode Character; "data", "hybrid", or "both".
#' @param ref_similarity Optional reference similarity matrix (e.g., PPM Ncor).
#' @param ref_db Optional database name to compute reference similarity when ref_similarity is NULL.
#' @param cache_dir Optional cache directory (default: out_dir/cache).
#' @param alpha Numeric in [0,1]; hybrid similarity weight for ref_similarity.
#' @param min_peak_support Numeric; minimum weighted peak support per motif.
#' @param sim_method Character; "jaccard" or "cosine".
#' @param linkage Linkage method for hclust().
#' @param target_clusters Target number of clusters for cutree.
#' @param run_mode Character; "full" runs full clustering, "pre" runs minimal K selection only.
#' @param qc_mode Character; "fast" skips slow QC plots (default "fast").
#' @param motif_db Optional motif table (must include motif + gene_symbol); if provided,
#'   adds sub_cluster and sub_cluster_name columns from clustering.
#' @param cluster_source Character; "hybrid" or "data" (default "hybrid").
#' @param cores Integer number of workers (default 1L).
#' @param weight_by_group_size Logical; weight by id_map$group_size if present.
#' @param verbose Logical; emit concise progress messages.
#' @return List with similarity matrices and cluster tables.
#' @export
run_fp_motif_clustering <- function(
  fp_aligned,
  out_dir = NULL,
  base_dir = NULL,
  mode = c("both", "data", "hybrid"),
  ref_similarity = NULL,
  ref_db = NULL,
  cache_dir = NULL,
  alpha = 0.5,
  min_peak_support = 100,
  sim_method = c("jaccard", "cosine"),
  linkage = "average",
  target_clusters = NULL,
  k_grid = NULL,
  run_mode = c("full", "pre"),
  qc_mode = c("none", "fast", "full"),
  motif_db = NULL,
  cluster_source = c("hybrid", "data"),
  cores = NULL,
  save_motif_db = FALSE,
  motif_db_out_dir = NULL,
  motif_db_prefix = "motif_db",
  weight_by_group_size = TRUE,
  verbose = TRUE
) {
  stopifnot(is.list(fp_aligned))
  mode <- match.arg(mode)
  run_mode <- match.arg(run_mode)
  qc_mode <- match.arg(qc_mode)
  do_qc <- !identical(qc_mode, "none")
  cluster_source <- match.arg(cluster_source)
  if (is.null(out_dir)) {
    if (is.null(base_dir) || !nzchar(base_dir)) {
      cli::cli_abort("Provide out_dir or base_dir.")
    }
    out_dir <- file.path(base_dir, "tf_motif_clustering")
  }
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  sim_method <- match.arg(sim_method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cache_dir)) cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cores) || !is.finite(cores) || cores < 1L) {
    cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
    if (!is.finite(cores) || cores < 1L) {
      cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
    }
  }
  cores <- as.integer(cores)
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(cores)
    if (isTRUE(verbose)) {
      .log_inform("data.table threads set to {data.table::getDTthreads()}.")
    }
  }
  if (isTRUE(verbose)) .log_inform("Mode: {mode}.")
  if (isTRUE(verbose)) .log_inform("Run mode: {run_mode}.")

  sim_hybrid <- NULL
  cl_hybrid <- NULL
  target_clusters_h <- NULL
  k_df_h <- NULL

  if (isTRUE(verbose)) {
    .log_inform("Running footprint-driven motif clustering.")
  }

  fp_annotation <- fp_aligned$fp_annotation
  id_map <- fp_aligned$id_map
  if (!is.data.frame(fp_annotation)) stop("fp_aligned$fp_annotation is missing.")
  total_motifs_raw <- length(unique(fp_annotation$motifs))
  total_peaks_raw <- length(unique(fp_annotation$fp_peak))

  tf_map <- NULL
  tf_map_db <- NULL
  if (mode %in% c("hybrid", "both")) {
    if (is.null(ref_similarity) && !is.null(ref_db)) {
      if (isTRUE(verbose)) .log_inform("Reference similarity: start ({ref_db}).")
      ref_similarity <- .get_ref_similarity(ref_db, cache_dir = cache_dir, cores = cores, verbose = verbose)
      if (isTRUE(verbose)) .log_inform("Reference similarity: done.")
    }
  }

  if (!is.null(ref_similarity) && mode %in% c("hybrid", "both")) {
    norm <- .normalize_similarity_ids(ref_similarity)
    ref_similarity <- norm$sim
    tf_map <- norm$tf_map
    norm_ids <- .extract_motif_key(fp_annotation$motifs, db = ref_db)
    if (sum(!is.na(norm_ids)) > 0) {
      changed <- sum(!is.na(norm_ids) & norm_ids != fp_annotation$motifs)
      if (isTRUE(verbose) && changed > 0) {
        .log_inform("Normalized motif IDs to JASPAR matrix IDs for {changed} rows.")
      }
      fp_annotation$motifs <- ifelse(!is.na(norm_ids), norm_ids, fp_annotation$motifs)
    }
  }

  motif_id_tag <- if (!is.null(ref_similarity) && mode %in% c("hybrid", "both")) "norm" else "raw"
  stats_cache_path <- file.path(
    cache_dir,
    sprintf(
      "motif_stats_%s_minp%s_w%s_m%s_p%s.rds",
      motif_id_tag,
      min_peak_support,
      if (isTRUE(weight_by_group_size)) "TRUE" else "FALSE",
      format(length(unique(fp_annotation$motifs)), scientific = FALSE),
      format(length(unique(fp_annotation$fp_peak)), scientific = FALSE)
    )
  )
  stats <- NULL
  if (file.exists(stats_cache_path)) {
    stats <- tryCatch(readRDS(stats_cache_path), error = function(e) NULL)
    if (isTRUE(verbose) && !is.null(stats)) {
      .log_inform("Using cached motif stats: {stats_cache_path}")
    }
  }
  if (is.null(stats) || !is.list(stats) || is.null(stats$motif_counts) || is.null(stats$motifs_keep)) {
    stats <- build_motif_peak_stats(
      fp_annotation = fp_annotation,
      id_map = id_map,
      min_peak_support = min_peak_support,
      weight_by_group_size = weight_by_group_size,
      verbose = verbose
    )
    if (!is.null(stats) && is.list(stats)) {
      saveRDS(stats, stats_cache_path)
      if (isTRUE(verbose)) {
        .log_inform("Saved motif stats cache: {stats_cache_path}")
      }
    }
  }
  motif_ids <- stats$motifs_keep
  if (is.data.frame(motif_db) && "motif" %in% names(motif_db)) {
    gene_col <- NULL
    if ("gene_symbol" %in% names(motif_db)) gene_col <- "gene_symbol"
    if (is.null(gene_col) && "HGNC" %in% names(motif_db)) gene_col <- "HGNC"
    if (!is.null(gene_col)) {
      motif_ids_db <- .extract_motif_key(motif_db$motif, db = ref_db)
      motif_ids_db[is.na(motif_ids_db)] <- motif_db$motif[is.na(motif_ids_db)]
      gene_vals <- motif_db[[gene_col]]
      ok <- !is.na(motif_ids_db) & nzchar(motif_ids_db) & !is.na(gene_vals) & nzchar(gene_vals)
      if (any(ok)) {
        motif_ids_ok <- motif_ids_db[ok]
        gene_ok <- gene_vals[ok]
        keep <- !duplicated(motif_ids_ok)
        tf_map_db <- stats::setNames(gene_ok[keep], motif_ids_ok[keep])
      }
    }
  }
  resolve_tf_names <- function(ids) {
    nm <- rep(NA_character_, length(ids))
    if (!is.null(tf_map_db)) {
      nm <- tf_map_db[ids]
    }
    if (!is.null(tf_map)) {
      nm2 <- tf_map[ids]
      nm[is.na(nm)] <- nm2[is.na(nm)]
    }
    nm[is.na(nm)] <- ids[is.na(nm)]
    nm
  }
  k_df <- NULL
  k_info <- NULL
  k_df_h <- NULL
  k_info_h <- NULL
  target_clusters_h <- NULL

  k_grid_id <- if (is.null(k_grid) || !length(k_grid)) {
    "auto5"
  } else {
    paste0("custom_", length(k_grid), "_", min(k_grid), "_", max(k_grid))
  }
  sim_cache_key <- .make_cache_key("sim_data", sim_method, paste0("minp", min_peak_support),
                                   paste0("w", weight_by_group_size), k_grid_id)
  sim_cache_path <- file.path(cache_dir, paste0(sim_cache_key, ".rds"))
  if (file.exists(sim_cache_path)) {
    if (isTRUE(verbose)) .log_inform("Using cached data similarity: {sim_cache_path}")
    sim_data <- readRDS(sim_cache_path)
  } else {
    if (isTRUE(verbose)) .log_inform("Compute motif similarity ({sim_method}) start.")
    sim_data <- compute_motif_similarity(stats, method = sim_method, cores = cores)
    if (isTRUE(verbose)) .log_inform("Compute motif similarity done.")
    saveRDS(sim_data, sim_cache_path)
  }
  total_weight <- sum(stats$motif_counts$weight, na.rm = TRUE)
  dist_data <- stats::as.dist(1 - sim_data)
  hc_data <- stats::hclust(dist_data, method = linkage)
  target_missing <- is.null(target_clusters)
  if (mode %in% c("data", "both")) {
    k_scores_cache <- file.path(
      cache_dir,
      paste0(.make_cache_key("k_scores_data", sim_cache_key), ".csv")
    )
    if (isTRUE(target_missing)) {
      if (file.exists(k_scores_cache)) {
        k_df <- readr::read_csv(k_scores_cache, show_col_types = FALSE)
        k_info <- .k_curve_with_second_derivative(k_df)
        target_clusters <- .select_best_k_from_scores(k_info$curve, k_info$d2)
      } else {
        k_pick <- .choose_best_k(hc_data, dist_data, k_grid = k_grid, verbose = verbose)
        k_df <- k_pick$scores
        k_info <- list(curve = k_pick$scores, d2 = k_pick$d2)
        target_clusters <- k_pick$k
        if (!is.null(k_df)) {
          readr::write_csv(k_df, k_scores_cache)
        }
      }
      if (!is.null(k_df)) {
        readr::write_csv(k_df, file.path(out_dir, "qc_best_k_scores_data.csv"))
        if (do_qc) {
          .plot_best_k_curve(
            k_df = k_info$curve,
            d2_df = k_info$d2,
            out_path = file.path(out_dir, "qc_best_k_curve_data.pdf"),
            title = "Best K selection curve (data-only)",
            line_color = "#2b8cbe"
          )
        }
      }
      if (isTRUE(verbose)) {
        .log_inform("Best K (data-only) selected: {target_clusters}.")
      }
    }
    cl_data <- if (run_mode == "full") stats::cutree(hc_data, k = target_clusters) else NULL
  } else {
    cl_data <- NULL
  }

  if (mode %in% c("data", "both") && run_mode == "full") {
    motif_ids <- names(cl_data)
    cluster_ids <- split(motif_ids, cl_data)
    cl_map <- data.frame(motif = motif_ids, cluster_id = as.integer(cl_data), stringsAsFactors = FALSE)
    peak_summary <- .cluster_peak_summary(fp_annotation, cl_map, top_n = 5L)
    cluster_df <- data.frame(
      cluster = sprintf("cluster_%03d", seq_along(cluster_ids)),
      cluster_id = as.integer(names(cluster_ids)),
      id = vapply(cluster_ids, function(x) paste(x, collapse = ","), character(1)),
      tf_name = vapply(cluster_ids, function(x) {
        nm <- resolve_tf_names(x)
        paste(unique(nm), collapse = ",")
      }, character(1)),
      n_motifs = vapply(cluster_ids, length, integer(1)),
      pct_motifs = vapply(cluster_ids, length, integer(1)) / length(motif_ids),
      peak_support = vapply(cluster_ids, function(x) {
        sum(stats$motif_counts$weight[match(x, stats$motif_counts$motif)], na.rm = TRUE)
      }, numeric(1)),
      pct_peak_support = vapply(cluster_ids, function(x) {
        sum(stats$motif_counts$weight[match(x, stats$motif_counts$motif)], na.rm = TRUE) / total_weight
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
    cluster_df$pct_motifs_of_kept <- cluster_df$pct_motifs
    cluster_df$pct_peak_support_of_total <- cluster_df$pct_peak_support
    if (!is.null(peak_summary)) {
      cluster_df <- merge(cluster_df, peak_summary, by = "cluster_id", all.x = TRUE)
    }
    cluster_df$cluster_name <- vapply(
      seq_len(nrow(cluster_df)),
      function(i) {
        .make_cluster_name(
          motif_ids = strsplit(cluster_df$id[i], ",", fixed = TRUE)[[1]],
          tf_names = strsplit(cluster_df$tf_name[i], ",", fixed = TRUE)[[1]],
          mode = "rich"
        )
      },
      character(1)
    )
    if (anyDuplicated(cluster_df$cluster_name)) {
      dup_ids <- which(duplicated(cluster_df$cluster_name) | duplicated(cluster_df$cluster_name, fromLast = TRUE))
      cluster_df$cluster_name[dup_ids] <- paste0(
        cluster_df$cluster_name[dup_ids],
        "_",
        cluster_df$cluster_id[dup_ids]
      )
    }
    cluster_df <- cluster_df[order(-cluster_df$n_motifs, -cluster_df$peak_support), , drop = FALSE]
  } else {
    cluster_df <- NULL
    peak_summary <- NULL
  }

  if (run_mode == "full") {
    motif_counts_out <- stats$motif_counts
    motif_counts_out$tf_name <- resolve_tf_names(motif_counts_out$motif)
    readr::write_csv(motif_counts_out, file.path(out_dir, "motif_peak_support.csv"))
    if (!is.null(peak_summary)) {
      readr::write_csv(peak_summary, file.path(out_dir, "cluster_peak_summary_data.csv"))
    }
    if (do_qc && mode %in% c("data", "both")) {
      if (isTRUE(verbose)) .log_inform("QC plots (data-only): start.")
      k_tag <- if (!is.null(target_clusters)) paste0("_K", target_clusters) else ""
      heatmap_data_path <- file.path(out_dir, sprintf("qc_similarity_heatmap_data%s.pdf", k_tag))
      if (!file.exists(heatmap_data_path)) {
        .save_similarity_heatmap(
          sim_mat = sim_data,
          clusters = cl_data,
          out_path = heatmap_data_path,
          title = "Motif similarity heatmap (footprint co-occurrence)",
          show_labels = FALSE
        )
      } else if (isTRUE(verbose)) {
        .log_inform("QC plots (data-only): skipping existing {heatmap_data_path}.")
      }
      cluster_labels <- NULL
      if (!is.null(cluster_df) && all(c("cluster_id", "cluster_name") %in% names(cluster_df))) {
        cluster_labels <- stats::setNames(cluster_df$cluster_name, cluster_df$cluster_id)
      }
      cluster_sim <- .cluster_similarity_matrix(sim_data, cl_data, cluster_labels = cluster_labels)
      if (!is.null(cluster_sim)) {
        cluster_annot <- stats::setNames(rep("cluster", nrow(cluster_sim)), rownames(cluster_sim))
        heatmap_cluster_path <- file.path(out_dir, sprintf("qc_similarity_heatmap_clusters_data%s.pdf", k_tag))
        if (!file.exists(heatmap_cluster_path)) {
          .save_similarity_heatmap(
            sim_mat = cluster_sim,
            clusters = cluster_annot,
            out_path = heatmap_cluster_path,
            title = "Cluster similarity heatmap (footprint co-occurrence)",
            show_labels = FALSE
          )
        } else if (isTRUE(verbose)) {
          .log_inform("QC plots (data-only): skipping existing {heatmap_cluster_path}.")
        }
      }
      .save_silhouette_plot(
        dist_obj = dist_data,
        clusters = cl_data,
        out_path = file.path(out_dir, "qc_silhouette_profile_data.pdf"),
        title = paste0("Silhouette profile (data-only, K=", target_clusters, ")")
      )
      .save_embedding_plot(
        dist_obj = dist_data,
        clusters = cl_data,
        out_path = file.path(out_dir, "qc_mds_scatter_data.pdf"),
        title = paste0("MDS scatter (data-only, K=", target_clusters, ")")
      )
      if (isTRUE(verbose)) .log_inform("QC plots (data-only): done.")
    }
    readme_path <- file.path(out_dir, "README_motif_clustering.txt")
    readme_lines <- c(
    "Motif clustering outputs (footprint-driven + optional hybrid)",
    "",
    "Overview",
    "- mode=data: clustering from footprint co-occurrence similarity only.",
    "- mode=hybrid: clustering from hybrid similarity only.",
    "- mode=both: compute and report both data and hybrid outputs.",
    "- hybrid: similarity = alpha * reference PPM similarity + (1-alpha) * footprint similarity.",
    "- reference similarity is cached under: <out_dir>/cache/ref_similarity_<db>.rds when ref_db is used.",
    "- cluster_name is built from TF-family heuristics (no external cluster labels).",
    "- cluster_name supports compact/rich modes; rich mode is default.",
    "",
    "Main tables",
    "1) clusters_data_only.csv / clusters_hybrid.csv",
    "   - cluster: formatted cluster label (cluster_###).",
    "   - cluster_id: numeric cluster id from cutree().",
    "   - cluster_name: compact, human-readable label (unique; derived from TF names).",
    "   - id: motif IDs in the cluster (JASPAR matrix IDs).",
    "   - tf_name: TF names mapped from motif_db (or reference similarity); falls back to motif ID when missing.",
    "   - n_motifs: number of motifs in the cluster.",
    "   - pct_motifs / pct_motifs_of_kept: fraction of kept motifs in this cluster.",
    "   - peak_support: total weighted aligned-peak support across motifs in the cluster.",
    "   - pct_peak_support / pct_peak_support_of_total: fraction of total peak_support contributed by this cluster.",
    "   - n_peaks: number of aligned peaks where any motif in the cluster appears.",
    "   - shared_peaks: number of aligned peaks where >=2 motifs from this cluster co-occur.",
    "   - pct_shared_peaks: shared_peaks / n_peaks.",
    "   - avg_motifs_per_peak: average number of cluster motifs per aligned peak.",
    "   - top_shared_peaks: up to 5 example peaks with highest within-cluster motif overlap.",
    "",
    "2) motif_peak_support.csv",
    "   - motif: motif ID (JASPAR matrix ID).",
    "   - weight: weighted count of aligned peaks containing the motif (uses group_size if enabled).",
    "   - tf_name: mapped TF name when available.",
    "",
    "3) cluster_peak_summary_data.csv / cluster_peak_summary_hybrid.csv",
    "   - cluster_id: numeric cluster id.",
    "   - n_peaks: number of aligned peaks with any motif from this cluster.",
    "   - shared_peaks: number of aligned peaks with >=2 motifs from this cluster.",
    "   - pct_shared_peaks: shared_peaks / n_peaks.",
    "   - avg_motifs_per_peak: mean number of cluster motifs per aligned peak.",
    "   - top_shared_peaks: up to 5 example peaks with highest within-cluster overlap.",
    "",
    "4) clustering_summary_data.csv / clustering_summary_hybrid.csv",
    "   - total_motifs_raw: total motifs before filtering.",
    "   - motifs_kept: motifs retained after min_peak_support filter.",
    "   - pct_motifs_kept: motifs_kept / total_motifs_raw.",
    "   - total_aligned_peaks: number of aligned peaks in fp_annotation.",
    "   - total_peak_support_weight: sum of motif peak support across kept motifs.",
    "   - min_peak_support: filter threshold.",
    "   - sim_method: footprint similarity method (jaccard or cosine).",
    "   - weight_by_group_size: whether aligned-peak weights were applied.",
    "   - target_clusters_data / target_clusters_hybrid: chosen K for clustering.",
    "   - selected_by_best_k: TRUE if K was chosen automatically by silhouette.",
    "   - alpha (hybrid only): blend weight for reference similarity.",
    "",
    "Similarity matrices",
    "- similarity_data_only.rds: footprint similarity matrix (motif x motif).",
    "- similarity_hybrid.rds: hybrid similarity matrix (if reference provided).",
    "",
    "Figures (PDF)",
    "- qc_motifs_per_peak.pdf: distribution of how many motifs annotate each aligned peak.",
    "- qc_motif_support.pdf: histogram of motif peak support (weighted).",
    "- qc_motif_support_cdf.pdf: cumulative distribution of motif peak support.",
    "- qc_similarity_distribution.pdf: distribution of footprint-based similarity scores.",
    "- qc_cluster_sizes_data.pdf: cluster size distribution (data-only).",
    "- qc_cluster_sizes_top20_data.pdf: largest 20 clusters (data-only).",
    "- qc_best_k_curve_data.pdf: silhouette vs K and second derivative; top 5 K highlighted (data-only).",
    "- qc_similarity_ref_vs_data.pdf: PPM reference similarity vs footprint similarity (with Spearman/Pearson/Kendall; qc_mode=full only).",
    "- qc_similarity_heatmap_data_K<k>.pdf: heatmap of footprint similarity; clustered by motif similarity.",
    "- qc_similarity_heatmap_clusters_data_K<k>.pdf: heatmap of cluster-by-cluster footprint similarity.",
    "- qc_silhouette_profile_data.pdf: silhouette profile for the selected K (data-only).",
    "- qc_mds_scatter_data.pdf: 2D MDS scatter of motifs colored by cluster (data-only).",
    "- qc_cluster_sizes_hybrid.pdf: cluster size distribution (hybrid).",
    "- qc_cluster_sizes_top20_hybrid.pdf: largest 20 clusters (hybrid).",
    "- qc_best_k_curve_hybrid.pdf: silhouette vs K and second derivative; top 5 K highlighted (hybrid).",
    "- qc_similarity_heatmap_hybrid_K<k>.pdf: heatmap of hybrid similarity; clustered by motif similarity.",
    "- qc_similarity_heatmap_clusters_hybrid_K<k>.pdf: heatmap of cluster-by-cluster hybrid similarity.",
    "- qc_silhouette_profile_hybrid.pdf: silhouette profile for the selected K (hybrid).",
    "- qc_mds_scatter_hybrid.pdf: 2D MDS scatter of motifs colored by cluster (hybrid).",
    "",
    "Notes",
    "- Peak support uses aligned peak groups; higher values indicate motifs that recur across more aligned peaks.",
    "- Shared peaks emphasize co-occurrence (two or more motifs in the same aligned peak).",
    "- Hybrid clusters reflect both sequence similarity and data-driven co-occurrence."
  )
    if (isTRUE(verbose)) .log_inform("README: writing.")
    writeLines(readme_lines, readme_path)
    if (isTRUE(verbose)) .log_inform("README: done.")
    if (!is.null(cluster_df)) {
      readr::write_csv(cluster_df, file.path(out_dir, "clusters_data_only.csv"))
    }
    saveRDS(sim_data, file.path(out_dir, "similarity_data_only.rds"))
    if (mode %in% c("data", "both")) {
      summary_df <- data.frame(
        total_motifs_raw = total_motifs_raw,
        motifs_kept = length(motif_ids),
        pct_motifs_kept = length(motif_ids) / max(1L, total_motifs_raw),
        total_aligned_peaks = total_peaks_raw,
        total_peak_support_weight = total_weight,
        min_peak_support = min_peak_support,
        sim_method = sim_method,
        weight_by_group_size = isTRUE(weight_by_group_size),
        target_clusters_data = target_clusters,
        selected_by_best_k = isTRUE(target_missing),
        stringsAsFactors = FALSE
      )
      readr::write_csv(summary_df, file.path(out_dir, "clustering_summary_data.csv"))
    }
  }

  if (!is.null(ref_similarity) && mode %in% c("hybrid", "both")) {
    ref_ids <- rownames(ref_similarity)
    if (is.null(ref_ids)) {
      cli::cli_warn("ref_similarity has no rownames; skipping hybrid similarity.")
      ref_similarity <- NULL
    } else {
      missing_in_ref <- setdiff(motif_ids, ref_ids)
      if (length(missing_in_ref)) {
        cli::cli_warn("ref_similarity missing {length(missing_in_ref)} motif(s); hybrid uses intersection only.")
      }
      keep_ids <- intersect(motif_ids, ref_ids)
      if (length(keep_ids) < 2L) {
        cli::cli_warn("Too few overlapping motifs for hybrid similarity; skipping.")
        ref_similarity <- NULL
      } else {
        ref_similarity <- ref_similarity[keep_ids, keep_ids, drop = FALSE]
        sim_data <- sim_data[keep_ids, keep_ids, drop = FALSE]
        motif_ids <- keep_ids
      }
    }
  }

  if (!is.null(ref_similarity)) {
    alpha <- min(1, max(0, alpha))
    sim_hybrid_cache_key <- .make_cache_key(
      "sim_hybrid",
      sim_method,
      paste0("minp", min_peak_support),
      paste0("w", weight_by_group_size),
      paste0("alpha", alpha),
      paste0("ref", ref_db),
      k_grid_id
    )
    sim_hybrid_cache_path <- file.path(cache_dir, paste0(sim_hybrid_cache_key, ".rds"))
    if (file.exists(sim_hybrid_cache_path)) {
      if (isTRUE(verbose)) .log_inform("Using cached hybrid similarity: {sim_hybrid_cache_path}")
      sim_hybrid <- readRDS(sim_hybrid_cache_path)
      has_ids <- !is.null(rownames(sim_hybrid)) && all(motif_ids %in% rownames(sim_hybrid))
      if (!has_ids) {
        sim_hybrid <- NULL
      } else {
        sim_hybrid <- sim_hybrid[motif_ids, motif_ids, drop = FALSE]
      }
    }
    if (is.null(sim_hybrid)) {
      sim_hybrid <- alpha * ref_similarity + (1 - alpha) * sim_data
      saveRDS(sim_hybrid, sim_hybrid_cache_path)
    }

    dist_hybrid <- stats::as.dist(1 - sim_hybrid)
    hc_hybrid <- stats::hclust(dist_hybrid, method = linkage)
    k_scores_cache_h <- file.path(
      cache_dir,
      paste0(.make_cache_key("k_scores_hybrid", sim_hybrid_cache_key), ".csv")
    )
    if (isTRUE(target_missing)) {
      if (file.exists(k_scores_cache_h)) {
        k_df_h <- readr::read_csv(k_scores_cache_h, show_col_types = FALSE)
        k_info_h <- .k_curve_with_second_derivative(k_df_h)
        target_clusters_h <- .select_best_k_from_scores(k_info_h$curve, k_info_h$d2)
      } else {
        k_pick_h <- .choose_best_k(hc_hybrid, dist_hybrid, k_grid = k_grid, verbose = verbose)
        k_df_h <- k_pick_h$scores
        k_info_h <- list(curve = k_pick_h$scores, d2 = k_pick_h$d2)
        target_clusters_h <- k_pick_h$k
        if (!is.null(k_df_h)) {
          readr::write_csv(k_df_h, k_scores_cache_h)
        }
      }
      if (!is.null(k_df_h)) {
        readr::write_csv(k_df_h, file.path(out_dir, "qc_best_k_scores_hybrid.csv"))
        if (do_qc) {
          .plot_best_k_curve(
            k_df = k_info_h$curve,
            d2_df = k_info_h$d2,
            out_path = file.path(out_dir, "qc_best_k_curve_hybrid.pdf"),
            title = "Best K selection curve (hybrid)",
            line_color = "#e6550d"
          )
        }
      }
      if (isTRUE(verbose)) {
        .log_inform("Best K (hybrid) selected: {target_clusters_h}.")
      }
    } else {
      target_clusters_h <- target_clusters
    }
    cl_hybrid <- if (run_mode == "full") stats::cutree(hc_hybrid, k = target_clusters_h) else NULL

    if (run_mode == "full") {
      cluster_ids_h <- split(names(cl_hybrid), cl_hybrid)
      cl_map_h <- data.frame(motif = names(cl_hybrid), cluster_id = as.integer(cl_hybrid), stringsAsFactors = FALSE)
      peak_summary_h <- .cluster_peak_summary(fp_annotation, cl_map_h, top_n = 5L)
      cluster_df_h <- data.frame(
      cluster = sprintf("cluster_%03d", seq_along(cluster_ids_h)),
      cluster_id = as.integer(names(cluster_ids_h)),
      id = vapply(cluster_ids_h, function(x) paste(x, collapse = ","), character(1)),
      tf_name = vapply(cluster_ids_h, function(x) {
          nm <- resolve_tf_names(x)
          paste(unique(nm), collapse = ",")
        }, character(1)),
        n_motifs = vapply(cluster_ids_h, length, integer(1)),
        pct_motifs = vapply(cluster_ids_h, length, integer(1)) / length(names(cl_hybrid)),
        peak_support = vapply(cluster_ids_h, function(x) {
          sum(stats$motif_counts$weight[match(x, stats$motif_counts$motif)], na.rm = TRUE)
        }, numeric(1)),
        pct_peak_support = vapply(cluster_ids_h, function(x) {
          sum(stats$motif_counts$weight[match(x, stats$motif_counts$motif)], na.rm = TRUE) / total_weight
        }, numeric(1)),
        stringsAsFactors = FALSE
      )
      cluster_df_h$pct_motifs_of_kept <- cluster_df_h$pct_motifs
      cluster_df_h$pct_peak_support_of_total <- cluster_df_h$pct_peak_support
      if (!is.null(peak_summary_h)) {
        cluster_df_h <- merge(cluster_df_h, peak_summary_h, by = "cluster_id", all.x = TRUE)
      }
      cluster_df_h$cluster_name <- vapply(
        seq_len(nrow(cluster_df_h)),
        function(i) {
          .make_cluster_name(
            motif_ids = strsplit(cluster_df_h$id[i], ",", fixed = TRUE)[[1]],
            tf_names = strsplit(cluster_df_h$tf_name[i], ",", fixed = TRUE)[[1]],
            mode = "rich"
          )
        },
        character(1)
      )
      if (anyDuplicated(cluster_df_h$cluster_name)) {
        dup_ids <- which(duplicated(cluster_df_h$cluster_name) | duplicated(cluster_df_h$cluster_name, fromLast = TRUE))
        cluster_df_h$cluster_name[dup_ids] <- paste0(
          cluster_df_h$cluster_name[dup_ids],
          "_",
          cluster_df_h$cluster_id[dup_ids]
        )
      }
      cluster_df_h <- cluster_df_h[order(-cluster_df_h$n_motifs, -cluster_df_h$peak_support), , drop = FALSE]

      if (isTRUE(verbose)) .log_inform("Hybrid clustering: start.")
      readr::write_csv(cluster_df_h, file.path(out_dir, "clusters_hybrid.csv"))
      if (!is.null(peak_summary_h)) {
        readr::write_csv(peak_summary_h, file.path(out_dir, "cluster_peak_summary_hybrid.csv"))
      }
      if (do_qc) {
        if (isTRUE(verbose)) .log_inform("QC plots (hybrid): start.")
        k_tag <- if (!is.null(target_clusters_h)) paste0("_K", target_clusters_h) else ""
        heatmap_hybrid_path <- file.path(out_dir, sprintf("qc_similarity_heatmap_hybrid%s.pdf", k_tag))
        if (!file.exists(heatmap_hybrid_path)) {
          .save_similarity_heatmap(
            sim_mat = sim_hybrid,
            clusters = cl_hybrid,
            out_path = heatmap_hybrid_path,
            title = "Motif similarity heatmap (hybrid reference + footprint)",
            show_labels = FALSE
          )
        } else if (isTRUE(verbose)) {
          .log_inform("QC plots (hybrid): skipping existing {heatmap_hybrid_path}.")
        }
        cluster_labels_h <- NULL
        if (!is.null(cluster_df_h) && all(c("cluster_id", "cluster_name") %in% names(cluster_df_h))) {
          cluster_labels_h <- stats::setNames(cluster_df_h$cluster_name, cluster_df_h$cluster_id)
        }
        cluster_sim_h <- .cluster_similarity_matrix(sim_hybrid, cl_hybrid, cluster_labels = cluster_labels_h)
        if (!is.null(cluster_sim_h)) {
          cluster_annot_h <- stats::setNames(rep("cluster", nrow(cluster_sim_h)), rownames(cluster_sim_h))
          heatmap_cluster_h_path <- file.path(out_dir, sprintf("qc_similarity_heatmap_clusters_hybrid%s.pdf", k_tag))
          if (!file.exists(heatmap_cluster_h_path)) {
            .save_similarity_heatmap(
              sim_mat = cluster_sim_h,
              clusters = cluster_annot_h,
              out_path = heatmap_cluster_h_path,
              title = "Cluster similarity heatmap (hybrid reference + footprint)",
              show_labels = FALSE
            )
          } else if (isTRUE(verbose)) {
            .log_inform("QC plots (hybrid): skipping existing {heatmap_cluster_h_path}.")
          }
        }
        .save_silhouette_plot(
          dist_obj = dist_hybrid,
          clusters = cl_hybrid,
          out_path = file.path(out_dir, "qc_silhouette_profile_hybrid.pdf"),
          title = paste0("Silhouette profile (hybrid, K=", target_clusters_h, ")")
        )
        .save_embedding_plot(
          dist_obj = dist_hybrid,
          clusters = cl_hybrid,
          out_path = file.path(out_dir, "qc_mds_scatter_hybrid.pdf"),
          title = paste0("MDS scatter (hybrid, K=", target_clusters_h, ")")
        )
        if (isTRUE(verbose)) .log_inform("QC plots (hybrid): done.")
      }
      if (isTRUE(verbose)) .log_inform("Hybrid clustering: done.")
      saveRDS(sim_hybrid, file.path(out_dir, "similarity_hybrid.rds"))
      summary_h <- data.frame(
        total_motifs_raw = total_motifs_raw,
        motifs_kept = length(names(cl_hybrid)),
        pct_motifs_kept = length(names(cl_hybrid)) / max(1L, total_motifs_raw),
        total_aligned_peaks = total_peaks_raw,
        total_peak_support_weight = total_weight,
        min_peak_support = min_peak_support,
        sim_method = sim_method,
        weight_by_group_size = isTRUE(weight_by_group_size),
        target_clusters_hybrid = target_clusters_h,
        selected_by_best_k = isTRUE(target_missing),
        alpha = alpha,
        stringsAsFactors = FALSE
      )
      readr::write_csv(summary_h, file.path(out_dir, "clustering_summary_hybrid.csv"))
    } else {
      cluster_df_h <- NULL
    }
  } else {
    sim_hybrid <- NULL
    cl_hybrid <- NULL
  }

  if (run_mode == "pre") {
    if (isTRUE(verbose)) {
      .log_inform("Pre-run mode: skipping full clustering outputs.")
    }
    return(list(
      target_clusters_data = if (mode %in% c("data", "both")) target_clusters else NULL,
      target_clusters_hybrid = if (!is.null(target_clusters_h)) target_clusters_h else NULL,
      k_scores_data = k_df,
      k_scores_hybrid = k_df_h,
      similarity_data = sim_data,
      similarity_hybrid = sim_hybrid
    ))
  }

  if (!do_qc) {
    if (isTRUE(verbose)) .log_inform("QC plots disabled (qc_mode = 'none').")
  }

  # QC plots
  if (do_qc && isTRUE(verbose)) .log_inform("QC plots (global distributions): start.")
  motifs_per_peak <- if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(fp_annotation[, c("fp_peak", "motifs")])
    dt <- unique(dt)
    dt[, .N, by = fp_peak]$N
  } else {
    tapply(fp_annotation$motifs, fp_annotation$fp_peak, function(x) length(unique(x)))
  }

  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_motifs_per_peak.pdf"),
    function(ggplot2 = TRUE) {
      if (ggplot2) {
        df <- data.frame(n_motifs = motifs_per_peak)
        ggplot2::ggplot(df, ggplot2::aes(x = n_motifs)) +
          ggplot2::geom_histogram(bins = 60, fill = "#2b8cbe", color = "white", linewidth = 0.2) +
          ggplot2::labs(
            title = "Motifs per aligned peak (footprint annotation)",
            x = "Motif count per aligned peak",
            y = "Aligned peak count",
            caption = paste0("Motifs de-duplicated within peak; weighted=", isTRUE(weight_by_group_size))
          ) +
          .theme_pub()
      } else {
        hist(motifs_per_peak, breaks = 60,
             main = "Motifs per aligned peak (footprint annotation)",
             xlab = "Motif count per aligned peak", ylab = "Aligned peak count")
      }
    }
  )

  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_motif_support.pdf"),
    function(ggplot2 = TRUE) {
      df <- stats$motif_counts
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = weight)) +
          ggplot2::geom_histogram(bins = 60, fill = "#1b9e77", color = "white", linewidth = 0.2) +
          ggplot2::labs(
            title = "Motif peak support (weighted aligned peaks)",
            x = "Weighted aligned peak count per motif",
            y = "Motif count",
            caption = paste0("Minimum support filter: ", min_peak_support)
          ) +
          .theme_pub()
      } else {
        hist(df$weight, breaks = 60,
             main = "Motif peak support (weighted aligned peaks)",
             xlab = "Weighted aligned peak count per motif", ylab = "Motif count")
      }
    }
  )

  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_motif_support_cdf.pdf"),
    function(ggplot2 = TRUE) {
      df <- stats$motif_counts
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = weight)) +
          ggplot2::stat_ecdf(geom = "step", color = "#1b9e77", linewidth = 0.8) +
          ggplot2::labs(
            title = "Motif peak support cumulative distribution",
            x = "Weighted aligned peak count per motif",
            y = "Cumulative fraction of motifs",
            caption = paste0("Minimum support filter: ", min_peak_support)
          ) +
          .theme_pub()
      } else {
        plot(stats::ecdf(df$weight), main = "Motif peak support cumulative distribution",
             xlab = "Weighted aligned peak count per motif",
             ylab = "Cumulative fraction of motifs")
      }
    }
  )

  sim_vals <- sim_data[upper.tri(sim_data)]
  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_similarity_distribution.pdf"),
    function(ggplot2 = TRUE) {
      df <- data.frame(sim = sim_vals)
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = sim)) +
          ggplot2::geom_histogram(bins = 60, fill = "#e08214", color = "white", linewidth = 0.2) +
          ggplot2::labs(
            title = "Motif similarity distribution (footprint co-occurrence)",
            x = "Similarity score",
            y = "Motif pair count",
            caption = paste0("Similarity method: ", sim_method, "; min_peak_support=", min_peak_support)
          ) +
          .theme_pub()
      } else {
        hist(sim_vals, breaks = 60,
             main = "Motif similarity distribution (footprint co-occurrence)",
             xlab = "Similarity score", ylab = "Motif pair count")
      }
    }
  )

  cluster_sizes <- sort(table(cl_data), decreasing = TRUE)
  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_cluster_sizes_data.pdf"),
    function(ggplot2 = TRUE) {
      df <- data.frame(size = as.integer(cluster_sizes))
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = size)) +
          ggplot2::geom_histogram(bins = 60, fill = "#756bb1", color = "white", linewidth = 0.2) +
          ggplot2::labs(
            title = "Cluster size distribution (footprint similarity only)",
            x = "Motifs per cluster",
            y = "Cluster count",
            caption = paste0("Target clusters: ", target_clusters)
          ) +
          .theme_pub()
      } else {
        hist(df$size, breaks = 60,
             main = "Cluster size distribution (footprint similarity only)",
             xlab = "Motifs per cluster", ylab = "Cluster count")
      }
    }
  )

  if (do_qc) .save_plot_safe(
    file.path(out_dir, "qc_cluster_sizes_top20_data.pdf"),
    function(ggplot2 = TRUE) {
      top_sizes <- sort(as.integer(cluster_sizes), decreasing = TRUE)[seq_len(min(20L, length(cluster_sizes)))]
      df <- data.frame(rank = seq_along(top_sizes), size = top_sizes)
      if (ggplot2) {
        ggplot2::ggplot(df, ggplot2::aes(x = factor(rank), y = size)) +
          ggplot2::geom_col(fill = "#756bb1") +
          ggplot2::labs(
            title = "Top 20 cluster sizes (footprint similarity only)",
            x = "Cluster rank (largest to smallest)",
            y = "Motifs per cluster",
            caption = paste0("Target clusters: ", target_clusters)
          ) +
          .theme_pub()
      } else {
        barplot(top_sizes, main = "Top 20 cluster sizes (footprint similarity only)",
                xlab = "Cluster rank (largest to smallest)",
                ylab = "Motifs per cluster")
      }
    }
  )

  if (do_qc) {
    if (file.exists(file.path(out_dir, "qc_best_k_scores_data.csv"))) {
      k_df <- readr::read_csv(file.path(out_dir, "qc_best_k_scores_data.csv"), show_col_types = FALSE)
      k_info <- .k_curve_with_second_derivative(k_df)
      .plot_best_k_curve(
        k_df = k_info$curve,
        d2_df = k_info$d2,
        out_path = file.path(out_dir, "qc_best_k_curve_data.pdf"),
        title = "Best K selection curve (data-only)",
        line_color = "#2b8cbe"
      )
    }

    if (!is.null(ref_similarity) && !is.null(sim_hybrid)) {
      if (qc_mode == "full") {
        ref_vals <- ref_similarity[upper.tri(ref_similarity)]
        .save_plot_safe(
          file.path(out_dir, "qc_similarity_ref_vs_data.pdf"),
          function(ggplot2 = TRUE) {
            df <- data.frame(ref = ref_vals, data = sim_vals)
            corr_s <- stats::cor(df$ref, df$data, method = "spearman", use = "complete.obs")
            corr_p <- stats::cor(df$ref, df$data, method = "pearson", use = "complete.obs")
            corr_k <- stats::cor(df$ref, df$data, method = "kendall", use = "complete.obs")
            if (ggplot2) {
              ggplot2::ggplot(df, ggplot2::aes(x = ref, y = data)) +
                ggplot2::geom_point(alpha = 0.06, size = 0.2, color = "#3c3c3c") +
                ggplot2::labs(
                  title = "Similarity: PPM reference vs footprint co-occurrence",
                  x = "Reference similarity (PPM length-normalized correlation)",
                  y = "Footprint similarity",
                  caption = paste0(
                    "Hybrid alpha=", alpha, "; similarity method=", sim_method,
                    "; Spearman r=", sprintf("%.3f", corr_s),
                    "; Pearson r=", sprintf("%.3f", corr_p),
                    "; Kendall tau=", sprintf("%.3f", corr_k)
                  )
                ) +
                .theme_pub()
            } else {
              plot(ref_vals, sim_vals, pch = 16, cex = 0.3,
                   main = "Similarity: PPM reference vs footprint co-occurrence",
                   xlab = "Reference similarity (PPM length-normalized correlation)",
                   ylab = "Footprint similarity")
              legend("topleft", legend = c(
                paste0("Spearman r=", sprintf("%.3f", corr_s)),
                paste0("Pearson r=", sprintf("%.3f", corr_p)),
                paste0("Kendall tau=", sprintf("%.3f", corr_k))
              ), bty = "n")
            }
          }
        )
      } else if (isTRUE(verbose)) {
        .log_inform("QC plots (hybrid): skipping qc_similarity_ref_vs_data.pdf (qc_mode=fast).")
      }

      cluster_sizes_h <- sort(table(cl_hybrid), decreasing = TRUE)
      .save_plot_safe(
        file.path(out_dir, "qc_cluster_sizes_hybrid.pdf"),
        function(ggplot2 = TRUE) {
          df <- data.frame(size = as.integer(cluster_sizes_h))
          if (ggplot2) {
            ggplot2::ggplot(df, ggplot2::aes(x = size)) +
              ggplot2::geom_histogram(bins = 60, fill = "#e6550d", color = "white", linewidth = 0.2) +
              ggplot2::labs(
                title = "Cluster size distribution (hybrid reference + footprint)",
                x = "Motifs per cluster",
                y = "Cluster count",
                caption = paste0("Hybrid alpha=", alpha, "; target clusters=", target_clusters)
              ) +
              .theme_pub()
          } else {
            hist(df$size, breaks = 60,
                 main = "Cluster size distribution (hybrid reference + footprint)",
                 xlab = "Motifs per cluster", ylab = "Cluster count")
          }
        }
      )

      .save_plot_safe(
        file.path(out_dir, "qc_cluster_sizes_top20_hybrid.pdf"),
        function(ggplot2 = TRUE) {
          top_sizes <- sort(as.integer(cluster_sizes_h), decreasing = TRUE)[seq_len(min(20L, length(cluster_sizes_h)))]
          df <- data.frame(rank = seq_along(top_sizes), size = top_sizes)
          if (ggplot2) {
            ggplot2::ggplot(df, ggplot2::aes(x = factor(rank), y = size)) +
              ggplot2::geom_col(fill = "#e6550d") +
              ggplot2::labs(
                title = "Top 20 cluster sizes (hybrid reference + footprint)",
                x = "Cluster rank (largest to smallest)",
                y = "Motifs per cluster",
                caption = paste0("Hybrid alpha=", alpha)
              ) +
              .theme_pub()
          } else {
            barplot(top_sizes, main = "Top 20 cluster sizes (hybrid reference + footprint)",
                    xlab = "Cluster rank (largest to smallest)",
                    ylab = "Motifs per cluster")
          }
        }
      )

      if (file.exists(file.path(out_dir, "qc_best_k_scores_hybrid.csv"))) {
        k_df <- readr::read_csv(file.path(out_dir, "qc_best_k_scores_hybrid.csv"), show_col_types = FALSE)
        k_info <- .k_curve_with_second_derivative(k_df)
        .plot_best_k_curve(
          k_df = k_info$curve,
          d2_df = k_info$d2,
          out_path = file.path(out_dir, "qc_best_k_curve_hybrid.pdf"),
          title = "Best K selection curve (hybrid)",
          line_color = "#e6550d"
        )
      }
    }
  }

  cluster_df_use <- if (cluster_source == "hybrid") {
    get0("cluster_df_h", inherits = FALSE)
  } else {
    cluster_df
  }
  if (cluster_source == "hybrid" && is.null(cluster_df_use)) {
    cluster_df_use <- cluster_df
  }
  motif_db_updated <- .apply_cluster_to_motif_db(motif_db, cluster_df_use, db = ref_db)

  if (isTRUE(save_motif_db) && !is.null(motif_db_updated)) {
    if (is.null(motif_db_out_dir)) motif_db_out_dir <- out_dir
    if (!dir.exists(motif_db_out_dir)) {
      dir.create(motif_db_out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    db_tag <- if (!is.null(ref_db)) ref_db else "motifdb"
    k_tag <- if (!is.null(target_clusters)) target_clusters else "auto"
    out_path <- file.path(
      motif_db_out_dir,
      sprintf("%s_%s_%s_K%s.csv", motif_db_prefix, db_tag, mode, k_tag)
    )
    readr::write_csv(motif_db_updated, out_path)
    if (isTRUE(verbose)) .log_inform("Saved motif DB: {out_path}")
  }

  list(
    motif_stats = stats,
    similarity_data = sim_data,
    clusters_data = cl_data,
    clusters_table_data = cluster_df,
    target_clusters_data = target_clusters,
    similarity_hybrid = sim_hybrid,
    clusters_hybrid = cl_hybrid,
    target_clusters_hybrid = if (exists("target_clusters_h", inherits = FALSE)) target_clusters_h else NULL,
    motif_db = motif_db_updated
  )
}

#' Pre-run motif clustering only if best-K grid outputs are missing
#'
#' @param fp_aligned Output from align_footprints().
#' @param out_dir,base_dir Output directory or base dir (tf_motif_clustering).
#' @param mode One of "both", "data", "hybrid".
#' @param ... Passed to run_fp_motif_clustering().
#' @return List with `skipped` flag and, if run, the pre-run result.
#' @noRd
run_fp_motif_clustering_pre_if_needed <- function(
  fp_aligned,
  out_dir = NULL,
  base_dir = NULL,
  mode = c("both", "data", "hybrid"),
  ...
) {
  mode <- match.arg(mode)
  if (is.null(out_dir)) {
    if (is.null(base_dir) || !nzchar(base_dir)) {
      cli::cli_abort("Provide out_dir or base_dir.")
    }
    out_dir <- file.path(base_dir, "tf_motif_clustering")
  }
  stopifnot(is.character(out_dir), length(out_dir) == 1L)
  data_ok <- file.exists(file.path(out_dir, "qc_best_k_scores_data.csv"))
  hybrid_ok <- file.exists(file.path(out_dir, "qc_best_k_scores_hybrid.csv"))
  if (mode == "data" && data_ok) {
    if (exists(".log_inform")) .log_inform("Pre-run: cached K grid found (data). Skipping.")
    return(list(skipped = TRUE, result = NULL))
  }
  if (mode == "hybrid" && hybrid_ok) {
    if (exists(".log_inform")) .log_inform("Pre-run: cached K grid found (hybrid). Skipping.")
    return(list(skipped = TRUE, result = NULL))
  }
  if (mode == "both" && data_ok && hybrid_ok) {
    if (exists(".log_inform")) .log_inform("Pre-run: cached K grids found (data + hybrid). Skipping.")
    return(list(skipped = TRUE, result = NULL))
  }
  res <- run_fp_motif_clustering(
    fp_aligned = fp_aligned,
    out_dir = out_dir,
    base_dir = base_dir,
    mode = mode,
    run_mode = "pre",
    ...
  )
  list(skipped = FALSE, result = res)
}
