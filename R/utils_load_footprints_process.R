#' Trim a single leading "_" and rename files accordingly
#' @export
fp_manifest_trim <- function(df, rename_files = TRUE, verbose = TRUE) {
  req <- c("motif", "score", "bound", "annot")
  miss <- setdiff(req, names(df)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

  strip1 <- function(s) sub("^_", "", s)
  fix_path <- function(p) {
    ifelse(is.na(p) | p == "", p, file.path(dirname(p), strip1(basename(p))))
  }
  do_rename <- function(old, new) {
    idx <- !is.na(old) & old != new
    if (!any(idx) || !rename_files) return(invisible(NULL))
    ok_src  <- file.exists(old[idx])
    tgt_ok  <- !file.exists(new[idx])
    todo    <- which(idx)[ok_src & tgt_ok]
    res <- if (length(todo)) file.rename(old[todo], new[todo]) else logical(0)
    if (verbose) {
      if (any(idx & !ok_src)) warning("Source missing: ", paste(old[idx & !ok_src], collapse = "; "))
      if (any(idx & !tgt_ok)) warning("Target exists: ",  paste(new[idx & !tgt_ok], collapse = "; "))
      if (length(res) && !all(res)) warning("Some renames failed.")
    }
  }

  df$motif <- strip1(df$motif)

  for (col in c("score","bound","annot")) {
    old <- df[[col]]
    new <- fix_path(old)
    do_rename(old, new)
    df[[col]] <- new
  }
  df
}

#' Trim leading "_" in all annotation CSVs listed in a manifest
#'
#' Reads each path in `manifest[[annot_col]]`, removes exactly one leading
#' underscore from the `motifs` column, and overwrites the file.
#'
#' @param manifest Tibble/data.frame with column `annot` (paths to CSVs).
#' @param annot_col Column name holding annotation CSV paths. Default: "annot".
#' @param motif_col Column name inside each CSV to fix. Default: "motifs".
#' @param n_workers Integer; if >1 and future.apply is available, process in parallel.
#' @param verbose Logical; print progress.
#' @return Tibble with path, n_rows, n_fixed (rows changed), and status.
#' @export
fp_manifest_trim_annots <- function(manifest,
                                 annot_col = "annot",
                                 motif_col = "motifs",
                                 n_workers = 1L,
                                 verbose = TRUE) {
  if (!is.data.frame(manifest) || !annot_col %in% names(manifest)) {
    cli::cli_abort("`manifest` must contain column {.val {annot_col}}.")
  }

  paths <- unique(stats::na.omit(manifest[[annot_col]]))
  if (!length(paths)) {
    cli::cli_abort("No valid paths found in {.val {annot_col}}.")
  }

  # worker fn
  fix_one <- function(p) {
    if (!file.exists(p)) {
      return(tibble::tibble(path = p, n_rows = NA_integer_, n_fixed = NA_integer_, status = "missing"))
    }
    df <- tryCatch(readr::read_csv(p, show_col_types = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) {
      return(tibble::tibble(path = p, n_rows = NA_integer_, n_fixed = NA_integer_, status = "read_error"))
    }
    if (!motif_col %in% names(df)) {
      return(tibble::tibble(path = p, n_rows = nrow(df), n_fixed = NA_integer_, status = "no_motif_col"))
    }
    old <- df[[motif_col]]
    new <- sub("^_", "", old)
    n_fixed <- sum(!is.na(old) & !is.na(new) & old != new)
    df[[motif_col]] <- new

    # overwrite
    ok <- tryCatch({ readr::write_csv(df, p); TRUE }, error = function(e) FALSE)
    tibble::tibble(path = p, n_rows = nrow(df), n_fixed = n_fixed, status = if (ok) "ok" else "write_error")
  }

  if (verbose) cli::cli_inform("Fixing {length(paths)} annotation file{?s}...")

  run_parallel <- n_workers > 1L &&
    requireNamespace("future.apply", quietly = TRUE)
  if (run_parallel) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_workers)
    res <- future.apply::future_lapply(paths, fix_one, future.seed = FALSE)
  } else {
    res <- lapply(paths, fix_one)
  }

  out <- tibble::as_tibble(dplyr::bind_rows(res))
  if (verbose) {
    ok_n <- sum(out$status == "ok", na.rm = TRUE)
    cli::cli_inform("Done: {ok_n}/{nrow(out)} OK. Changed rows (total): {sum(out$n_fixed, na.rm = TRUE)}.")
  }
  out
}


#' @export
align_footprints <- function(
    fp_filtered_manifest,
    mid_slop         = 10L,    # midpoint tolerance (bp)
    round_digits     = 1L,     # round scores before comparing vectors
    score_match_pct  = 0.8,    # fraction of samples that must match (<=1, >=0)
    verbose          = TRUE,
    threads          = max(1L, parallel::detectCores(TRUE))
) {
  stopifnot(is.data.frame(fp_filtered_manifest),
            all(c("motif","score","bound","annot") %in% names(fp_filtered_manifest)))
  # allow either n_rows or n_peaks
  count_col <- if ("n_rows" %in% names(fp_filtered_manifest)) {
    "n_rows"
  } else if ("n_peaks" %in% names(fp_filtered_manifest)) {
    "n_peaks"
  } else {
    cli::cli_abort("Manifest must contain either 'n_rows' or 'n_peaks'.")
  }

  stopifnot(is.numeric(score_match_pct), length(score_match_pct) == 1L,
            score_match_pct > 0 && score_match_pct <= 1)

  # Use base filtering to avoid dplyr data mask quirks
  use <- fp_filtered_manifest[
    !is.na(fp_filtered_manifest[[count_col]]) & fp_filtered_manifest[[count_col]] > 0 &
      !is.na(fp_filtered_manifest$score) & file.exists(fp_filtered_manifest$score) &
      !is.na(fp_filtered_manifest$bound) & file.exists(fp_filtered_manifest$bound) &
      !is.na(fp_filtered_manifest$annot) & file.exists(fp_filtered_manifest$annot),
    , drop = FALSE
  ]

  if (!nrow(use)) {
    if (verbose) cli::cli_inform("No non-empty motifs with existing files.")
    return(list(
      fp_score = tibble::tibble(),
      fp_bound = tibble::tibble(),
      fp_annotation = tibble::tibble(),
      id_map = tibble::tibble()
    ))
  }

  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(threads)

  fread_fast <- function(p) data.table::fread(p, nThread = threads, showProgress = FALSE)

  # ---- readers ----
  read_score <- function(p) {
    dt <- fread_fast(p)
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    ids <- setdiff(names(dt), "peak_ID")
    for (j in ids) data.table::set(dt, j = j, value = as.numeric(dt[[j]]))
    dt[]
  }
  read_bound <- function(p) {
    dt <- fread_fast(p)
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    ids <- setdiff(names(dt), "peak_ID")
    for (j in ids) data.table::set(dt, j = j, value = as.integer(dt[[j]]))
    dt[]
  }
  read_annot <- function(p) {
    dt <- fread_fast(p)
    if (!"fp_peak"  %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID",   "fp_peak")
    if (!"atac_peak"%in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"   %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    dt[]
  }

  if (verbose) cli::cli_inform("Loading & combining ({nrow(use)} motifs)...")
  score_dt <- data.table::rbindlist(lapply(use$score,  read_score), use.names = TRUE, fill = TRUE)
  bound_dt <- data.table::rbindlist(lapply(use$bound,  read_bound), use.names = TRUE, fill = TRUE)
  annot_dt <- data.table::rbindlist(lapply(use$annot,  read_annot), use.names = TRUE, fill = TRUE)
  # Keep annot_dt EXACT as loaded (row count + order).

  # De-dupe score/bound by old peak id defensively
  if ("peak_ID" %in% names(score_dt)) score_dt <- score_dt[!duplicated(peak_ID)]
  if ("peak_ID" %in% names(bound_dt)) bound_dt <- bound_dt[!duplicated(peak_ID)]

  # ---- assert 1:1 fp_peak -> atac_peak ----
  map_fp_atac <- unique(annot_dt[, .(fp_peak, atac_peak)])
  if (anyNA(map_fp_atac$atac_peak)) stop("Found NA atac_peak in annotation; this function assumes all exist.")
  dup <- map_fp_atac[, .N, by = fp_peak][N > 1L]
  if (nrow(dup)) stop("Some fp_peak map to multiple atac_peak; expected one-to-one mapping.")

  # ---- parse genomic coords from score_dt$peak_ID ----
  parse_peak <- function(x) {
    chr <- sub(":.*$", "", x, perl = TRUE)
    se  <- sub("^[^:]+:", "", x, perl = TRUE)
    st  <- as.integer(sub("-.*$", "", se,  perl = TRUE))
    en  <- as.integer(sub("^.*-", "", se,  perl = TRUE))
    list(chr = chr, start = st, end = en)
  }
  pp <- parse_peak(score_dt$peak_ID)
  score_dt[, `:=`(chr = pp$chr, start = pp$start, end = pp$end)]
  score_dt[, mid := as.integer((start + end) %/% 2)]

  # ---- add atac_peak to score_dt (non-destructive) ----
  data.table::setkey(map_fp_atac, fp_peak)
  score_dt[map_fp_atac, atac_peak := i.atac_peak, on = .(peak_ID = fp_peak)]
  stopifnot(!anyNA(score_dt$atac_peak))

  # ---- signature of rounded scores (temp cols) ----
  sample_cols <- setdiff(names(score_dt), c("peak_ID","atac_peak","chr","start","end","mid"))
  if (!length(sample_cols)) stop("No sample columns detected in fp_score.")
  mult <- 10^round_digits
  tmp_cols <- paste0("__sig__", seq_along(sample_cols))
  score_dt[, (tmp_cols) := lapply(.SD, function(x) as.integer(round(x * mult))), .SDcols = sample_cols]

  # ---- assign sig_id (EXACT vs FUZZY) within (atac_peak, chr) ----
  if (score_match_pct >= 0.9999) {
    data.table::setkeyv(score_dt, c("atac_peak","chr", tmp_cols, "mid"))
    score_dt[, sig_id := data.table::rleidv(.SD), by = .(atac_peak, chr), .SDcols = tmp_cols]
  } else {
    .eq_counts <- function(A, v) {
      B <- matrix(v, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
      E <- (A == B) | (is.na(A) & is.na(B))
      rowSums(E)
    }
    k_req <- ceiling(score_match_pct * length(tmp_cols))

    # One task per atac_peak (no splitting within a peak)
    score_dt[, grp__ := .GRP, by = .(atac_peak)]
    splits <- split(seq_len(nrow(score_dt)), score_dt$grp__, drop = TRUE)

    # Precompute signature matrix once to reduce serialization
    S <- as.matrix(score_dt[, ..tmp_cols])

    # Component assigner (unchanged logic)
    assign_components <- function(M, k_req) {
      n <- nrow(M)
      if (n == 0L) return(integer(0))
      ids <- integer(n)   # 0 = unassigned
      comp <- 0L
      for (i in seq_len(n)) {
        if (ids[i] != 0L) next
        comp <- comp + 1L
        ids[i] <- comp
        q <- i
        while (length(q)) {
          a <- q[1]; q <- q[-1]
          unassigned <- which(ids == 0L)
          if (!length(unassigned)) break
          eqc <- .eq_counts(M[unassigned, , drop = FALSE], M[a, ])
          to_add <- unassigned[eqc >= k_req]
          if (length(to_add)) {
            ids[to_add] <- comp
            q <- c(q, to_add)
          }
        }
      }
      ids
    }

    # Parallelize across atac_peak groups (no intra-peak splitting)
    has_future <- requireNamespace("future.apply", quietly = TRUE)
    apply_fun <- if (has_future) future.apply::future_lapply else lapply
    if (has_future) {
      old_plan <- future::plan()
      on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)
      workers <- max(1L, min(threads, length(splits)))
      # Prefer multicore on Unix when possible; fallback to multisession otherwise
      if (.Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") == "") {
        future::plan(future::multicore, workers = workers)
      } else {
        future::plan(future::multisession, workers = workers)
      }
    }

    # Largest groups first so big peaks spread across workers
    ord <- order(vapply(splits, length, integer(1L)), decreasing = TRUE)
    splits <- splits[ord]

    comp_list <- apply_fun(
      X = splits,
      FUN = function(idx) {
        M <- S[idx, , drop = FALSE]
        assign_components(M, k_req)
      },
      future.scheduling = 2,
      future.seed = TRUE
    )

    # Reassemble into a single sig_id vector, making ids unique across peaks
    sig_id_vec <- integer(nrow(score_dt))
    next_id <- 1L
    for (g in seq_along(splits)) {
      ids <- comp_list[[g]]
      sig_id_vec[splits[[g]]] <- ids + (next_id - 1L)
      next_id <- next_id + max(ids, 0L)
    }

    score_dt[, sig_id := sig_id_vec]
    score_dt[, grp__ := NULL]
  }

  # ---- connected components on mid±slop within (atac_peak, chr, sig_id) ----
  score_dt[, `:=`(s_exp = mid - mid_slop, e_exp = mid + mid_slop)]
  data.table::setorder(score_dt, atac_peak, chr, sig_id, s_exp, e_exp)
  score_dt[, comp_id := {
    e <- e_exp; s <- s_exp
    if (.N <= 1L) 1L else {
      prev <- c(-Inf, cummax(e[-.N]))
      cumsum(s > prev)
    }
  }, by = .(atac_peak, chr, sig_id)]

  # ---- compute merged coords & new IDs IN-PLACE ----
  score_dt[, `:=`(
    new_start  = min(start),
    new_end    = max(end),
    group_size = .N
  ), by = .(atac_peak, chr, sig_id, comp_id)]
  score_dt[, new_peak_ID := paste0(chr, ":", new_start, "-", new_end)]

  # ---- unique old->new mapping (per old peak) ----
  map_old_new <- unique(score_dt[, .(peak_ID_old = peak_ID, atac_peak, new_peak_ID, group_size)])
  data.table::setkey(map_old_new, peak_ID_old)

  # ---- rewrite IDs in place ----
  score_dt[map_old_new, peak_ID := i.new_peak_ID, on = .(peak_ID = peak_ID_old)]
  score_dt <- score_dt[, c("peak_ID", sample_cols), with = FALSE][!duplicated(peak_ID)]

  if (!"peak_ID" %in% names(bound_dt) && "fp_peak" %in% names(bound_dt))
    data.table::setnames(bound_dt, "fp_peak", "peak_ID")
  bound_dt[map_old_new, peak_ID := i.new_peak_ID, on = .(peak_ID = peak_ID_old)]
  bound_dt <- bound_dt[!duplicated(peak_ID)]

  annot_dt[, fp_peak_bak := fp_peak]
  annot_dt[map_old_new, fp_peak := i.new_peak_ID, on = .(fp_peak = peak_ID_old)]

  # ===== ALIGN EVERYTHING TO annot_dt's ROW ORDER =====
  i_ids <- data.table::data.table(peak_ID = annot_dt$fp_peak)
  fp_score_out <- merge(i_ids, score_dt, by = "peak_ID", all.x = TRUE, sort = FALSE)
  fp_bound_out <- merge(i_ids, bound_dt, by = "peak_ID", all.x = TRUE, sort = FALSE)

  gs_new <- map_old_new[, .(group_size = max(group_size, na.rm = TRUE)), by = .(new_peak_ID, atac_peak)]
  data.table::setnames(gs_new, "new_peak_ID", "peak_ID")
  id_map_aligned <- data.table::data.table(
    peak_ID      = annot_dt$fp_peak,
    fp_peak_bak  = annot_dt$fp_peak_bak,
    atac_peak    = annot_dt$atac_peak
  )[gs_new, group_size := i.group_size, on = .(peak_ID)]

  # Final sanity
  stopifnot(
    nrow(fp_score_out)      == nrow(annot_dt),
    nrow(fp_bound_out)      == nrow(annot_dt),
    nrow(id_map_aligned)    == nrow(annot_dt),
    identical(fp_score_out$peak_ID, annot_dt$fp_peak),
    identical(fp_bound_out$peak_ID, annot_dt$fp_peak),
    identical(id_map_aligned$peak_ID, annot_dt$fp_peak)
  )

  list(
    fp_score      = tibble::as_tibble(fp_score_out),
    fp_bound      = tibble::as_tibble(fp_bound_out),
    fp_annotation = tibble::as_tibble(annot_dt),
    id_map        = tibble::as_tibble(id_map_aligned)
  )
}



# Quantile normalize *by unique peak_ID*, then broadcast to all rows.
# - Keeps the input row order and N exactly the same.
# - Unique weighting prevents duplicated peak_ID from over-influencing the reference.
# - Tie-aware assignment within each column (avg of ref values across the tie block).
#' @export
qn_footprints <- function(fp_tbl, id_col = "peak_ID", tie_method = c("average","first")) {
  tie_method <- match.arg(tie_method)
  stopifnot(id_col %in% names(fp_tbl))

  samp <- setdiff(names(fp_tbl), id_col)
  if (!length(samp)) return(fp_tbl)

  # Original order mapping
  ids_all <- fp_tbl[[id_col]]
  # First occurrence of each peak_ID (we assume duplicates are identical rows from the alignment step)
  uniq_first <- !duplicated(ids_all)
  ids_unique <- ids_all[uniq_first]
  # Group index for broadcasting back
  grp <- match(ids_all, ids_unique)

  # Unique matrix (one row per unique peak)
  X <- as.matrix(fp_tbl[uniq_first, samp, drop = FALSE])
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  if (n == 0L || p == 0L) return(fp_tbl)

  # Build reference profile without constructing a full sorted matrix
  ref_sum   <- numeric(n)
  ref_count <- integer(n)

  # We'll also store per-column orders & non-NA counts for reuse
  ord_list <- vector("list", p)
  m_list   <- integer(p)

  for (j in seq_len(p)) {
    v <- X[, j]
    o <- order(v, na.last = TRUE)
    # number of non-NA values
    m <- sum(!is.na(v))
    ord_list[[j]] <- o
    m_list[j]     <- m
    if (m > 0L) {
      # Add to reference accumulators (first m positions)
      # we only need v sorted for the count (ref is just average of values at each rank)
      ref_sum[seq_len(m)]   <- ref_sum[seq_len(m)]   + v[o][seq_len(m)]
      ref_count[seq_len(m)] <- ref_count[seq_len(m)] + 1L
    }
  }

  # Mean across columns for each rank (ignore ranks with 0 count)
  ref <- ref_sum
  nz  <- ref_count > 0L
  ref[nz] <- ref_sum[nz] / ref_count[nz]
  ref[!nz] <- NA_real_

  # Map each column to the reference (tie-aware if requested)
  Xqn <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    v <- X[, j]
    o <- ord_list[[j]]
    m <- m_list[j]
    if (m == 0L) next

    if (tie_method == "first") {
      # simple, fast: one-to-one rank mapping
      Xqn[o[seq_len(m)], j] <- ref[seq_len(m)]
    } else {
      # tie-aware: average the ref within each tie block
      s <- v[o[seq_len(m)]]
      r <- rle(s)
      ends   <- cumsum(r$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)
      # assign the mean(ref[start:end]) to all positions in that tie block
      for (k in seq_along(ends)) {
        idx <- starts[k]:ends[k]
        Xqn[o[idx], j] <- mean(ref[idx], na.rm = TRUE)
      }
    }
  }

  # Broadcast normalized unique values back to all rows WITHOUT materializing a big matrix
  out <- fp_tbl
  for (j in seq_len(p)) {
    out[[samp[j]]] <- Xqn[, j][grp]
  }
  out
}



# Save aligned/normalized FP results per motif and return a manifest
#' @export
save_footprints <- function(
    fp_aligned_normalized,
    out_dir  = NULL,
    threads  = max(1L, data.table::getDTthreads()),
    verbose  = TRUE
) {
  # pull tables
  score_dt <- data.table::as.data.table(fp_aligned_normalized$fp_score)
  bound_dt <- data.table::as.data.table(fp_aligned_normalized$fp_bound)
  annot_dt <- data.table::as.data.table(fp_aligned_normalized$fp_annotation)

  # sanity checks: same N and order
  stopifnot(nrow(score_dt) == nrow(bound_dt),
            nrow(bound_dt) == nrow(annot_dt),
            identical(score_dt$peak_ID, bound_dt$peak_ID),
            identical(score_dt$peak_ID, annot_dt$fp_peak),
            "motifs" %in% names(annot_dt))

  # ensure output dir
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  data.table::setDTthreads(threads)

  motifs <- unique(annot_dt$motifs)
  n_m <- length(motifs)

  # pre-allocate manifest fields
  n_rows_vec   <- integer(n_m)
  score_paths  <- character(n_m)
  bound_paths  <- character(n_m)
  annot_paths  <- character(n_m)

  if (verbose) cli::cli_inform("Writing per-motif CSVs to {out_dir} ({n_m} motif{?s}) ...")

  for (i in seq_len(n_m)) {
    m <- motifs[i]
    idx <- which(annot_dt$motifs == m)

    # slice using the aligned row order
    sc <- score_dt[idx]
    bo <- bound_dt[idx]
    an <- annot_dt[idx]

    # file paths (match prior naming convention)
    base        <- file.path(out_dir, m)
    f_score     <- paste0(base, "_score.csv")
    f_bound     <- paste0(base, "_bound.csv")
    f_annot     <- paste0(base, "_annot.csv")

    # write
    data.table::fwrite(sc, f_score, nThread = threads, bom = FALSE)
    data.table::fwrite(bo, f_bound, nThread = threads, bom = FALSE)
    data.table::fwrite(an, f_annot, nThread = threads, bom = FALSE)

    # manifest rows
    n_rows_vec[i]  <- length(idx)
    score_paths[i] <- f_score
    bound_paths[i] <- f_bound
    annot_paths[i] <- f_annot
  }

  if (verbose) cli::cli_inform("Done. Creating manifest tibble ...")

  manifest <- tibble::tibble(
    motif = motifs,
    n_rows = n_rows_vec,
    score  = score_paths,
    bound  = bound_paths,
    annot  = annot_paths
  )

  manifest
}

# per-motif filtering
# BEGIN EDIT: per-motif processor + parallel driver (uses update_fp_bound_by_overlap)

#' Process ONE motif from fp_manifest (strict set) and write filtered CSVs
#'
#' Reads the three per-motif files (score/bound/annotation), builds a strict GRN set
#' using provided `build_args`, applies the standard filtering/update steps, and writes
#' filtered outputs to `out_dir` as "<motif>_{score,bound,annotation}.csv".
#'
#' Required columns in `entry`: motif, score, bound, annot (absolute or relative paths).
#' `update_fp_bound_by_overlap()` is used as-is.
#'
#' @param entry Single-row tibble/data.frame from fp_manifest.
#' @param out_dir Output directory (e.g., "inst/extdata/fp_filtered_jaspar2024").
#' @param build_args Named list of args passed to `build_grn_set()` in addition to fp_*.
#'                   Should include: atac_score, atac_overlap, rna, metadata, tf_list,
#'                   motif_db, label_col, expected_n.
#' @param threshold_tf_expr Numeric; forwarded to `add_tf_expr_flags()`.
#' @param skip_existing If TRUE and all three outputs exist, skip work.
#' @param verbose Print progress.
#' @return Tibble with motif, n_rows (filtered fp_bound rows), and output paths.
process_one_motif_from_manifest <- function(entry,
                                            out_dir,
                                            build_args,
                                            threshold_tf_expr = 10,
                                            skip_existing = TRUE,
                                            verbose = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  if (!all(need_cols %in% names(entry))) {
    cli::cli_abort("`entry` must contain columns: {cli::fmt_cols(setdiff(need_cols, names(entry)))}")
  }
  m <- as.character(entry$motif)[1]
  f_score <- as.character(entry$score)[1]
  f_bound <- as.character(entry$bound)[1]
  f_annot <- as.character(entry$annot)[1]

  if (!file.exists(f_score) || !file.exists(f_bound) || !file.exists(f_annot)) {
    cli::cli_abort("Missing file(s) for motif {m}:\n  score={f_score}\n  bound={f_bound}\n  annot={f_annot}")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_score <- file.path(out_dir, paste0(m, "_score.csv"))
  out_bound <- file.path(out_dir, paste0(m, "_bound.csv"))
  out_annot <- file.path(out_dir, paste0(m, "_annotation.csv"))

  if (skip_existing && file.exists(out_score) && file.exists(out_bound) && file.exists(out_annot)) {
    return(tibble::tibble(motif = m, n_rows = NA_integer_,
                          score = out_score, bound = out_bound, annot = out_annot))
  }

  if (verbose) cli::cli_inform("Processing motif {m} ...")

  fp_score <- readr::read_csv(f_score, show_col_types = FALSE)
  fp_bound <- readr::read_csv(f_bound, show_col_types = FALSE)
  fp_annot <- readr::read_csv(f_annot, show_col_types = FALSE)

  if (!all(c("fp_peak","atac_peak","motifs") %in% names(fp_annot))) {
    if (all(c("peak_ID","peak_ATAC","TFBS_name") %in% names(fp_annot))) {
      fp_annot <- dplyr::rename(fp_annot, fp_peak = "peak_ID", atac_peak = "peak_ATAC", motifs = "TFBS_name")
    } else {
      cli::cli_abort("Annotation for {m} must have fp_peak/atac_peak/motifs (or peak_ID/peak_ATAC/TFBS_name).")
    }
  }

  ss <- do.call(build_grn_set, c(list(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annot
  ), build_args))

  # --- EARLY-EXIT GUARDS (no writes if empty) -----------------
  empty_row <- function() tibble::tibble(motif = m, n_rows = 0L,
                                         score = NA_character_,
                                         bound = NA_character_,
                                         annot = NA_character_)

  ss <- filter_by_min_bound(ss, min_bound = 1L)
  if (nrow(ss$fp_bound) == 0L) { if (verbose) cli::cli_inform("Motif {m}: empty after min-bound filter; skipping."); return(empty_row()) }

  ss <- update_fp_bound_by_overlap(ss)
  ss <- filter_by_min_bound(ss, min_bound = 1L)
  if (nrow(ss$fp_bound) == 0L) { if (verbose) cli::cli_inform("Motif {m}: empty after overlap gating; skipping."); return(empty_row()) }

  ss <- add_tf_expr_flags(ss, threshold = threshold_tf_expr)
  ss2 <- update_fp_bound_by_tf_expr(ss, group_size = 1L)
  ss2 <- filter_by_min_bound(ss2, min_bound = 1L)
  if (nrow(ss2$fp_bound) == 0L) { if (verbose) cli::cli_inform("Motif {m}: empty after TF expression gating; skipping."); return(empty_row()) }
  # ------------------------------------------------------------

  # Only compute correlations if we still have rows
  # ss2 <- annotate_fp_tf_corr_one_motif(ss2, cor_method = "pearson", min_non_na = 5L) # do not calculate correlation for now, will calculate after quantile normalization

  out_score_tbl <- ss2$fp_score
  out_bound_tbl <- ss2$fp_bound
  out_annot_tbl <- ss2$fp_annotation

  if (!"peak_ID" %in% names(out_score_tbl) && "fp_peak" %in% names(out_score_tbl)) {
    out_score_tbl <- dplyr::rename(out_score_tbl, peak_ID = "fp_peak")
  }
  if (!"peak_ID" %in% names(out_bound_tbl) && "fp_peak" %in% names(out_bound_tbl)) {
    out_bound_tbl <- dplyr::rename(out_bound_tbl, peak_ID = "fp_peak")
  }
  if (!all(c("fp_peak","atac_peak","motifs") %in% names(out_annot_tbl))) {
    out_annot_tbl <- out_annot_tbl |>
      dplyr::rename_with(~"fp_peak",  dplyr::any_of("peak_ID")) |>
      dplyr::rename_with(~"atac_peak",dplyr::any_of("peak_ATAC")) |>
      dplyr::rename_with(~"motifs",   dplyr::any_of("TFBS_name"))
  }

  readr::write_csv(out_score_tbl, out_score)
  readr::write_csv(out_bound_tbl, out_bound)
  readr::write_csv(out_annot_tbl, out_annot)

  tibble::tibble(motif = m,
                 n_rows = nrow(out_bound_tbl),
                 score  = out_score,
                 bound  = out_bound,
                 annot  = out_annot)
}

#' Parallel driver: map `process_one_motif_from_manifest()` over fp_manifest
#'
#' @param fp_manifest Tibble with columns motif, score, bound, annot.
#' @param out_dir Output directory for filtered CSVs.
#' @param build_args Named list passed to `build_grn_set()` (see above).
#' @param motif_ids Optional subset of motifs to process.
#' @param n_workers Parallel workers (motif-level). Default: all logical cores.
#' @param set_plan Whether to set/reset a future plan automatically.
#' @param skip_existing If TRUE, skip motifs whose three outputs already exist.
#' @param threshold_tf_expr Numeric threshold for `add_tf_expr_flags()`.
#' @param verbose Print progress.
#' @return Tibble manifest of processed motifs.
#'
#' @export
filter_footprints <- function(fp_manifest,
                                       out_dir,
                                       build_args,
                                       motif_ids = NULL,
                                       n_workers = max(1L, parallel::detectCores(TRUE)),
                                       set_plan = TRUE,
                                       skip_existing = TRUE,
                                       threshold_tf_expr = 10,
                                       verbose = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  if (!isTRUE(is.data.frame(fp_manifest)) || !all(need_cols %in% names(fp_manifest))) {
    cli::cli_abort("`fp_manifest` must be a data.frame with columns: {cli::fmt_cols(need_cols)}")
  }

  rows <- fp_manifest
  if (!is.null(motif_ids)) {
    rows <- dplyr::semi_join(rows, tibble::tibble(motif = motif_ids), by = "motif")
    if (!nrow(rows)) cli::cli_abort("No matching motifs found in `fp_manifest` for the requested subset.")
  }

  if (verbose) cli::cli_inform("Submitting {nrow(rows)} motif(s) to {n_workers} worker(s)...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1") ||
      (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !is_rstudio &&
        isTRUE(requireNamespace("parallelly", quietly = TRUE)) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore,   workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
  }

  res <- future.apply::future_lapply(seq_len(nrow(rows)), function(i) {
    entry <- rows[i, , drop = FALSE]
    process_one_motif_from_manifest(
      entry               = entry,
      out_dir             = out_dir,
      build_args          = build_args,
      threshold_tf_expr   = threshold_tf_expr,
      skip_existing       = skip_existing,
      verbose             = FALSE
    )
  }, future.seed = FALSE)

  tibble::as_tibble(dplyr::bind_rows(res))
}




#' Combine FP correlation tables (split dimers), run BH FDR once across all,
#' sanity-check p vs p_adj, filter, and align score/bound to the kept annotations.
#'
#' Robust, minimal pipeline:
#' 1) Normalizes columns when reading each file.
#' 2) Splits comma-separated TF/r/p columns row-wise (handles dimers) with safe padding.
#' 3) Runs **one** Benjamini–Hochberg adjustment across the combined annotations.
#' 4) Enforces `p_adj >= p` and reports monotonicity stats.
#' 5) Filters by `r_thr` (positive) and `p_thr` (FDR).
#' 6) Reads score/bound only for peaks that pass, de-dups on `peak_ID`, and aligns order.
#'
#' @param fp_filtered_manifest Tibble/data.frame with columns:
#'   `motif, n_rows, score, bound, annot` (CSV paths).
#' @param p_thr Numeric FDR threshold (default 0.05).
#' @param r_thr Numeric correlation threshold (keep r > r_thr; default 0.3).
#' @param verbose Logical; print progress (default TRUE).
#' @param threads Integer; data.table threads (default = all detected cores).
#' @param output_bed Optional directory path. If provided, the function will write per-TF
#'   BED-like files: `<TF>_all.bed`, `<TF>_bound.bed`, `<TF>_unbound.bed` (no headers),
#'   and `<TF>_overview.txt` (headered; same columns as `_all.bed` plus a `_bound` 0/1 column).
#'   Columns: TFBS_chr, TFBS_start, TFBS_end, TFBS_name, peak_chr, peak_start, peak_end,
#'   TF, corr_fp_tf_r, corr_fp_tf_p, corr_fp_tf_p_adj; bound/unbound defined by r/p filters.
#'
#' @return A list of tibbles: `fp_score`, `fp_bound`, `fp_annotation`
#'         (all aligned to the same row order).
#'
#' @examples
#' # res <- combine_filtered_fp_tables_split_dimers_simple(fp_corr_manifest)
combine_filtered_fp_tables_split_dimers_simple <- function(
    fp_filtered_manifest,
    p_thr   = 0.05,
    r_thr   = 0.3,
    verbose = TRUE,
    threads = max(1L, parallel::detectCores(TRUE)),
    # BEGIN EDIT: add optional output directory for BED exports
    output_bed = NULL
    # END EDIT
) {
  # ---- checks ----
  if (!is.data.frame(fp_filtered_manifest) ||
      !all(c("motif","n_rows","score","bound","annot") %in% names(fp_filtered_manifest))) {
    cli::cli_abort("`fp_filtered_manifest` must have columns: motif, n_rows, score, bound, annot.")
  }

  use <- fp_filtered_manifest |>
    dplyr::filter(
      !is.na(.data$n_rows), .data$n_rows > 0,
      !is.na(.data$score), file.exists(.data$score),
      !is.na(.data$bound), file.exists(.data$bound),
      !is.na(.data$annot), file.exists(.data$annot)
    )

  if (!nrow(use)) {
    if (verbose) cli::cli_inform("No non-empty motifs with existing files; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }
  if (verbose) cli::cli_inform("Combining {nrow(use)} motif{?s} ...")

  # ---- threading ----
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(threads))

  fread_fast <- function(path) data.table::fread(path, nThread = threads, showProgress = FALSE)

  # ---- helpers ----
  split_commas <- function(x) {
    x <- as.character(x)
    lapply(x, function(s) {
      if (is.na(s) || !nzchar(s)) character(0)
      else strsplit(s, "\\s*,\\s*", perl = TRUE)[[1]]
    })
  }
  pad_to <- function(lst, lens, fill) {
    out <- vector("list", length(lst))
    for (i in seq_along(lst)) {
      v <- lst[[i]]
      # recycle singleton
      if (length(v) == 1L && lens[i] > 1L) v <- rep(v, lens[i])
      # right-pad
      if (length(v) < lens[i]) v <- c(v, rep(fill, lens[i] - length(v)))
      out[[i]] <- v
    }
    unlist(out, use.names = FALSE)
  }
  safe_num <- function(x) {
    x <- as.character(x)
    x[!nzchar(x)] <- NA_character_   # treat "" as NA
    suppressWarnings(as.numeric(x))  # quiet "NAs introduced by coercion"
  }

  # ---- read one annotation, normalize names, split dimers ----
  read_annot_split <- function(p) {
    dt <- fread_fast(p)

    # Normalize names to: fp_peak / atac_peak / motifs / tfs / corr_fp_tf_r / corr_fp_tf_p
    if (!"fp_peak"   %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID",   "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"    %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    if (!"tfs"       %in% names(dt)) dt[["tfs"]] <- NA_character_
    if (!"corr_fp_tf_r" %in% names(dt)) dt[["corr_fp_tf_r"]] <- NA_character_
    if (!"corr_fp_tf_p" %in% names(dt)) dt[["corr_fp_tf_p"]] <- NA_character_

    tfs_list <- split_commas(dt[["tfs"]])
    r_list   <- lapply(split_commas(dt[["corr_fp_tf_r"]]), safe_num)
    p_list   <- lapply(split_commas(dt[["corr_fp_tf_p"]]), safe_num)

    # how many repeats per row
    nrep <- pmax(lengths(tfs_list), lengths(r_list), lengths(p_list))
    nrep[nrep == 0L] <- 1L
    idx <- rep.int(seq_len(nrow(dt)), nrep)

    data.table::data.table(
      fp_peak      = dt$fp_peak[idx],
      atac_peak    = dt$atac_peak[idx],
      motifs       = dt$motifs[idx],
      tfs          = pad_to(tfs_list, nrep, NA_character_),
      corr_fp_tf_r = safe_num(pad_to(r_list, nrep, NA_real_)),
      corr_fp_tf_p = safe_num(pad_to(p_list, nrep, NA_real_))
    )
  }

  # ---- 1) combine all annotations ----
  annot_dt <- data.table::rbindlist(lapply(use$annot, read_annot_split),
                                    use.names = TRUE, fill = TRUE)

  # Optional quick NA summary
  if (verbose) {
    n_p_na <- sum(!is.finite(annot_dt$corr_fp_tf_p))
    n_r_na <- sum(!is.finite(annot_dt$corr_fp_tf_r))
    cli::cli_inform("Pre-filter NA counts: p={n_p_na}, r={n_r_na}")
  }

  # ---- sanity: p in [0,1] (raw) ----
  bad_gt1 <- sum(annot_dt$corr_fp_tf_p > 1, na.rm = TRUE)
  bad_lt0 <- sum(annot_dt$corr_fp_tf_p < 0, na.rm = TRUE)
  if (bad_gt1 || bad_lt0) {
    cli::cli_abort("[raw p] Found {bad_gt1} values > 1 and {bad_lt0} values < 0. Please fix upstream.")
  }

  # ---- pre-BH: dedup (fp_peak, tfs) by keeping the FIRST occurrence ----
  if (verbose) cli::cli_inform("Pre-BH: rows before pair-dedup: {nrow(annot_dt)}")

  keep_idx <- !duplicated(annot_dt[, .(fp_peak, tfs)])
  n_drop   <- sum(!keep_idx)

  annot_dt <- annot_dt[keep_idx]

  if (verbose) {
    cli::cli_inform("Pre-BH: dropped {n_drop} duplicate (fp_peak, tfs) pairs (kept first).")
    cli::cli_inform("Pre-BH: rows after pair-dedup: {nrow(annot_dt)}")
  }

  # ---- BH adjust across ALL rows at once (finite p only) ----
  good_p <- is.finite(annot_dt$corr_fp_tf_p)
  annot_dt[, corr_fp_tf_p_adj := NA_real_]
  if (any(good_p)) {
    annot_dt[good_p, corr_fp_tf_p_adj := stats::p.adjust(annot_dt$corr_fp_tf_p[good_p], method = "BH")]
    # enforce theoretical constraint (defensive against FP/rounding)
    annot_dt[good_p, corr_fp_tf_p_adj := pmax(corr_fp_tf_p_adj, corr_fp_tf_p)]
  }

  # ---- post-BH sanity stats ----
  cmp <- data.table::fcase(
    !is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj), NA_character_,
    annot_dt$corr_fp_tf_p_adj >  annot_dt$corr_fp_tf_p, ">",
    annot_dt$corr_fp_tf_p_adj == annot_dt$corr_fp_tf_p, "==",
    annot_dt$corr_fp_tf_p_adj <  annot_dt$corr_fp_tf_p, "<"
  )
  if (verbose) {
    cli::cli_inform(c(
      "BH sanity (p_adj vs p):",
      "  >  = {sum(cmp == '>', na.rm = TRUE)}",
      "  == = {sum(cmp == '==', na.rm = TRUE)}",
      "  <  = {sum(cmp == '<', na.rm = TRUE)}  (after pmax() this should be 0)",
      "  NA = {sum(!is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj))}",
      "  total rows = {nrow(annot_dt)}"
    ))
  }

  # BEGIN EDIT: snapshot pre-filter annotations for BED exports
  annot_pre_filter <- data.table::copy(annot_dt)
  # END EDIT

  # ---- 2) filter: positive r and FDR ----
  keep <- is.finite(annot_dt$corr_fp_tf_r) &
    is.finite(annot_dt$corr_fp_tf_p_adj) &
    (annot_dt$corr_fp_tf_r > r_thr) &
    (annot_dt$corr_fp_tf_p_adj < p_thr)
  annot_dt <- annot_dt[keep]
  if (verbose) {
    cli::cli_inform("After corr/FDR filter: {nrow(annot_dt)} rows (unique peaks: {length(unique(annot_dt$fp_peak))}).")
  }
  if (nrow(annot_dt) == 0L) {
    # BEGIN EDIT: still optionally create empty per-TF files if requested
    if (!is.null(output_bed)) {
      if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)
      if (verbose) cli::cli_inform("No rows passed; requested BED export directory initialized at {output_bed}.")
    }
    # END EDIT
    if (verbose) cli::cli_inform("Nothing passed thresholds; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }

  peaks_keep <- unique(annot_dt$fp_peak)

  # ---- helpers to normalize score/bound ----
  norm_score <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.numeric(dt[[j]])))
    dt[]
  }
  norm_bound <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.integer(safe_num(dt[[j]]))))
    dt[]
  }

  read_score_subset <- function(p) {
    dt <- norm_score(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }
  read_bound_subset <- function(p) {
    dt <- norm_bound(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }

  score_dt <- data.table::rbindlist(lapply(use$score, read_score_subset), use.names = TRUE, fill = TRUE)
  bound_dt <- data.table::rbindlist(lapply(use$bound, read_bound_subset), use.names = TRUE, fill = TRUE)

  # de-dup by first occurrence
  if (nrow(score_dt)) score_dt <- score_dt[!duplicated(peak_ID)]
  if (nrow(bound_dt)) bound_dt <- bound_dt[!duplicated(peak_ID)]

  # ---- 3) align score/bound EXACTLY to annotation order ----
  key_vec <- annot_dt$fp_peak
  key_dt  <- data.table::data.table(peak_ID = key_vec)

  fp_score_dt <- score_dt[key_dt, on = "peak_ID"]  # preserves key_dt order
  fp_bound_dt <- bound_dt[key_dt, on  = "peak_ID"]

  # final sanity: 1:1 and same order
  stopifnot(
    nrow(fp_score_dt) == nrow(annot_dt),
    nrow(fp_bound_dt) == nrow(annot_dt),
    identical(fp_score_dt$peak_ID, key_vec),
    identical(fp_bound_dt$peak_ID, key_vec)
  )

  # BEGIN EDIT: optional per-TF BED & overview export
  if (!is.null(output_bed)) {
    if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)

    .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
    .mk_table <- function(dt) {
      if (!nrow(dt)) return(data.table::data.table())
      dt2 <- dt[!is.na(tfs) & nzchar(tfs) & !is.na(fp_peak) & !is.na(atac_peak)]
      if (!nrow(dt2)) return(data.table::data.table())
      fp_s <- data.table::tstrsplit(dt2$fp_peak, "[:-]", perl = TRUE)
      at_s <- data.table::tstrsplit(dt2$atac_peak, "[:-]", perl = TRUE)
      data.table::data.table(
        TFBS_chr         = fp_s[[1]],
        TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
        TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
        TFBS_name        = dt2$motifs,
        peak_chr         = at_s[[1]],
        peak_start       = suppressWarnings(as.integer(at_s[[2]])),
        peak_end         = suppressWarnings(as.integer(at_s[[3]])),
        TF               = dt2$tfs,
        corr_fp_tf_r     = dt2$corr_fp_tf_r,
        corr_fp_tf_p     = dt2$corr_fp_tf_p,
        corr_fp_tf_p_adj = dt2$corr_fp_tf_p_adj
      )
    }

    all_dt    <- .mk_table(annot_pre_filter)
    bound_dt2 <- .mk_table(annot_pre_filter[keep])

    if (nrow(all_dt))   all_dt[,   .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]
    if (nrow(bound_dt2)) bound_dt2[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]

    tfs_vec <- sort(unique(all_dt$TF))
    for (tf in tfs_vec) {
      tf_slug <- .slug(tf)
      f_all      <- file.path(output_bed, paste0(tf_slug, "_all.bed"))
      f_bound    <- file.path(output_bed, paste0(tf_slug, "_bound.bed"))
      f_unbound  <- file.path(output_bed, paste0(tf_slug, "_unbound.bed"))
      f_overview <- file.path(output_bed, paste0(tf_slug, "_overview.txt"))

      sub_all   <- all_dt[TF == tf]
      sub_bound <- bound_dt2[TF == tf]

      # --- FIX 1: use %in% (not %chin%) ---
      ids_bound   <- if (nrow(sub_bound)) sub_bound$.id else character(0)
      sub_unbound <- sub_all[!(.id %in% ids_bound)]

      # --- FIX 2: drop ".id" (not "._id"); keep=FALSE selection with data.table ---
      data.table::fwrite(sub_all[,   !".id", with = FALSE], f_all,     sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_bound[, !".id", with = FALSE], f_bound,   sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_unbound[, !".id", with = FALSE], f_unbound, sep = "\t", col.names = FALSE, quote = FALSE)

      # headered overview with _bound flag
      if (nrow(sub_all)) {
        sub_all[, `_bound` := as.integer(.id %in% ids_bound)]   # <-- FIX 1 applied here too
        data.table::fwrite(sub_all[, !".id", with = FALSE], f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      } else {
        hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
                 "peak_chr","peak_start","peak_end","TF",
                 "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
        # create a 0-row table with the right header
        empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode="character", length=0)), length(hdr)), hdr))
        data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      }
    }
    if (verbose) cli::cli_inform("Per-TF BED exports written to: {output_bed}")
  }
  # END EDIT


  # ---- return ----
  list(
    fp_score      = tibble::as_tibble(fp_score_dt),
    fp_bound      = tibble::as_tibble(fp_bound_dt),
    fp_annotation = tibble::as_tibble(annot_dt)
  )
}
