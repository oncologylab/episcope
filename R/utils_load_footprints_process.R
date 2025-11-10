#' Build a GRN-aligned set from pre-cleaned inputs (general, package-ready)
#'
#' @param fp_score      tibble with 'peak_ID' + ATAC id columns (numeric)
#' @param fp_bound      tibble with 'peak_ID' + ATAC id columns (0/1 or logical)
#' @param fp_annotation tibble with motif hits (passed through)
#' @param atac_score    tibble with 'atac_peak' + ATAC id columns (numeric)
#' @param atac_overlap  tibble with 'atac_peak' + ATAC id columns (0/1 or logical)
#' @param rna           tibble with 'ensembl_gene_id','HGNC' + RNA columns (already named by ATAC ids)
#' @param metadata      tibble with at least an 'id' column; order defines column order
#' @param tf_list       character vector of TF symbols (passed through)
#' @param motif_db      optional tibble (passed through)
#' @param label_col     optional metadata column to carry along as sample labels (e.g. "strict_match_rna" or "cell_stress_type")
#' @param expected_n    optional integer for info message
#' @return list(fp_score, fp_bound, atac_score, atac_overlap, rna, fp_annotation, tf_list, motif_db,
#'              sample_metadata_used, sample_labels, dropped_ids)
build_grn_set <- function(
    fp_score,
    fp_bound,
    fp_annotation,
    atac_score,
    atac_overlap,
    rna,
    metadata,
    tf_list,
    motif_db   = NULL,
    label_col  = NULL,
    expected_n = NULL
) {
  # ---- checks
  if (!all(c("peak_ID") %in% names(fp_score))) cli::cli_abort("`fp_score` needs 'peak_ID'.")
  if (!all(c("peak_ID") %in% names(fp_bound))) cli::cli_abort("`fp_bound` needs 'peak_ID'.")
  if (!"atac_peak" %in% names(atac_score))     cli::cli_abort("`atac_score` needs 'atac_peak'.")
  if (!"atac_peak" %in% names(atac_overlap))   cli::cli_abort("`atac_overlap` needs 'atac_peak'.")
  if (!all(c("ensembl_gene_id","HGNC") %in% names(rna)))
    cli::cli_abort("`rna` needs 'ensembl_gene_id' and 'HGNC'.")
  if (!"id" %in% names(metadata))
    cli::cli_abort("`metadata` needs an 'id' column.")
  if (!is.null(label_col) && !label_col %in% names(metadata))
    cli::cli_abort("label_col {.val {label_col}} not found in `metadata`.")
  if (!is.character(tf_list)) cli::cli_abort("`tf_list` must be a character vector.")

  # ---- IDs (preserve metadata order; drop NA/dups)
  ids <- metadata$id
  ids <- ids[!is.na(ids)]
  ids <- ids[!duplicated(ids)]
  if (!length(ids)) cli::cli_abort("No usable sample ids found in `metadata$id`.")

  # ---- availability & subsetting
  miss_fp_score <- setdiff(ids, names(fp_score))
  miss_fp_bound <- setdiff(ids, names(fp_bound))
  miss_atac_sc  <- setdiff(ids, names(atac_score))
  miss_atac_ol  <- setdiff(ids, names(atac_overlap))
  miss_rna      <- setdiff(ids, names(rna))

  if (length(miss_fp_score)) message("Dropping ", length(miss_fp_score), " id(s) missing in `fp_score`: ", paste(miss_fp_score, collapse=", "))
  if (length(miss_fp_bound)) message("Dropping ", length(miss_fp_bound), " id(s) missing in `fp_bound`: ", paste(miss_fp_bound, collapse=", "))
  if (length(miss_atac_sc))  message("`atac_score` missing ", length(miss_atac_sc), " id(s): ", paste(miss_atac_sc, collapse=", "))
  if (length(miss_atac_ol))  message("`atac_overlap` missing ", length(miss_atac_ol), " id(s): ", paste(miss_atac_ol, collapse=", "))
  if (length(miss_rna))      message("`rna` missing ", length(miss_rna), " id(s): ", paste(miss_rna, collapse=", "))

  keep_ids <- ids[ids %in% names(fp_score) &
                    ids %in% names(fp_bound) &
                    ids %in% names(rna)]

  if (!length(keep_ids)) cli::cli_abort("No overlapping ids across fp_score/fp_bound/rna after checks.")

  if (!is.null(expected_n) && length(keep_ids) != expected_n) {
    message("Note: aligned sample count = ", length(keep_ids),
            " (expected ", expected_n, "). Proceeding.")
  }


  # BEGIN EDIT: subset atac_* rows to peaks present in filtered fp_annotation
  # rows we keep for ATAC matrices (from filtered fp_annotation)
  peaks_keep_atac <- unique(fp_annotation$atac_peak)
  peaks_keep_atac <- peaks_keep_atac[!is.na(peaks_keep_atac)]

  if (length(peaks_keep_atac) == 0L) {
    cli::cli_abort("No usable 'atac_peak' values found in `fp_annotation` after filtering.")
  }

  # subset rows in atac_* to only those peaks (do this BEFORE column selection)
  atac_score    <- dplyr::semi_join(atac_score,    tibble::tibble(atac_peak = peaks_keep_atac), by = "atac_peak")
  atac_overlap  <- dplyr::semi_join(atac_overlap,  tibble::tibble(atac_peak = peaks_keep_atac), by = "atac_peak")

  # optional: brief info
  message(sprintf("ATAC rows kept by fp_annotation peaks: score=%s, overlap=%s",
                  nrow(atac_score), nrow(atac_overlap)))
  # END EDIT


  # ---- ordered subsets
  fp_score_sub    <- dplyr::select(fp_score, "peak_ID", dplyr::all_of(keep_ids))
  fp_bound_sub    <- dplyr::select(fp_bound, "peak_ID", dplyr::all_of(keep_ids))
  atac_score_sub  <- dplyr::select(atac_score, "atac_peak", dplyr::all_of(intersect(keep_ids, names(atac_score))))
  atac_overlap_sub<- dplyr::select(atac_overlap, "atac_peak", dplyr::all_of(intersect(keep_ids, names(atac_overlap))))
  rna_sub         <- dplyr::select(rna, "ensembl_gene_id","HGNC", dplyr::all_of(keep_ids))

  # ---- metadata used / labels
  meta_used <- dplyr::semi_join(metadata, tibble::tibble(id = keep_ids), by = "id")
  # keep metadata order
  meta_used <- meta_used[match(keep_ids, meta_used$id), , drop = FALSE]
  sample_labels <- if (!is.null(label_col)) meta_used[[label_col]] else NULL

  list(
    fp_score             = fp_score_sub,
    fp_bound             = fp_bound_sub,
    atac_score           = atac_score_sub,
    atac_overlap         = atac_overlap_sub,
    rna                  = rna_sub,
    fp_annotation        = fp_annotation,
    tf_list              = tf_list,
    motif_db             = motif_db,
    sample_metadata_used = meta_used,
    sample_labels        = sample_labels,
    dropped_ids          = setdiff(ids, keep_ids)
  )
}

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

# ── 2) Keep rows with ≥ N bound samples (general, reusable) ────────────────────
# Counts 1s per row in fp_bound across selected sample columns and filters via (1).
#
# - min_bound: integer threshold (your "group size"; default 1)
# - samples: optional character vector of sample IDs to consider; default = all ids
# - na_as_unbound: treat NA as 0 (TRUE) or drop rows with NA (FALSE)
#' @export
filter_by_min_bound <- function(
    grn_set,
    min_bound     = 1L,
    samples       = NULL,
    na_as_unbound = TRUE,
    verbose       = TRUE
) {
  stopifnot(is.list(grn_set), "fp_bound" %in% names(grn_set))
  fb <- grn_set$fp_bound
  if (!"peak_ID" %in% names(fb))
    stop("`fp_bound` must contain a 'peak_ID' column.")

  # sample columns = all except key
  all_ids <- setdiff(names(fb), "peak_ID")
  if (is.null(samples)) {
    use_ids <- all_ids
  } else {
    miss <- setdiff(samples, all_ids)
    if (length(miss)) {
      stop("Some `samples` not found in fp_bound: ", paste(miss, collapse=", "))
    }
    use_ids <- samples
  }
  if (!length(use_ids)) stop("No sample columns selected to evaluate binding.")

  mat <- as.data.frame(fb[use_ids], stringsAsFactors = FALSE)

  # coerce to numeric 0/1
  mat[] <- lapply(mat, function(x) {
    if (is.logical(x)) as.integer(x)
    else if (is.numeric(x)) as.integer(x)
    else suppressWarnings(as.integer(x))
  })

  if (na_as_unbound) {
    mat[is.na(as.matrix(mat))] <- 0L
    keep <- rowSums(as.matrix(mat)) >= as.integer(min_bound)
  } else {
    # drop rows where any NA among selected samples
    any_na <- apply(as.matrix(mat), 1L, function(r) any(is.na(r)))
    keep   <- !any_na & rowSums(as.matrix(mat), na.rm = TRUE) >= as.integer(min_bound)
  }

  if (verbose) {
    message("Rows meeting ≥", min_bound, " bound samples among ",
            length(use_ids), " sample(s): ", sum(keep), "/", nrow(fb))
  }

  # delegate to the synced triplet filter to keep consistency
  filter_fp_rows(
    grn_set   = grn_set,
    peaks     = fb$peak_ID[keep],
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose   = verbose
  )
}


# Minimal update: fp_bound <- fp_bound & atac_overlap[ mapped atac_peak ]
update_fp_bound_by_overlap <- function(grn_set,
                                       key_fp   = "peak_ID",
                                       key_fpa  = "fp_peak",
                                       key_atac = "atac_peak",
                                       na_as_unbound = TRUE) {
  stopifnot(is.list(grn_set),
            all(c("fp_bound","fp_annotation","atac_overlap") %in% names(grn_set)))

  fb  <- grn_set$fp_bound
  fa  <- grn_set$fp_annotation
  aol <- grn_set$atac_overlap

  # sample IDs common to fp_bound and atac_overlap
  ids_fb  <- setdiff(names(fb),  key_fp)
  ids_ol  <- setdiff(names(aol), key_atac)
  ids_use <- intersect(ids_fb, ids_ol)
  if (!length(ids_use)) stop("No common sample IDs between fp_bound and atac_overlap.")

  # one-to-one map: fp_peak -> atac_peak (drop motif column)
  map <- dplyr::distinct(fa, !!rlang::sym(key_fpa), !!rlang::sym(key_atac))
  dup <- map |>
    dplyr::count(!!rlang::sym(key_fpa), name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup)) stop("Each fp_peak must map to a single atac_peak.")

  # per-fp_peak overlap flags, aligned to fp_bound row order
  ol_fp <- tibble::tibble(!!key_fpa := fb[[key_fp]]) |>
    dplyr::left_join(map, by = setNames(key_fpa, key_fpa)) |>
    dplyr::left_join(
      dplyr::select(aol, !!key_atac, dplyr::all_of(ids_use)),
      by = setNames(key_atac, key_atac)
    )

  # matrices for pointwise logical AND (treat NA as 0 by default)
  fb_mat  <- as.matrix(fb[ids_use])
  ol_mat  <- as.matrix(ol_fp[ids_use])
  storage.mode(fb_mat) <- "integer"
  storage.mode(ol_mat) <- "integer"
  if (na_as_unbound) {
    fb_mat[is.na(fb_mat)] <- 0L
    ol_mat[is.na(ol_mat)] <- 0L
  }

  new_mat <- (fb_mat > 0L) & (ol_mat > 0L)

  # rebuild fp_bound with identical shape/order; keep any extra (non-overlap) columns unchanged
  out <- dplyr::bind_cols(
    fb[key_fp],
    tibble::as_tibble(matrix(as.integer(new_mat), nrow(fb), length(ids_use),
                             dimnames = list(NULL, ids_use)))
  )
  if (length(extra <- setdiff(ids_fb, ids_use))) {
    out <- dplyr::bind_cols(out, fb[extra])
    out <- out[, c(key_fp, setdiff(names(fb), key_fp)), drop = FALSE]
  }

  # final guards
  if (!identical(fb[[key_fp]], out[[key_fp]])) stop("Row order changed unexpectedly.")
  if (!setequal(names(fb), names(out)))        stop("Column set changed unexpectedly.")

  grn_set$fp_bound <- out
  grn_set
}

# Subset RNA to TFs and create a per-sample 0/1 flag
add_tf_expr_flags <- function(grn_set, threshold = 0, gt = TRUE) {
  stopifnot(is.list(grn_set),
            all(c("rna","tf_list") %in% names(grn_set)))
  rna <- grn_set$rna
  ids <- setdiff(names(rna), c("ensembl_gene_id","HGNC"))
  if (!length(ids)) stop("No RNA sample columns found.")

  rna_tf <- dplyr::filter(rna, .data$HGNC %in% grn_set$tf_list)

  flag_mat <- as.matrix(rna_tf[ids])
  if (gt) {
    flag_mat <- (flag_mat >  threshold)
  } else {
    flag_mat <- (flag_mat >= threshold)
  }
  storage.mode(flag_mat) <- "integer"

  grn_set$rna_tf       <- rna_tf
  grn_set$tf_expr_flag <- dplyr::bind_cols(
    rna_tf[c("ensembl_gene_id","HGNC")],
    tibble::as_tibble(flag_mat, .name_repair = "minimal")
  )
  grn_set
}


# Gate fp_bound by TF expression, handling motif dimers (A::B) properly
# Requires grn_set$tf_expr_flag (from add_tf_expr_flags) and grn_set$motif_db
update_fp_bound_by_tf_expr <- function(grn_set, group_size = 1L) {
  stopifnot(is.list(grn_set),
            all(c("fp_bound","fp_annotation","tf_expr_flag","motif_db") %in% names(grn_set)))

  fb   <- grn_set$fp_bound
  fpa  <- grn_set$fp_annotation
  tff  <- grn_set$tf_expr_flag

  # shared sample IDs (preserve fp_bound order)
  ids_fb  <- setdiff(names(fb), "peak_ID")
  ids_tf  <- setdiff(names(tff), c("ensembl_gene_id","HGNC"))
  ids     <- unique(intersect(ids_fb, ids_tf))
  if (!length(ids)) stop("No common sample IDs between fp_bound and tf_expr_flag.")

  # --- expand motifs -> TFs (handle dimers A::B) and join TF flags -------------
  mdb_exp <- grn_set$motif_db |>
    tidyr::separate_rows("HGNC", sep = "::") |>
    dplyr::mutate(HGNC = stringr::str_trim(.data$HGNC)) |>
    dplyr::filter(.data$HGNC != "") |>
    dplyr::distinct(.data$motif, .data$HGNC)

  expanded <- fpa |>
    tidyr::separate_rows("motifs", sep = "\\s*,\\s*") |>
    dplyr::rename(motif = .data$motifs) |>
    dplyr::left_join(mdb_exp,           by = "motif") |>
    dplyr::left_join(dplyr::select(tff, "HGNC", dplyr::all_of(ids)), by = "HGNC")

  # --- per (fp_peak, motif-mapped TF): keep rows with >= group_size positives ---
  keep_row <- expanded |>
    dplyr::transmute(
      n_ok = rowSums(dplyr::across(dplyr::all_of(ids), ~ as.integer(.x == 1L)), na.rm = TRUE),
      ok   = .data$n_ok >= group_size
    ) |>
    dplyr::pull("ok")

  expanded$HGNC[!keep_row] <- NA_character_
  if (any(!keep_row)) expanded[!keep_row, ids] <- NA_integer_

  # --- collapse to one row per (fp_peak, atac_peak): ANY across TFs per sample ---
  data.table::setDT(expanded)
  collapsed_dt <- expanded[
    , c(
      as.list(vapply(.SD,
                     function(x) if (all(is.na(x))) NA_integer_
                     else as.integer(any(x == 1L, na.rm = TRUE)),
                     integer(1L))),
      list(tfs = { u <- unique(stats::na.omit(HGNC))
      if (length(u)) paste(u, collapse = ",") else NA_character_ })
    ),
    by = .(fp_peak, atac_peak),
    .SDcols = ids
  ]
  collapsed <- tibble::as_tibble(collapsed_dt)

  # peaks that still have at least one TF after gating
  keep_peaks <- collapsed |>
    dplyr::filter(!is.na(.data$tfs)) |>
    dplyr::pull("fp_peak")

  # --- AND with fp_bound (pointwise), preserving shape/order --------------------
  fb_keep <- dplyr::semi_join(fb, tibble::tibble(peak_ID = keep_peaks), by = "peak_ID")

  coll_aligned <- fb_keep["peak_ID"] |>
    dplyr::left_join(
      dplyr::select(collapsed, "fp_peak", dplyr::all_of(ids)),
      by = c("peak_ID" = "fp_peak")
    )

  fb_mat <- as.matrix(fb_keep[, ids, drop = FALSE]); storage.mode(fb_mat) <- "integer"
  tf_mat <- as.matrix(coll_aligned[, ids, drop = FALSE]); storage.mode(tf_mat) <- "integer"

  fb_new <- (fb_mat > 0L) & (tf_mat > 0L)                # logical matrix, same dims as fb_mat
  # SAFER rebuild: derive column names from the actual matrix we used
  out_mat <- matrix(as.integer(fb_new), nrow = nrow(fb_new), ncol = ncol(fb_new))
  colnames(out_mat) <- colnames(fb_mat)                  # <- avoids dimnames-length mismatch

  fb_out <- dplyr::bind_cols(
    fb_keep["peak_ID"],
    tibble::as_tibble(out_mat)
  )

  # add back any extra fp_bound columns (not part of ids), in original order
  extra <- setdiff(names(fb), c("peak_ID", ids))
  if (length(extra)) {
    fb_out <- dplyr::bind_cols(fb_out, fb_keep[extra]) |>
      dplyr::select("peak_ID", dplyr::all_of(setdiff(names(fb), "peak_ID")))
  }

  # --- update set (drop rows consistently; replace fp_bound; keep TF summary) ---
  grn_set <- filter_fp_rows(grn_set, peaks = keep_peaks)
  grn_set$fp_bound <- fb_out
  grn_set$fp_tfs   <- collapsed |>
    dplyr::filter(.data$fp_peak %in% keep_peaks) |>
    dplyr::select("fp_peak", "atac_peak", "tfs")

  grn_set
}

filter_fp_rows <- function(
    grn_set,
    peaks     = NULL,
    predicate = NULL,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose   = TRUE
) {
  stopifnot(is.list(grn_set),
            all(c("fp_score","fp_bound","fp_annotation") %in% names(grn_set)))

  fp_score <- grn_set$fp_score
  fp_bound <- grn_set$fp_bound
  fp_annot <- grn_set$fp_annotation

  s_ids <- fp_score[[score_key]]
  b_ids <- fp_bound[[bound_key]]
  a_ids <- fp_annot[[annot_key]]

  if (!is.character(s_ids) || !is.character(b_ids) || !is.character(a_ids)) {
    stop("Keys must be character; coerce before calling if needed.")
  }

  # --- NEW: if all three are already empty, just handle requests sanely and return
  if (length(s_ids) == 0L && length(b_ids) == 0L && length(a_ids) == 0L) {
    # Nothing to reconcile; respect peaks/predicate but result will remain empty.
    if (!is.null(predicate)) {
      keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
      if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
        stop("`predicate` must return a logical vector of length nrow(fp_score).")
      }
    }
    # If user passed non-empty peaks, warn that none can match.
    if (!is.null(peaks) && length(peaks) > 0L) {
      if (requireNamespace("cli", quietly = TRUE)) {
        cli::cli_warn("Requested {length(peaks)} peak ID(s), but the current set is empty; returning empty triplet.")
      } else {
        warning("Requested peaks in an empty set; returning empty triplet.")
      }
    }
    if (verbose) message("Filtered to 0 peaks (from 0).")
    grn_set$fp_score      <- fp_score[0, , drop = FALSE]
    grn_set$fp_bound      <- fp_bound[0, , drop = FALSE]
    grn_set$fp_annotation <- fp_annot[0, , drop = FALSE]
    return(grn_set)
  }

  # --- Reconcile key universe FIRST (align all three to common IDs; keep annot order)
  common_ids <- Reduce(intersect, list(unique(s_ids), unique(b_ids), unique(a_ids)))
  if (length(common_ids) == 0L) {
    stop("No common keys across fp_score/fp_bound/fp_annotation.")
  }
  if (!(setequal(s_ids, b_ids) && setequal(s_ids, a_ids))) {
    n_drop <- c(
      score_only = sum(!(s_ids %in% common_ids)),
      bound_only = sum(!(b_ids %in% common_ids)),
      annot_only = sum(!(a_ids %in% common_ids))
    )
    if (requireNamespace("cli", quietly = TRUE)) {
      cli::cli_warn(c(
        "!" = "Key sets differ across tables; reconciling to their intersection ({length(common_ids)} rows).",
        "i" = "Dropped: score-only={n_drop['score_only']}, bound-only={n_drop['bound_only']}, annot-only={n_drop['annot_only']}."
      ))
    } else {
      warning(sprintf(
        "Key sets differ; aligning to intersection (%d rows). Dropped score-only=%d, bound-only=%d, annot-only=%d.",
        length(common_ids), n_drop["score_only"], n_drop["bound_only"], n_drop["annot_only"]
      ))
    }

    # anchor order to annotation (deterministic)
    a_keep <- fp_annot[[annot_key]] %in% common_ids
    a_ord  <- fp_annot[[annot_key]][a_keep]

    s_keep <- fp_score[[score_key]] %in% common_ids
    b_keep <- fp_bound[[bound_key]] %in% common_ids

    fp_annot <- fp_annot[a_keep, , drop = FALSE]

    fp_score <- fp_score[s_keep, , drop = FALSE]
    idx_s    <- match(a_ord, fp_score[[score_key]])
    fp_score <- fp_score[idx_s, , drop = FALSE]

    fp_bound <- fp_bound[b_keep, , drop = FALSE]
    idx_b    <- match(a_ord, fp_bound[[bound_key]])
    fp_bound <- fp_bound[idx_b, , drop = FALSE]

    s_ids <- fp_score[[score_key]]
    b_ids <- fp_bound[[bound_key]]
    a_ids <- fp_annot[[annot_key]]

    stopifnot(identical(s_ids, a_ids), identical(b_ids, a_ids))
  }

  # --- Determine peaks to keep
  if (!is.null(predicate)) {
    keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
    if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
      stop("`predicate` must return a logical vector of length nrow(fp_score).")
    }
    keep_ids <- s_ids[keep_lgl %in% TRUE]
  } else if (!is.null(peaks)) {
    if (!all(peaks %in% s_ids)) {
      miss <- setdiff(peaks, s_ids)
      if (length(miss)) {
        if (requireNamespace("cli", quietly = TRUE)) {
          cli::cli_warn("Some requested `peaks` not present after reconciliation; dropping {length(miss)} missing id(s).")
        } else {
          warning(sprintf("Some requested `peaks` are missing after reconciliation; dropping %d id(s).", length(miss)))
        }
      }
    }
    peaks_uniq <- peaks[!duplicated(peaks)]
    keep_ids   <- peaks_uniq[peaks_uniq %in% s_ids]
  } else {
    if (verbose) message("No filter applied; returning original set.")
    grn_set$fp_score      <- fp_score
    grn_set$fp_bound      <- fp_bound
    grn_set$fp_annotation <- fp_annot
    return(grn_set)
  }

  if (!length(keep_ids)) {
    if (verbose) message("Filtered to 0 peaks (from ", length(s_ids), ").")
    fp_score2 <- fp_score[0, , drop = FALSE]
    fp_bound2 <- fp_bound[0, , drop = FALSE]
    fp_annot2 <- fp_annot[0, , drop = FALSE]
  } else {
    keep_s <- fp_score[[score_key]] %in% keep_ids
    keep_b <- fp_bound[[bound_key]] %in% keep_ids
    keep_a <- fp_annot[[annot_key]] %in% keep_ids

    fp_score2 <- fp_score[keep_s, , drop = FALSE]
    fp_bound2 <- fp_bound[keep_b, , drop = FALSE]
    fp_annot2 <- fp_annot[keep_a, , drop = FALSE]

    ord_a <- fp_annot2[[annot_key]]
    if (!is.null(peaks)) {
      target <- keep_ids
      ord_idx <- match(ord_a, target)
      o <- order(ord_idx, seq_along(ord_idx), na.last = TRUE)
      fp_annot2 <- fp_annot2[o, , drop = FALSE]
      m_s <- match(fp_annot2[[annot_key]], fp_score2[[score_key]])
      m_b <- match(fp_annot2[[annot_key]], fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    } else {
      m_s <- match(ord_a, fp_score2[[score_key]])
      m_b <- match(ord_a, fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    }
  }

  s2 <- fp_score2[[score_key]]
  b2 <- fp_bound2[[bound_key]]
  a2 <- fp_annot2[[annot_key]]
  if (!(identical(s2, a2) && identical(b2, a2))) {
    ks <- head(s2[!(s2 %in% a2)], 3)
    kb <- head(b2[!(b2 %in% a2)], 3)
    stop(paste0(
      "Post-check failed: filtered key sets are not identical across the three tables. ",
      "Examples — score-only: ", paste(ks, collapse = ", "),
      "; bound-only: ", paste(kb, collapse = ", ")
    ))
  }

  if (verbose) {
    message("Filtered to ", nrow(fp_annot2), " peaks (from ", length(s_ids), ").")
  }

  grn_set$fp_score       <- fp_score2
  grn_set$fp_bound       <- fp_bound2
  grn_set$fp_annotation  <- fp_annot2
  grn_set
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



#' Build ATAC score/overlap tibbles; sort by chr/start/end first (if enabled).
#'
#' @param atac_data data.frame/tibble with coords, metadata, numeric samples, and overlap cols.
#' @param chr_col,start_col,end_col coordinate column names (defaults: "Chr","Start","End").
#' @param overlap_prefix prefix for overlap columns (default "Overlap_").
#' @param metadata_cols metadata columns to exclude from samples (besides coords).
#' @param sort_peaks logical: sort rows by chr/start/end before processing (default TRUE).
#' @return list(score = atac_score, overlap = atac_overlap)
#' @export
load_atac <- function(
    atac_data,
    chr_col        = "Chr",
    start_col      = "Start",
    end_col        = "End",
    overlap_prefix = "Overlap_",
    metadata_cols  = c("PeakID","Distance.to.TSS","Nearest.PromoterID","Entrez.ID",
                       "Nearest.Unigene","Nearest.Refseq","Nearest.Ensembl",
                       "Gene.Name","Gene.Alias","Gene.Description","Gene.Type"),
    sort_peaks     = TRUE
) {
  if (!is.data.frame(atac_data)) cli::cli_abort("`atac_data` must be a data.frame/tibble.")
  need <- c(chr_col, start_col, end_col)
  if (!all(need %in% names(atac_data)))
    cli::cli_abort("Missing coord columns: {.val {setdiff(need, names(atac_data))}}.")

  df <- atac_data

  # ---- sort by chr/start/end first ----
  if (isTRUE(sort_peaks)) {
    chr <- tolower(gsub("^chr", "", as.character(df[[chr_col]])))
    chr_num <- suppressWarnings(as.integer(chr))
    chr_num[is.na(chr_num) & chr %in% c("x")]  <- 23L
    chr_num[is.na(chr_num) & chr %in% c("y")]  <- 24L
    chr_num[is.na(chr_num) & chr %in% c("m","mt","dmel_mito","mito")] <- 25L
    ord <- order(is.na(chr_num), chr_num, chr,
                 as.integer(round(df[[start_col]])),
                 as.integer(round(df[[end_col]])),
                 na.last = TRUE, method = "radix")
    df <- df[ord, , drop = FALSE]
  }

  # unified peak key
  atac_peak <- paste0(
    df[[chr_col]], ":", as.integer(round(df[[start_col]])),
    "-", as.integer(round(df[[end_col]]))
  )

  # detect columns
  overlap_cols <- grep(paste0("^", overlap_prefix), names(df), value = TRUE)
  meta_set     <- unique(c(metadata_cols, chr_col, start_col, end_col))
  candidates   <- setdiff(names(df), c(meta_set, overlap_cols))
  sample_cols  <- candidates[vapply(df[candidates], is.numeric, logical(1))]
  if (!length(sample_cols))
    cli::cli_abort("No numeric sample columns found after excluding metadata and {.val {overlap_prefix}}*.")

  # outputs
  atac_score <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, sample_cols, drop = FALSE])

  atac_overlap <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, overlap_cols, drop = FALSE])

  # strip prefix from overlap column names
  if (length(overlap_cols)) {
    new_names <- sub(paste0("^", overlap_prefix), "", names(atac_overlap))
    new_names[1L] <- "atac_peak"  # keep key intact
    # guard against accidental duplicates
    names(atac_overlap) <- make.unique(new_names, sep = "_")
  }

  list(score = dplyr::as_tibble(atac_score),
       overlap = dplyr::as_tibble(atac_overlap))
}

#' Clean HGNC symbols (drop any whitespace suffix), prints a compact diff.
#' @param df data.frame with column HGNC
#' @param label tag used in printed message
#' @return data.frame with cleaned HGNC column
clean_hgnc <- function(df, label = "rna") {
  df <- dplyr::mutate(df, HGNC = as.character(.data$HGNC))
  after   <- ifelse(is.na(df$HGNC), NA_character_,
                    stringr::str_replace(df$HGNC, "\\s+.*$", ""))
  changed <- which(!is.na(df$HGNC) & df$HGNC != after)
  if (length(changed)) {
    cat("=== HGNC cleanup:", label, "===\n")
    cat(paste(unique(paste0(df$HGNC[changed], " -> ", after[changed]))), sep = "\n")
  }
  df$HGNC <- after
  df
}


#' Filter RNA-seq by TF/non-TF thresholds with "at least N samples" logic.
#'
#' Keeps rows where >= `min_samples` expression columns meet the threshold:
#'   - non-TFs use `gene_min`
#'   - TFs     use `tf_min`
#'
#' @param rna data.frame/tibble: genes x samples (numeric sample columns).
#' @param tf_list character: TF symbols.
#' @param hgnc_col character: gene symbol column name (default "HGNC").
#' @param gene_min numeric: non-TF threshold (default 5).
#' @param tf_min numeric: TF threshold (default 10).
#' @param min_samples integer: require at least this many samples pass (default 1).
#' @param id_cols character: non-expression columns to exclude when testing.
#'   Defaults to c("ensembl_gene_id", `hgnc_col`).
#' @return tibble with filtered rows (original columns preserved).
filter_rna_expr <- function(
    rna,
    tf_list,
    hgnc_col   = "HGNC",
    gene_min   = 5,
    tf_min     = 10,
    min_samples = 1L,
    id_cols    = c("ensembl_gene_id", hgnc_col)
) {
  if (!is.data.frame(rna)) cli::cli_abort("`rna` must be a data.frame/tibble.")
  if (!is.character(tf_list)) cli::cli_abort("`tf_list` must be a character vector.")
  if (!hgnc_col %in% names(rna)) cli::cli_abort("Column {.val {hgnc_col}} not found in `rna`.")
  if (!is.numeric(min_samples) || length(min_samples) != 1 || min_samples < 1)
    cli::cli_abort("`min_samples` must be a positive integer.")

  id_cols   <- unique(id_cols)
  expr_cols <- setdiff(names(rna), id_cols)
  num_cols  <- expr_cols[vapply(rna[expr_cols], is.numeric, logical(1))]
  if (length(num_cols) == 0) cli::cli_abort("No numeric expression columns found after excluding ID columns.")

  # Build logical pass matrices and per-row counts
  M <- as.data.frame(rna[num_cols], stringsAsFactors = FALSE)
  pass_gene <- as.data.frame(lapply(M, function(x) as.numeric(x) >= gene_min), check.names = FALSE)
  pass_tf   <- as.data.frame(lapply(M, function(x) as.numeric(x) >= tf_min),   check.names = FALSE)
  cnt_gene  <- rowSums(pass_gene, na.rm = TRUE)
  cnt_tf    <- rowSums(pass_tf,   na.rm = TRUE)

  is_tf <- rna[[hgnc_col]] %in% tf_list
  keep  <- (!is_tf & cnt_gene >= min_samples) | (is_tf & cnt_tf >= min_samples)

  dplyr::as_tibble(rna[keep, , drop = FALSE])
}


# Correlate footprint score vs TF RNA (per motif)
annotate_fp_tf_corr_one_motif <- function(grn_set,
                                          cor_method = c("pearson","spearman"),
                                          min_non_na = 5L,
                                          tfs_col    = "tfs",
                                          digits     = 6,
                                          dedup_fun  = c("mean","median","first")) {
  cor_method <- match.arg(cor_method)
  dedup_fun  <- match.arg(dedup_fun)

  stopifnot(
    is.list(grn_set),
    all(c("fp_score","rna","fp_tfs","fp_annotation") %in% names(grn_set)),
    "peak_ID" %in% names(grn_set$fp_score),
    "HGNC"    %in% names(grn_set$rna),
    "fp_peak" %in% names(grn_set$fp_tfs),
    tfs_col %in% names(grn_set$fp_tfs) || (tfs_col == "tfs")
  )

  fp_tbl  <- grn_set$fp_score
  rna_tbl <- grn_set$rna
  ft_raw  <- grn_set$fp_tfs

  # collapse to one row per fp_peak with comma-separated TFs
  if (tfs_col == "tfs") {
    ft2 <- ft_raw[, c("fp_peak","tfs"), drop = FALSE]
  } else {
    ft2 <- ft_raw |>
      dplyr::filter(!is.na(.data[[tfs_col]]), .data[[tfs_col]] != "") |>
      dplyr::group_by(fp_peak) |>
      dplyr::summarise(tfs = paste(unique(.data[[tfs_col]]), collapse = ","), .groups = "drop")
  }

  # shared sample columns
  fp_ids  <- setdiff(names(fp_tbl), "peak_ID")
  rna_ids <- setdiff(names(rna_tbl), c("ensembl_gene_id","HGNC"))
  sample_cols <- intersect(fp_ids, rna_ids)
  if (length(sample_cols) < 2L) stop("Need >= 2 shared samples between fp_score and RNA.")

  # aggregate duplicate peak_IDs if needed
  X0 <- fp_tbl[, c("peak_ID", sample_cols), drop = FALSE]
  if (anyDuplicated(X0$peak_ID)) {
    aggfun <- switch(
      dedup_fun,
      mean   = function(z){ m <- mean(z, na.rm = TRUE); if (is.nan(m)) NA_real_ else m },
      median = function(z){ m <- stats::median(z, na.rm = TRUE); if (is.nan(m)) NA_real_ else m },
      first  = function(z){ i <- which(!is.na(z))[1]; if (length(i)) z[i] else NA_real_ }
    )
    Xuni <- X0 |>
      dplyr::group_by(peak_ID) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), aggfun), .groups = "drop")
  } else {
    Xuni <- X0
  }

  # convert to base df before rownames to avoid tibble warning
  Xuni_df <- as.data.frame(Xuni, check.names = FALSE)
  rownames(Xuni_df) <- Xuni_df$peak_ID
  Xmat <- as.matrix(Xuni_df[, sample_cols, drop = FALSE]); storage.mode(Xmat) <- "double"
  mXi  <- is.finite(Xmat)

  # list unique TF symbols
  tf_unique <- unique(unlist(strsplit(paste(ft2$tfs, collapse = ","), "\\s*,\\s*")))
  tf_unique <- tf_unique[nzchar(tf_unique)]
  if (!length(tf_unique)) {
    grn_set$fp_tfs <- dplyr::mutate(ft2, fp_tf_corr_r = NA_character_, fp_tf_corr_p = NA_character_)
    grn_set$fp_annotation <- grn_set$fp_annotation |>
      dplyr::left_join(dplyr::rename(grn_set$fp_tfs,
                                     corr_fp_tf_r = fp_tf_corr_r,
                                     corr_fp_tf_p = fp_tf_corr_p),
                       by = "fp_peak")
    grn_set$fp_tf_corr_diag <- tibble::tibble()
    return(grn_set)
  }

  # RNA matrix for those TFs
  Rdf <- as.data.frame(rna_tbl[rna_tbl$HGNC %in% tf_unique, c("HGNC", sample_cols), drop = FALSE],
                       check.names = FALSE)
  Rdf <- Rdf[!duplicated(Rdf$HGNC), , drop = FALSE]
  rownames(Rdf) <- Rdf$HGNC
  Rmat <- as.matrix(Rdf[, sample_cols, drop = FALSE]); storage.mode(Rmat) <- "double"

  # diagnostics
  diag_df <- tibble::tibble(
    HGNC      = tf_unique,
    present   = tf_unique %in% rownames(Rdf),
    n_finite  = NA_integer_,
    sd_expr   = NA_real_,
    n_ok_peaks= NA_integer_
  )

  # correlation maps
  r_map <- p_map <- setNames(vector("list", length(tf_unique)), tf_unique)

  for (i in seq_along(tf_unique)) {
    tf <- tf_unique[i]
    if (!tf %in% rownames(Rdf)) {
      r_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      p_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      diag_df$n_finite[i]  <- 0L
      diag_df$sd_expr[i]   <- NA_real_
      diag_df$n_ok_peaks[i] <- 0L
      next
    }
    y  <- Rmat[tf, , drop = TRUE]
    my <- is.finite(y)
    diag_df$n_finite[i] <- sum(my)
    diag_df$sd_expr[i]  <- stats::sd(y[my])

    if (sum(my) < 2L || isTRUE(diag_df$sd_expr[i] == 0)) {
      r_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      p_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      diag_df$n_ok_peaks[i] <- 0L
      next
    }

    r_v  <- suppressWarnings(as.vector(stats::cor(t(Xmat), y,
                                                  use = "pairwise.complete.obs",
                                                  method = cor_method)))
    n_vec <- as.integer(mXi %*% as.integer(my))
    ok    <- is.finite(r_v) & (n_vec >= min_non_na) & (n_vec > 2L)

    t_stat <- rep(NA_real_, length(r_v))
    t_stat[ok] <- r_v[ok] * sqrt((n_vec[ok]-2) / pmax(1e-12, 1 - r_v[ok]^2))
    p_v <- rep(NA_real_, length(r_v))
    p_v[ok] <- 2 * stats::pt(abs(t_stat[ok]), df = n_vec[ok]-2, lower.tail = FALSE)

    names(r_v) <- rownames(Xmat)
    names(p_v) <- rownames(Xmat)
    r_map[[tf]] <- r_v
    p_map[[tf]] <- p_v
    diag_df$n_ok_peaks[i] <- sum(ok, na.rm = TRUE)
  }

  # format helper
  fmt <- function(v) {
    out <- sprintf(paste0("%.", digits, "g"), v)
    out[is.na(v)] <- NA_character_
    out
  }

  # build comma-separated strings per fp_peak
  assemble <- function(peak, tf_str, which = c("r","p")) {
    which <- match.arg(which)
    if (is.na(tf_str) || !nzchar(tf_str)) return(NA_character_)
    tfs <- strsplit(tf_str, "\\s*,\\s*")[[1]]
    tfs <- tfs[nzchar(tfs)]
    if (!length(tfs)) return(NA_character_)
    vals <- vapply(
      tfs,
      function(tf) if (which == "r") r_map[[tf]][peak] else p_map[[tf]][peak],
      numeric(1)
    )
    paste(fmt(vals), collapse = ",")
  }

  corr_r <- mapply(assemble, ft2$fp_peak, ft2$tfs, MoreArgs = list(which = "r"), USE.NAMES = FALSE)
  corr_p <- mapply(assemble, ft2$fp_peak, ft2$tfs, MoreArgs = list(which = "p"), USE.NAMES = FALSE)

  # write back
  grn_set$fp_tfs <- ft2 |>
    dplyr::mutate(fp_tf_corr_r = as.character(corr_r),
                  fp_tf_corr_p = as.character(corr_p))

  grn_set$fp_annotation <- grn_set$fp_annotation |>
    dplyr::left_join(dplyr::rename(grn_set$fp_tfs,
                                   corr_fp_tf_r = fp_tf_corr_r,
                                   corr_fp_tf_p = fp_tf_corr_p),
                     by = "fp_peak")

  grn_set$fp_tf_corr_diag <- diag_df
  grn_set
}

# Helper: build fp_tfs if missing
.derive_fp_tfs <- function(motif_name, fp_annot, build_args) {
  # Expect fp_annot with fp_peak + motifs; build_args$motif_db has motif + HGNC
  if (!is.null(build_args$motif_db) &&
      is.data.frame(build_args$motif_db) &&
      all(c("motif","HGNC") %in% names(build_args$motif_db)) &&
      "motifs" %in% names(fp_annot)) {
    ft <- fp_annot |>
      dplyr::distinct(fp_peak, motifs) |>
      dplyr::left_join(build_args$motif_db |>
                         dplyr::select(motif, HGNC),
                       by = c("motifs" = "motif")) |>
      tidyr::separate_rows(HGNC, sep = "::", convert = FALSE) |>
      dplyr::filter(!is.na(HGNC), HGNC != "") |>
      dplyr::group_by(fp_peak) |>
      dplyr::summarise(tfs = paste(unique(HGNC), collapse = ","), .groups = "drop")
    if (nrow(ft)) return(ft)
  }

  # Fallback: guess TF symbol(s) from motif prefix
  guess <- toupper(sub("_.*$", "", motif_name))
  tibble::tibble(fp_peak = fp_annot$fp_peak,
                 tfs     = guess)
}

# Process ONE motif (load, annotate, write)
process_one_motif_corr <- function(entry,
                                   out_dir,
                                   build_args,
                                   cor_method    = "pearson",
                                   min_non_na    = 5L,
                                   skip_existing = TRUE,
                                   verbose       = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  stopifnot(all(need_cols %in% names(entry)))

  m       <- as.character(entry$motif)[1]
  f_score <- as.character(entry$score)[1]
  f_bound <- as.character(entry$bound)[1]
  f_annot <- as.character(entry$annot)[1]

  # Early graceful skip if any path is NA/missing
  if (is.na(f_score) || is.na(f_bound) || is.na(f_annot) ||
      !file.exists(f_score) || !file.exists(f_bound) || !file.exists(f_annot)) {
    if (verbose) message("Skipping motif ", m, " (missing input file paths).")
    return(tibble::tibble(motif = m,
                          n_rows = 0L,
                          score  = NA_character_,
                          bound  = NA_character_,
                          annot  = NA_character_))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_score <- file.path(out_dir, paste0(m, "_score.csv"))
  out_bound <- file.path(out_dir, paste0(m, "_bound.csv"))
  # NOTE: keep the annot output suffix consistent with your earlier pipeline.
  # You previously used "_annotation.csv" when building the filtered set.
  # If you want consistency, prefer "_annotation.csv" here too:
  out_annot <- file.path(out_dir, paste0(m, "_annotation.csv"))

  if (skip_existing && file.exists(out_score) && file.exists(out_bound) && file.exists(out_annot)) {
    if (verbose) message("Skipping ", m, " (already exists).")
    return(tibble::tibble(motif = m, n_rows = NA_integer_,
                          score = out_score, bound = out_bound, annot = out_annot))
  }

  if (verbose) message("Correlating motif ", m, " ...")

  fp_score <- readr::read_csv(f_score, show_col_types = FALSE)
  fp_bound <- readr::read_csv(f_bound, show_col_types = FALSE)
  fp_annot <- readr::read_csv(f_annot, show_col_types = FALSE)

  # Normalize common column names
  if (!"peak_ID" %in% names(fp_score) && "fp_peak" %in% names(fp_score)) {
    fp_score <- dplyr::rename(fp_score, peak_ID = "fp_peak")
  }
  if (!"peak_ID" %in% names(fp_bound) && "fp_peak" %in% names(fp_bound)) {
    fp_bound <- dplyr::rename(fp_bound, peak_ID = "fp_peak")
  }
  if (!all(c("fp_peak","atac_peak","motifs") %in% names(fp_annot))) {
    fp_annot <- fp_annot |>
      dplyr::rename_with(~"fp_peak",   dplyr::any_of("peak_ID")) |>
      dplyr::rename_with(~"atac_peak", dplyr::any_of("peak_ATAC")) |>
      dplyr::rename_with(~"motifs",    dplyr::any_of("TFBS_name"))
  }

  ss <- do.call(build_grn_set, c(list(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annot
  ), build_args))

  if (!("rna" %in% names(ss))) ss$rna <- build_args$rna

  if (!("fp_tfs" %in% names(ss)) || !NROW(ss$fp_tfs)) {
    ss$fp_tfs <- .derive_fp_tfs(motif_name = m, fp_annot = ss$fp_annotation, build_args = build_args)
  } else if (!"tfs" %in% names(ss$fp_tfs)) {
    if ("HGNC" %in% names(ss$fp_tfs)) {
      ss$fp_tfs <- ss$fp_tfs |>
        dplyr::filter(!is.na(HGNC), HGNC != "") |>
        dplyr::group_by(fp_peak) |>
        dplyr::summarise(tfs = paste(unique(HGNC), collapse = ","), .groups = "drop")
    } else {
      ss$fp_tfs <- .derive_fp_tfs(motif_name = m, fp_annot = ss$fp_annotation, build_args = build_args)
    }
  }

  have_all <- all(c("fp_score","fp_annotation","rna","fp_tfs") %in% names(ss)) &&
    is.data.frame(ss$fp_score) && is.data.frame(ss$fp_annotation) &&
    NROW(ss$rna) > 0L && NROW(ss$fp_tfs) > 0L
  if (have_all) {
    ss <- annotate_fp_tf_corr_one_motif(ss,
                                        cor_method = cor_method,
                                        min_non_na = min_non_na,
                                        tfs_col    = "tfs")
  } else if (verbose) {
    message("Skipping correlation for ", m, ": missing one of fp_score/fp_annotation/rna/fp_tfs.")
  }

  readr::write_csv(ss$fp_score,      out_score)
  readr::write_csv(ss$fp_bound,      out_bound)
  readr::write_csv(ss$fp_annotation, out_annot)

  tibble::tibble(motif = m,
                 n_rows = nrow(ss$fp_bound),
                 score  = out_score,
                 bound  = out_bound,
                 annot  = out_annot)
}


#' Parallel driver to correlate TF footprints with expression (manifest-aware)
#'
#' Runs motif-by-motif correlation and returns a manifest. When `skip_existing=TRUE`
#' and existing result CSVs are detected for a motif, this function now *backfills*
#' the `n_rows` column by counting rows from the existing `score` CSV so you don't
#' end up with `NA` in the output manifest.
#'
#' @param fp_manifest data.frame/tibble with columns: motif, score, bound, annot
#'        (paths to per-motif CSVs or inputs needed by `process_one_motif_corr`)
#' @param out_dir output directory for new results
#' @param build_args list forwarded to `process_one_motif_corr()`
#' @param motif_ids optional character vector to subset motifs
#' @param n_workers integer number of workers (default: detectCores)
#' @param set_plan logical; if TRUE set a future plan inside the function
#' @param skip_existing logical; if TRUE, skip recomputation when outputs exist
#' @param cor_method character; correlation method, e.g. "pearson"
#' @param min_non_na integer; minimum non-NA pairs for correlation
#' @param verbose logical; emit progress messages
#'
#' @return tibble with columns including motif, n_rows, score, bound, annot
#' @export
tf_corr_footprints <- function(fp_manifest,
                               out_dir,
                               build_args,
                               motif_ids     = NULL,
                               n_workers     = max(1L, parallel::detectCores(TRUE)),
                               set_plan      = TRUE,
                               skip_existing = TRUE,
                               cor_method    = "pearson",
                               min_non_na    = 5L,
                               verbose       = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  stopifnot(is.data.frame(fp_manifest), all(need_cols %in% names(fp_manifest)))

  rows <- fp_manifest
  if (!is.null(motif_ids)) {
    rows <- dplyr::semi_join(rows, tibble::tibble(motif = motif_ids), by = "motif")
    if (!nrow(rows)) stop("No matching motifs found in `fp_manifest` for the requested subset.")
  }

  # Keep only rows with non-NA paths AND files that exist
  rows <- rows |>
    dplyr::filter(!is.na(score), !is.na(bound), !is.na(annot)) |>
    dplyr::mutate(
      .ok = file.exists(score) & file.exists(bound) & file.exists(annot)
    ) |>
    dplyr::filter(.ok) |>
    dplyr::select(-.ok)

  if (!nrow(rows)) stop("After filtering, no motifs have available input files.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (verbose) message("Submitting ", nrow(rows), " motif(s) to ", n_workers, " worker(s) for correlation...")

  # Choose a plan if requested
  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !identical(Sys.getenv("RSTUDIO"), "1") &&
        requireNamespace("parallelly", quietly = TRUE) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore, workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
  }

  # Helper: fast row count for CSV (header counted by fread, we want nrow)
  fast_count_csv <- function(path) {
    # Prefer data.table::fread if available; otherwise fallback to readr
    if (requireNamespace("data.table", quietly = TRUE)) {
      dt <- tryCatch(
        data.table::fread(path, nThread = 1L, showProgress = FALSE),
        error = function(e) NULL
      )
      return(if (is.null(dt)) NA_integer_ else nrow(dt))
    } else if (requireNamespace("readr", quietly = TRUE)) {
      df <- tryCatch(
        readr::read_csv(path, progress = FALSE, show_col_types = FALSE),
        error = function(e) NULL
      )
      return(if (is.null(df)) NA_integer_ else nrow(df))
    } else {
      # Very small fallback using base; may be slower on large files
      con <- NULL
      out <- NA_integer_
      try({
        con <- file(path, open = "r")
        # Read & discard header
        header <- readLines(con, n = 1L, warn = FALSE)
        cnt <- 0L
        chunk <- character()
        repeat {
          chunk <- readLines(con, n = 100000L, warn = FALSE)
          if (!length(chunk)) break
          cnt <- cnt + length(chunk)
        }
        out <- cnt
      }, silent = TRUE)
      if (!is.null(con)) close(con)
      return(out)
    }
  }

  # Run per motif
  res <- future.apply::future_lapply(seq_len(nrow(rows)), function(i) {
    entry <- rows[i, , drop = FALSE]
    process_one_motif_corr(entry = entry,
                           out_dir = out_dir,
                           build_args = build_args,
                           cor_method = cor_method,
                           min_non_na = min_non_na,
                           skip_existing = skip_existing,
                           verbose = FALSE)
  }, future.seed = FALSE)

  out <- tibble::as_tibble(dplyr::bind_rows(res))

  # Backfill n_rows when skip_existing caused a bypass and left NA
  # Criteria: n_rows is NA, score path is non-NA and exists
  # We count rows from the score CSV as authoritative per-motif row count.
  if ("n_rows" %in% names(out)) {
    need_fill <- which(is.na(out$n_rows) & !is.na(out$score) & file.exists(out$score))
    if (length(need_fill)) {
      if (verbose) message("Backfilling n_rows for ", length(need_fill), " motif(s) from existing score CSVs...")
      filled <- vapply(need_fill, function(idx) fast_count_csv(out$score[[idx]]), integer(1))
      out$n_rows[need_fill] <- filled
    }
  } else {
    # If the downstream returns no n_rows column, create and fill it where possible
    out$n_rows <- NA_integer_
    need_fill <- which(!is.na(out$score) & file.exists(out$score))
    if (length(need_fill)) {
      if (verbose) message("Populating n_rows from existing score CSVs for ", length(need_fill), " motif(s)...")
      filled <- vapply(need_fill, function(idx) fast_count_csv(out$score[[idx]]), integer(1))
      out$n_rows[need_fill] <- filled
    }
    # Ensure conventional column order: motif, n_rows, then the paths
    front_cols <- c("motif", "n_rows")
    rest_cols  <- setdiff(names(out), front_cols)
    out <- out[, c(front_cols, rest_cols), drop = FALSE]
  }

  out
}


#' Correlate one TF with ALL footprints and export per-TF BEDs
#'
#' Minimal, fast path to:
#' 1) load a single TF expression vector,
#' 2) load *all* footprint score/annotation files from a manifest,
#' 3) compute row-wise correlations (per TFBS) vs TF expression,
#' 4) BH-adjust p-values,
#' 5) write 3 BED files (`*_all.bed`, `*_bound.bed`, `*_unbound.bed`) and one
#'    tab-delimited overview (`*_overview.txt`) for the TF.
#'
#' Assumptions:
#' - `fp_manifest` has columns `motif, n_rows, score, bound, annot`.
#' - Each `score` CSV has a `peak_ID` (or `fp_peak`) column and per-sample columns.
#' - Each `annot` CSV has `fp_peak` (or `peak_ID`), `atac_peak` (or `peak_ATAC`),
#'   and `motifs` (or `TFBS_name`) columns linking TFBS→ATAC peak and name.
#'
#' @param fp_manifest data.frame with columns motif,n_rows,score,bound,annot.
#' @param tf_name     Character scalar; label written to outputs (e.g., "HNF4A").
#' @param tf_expr     Named numeric vector (preferred) of TF expression with names
#'                    matching the footprint score column names. Alternatively, a
#'                    one-row numeric data.frame with column names as sample IDs.
#' @param out_dir     Output directory to write files.
#' @param cor_method  "pearson" (default) or "spearman".
#' @param min_non_na  Minimum paired non-NA values required to test (default 5).
#' @param r_thr       Correlation threshold for "bound" (default 0.30).
#' @param fdr_thr     BH-FDR threshold for "bound" (default 0.05).
#' @param threads     Integer threads for data.table fread (default: detectCores()).
#' @param verbose     Logical; progress messages (default TRUE).
#'
#' @return (Invisibly) a list with:
#' \describe{
#'   \item{annot}{Tibble of TFBS annotations with r, p, p_adj.}
#'   \item{files}{Named character vector of written file paths.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' tf_vec <- c(S1=10, S2=8, S3=12)  # names must match columns in score CSVs
#' res <- tf_corr_footprints_all_tfbs(
#'   fp_manifest = fp_aligned_normalized_filtered_manifest,
#'   tf_name  = "HNF4A",
#'   tf_expr  = tf_vec,
#'   out_dir  = file.path(tempdir(), "hnf4a_out")
#' )
#' res$files
#' }
tf_corr_footprints_all_tfbs <- function(
    fp_manifest,
    tf_name,
    tf_expr,
    out_dir,
    cor_method  = "pearson",
    min_non_na  = 5L,
    r_thr       = 0.30,
    fdr_thr     = 0.05,
    threads     = max(1L, parallel::detectCores(TRUE)),
    verbose     = TRUE
) {
  # ---- checks ----
  need_cols <- c("motif","n_rows","score","bound","annot")
  if (!is.data.frame(fp_manifest) || !all(need_cols %in% names(fp_manifest))) {
    cli::cli_abort("`fp_manifest` must be a data.frame with columns: {need_cols}.")
  }
  if (!is.character(tf_name) || length(tf_name) != 1L || !nzchar(tf_name)) {
    cli::cli_abort("`tf_name` must be a non-empty character scalar.")
  }
  # Coerce tf_expr
  if (is.data.frame(tf_expr)) {
    if (nrow(tf_expr) != 1L) cli::cli_abort("`tf_expr` data.frame must have exactly one row.")
    tf_expr <- as.numeric(tf_expr[1, , drop = TRUE])
    names(tf_expr) <- colnames(tf_expr)
  }
  if (!is.numeric(tf_expr) || is.null(names(tf_expr)) || !length(tf_expr)) {
    cli::cli_abort("`tf_expr` must be a named numeric vector (names = sample IDs).")
  }

  # keep rows with non-empty files & n_rows>0 (score+annot are essential)
  use <- fp_manifest[
    is.finite(fp_manifest$n_rows) & fp_manifest$n_rows > 0 &
      !is.na(fp_manifest$score) & file.exists(fp_manifest$score) &
      !is.na(fp_manifest$annot) & file.exists(fp_manifest$annot),
    , drop = FALSE
  ]
  if (!nrow(use)) {
    cli::cli_abort("No non-empty motifs with existing {`score`} and {`annot`} files.")
  }

  # ---- IO helpers ----
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fread_fast <- function(path) data.table::fread(path, nThread = as.integer(threads), showProgress = FALSE)

  # threads scope
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(threads))

  # ---- read & normalize annotations (fp_peak, atac_peak, motifs) ----
  norm_annot <- function(dt) {
    if (!"fp_peak"   %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID", "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"    %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    keep <- intersect(c("fp_peak","atac_peak","motifs"), names(dt))
    dt <- dt[, ..keep]
    # drop rows without fp_peak
    dt <- dt[!is.na(fp_peak) & nzchar(fp_peak)]
    dt
  }
  annot_dt <- data.table::rbindlist(lapply(use$annot, function(p) norm_annot(fread_fast(p))),
                                    use.names = TRUE, fill = TRUE)
  if (!nrow(annot_dt) || !"fp_peak" %in% names(annot_dt)) {
    cli::cli_abort("Annotation inputs lack required `fp_peak` column after normalization.")
  }
  # Deduplicate TFBS rows by first occurrence
  annot_dt <- annot_dt[!duplicated(fp_peak)]

  # ---- read & normalize scores (peak_ID + sample columns) ----
  norm_score <- function(dt) {
    if (!"peak_ID" %in% names(dt) && "fp_peak" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    if (!"peak_ID" %in% names(dt)) cli::cli_abort("Score file is missing `peak_ID` (or `fp_peak`) column.")
    # keep only peak_ID + intersecting sample columns (matching tf_expr names)
    sample_cols <- intersect(setdiff(names(dt), "peak_ID"), names(tf_expr))
    if (!length(sample_cols)) cli::cli_abort("No overlapping sample columns between score file and `tf_expr` names.")
    # coerce sample columns to numeric
    for (j in sample_cols) data.table::set(dt, j = j, value = suppressWarnings(as.numeric(dt[[j]])))
    dt[, c("peak_ID", sample_cols), with = FALSE]
  }
  score_dt <- data.table::rbindlist(lapply(use$score, function(p) norm_score(fread_fast(p))),
                                    use.names = TRUE, fill = TRUE)
  # de-dup by first occurrence and keep only peaks present in annotations
  score_dt <- score_dt[!duplicated(peak_ID)]
  score_dt <- score_dt[peak_ID %chin% annot_dt$fp_peak]

  if (!nrow(score_dt)) {
    cli::cli_abort("After alignment, no overlapping TFBS between scores and annotations.")
  }

  # align sample columns to tf_expr names & compute correlation per row
  sample_cols <- intersect(setdiff(names(score_dt), "peak_ID"), names(tf_expr))
  if (length(sample_cols) < min_non_na) {
    cli::cli_abort("Not enough overlapping samples (need >= {min_non_na}, got {length(sample_cols)}).")
  }
  # order tf_expr to score columns
  y <- suppressWarnings(as.numeric(tf_expr[sample_cols]))
  names(y) <- sample_cols

  X <- as.matrix(score_dt[, ..sample_cols])

  row_stat <- function(v) {
    ok <- is.finite(v) & is.finite(y)
    n  <- sum(ok)
    if (n < min_non_na) return(c(NA_real_, NA_real_, as.integer(n)))
    r <- suppressWarnings(stats::cor(v[ok], y[ok], method = cor_method))
    if (!is.finite(r)) return(c(NA_real_, NA_real_, as.integer(n)))
    tval <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
    pval <- 2 * stats::pt(-abs(tval), df = n - 2)
    c(r, pval, as.integer(n))
  }
  mat <- t(apply(X, 1L, row_stat))
  score_dt[, `:=`(
    corr_fp_tf_r = mat[, 1],
    corr_fp_tf_p = mat[, 2],
    pairs        = as.integer(mat[, 3])
  )]

  # BH across all finite p
  good_p <- is.finite(score_dt$corr_fp_tf_p)
  score_dt[, corr_fp_tf_p_adj := NA_real_]
  if (any(good_p)) {
    score_dt[good_p, corr_fp_tf_p_adj := stats::p.adjust(score_dt$corr_fp_tf_p[good_p], method = "BH")]
    score_dt[good_p, corr_fp_tf_p_adj := pmax(corr_fp_tf_p_adj, corr_fp_tf_p)]
  }

  # ---- join annotation + stats & build export tables ----
  merged <- score_dt[, .(fp_peak = peak_ID, corr_fp_tf_r, corr_fp_tf_p, corr_fp_tf_p_adj)][
    annot_dt, on = "fp_peak"
  ]
  merged[, TF := tf_name]

  # split coordinates
  mk_table <- function(dt) {
    if (!nrow(dt)) return(data.table::data.table())
    fp_s <- data.table::tstrsplit(dt$fp_peak, "[:-]", perl = TRUE)
    at_s <- data.table::tstrsplit(dt$atac_peak, "[:-]", perl = TRUE)
    data.table::data.table(
      TFBS_chr         = fp_s[[1]],
      TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
      TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
      TFBS_name        = dt$motifs,
      peak_chr         = at_s[[1]],
      peak_start       = suppressWarnings(as.integer(at_s[[2]])),
      peak_end         = suppressWarnings(as.integer(at_s[[3]])),
      TF               = dt$TF,
      corr_fp_tf_r     = dt$corr_fp_tf_r,
      corr_fp_tf_p     = dt$corr_fp_tf_p,
      corr_fp_tf_p_adj = dt$corr_fp_tf_p_adj
    )
  }
  all_dt <- mk_table(merged)

  # bound vs unbound
  is_bound <- is.finite(all_dt$corr_fp_tf_r) &
    is.finite(all_dt$corr_fp_tf_p_adj) &
    (all_dt$corr_fp_tf_r > r_thr) &
    (all_dt$corr_fp_tf_p_adj < fdr_thr)
  bound_dt   <- all_dt[is_bound]
  unbound_dt <- all_dt[!is_bound]

  # ---- write outputs ----
  .slug <- function(x) {
    out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    out <- gsub("^_+|_+$", "", out)
    if (!nzchar(out)) out <- "TF"
    out
  }
  tf_slug <- .slug(tf_name)
  f_all      <- file.path(out_dir, paste0(tf_slug, "_all.bed"))
  f_bound    <- file.path(out_dir, paste0(tf_slug, "_bound.bed"))
  f_unbound  <- file.path(out_dir, paste0(tf_slug, "_unbound.bed"))
  f_overview <- file.path(out_dir, paste0(tf_slug, "_overview.txt"))

  # BEDs: no header
  data.table::fwrite(all_dt,     f_all,     sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")
  data.table::fwrite(bound_dt,   f_bound,   sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")
  data.table::fwrite(unbound_dt, f_unbound, sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")

  # overview: header + _bound flag
  if (nrow(all_dt)) {
    all_dt[, `_bound` := as.integer(is_bound)]
    data.table::fwrite(all_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE, na = "NA")
  } else {
    hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
             "peak_chr","peak_start","peak_end","TF",
             "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
    empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode="character", length=0)), length(hdr)), hdr))
    data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
  }

  if (verbose) {
    cli::cli_inform(c(
      "TF: {tf_name}",
      "All TFBS: {nrow(all_dt)}; Bound (r>{r_thr}, FDR<{fdr_thr}): {nrow(bound_dt)}; Unbound: {nrow(unbound_dt)}",
      "Wrote:",
      "  - {f_all}",
      "  - {f_bound}",
      "  - {f_unbound}",
      "  - {f_overview}"
    ))
  }

  invisible(list(
    annot = tibble::as_tibble(merged),
    files = c(all = f_all, bound = f_bound, unbound = f_unbound, overview = f_overview)
  ))
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
#' @export
#' @examples
#' # res <- tf_corr_footprints_filter(fp_corr_manifest)
tf_corr_footprints_filter <- function(
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

