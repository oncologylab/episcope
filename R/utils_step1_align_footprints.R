#' Align footprints by peak similarity
#'
#' Aligns footprint peaks across motifs and samples by grouping peaks that match
#' on rounded score signatures within each ATAC peak and chromosome, then merges
#' nearby peaks by midpoint proximity.
#'
#' @param fp_filtered_manifest Tibble/data.frame with columns
#'   \code{motif, score, bound, annot} and either \code{n_rows} or \code{n_peaks}.
#' @param mid_slop Integer; midpoint tolerance in bp for merging peaks.
#' @param round_digits Integer; number of digits for rounding scores before matching.
#' @param score_match_pct Numeric; fraction of score columns that must match.
#' @param verbose Logical; print progress.
#' @param threads Integer; data.table threads and parallel workers for fuzzy matching.
#' @param cache_dir Optional directory for cached outputs.
#' @param cache_tag Optional cache tag (used to build cache file names).
#' @param use_cache Logical; read cached outputs if available.
#' @param write_cache Logical; write cached outputs.
#' @param log_file Optional log file path for per-chromosome summaries.
#' @param parallel_by Parallelization scope. Use \code{"atac_peak"} to parallelize
#'   within each chromosome (default), \code{"chromosome"} to parallelize across
#'   chromosomes, or \code{"none"} for sequential execution.
#' @param output_mode Output structure. Use \code{"full"} (default) or
#'   \code{"distinct"}; both return the 4-slot list with aligned, unique rows.
#'   In \code{"distinct"}, \code{fp_score} has one row per aligned \code{peak_ID}
#'   (no motif duplication) and \code{fp_bound} matches that row set.
#' @param save_prealign_score Logical; if TRUE, write a merged footprint score
#'   table immediately after loading motif score files and before any alignment.
#' @param prealign_score_path Optional CSV path for the pre-alignment score
#'   output. If NULL and \code{cache_dir} is set, defaults to
#'   \code{fp_scores_prealign_<cache_tag>.csv} in \code{cache_dir}.
#' @param prealign_output_mode Output mode for pre-alignment score export:
#'   \code{"distinct"} (default; one row per \code{peak_ID}) or \code{"all"}
#'   (all rows as loaded from motif files).
#' @param prealign_only Logical; if TRUE and \code{save_prealign_score = TRUE},
#'   return immediately after writing the pre-alignment score table (skip alignment).
#'
#' @return List with tibbles: \code{fp_score}, \code{fp_bound},
#'   \code{fp_annotation}, and \code{id_map}.
#' @export
#'
#' @examples
#' \dontrun{
#' fp_aligned <- align_footprints(fp_manifest, cache_dir = "cache", cache_tag = "JASPAR2024")
#' }
align_footprints <- function(
    fp_filtered_manifest,
    mid_slop         = 10L,    # midpoint tolerance (bp)
    round_digits     = 1L,     # round scores before comparing vectors
    score_match_pct  = 0.8,    # fraction of samples that must match (<=1, >=0)
    verbose          = TRUE,
    threads          = max(1L, parallel::detectCores(TRUE)),
    cache_dir        = NULL,
    cache_tag        = NULL,
    use_cache        = TRUE,
    write_cache      = TRUE,
    log_file         = NULL,
    parallel_by      = c("atac_peak", "chromosome", "none"),
    output_mode      = c("full", "distinct"),
    save_prealign_score = FALSE,
    prealign_score_path = NULL,
    prealign_output_mode = c("distinct", "all"),
    prealign_only = FALSE
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
  prealign_output_mode <- match.arg(prealign_output_mode)

  prealign_path_use <- prealign_score_path
  if (is.null(prealign_path_use) && !is.null(cache_dir)) {
    tag <- if (!is.null(cache_tag) && nzchar(cache_tag)) cache_tag else "cache"
    prealign_path_use <- file.path(cache_dir, sprintf("fp_scores_prealign_%s.csv", tag))
  }

  cache_paths <- NULL
  if (!is.null(cache_dir) && !is.null(cache_tag)) {
    cache_paths <- list(
      fp_bound = file.path(cache_dir, sprintf("fp_bounds_%s.csv", cache_tag)),
      fp_score = file.path(cache_dir, sprintf("fp_scores_%s.csv", cache_tag)),
      fp_annotation = file.path(cache_dir, sprintf("fp_annotation_%s.csv", cache_tag)),
      id_map = file.path(cache_dir, sprintf("fp_id_map_%s.csv", cache_tag))
    )
  }

  if (!is.null(cache_paths) && isTRUE(use_cache) &&
      all(file.exists(unlist(cache_paths[c("fp_bound", "fp_score", "fp_annotation")])))) {
    if (verbose) .log_inform("Using cached aligned footprints from {.path {cache_dir}} (tag = {cache_tag}).")
    fp_bound <- readr::read_csv(cache_paths$fp_bound, show_col_types = FALSE)
    fp_score <- readr::read_csv(cache_paths$fp_score, show_col_types = FALSE)
    fp_annotation <- readr::read_csv(cache_paths$fp_annotation, show_col_types = FALSE)
    id_map <- if (file.exists(cache_paths$id_map)) {
      readr::read_csv(cache_paths$id_map, show_col_types = FALSE)
    } else if (all(c("fp_peak", "atac_peak") %in% names(fp_annotation))) {
      tibble::tibble(
        peak_ID = fp_annotation$fp_peak,
        fp_peak_bak = fp_annotation$fp_peak,
        atac_peak = fp_annotation$atac_peak,
        group_size = NA_real_
      )
    } else {
      tibble::tibble()
    }
    out_cached <- list(
      fp_score = tibble::as_tibble(fp_score),
      fp_bound = tibble::as_tibble(fp_bound),
      fp_annotation = tibble::as_tibble(fp_annotation),
      id_map = tibble::as_tibble(id_map)
    )
    if (output_mode == "distinct") {
      out_cached$fp_score <- out_cached$fp_score[!duplicated(out_cached$fp_score$peak_ID), , drop = FALSE]
      out_cached$fp_bound <- out_cached$fp_bound[!duplicated(out_cached$fp_bound$peak_ID), , drop = FALSE]
    }
    if (!identical(out_cached$fp_score$peak_ID, out_cached$fp_bound$peak_ID)) {
      cli::cli_abort("fp_score and fp_bound peak_ID rows are not in the same order.")
    }
    if (isTRUE(save_prealign_score)) {
      if (!is.null(prealign_path_use) && file.exists(prealign_path_use)) {
        if (isTRUE(verbose)) {
          .log_inform("Using existing pre-alignment score table: {.path {prealign_path_use}}.")
        }
      } else {
        .log_warn("`save_prealign_score = TRUE` but aligned cache was used; set `use_cache = FALSE` to regenerate pre-alignment score output.")
      }
    }
    return(out_cached)
  }

  # Use base filtering to avoid dplyr data mask quirks
  use <- fp_filtered_manifest[
    !is.na(fp_filtered_manifest[[count_col]]) & fp_filtered_manifest[[count_col]] > 0 &
      !is.na(fp_filtered_manifest$score) & file.exists(fp_filtered_manifest$score) &
      !is.na(fp_filtered_manifest$bound) & file.exists(fp_filtered_manifest$bound) &
      !is.na(fp_filtered_manifest$annot) & file.exists(fp_filtered_manifest$annot),
    , drop = FALSE
  ]

  if (!nrow(use)) {
    if (verbose) .log_inform("No non-empty motifs with existing files.")
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

  if (verbose) .log_inform("Loading & combining ({nrow(use)} motifs)...")
  score_dt <- data.table::rbindlist(lapply(use$score,  read_score), use.names = TRUE, fill = TRUE)
  bound_dt <- data.table::rbindlist(lapply(use$bound,  read_bound), use.names = TRUE, fill = TRUE)
  annot_dt <- data.table::rbindlist(lapply(use$annot,  read_annot), use.names = TRUE, fill = TRUE)
  # Keep annot_dt EXACT as loaded (row count + order).

  if (isTRUE(save_prealign_score)) {
    if (is.null(prealign_path_use) || !nzchar(prealign_path_use)) {
      .log_warn("`save_prealign_score = TRUE` but no writable output path was resolved; skipping pre-alignment export.")
    } else {
      prealign_out <- data.table::copy(score_dt)
      if (identical(prealign_output_mode, "distinct") && "peak_ID" %in% names(prealign_out)) {
        prealign_out <- prealign_out[!duplicated(peak_ID)]
      }
      dir.create(dirname(prealign_path_use), recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(tibble::as_tibble(prealign_out), prealign_path_use)
      if (isTRUE(verbose)) {
        .log_inform("Saved pre-alignment score table ({prealign_output_mode}): {.path {prealign_path_use}}.")
      }
      if (isTRUE(prealign_only)) {
        if (isTRUE(verbose)) {
          .log_inform("`prealign_only = TRUE`: skipping alignment and returning pre-alignment table only.")
        }
        return(list(
          fp_score = tibble::as_tibble(prealign_out),
          fp_bound = tibble::tibble(),
          fp_annotation = tibble::tibble(),
          id_map = tibble::tibble()
        ))
      }
    }
  }

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

  # ---- signature of rounded scores (per chromosome) ----
  sample_cols <- setdiff(names(score_dt), c("peak_ID","atac_peak","chr","start","end","mid"))
  if (!length(sample_cols)) stop("No sample columns detected in fp_score.")
  mult <- 10^round_digits
  sig_id_vec <- integer(nrow(score_dt))
  chr_groups <- split(seq_len(nrow(score_dt)), score_dt$chr, drop = TRUE)

  log_file_use <- log_file
  if (is.null(log_file_use) && !is.null(cache_dir)) {
    tag <- if (!is.null(cache_tag) && nzchar(cache_tag)) cache_tag else "cache"
    log_file_use <- file.path(cache_dir, sprintf("align_footprints_%s_log.txt", tag))
  }
  if (!is.null(log_file_use)) {
    dir.create(dirname(log_file_use), recursive = TRUE, showWarnings = FALSE)
  }
  log_line <- function(msg) {
    if (is.null(log_file_use)) return(invisible(NULL))
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write(paste(stamp, msg, sep = "\t"), file = log_file_use, append = TRUE)
    invisible(NULL)
  }

  .eq_counts <- function(A, v) {
    B <- matrix(v, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
    E <- (A == B) | (is.na(A) & is.na(B))
    rowSums(E)
  }
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

  parallel_by <- match.arg(parallel_by)
  output_mode <- match.arg(output_mode)
  use_parallel <- .Platform$OS.type == "unix" && threads > 1L
  use_inner_parallel <- use_parallel && parallel_by == "atac_peak"
  if (verbose && !is.null(log_file_use)) {
    .log_inform("Per-chromosome logs will be written to {.path {log_file_use}}.")
  }

  process_chr <- function(chr, idx_chr) {
    if (!length(idx_chr)) {
      return(list(idx_chr = idx_chr, sig_local = integer(0)))
    }

    S_sig <- as.matrix(score_dt[idx_chr, ..sample_cols])
    S_sig <- round(S_sig * mult)
    storage.mode(S_sig) <- "integer"

    atac_chr <- score_dt$atac_peak[idx_chr]
    splits <- split(seq_len(length(idx_chr)), atac_chr, drop = TRUE)
    log_line(sprintf("chr=%s n_rows=%d n_atac_peaks=%d", chr, length(idx_chr), length(splits)))

    if (score_match_pct >= 0.9999) {
      tmp_dt <- data.table::as.data.table(S_sig)
      tmp_dt[, atac_peak := atac_chr]
      sig_cols <- names(tmp_dt)[seq_len(ncol(tmp_dt) - 1L)]
      tmp_dt[, sig_id := data.table::rleidv(.SD), by = .(atac_peak), .SDcols = sig_cols]
      sig_local <- tmp_dt$sig_id
    } else {
      k_req <- ceiling(score_match_pct * ncol(S_sig))
      ord <- order(vapply(splits, length, integer(1L)), decreasing = TRUE)
      splits <- splits[ord]

      if (use_inner_parallel) {
        n_workers <- max(1L, min(threads, length(splits)))
        if (length(splits) == 1L) {
          comp_list <- list(assign_components(S_sig[splits[[1L]], , drop = FALSE], k_req))
        } else {
          batch_ids <- parallel::splitIndices(length(splits), n_workers)
          comp_batches <- parallel::mclapply(
            batch_ids,
            function(idx) {
              lapply(idx, function(i) assign_components(S_sig[splits[[i]], , drop = FALSE], k_req))
            },
            mc.cores = n_workers,
            mc.preschedule = FALSE
          )
          comp_list <- vector("list", length(splits))
          for (b in seq_along(batch_ids)) {
            idx <- batch_ids[[b]]
            vals <- comp_batches[[b]]
            for (j in seq_along(idx)) comp_list[[idx[[j]]]] <- vals[[j]]
          }
        }
      } else {
        comp_list <- lapply(
          splits,
          function(ii) assign_components(S_sig[ii, , drop = FALSE], k_req)
        )
      }

      sig_local <- integer(length(idx_chr))
      next_id <- 1L
      for (g in seq_along(splits)) {
        ids <- comp_list[[g]]
        sig_local[splits[[g]]] <- ids + (next_id - 1L)
        next_id <- next_id + max(ids, 0L)
      }
    }
    list(idx_chr = idx_chr, sig_local = sig_local)
  }

  if (parallel_by == "chromosome" && use_parallel) {
    chr_names <- names(chr_groups)
    chr_results <- parallel::mclapply(
      chr_names,
      function(chr) process_chr(chr, chr_groups[[chr]]),
      mc.cores = max(1L, min(threads, length(chr_names))),
      mc.preschedule = FALSE
    )
    for (res in chr_results) {
      if (length(res$idx_chr)) sig_id_vec[res$idx_chr] <- res$sig_local
    }
  } else {
    for (chr in names(chr_groups)) {
      res <- process_chr(chr, chr_groups[[chr]])
      if (length(res$idx_chr)) sig_id_vec[res$idx_chr] <- res$sig_local
    }
  }

  score_dt[, sig_id := sig_id_vec]

  if (verbose && !is.null(log_file_use)) {
    .log_inform("Alignment log saved to {.path {log_file_use}}.")
  }

  # ---- connected components on mid-slop within (atac_peak, chr, sig_id) ----
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

  out <- list(
    fp_score      = tibble::as_tibble(fp_score_out),
    fp_bound      = tibble::as_tibble(fp_bound_out),
    fp_annotation = tibble::as_tibble(annot_dt),
    id_map        = tibble::as_tibble(id_map_aligned)
  )

  if (output_mode == "distinct") {
    out$fp_score <- out$fp_score[!duplicated(out$fp_score$peak_ID), , drop = FALSE]
    out$fp_bound <- out$fp_bound[!duplicated(out$fp_bound$peak_ID), , drop = FALSE]
  }
  if (!identical(out$fp_score$peak_ID, out$fp_bound$peak_ID)) {
    cli::cli_abort("fp_score and fp_bound peak_ID rows are not in the same order.")
  }

  if (!is.null(cache_paths) && isTRUE(write_cache)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(out$fp_bound, cache_paths$fp_bound)
    readr::write_csv(out$fp_score, cache_paths$fp_score)
    readr::write_csv(out$fp_annotation, cache_paths$fp_annotation)
    readr::write_csv(out$id_map, cache_paths$id_map)
  }

  out
}
