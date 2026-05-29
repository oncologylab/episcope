# File: utils_step1_align_footprints.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 footprint alignment helpers.
#
# Inputs:
# - raw footprint score, bound, and annotation tables
# - alignment tolerances and cache settings
#
# Outputs:
# - aligned and optionally cached footprint objects for downstream Step 1 use
#
# Notes:
# - Keep alignment logic isolated from correlation and plotting code.

#' Module 1 footprint alignment helpers
#'
#' @noRd
NULL

.fp_eq_counts <- function(a, v) {
  b <- matrix(v, nrow = nrow(a), ncol = ncol(a), byrow = TRUE)
  e <- (a == b) | (is.na(a) & is.na(b))
  rowSums(e)
}

.assign_fp_score_components_base <- function(m, k_req) {
  n <- nrow(m)
  if (n == 0L) return(integer(0))
  ids <- integer(n)
  comp <- 0L
  for (i in seq_len(n)) {
    if (ids[i] != 0L) next
    comp <- comp + 1L
    ids[i] <- comp
    q <- i
    while (length(q)) {
      a <- q[1L]
      q <- q[-1L]
      unassigned <- which(ids == 0L)
      if (!length(unassigned)) break
      eqc <- .fp_eq_counts(m[unassigned, , drop = FALSE], m[a, ])
      to_add <- unassigned[eqc >= k_req]
      if (length(to_add)) {
        ids[to_add] <- comp
        q <- c(q, to_add)
      }
    }
  }
  ids
}

.assign_fp_score_components <- function(m,
                                        k_req,
                                        compress_duplicates = TRUE,
                                        compress_min_n = 64L) {
  n <- nrow(m)
  if (n == 0L) return(integer(0))
  if (!isTRUE(compress_duplicates) || n < as.integer(compress_min_n)) {
    return(.assign_fp_score_components_base(m, k_req = k_req))
  }
  keys <- do.call(paste, c(as.data.frame(m, optional = TRUE), sep = "\r"))
  unique_pos <- !duplicated(keys)
  if (all(unique_pos)) {
    return(.assign_fp_score_components_base(m, k_req = k_req))
  }
  unique_keys <- keys[unique_pos]
  key_map <- match(keys, unique_keys)
  ids_unique <- .assign_fp_score_components_base(m[unique_pos, , drop = FALSE], k_req = k_req)
  ids_unique[key_map]
}

# Alignment cache writes compact score/bound matrices plus fp_sites.
# fp_sites preserves all motif labels and original pre-alignment FP peaks per aligned FP.
align_footprints <- function(
    fp_filtered_manifest,
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 0.8,
    verbose = TRUE,
    threads = .available_cores(logical = TRUE),
    cache_dir = NULL,
    cache_tag = NULL,
    use_cache = TRUE,
    write_cache = TRUE,
    log_file = NULL,
    parallel_by = c("none", "chromosome", "atac_peak"),
    output_mode = c("full", "distinct"),
    return_data = TRUE,
    save_prealign_score = FALSE,
    prealign_score_path = NULL,
    prealign_output_mode = c("distinct", "all"),
    prealign_only = FALSE,
    return_id_map = TRUE,
    write_id_map = return_id_map,
    write_fp_sites = TRUE,
    include_source_fp_peaks = TRUE,
    cache_format = c("csv", "parquet", "both")
) {
  stopifnot(
    is.data.frame(fp_filtered_manifest),
    all(c("motif", "score", "bound", "annot") %in% names(fp_filtered_manifest))
  )
  count_col <- if ("n_rows" %in% names(fp_filtered_manifest)) {
    "n_rows"
  } else if ("n_peaks" %in% names(fp_filtered_manifest)) {
    "n_peaks"
  } else {
    .log_abort("Manifest must contain either 'n_rows' or 'n_peaks'.")
  }

  stopifnot(is.numeric(score_match_pct), length(score_match_pct) == 1L, score_match_pct > 0, score_match_pct <= 1)
  prealign_output_mode <- match.arg(prealign_output_mode)
  cache_format <- match.arg(cache_format)
  parallel_by <- match.arg(parallel_by)
  output_mode <- match.arg(output_mode)
  stopifnot(is.logical(return_data), length(return_data) == 1L, !is.na(return_data))
  stopifnot(is.logical(write_fp_sites), length(write_fp_sites) == 1L, !is.na(write_fp_sites))
  stopifnot(is.logical(include_source_fp_peaks), length(include_source_fp_peaks) == 1L, !is.na(include_source_fp_peaks))
  cache_format_primary <- if (identical(cache_format, "both")) "csv" else cache_format

  cache_summary <- function(paths, counts, format) {
    list(
      paths = paths,
      counts = counts,
      format = format,
      output_mode = output_mode,
      cache_schema = "compact",
      id_map_written = isTRUE(write_id_map),
      fp_sites_written = isTRUE(write_fp_sites)
    )
  }

  cache_counts_from_paths <- function(paths, format) {
    if (!identical(format, "csv")) {
      return(list(fp_score = NA_integer_, fp_bound = NA_integer_, fp_annotation = NA_integer_, id_map = NA_integer_, fp_sites = NA_integer_))
    }
    count_path <- function(path, missing = NA_integer_) {
      if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
        return(missing)
      }
      .count_rows_fast(path)
    }
    list(
      fp_score = count_path(paths$fp_score),
      fp_bound = count_path(paths$fp_bound),
      fp_annotation = count_path(paths$fp_annotation, missing = 0L),
      id_map = count_path(paths$id_map, missing = 0L),
      fp_sites = count_path(paths$fp_sites, missing = 0L)
    )
  }

  prealign_path_use <- prealign_score_path
  if (is.null(prealign_path_use) && !is.null(cache_dir)) {
    tag <- if (!is.null(cache_tag) && nzchar(cache_tag)) cache_tag else "cache"
    prealign_path_use <- file.path(cache_dir, sprintf("fp_scores_prealign_%s.csv", tag))
  }

  cache_paths <- NULL
  if (!is.null(cache_dir) && !is.null(cache_tag)) {
    cache_paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, "csv")
  }

  if (!is.null(cache_paths) && isTRUE(use_cache)) {
    cache_read_format <- if (identical(cache_format, "csv")) "csv" else "auto"
    cache_info <- .aligned_fp_cache_choose_format(cache_dir, cache_tag, cache_read_format)
    if (!isTRUE(cache_info$exists)) {
      cache_info <- .aligned_fp_cache_choose_format(cache_dir, cache_tag, "csv")
    }
    if (isTRUE(cache_info$exists)) {
      if (!isTRUE(return_data)) {
        paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, cache_info$format)
        return(cache_summary(paths, cache_counts_from_paths(paths, cache_info$format), cache_info$format))
      }
      return(load_fp_aligned_from_cache(
        cache_dir = cache_dir,
        cache_tag = cache_tag,
        output_mode = output_mode,
        load_id_map = return_id_map,
        cache_format = cache_info$format,
        verbose = verbose
      ))
    }
  }


  use <- fp_filtered_manifest[
    !is.na(fp_filtered_manifest[[count_col]]) & fp_filtered_manifest[[count_col]] > 0 &
      !is.na(fp_filtered_manifest$score) & file.exists(fp_filtered_manifest$score) &
      !is.na(fp_filtered_manifest$bound) & file.exists(fp_filtered_manifest$bound) &
      !is.na(fp_filtered_manifest$annot) & file.exists(fp_filtered_manifest$annot),
    , drop = FALSE
  ]

  if (!nrow(use)) {
    if (isTRUE(verbose)) .log_inform("No non-empty motifs with existing files.")
    if (!isTRUE(return_data)) {
      paths <- if (!is.null(cache_dir) && !is.null(cache_tag)) .aligned_fp_cache_paths(cache_dir, cache_tag, cache_format_primary) else list()
      return(cache_summary(paths, list(fp_score = 0L, fp_bound = 0L, fp_annotation = 0L, id_map = 0L, fp_sites = 0L), cache_format_primary))
    }
    return(list(fp_score = tibble::tibble(), fp_bound = tibble::tibble(), fp_annotation = tibble::tibble(), fp_sites = tibble::tibble(), id_map = tibble::tibble(), cache_schema = "compact"))
  }

  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(threads)

  fread_fast <- function(p) data.table::fread(p, nThread = threads, showProgress = FALSE)
  read_score <- function(p) {
    dt <- fread_fast(p)
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    ids <- setdiff(names(dt), "peak_ID")
    for (j in ids) {
      if (!is.numeric(dt[[j]])) data.table::set(dt, j = j, value = as.numeric(dt[[j]]))
    }
    dt[]
  }
  read_bound <- function(p) {
    dt <- fread_fast(p)
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    ids <- setdiff(names(dt), "peak_ID")
    for (j in ids) {
      if (!is.integer(dt[[j]])) data.table::set(dt, j = j, value = as.integer(dt[[j]]))
    }
    dt[]
  }
  read_annot <- function(p) {
    dt <- fread_fast(p)
    if (!"fp_peak" %in% names(dt) && "peak_ID" %in% names(dt)) data.table::setnames(dt, "peak_ID", "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs" %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    dt[]
  }

  if (isTRUE(verbose)) {
    .log_inform("Loading and combining {nrow(use)} motifs.")
    .log_inform("Alignment threads: {threads}; parallelization: {match.arg(parallel_by)}.")
  }
  score_dt <- data.table::rbindlist(lapply(use$score, read_score), use.names = TRUE, fill = TRUE)
  annot_dt <- data.table::rbindlist(lapply(use$annot, read_annot), use.names = TRUE, fill = TRUE)

  if (isTRUE(save_prealign_score) && !is.null(prealign_path_use) && nzchar(prealign_path_use)) {
    prealign_out <- data.table::copy(score_dt)
    if (identical(prealign_output_mode, "distinct") && "peak_ID" %in% names(prealign_out)) {
      prealign_out <- prealign_out[!duplicated(peak_ID)]
    }
    dir.create(dirname(prealign_path_use), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(tibble::as_tibble(prealign_out), prealign_path_use)
    if (isTRUE(prealign_only)) {
      return(list(fp_score = tibble::as_tibble(prealign_out), fp_bound = tibble::tibble(), fp_annotation = tibble::tibble(), id_map = tibble::tibble()))
    }
  }

  if ("peak_ID" %in% names(score_dt)) score_dt <- score_dt[!duplicated(peak_ID)]
  map_fp_atac <- unique(annot_dt[, .(fp_peak, atac_peak)])
  if (anyNA(map_fp_atac$atac_peak)) .log_abort("Found NA atac_peak in annotation.")
  dup <- map_fp_atac[, .N, by = fp_peak][N > 1L]
  if (nrow(dup)) .log_abort("Some fp_peak values map to multiple atac_peak values.")

  parse_peak <- function(x) {
    chr <- sub(":.*$", "", x, perl = TRUE)
    se <- sub("^[^:]+:", "", x, perl = TRUE)
    st <- as.integer(sub("-.*$", "", se, perl = TRUE))
    en <- as.integer(sub("^.*-", "", se, perl = TRUE))
    list(chr = chr, start = st, end = en)
  }
  pp <- parse_peak(score_dt$peak_ID)
  score_dt[, `:=`(chr = pp$chr, start = pp$start, end = pp$end)]
  score_dt[, mid := as.integer((start + end) %/% 2)]

  data.table::setkey(map_fp_atac, fp_peak)
  score_dt[map_fp_atac, atac_peak := i.atac_peak, on = .(peak_ID = fp_peak)]
  if (anyNA(score_dt$atac_peak)) .log_abort("Failed to map some footprint peaks to atac_peak.")

  sample_cols <- setdiff(names(score_dt), c("peak_ID", "atac_peak", "chr", "start", "end", "mid"))
  if (!length(sample_cols)) .log_abort("No sample columns detected in fp_score.")
  mult <- 10^round_digits
  sig_id_vec <- integer(nrow(score_dt))
  chr_groups <- split(seq_len(nrow(score_dt)), score_dt$chr, drop = TRUE)

  log_file_use <- log_file
  if (is.null(log_file_use) && !is.null(cache_dir)) {
    tag <- if (!is.null(cache_tag) && nzchar(cache_tag)) cache_tag else "cache"
    log_file_use <- file.path(cache_dir, sprintf("align_footprints_%s_log.txt", tag))
  }
  if (!is.null(log_file_use)) dir.create(dirname(log_file_use), recursive = TRUE, showWarnings = FALSE)
  log_line <- function(msg) {
    if (is.null(log_file_use)) return(invisible(NULL))
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write(paste(stamp, msg, sep = "\t"), file = log_file_use, append = TRUE)
    invisible(NULL)
  }

  use_parallel <- .Platform$OS.type == "unix" && threads > 1L
  use_inner_parallel <- use_parallel && identical(parallel_by, "atac_peak")

  process_chr <- function(chr, idx_chr) {
    if (!length(idx_chr)) return(list(idx_chr = idx_chr, sig_local = integer(0)))
    s_sig <- as.matrix(score_dt[idx_chr, ..sample_cols])
    s_sig <- round(s_sig * mult)
    storage.mode(s_sig) <- "integer"
    atac_chr <- score_dt$atac_peak[idx_chr]
    splits <- split(seq_len(length(idx_chr)), atac_chr, drop = TRUE)
    log_line(sprintf("chr=%s n_rows=%d n_atac_peaks=%d", chr, length(idx_chr), length(splits)))

    if (score_match_pct >= 0.9999) {
      tmp_dt <- data.table::as.data.table(s_sig)
      tmp_dt[, atac_peak := atac_chr]
      sig_cols <- names(tmp_dt)[seq_len(ncol(tmp_dt) - 1L)]
      tmp_dt[, sig_id := data.table::rleidv(.SD), by = .(atac_peak), .SDcols = sig_cols]
      sig_local <- tmp_dt$sig_id
    } else {
      k_req <- ceiling(score_match_pct * ncol(s_sig))
      ord <- order(vapply(splits, length, integer(1L)), decreasing = TRUE)
      splits <- splits[ord]
      if (use_inner_parallel) {
        n_workers <- max(1L, min(threads, length(splits)))
        if (length(splits) == 1L) {
          comp_list <- list(.assign_fp_score_components(s_sig[splits[[1L]], , drop = FALSE], k_req))
        } else {
          batch_ids <- parallel::splitIndices(length(splits), n_workers)
          comp_batches <- parallel::mclapply(
            batch_ids,
            function(idx) lapply(idx, function(i) .assign_fp_score_components(s_sig[splits[[i]], , drop = FALSE], k_req)),
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
        comp_list <- lapply(splits, function(ii) .assign_fp_score_components(s_sig[ii, , drop = FALSE], k_req))
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

  if (identical(parallel_by, "chromosome") && use_parallel) {
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
  score_dt[, `:=`(s_exp = mid - mid_slop, e_exp = mid + mid_slop)]
  data.table::setorder(score_dt, atac_peak, chr, sig_id, s_exp, e_exp)
  score_dt[, comp_id := {
    e <- e_exp
    s <- s_exp
    if (.N <= 1L) 1L else {
      prev <- c(-Inf, cummax(e[-.N]))
      cumsum(s > prev)
    }
  }, by = .(atac_peak, chr, sig_id)]

  score_dt[, `:=`(new_start = min(start), new_end = max(end), group_size = .N), by = .(atac_peak, chr, sig_id, comp_id)]
  score_dt[, new_peak_ID := paste0(chr, ":", new_start, "-", new_end)]

  map_old_new <- unique(score_dt[, .(peak_ID_old = peak_ID, atac_peak, new_peak_ID, group_size)])
  data.table::setkey(map_old_new, peak_ID_old)
  score_dt[map_old_new, peak_ID := i.new_peak_ID, on = .(peak_ID = peak_ID_old)]
  score_dt <- score_dt[, c("peak_ID", sample_cols), with = FALSE][!duplicated(peak_ID)]

  rm(sig_id_vec, chr_groups)
  gc(verbose = FALSE)

  if (isTRUE(verbose)) .log_inform("Loading and remapping footprint bound matrices.")
  bound_dt <- data.table::rbindlist(lapply(use$bound, read_bound), use.names = TRUE, fill = TRUE)
  if (!"peak_ID" %in% names(bound_dt) && "fp_peak" %in% names(bound_dt)) {
    data.table::setnames(bound_dt, "fp_peak", "peak_ID")
  }
  if ("peak_ID" %in% names(bound_dt)) bound_dt <- bound_dt[!duplicated(peak_ID)]
  bound_dt[map_old_new, peak_ID := i.new_peak_ID, on = .(peak_ID = peak_ID_old)]
  bound_dt <- bound_dt[!duplicated(peak_ID)]

  annot_dt[, source_fp_peak := fp_peak]
  annot_dt[map_old_new, fp_peak := i.new_peak_ID, on = .(fp_peak = peak_ID_old)]

  i_ids <- data.table::data.table(peak_ID = annot_dt$fp_peak)
  if (identical(output_mode, "distinct")) {
    i_ids <- unique(i_ids)
  }
  fp_score_out <- merge(i_ids, score_dt, by = "peak_ID", all.x = TRUE, sort = FALSE)
  fp_bound_out <- merge(i_ids, bound_dt, by = "peak_ID", all.x = TRUE, sort = FALSE)

  gs_new <- map_old_new[, .(group_size = max(group_size, na.rm = TRUE)), by = .(new_peak_ID, atac_peak)]
  data.table::setnames(gs_new, "new_peak_ID", "peak_ID")
  id_map_aligned <- if (isTRUE(return_id_map) || isTRUE(write_id_map)) {
    data.table::data.table(
      peak_ID = annot_dt$fp_peak,
      source_fp_peak = annot_dt$source_fp_peak,
      atac_peak = annot_dt$atac_peak
    )[gs_new, group_size := i.group_size, on = .(peak_ID)]
  } else {
    data.table::data.table()
  }

  annot_for_sites <- annot_dt
  annot_dt <- unique(annot_dt, by = c("fp_peak", "atac_peak", "motifs"))

  fp_sites <- data.table::data.table()
  source_fp_peak <- NULL
  if (isTRUE(write_fp_sites)) {
    site_long <- unique(annot_for_sites[, .(
      peak_ID = fp_peak,
      atac_peak = atac_peak,
      motif = motifs,
      source_fp_peak = source_fp_peak
    )])
    fp_sites <- site_long[, .(
      motifs_all = paste(sort(unique(as.character(motif[!is.na(motif)]))), collapse = ";"),
      n_motifs = data.table::uniqueN(motif, na.rm = TRUE),
      source_fp_peaks = if (isTRUE(include_source_fp_peaks)) paste(sort(unique(as.character(source_fp_peak[!is.na(source_fp_peak)]))), collapse = ";") else NA_character_,
      n_source_fp_peaks = data.table::uniqueN(source_fp_peak, na.rm = TRUE)
    ), by = .(peak_ID, atac_peak)]
    pp_sites <- parse_peak(fp_sites$peak_ID)
    fp_sites[, `:=`(chr = pp_sites$chr, start = pp_sites$start, end = pp_sites$end)]
    fp_sites[gs_new, group_size := i.group_size, on = .(peak_ID, atac_peak)]
    data.table::setcolorder(fp_sites, c(
      "peak_ID", "chr", "start", "end", "atac_peak", "group_size",
      "n_motifs", "motifs_all", "n_source_fp_peaks", "source_fp_peaks"
    ))
  }

  if (identical(output_mode, "distinct")) {
    fp_score_out <- fp_score_out[!duplicated(peak_ID)]
    fp_bound_out <- fp_bound_out[!duplicated(peak_ID)]
  }
  if (!identical(fp_score_out$peak_ID, fp_bound_out$peak_ID)) {
    .log_abort("fp_score and fp_bound peak_ID rows are not in the same order.")
  }

  fp_score_out <- fp_score_out[!duplicated(peak_ID)]
  fp_bound_out <- fp_bound_out[!duplicated(peak_ID)]

  written_paths <- list()
  written_format <- cache_format_primary
  if (!is.null(cache_paths) && isTRUE(write_cache)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    formats <- if (identical(cache_format, "both")) c("csv", "parquet") else cache_format
    for (fmt in formats) {
      paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, fmt)
      .aligned_fp_write_table(tibble::as_tibble(fp_bound_out), paths$fp_bound, fmt)
      .aligned_fp_write_table(tibble::as_tibble(fp_score_out), paths$fp_score, fmt)
      if (isTRUE(write_id_map)) {
        .aligned_fp_write_table(tibble::as_tibble(id_map_aligned), paths$id_map, fmt)
      }
      if (isTRUE(write_fp_sites)) {
        .aligned_fp_write_table(tibble::as_tibble(fp_sites), paths$fp_sites, fmt)
      }
      if (identical(fmt, written_format)) written_paths <- paths
    }
  }

  counts <- list(
    fp_score = nrow(fp_score_out),
    fp_bound = nrow(fp_bound_out),
    fp_annotation = 0L,
    id_map = if (isTRUE(return_id_map) || isTRUE(write_id_map)) nrow(id_map_aligned) else 0L,
    fp_sites = if (isTRUE(write_fp_sites)) nrow(fp_sites) else 0L
  )
  if (!isTRUE(return_data)) {
    paths <- if (length(written_paths)) written_paths else if (!is.null(cache_dir) && !is.null(cache_tag)) .aligned_fp_cache_paths(cache_dir, cache_tag, written_format) else list()
    return(cache_summary(paths, counts, written_format))
  }

  list(
    fp_score = tibble::as_tibble(fp_score_out),
    fp_bound = tibble::as_tibble(fp_bound_out),
    fp_annotation = tibble::as_tibble(annot_dt),
    id_map = if (isTRUE(return_id_map)) tibble::as_tibble(id_map_aligned) else tibble::tibble(),
    fp_sites = if (isTRUE(write_fp_sites)) tibble::as_tibble(fp_sites) else tibble::tibble(),
    cache_schema = "compact"
  )
}
