# File: utils_step1_io.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Provide Module 1 input loading, cache loading, and small Step 1 IO helpers.
#
# Inputs:
# - footprint directories and cache directories
# - ATAC, RNA, and metadata file paths
# - Module 1 config values
#
# Outputs:
# - Step 1 input tables and cached aligned footprint objects
#
# Notes:
# - Keep Step 1 IO separate from Step 1 model logic.

#' Module 1 IO helpers
#'
#' @noRd
NULL

.available_cores <- function(logical = TRUE) {
  cores <- suppressWarnings(parallel::detectCores(logical = logical))
  cores <- suppressWarnings(as.integer(cores))
  if (!is.finite(cores) || is.na(cores) || cores < 1L) {
    return(1L)
  }
  max(1L, cores)
}

.arrow_codec_available <- function(codec) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    return(FALSE)
  }
  isTRUE(tryCatch(arrow::codec_is_available(codec), error = function(...) FALSE))
}

.write_parquet_table <- function(x, path, compression = "zstd") {
  codec <- if (.arrow_codec_available(compression)) compression else "uncompressed"
  arrow::write_parquet(x, path, compression = codec)
  invisible(codec)
}

.aligned_fp_cache_paths <- function(cache_dir, cache_tag, format = c("csv", "parquet")) {
  format <- match.arg(format)
  ext <- if (identical(format, "parquet")) "parquet" else "csv"
  list(
    fp_bound = file.path(cache_dir, sprintf("fp_bounds_%s.%s", cache_tag, ext)),
    fp_score = file.path(cache_dir, sprintf("fp_scores_%s.%s", cache_tag, ext)),
    id_map = file.path(cache_dir, sprintf("fp_id_map_%s.%s", cache_tag, ext)),
    fp_sites = file.path(cache_dir, sprintf("fp_sites_%s.%s", cache_tag, ext))
  )
}

.aligned_fp_cache_required_exists <- function(cache_paths) {
  all(file.exists(unlist(cache_paths[c("fp_bound", "fp_score", "fp_sites")], use.names = FALSE)))
}

.aligned_fp_cache_schema <- function(cache_paths) {
  if (.aligned_fp_cache_required_exists(cache_paths)) "compact" else "missing"
}

.aligned_fp_cache_choose_format <- function(cache_dir,
                                            cache_tag,
                                            cache_format = c("auto", "parquet", "csv"),
                                            require_arrow = FALSE) {
  cache_format <- match.arg(cache_format)
  csv_paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, "csv")
  parquet_paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, "parquet")
  has_csv <- .aligned_fp_cache_required_exists(csv_paths)
  has_parquet <- .aligned_fp_cache_required_exists(parquet_paths)
  has_arrow <- requireNamespace("arrow", quietly = TRUE)

  if (identical(cache_format, "parquet")) {
    if (!has_arrow && isTRUE(require_arrow)) .log_abort("Package {.pkg arrow} is required to read Parquet aligned footprint cache.")
    return(list(format = "parquet", paths = parquet_paths, exists = has_parquet, arrow = has_arrow))
  }
  if (identical(cache_format, "csv")) {
    return(list(format = "csv", paths = csv_paths, exists = has_csv, arrow = has_arrow))
  }
  if (has_parquet && has_arrow) {
    return(list(format = "parquet", paths = parquet_paths, exists = TRUE, arrow = TRUE))
  }
  list(format = "csv", paths = csv_paths, exists = has_csv, arrow = has_arrow)
}

.aligned_fp_read_table <- function(path, format) {
  if (identical(format, "parquet")) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      .log_abort("Package {.pkg arrow} is required to read {.path {path}}.")
    }
    return(data.table::as.data.table(arrow::read_parquet(path)))
  }
  data.table::fread(path, showProgress = FALSE)
}

.aligned_fp_write_table <- function(x, path, format) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (identical(format, "parquet")) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      .log_warn("Package {.pkg arrow} is not installed; skipping Parquet cache {.path {path}}.")
      return(invisible(FALSE))
    }
    .write_parquet_table(x, path)
    return(invisible(TRUE))
  }
  data.table::fwrite(data.table::as.data.table(x), path)
  invisible(TRUE)
}

#' Load aligned footprint outputs from cache
#'
#' Reads cached aligned footprint tables produced by `align_footprints()`.
#'
#' @param cache_dir Directory containing cached aligned footprint files.
#' @param cache_tag Cache tag used in aligned footprint file names.
#' @param output_mode Output mode. One of `"full"` or `"distinct"`.
#' @param load_id_map Logical; if `TRUE`, load the optional footprint ID map.
#' @param cache_format Cache format. The default `"csv"` reads CSV caches.
#'   Use `"parquet"` to read Parquet caches, or `"auto"` to prefer Parquet
#'   when available and fall back to CSV.
#' @param verbose Logical; if `TRUE`, emit concise status messages.
#'
#' @return A list with `fp_score`, `fp_bound`, `fp_annotation`, and `id_map`.
#' @noRd
load_fp_aligned_from_cache <- function(
    cache_dir,
    cache_tag,
    output_mode = c("full", "distinct"),
    load_id_map = FALSE,
    cache_format = c("csv", "parquet", "auto"),
    expand_compact_annotation = TRUE,
    verbose = TRUE
) {
  .assert_pkg("data.table")

  output_mode <- match.arg(output_mode)
  cache_format <- match.arg(cache_format)
  stopifnot(is.logical(expand_compact_annotation), length(expand_compact_annotation) == 1L, !is.na(expand_compact_annotation))
  cache_info <- .aligned_fp_cache_choose_format(cache_dir, cache_tag, cache_format)
  cache_paths <- cache_info$paths

  if (!isTRUE(cache_info$exists)) {
    .log_abort("Cached aligned footprints not found for tag={cache_tag} in {cache_dir}.")
  }

  if (isTRUE(verbose)) {
    .log_inform("Using {cache_info$format} cached aligned footprints from {.path {cache_dir}} (tag = {cache_tag}).")
  }

  schema_detected <- .aligned_fp_cache_schema(cache_paths)

  fp_bound <- .aligned_fp_read_table(cache_paths$fp_bound, cache_info$format)
  fp_score <- .aligned_fp_read_table(cache_paths$fp_score, cache_info$format)
  fp_sites <- if (file.exists(cache_paths$fp_sites)) {
    .aligned_fp_read_table(cache_paths$fp_sites, cache_info$format)
  } else {
    data.table::data.table()
  }

  if (isTRUE(expand_compact_annotation) && nrow(fp_sites) > 0L && all(c("peak_ID", "atac_peak", "motifs_all") %in% names(fp_sites))) {
    motif_parts <- strsplit(as.character(fp_sites$motifs_all), ";", fixed = TRUE)
    motif_lens <- lengths(motif_parts)
    keep <- motif_lens > 0L
    fp_annotation <- data.table::data.table(
      fp_peak = rep(fp_sites$peak_ID[keep], motif_lens[keep]),
      atac_peak = rep(fp_sites$atac_peak[keep], motif_lens[keep]),
      motifs = unlist(motif_parts[keep], use.names = FALSE)
    )
  } else {
    fp_annotation <- data.table::data.table()
  }
  if (identical(cache_tag, "HOCOMOCOv13") && "motifs" %in% names(fp_annotation)) {
    fp_annotation$motifs <- .normalize_motif_id(fp_annotation$motifs, db_name = cache_tag)
  }
  id_map <- if (isTRUE(load_id_map) && file.exists(cache_paths$id_map)) {
    .aligned_fp_read_table(cache_paths$id_map, cache_info$format)
  } else {
    data.table::data.table()
  }

  out_cached <- list(
    fp_score = tibble::as_tibble(fp_score),
    fp_bound = tibble::as_tibble(fp_bound),
    fp_annotation = tibble::as_tibble(fp_annotation),
    fp_sites = tibble::as_tibble(fp_sites),
    id_map = tibble::as_tibble(id_map),
    cache_schema = schema_detected
  )

  if (identical(output_mode, "distinct")) {
    out_cached$fp_score <- out_cached$fp_score[!duplicated(out_cached$fp_score$peak_ID), , drop = FALSE]
    out_cached$fp_bound <- out_cached$fp_bound[!duplicated(out_cached$fp_bound$peak_ID), , drop = FALSE]
  }

  if (!identical(out_cached$fp_score$peak_ID, out_cached$fp_bound$peak_ID)) {
    .log_abort("`fp_score` and `fp_bound` peak_ID rows are not in the same order.")
  }

  out_cached
}

load_atac <- function(
    atac_data,
    chr_col = "Chr",
    start_col = "Start",
    end_col = "End",
    overlap_prefix = "Overlap_",
    metadata_cols = c(
      "PeakID", "Distance.to.TSS", "Nearest.PromoterID", "Entrez.ID",
      "Nearest.Unigene", "Nearest.Refseq", "Nearest.Ensembl",
      "Gene.Name", "Gene.Alias", "Gene.Description", "Gene.Type"
    ),
    sort_peaks = TRUE
) {
  if (!is.data.frame(atac_data)) .log_abort("`atac_data` must be a data.frame or tibble.")
  need <- c(chr_col, start_col, end_col)
  if (!all(need %in% names(atac_data))) {
    .log_abort("Missing coord columns: {.val {setdiff(need, names(atac_data))}}.")
  }

  df <- atac_data
  if (isTRUE(sort_peaks)) {
    chr <- tolower(gsub("^chr", "", as.character(df[[chr_col]])))
    chr_num <- suppressWarnings(as.integer(chr))
    chr_num[is.na(chr_num) & chr %in% c("x")] <- 23L
    chr_num[is.na(chr_num) & chr %in% c("y")] <- 24L
    chr_num[is.na(chr_num) & chr %in% c("m", "mt", "dmel_mito", "mito")] <- 25L
    ord <- order(
      is.na(chr_num),
      chr_num,
      chr,
      as.integer(round(df[[start_col]])),
      as.integer(round(df[[end_col]])),
      na.last = TRUE,
      method = "radix"
    )
    df <- df[ord, , drop = FALSE]
  }

  atac_peak <- paste0(
    df[[chr_col]],
    ":",
    as.integer(round(df[[start_col]])),
    "-",
    as.integer(round(df[[end_col]]))
  )

  overlap_cols <- grep(paste0("^", overlap_prefix), names(df), value = TRUE)
  meta_set <- unique(c(metadata_cols, chr_col, start_col, end_col))
  candidates <- setdiff(names(df), c(meta_set, overlap_cols))
  sample_cols <- candidates[vapply(df[candidates], is.numeric, logical(1))]
  if (!length(sample_cols)) {
    .log_abort("No numeric sample columns found after excluding metadata and {.val {overlap_prefix}}*.")
  }

  atac_score <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, sample_cols, drop = FALSE])

  atac_overlap <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, overlap_cols, drop = FALSE])

  if (length(overlap_cols)) {
    new_names <- sub(paste0("^", overlap_prefix), "", names(atac_overlap))
    new_names[1L] <- "atac_peak"
    names(atac_overlap) <- make.unique(new_names, sep = "_")
  }

  list(
    score = dplyr::as_tibble(atac_score),
    overlap = dplyr::as_tibble(atac_overlap)
  )
}

.filter_any_bound <- function(wide_tbl) {
  bd <- grep("_ATAC_bound$", names(wide_tbl), value = TRUE)
  if (!length(bd)) {
    return(tibble::as_tibble(wide_tbl[0, , drop = FALSE]))
  }
  keep <- rowSums(as.data.frame(wide_tbl[bd]), na.rm = TRUE) > 0
  tibble::as_tibble(wide_tbl[keep, , drop = FALSE])
}

.select_fp_scores <- function(x) {
  if (!"peak_ID" %in% names(x)) .log_abort("Input must contain a `peak_ID` column.")
  score_cols <- grep("_ATAC_score$", names(x), value = TRUE)
  if (!length(score_cols)) .log_abort("No columns ending with `_ATAC_score` were found.")
  out <- x[, c("peak_ID", score_cols), drop = FALSE]
  colnames(out) <- sub("_ATAC_score$", "", colnames(out))
  tibble::as_tibble(out)
}

.select_fp_bounds <- function(x) {
  if (!"peak_ID" %in% names(x)) .log_abort("Input must contain a `peak_ID` column.")
  bound_cols <- grep("_ATAC_bound$", names(x), value = TRUE)
  if (!length(bound_cols)) .log_abort("No columns ending with `_ATAC_bound` were found.")
  out <- x[, c("peak_ID", bound_cols), drop = FALSE]
  colnames(out) <- sub("_ATAC_bound$", "", colnames(out))
  tibble::as_tibble(out)
}

.select_fp_annots <- function(x) {
  if (!"peak_ID" %in% names(x)) .log_abort("Input must contain a `peak_ID` column.")
  if (!"peak_ATAC" %in% names(x)) .log_abort("Input must contain a `peak_ATAC` column.")
  if (!"TFBS_name" %in% names(x)) .log_abort("Input must contain a `TFBS_name` column.")
  out <- x[, c("peak_ID", "peak_ATAC", "TFBS_name"), drop = FALSE]
  out <- data.table::as.data.table(out)
  data.table::setnames(out, c("peak_ID", "peak_ATAC", "TFBS_name"), c("fp_peak", "atac_peak", "motifs"))
  tibble::as_tibble(out)
}


.normalize_motif_id <- function(x, db_name = NULL) {
  if (!is.character(x)) return(x)
  db_name <- if (!is.null(db_name) && nzchar(db_name)) as.character(db_name) else ""
  if (identical(db_name, "HOCOMOCOv13")) {
    return(sub("^_+", "", x))
  }
  x
}

discover_overview_motifs <- function(
    root_dir,
    db_name,
    sample_ids = NULL,
    n_samples = NULL,
    verbose = TRUE
) {
  if (!dir.exists(root_dir)) .log_abort("Root directory not found: {root_dir}")

  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) .log_abort("No subfolders were found under {root_dir}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) .log_abort("None of the requested samples exist under {root_dir}.")
  } else {
    sample_ids <- all_samples
  }
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) .log_abort("`n_samples` must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }

  if (isTRUE(verbose)) {
    .log_inform("Indexing motif overview files across {length(sample_ids)} sample(s) (DB = {db_name}).")
  }

  motifs_union <- unique(unlist(lapply(sample_ids, function(sid) {
    db_dir <- file.path(root_dir, sid, db_name)
    if (!dir.exists(db_dir)) return(character(0))
    subdirs <- base::list.dirs(db_dir, recursive = FALSE, full.names = FALSE)
    if (!length(subdirs)) return(character(0))
    has_file <- vapply(subdirs, function(m) {
      file.exists(file.path(db_dir, m, paste0(m, "_overview.txt")))
    }, logical(1))
    subdirs[has_file]
  })))

  if (!length(motifs_union)) .log_abort("No motif folders with overview files were found.")
  if (isTRUE(verbose)) .log_inform("Discovered {length(motifs_union)} motif(s) with overview files.")
  motifs_union
}

load_one_motif_wide <- function(
    root_dir,
    db_name,
    sample_ids,
    motif_id,
    verbose = TRUE
) {
  prev_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(prev_threads), add = TRUE)
  data.table::setDTthreads(1L)

  read_one_overview_dt <- function(file, sid) {
    fixed <- c(
      "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
      "peak_chr", "peak_start", "peak_end"
    )
    sc <- paste0(sid, "_ATAC_score")
    bd <- paste0(sid, "_ATAC_bound")
    need <- c(fixed, sc, bd)

    dt <- tryCatch(
      data.table::fread(file, select = need, showProgress = FALSE),
      error = function(e) NULL
    )
    if (is.null(dt)) {
      hdr <- data.table::fread(file, nrows = 0L, showProgress = FALSE)
      nm <- names(hdr)

      sc <- nm[grepl(paste0("^", sid, "_ATAC_score$"), nm)]
      bd <- nm[grepl(paste0("^", sid, "_ATAC_bound$"), nm)]
      if (!length(sc)) sc <- nm[grepl("_ATAC_score$", nm)]
      if (!length(bd)) bd <- nm[grepl("_ATAC_bound$", nm)]
      if (!length(sc)) sc <- nm[grepl(paste0("^", sid, ".*footprints_score$"), nm)]
      if (!length(bd)) bd <- nm[grepl(paste0("^", sid, ".*footprints_bound$"), nm)]
      if (!length(sc)) sc <- nm[grepl("footprints_score$", nm)]
      if (!length(bd)) bd <- nm[grepl("footprints_bound$", nm)]
      if (length(sc) != 1L || length(bd) != 1L) return(NULL)

      need <- c(fixed, sc, bd)
      if (!all(need %in% nm)) return(NULL)
      dt <- data.table::fread(file, select = need, showProgress = FALSE)
    }

    data.table::setnames(dt, old = sc, new = "ATAC_score")
    data.table::setnames(dt, old = bd, new = "ATAC_bound")
    dt[, peak_ID := paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end)]
    dt[, peak_ATAC := paste0(peak_chr, ":", peak_start, "-", peak_end)]
    dt[, `:=`(sample_id = sid)]
    dt[, .(sample_id, TFBS_name, peak_ID, peak_ATAC, ATAC_score, ATAC_bound)]
  }

  chunks <- list()
  samples_present <- character(0L)
  for (sid in sample_ids) {
    f <- file.path(root_dir, sid, db_name, motif_id, paste0(motif_id, "_overview.txt"))
    if (!file.exists(f)) next
    dtj <- read_one_overview_dt(f, sid)
    if (is.null(dtj) || !nrow(dtj)) next
    motif_id_norm <- .normalize_motif_id(motif_id, db_name = db_name)
    dtj[, TFBS_name := .normalize_motif_id(TFBS_name, db_name = db_name)]
    dtj <- dtj[TFBS_name == motif_id_norm]
    if (!nrow(dtj)) next
    dtj <- unique(dtj, by = c("TFBS_name", "peak_ID"))

    if (!sid %in% samples_present) samples_present <- c(samples_present, sid)
    chunks[[sid]] <- dtj
  }

  chunks <- Filter(Negate(is.null), chunks)
  if (!length(chunks)) {
    if (isTRUE(verbose)) .log_inform("Motif {motif_id}: no per-sample overview rows found.")
    return(tibble::as_tibble(data.frame()))
  }

  long <- data.table::rbindlist(chunks, use.names = TRUE)
  if (!nrow(long)) {
    if (isTRUE(verbose)) .log_inform("Motif {motif_id}: no per-sample overview rows found.")
    return(tibble::as_tibble(data.frame()))
  }
  key_cols <- c("peak_ID", "peak_ATAC", "TFBS_name")
  long <- unique(long, by = c(key_cols, "sample_id"))
  score_wide <- data.table::dcast(
    long,
    peak_ID + peak_ATAC + TFBS_name ~ sample_id,
    value.var = "ATAC_score"
  )
  bound_wide <- data.table::dcast(
    long,
    peak_ID + peak_ATAC + TFBS_name ~ sample_id,
    value.var = "ATAC_bound"
  )
  sample_score_cols <- intersect(samples_present, names(score_wide))
  sample_bound_cols <- intersect(samples_present, names(bound_wide))
  data.table::setnames(score_wide, sample_score_cols, paste0(sample_score_cols, "_ATAC_score"))
  data.table::setnames(bound_wide, sample_bound_cols, paste0(sample_bound_cols, "_ATAC_bound"))
  acc <- merge(score_wide, bound_wide, by = key_cols, all = TRUE, sort = FALSE)
  desired <- c(key_cols, paste0(samples_present, "_ATAC_score"), paste0(samples_present, "_ATAC_bound"))
  acc <- acc[, desired, with = FALSE]
  tibble::as_tibble(acc)
}

.count_rows_fast <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  tryCatch({
    dt <- data.table::fread(path, select = 1L, showProgress = FALSE)
    as.integer(nrow(dt))
  }, error = function(e) NA_integer_)
}

.fp_cache_triplet_status <- function(score_path, bound_path, annot_path) {
  paths <- c(score = score_path, bound = bound_path, annot = annot_path)
  exists <- file.exists(paths)
  if (!all(exists)) {
    reason <- if (!any(exists)) {
      "not_found"
    } else {
      paste0("partial_missing:", paste0(names(paths)[!exists], collapse = ","))
    }
    return(list(valid = FALSE, n_rows = NA_integer_, reason = reason))
  }
  n_score <- .count_rows_fast(score_path)
  n_bound <- .count_rows_fast(bound_path)
  n_annot <- .count_rows_fast(annot_path)
  counts <- c(score = n_score, bound = n_bound, annot = n_annot)
  if (anyNA(counts)) {
    return(list(valid = FALSE, n_rows = NA_integer_, reason = "unreadable"))
  }
  if (!identical(as.integer(n_score), as.integer(n_bound)) || !identical(as.integer(n_score), as.integer(n_annot))) {
    return(list(valid = FALSE, n_rows = as.integer(n_score), reason = sprintf("row_mismatch:score=%d,bound=%d,annot=%d", n_score, n_bound, n_annot)))
  }
  list(valid = TRUE, n_rows = as.integer(n_score), reason = "complete")
}

load_footprints <- function(
    root_dir,
    db_name,
    out_dir = NULL,
    sample_ids = NULL,
    n_samples = NULL,
    motif_ids = NULL,
    n_motifs = NULL,
    n_workers = .available_cores(logical = TRUE),
    set_plan = TRUE,
    skip_existing = TRUE,
    verbose = TRUE,
    log_file = NULL,
    memory_check = TRUE,
    memory_check_n_motifs = 3L,
    memory_check_factor = 2.5,
    memory_check_fraction = 0.75,
    ask_on_risk = TRUE
) {
  .assert_pkg("data.table")
  .assert_pkg("future")
  .assert_pkg("future.apply")
  .assert_pkg("readr")

  can_use_manifest_before_scan <- isTRUE(skip_existing) &&
    !is.null(out_dir) &&
    is.null(sample_ids) &&
    is.null(n_samples) &&
    is.null(motif_ids) &&
    is.null(n_motifs)
  if (can_use_manifest_before_scan) {
    manifest_path <- file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
    if (file.exists(manifest_path)) {
      man <- tryCatch(readr::read_csv(manifest_path, show_col_types = FALSE), error = function(e) NULL)
      required_cols <- c("motif", "n_peaks", "score", "bound", "annot")
      if (!is.null(man) &&
          all(required_cols %in% names(man)) &&
          nrow(man) > 0 &&
          !any(is.na(man$n_peaks))) {
        cache_ok <- vapply(seq_len(nrow(man)), function(i) {
          status <- .fp_cache_triplet_status(man$score[[i]], man$bound[[i]], man$annot[[i]])
          isTRUE(status$valid) && identical(as.integer(status$n_rows), as.integer(man$n_peaks[[i]]))
        }, logical(1))
        if (all(cache_ok)) {
          if (isTRUE(verbose)) .log_inform("Found existing manifest at {.path {manifest_path}}; returning cached summary.")
          attr(man, "from_cache") <- TRUE
          return(man)
        }
      }
    }
  }

  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) .log_abort("No subfolders under {root_dir}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) .log_abort("None of the requested samples exist under {root_dir}.")
  } else {
    sample_ids <- all_samples
  }
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) .log_abort("`n_samples` must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }

  motifs_union <- discover_overview_motifs(root_dir, db_name, sample_ids, NULL, verbose = verbose)
  if (!is.null(motif_ids)) {
    motifs_union <- intersect(motifs_union, motif_ids)
    if (!length(motifs_union)) .log_abort("Requested motif_ids were not found under the selected samples.")
  }
  if (!is.null(n_motifs)) {
    if (!is.numeric(n_motifs) || n_motifs < 1) .log_abort("`n_motifs` must be a positive integer.")
    motifs_union <- utils::head(motifs_union, n_motifs)
  }

  if (isTRUE(memory_check) && is.finite(n_workers) && n_workers > 1L) {
    .available_mem_bytes <- function() {
      os <- Sys.info()[["sysname"]]
      if (identical(os, "Linux") && file.exists("/proc/meminfo")) {
        meminfo <- readLines("/proc/meminfo", warn = FALSE)
        line <- meminfo[grepl("^MemAvailable:", meminfo)]
        if (length(line)) {
          kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", line[1])))
          if (is.finite(kb)) return(kb * 1024)
        }
      }
      if (identical(os, "Windows")) {
        lim <- suppressWarnings(utils::memory.limit())
        if (is.finite(lim) && lim > 0) return(as.numeric(lim) * 1024^2)
      }
      NA_real_
    }
    .fmt_bytes <- function(x) {
      if (!is.finite(x)) return("NA")
      units <- c("B", "KB", "MB", "GB", "TB")
      p <- max(0, min(length(units) - 1L, floor(log(x, 1024))))
      sprintf("%.2f %s", x / (1024^p), units[p + 1L])
    }
    mem_avail <- .available_mem_bytes()
    if (is.finite(mem_avail)) {
      n_sample <- max(1L, min(as.integer(memory_check_n_motifs), length(motifs_union)))
      motifs_sample <- motifs_union[seq_len(n_sample)]
      sizes <- vapply(motifs_sample, function(m) {
        paths <- file.path(root_dir, sample_ids, db_name, m, paste0(m, "_overview.txt"))
        info <- file.info(paths)
        sum(info$size[is.finite(info$size)], na.rm = TRUE)
      }, numeric(1))
      sizes <- sizes[is.finite(sizes) & sizes > 0]
      if (length(sizes)) {
        est_per_worker <- stats::median(sizes) * memory_check_factor
        est_total <- est_per_worker * as.numeric(n_workers)
        budget <- mem_avail * memory_check_fraction
        if (est_total > budget) {
          msg <- paste(
            "Potential OOM risk detected.",
            sprintf("Median total size across sampled motifs = %s.", .fmt_bytes(stats::median(sizes))),
            sprintf("Estimated peak RAM = %s (per worker ~ %s; workers = %d).", .fmt_bytes(est_total), .fmt_bytes(est_per_worker), as.integer(n_workers)),
            sprintf("Available RAM = %s; safety budget = %s.", .fmt_bytes(mem_avail), .fmt_bytes(budget)),
            "Consider reducing n_workers or using skip_existing = TRUE."
          )
          .log_warn(msg)
          if (isTRUE(ask_on_risk) && interactive()) {
            ans <- readline("OOM risk detected. Continue? [y/N]: ")
            if (!nzchar(ans) || !tolower(substr(ans, 1, 1)) %in% "y") {
              .log_abort("Aborted by user due to OOM risk.")
            }
          }
        }
      }
    }
  }

  manifest_path <- NULL
  if (!is.null(out_dir)) {
    manifest_path <- file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
    if (skip_existing && file.exists(manifest_path)) {
      manifest_ok <- FALSE
      man <- tryCatch(readr::read_csv(manifest_path, show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(man)) {
        required_cols <- c("motif", "n_peaks", "score", "bound", "annot")
        if (all(required_cols %in% names(man))) {
          motifs_in_manifest <- unique(man$motif)
          coverage_ok <- all(motifs_union %in% motifs_in_manifest)
          if (coverage_ok) {
            man <- man[man$motif %in% motifs_union, , drop = FALSE]
            cache_ok <- vapply(seq_len(nrow(man)), function(i) {
              status <- .fp_cache_triplet_status(man$score[[i]], man$bound[[i]], man$annot[[i]])
              isTRUE(status$valid) && identical(as.integer(status$n_rows), as.integer(man$n_peaks[[i]]))
            }, logical(1))
            if (nrow(man) > 0 && !any(is.na(man$n_peaks)) && all(cache_ok)) {
              manifest_ok <- TRUE
            }
          }
        }
      }
      if (manifest_ok) {
        if (isTRUE(verbose)) .log_inform("Found existing manifest at {.path {manifest_path}}; returning cached summary.")
        attr(man, "from_cache") <- TRUE
        return(man)
      }
    }
  }

  if (isTRUE(verbose)) {
    .log_inform("Streaming {length(motifs_union)} motif(s) with up to {length(sample_ids)} sample(s) each.")
    .log_inform("Writing outputs to {.path {out_dir}}")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(set_plan)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1") ||
      (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !is_rstudio &&
        isTRUE(requireNamespace("parallelly", quietly = TRUE)) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore, workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
  }

  have_progressr <- requireNamespace("progressr", quietly = TRUE)
  use_parallel <- isTRUE(n_workers > 1L)
  log_file_use <- log_file
  if (use_parallel && is.null(log_file_use) && !is.null(out_dir)) {
    log_file_use <- file.path(dirname(out_dir), paste0(basename(out_dir), "_load_footprints_log.txt"))
  }
  if (!is.null(log_file_use)) dir.create(dirname(log_file_use), recursive = TRUE, showWarnings = FALSE)
  log_line <- function(msg) {
    if (is.null(log_file_use)) return(invisible(NULL))
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    write(paste(stamp, msg, sep = "\t"), file = log_file_use, append = TRUE)
    invisible(NULL)
  }

  process_one_motif <- function(m) {
    base <- .normalize_motif_id(m, db_name = db_name)
    f_score <- file.path(out_dir, paste0(base, "_score.csv"))
    f_bound <- file.path(out_dir, paste0(base, "_bound.csv"))
    f_annot <- file.path(out_dir, paste0(base, "_annotation.csv"))

    if (isTRUE(skip_existing)) {
      cache_status <- .fp_cache_triplet_status(f_score, f_bound, f_annot)
      if (isTRUE(cache_status$valid)) {
        if (!is.null(log_file_use)) {
          log_line(sprintf("Motif %s: using cached outputs (n_peaks=%s).", m, cache_status$n_rows))
        }
        return(tibble::tibble(motif = base, n_peaks = cache_status$n_rows, score = f_score, bound = f_bound, annot = f_annot))
      } else if (!is.null(log_file_use)) {
        if (identical(cache_status$reason, "not_found")) {
          log_line(sprintf("Motif %s: preparing cache outputs for this run.", m))
        } else {
          log_line(sprintf("Motif %s: refreshing partial cache (%s).", m, cache_status$reason))
        }
      }
    }

    wd <- load_one_motif_wide(root_dir, db_name, sample_ids, m, verbose = FALSE)
    if (!nrow(wd)) {
      data.table::fwrite(data.table::data.table(peak_ID = character()), f_score)
      data.table::fwrite(data.table::data.table(peak_ID = character()), f_bound)
      data.table::fwrite(data.table::data.table(fp_peak = character(), atac_peak = character(), motifs = character()), f_annot)
      return(tibble::tibble(motif = base, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    anyb <- .filter_any_bound(wd)
    if (!nrow(anyb)) {
      data.table::fwrite(data.table::as.data.table(.select_fp_scores(wd)[0, ]), f_score)
      data.table::fwrite(data.table::as.data.table(.select_fp_bounds(wd)[0, ]), f_bound)
      data.table::fwrite(data.table::as.data.table(.select_fp_annots(wd)[0, ]), f_annot)
      return(tibble::tibble(motif = base, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    score_tbl <- .select_fp_scores(anyb)
    bound_tbl <- .select_fp_bounds(anyb)
    annot_tbl <- .select_fp_annots(anyb)
    n_out <- nrow(score_tbl)
    data.table::fwrite(data.table::as.data.table(score_tbl), f_score)
    data.table::fwrite(data.table::as.data.table(bound_tbl), f_bound)
    data.table::fwrite(data.table::as.data.table(annot_tbl), f_annot)

    rm(wd, anyb, score_tbl, bound_tbl, annot_tbl)
    gc()
    tibble::tibble(motif = base, n_peaks = n_out, score = f_score, bound = f_bound, annot = f_annot)
  }

  if (isTRUE(have_progressr)) {
    results <- progressr::with_progress({
      p <- progressr::progressor(along = motifs_union)
      future.apply::future_lapply(motifs_union, function(m) {
        if (isTRUE(verbose)) p(message = sprintf("Processing motif %s", m)) else p()
        process_one_motif(m)
      }, future.seed = FALSE)
    })
  } else {
    results <- future.apply::future_lapply(motifs_union, process_one_motif, future.seed = FALSE)
  }

  out <- dplyr::bind_rows(results)
  if (!is.null(out_dir)) {
    manifest_path <- file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
    data.table::fwrite(data.table::as.data.table(out), manifest_path)
  }
  attr(out, "from_cache") <- FALSE
  out
}
