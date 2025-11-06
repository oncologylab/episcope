# utils_tobias_overview_loader_simple.R — fast, simple loader to wide table
# Author: Yaoxiang Li
# Updated: 2025-10-14
#
# Goal
# ----
# Read all "<motif>/<motif>_overview.txt" under each sample's motif DB folder and
# return ONE de-duplicated wide tibble:
#   columns: peak_ID, TFBS_name, <sample>_ATAC_score, <sample>_ATAC_bound (for each sample)
#
# Key features
# ------------
# - Limit how many samples to process with `n_samples` (e.g., 2 for a quick run).
# - Simple and robust: minimal logic, explicit namespacing, de-duplicates rows.
# - Keeps only the requested columns to keep memory down and I/O fast.

#' Load TOBIAS *_overview.txt across samples to a wide table (optionally limit samples)
#'
#' @param root_dir Character. Path with per-sample folders (e.g., "Y:/.../TOBIAS_merged_peaks").
#' @param db_name  Character. Motif DB folder name inside each sample (e.g., "JASPAR2024", "HOCOMOCOv13").
#' @param sample_ids Optional character vector of sample folder names to include. Default: discover all.
#' @param n_samples Optional integer; if set, only the first `n_samples` (after filtering) are processed.
#' @param out_tsv Optional path to write the wide table as TSV. Default NULL (do not write).
#' @param verbose Logical; print progress messages. Default TRUE.
#'
#' @return A tibble with columns:
#'   - peak_ID (chr:start-end from TFBS_chr/start/end)
#'   - TFBS_name
#'   - <sample>_ATAC_score, <sample>_ATAC_bound for each processed sample
#'
#' @examples
#' \dontrun{
#' wide <- load_tobias_overviews_wide(
#'   root_dir  = "Y:/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#'   db_name   = "JASPAR2024",
#'   n_samples = 2,       # quick test on 2 samples
#'   verbose   = TRUE
#' )
#' }
# Minimal change: also keep peak_chr/peak_start/peak_end as "peak_ATAC" (chr:start-end)
# Updated functions only — still Windows-friendly (future) and returns TIBBLES.

#' Load TOBIAS *_overview.txt to a wide tibble (future-parallel) + keep peak_ATAC
load_tobias_overviews_wide <- function(root_dir,
                                       db_name,
                                       sample_ids = NULL,
                                       n_samples = NULL,
                                       n_workers = max(1L, parallel::detectCores(logical = TRUE)),
                                       set_plan = TRUE,
                                       out_tsv = NULL,
                                       verbose = TRUE) {
  if (!dir.exists(root_dir)) cli::cli_abort("Root directory not found: {.path {root_dir}}")
  if (!is.character(db_name) || length(db_name) != 1L || !nzchar(db_name)) cli::cli_abort("{.arg db_name} must be a non-empty single string.")

  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) cli::cli_abort("No subfolders under {.path {root_dir}}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) cli::cli_abort("None of the requested samples exist under {.path {root_dir}}.")
  } else sample_ids <- all_samples
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) cli::cli_abort("{.arg n_samples} must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }
  if (verbose) cli::cli_inform("Scanning {length(sample_ids)} sample(s) under {.path {root_dir}} (DB = {db_name})...")

  # ---- per-file reader (inside future workers) ----
  read_one_overview_dt <- function(file, sid) {
    data.table::setDTthreads(1L)  # avoid oversubscription in workers

    hdr <- data.table::fread(file, nrows = 0L, showProgress = FALSE)
    nm  <- names(hdr)

    sc <- nm[grepl(paste0("^", sid, "_ATAC_score$"), nm)]; if (!length(sc)) sc <- nm[grepl("_ATAC_score$", nm)]
    bd <- nm[grepl(paste0("^", sid, "_ATAC_bound$"), nm)]; if (!length(bd)) bd <- nm[grepl("_ATAC_bound$", nm)]
    if (length(sc) != 1L || length(bd) != 1L) return(NULL)

    # KEEP peak_chr/peak_start/peak_end to build peak_ATAC
    need <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
              "peak_chr","peak_start","peak_end", sc, bd)
    if (!all(need %in% nm)) return(NULL)

    dt <- data.table::fread(file, select = match(need, nm), showProgress = FALSE)
    data.table::setnames(dt, old = sc, new = "ATAC_score")
    data.table::setnames(dt, old = bd, new = "ATAC_bound")

    dt[, peak_ID  := paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end)]
    dt[, peak_ATAC:= paste0(peak_chr,  ":", peak_start, "-", peak_end)]
    dt[, `:=`(sample_id = sid)]
    dt[, .(sample_id, TFBS_name, peak_ID, peak_ATAC, ATAC_score, ATAC_bound)]
  }

  load_one_sample <- function(sid) {
    db_dir <- file.path(root_dir, sid, db_name)
    if (!dir.exists(db_dir)) { if (verbose) cli::cli_inform("Sample {sid}: missing {.path {db_dir}} — skipping."); return(NULL) }
    files <- list.files(db_dir, pattern = "_overview\\.txt$", full.names = TRUE, recursive = TRUE)
    if (!length(files)) { if (verbose) cli::cli_inform("Sample {sid}: no overview files — skipping."); return(NULL) }

    parts <- future.apply::future_lapply(files, function(f) read_one_overview_dt(f, sid), future.seed = FALSE)
    dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
    if (!nrow(dt)) return(NULL)
    dt <- unique(dt, by = c("sample_id", "TFBS_name", "peak_ID"))
    if (verbose) cli::cli_inform("Sample {sid}: loaded {nrow(dt)} rows from {length(files)} overview file(s).")
    dt
  }

  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)

    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1") ||
      (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())

    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !is_rstudio &&
        isTRUE(requireNamespace("parallelly", quietly = TRUE)) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore,   workers = n_workers)  # fast on Linux when not in RStudio
    } else {
      future::plan(future::multisession, workers = n_workers)  # safe everywhere incl. RStudio/Windows
    }
  }


  long_list <- lapply(sample_ids, load_one_sample)
  long_dt   <- data.table::rbindlist(long_list, use.names = TRUE, fill = TRUE)
  if (!nrow(long_dt)) cli::cli_abort("No rows were read. Check paths/db name.")
  long_dt <- unique(long_dt, by = c("sample_id", "TFBS_name", "peak_ID"))

  # ---- pivot to WIDE; include peak_ATAC as an ID ----
  wide_dt <- data.table::dcast(
    long_dt,
    peak_ID + peak_ATAC + TFBS_name ~ sample_id,
    value.var = c("ATAC_score", "ATAC_bound")
  )

  # rename ATAC_score_cy249 -> cy249_ATAC_score (and bound likewise)
  sc_old <- grep("^ATAC_score_", names(wide_dt), value = TRUE)
  bd_old <- grep("^ATAC_bound_", names(wide_dt), value = TRUE)
  sc_new <- sub("^ATAC_score_", "", sc_old)
  bd_new <- sub("^ATAC_bound_", "", bd_old)
  data.table::setnames(wide_dt, sc_old, paste0(sc_new, "_ATAC_score"))
  data.table::setnames(wide_dt, bd_old, paste0(bd_new, "_ATAC_bound"))

  samples_present <- unique(c(sc_new, bd_new))
  desired <- c("peak_ID", "peak_ATAC", "TFBS_name",
               paste0(samples_present, "_ATAC_score"),
               paste0(samples_present, "_ATAC_bound"))
  have <- intersect(desired, names(wide_dt))
  wide_dt <- wide_dt[, ..have]

  if (!is.null(out_tsv)) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(tibble::as_tibble(wide_dt), out_tsv)
    if (verbose) cli::cli_inform("Wrote combined table to {.path {out_tsv}}")
  }
  tibble::as_tibble(wide_dt)
}

# ---- helpers (unchanged API) now preserve/return peak_ATAC too ----

filter_any_bound <- function(wide_tbl) {
  bd <- grep("_ATAC_bound$", names(wide_tbl), value = TRUE)
  if (!length(bd)) return(tibble::as_tibble(wide_tbl[0, ]))
  tibble::as_tibble(wide_tbl[rowSums(as.data.frame(wide_tbl[bd])) > 0, ])
}

collapse_by_peak_fast <- function(any_bound_tbl) {
  dt <- data.table::as.data.table(any_bound_tbl)
  score_cols <- grep("_ATAC_score$", names(dt), value = TRUE)
  bound_cols <- grep("_ATAC_bound$", names(dt), value = TRUE)

  # keep first peak_ATAC + first sample columns per peak_ID
  first_dt <- dt[, .SD[1L], by = peak_ID, .SDcols = c("peak_ATAC", score_cols, bound_cols)]
  name_dt  <- dt[, .(TFBS_name = paste(sort(unique(TFBS_name)), collapse = ", ")), by = peak_ID]
  out_dt   <- first_dt[name_dt, on = "peak_ID"]
  data.table::setcolorder(out_dt, c("peak_ID", "peak_ATAC", "TFBS_name", score_cols, bound_cols))
  tibble::as_tibble(out_dt)
}

# res <- load_tobias_overviews_wide(
#   root_dir = "Y:/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#   db_name  = "JASPAR2024",
#   n_samples = 2,       # quick test on 2 samples
#   verbose  = TRUE
# )


# rows where any *_ATAC_bound == 1
# res_any_bound <- filter_any_bound(res)

.select_fp_scores <- function(x) {
  # ---- checks ----
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  score_cols <- names(x)[tidyselect::eval_select(rlang::expr(tidyselect::ends_with("_ATAC_score")), x)]
  if (length(score_cols) == 0) {
    cli::cli_abort("No columns ending with {.val _ATAC_score} were found.")
  }

  # ---- select + rename ----
  out <- dplyr::select(x, "peak_ID", tidyselect::ends_with("_ATAC_score"))
  out <- dplyr::rename_with(out, ~ sub("_ATAC_score$", "", .x), tidyselect::ends_with("_ATAC_score"))

  # Return tibble
  tibble::as_tibble(out)
}
.select_fp_bounds <- function(x) {
  # ---- checks ----
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  score_cols <- names(x)[tidyselect::eval_select(rlang::expr(tidyselect::ends_with("_ATAC_bound")), x)]
  if (length(score_cols) == 0) {
    cli::cli_abort("No columns ending with {.val _ATAC_bound} were found.")
  }

  # ---- select + rename ----
  out <- dplyr::select(x, "peak_ID", tidyselect::ends_with("_ATAC_bound"))
  out <- dplyr::rename_with(out, ~ sub("_ATAC_bound$", "", .x), tidyselect::ends_with("_ATAC_bound"))

  # Return tibble
  tibble::as_tibble(out)
}
.select_fp_annots <- function(x) {
  # ---- checks ----
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  if (!"peak_ATAC" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  if (!"TFBS_name" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  out <- x |>
    dplyr::select(peak_ID, peak_ATAC, TFBS_name) |>
    dplyr::rename(fp_peak = peak_ID, atac_peak = peak_ATAC, motifs = TFBS_name)

  # Return tibble
  tibble::as_tibble(out)
}

# .select_fp_scores(res_any_bound)
# .select_fp_bounds(res_any_bound)
# .select_fp_annots(res_any_bound)

# Streaming per-motif writer (low memory, parallel-safe) -----------------

#' Discover motif folders that have "<motif>_overview.txt" in at least one sample
#'
#' @param root_dir Path containing per-sample folders.
#' @param db_name  Motif DB folder name inside each sample.
#' @param sample_ids Optional subset of sample folder names. Default: discover all.
#' @param n_samples Optional limit on number of samples (after filtering).
#' @param verbose   Print progress.
#' @return Character vector of motif folder names.
discover_overview_motifs <- function(root_dir,
                                     db_name,
                                     sample_ids = NULL,
                                     n_samples = NULL,
                                     verbose = TRUE) {
  if (!dir.exists(root_dir)) cli::cli_abort("Root directory not found: {.path {root_dir}}")

  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) cli::cli_abort("No subfolders under {.path {root_dir}}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) cli::cli_abort("None of the requested samples exist under {.path {root_dir}}.")
  } else sample_ids <- all_samples
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) cli::cli_abort("{.arg n_samples} must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }

  if (verbose) cli::cli_inform("Indexing *_overview.txt across {length(sample_ids)} sample(s) (DB = {db_name})...")

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

  if (!length(motifs_union)) cli::cli_abort("No motif folders with *_overview.txt found.")

  if (verbose) cli::cli_inform("Discovered {length(motifs_union)} motif(s) with overview files.")
  motifs_union
}

#' Load ONE motif across ALL samples into a wide tibble (scores/bounds)
#'
#' This is intentionally single-threaded inside to minimize oversubscription.
#' It reads only "<sid>/<db>/<motif>/<motif>_overview.txt" files.
#'
#' @param root_dir,db_name,sample_ids As above.
#' @param motif_id Single motif folder name (e.g., "AhrArnt_MA0006.2").
#' @param verbose Print progress.
#' @return Tibble with columns: peak_ID, peak_ATAC, TFBS_name,
#'         <sample>_ATAC_score, <sample>_ATAC_bound. Returns 0-row tibble if none.
load_one_motif_wide <- function(root_dir,
                                db_name,
                                sample_ids,
                                motif_id,
                                verbose = TRUE) {
  data.table::setDTthreads(1L)

  # tiny reader for this motif+sample
  read_one_overview_dt <- function(file, sid) {
    hdr <- data.table::fread(file, nrows = 0L, showProgress = FALSE)
    nm  <- names(hdr)

    sc <- nm[grepl(paste0("^", sid, "_ATAC_score$"), nm)]; if (!length(sc)) sc <- nm[grepl("_ATAC_score$", nm)]
    bd <- nm[grepl(paste0("^", sid, "_ATAC_bound$"), nm)]; if (!length(bd)) bd <- nm[grepl("_ATAC_bound$", nm)]
    if (length(sc) != 1L || length(bd) != 1L) return(NULL)

    need <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
              "peak_chr","peak_start","peak_end", sc, bd)
    if (!all(need %in% nm)) return(NULL)

    dt <- data.table::fread(file, select = match(need, nm), showProgress = FALSE)
    data.table::setnames(dt, old = sc, new = "ATAC_score")
    data.table::setnames(dt, old = bd, new = "ATAC_bound")

    dt[, peak_ID  := paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end)]
    dt[, peak_ATAC:= paste0(peak_chr,  ":", peak_start, "-", peak_end)]
    dt[, `:=`(sample_id = sid)]
    dt[, .(sample_id, TFBS_name, peak_ID, peak_ATAC, ATAC_score, ATAC_bound)]
  }

  parts <- vector("list", length(sample_ids))
  kept  <- 0L
  for (j in seq_along(sample_ids)) {
    sid <- sample_ids[[j]]
    f <- file.path(root_dir, sid, db_name, motif_id, paste0(motif_id, "_overview.txt"))
    if (!file.exists(f)) next
    dtj <- read_one_overview_dt(f, sid)
    if (!is.null(dtj) && nrow(dtj)) {
      parts[[j]] <- dtj
      kept <- kept + 1L
    }
  }
  if (kept == 0L) {
    if (verbose) cli::cli_inform("Motif {motif_id}: no per-sample overview rows found — skipping.")
    return(tibble::as_tibble(data.frame()))
  }

  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) return(tibble::as_tibble(data.frame()))
  # enforce exact motif match (defensive)
  dt <- dt[TFBS_name == motif_id]
  dt <- unique(dt, by = c("sample_id", "TFBS_name", "peak_ID"))

  # cast wide
  wd <- data.table::dcast(
    dt,
    peak_ID + peak_ATAC + TFBS_name ~ sample_id,
    value.var = c("ATAC_score", "ATAC_bound")
  )

  # rename ATAC_score_cyXXX -> cyXXX_ATAC_score (likewise bound)
  sc_old <- grep("^ATAC_score_", names(wd), value = TRUE)
  bd_old <- grep("^ATAC_bound_", names(wd), value = TRUE)
  sc_new <- sub("^ATAC_score_", "", sc_old)
  bd_new <- sub("^ATAC_bound_", "", bd_old)
  data.table::setnames(wd, sc_old, paste0(sc_new, "_ATAC_score"))
  data.table::setnames(wd, bd_old, paste0(bd_new, "_ATAC_bound"))

  samples_present <- unique(c(sc_new, bd_new))
  desired <- c("peak_ID", "peak_ATAC", "TFBS_name",
               paste0(samples_present, "_ATAC_score"),
               paste0(samples_present, "_ATAC_bound"))
  have <- intersect(desired, names(wd))
  wd <- wd[, ..have]

  tibble::as_tibble(wd)
}

#' Stream motifs: load -> filter(any bound) -> write per motif CSVs (parallel)
#'
#' This function avoids building a giant all-motifs table. It:
#'  1) Discovers motifs (or uses user-provided list)
#'  2) In parallel (per motif): loads one motif across all samples,
#'     filters to rows with any *_ATAC_bound == 1, then calls
#'     `write_fp_tables_by_motif()` to emit "<motif>_{score,bound,annotation}.csv"
#'     to `out_dir`.
#'
#' @param root_dir,db_name Paths as above.
#' @param out_dir Output directory (e.g., "inst/extdata/fp_scores_jaspar2024").
#' @param sample_ids,n_samples Subset/limit samples.
#' @param motif_ids,n_motifs   Subset/limit motifs (folder names).
#' @param n_workers Parallel workers (motif-level). Default: all logical cores.
#' @param set_plan  Whether to set a future plan. Default TRUE.
#' @param skip_existing If TRUE, skip a motif when all three output files already exist.
#' @param verbose   Print progress.
#' @return A tibble summarizing outputs: motif, n_peaks, score, bound, annot.
stream_write_overviews_by_motif <- function(root_dir,
                                            db_name,
                                            out_dir = file.path("inst","extdata","fp_scores_jaspar2024"),
                                            sample_ids = NULL,
                                            n_samples = NULL,
                                            motif_ids = NULL,
                                            n_motifs = NULL,
                                            n_workers = max(1L, parallel::detectCores(logical = TRUE)),
                                            set_plan = TRUE,
                                            skip_existing = TRUE,
                                            verbose = TRUE) {
  # discover samples (reusing logic from discover_overview_motifs)
  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) cli::cli_abort("No subfolders under {.path {root_dir}}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) cli::cli_abort("None of the requested samples exist under {.path {root_dir}}.")
  } else sample_ids <- all_samples
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) cli::cli_abort("{.arg n_samples} must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }

  motifs_union <- discover_overview_motifs(root_dir, db_name, sample_ids, NULL, verbose = verbose)
  if (!is.null(motif_ids)) {
    motifs_union <- intersect(motifs_union, motif_ids)
    if (!length(motifs_union)) cli::cli_abort("Requested motif_ids not found under selected samples.")
  }
  if (!is.null(n_motifs)) {
    if (!is.numeric(n_motifs) || n_motifs < 1) cli::cli_abort("{.arg n_motifs} must be a positive integer.")
    motifs_union <- utils::head(motifs_union, n_motifs)
  }

  if (verbose) {
    cli::cli_inform("Streaming {length(motifs_union)} motif(s) with up to {length(sample_ids)} sample(s) each...")
    cli::cli_inform("Writing outputs to {.path {out_dir}}")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # filename sanitizer consistent with your previous writer
  sanitize_fname <- function(s) chartr("/\\", "__", s)

  # set future plan if requested
  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)

    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1") ||
      (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())

    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !is_rstudio &&
        isTRUE(requireNamespace("parallelly", quietly = TRUE)) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore,   workers = n_workers)  # fast on Linux non-RStudio
    } else {
      future::plan(future::multisession, workers = n_workers) # safe everywhere
    }
  }

  # motif-level parallel loop
  results <- future.apply::future_lapply(motifs_union, function(m) {
    base <- sanitize_fname(m)
    f_score <- file.path(out_dir, paste0(base, "_score.csv"))
    f_bound <- file.path(out_dir, paste0(base, "_bound.csv"))
    f_annot <- file.path(out_dir, paste0(base, "_annotation.csv"))

    # optionally skip finished motifs
    if (skip_existing && file.exists(f_score) && file.exists(f_bound) && file.exists(f_annot)) {
      return(tibble::tibble(motif = m, n_peaks = NA_integer_, score = f_score, bound = f_bound, annot = f_annot))
    }

    # 1) load single motif (wide)
    wd <- load_one_motif_wide(root_dir, db_name, sample_ids, m, verbose = FALSE)
    if (!nrow(wd)) {
      # write empty CSVs to mark "done but empty"
      readr::write_csv(tibble::tibble(peak_ID=character()), f_score)
      readr::write_csv(tibble::tibble(peak_ID=character()), f_bound)
      readr::write_csv(tibble::tibble(fp_peak=character(), atac_peak=character(), motifs=character()), f_annot)
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    # 2) filter rows where any *_ATAC_bound == 1
    anyb <- filter_any_bound(wd)
    if (!nrow(anyb)) {
      readr::write_csv(.select_fp_scores(wd)[0, ], f_score)
      readr::write_csv(.select_fp_bounds(wd)[0, ], f_bound)
      readr::write_csv(.select_fp_annots(wd)[0, ], f_annot)
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    # 3) write the three CSVs via your existing helpers (single motif in table)
    #    We inline a minimal version to avoid scanning all motifs again.
    readr::write_csv(.select_fp_scores(anyb), f_score)
    readr::write_csv(.select_fp_bounds(anyb), f_bound)
    readr::write_csv({
      # ensure columns + rename like write_fp_tables_by_motif
      x <- anyb |> dplyr::select(peak_ID, peak_ATAC, TFBS_name) |>
        dplyr::rename(fp_peak = peak_ID, atac_peak = peak_ATAC, motifs = TFBS_name)
      tibble::as_tibble(x)
    }, f_annot)

    # 4) free memory aggressively within worker
    n_out <- nrow(.select_fp_scores(anyb))
    rm(wd, anyb); gc()

    tibble::tibble(motif = m, n_peaks = n_out, score = f_score, bound = f_bound, annot = f_annot)
  }, future.seed = FALSE)

  out <- dplyr::bind_rows(results)
  if (verbose) cli::cli_inform("Done. Wrote {sum(!is.na(out$n_peaks))} motif result(s).")
  out
}

# fp_manifest <- stream_write_overviews_by_motif(
#   root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#   db_name    = "JASPAR2024",
#   out_dir    = file.path("inst","extdata","fp_jaspar2024"),
#   n_workers  = 36,
#   verbose    = TRUE
# )
# readr::write_csv(fp_manifest, "inst/extdata/fp_jaspar2024_manifest.csv")
# fp_manifest <- readr::read_csv("inst/extdata/fp_jaspar2024_manifest.csv")
# fp_manifest <- stream_write_overviews_by_motif(
#   root_dir   = "/data/homes/cy232/cutntag/humanPDAC/Nutrients_basal/TOBIAS_merged_peaks",
#   db_name    = "HOCOMOCOv13",
#   out_dir    = file.path("inst","extdata","fp_hocomocov13"),
#   n_workers  = 36,
#   verbose    = TRUE
# )
# readr::write_csv(fp_manifest, "inst/extdata/fp_hocomocov13_manifest.csv")







# res_collapsed <- collapse_by_peak_fast(res_any_bound)

# readr::write_csv(.select_fp_scores(res_collapsed), "inst/extdata/fp_scores_hocomocov13.csv")
# readr::write_csv(.select_fp_bounds(res_collapsed), "inst/extdata/fp_bounds_hocomocov13.csv")
#
# readr::write_csv(res_collapsed |>
#                    dplyr::select(peak_ID, peak_ATAC, TFBS_name) |>
#                    dplyr::rename(fp_peak = peak_ID, atac_peak = peak_ATAC, motifs = TFBS_name), "inst/extdata/fp_annotation_hocomocov13.csv")


# readr::write_csv(.select_fp_scores(res_collapsed), "inst/extdata/fp_scores_jaspar2024.csv")
# readr::write_csv(.select_fp_bounds(res_collapsed), "inst/extdata/fp_bounds_jaspar2024.csv")
#
# readr::write_csv(res_collapsed |>
#                    dplyr::select(peak_ID, peak_ATAC, TFBS_name) |>
#                    dplyr::rename(fp_peak = peak_ID, atac_peak = peak_ATAC, motifs = TFBS_name), "inst/extdata/fp_annotation_jaspar2024.csv")
#
