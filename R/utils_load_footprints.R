.filter_any_bound <- function(wide_tbl) {
  bd <- grep("_ATAC_bound$", names(wide_tbl), value = TRUE)
  if (!length(bd)) {
    return(tibble::as_tibble(wide_tbl[0, , drop = FALSE]))
  }
  keep <- rowSums(as.data.frame(wide_tbl[bd]), na.rm = TRUE) > 0
  tibble::as_tibble(wide_tbl[keep, , drop = FALSE])
}
.select_fp_scores <- function(x) {
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  score_cols <- grep("_ATAC_score$", names(x), value = TRUE)
  if (length(score_cols) == 0) {
    cli::cli_abort("No columns ending with {.val _ATAC_score} were found.")
  }
  out <- x[, c("peak_ID", score_cols), drop = FALSE]
  # rename: *_ATAC_score -> sample id
  colnames(out) <- sub("_ATAC_score$", "", colnames(out))
  tibble::as_tibble(out)
}
.select_fp_bounds <- function(x) {
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  bound_cols <- grep("_ATAC_bound$", names(x), value = TRUE)
  if (length(bound_cols) == 0) {
    cli::cli_abort("No columns ending with {.val _ATAC_bound} were found.")
  }
  out <- x[, c("peak_ID", bound_cols), drop = FALSE]
  colnames(out) <- sub("_ATAC_bound$", "", colnames(out))
  tibble::as_tibble(out)
}
.select_fp_annots <- function(x) {
  if (!"peak_ID" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ID} column.")
  }
  if (!"peak_ATAC" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field peak_ATAC} column.")
  }
  if (!"TFBS_name" %in% names(x)) {
    cli::cli_abort("Input must contain a {.field TFBS_name} column.")
  }
  out <- x[, c("peak_ID", "peak_ATAC", "TFBS_name"), drop = FALSE]
  out <- tibble::as_tibble(out)
  out <- dplyr::rename(out, fp_peak = peak_ID, atac_peak = peak_ATAC, motifs = TFBS_name)
  tibble::as_tibble(out)
}

# Streaming per-motif writer
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
  } else {
    sample_ids <- all_samples
  }
  if (!is.null(n_samples)) {
    if (!is.numeric(n_samples) || n_samples < 1) cli::cli_abort("{.arg n_samples} must be a positive integer.")
    sample_ids <- utils::head(sample_ids, n_samples)
  }

  if (verbose) cli::cli_inform("Indexing *_overview.txt across {length(sample_ids)} sample(s) (DB = {db_name})...")

  motifs_union <- unique(unlist(lapply(sample_ids, function(sid) {
    db_dir <- file.path(root_dir, sid, db_name)
    if (!dir.exists(db_dir)) {
      return(character(0))
    }
    subdirs <- base::list.dirs(db_dir, recursive = FALSE, full.names = FALSE)
    if (!length(subdirs)) {
      return(character(0))
    }
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
  # Save/restore global data.table threads
  prev_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(prev_threads), add = TRUE)
  data.table::setDTthreads(1L)

  read_one_overview_dt <- function(file, sid) {
    hdr <- data.table::fread(file, nrows = 0L, showProgress = FALSE)
    nm <- names(hdr)

    sc <- nm[grepl(paste0("^", sid, "_ATAC_score$"), nm)]
    if (!length(sc)) sc <- nm[grepl("_ATAC_score$", nm)]
    bd <- nm[grepl(paste0("^", sid, "_ATAC_bound$"), nm)]
    if (!length(bd)) bd <- nm[grepl("_ATAC_bound$", nm)]
    if (length(sc) != 1L || length(bd) != 1L) {
      return(NULL)
    }

    need <- c(
      "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
      "peak_chr", "peak_start", "peak_end", sc, bd
    )
    if (!all(need %in% nm)) {
      return(NULL)
    }

    dt <- data.table::fread(file, select = match(need, nm), showProgress = FALSE)
    data.table::setnames(dt, old = sc, new = "ATAC_score")
    data.table::setnames(dt, old = bd, new = "ATAC_bound")

    dt[, peak_ID := paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end)]
    dt[, peak_ATAC := paste0(peak_chr, ":", peak_start, "-", peak_end)]
    dt[, `:=`(sample_id = sid)]
    dt[, .(sample_id, TFBS_name, peak_ID, peak_ATAC, ATAC_score, ATAC_bound)]
  }

  parts <- vector("list", length(sample_ids))
  kept <- 0L
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
  if (!nrow(dt)) {
    return(tibble::as_tibble(data.frame()))
  }
  dt <- dt[TFBS_name == motif_id]
  dt <- unique(dt, by = c("sample_id", "TFBS_name", "peak_ID"))

  wd <- data.table::dcast(
    dt,
    peak_ID + peak_ATAC + TFBS_name ~ sample_id,
    value.var = c("ATAC_score", "ATAC_bound")
  )

  sc_old <- grep("^ATAC_score_", names(wd), value = TRUE)
  bd_old <- grep("^ATAC_bound_", names(wd), value = TRUE)
  sc_new <- sub("^ATAC_score_", "", sc_old)
  bd_new <- sub("^ATAC_bound_", "", bd_old)
  data.table::setnames(wd, sc_old, paste0(sc_new, "_ATAC_score"))
  data.table::setnames(wd, bd_old, paste0(bd_new, "_ATAC_bound"))

  samples_present <- unique(c(sc_new, bd_new))
  desired <- c(
    "peak_ID", "peak_ATAC", "TFBS_name",
    paste0(samples_present, "_ATAC_score"),
    paste0(samples_present, "_ATAC_bound")
  )
  have <- intersect(desired, names(wd))
  wd <- wd[, ..have]

  tibble::as_tibble(wd)
}

.count_rows_fast <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  # Read a single column only; fast even for large CSVs
  out <- tryCatch({
    dt <- data.table::fread(path, select = 1L, showProgress = FALSE)
    as.integer(nrow(dt))
  }, error = function(e) NA_integer_)
  out
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
#'
#' @export
load_footprints <- function(root_dir,
                            db_name,
                            out_dir = NULL,
                            sample_ids = NULL,
                            n_samples = NULL,
                            motif_ids = NULL,
                            n_motifs = NULL,
                            n_workers = max(1L, parallel::detectCores(logical = TRUE)),
                            set_plan = TRUE,
                            skip_existing = TRUE,
                            verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg data.table} is required.")
  }

  all_samples <- base::list.dirs(root_dir, recursive = FALSE, full.names = FALSE)
  if (!length(all_samples)) cli::cli_abort("No subfolders under {.path {root_dir}}.")
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, all_samples)
    if (!length(sample_ids)) cli::cli_abort("None of the requested samples exist under {.path {root_dir}}.")
  } else {
    sample_ids <- all_samples
  }
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

  # Stronger, portable filename sanitizer
  # sanitize_fname <- function(s) {
  #   out <- gsub('[<>:"/\\\\|?*]', "_", s, perl = TRUE)
  #   out <- gsub("[[:cntrl:]]", "_", out, perl = TRUE)
  #   out <- gsub("_+", "_", out, perl = TRUE)
  #   out <- trimws(out)
  #   out <- sub("[\\.+\\s]+$", "", out, perl = TRUE)
  #   if (nchar(out) > 240L) out <- substr(out, 1L, 240L)
  #   out
  # }

  if (set_plan) {
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

  results <- future.apply::future_lapply(motifs_union, function(m) {
    # base <- sanitize_fname(m)
    base <- m

    f_score <- file.path(out_dir, paste0(base, "_score.csv"))
    f_bound <- file.path(out_dir, paste0(base, "_bound.csv"))
    f_annot <- file.path(out_dir, paste0(base, "_annotation.csv"))

    if (skip_existing && file.exists(f_score) && file.exists(f_bound) && file.exists(f_annot)) {
      n_out <- .count_rows_fast(f_score)
      return(tibble::tibble(
        motif = m,
        n_peaks = n_out,
        score = f_score,
        bound = f_bound,
        annot = f_annot
      ))
    }


    wd <- load_one_motif_wide(root_dir, db_name, sample_ids, m, verbose = FALSE)
    if (!nrow(wd)) {
      readr::write_csv(tibble::tibble(peak_ID = character()), f_score)
      readr::write_csv(tibble::tibble(peak_ID = character()), f_bound)
      readr::write_csv(tibble::tibble(fp_peak = character(), atac_peak = character(), motifs = character()), f_annot)
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    anyb <- .filter_any_bound(wd)
    if (!nrow(anyb)) {
      readr::write_csv(.select_fp_scores(wd)[0, ], f_score)
      readr::write_csv(.select_fp_bounds(wd)[0, ], f_bound)
      readr::write_csv(.select_fp_annots(wd)[0, ], f_annot)
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    readr::write_csv(.select_fp_scores(anyb), f_score)
    readr::write_csv(.select_fp_bounds(anyb), f_bound)
    readr::write_csv(.select_fp_annots(anyb), f_annot)

    n_out <- nrow(.select_fp_scores(anyb))
    rm(wd, anyb)
    gc()

    tibble::tibble(motif = m, n_peaks = n_out, score = f_score, bound = f_bound, annot = f_annot)
  }, future.seed = FALSE)

  out <- dplyr::bind_rows(results)
  if (verbose) cli::cli_inform("Done. Wrote {sum(!is.na(out$n_peaks))} motif result(s).")
  out
}
