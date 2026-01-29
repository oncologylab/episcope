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

  if (verbose) .log_inform("Indexing *_overview.txt across {length(sample_ids)} sample(s) (DB = {db_name})...")

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

  if (verbose) .log_inform("Discovered {length(motifs_union)} motif(s) with overview files.")
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
    nm  <- names(hdr)

    # BEGIN EDIT: robust score/bound column detection (old + new naming)
    sc <- nm[grepl(paste0("^", sid, "_ATAC_score$"), nm)]
    bd <- nm[grepl(paste0("^", sid, "_ATAC_bound$"), nm)]

    # fallback: any *_ATAC_score / *_ATAC_bound
    if (!length(sc)) sc <- nm[grepl("_ATAC_score$", nm)]
    if (!length(bd)) bd <- nm[grepl("_ATAC_bound$", nm)]

    # new TOBIAS / fptools style: SID.*footprints_score / SID.*footprints_bound
    if (!length(sc)) sc <- nm[grepl(paste0("^", sid, ".*footprints_score$"), nm)]
    if (!length(bd)) bd <- nm[grepl(paste0("^", sid, ".*footprints_bound$"), nm)]

    # final generic fallback: any footprints_score / footprints_bound
    if (!length(sc)) sc <- nm[grepl("footprints_score$", nm)]
    if (!length(bd)) bd <- nm[grepl("footprints_bound$", nm)]

    if (length(sc) != 1L || length(bd) != 1L) {
      return(NULL)
    }
    # END EDIT

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

    dt[, peak_ID   := paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end)]
    dt[, peak_ATAC := paste0(peak_chr,  ":", peak_start,  "-", peak_end)]
    dt[, `:=`(sample_id = sid)]
    dt[, .(sample_id, TFBS_name, peak_ID, peak_ATAC, ATAC_score, ATAC_bound)]
  }

  acc <- NULL
  samples_present <- character(0)
  for (j in seq_along(sample_ids)) {
    sid <- sample_ids[[j]]
    f <- file.path(root_dir, sid, db_name, motif_id, paste0(motif_id, "_overview.txt"))
    if (!file.exists(f)) next
    dtj <- read_one_overview_dt(f, sid)
    if (is.null(dtj) || !nrow(dtj)) next

    dtj <- dtj[TFBS_name == motif_id]
    if (!nrow(dtj)) next

    dtj <- unique(dtj, by = c("TFBS_name", "peak_ID"))
    dtj[, sample_id := NULL]

    sc_name <- paste0(sid, "_ATAC_score")
    bd_name <- paste0(sid, "_ATAC_bound")
    data.table::setnames(dtj, c("ATAC_score", "ATAC_bound"), c(sc_name, bd_name))

    if (!sid %in% samples_present) {
      samples_present <- c(samples_present, sid)
    }

    if (is.null(acc)) {
      acc <- dtj
    } else {
      acc <- merge(
        acc,
        dtj,
        by = c("peak_ID", "peak_ATAC", "TFBS_name"),
        all = TRUE,
        sort = FALSE
      )
    }
  }

  if (is.null(acc) || !nrow(acc)) {
    if (verbose) .log_inform("Motif {motif_id}: no per-sample overview rows found - skipping.")
    return(tibble::as_tibble(data.frame()))
  }

  desired <- c(
    "peak_ID", "peak_ATAC", "TFBS_name",
    paste0(samples_present, "_ATAC_score"),
    paste0(samples_present, "_ATAC_bound")
  )
  have <- intersect(desired, names(acc))
  acc <- acc[, ..have]

  tibble::as_tibble(acc)
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
#' @param log_file Optional log file path for per-motif messages (useful with parallel runs).
#' @param memory_check Logical; run a conservative RAM-risk preflight check.
#' @param memory_check_n_motifs Integer; number of motifs to sample for RAM estimate.
#' @param memory_check_factor Numeric; multiplier for in-memory overhead (conservative).
#' @param memory_check_fraction Numeric; fraction of available RAM allowed for estimate.
#' @param ask_on_risk Logical; if TRUE and interactive, prompt on OOM risk.
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
                            n_workers = max(1L, parallel::detectCores(TRUE)),
                            set_plan = TRUE,
                            skip_existing = TRUE,
                            verbose = TRUE,
                            log_file = NULL,
                            memory_check = TRUE,
                            memory_check_n_motifs = 3L,
                            memory_check_factor = 2.5,
                            memory_check_fraction = 0.75,
                            ask_on_risk = TRUE) {
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

  # ---- Optional memory-risk preflight (conservative heuristic)
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
        if (requireNamespace("utils", quietly = TRUE)) {
          lim <- suppressWarnings(utils::memory.limit())
          if (is.finite(lim) && lim > 0) return(as.numeric(lim) * 1024^2)
        }
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
            "Potential OOM risk detected (conservative heuristic).",
            sprintf("Samples=%d motifs; median total size across samples = %s; overhead factor = %.2f.",
                    n_sample, .fmt_bytes(stats::median(sizes)), memory_check_factor),
            sprintf("Estimated peak RAM = %s (per worker ~ %s; workers = %d).",
                    .fmt_bytes(est_total), .fmt_bytes(est_per_worker), as.integer(n_workers)),
            sprintf("Available RAM = %s; safety budget = %s (%.0f%%).",
                    .fmt_bytes(mem_avail), .fmt_bytes(budget), memory_check_fraction * 100),
            "Consider reducing n_workers or using skip_existing = TRUE.",
            sep = " "
          )
          cli::cli_warn(msg)
          if (isTRUE(ask_on_risk) && interactive()) {
            ans <- readline("OOM risk detected. Continue? [y/N]: ")
            if (!nzchar(ans) || !tolower(substr(ans, 1, 1)) %in% "y") {
              cli::cli_abort("Aborted by user due to OOM risk.")
            }
          }
        }
      }
    } else if (verbose) {
      .log_inform("RAM preflight check skipped: available memory could not be determined.")
    }
  }

  # BEGIN EDIT: optional early return from an existing, valid manifest
  manifest_path <- NULL
  if (!is.null(out_dir)) {
    manifest_path <- file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
    if (skip_existing && file.exists(manifest_path)) {
      manifest_ok <- FALSE
      man <- tryCatch(
        readr::read_csv(manifest_path, show_col_types = FALSE),
        error = function(e) NULL
      )
      if (!is.null(man)) {
        required_cols <- c("motif", "n_peaks", "score", "bound", "annot")
        if (all(required_cols %in% names(man))) {
          motifs_in_manifest <- unique(man$motif)
          coverage_ok <- all(motifs_union %in% motifs_in_manifest)
          if (coverage_ok) {
            # Restrict to the motifs requested in this call
            man <- man[man$motif %in% motifs_union, , drop = FALSE]
            if (nrow(man) > 0 &&
                !any(is.na(man$n_peaks)) &&
                all(file.exists(man$score)) &&
                all(file.exists(man$bound)) &&
                all(file.exists(man$annot))) {
              manifest_ok <- TRUE
            }
          }
        }
      }
      if (manifest_ok) {
        if (verbose) {
          .log_inform(
            "Found existing manifest at {.path {manifest_path}}; returning cached summary. Set {.code skip_existing = FALSE} to rebuild."
          )
        }
        attr(man, "from_cache") <- TRUE
        return(man)
      } else if (verbose) {
        .log_inform(
          "Existing manifest at {.path {manifest_path}} is incomplete or out-of-date; rebuilding summary from cached motif CSVs."
        )
      }
    }
  }
  # END EDIT

  if (verbose) {
    .log_inform("Streaming {length(motifs_union)} motif(s) with up to {length(sample_ids)} sample(s) each...")
    .log_inform("Writing outputs to {.path {out_dir}}")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

  ## Optional progressr support (future-friendly)
  have_progressr <- requireNamespace("progressr", quietly = TRUE)
  use_parallel <- isTRUE(n_workers > 1L)

  log_file_use <- log_file
  if (use_parallel && is.null(log_file_use) && !is.null(out_dir)) {
    log_file_use <- file.path(dirname(out_dir), paste0(basename(out_dir), "_load_footprints_log.txt"))
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
  if (use_parallel && !is.null(log_file_use) && verbose) {
    .log_inform("Per-motif logs will be written to {.path {log_file_use}}.")
  }

  process_one_motif <- function(m) {
    # base <- sanitize_fname(m)
    base <- m

    f_score <- file.path(out_dir, paste0(base, "_score.csv"))
    f_bound <- file.path(out_dir, paste0(base, "_bound.csv"))
    f_annot <- file.path(out_dir, paste0(base, "_annotation.csv"))

    if (!have_progressr && verbose && is.null(log_file_use)) {
      .log_inform("Processing motif {m} ...")
    }

    if (skip_existing && file.exists(f_score) && file.exists(f_bound) && file.exists(f_annot)) {
      n_out <- .count_rows_fast(f_score)
      if (!is.null(log_file_use)) {
        log_line(sprintf("Motif %s: outputs exist; skipped (n_peaks=%s).", m, n_out))
      } else if (verbose) {
        .log_inform(
          "Motif {m}: all three output files already exist ({.path {basename(f_score)}}, {.path {basename(f_bound)}}, {.path {basename(f_annot)}}); skipping. Set {.code skip_existing = FALSE} to overwrite."
        )
      }
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
      if (!is.null(log_file_use)) {
        log_line(sprintf("Motif %s: no peaks found; wrote empty tables.", m))
      } else if (verbose) {
        .log_inform("Motif {m}: no peaks found; wrote empty tables.")
      }
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    anyb <- .filter_any_bound(wd)
    if (!nrow(anyb)) {
      readr::write_csv(.select_fp_scores(wd)[0, ], f_score)
      readr::write_csv(.select_fp_bounds(wd)[0, ], f_bound)
      readr::write_csv(.select_fp_annots(wd)[0, ], f_annot)
      if (!is.null(log_file_use)) {
        log_line(sprintf("Motif %s: no bound peaks after filtering; wrote empty tables.", m))
      } else if (verbose) {
        .log_inform("Motif {m}: no bound peaks after filtering; wrote empty tables.")
      }
      return(tibble::tibble(motif = m, n_peaks = 0L, score = f_score, bound = f_bound, annot = f_annot))
    }

    score_tbl <- .select_fp_scores(anyb)
    bound_tbl <- .select_fp_bounds(anyb)
    annot_tbl <- .select_fp_annots(anyb)

    readr::write_csv(score_tbl, f_score)
    readr::write_csv(bound_tbl, f_bound)
    readr::write_csv(annot_tbl, f_annot)

    n_out <- nrow(score_tbl)
    if (!is.null(log_file_use)) {
      log_line(sprintf("Motif %s: wrote %d peak(s) to CSVs.", m, n_out))
    } else if (verbose) {
      .log_inform("Motif {m}: wrote {n_out} peak(s) to CSVs.")
    }

    rm(wd, anyb, score_tbl, bound_tbl, annot_tbl)
    gc()

    tibble::tibble(motif = m, n_peaks = n_out, score = f_score, bound = f_bound, annot = f_annot)
  }

  if (have_progressr) {
    results <- progressr::with_progress({
      p <- progressr::progressor(along = motifs_union)
      future.apply::future_lapply(motifs_union, function(m) {
        if (verbose) {
          p(message = sprintf("Processing motif %s", m))
        } else {
          p()
        }
        process_one_motif(m)
      }, future.seed = FALSE)
    })
  } else {
    results <- future.apply::future_lapply(motifs_union, process_one_motif, future.seed = FALSE)
  }

  out <- dplyr::bind_rows(results)

  # BEGIN EDIT: always (re)write manifest inferred from out_dir
  if (!is.null(out_dir)) {
    manifest_path <- file.path(dirname(out_dir), paste0(basename(out_dir), "_manifest.csv"))
    readr::write_csv(out, manifest_path)
    if (verbose) {
      .log_inform(
        "Saved manifest with {nrow(out)} motif(s) to {.path {manifest_path}}."
      )
    }
  }
  attr(out, "from_cache") <- FALSE
  # END EDIT

  if (verbose) .log_inform("Done. Wrote {sum(!is.na(out$n_peaks))} motif result(s).")
  out
}
