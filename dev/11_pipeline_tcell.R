library(episcope)
load_config("episcope_grn.yaml")
# source("R/utils_connect_tf_enhancers_to_target_genes.R")
# source("R/utils_load_footprints_process.R")

# ──────────────────────────────────────────────────────────────────────────────
# Turn modules ON/OFF
# ──────────────────────────────────────────────────────────────────────────────
do_load_footprints_preprocess    <- TRUE
do_tf_binding_sites_prediction   <- TRUE
do_tf_to_target_genes_prediction <- FALSE
do_build_grn                     <- FALSE
do_differential_grn              <- FALSE
do_topic_subplot_single          <- FALSE
do_topic_subplot_bulk            <- FALSE
do_tf_centric_subplot_single     <- FALSE
do_tf_centric_subplot_bulk       <- FALSE

base_dir <- "/data/homes/yl814/episcope_test/GSE192390_T_Cells"

# Step 0. Load footprint data and preprocess------------------------------------
## load_footprints
## fp_manifest_trim
## fp_manifest_trim_annots
## align_footprints
## qn_footprints
## save_footprints
if (do_load_footprints_preprocess == TRUE) {
  sample_metadata <- readr::read_tsv(file.path(base_dir, "GSE192390_CD8_T_Cells_samples.txt"), na = "NA")
  sample_metadata$id <- sample_metadata$ID

  fp_manifest <- load_footprints(
    root_dir   = "/data/homes/yl814/tbio_cy232/GSE192390_TOBIAS",
    sample_ids = sample_metadata$ID,
    db_name    = "JASPAR2024",
    out_dir    = file.path(base_dir, "cache", sprintf("fp_%s", db)),
    n_workers  = 20
  )


  readr::write_csv(fp_manifest, file.path(base_dir, "cache", sprintf("fp_%s_manifest.csv", db)))
  fp_manifest <- readr::read_csv(file.path(base_dir, "cache", sprintf("fp_%s_manifest.csv", db)))

  # (Optional, when using HOCOMOCO database)
  if (db == "hocomocov13"){
    fp_manifest <- fp_manifest_trim(fp_manifest) # rename files on disk by default
    # Overwrite every annotation CSV referenced by the manifest:
    summary_tbl <- fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = TRUE)

    # Inspect what changed:
    dplyr::count(summary_tbl, status)
    sum(summary_tbl$n_fixed, na.rm = TRUE)
  }

  options(future.globals.maxSize = 32 * 1024^3)
  # Align the peaks based on the peak similarity
  fp_aligned <- align_footprints(fp_manifest,
                                 mid_slop        = 10L, # midpoint tolerance (bp)
                                 round_digits    = 1L, # round scores before comparing vectors
                                 score_match_pct = 0.8) # fraction of samples that must match (<=1, >=0)

  readr::write_csv(fp_aligned$fp_bound, file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_bounds_%s.csv", db)))
  readr::write_csv(fp_aligned$fp_score, file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_scores_%s.csv", db)))
  readr::write_csv(fp_aligned$fp_annotation, file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_annotation_%s.csv", db)))


  fp_score_raw <- fp_aligned$fp_score
  fp_aligned_normalized <- fp_aligned
  fp_aligned_normalized$fp_score <- qn_footprints(fp_aligned_normalized$fp_score, id_col = "peak_ID")
  readr::write_csv(fp_aligned_normalized$fp_score, file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_scores_qn_%s.csv", db)))

  fp_aligned_normalized_manifest <- save_footprints(fp_aligned_normalized, out_dir = file.path(base_dir, "cache", sprintf("fp_aligned_normalized_%s", db)))
  readr::write_csv(fp_aligned_normalized_manifest, file.path(base_dir, "cache", sprintf("fp_aligned_normalized_manifest_%s.csv", db)))
  fp_aligned_normalized_manifest <- readr::read_csv(file.path(base_dir, "cache", sprintf("fp_aligned_normalized_manifest_%s.csv", db)))

  if (db == "jaspar2024") {
    motif_db <- readr::read_tsv("inst/extdata/genome/JASPAR2024_mm.txt") |> dplyr::filter(!is.na(HGNC)) # mouse
  } else if (db == "hocomocov13") {
    motif_db <- readr::read_tsv("inst/extdata/genome/HOCOMOCOv13.txt")
  }

  tf_list <- motif_db |>
    tidyr::separate_rows(HGNC, sep = "::") |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::distinct(HGNC) |>
    dplyr::pull(HGNC)

  # Load RNA data and ATAC data
  atac_data <- readr::read_tsv(file.path(base_dir, "GSE192390_CD8_ATAC.master_table.narrowpeaks.10mil.txt"))

  atac_out <- load_atac(atac_data, sort_peaks = TRUE)
  atac_score <- atac_out$score
  atac_overlap <- atac_out$overlap

  rna <- readr::read_tsv(file.path(base_dir, "GSE192389_TCell_InVitro_InVivo_Normalized_Counts.tsv"))
  rna <- clean_hgnc(rna) # clean "HGNC" column
  # Filter rna expression genes/tfs needs to reach the threshold in at least 1 sample (group size)
  rna <- filter_rna_expr(rna, tf_list, hgnc_col = "HGNC", gene_min = threshold_gene_expr, tf_min = threshold_tf_expr, min_samples = 1L)

  # Build RNA and ATAC data list
  rna_atac_build_args <- list(
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = rna,
    metadata      = sample_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db, # tibble
    label_col     = "Sample",
    expected_n    = 20
  )

  # ATAC peak and/or RNA expression based, footprint filtering
  fp_aligned_normalized_filtered_manifest <- filter_footprints(
    fp_manifest         = fp_aligned_normalized_manifest,
    out_dir             = file.path(base_dir, "cache", sprintf("fp_aligned_normalized_filtered_%s", db)),
    build_args          = rna_atac_build_args,
    n_workers           = 10,
    skip_existing       = TRUE,
    tf_expr_filter      = FALSE,
    verbose             = TRUE
  )

  readr::write_csv(fp_aligned_normalized_filtered_manifest, file.path(base_dir, "cache", sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))
  fp_aligned_normalized_filtered_manifest <-  readr::read_csv(file.path(base_dir, "cache", sprintf("fp_aligned_normalized_filtered_manifest_%s.csv", db)))
}


# Step 1. Predict TF binding sites ----------------------------------------

if (do_tf_binding_sites_prediction == TRUE) {
  # correlate with all TFBS
  tfs <- tf_list[tf_list %in% unique(rna$HGNC)]

  # Build a *named list* of TF expression vectors
  tf_expr_list <- lapply(tfs, function(tf) {
    rna |>
      dplyr::filter(HGNC == tf) |>
      dplyr::select(-ensembl_gene_id, -HGNC) |>
      simplify2array()
  })
  names(tf_expr_list) <- tfs

  progressr::handlers("progress")

  res_list <- tf_corr_footprints_all_tfbs_multi(
    fp_manifest = fp_aligned_normalized_filtered_manifest,
    tf_expr_list = tf_expr_list,
    out_dir   = file.path(base_dir,
                          "predict_tf_binding_sites",
                          "binding_probability_overviews"),
    motif_db  = motif_db,
    cor_method = "pearson",
    min_non_na = 5L,
    r_thr      = 0.30,
    fdr_thr    = 0.05,
    threads    = 4L,            # data.table threads per worker
    write_bed  = FALSE,         # TRUE for *_all/_bound/_unbound.bed
    n_workers  = 36L,            # number of TFs processed in parallel
    set_plan   = TRUE,
    verbose    = TRUE
  )



  # TF expr to footprints
  # fp_corr_manifest <- tf_corr_footprints(
  #   fp_manifest = fp_aligned_normalized_filtered_manifest,
  #   out_dir     = file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_%s", db)),
  #   build_args  = rna_atac_build_args,
  #   n_workers   = 15,
  #   cor_method  = "pearson",
  #   min_non_na  = 5L
  # )


  # readr::write_csv(fp_corr_manifest, file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))
  # fp_corr_manifest <- readr::read_csv(file.path(base_dir, sprintf("fp_aligned_normalized_filtered_corr_manifest_%s.csv", db)))

  # tf_corr_qc(fp_corr_manifest)
  # threshold_fp_tf_corr_p <- 0.3
  # combined <- tf_corr_footprints_filter(fp_corr_manifest, p_thr = threshold_fp_tf_corr_p, r_thr = threshold_fp_tf_corr_r, output_bed = file.path(base_dir, sprintf("fp_predicted_tfbs_%s", db)))

  # fp_score <- combined$fp_score
  # fp_bound <- combined$fp_bound
  # fp_annotation <- combined$fp_annotation

  # readr::write_csv(fp_score, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_score_strict_tf_filtered_corr_", db, ".csv"))
  # readr::write_csv(fp_bound, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_bound_strict_tf_filtered_corr_", db, ".csv"))
  # readr::write_csv(fp_annotation, paste0("/data/homes/yl814/episcope_test/GSE87218_ATAC/fp_annotation_strict_tf_filtered_corr_", db, ".csv"))

  # unique(fp_annotation$fp_peak)

  # fp_score <- readr::read_csv(file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db)))
  # fp_bound <- readr::read_csv(file.path(base_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db)))
  # fp_annotation <- readr::read_csv(file.path(base_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db)))


  # grn_set
  # grn_set <- build_grn_set(
  #   fp_score      = fp_score,
  #   fp_bound      = fp_bound,
  #   fp_annotation = fp_annotation,
  #   atac_score    = atac_score,
  #   atac_overlap  = atac_overlap,
  #   rna           = rna,
  #   metadata      = sample_metadata,
  #   tf_list       = tf_list,
  #   motif_db      = motif_db,
  #   label_col     = "Sample",
  #   expected_n    = 13
  # )
  # grn_set <- filter_by_min_bound(grn_set, min_bound = 1L)
}

# (Optional correlate with all tfbs, not only canonical tfbs)
# tf_vec <- strict_rna |> dplyr::filter(HGNC == "ZEB1") |> dplyr::select(-ensembl_gene_id, -HGNC) |> simplify2array()
# res <- tf_corr_footprints_all_tfbs(
#   fp_manifest = fp_aligned_normalized_filtered_manifest,
#   tf_name  = "ZEB1",
#   tf_expr  = tf_vec,
#   out_dir  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
# )
# res$files


#' Collect per-TF peak_ID sets from *_overview.txt files
#'
#' Parallel loader for per-TF overview files produced by
#' `tf_corr_footprints_all_tfbs()` / `_multi()`.
#'
#' For each overview file, we:
#'   * read the overview table
#'   * (optionally) mark canonical rows if `canonical_only = TRUE`
#'   * filter by correlation thresholds
#'   * construct `peak_ID = TFBS_chr:TFBS_start-TFBS_end`
#'   * return a named list: one character vector of peak_ID per TF.
#'
#' Expected columns in each overview file:
#'   * TFBS_chr, TFBS_start, TFBS_end
#'   * corr_fp_tf_r, corr_fp_tf_p
#'   * TF
#'   * optionally: corr_fp_tf_p_adj, canonical (0/1), TFBS_canonical (semicolon-separated)
#'
#' @param overview_dir Directory containing per-TF `*_overview.txt` files.
#' @param threshold_fp_tf_corr_r Numeric, minimum r (default 0.3).
#' @param threshold_fp_tf_corr_p Numeric, maximum p (default 0.05).
#' @param pattern Regex to detect overview files. Default: `"_overview\\.txt$"`.
#' @param canonical_only Logical; if `TRUE`, keep only rows that are
#'   canonical for that TF (see Details). Default: `FALSE`.
#'   Canonical rows are defined as:
#'   * If column `canonical` exists: `canonical == 1`.
#'   * Else if `TFBS_canonical` exists: rows where the TF symbol appears
#'     as a standalone token or part of a dimer token (e.g. Ahr in Ahr::Arnt).
#'   * Else: all rows are treated as non-canonical and a warning is emitted.
#' @param use_p_adj Logical; if `TRUE` (default), use the adjusted p-value
#'   column `corr_fp_tf_p_adj` when available, otherwise fall back to the
#'   raw p-value column `corr_fp_tf_p`. If `FALSE`, always use raw p.
#' @param workers Integer number of parallel workers. If `NULL` (default),
#'   uses `max(1, floor(parallel::detectCores() / 2))`.
#' @param verbose Logical; if `TRUE` (default), print progress via cli.
#'
#' @return Named list `tf_peak_ids`, where each element is a character
#'   vector of `peak_ID` strings for that TF.
#' @export
collect_tf_peak_ids_by_corr <- function(overview_dir,
                                        threshold_fp_tf_corr_r = 0.30,
                                        threshold_fp_tf_corr_p = 0.05,
                                        pattern                 = "_overview\\.txt$",
                                        canonical_only          = FALSE,
                                        use_p_adj               = TRUE,
                                        workers                 = NULL,
                                        verbose                 = TRUE) {
  # ---------------------------------------------------------------------------
  # Basic checks & file discovery
  # ---------------------------------------------------------------------------
  if (missing(overview_dir) || is.null(overview_dir)) {
    cli::cli_abort("{.arg overview_dir} must be a directory containing *_overview.txt files.")
  }

  if (!dir.exists(overview_dir)) {
    cli::cli_abort("Directory {.path {overview_dir}} does not exist.")
  }

  overview_dir <- normalizePath(overview_dir, mustWork = TRUE)

  overview_files <- list.files(
    overview_dir,
    pattern    = pattern,
    full.names = TRUE
  )

  if (!length(overview_files)) {
    cli::cli_warn(
      c(
        "No overview files found in {.path {overview_dir}}.",
        "i" = "Expected files matching pattern {.code {pattern}}."
      )
    )
    return(structure(list(), names = character(0)))
  }

  if (verbose) {
    cli::cli_inform(
      "Found {.val {length(overview_files)}} overview file{?s} in {.path {overview_dir}}."
    )
  }

  # ---------------------------------------------------------------------------
  # Resolve workers
  # ---------------------------------------------------------------------------
  if (is.null(workers)) {
    n_cores <- parallel::detectCores()
    if (is.na(n_cores) || n_cores < 1L) {
      workers <- 1L
    } else {
      workers <- max(1L, floor(n_cores / 2L))
    }
  } else {
    workers <- as.integer(workers)
    if (is.na(workers) || workers < 1L) {
      cli::cli_abort("{.arg workers} must be a positive integer.")
    }
  }

  # ---------------------------------------------------------------------------
  # Helper: process ONE overview file
  # ---------------------------------------------------------------------------
  process_one_overview <- function(path,
                                   canonical_only_local = FALSE,
                                   use_p_adj_local      = TRUE,
                                   verbose_local        = FALSE) {
    if (!file.exists(path)) {
      cli::cli_warn("Overview file not found: {.path {path}}")
      return(list(tf = NA_character_, peaks = character(0)))
    }

    dat <- readr::read_tsv(
      path,
      col_types      = readr::cols(.default = "c"),
      show_col_types = FALSE
    )

    req_cols <- c("TFBS_chr", "TFBS_start", "TFBS_end",
                  "corr_fp_tf_r", "corr_fp_tf_p", "TF")
    miss     <- setdiff(req_cols, colnames(dat))
    if (length(miss) > 0L) {
      cli::cli_warn(
        c(
          "Skipping overview file with missing required columns.",
          "file" = path,
          "x"    = sprintf("Missing: %s", paste(miss, collapse = ", "))
        )
      )
      return(list(tf = NA_character_, peaks = character(0)))
    }

    # TF name -----------------------------------------------------------------
    tf_vals <- unique(as.character(dat[["TF"]]))
    tf_vals <- tf_vals[!is.na(tf_vals) & tf_vals != ""]
    if (!length(tf_vals)) {
      cli::cli_warn(
        c(
          "No valid TF name found in overview file.",
          "file" = path
        )
      )
      return(list(tf = NA_character_, peaks = character(0)))
    }
    if (length(tf_vals) > 1L) {
      cli::cli_warn(
        c(
          "Overview file has multiple TF values; using the first.",
          "file" = path,
          "i"    = sprintf("TF values: %s", paste(head(tf_vals, 5L), collapse = ", "))
        )
      )
    }
    tf_name <- tf_vals[[1]]

    if (verbose_local) {
      cli::cli_inform("Reading TF {.val {tf_name}} from {.path {basename(path)}}")
    }

    # Numeric r / p -----------------------------------------------------------
    r_num <- suppressWarnings(as.numeric(dat[["corr_fp_tf_r"]]))

    # choose p column: adjusted first if requested and available
    if (use_p_adj_local && "corr_fp_tf_p_adj" %in% colnames(dat)) {
      p_num <- suppressWarnings(as.numeric(dat[["corr_fp_tf_p_adj"]]))
      p_col_used <- "corr_fp_tf_p_adj"
    } else {
      if (use_p_adj_local && !"corr_fp_tf_p_adj" %in% colnames(dat)) {
        cli::cli_warn(
          c(
            "Requested use_p_adj = TRUE but no {.field corr_fp_tf_p_adj} column.",
            "file" = path,
            "i"    = "Falling back to raw {.field corr_fp_tf_p}."
          )
        )
      }
      p_num <- suppressWarnings(as.numeric(dat[["corr_fp_tf_p"]]))
      p_col_used <- "corr_fp_tf_p"
    }

    ok_rp <- !is.na(r_num) & !is.na(p_num)

    if (!any(ok_rp)) {
      cli::cli_warn(
        c(
          "No usable r/p values for TF {.val {tf_name}}.",
          "file" = path
        )
      )
      return(list(tf = tf_name, peaks = character(0)))
    }

    # Canonical flag (if requested) ------------------------------------------
    canonical_flag <- rep(TRUE, length(r_num))

    if (canonical_only_local) {
      if ("canonical" %in% colnames(dat)) {
        canon_vec <- suppressWarnings(as.integer(dat[["canonical"]]))
        canonical_flag <- !is.na(canon_vec) & (canon_vec == 1L)
      } else if ("TFBS_canonical" %in% colnames(dat)) {
        tf_bs_canon <- dat[["TFBS_canonical"]]
        tf_bs_canon[is.na(tf_bs_canon)] <- ""

        # Escape TF for regex
        tf_esc <- gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", tf_name)

        # Match TF as a token or part of a dimer token
        pattern <- paste0(
          "(^|;)", tf_esc, "(;|$|::)",      # token equals TF or TF followed by "::"
          "|(::", tf_esc, ")(;|$)"          # or token ends with "::TF"
        )

        canonical_flag <- grepl(pattern, tf_bs_canon, perl = TRUE)
      } else {
        cli::cli_warn(
          c(
            "canonical_only = TRUE but no {.field canonical} or {.field TFBS_canonical} column.",
            "file" = path,
            "i"    = "Treating all rows as non-canonical; this TF will contribute 0 peaks."
          )
        )
        canonical_flag <- rep(FALSE, length(r_num))
      }
    }

    # Final keep mask ---------------------------------------------------------
    keep <- ok_rp &
      (r_num >= threshold_fp_tf_corr_r) &
      (p_num <= threshold_fp_tf_corr_p) &
      canonical_flag

    if (!any(keep)) {
      if (verbose_local) {
        cli::cli_inform(
          "TF {.val {tf_name}}: no peaks passed r >= {threshold_fp_tf_corr_r}, {p_col_used} <= {threshold_fp_tf_corr_p}{if (canonical_only_local) ' and canonical == 1' else ''}."
        )
      }
      return(list(tf = tf_name, peaks = character(0)))
    }

    dat_sub <- dat[keep, , drop = FALSE]

    # Construct peak_ID = TFBS_chr:TFBS_start-TFBS_end -----------------------
    chr   <- as.character(dat_sub[["TFBS_chr"]])
    start <- suppressWarnings(as.integer(dat_sub[["TFBS_start"]]))
    end   <- suppressWarnings(as.integer(dat_sub[["TFBS_end"]]))

    ok_pos <- !is.na(chr) & nzchar(chr) & !is.na(start) & !is.na(end)

    if (!any(ok_pos)) {
      cli::cli_warn(
        c(
          "TF {.val {tf_name}}: no valid genomic coordinates in kept rows.",
          "file" = path
        )
      )
      return(list(tf = tf_name, peaks = character(0)))
    }

    chr   <- chr[ok_pos]
    start <- start[ok_pos]
    end   <- end[ok_pos]

    peak_ids <- unique(paste0(chr, ":", start, "-", end))

    list(tf = tf_name, peaks = peak_ids)
  }

  # ---------------------------------------------------------------------------
  # Parallel over overview files
  # ---------------------------------------------------------------------------
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future::multisession, workers = workers)

  if (verbose) {
    cli::cli_inform(
      "Collecting peak_IDs with {.val {workers}} worker{?s} (canonical_only = {.val {canonical_only}}, use_p_adj = {.val {use_p_adj}})."
    )
  }

  idxs <- seq_along(overview_files)

  res_list <- future.apply::future_lapply(
    idxs,
    function(i) {
      process_one_overview(
        path                 = overview_files[[i]],
        canonical_only_local = canonical_only,
        use_p_adj_local      = use_p_adj,
        verbose_local        = verbose
      )
    },
    future.seed = TRUE
  )

  # ---------------------------------------------------------------------------
  # Assemble named list
  # ---------------------------------------------------------------------------
  tf_names  <- vapply(res_list, function(x) x$tf,   FUN.VALUE = character(1L))
  peak_vecs <- lapply(res_list, function(x) x$peaks)

  # Drop any NA TFs
  keep_idx  <- !is.na(tf_names) & nzchar(tf_names)
  tf_names  <- tf_names[keep_idx]
  peak_vecs <- peak_vecs[keep_idx]

  # If duplicate TFs somehow appear, merge their peak sets --------------------
  if (length(tf_names)) {
    uniq_tfs <- unique(tf_names)
    out <- vector("list", length(uniq_tfs))
    names(out) <- uniq_tfs
    for (tf in uniq_tfs) {
      idx <- which(tf_names == tf)
      out[[tf]] <- unique(unlist(peak_vecs[idx], use.names = FALSE))
    }
  } else {
    out <- list()
  }

  if (verbose) {
    cli::cli_inform(
      "Collected peak_ID sets for {.val {length(out)}} TF{?s}."
    )
  }

  out
}



#' Subset footprint tables by TF–peak correlations
#'
#' Given full footprint score / bound / annotation tables and a list of
#' per-TF peak_ID vectors (e.g. from collect_tf_peak_ids_by_corr()),
#' keep only rows whose peak IDs are in the union across all TFs.
#'
#' @param fp_score      Tibble/data.frame with a peak column named
#'                      either "peak_ID" or "fp_peak".
#' @param fp_bound      Tibble/data.frame with a peak column named
#'                      either "peak_ID" or "fp_peak".
#' @param fp_annotation Tibble/data.frame with a peak column named
#'                      either "fp_peak" or "peak_ID".
#' @param tf_peak_ids   Named list of character vectors; each element
#'                      is a set of peak_ID strings for one TF.
#' @param verbose       Logical; if TRUE, print a brief summary.
#'
#' @return A list with elements:
#'   - fp_score      : subset of fp_score
#'   - fp_bound      : subset of fp_bound
#'   - fp_annotation : subset of fp_annotation
subset_fp_by_tf_corr <- function(fp_score,
                                 fp_bound,
                                 fp_annotation,
                                 tf_peak_ids,
                                 verbose = TRUE) {
  # Union of all TF peak_IDs -----------------------------------------------
  if (is.null(tf_peak_ids) || !length(tf_peak_ids)) {
    cli::cli_abort("`tf_peak_ids` must be a non-empty list of peak_ID vectors.")
  }

  peak_union <- unique(unlist(tf_peak_ids, use.names = FALSE))
  peak_union <- peak_union[!is.na(peak_union) & nzchar(peak_union)]

  if (!length(peak_union)) {
    cli::cli_abort("No valid peak IDs found in `tf_peak_ids`.")
  }

  # Helper: pick ID column name -------------------------------------------
  pick_id_col <- function(x, score = FALSE, annot = FALSE) {
    nm <- names(x)

    if (score || !annot) {
      # For score/bound: usually peak_ID or fp_peak
      if ("peak_ID" %in% nm) return("peak_ID")
      if ("fp_peak"  %in% nm) return("fp_peak")
    }

    if (annot) {
      if ("fp_peak"  %in% nm) return("fp_peak")
      if ("peak_ID"  %in% nm) return("peak_ID")
    }

    cli::cli_abort(
      "Could not find a peak ID column in table with columns: {cli::fmt_cols(nm)}"
    )
  }

  score_id_col <- pick_id_col(fp_score, score = TRUE, annot = FALSE)
  bound_id_col <- pick_id_col(fp_bound, score = FALSE, annot = FALSE)
  annot_id_col <- pick_id_col(fp_annotation, score = FALSE, annot = TRUE)

  # Subset each table ------------------------------------------------------
  fp_score_sub <- fp_score[fp_score[[score_id_col]] %in% peak_union, , drop = FALSE]
  fp_bound_sub <- fp_bound[fp_bound[[bound_id_col]] %in% peak_union, , drop = FALSE]
  fp_annot_sub <- fp_annotation[fp_annotation[[annot_id_col]] %in% peak_union, , drop = FALSE]

  if (verbose) {
    cli::cli_inform(c(
      "Subsetted footprint tables by TF-correlated peaks:",
      "  Total unique peak_IDs in union: {length(peak_union)}",
      "  fp_score: {nrow(fp_score_sub)} / {nrow(fp_score)} rows kept",
      "  fp_bound: {nrow(fp_bound_sub)} / {nrow(fp_bound)} rows kept",
      "  fp_annotation: {nrow(fp_annot_sub)} / {nrow(fp_annotation)} rows kept"
    ))
  }

  list(
    fp_score      = fp_score_sub,
    fp_bound      = fp_bound_sub,
    fp_annotation = fp_annot_sub
  )
}




# Step 2. Connect TF-occupied enhancers to target genes -------------------
if (do_tf_to_target_genes_prediction == TRUE) {

  overview_dir <- file.path(
    base_dir,
    "predict_tf_binding_sites",
    "binding_probability_overviews"
  )



  tf_peak_ids_full <- collect_tf_peak_ids_by_corr(
    overview_dir              = file.path(base_dir, "predict_tf_binding_sites", "binding_probability_overviews"),
    threshold_fp_tf_corr_r    = threshold_fp_tf_corr_r,
    threshold_fp_tf_corr_p    = threshold_fp_tf_corr_p,
    canonical_only            = FALSE,
    use_p_adj                 = TRUE,
    workers                   = 36,
    verbose                   = TRUE
  )

  tf_peak_ids_canonical <- collect_tf_peak_ids_by_corr(
    overview_dir              = file.path(base_dir, "predict_tf_binding_sites", "binding_probability_overviews"),
    threshold_fp_tf_corr_r    = threshold_fp_tf_corr_r,
    threshold_fp_tf_corr_p    = threshold_fp_tf_corr_p,
    canonical_only            = TRUE,
    use_p_adj                 = TRUE,
    workers                   = 36,
    verbose                   = TRUE
  )




  fp_score      <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_scores_qn_%s.csv", db)))
  fp_bound      <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_bounds_%s.csv", db)))
  fp_annotation <- readr::read_csv(file.path(base_dir, "predict_tf_binding_sites", sprintf("fp_annotation_%s.csv", db)))


  # footprints significant in at least one TF
  subset_res <- subset_fp_by_tf_corr(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annotation,
    tf_peak_ids   = tf_peak_ids_canonical,
    verbose       = TRUE
  )

  fp_score_sub      <- subset_res$fp_score
  fp_bound_sub      <- subset_res$fp_bound
  fp_annotation_sub <- subset_res$fp_annotation


  # pivot longer the fp_annotation by motifs -> tfs -> tf::tf dimers

  # grn_set
  grn_set <- build_grn_set(
    fp_score      = subset_res$fp_score,
    fp_bound      = subset_res$fp_bound,
    fp_annotation = subset_res$fp_annotation,
    atac_score    = atac_score,
    atac_overlap  = atac_overlap,
    rna           = rna,
    metadata      = sample_metadata,
    tf_list       = tf_list,
    motif_db      = motif_db,
    label_col     = "Sample",
    expected_n    = 20
  )
  grn_set <- filter_by_min_bound(grn_set, min_bound = 1L)

  # gene_annot_ref_hg38 <- episcope_build_gene_annot("hg38")
  gene_annot_ref_mm10 <- episcope_build_gene_annot("mm10")

  # 30kb/25kb/20kb based
  gh_std <- episcope_make_windowed_gh(
    peaks      = grn_set$atac_score,
    genes      = unique(c(grn_set$rna$HGNC, grn_set$rna$ensembl_gene_id)),
    flank_bp   = 100000, # 100kb
    mode       = "TSS",
    gene_annot = gene_annot_ref_mm10,
    id_col     = "HGNC" # or "ensembl_gene_id"
  )

  # atac_peak to target gene corr




  GeneHancer_AnnotSV_tissues_v5.24 <- readr::read_tsv("/data/homes/yl814/genome/GeneHancer_v5.24/GeneHancer_AnnotSV_tissues_v5.24.txt")
  unique(GeneHancer_AnnotSV_tissues_v5.24$tissue)


  readr::write_csv(GeneHancer_AnnotSV_tissues_v5.24 |> dplyr::distinct(tissue), "/data/homes/yl814/genome/GeneHancer_v5.24/GeneHancer_tissues_cleaned_v5.24.csv")
  # GeneHancer based
  # gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_full.csv")) # GeneHancer_v5.24_elite_panc.csv
  gh_std   <- load_genehancer_panc(file.path("inst","extdata","GeneHancer_v5.24_full.csv"))
  options(future.globals.maxSize = 64 * 1024^3)

  # Run atac gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
  # threshold_atac_gene_corr_p <- 0.3
  # atac_res <- correlate_atac_to_genes(grn_set = grn_set, gh_tbl = gh_std, gene_mode = "both", fdr = threshold_atac_gene_corr_p, r_abs_min = threshold_atac_gene_corr_abs_r)

  # shrink to only correlated peaks
  # grn_set_filtered <- filter_grn_by_corr(grn_set, atac_res$atac_gene_corr_kept)

  # Run fp gene correlations (adjusted p-value < 0.05, |r| >= 0.3)
  # threshold_fp_gene_corr_p <- 0.3
  # fp_res <- correlate_fp_to_genes(
  #   grn_set              = grn_set_filtered,
  #   atac_gene_corr_kept  = atac_res$atac_gene_corr_kept,
  #   fdr                  = threshold_fp_gene_corr_p,
  #   r_abs_min            = threshold_fp_gene_corr_abs_r,
  #   workers              = 15
  # )

  atac_res <- correlate_atac_to_genes(#TODO Update the function name as correlate_atac_tf_to_genes
    grn_set  = grn_set,
    gh_tbl   = gh_std,
    gene_mode = "both",
    fdr       = threshold_atac_gene_corr_p,
    r_abs_min = threshold_atac_gene_corr_abs_r,
    cache_dir = file.path(base_dir, "cache"),
    cache_tag = "atac_to_gene_links_100kb",
    workers   = 36,
    cache_verbose = TRUE
  )





  fp_res <- correlate_fp_to_genes(#TODO Update the function name as correlate_fp_tf_to_genes
    grn_set             = grn_set,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,  # use `atac_gene_corr_full` to not filter by any atac peaks correlation
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 36,
    cache_dir           = file.path(base_dir, "cache"),
    cache_tag           = "fp_to_gene_links_100kb",
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )

  rna_res <- correlate_rna_tf_to_genes(
    grn_set             = grn_set,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,  # use `atac_gene_corr_full` to not filter by any atac peaks correlation
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 36,
    cache_dir           = file.path(base_dir, "cache"),
    cache_tag           = "rna_to_gene_links_100kb",
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )

  lmm_res <- lmm_tf_to_genes(
    grn_set             = grn_set,
    atac_gene_corr_kept = atac_res$atac_gene_corr_full,  # use `atac_gene_corr_full` to not filter by any atac peaks correlation
    fdr                 = threshold_fp_gene_corr_p,
    r_abs_min           = threshold_fp_gene_corr_abs_r,
    method              = "pearson",
    workers             = 36,
    cache_dir           = file.path(base_dir, "cache"),
    cache_tag           = "lmm_tf_to_gene_links_100kb",
    cache_chunk_size    = 5000L,
    cache_verbose       = TRUE
  )


  readr::write_csv(atac_res$atac_gene_corr_full, file.path(base_dir, sprintf("atac_gene_corr_full_100kb_%s.csv", db)))
  readr::write_csv(fp_res$fp_gene_corr_full, file.path(base_dir, sprintf("fp_gene_corr_full_100kb_%s.csv", db)))

}


if (FALSE) {
  # match_mode <- "strict"  # strict lenient
  regulated_genes <- 1.5 # 1.5 2
  delta_link <- 1 # 2 1
  regulation_priors <- "300kb" # "genehancer"
  lighting_folder <- file.path(base_dir, paste("lighting", "fp_tf_corr_FDR", threshold_fp_tf_corr_p, regulation_priors, db, "regulated_genes", regulated_genes, "delta_link", delta_link, sep = "_"))
  print(lighting_folder)


  # Step 3. Build basal GRN & identify active regulatory edges per c --------
  # basal
  basal <- make_basal_links(
    fp_gene_corr_kept = fp_res$fp_gene_corr_kept,
    fp_annotation     = grn_set$fp_annotation,
    out_dir           = lighting_folder,
    prefix            = "lighting"
  )


  light_by_condition(
    ds = grn_set,
    basal_links = basal,
    out_dir = lighting_folder,
    prefix = "lighting",
    label_col = "Sample",
    link_score_threshold = link_score_threshold,
    fp_score_threshold = fp_score_threshold,
    tf_expr_threshold = threshold_tf_expr,
    use_parallel = TRUE,
    workers = 4
  )

  # Step 4. Perform differential GRN analysis & identify master TFs ---------

  specs <- build_cellwise_contrasts_from_index(
    index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
    out_dir = lighting_folder,
    prefix = "lighting",
    ctrl_tag = "Ctrl",
    clean_names = FALSE
  )
  str(specs)


  # Adding comparisons vs. culture:
  # --- Add treatment-vs-culture matched-time contrasts -------------------------

  # map the time token in names like "BG_1h_post" -> target culture name
  time_to_culture <- c(
    "1h_post"        = "culture_1h_post",
    "4h_post"        = "culture_4h_post",
    "24h_post"       = "culture_24h_post",
    "24h_culture_5d" = "culture_6d"
  )

  # lookup: condition name -> its per-condition source path (from existing specs)
  src_lookup <- setNames(specs$cond1_source, specs$cond1_name)

  # base output folder (reuse what specs already uses)
  out_base <- dirname(specs$out_file[1])

  # pick treatments that need culture controls
  treat_rows <- specs[grepl("^(BG|LPS)_", specs$cond1_name) & specs$cond2_name == "Ctrl", , drop = FALSE]

  # derive the matched culture name for each treatment
  time_token <- sub("^[^_]+_", "", treat_rows$cond1_name) # e.g., "1h_post", "24h_culture_5d"
  matched_culture <- unname(time_to_culture[time_token])

  # keep only rows that have a valid mapping and available culture sources
  ok <- !is.na(matched_culture) & matched_culture %in% names(src_lookup) & treat_rows$cond1_name %in% names(src_lookup)

  new_specs <- tibble::tibble(
    cond1_label  = treat_rows$cond1_label[ok],
    cond2_label  = matched_culture[ok],
    cond1_source = unname(src_lookup[treat_rows$cond1_name[ok]]),
    cond2_source = unname(src_lookup[matched_culture[ok]]),
    cond1_name   = treat_rows$cond1_name[ok],
    cond2_name   = matched_culture[ok],
    out_file     = file.path(out_base, paste0(treat_rows$cond1_name[ok], "_vs_", matched_culture[ok], "_delta_links.csv"))
  )

  # append, de-duplicate if needed
  specs2 <- rbind(specs, new_specs)
  specs2 <- specs2[!duplicated(specs2[, c("cond1_name", "cond2_name")]), , drop = FALSE]

  specs <- specs2
  rm(specs2, new_specs, treat_rows, time_token, matched_culture, ok, src_lookup, out_base, time_to_culture)

  # quick check
  specs[, c("cond1_name", "cond2_name", "out_file")]


  # run_links_deltas_driver(
  #   out_dir = lighting_folder,
  #   prefix  = "lighting",
  #   index_csv = file.path(lighting_folder, "lighting_per_condition_index.csv"),
  #   ctrl_tag = "Ctrl",
  #   clean_names = FALSE,
  #   parallel = TRUE
  # )

  run_links_deltas_driver(
    specs       = specs,
    clean_names = FALSE,
    parallel    = TRUE
  )

  delta_csvs <- list.files(lighting_folder, "_delta_links.csv", full.names = TRUE)
  de_gene_log2_abs_min <- if (regulated_genes == 1.5) 0.585 else if (regulated_genes == 2) 1 else NA_real_

  bulk <- episcope::filter_links_deltas_bulk(
    delta_csvs,
    gene_expr_min = threshold_gene_expr,
    tf_expr_min = threshold_tf_expr,
    fp_min = threshold_fp_score,
    link_min = threshold_link_score,
    abs_delta_min = delta_link,
    apply_de_gene = TRUE,
    de_gene_log2_abs_min = de_gene_log2_abs_min,
    enforce_link_expr_sign = TRUE,
    expr_dir_col = "log2FC_gene_expr",
    workers = 20
  )

  bulk$filtered_paths
  bulk$filtered_dirs
  bulk$manifest_path

  # LDA across all filtered CSVs produced above
  filtered_csvs <- bulk$filtered_paths

  lda_out <- episcope::annotate_links_with_lda_topic_bulk(
    filtered_csvs,
    K = 20, # topics
    which = "abs", # |delta_link_score|
    min_df = 2, # drop ultra-rare terms
    gamma_cutoff = 0.2, # keep all topics with gamma >= 0.2 per TF
    parallel = TRUE, workers = 20, plan = "multisession", seed = 1
  )

  lda_out$assigned_paths # vector of "*_filtered_lda_K20.csv"
  lda_out$assigned_dirs # unique directories
  lda_out$manifest_path # CSV mapping filtered -> annotated

  # Summarize all annotated CSVs produced by LDA
  annotated_csvs <- lda_out$assigned_paths

  sum_out <- episcope::summarize_lda_annotations_bulk(
    annotated_csvs,
    edge_filter_min = delta_link, # keep edges used to call delta links
    parallel = TRUE, plan = "multisession", workers = 20,
    use_tf_r_weight = TRUE # optional: weight deltas by TF r column if present
  )

  sum_out$summary_paths # "*_filtered_lda_K20_summary.csv"
  sum_out$summary_dirs
  sum_out$manifest_path


  # Using outputs from summarize_lda_annotations_bulk():
  summary_csvs <- sum_out$summary_paths

  # Generate all bubble & hub PDFs (overall / activate / repress + no-cancel variants)
  episcope::plot_from_summary_bulk(
    summary_csvs,
    top_tf_per_topic = 20,
    parallel = TRUE,
    workers = 20,
    plan = "multisession",
    color_sigma = 2
  )


  # TODO TF HITS score summary plot (waterfall)

  # Step 5. Generate interactive Topic & TF regulatory hub subnetwork
  # ---- Step 5A. Single Δ-topic from a known pair ----

  comp_csv <- file.path(lighting_folder, "AsPC1_10_Gln.Arg_vs_AsPC1_10_FBS_delta_links_filtered_lda_K20.csv")
  summary_csv <- file.path(lighting_folder, "AsPC1_10_Gln.Arg_vs_AsPC1_10_FBS_delta_links_filtered_lda_K20_summary.csv")

  # Render Topic 7 with visual rules (Δ = stress − control)
  episcope::render_link_network_delta_topic_simple(
    comp_csv = comp_csv,
    summary_csv = summary_csv,
    topic_id = 7,
    edge_filter_min = 1,
    gene_fc_thresh = 1.5,
    de_reference = "str_over_ctrl",
    top_n_tfs_per_topic = 20,
    out_html = file.path(
      dirname(comp_csv),
      sprintf(
        "%s_topic-%d_subnetwork_delta.html",
        tools::file_path_sans_ext(basename(comp_csv)), 7
      )
    ),
    verbose = TRUE
  )


  # ---- Step 5B. Bulk Δ-topic across many summaries (auto-pairs) ----
  summary_csvs <- list.files(
    lighting_folder,
    pattern = "_filtered_lda_K20_summary\\.csv$",
    full.names = TRUE
  )

  for (sum_path in summary_csvs) {
    # Try to find the matching annotated comparison file
    cand1 <- sub("_summary\\.csv$", ".csv", sum_path)
    cand2 <- sub("_filtered_lda_K\\d+_summary\\.csv$", "_delta_links.csv", sum_path)
    comp_csv <- if (file.exists(cand1)) cand1 else if (file.exists(cand2)) cand2 else NA_character_
    if (!isTRUE(file.exists(comp_csv))) {
      cli::cli_inform("Skip (no comp CSV found for): {sum_path}")
      next
    }

    S <- readr::read_csv(sum_path, show_col_types = FALSE)
    topic_col <- if ("topic" %in% names(S)) "topic" else if ("main_topic" %in% names(S)) "main_topic" else NA_character_
    if (!isTRUE(topic_col %in% names(S))) next
    if ("topic_rank" %in% names(S)) S <- S[order(S$topic_rank), , drop = FALSE]
    topics <- head(unique(as.integer(S[[topic_col]])), 20)

    for (t in topics) {
      episcope::render_link_network_delta_topic_simple(
        comp_csv = comp_csv,
        summary_csv = sum_path,
        topic_id = t,
        edge_filter_min = 1,
        gene_fc_thresh = 1.5,
        de_reference = "str_over_ctrl",
        top_n_tfs_per_topic = 20,
        out_html = NULL,
        verbose = TRUE
      )
    }
  }

  # 5C TF centric plots
  render_tf_hub_delta_network(
    comp_csv = file.path(
      lighting_folder,
      "BG_24h_culture_5d_vs_culture_6d_delta_links_filtered_lda_K20.csv"
    ),
    input_tf = "TBX21",
    edge_filter_min = 1,
    edge_filter_on = "either",
    gene_fc_thresh = 1.5,
    de_reference = "str_over_ctrl",
    # ring_tf_direct_only = TRUE,   # only direct TFs form the ring
    motif_db = "jaspar2024",
    verbose = TRUE
  )

  comp_csvs <- list.files(lighting_folder, pattern = "_delta_links_filtered_lda_K20\\.csv$", full.names = TRUE)
  tf_list <- c(
    "TBX21", "TCF3", "ZNF610", "SOX4", "ZNF85", "SP4", "REL", "MTF1", "STAT4", "JUNB",
    "MITF", "IRF1", "IRF8", "NRF1", "NFKB1", "NFKB2", "USF2", "STAT2", "STAT5A", "EGR2", "KLF6"
  )

  for (f in comp_csvs) {
    for (tf in tf_list) {
      res <- try(
        {
          render_tf_hub_delta_network(
            comp_csv = f,
            input_tf = tf,
            edge_filter_min = 1,
            edge_filter_on = "either",
            gene_fc_thresh = 1.5,
            de_reference = "str_over_ctrl",
            motif_db = "jaspar2024",
            verbose = TRUE
          )
          TRUE
        },
        silent = TRUE
      )
    }
  }

}






