# File: utils_step1_direct_bound_filter.R
# Author: Yaoxiang Li
# Created: 2026-04-01
# Updated: 2026-04-01
#
# Purpose:
# Provide a post-Step-1 direct-binding filter that derives condition-specific
# direct-bound peak BED files from Module 1 TFBS overview outputs.
#
# Inputs:
# - Step 1 output directory containing overview files and numbered prep outputs
# - condition names to filter
# - motif database label and reference genome
# - correlation cutoff and optional TF subset
#
# Outputs:
# - one direct-bound BED file per requested condition
# - one summary TSV per requested condition
#
# Notes:
# - This is a post-Step-1 filter layer. It should not change the core Step 1
#   correlation engine.
# - Peaks are kept when they have at least one direct TF-binding assignment for
#   the requested condition and exceed the supplied correlation cutoff.

#' Module 1 direct-bound post-filter helpers
#'
#' @noRd
NULL

.step1_overview_dir <- function(step1_out_dir, db, mode = "all") {
  file.path(step1_out_dir, sprintf("06_fp_predicted_tfbs_%s%s", db, .step1_mode_suffix(mode)))
}

.step1_direct_filter_summary_path <- function(out_dir, condition, r_cutoff) {
  file.path(
    out_dir,
    sprintf(
      "%s_direct_filter_r_%s_summary.tsv",
      condition,
      formatC(r_cutoff, format = "f", digits = 2)
    )
  )
}

.step1_direct_filter_bed_path <- function(out_dir, condition, r_cutoff) {
  file.path(
    out_dir,
    sprintf(
      "%s_direct_filter_r_%s.bed",
      condition,
      formatC(r_cutoff, format = "f", digits = 2)
    )
  )
}

#' Write condition-specific direct-bound peak BED files from Module 1 outputs
#'
#' Apply a post-Step-1 direct-binding filter using Module 1 overview files,
#' footprint-bound condition flags, gene-expression flags, and the packaged
#' motif cluster map for the requested genome. For each requested condition,
#' this writes a BED file containing peaks that pass the correlation cutoff and
#' have at least one direct TF-binding assignment in that condition.
#'
#' @param step1_out_dir Module 1 output directory, typically
#'   `predict_tf_binding_sites/`.
#' @param conditions Character vector of condition names, for example
#'   `c("AsPC1_10_FBS", "HPAFII_10_FBS")`.
#' @param db Motif database label used in the Step 1 outputs, default
#'   `"JASPAR2024"`.
#' @param genome Reference genome used to select the packaged cluster map, for
#'   example `"hg38"` or `"mm10"`.
#' @param r_cutoff Minimum absolute correlation required to keep a TFBS.
#' @param out_dir Output directory. Defaults to
#'   `file.path(step1_out_dir, "direct_bound_filter")`.
#' @param tf_subset Optional TF subset.
#' @param n_cores Number of cores used for per-TF filtering.
#' @param verbose Logical; if `TRUE`, emit concise progress messages.
#'
#' @return A named character vector of written BED-file paths.
#' @examples
#' \dontrun{
#' write_direct_bound_filter_bed(
#'   step1_out_dir = file.path("predict_tf_binding_sites"),
#'   conditions = c("AsPC1_10_FBS"),
#'   db = "JASPAR2024",
#'   genome = "hg38",
#'   r_cutoff = 0.30
#' )
#' }
#' @keywords internal
write_direct_bound_filter_bed <- function(
    step1_out_dir,
    conditions,
    db = "JASPAR2024",
    genome,
    r_cutoff = 0.30,
    out_dir = NULL,
    tf_subset = NULL,
    n_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
    verbose = TRUE
) {
  if (!is.character(step1_out_dir) || !nzchar(step1_out_dir)) {
    .log_abort("`step1_out_dir` must be a non-empty path.")
  }
  if (!is.character(conditions) || !length(conditions)) {
    .log_abort("`conditions` must be a non-empty character vector.")
  }
  conditions <- unique(conditions[nzchar(conditions)])
  if (!length(conditions)) {
    .log_abort("`conditions` must include at least one non-empty value.")
  }
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }
  if (!is.character(genome) || !nzchar(genome)) {
    .log_abort("`genome` must be a non-empty string.")
  }
  if (!is.finite(r_cutoff)) {
    .log_abort("`r_cutoff` must be a finite numeric scalar.")
  }
  if (is.null(out_dir)) {
    out_dir <- file.path(step1_out_dir, "direct_bound_filter")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  overview_dir <- .step1_overview_dir(step1_out_dir, db = db, mode = "all")
  if (!dir.exists(overview_dir)) {
    .log_abort("Step 1 overview directory not found: {overview_dir}")
  }

  fp_annotation_file <- file.path(step1_out_dir, sprintf("03_fp_annotation_%s.csv", db))
  fp_bound_condition_file <- file.path(step1_out_dir, sprintf("03_fp_bound_condition_%s.csv", db))
  gene_expr_flag_file <- file.path(step1_out_dir, sprintf("05_gene_expr_flag_%s.csv", db))
  cluster_map_file <- resolve_motif_db_path(db = db, ref_genome = genome)

  for (path in c(fp_annotation_file, fp_bound_condition_file, gene_expr_flag_file, cluster_map_file)) {
    if (!file.exists(path)) .log_abort("Required input file not found: {path}")
  }

  overview_files <- sort(list.files(overview_dir, pattern = "_overview\\.txt$", full.names = TRUE))
  if (!length(overview_files)) {
    .log_abort("No *_overview.txt files found under {overview_dir}")
  }

  overview_tf_names <- sub("_overview\\.txt$", "", basename(overview_files))
  if (is.character(tf_subset) && length(tf_subset)) {
    tf_subset <- unique(tf_subset[nzchar(tf_subset)])
    overview_keep <- overview_tf_names %in% tf_subset
    overview_files <- overview_files[overview_keep]
    overview_tf_names <- overview_tf_names[overview_keep]
  }
  if (!length(overview_files)) {
    .log_abort("No overview files remain after applying `tf_subset`.")
  }

  fp_annotation_dt <- data.table::fread(
    fp_annotation_file,
    sep = ",",
    select = c("fp_peak", "motifs"),
    showProgress = FALSE
  )
  if (identical(db, "HOCOMOCOv13") && "motifs" %in% names(fp_annotation_dt)) {
    fp_annotation_dt[, motifs := sub("^_+", "", motifs)]
  }
  cluster_map_dt <- data.table::fread(cluster_map_file, sep = "\t", showProgress = FALSE)
  gene_col <- if ("gene_symbol" %in% names(cluster_map_dt)) "gene_symbol" else if ("HGNC" %in% names(cluster_map_dt)) "HGNC" else NULL
  if (is.null(gene_col) || !"motif" %in% names(cluster_map_dt)) {
    .log_abort("Cluster map file must contain motif and gene-symbol columns.")
  }
  cluster_map_dt <- unique(cluster_map_dt[, .(motif, TFBS_name = get(gene_col))])
  if (identical(db, "HOCOMOCOv13") && "motif" %in% names(cluster_map_dt)) {
    cluster_map_dt[, motif := sub("^_+", "", motif)]
  }
  cluster_map_dt <- cluster_map_dt[, .(TFBS_name = unlist(strsplit(TFBS_name, "::", fixed = TRUE))), by = motif]
  cluster_map_dt[, TFBS_name := trimws(TFBS_name)]
  cluster_map_dt <- unique(cluster_map_dt[nzchar(TFBS_name)])

  fp_annotation_dt <- merge(
    fp_annotation_dt,
    cluster_map_dt,
    by.x = "motifs",
    by.y = "motif",
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )[, .(fp_peak, TFBS_name)]
  fp_annotation_dt <- unique(fp_annotation_dt)

  fp_peak_universe_dt <- unique(fp_annotation_dt[, .(fp_peak)])
  fp_peak_universe_dt[, fp_peak_idx := .I]
  fp_peak_idx_dt <- data.table::copy(fp_peak_universe_dt)
  data.table::setkey(fp_annotation_dt, fp_peak)
  data.table::setkey(fp_peak_idx_dt, fp_peak)

  fp_annotation_idx_dt <- fp_peak_idx_dt[
    fp_annotation_dt,
    on = "fp_peak",
    nomatch = 0L
  ][, .(TFBS_name, fp_peak_idx)]
  fp_annotation_idx_dt <- unique(fp_annotation_idx_dt)
  fp_peak_idx_by_tf <- split(fp_annotation_idx_dt$fp_peak_idx, fp_annotation_idx_dt$TFBS_name)

  fp_bound_condition_dt <- data.table::fread(fp_bound_condition_file, sep = ",", showProgress = FALSE)
  gene_expr_flag_dt <- data.table::fread(gene_expr_flag_file, sep = ",", showProgress = FALSE)
  gene_symbol_col <- if ("HGNC" %in% names(gene_expr_flag_dt)) "HGNC" else "gene_symbol"
  if (!gene_symbol_col %in% names(gene_expr_flag_dt)) {
    .log_abort("Gene-expression flag file must include HGNC or gene_symbol.")
  }

  overview_direct_cols <- c(
    "TFBS_chr", "TFBS_start", "TFBS_end",
    "pearson_corr_fp_tf_r",
    "spearman_corr_fp_tf_r"
  )

  prepare_one_direct_filter_dt <- function(overview_path, bound_peak_ids_condition, tf_name, r_cutoff) {
    direct_bound_col <- paste0(tf_name, "_direct_bound")
    direct_bound_vec <- integer(nrow(fp_peak_universe_dt))
    overview_dt <- data.table::fread(
      overview_path,
      sep = "\t",
      select = overview_direct_cols,
      showProgress = FALSE
    )

    corr_r_abs_max <- pmax(
      abs(overview_dt[["pearson_corr_fp_tf_r"]]),
      abs(overview_dt[["spearman_corr_fp_tf_r"]]),
      na.rm = TRUE
    )
    keep_idx <- which(is.finite(corr_r_abs_max) & corr_r_abs_max > r_cutoff)

    if (!length(keep_idx)) {
      return(list(
        col_name = direct_bound_col,
        values = direct_bound_vec,
        summary = data.table::data.table(
          tf = tf_name,
          n_overview_rows = nrow(overview_dt),
          n_keep_r = 0L,
          n_keep_r_and_bound = 0L,
          n_direct_positive = 0L
        )
      ))
    }

    keep_fp_peak <- unique(sprintf(
      "%s:%s-%s",
      overview_dt[["TFBS_chr"]][keep_idx],
      overview_dt[["TFBS_start"]][keep_idx],
      overview_dt[["TFBS_end"]][keep_idx]
    ))
    keep_fp_peak <- intersect(keep_fp_peak, bound_peak_ids_condition)

    if (!length(keep_fp_peak)) {
      return(list(
        col_name = direct_bound_col,
        values = direct_bound_vec,
        summary = data.table::data.table(
          tf = tf_name,
          n_overview_rows = nrow(overview_dt),
          n_keep_r = length(keep_idx),
          n_keep_r_and_bound = 0L,
          n_direct_positive = 0L
        )
      ))
    }

    keep_fp_peak_idx <- fp_peak_idx_dt[
      .(fp_peak = keep_fp_peak),
      on = "fp_peak",
      nomatch = 0L
    ][["fp_peak_idx"]]

    tf_direct_idx <- fp_peak_idx_by_tf[[tf_name]]
    if (length(tf_direct_idx) && length(keep_fp_peak_idx)) {
      direct_bound_vec[intersect(keep_fp_peak_idx, tf_direct_idx)] <- 1L
    }

    list(
      col_name = direct_bound_col,
      values = direct_bound_vec,
      summary = data.table::data.table(
        tf = tf_name,
        n_overview_rows = nrow(overview_dt),
        n_keep_r = length(keep_idx),
        n_keep_r_and_bound = length(keep_fp_peak),
        n_direct_positive = sum(direct_bound_vec)
      )
    )
  }

  out_paths <- setNames(character(length(conditions)), conditions)
  for (condition in conditions) {
    condition_bound_col <- if (condition %in% names(fp_bound_condition_dt)) {
      condition
    } else {
      paste0(condition, "_bound")
    }
    if (!condition_bound_col %in% names(fp_bound_condition_dt)) {
      .log_abort("Condition-bound column not found: {condition_bound_col}")
    }
    if (!condition %in% names(gene_expr_flag_dt)) {
      .log_abort("Gene-expression condition column not found: {condition}")
    }

    bound_peak_ids_condition <- unique(fp_bound_condition_dt[get(condition_bound_col) == 1, peak_ID])
    bound_peak_ids_condition <- bound_peak_ids_condition[!is.na(bound_peak_ids_condition) & nzchar(bound_peak_ids_condition)]
    expressed_tf_vec <- unique(gene_expr_flag_dt[get(condition) == 1, get(gene_symbol_col)])
    expressed_tf_vec <- trimws(as.character(expressed_tf_vec))
    expressed_tf_vec <- expressed_tf_vec[!is.na(expressed_tf_vec) & nzchar(expressed_tf_vec)]

    keep_overview <- overview_tf_names %in% expressed_tf_vec
    overview_files_condition <- overview_files[keep_overview]
    overview_tf_names_condition <- overview_tf_names[keep_overview]

    if (!length(overview_files_condition)) {
      .log_warn("No overview files remain after filtering by expressed TFs for condition {.val {condition}}.")
      bed_path <- .step1_direct_filter_bed_path(out_dir, condition, r_cutoff)
      data.table::fwrite(data.table::data.table(V1 = character(), V2 = integer(), V3 = integer()), bed_path, sep = "\t", col.names = FALSE)
      out_paths[[condition]] <- bed_path
      next
    }

    if (isTRUE(verbose)) {
      .log_inform("Direct-bound filter for {.val {condition}}: {length(overview_files_condition)} expressed TF(s), {length(bound_peak_ids_condition)} condition-bound peak(s).")
    }

    direct_filter_prepare_list <- if (.Platform$OS.type == "unix" && n_cores > 1L) {
      parallel::mclapply(
        seq_along(overview_files_condition),
        function(i) {
          prepare_one_direct_filter_dt(
            overview_path = overview_files_condition[[i]],
            bound_peak_ids_condition = bound_peak_ids_condition,
            tf_name = overview_tf_names_condition[[i]],
            r_cutoff = r_cutoff
          )
        },
        mc.cores = n_cores
      )
    } else {
      lapply(
        seq_along(overview_files_condition),
        function(i) {
          prepare_one_direct_filter_dt(
            overview_path = overview_files_condition[[i]],
            bound_peak_ids_condition = bound_peak_ids_condition,
            tf_name = overview_tf_names_condition[[i]],
            r_cutoff = r_cutoff
          )
        }
      )
    }

    direct_filter_prepare_dt <- fp_peak_universe_dt[, .(fp_peak)]
    for (res in direct_filter_prepare_list) {
      direct_filter_prepare_dt[, (res$col_name) := res$values]
    }
    direct_filter_summary_dt <- data.table::rbindlist(lapply(direct_filter_prepare_list, `[[`, "summary"), use.names = TRUE, fill = TRUE)
    direct_bound_cols <- setdiff(names(direct_filter_prepare_dt), "fp_peak")
    direct_filter_prepare_dt[
      ,
      direct_bound := as.integer(rowSums(.SD) > 0L),
      .SDcols = direct_bound_cols
    ]

    direct_filter_bed_dt <- data.table::copy(direct_filter_prepare_dt[direct_bound == 1L, .(fp_peak)])
    if (nrow(direct_filter_bed_dt)) {
      direct_filter_bed_dt[, c("V1", "V2", "V3") := data.table::tstrsplit(fp_peak, "[:-]", fixed = FALSE)]
      direct_filter_bed_dt <- direct_filter_bed_dt[, .(
        V1,
        V2 = as.integer(V2),
        V3 = as.integer(V3)
      )]
    } else {
      direct_filter_bed_dt <- data.table::data.table(
        V1 = character(),
        V2 = integer(),
        V3 = integer()
      )
    }

    bed_path <- .step1_direct_filter_bed_path(out_dir, condition, r_cutoff)
    summary_path <- .step1_direct_filter_summary_path(out_dir, condition, r_cutoff)
    data.table::fwrite(direct_filter_bed_dt, bed_path, sep = "\t", col.names = FALSE)
    data.table::fwrite(direct_filter_summary_dt, summary_path, sep = "\t")
    out_paths[[condition]] <- bed_path

    if (isTRUE(verbose)) {
      .log_inform("Direct-bound BED written for {.val {condition}}: {bed_path}")
    }
  }

  invisible(out_paths)
}

utils::globalVariables(c("V1", "V2", "V3", "direct_bound", "fp_peak_idx", "motif"))
