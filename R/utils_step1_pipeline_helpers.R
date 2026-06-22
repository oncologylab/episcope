# File: utils_step1_pipeline_helpers.R
# Author: Yaoxiang Li
# Created: 2026-03-31
# Updated: 2026-03-31
#
# Purpose:
# Define the user-facing Module 1 wrappers for preparing multi-omic data and
# predicting TF binding sites for the rebuilt package path.
#
# Inputs:
# - YAML config path or explicit ATAC, RNA, metadata, and aligned-footprint
#   inputs
# - footprint preprocessing settings, cache locations, and sample-scope rules
# - Step 1 thresholds and optional motif metadata
#
# Outputs:
# - prepared Step 1 multi-omic object with aligned score, bound, annotation,
#   condition-level, and expression-flag components
# - compact prepared multiomic RDS cache under predict_tf_binding_sites/
# - on-disk RDS cache for fast reruns when requested
#
# Notes:
# - Keep the public API small and stable.
# - Move implementation into focused internal Step 1 files as migration
#   proceeds.

#' Module 1 pipeline helpers
#'
#' @noRd
NULL

.write_fp_score_qn_csv <- function(grn_set, out_dir, db, verbose = TRUE) {
  if (!is.list(grn_set) || !is.data.frame(grn_set$fp_score_condition_qn)) {
    .log_abort("`grn_set$fp_score_condition_qn` is required to write normalized footprint scores.")
  }
  db_use <- if (is.character(db) && length(db) >= 1L && nzchar(db[[1L]])) db[[1L]] else "db"
  out_path <- file.path(out_dir, sprintf("01_fp_scores_qn_%s.csv", db_use))
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(grn_set$fp_score_condition_qn), out_path)
  if (isTRUE(verbose)) {
    .log_inform("Saved quantile-normalized footprint scores: {.path {out_path}}")
  }
  invisible(out_path)
}

#' Load and prepare the Module 1 multi-omic object
#'
#' Build the rebuilt Module 1 data object from cached aligned footprints or from
#' raw footprint overview files plus ATAC, RNA, and sample metadata inputs. The
#' returned object is the canonical input for downstream Step 1 TFBS
#' correlation.
#'
#' @param config Optional YAML config path.
#' @param genome Optional genome string used to override the config value.
#' @param gene_symbol_col Gene-symbol column in the RNA table.
#' @param fp_aligned Optional pre-aligned footprint object.
#' @param do_preprocess Logical; if `TRUE`, load and align raw footprints before
#'   building the object. If `FALSE`, use cached aligned footprints.
#' @param do_motif_clustering Logical; if `TRUE`, run motif clustering during
#'   preprocessing when available.
#' @param trim_hocomoco Logical; trim HOCOMOCO manifests when the trimming helper
#'   is available.
#' @param fp_root_dir Optional root directory for raw footprint overview files.
#' @param fp_cache_dir Cache directory for aligned footprint files.
#' @param fp_cache_tag Cache tag, typically the motif database name.
#' @param footprint_sample_scope Footprint sample selection rule.
#' @param mid_slop,round_digits,score_match_pct Alignment parameters passed to
#'   `align_footprints()`.
#' @param output_mode Output mode for aligned footprints. One of `"full"` or
#'   `"distinct"`.
#' @param write_outputs Logical; if `TRUE`, save the prepared object as an RDS
#'   cache under `predict_tf_binding_sites/`.
#' @param write_fp_score_qn_csv Logical; if `TRUE` and `write_outputs = TRUE`,
#'   also save quantile-normalized footprint scores as
#'   `01_fp_scores_qn_<db>.csv` under the Module 1 output directory.
#' @param atac_data,rna_tbl,metadata Optional in-memory input tables.
#' @param atac_data_path,rna_path,metadata_path Optional explicit file paths for
#'   the input tables.
#' @param step1_out_dir_name Output folder name under `base_dir`.
#' @param label_col Metadata column used to aggregate matched conditions.
#' @param expected_n Optional expected matched sample count.
#' @param tf_list Optional TF allowlist for downstream correlation.
#' @param motif_db Optional motif metadata table.
#' @param threshold_gene_expr Expression threshold for Step 1 expression flags.
#' @param threshold_fp_score Footprint-score threshold for Step 1 bound flags.
#' @param use_parallel Logical; if `TRUE`, allow parallel work in supported
#'   helpers.
#' @param verbose Logical; if `TRUE`, emit concise progress messages.
#' @param time_log Logical; if TRUE, emit elapsed-time messages.
#'
#' @return A rebuilt Module 1 multi-omic object.
#' @examples
#' \donttest{
#' config_path <- "dev/config/pdac_nutrient_stress_strict_jaspar2024_demo.yaml"
#' if (file.exists(config_path)) {
#'   omics_data <- load_prep_multiomic_data(
#'     config = config_path,
#'     genome = "hg38",
#'     label_col = "strict_match_rna",
#'     do_preprocess = FALSE,
#'     verbose = TRUE
#'   )
#' }
#' }
#' @export
load_prep_multiomic_data <- function(
    config = NULL,
    genome = NULL,
    gene_symbol_col = "HGNC",
    fp_aligned = NULL,
    do_preprocess = FALSE,
    do_motif_clustering = FALSE,
    trim_hocomoco = FALSE,
    fp_root_dir = NULL,
    fp_cache_dir = NULL,
    fp_cache_tag = NULL,
    footprint_sample_scope = "metadata",
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 0.8,
    output_mode = c("full", "distinct"),
    write_outputs = FALSE,
    write_fp_score_qn_csv = TRUE,
    atac_data = NULL,
    rna_tbl = NULL,
    metadata = NULL,
    atac_data_path = NULL,
    rna_path = NULL,
    metadata_path = NULL,
    step1_out_dir_name = "predict_tf_binding_sites",
    label_col,
    expected_n = NULL,
    tf_list = NULL,
    motif_db = NULL,
    threshold_gene_expr = NULL,
    threshold_fp_score = NULL,
    use_parallel = TRUE,
    verbose = TRUE,
    time_log = verbose
) {
  start_time <- Sys.time()
  if (isTRUE(time_log)) {
    on.exit(.log_inform("Module 1 multiomic preparation completed in {round(as.numeric(difftime(Sys.time(), start_time, units = 'secs')), 1)} sec."), add = TRUE)
  }
  output_mode <- match.arg(output_mode)

  if (!is.null(config)) {
    if (is.character(config) && length(config) == 1L && file.exists(config)) {
      load_config(config)
    } else {
      .log_abort("`config` must be a path to a YAML file.")
    }
  }
  if (!is.null(genome) && nzchar(genome)) {
    .cfg_set("ref_genome", genome)
  }

  get_cfg <- function(name, default = NULL) .cfg_get(name, default = default)
  resolve_path <- function(arg, candidates) {
    if (!is.null(arg) && nzchar(arg)) return(arg)
    for (nm in candidates) {
      val <- get_cfg(nm)
      if (is.character(val) && length(val) && nzchar(val[1])) return(val[1])
    }
    NULL
  }

  step1_out_dir_name_use <- if (!is.null(step1_out_dir_name) && nzchar(step1_out_dir_name)) {
    as.character(step1_out_dir_name)
  } else {
    "predict_tf_binding_sites"
  }

  if (is.null(fp_cache_dir)) {
    base_dir_cfg <- get_cfg("base_dir")
    if (is.character(base_dir_cfg) && nzchar(base_dir_cfg)) {
      fp_cache_dir <- file.path(base_dir_cfg, "cache")
    }
  }
  if (is.null(fp_cache_tag)) {
    fp_cache_tag <- get_cfg("db")
  }
  if (is.null(threshold_gene_expr)) {
    threshold_gene_expr <- get_cfg("threshold_gene_expr")
  }
  if (is.null(threshold_fp_score)) {
    threshold_fp_score <- get_cfg("threshold_fp_score")
  }
  if (is.null(tf_list)) {
    tf_list <- get_cfg("tf_list")
  }
  if (is.null(motif_db)) {
    motif_db <- get_cfg("motif_db")
  }

  compact_rds_path <- NULL
  base_dir_cfg <- get_cfg("base_dir")
  if (is.character(base_dir_cfg) && nzchar(base_dir_cfg)) {
    compact_rds_path <- file.path(
      base_dir_cfg,
      step1_out_dir_name_use,
      sprintf("01_multiomic_data_object_%s.rds", get_cfg("db"))
    )
  }

  if (is.null(fp_aligned)) {
    if (isTRUE(do_preprocess)) {
      if (!exists("load_footprints", mode = "function") || !exists("align_footprints", mode = "function")) {
        .log_abort("Preprocessing requires `load_footprints()` and `align_footprints()`, which are not yet available in the rebuilt Step 1 path.")
      }
      if (is.null(fp_root_dir)) fp_root_dir <- get_cfg("fp_root_dir")
      if (is.null(fp_root_dir) || !nzchar(fp_root_dir)) {
        .log_abort("`fp_root_dir` is not set. Provide it or set it in config.")
      }
      db_use <- get_cfg("db")
      if (!is.character(db_use) || !nzchar(db_use)) {
        .log_abort("`db` must be set in config when preprocessing footprints.")
      }
      fp_scope_cfg <- get_cfg("footprint_sample_scope")
      fp_scope <- if (is.character(fp_scope_cfg) && length(fp_scope_cfg) >= 1L && nzchar(fp_scope_cfg[[1]])) {
        as.character(fp_scope_cfg[[1]])
      } else {
        as.character(footprint_sample_scope)
      }
      fp_scope <- tolower(fp_scope)
      fp_sample_ids <- NULL
      if (identical(fp_scope, "metadata")) {
        meta_for_fp <- metadata
        if (is.null(meta_for_fp)) {
          meta_path_for_fp <- resolve_path(metadata_path, c("sample_metadata", "metadata_path", "sample_metadata_path"))
          if (!is.null(meta_path_for_fp) && file.exists(meta_path_for_fp)) {
            meta_for_fp <- readr::read_csv(meta_path_for_fp, show_col_types = FALSE)
          }
        }
        if (is.data.frame(meta_for_fp) && "id" %in% names(meta_for_fp)) {
          fp_id_col <- if ("id_fp" %in% names(meta_for_fp)) "id_fp" else "id"
          fp_sample_ids <- unique(as.character(meta_for_fp[[fp_id_col]]))
          fp_sample_ids <- fp_sample_ids[!is.na(fp_sample_ids) & nzchar(fp_sample_ids)]
        }
      }
      fp_manifest <- load_footprints(
        root_dir = fp_root_dir,
        db_name = db_use,
        out_dir = file.path(fp_cache_dir, paste0("fp_", db_use)),
        sample_ids = fp_sample_ids
      )
      fp_manifest_trim_fn <- get0("fp_manifest_trim", mode = "function")
      if (isTRUE(trim_hocomoco) && identical(db_use, "HOCOMOCOv13") && is.function(fp_manifest_trim_fn)) {
        fp_manifest <- fp_manifest_trim_fn(fp_manifest)
      }
      fp_aligned <- align_footprints(
        fp_manifest,
        mid_slop = mid_slop,
        round_digits = round_digits,
        score_match_pct = score_match_pct,
        cache_dir = fp_cache_dir,
        cache_tag = db_use,
        output_mode = output_mode
      )
      run_fp_motif_clustering_fn <- get0("run_fp_motif_clustering", mode = "function")
      if (isTRUE(do_motif_clustering) && is.function(run_fp_motif_clustering_fn)) {
        cluster_out_dir <- file.path(
          get_cfg("base_dir"),
          step1_out_dir_name_use,
          "04_tf_motif_clustering"
        )
        run_fp_motif_clustering_fn(
          fp_aligned = fp_aligned,
          out_dir = cluster_out_dir,
          base_dir = get_cfg("base_dir"),
          ref_db = db_use,
          motif_db = motif_db,
          mode = "data",
          target_clusters = 220,
          qc_mode = "fast",
          save_motif_db = TRUE
        )
      }
      if (exists("plot_fp_merge_summary", mode = "function")) {
        plot_fp_merge_summary(
          fp_aligned,
          out_dir = file.path(get_cfg("base_dir"), step1_out_dir_name_use),
          db = db_use,
          verbose = verbose
        )
      }
    } else {
      fp_aligned <- load_fp_aligned_from_cache(
        cache_dir = fp_cache_dir,
        cache_tag = fp_cache_tag,
        output_mode = output_mode,
        verbose = verbose
      )
    }
  }

  if (is.null(metadata)) {
    path <- resolve_path(metadata_path, c("sample_metadata", "metadata_path", "sample_metadata_path"))
    if (is.null(path) || !file.exists(path)) {
      .log_abort("Missing `metadata` and no valid configured path.")
    }
    metadata <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (is.null(atac_data)) {
    path <- resolve_path(atac_data_path, c("atac_master", "atac_data_path", "atac_master_path"))
    if (is.null(path) || !file.exists(path)) {
      .log_abort("Missing `atac_data` and no valid configured path.")
    }
    atac_data <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (is.null(rna_tbl)) {
    path <- resolve_path(rna_path, c("rna_mapped", "rna_mapped_path", "rna_path"))
    if (is.null(path) || !file.exists(path)) {
      .log_abort("Missing `rna_tbl` and no valid configured path.")
    }
    rna_tbl <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (!is.null(gene_symbol_col) && gene_symbol_col != "HGNC") {
    if (gene_symbol_col %in% names(rna_tbl) && !"HGNC" %in% names(rna_tbl)) {
      rna_tbl <- dplyr::rename(rna_tbl, HGNC = dplyr::all_of(gene_symbol_col))
    }
  }

  atac_out <- load_atac(atac_data, sort_peaks = TRUE)

  grn_set <- build_grn_set(
    fp_score = fp_aligned$fp_score,
    fp_bound = fp_aligned$fp_bound,
    fp_annotation = fp_aligned$fp_annotation,
    atac_score = atac_out$score,
    atac_overlap = atac_out$overlap,
    rna = rna_tbl,
    metadata = metadata,
    tf_list = tf_list,
    motif_db = motif_db,
    label_col = label_col,
    expected_n = expected_n
  )

  grn_set <- grn_add_fp_score_condition(grn_set, label_col = label_col)
  grn_set <- grn_add_fp_bound_condition(
    grn_set,
    label_col = label_col,
    threshold_fp_score = threshold_fp_score
  )
  grn_set <- grn_filter_fp_bound_condition(
    grn_set,
    min_bound = 1L,
    use_parallel = use_parallel
  )
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  grn_set <- grn_add_rna_expressed(
    grn_set,
    label_col = label_col,
    threshold_gene_expr = threshold_gene_expr
  )
  grn_set <- grn_add_rna_condition(grn_set, label_col = label_col)

  compact <- as_multiomic_object(
    grn_set,
    project = list(db = get_cfg("db"), ref_genome = get_cfg("ref_genome"), label_col = label_col),
    label_col = label_col,
    verbose = FALSE
  )

  if (isTRUE(write_outputs)) {
    out_dir <- file.path(base_dir_cfg, step1_out_dir_name_use)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (exists("plot_fp_norm_bound_qc", mode = "function")) {
      plot_fp_norm_bound_qc(
        omics_data = grn_set,
        out_dir = out_dir,
        db = get_cfg("db"),
        threshold_fp_score = threshold_fp_score,
        verbose = verbose
      )
    }
    if (exists("plot_gene_expr_qc", mode = "function")) {
      plot_gene_expr_qc(
        omics_data = grn_set,
        out_dir = out_dir,
        db = get_cfg("db"),
        threshold_gene_expr = threshold_gene_expr,
        verbose = verbose
      )
    }
    if (isTRUE(write_fp_score_qn_csv)) {
      .write_fp_score_qn_csv(
        grn_set = grn_set,
        out_dir = out_dir,
        db = get_cfg("db"),
        verbose = verbose
      )
    }
    if (is.character(compact_rds_path) && nzchar(compact_rds_path)) {
      dir.create(dirname(compact_rds_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(compact, compact_rds_path)
      if (isTRUE(verbose)) .log_inform("Saved Module 1 multiomic object: {.path {compact_rds_path}}")
    }
  }

  compact
}
