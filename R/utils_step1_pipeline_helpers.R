#' Module 1 pipeline helpers
#'
#' Lightweight helpers for assembling multi-omic inputs in Module 1.
#'
#' @name module1_pipeline_helpers
#' @rdname module1_pipeline_helpers
#' @noRd
NULL

#' Load aligned footprint outputs from cache
#'
#' Reads cached aligned footprint CSVs produced by [align_footprints()].
#'
#' @param cache_dir Directory containing cached footprint files.
#' @param cache_tag Cache tag used in file names (e.g. "JASPAR2024").
#' @param output_mode Output mode passed through (see [align_footprints()]).
#' @param verbose Emit status messages.
#'
#' @return List with tibbles: \code{fp_score}, \code{fp_bound},
#'   \code{fp_annotation}, and \code{id_map}.
#' @export
load_fp_aligned_from_cache <- function(cache_dir,
                                       cache_tag,
                                       output_mode = c("full", "distinct"),
                                       verbose = TRUE) {
  .assert_pkg("readr")
  output_mode <- match.arg(output_mode)
  cache_paths <- list(
    fp_bound = file.path(cache_dir, sprintf("fp_bounds_%s.csv", cache_tag)),
    fp_score = file.path(cache_dir, sprintf("fp_scores_%s.csv", cache_tag)),
    fp_annotation = file.path(cache_dir, sprintf("fp_annotation_%s.csv", cache_tag)),
    id_map = file.path(cache_dir, sprintf("fp_id_map_%s.csv", cache_tag))
  )
  if (!all(file.exists(unlist(cache_paths[c("fp_bound", "fp_score", "fp_annotation")])))) {
    .log_abort("Cached aligned footprints not found for tag={cache_tag} in {cache_dir}.")
  }
  if (isTRUE(verbose)) .log_inform("Using cached aligned footprints from {.path {cache_dir}} (tag = {cache_tag}).")
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
  out_cached
}

#' Load and assemble a combined multi-omic data object
#'
#' Runs ATAC preprocessing, builds a GRN-aligned object, and adds
#' binding/normalization summaries followed by RNA expression flags.
#'
#' @param config Optional YAML config path; when provided, loads settings via
#'   [load_config()] and uses configured paths for input data and caches.
#' @param genome Optional genome string (e.g. "hg38", "mm10"). When provided,
#'   sets \code{ref_genome} for motif resolution.
#' @param gene_symbol_col Column in RNA table to treat as gene symbols. If not
#'   \code{"HGNC"}, it will be renamed to \code{"HGNC"} before processing.
#' @param fp_aligned Aligned footprint object returned by align_footprints().
#' @param do_preprocess If TRUE and \code{fp_aligned} is NULL, run footprint
#'   preprocessing (load + align) using config paths.
#' @param do_motif_clustering If TRUE, run motif clustering on aligned footprints.
#' @param trim_hocomoco If TRUE, trims HOCOMOCOv13 footprint annotations after
#'   loading the manifest (ignored for other databases).
#' @param fp_root_dir Optional footprint root directory (overrides config).
#' @param fp_cache_dir Optional cache directory for footprint outputs.
#' @param fp_cache_tag Optional cache tag used for aligned footprint files.
#' @param mid_slop Integer; midpoint tolerance in bp for aligning footprints.
#' @param round_digits Integer; digits for rounding scores before matching.
#' @param score_match_pct Numeric; fraction of score columns that must match.
#' @param output_mode Output mode passed to [align_footprints()].
#' @param write_outputs If TRUE, write GRN outputs and QC PDFs.
#' @param atac_data ATAC master table (tibble/data.frame). If NULL, load from
#'   \code{atac_data_path} or YAML config (\code{atac_master} preferred).
#' @param rna_tbl RNA table with columns already matched to sample IDs. If NULL,
#'   load from \code{rna_path} or YAML config (\code{rna_mapped} preferred).
#' @param metadata Sample metadata with id column. If NULL, load from
#'   \code{metadata_path} or YAML config (\code{sample_metadata} preferred).
#' @param atac_data_path Optional CSV path for ATAC master table.
#' @param rna_path Optional CSV path for RNA (mapped) table.
#' @param metadata_path Optional CSV path for sample metadata.
#' @param fp_cache_dir Optional cache directory for aligned footprint files.
#' @param fp_cache_tag Optional cache tag used by [align_footprints()] when
#'   writing cached outputs.
#' @param label_col Metadata column used for condition labels.
#' @param expected_n Optional expected number of matched samples.
#' @param tf_list Character vector of TFs.
#' @param motif_db Motif database tibble.
#' @param threshold_gene_expr Expression threshold for gene flags.
#' @param threshold_fp_score Footprint score threshold for bound flags.
#' @param use_parallel Use parallel filtering when available.
#' @param verbose Emit status messages.
#'
#' @return Multi-omic data list produced by build_grn_set(), updated with condition matrices.
#' @export
load_multiomic_data <- function(
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
    mid_slop = 10L,
    round_digits = 1L,
    score_match_pct = 0.8,
    output_mode = c("full", "distinct"),
    write_outputs = FALSE,
    atac_data = NULL,
    rna_tbl = NULL,
    metadata = NULL,
    atac_data_path = NULL,
    rna_path = NULL,
    metadata_path = NULL,
    label_col,
    expected_n,
    tf_list,
    motif_db,
    threshold_gene_expr,
    threshold_fp_score,
    use_parallel = TRUE,
    verbose = TRUE
) {
  output_mode <- match.arg(output_mode)
  if (!is.null(config)) {
    if (is.character(config) && length(config) == 1L && file.exists(config)) {
      load_config(config)
    } else {
      .log_abort("`config` must be a path to a YAML file.")
    }
  }
  if (!is.null(genome) && nzchar(genome)) {
    assign("ref_genome", genome, envir = .GlobalEnv)
  }

  if (isTRUE(verbose)) {
    .log_inform("Loading input data and creating a combined multi-omic data object.")
  }

  get_cfg <- function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      return(get(name, envir = .GlobalEnv))
    }
    NULL
  }

  resolve_path <- function(arg, candidates) {
    if (!is.null(arg) && nzchar(arg)) return(arg)
    for (nm in candidates) {
      val <- get_cfg(nm)
      if (is.character(val) && length(val) && nzchar(val[1])) return(val[1])
    }
    NULL
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

  if (is.null(fp_aligned)) {
    if (isTRUE(do_preprocess)) {
      if (is.null(fp_root_dir)) fp_root_dir <- get_cfg("fp_root_dir")
      if (is.null(fp_root_dir) || !nzchar(fp_root_dir)) {
        .log_abort("`fp_root_dir` not set. Provide it or set in config.")
      }
      db_use <- get_cfg("db")
      if (!is.character(db_use) || !nzchar(db_use)) {
        .log_abort("`db` must be set in config when preprocessing footprints.")
      }
      if (isTRUE(verbose)) {
        .log_inform("Loading footprint files for {db_use}.")
      }
      fp_manifest <- load_footprints(
        root_dir = fp_root_dir,
        db_name = db_use,
        out_dir = file.path(fp_cache_dir, paste0("fp_", db_use))
      )
      if (isTRUE(trim_hocomoco) && identical(db_use, "HOCOMOCOv13") && !isTRUE(attr(fp_manifest, "from_cache"))) {
        if (isTRUE(verbose)) {
          .log_inform("Trimming HOCOMOCOv13 footprint manifest annotations.")
        }
        fp_manifest <- fp_manifest_trim(fp_manifest)
        fp_manifest_trim_annots(fp_manifest, n_workers = 18, verbose = verbose)
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
      plot_fp_merge_summary(
        fp_aligned,
        out_dir = file.path(get_cfg("base_dir"), "predict_tf_binding_sites"),
        db = db_use,
        verbose = verbose
      )

      if (isTRUE(do_motif_clustering)) {
        run_fp_motif_clustering_pre_if_needed(
          fp_aligned = fp_aligned,
          base_dir = get_cfg("base_dir"),
          ref_db = db_use,
          motif_db = motif_db
        )
        run_fp_motif_clustering(
          fp_aligned = fp_aligned,
          base_dir = get_cfg("base_dir"),
          ref_db = db_use,
          motif_db = motif_db,
          mode = "data",
          target_clusters = 220,
          qc_mode = "fast",
          save_motif_db = TRUE
        )
        run_fp_motif_clustering(
          fp_aligned = fp_aligned,
          base_dir = get_cfg("base_dir"),
          ref_db = db_use,
          motif_db = motif_db,
          mode = "hybrid",
          target_clusters = 165,
          qc_mode = "fast",
          save_motif_db = TRUE
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
    if (is.null(path)) {
      .log_abort("Missing `metadata` and no configured path. Set `sample_metadata` (preferred) in YAML or pass `metadata_path`.")
    }
    if (!file.exists(path)) .log_abort("Metadata file not found: {path}")
    if (isTRUE(verbose)) .log_inform("Loading metadata: {path}")
    metadata <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (is.null(atac_data)) {
    path <- resolve_path(atac_data_path, c("atac_master", "atac_data_path", "atac_master_path"))
    if (is.null(path)) {
      .log_abort("Missing `atac_data` and no configured path. Set `atac_master` (preferred) in YAML or pass `atac_data_path`.")
    }
    if (!file.exists(path)) .log_abort("ATAC master file not found: {path}")
    if (isTRUE(verbose)) .log_inform("Loading ATAC master table: {path}")
    atac_data <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (is.null(rna_tbl)) {
    path <- resolve_path(rna_path, c("rna_mapped", "rna_mapped_path", "rna_path"))
    if (is.null(path)) {
      .log_abort("Missing `rna_tbl` and no configured path. Set `rna_mapped` (preferred) in YAML or pass `rna_path`.")
    }
    if (!file.exists(path)) .log_abort("RNA table not found: {path}")
    if (isTRUE(verbose)) .log_inform("Loading RNA table: {path}")
    rna_tbl <- readr::read_csv(path, show_col_types = FALSE)
  }

  if (!is.null(gene_symbol_col) && gene_symbol_col != "HGNC") {
    if (gene_symbol_col %in% names(rna_tbl) && !"HGNC" %in% names(rna_tbl)) {
      rna_tbl <- dplyr::rename(rna_tbl, HGNC = dplyr::all_of(gene_symbol_col))
    }
  }

  if (is.null(fp_aligned)) {
    cache_dir <- fp_cache_dir
    if (is.null(cache_dir)) {
      base_dir_cfg <- get_cfg("base_dir")
      if (is.character(base_dir_cfg) && nzchar(base_dir_cfg)) {
        cache_dir <- file.path(base_dir_cfg, "cache")
      }
    }
    cache_tag <- fp_cache_tag
    if (is.null(cache_tag)) {
      cache_tag <- get_cfg("db")
    }
    if (is.null(cache_dir) || is.null(cache_tag)) {
      .log_abort("`fp_aligned` is NULL and no cache_dir/cache_tag could be resolved. Provide fp_aligned or set base_dir/db in config.")
    }
    fp_aligned <- load_fp_aligned_from_cache(cache_dir = cache_dir, cache_tag = cache_tag, output_mode = "distinct", verbose = verbose)
  }

  atac_out <- load_atac(atac_data, sort_peaks = TRUE)

  grn_set <- build_grn_set(
    fp_score      = fp_aligned$fp_score,
    fp_bound      = fp_aligned$fp_bound,
    fp_annotation = fp_aligned$fp_annotation,
    atac_score    = atac_out$score,
    atac_overlap  = atac_out$overlap,
    rna           = rna_tbl,
    metadata      = metadata,
    tf_list       = tf_list,
    motif_db      = motif_db,
    label_col     = label_col,
    expected_n    = expected_n
  )

  if (isTRUE(verbose)) {
    .log_inform("Generating condition-specific footprint bound matrix.")
  }
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
  if (isTRUE(verbose)) {
    .log_inform("Quantile-normalizing footprint scores across conditions.")
  }
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  if (isTRUE(verbose)) {
    .log_inform("Generating condition-specific gene expression flags.")
  }
  grn_set <- grn_add_rna_expressed(
    grn_set,
    label_col = label_col,
    threshold_gene_expr = threshold_gene_expr
  )
  if (isTRUE(verbose)) {
    .log_inform("Combined multi-omic data object ready for downstream analysis.")
  }

  if (isTRUE(write_outputs)) {
    step1_out_dir <- file.path(get_cfg("base_dir"), "predict_tf_binding_sites")
    write_grn_outputs(
      grn_set,
      out_dir = step1_out_dir,
      db = get_cfg("db"),
      qn_base_dir = file.path(step1_out_dir, "cache")
    )
    plot_fp_norm_bound_qc(
      omics_data = grn_set,
      out_dir = step1_out_dir,
      db = get_cfg("db"),
      threshold_fp_score = threshold_fp_score,
      max_points = 100000L,
      verbose = verbose
    )
    plot_gene_expr_qc(
      omics_data = grn_set,
      out_dir = step1_out_dir,
      db = get_cfg("db"),
      threshold_gene_expr = threshold_gene_expr,
      verbose = verbose
    )
    saveRDS(
      grn_set,
      file.path(step1_out_dir, sprintf("01_multiomic_data_object_%s.rds", get_cfg("db")))
    )
  }

  grn_set
}
