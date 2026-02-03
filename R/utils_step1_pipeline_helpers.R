#' Module 1 pipeline helpers
#'
#' Lightweight helpers for assembling multi-omic inputs in Module 1.
#'
#' @name module1_pipeline_helpers
#' @rdname module1_pipeline_helpers
#' @noRd
NULL

#' Load and assemble a combined multi-omic data object
#'
#' Runs ATAC preprocessing, builds a GRN-aligned object, and adds
#' binding/normalization summaries followed by RNA expression flags.
#'
#' @param fp_aligned Aligned footprint object returned by align_footprints().
#' @param atac_data ATAC master table (tibble/data.frame). If NULL, load from
#'   \code{atac_data_path} or YAML config (\code{atac_master} preferred).
#' @param rna_tbl RNA table with columns already matched to sample IDs. If NULL,
#'   load from \code{rna_path} or YAML config (\code{rna_mapped} preferred).
#' @param metadata Sample metadata with id column. If NULL, load from
#'   \code{metadata_path} or YAML config (\code{sample_metadata} preferred).
#' @param atac_data_path Optional CSV path for ATAC master table.
#' @param rna_path Optional CSV path for RNA (mapped) table.
#' @param metadata_path Optional CSV path for sample metadata.
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
    fp_aligned,
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
  if (isTRUE(verbose)) {
    .log_inform("Process: Load input data and create a combined multiomic data object.")
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
    .log_inform("2.2 Generate a binary, condition-specific footprint binding flag matrix (bound=1; unbound=0).")
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
    .log_inform("2.3 Perform cross-condition quantile normalization on bound footprints.")
  }
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID")
  if (isTRUE(verbose)) {
    .log_inform("2.4 Generate a binary, condition-specific gene expression flag matrix.")
  }
  grn_set <- grn_add_rna_expressed(
    grn_set,
    label_col = label_col,
    threshold_gene_expr = threshold_gene_expr
  )
  if (isTRUE(verbose)) {
    .log_inform("Output: Combined multi-omic data object ready for downstream analysis.")
  }

  grn_set
}
