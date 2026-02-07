#' Build a GRN-aligned set from pre-cleaned inputs (general, package-ready)
#'
#' @param fp_score      tibble with 'peak_ID' + ATAC id columns (numeric)
#' @param fp_bound      tibble with 'peak_ID' + ATAC id columns (0/1 or logical)
#' @param fp_annotation tibble with motif hits (passed through)
#' @param atac_score    tibble with 'atac_peak' + ATAC id columns (numeric)
#' @param atac_overlap  tibble with 'atac_peak' + ATAC id columns (0/1 or logical)
#' @param rna           tibble with 'ensembl_gene_id','HGNC' + RNA columns (already named by ATAC ids)
#' @param metadata      tibble with at least an 'id' column; order defines column order
#' @param tf_list       character vector of TF symbols (passed through)
#' @param motif_db      optional tibble (passed through)
#' @param label_col     optional metadata column to carry along as sample labels (e.g. "strict_match_rna" or "cell_stress_type")
#' @param expected_n    optional integer for info message
#' @param filter_atac_by_fp_annotation Logical; if TRUE, subset atac_score/atac_overlap
#'   rows to peaks present in \code{fp_annotation}. Default TRUE.
#' @return list(fp_score, fp_bound, atac_score, atac_overlap, rna, fp_annotation, tf_list, motif_db,
#'              sample_metadata_used, sample_labels, dropped_ids)
#' @export
build_grn_set <- function(
    fp_score,
    fp_bound,
    fp_annotation,
    atac_score,
    atac_overlap,
    rna,
    metadata,
    tf_list,
    motif_db   = NULL,
    label_col  = NULL,
    expected_n = NULL,
    filter_atac_by_fp_annotation = TRUE
) {
  # ---- checks
  if (!all(c("peak_ID") %in% names(fp_score))) .log_abort("`fp_score` needs 'peak_ID'.")
  if (!all(c("peak_ID") %in% names(fp_bound))) .log_abort("`fp_bound` needs 'peak_ID'.")
  if (!"atac_peak" %in% names(atac_score))     .log_abort("`atac_score` needs 'atac_peak'.")
  if (!"atac_peak" %in% names(atac_overlap))   .log_abort("`atac_overlap` needs 'atac_peak'.")
  if (!all(c("ensembl_gene_id","HGNC") %in% names(rna)))
    .log_abort("`rna` needs 'ensembl_gene_id' and 'HGNC'.")
  if (!"id" %in% names(metadata))
    .log_abort("`metadata` needs an 'id' column.")
  if (!is.null(label_col) && !label_col %in% names(metadata))
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  if (!is.character(tf_list)) .log_abort("`tf_list` must be a character vector.")

  # ---- IDs (preserve metadata order; drop NA/dups)
  ids <- metadata$id
  ids <- ids[!is.na(ids)]
  ids <- ids[!duplicated(ids)]
  if (!length(ids)) .log_abort("No usable sample ids found in `metadata$id`.")

  # ---- availability & subsetting
  miss_fp_score <- setdiff(ids, names(fp_score))
  miss_fp_bound <- setdiff(ids, names(fp_bound))
  miss_atac_sc  <- setdiff(ids, names(atac_score))
  miss_atac_ol  <- setdiff(ids, names(atac_overlap))
  miss_rna      <- setdiff(ids, names(rna))

  if (length(miss_fp_score)) message("Dropping ", length(miss_fp_score), " id(s) missing in `fp_score`: ", paste(miss_fp_score, collapse=", "))
  if (length(miss_fp_bound)) message("Dropping ", length(miss_fp_bound), " id(s) missing in `fp_bound`: ", paste(miss_fp_bound, collapse=", "))
  if (length(miss_atac_sc))  message("`atac_score` missing ", length(miss_atac_sc), " id(s): ", paste(miss_atac_sc, collapse=", "))
  if (length(miss_atac_ol))  message("`atac_overlap` missing ", length(miss_atac_ol), " id(s): ", paste(miss_atac_ol, collapse=", "))
  if (length(miss_rna))      message("`rna` missing ", length(miss_rna), " id(s): ", paste(miss_rna, collapse=", "))

  keep_ids <- ids[ids %in% names(fp_score) &
                    ids %in% names(fp_bound) &
                    ids %in% names(rna)]

  if (!length(keep_ids)) .log_abort("No overlapping ids across fp_score/fp_bound/rna after checks.")

  if (!is.null(expected_n) && length(keep_ids) != expected_n) {
    message("Note: aligned sample count = ", length(keep_ids),
            " (expected ", expected_n, "). Proceeding.")
  }


  if (isTRUE(filter_atac_by_fp_annotation)) {
    # subset atac_* rows to peaks present in filtered fp_annotation
    peaks_keep_atac <- unique(fp_annotation$atac_peak)
    peaks_keep_atac <- peaks_keep_atac[!is.na(peaks_keep_atac)]

    if (length(peaks_keep_atac) == 0L) {
      .log_abort("No usable 'atac_peak' values found in `fp_annotation` after filtering.")
    }

    # subset rows in atac_* to only those peaks (do this BEFORE column selection)
    atac_score   <- dplyr::semi_join(atac_score,   tibble::tibble(atac_peak = peaks_keep_atac), by = "atac_peak")
    atac_overlap <- dplyr::semi_join(atac_overlap, tibble::tibble(atac_peak = peaks_keep_atac), by = "atac_peak")

    # optional: brief info
    message(sprintf("ATAC rows kept by fp_annotation peaks: score=%s, overlap=%s",
                    nrow(atac_score), nrow(atac_overlap)))
  }


  # ---- ordered subsets
  fp_score_sub    <- dplyr::select(fp_score, "peak_ID", dplyr::all_of(keep_ids))
  fp_bound_sub    <- dplyr::select(fp_bound, "peak_ID", dplyr::all_of(keep_ids))
  atac_score_sub  <- dplyr::select(atac_score, "atac_peak", dplyr::all_of(intersect(keep_ids, names(atac_score))))
  atac_overlap_sub<- dplyr::select(atac_overlap, "atac_peak", dplyr::all_of(intersect(keep_ids, names(atac_overlap))))
  rna_sub         <- dplyr::select(rna, "ensembl_gene_id","HGNC", dplyr::all_of(keep_ids))

  # ---- metadata used / labels
  meta_used <- dplyr::semi_join(metadata, tibble::tibble(id = keep_ids), by = "id")
  # keep metadata order
  meta_used <- meta_used[match(keep_ids, meta_used$id), , drop = FALSE]
  sample_labels <- if (!is.null(label_col)) meta_used[[label_col]] else NULL

  grn_set <- list(
    fp_score             = fp_score_sub,
    fp_bound             = fp_bound_sub,
    atac_score           = atac_score_sub,
    atac_overlap         = atac_overlap_sub,
    rna                  = rna_sub,
    fp_annotation        = fp_annotation,
    tf_list              = tf_list,
    motif_db             = motif_db,
    sample_metadata_used = meta_used,
    sample_labels        = sample_labels,
    dropped_ids          = setdiff(ids, keep_ids)
  )

  grn_set <- grn_status_init(grn_set)
  grn_status_set(grn_set, "built")
}

grn_status_init <- function(grn_set) {
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (is.null(grn_set$status) || !is.list(grn_set$status)) {
    grn_set$status <- list()
  }
  grn_set
}

grn_status_is <- function(grn_set, key) {
  is.list(grn_set) && is.list(grn_set$status) && isTRUE(grn_set$status[[key]])
}

grn_status_set <- function(grn_set, key, value = TRUE) {
  grn_set <- grn_status_init(grn_set)
  grn_set$status[[key]] <- isTRUE(value)
  grn_set
}

assert_fp_alignment <- function(
    grn_set,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak"
) {
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (!all(c("fp_score", "fp_bound", "fp_annotation") %in% names(grn_set))) {
    .log_abort("`grn_set` must include fp_score, fp_bound, and fp_annotation.")
  }
  score <- grn_set$fp_score
  bound <- grn_set$fp_bound
  annot <- grn_set$fp_annotation
  if (!is.data.frame(score) || !is.data.frame(bound) || !is.data.frame(annot)) {
    .log_abort("fp_score, fp_bound, and fp_annotation must be data frames.")
  }
  if (!all(c(score_key) %in% names(score))) {
    .log_abort("fp_score is missing column {.val {score_key}}.")
  }
  if (!all(c(bound_key) %in% names(bound))) {
    .log_abort("fp_bound is missing column {.val {bound_key}}.")
  }
  if (!all(c(annot_key) %in% names(annot))) {
    .log_abort("fp_annotation is missing column {.val {annot_key}}.")
  }
  if (nrow(score) != nrow(bound) || nrow(score) != nrow(annot)) {
    .log_abort(
      "Row counts differ: fp_score={nrow(score)}, fp_bound={nrow(bound)}, fp_annotation={nrow(annot)}."
    )
  }
  score_ids <- score[[score_key]]
  bound_ids <- bound[[bound_key]]
  annot_ids <- annot[[annot_key]]
  if (!identical(score_ids, bound_ids)) {
    n_bad <- sum(score_ids != bound_ids, na.rm = TRUE)
    .log_abort("fp_score and fp_bound order mismatch (n_mismatch={n_bad}).")
  }
  if (!identical(score_ids, annot_ids)) {
    n_bad <- sum(score_ids != annot_ids, na.rm = TRUE)
    .log_abort("fp_score and fp_annotation order mismatch (n_mismatch={n_bad}).")
  }
  invisible(TRUE)
}

#' Load motif database and derive TF list
#'
#' @param db Character; supported values: "JASPAR2024", "HOCOMOCOv13".
#' @param ref_genome Character; "hg38" or "mm10" (used for JASPAR2024 filenames).
#' @return List with `motif_db` tibble and `tf_list` character vector.
#' @export
init_motif_db <- function(db, ref_genome = NULL) {
  if (!is.character(db) || length(db) != 1L || !nzchar(db)) {
    .log_abort("`db` must be a non-empty character scalar.")
  }
  motif_path <- resolve_motif_db_path(db, ref_genome = ref_genome)
  motif_db <- readr::read_tsv(motif_path, show_col_types = FALSE)
  if ("HGNC" %in% names(motif_db) && !("gene_symbol" %in% names(motif_db))) {
    motif_db <- dplyr::rename(motif_db, gene_symbol = HGNC)
  }
  if (!all(c("motif", "gene_symbol") %in% names(motif_db))) {
    .log_abort("Motif DB must include columns: motif, gene_symbol.")
  }
  motif_db <- motif_db |>
    dplyr::select(motif, gene_symbol)

  tf_list <- motif_db |>
    tidyr::separate_rows(gene_symbol, sep = "::") |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
    dplyr::distinct(gene_symbol) |>
    dplyr::pull(gene_symbol)
  list(motif_db = motif_db, tf_list = tf_list)
}

.motif_gene_col <- function(motif_db) {
  if (!is.data.frame(motif_db)) return(NULL)
  if ("gene_symbol" %in% names(motif_db)) return("gene_symbol")
  if ("HGNC" %in% names(motif_db)) return("HGNC")
  NULL
}

#' Resolve motif DB path by db and reference genome
#'
#' @param db Character; "JASPAR2024" or "HOCOMOCOv13".
#' @param ref_genome Character; "hg38" or "mm10".
#' @return Character path to motif DB file.
#' @export
resolve_motif_db_path <- function(db, ref_genome = NULL) {
  if (is.null(ref_genome) || !nzchar(ref_genome)) {
    ref_genome <- get0("ref_genome", envir = .GlobalEnv, inherits = FALSE)
  }
  if (is.null(ref_genome) || !nzchar(ref_genome)) ref_genome <- "hg38"
  ref_genome <- tolower(as.character(ref_genome))
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty character scalar.")
  }
  db <- as.character(db)
  if (db == "JASPAR2024") {
    if (ref_genome == "mm10") {
      cand <- c("inst/extdata/genome/JASPAR2024_mm10.txt",
                "inst/extdata/genome/JASPAR2024_mm.txt",
                "inst/extdata/genome/JASPAR2024.txt")
    } else {
      cand <- c("inst/extdata/genome/JASPAR2024_hg38.txt",
                "inst/extdata/genome/JASPAR2024.txt")
    }
    path <- cand[file.exists(cand)][1]
  } else if (db == "HOCOMOCOv13") {
    if (ref_genome == "mm10") {
      cand <- c("inst/extdata/genome/HOCOMOCOv13_mm10.txt",
                "inst/extdata/genome/HOCOMOCOv13.txt")
    } else {
      cand <- c("inst/extdata/genome/HOCOMOCOv13_hg38.txt",
                "inst/extdata/genome/HOCOMOCOv13.txt")
    }
    path <- cand[file.exists(cand)][1]
  } else {
    .log_abort("Unsupported db {.val {db}}. Expected 'JASPAR2024' or 'HOCOMOCOv13'.")
  }
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    .log_abort("Motif DB file not found for db={db}, ref_genome={ref_genome}.")
  }
  path
}

#' Build per-condition gene expression flags (1/0)
#'
#' @param rna_tbl Tibble with ensembl_gene_id, HGNC + sample columns.
#' @param metadata Tibble with at least id and label_col.
#' @param label_col Metadata column holding condition labels.
#' @param threshold_gene_expr Numeric threshold for gene expression.
#' @param min_samples Integer; minimum number of samples per condition required to flag.
#' @return Tibble with ensembl_gene_id, HGNC + condition columns (0/1).
#' @export
make_rna_expressed <- function(
    rna_tbl,
    metadata,
    label_col,
    threshold_gene_expr = 10,
    min_samples = 1L
) {
  if (!is.data.frame(rna_tbl)) .log_abort("`rna_tbl` must be a data.frame.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl))) {
    .log_abort("`rna_tbl` must include ensembl_gene_id and HGNC.")
  }
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) {
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  }
  if (!is.numeric(threshold_gene_expr) || length(threshold_gene_expr) != 1L) {
    .log_abort("`threshold_gene_expr` must be a single numeric value.")
  }
  min_samples <- as.integer(min_samples)
  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 1L) {
    .log_abort("`min_samples` must be a positive integer.")
  }

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(rna_tbl))
  if (!length(ids_present)) {
    .log_abort("No overlap between metadata ids and RNA columns.")
  }
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  out <- tibble::tibble(
    ensembl_gene_id = rna_tbl$ensembl_gene_id,
    HGNC = rna_tbl$HGNC
  )

  for (cond in cond_levels) {
    ids_cond <- meta_use$id[meta_use[[label_col]] == cond]
    ids_cond <- intersect(ids_cond, names(rna_tbl))
    if (!length(ids_cond)) {
      .log_warn("No RNA columns found for condition {.val {cond}}.")
      out[[cond]] <- integer(nrow(rna_tbl))
      next
    }
    mat <- as.matrix(rna_tbl[, ids_cond, drop = FALSE])
    flag <- rowSums(mat >= threshold_gene_expr, na.rm = TRUE) >= min_samples
    out[[cond]] <- as.integer(flag)
  }
  out
}

#' Aggregate RNA expression by condition
#'
#' @param rna_tbl Tibble with ensembl_gene_id, HGNC + sample columns.
#' @param metadata Tibble with id and label_col.
#' @param label_col Metadata column holding condition labels.
#' @param agg_fun Aggregate function: "mean" or "median".
#' @return Tibble with ensembl_gene_id, HGNC + condition columns.
#'   If each condition has a single sample, values are taken directly
#'   from that sample (no aggregation).
#' @export
make_rna_condition <- function(
    rna_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(rna_tbl)) .log_abort("`rna_tbl` must be a data.frame.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl))) {
    .log_abort("`rna_tbl` must include ensembl_gene_id and HGNC.")
  }
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) {
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  }
  agg_fun <- match.arg(agg_fun)

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(rna_tbl))
  if (!length(ids_present)) {
    .log_abort("No overlap between metadata ids and RNA columns.")
  }
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  out <- tibble::tibble(
    ensembl_gene_id = rna_tbl$ensembl_gene_id,
    HGNC = rna_tbl$HGNC
  )

  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, names(rna_tbl))
      if (!length(ids_cond)) {
        .log_warn("No RNA columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      out[[cond]] <- rna_tbl[[ids_cond[[1L]]]]
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, names(rna_tbl))
      if (!length(ids_cond)) {
        .log_warn("No RNA columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      mat <- as.matrix(rna_tbl[, ids_cond, drop = FALSE])
      if (agg_fun == "median") {
        out[[cond]] <- apply(mat, 1L, stats::median, na.rm = TRUE)
      } else {
        out[[cond]] <- rowMeans(mat, na.rm = TRUE)
      }
    }
  }
  out
}

#' Aggregate fp_score by condition
#'
#' @param fp_score_tbl Tibble with peak_ID + sample columns.
#' @param metadata Tibble with id and label_col.
#' @param label_col Metadata column holding condition labels.
#' @param agg_fun Aggregate function: "mean" or "median".
#' @return Tibble with peak_ID + condition columns.
#'   If each condition has a single sample, values are taken directly
#'   from that sample (no aggregation).
#' @export
make_fp_score_condition <- function(
    fp_score_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(fp_score_tbl)) .log_abort("`fp_score_tbl` must be a data.frame.")
  if (!"peak_ID" %in% names(fp_score_tbl)) .log_abort("`fp_score_tbl` needs a peak_ID column.")
  if (!is.data.frame(metadata)) .log_abort("`metadata` must be a data.frame.")
  if (!"id" %in% names(metadata)) .log_abort("`metadata` must include an id column.")
  if (!label_col %in% names(metadata)) {
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  }
  agg_fun <- match.arg(agg_fun)

  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_present <- intersect(meta_use$id, names(fp_score_tbl))
  if (!length(ids_present)) {
    .log_abort("No overlap between metadata ids and fp_score columns.")
  }
  meta_use <- meta_use[meta_use$id %in% ids_present, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")

  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  out <- tibble::tibble(peak_ID = fp_score_tbl$peak_ID)
  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, names(fp_score_tbl))
      if (!length(ids_cond)) {
        .log_warn("No fp_score columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      out[[cond]] <- fp_score_tbl[[ids_cond[[1L]]]]
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, names(fp_score_tbl))
      if (!length(ids_cond)) {
        .log_warn("No fp_score columns found for condition {.val {cond}}.")
        out[[cond]] <- NA_real_
        next
      }
      mat <- as.matrix(fp_score_tbl[, ids_cond, drop = FALSE])
      if (agg_fun == "median") {
        out[[cond]] <- apply(mat, 1L, stats::median, na.rm = TRUE)
      } else {
        out[[cond]] <- rowMeans(mat, na.rm = TRUE)
      }
    }
  }
  out
}

#' Aggregate ATAC score by condition
#'
#' @param atac_score_tbl Tibble with atac_peak + sample columns.
#' @param metadata Tibble with id and label_col.
#' @param label_col Metadata column holding condition labels.
#' @param agg_fun Aggregate function: "mean" or "median".
#' @return Tibble with atac_peak + condition columns.
#'   If each condition has a single sample, values are taken directly
#'   from that sample (no aggregation).
#' @export
make_atac_score_condition <- function(
    atac_score_tbl,
    metadata,
    label_col,
    agg_fun = c("mean", "median")
) {
  if (!is.data.frame(atac_score_tbl)) .log_abort("`atac_score_tbl` must be a data.frame.")
  if (!"atac_peak" %in% names(atac_score_tbl)) .log_abort("`atac_score_tbl` needs an atac_peak column.")
  atac_tmp <- dplyr::rename(atac_score_tbl, peak_ID = atac_peak)
  out <- make_fp_score_condition(
    fp_score_tbl = atac_tmp,
    metadata = metadata,
    label_col = label_col,
    agg_fun = agg_fun
  )
  dplyr::rename(out, atac_peak = peak_ID)
}

#' Build per-condition FP bound flags (1/0) with score and ATAC overlap gating
#'
#' @param fp_bound_tbl Tibble with peak_ID + sample columns (0/1).
#' @param fp_score_tbl Tibble with peak_ID + sample columns (numeric).
#' @param atac_overlap_tbl Tibble with atac_peak + sample columns (0/1).
#' @param fp_annotation_tbl Tibble with fp_peak and atac_peak columns.
#' @param metadata Tibble with id and label_col.
#' @param label_col Metadata column holding condition labels.
#' @param threshold_fp_score Numeric; require fp_score >= threshold to keep bound.
#' @param min_samples Integer; minimum number of samples per condition required to flag.
#' @return Tibble with peak_ID + condition columns (0/1). If each condition has
#'   a single sample, values are taken directly from that sample (no aggregation).
#' @export
make_fp_bound_condition <- function(
    fp_bound_tbl,
    fp_score_tbl,
    atac_overlap_tbl,
    fp_annotation_tbl,
    metadata,
    label_col,
    threshold_fp_score = 2,
    min_samples = 1L
) {
  if (!all(c("peak_ID") %in% names(fp_bound_tbl))) .log_abort("`fp_bound_tbl` needs peak_ID.")
  if (!all(c("peak_ID") %in% names(fp_score_tbl))) .log_abort("`fp_score_tbl` needs peak_ID.")
  if (!all(c("atac_peak") %in% names(atac_overlap_tbl))) .log_abort("`atac_overlap_tbl` needs atac_peak.")
  if (!all(c("fp_peak", "atac_peak") %in% names(fp_annotation_tbl))) {
    .log_abort("`fp_annotation_tbl` needs fp_peak and atac_peak.")
  }
  if (!is.data.frame(metadata) || !"id" %in% names(metadata)) {
    .log_abort("`metadata` must include an id column.")
  }
  if (!label_col %in% names(metadata)) {
    .log_abort("label_col {.val {label_col}} not found in `metadata`.")
  }
  min_samples <- as.integer(min_samples)
  if (!is.numeric(min_samples) || length(min_samples) != 1L || min_samples < 1L) {
    .log_abort("`min_samples` must be a positive integer.")
  }
  if (!is.numeric(threshold_fp_score) || length(threshold_fp_score) != 1L) {
    .log_abort("`threshold_fp_score` must be a single numeric value.")
  }

  ids_bound <- setdiff(names(fp_bound_tbl), "peak_ID")
  ids_score <- setdiff(names(fp_score_tbl), "peak_ID")
  ids_atac  <- setdiff(names(atac_overlap_tbl), "atac_peak")
  meta_use <- metadata[, c("id", label_col), drop = FALSE]
  meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
  meta_use <- meta_use[!duplicated(meta_use$id), , drop = FALSE]
  ids_use <- intersect(meta_use$id, intersect(ids_bound, intersect(ids_score, ids_atac)))
  if (!length(ids_use)) {
    .log_abort("No overlapping ids across fp_bound/fp_score/atac_overlap and metadata.")
  }
  meta_use <- meta_use[meta_use$id %in% ids_use, , drop = FALSE]
  cond_levels <- unique(meta_use[[label_col]])
  if (!length(cond_levels)) .log_abort("No condition labels found in metadata.")
  meta_use[[label_col]] <- factor(meta_use[[label_col]], levels = cond_levels)
  cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)
  cond_sizes <- vapply(cond_groups, length, integer(1L))
  single_sample <- all(cond_sizes == 1L)

  # align fp_score rows to fp_bound by peak_ID if needed
  if (!identical(fp_bound_tbl$peak_ID, fp_score_tbl$peak_ID)) {
    idx <- match(fp_bound_tbl$peak_ID, fp_score_tbl$peak_ID)
    if (anyNA(idx)) .log_abort("fp_score rows do not align with fp_bound by peak_ID.")
    score_tbl <- fp_score_tbl[idx, , drop = FALSE]
  } else {
    score_tbl <- fp_score_tbl
  }

  bound_mat <- as.matrix(fp_bound_tbl[, ids_use, drop = FALSE])
  score_mat <- as.matrix(score_tbl[, ids_use, drop = FALSE])
  bound_mat[is.na(bound_mat)] <- 0L
  score_mat[is.na(score_mat)] <- NA_real_
  bound_mat <- bound_mat > 0 & score_mat >= threshold_fp_score
  storage.mode(bound_mat) <- "integer"

  # map fp_peak -> atac_peak and align atac_overlap to fp_bound
  map_fp_atac <- dplyr::distinct(fp_annotation_tbl, .data$fp_peak, .data$atac_peak)
  dup <- map_fp_atac |>
    dplyr::count(.data$fp_peak, name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup)) .log_abort("Each fp_peak must map to a single atac_peak.")

  idx_fp <- match(fp_bound_tbl$peak_ID, map_fp_atac$fp_peak)
  if (anyNA(idx_fp)) .log_abort("Some fp_bound peak_ID not found in fp_annotation fp_peak.")
  atac_peaks <- map_fp_atac$atac_peak[idx_fp]
  idx_atac <- match(atac_peaks, atac_overlap_tbl$atac_peak)
  if (anyNA(idx_atac)) .log_abort("Mapped atac_peak missing in atac_overlap.")
  atac_tbl <- atac_overlap_tbl[idx_atac, ids_use, drop = FALSE]
  atac_mat <- as.matrix(atac_tbl)
  atac_mat[is.na(atac_mat)] <- 0L
  atac_mat <- atac_mat > 0

  out <- tibble::tibble(peak_ID = fp_bound_tbl$peak_ID)
  if (single_sample) {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, ids_use)
      if (!length(ids_cond)) {
        .log_warn("No fp_bound columns found for condition {.val {cond}}.")
        out[[cond]] <- integer(nrow(fp_bound_tbl))
        next
      }
      id <- ids_cond[[1L]]
      out[[cond]] <- as.integer(bound_mat[, id, drop = FALSE] > 0 & atac_mat[, id, drop = FALSE] > 0)
    }
  } else {
    for (cond in cond_levels) {
      ids_cond <- cond_groups[[cond]]
      ids_cond <- intersect(ids_cond, ids_use)
      if (!length(ids_cond)) {
        .log_warn("No fp_bound columns found for condition {.val {cond}}.")
        out[[cond]] <- integer(nrow(fp_bound_tbl))
        next
      }
      b_ok <- rowSums(bound_mat[, ids_cond, drop = FALSE], na.rm = TRUE) >= min_samples
      a_ok <- rowSums(atac_mat[, ids_cond, drop = FALSE], na.rm = TRUE) > 0L
      out[[cond]] <- as.integer(b_ok & a_ok)
    }
  }
  out
}

#' Select peaks bound in at least N conditions
#'
#' @param fp_bound_tbl Tibble with peak_ID + condition/sample columns (0/1 or logical).
#' @param min_bound Integer threshold for the number of bound conditions.
#' @param samples Optional character vector of columns to consider (default: all).
#' @param na_as_unbound Logical; treat NA as unbound (0) if TRUE.
#' @param use_parallel Logical; if TRUE, use multicore chunking when available.
#' @param n_workers Integer; number of workers for multicore (default: getOption("mc.cores", 1L)).
#' @return Character vector of peak_ID values meeting the threshold.
#' @export
get_min_bound_peaks <- function(
    fp_bound_tbl,
    min_bound = 1L,
    samples = NULL,
    na_as_unbound = TRUE,
    use_parallel = FALSE,
    n_workers = getOption("mc.cores", 1L)
) {
  if (!is.data.frame(fp_bound_tbl)) .log_abort("`fp_bound_tbl` must be a data.frame.")
  if (!"peak_ID" %in% names(fp_bound_tbl)) .log_abort("`fp_bound_tbl` needs a peak_ID column.")
  min_bound <- as.integer(min_bound)
  if (length(min_bound) != 1L || is.na(min_bound) || min_bound < 1L) {
    .log_abort("`min_bound` must be a positive integer.")
  }

  all_ids <- setdiff(names(fp_bound_tbl), "peak_ID")
  if (is.null(samples)) {
    use_ids <- all_ids
  } else {
    miss <- setdiff(samples, all_ids)
    if (length(miss)) {
      .log_abort("Some `samples` not found in fp_bound: {paste(miss, collapse = ', ')}.")
    }
    use_ids <- samples
  }
  if (!length(use_ids)) .log_abort("No sample columns selected to evaluate binding.")

  mat_raw <- as.matrix(fp_bound_tbl[, use_ids, drop = FALSE])
  if (isTRUE(na_as_unbound)) {
    mat_raw[is.na(mat_raw)] <- 0L
  }
  mat <- mat_raw > 0

  if (isTRUE(use_parallel) && .Platform$OS.type == "windows") {
    .log_warn("Multicore is not supported on Windows; using single-core rowSums().")
    use_parallel <- FALSE
  }

  if (isTRUE(use_parallel) && n_workers > 1L) {
    n_rows <- nrow(mat)
    if (!n_rows) return(character(0))
    chunk_size <- ceiling(n_rows / n_workers)
    chunks <- split(seq_len(n_rows), ceiling(seq_len(n_rows) / chunk_size))
    sums <- unlist(
      parallel::mclapply(
        chunks,
        function(idx) rowSums(mat[idx, , drop = FALSE]),
        mc.cores = n_workers
      ),
      use.names = FALSE
    )
  } else {
    sums <- rowSums(mat)
  }

  if (!isTRUE(na_as_unbound)) {
    any_na <- apply(is.na(mat_raw), 1L, any)
    keep <- !any_na & sums >= min_bound
  } else {
    keep <- sums >= min_bound
  }

  fp_bound_tbl$peak_ID[keep]
}

#' GRN set helpers for condition summaries and outputs
#'
#' These helpers add condition-level summaries (RNA and footprint),
#' compute filtered TF-annotation correlations, and write out
#' standard GRN outputs.
#'
#' @param grn_set A list returned by \code{build_grn_set()}.
#' @param label_col Metadata column used to define conditions.
#' @param threshold_gene_expr Expression threshold for gene activity.
#' @param min_samples Minimum samples meeting a threshold.
#' @param agg_fun Aggregation function for condition summaries.
#' @param threshold_fp_score Footprint score threshold for bound status.
#' @param min_bound Minimum number of bound samples per peak to keep.
#' @param use_parallel Use parallel computation when available.
#' @param id_col Identifier column for quantile normalization.
#' @param method Correlation method.
#' @param cores Number of CPU cores to use.
#' @param chunk_size Row chunk size for parallel processing.
#' @param min_non_na Minimum non-missing values required for correlations.
#' @param r_thr Correlation threshold.
#' @param p_thr Adjusted P-value threshold.
#' @param set_active Whether to set active flags after filtering.
#' @param output_bed Optional BED output path.
#' @param output_bed_condition Optional per-condition BED output path.
#' @param label_col Metadata column used for per-condition outputs.
#' @param out_dir Output directory for writing artifacts.
#' @param db Motif database label.
#' @param qn_base_dir Base directory for quantile-normalized outputs.
#' @param force Recompute even if already present.
#' @param verbose Emit progress messages.
#' @return Updated \code{grn_set} or \code{invisible(TRUE)} for writer helpers.
#' @rdname grn_set_helpers
#' @export
grn_add_rna_expressed <- function(
    grn_set,
    label_col,
    threshold_gene_expr = 10,
    min_samples = 1L,
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "rna_expressed") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$rna)) .log_abort("`grn_set$rna` is missing or invalid.")
  grn_set$rna_expressed <- make_rna_expressed(
    rna_tbl = grn_set$rna,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    threshold_gene_expr = threshold_gene_expr,
    min_samples = min_samples
  )
  if (isTRUE(verbose)) {
    n_genes <- nrow(grn_set$rna_expressed)
    n_cond <- max(ncol(grn_set$rna_expressed) - 2L, 0L)
    n_flag <- sum(as.matrix(grn_set$rna_expressed[, -(1:2), drop = FALSE]) > 0L, na.rm = TRUE)
    .log_inform(
      "RNA expressed flags: {n_genes} genes × {n_cond} conditions; {n_flag} total positives."
    )
  }
  grn_status_set(grn_set, "rna_expressed")
}

#' @rdname grn_set_helpers
#' @export
grn_add_rna_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "rna_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$rna)) .log_abort("`grn_set$rna` is missing or invalid.")
  grn_set$rna_condition <- make_rna_condition(
    rna_tbl = grn_set$rna,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  if (isTRUE(verbose)) {
    n_genes <- nrow(grn_set$rna_condition)
    n_cond <- max(ncol(grn_set$rna_condition) - 2L, 0L)
    .log_inform("RNA by condition: {n_genes} genes x {n_cond} conditions.")
  }
  grn_status_set(grn_set, "rna_condition")
}

#' @rdname grn_set_helpers
#' @export
grn_add_fp_score_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_score_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_score)) .log_abort("`grn_set$fp_score` is missing or invalid.")
  grn_set$fp_score_condition <- make_fp_score_condition(
    fp_score_tbl = grn_set$fp_score,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  if (isTRUE(verbose)) {
    n_peaks <- nrow(grn_set$fp_score_condition)
    n_cond <- max(ncol(grn_set$fp_score_condition) - 1L, 0L)
    .log_inform(
      "FP score by condition: {n_peaks} peaks × {n_cond} conditions."
    )
  }
  grn_status_set(grn_set, "fp_score_condition")
}

#' @rdname grn_set_helpers
#' @export
grn_add_atac_score_condition <- function(
    grn_set,
    label_col,
    agg_fun = c("mean", "median"),
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "atac_score_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$atac_score)) .log_abort("`grn_set$atac_score` is missing or invalid.")
  grn_set$atac_score_condition <- make_atac_score_condition(
    atac_score_tbl = grn_set$atac_score,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    agg_fun = agg_fun
  )
  if (isTRUE(verbose)) {
    n_peaks <- nrow(grn_set$atac_score_condition)
    n_cond <- max(ncol(grn_set$atac_score_condition) - 1L, 0L)
    .log_inform("ATAC score by condition: {n_peaks} peaks × {n_cond} conditions.")
  }
  grn_status_set(grn_set, "atac_score_condition")
}

#' @rdname grn_set_helpers
#' @export
grn_add_fp_bound_condition <- function(
    grn_set,
    label_col,
    threshold_fp_score = 2,
    min_samples = 1L,
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_bound_condition") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_bound)) .log_abort("`grn_set$fp_bound` is missing or invalid.")
  grn_set$fp_bound_condition <- make_fp_bound_condition(
    fp_bound_tbl = grn_set$fp_bound,
    fp_score_tbl = grn_set$fp_score,
    atac_overlap_tbl = grn_set$atac_overlap,
    fp_annotation_tbl = grn_set$fp_annotation,
    metadata = grn_set$sample_metadata_used,
    label_col = label_col,
    threshold_fp_score = threshold_fp_score,
    min_samples = min_samples
  )
  if (isTRUE(verbose)) {
    n_peaks <- nrow(grn_set$fp_bound_condition)
    n_cond <- max(ncol(grn_set$fp_bound_condition) - 1L, 0L)
    mat <- as.matrix(grn_set$fp_bound_condition[, -1, drop = FALSE])
    total_bound <- sum(mat > 0L, na.rm = TRUE)
    .log_inform(
      "FP bound by condition: {n_peaks} peaks × {n_cond} conditions; {total_bound} total bound calls."
    )
  }
  grn_status_set(grn_set, "fp_bound_condition")
}

#' @rdname grn_set_helpers
#' @export
grn_filter_fp_bound_condition <- function(
    grn_set,
    min_bound = 1L,
    use_parallel = TRUE,
    force = FALSE,
    verbose = FALSE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_filtered") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_bound_condition)) {
    .log_abort("`grn_set$fp_bound_condition` is missing; run grn_add_fp_bound_condition() first.")
  }
  n_before <- nrow(grn_set$fp_bound_condition)
  keep_peaks <- get_min_bound_peaks(
    grn_set$fp_bound_condition,
    min_bound = min_bound,
    use_parallel = use_parallel
  )
  grn_set <- filter_fp_rows(
    grn_set = grn_set,
    peaks = keep_peaks,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose = FALSE
  )
  if (is.data.frame(grn_set$fp_score_condition)) {
    ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition$peak_ID)
    if (anyNA(ord)) .log_abort("fp_score_condition missing peaks after filtering.")
    grn_set$fp_score_condition <- grn_set$fp_score_condition[ord, , drop = FALSE]
  }
  if (is.data.frame(grn_set$fp_bound_condition)) {
    ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_bound_condition$peak_ID)
    if (anyNA(ord)) .log_abort("fp_bound_condition missing peaks after filtering.")
    grn_set$fp_bound_condition <- grn_set$fp_bound_condition[ord, , drop = FALSE]
  }
  if (isTRUE(verbose)) {
    n_after <- nrow(grn_set$fp_bound_condition)
    .log_inform("FP filter: kept {n_after}/{n_before} peaks (min_bound = {min_bound}).")
  }
  grn_status_set(grn_set, "fp_filtered")
}

#' @rdname grn_set_helpers
#' @export
grn_add_fp_score_qn <- function(
    grn_set,
    id_col = "peak_ID",
    force = FALSE,
    verbose = TRUE
) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_score_qn") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_score_condition)) {
    .log_abort("`grn_set$fp_score_condition` is missing; run grn_add_fp_score_condition() first.")
  }
  grn_set$fp_score_condition_qn <- qn_footprints(grn_set$fp_score_condition, id_col = id_col)
  if (isTRUE(verbose)) {
    n_peaks <- nrow(grn_set$fp_score_condition_qn)
    n_cond <- max(ncol(grn_set$fp_score_condition_qn) - 1L, 0L)
    .log_inform("FP score QN: {n_peaks} peaks × {n_cond} conditions.")
  }
  grn_status_set(grn_set, "fp_score_qn")
}

#' @rdname grn_set_helpers
#' @export
grn_add_fp_tfs <- function(grn_set, force = FALSE, verbose = TRUE) {
  grn_set <- grn_status_init(grn_set)
  if (grn_status_is(grn_set, "fp_tfs") && !isTRUE(force)) return(grn_set)
  if (!is.data.frame(grn_set$fp_annotation)) {
    .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  }
  annot <- grn_set$fp_annotation
  if (!"fp_peak" %in% names(annot)) .log_abort("`fp_annotation` needs fp_peak.")
  if (!"motifs" %in% names(annot)) .log_abort("`fp_annotation` needs motifs.")

  gene_col <- .motif_gene_col(grn_set$motif_db)
  has_db <- is.data.frame(grn_set$motif_db) &&
    !is.null(gene_col) &&
    all(c("motif", gene_col) %in% names(grn_set$motif_db))

  if (has_db) {
    annot <- annot |>
      dplyr::left_join(
        grn_set$motif_db |>
          dplyr::select(motif, gene_symbol = dplyr::all_of(gene_col)),
        by = c("motifs" = "motif")
      )
    annot$tfs <- ifelse(
      !is.na(annot$gene_symbol) & nzchar(annot$gene_symbol),
      gsub("\\s*::\\s*", ",", annot$gene_symbol),
      toupper(sub("_.*$", "", annot$motifs))
    )
    annot$gene_symbol <- NULL
  } else {
    annot$tfs <- toupper(sub("_.*$", "", annot$motifs))
  }

  grn_set$fp_annotation <- annot

  fp_tfs <- annot |>
    dplyr::select(fp_peak, atac_peak, tfs) |>
    tidyr::separate_rows(tfs, sep = "\\s*,\\s*|\\s*::\\s*") |>
    dplyr::filter(!is.na(tfs), tfs != "") |>
    dplyr::distinct(fp_peak, atac_peak, tfs) |>
    dplyr::group_by(fp_peak, atac_peak) |>
    dplyr::summarise(tfs = paste(unique(tfs), collapse = ","), .groups = "drop")

  grn_set$fp_tfs <- fp_tfs

  if (isTRUE(verbose)) {
    n_rows <- nrow(annot)
    n_peaks <- length(unique(annot$fp_peak))
    n_tfs <- if (nrow(fp_tfs)) {
      length(unique(unlist(strsplit(fp_tfs$tfs, "\\s*,\\s*"))))
    } else {
      0L
    }
    .log_inform("FP TF mapping: {n_rows} annotations; {n_peaks} peaks; {n_tfs} TFs.")
  }

  grn_status_set(grn_set, "fp_tfs")
}

.expand_fp_annotation_tfs <- function(fp_annotation, tfs_col = "tfs") {
  if (!is.data.frame(fp_annotation)) .log_abort("`fp_annotation` must be a data.frame.")
  if (!tfs_col %in% names(fp_annotation)) {
    .log_abort("`fp_annotation` is missing column {.val {tfs_col}}.")
  }
  tfs_sym <- rlang::sym(tfs_col)
  fp_annotation |>
    tidyr::separate_rows(!!tfs_sym, sep = "\\s*,\\s*|\\s*::\\s*") |>
    dplyr::filter(!is.na(!!tfs_sym), !!tfs_sym != "")
}

.fp_tf_corr_add_p_adj <- function(fp_annotation, p_col = "corr_fp_tf_p") {
  if (!is.data.frame(fp_annotation)) .log_abort("`fp_annotation` must be a data.frame.")
  need <- c("fp_peak", "tfs", p_col)
  if (!all(need %in% names(fp_annotation))) {
    .log_abort("`fp_annotation` must include columns: {paste(need, collapse = ', ')}.")
  }

  annot <- fp_annotation |>
    dplyr::distinct(fp_peak, tfs, .keep_all = TRUE)

  p_vals <- annot[[p_col]]
  annot$corr_fp_tf_p_adj <- NA_real_
  good_p <- is.finite(p_vals)
  if (any(good_p)) {
    adj <- stats::p.adjust(p_vals[good_p], method = "BH")
    annot$corr_fp_tf_p_adj[good_p] <- pmax(adj, p_vals[good_p])
  }
  annot
}

#' @rdname grn_set_helpers
#' @export
grn_add_fp_tf_corr <- function(
    grn_set,
    method = c("pearson", "spearman", "kendall"),
    mode = c("canonical", "all"),
    tf_subset = NULL,
    cores = 10L,
    chunk_size = 5000L,
    min_non_na = 5L,
    force = FALSE,
    verbose = TRUE
) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  grn_set <- grn_status_init(grn_set)
  status_key <- paste0("fp_tf_corr_", method, "_", mode)
  if (grn_status_is(grn_set, status_key) && !isTRUE(force)) return(grn_set)

  if (!is.data.frame(grn_set$fp_score_condition_qn)) {
    .log_abort("`grn_set$fp_score_condition_qn` is missing; run grn_add_fp_score_qn() first.")
  }
  if (!is.data.frame(grn_set$rna_condition)) {
    .log_abort("`grn_set$rna_condition` is missing; run grn_add_rna_condition() first.")
  }
  if (!is.data.frame(grn_set$fp_annotation)) {
    .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  }

  if (!"tfs" %in% names(grn_set$fp_annotation)) {
    grn_set <- grn_add_fp_tfs(grn_set, force = TRUE, verbose = verbose)
  }

  tf_subset <- if (is.null(tf_subset)) NULL else unique(as.character(tf_subset))
  if (identical(mode, "canonical")) {
    ann_exp <- .expand_fp_annotation_tfs(grn_set$fp_annotation, tfs_col = "tfs")
    if (!is.null(tf_subset) && length(tf_subset)) {
      ann_exp <- ann_exp[ann_exp$tfs %in% tf_subset, , drop = FALSE]
    }
  } else {
    rna_cond <- grn_set$rna_condition
    tf_all <- unique(rna_cond$HGNC)
    tf_all <- tf_all[!is.na(tf_all) & nzchar(tf_all)]
    if (!is.null(tf_subset) && length(tf_subset)) {
      tf_all <- intersect(tf_all, tf_subset)
    }
    if (!length(tf_all)) {
      ann_name <- paste0("fp_annotation_", method)
      ann_exp <- grn_set$fp_annotation[0, , drop = FALSE]
      if (!"tfs" %in% names(ann_exp)) ann_exp$tfs <- character(0)
      grn_set[[ann_name]] <- ann_exp
      grn_status_set(grn_set, status_key)
      return(grn_set)
    }
    base <- grn_set$fp_annotation[, c("fp_peak", "atac_peak", "motifs"), drop = FALSE]
    base <- base[!duplicated(base[, c("fp_peak", "motifs")]), , drop = FALSE]
    base <- base[!is.na(base$motifs) & nzchar(base$motifs), , drop = FALSE]
    n_peaks <- nrow(base)
    ann_exp <- data.frame(
      fp_peak = rep(base$fp_peak, each = length(tf_all)),
      atac_peak = rep(base$atac_peak, each = length(tf_all)),
      tfs = rep(tf_all, times = n_peaks),
      motifs = rep(base$motifs, each = length(tf_all)),
      stringsAsFactors = FALSE
    )
  }
  if (!nrow(ann_exp)) {
    ann_name <- paste0("fp_annotation_", method)
    grn_set[[ann_name]] <- ann_exp
    grn_status_set(grn_set, status_key)
    return(grn_set)
  }

  fp_score_qn <- grn_set$fp_score_condition_qn
  if (anyDuplicated(fp_score_qn$peak_ID)) {
    fp_score_qn <- fp_score_qn[!duplicated(fp_score_qn$peak_ID), , drop = FALSE]
  }
  rna_cond <- grn_set$rna_condition
  if (anyDuplicated(rna_cond$HGNC)) {
    rna_cond <- rna_cond[!duplicated(rna_cond$HGNC), , drop = FALSE]
  }

  tmp <- list(
    fp_annotation = ann_exp,
    fp_score = fp_score_qn,
    rna = rna_cond
  )

  ann_corr <- make_fp_annotation_corr(
    tmp,
    method = method,
    cores = cores,
    chunk_size = chunk_size,
    min_non_na = min_non_na
  )
  ann_corr <- .fp_tf_corr_add_p_adj(ann_corr, p_col = "corr_fp_tf_p")

  ann_name <- paste0("fp_annotation_", method)
  grn_set[[ann_name]] <- ann_corr

  if (isTRUE(verbose)) {
    n_rows <- nrow(ann_corr)
    n_peaks <- length(unique(ann_corr$fp_peak))
    n_tfs <- length(unique(ann_corr$tfs))
    .log_inform(
      "FP-TF correlation ({method}): {n_rows} rows; {n_peaks} peaks; {n_tfs} TFs."
    )
  }

  grn_status_set(grn_set, status_key)
}

#' @rdname grn_set_helpers
#' @export
grn_filter_fp_tf_corr <- function(
    grn_set,
    method = c("pearson", "spearman", "kendall"),
    mode = c("canonical", "all"),
    r_thr = 0.3,
    p_thr = 0.05,
    set_active = TRUE,
    output_bed = NULL,
    output_bed_condition = NULL,
    label_col = NULL,
    verbose = TRUE
) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  ann_name <- paste0("fp_annotation_", method)
  if (!is.data.frame(grn_set[[ann_name]])) {
    .log_abort("`grn_set` is missing {ann_name}; run grn_add_fp_tf_corr() first.")
  }

  ann_all <- grn_set[[ann_name]]
  if (!"corr_fp_tf_p_adj" %in% names(ann_all)) {
    ann_all <- .fp_tf_corr_add_p_adj(ann_all, p_col = "corr_fp_tf_p")
  }
  grn_set[[ann_name]] <- ann_all

  keep <- is.finite(ann_all$corr_fp_tf_r) & is.finite(ann_all$corr_fp_tf_p_adj)
  keep <- keep & (ann_all$corr_fp_tf_r > r_thr) & (ann_all$corr_fp_tf_p_adj < p_thr)
  ann_filt <- ann_all[keep, , drop = FALSE]

  grn_set[[paste0(ann_name, "_filtered")]] <- ann_filt
  grn_status_set(grn_set, paste0("fp_tf_corr_filtered_", method))

  if (isTRUE(set_active)) {
    grn_set$fp_annotation <- ann_filt
    if (!nrow(ann_filt)) {
      grn_set$fp_score <- grn_set$fp_score[0, , drop = FALSE]
      grn_set$fp_bound <- grn_set$fp_bound[0, , drop = FALSE]
      if (is.data.frame(grn_set$fp_score_condition)) {
        grn_set$fp_score_condition <- grn_set$fp_score_condition[0, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_bound_condition)) {
        grn_set$fp_bound_condition <- grn_set$fp_bound_condition[0, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_score_condition_qn)) {
        grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[0, , drop = FALSE]
      }
    } else {
      peaks_keep <- unique(ann_filt$fp_peak)
      grn_set <- filter_fp_rows(
        grn_set = grn_set,
        peaks = peaks_keep,
        score_key = "peak_ID",
        bound_key = "peak_ID",
        annot_key = "fp_peak",
        verbose = FALSE,
        warn_on_mismatch = FALSE
      )
      if (is.data.frame(grn_set$fp_score_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition$peak_ID)
        if (anyNA(ord)) .log_abort("fp_score_condition missing peaks after TF-corr filter.")
        grn_set$fp_score_condition <- grn_set$fp_score_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_bound_condition)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_bound_condition$peak_ID)
        if (anyNA(ord)) .log_abort("fp_bound_condition missing peaks after TF-corr filter.")
        grn_set$fp_bound_condition <- grn_set$fp_bound_condition[ord, , drop = FALSE]
      }
      if (is.data.frame(grn_set$fp_score_condition_qn)) {
        ord <- match(grn_set$fp_score$peak_ID, grn_set$fp_score_condition_qn$peak_ID)
        if (anyNA(ord)) .log_abort("fp_score_condition_qn missing peaks after TF-corr filter.")
        grn_set$fp_score_condition_qn <- grn_set$fp_score_condition_qn[ord, , drop = FALSE]
      }
    }
  }

  if (isTRUE(verbose)) {
    n_all <- nrow(ann_all)
    n_keep <- nrow(ann_filt)
    .log_inform(
      "FP-TF filter ({method}): kept {n_keep}/{n_all} rows (r>{r_thr}, p_adj<{p_thr})."
    )
  }

  if (!is.null(output_bed) || !is.null(output_bed_condition)) {
    .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
    .mk_table <- function(dt) {
      if (!nrow(dt)) return(data.table::data.table())
      fp_col <- if ("fp_peak" %in% names(dt)) "fp_peak" else "peak_ID"
      if (!fp_col %in% names(dt)) .log_abort("`fp_annotation` is missing fp_peak/peak_ID.")
      if (!"atac_peak" %in% names(dt)) .log_abort("`fp_annotation` is missing atac_peak.")
      if (!"tfs" %in% names(dt)) .log_abort("`fp_annotation` is missing tfs.")
      dt <- data.table::as.data.table(dt)
      dt2 <- dt[!is.na(get(fp_col)) & !is.na(atac_peak) &
                  !is.na(tfs) & nzchar(tfs)]
      if (!nrow(dt2)) return(data.table::data.table())
      fp_s <- data.table::tstrsplit(dt2[[fp_col]], "[:-]", perl = TRUE)
      at_s <- data.table::tstrsplit(dt2[["atac_peak"]], "[:-]", perl = TRUE)
      data.table::data.table(
        fp_key          = dt2[[fp_col]],
        TFBS_chr         = fp_s[[1]],
        TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
        TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
        TFBS_name        = dt2$motifs,
        peak_chr         = at_s[[1]],
        peak_start       = suppressWarnings(as.integer(at_s[[2]])),
        peak_end         = suppressWarnings(as.integer(at_s[[3]])),
        atac_peak        = dt2$atac_peak,
        TF               = dt2$tfs,
        corr_fp_tf_r     = dt2$corr_fp_tf_r,
        corr_fp_tf_p     = dt2$corr_fp_tf_p,
        corr_fp_tf_p_adj = dt2$corr_fp_tf_p_adj
      )
    }

    all_dt <- .mk_table(ann_all)
    bound_dt <- .mk_table(ann_filt)
  }

  if (!is.null(output_bed)) {
    if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)

    if (nrow(all_dt)) all_dt[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]
    if (nrow(bound_dt)) bound_dt[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]

    cond_cols <- character(0)
    fp_bound_mat <- NULL
    tf_expr_mat <- NULL
    if (!is.null(label_col) &&
        is.data.frame(grn_set$fp_bound_condition) &&
        is.data.frame(grn_set$rna_expressed)) {
      cond_cols <- intersect(
        setdiff(names(grn_set$fp_bound_condition), "peak_ID"),
        setdiff(names(grn_set$rna_expressed), c("ensembl_gene_id", "HGNC"))
      )
      if (length(cond_cols)) {
        fp_bound_tbl <- grn_set$fp_bound_condition[, c("peak_ID", cond_cols), drop = FALSE]
        fp_bound_tbl$peak_ID <- as.character(fp_bound_tbl$peak_ID)
        fp_bound_mat <- as.matrix(fp_bound_tbl[, cond_cols, drop = FALSE])
        storage.mode(fp_bound_mat) <- "integer"
        rownames(fp_bound_mat) <- fp_bound_tbl$peak_ID

        tf_expr_tbl <- grn_set$rna_expressed[, c("HGNC", cond_cols), drop = FALSE]
        tf_expr_tbl$HGNC <- toupper(tf_expr_tbl$HGNC)
        tf_expr_tbl <- tf_expr_tbl[!is.na(tf_expr_tbl$HGNC) & tf_expr_tbl$HGNC != "", , drop = FALSE]
        tf_expr_tbl <- tf_expr_tbl[!duplicated(tf_expr_tbl$HGNC), , drop = FALSE]
        tf_expr_mat <- as.matrix(tf_expr_tbl[, cond_cols, drop = FALSE])
        storage.mode(tf_expr_mat) <- "integer"
        rownames(tf_expr_mat) <- tf_expr_tbl$HGNC
      }
    }

    tfs_vec <- sort(unique(all_dt$TF))
    for (tf in tfs_vec) {
      tf_slug <- .slug(tf)
      f_all <- file.path(output_bed, paste0(tf_slug, "_all.bed"))
      f_bound <- file.path(output_bed, paste0(tf_slug, "_bound.bed"))
      f_unbound <- file.path(output_bed, paste0(tf_slug, "_unbound.bed"))
      f_overview <- file.path(output_bed, paste0(tf_slug, "_overview.txt"))

      sub_all <- all_dt[TF == tf]
      sub_bound <- bound_dt[TF == tf]

      ids_bound <- if (nrow(sub_bound)) sub_bound$.id else character(0)
      sub_unbound <- sub_all[!(.id %in% ids_bound)]

      data.table::fwrite(sub_all[, !c(".id","atac_peak","fp_key"), with = FALSE], f_all, sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_bound[, !c(".id","atac_peak","fp_key"), with = FALSE], f_bound, sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_unbound[, !c(".id","atac_peak","fp_key"), with = FALSE], f_unbound, sep = "\t", col.names = FALSE, quote = FALSE)

      if (nrow(sub_all)) {
        sub_all[, `_bound` := as.integer(.id %in% ids_bound)]
        if (length(cond_cols) && !is.null(fp_bound_mat) && !is.null(tf_expr_mat)) {
          tf_key <- toupper(tf)
          tf_idx <- match(tf_key, rownames(tf_expr_mat))
          if (is.na(tf_idx)) {
            tf_ok <- rep(0L, length(cond_cols))
          } else {
            tf_ok <- as.integer(tf_expr_mat[tf_idx, ] > 0L)
          }
          fp_idx <- match(sub_all$fp_key, rownames(fp_bound_mat))
          fp_ok <- matrix(0L, nrow = nrow(sub_all), ncol = length(cond_cols))
          ok_idx <- !is.na(fp_idx)
          if (any(ok_idx)) {
            fp_ok[ok_idx, ] <- fp_bound_mat[fp_idx[ok_idx], , drop = FALSE]
          }
          bound_vec <- as.integer(sub_all$`_bound` > 0L)
          fp_ok <- fp_ok * bound_vec
          fp_ok <- sweep(fp_ok, 2, tf_ok, `*`)
          colnames(fp_ok) <- paste0(cond_cols, "_bound")
          sub_all[, (colnames(fp_ok)) := data.table::as.data.table(fp_ok)]
        }
        data.table::fwrite(sub_all[, !c(".id","atac_peak","fp_key"), with = FALSE], f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      } else {
        hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
                 "peak_chr","peak_start","peak_end","TF",
                 "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
        if (length(cond_cols)) {
          hdr <- c(hdr, paste0(cond_cols, "_bound"))
        }
        empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode = "character", length = 0)), length(hdr)), hdr))
        data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      }
    }
    if (isTRUE(verbose)) .log_inform("Per-TF BED exports written to: {output_bed}")
  }

  if (!is.null(output_bed_condition)) {
    if (!dir.exists(output_bed_condition)) dir.create(output_bed_condition, recursive = TRUE, showWarnings = FALSE)
    if (is.null(label_col)) {
      .log_abort("`label_col` is required when `output_bed_condition` is set.")
    }
    if (!is.data.frame(grn_set$atac_overlap)) {
      .log_abort("`grn_set$atac_overlap` is missing; cannot build condition-specific BEDs.")
    }
    meta_use <- grn_set$sample_metadata_used[, c("id", label_col), drop = FALSE]
    meta_use <- meta_use[!is.na(meta_use$id) & !is.na(meta_use[[label_col]]), , drop = FALSE]
    if (!nrow(meta_use)) .log_abort("No condition labels found in metadata for label_col {label_col}.")

    ov <- grn_set$atac_overlap
    if (!"atac_peak" %in% names(ov)) .log_abort("`grn_set$atac_overlap` must contain atac_peak.")
    ov_ids <- setdiff(names(ov), "atac_peak")
    cond_groups <- split(meta_use$id, meta_use[[label_col]], drop = TRUE)

    .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
    for (cond in names(cond_groups)) {
      ids <- intersect(cond_groups[[cond]], ov_ids)
      if (!length(ids)) next
      ov_mat <- as.matrix(ov[, ids, drop = FALSE])
      storage.mode(ov_mat) <- "integer"
      keep_atac <- ov$atac_peak[rowSums(ov_mat > 0L, na.rm = TRUE) > 0L]
      if (!length(keep_atac)) next

      cond_dir <- file.path(output_bed_condition, .slug(cond))
      if (!dir.exists(cond_dir)) dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)

      all_cond <- all_dt[atac_peak %in% keep_atac]
      bound_cond <- bound_dt[atac_peak %in% keep_atac]
      if (nrow(all_cond)) all_cond[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]
      if (nrow(bound_cond)) bound_cond[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]

      tfs_vec <- sort(unique(all_cond$TF))
      for (tf in tfs_vec) {
        tf_slug <- .slug(tf)
        f_all <- file.path(cond_dir, paste0(tf_slug, "_all.bed"))
        f_bound <- file.path(cond_dir, paste0(tf_slug, "_bound.bed"))
        f_unbound <- file.path(cond_dir, paste0(tf_slug, "_unbound.bed"))

        sub_all <- all_cond[TF == tf]
        sub_bound <- bound_cond[TF == tf]

        ids_bound <- if (nrow(sub_bound)) sub_bound$.id else character(0)
        sub_unbound <- sub_all[!(.id %in% ids_bound)]

        data.table::fwrite(sub_all[, !c(".id","atac_peak","fp_key"), with = FALSE], f_all, sep = "\t", col.names = FALSE, quote = FALSE)
        data.table::fwrite(sub_bound[, !c(".id","atac_peak","fp_key"), with = FALSE], f_bound, sep = "\t", col.names = FALSE, quote = FALSE)
        data.table::fwrite(sub_unbound[, !c(".id","atac_peak","fp_key"), with = FALSE], f_unbound, sep = "\t", col.names = FALSE, quote = FALSE)
      }
    }
    if (isTRUE(verbose)) .log_inform("Per-condition TF BED exports written to: {output_bed_condition}")
  }

  grn_set
}

#' Correlate TF expression to footprint scores (Pearson + Spearman)
#'
#' Runs Pearson and Spearman correlations, applies filters, and optionally
#' generates QC outputs and per-TF overview tables. Designed as a simple
#' user-facing wrapper for Module 1.
#'
#' @param omics_data A multi-omic data list from [load_multiomic_data()].
#' @param grn_set (Deprecated) Use `omics_data`.
#' @param config Optional YAML config path used to load inputs when
#'   \code{omics_data} is NULL.
#' @param genome Optional genome string (e.g. "hg38", "mm10").
#' @param gene_symbol_col Column in RNA table to treat as gene symbols.
#' @param fp_aligned Optional aligned footprint object returned by
#'   [align_footprints()]. When NULL, cached aligned footprints are loaded.
#' @param do_preprocess If TRUE and \code{fp_aligned} is NULL, runs footprint
#'   preprocessing (load + align) using config paths.
#' @param do_motif_clustering If TRUE, runs motif clustering before correlation.
#' @param fp_root_dir Optional root directory containing footprint outputs.
#' @param fp_cache_dir Optional cache directory for aligned footprint files.
#' @param fp_cache_tag Optional cache tag used for aligned footprint files.
#' @param mode Correlation mode: "canonical" or "all".
#' @param label_col Metadata label column used for conditions.
#' @param out_dir Output directory for QC and overviews. If NULL, defaults to
#'   \code{file.path(base_dir, "predict_tf_binding_sites")}.
#' @param db Database tag used in output filenames.
#' @param r_thr Correlation threshold.
#' @param p_thr Adjusted P-value threshold.
#' @param tf_subset Optional TF subset.
#' @param cores_pearson Cores for Pearson correlation.
#' @param cores_spearman Cores for Spearman correlation.
#' @param chunk_size Chunk size for correlation.
#' @param min_non_na Minimum non-NA values.
#' @param qc If TRUE, generate QC PDF and overview tables.
#' @param write_bed If TRUE, emit per-TF BEDs (off by default).
#' @param use_cache If TRUE, reuse cached correlation tables.
#' @param cache_dir Optional cache directory (defaults under out_dir).
#' @param verbose Emit status messages.
#'
#' @return Updated grn_set with Pearson + Spearman annotations.
#' @export
correlate_tf_to_fp <- function(
    omics_data = NULL,
    grn_set = NULL,
    config = NULL,
    genome = NULL,
    gene_symbol_col = "HGNC",
    fp_aligned = NULL,
    do_preprocess = FALSE,
    do_motif_clustering = FALSE,
    fp_root_dir = NULL,
    fp_cache_dir = NULL,
    fp_cache_tag = NULL,
    mode = c("canonical", "all"),
    label_col,
    out_dir = NULL,
    db,
    r_thr = 0.3,
    p_thr = 0.05,
    tf_subset = NULL,
    cores_pearson = 20L,
    cores_spearman = 36L,
    chunk_size = 5000L,
    min_non_na = 5L,
    qc = TRUE,
    write_bed = FALSE,
    write_outputs = TRUE,
    use_cache = TRUE,
    cache_dir = NULL,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  if (is.null(omics_data)) omics_data <- grn_set
  # Default config lookup helper (reads from globals set by load_config)
  cfg_val <- function(name) {
    if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      return(get(name, envir = .GlobalEnv))
    }
    NULL
  }
  if (is.null(omics_data) &&
      is.null(config) &&
      is.null(fp_aligned) &&
      !isTRUE(do_preprocess)) {
    .log_abort("`omics_data` not found. Provide it or set `config`/`do_preprocess`.")
  }
  if (is.null(omics_data) && (!is.null(config) || !is.null(fp_aligned) || isTRUE(do_preprocess))) {
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
    cfg_val <- function(name) {
      if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
        return(get(name, envir = .GlobalEnv))
      }
      NULL
    }
    if (is.null(fp_cache_dir)) {
      base_dir_cfg <- cfg_val("base_dir")
      if (is.character(base_dir_cfg) && nzchar(base_dir_cfg)) {
        fp_cache_dir <- file.path(base_dir_cfg, "cache")
      }
    }
    if (is.null(fp_cache_tag)) {
      fp_cache_tag <- cfg_val("db")
    }

    if (is.null(fp_aligned)) {
      if (isTRUE(do_preprocess)) {
        if (is.null(fp_root_dir)) fp_root_dir <- cfg_val("fp_root_dir")
        if (is.null(fp_root_dir) || !nzchar(fp_root_dir)) {
          .log_abort("`fp_root_dir` not set. Provide it or set in config.")
        }
        if (is.null(db) || !nzchar(db)) {
          db <- cfg_val("db")
        }
        if (is.null(db) || !nzchar(db)) .log_abort("`db` must be provided or set in config.")
        fp_manifest <- load_footprints(
          root_dir = fp_root_dir,
          db_name = db,
          out_dir = file.path(fp_cache_dir, paste0("fp_", db))
        )
        fp_aligned <- align_footprints(
          fp_manifest,
          mid_slop = 10L,
          round_digits = 1L,
          score_match_pct = 0.8,
          cache_dir = fp_cache_dir,
          cache_tag = db,
          output_mode = "distinct"
        )
      } else {
        fp_aligned <- load_fp_aligned_from_cache(
          cache_dir = fp_cache_dir,
          cache_tag = fp_cache_tag,
          output_mode = "distinct",
          verbose = verbose
        )
      }
    }

    if (isTRUE(do_motif_clustering)) {
      run_fp_motif_clustering_pre_if_needed(
        fp_aligned = fp_aligned,
        base_dir = cfg_val("base_dir"),
        ref_db = db,
        motif_db = cfg_val("motif_db")
      )
      run_fp_motif_clustering(
        fp_aligned = fp_aligned,
        base_dir = cfg_val("base_dir"),
        ref_db = db,
        motif_db = cfg_val("motif_db"),
        mode = "data",
        target_clusters = 220,
        qc_mode = "fast",
        save_motif_db = TRUE
      )
      run_fp_motif_clustering(
        fp_aligned = fp_aligned,
        base_dir = cfg_val("base_dir"),
        ref_db = db,
        motif_db = cfg_val("motif_db"),
        mode = "hybrid",
        target_clusters = 165,
        qc_mode = "fast",
        save_motif_db = TRUE
      )
    }

    omics_data <- load_multiomic_data(
      config = config,
      genome = genome,
      gene_symbol_col = gene_symbol_col,
      fp_aligned = fp_aligned,
      label_col = label_col,
      expected_n = cfg_val("expected_n"),
      tf_list = cfg_val("tf_list"),
      motif_db = cfg_val("motif_db"),
      threshold_gene_expr = cfg_val("threshold_gene_expr"),
      threshold_fp_score = cfg_val("threshold_fp_score"),
      use_parallel = TRUE,
      verbose = verbose
    )
  }
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.character(label_col) || !nzchar(label_col)) {
    .log_abort("`label_col` must be a non-empty string.")
  }
  if (is.null(out_dir) || !is.character(out_dir) || !nzchar(out_dir)) {
    base_dir_cfg <- cfg_val("base_dir")
    if (!is.character(base_dir_cfg) || !nzchar(base_dir_cfg)) {
      .log_abort("`out_dir` not set and `base_dir` not available. Provide `out_dir` or load config.")
    }
    out_dir <- file.path(base_dir_cfg, "predict_tf_binding_sites")
  }
  if (!is.character(db) || !nzchar(db)) {
    db <- cfg_val("db")
  }
  if (!is.character(db) || !nzchar(db)) {
    .log_abort("`db` must be a non-empty string.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cache_dir)) {
    cache_dir <- file.path(out_dir, "cache", "fp_tf_corr")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  grn_set <- omics_data
  if (isTRUE(verbose)) {
    .log_inform("Computing TF expression vs footprint score correlations (Pearson + Spearman).")
  }
  grn_set <- grn_add_rna_condition(grn_set, label_col = label_col, verbose = verbose)
  grn_set <- grn_add_fp_tfs(grn_set, verbose = verbose)
  grn_set <- grn_add_fp_score_qn(grn_set, id_col = "peak_ID", verbose = verbose)

  .load_cache <- function(method) {
    cache_path <- file.path(cache_dir, sprintf("fp_tf_corr_%s_%s.rds", method, mode))
    if (isTRUE(use_cache) && file.exists(cache_path)) {
      if (isTRUE(verbose)) .log_inform("Using cached {method} TF correlations: {cache_path}")
      readRDS(cache_path)
    } else {
      NULL
    }
  }
  .save_cache <- function(method, tbl) {
    cache_path <- file.path(cache_dir, sprintf("fp_tf_corr_%s_%s.rds", method, mode))
    if (isTRUE(use_cache) && is.data.frame(tbl)) {
      saveRDS(tbl, cache_path)
    }
  }

  grn_set_base <- grn_set

  # Pearson (set_active = TRUE)
  ann_p <- .load_cache("pearson")
  if (is.data.frame(ann_p)) {
    grn_set$fp_annotation_pearson <- ann_p
    grn_status_set(grn_set, paste0("fp_tf_corr_pearson_", mode))
  } else {
    grn_set <- grn_add_fp_tf_corr(
      grn_set,
      method = "pearson",
      mode = mode,
      tf_subset = tf_subset,
      cores = cores_pearson,
      chunk_size = chunk_size,
      min_non_na = min_non_na,
      verbose = verbose
    )
    .save_cache("pearson", grn_set$fp_annotation_pearson)
  }
  out_bed <- if (isTRUE(write_bed)) file.path(out_dir, sprintf("fp_predicted_tfbs_%s", db)) else NULL
  out_bed_cond <- if (isTRUE(write_bed)) file.path(out_dir, sprintf("fp_predicted_tfbs_%s_by_condition", db)) else NULL
  grn_set <- grn_filter_fp_tf_corr(
    grn_set,
    method = "pearson",
    mode = mode,
    r_thr = r_thr,
    p_thr = p_thr,
    set_active = TRUE,
    output_bed = out_bed,
    output_bed_condition = out_bed_cond,
    label_col = label_col,
    verbose = verbose
  )
  if (isTRUE(write_outputs)) {
    write_grn_tf_corr_outputs(grn_set, out_dir = out_dir, db = db)
  }

  # Spearman on full (unfiltered) set; preserve pearson-filtered fp_annotation in grn_set
  grn_set_spearman <- grn_set_base
  ann_s <- .load_cache("spearman")
  if (is.data.frame(ann_s)) {
    grn_set_spearman$fp_annotation_spearman <- ann_s
    grn_status_set(grn_set_spearman, paste0("fp_tf_corr_spearman_", mode))
  } else {
    grn_set_spearman <- grn_add_fp_tf_corr(
      grn_set_spearman,
      method = "spearman",
      mode = mode,
      tf_subset = tf_subset,
      cores = cores_spearman,
      chunk_size = chunk_size,
      min_non_na = min_non_na,
      verbose = verbose
    )
    .save_cache("spearman", grn_set_spearman$fp_annotation_spearman)
  }
  grn_set_spearman <- grn_filter_fp_tf_corr(
    grn_set_spearman,
    method = "spearman",
    mode = mode,
    r_thr = r_thr,
    p_thr = p_thr,
    set_active = FALSE,
    verbose = verbose
  )

  grn_set$fp_annotation_spearman <- grn_set_spearman$fp_annotation_spearman
  grn_set$fp_annotation_spearman_filtered <- grn_set_spearman$fp_annotation_spearman_filtered

  if (isTRUE(qc)) {
    plot_tf_corr_stats_pdf(
      ann_pearson = grn_set$fp_annotation_pearson_filtered,
      ann_spearman = grn_set$fp_annotation_spearman_filtered,
      out_dir = out_dir,
      db = db,
      mode = mode,
      r_thr = r_thr,
      p_thr = p_thr,
      ann_pearson_all = grn_set$fp_annotation_pearson,
      ann_spearman_all = grn_set$fp_annotation_spearman,
      verbose = verbose
    )
    overview_dir <- file.path(out_dir, sprintf("06_fp_predicted_tfbs_%s_%s", db, mode))
    write_tf_tfbs_overviews(
      omics_data = grn_set,
      ann_pearson = grn_set$fp_annotation_pearson,
      ann_spearman = grn_set$fp_annotation_spearman,
      out_dir = overview_dir,
      db = db,
      label_col = label_col,
      r_thr = r_thr,
      p_thr = p_thr,
      verbose = verbose
    )
  }

  grn_set
}

#' @rdname grn_set_helpers
#' @export
write_grn_tf_corr_outputs <- function(grn_set, out_dir, db) {
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(out_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.data.frame(grn_set$fp_score)) .log_abort("`grn_set$fp_score` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_bound)) .log_abort("`grn_set$fp_bound` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_annotation)) .log_abort("`grn_set$fp_annotation` is missing or invalid.")

  readr::write_csv(
    grn_set$fp_score,
    file.path(cache_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db))
  )
  readr::write_csv(
    grn_set$fp_bound,
    file.path(cache_dir, sprintf("fp_bound_strict_tf_filtered_corr_%s.csv", db))
  )
  readr::write_csv(
    grn_set$fp_annotation,
    file.path(cache_dir, sprintf("fp_annotation_strict_tf_filtered_corr_%s.csv", db))
  )

  invisible(out_dir)
}

#' Write per-TF overview tables with Pearson + Spearman correlations
#'
#' Generates per-TF `<TF>_overview.txt` files and a summary table of
#' predicted binding site counts per TF per condition. By default, bound
#' calls are flagged when either Pearson or Spearman passes the thresholds.
#'
#' @param omics_data A multi-omic data list containing fp_bound_condition and rna_expressed.
#' @param grn_set (Deprecated) Use `omics_data`.
#' @param ann_pearson Annotation table with Pearson correlations (full).
#' @param ann_spearman Annotation table with Spearman correlations (full).
#' @param out_dir Directory to write per-TF overview files.
#' @param db Database tag used in the summary filename.
#' @param label_col Label column (used for condition names if needed).
#' @param r_thr Correlation threshold.
#' @param p_thr Adjusted P-value threshold.
#' @param verbose Emit status messages.
#'
#' @return Invisibly returns the output directory.
#' @export
write_tf_tfbs_overviews <- function(
    omics_data = NULL,
    grn_set = NULL,
    ann_pearson,
    ann_spearman,
    out_dir,
    db,
    label_col = NULL,
    r_thr = 0.3,
    p_thr = 0.05,
    verbose = TRUE
) {
  if (is.null(omics_data)) omics_data <- grn_set
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  if (!is.data.frame(ann_pearson) && !is.data.frame(ann_spearman)) {
    .log_abort("Both Pearson and Spearman annotation tables are missing.")
  }
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  grn_set <- omics_data
  if (!is.data.frame(grn_set$fp_bound_condition)) {
    .log_abort("`omics_data$fp_bound_condition` is missing; run grn_add_fp_bound_condition() first.")
  }
  if (!is.data.frame(grn_set$rna_expressed)) {
    .log_abort("`omics_data$rna_expressed` is missing; run grn_add_rna_expressed() first.")
  }

  r_thr_use <- if (is.finite(r_thr)) r_thr else -Inf
  p_thr_use <- if (is.finite(p_thr)) p_thr else 1

  normalize_ann <- function(dt, prefix) {
    if (!is.data.frame(dt) || !nrow(dt)) {
      return(data.table::data.table())
    }
    if (!"corr_fp_tf_p_adj" %in% names(dt)) {
      dt <- .fp_tf_corr_add_p_adj(dt, p_col = "corr_fp_tf_p")
    }
    fp_col <- if ("fp_peak" %in% names(dt)) "fp_peak" else if ("peak_ID" %in% names(dt)) "peak_ID" else NA_character_
    if (is.na(fp_col)) .log_abort("Annotation table missing fp_peak/peak_ID.")
    tf_col <- if ("tfs" %in% names(dt)) "tfs" else if ("TF" %in% names(dt)) "TF" else NA_character_
    if (is.na(tf_col)) .log_abort("Annotation table missing tfs/TF.")
    if (!"atac_peak" %in% names(dt)) .log_abort("Annotation table missing atac_peak.")
    if (!"motifs" %in% names(dt)) .log_abort("Annotation table missing motifs.")

    dt <- data.table::as.data.table(dt)
    out <- dt[, .(
      fp_peak = as.character(get(fp_col)),
      atac_peak = as.character(atac_peak),
      motifs = as.character(motifs),
      tfs = as.character(get(tf_col)),
      corr_fp_tf_r = corr_fp_tf_r,
      corr_fp_tf_p = corr_fp_tf_p,
      corr_fp_tf_p_adj = corr_fp_tf_p_adj
    )]
    data.table::setnames(
      out,
      c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj"),
      paste0(prefix, c("corr_fp_tf_r", "corr_fp_tf_p", "corr_fp_tf_p_adj"))
    )
    out
  }

  ann_p <- normalize_ann(ann_pearson, "pearson_")
  ann_s <- normalize_ann(ann_spearman, "spearman_")
  if (!nrow(ann_p) && !nrow(ann_s)) {
    .log_abort("No rows in Pearson or Spearman annotations.")
  }

  all_dt <- if (nrow(ann_p) && nrow(ann_s)) {
    data.table::as.data.table(
      base::merge(
        ann_p,
        ann_s,
        by = c("fp_peak", "atac_peak", "motifs", "tfs"),
        all = TRUE
      )
    )
  } else if (nrow(ann_p)) {
    ann_p
  } else {
    ann_s
  }

  # overall bound flag (either method passes)
  pearson_ok <- with(all_dt, is.finite(pearson_corr_fp_tf_r) &
                       is.finite(pearson_corr_fp_tf_p_adj) &
                       pearson_corr_fp_tf_r > r_thr_use &
                       pearson_corr_fp_tf_p_adj < p_thr_use)
  spearman_ok <- with(all_dt, is.finite(spearman_corr_fp_tf_r) &
                        is.finite(spearman_corr_fp_tf_p_adj) &
                        spearman_corr_fp_tf_r > r_thr_use &
                        spearman_corr_fp_tf_p_adj < p_thr_use)
  all_dt[, `_bound` := as.integer(pearson_ok | spearman_ok)]

  fp_s <- data.table::tstrsplit(all_dt[["fp_peak"]], "[:-]", perl = TRUE)
  at_s <- data.table::tstrsplit(all_dt[["atac_peak"]], "[:-]", perl = TRUE)
  all_dt[, `:=`(
    TFBS_chr = fp_s[[1]],
    TFBS_start = suppressWarnings(as.integer(fp_s[[2]])),
    TFBS_end = suppressWarnings(as.integer(fp_s[[3]])),
    TFBS_name = motifs,
    peak_chr = at_s[[1]],
    peak_start = suppressWarnings(as.integer(at_s[[2]])),
    peak_end = suppressWarnings(as.integer(at_s[[3]])),
    TF = tfs
  )]

  cond_cols <- intersect(
    setdiff(names(grn_set$fp_bound_condition), "peak_ID"),
    setdiff(names(grn_set$rna_expressed), c("ensembl_gene_id", "HGNC"))
  )
  if (length(cond_cols)) {
    fp_bound_tbl <- grn_set$fp_bound_condition[, c("peak_ID", cond_cols), drop = FALSE]
    fp_bound_tbl$peak_ID <- as.character(fp_bound_tbl$peak_ID)
    fp_bound_mat <- as.matrix(fp_bound_tbl[, cond_cols, drop = FALSE])
    storage.mode(fp_bound_mat) <- "integer"
    rownames(fp_bound_mat) <- fp_bound_tbl$peak_ID

    tf_expr_tbl <- grn_set$rna_expressed[, c("HGNC", cond_cols), drop = FALSE]
    tf_expr_tbl$HGNC <- toupper(tf_expr_tbl$HGNC)
    tf_expr_tbl <- tf_expr_tbl[!is.na(tf_expr_tbl$HGNC) & tf_expr_tbl$HGNC != "", , drop = FALSE]
    tf_expr_tbl <- tf_expr_tbl[!duplicated(tf_expr_tbl$HGNC), , drop = FALSE]
    tf_expr_mat <- as.matrix(tf_expr_tbl[, cond_cols, drop = FALSE])
    storage.mode(tf_expr_mat) <- "integer"
    rownames(tf_expr_mat) <- tf_expr_tbl$HGNC

    fp_idx <- match(all_dt$fp_peak, rownames(fp_bound_mat))
    fp_ok <- matrix(0L, nrow = nrow(all_dt), ncol = length(cond_cols))
    ok_idx <- !is.na(fp_idx)
    if (any(ok_idx)) {
      fp_ok[ok_idx, ] <- fp_bound_mat[fp_idx[ok_idx], , drop = FALSE]
    }
    tf_key <- toupper(all_dt$TF)
    tf_idx <- match(tf_key, rownames(tf_expr_mat))
    tf_ok <- matrix(0L, nrow = nrow(all_dt), ncol = length(cond_cols))
    ok_tf <- !is.na(tf_idx)
    if (any(ok_tf)) {
      tf_ok[ok_tf, ] <- tf_expr_mat[tf_idx[ok_tf], , drop = FALSE]
    }
    bound_vec <- as.integer(all_dt$`_bound` > 0L)
    fp_ok <- fp_ok * bound_vec
    fp_ok <- fp_ok * tf_ok
    colnames(fp_ok) <- paste0(cond_cols, "_bound")
    all_dt[, (colnames(fp_ok)) := data.table::as.data.table(fp_ok)]
  }

  keep_cols <- c(
    "TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name",
    "peak_chr", "peak_start", "peak_end", "TF",
    "pearson_corr_fp_tf_r", "pearson_corr_fp_tf_p", "pearson_corr_fp_tf_p_adj",
    "spearman_corr_fp_tf_r", "spearman_corr_fp_tf_p", "spearman_corr_fp_tf_p_adj",
    "_bound"
  )
  if (length(cond_cols)) {
    keep_cols <- c(keep_cols, paste0(cond_cols, "_bound"))
  }

  .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
  tfs_vec <- sort(unique(all_dt$TF))
  for (tf in tfs_vec) {
    tf_slug <- .slug(tf)
    f_overview <- file.path(out_dir, paste0(tf_slug, "_overview.txt"))
    sub_all <- all_dt[TF == tf]
    if (nrow(sub_all)) {
      data.table::fwrite(sub_all[, ..keep_cols], f_overview, sep = "\t", col.names = TRUE, quote = FALSE, na = "NA")
    } else {
      empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode = "character", length = 0)), length(keep_cols)), keep_cols))
      data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
    }
  }

  if (length(cond_cols)) {
    bound_cols <- paste0(cond_cols, "_bound")
    summary_tbl <- all_dt[, lapply(.SD, function(x) sum(as.integer(x) > 0L, na.rm = TRUE)),
                          by = TF, .SDcols = bound_cols]
    data.table::setnames(summary_tbl, bound_cols, cond_cols)
    summary_path <- file.path(dirname(out_dir), sprintf("06_tf_binding_site_counts_%s.csv", db))
    data.table::fwrite(summary_tbl, summary_path, sep = ",", col.names = TRUE)
  }

  if (isTRUE(verbose)) {
    .log_inform("Per-TF overviews written to: {out_dir}")
  }
  invisible(out_dir)
}

#' Export per-condition TFBS BED files from an omics object
#'
#' Convenience wrapper around [correlate_tf_to_fp()] that writes
#' per-condition TFBS BED files and returns the updated object.
#'
#' @param omics_data A multi-omic data list (output of [load_multiomic_data()]).
#' @param out_dir Output directory for TFBS BED files.
#' @param db Motif database tag used in output file names.
#' @param label_col Metadata column used for condition labels.
#' @param mode Correlation mode: \code{"canonical"} or \code{"all"}.
#' @param r_thr Correlation threshold.
#' @param p_thr Adjusted P-value threshold.
#' @param cores_pearson Number of cores for Pearson.
#' @param cores_spearman Number of cores for Spearman.
#' @param chunk_size Chunk size for correlation.
#' @param min_non_na Minimum non-NA observations per TF/peak.
#' @param use_cache Reuse cached correlations when available.
#' @param verbose Emit status messages.
#'
#' @return Updated multi-omic data list.
#' @export
output_per_condition_tfbs_bed_files <- function(
    omics_data,
    out_dir,
    db,
    label_col,
    mode = c("canonical", "all"),
    r_thr = 0.3,
    p_thr = 0.05,
    cores_pearson = 20L,
    cores_spearman = 36L,
    chunk_size = 5000L,
    min_non_na = 5L,
    use_cache = TRUE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  correlate_tf_to_fp(
    omics_data = omics_data,
    mode = mode,
    out_dir = out_dir,
    label_col = label_col,
    r_thr = r_thr,
    p_thr = p_thr,
    db = db,
    cores_pearson = cores_pearson,
    cores_spearman = cores_spearman,
    chunk_size = chunk_size,
    min_non_na = min_non_na,
    qc = FALSE,
    write_bed = TRUE,
    write_outputs = FALSE,
    use_cache = use_cache,
    verbose = verbose
  )
}

#' Experimental benchmarking of TFBS predictions
#'
#' Runs TFBS correlation in \code{"all"} mode and writes outputs for inspection.
#' Intended for exploratory benchmarking workflows.
#'
#' @param omics_data A multi-omic data list (output of [load_multiomic_data()]).
#' @param out_dir Output directory for benchmarking artifacts.
#' @param db Motif database tag used in output file names.
#' @param label_col Metadata column used for condition labels.
#' @param r_thr Correlation threshold.
#' @param p_thr Adjusted P-value threshold.
#' @param cores_pearson Number of cores for Pearson.
#' @param cores_spearman Number of cores for Spearman.
#' @param chunk_size Chunk size for correlation.
#' @param min_non_na Minimum non-NA observations per TF/peak.
#' @param use_cache Reuse cached correlations when available.
#' @param verbose Emit status messages.
#'
#' @return Updated multi-omic data list.
#' @export
experimental_benchmarking_of_tfbs_predictions <- function(
    omics_data,
    out_dir,
    db,
    label_col,
    r_thr = 0.3,
    p_thr = 0.05,
    cores_pearson = 20L,
    cores_spearman = 36L,
    chunk_size = 5000L,
    min_non_na = 5L,
    use_cache = TRUE,
    verbose = TRUE
) {
  correlate_tf_to_fp(
    omics_data = omics_data,
    mode = "all",
    out_dir = out_dir,
    label_col = label_col,
    r_thr = r_thr,
    p_thr = p_thr,
    db = db,
    cores_pearson = cores_pearson,
    cores_spearman = cores_spearman,
    chunk_size = chunk_size,
    min_non_na = min_non_na,
    qc = TRUE,
    write_bed = TRUE,
    write_outputs = TRUE,
    use_cache = use_cache,
    verbose = verbose
  )
}

#' @rdname grn_set_helpers
#' @export
write_grn_outputs <- function(grn_set, out_dir, db, qn_base_dir = NULL) {
  if (!is.list(grn_set)) .log_abort("`grn_set` must be a list.")
  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!is.character(db) || !nzchar(db)) .log_abort("`db` must be a non-empty string.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.data.frame(grn_set$fp_score)) .log_abort("`grn_set$fp_score` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_annotation)) .log_abort("`grn_set$fp_annotation` is missing or invalid.")
  if (!is.data.frame(grn_set$rna_expressed)) .log_abort("`grn_set$rna_expressed` is missing or invalid.")
  if (!is.data.frame(grn_set$fp_bound_condition)) .log_abort("`grn_set$fp_bound_condition` is missing or invalid.")

  fp_score_out <- grn_set$fp_score
  if (is.data.frame(fp_score_out) && "peak_ID" %in% names(fp_score_out)) {
    fp_score_out <- dplyr::distinct(fp_score_out, .data$peak_ID, .keep_all = TRUE)
  }
  readr::write_csv(fp_score_out, file.path(out_dir, sprintf("02_fp_score_raw_%s.csv", db)))
  readr::write_csv(grn_set$fp_annotation, file.path(out_dir, sprintf("03_fp_annotation_%s.csv", db)))
  readr::write_csv(grn_set$rna_expressed, file.path(out_dir, sprintf("05_gene_expr_flag_%s.csv", db)))
  fp_bound_out <- grn_set$fp_bound_condition
  if (is.data.frame(fp_bound_out) && "peak_ID" %in% names(fp_bound_out)) {
    fp_bound_out <- dplyr::distinct(fp_bound_out, .data$peak_ID, .keep_all = TRUE)
  }

  if (is.data.frame(grn_set$fp_score_condition_qn)) {
    fp_qn_out <- grn_set$fp_score_condition_qn
    if (is.data.frame(fp_qn_out) && "peak_ID" %in% names(fp_qn_out)) {
      fp_qn_out <- dplyr::distinct(fp_qn_out, .data$peak_ID, .keep_all = TRUE)
    }
    if (is.data.frame(fp_qn_out) && "peak_ID" %in% names(fp_qn_out) &&
        is.data.frame(fp_bound_out) && "peak_ID" %in% names(fp_bound_out)) {
      idx_bound <- match(fp_qn_out$peak_ID, fp_bound_out$peak_ID)
      if (anyNA(idx_bound)) {
        .log_warn("Some fp_bound_condition peaks are missing in fp_score_qn; dropping unmatched rows.")
        idx_keep <- which(!is.na(idx_bound))
        idx_bound <- idx_bound[idx_keep]
      }
      fp_bound_out <- fp_bound_out[idx_bound, , drop = FALSE]
    }
    readr::write_csv(
      fp_bound_out,
      file.path(out_dir, sprintf("03_fp_bound_condition_%s.csv", db))
    )
    if (!is.null(qn_base_dir)) {
      if (!is.character(qn_base_dir) || !nzchar(qn_base_dir)) {
        .log_abort("`qn_base_dir` must be a non-empty path when provided.")
      }
      dir.create(qn_base_dir, recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(
        fp_qn_out,
        file.path(qn_base_dir, sprintf("fp_scores_qn_%s.csv", db))
      )
    }
    readr::write_csv(
      fp_qn_out,
      file.path(out_dir, sprintf("03_fp_score_qn_%s.csv", db))
    )
  }

  invisible(out_dir)
}

#' Trim a single leading "_" and rename files accordingly
#' @export
fp_manifest_trim <- function(df, rename_files = TRUE, verbose = TRUE) {
  req <- c("motif", "score", "bound", "annot")
  miss <- setdiff(req, names(df)); if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "))

  strip1 <- function(s) sub("^_", "", s)
  fix_path <- function(p) {
    ifelse(is.na(p) | p == "", p, file.path(dirname(p), strip1(basename(p))))
  }
  do_rename <- function(old, new) {
    idx <- !is.na(old) & old != new
    if (!any(idx) || !rename_files) return(invisible(NULL))
    ok_src  <- file.exists(old[idx])
    tgt_ok  <- !file.exists(new[idx])
    todo    <- which(idx)[ok_src & tgt_ok]
    res <- if (length(todo)) file.rename(old[todo], new[todo]) else logical(0)
    if (verbose) {
      if (any(idx & !ok_src)) warning("Source missing: ", paste(old[idx & !ok_src], collapse = "; "))
      if (any(idx & !tgt_ok)) warning("Target exists: ",  paste(new[idx & !tgt_ok], collapse = "; "))
      if (length(res) && !all(res)) warning("Some renames failed.")
    }
  }

  df$motif <- strip1(df$motif)

  for (col in c("score","bound","annot")) {
    old <- df[[col]]
    new <- fix_path(old)
    do_rename(old, new)
    df[[col]] <- new
  }
  df
}

#' Trim leading "_" in all annotation CSVs listed in a manifest
#'
#' Reads each path in `manifest[[annot_col]]`, removes exactly one leading
#' underscore from the `motifs` column, and overwrites the file.
#'
#' @param manifest Tibble/data.frame with column `annot` (paths to CSVs).
#' @param annot_col Column name holding annotation CSV paths. Default: "annot".
#' @param motif_col Column name inside each CSV to fix. Default: "motifs".
#' @param n_workers Integer; if >1 and future.apply is available, process in parallel.
#' @param verbose Logical; print progress.
#' @return Tibble with path, n_rows, n_fixed (rows changed), and status.
#' @export
fp_manifest_trim_annots <- function(manifest,
                                 annot_col = "annot",
                                 motif_col = "motifs",
                                 n_workers = 1L,
                                 verbose = TRUE) {
  if (!is.data.frame(manifest) || !annot_col %in% names(manifest)) {
    .log_abort("`manifest` must contain column {.val {annot_col}}.")
  }

  paths <- unique(stats::na.omit(manifest[[annot_col]]))
  if (!length(paths)) {
    .log_abort("No valid paths found in {.val {annot_col}}.")
  }

  # worker fn
  fix_one <- function(p) {
    if (!file.exists(p)) {
      return(tibble::tibble(path = p, n_rows = NA_integer_, n_fixed = NA_integer_, status = "missing"))
    }
    df <- tryCatch(readr::read_csv(p, show_col_types = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) {
      return(tibble::tibble(path = p, n_rows = NA_integer_, n_fixed = NA_integer_, status = "read_error"))
    }
    if (!motif_col %in% names(df)) {
      return(tibble::tibble(path = p, n_rows = nrow(df), n_fixed = NA_integer_, status = "no_motif_col"))
    }
    old <- df[[motif_col]]
    new <- sub("^_", "", old)
    n_fixed <- sum(!is.na(old) & !is.na(new) & old != new)
    df[[motif_col]] <- new

    # overwrite
    ok <- tryCatch({ readr::write_csv(df, p); TRUE }, error = function(e) FALSE)
    tibble::tibble(path = p, n_rows = nrow(df), n_fixed = n_fixed, status = if (ok) "ok" else "write_error")
  }

  if (verbose) .log_inform("Fixing {length(paths)} annotation file{?s}...")

  run_parallel <- n_workers > 1L &&
    requireNamespace("future.apply", quietly = TRUE)
  if (run_parallel) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_workers)
    res <- future.apply::future_lapply(paths, fix_one, future.seed = FALSE)
  } else {
    res <- lapply(paths, fix_one)
  }

  out <- tibble::as_tibble(dplyr::bind_rows(res))
  if (verbose) {
    ok_n <- sum(out$status == "ok", na.rm = TRUE)
    .log_inform("Done: {ok_n}/{nrow(out)} OK. Changed rows (total): {sum(out$n_fixed, na.rm = TRUE)}.")
  }
  out
}

# Quantile normalize *by unique peak_ID*, then broadcast to all rows.
# - Keeps the input row order and N exactly the same.
# - Unique weighting prevents duplicated peak_ID from over-influencing the reference.
# - Tie-aware assignment within each column (avg of ref values across the tie block).
#' @export
qn_footprints <- function(fp_tbl, id_col = "peak_ID", tie_method = c("average","first")) {
  tie_method <- match.arg(tie_method)
  stopifnot(id_col %in% names(fp_tbl))

  samp <- setdiff(names(fp_tbl), id_col)
  if (!length(samp)) return(fp_tbl)

  # Original order mapping
  ids_all <- fp_tbl[[id_col]]
  # First occurrence of each peak_ID (we assume duplicates are identical rows from the alignment step)
  uniq_first <- !duplicated(ids_all)
  ids_unique <- ids_all[uniq_first]
  # Group index for broadcasting back
  grp <- match(ids_all, ids_unique)

  # Unique matrix (one row per unique peak)
  X <- as.matrix(fp_tbl[uniq_first, samp, drop = FALSE])
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  if (n == 0L || p == 0L) return(fp_tbl)

  # Build reference profile without constructing a full sorted matrix
  ref_sum   <- numeric(n)
  ref_count <- integer(n)

  # We'll also store per-column orders & non-NA counts for reuse
  ord_list <- vector("list", p)
  m_list   <- integer(p)

  for (j in seq_len(p)) {
    v <- X[, j]
    o <- order(v, na.last = TRUE)
    # number of non-NA values
    m <- sum(!is.na(v))
    ord_list[[j]] <- o
    m_list[j]     <- m
    if (m > 0L) {
      # Add to reference accumulators (first m positions)
      # we only need v sorted for the count (ref is just average of values at each rank)
      ref_sum[seq_len(m)]   <- ref_sum[seq_len(m)]   + v[o][seq_len(m)]
      ref_count[seq_len(m)] <- ref_count[seq_len(m)] + 1L
    }
  }

  # Mean across columns for each rank (ignore ranks with 0 count)
  ref <- ref_sum
  nz  <- ref_count > 0L
  ref[nz] <- ref_sum[nz] / ref_count[nz]
  ref[!nz] <- NA_real_

  # Map each column to the reference (tie-aware if requested)
  Xqn <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) {
    v <- X[, j]
    o <- ord_list[[j]]
    m <- m_list[j]
    if (m == 0L) next

    if (tie_method == "first") {
      # simple, fast: one-to-one rank mapping
      Xqn[o[seq_len(m)], j] <- ref[seq_len(m)]
    } else {
      # tie-aware: average the ref within each tie block
      s <- v[o[seq_len(m)]]
      r <- rle(s)
      ends   <- cumsum(r$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)
      # assign the mean(ref[start:end]) to all positions in that tie block
      for (k in seq_along(ends)) {
        idx <- starts[k]:ends[k]
        Xqn[o[idx], j] <- mean(ref[idx], na.rm = TRUE)
      }
    }
  }

  # Broadcast normalized unique values back to all rows WITHOUT materializing a big matrix
  out <- fp_tbl
  for (j in seq_len(p)) {
    out[[samp[j]]] <- Xqn[, j][grp]
  }
  out
}



# Save aligned/normalized FP results per motif and return a manifest
#' @export
save_footprints <- function(
    fp_aligned_normalized,
    out_dir  = NULL,
    threads  = max(1L, data.table::getDTthreads()),
    verbose  = TRUE
) {
  # pull tables
  score_dt <- data.table::as.data.table(fp_aligned_normalized$fp_score)
  bound_dt <- data.table::as.data.table(fp_aligned_normalized$fp_bound)
  annot_dt <- data.table::as.data.table(fp_aligned_normalized$fp_annotation)

  # sanity checks: same N and order
  stopifnot(nrow(score_dt) == nrow(bound_dt),
            nrow(bound_dt) == nrow(annot_dt),
            identical(score_dt$peak_ID, bound_dt$peak_ID),
            identical(score_dt$peak_ID, annot_dt$fp_peak),
            "motifs" %in% names(annot_dt))

  # ensure output dir
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  data.table::setDTthreads(threads)

  motifs <- unique(annot_dt$motifs)
  n_m <- length(motifs)

  # pre-allocate manifest fields
  n_rows_vec   <- integer(n_m)
  score_paths  <- character(n_m)
  bound_paths  <- character(n_m)
  annot_paths  <- character(n_m)

  if (verbose) .log_inform("Writing per-motif CSVs to {out_dir} ({n_m} motif{?s}) ...")

  for (i in seq_len(n_m)) {
    m <- motifs[i]
    idx <- which(annot_dt$motifs == m)

    # slice using the aligned row order
    sc <- score_dt[idx]
    bo <- bound_dt[idx]
    an <- annot_dt[idx]

    # file paths (match prior naming convention)
    base        <- file.path(out_dir, m)
    f_score     <- paste0(base, "_score.csv")
    f_bound     <- paste0(base, "_bound.csv")
    f_annot     <- paste0(base, "_annot.csv")

    # write
    data.table::fwrite(sc, f_score, nThread = threads, bom = FALSE)
    data.table::fwrite(bo, f_bound, nThread = threads, bom = FALSE)
    data.table::fwrite(an, f_annot, nThread = threads, bom = FALSE)

    # manifest rows
    n_rows_vec[i]  <- length(idx)
    score_paths[i] <- f_score
    bound_paths[i] <- f_bound
    annot_paths[i] <- f_annot
  }

  if (verbose) .log_inform("Done. Creating manifest tibble ...")

  manifest <- tibble::tibble(
    motif = motifs,
    n_rows = n_rows_vec,
    score  = score_paths,
    bound  = bound_paths,
    annot  = annot_paths
  )

  manifest
}

#
#' Filter GRN peaks by minimum number of bound samples
#'
#' Keeps rows where the count of bound samples in \code{fp_bound} meets
#' \code{min_bound}, and applies the same peak filter to the aligned
#' fp_score/fp_bound/fp_annotation tables.
#'
#' @param grn_set List containing \code{fp_score}, \code{fp_bound},
#'   and \code{fp_annotation}.
#' @param min_bound Integer threshold (default 1).
#' @param samples Optional character vector of sample IDs to consider;
#'   default uses all sample columns in \code{fp_bound}.
#' @param na_as_unbound Logical; treat NA as unbound (0) if TRUE, otherwise
#'   drop rows with any NA among the selected samples.
#' @param verbose Logical; print progress messages.
#' @param check_align Logical; if TRUE, assert alignment before filtering.
#'
#' @return A filtered \code{grn_set} list with aligned fp tables.
#' @export
filter_by_min_bound <- function(
    grn_set,
    min_bound     = 1L,
    samples       = NULL,
    na_as_unbound = TRUE,
    verbose       = TRUE,
    check_align   = FALSE
) {
  stopifnot(is.list(grn_set), "fp_bound" %in% names(grn_set))
  if (isTRUE(check_align)) {
    assert_fp_alignment(grn_set)
  }
  fb <- grn_set$fp_bound
  if (!"peak_ID" %in% names(fb))
    stop("`fp_bound` must contain a 'peak_ID' column.")

  # sample columns = all except key
  all_ids <- setdiff(names(fb), "peak_ID")
  if (is.null(samples)) {
    use_ids <- all_ids
  } else {
    miss <- setdiff(samples, all_ids)
    if (length(miss)) {
      stop("Some `samples` not found in fp_bound: ", paste(miss, collapse=", "))
    }
    use_ids <- samples
  }
  if (!length(use_ids)) stop("No sample columns selected to evaluate binding.")

  mat <- as.data.frame(fb[use_ids], stringsAsFactors = FALSE)

  # coerce to numeric 0/1
  mat[] <- lapply(mat, function(x) {
    if (is.logical(x)) as.integer(x)
    else if (is.numeric(x)) as.integer(x)
    else suppressWarnings(as.integer(x))
  })

  if (na_as_unbound) {
    mat[is.na(as.matrix(mat))] <- 0L
    keep <- rowSums(as.matrix(mat)) >= as.integer(min_bound)
  } else {
    # drop rows where any NA among selected samples
    any_na <- apply(as.matrix(mat), 1L, function(r) any(is.na(r)))
    keep   <- !any_na & rowSums(as.matrix(mat), na.rm = TRUE) >= as.integer(min_bound)
  }

  if (verbose) {
    message("Rows meeting ΓëÑ", min_bound, " bound samples among ",
            length(use_ids), " sample(s): ", sum(keep), "/", nrow(fb))
  }

  # delegate to the synced triplet filter to keep consistency
  filter_fp_rows(
    grn_set   = grn_set,
    peaks     = fb$peak_ID[keep],
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose   = verbose
  )
}


# Minimal update: fp_bound <- fp_bound & atac_overlap[ mapped atac_peak ]
update_fp_bound_by_overlap <- function(grn_set,
                                       key_fp   = "peak_ID",
                                       key_fpa  = "fp_peak",
                                       key_atac = "atac_peak",
                                       na_as_unbound = TRUE) {
  stopifnot(is.list(grn_set),
            all(c("fp_bound","fp_annotation","atac_overlap") %in% names(grn_set)))

  fb  <- grn_set$fp_bound
  fa  <- grn_set$fp_annotation
  aol <- grn_set$atac_overlap

  # sample IDs common to fp_bound and atac_overlap
  ids_fb  <- setdiff(names(fb),  key_fp)
  ids_ol  <- setdiff(names(aol), key_atac)
  ids_use <- intersect(ids_fb, ids_ol)
  if (!length(ids_use)) stop("No common sample IDs between fp_bound and atac_overlap.")

  # one-to-one map: fp_peak -> atac_peak (drop motif column)
  map <- dplyr::distinct(fa, !!rlang::sym(key_fpa), !!rlang::sym(key_atac))
  dup <- map |>
    dplyr::count(!!rlang::sym(key_fpa), name = "n") |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup)) stop("Each fp_peak must map to a single atac_peak.")

  # per-fp_peak overlap flags, aligned to fp_bound row order
  ol_fp <- tibble::tibble(!!key_fpa := fb[[key_fp]]) |>
    dplyr::left_join(map, by = setNames(key_fpa, key_fpa)) |>
    dplyr::left_join(
      dplyr::select(aol, !!key_atac, dplyr::all_of(ids_use)),
      by = setNames(key_atac, key_atac)
    )

  # matrices for pointwise logical AND (treat NA as 0 by default)
  fb_mat  <- as.matrix(fb[ids_use])
  ol_mat  <- as.matrix(ol_fp[ids_use])
  storage.mode(fb_mat) <- "integer"
  storage.mode(ol_mat) <- "integer"
  if (na_as_unbound) {
    fb_mat[is.na(fb_mat)] <- 0L
    ol_mat[is.na(ol_mat)] <- 0L
  }

  new_mat <- (fb_mat > 0L) & (ol_mat > 0L)

  # rebuild fp_bound with identical shape/order; keep any extra (non-overlap) columns unchanged
  out <- dplyr::bind_cols(
    fb[key_fp],
    tibble::as_tibble(matrix(as.integer(new_mat), nrow(fb), length(ids_use),
                             dimnames = list(NULL, ids_use)))
  )
  if (length(extra <- setdiff(ids_fb, ids_use))) {
    out <- dplyr::bind_cols(out, fb[extra])
    out <- out[, c(key_fp, setdiff(names(fb), key_fp)), drop = FALSE]
  }

  # final guards
  if (!identical(fb[[key_fp]], out[[key_fp]])) stop("Row order changed unexpectedly.")
  if (!setequal(names(fb), names(out)))        stop("Column set changed unexpectedly.")

  grn_set$fp_bound <- out
  grn_set
}

# Subset RNA to TFs and create a per-sample 0/1 flag
add_tf_expr_flags <- function(grn_set, threshold = 0, gt = TRUE) {
  stopifnot(is.list(grn_set),
            all(c("rna","tf_list") %in% names(grn_set)))
  rna <- grn_set$rna
  ids <- setdiff(names(rna), c("ensembl_gene_id","HGNC"))
  if (!length(ids)) stop("No RNA sample columns found.")

  rna_tf <- dplyr::filter(rna, .data$HGNC %in% grn_set$tf_list)

  flag_mat <- as.matrix(rna_tf[ids])
  if (gt) {
    flag_mat <- (flag_mat >  threshold)
  } else {
    flag_mat <- (flag_mat >= threshold)
  }
  storage.mode(flag_mat) <- "integer"

  grn_set$rna_tf       <- rna_tf
  grn_set$tf_expr_flag <- dplyr::bind_cols(
    rna_tf[c("ensembl_gene_id","HGNC")],
    tibble::as_tibble(flag_mat, .name_repair = "minimal")
  )
  grn_set
}


# Gate fp_bound by TF expression, handling motif dimers (A::B) properly
# Requires grn_set$tf_expr_flag (from add_tf_expr_flags) and grn_set$motif_db
update_fp_bound_by_tf_expr <- function(grn_set, group_size = 1L) {
  stopifnot(is.list(grn_set),
            all(c("fp_bound","fp_annotation","tf_expr_flag","motif_db") %in% names(grn_set)))

  fb   <- grn_set$fp_bound
  fpa  <- grn_set$fp_annotation
  tff  <- grn_set$tf_expr_flag

  # shared sample IDs (preserve fp_bound order)
  ids_fb  <- setdiff(names(fb), "peak_ID")
  ids_tf  <- setdiff(names(tff), c("ensembl_gene_id","HGNC"))
  ids     <- unique(intersect(ids_fb, ids_tf))
  if (!length(ids)) stop("No common sample IDs between fp_bound and tf_expr_flag.")

  # --- expand motifs -> TFs (handle dimers A::B) and join TF flags -------------
  gene_col <- .motif_gene_col(grn_set$motif_db)
  if (is.null(gene_col)) {
    stop("motif_db must include gene_symbol (or HGNC) for TF expression gating.")
  }
  mdb_exp <- grn_set$motif_db |>
    tidyr::separate_rows(dplyr::all_of(gene_col), sep = "::") |>
    dplyr::mutate(gene_symbol = stringr::str_trim(.data[[gene_col]])) |>
    dplyr::filter(.data$gene_symbol != "") |>
    dplyr::distinct(.data$motif, .data$gene_symbol)

  expanded <- fpa |>
    tidyr::separate_rows("motifs", sep = "\\s*,\\s*") |>
    dplyr::rename(motif = .data$motifs) |>
    dplyr::left_join(mdb_exp,           by = "motif") |>
    dplyr::left_join(dplyr::select(tff, "HGNC", dplyr::all_of(ids)),
                     by = c("gene_symbol" = "HGNC"))

  keep_row <- expanded |>
    dplyr::transmute(
      n_ok = rowSums(dplyr::across(dplyr::all_of(ids), ~ as.integer(.x == 1L)), na.rm = TRUE),
      ok   = .data$n_ok >= group_size
    ) |>
    dplyr::pull("ok")

  expanded$gene_symbol[!keep_row] <- NA_character_
  if (any(!keep_row)) expanded[!keep_row, ids] <- NA_integer_

  # --- collapse to one row per (fp_peak, atac_peak): ANY across TFs per sample ---
  data.table::setDT(expanded)
  collapsed_dt <- expanded[
    , c(
      as.list(vapply(.SD,
                     function(x) if (all(is.na(x))) NA_integer_
                     else as.integer(any(x == 1L, na.rm = TRUE)),
                     integer(1L))),
      list(tfs = { u <- unique(stats::na.omit(gene_symbol))
      if (length(u)) paste(u, collapse = ",") else NA_character_ })
    ),
    by = .(fp_peak, atac_peak),
    .SDcols = ids
  ]
  collapsed <- tibble::as_tibble(collapsed_dt)

  # peaks that still have at least one TF after gating
  keep_peaks <- collapsed |>
    dplyr::filter(!is.na(.data$tfs)) |>
    dplyr::pull("fp_peak")

  # --- AND with fp_bound (pointwise), preserving shape/order --------------------
  fb_keep <- dplyr::semi_join(fb, tibble::tibble(peak_ID = keep_peaks), by = "peak_ID")

  coll_aligned <- fb_keep["peak_ID"] |>
    dplyr::left_join(
      dplyr::select(collapsed, "fp_peak", dplyr::all_of(ids)),
      by = c("peak_ID" = "fp_peak")
    )

  fb_mat <- as.matrix(fb_keep[, ids, drop = FALSE]); storage.mode(fb_mat) <- "integer"
  tf_mat <- as.matrix(coll_aligned[, ids, drop = FALSE]); storage.mode(tf_mat) <- "integer"

  fb_new <- (fb_mat > 0L) & (tf_mat > 0L)                # logical matrix, same dims as fb_mat
  # SAFER rebuild: derive column names from the actual matrix we used
  out_mat <- matrix(as.integer(fb_new), nrow = nrow(fb_new), ncol = ncol(fb_new))
  colnames(out_mat) <- colnames(fb_mat)                  # <- avoids dimnames-length mismatch

  fb_out <- dplyr::bind_cols(
    fb_keep["peak_ID"],
    tibble::as_tibble(out_mat)
  )

  # add back any extra fp_bound columns (not part of ids), in original order
  extra <- setdiff(names(fb), c("peak_ID", ids))
  if (length(extra)) {
    fb_out <- dplyr::bind_cols(fb_out, fb_keep[extra]) |>
      dplyr::select("peak_ID", dplyr::all_of(setdiff(names(fb), "peak_ID")))
  }

  # --- update set (drop rows consistently; replace fp_bound; keep TF summary) ---
  grn_set <- filter_fp_rows(grn_set, peaks = keep_peaks)
  grn_set$fp_bound <- fb_out
  grn_set$fp_tfs   <- collapsed |>
    dplyr::filter(.data$fp_peak %in% keep_peaks) |>
    dplyr::select("fp_peak", "atac_peak", "tfs")

  grn_set
}

filter_fp_rows <- function(
    grn_set,
    peaks     = NULL,
    predicate = NULL,
    score_key = "peak_ID",
    bound_key = "peak_ID",
    annot_key = "fp_peak",
    verbose   = TRUE,
    warn_on_mismatch = TRUE
) {
  stopifnot(is.list(grn_set),
            all(c("fp_score","fp_bound","fp_annotation") %in% names(grn_set)))

  fp_score <- grn_set$fp_score
  fp_bound <- grn_set$fp_bound
  fp_annot <- grn_set$fp_annotation

  s_ids <- fp_score[[score_key]]
  b_ids <- fp_bound[[bound_key]]
  a_ids <- fp_annot[[annot_key]]

  if (!is.character(s_ids) || !is.character(b_ids) || !is.character(a_ids)) {
    stop("Keys must be character; coerce before calling if needed.")
  }

  # --- NEW: if all three are already empty, just handle requests sanely and return
  if (length(s_ids) == 0L && length(b_ids) == 0L && length(a_ids) == 0L) {
    # Nothing to reconcile; respect peaks/predicate but result will remain empty.
    if (!is.null(predicate)) {
      keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
      if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
        stop("`predicate` must return a logical vector of length nrow(fp_score).")
      }
    }
    # If user passed non-empty peaks, warn that none can match.
    if (!is.null(peaks) && length(peaks) > 0L) {
      if (requireNamespace("cli", quietly = TRUE)) {
        .log_warn("Requested {length(peaks)} peak ID(s), but the current set is empty; returning empty triplet.")
      } else {
        warning("Requested peaks in an empty set; returning empty triplet.")
      }
    }
    if (verbose) message("Filtered to 0 peaks (from 0).")
    grn_set$fp_score      <- fp_score[0, , drop = FALSE]
    grn_set$fp_bound      <- fp_bound[0, , drop = FALSE]
    grn_set$fp_annotation <- fp_annot[0, , drop = FALSE]
    return(grn_set)
  }

  # --- Reconcile key universe FIRST (align all three to common IDs; keep annot order)
  common_ids <- Reduce(intersect, list(unique(s_ids), unique(b_ids), unique(a_ids)))
  if (length(common_ids) == 0L) {
    stop("No common keys across fp_score/fp_bound/fp_annotation.")
  }
  if (!(setequal(s_ids, b_ids) && setequal(s_ids, a_ids))) {
    n_drop <- c(
      score_only = sum(!(s_ids %in% common_ids)),
      bound_only = sum(!(b_ids %in% common_ids)),
      annot_only = sum(!(a_ids %in% common_ids))
    )
    if (isTRUE(warn_on_mismatch)) {
      if (requireNamespace("cli", quietly = TRUE)) {
        .log_warn(c(
          "!" = "Key sets differ across tables; reconciling to their intersection ({length(common_ids)} rows).",
          "i" = "Dropped: score-only={n_drop['score_only']}, bound-only={n_drop['bound_only']}, annot-only={n_drop['annot_only']}."
        ))
      } else {
        warning(sprintf(
          "Key sets differ; aligning to intersection (%d rows). Dropped score-only=%d, bound-only=%d, annot-only=%d.",
          length(common_ids), n_drop["score_only"], n_drop["bound_only"], n_drop["annot_only"]
        ))
      }
    }

    # anchor order to annotation (deterministic)
    a_keep <- fp_annot[[annot_key]] %in% common_ids
    a_ord  <- fp_annot[[annot_key]][a_keep]

    s_keep <- fp_score[[score_key]] %in% common_ids
    b_keep <- fp_bound[[bound_key]] %in% common_ids

    fp_annot <- fp_annot[a_keep, , drop = FALSE]

    fp_score <- fp_score[s_keep, , drop = FALSE]
    idx_s    <- match(a_ord, fp_score[[score_key]])
    fp_score <- fp_score[idx_s, , drop = FALSE]

    fp_bound <- fp_bound[b_keep, , drop = FALSE]
    idx_b    <- match(a_ord, fp_bound[[bound_key]])
    fp_bound <- fp_bound[idx_b, , drop = FALSE]

    s_ids <- fp_score[[score_key]]
    b_ids <- fp_bound[[bound_key]]
    a_ids <- fp_annot[[annot_key]]

    stopifnot(identical(s_ids, a_ids), identical(b_ids, a_ids))
  }

  # --- Determine peaks to keep
  if (!is.null(predicate)) {
    keep_lgl <- predicate(fp_score, fp_bound, fp_annot)
    if (!is.logical(keep_lgl) || length(keep_lgl) != nrow(fp_score)) {
      stop("`predicate` must return a logical vector of length nrow(fp_score).")
    }
    keep_ids <- s_ids[keep_lgl %in% TRUE]
  } else if (!is.null(peaks)) {
    if (!all(peaks %in% s_ids)) {
      miss <- setdiff(peaks, s_ids)
      if (length(miss)) {
        if (isTRUE(warn_on_mismatch)) {
          if (requireNamespace("cli", quietly = TRUE)) {
            .log_warn("Some requested `peaks` not present after reconciliation; dropping {length(miss)} missing id(s).")
          } else {
            warning(sprintf("Some requested `peaks` are missing after reconciliation; dropping %d id(s).", length(miss)))
          }
        }
      }
    }
    peaks_uniq <- peaks[!duplicated(peaks)]
    keep_ids   <- peaks_uniq[peaks_uniq %in% s_ids]
  } else {
    if (verbose) message("No filter applied; returning original set.")
    grn_set$fp_score      <- fp_score
    grn_set$fp_bound      <- fp_bound
    grn_set$fp_annotation <- fp_annot
    return(grn_set)
  }

  if (!length(keep_ids)) {
    if (verbose) message("Filtered to 0 peaks (from ", length(s_ids), ").")
    fp_score2 <- fp_score[0, , drop = FALSE]
    fp_bound2 <- fp_bound[0, , drop = FALSE]
    fp_annot2 <- fp_annot[0, , drop = FALSE]
  } else {
    keep_s <- fp_score[[score_key]] %in% keep_ids
    keep_b <- fp_bound[[bound_key]] %in% keep_ids
    keep_a <- fp_annot[[annot_key]] %in% keep_ids

    fp_score2 <- fp_score[keep_s, , drop = FALSE]
    fp_bound2 <- fp_bound[keep_b, , drop = FALSE]
    fp_annot2 <- fp_annot[keep_a, , drop = FALSE]

    ord_a <- fp_annot2[[annot_key]]
    if (!is.null(peaks)) {
      target <- keep_ids
      ord_idx <- match(ord_a, target)
      o <- order(ord_idx, seq_along(ord_idx), na.last = TRUE)
      fp_annot2 <- fp_annot2[o, , drop = FALSE]
      m_s <- match(fp_annot2[[annot_key]], fp_score2[[score_key]])
      m_b <- match(fp_annot2[[annot_key]], fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    } else {
      m_s <- match(ord_a, fp_score2[[score_key]])
      m_b <- match(ord_a, fp_bound2[[bound_key]])
      fp_score2 <- fp_score2[m_s, , drop = FALSE]
      fp_bound2 <- fp_bound2[m_b, , drop = FALSE]
    }
  }

  s2 <- fp_score2[[score_key]]
  b2 <- fp_bound2[[bound_key]]
  a2 <- fp_annot2[[annot_key]]
  if (!(identical(s2, a2) && identical(b2, a2))) {
    ks <- head(s2[!(s2 %in% a2)], 3)
    kb <- head(b2[!(b2 %in% a2)], 3)
    stop(paste0(
      "Post-check failed: filtered key sets are not identical across the three tables. ",
      "Examples ΓÇö score-only: ", paste(ks, collapse = ", "),
      "; bound-only: ", paste(kb, collapse = ", ")
    ))
  }

  if (verbose) {
    message("Filtered to ", nrow(fp_annot2), " peaks (from ", length(s_ids), ").")
  }

  grn_set$fp_score       <- fp_score2
  grn_set$fp_bound       <- fp_bound2
  grn_set$fp_annotation  <- fp_annot2
  grn_set
}
# per-motif filtering
# BEGIN EDIT: per-motif processor + parallel driver (uses update_fp_bound_by_overlap)

#' Process ONE motif from fp_manifest (strict set) and write filtered CSVs
#'
#' Reads the three per-motif files (score/bound/annotation), builds a strict GRN set
#' using provided `build_args`, applies the standard filtering/update steps, and writes
#' filtered outputs to `out_dir` as "<motif>_{score,bound,annotation}.csv".
#'
#' Required columns in `entry`: motif, score, bound, annot (absolute or relative paths).
#' `update_fp_bound_by_overlap()` is used as-is.
#'
#' @param entry Single-row tibble/data.frame from fp_manifest.
#' @param out_dir Output directory (e.g., "inst/extdata/fp_filtered_jaspar2024").
#' @param build_args Named list of args passed to `build_grn_set()` in addition to fp_*.
#'                   Should include: atac_score, atac_overlap, rna, metadata, tf_list,
#'                   motif_db, label_col, expected_n.
#' @param threshold_tf_expr Numeric; forwarded to `add_tf_expr_flags()`.
#' @param skip_existing If TRUE and all three outputs exist, skip work.
#' @param verbose Print progress.
#' @return Tibble with motif, n_rows (filtered fp_bound rows), and output paths.
process_one_motif_from_manifest <- function(entry,
                                            out_dir,
                                            build_args,
                                            threshold_tf_expr = 10,
                                            skip_existing = TRUE,
                                            verbose = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  if (!all(need_cols %in% names(entry))) {
    .log_abort("`entry` must contain columns: {cli::fmt_cols(setdiff(need_cols, names(entry)))}")
  }
  m <- as.character(entry$motif)[1]
  f_score <- as.character(entry$score)[1]
  f_bound <- as.character(entry$bound)[1]
  f_annot <- as.character(entry$annot)[1]

  if (!file.exists(f_score) || !file.exists(f_bound) || !file.exists(f_annot)) {
    .log_abort("Missing file(s) for motif {m}:\n  score={f_score}\n  bound={f_bound}\n  annot={f_annot}")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_score <- file.path(out_dir, paste0(m, "_score.csv"))
  out_bound <- file.path(out_dir, paste0(m, "_bound.csv"))
  out_annot <- file.path(out_dir, paste0(m, "_annotation.csv"))

  if (skip_existing && file.exists(out_score) && file.exists(out_bound) && file.exists(out_annot)) {
    return(tibble::tibble(motif = m, n_rows = NA_integer_,
                          score = out_score, bound = out_bound, annot = out_annot))
  }

  if (verbose) .log_inform("Processing motif {m} ...")

  fp_score <- readr::read_csv(f_score, show_col_types = FALSE)
  fp_bound <- readr::read_csv(f_bound, show_col_types = FALSE)
  fp_annot <- readr::read_csv(f_annot, show_col_types = FALSE)

  if (!all(c("fp_peak","atac_peak","motifs") %in% names(fp_annot))) {
    if (all(c("peak_ID","peak_ATAC","TFBS_name") %in% names(fp_annot))) {
      fp_annot <- dplyr::rename(fp_annot, fp_peak = "peak_ID", atac_peak = "peak_ATAC", motifs = "TFBS_name")
    } else {
      .log_abort("Annotation for {m} must have fp_peak/atac_peak/motifs (or peak_ID/peak_ATAC/TFBS_name).")
    }
  }

  ss <- do.call(build_grn_set, c(list(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annot
  ), build_args))

  # --- EARLY-EXIT GUARDS (no writes if empty) -----------------
  empty_row <- function() tibble::tibble(motif = m, n_rows = 0L,
                                         score = NA_character_,
                                         bound = NA_character_,
                                         annot = NA_character_)

  ss <- filter_by_min_bound(ss, min_bound = 1L)
  if (nrow(ss$fp_bound) == 0L) { if (verbose) .log_inform("Motif {m}: empty after min-bound filter; skipping."); return(empty_row()) }

  ss <- update_fp_bound_by_overlap(ss)
  ss <- filter_by_min_bound(ss, min_bound = 1L)
  if (nrow(ss$fp_bound) == 0L) { if (verbose) .log_inform("Motif {m}: empty after overlap gating; skipping."); return(empty_row()) }

  ss <- add_tf_expr_flags(ss, threshold = threshold_tf_expr)
  ss2 <- update_fp_bound_by_tf_expr(ss, group_size = 1L)
  ss2 <- filter_by_min_bound(ss2, min_bound = 1L)
  if (nrow(ss2$fp_bound) == 0L) { if (verbose) .log_inform("Motif {m}: empty after TF expression gating; skipping."); return(empty_row()) }
  # ------------------------------------------------------------

  # Only compute correlations if we still have rows
  # ss2 <- annotate_fp_tf_corr_one_motif(ss2, cor_method = "pearson", min_non_na = 5L) # do not calculate correlation for now, will calculate after quantile normalization

  out_score_tbl <- ss2$fp_score
  out_bound_tbl <- ss2$fp_bound
  out_annot_tbl <- ss2$fp_annotation

  if (!"peak_ID" %in% names(out_score_tbl) && "fp_peak" %in% names(out_score_tbl)) {
    out_score_tbl <- dplyr::rename(out_score_tbl, peak_ID = "fp_peak")
  }
  if (!"peak_ID" %in% names(out_bound_tbl) && "fp_peak" %in% names(out_bound_tbl)) {
    out_bound_tbl <- dplyr::rename(out_bound_tbl, peak_ID = "fp_peak")
  }
  if (!all(c("fp_peak","atac_peak","motifs") %in% names(out_annot_tbl))) {
    out_annot_tbl <- out_annot_tbl |>
      dplyr::rename_with(~"fp_peak",  dplyr::any_of("peak_ID")) |>
      dplyr::rename_with(~"atac_peak",dplyr::any_of("peak_ATAC")) |>
      dplyr::rename_with(~"motifs",   dplyr::any_of("TFBS_name"))
  }

  readr::write_csv(out_score_tbl, out_score)
  readr::write_csv(out_bound_tbl, out_bound)
  readr::write_csv(out_annot_tbl, out_annot)

  tibble::tibble(motif = m,
                 n_rows = nrow(out_bound_tbl),
                 score  = out_score,
                 bound  = out_bound,
                 annot  = out_annot)
}

#' Parallel driver: map `process_one_motif_from_manifest()` over fp_manifest
#'
#' @param fp_manifest Tibble with columns motif, score, bound, annot.
#' @param out_dir Output directory for filtered CSVs.
#' @param build_args Named list passed to `build_grn_set()` (see above).
#' @param motif_ids Optional subset of motifs to process.
#' @param n_workers Parallel workers (motif-level). Default: all logical cores.
#' @param set_plan Whether to set/reset a future plan automatically.
#' @param skip_existing If TRUE, skip motifs whose three outputs already exist.
#' @param threshold_tf_expr Numeric threshold for `add_tf_expr_flags()`.
#' @param verbose Print progress.
#' @return Tibble manifest of processed motifs.
#'
#' @export
filter_footprints <- function(fp_manifest,
                                       out_dir,
                                       build_args,
                                       motif_ids = NULL,
                                       n_workers = max(1L, parallel::detectCores(TRUE)),
                                       set_plan = TRUE,
                                       skip_existing = TRUE,
                                       threshold_tf_expr = 10,
                                       verbose = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  if (!isTRUE(is.data.frame(fp_manifest)) || !all(need_cols %in% names(fp_manifest))) {
    .log_abort("`fp_manifest` must be a data.frame with columns: {cli::fmt_cols(need_cols)}")
  }

  rows <- fp_manifest
  if (!is.null(motif_ids)) {
    rows <- dplyr::semi_join(rows, tibble::tibble(motif = motif_ids), by = "motif")
    if (!nrow(rows)) .log_abort("No matching motifs found in `fp_manifest` for the requested subset.")
  }

  if (verbose) .log_inform("Submitting {nrow(rows)} motif(s) to {n_workers} worker(s)...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    is_rstudio <- identical(Sys.getenv("RSTUDIO"), "1") ||
      (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !is_rstudio &&
        isTRUE(requireNamespace("parallelly", quietly = TRUE)) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore,   workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
  }

  res <- future.apply::future_lapply(seq_len(nrow(rows)), function(i) {
    entry <- rows[i, , drop = FALSE]
    process_one_motif_from_manifest(
      entry               = entry,
      out_dir             = out_dir,
      build_args          = build_args,
      threshold_tf_expr   = threshold_tf_expr,
      skip_existing       = skip_existing,
      verbose             = FALSE
    )
  }, future.seed = FALSE)

  tibble::as_tibble(dplyr::bind_rows(res))
}




#' Combine FP correlation tables (split dimers), run BH FDR once across all,
#' sanity-check p vs p_adj, filter, and align score/bound to the kept annotations.
#'
#' Robust, minimal pipeline:
#' 1) Normalizes columns when reading each file.
#' 2) Splits comma-separated TF/r/p columns row-wise (handles dimers) with safe padding.
#' 3) Runs **one** BenjaminiΓÇôHochberg adjustment across the combined annotations.
#' 4) Enforces `p_adj >= p` and reports monotonicity stats.
#' 5) Filters by `r_thr` (positive) and `p_thr` (FDR).
#' 6) Reads score/bound only for peaks that pass, de-dups on `peak_ID`, and aligns order.
#'
#' @param fp_filtered_manifest Tibble/data.frame with columns:
#'   `motif, n_rows, score, bound, annot` (CSV paths).
#' @param p_thr Numeric FDR threshold (default 0.05).
#' @param r_thr Numeric correlation threshold (keep r > r_thr; default 0.3).
#' @param verbose Logical; print progress (default TRUE).
#' @param threads Integer; data.table threads (default = all detected cores).
#' @param output_bed Optional directory path. If provided, the function will write per-TF
#'   BED-like files: `<TF>_all.bed`, `<TF>_bound.bed`, `<TF>_unbound.bed` (no headers),
#'   and `<TF>_overview.txt` (headered; same columns as `_all.bed` plus a `_bound` 0/1 column).
#'   Columns: TFBS_chr, TFBS_start, TFBS_end, TFBS_name, peak_chr, peak_start, peak_end,
#'   TF, corr_fp_tf_r, corr_fp_tf_p, corr_fp_tf_p_adj; bound/unbound defined by r/p filters.
#'
#' @return A list of tibbles: `fp_score`, `fp_bound`, `fp_annotation`
#'         (all aligned to the same row order).
#'
#' @examples
#' # res <- combine_filtered_fp_tables_split_dimers_simple(fp_corr_manifest)
combine_filtered_fp_tables_split_dimers_simple <- function(
    fp_filtered_manifest,
    p_thr   = 0.05,
    r_thr   = 0.3,
    verbose = TRUE,
    threads = max(1L, parallel::detectCores(TRUE)),
    # BEGIN EDIT: add optional output directory for BED exports
    output_bed = NULL
    # END EDIT
) {
  # ---- checks ----
  if (!is.data.frame(fp_filtered_manifest) ||
      !all(c("motif","n_rows","score","bound","annot") %in% names(fp_filtered_manifest))) {
    .log_abort("`fp_filtered_manifest` must have columns: motif, n_rows, score, bound, annot.")
  }

  use <- fp_filtered_manifest |>
    dplyr::filter(
      !is.na(.data$n_rows), .data$n_rows > 0,
      !is.na(.data$score), file.exists(.data$score),
      !is.na(.data$bound), file.exists(.data$bound),
      !is.na(.data$annot), file.exists(.data$annot)
    )

  if (!nrow(use)) {
    if (verbose) .log_inform("No non-empty motifs with existing files; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }
  if (verbose) .log_inform("Combining {nrow(use)} motif{?s} ...")

  # ---- threading ----
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(threads))

  fread_fast <- function(path) data.table::fread(path, nThread = threads, showProgress = FALSE)

  # ---- helpers ----
  split_commas <- function(x) {
    x <- as.character(x)
    lapply(x, function(s) {
      if (is.na(s) || !nzchar(s)) character(0)
      else strsplit(s, "\\s*,\\s*", perl = TRUE)[[1]]
    })
  }
  pad_to <- function(lst, lens, fill) {
    out <- vector("list", length(lst))
    for (i in seq_along(lst)) {
      v <- lst[[i]]
      # recycle singleton
      if (length(v) == 1L && lens[i] > 1L) v <- rep(v, lens[i])
      # right-pad
      if (length(v) < lens[i]) v <- c(v, rep(fill, lens[i] - length(v)))
      out[[i]] <- v
    }
    unlist(out, use.names = FALSE)
  }
  safe_num <- function(x) {
    x <- as.character(x)
    x[!nzchar(x)] <- NA_character_   # treat "" as NA
    suppressWarnings(as.numeric(x))  # quiet "NAs introduced by coercion"
  }

  # ---- read one annotation, normalize names, split dimers ----
  read_annot_split <- function(p) {
    dt <- fread_fast(p)

    # Normalize names to: fp_peak / atac_peak / motifs / tfs / corr_fp_tf_r / corr_fp_tf_p
    if (!"fp_peak"   %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID",   "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"    %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    if (!"tfs"       %in% names(dt)) dt[["tfs"]] <- NA_character_
    if (!"corr_fp_tf_r" %in% names(dt)) dt[["corr_fp_tf_r"]] <- NA_character_
    if (!"corr_fp_tf_p" %in% names(dt)) dt[["corr_fp_tf_p"]] <- NA_character_

    tfs_list <- split_commas(dt[["tfs"]])
    r_list   <- lapply(split_commas(dt[["corr_fp_tf_r"]]), safe_num)
    p_list   <- lapply(split_commas(dt[["corr_fp_tf_p"]]), safe_num)

    # how many repeats per row
    nrep <- pmax(lengths(tfs_list), lengths(r_list), lengths(p_list))
    nrep[nrep == 0L] <- 1L
    idx <- rep.int(seq_len(nrow(dt)), nrep)

    data.table::data.table(
      fp_peak      = dt$fp_peak[idx],
      atac_peak    = dt$atac_peak[idx],
      motifs       = dt$motifs[idx],
      tfs          = pad_to(tfs_list, nrep, NA_character_),
      corr_fp_tf_r = safe_num(pad_to(r_list, nrep, NA_real_)),
      corr_fp_tf_p = safe_num(pad_to(p_list, nrep, NA_real_))
    )
  }

  # ---- 1) combine all annotations ----
  annot_dt <- data.table::rbindlist(lapply(use$annot, read_annot_split),
                                    use.names = TRUE, fill = TRUE)

  # Optional quick NA summary
  if (verbose) {
    n_p_na <- sum(!is.finite(annot_dt$corr_fp_tf_p))
    n_r_na <- sum(!is.finite(annot_dt$corr_fp_tf_r))
    .log_inform("Pre-filter NA counts: p={n_p_na}, r={n_r_na}")
  }

  # ---- sanity: p in [0,1] (raw) ----
  bad_gt1 <- sum(annot_dt$corr_fp_tf_p > 1, na.rm = TRUE)
  bad_lt0 <- sum(annot_dt$corr_fp_tf_p < 0, na.rm = TRUE)
  if (bad_gt1 || bad_lt0) {
    .log_abort("[raw p] Found {bad_gt1} values > 1 and {bad_lt0} values < 0. Please fix upstream.")
  }

  # ---- pre-BH: dedup (fp_peak, tfs) by keeping the FIRST occurrence ----
  if (verbose) .log_inform("Pre-BH: rows before pair-dedup: {nrow(annot_dt)}")

  keep_idx <- !duplicated(annot_dt[, .(fp_peak, tfs)])
  n_drop   <- sum(!keep_idx)

  annot_dt <- annot_dt[keep_idx]

  if (verbose) {
    .log_inform("Pre-BH: dropped {n_drop} duplicate (fp_peak, tfs) pairs (kept first).")
    .log_inform("Pre-BH: rows after pair-dedup: {nrow(annot_dt)}")
  }

  # ---- BH adjust across ALL rows at once (finite p only) ----
  good_p <- is.finite(annot_dt$corr_fp_tf_p)
  annot_dt[, corr_fp_tf_p_adj := NA_real_]
  if (any(good_p)) {
    annot_dt[good_p, corr_fp_tf_p_adj := stats::p.adjust(annot_dt$corr_fp_tf_p[good_p], method = "BH")]
    # enforce theoretical constraint (defensive against FP/rounding)
    annot_dt[good_p, corr_fp_tf_p_adj := pmax(corr_fp_tf_p_adj, corr_fp_tf_p)]
  }

  # ---- post-BH sanity stats ----
  cmp <- data.table::fcase(
    !is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj), NA_character_,
    annot_dt$corr_fp_tf_p_adj >  annot_dt$corr_fp_tf_p, ">",
    annot_dt$corr_fp_tf_p_adj == annot_dt$corr_fp_tf_p, "==",
    annot_dt$corr_fp_tf_p_adj <  annot_dt$corr_fp_tf_p, "<"
  )
  if (verbose) {
    .log_inform(c(
      "BH sanity (p_adj vs p):",
      "  >  = {sum(cmp == '>', na.rm = TRUE)}",
      "  == = {sum(cmp == '==', na.rm = TRUE)}",
      "  <  = {sum(cmp == '<', na.rm = TRUE)}  (after pmax() this should be 0)",
      "  NA = {sum(!is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj))}",
      "  total rows = {nrow(annot_dt)}"
    ))
  }

  # BEGIN EDIT: snapshot pre-filter annotations for BED exports
  annot_pre_filter <- data.table::copy(annot_dt)
  # END EDIT

  # ---- 2) filter: positive r and FDR ----
  keep <- is.finite(annot_dt$corr_fp_tf_r) &
    is.finite(annot_dt$corr_fp_tf_p_adj) &
    (annot_dt$corr_fp_tf_r > r_thr) &
    (annot_dt$corr_fp_tf_p_adj < p_thr)
  annot_dt <- annot_dt[keep]
  if (verbose) {
    .log_inform("After corr/FDR filter: {nrow(annot_dt)} rows (unique peaks: {length(unique(annot_dt$fp_peak))}).")
  }
  if (nrow(annot_dt) == 0L) {
    # BEGIN EDIT: still optionally create empty per-TF files if requested
    if (!is.null(output_bed)) {
      if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)
      if (verbose) .log_inform("No rows passed; requested BED export directory initialized at {output_bed}.")
    }
    # END EDIT
    if (verbose) .log_inform("Nothing passed thresholds; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }

  peaks_keep <- unique(annot_dt$fp_peak)

  # ---- helpers to normalize score/bound ----
  norm_score <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.numeric(dt[[j]])))
    dt[]
  }
  norm_bound <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.integer(safe_num(dt[[j]]))))
    dt[]
  }

  read_score_subset <- function(p) {
    dt <- norm_score(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }
  read_bound_subset <- function(p) {
    dt <- norm_bound(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }

  score_dt <- data.table::rbindlist(lapply(use$score, read_score_subset), use.names = TRUE, fill = TRUE)
  bound_dt <- data.table::rbindlist(lapply(use$bound, read_bound_subset), use.names = TRUE, fill = TRUE)

  # de-dup by first occurrence
  if (nrow(score_dt)) score_dt <- score_dt[!duplicated(peak_ID)]
  if (nrow(bound_dt)) bound_dt <- bound_dt[!duplicated(peak_ID)]

  # ---- 3) align score/bound EXACTLY to annotation order ----
  key_vec <- annot_dt$fp_peak
  key_dt  <- data.table::data.table(peak_ID = key_vec)

  fp_score_dt <- score_dt[key_dt, on = "peak_ID"]  # preserves key_dt order
  fp_bound_dt <- bound_dt[key_dt, on  = "peak_ID"]

  # final sanity: 1:1 and same order
  stopifnot(
    nrow(fp_score_dt) == nrow(annot_dt),
    nrow(fp_bound_dt) == nrow(annot_dt),
    identical(fp_score_dt$peak_ID, key_vec),
    identical(fp_bound_dt$peak_ID, key_vec)
  )

  # BEGIN EDIT: optional per-TF BED & overview export
  if (!is.null(output_bed)) {
    if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)

    .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
    .mk_table <- function(dt) {
      if (!nrow(dt)) return(data.table::data.table())
      dt2 <- dt[!is.na(tfs) & nzchar(tfs) & !is.na(fp_peak) & !is.na(atac_peak)]
      if (!nrow(dt2)) return(data.table::data.table())
      fp_s <- data.table::tstrsplit(dt2$fp_peak, "[:-]", perl = TRUE)
      at_s <- data.table::tstrsplit(dt2$atac_peak, "[:-]", perl = TRUE)
      data.table::data.table(
        TFBS_chr         = fp_s[[1]],
        TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
        TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
        TFBS_name        = dt2$motifs,
        peak_chr         = at_s[[1]],
        peak_start       = suppressWarnings(as.integer(at_s[[2]])),
        peak_end         = suppressWarnings(as.integer(at_s[[3]])),
        TF               = dt2$tfs,
        corr_fp_tf_r     = dt2$corr_fp_tf_r,
        corr_fp_tf_p     = dt2$corr_fp_tf_p,
        corr_fp_tf_p_adj = dt2$corr_fp_tf_p_adj
      )
    }

    all_dt    <- .mk_table(annot_pre_filter)
    bound_dt2 <- .mk_table(annot_pre_filter[keep])

    if (nrow(all_dt))   all_dt[,   .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]
    if (nrow(bound_dt2)) bound_dt2[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]

    tfs_vec <- sort(unique(all_dt$TF))
    for (tf in tfs_vec) {
      tf_slug <- .slug(tf)
      f_all      <- file.path(output_bed, paste0(tf_slug, "_all.bed"))
      f_bound    <- file.path(output_bed, paste0(tf_slug, "_bound.bed"))
      f_unbound  <- file.path(output_bed, paste0(tf_slug, "_unbound.bed"))
      f_overview <- file.path(output_bed, paste0(tf_slug, "_overview.txt"))

      sub_all   <- all_dt[TF == tf]
      sub_bound <- bound_dt2[TF == tf]

      # --- FIX 1: use %in% (not %chin%) ---
      ids_bound   <- if (nrow(sub_bound)) sub_bound$.id else character(0)
      sub_unbound <- sub_all[!(.id %in% ids_bound)]

      # --- FIX 2: drop ".id" (not "._id"); keep=FALSE selection with data.table ---
      data.table::fwrite(sub_all[,   !".id", with = FALSE], f_all,     sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_bound[, !".id", with = FALSE], f_bound,   sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_unbound[, !".id", with = FALSE], f_unbound, sep = "\t", col.names = FALSE, quote = FALSE)

      # headered overview with _bound flag
      if (nrow(sub_all)) {
        sub_all[, `_bound` := as.integer(.id %in% ids_bound)]   # <-- FIX 1 applied here too
        data.table::fwrite(sub_all[, !".id", with = FALSE], f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      } else {
        hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
                 "peak_chr","peak_start","peak_end","TF",
                 "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
        # create a 0-row table with the right header
        empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode="character", length=0)), length(hdr)), hdr))
        data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      }
    }
    if (verbose) .log_inform("Per-TF BED exports written to: {output_bed}")
  }
  # END EDIT


  # ---- return ----
  list(
    fp_score      = tibble::as_tibble(fp_score_dt),
    fp_bound      = tibble::as_tibble(fp_bound_dt),
    fp_annotation = tibble::as_tibble(annot_dt)
  )
}



#' Build ATAC score/overlap tibbles; sort by chr/start/end first (if enabled).
#'
#' @param atac_data data.frame/tibble with coords, metadata, numeric samples, and overlap cols.
#' @param chr_col,start_col,end_col coordinate column names (defaults: "Chr","Start","End").
#' @param overlap_prefix prefix for overlap columns (default "Overlap_").
#' @param metadata_cols metadata columns to exclude from samples (besides coords).
#' @param sort_peaks logical: sort rows by chr/start/end before processing (default TRUE).
#' @return list(score = atac_score, overlap = atac_overlap)
#' @export
load_atac <- function(
    atac_data,
    chr_col        = "Chr",
    start_col      = "Start",
    end_col        = "End",
    overlap_prefix = "Overlap_",
    metadata_cols  = c("PeakID","Distance.to.TSS","Nearest.PromoterID","Entrez.ID",
                       "Nearest.Unigene","Nearest.Refseq","Nearest.Ensembl",
                       "Gene.Name","Gene.Alias","Gene.Description","Gene.Type"),
    sort_peaks     = TRUE
) {
  if (!is.data.frame(atac_data)) .log_abort("`atac_data` must be a data.frame/tibble.")
  need <- c(chr_col, start_col, end_col)
  if (!all(need %in% names(atac_data)))
    .log_abort("Missing coord columns: {.val {setdiff(need, names(atac_data))}}.")

  df <- atac_data

  # ---- sort by chr/start/end first ----
  if (isTRUE(sort_peaks)) {
    chr <- tolower(gsub("^chr", "", as.character(df[[chr_col]])))
    chr_num <- suppressWarnings(as.integer(chr))
    chr_num[is.na(chr_num) & chr %in% c("x")]  <- 23L
    chr_num[is.na(chr_num) & chr %in% c("y")]  <- 24L
    chr_num[is.na(chr_num) & chr %in% c("m","mt","dmel_mito","mito")] <- 25L
    ord <- order(is.na(chr_num), chr_num, chr,
                 as.integer(round(df[[start_col]])),
                 as.integer(round(df[[end_col]])),
                 na.last = TRUE, method = "radix")
    df <- df[ord, , drop = FALSE]
  }

  # unified peak key
  atac_peak <- paste0(
    df[[chr_col]], ":", as.integer(round(df[[start_col]])),
    "-", as.integer(round(df[[end_col]]))
  )

  # detect columns
  overlap_cols <- grep(paste0("^", overlap_prefix), names(df), value = TRUE)
  meta_set     <- unique(c(metadata_cols, chr_col, start_col, end_col))
  candidates   <- setdiff(names(df), c(meta_set, overlap_cols))
  sample_cols  <- candidates[vapply(df[candidates], is.numeric, logical(1))]
  if (!length(sample_cols))
    .log_abort("No numeric sample columns found after excluding metadata and {.val {overlap_prefix}}*.")

  # outputs
  atac_score <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, sample_cols, drop = FALSE])

  atac_overlap <- dplyr::tibble(atac_peak = atac_peak) |>
    dplyr::bind_cols(df[, overlap_cols, drop = FALSE])

  # strip prefix from overlap column names
  if (length(overlap_cols)) {
    new_names <- sub(paste0("^", overlap_prefix), "", names(atac_overlap))
    new_names[1L] <- "atac_peak"  # keep key intact
    # guard against accidental duplicates
    names(atac_overlap) <- make.unique(new_names, sep = "_")
  }

  list(score = dplyr::as_tibble(atac_score),
       overlap = dplyr::as_tibble(atac_overlap))
}

#' Clean HGNC symbols (drop any whitespace suffix), prints a compact diff.
#' @param df data.frame with column HGNC
#' @param label tag used in printed message
#' @return data.frame with cleaned HGNC column
#' @export
clean_hgnc <- function(df, label = "rna") {
  df <- dplyr::mutate(df, HGNC = as.character(.data$HGNC))
  after   <- ifelse(is.na(df$HGNC), NA_character_,
                    stringr::str_replace(df$HGNC, "\\s+.*$", ""))
  changed <- which(!is.na(df$HGNC) & df$HGNC != after)
  if (length(changed)) {
    cat("=== HGNC cleanup:", label, "===\n")
    cat(paste(unique(paste0(df$HGNC[changed], " -> ", after[changed]))), sep = "\n")
  }
  df$HGNC <- after
  df
}


#' Filter RNA-seq by TF/non-TF thresholds with "at least N samples" logic.
#'
#' Keeps rows where >= `min_samples` expression columns meet the threshold:
#'   - non-TFs use `gene_min`
#'   - TFs     use `tf_min`
#'
#' @param rna data.frame/tibble: genes x samples (numeric sample columns).
#' @param tf_list character: TF symbols.
#' @param hgnc_col character: gene symbol column name (default "HGNC").
#' @param gene_min numeric: non-TF threshold (default 5).
#' @param tf_min numeric: TF threshold (default 10).
#' @param min_samples integer: require at least this many samples pass (default 1).
#' @param id_cols character: non-expression columns to exclude when testing.
#'   Defaults to c("ensembl_gene_id", `hgnc_col`).
#' @return tibble with filtered rows (original columns preserved).
#' @export
filter_rna_expr <- function(
    rna,
    tf_list,
    hgnc_col   = "HGNC",
    gene_min   = 5,
    tf_min     = 10,
    min_samples = 1L,
    id_cols    = c("ensembl_gene_id", hgnc_col)
) {
  if (!is.data.frame(rna)) .log_abort("`rna` must be a data.frame/tibble.")
  if (!is.character(tf_list)) .log_abort("`tf_list` must be a character vector.")
  if (!hgnc_col %in% names(rna)) .log_abort("Column {.val {hgnc_col}} not found in `rna`.")
  if (!is.numeric(min_samples) || length(min_samples) != 1 || min_samples < 1)
    .log_abort("`min_samples` must be a positive integer.")

  id_cols   <- unique(id_cols)
  expr_cols <- setdiff(names(rna), id_cols)
  num_cols  <- expr_cols[vapply(rna[expr_cols], is.numeric, logical(1))]
  if (length(num_cols) == 0) .log_abort("No numeric expression columns found after excluding ID columns.")

  # Build logical pass matrices and per-row counts
  M <- as.data.frame(rna[num_cols], stringsAsFactors = FALSE)
  pass_gene <- as.data.frame(lapply(M, function(x) as.numeric(x) >= gene_min), check.names = FALSE)
  pass_tf   <- as.data.frame(lapply(M, function(x) as.numeric(x) >= tf_min),   check.names = FALSE)
  cnt_gene  <- rowSums(pass_gene, na.rm = TRUE)
  cnt_tf    <- rowSums(pass_tf,   na.rm = TRUE)

  is_tf <- rna[[hgnc_col]] %in% tf_list
  keep  <- (!is_tf & cnt_gene >= min_samples) | (is_tf & cnt_tf >= min_samples)

  dplyr::as_tibble(rna[keep, , drop = FALSE])
}


# Correlate footprint score vs TF RNA (per motif)
annotate_fp_tf_corr_one_motif <- function(grn_set,
                                          cor_method = c("pearson","spearman"),
                                          min_non_na = 5L,
                                          tfs_col    = "tfs",
                                          digits     = 6,
                                          dedup_fun  = c("mean","median","first")) {
  cor_method <- match.arg(cor_method)
  dedup_fun  <- match.arg(dedup_fun)

  stopifnot(
    is.list(grn_set),
    all(c("fp_score","rna","fp_tfs","fp_annotation") %in% names(grn_set)),
    "peak_ID" %in% names(grn_set$fp_score),
    "HGNC"    %in% names(grn_set$rna),
    "fp_peak" %in% names(grn_set$fp_tfs),
    tfs_col %in% names(grn_set$fp_tfs) || (tfs_col == "tfs")
  )

  fp_tbl  <- grn_set$fp_score
  rna_tbl <- grn_set$rna
  ft_raw  <- grn_set$fp_tfs

  # collapse to one row per fp_peak with comma-separated TFs
  if (tfs_col == "tfs") {
    ft2 <- ft_raw[, c("fp_peak","tfs"), drop = FALSE]
  } else {
    ft2 <- ft_raw |>
      dplyr::filter(!is.na(.data[[tfs_col]]), .data[[tfs_col]] != "") |>
      dplyr::group_by(fp_peak) |>
      dplyr::summarise(tfs = paste(unique(.data[[tfs_col]]), collapse = ","), .groups = "drop")
  }

  # shared sample columns
  fp_ids  <- setdiff(names(fp_tbl), "peak_ID")
  rna_ids <- setdiff(names(rna_tbl), c("ensembl_gene_id","HGNC"))
  sample_cols <- intersect(fp_ids, rna_ids)
  if (length(sample_cols) < 2L) stop("Need >= 2 shared samples between fp_score and RNA.")

  # aggregate duplicate peak_IDs if needed
  X0 <- fp_tbl[, c("peak_ID", sample_cols), drop = FALSE]
  if (anyDuplicated(X0$peak_ID)) {
    aggfun <- switch(
      dedup_fun,
      mean   = function(z){ m <- mean(z, na.rm = TRUE); if (is.nan(m)) NA_real_ else m },
      median = function(z){ m <- stats::median(z, na.rm = TRUE); if (is.nan(m)) NA_real_ else m },
      first  = function(z){ i <- which(!is.na(z))[1]; if (length(i)) z[i] else NA_real_ }
    )
    Xuni <- X0 |>
      dplyr::group_by(peak_ID) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(sample_cols), aggfun), .groups = "drop")
  } else {
    Xuni <- X0
  }

  # convert to base df before rownames to avoid tibble warning
  Xuni_df <- as.data.frame(Xuni, check.names = FALSE)
  rownames(Xuni_df) <- Xuni_df$peak_ID
  Xmat <- as.matrix(Xuni_df[, sample_cols, drop = FALSE]); storage.mode(Xmat) <- "double"
  mXi  <- is.finite(Xmat)

  # list unique TF symbols
  tf_unique <- unique(unlist(strsplit(paste(ft2$tfs, collapse = ","), "\\s*,\\s*")))
  tf_unique <- tf_unique[nzchar(tf_unique)]
  if (!length(tf_unique)) {
    grn_set$fp_tfs <- dplyr::mutate(ft2, fp_tf_corr_r = NA_character_, fp_tf_corr_p = NA_character_)
    grn_set$fp_annotation <- grn_set$fp_annotation |>
      dplyr::left_join(dplyr::rename(grn_set$fp_tfs,
                                     corr_fp_tf_r = fp_tf_corr_r,
                                     corr_fp_tf_p = fp_tf_corr_p),
                       by = "fp_peak")
    grn_set$fp_tf_corr_diag <- tibble::tibble()
    return(grn_set)
  }

  # RNA matrix for those TFs
  Rdf <- as.data.frame(rna_tbl[rna_tbl$HGNC %in% tf_unique, c("HGNC", sample_cols), drop = FALSE],
                       check.names = FALSE)
  Rdf <- Rdf[!duplicated(Rdf$HGNC), , drop = FALSE]
  rownames(Rdf) <- Rdf$HGNC
  Rmat <- as.matrix(Rdf[, sample_cols, drop = FALSE]); storage.mode(Rmat) <- "double"

  # diagnostics
  diag_df <- tibble::tibble(
    HGNC      = tf_unique,
    present   = tf_unique %in% rownames(Rdf),
    n_finite  = NA_integer_,
    sd_expr   = NA_real_,
    n_ok_peaks= NA_integer_
  )

  # correlation maps
  r_map <- p_map <- setNames(vector("list", length(tf_unique)), tf_unique)

  for (i in seq_along(tf_unique)) {
    tf <- tf_unique[i]
    if (!tf %in% rownames(Rdf)) {
      r_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      p_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      diag_df$n_finite[i]  <- 0L
      diag_df$sd_expr[i]   <- NA_real_
      diag_df$n_ok_peaks[i] <- 0L
      next
    }
    y  <- Rmat[tf, , drop = TRUE]
    my <- is.finite(y)
    diag_df$n_finite[i] <- sum(my)
    diag_df$sd_expr[i]  <- stats::sd(y[my])

    if (sum(my) < 2L || isTRUE(diag_df$sd_expr[i] == 0)) {
      r_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      p_map[[tf]] <- setNames(rep(NA_real_, nrow(Xmat)), rownames(Xmat))
      diag_df$n_ok_peaks[i] <- 0L
      next
    }

    r_v  <- suppressWarnings(as.vector(stats::cor(t(Xmat), y,
                                                  use = "pairwise.complete.obs",
                                                  method = cor_method)))
    n_vec <- as.integer(mXi %*% as.integer(my))
    ok    <- is.finite(r_v) & (n_vec >= min_non_na) & (n_vec > 2L)

    t_stat <- rep(NA_real_, length(r_v))
    t_stat[ok] <- r_v[ok] * sqrt((n_vec[ok]-2) / pmax(1e-12, 1 - r_v[ok]^2))
    p_v <- rep(NA_real_, length(r_v))
    p_v[ok] <- 2 * stats::pt(abs(t_stat[ok]), df = n_vec[ok]-2, lower.tail = FALSE)

    names(r_v) <- rownames(Xmat)
    names(p_v) <- rownames(Xmat)
    r_map[[tf]] <- r_v
    p_map[[tf]] <- p_v
    diag_df$n_ok_peaks[i] <- sum(ok, na.rm = TRUE)
  }

  # format helper
  fmt <- function(v) {
    out <- sprintf(paste0("%.", digits, "g"), v)
    out[is.na(v)] <- NA_character_
    out
  }

  # build comma-separated strings per fp_peak
  assemble <- function(peak, tf_str, which = c("r","p")) {
    which <- match.arg(which)
    if (is.na(tf_str) || !nzchar(tf_str)) return(NA_character_)
    tfs <- strsplit(tf_str, "\\s*,\\s*")[[1]]
    tfs <- tfs[nzchar(tfs)]
    if (!length(tfs)) return(NA_character_)
    vals <- vapply(
      tfs,
      function(tf) if (which == "r") r_map[[tf]][peak] else p_map[[tf]][peak],
      numeric(1)
    )
    paste(fmt(vals), collapse = ",")
  }

  corr_r <- mapply(assemble, ft2$fp_peak, ft2$tfs, MoreArgs = list(which = "r"), USE.NAMES = FALSE)
  corr_p <- mapply(assemble, ft2$fp_peak, ft2$tfs, MoreArgs = list(which = "p"), USE.NAMES = FALSE)

  # write back
  grn_set$fp_tfs <- ft2 |>
    dplyr::mutate(fp_tf_corr_r = as.character(corr_r),
                  fp_tf_corr_p = as.character(corr_p))

  grn_set$fp_annotation <- grn_set$fp_annotation |>
    dplyr::left_join(dplyr::rename(grn_set$fp_tfs,
                                   corr_fp_tf_r = fp_tf_corr_r,
                                   corr_fp_tf_p = fp_tf_corr_p),
                     by = "fp_peak")

  grn_set$fp_tf_corr_diag <- diag_df
  grn_set
}

# Helper: build fp_tfs if missing
.derive_fp_tfs <- function(motif_name, fp_annot, build_args) {
  # Expect fp_annot with fp_peak + motifs; build_args$motif_db has motif + gene_symbol
  gene_col <- .motif_gene_col(build_args$motif_db)
  if (!is.null(build_args$motif_db) &&
      is.data.frame(build_args$motif_db) &&
      !is.null(gene_col) &&
      all(c("motif", gene_col) %in% names(build_args$motif_db)) &&
      "motifs" %in% names(fp_annot)) {
    ft <- fp_annot |>
      dplyr::distinct(fp_peak, motifs) |>
      dplyr::left_join(build_args$motif_db |>
                         dplyr::select(motif, gene_symbol = dplyr::all_of(gene_col)),
                       by = c("motifs" = "motif")) |>
      tidyr::separate_rows(gene_symbol, sep = "::", convert = FALSE) |>
      dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
      dplyr::group_by(fp_peak) |>
      dplyr::summarise(tfs = paste(unique(gene_symbol), collapse = ","), .groups = "drop")
    if (nrow(ft)) return(ft)
  }

  # Fallback: guess TF symbol(s) from motif prefix
  guess <- toupper(sub("_.*$", "", motif_name))
  tibble::tibble(fp_peak = fp_annot$fp_peak,
                 tfs     = guess)
}

# Process ONE motif (load, annotate, write)
process_one_motif_corr <- function(entry,
                                   out_dir,
                                   build_args,
                                   cor_method    = "pearson",
                                   min_non_na    = 5L,
                                   skip_existing = TRUE,
                                   verbose       = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  stopifnot(all(need_cols %in% names(entry)))

  m       <- as.character(entry$motif)[1]
  f_score <- as.character(entry$score)[1]
  f_bound <- as.character(entry$bound)[1]
  f_annot <- as.character(entry$annot)[1]

  # Early graceful skip if any path is NA/missing
  if (is.na(f_score) || is.na(f_bound) || is.na(f_annot) ||
      !file.exists(f_score) || !file.exists(f_bound) || !file.exists(f_annot)) {
    if (verbose) message("Skipping motif ", m, " (missing input file paths).")
    return(tibble::tibble(motif = m,
                          n_rows = 0L,
                          score  = NA_character_,
                          bound  = NA_character_,
                          annot  = NA_character_))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_score <- file.path(out_dir, paste0(m, "_score.csv"))
  out_bound <- file.path(out_dir, paste0(m, "_bound.csv"))
  # NOTE: keep the annot output suffix consistent with your earlier pipeline.
  # You previously used "_annotation.csv" when building the filtered set.
  # If you want consistency, prefer "_annotation.csv" here too:
  out_annot <- file.path(out_dir, paste0(m, "_annotation.csv"))

  if (skip_existing && file.exists(out_score) && file.exists(out_bound) && file.exists(out_annot)) {
    if (verbose) message("Skipping ", m, " (already exists).")
    return(tibble::tibble(motif = m, n_rows = NA_integer_,
                          score = out_score, bound = out_bound, annot = out_annot))
  }

  if (verbose) message("Correlating motif ", m, " ...")

  fp_score <- readr::read_csv(f_score, show_col_types = FALSE)
  fp_bound <- readr::read_csv(f_bound, show_col_types = FALSE)
  fp_annot <- readr::read_csv(f_annot, show_col_types = FALSE)

  # Normalize common column names
  if (!"peak_ID" %in% names(fp_score) && "fp_peak" %in% names(fp_score)) {
    fp_score <- dplyr::rename(fp_score, peak_ID = "fp_peak")
  }
  if (!"peak_ID" %in% names(fp_bound) && "fp_peak" %in% names(fp_bound)) {
    fp_bound <- dplyr::rename(fp_bound, peak_ID = "fp_peak")
  }
  if (!all(c("fp_peak","atac_peak","motifs") %in% names(fp_annot))) {
    fp_annot <- fp_annot |>
      dplyr::rename_with(~"fp_peak",   dplyr::any_of("peak_ID")) |>
      dplyr::rename_with(~"atac_peak", dplyr::any_of("peak_ATAC")) |>
      dplyr::rename_with(~"motifs",    dplyr::any_of("TFBS_name"))
  }

  ss <- do.call(build_grn_set, c(list(
    fp_score      = fp_score,
    fp_bound      = fp_bound,
    fp_annotation = fp_annot
  ), build_args))

  if (!("rna" %in% names(ss))) ss$rna <- build_args$rna

  if (!("fp_tfs" %in% names(ss)) || !NROW(ss$fp_tfs)) {
    ss$fp_tfs <- .derive_fp_tfs(motif_name = m, fp_annot = ss$fp_annotation, build_args = build_args)
  } else if (!"tfs" %in% names(ss$fp_tfs)) {
    if ("HGNC" %in% names(ss$fp_tfs)) {
      ss$fp_tfs <- ss$fp_tfs |>
        dplyr::filter(!is.na(HGNC), HGNC != "") |>
        dplyr::group_by(fp_peak) |>
        dplyr::summarise(tfs = paste(unique(HGNC), collapse = ","), .groups = "drop")
    } else {
      ss$fp_tfs <- .derive_fp_tfs(motif_name = m, fp_annot = ss$fp_annotation, build_args = build_args)
    }
  }

  have_all <- all(c("fp_score","fp_annotation","rna","fp_tfs") %in% names(ss)) &&
    is.data.frame(ss$fp_score) && is.data.frame(ss$fp_annotation) &&
    NROW(ss$rna) > 0L && NROW(ss$fp_tfs) > 0L
  if (have_all) {
    ss <- annotate_fp_tf_corr_one_motif(ss,
                                        cor_method = cor_method,
                                        min_non_na = min_non_na,
                                        tfs_col    = "tfs")
  } else if (verbose) {
    message("Skipping correlation for ", m, ": missing one of fp_score/fp_annotation/rna/fp_tfs.")
  }

  readr::write_csv(ss$fp_score,      out_score)
  readr::write_csv(ss$fp_bound,      out_bound)
  readr::write_csv(ss$fp_annotation, out_annot)

  tibble::tibble(motif = m,
                 n_rows = nrow(ss$fp_bound),
                 score  = out_score,
                 bound  = out_bound,
                 annot  = out_annot)
}


#' Parallel driver to correlate TF footprints with expression (manifest-aware)
#'
#' Runs motif-by-motif correlation and returns a manifest. When `skip_existing=TRUE`
#' and existing result CSVs are detected for a motif, this function now *backfills*
#' the `n_rows` column by counting rows from the existing `score` CSV so you don't
#' end up with `NA` in the output manifest.
#'
#' @param fp_manifest data.frame/tibble with columns: motif, score, bound, annot
#'        (paths to per-motif CSVs or inputs needed by `process_one_motif_corr`)
#' @param out_dir output directory for new results
#' @param build_args list forwarded to `process_one_motif_corr()`
#' @param motif_ids optional character vector to subset motifs
#' @param n_workers integer number of workers (default: detectCores)
#' @param set_plan logical; if TRUE set a future plan inside the function
#' @param skip_existing logical; if TRUE, skip recomputation when outputs exist
#' @param cor_method character; correlation method, e.g. "pearson"
#' @param min_non_na integer; minimum non-NA pairs for correlation
#' @param verbose logical; emit progress messages
#'
#' @return tibble with columns including motif, n_rows, score, bound, annot
#' @export
tf_corr_footprints <- function(fp_manifest,
                               out_dir,
                               build_args,
                               motif_ids     = NULL,
                               n_workers     = max(1L, parallel::detectCores(TRUE)),
                               set_plan      = TRUE,
                               skip_existing = TRUE,
                               cor_method    = "pearson",
                               min_non_na    = 5L,
                               verbose       = TRUE) {
  need_cols <- c("motif","score","bound","annot")
  stopifnot(is.data.frame(fp_manifest), all(need_cols %in% names(fp_manifest)))

  rows <- fp_manifest
  if (!is.null(motif_ids)) {
    rows <- dplyr::semi_join(rows, tibble::tibble(motif = motif_ids), by = "motif")
    if (!nrow(rows)) stop("No matching motifs found in `fp_manifest` for the requested subset.")
  }

  # Keep only rows with non-NA paths AND files that exist
  rows <- rows |>
    dplyr::filter(!is.na(score), !is.na(bound), !is.na(annot)) |>
    dplyr::mutate(
      .ok = file.exists(score) & file.exists(bound) & file.exists(annot)
    ) |>
    dplyr::filter(.ok) |>
    dplyr::select(-.ok)

  if (!nrow(rows)) stop("After filtering, no motifs have available input files.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (verbose) message("Submitting ", nrow(rows), " motif(s) to ", n_workers, " worker(s) for correlation...")

  # Choose a plan if requested
  if (set_plan) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    if (.Platform$OS.type == "unix" &&
        n_workers > 1L &&
        !identical(Sys.getenv("RSTUDIO"), "1") &&
        requireNamespace("parallelly", quietly = TRUE) &&
        parallelly::supportsMulticore()) {
      future::plan(future::multicore, workers = n_workers)
    } else {
      future::plan(future::multisession, workers = n_workers)
    }
  }

  # Helper: fast row count for CSV (header counted by fread, we want nrow)
  fast_count_csv <- function(path) {
    # Prefer data.table::fread if available; otherwise fallback to readr
    if (requireNamespace("data.table", quietly = TRUE)) {
      dt <- tryCatch(
        data.table::fread(path, nThread = 1L, showProgress = FALSE),
        error = function(e) NULL
      )
      return(if (is.null(dt)) NA_integer_ else nrow(dt))
    } else if (requireNamespace("readr", quietly = TRUE)) {
      df <- tryCatch(
        readr::read_csv(path, progress = FALSE, show_col_types = FALSE),
        error = function(e) NULL
      )
      return(if (is.null(df)) NA_integer_ else nrow(df))
    } else {
      # Very small fallback using base; may be slower on large files
      con <- NULL
      out <- NA_integer_
      try({
        con <- file(path, open = "r")
        # Read & discard header
        header <- readLines(con, n = 1L, warn = FALSE)
        cnt <- 0L
        chunk <- character()
        repeat {
          chunk <- readLines(con, n = 100000L, warn = FALSE)
          if (!length(chunk)) break
          cnt <- cnt + length(chunk)
        }
        out <- cnt
      }, silent = TRUE)
      if (!is.null(con)) close(con)
      return(out)
    }
  }

  # Run per motif
  res <- future.apply::future_lapply(seq_len(nrow(rows)), function(i) {
    entry <- rows[i, , drop = FALSE]
    process_one_motif_corr(entry = entry,
                           out_dir = out_dir,
                           build_args = build_args,
                           cor_method = cor_method,
                           min_non_na = min_non_na,
                           skip_existing = skip_existing,
                           verbose = FALSE)
  }, future.seed = FALSE)

  out <- tibble::as_tibble(dplyr::bind_rows(res))

  # Backfill n_rows when skip_existing caused a bypass and left NA
  # Criteria: n_rows is NA, score path is non-NA and exists
  # We count rows from the score CSV as authoritative per-motif row count.
  if ("n_rows" %in% names(out)) {
    need_fill <- which(is.na(out$n_rows) & !is.na(out$score) & file.exists(out$score))
    if (length(need_fill)) {
      if (verbose) message("Backfilling n_rows for ", length(need_fill), " motif(s) from existing score CSVs...")
      filled <- vapply(need_fill, function(idx) fast_count_csv(out$score[[idx]]), integer(1))
      out$n_rows[need_fill] <- filled
    }
  } else {
    # If the downstream returns no n_rows column, create and fill it where possible
    out$n_rows <- NA_integer_
    need_fill <- which(!is.na(out$score) & file.exists(out$score))
    if (length(need_fill)) {
      if (verbose) message("Populating n_rows from existing score CSVs for ", length(need_fill), " motif(s)...")
      filled <- vapply(need_fill, function(idx) fast_count_csv(out$score[[idx]]), integer(1))
      out$n_rows[need_fill] <- filled
    }
    # Ensure conventional column order: motif, n_rows, then the paths
    front_cols <- c("motif", "n_rows")
    rest_cols  <- setdiff(names(out), front_cols)
    out <- out[, c(front_cols, rest_cols), drop = FALSE]
  }

  out
}


#' Correlate one TF with ALL footprints and export per-TF BEDs
#'
#' Minimal, fast path to:
#' 1) load a single TF expression vector,
#' 2) load *all* footprint score/annotation files from a manifest,
#' 3) compute row-wise correlations (per TFBS) vs TF expression,
#' 4) BH-adjust p-values,
#' 5) write 3 BED files (`*_all.bed`, `*_bound.bed`, `*_unbound.bed`) and one
#'    tab-delimited overview (`*_overview.txt`) for the TF.
#'
#' Assumptions:
#' - `fp_manifest` has columns `motif, n_rows, score, bound, annot`.
#' - Each `score` CSV has a `peak_ID` (or `fp_peak`) column and per-sample columns.
#' - Each `annot` CSV has `fp_peak` (or `peak_ID`), `atac_peak` (or `peak_ATAC`),
#'   and `motifs` (or `TFBS_name`) columns linking TFBSΓåÆATAC peak and name.
#'
#' @param fp_manifest data.frame with columns motif,n_rows,score,bound,annot.
#' @param tf_name     Character scalar; label written to outputs (e.g., "HNF4A").
#' @param tf_expr     Named numeric vector (preferred) of TF expression with names
#'                    matching the footprint score column names. Alternatively, a
#'                    one-row numeric data.frame with column names as sample IDs.
#' @param out_dir     Output directory to write files.
#' @param cor_method  "pearson" (default) or "spearman".
#' @param min_non_na  Minimum paired non-NA values required to test (default 5).
#' @param r_thr       Correlation threshold for "bound" (default 0.30).
#' @param fdr_thr     BH-FDR threshold for "bound" (default 0.05).
#' @param threads     Integer threads for data.table fread (default: detectCores()).
#' @param verbose     Logical; progress messages (default TRUE).
#'
#' @return (Invisibly) a list with:
#' \describe{
#'   \item{annot}{Tibble of TFBS annotations with r, p, p_adj.}
#'   \item{files}{Named character vector of written file paths.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' tf_vec <- c(S1=10, S2=8, S3=12)  # names must match columns in score CSVs
#' res <- tf_corr_footprints_all_tfbs(
#'   fp_manifest = fp_aligned_normalized_filtered_manifest,
#'   tf_name  = "HNF4A",
#'   tf_expr  = tf_vec,
#'   out_dir  = file.path(tempdir(), "hnf4a_out")
#' )
#' res$files
#' }
tf_corr_footprints_all_tfbs <- function(
    fp_manifest,
    tf_name,
    tf_expr,
    out_dir,
    cor_method  = "pearson",
    min_non_na  = 5L,
    r_thr       = 0.30,
    fdr_thr     = 0.05,
    threads     = max(1L, parallel::detectCores(TRUE)),
    verbose     = TRUE
) {
  # ---- checks ----
  need_cols <- c("motif","n_rows","score","bound","annot")
  if (!is.data.frame(fp_manifest) || !all(need_cols %in% names(fp_manifest))) {
    .log_abort("`fp_manifest` must be a data.frame with columns: {need_cols}.")
  }
  if (!is.character(tf_name) || length(tf_name) != 1L || !nzchar(tf_name)) {
    .log_abort("`tf_name` must be a non-empty character scalar.")
  }
  # Coerce tf_expr
  if (is.data.frame(tf_expr)) {
    if (nrow(tf_expr) != 1L) .log_abort("`tf_expr` data.frame must have exactly one row.")
    tf_expr <- as.numeric(tf_expr[1, , drop = TRUE])
    names(tf_expr) <- colnames(tf_expr)
  }
  if (!is.numeric(tf_expr) || is.null(names(tf_expr)) || !length(tf_expr)) {
    .log_abort("`tf_expr` must be a named numeric vector (names = sample IDs).")
  }

  use <- fp_manifest[
    is.finite(fp_manifest$n_rows) & fp_manifest$n_rows > 0 &
      !is.na(fp_manifest$score) & file.exists(fp_manifest$score) &
      !is.na(fp_manifest$annot) & file.exists(fp_manifest$annot),
    , drop = FALSE
  ]
  if (!nrow(use)) {
    .log_abort("No non-empty motifs with existing {`score`} and {`annot`} files.")
  }

  # ---- IO helpers ----
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fread_fast <- function(path) data.table::fread(path, nThread = as.integer(threads), showProgress = FALSE)

  # threads scope
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(threads))

  # ---- read & normalize annotations (fp_peak, atac_peak, motifs) ----
  norm_annot <- function(dt) {
    if (!"fp_peak"   %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID", "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"    %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    keep <- intersect(c("fp_peak","atac_peak","motifs"), names(dt))
    dt <- dt[, ..keep]
    # drop rows without fp_peak
    dt <- dt[!is.na(fp_peak) & nzchar(fp_peak)]
    dt
  }
  annot_dt <- data.table::rbindlist(lapply(use$annot, function(p) norm_annot(fread_fast(p))),
                                    use.names = TRUE, fill = TRUE)
  if (!nrow(annot_dt) || !"fp_peak" %in% names(annot_dt)) {
    .log_abort("Annotation inputs lack required `fp_peak` column after normalization.")
  }
  # Deduplicate TFBS rows by first occurrence
  annot_dt <- annot_dt[!duplicated(fp_peak)]

  # ---- read & normalize scores (peak_ID + sample columns) ----
  norm_score <- function(dt) {
    if (!"peak_ID" %in% names(dt) && "fp_peak" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    if (!"peak_ID" %in% names(dt)) .log_abort("Score file is missing `peak_ID` (or `fp_peak`) column.")
    # keep only peak_ID + intersecting sample columns (matching tf_expr names)
    sample_cols <- intersect(setdiff(names(dt), "peak_ID"), names(tf_expr))
    if (!length(sample_cols)) .log_abort("No overlapping sample columns between score file and `tf_expr` names.")
    # coerce sample columns to numeric
    for (j in sample_cols) data.table::set(dt, j = j, value = suppressWarnings(as.numeric(dt[[j]])))
    dt[, c("peak_ID", sample_cols), with = FALSE]
  }
  score_dt <- data.table::rbindlist(lapply(use$score, function(p) norm_score(fread_fast(p))),
                                    use.names = TRUE, fill = TRUE)
  # de-dup by first occurrence and keep only peaks present in annotations
  score_dt <- score_dt[!duplicated(peak_ID)]
  score_dt <- score_dt[peak_ID %chin% annot_dt$fp_peak]

  if (!nrow(score_dt)) {
    .log_abort("After alignment, no overlapping TFBS between scores and annotations.")
  }

  # align sample columns to tf_expr names & compute correlation per row
  sample_cols <- intersect(setdiff(names(score_dt), "peak_ID"), names(tf_expr))
  if (length(sample_cols) < min_non_na) {
    .log_abort("Not enough overlapping samples (need >= {min_non_na}, got {length(sample_cols)}).")
  }
  # order tf_expr to score columns
  y <- suppressWarnings(as.numeric(tf_expr[sample_cols]))
  names(y) <- sample_cols

  X <- as.matrix(score_dt[, ..sample_cols])

  row_stat <- function(v) {
    ok <- is.finite(v) & is.finite(y)
    n  <- sum(ok)
    if (n < min_non_na) return(c(NA_real_, NA_real_, as.integer(n)))
    r <- suppressWarnings(stats::cor(v[ok], y[ok], method = cor_method))
    if (!is.finite(r)) return(c(NA_real_, NA_real_, as.integer(n)))
    tval <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))
    pval <- 2 * stats::pt(-abs(tval), df = n - 2)
    c(r, pval, as.integer(n))
  }
  mat <- t(apply(X, 1L, row_stat))
  score_dt[, `:=`(
    corr_fp_tf_r = mat[, 1],
    corr_fp_tf_p = mat[, 2],
    pairs        = as.integer(mat[, 3])
  )]

  # BH across all finite p
  good_p <- is.finite(score_dt$corr_fp_tf_p)
  score_dt[, corr_fp_tf_p_adj := NA_real_]
  if (any(good_p)) {
    score_dt[good_p, corr_fp_tf_p_adj := stats::p.adjust(score_dt$corr_fp_tf_p[good_p], method = "BH")]
    score_dt[good_p, corr_fp_tf_p_adj := pmax(corr_fp_tf_p_adj, corr_fp_tf_p)]
  }

  # ---- join annotation + stats & build export tables ----
  merged <- score_dt[, .(fp_peak = peak_ID, corr_fp_tf_r, corr_fp_tf_p, corr_fp_tf_p_adj)][
    annot_dt, on = "fp_peak"
  ]
  merged[, TF := tf_name]

  # split coordinates
  mk_table <- function(dt) {
    if (!nrow(dt)) return(data.table::data.table())
    fp_s <- data.table::tstrsplit(dt$fp_peak, "[:-]", perl = TRUE)
    at_s <- data.table::tstrsplit(dt$atac_peak, "[:-]", perl = TRUE)
    data.table::data.table(
      TFBS_chr         = fp_s[[1]],
      TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
      TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
      TFBS_name        = dt$motifs,
      peak_chr         = at_s[[1]],
      peak_start       = suppressWarnings(as.integer(at_s[[2]])),
      peak_end         = suppressWarnings(as.integer(at_s[[3]])),
      TF               = dt$TF,
      corr_fp_tf_r     = dt$corr_fp_tf_r,
      corr_fp_tf_p     = dt$corr_fp_tf_p,
      corr_fp_tf_p_adj = dt$corr_fp_tf_p_adj
    )
  }
  all_dt <- mk_table(merged)

  # bound vs unbound
  is_bound <- is.finite(all_dt$corr_fp_tf_r) &
    is.finite(all_dt$corr_fp_tf_p_adj) &
    (all_dt$corr_fp_tf_r > r_thr) &
    (all_dt$corr_fp_tf_p_adj < fdr_thr)
  bound_dt   <- all_dt[is_bound]
  unbound_dt <- all_dt[!is_bound]

  # ---- write outputs ----
  .slug <- function(x) {
    out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    out <- gsub("^_+|_+$", "", out)
    if (!nzchar(out)) out <- "TF"
    out
  }
  tf_slug <- .slug(tf_name)
  f_all      <- file.path(out_dir, paste0(tf_slug, "_all.bed"))
  f_bound    <- file.path(out_dir, paste0(tf_slug, "_bound.bed"))
  f_unbound  <- file.path(out_dir, paste0(tf_slug, "_unbound.bed"))
  f_overview <- file.path(out_dir, paste0(tf_slug, "_overview.txt"))

  # BEDs: no header
  data.table::fwrite(all_dt,     f_all,     sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")
  data.table::fwrite(bound_dt,   f_bound,   sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")
  data.table::fwrite(unbound_dt, f_unbound, sep = "\t", col.names = FALSE, quote = FALSE, na = "NA")

  # overview: header + _bound flag
  if (nrow(all_dt)) {
    all_dt[, `_bound` := as.integer(is_bound)]
    data.table::fwrite(all_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE, na = "NA")
  } else {
    hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
             "peak_chr","peak_start","peak_end","TF",
             "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
    empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode="character", length=0)), length(hdr)), hdr))
    data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
  }

  if (verbose) {
    .log_inform(c(
      "TF: {tf_name}",
      "All TFBS: {nrow(all_dt)}; Bound (r>{r_thr}, FDR<{fdr_thr}): {nrow(bound_dt)}; Unbound: {nrow(unbound_dt)}",
      "Wrote:",
      "  - {f_all}",
      "  - {f_bound}",
      "  - {f_unbound}",
      "  - {f_overview}"
    ))
  }

  invisible(list(
    annot = tibble::as_tibble(merged),
    files = c(all = f_all, bound = f_bound, unbound = f_unbound, overview = f_overview)
  ))
}

#' Combine FP correlation tables (split dimers), run BH FDR once across all,
#' sanity-check p vs p_adj, filter, and align score/bound to the kept annotations.
#'
#' Robust, minimal pipeline:
#' 1) Normalizes columns when reading each file.
#' 2) Splits comma-separated TF/r/p columns row-wise (handles dimers) with safe padding.
#' 3) Runs **one** BenjaminiΓÇôHochberg adjustment across the combined annotations.
#' 4) Enforces `p_adj >= p` and reports monotonicity stats.
#' 5) Filters by `r_thr` (positive) and `p_thr` (FDR).
#' 6) Reads score/bound only for peaks that pass, de-dups on `peak_ID`, and aligns order.
#'
#' @param fp_filtered_manifest Tibble/data.frame with columns:
#'   `motif, n_rows, score, bound, annot` (CSV paths).
#' @param p_thr Numeric FDR threshold (default 0.05).
#' @param r_thr Numeric correlation threshold (keep r > r_thr; default 0.3).
#' @param verbose Logical; print progress (default TRUE).
#' @param threads Integer; data.table threads (default = all detected cores).
#' @param output_bed Optional directory path. If provided, the function will write per-TF
#'   BED-like files: `<TF>_all.bed`, `<TF>_bound.bed`, `<TF>_unbound.bed` (no headers),
#'   and `<TF>_overview.txt` (headered; same columns as `_all.bed` plus a `_bound` 0/1 column).
#'   Columns: TFBS_chr, TFBS_start, TFBS_end, TFBS_name, peak_chr, peak_start, peak_end,
#'   TF, corr_fp_tf_r, corr_fp_tf_p, corr_fp_tf_p_adj; bound/unbound defined by r/p filters.
#'
#' @return A list of tibbles: `fp_score`, `fp_bound`, `fp_annotation`
#'         (all aligned to the same row order).
#' @export
#' @examples
#' # res <- tf_corr_footprints_filter(fp_corr_manifest)
tf_corr_footprints_filter <- function(
    fp_filtered_manifest,
    p_thr   = 0.05,
    r_thr   = 0.3,
    verbose = TRUE,
    threads = max(1L, parallel::detectCores(TRUE)),
    # BEGIN EDIT: add optional output directory for BED exports
    output_bed = NULL
    # END EDIT
) {
  # ---- checks ----
  if (!is.data.frame(fp_filtered_manifest) ||
      !all(c("motif","n_rows","score","bound","annot") %in% names(fp_filtered_manifest))) {
    .log_abort("`fp_filtered_manifest` must have columns: motif, n_rows, score, bound, annot.")
  }

  use <- fp_filtered_manifest |>
    dplyr::filter(
      !is.na(.data$n_rows), .data$n_rows > 0,
      !is.na(.data$score), file.exists(.data$score),
      !is.na(.data$bound), file.exists(.data$bound),
      !is.na(.data$annot), file.exists(.data$annot)
    )

  if (!nrow(use)) {
    if (verbose) .log_inform("No non-empty motifs with existing files; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }
  if (verbose) .log_inform("Combining {nrow(use)} motif{?s} ...")

  # ---- threading ----
  old_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_threads), add = TRUE)
  data.table::setDTthreads(as.integer(threads))

  fread_fast <- function(path) data.table::fread(path, nThread = threads, showProgress = FALSE)

  # ---- helpers ----
  split_commas <- function(x) {
    x <- as.character(x)
    lapply(x, function(s) {
      if (is.na(s) || !nzchar(s)) character(0)
      else strsplit(s, "\\s*,\\s*", perl = TRUE)[[1]]
    })
  }
  pad_to <- function(lst, lens, fill) {
    out <- vector("list", length(lst))
    for (i in seq_along(lst)) {
      v <- lst[[i]]
      # recycle singleton
      if (length(v) == 1L && lens[i] > 1L) v <- rep(v, lens[i])
      # right-pad
      if (length(v) < lens[i]) v <- c(v, rep(fill, lens[i] - length(v)))
      out[[i]] <- v
    }
    unlist(out, use.names = FALSE)
  }
  safe_num <- function(x) {
    x <- as.character(x)
    x[!nzchar(x)] <- NA_character_   # treat "" as NA
    suppressWarnings(as.numeric(x))  # quiet "NAs introduced by coercion"
  }

  # ---- read one annotation, normalize names, split dimers ----
  read_annot_split <- function(p) {
    dt <- fread_fast(p)

    # Normalize names to: fp_peak / atac_peak / motifs / tfs / corr_fp_tf_r / corr_fp_tf_p
    if (!"fp_peak"   %in% names(dt) && "peak_ID"   %in% names(dt)) data.table::setnames(dt, "peak_ID",   "fp_peak")
    if (!"atac_peak" %in% names(dt) && "peak_ATAC" %in% names(dt)) data.table::setnames(dt, "peak_ATAC", "atac_peak")
    if (!"motifs"    %in% names(dt) && "TFBS_name" %in% names(dt)) data.table::setnames(dt, "TFBS_name", "motifs")
    if (!"tfs"       %in% names(dt)) dt[["tfs"]] <- NA_character_
    if (!"corr_fp_tf_r" %in% names(dt)) dt[["corr_fp_tf_r"]] <- NA_character_
    if (!"corr_fp_tf_p" %in% names(dt)) dt[["corr_fp_tf_p"]] <- NA_character_

    tfs_list <- split_commas(dt[["tfs"]])
    r_list   <- lapply(split_commas(dt[["corr_fp_tf_r"]]), safe_num)
    p_list   <- lapply(split_commas(dt[["corr_fp_tf_p"]]), safe_num)

    # how many repeats per row
    nrep <- pmax(lengths(tfs_list), lengths(r_list), lengths(p_list))
    nrep[nrep == 0L] <- 1L
    idx <- rep.int(seq_len(nrow(dt)), nrep)

    data.table::data.table(
      fp_peak      = dt$fp_peak[idx],
      atac_peak    = dt$atac_peak[idx],
      motifs       = dt$motifs[idx],
      tfs          = pad_to(tfs_list, nrep, NA_character_),
      corr_fp_tf_r = safe_num(pad_to(r_list, nrep, NA_real_)),
      corr_fp_tf_p = safe_num(pad_to(p_list, nrep, NA_real_))
    )
  }

  # ---- 1) combine all annotations ----
  annot_dt <- data.table::rbindlist(lapply(use$annot, read_annot_split),
                                    use.names = TRUE, fill = TRUE)

  # Optional quick NA summary
  if (verbose) {
    n_p_na <- sum(!is.finite(annot_dt$corr_fp_tf_p))
    n_r_na <- sum(!is.finite(annot_dt$corr_fp_tf_r))
    .log_inform("Pre-filter NA counts: p={n_p_na}, r={n_r_na}")
  }

  # ---- sanity: p in [0,1] (raw) ----
  bad_gt1 <- sum(annot_dt$corr_fp_tf_p > 1, na.rm = TRUE)
  bad_lt0 <- sum(annot_dt$corr_fp_tf_p < 0, na.rm = TRUE)
  if (bad_gt1 || bad_lt0) {
    .log_abort("[raw p] Found {bad_gt1} values > 1 and {bad_lt0} values < 0. Please fix upstream.")
  }

  # ---- pre-BH: dedup (fp_peak, tfs) by keeping the FIRST occurrence ----
  if (verbose) .log_inform("Pre-BH: rows before pair-dedup: {nrow(annot_dt)}")

  keep_idx <- !duplicated(annot_dt[, .(fp_peak, tfs)])
  n_drop   <- sum(!keep_idx)

  annot_dt <- annot_dt[keep_idx]

  if (verbose) {
    .log_inform("Pre-BH: dropped {n_drop} duplicate (fp_peak, tfs) pairs (kept first).")
    .log_inform("Pre-BH: rows after pair-dedup: {nrow(annot_dt)}")
  }

  # ---- BH adjust across ALL rows at once (finite p only) ----
  good_p <- is.finite(annot_dt$corr_fp_tf_p)
  annot_dt[, corr_fp_tf_p_adj := NA_real_]
  if (any(good_p)) {
    annot_dt[good_p, corr_fp_tf_p_adj := stats::p.adjust(annot_dt$corr_fp_tf_p[good_p], method = "BH")]
    # enforce theoretical constraint (defensive against FP/rounding)
    annot_dt[good_p, corr_fp_tf_p_adj := pmax(corr_fp_tf_p_adj, corr_fp_tf_p)]
  }

  # ---- post-BH sanity stats ----
  cmp <- data.table::fcase(
    !is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj), NA_character_,
    annot_dt$corr_fp_tf_p_adj >  annot_dt$corr_fp_tf_p, ">",
    annot_dt$corr_fp_tf_p_adj == annot_dt$corr_fp_tf_p, "==",
    annot_dt$corr_fp_tf_p_adj <  annot_dt$corr_fp_tf_p, "<"
  )
  if (verbose) {
    .log_inform(c(
      "BH sanity (p_adj vs p):",
      "  >  = {sum(cmp == '>', na.rm = TRUE)}",
      "  == = {sum(cmp == '==', na.rm = TRUE)}",
      "  <  = {sum(cmp == '<', na.rm = TRUE)}  (after pmax() this should be 0)",
      "  NA = {sum(!is.finite(annot_dt$corr_fp_tf_p) | !is.finite(annot_dt$corr_fp_tf_p_adj))}",
      "  total rows = {nrow(annot_dt)}"
    ))
  }

  # BEGIN EDIT: snapshot pre-filter annotations for BED exports
  annot_pre_filter <- data.table::copy(annot_dt)
  # END EDIT

  # ---- 2) filter: positive r and FDR ----
  keep <- is.finite(annot_dt$corr_fp_tf_r) &
    is.finite(annot_dt$corr_fp_tf_p_adj) &
    (annot_dt$corr_fp_tf_r > r_thr) &
    (annot_dt$corr_fp_tf_p_adj < p_thr)
  annot_dt <- annot_dt[keep]
  if (verbose) {
    .log_inform("After corr/FDR filter: {nrow(annot_dt)} rows (unique peaks: {length(unique(annot_dt$fp_peak))}).")
  }
  if (nrow(annot_dt) == 0L) {
    # BEGIN EDIT: still optionally create empty per-TF files if requested
    if (!is.null(output_bed)) {
      if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)
      if (verbose) .log_inform("No rows passed; requested BED export directory initialized at {output_bed}.")
    }
    # END EDIT
    if (verbose) .log_inform("Nothing passed thresholds; returning empties.")
    return(list(
      fp_score      = tibble::tibble(),
      fp_bound      = tibble::tibble(),
      fp_annotation = tibble::tibble()
    ))
  }

  peaks_keep <- unique(annot_dt$fp_peak)

  # ---- helpers to normalize score/bound ----
  norm_score <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.numeric(dt[[j]])))
    dt[]
  }
  norm_bound <- function(dt) {
    if ("fp_peak" %in% names(dt) && !"peak_ID" %in% names(dt)) data.table::setnames(dt, "fp_peak", "peak_ID")
    idc <- setdiff(names(dt), "peak_ID")
    for (j in idc) data.table::set(dt, j = j, value = suppressWarnings(as.integer(safe_num(dt[[j]]))))
    dt[]
  }

  read_score_subset <- function(p) {
    dt <- norm_score(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }
  read_bound_subset <- function(p) {
    dt <- norm_bound(fread_fast(p))
    dt[peak_ID %chin% peaks_keep]
  }

  score_dt <- data.table::rbindlist(lapply(use$score, read_score_subset), use.names = TRUE, fill = TRUE)
  bound_dt <- data.table::rbindlist(lapply(use$bound, read_bound_subset), use.names = TRUE, fill = TRUE)

  # de-dup by first occurrence
  if (nrow(score_dt)) score_dt <- score_dt[!duplicated(peak_ID)]
  if (nrow(bound_dt)) bound_dt <- bound_dt[!duplicated(peak_ID)]

  # ---- 3) align score/bound EXACTLY to annotation order ----
  key_vec <- annot_dt$fp_peak
  key_dt  <- data.table::data.table(peak_ID = key_vec)

  fp_score_dt <- score_dt[key_dt, on = "peak_ID"]  # preserves key_dt order
  fp_bound_dt <- bound_dt[key_dt, on  = "peak_ID"]

  # final sanity: 1:1 and same order
  stopifnot(
    nrow(fp_score_dt) == nrow(annot_dt),
    nrow(fp_bound_dt) == nrow(annot_dt),
    identical(fp_score_dt$peak_ID, key_vec),
    identical(fp_bound_dt$peak_ID, key_vec)
  )

  # BEGIN EDIT: optional per-TF BED & overview export
  if (!is.null(output_bed)) {
    if (!dir.exists(output_bed)) dir.create(output_bed, recursive = TRUE, showWarnings = FALSE)

    .slug <- function(x) gsub("[^A-Za-z0-9]+", "_", as.character(x))
    .mk_table <- function(dt) {
      if (!nrow(dt)) return(data.table::data.table())
      dt2 <- dt[!is.na(tfs) & nzchar(tfs) & !is.na(fp_peak) & !is.na(atac_peak)]
      if (!nrow(dt2)) return(data.table::data.table())
      fp_s <- data.table::tstrsplit(dt2$fp_peak, "[:-]", perl = TRUE)
      at_s <- data.table::tstrsplit(dt2$atac_peak, "[:-]", perl = TRUE)
      data.table::data.table(
        TFBS_chr         = fp_s[[1]],
        TFBS_start       = suppressWarnings(as.integer(fp_s[[2]])),
        TFBS_end         = suppressWarnings(as.integer(fp_s[[3]])),
        TFBS_name        = dt2$motifs,
        peak_chr         = at_s[[1]],
        peak_start       = suppressWarnings(as.integer(at_s[[2]])),
        peak_end         = suppressWarnings(as.integer(at_s[[3]])),
        TF               = dt2$tfs,
        corr_fp_tf_r     = dt2$corr_fp_tf_r,
        corr_fp_tf_p     = dt2$corr_fp_tf_p,
        corr_fp_tf_p_adj = dt2$corr_fp_tf_p_adj
      )
    }

    all_dt    <- .mk_table(annot_pre_filter)
    bound_dt2 <- .mk_table(annot_pre_filter[keep])

    if (nrow(all_dt))   all_dt[,   .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]
    if (nrow(bound_dt2)) bound_dt2[, .id := paste(TF, TFBS_chr, TFBS_start, TFBS_end, sep = ";")]

    tfs_vec <- sort(unique(all_dt$TF))
    for (tf in tfs_vec) {
      tf_slug <- .slug(tf)
      f_all      <- file.path(output_bed, paste0(tf_slug, "_all.bed"))
      f_bound    <- file.path(output_bed, paste0(tf_slug, "_bound.bed"))
      f_unbound  <- file.path(output_bed, paste0(tf_slug, "_unbound.bed"))
      f_overview <- file.path(output_bed, paste0(tf_slug, "_overview.txt"))

      sub_all   <- all_dt[TF == tf]
      sub_bound <- bound_dt2[TF == tf]

      # --- FIX 1: use %in% (not %chin%) ---
      ids_bound   <- if (nrow(sub_bound)) sub_bound$.id else character(0)
      sub_unbound <- sub_all[!(.id %in% ids_bound)]

      # --- FIX 2: drop ".id" (not "._id"); keep=FALSE selection with data.table ---
      data.table::fwrite(sub_all[,   !".id", with = FALSE], f_all,     sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_bound[, !".id", with = FALSE], f_bound,   sep = "\t", col.names = FALSE, quote = FALSE)
      data.table::fwrite(sub_unbound[, !".id", with = FALSE], f_unbound, sep = "\t", col.names = FALSE, quote = FALSE)

      # headered overview with _bound flag
      if (nrow(sub_all)) {
        sub_all[, `_bound` := as.integer(.id %in% ids_bound)]   # <-- FIX 1 applied here too
        data.table::fwrite(sub_all[, !".id", with = FALSE], f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      } else {
        hdr <- c("TFBS_chr","TFBS_start","TFBS_end","TFBS_name",
                 "peak_chr","peak_start","peak_end","TF",
                 "corr_fp_tf_r","corr_fp_tf_p","corr_fp_tf_p_adj","_bound")
        # create a 0-row table with the right header
        empty_dt <- data.table::as.data.table(setNames(rep(list(vector(mode="character", length=0)), length(hdr)), hdr))
        data.table::fwrite(empty_dt, f_overview, sep = "\t", col.names = TRUE, quote = FALSE)
      }
    }
    if (verbose) .log_inform("Per-TF BED exports written to: {output_bed}")
  }
  # END EDIT


  # ---- return ----
  list(
    fp_score      = tibble::as_tibble(fp_score_dt),
    fp_bound      = tibble::as_tibble(fp_bound_dt),
    fp_annotation = tibble::as_tibble(annot_dt)
  )
}

