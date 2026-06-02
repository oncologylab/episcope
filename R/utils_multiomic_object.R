# File: utils_multiomic_object.R
# Purpose: CraftGRN multiomic object helpers for Module 1.

.multiomic_object_schema <- function() "craftgrn_multiomic_v1"

#' @noRd
is_multiomic_object <- function(x) {
  is.list(x) &&
    identical(x$schema, .multiomic_object_schema()) &&
    is.list(x$features) &&
    is.list(x$matrices)
}


.compact_condition_columns <- function(x, id_cols) {
  setdiff(names(x), id_cols)
}

.compact_matrix_from_table <- function(x, id_col, storage = c("numeric", "logical")) {
  storage <- match.arg(storage)
  if (!is.data.frame(x)) .log_abort("Expected a data.frame for matrix conversion.")
  if (!id_col %in% names(x)) .log_abort("Missing required id column: {id_col}")
  ids <- as.character(x[[id_col]])
  keep <- !duplicated(ids)
  x <- x[keep, , drop = FALSE]
  ids <- ids[keep]
  cols <- setdiff(names(x), id_col)
  mat <- data.matrix(x[, cols, drop = FALSE])
  rownames(mat) <- ids
  if (identical(storage, "logical")) {
    mat <- mat > 0
    mat[is.na(mat)] <- FALSE
  }
  mat
}

.compact_gene_matrix_from_table <- function(x, storage = c("numeric", "logical")) {
  storage <- match.arg(storage)
  if (!is.data.frame(x)) .log_abort("Expected a data.frame for gene matrix conversion.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(x))) {
    .log_abort("Gene tables must include ensembl_gene_id and HGNC.")
  }
  gene_id <- as.character(x$HGNC)
  missing_gene <- is.na(gene_id) | !nzchar(gene_id)
  gene_id[missing_gene] <- as.character(x$ensembl_gene_id[missing_gene])
  keep <- !duplicated(gene_id)
  x <- x[keep, , drop = FALSE]
  gene_id <- gene_id[keep]
  cols <- setdiff(names(x), c("ensembl_gene_id", "HGNC"))
  mat <- data.matrix(x[, cols, drop = FALSE])
  rownames(mat) <- gene_id
  if (identical(storage, "logical")) {
    mat <- mat > 0
    mat[is.na(mat)] <- FALSE
  }
  mat
}

.compact_df_from_matrix <- function(mat, id_col) {
  mat <- as.matrix(mat)
  ids <- rownames(mat)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(mat)))
  out <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  out <- cbind(stats::setNames(data.frame(ids, stringsAsFactors = FALSE), id_col), out)
  tibble::as_tibble(out)
}

.compact_gene_df_from_matrix <- function(mat, genes, logical_out = FALSE) {
  mat <- as.matrix(mat)
  if (!is.data.frame(genes)) .log_abort("CraftGRN multiomic object is missing gene features.")
  gene_ids <- rownames(mat)
  if (is.null(gene_ids)) gene_ids <- as.character(seq_len(nrow(mat)))
  idx <- match(gene_ids, genes$gene_id)
  ensembl <- if ("ensembl_gene_id" %in% names(genes)) genes$ensembl_gene_id[idx] else gene_ids
  hgnc <- if ("hgnc_symbol" %in% names(genes)) genes$hgnc_symbol[idx] else gene_ids
  vals <- if (isTRUE(logical_out)) matrix(as.integer(mat > 0), nrow = nrow(mat), dimnames = dimnames(mat)) else mat
  out <- as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
  tibble::as_tibble(cbind(
    data.frame(ensembl_gene_id = ensembl, HGNC = hgnc, stringsAsFactors = FALSE),
    out
  ))
}

.compact_fp_features <- function(omics_data) {
  fp_ids <- rownames(.compact_matrix_from_table(omics_data$fp_score_condition_qn, "peak_ID"))
  ann <- omics_data$fp_annotation
  if (is.data.frame(ann) && all(c("fp_peak", "atac_peak") %in% names(ann))) {
    ann_base <- ann[!duplicated(ann$fp_peak), c("fp_peak", "atac_peak"), drop = FALSE]
    atac_peak <- ann_base$atac_peak[match(fp_ids, ann_base$fp_peak)]
  } else {
    atac_peak <- rep(NA_character_, length(fp_ids))
  }
  coords <- .module1_parse_fp_coordinates(fp_ids)
  tibble::tibble(
    fp_id = fp_ids,
    chr = coords$chr,
    start = coords$start,
    end = coords$end,
    atac_peak = as.character(atac_peak),
    group_size = NA_integer_
  )
}

.compact_fp_motifs <- function(omics_data) {
  empty <- tibble::tibble(fp_id = character(), motif = character(), tf = character())
  ann <- omics_data$fp_annotation
  if (!is.data.frame(ann)) return(empty)

  if (all(c("fp_peak", "motifs", "tfs") %in% names(ann))) {
    dt <- data.table::as.data.table(ann[, c("fp_peak", "motifs", "tfs"), drop = FALSE])
    data.table::setnames(dt, c("fp_peak", "motifs", "tfs"), c("fp_id", "motif", "tf"))
  } else if (all(c("fp_peak", "motifs") %in% names(ann))) {
    dt <- data.table::as.data.table(ann[, c("fp_peak", "motifs"), drop = FALSE])
    data.table::setnames(dt, c("fp_peak", "motifs"), c("fp_id", "motif"))
    gene_col <- .module1_motif_gene_col(omics_data$motif_db)
    if (is.data.frame(omics_data$motif_db) && !is.null(gene_col) && all(c("motif", gene_col) %in% names(omics_data$motif_db))) {
      map <- data.table::as.data.table(omics_data$motif_db[, c("motif", gene_col), drop = FALSE])
      data.table::setnames(map, c("motif", gene_col), c("motif", "tf"))
      map[, motif := sub("^_+", "", as.character(motif))]
      map <- unique(map[!is.na(motif) & nzchar(motif) & !is.na(tf) & nzchar(tf)])
      dt[, motif := sub("^_+", "", as.character(motif))]
      dt <- map[dt, on = "motif"]
    }
    if (!"tf" %in% names(dt) || all(is.na(dt$tf) | !nzchar(dt$tf))) {
      dt[, tf := toupper(sub("\\..*$", "", sub("^_+", "", as.character(motif))))]
      dt[, tf := gsub("[^A-Za-z0-9]+", "", tf)]
    }
  } else if (is.data.frame(omics_data$fp_tfs) && all(c("fp_peak", "tfs") %in% names(omics_data$fp_tfs))) {
    dt <- data.table::as.data.table(omics_data$fp_tfs[, c("fp_peak", "tfs"), drop = FALSE])
    data.table::setnames(dt, c("fp_peak", "tfs"), c("fp_id", "tf"))
    dt[, motif := NA_character_]
  } else {
    return(empty)
  }

  dt[, tf := gsub("\\s*::\\s*", ",", as.character(tf))]
  dt[, tf := strsplit(tf, ",", fixed = TRUE)]
  dt <- dt[, .(tf = unlist(tf, use.names = FALSE)), by = .(fp_id, motif)]
  dt[, tf := trimws(as.character(tf))]
  dt <- unique(dt[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf)])
  tibble::as_tibble(dt[, .(fp_id, motif, tf)])
}
.compact_gene_features <- function(omics_data) {
  rna <- omics_data$rna_condition
  if (!is.data.frame(rna)) .log_abort("`omics_data$rna_condition` is required.")
  gene_id <- as.character(rna$HGNC)
  missing_gene <- is.na(gene_id) | !nzchar(gene_id)
  gene_id[missing_gene] <- as.character(rna$ensembl_gene_id[missing_gene])
  keep <- !duplicated(gene_id)
  rna <- rna[keep, , drop = FALSE]
  gene_id <- gene_id[keep]
  tf <- if (is.null(omics_data$tf_list)) character() else unique(as.character(omics_data$tf_list))
  tibble::tibble(
    gene_id = gene_id,
    ensembl_gene_id = as.character(rna$ensembl_gene_id),
    hgnc_symbol = as.character(rna$HGNC),
    is_tf = gene_id %in% tf | as.character(rna$HGNC) %in% tf
  )
}

.compact_atac_features <- function(omics_data) {
  if (!is.data.frame(omics_data$atac_score) || !"atac_peak" %in% names(omics_data$atac_score)) {
    return(tibble::tibble(atac_peak = character(), chr = character(), start = integer(), end = integer()))
  }
  ids <- as.character(omics_data$atac_score$atac_peak)
  coords <- .module1_parse_fp_coordinates(ids)
  tibble::tibble(atac_peak = ids, chr = coords$chr, start = coords$start, end = coords$end)
}

.compact_samples <- function(omics_data, label_col = NULL) {
  meta <- omics_data$sample_metadata_used
  cond <- colnames(.compact_matrix_from_table(omics_data$fp_score_condition_qn, "peak_ID"))
  if (!is.data.frame(meta)) {
    return(tibble::tibble(condition_id = cond))
  }
  out <- tibble::as_tibble(meta)
  if ("id" %in% names(out) && !"sample_id" %in% names(out)) {
    out$sample_id <- as.character(out$id)
  }
  if (!is.null(label_col) && label_col %in% names(out)) {
    out$condition_id <- as.character(out[[label_col]])
  } else if (!"condition_id" %in% names(out)) {
    out$condition_id <- as.character(cond[seq_len(min(nrow(out), length(cond)))])
  }
  out
}

#' @noRd
as_multiomic_object <- function(omics_data,
                              project = list(),
                              paths = list(),
                              label_col = NULL,
                              verbose = TRUE) {
  if (is_multiomic_object(omics_data)) {
    validate_multiomic_object(omics_data)
    return(omics_data)
  }
  if (!is.list(omics_data)) .log_abort("`omics_data` must be a list.")
  required <- c("fp_score_condition_qn", "fp_bound_condition", "rna_condition", "rna_expressed")
  missing <- required[!vapply(required, function(nm) is.data.frame(omics_data[[nm]]), logical(1))]
  if (length(missing)) .log_abort("Cannot build CraftGRN multiomic object; missing: {paste(missing, collapse = ', ')}.")

  fp_score <- .compact_matrix_from_table(omics_data$fp_score_condition_qn, "peak_ID", "numeric")
  fp_bound <- .compact_matrix_from_table(omics_data$fp_bound_condition, "peak_ID", "logical")
  gene_expr <- .compact_gene_matrix_from_table(omics_data$rna_condition, "numeric")
  gene_on <- .compact_gene_matrix_from_table(omics_data$rna_expressed, "logical")

  object <- list(
    schema = .multiomic_object_schema(),
    project = project,
    samples = .compact_samples(omics_data, label_col = label_col),
    features = list(
      fp = .compact_fp_features(omics_data),
      fp_motif = .compact_fp_motifs(omics_data),
      atac = .compact_atac_features(omics_data),
      gene = .compact_gene_features(omics_data)
    ),
    matrices = list(
      fp_score = fp_score,
      fp_bound = fp_bound,
      gene_expr = gene_expr,
      gene_on = gene_on,
      atac_score = if (is.data.frame(omics_data$atac_score)) .compact_matrix_from_table(omics_data$atac_score, "atac_peak", "numeric") else NULL,
      atac_open = if (is.data.frame(omics_data$atac_overlap)) .compact_matrix_from_table(omics_data$atac_overlap, "atac_peak", "logical") else NULL
    ),
    refs = list(
      motif_db = omics_data$motif_db,
      tf = omics_data$tf_list
    ),
    qc = list(),
    paths = paths
  )
  class(object) <- c("craftgrn_multiomic", "list")
  validate_multiomic_object(object)
  if (isTRUE(verbose)) .log_inform("Built CraftGRN multiomic object.")
  object
}

#' @noRd
validate_multiomic_object <- function(omics_data) {
  if (!is_multiomic_object(omics_data)) .log_abort("Object is not a CraftGRN multiomic object.")
  mats <- omics_data$matrices
  feats <- omics_data$features
  if (!is.matrix(mats$fp_score)) .log_abort("`matrices$fp_score` must be a matrix.")
  if (!is.matrix(mats$fp_bound)) .log_abort("`matrices$fp_bound` must be a matrix.")
  if (!is.matrix(mats$gene_expr)) .log_abort("`matrices$gene_expr` must be a matrix.")
  if (!is.matrix(mats$gene_on)) .log_abort("`matrices$gene_on` must be a matrix.")
  if (!is.data.frame(feats$fp) || !"fp_id" %in% names(feats$fp)) .log_abort("`features$fp` must include fp_id.")
  if (!is.data.frame(feats$gene) || !"gene_id" %in% names(feats$gene)) .log_abort("`features$gene` must include gene_id.")
  if (!identical(rownames(mats$fp_score), as.character(feats$fp$fp_id))) {
    .log_abort("FP score rownames must match `features$fp$fp_id`.")
  }
  if (!identical(rownames(mats$fp_bound), as.character(feats$fp$fp_id))) {
    .log_abort("FP bound rownames must match `features$fp$fp_id`.")
  }
  if (!identical(rownames(mats$gene_expr), as.character(feats$gene$gene_id))) {
    .log_abort("Gene expression rownames must match `features$gene$gene_id`.")
  }
  if (!identical(rownames(mats$gene_on), as.character(feats$gene$gene_id))) {
    .log_abort("Gene-on rownames must match `features$gene$gene_id`.")
  }
  if (!identical(colnames(mats$fp_score), colnames(mats$fp_bound))) {
    .log_abort("FP score and FP bound condition columns must match.")
  }
  if (!identical(colnames(mats$gene_expr), colnames(mats$gene_on))) {
    .log_abort("Gene expression and gene-on condition columns must match.")
  }
  invisible(omics_data)
}



.as_module1_analysis_data <- function(omics_data, verbose = TRUE) {
  if (!is_multiomic_object(omics_data)) {
    .log_abort("`omics_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  }
  validate_multiomic_object(omics_data)
  fp <- data.table::as.data.table(omics_data$features$fp)
  genes <- data.table::as.data.table(omics_data$features$gene)
  motif <- data.table::as.data.table(omics_data$features$fp_motif)
  if (!is.data.frame(motif) || !all(c("fp_id", "motif", "tf") %in% names(motif))) {
    motif <- data.table::data.table(fp_id = character(), motif = character(), tf = character())
  }
  fp_keep <- fp[, .(fp_id = as.character(fp_id), atac_peak = as.character(atac_peak))]
  data.table::setkey(motif, fp_id)
  data.table::setkey(fp_keep, fp_id)
  ann <- fp_keep[motif, on = "fp_id"]
  fp_annotation <- data.table::data.table(
    fp_peak = as.character(ann$fp_id),
    atac_peak = as.character(ann$atac_peak),
    motifs = as.character(ann$motif),
    tfs = as.character(ann$tf)
  )
  out <- list(
    fp_score_condition_qn = .compact_df_from_matrix(omics_data$matrices$fp_score, "peak_ID"),
    fp_bound_condition = .compact_df_from_matrix(omics_data$matrices$fp_bound, "peak_ID"),
    fp_annotation = tibble::as_tibble(fp_annotation),
    rna_condition = .compact_gene_df_from_matrix(omics_data$matrices$gene_expr, genes),
    rna_expressed = .compact_gene_df_from_matrix(omics_data$matrices$gene_on, genes, logical_out = TRUE),
    tf_list = omics_data$refs$tf,
    motif_db = omics_data$refs$motif_db,
    sample_metadata_used = omics_data$samples,
    sample_labels = colnames(omics_data$matrices$fp_score),
    status = list(module1_analysis_view = TRUE)
  )
  if (isTRUE(verbose)) .log_inform("Prepared Module 1 analysis view from the CraftGRN multiomic object.")
  out
}
