# File: utils_step1_predicted_tfbs.R
# Purpose: Internal Module 1 helpers for compact TFBS prediction tables.

.module1_parse_fp_coordinates <- function(fp_id) {
  fp_id <- as.character(fp_id)
  chr <- sub(":.*$", "", fp_id, perl = TRUE)
  pos <- sub("^[^:]+:", "", fp_id, perl = TRUE)
  start <- suppressWarnings(as.integer(sub("-.*$", "", pos, perl = TRUE)))
  end <- suppressWarnings(as.integer(sub("^.*-", "", pos, perl = TRUE)))
  tibble::tibble(chr = chr, start = start, end = end)
}

.module1_best_corr <- function(pearson_r, spearman_r) {
  pearson_r <- as.numeric(pearson_r)
  spearman_r <- as.numeric(spearman_r)
  p_use <- is.finite(pearson_r)
  s_use <- is.finite(spearman_r)
  best_r <- rep(NA_real_, length(pearson_r))
  best_method <- rep(NA_character_, length(pearson_r))

  use_p <- p_use & (!s_use | pearson_r >= spearman_r)
  use_s <- s_use & (!p_use | spearman_r > pearson_r)
  best_r[use_p] <- pearson_r[use_p]
  best_method[use_p] <- "pearson"
  best_r[use_s] <- spearman_r[use_s]
  best_method[use_s] <- "spearman"

  list(best_r = best_r, best_method = best_method)
}

.module1_tfbs_cutoffs <- function(r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  r_cutoff <- as.numeric(r_cutoff)[[1L]]
  if (!is.finite(r_cutoff)) .log_abort("`r_cutoff` must be finite.")
  p_cutoff <- if (is.null(p_cutoff) || length(p_cutoff) == 0L || is.na(p_cutoff[[1L]])) Inf else as.numeric(p_cutoff)[[1L]]
  fdr_cutoff <- if (is.null(fdr_cutoff) || length(fdr_cutoff) == 0L || is.na(fdr_cutoff[[1L]])) Inf else as.numeric(fdr_cutoff)[[1L]]
  if (!is.finite(p_cutoff) && !identical(p_cutoff, Inf)) .log_abort("`p_cutoff` must be finite or NULL.")
  if (!is.finite(fdr_cutoff) && !identical(fdr_cutoff, Inf)) .log_abort("`fdr_cutoff` must be finite or NULL.")
  list(r = r_cutoff, p = p_cutoff, fdr = fdr_cutoff)
}

.module1_apply_tfbs_cutoffs <- function(x, cutoffs) {
  if (!is.data.frame(x)) .log_abort("`x` must be a data.frame.")
  if (!all(c("best_r", "best_method") %in% names(x))) {
    .log_abort("`x` must include best_r and best_method.")
  }
  pearson_best <- rep(FALSE, nrow(x))
  if ("pearson_p" %in% names(x)) {
    pearson_best <- is.finite(x$pearson_r) & (x$best_method == "pearson")
  }
  spearman_best <- rep(FALSE, nrow(x))
  if ("spearman_p" %in% names(x)) {
    spearman_best <- is.finite(x$spearman_r) & (x$best_method == "spearman")
  }
  best_p <- rep(NA_real_, nrow(x))
  best_fdr <- rep(NA_real_, nrow(x))
  if ("pearson_p" %in% names(x)) best_p[pearson_best] <- x$pearson_p[pearson_best]
  if ("spearman_p" %in% names(x)) best_p[spearman_best] <- x$spearman_p[spearman_best]
  if ("pearson_p_adj" %in% names(x)) best_fdr[pearson_best] <- x$pearson_p_adj[pearson_best]
  if ("spearman_p_adj" %in% names(x)) best_fdr[spearman_best] <- x$spearman_p_adj[spearman_best]
  x$best_p <- best_p
  x$best_fdr <- best_fdr
  x$pass_r <- is.finite(x$best_r) & x$best_r >= cutoffs$r
  x$pass_p <- is.infinite(cutoffs$p) | (is.finite(x$best_p) & x$best_p <= cutoffs$p)
  x$pass_fdr <- is.infinite(cutoffs$fdr) | (is.finite(x$best_fdr) & x$best_fdr <= cutoffs$fdr)
  x$pass <- x$pass_r & x$pass_p & x$pass_fdr
  tibble::as_tibble(x)
}

.module1_merge_correlation_stats <- function(pearson_stats, spearman_stats, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  keys <- c("fp_id", "atac_peak", "tf")
  if (!is.data.frame(pearson_stats)) .log_abort("`pearson_stats` must be a data.frame.")
  if (!is.data.frame(spearman_stats)) .log_abort("`spearman_stats` must be a data.frame.")
  if (!all(keys %in% names(pearson_stats))) .log_abort("`pearson_stats` is missing required key columns.")
  if (!all(keys %in% names(spearman_stats))) .log_abort("`spearman_stats` is missing required key columns.")

  collapse_stats <- function(x, value_cols) {
    keep <- intersect(unique(c(keys, value_cols, "motifs")), names(x))
    dt <- data.table::as.data.table(x[, keep, drop = FALSE])
    for (nm in keys) dt[, (nm) := as.character(get(nm))]
    value_cols <- intersect(value_cols, names(dt))
    if (!"motifs" %in% names(dt)) dt[, motifs := NA_character_]
    first_non_missing <- function(v) {
      idx <- which(!is.na(v))
      if (length(idx)) v[[idx[[1L]]]] else v[NA_integer_][[1L]]
    }
    out <- dt[
      ,
      c(
        lapply(.SD, first_non_missing),
        list(motifs = paste(sort(unique(as.character(motifs[!is.na(motifs) & nzchar(motifs)]))), collapse = ";"))
      ),
      by = keys,
      .SDcols = value_cols
    ]
    tibble::as_tibble(out)
  }

  pearson_tbl <- collapse_stats(pearson_stats, c("pearson_r", "pearson_p", "pearson_p_adj"))
  spearman_tbl <- collapse_stats(spearman_stats, c("spearman_r", "spearman_p", "spearman_p_adj"))

  merged <- dplyr::full_join(
    pearson_tbl,
    spearman_tbl,
    by = keys,
    suffix = c(".pearson", ".spearman")
  )

  if (!"pearson_r" %in% names(merged)) merged$pearson_r <- NA_real_
  if (!"spearman_r" %in% names(merged)) merged$spearman_r <- NA_real_
  best <- .module1_best_corr(merged$pearson_r, merged$spearman_r)
  merged$best_r <- best$best_r
  merged$best_method <- best$best_method
  .module1_apply_tfbs_cutoffs(
    merged,
    .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
  )
}

.module1_select_high_confidence_footprints <- function(motif_supported_correlations, r_cutoff = 0.3, p_cutoff = NULL, fdr_cutoff = NULL) {
  if (!is.data.frame(motif_supported_correlations)) {
    .log_abort("`motif_supported_correlations` must be a data.frame.")
  }
  need <- c("fp_id", "atac_peak")
  if (!all(need %in% names(motif_supported_correlations))) {
    .log_abort("`motif_supported_correlations` must include fp_id and atac_peak.")
  }
  if (!"best_r" %in% names(motif_supported_correlations)) {
    if (!all(c("pearson_r", "spearman_r") %in% names(motif_supported_correlations))) {
      .log_abort("`motif_supported_correlations` must include best_r or Pearson/Spearman columns.")
    }
    best <- .module1_best_corr(
      motif_supported_correlations$pearson_r,
      motif_supported_correlations$spearman_r
    )
    motif_supported_correlations$best_r <- best$best_r
    motif_supported_correlations$best_method <- best$best_method
  }

  if (!"pass" %in% names(motif_supported_correlations)) {
    motif_supported_correlations <- .module1_apply_tfbs_cutoffs(
      motif_supported_correlations,
      .module1_tfbs_cutoffs(r_cutoff = r_cutoff, p_cutoff = p_cutoff, fdr_cutoff = fdr_cutoff)
    )
  }
  keep <- motif_supported_correlations$pass %in% TRUE
  passed <- motif_supported_correlations[keep, , drop = FALSE]
  if (!nrow(passed)) {
    return(tibble::tibble(
      fp_id = character(), chr = character(), start = integer(), end = integer(),
      atac_peak = character(), n_canonical_bound_tfs = integer(),
      canonical_bound_tfs = character(), canonical_bound_motifs = character(),
      has_canonical_bound = logical()
    ))
  }
  motif_cols <- intersect(c("motifs", "motifs.pearson", "motifs.spearman"), names(passed))
  if (length(motif_cols)) {
    passed$motif_use <- passed[[motif_cols[[1L]]]]
  } else {
    passed$motif_use <- NA_character_
  }
  dt <- data.table::as.data.table(passed)
  out <- dt[, .(
    atac_peak = as.character(atac_peak[[1L]]),
    n_canonical_bound_tfs = data.table::uniqueN(tf),
    canonical_bound_tfs = paste(sort(unique(as.character(tf))), collapse = ";"),
    canonical_bound_motifs = paste(sort(unique(as.character(motif_use[!is.na(motif_use) & nzchar(motif_use)]))), collapse = ";"),
    has_canonical_bound = TRUE
  ), by = .(fp_id)]
  coords <- .module1_parse_fp_coordinates(out$fp_id)
  out_tbl <- tibble::as_tibble(out)
  tibble::as_tibble(cbind(out_tbl["fp_id"], coords, out_tbl[c("atac_peak", "n_canonical_bound_tfs", "canonical_bound_tfs", "canonical_bound_motifs", "has_canonical_bound")]))
}

.module1_condition_support <- function(fp_id, tf, omics_data) {
  if (!is.list(omics_data)) return(rep(NA_integer_, length(fp_id)))
  if (!is.data.frame(omics_data$fp_bound_condition) || !is.data.frame(omics_data$rna_expressed)) {
    return(rep(NA_integer_, length(fp_id)))
  }
  fp_bound <- omics_data$fp_bound_condition
  rna_expr <- omics_data$rna_expressed
  cond_cols <- intersect(
    setdiff(names(fp_bound), "peak_ID"),
    setdiff(names(rna_expr), c("ensembl_gene_id", "HGNC"))
  )
  if (!length(cond_cols)) return(rep(NA_integer_, length(fp_id)))

  fp_idx <- match(fp_id, fp_bound$peak_ID)
  rna_expr <- rna_expr[!duplicated(toupper(rna_expr$HGNC)), , drop = FALSE]
  tf_idx <- match(toupper(tf), toupper(rna_expr$HGNC))
  out <- rep(0L, length(fp_id))
  valid <- !is.na(fp_idx) & !is.na(tf_idx)
  if (!any(valid)) return(out)

  fp_mat <- as.matrix(fp_bound[, cond_cols, drop = FALSE]) > 0L
  tf_mat <- as.matrix(rna_expr[, cond_cols, drop = FALSE]) > 0L
  fp_mat[is.na(fp_mat)] <- FALSE
  tf_mat[is.na(tf_mat)] <- FALSE

  valid_idx <- which(valid)
  chunk_size <- 1000000L
  chunks <- split(valid_idx, ceiling(seq_along(valid_idx) / chunk_size))
  for (idx in chunks) {
    out[idx] <- rowSums(fp_mat[fp_idx[idx], , drop = FALSE] & tf_mat[tf_idx[idx], , drop = FALSE])
  }
  out
}

.module1_build_prediction_stats <- function(correlation_stats, high_confidence_footprints, omics_data = NULL) {
  if (!is.data.frame(correlation_stats)) .log_abort("`correlation_stats` must be a data.frame.")
  if (!is.data.frame(high_confidence_footprints)) .log_abort("`high_confidence_footprints` must be a data.frame.")
  need_stats <- c("fp_id", "atac_peak", "tf", "best_r", "best_method", "pass")
  need_high <- c("fp_id", "chr", "start", "end", "atac_peak")
  if (!all(need_stats %in% names(correlation_stats))) .log_abort("`correlation_stats` is missing required columns.")
  if (!all(need_high %in% names(high_confidence_footprints))) .log_abort("`high_confidence_footprints` is missing required columns.")

  high <- high_confidence_footprints[!duplicated(high_confidence_footprints$fp_id), need_high, drop = FALSE]
  stats_pass <- correlation_stats[correlation_stats$pass %in% TRUE, , drop = FALSE]
  if (!nrow(stats_pass)) {
    return(tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      best_r = numeric(),
      best_method = character(),
      condition_support = integer()
    ))
  }

  high_idx <- match(stats_pass$fp_id, high$fp_id)
  keep <- !is.na(high_idx)
  stats_pass <- stats_pass[keep, , drop = FALSE]
  high_idx <- high_idx[keep]
  if (!nrow(stats_pass)) {
    return(tibble::tibble(
      fp_id = character(),
      chr = character(),
      start = integer(),
      end = integer(),
      atac_peak = character(),
      tf = character(),
      best_r = numeric(),
      best_method = character(),
      condition_support = integer()
    ))
  }

  out <- data.table::as.data.table(stats_pass)
  out[, chr := as.character(high$chr[high_idx])]
  out[, start := as.integer(high$start[high_idx])]
  out[, end := as.integer(high$end[high_idx])]
  out[, condition_support := .module1_condition_support(fp_id, tf, omics_data)]
  keep_cols <- c(
    "fp_id", "chr", "start", "end", "atac_peak", "tf",
    intersect(c("motifs", "pearson_r", "spearman_r", "pearson_p", "spearman_p"), names(out)),
    "best_r", "best_method", "condition_support"
  )
  out <- out[, ..keep_cols]
  data.table::setorder(out, fp_id, tf)
  tibble::as_tibble(out)
}

# ---- Predicted TFBS handoff for Module 2 --------------------------------

.build_predicted_tfbs_table <- function(prediction_stats, include_support = TRUE, id_offset = 0L) {
  if (!is.data.frame(prediction_stats)) .log_abort("prediction_stats must be a data.frame.")
  need <- c("fp_id", "chr", "start", "end", "atac_peak", "tf")
  missing <- setdiff(need, names(prediction_stats))
  if (length(missing)) .log_abort("prediction_stats is missing required predicted TFBS columns: {paste(missing, collapse = ', ')}.")
  keep <- need
  if (isTRUE(include_support)) keep <- c(keep, intersect(c("condition_support"), names(prediction_stats)))
  out <- data.table::as.data.table(prediction_stats[, keep, drop = FALSE])
  out <- unique(out[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf)])
  out[, tfbs_id := sprintf("tfbs_%08d", as.integer(id_offset) + seq_len(.N))]
  data.table::setcolorder(out, c("tfbs_id", setdiff(names(out), "tfbs_id")))
  data.table::setorder(out, tf, fp_id)
  tibble::as_tibble(out)
}

#' Build the compact predicted TFBS handoff for Module 2
#'
#' @param prediction_stats Module 1 TFBS prediction statistic table.
#' @param include_support Include compact condition support when available.
#' @return A tibble with one row per predicted FP-TF binding event.
#' @noRd
build_predicted_tfbs <- function(prediction_stats, include_support = TRUE) {
  .build_predicted_tfbs_table(prediction_stats, include_support = include_support)
}

#' Output predicted TFBS
#'
#' @param prediction_stats Module 1 TFBS prediction statistic table.
#' @param out_dir Optional output directory. If supplied, a predicted TFBS table
#'   and manifest are written for Module 2.
#' @param output_format Output format: auto, parquet, or csv.
#' @param include_support Include compact condition support when available.
#' @return A predicted TFBS tibble when `out_dir` is NULL; otherwise a list with
#'   output paths and row counts.
#' @export
output_predicted_tfbs <- function(prediction_stats, out_dir = NULL, output_format = c("auto", "parquet", "csv"), include_support = TRUE) {
  predicted_tfbs <- build_predicted_tfbs(prediction_stats, include_support = include_support)
  if (is.null(out_dir) || !nzchar(out_dir)) {
    return(predicted_tfbs)
  }
  .write_predicted_tfbs_table(predicted_tfbs, out_dir = out_dir, output_format = output_format)
}

.predicted_tfbs_output_format <- function(output_format = c("auto", "parquet", "csv")) {
  output_format <- match.arg(output_format)
  if (identical(output_format, "auto")) {
    if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "csv"
  } else {
    output_format
  }
}

.write_predicted_tfbs_table <- function(predicted_tfbs, out_dir, output_format = c("auto", "parquet", "csv"), filename = NULL) {
  output_format <- .predicted_tfbs_output_format(output_format)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(output_format, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
    path <- file.path(out_dir, if (is.null(filename)) "module1_predicted_tfbs.parquet" else filename)
    .write_parquet_table(predicted_tfbs, path)
  } else {
    output_format <- "csv"
    path <- file.path(out_dir, if (is.null(filename)) "module1_predicted_tfbs.csv.gz" else filename)
    readr::write_csv(predicted_tfbs, path)
  }
  manifest <- tibble::tibble(table = "module1_predicted_tfbs", path = path, format = output_format, n_rows = nrow(predicted_tfbs))
  manifest_path <- file.path(out_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(manifest, manifest_path)
  list(path = path, manifest = manifest_path, format = output_format, n_rows = nrow(predicted_tfbs))
}

.write_predicted_tfbs_from_prediction_stats_manifest <- function(prediction_stats_manifest, out_dir, output_format = c("auto", "parquet", "csv")) {
  if (!is.data.frame(prediction_stats_manifest)) .log_abort("prediction_stats_manifest must be a data.frame.")
  if (!all(c("path", "format") %in% names(prediction_stats_manifest))) .log_abort("prediction_stats_manifest must contain path and format columns.")
  output_format <- .predicted_tfbs_output_format(output_format)
  chunk_dir <- file.path(out_dir, "module1_predicted_tfbs_chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  offset <- 0L
  manifest <- vector("list", nrow(prediction_stats_manifest))
  for (i in seq_len(nrow(prediction_stats_manifest))) {
    stats_path <- as.character(prediction_stats_manifest$path[[i]])
    stats_format <- as.character(prediction_stats_manifest$format[[i]])
    if (!file.exists(stats_path)) .log_abort("TFBS prediction statistic chunk not found: {stats_path}")
    if (identical(stats_format, "parquet") || grepl("\\.parquet$", stats_path, ignore.case = TRUE)) {
      if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet TFBS prediction statistic chunks.")
      stats_i <- tibble::as_tibble(arrow::read_parquet(
        stats_path,
        col_select = c("fp_id", "chr", "start", "end", "atac_peak", "tf", "condition_support")
      ))
    } else {
      stats_i <- tibble::as_tibble(readr::read_csv(stats_path, show_col_types = FALSE))
    }
    pred_i <- .build_predicted_tfbs_table(stats_i, include_support = TRUE, id_offset = offset)
    offset <- offset + nrow(pred_i)
    if (identical(output_format, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
      out_path <- file.path(chunk_dir, sprintf("module1_predicted_tfbs_chunk_%04d.parquet", i))
      .write_parquet_table(pred_i, out_path)
      fmt <- "parquet"
    } else {
      fmt <- "csv"
      out_path <- file.path(chunk_dir, sprintf("module1_predicted_tfbs_chunk_%04d.csv.gz", i))
      readr::write_csv(pred_i, out_path)
    }
    manifest[[i]] <- tibble::tibble(chunk_id = i, path = out_path, format = fmt, n_rows = nrow(pred_i))
    rm(stats_i, pred_i)
  }
  manifest <- dplyr::bind_rows(manifest)
  manifest_path <- file.path(out_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(manifest, manifest_path)
  list(path = manifest_path, manifest = manifest_path, chunks = chunk_dir, format = output_format, n_rows = sum(manifest$n_rows))
}

#' Load TFBS predicted from Module 1
#'
#' @param path Path to a predicted TFBS manifest, Parquet file, or CSV file.
#' @return A tibble.
#' @export
load_predicted_tfbs <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) .log_abort("path must be a non-empty file path.")
  if (!file.exists(path)) .log_abort("Predicted TFBS path not found: {path}")
  if (grepl("_manifest\\.csv$", basename(path))) {
    man <- data.table::fread(path, showProgress = FALSE)
    if (!all(c("path", "format") %in% names(man)) || !nrow(man)) .log_abort("Invalid predicted TFBS manifest: {path}")
    read_one <- function(path_i, format_i) {
      if (identical(format_i, "parquet") || grepl("\\.parquet$", path_i, ignore.case = TRUE)) {
        if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet predicted TFBS output.")
        return(tibble::as_tibble(arrow::read_parquet(path_i)))
      }
      tibble::as_tibble(readr::read_csv(path_i, show_col_types = FALSE))
    }
    return(dplyr::bind_rows(lapply(seq_len(nrow(man)), function(i) read_one(as.character(man$path[[i]]), as.character(man$format[[i]])))))
  }
  if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet predicted TFBS output.")
    return(tibble::as_tibble(arrow::read_parquet(path)))
  }
  tibble::as_tibble(readr::read_csv(path, show_col_types = FALSE))
}

#' Export predicted TFBS as BED files
#'
#' @param predicted_tfbs Compact predicted TFBS table or path.
#' @param out_file BED output path. Required when split_by is none.
#' @param out_dir Directory for split BED outputs.
#' @param tf Optional TF subset.
#' @param split_by One of none or tf.
#' @return Output path or manifest tibble, invisibly.
#' @export
export_predicted_tfbs_bed <- function(predicted_tfbs, out_file = NULL, out_dir = NULL, tf = NULL, split_by = c("none", "tf")) {
  split_by <- match.arg(split_by)
  if (is.character(predicted_tfbs) && length(predicted_tfbs) == 1L && file.exists(predicted_tfbs)) predicted_tfbs <- load_predicted_tfbs(predicted_tfbs)
  if (!is.data.frame(predicted_tfbs)) .log_abort("predicted_tfbs must be a data.frame or path.")
  need <- c("chr", "start", "end", "tf", "fp_id")
  if (!all(need %in% names(predicted_tfbs))) .log_abort("predicted_tfbs is missing required BED columns.")
  dt <- data.table::as.data.table(predicted_tfbs)
  if (!is.null(tf)) {
    tf_filter <- unique(as.character(tf))
    tf_filter <- tf_filter[nzchar(tf_filter)]
    dt <- dt[as.character(dt[["tf"]]) %in% tf_filter]
  }
  make_bed <- function(x) data.table::data.table(chrom = as.character(x$chr), chromStart = as.integer(x$start), chromEnd = as.integer(x$end), name = paste(as.character(x$tf), as.character(x$fp_id), sep = "|"), score = 1000L, strand = ".")
  if (identical(split_by, "none")) {
    if (is.null(out_file) || !nzchar(out_file)) .log_abort("out_file is required when split_by is none.")
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(make_bed(dt), out_file, sep = "\t", col.names = FALSE)
    return(invisible(out_file))
  }
  if (is.null(out_dir) || !nzchar(out_dir)) .log_abort("out_dir is required when split_by is tf.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tfs <- sort(unique(as.character(dt$tf)))
  manifest <- lapply(tfs, function(tf_i) {
    x <- dt[as.character(dt[["tf"]]) == tf_i]
    path <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]+", "_", tf_i), ".bed"))
    data.table::fwrite(make_bed(x), path, sep = "\t", col.names = FALSE)
    data.frame(tf = tf_i, path = path, n_rows = nrow(x), stringsAsFactors = FALSE)
  })
  manifest <- dplyr::bind_rows(manifest)
  readr::write_csv(manifest, file.path(out_dir, "predicted_tfbs_bed_manifest.csv"))
  invisible(tibble::as_tibble(manifest))
}
