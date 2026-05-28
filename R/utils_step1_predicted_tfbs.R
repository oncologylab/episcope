# File: utils_step1_predicted_tfbs.R
# Purpose: Compact Module 1 predicted TFBS handoff helpers.

.build_predicted_tfbs_table <- function(tfbs_links, include_support = TRUE, id_offset = 0L) {
  if (!is.data.frame(tfbs_links)) .log_abort("tfbs_links must be a data.frame.")
  need <- c("fp_id", "chr", "start", "end", "atac_peak", "tf")
  missing <- setdiff(need, names(tfbs_links))
  if (length(missing)) .log_abort("tfbs_links is missing required predicted TFBS columns: {paste(missing, collapse = ', ')}.")
  keep <- need
  if (isTRUE(include_support)) keep <- c(keep, intersect(c("condition_support", "sample_support"), names(tfbs_links)))
  out <- data.table::as.data.table(tfbs_links[, keep, drop = FALSE])
  out <- unique(out[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf)])
  out[, tfbs_id := sprintf("tfbs_%08d", as.integer(id_offset) + seq_len(.N))]
  data.table::setcolorder(out, c("tfbs_id", setdiff(names(out), "tfbs_id")))
  data.table::setorder(out, tf, fp_id)
  tibble::as_tibble(out)
}

#' Build the compact predicted TFBS handoff for Module 2
#'
#' @param tfbs_links Module 1 TFBS link table.
#' @param include_support Include compact condition support when available.
#' @return A tibble with one row per predicted FP-TF binding event.
#' @export
build_predicted_tfbs <- function(tfbs_links, include_support = TRUE) {
  .build_predicted_tfbs_table(tfbs_links, include_support = include_support)
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
    arrow::write_parquet(predicted_tfbs, path, compression = "zstd")
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

.write_predicted_tfbs_from_link_manifest <- function(link_manifest, out_dir, output_format = c("auto", "parquet", "csv")) {
  if (!is.data.frame(link_manifest)) .log_abort("link_manifest must be a data.frame.")
  if (!all(c("path", "format") %in% names(link_manifest))) .log_abort("link_manifest must contain path and format columns.")
  output_format <- .predicted_tfbs_output_format(output_format)
  chunk_dir <- file.path(out_dir, "module1_predicted_tfbs_chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
  offset <- 0L
  manifest <- vector("list", nrow(link_manifest))
  for (i in seq_len(nrow(link_manifest))) {
    link_path <- as.character(link_manifest$path[[i]])
    link_format <- as.character(link_manifest$format[[i]])
    if (!file.exists(link_path)) .log_abort("TFBS link chunk not found: {link_path}")
    if (identical(link_format, "parquet") || grepl("\\.parquet$", link_path, ignore.case = TRUE)) {
      if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet TFBS link chunks.")
      links_i <- tibble::as_tibble(arrow::read_parquet(link_path))
    } else {
      links_i <- tibble::as_tibble(readr::read_csv(link_path, show_col_types = FALSE))
    }
    pred_i <- .build_predicted_tfbs_table(links_i, include_support = TRUE, id_offset = offset)
    offset <- offset + nrow(pred_i)
    if (identical(output_format, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
      out_path <- file.path(chunk_dir, sprintf("module1_predicted_tfbs_chunk_%04d.parquet", i))
      arrow::write_parquet(pred_i, out_path, compression = "zstd")
      fmt <- "parquet"
    } else {
      fmt <- "csv"
      out_path <- file.path(chunk_dir, sprintf("module1_predicted_tfbs_chunk_%04d.csv.gz", i))
      readr::write_csv(pred_i, out_path)
    }
    manifest[[i]] <- tibble::tibble(chunk_id = i, path = out_path, format = fmt, n_rows = nrow(pred_i))
    rm(links_i, pred_i)
  }
  manifest <- dplyr::bind_rows(manifest)
  manifest_path <- file.path(out_dir, "module1_predicted_tfbs_manifest.csv")
  readr::write_csv(manifest, manifest_path)
  list(path = manifest_path, manifest = manifest_path, chunks = chunk_dir, format = output_format, n_rows = sum(manifest$n_rows))
}

#' Load compact predicted TFBS output
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
  if (!is.null(tf)) dt <- dt[tf %in% as.character(tf)]
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
    x <- dt[tf == tf_i]
    path <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]+", "_", tf_i), ".bed"))
    data.table::fwrite(make_bed(x), path, sep = "\t", col.names = FALSE)
    data.frame(tf = tf_i, path = path, n_rows = nrow(x), stringsAsFactors = FALSE)
  })
  manifest <- dplyr::bind_rows(manifest)
  readr::write_csv(manifest, file.path(out_dir, "predicted_tfbs_bed_manifest.csv"))
  invisible(tibble::as_tibble(manifest))
}
