# Internal adapter for fp-tools compact motif-site caches.

.fp_tools_compact_paths <- function(root_dir, sample_id) {
  base <- file.path(root_dir, sample_id, "match_motifs")
  list(
    sites = file.path(base, "cache", "motif_sites.tsv.gz"),
    thresholds = file.path(base, "cache", "thresholds.json"),
    motif_results = file.path(base, "motif_matches_results.txt")
  )
}

.fp_tools_compact_threshold <- function(path) {
  if (!file.exists(path)) .log_abort("Missing fp-tools threshold cache: {.path {path}}.")
  parsed <- jsonlite::read_json(path, simplifyVector = TRUE)
  values <- suppressWarnings(as.numeric(unlist(parsed$thresholds, use.names = FALSE)))
  values <- values[is.finite(values)]
  if (length(values) != 1L) {
    .log_abort("Expected one finite footprint threshold in {.path {path}}.")
  }
  values[[1L]]
}

.fp_tools_compact_sites <- function(path, include_motifs = FALSE) {
  if (!file.exists(path)) .log_abort("Missing fp-tools motif-site cache: {.path {path}}.")
  x <- data.table::fread(
    path,
    select = c(
      "motif", "TFBS_chr", "TFBS_start", "TFBS_end",
      "peak_chr", "peak_start", "peak_end", "score"
    ),
    showProgress = FALSE
  )
  x[, `:=`(
    peak_ID = paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end),
    atac_peak = paste0(peak_chr, ":", peak_start, "-", peak_end),
    score = as.numeric(score)
  )]
  grouped <- if (isTRUE(include_motifs)) {
    x[, .(
      score = mean(score, na.rm = TRUE),
      motifs_all = paste(sort(unique(as.character(motif))), collapse = ";")
    ), by = .(peak_ID, atac_peak)]
  } else {
    x[, .(score = mean(score, na.rm = TRUE)), by = .(peak_ID, atac_peak)]
  }
  if (anyDuplicated(grouped$peak_ID)) {
    .log_abort("fp-tools cache maps one TFBS coordinate to multiple ATAC peaks: {.path {path}}.")
  }
  grouped
}

.fp_tools_compact_header_signature <- function(path, n = 100L) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, n = as.integer(n) + 1L, warn = FALSE)
  if (length(lines) < 2L) .log_abort("Empty fp-tools motif-site cache: {.path {path}}.")
  lines <- sub("\t[^\t]*$", "", lines)
  digest::digest(lines, algo = "xxhash64", serialize = FALSE)
}

.fp_tools_compact_reference <- function(path) {
  x <- data.table::fread(
    path,
    select = c(
      "motif", "TFBS_chr", "TFBS_start", "TFBS_end",
      "peak_chr", "peak_start", "peak_end", "score"
    ),
    showProgress = FALSE
  )
  x[, `:=`(
    peak_ID = paste0(TFBS_chr, ":", TFBS_start, "-", TFBS_end),
    atac_peak = paste0(peak_chr, ":", peak_start, "-", peak_end),
    score = as.numeric(score)
  )]
  grouped <- x[, .(
    score = mean(score, na.rm = TRUE),
    motifs_all = paste(sort(unique(as.character(motif))), collapse = ";")
  ), by = .(peak_ID, atac_peak)]
  if (anyDuplicated(grouped$peak_ID)) {
    .log_abort("fp-tools cache maps one TFBS coordinate to multiple ATAC peaks: {.path {path}}.")
  }
  group_index <- match(
    paste(x$peak_ID, x$atac_peak, sep = "|"),
    paste(grouped$peak_ID, grouped$atac_peak, sep = "|")
  )
  list(
    grouped = grouped,
    group_index = group_index,
    group_size = tabulate(group_index, nbins = nrow(grouped)),
    n_raw_rows = nrow(x),
    signature = .fp_tools_compact_header_signature(path)
  )
}

.fp_tools_compact_score <- function(path, reference) {
  if (!identical(.fp_tools_compact_header_signature(path), reference$signature)) {
    .log_abort("fp-tools motif-site row signature differs: {.path {path}}.")
  }
  values <- data.table::fread(path, select = "score", showProgress = FALSE)$score
  if (length(values) != reference$n_raw_rows) {
    .log_abort("fp-tools motif-site row count differs: {.path {path}}.")
  }
  values <- as.numeric(values)
  values[!is.finite(values)] <- 0
  as.numeric(rowsum(values, reference$group_index, reorder = FALSE)) / reference$group_size
}

.fp_tools_compact_motif_db <- function(path) {
  if (!file.exists(path)) .log_abort("Missing fp-tools motif result table: {.path {path}}.")
  x <- data.table::fread(path, showProgress = FALSE)
  required <- c("output_prefix", "name", "motif_id")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    .log_abort("fp-tools motif result table is missing: {paste(missing, collapse = ', ')}.")
  }
  unique(x[, .(
    motif = as.character(output_prefix),
    gene_symbol = as.character(name),
    motif_id = as.character(motif_id)
  )])
}

.load_fp_tools_compact_footprints <- function(root_dir,
                                               sample_ids,
                                               cache_dir = NULL,
                                               cache_tag = "fp_tools_compact",
                                               verbose = TRUE) {
  sample_ids <- unique(as.character(sample_ids))
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  if (!length(sample_ids)) .log_abort("No sample IDs were supplied for fp-tools compact import.")
  if (!dir.exists(root_dir)) .log_abort("fp-tools sample root does not exist: {.path {root_dir}}.")

  if (!is.null(cache_dir)) {
    cached <- .aligned_fp_cache_choose_format(cache_dir, cache_tag, "auto")
    if (isTRUE(cached$exists)) {
      out <- load_fp_aligned_from_cache(
        cache_dir = cache_dir,
        cache_tag = cache_tag,
        cache_format = "auto",
        verbose = verbose
      )
      motif_path <- file.path(cache_dir, sprintf("motif_db_%s.csv", cache_tag))
      if (file.exists(motif_path)) out$motif_db <- tibble::as_tibble(data.table::fread(motif_path))
      return(out)
    }
  }

  paths <- lapply(sample_ids, function(id) .fp_tools_compact_paths(root_dir, id))
  names(paths) <- sample_ids
  missing <- unlist(lapply(paths, function(x) unlist(x)[!file.exists(unlist(x))]), use.names = FALSE)
  if (length(missing)) {
    .log_abort("Missing fp-tools compact input files, beginning with: {.path {missing[[1L]]}}.")
  }
  if (isTRUE(verbose)) {
    .log_inform("Importing fp-tools compact footprints for {length(sample_ids)} sample(s).")
  }

  reference <- .fp_tools_compact_reference(paths[[1L]]$sites)
  first <- reference$grouped
  .resource_preflight(
    estimated_bytes = nrow(first) * (length(sample_ids) * 16 + 512),
    stage = "fp-tools compact footprint import",
    verbose = verbose
  )
  score <- data.table::data.table(peak_ID = first$peak_ID)
  bound <- data.table::data.table(peak_ID = first$peak_ID)
  for (sample_id in sample_ids) {
    values <- if (identical(sample_id, sample_ids[[1L]])) {
      as.numeric(first$score)
    } else {
      .fp_tools_compact_score(paths[[sample_id]]$sites, reference)
    }
    score[[sample_id]] <- values
    bound[[sample_id]] <- values >= .fp_tools_compact_threshold(paths[[sample_id]]$thresholds)
  }
  fp_sites <- first[, .(
    peak_ID,
    atac_peak,
    motifs_all,
    source_fp_peaks = peak_ID,
    n_source_fp_peaks = 1L,
    chr = sub(":.*$", "", peak_ID),
    start = as.integer(sub("-.*$", "", sub("^[^:]+:", "", peak_ID))),
    end = as.integer(sub("^.*-", "", peak_ID))
  )]
  motif_db <- .fp_tools_compact_motif_db(paths[[1L]]$motif_results)
  motif_parts <- strsplit(fp_sites$motifs_all, ";", fixed = TRUE)
  fp_annotation <- data.table::data.table(
    fp_peak = rep(fp_sites$peak_ID, lengths(motif_parts)),
    atac_peak = rep(fp_sites$atac_peak, lengths(motif_parts)),
    motifs = unlist(motif_parts, use.names = FALSE)
  )
  out <- list(
    fp_score = tibble::as_tibble(score),
    fp_bound = tibble::as_tibble(bound),
    fp_annotation = tibble::as_tibble(fp_annotation),
    fp_sites = tibble::as_tibble(fp_sites),
    id_map = tibble::tibble(),
    motif_db = tibble::as_tibble(motif_db),
    cache_schema = "fp_tools_compact"
  )
  rm(reference)
  if (!is.null(cache_dir)) {
    format <- if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "csv"
    cache_paths <- .aligned_fp_cache_paths(cache_dir, cache_tag, format)
    .aligned_fp_write_table(score, cache_paths$fp_score, format)
    .aligned_fp_write_table(bound, cache_paths$fp_bound, format)
    .aligned_fp_write_table(fp_sites, cache_paths$fp_sites, format)
    data.table::fwrite(motif_db, file.path(cache_dir, sprintf("motif_db_%s.csv", cache_tag)))
  }
  out
}
