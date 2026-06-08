# File: utils_step2_link_tf_targets.R
# Purpose: Compact relational Module 2 TF-target linking.

.module2_default_cores <- function(n_cores = NULL) {
  if (is.null(n_cores) || length(n_cores) == 0L || is.na(n_cores[[1L]])) {
    return(max(1L, parallel::detectCores(logical = TRUE)))
  }
  max(1L, as.integer(n_cores[[1L]]))
}

.module2_cfg <- function(project_config = NULL) {
  cfg <- list()
  if (is.character(project_config) && length(project_config) == 1L && nzchar(project_config) && file.exists(project_config)) {
    cfg <- yaml::read_yaml(project_config)
  } else if (is.list(project_config)) {
    cfg <- project_config
  }
  if (is.list(cfg$module2)) cfg <- utils::modifyList(cfg, cfg$module2)
  cfg
}

.module2_cfg_value <- function(cfg, name, default = NULL) {
  if (is.list(cfg) && !is.null(cfg[[name]])) cfg[[name]] else default
}

.module2_corr_cutoffs <- function(cfg, prefix, r_default = 0.3) {
  list(
    r = as.numeric(.module2_cfg_value(cfg, paste0("threshold_", prefix, "_corr_r"), r_default))[[1L]],
    p = .module2_cfg_value(cfg, paste0("threshold_", prefix, "_corr_p"), NULL),
    fdr = .module2_cfg_value(cfg, paste0("threshold_", prefix, "_corr_fdr"), NULL)
  )
}

.module2_empty_corr <- function(id1_col, id2_col) {
  out <- data.frame(
    id1 = character(),
    id2 = character(),
    pearson_r = numeric(),
    pearson_p = numeric(),
    pearson_fdr = numeric(),
    spearman_r = numeric(),
    spearman_p = numeric(),
    spearman_fdr = numeric(),
    best_r = numeric(),
    best_method = character(),
    best_p = numeric(),
    best_fdr = numeric(),
    pass = logical(),
    stringsAsFactors = FALSE
  )
  names(out)[1:2] <- c(id1_col, id2_col)
  tibble::as_tibble(out)
}

.module2_apply_corr_cutoffs <- function(x, cutoffs) {
  if (!nrow(x)) return(x)
  best <- .module1_best_corr(x$pearson_r, x$spearman_r)
  x$best_r <- best$best_r
  x$best_method <- best$best_method
  use_p <- !is.na(x$best_method) & x$best_method == "pearson"
  use_s <- !is.na(x$best_method) & x$best_method == "spearman"
  x$best_p <- NA_real_
  x$best_fdr <- NA_real_
  x$best_p[use_p] <- x$pearson_p[use_p]
  x$best_p[use_s] <- x$spearman_p[use_s]
  x$best_fdr[use_p] <- x$pearson_fdr[use_p]
  x$best_fdr[use_s] <- x$spearman_fdr[use_s]
  r_cut <- as.numeric(cutoffs$r)[[1L]]
  if (!is.finite(r_cut)) r_cut <- 0.3
  p_cut <- cutoffs$p
  fdr_cut <- cutoffs$fdr
  p_ok <- is.null(p_cut) || length(p_cut) == 0L || is.na(p_cut[[1L]])
  fdr_ok <- is.null(fdr_cut) || length(fdr_cut) == 0L || is.na(fdr_cut[[1L]])
  x$pass <- is.finite(x$best_r) & x$best_r >= r_cut
  if (!p_ok) x$pass <- x$pass & is.finite(x$best_p) & x$best_p <= as.numeric(p_cut)[[1L]]
  if (!fdr_ok) x$pass <- x$pass & is.finite(x$best_fdr) & x$best_fdr <= as.numeric(fdr_cut)[[1L]]
  tibble::as_tibble(x)
}

.module2_pair_correlations <- function(x_mat, y_mat, pairs, x_col, y_col, cutoffs, n_cores = NULL) {
  if (!is.data.frame(pairs) || !nrow(pairs)) return(.module2_empty_corr(x_col, y_col))
  x_ids <- rownames(x_mat)
  y_ids <- rownames(y_mat)
  xi <- match(as.character(pairs[[x_col]]), x_ids)
  yi <- match(as.character(pairs[[y_col]]), y_ids)
  valid <- !is.na(xi) & !is.na(yi)
  pairs <- pairs[valid, c(x_col, y_col), drop = FALSE]
  xi <- xi[valid]
  yi <- yi[valid]
  if (!nrow(pairs)) return(.module2_empty_corr(x_col, y_col))
  ux <- unique(xi)
  uy <- unique(yi)
  x_sub <- as.matrix(x_mat[ux, , drop = FALSE])
  y_sub <- as.matrix(y_mat[uy, , drop = FALSE])
  x_rank <- .module1_rank_matrix_rows(x_sub)
  y_rank <- .module1_rank_matrix_rows(y_sub)
  cpp_fun <- get0(".sparse_pair_correlations_cpp", envir = asNamespace("craftgrn"), mode = "function")
  if (!is.function(cpp_fun)) .log_abort("C++ sparse pair correlation backend is unavailable.")
  stats <- tibble::as_tibble(cpp_fun(x_sub, y_sub, x_rank, y_rank, as.integer(match(xi, ux) - 1L), as.integer(match(yi, uy) - 1L), .module2_default_cores(n_cores)))
  out <- tibble::as_tibble(cbind(pairs, stats))
  out$pearson_fdr <- stats::p.adjust(out$pearson_p, method = "BH")
  out$spearman_fdr <- stats::p.adjust(out$spearman_p, method = "BH")
  .module2_apply_corr_cutoffs(out, cutoffs)
}

.module2_write_table <- function(x, out_dir, name, output_format = c("auto", "parquet", "csv")) {
  fmt <- .predicted_tfbs_output_format(output_format)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(fmt, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
    path <- file.path(out_dir, paste0(name, ".parquet"))
    .write_parquet_table(x, path)
  } else {
    fmt <- "csv"
    path <- file.path(out_dir, paste0(name, ".csv.gz"))
    readr::write_csv(x, path)
  }
  tibble::tibble(table = name, path = path, format = fmt, n_rows = nrow(x))
}

.module2_manifest_table <- function(x, table_name) {
  x$table <- table_name
  x
}

.module2_write_run_summary <- function(x, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, "module2_run_summary.csv")
  readr::write_csv(x, path)
  tibble::tibble(table = "module2_qc_summary", path = path, format = "csv", n_rows = nrow(x))
}

.module2_normalize_gene_tss <- function(gene_tss) {
  if (!is.data.frame(gene_tss)) .log_abort("gene_tss must be a data.frame.")
  nms <- names(gene_tss)
  gene_col <- intersect(c("target_gene", "gene_id", "gene", "HGNC", "hgnc_symbol", "ensembl_gene_id"), nms)[1L]
  chr_col <- intersect(c("target_chr", "chr", "chrom", "seqnames"), nms)[1L]
  tss_col <- intersect(c("target_tss", "tss", "gene_tss"), nms)[1L]
  strand_col <- intersect(c("target_strand", "strand", "gene_strand"), nms)[1L]
  if (is.na(gene_col) || is.na(chr_col) || is.na(tss_col)) .log_abort("gene_tss must include gene, chromosome, and TSS columns.")
  out <- tibble::tibble(target_gene = as.character(gene_tss[[gene_col]]), target_chr = as.character(gene_tss[[chr_col]]), target_tss = as.integer(gene_tss[[tss_col]]), target_strand = if (is.na(strand_col)) "." else as.character(gene_tss[[strand_col]]))
  out <- out[!is.na(out$target_gene) & nzchar(out$target_gene) & !is.na(out$target_chr) & is.finite(out$target_tss), , drop = FALSE]
  out[!duplicated(out$target_gene), , drop = FALSE]
}

.module2_normalize_prior <- function(regulatory_prior = NULL, predicted_tfbs = NULL) {
  if (is.null(regulatory_prior)) return(tibble::tibble(fp_id = character(), target_gene = character(), prior_id = character(), prior_source = character(), prior_score = numeric(), prior_status = character()))
  if (!is.data.frame(regulatory_prior)) .log_abort("regulatory_prior must be a data.frame when supplied.")
  nms <- names(regulatory_prior)
  gene_col <- intersect(c("target_gene", "gene_key", "gene", "HGNC", "connected_gene", "ensembl_gene_id"), nms)[1L]
  fp_col <- intersect(c("fp_id", "peak_ID"), nms)[1L]
  atac_col <- intersect(c("atac_peak", "peak"), nms)[1L]
  if (is.na(gene_col)) .log_abort("regulatory_prior must include a target gene column.")
  prior <- tibble::tibble(target_gene = as.character(regulatory_prior[[gene_col]]))
  if (!is.na(fp_col)) {
    prior$fp_id <- as.character(regulatory_prior[[fp_col]])
  } else if (!is.na(atac_col) && is.data.frame(predicted_tfbs) && "atac_peak" %in% names(predicted_tfbs)) {
    map <- unique(predicted_tfbs[, c("fp_id", "atac_peak"), drop = FALSE])
    prior$atac_peak <- as.character(regulatory_prior[[atac_col]])
    prior <- dplyr::left_join(prior, map, by = "atac_peak")
  } else {
    prior$fp_id <- NA_character_
  }
  prior$prior_id <- if ("prior_id" %in% nms) as.character(regulatory_prior$prior_id) else sprintf("prior_%08d", seq_len(nrow(prior)))
  prior$prior_source <- if ("prior_source" %in% nms) as.character(regulatory_prior$prior_source) else "regulatory_prior"
  prior$prior_score <- if ("prior_score" %in% nms) suppressWarnings(as.numeric(regulatory_prior$prior_score)) else NA_real_
  prior$prior_status <- if ("prior_status" %in% nms) as.character(regulatory_prior$prior_status) else "supported"
  prior <- prior[!is.na(prior$fp_id) & nzchar(prior$fp_id) & !is.na(prior$target_gene) & nzchar(prior$target_gene), , drop = FALSE]
  unique(prior[, c("fp_id", "target_gene", "prior_id", "prior_source", "prior_score", "prior_status"), drop = FALSE])
}

.module2_empty_candidates <- function() {
  tibble::tibble(
    candidate_id = character(),
    fp_id = character(),
    target_gene = character(),
    chr = character(),
    start = integer(),
    end = integer(),
    atac_peak = character(),
    target_chr = character(),
    target_tss = integer(),
    target_strand = character(),
    distance_to_tss = numeric(),
    candidate_source = character(),
    within_tss_window = logical(),
    prior_supported = logical(),
    prior_id = character(),
    prior_source = character(),
    prior_score = numeric(),
    prior_status = character()
  )
}

.module2_build_candidates <- function(predicted_tfbs, tf_target_pass, gene_tss, regulatory_prior = NULL, max_distance_bp = 100000, id_offset = 0L) {
  empty <- .module2_empty_candidates()
  if (!nrow(predicted_tfbs) || !nrow(tf_target_pass)) {
    return(empty)
  }
  pred <- data.table::as.data.table(predicted_tfbs[, c("fp_id", "tf", "chr", "start", "end", "atac_peak"), drop = FALSE])
  pred <- unique(pred[!is.na(fp_id) & nzchar(fp_id) & !is.na(tf) & nzchar(tf) & !is.na(chr) & is.finite(start) & is.finite(end)])
  if (!nrow(pred)) return(empty)
  pred[, fp_center := as.integer(floor((as.integer(start) + as.integer(end)) / 2))]
  pred[, point_start := fp_center]
  pred[, point_end := fp_center]

  pass <- data.table::as.data.table(tf_target_pass[, c("tf", "target_gene"), drop = FALSE])
  pass <- unique(pass[!is.na(tf) & nzchar(tf) & !is.na(target_gene) & nzchar(target_gene)])
  if (!nrow(pass)) return(empty)

  gt_all <- data.table::as.data.table(gene_tss)
  gt_all <- unique(gt_all[!is.na(target_gene) & nzchar(target_gene) & !is.na(target_chr) & is.finite(target_tss), .(target_gene, target_chr, target_tss, target_strand)])
  cand_parts <- list()

  gt <- data.table::copy(gt_all[target_gene %in% pass$target_gene])
  if (nrow(gt)) {
    data.table::setnames(gt, c("target_chr"), c("chr"), skip_absent = TRUE)
    gt[, window_start := pmax(0L, as.integer(target_tss - as.numeric(max_distance_bp)))]
    gt[, window_end := as.integer(target_tss + as.numeric(max_distance_bp))]
    pass_gt <- pass[gt, on = "target_gene", allow.cartesian = TRUE, nomatch = 0L]

    if (nrow(pass_gt)) {
      fps <- pred[, .(fp_id, tf, chr, point_start, point_end, start, end, atac_peak)]
      wins <- pass_gt[, .(tf, target_gene, chr, window_start, window_end, target_tss, target_strand)]
      data.table::setkey(fps, tf, chr, point_start, point_end)
      data.table::setkey(wins, tf, chr, window_start, window_end)
      hit <- data.table::foverlaps(
        fps,
        wins,
        by.x = c("tf", "chr", "point_start", "point_end"),
        by.y = c("tf", "chr", "window_start", "window_end"),
        type = "within",
        nomatch = 0L
      )

      if (nrow(hit)) {
        hit[, distance_to_tss := as.numeric(point_start) - as.numeric(target_tss)]
        hit[target_strand == "-", distance_to_tss := -distance_to_tss]
        hit[, within_tss_window := TRUE]
        hit[, prior_supported := FALSE]
        hit[, `:=`(prior_id = NA_character_, prior_source = NA_character_, prior_score = NA_real_, prior_status = NA_character_)]
        hit[, candidate_source := "tss_window"]
        cand_parts[[length(cand_parts) + 1L]] <- unique(hit[, .(fp_id, target_gene, chr, start, end, atac_peak, target_chr = chr, target_tss, target_strand, distance_to_tss, candidate_source, within_tss_window, prior_supported, prior_id, prior_source, prior_score, prior_status)])
      }
    }
  }

  prior <- .module2_normalize_prior(regulatory_prior, predicted_tfbs = predicted_tfbs)
  if (nrow(prior)) {
    fp_meta <- unique(pred[, .(fp_id, tf, chr, start, end, atac_peak, fp_center)])
    prior_dt <- data.table::as.data.table(prior)
    prior_dt <- fp_meta[prior_dt, on = "fp_id", allow.cartesian = TRUE, nomatch = 0L]
    prior_dt <- pass[prior_dt, on = c("tf", "target_gene"), nomatch = 0L]
    gt_prior <- data.table::copy(gt_all)
    prior_dt <- gt_prior[prior_dt, on = "target_gene", nomatch = 0L]
    if (nrow(prior_dt)) {
      d_genomic <- as.numeric(prior_dt$fp_center) - as.numeric(prior_dt$target_tss)
      prior_dt[, distance_to_tss := ifelse(as.character(target_strand) == "-", -d_genomic, d_genomic)]
      prior_dt[, `:=`(
        candidate_source = "regulatory_prior",
        within_tss_window = is.finite(distance_to_tss) & abs(distance_to_tss) <= as.numeric(max_distance_bp),
        prior_supported = TRUE
      )]
      cand_parts[[length(cand_parts) + 1L]] <- unique(prior_dt[, .(fp_id, target_gene, chr, start, end, atac_peak, target_chr, target_tss, target_strand, distance_to_tss, candidate_source, within_tss_window, prior_supported, prior_id, prior_source, prior_score, prior_status)])
    }
  }

  if (!length(cand_parts)) return(empty)
  cand <- data.table::rbindlist(cand_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(cand)) return(empty)
  cand$prior_supported_sort <- as.integer(cand$prior_supported %in% TRUE)
  cand$within_tss_window_sort <- as.integer(cand$within_tss_window %in% TRUE)
  cand <- cand[order(
    cand$fp_id,
    cand$target_gene,
    -cand$prior_supported_sort,
    -cand$within_tss_window_sort
  )]
  cand$prior_supported_sort <- NULL
  cand$within_tss_window_sort <- NULL
  cand <- data.table::as.data.table(cand)
  cand <- cand[, {
    within_any <- any(.SD$within_tss_window %in% TRUE)
    prior_any <- any(.SD$prior_supported %in% TRUE)
    x <- .SD[1L]
    x[, `:=`(
      within_tss_window = within_any,
      prior_supported = prior_any
    )]
    x[, candidate_source := fifelse(
      within_tss_window %in% TRUE & prior_supported %in% TRUE,
      "both",
      fifelse(within_tss_window %in% TRUE, "tss_window", "regulatory_prior")
    )]
    x
  }, by = .(fp_id, target_gene)]
  cand[, candidate_id := sprintf("cand_%08d", as.integer(id_offset) + seq_len(.N))]
  data.table::setcolorder(cand, c("candidate_id", setdiff(names(cand), "candidate_id")))
  tibble::as_tibble(cand)
}

.module2_build_candidates_from_fp_meta <- function(fp_meta, tf_target_pass, gene_tss, regulatory_prior = NULL, max_distance_bp = 100000, id_offset = 0L) {
  atac_peak <- candidate_source <- chr <- distance_to_tss <- end <- fp_center <- fp_id <- point_end <- point_start <- prior_id <- prior_score <- prior_source <- prior_status <- prior_supported <- start <- target_gene <- target_strand <- target_tss <- within_tss_window <- NULL
  empty <- .module2_empty_candidates()
  if (!nrow(fp_meta) || !nrow(tf_target_pass)) {
    return(empty)
  }
  fp <- data.table::as.data.table(fp_meta[, c("fp_id", "chr", "start", "end", "atac_peak"), drop = FALSE])
  fp <- unique(fp[!is.na(fp_id) & nzchar(fp_id) & !is.na(chr) & is.finite(start) & is.finite(end)])
  if (!nrow(fp)) return(empty)
  fp[, fp_center := as.integer(floor((as.integer(start) + as.integer(end)) / 2))]
  fp[, point_start := fp_center]
  fp[, point_end := fp_center]

  pass <- data.table::as.data.table(tf_target_pass[, c("tf", "target_gene"), drop = FALSE])
  pass <- unique(pass[!is.na(target_gene) & nzchar(target_gene)])
  target_keep <- unique(as.character(pass$target_gene))
  if (!length(target_keep)) return(empty)

  gt_all <- data.table::as.data.table(gene_tss)
  gt_all <- unique(gt_all[!is.na(target_gene) & nzchar(target_gene) & !is.na(target_chr) & is.finite(target_tss), .(target_gene, target_chr, target_tss, target_strand)])
  cand_parts <- list()

  gt <- data.table::copy(gt_all[target_gene %in% target_keep])
  if (nrow(gt)) {
    data.table::setnames(gt, c("target_chr"), c("chr"), skip_absent = TRUE)
    gt[, window_start := pmax(0L, as.integer(target_tss - as.numeric(max_distance_bp)))]
    gt[, window_end := as.integer(target_tss + as.numeric(max_distance_bp))]
    fps <- unique(fp[, .(fp_id, chr, point_start, point_end, start, end, atac_peak)])
    wins <- gt[, .(target_gene, chr, window_start, window_end, target_tss, target_strand)]
    data.table::setkey(fps, chr, point_start, point_end)
    data.table::setkey(wins, chr, window_start, window_end)
    hit <- data.table::foverlaps(
      fps,
      wins,
      by.x = c("chr", "point_start", "point_end"),
      by.y = c("chr", "window_start", "window_end"),
      type = "within",
      nomatch = 0L
    )

    if (nrow(hit)) {
      hit[, distance_to_tss := as.numeric(point_start) - as.numeric(target_tss)]
      hit[target_strand == "-", distance_to_tss := -distance_to_tss]
      hit[, within_tss_window := TRUE]
      hit[, prior_supported := FALSE]
      hit[, `:=`(prior_id = NA_character_, prior_source = NA_character_, prior_score = NA_real_, prior_status = NA_character_)]
      hit[, candidate_source := "tss_window"]
      cand_parts[[length(cand_parts) + 1L]] <- unique(hit[, .(fp_id, target_gene, chr, start, end, atac_peak, target_chr = chr, target_tss, target_strand, distance_to_tss, candidate_source, within_tss_window, prior_supported, prior_id, prior_source, prior_score, prior_status)])
    }
  }

  prior <- .module2_normalize_prior(regulatory_prior, predicted_tfbs = fp_meta)
  if (nrow(prior)) {
    fp_meta_dt <- unique(fp[, .(fp_id, chr, start, end, atac_peak, fp_center)])
    prior_dt <- data.table::as.data.table(prior)
    prior_dt <- prior_dt[target_gene %in% target_keep]
    prior_dt <- fp_meta_dt[prior_dt, on = "fp_id", allow.cartesian = TRUE, nomatch = 0L]
    gt_prior <- data.table::copy(gt_all)
    prior_dt <- gt_prior[prior_dt, on = "target_gene", nomatch = 0L]
    if (nrow(prior_dt)) {
      d_genomic <- as.numeric(prior_dt$fp_center) - as.numeric(prior_dt$target_tss)
      prior_dt[, distance_to_tss := ifelse(as.character(target_strand) == "-", -d_genomic, d_genomic)]
      prior_dt[, `:=`(
        candidate_source = "regulatory_prior",
        within_tss_window = is.finite(distance_to_tss) & abs(distance_to_tss) <= as.numeric(max_distance_bp),
        prior_supported = TRUE
      )]
      cand_parts[[length(cand_parts) + 1L]] <- unique(prior_dt[, .(fp_id, target_gene, chr, start, end, atac_peak, target_chr, target_tss, target_strand, distance_to_tss, candidate_source, within_tss_window, prior_supported, prior_id, prior_source, prior_score, prior_status)])
    }
  }

  if (!length(cand_parts)) return(empty)
  cand <- data.table::rbindlist(cand_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(cand)) return(empty)
  cand$prior_supported_sort <- as.integer(cand$prior_supported %in% TRUE)
  cand$within_tss_window_sort <- as.integer(cand$within_tss_window %in% TRUE)
  cand <- cand[order(
    cand$fp_id,
    cand$target_gene,
    -cand$prior_supported_sort,
    -cand$within_tss_window_sort
  )]
  cand$prior_supported_sort <- NULL
  cand$within_tss_window_sort <- NULL
  cand <- data.table::as.data.table(cand)
  cand <- cand[, {
    within_any <- any(.SD$within_tss_window %in% TRUE)
    prior_any <- any(.SD$prior_supported %in% TRUE)
    x <- .SD[1L]
    x[, `:=`(
      within_tss_window = within_any,
      prior_supported = prior_any
    )]
    x[, candidate_source := fifelse(
      within_tss_window %in% TRUE & prior_supported %in% TRUE,
      "both",
      fifelse(within_tss_window %in% TRUE, "tss_window", "regulatory_prior")
    )]
    x
  }, by = .(fp_id, target_gene)]
  cand[, candidate_id := sprintf("cand_%08d", as.integer(id_offset) + seq_len(.N))]
  data.table::setcolorder(cand, c("candidate_id", setdiff(names(cand), "candidate_id")))
  tibble::as_tibble(cand)
}

.module2_predicted_manifest <- function(predicted_tfbs) {
  if (!is.character(predicted_tfbs) || length(predicted_tfbs) != 1L || !file.exists(predicted_tfbs)) return(NULL)
  if (!grepl("_manifest\\.csv$", basename(predicted_tfbs))) return(NULL)
  man <- tibble::as_tibble(data.table::fread(predicted_tfbs, showProgress = FALSE))
  if (!all(c("path", "format") %in% names(man))) return(NULL)
  attr(man, "manifest_path") <- predicted_tfbs
  man
}

.module2_read_predicted_chunk <- function(path, format = NULL, columns = NULL) {
  if (identical(format, "parquet") || grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet predicted TFBS chunks.")
    if (is.null(columns)) return(tibble::as_tibble(arrow::read_parquet(path)))
    return(tibble::as_tibble(arrow::read_parquet(path, col_select = columns)))
  }
  if (is.null(columns)) return(tibble::as_tibble(readr::read_csv(path, show_col_types = FALSE)))
  tibble::as_tibble(readr::read_csv(path, col_select = dplyr::all_of(columns), show_col_types = FALSE))
}

.module2_write_chunk <- function(x, out_dir, prefix, chunk_id, output_format = c("auto", "parquet", "csv")) {
  fmt <- .predicted_tfbs_output_format(output_format)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(fmt, "parquet") && requireNamespace("arrow", quietly = TRUE)) {
    path <- file.path(out_dir, sprintf("%s_%04d.parquet", prefix, as.integer(chunk_id)))
    .write_parquet_table(x, path)
  } else {
    fmt <- "csv"
    path <- file.path(out_dir, sprintf("%s_%04d.csv.gz", prefix, as.integer(chunk_id)))
    readr::write_csv(x, path)
  }
  tibble::tibble(chunk_id = as.integer(chunk_id), path = path, format = fmt, n_rows = nrow(x))
}

.module2_link_tf_targets_streamed <- function(multiomic_data, predicted_manifest, gene_tss, regulatory_prior = NULL, project_config = NULL, output_dir, max_distance_bp = NULL, n_cores = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE, write_qc_report = TRUE, qc_report_validate = FALSE) {
  output_format <- match.arg(output_format)
  stopifnot(is.logical(write_qc_report), length(write_qc_report) == 1L, !is.na(write_qc_report))
  stopifnot(is.logical(qc_report_validate), length(qc_report_validate) == 1L, !is.na(qc_report_validate))
  cfg <- .module2_cfg(project_config)
  if (is.null(max_distance_bp)) max_distance_bp <- as.numeric(.module2_cfg_value(cfg, "max_distance_bp", .module2_cfg_value(cfg, "link_window_bp", 100000)))[[1L]]
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  validate_multiomic_object(multiomic_data)
  gene_tss <- .module2_resolve_gene_tss(gene_tss, cfg = cfg, multiomic_data = multiomic_data, verbose = verbose)
  mats <- multiomic_data$matrices
  gene_on <- mats$gene_on
  gene_expr <- mats$gene_expr
  fp_bound <- mats$fp_bound
  fp_score <- mats$fp_score
  expressed_genes <- rownames(gene_on)[rowSums(gene_on > 0, na.rm = TRUE) > 0]
  bound_fps <- rownames(fp_bound)[rowSums(fp_bound > 0, na.rm = TRUE) > 0]
  target_genes <- intersect(expressed_genes, gene_tss$target_gene)

  tfs <- character()
  n_predicted_tfbs <- 0
  fp_meta_chunks <- vector("list", nrow(predicted_manifest))
  for (i in seq_len(nrow(predicted_manifest))) {
    pred_i <- .module2_read_predicted_chunk(as.character(predicted_manifest$path[[i]]), as.character(predicted_manifest$format[[i]]), columns = c("fp_id", "tf", "chr", "start", "end", "atac_peak"))
    pred_i <- pred_i[pred_i$fp_id %in% bound_fps & pred_i$tf %in% expressed_genes, , drop = FALSE]
    n_predicted_tfbs <- n_predicted_tfbs + nrow(pred_i)
    tfs <- union(tfs, as.character(pred_i$tf))
    if (nrow(pred_i)) {
      fp_meta_chunks[[i]] <- unique(pred_i[, c("fp_id", "chr", "start", "end", "atac_peak"), drop = FALSE])
    }
    rm(pred_i)
  }
  fp_meta_chunks <- fp_meta_chunks[!vapply(fp_meta_chunks, is.null, logical(1L))]
  fp_meta <- if (length(fp_meta_chunks)) {
    tibble::as_tibble(unique(data.table::rbindlist(fp_meta_chunks, use.names = TRUE, fill = TRUE)))
  } else {
    tibble::tibble(fp_id = character(), chr = character(), start = integer(), end = integer(), atac_peak = character())
  }
  tfs <- sort(tfs)
  if (isTRUE(verbose)) .log_inform("Module 2 inputs: {length(tfs)} TF(s), {length(target_genes)} target gene(s), {n_predicted_tfbs} predicted TFBS row(s), {nrow(fp_meta)} unique predicted FP(s), streamed from {nrow(predicted_manifest)} chunk(s).")
  tf_pairs <- tidyr::crossing(tf = tfs, target_gene = target_genes)
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: testing {nrow(tf_pairs)} pair(s).")
  tf_cutoffs <- .module2_corr_cutoffs(cfg, "tf_target", r_default = .module2_cfg_value(cfg, "threshold_rna_gene_corr_r", 0.3))
  if (is.null(tf_cutoffs$p)) tf_cutoffs$p <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_p", NULL)
  if (is.null(tf_cutoffs$fdr)) tf_cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_fdr", NULL)
  tf_target_corr <- .module2_pair_correlations(gene_expr, gene_expr, tf_pairs, "tf", "target_gene", tf_cutoffs, n_cores = n_cores)
  tf_target_pass <- tf_target_corr[tf_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: {nrow(tf_target_pass)} pair(s) passed.")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data_dir <- file.path(output_dir, "data")
  cand_dir <- file.path(data_dir, "candidates")
  corr_dir <- file.path(data_dir, "correlations")
  link_dir <- file.path(data_dir, "links")
  dir.create(file.path(output_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
  candidates <- .module2_build_candidates_from_fp_meta(fp_meta, tf_target_pass, gene_tss, regulatory_prior = regulatory_prior, max_distance_bp = max_distance_bp)
  cand_table <- .module2_manifest_table(.module2_write_table(candidates, cand_dir, "fp_target_candidates", output_format), "module2_fp_target_candidates")
  fp_pair_dt <- if (nrow(candidates)) {
    unique(data.table::as.data.table(candidates[, c("fp_id", "target_gene"), drop = FALSE]))
  } else {
    data.table::data.table(fp_id = character(), target_gene = character())
  }
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target candidates: {nrow(candidates)} unique candidate row(s) from {nrow(fp_meta)} predicted FP(s).")
  fp_pairs <- tibble::as_tibble(fp_pair_dt)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: testing {nrow(fp_pairs)} unique restricted pair(s).")
  fp_cutoffs <- .module2_corr_cutoffs(cfg, "fp_target", r_default = .module2_cfg_value(cfg, "threshold_fp_gene_corr_r", 0.3))
  if (is.null(fp_cutoffs$p)) fp_cutoffs$p <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_p", NULL)
  if (is.null(fp_cutoffs$fdr)) fp_cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_fdr", NULL)
  fp_target_corr <- .module2_pair_correlations(fp_score, gene_expr, fp_pairs, "fp_id", "target_gene", fp_cutoffs, n_cores = n_cores)
  fp_target_pass <- fp_target_corr[fp_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: {nrow(fp_target_pass)} pair(s) passed.")
  corr_table <- .module2_manifest_table(.module2_write_table(fp_target_corr, corr_dir, "fp_target_corr", output_format), "module2_fp_target_corr")

  fp_pass_dt <- data.table::as.data.table(fp_target_pass[, c("fp_id", "target_gene"), drop = FALSE])
  data.table::setkey(fp_pass_dt, fp_id, target_gene)
  tf_pass_dt <- unique(data.table::as.data.table(tf_target_pass[, c("tf", "target_gene"), drop = FALSE]))
  data.table::setkey(tf_pass_dt, tf, target_gene)
  cand_all <- candidates[, c("candidate_id", "fp_id", "target_gene"), drop = FALSE]
  cand_all_dt <- unique(data.table::as.data.table(cand_all[, c("candidate_id", "fp_id", "target_gene"), drop = FALSE]))
  data.table::setkey(cand_all_dt, fp_id)
  link_manifest <- list()
  n_links <- 0L
  link_offset <- 0L
  for (i in seq_len(nrow(predicted_manifest))) {
    pred_i <- .module2_read_predicted_chunk(as.character(predicted_manifest$path[[i]]), as.character(predicted_manifest$format[[i]]), columns = c("fp_id", "tf"))
    pred_i <- pred_i[pred_i$fp_id %in% bound_fps & pred_i$tf %in% tfs, , drop = FALSE]
    if (nrow(cand_all_dt) && nrow(pred_i)) {
      pred_dt <- unique(data.table::as.data.table(pred_i[, c("fp_id", "tf"), drop = FALSE]))
      cand_dt <- cand_all_dt[unique(pred_dt$fp_id), on = "fp_id", nomatch = 0L]
      link_i <- pred_dt[cand_dt, on = "fp_id", allow.cartesian = TRUE, nomatch = 0L]
      link_i <- tf_pass_dt[link_i, on = c("tf", "target_gene"), nomatch = 0L]
      link_i <- fp_pass_dt[link_i, on = c("fp_id", "target_gene"), nomatch = 0L]
      link_i <- unique(link_i)
      if (nrow(link_i)) {
        link_i[, link_id := sprintf("link_%08d", as.integer(link_offset) + seq_len(.N))]
        link_offset <- link_offset + nrow(link_i)
        link_i[, `:=`(tf_target_pass = TRUE, fp_target_pass = TRUE, module2_link_pass = TRUE)]
        link_i <- link_i[, .(link_id, tf, fp_id, target_gene, candidate_id, tf_target_pass, fp_target_pass, module2_link_pass)]
      }
    } else {
      link_i <- data.table::data.table(link_id = character(), tf = character(), fp_id = character(), target_gene = character(), candidate_id = character(), tf_target_pass = logical(), fp_target_pass = logical(), module2_link_pass = logical())
    }
    n_links <- n_links + nrow(link_i)
    link_manifest[[i]] <- .module2_write_chunk(tibble::as_tibble(link_i), link_dir, "module2_links", i, output_format)
    if (isTRUE(verbose)) .log_inform("Module 2 link chunks: {i}/{nrow(predicted_manifest)} done, {nrow(link_i)} link row(s) in this chunk.")
    rm(pred_i, link_i)
  }
  link_manifest <- dplyr::bind_rows(link_manifest)
  link_manifest_path <- file.path(link_dir, "module2_links_manifest.csv")
  readr::write_csv(link_manifest, link_manifest_path)
  tf_target_manifest <- .module2_manifest_table(.module2_write_table(tf_target_corr, corr_dir, "tf_target_corr", output_format), "module2_tf_target_corr")
  qc_summary <- tibble::tibble(metric = c("n_predicted_tfbs", "n_tfs", "n_target_genes", "n_tf_target_pairs_tested", "n_tf_target_pairs_pass", "n_fp_target_candidates", "n_fp_target_pairs_tested", "n_fp_target_pairs_pass", "n_module2_links", "n_active_link_conditions"), value = c(n_predicted_tfbs, length(tfs), length(target_genes), nrow(tf_pairs), nrow(tf_target_pass), nrow(candidates), nrow(fp_pairs), nrow(fp_target_pass), n_links, NA_real_))
  qc_manifest <- .module2_write_run_summary(qc_summary, output_dir)
  predicted_manifest_path <- attr(predicted_manifest, "manifest_path")
  if (is.null(predicted_manifest_path) || !nzchar(predicted_manifest_path)) predicted_manifest_path <- file.path(output_dir, "module1_predicted_tfbs_manifest.csv")
  manifest <- dplyr::bind_rows(
    tibble::tibble(table = "module1_predicted_tfbs", path = predicted_manifest_path, format = "manifest", n_rows = n_predicted_tfbs),
    cand_table,
    tf_target_manifest,
    corr_table,
    tibble::tibble(table = "module2_links", path = link_manifest_path, format = "manifest", n_rows = n_links),
    qc_manifest
  )
  manifest_path <- file.path(output_dir, "module2_manifest.csv")
  readr::write_csv(manifest, manifest_path)
  out <- list(predicted_tfbs = tibble::tibble(), candidates = tibble::tibble(), tf_target_corr = tf_target_corr, fp_target_corr = tibble::tibble(), links = tibble::tibble(), module2_links_manifest = link_manifest, condition_activity = tibble::tibble(), qc_summary = qc_summary, manifest = manifest, reports = list(manifest = manifest_path, links_manifest = link_manifest_path), parameters = list(max_distance_bp = max_distance_bp, n_cores = .module2_default_cores(n_cores), output_format = output_format, streamed = TRUE, write_qc_report = isTRUE(write_qc_report), qc_report_validate = isTRUE(qc_report_validate)))
  class(out) <- c("craftgrn_module2", "list")
  if (isTRUE(write_qc_report)) {
    out$reports$qc_html <- build_module2_qc_report(
      module2 = out,
      multiomic_data = multiomic_data,
      output_dir = output_dir,
      scan_large_tables = FALSE,
      validate_integrity = isTRUE(qc_report_validate),
      verbose = verbose
    )
  }
  out
}

#' Link predicted TF binding sites to target genes
#'
#' @param multiomic_data CraftGRN multiomic object returned by `load_prep_multiomic_data()`.
#' @param predicted_tfbs Compact Module 1 predicted TFBS table or path.
#' @param gene_tss Optional gene TSS annotation table or path. If `NULL`, the
#'   table is resolved from `project_config$gene_tss` or generated from the
#'   configured `ref_genome`.
#' @param regulatory_prior Optional generic FP-target regulatory prior.
#' @param project_config Optional project YAML path or list.
#' @param output_dir Optional output directory.
#' @param max_distance_bp Maximum signed distance to TSS for window candidates.
#' @param n_cores Number of CPU cores.
#' @param output_format Output format: auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#' @param write_qc_report Write a Module 2 HTML QC report when `output_dir` is supplied.
#' @param qc_report_validate Run relational integrity checks in the automatic QC report.
#' @return Compact Module 2 relational result list.
#' @noRd
link_tf_targets <- function(multiomic_data, predicted_tfbs, gene_tss = NULL, regulatory_prior = NULL, project_config = NULL, output_dir = NULL, max_distance_bp = NULL, n_cores = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE, write_qc_report = TRUE, qc_report_validate = FALSE) {
  output_format <- match.arg(output_format)
  stopifnot(is.logical(write_qc_report), length(write_qc_report) == 1L, !is.na(write_qc_report))
  stopifnot(is.logical(qc_report_validate), length(qc_report_validate) == 1L, !is.na(qc_report_validate))
  cfg <- .module2_cfg(project_config)
  if (is.null(max_distance_bp)) max_distance_bp <- as.numeric(.module2_cfg_value(cfg, "max_distance_bp", .module2_cfg_value(cfg, "link_window_bp", 100000)))[[1L]]
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object returned by load_prep_multiomic_data().")
  validate_multiomic_object(multiomic_data)
  predicted_manifest <- .module2_predicted_manifest(predicted_tfbs)
  if (!is.null(predicted_manifest)) {
    if (is.null(output_dir) || !nzchar(output_dir)) .log_abort("output_dir is required for streamed Module 2 runs.")
    return(.module2_link_tf_targets_streamed(
      multiomic_data = multiomic_data,
      predicted_manifest = predicted_manifest,
      gene_tss = gene_tss,
      regulatory_prior = regulatory_prior,
      project_config = project_config,
      output_dir = output_dir,
      max_distance_bp = max_distance_bp,
      n_cores = n_cores,
      output_format = output_format,
      verbose = verbose,
      write_qc_report = isTRUE(write_qc_report),
      qc_report_validate = isTRUE(qc_report_validate)
    ))
  }
  if (is.character(predicted_tfbs) && length(predicted_tfbs) == 1L && file.exists(predicted_tfbs)) predicted_tfbs <- load_predicted_tfbs(predicted_tfbs)
  predicted_tfbs <- build_predicted_tfbs(predicted_tfbs)
  gene_tss <- .module2_resolve_gene_tss(gene_tss, cfg = cfg, multiomic_data = multiomic_data, verbose = verbose)
  mats <- multiomic_data$matrices
  gene_on <- mats$gene_on
  gene_expr <- mats$gene_expr
  fp_bound <- mats$fp_bound
  fp_score <- mats$fp_score
  expressed_genes <- rownames(gene_on)[rowSums(gene_on > 0, na.rm = TRUE) > 0]
  bound_fps <- rownames(fp_bound)[rowSums(fp_bound > 0, na.rm = TRUE) > 0]
  predicted_tfbs <- predicted_tfbs[predicted_tfbs$fp_id %in% bound_fps & predicted_tfbs$tf %in% expressed_genes, , drop = FALSE]
  target_genes <- intersect(expressed_genes, gene_tss$target_gene)
  tfs <- sort(unique(as.character(predicted_tfbs$tf)))
  if (isTRUE(verbose)) .log_inform("Module 2 inputs: {length(tfs)} TF(s), {length(target_genes)} target gene(s), {nrow(predicted_tfbs)} predicted TFBS row(s).")
  tf_pairs <- tidyr::crossing(tf = tfs, target_gene = target_genes)
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: testing {nrow(tf_pairs)} pair(s).")
  tf_cutoffs <- .module2_corr_cutoffs(cfg, "tf_target", r_default = .module2_cfg_value(cfg, "threshold_rna_gene_corr_r", 0.3))
  if (is.null(tf_cutoffs$p)) tf_cutoffs$p <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_p", NULL)
  if (is.null(tf_cutoffs$fdr)) tf_cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_rna_gene_corr_fdr", NULL)
  tf_target_corr <- .module2_pair_correlations(gene_expr, gene_expr, tf_pairs, "tf", "target_gene", tf_cutoffs, n_cores = n_cores)
  tf_target_pass <- tf_target_corr[tf_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: {nrow(tf_target_pass)} pair(s) passed.")
  candidates <- .module2_build_candidates(predicted_tfbs, tf_target_pass, gene_tss, regulatory_prior = regulatory_prior, max_distance_bp = max_distance_bp)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target candidates after TF-target and TSS/prior filters: {nrow(candidates)} pair(s).")
  fp_pairs <- unique(candidates[, c("fp_id", "target_gene"), drop = FALSE])
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: testing {nrow(fp_pairs)} restricted pair(s).")
  fp_cutoffs <- .module2_corr_cutoffs(cfg, "fp_target", r_default = .module2_cfg_value(cfg, "threshold_fp_gene_corr_r", 0.3))
  if (is.null(fp_cutoffs$p)) fp_cutoffs$p <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_p", NULL)
  if (is.null(fp_cutoffs$fdr)) fp_cutoffs$fdr <- .module2_cfg_value(cfg, "threshold_fp_gene_corr_fdr", NULL)
  fp_target_corr <- .module2_pair_correlations(fp_score, gene_expr, fp_pairs, "fp_id", "target_gene", fp_cutoffs, n_cores = n_cores)
  fp_target_pass <- fp_target_corr[fp_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: {nrow(fp_target_pass)} pair(s) passed.")
  pred_dt <- unique(data.table::as.data.table(predicted_tfbs[, c("fp_id", "tf"), drop = FALSE]))
  cand_dt <- unique(data.table::as.data.table(candidates[, c("candidate_id", "fp_id", "target_gene"), drop = FALSE]))
  tf_pass_dt <- unique(data.table::as.data.table(tf_target_pass[, c("tf", "target_gene"), drop = FALSE]))
  fp_pass_dt <- unique(data.table::as.data.table(fp_target_pass[, c("fp_id", "target_gene"), drop = FALSE]))
  links <- pred_dt[cand_dt, on = "fp_id", allow.cartesian = TRUE, nomatch = 0L]
  links <- tf_pass_dt[links, on = c("tf", "target_gene"), nomatch = 0L]
  links <- fp_pass_dt[links, on = c("fp_id", "target_gene"), nomatch = 0L]
  links <- unique(links)
  links$link_id <- sprintf("link_%08d", seq_len(nrow(links)))
  links$tf_target_pass <- TRUE
  links$fp_target_pass <- TRUE
  links$module2_link_pass <- TRUE
  links <- tibble::as_tibble(links[, c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "tf_target_pass", "fp_target_pass", "module2_link_pass"), drop = FALSE])
  activity <- .module2_condition_activity(links, predicted_tfbs, multiomic_data)
  qc_summary <- tibble::tibble(metric = c("n_predicted_tfbs", "n_tfs", "n_target_genes", "n_tf_target_pairs_tested", "n_tf_target_pairs_pass", "n_fp_target_candidates", "n_fp_target_pairs_tested", "n_fp_target_pairs_pass", "n_module2_links", "n_active_link_conditions"), value = c(nrow(predicted_tfbs), length(tfs), length(target_genes), nrow(tf_pairs), nrow(tf_target_pass), nrow(candidates), nrow(fp_pairs), nrow(fp_target_pass), nrow(links), sum(activity$active %in% TRUE)))
  reports <- list()
  manifest <- tibble::tibble()
  if (!is.null(output_dir) && nzchar(output_dir)) {
    data_dir <- file.path(output_dir, "data")
    dir.create(file.path(output_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
    manifest <- dplyr::bind_rows(
      .module2_manifest_table(.module2_write_table(predicted_tfbs, file.path(data_dir, "predicted_tfbs"), "predicted_tfbs", output_format), "module1_predicted_tfbs"),
      .module2_manifest_table(.module2_write_table(candidates, file.path(data_dir, "candidates"), "fp_target_candidates", output_format), "module2_fp_target_candidates"),
      .module2_manifest_table(.module2_write_table(tf_target_corr, file.path(data_dir, "correlations"), "tf_target_corr", output_format), "module2_tf_target_corr"),
      .module2_manifest_table(.module2_write_table(fp_target_corr, file.path(data_dir, "correlations"), "fp_target_corr", output_format), "module2_fp_target_corr"),
      .module2_manifest_table(.module2_write_table(links, file.path(data_dir, "links"), "module2_links", output_format), "module2_links"),
      .module2_manifest_table(.module2_write_table(activity, file.path(data_dir, "condition_activity"), "condition_activity", output_format), "module2_condition_activity"),
      .module2_write_run_summary(qc_summary, output_dir)
    )
    manifest_path <- file.path(output_dir, "module2_manifest.csv")
    readr::write_csv(manifest, manifest_path)
    reports$manifest <- manifest_path
  }
  out <- list(predicted_tfbs = predicted_tfbs, candidates = candidates, tf_target_corr = tf_target_corr, fp_target_corr = fp_target_corr, links = links, condition_activity = activity, qc_summary = qc_summary, manifest = manifest, reports = reports, parameters = list(max_distance_bp = max_distance_bp, n_cores = .module2_default_cores(n_cores), output_format = output_format, write_qc_report = isTRUE(write_qc_report), qc_report_validate = isTRUE(qc_report_validate)))
  class(out) <- c("craftgrn_module2", "list")
  if (!is.null(output_dir) && nzchar(output_dir) && isTRUE(write_qc_report)) {
    out$reports$qc_html <- build_module2_qc_report(
      module2 = out,
      multiomic_data = multiomic_data,
      output_dir = output_dir,
      scan_large_tables = FALSE,
      validate_integrity = isTRUE(qc_report_validate),
      verbose = verbose
    )
  }
  out
}

#' Predict TF targets through TFBS-target and TF-target correlations
#'
#' @param multiomic_data CraftGRN multiomic object returned by `load_prep_multiomic_data()`.
#' @param predicted_tfbs Compact Module 1 predicted TFBS table or manifest path.
#' @param gene_tss Optional gene TSS annotation table or path. If `NULL`, the
#'   table is resolved from `project_config$gene_tss` or generated from the
#'   configured `ref_genome`.
#' @param regulatory_prior Optional generic FP-target regulatory prior.
#' @param project_config Optional project YAML path or list.
#' @param output_dir Optional output directory.
#' @param max_distance_bp Maximum signed distance to TSS for window candidates.
#' @param n_cores Number of CPU cores.
#' @param output_format Output format: auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#' @param write_qc_report Write a Module 2 HTML QC report when `output_dir` is supplied.
#' @param qc_report_validate Run relational integrity checks in the automatic QC report.
#' @return Compact Module 2 relational result list.
#' @export
predict_tf_targets <- function(multiomic_data, predicted_tfbs, gene_tss = NULL, regulatory_prior = NULL, project_config = NULL, output_dir = NULL, max_distance_bp = NULL, n_cores = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE, write_qc_report = TRUE, qc_report_validate = FALSE) {
  link_tf_targets(
    multiomic_data = multiomic_data,
    predicted_tfbs = predicted_tfbs,
    gene_tss = gene_tss,
    regulatory_prior = regulatory_prior,
    project_config = project_config,
    output_dir = output_dir,
    max_distance_bp = max_distance_bp,
    n_cores = n_cores,
    output_format = output_format,
    verbose = verbose,
    write_qc_report = write_qc_report,
    qc_report_validate = qc_report_validate
  )
}

.module2_condition_activity <- function(links, predicted_tfbs, multiomic_data) {
  if (!nrow(links)) return(tibble::tibble(link_id = character(), condition = character(), tf_expr = numeric(), target_expr = numeric(), fp_score = numeric(), fp_bound = logical(), atac_score = numeric(), active = logical()))
  mats <- multiomic_data$matrices
  cond <- colnames(mats$fp_score)
  dt <- data.table::as.data.table(links)
  pred <- unique(data.table::as.data.table(predicted_tfbs[, c("fp_id", "atac_peak"), drop = FALSE]))
  dt <- pred[dt, on = "fp_id"]
  rows <- vector("list", length(cond))
  for (i in seq_along(cond)) {
    cc <- cond[[i]]
    tf_expr <- mats$gene_expr[match(dt$tf, rownames(mats$gene_expr)), cc]
    target_expr <- mats$gene_expr[match(dt$target_gene, rownames(mats$gene_expr)), cc]
    fs <- mats$fp_score[match(dt$fp_id, rownames(mats$fp_score)), cc]
    fb <- mats$fp_bound[match(dt$fp_id, rownames(mats$fp_bound)), cc]
    atac_score <- rep(NA_real_, nrow(dt))
    if (is.matrix(mats$atac_score) && "atac_peak" %in% names(dt)) atac_score <- mats$atac_score[match(dt$atac_peak, rownames(mats$atac_score)), cc]
    active <- dt$module2_link_pass %in% TRUE & is.finite(tf_expr) & tf_expr > 0 & is.finite(target_expr) & target_expr > 0 & is.finite(fs) & fs > 0 & fb %in% TRUE
    rows[[i]] <- tibble::tibble(link_id = dt$link_id, condition = cc, tf_expr = as.numeric(tf_expr), target_expr = as.numeric(target_expr), fp_score = as.numeric(fs), fp_bound = as.logical(fb), atac_score = as.numeric(atac_score), active = as.logical(active))
  }
  dplyr::bind_rows(rows)
}

#' Load predicted links from Module 2
#'
#' @param path Module 2 output directory or module2_manifest.csv path.
#' @return A named list of tables.
#' @noRd
load_module2_links <- function(path) {
  if (dir.exists(path)) path <- file.path(path, "module2_manifest.csv")
  if (!file.exists(path)) .log_abort("Module 2 manifest not found: {path}")
  manifest <- readr::read_csv(path, show_col_types = FALSE)
  out <- list()
  for (i in seq_len(nrow(manifest))) {
    p <- as.character(manifest$path[[i]])
    nm <- as.character(manifest$table[[i]])
    fmt <- as.character(manifest$format[[i]])
    if (identical(fmt, "manifest")) {
      out[[paste0(nm, "_manifest")]] <- readr::read_csv(p, show_col_types = FALSE)
    } else if (grepl("\\.parquet$", p, ignore.case = TRUE)) {
      if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet Module 2 outputs.")
      out[[nm]] <- tibble::as_tibble(arrow::read_parquet(p))
    } else {
      out[[nm]] <- tibble::as_tibble(readr::read_csv(p, show_col_types = FALSE))
    }
  }
  out$manifest <- tibble::as_tibble(manifest)
  out
}

#' Load predicted links from Module 2
#'
#' @param path Module 2 output directory or module2_manifest.csv path.
#' @return A named list of Module 2 tables.
#' @export
load_predicted_links <- function(path) {
  load_module2_links(path)
}

#' Query Module 2 links
#'
#' @param module2 Module 2 result list or loaded output list.
#' @param tf Optional TF filter.
#' @param fp_id Optional FP filter.
#' @param target_gene Optional target-gene filter.
#' @param pass_only Keep only passing links.
#' @return A tibble of matching final links.
#' @noRd
query_module2_links <- function(module2, tf = NULL, fp_id = NULL, target_gene = NULL, pass_only = TRUE) {
  links <- if (is.data.frame(module2$links) && ncol(module2$links)) module2$links else module2$module2_links
  link_manifest <- module2$module2_links_manifest
  tf_filter <- if (is.null(tf)) NULL else as.character(tf)
  fp_filter <- if (is.null(fp_id)) NULL else as.character(fp_id)
  target_filter <- if (is.null(target_gene)) NULL else as.character(target_gene)
  filter_dt <- function(dt) {
    if (!is.null(tf_filter)) dt <- dt[tf %in% tf_filter]
    if (!is.null(fp_filter)) dt <- dt[fp_id %in% fp_filter]
    if (!is.null(target_filter)) dt <- dt[target_gene %in% target_filter]
    if (isTRUE(pass_only) && "module2_link_pass" %in% names(dt)) dt <- dt[module2_link_pass %in% TRUE]
    dt
  }
  if (is.data.frame(links) && all(c("path", "format") %in% names(links)) && !all(c("link_id", "tf", "fp_id", "target_gene") %in% names(links))) {
    link_manifest <- links
    links <- NULL
  }
  if (is.data.frame(links)) return(tibble::as_tibble(filter_dt(data.table::as.data.table(links))))
  if (is.data.frame(link_manifest)) {
    rows <- lapply(seq_len(nrow(link_manifest)), function(i) {
      dt <- data.table::as.data.table(.module2_read_predicted_chunk(as.character(link_manifest$path[[i]]), as.character(link_manifest$format[[i]])))
      filter_dt(dt)
    })
    return(tibble::as_tibble(data.table::rbindlist(rows, fill = TRUE)))
  }
  .log_abort("Module 2 links table not found.")
}

#' Query specific links by TF(s) and/or distance to TSS
#'
#' @param module2 Module 2 result list or loaded output list.
#' @param tf Optional TF filter.
#' @param fp_id Optional FP filter.
#' @param target_gene Optional target-gene filter.
#' @param max_distance_to_tss Optional maximum absolute distance to TSS.
#' @param pass_only Keep only passing links.
#' @return A tibble of matching final links.
#' @export
query_predicted_links <- function(module2, tf = NULL, fp_id = NULL, target_gene = NULL, max_distance_to_tss = NULL, pass_only = TRUE) {
  links <- query_module2_links(module2, tf = tf, fp_id = fp_id, target_gene = target_gene, pass_only = pass_only)
  if (!is.null(max_distance_to_tss)) {
    candidates <- if (is.data.frame(module2$candidates) && nrow(module2$candidates)) {
      module2$candidates
    } else if (is.data.frame(module2$module2_fp_target_candidates) && nrow(module2$module2_fp_target_candidates)) {
      module2$module2_fp_target_candidates
    } else {
      .module2_report_read_table(module2, "module2_fp_target_candidates", columns = c("candidate_id", "distance_to_tss"))
    }
    if (!"distance_to_tss" %in% names(candidates)) .log_abort("Module 2 candidate table does not contain distance_to_tss.")
    cand <- candidates[, intersect(c("candidate_id", "distance_to_tss"), names(candidates)), drop = FALSE]
    links <- dplyr::left_join(links, cand, by = "candidate_id")
    links <- links[is.finite(links$distance_to_tss) & abs(links$distance_to_tss) <= as.numeric(max_distance_to_tss)[[1L]], , drop = FALSE]
  }
  tibble::as_tibble(links)
}

#' Validate Module 2 links
#'
#' @param module2 Module 2 result list or loaded output list.
#' @return TRUE invisibly when valid.
#' @noRd
validate_module2_links <- function(module2) {
  links <- if (is.data.frame(module2$links) && ncol(module2$links)) module2$links else module2$module2_links
  if (!is.data.frame(links) && is.data.frame(module2$module2_links_manifest)) links <- module2$module2_links_manifest
  if (!is.data.frame(links)) .log_abort("Module 2 links table not found.")
  if (all(c("path", "format") %in% names(links)) && !all(c("link_id", "tf", "fp_id", "target_gene") %in% names(links))) {
    if (!all(file.exists(links$path))) .log_abort("Module 2 link manifest contains missing chunk files.")
    return(invisible(TRUE))
  }
  need <- c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass")
  missing <- setdiff(need, names(links))
  if (length(missing)) {
    missing_label <- paste(missing, collapse = ", ")
    .log_abort("Module 2 links missing required columns: {missing_label}.")
  }
  invisible(TRUE)
}

#' Perform sanity checks for predicted Module 2 links
#'
#' @param module2 Module 2 result list or loaded output list.
#' @return TRUE invisibly when valid.
#' @noRd
check_module2_links <- function(module2) {
  validate_module2_links(module2)
}

#' Perform sanity check for predicted links for Module 2 diagnostics
#'
#' @param module2 Module 2 result list or loaded output list.
#' @return TRUE invisibly when valid.
#' @export
check_predicted_links <- function(module2) {
  validate_module2_links(module2)
}

#' Export predicted TF-target links as BEDPE
#'
#' @param module2 Module 2 result list or loaded output list.
#' @param output_file BEDPE output file.
#' @param tf Optional TF subset.
#' @return Output path invisibly.
#' @export
export_tf_target_bedpe <- function(module2, output_file, tf = NULL) {
  links <- query_predicted_links(module2, tf = tf, pass_only = TRUE)
  candidates <- if (is.data.frame(module2$candidates)) module2$candidates else module2$module2_fp_target_candidates
  if (!is.data.frame(candidates)) .log_abort("Module 2 candidate table not found.")
  need <- c("candidate_id", "chr", "start", "end", "target_chr", "target_tss")
  if (!all(need %in% names(candidates))) .log_abort("Candidate table is missing BEDPE columns.")
  bed <- dplyr::left_join(links, candidates, by = "candidate_id")
  if (!nrow(bed)) {
    out <- data.frame(chrom1 = character(), start1 = integer(), end1 = integer(), chrom2 = character(), start2 = integer(), end2 = integer(), name = character(), score = integer(), strand1 = character(), strand2 = character())
  } else {
    out <- data.frame(chrom1 = as.character(bed$chr), start1 = as.integer(bed$start), end1 = as.integer(bed$end), chrom2 = as.character(bed$target_chr), start2 = pmax(0L, as.integer(bed$target_tss) - 1L), end2 = as.integer(bed$target_tss), name = paste(as.character(bed$tf), as.character(bed$fp_id), as.character(bed$target_gene), sep = "|"), score = 1000L, strand1 = ".", strand2 = as.character(bed$target_strand), stringsAsFactors = FALSE)
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(out), output_file, sep = "\t", col.names = FALSE)
  invisible(output_file)
}
