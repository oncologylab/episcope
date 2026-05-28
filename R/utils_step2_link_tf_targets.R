# File: utils_step2_link_tf_targets.R
# Purpose: Compact relational Module 2 TF-target linking.

.module2_default_cores <- function(n_cores = NULL) {
  if (is.null(n_cores) || length(n_cores) == 0L || is.na(n_cores[[1L]])) {
    return(max(1L, parallel::detectCores(logical = TRUE) - 1L))
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
  use_p <- x$best_method == "pearson"
  use_s <- x$best_method == "spearman"
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
    arrow::write_parquet(x, path, compression = "zstd")
  } else {
    fmt <- "csv"
    path <- file.path(out_dir, paste0(name, ".csv.gz"))
    readr::write_csv(x, path)
  }
  tibble::tibble(table = name, path = path, format = fmt, n_rows = nrow(x))
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

.module2_build_candidates <- function(predicted_tfbs, tf_target_pass, gene_tss, regulatory_prior = NULL, max_distance_bp = 100000) {
  base <- dplyr::inner_join(predicted_tfbs[, c("fp_id", "tf"), drop = FALSE], tf_target_pass[, c("tf", "target_gene"), drop = FALSE], by = "tf")
  fp_target <- unique(base[, c("fp_id", "target_gene"), drop = FALSE])
  fp_meta <- unique(predicted_tfbs[, c("fp_id", "chr", "start", "end", "atac_peak"), drop = FALSE])
  cand <- dplyr::left_join(fp_target, fp_meta, by = "fp_id")
  cand <- dplyr::left_join(cand, gene_tss, by = "target_gene")
  cand$fp_center <- as.integer(floor((as.integer(cand$start) + as.integer(cand$end)) / 2))
  d_genomic <- as.numeric(cand$fp_center) - as.numeric(cand$target_tss)
  cand$distance_to_tss <- ifelse(as.character(cand$target_strand) == "-", -d_genomic, d_genomic)
  prior <- .module2_normalize_prior(regulatory_prior, predicted_tfbs = predicted_tfbs)
  cand <- dplyr::left_join(cand, prior, by = c("fp_id", "target_gene"))
  cand$within_tss_window <- is.finite(cand$distance_to_tss) & abs(cand$distance_to_tss) <= as.numeric(max_distance_bp)
  cand$prior_supported <- !is.na(cand$prior_id) & nzchar(cand$prior_id)
  cand <- cand[cand$within_tss_window | cand$prior_supported, , drop = FALSE]
  cand$candidate_source <- ifelse(cand$within_tss_window & cand$prior_supported, "both", ifelse(cand$prior_supported, "regulatory_prior", "tss_window"))
  cand$candidate_id <- sprintf("cand_%08d", seq_len(nrow(cand)))
  keep <- c("candidate_id", "fp_id", "target_gene", "chr", "start", "end", "atac_peak", "target_chr", "target_tss", "target_strand", "distance_to_tss", "candidate_source", "within_tss_window", "prior_supported", "prior_id", "prior_source", "prior_score", "prior_status")
  tibble::as_tibble(cand[, intersect(keep, names(cand)), drop = FALSE])
}

#' Link predicted TF binding sites to target genes
#'
#' @param multiomic_data Compact multiomic object or compatible Module 1 object.
#' @param predicted_tfbs Compact Module 1 predicted TFBS table or path.
#' @param gene_tss Gene TSS annotation table.
#' @param regulatory_prior Optional generic FP-target regulatory prior.
#' @param project_config Optional project YAML path or list.
#' @param output_dir Optional output directory.
#' @param max_distance_bp Maximum signed distance to TSS for window candidates.
#' @param n_cores Number of CPU cores.
#' @param output_format Output format: auto, parquet, or csv.
#' @param verbose Emit concise progress messages.
#' @return Compact Module 2 relational result list.
#' @export
link_tf_targets <- function(multiomic_data, predicted_tfbs, gene_tss, regulatory_prior = NULL, project_config = NULL, output_dir = NULL, max_distance_bp = NULL, n_cores = NULL, output_format = c("auto", "parquet", "csv"), verbose = TRUE) {
  output_format <- match.arg(output_format)
  cfg <- .module2_cfg(project_config)
  if (is.null(max_distance_bp)) max_distance_bp <- as.numeric(.module2_cfg_value(cfg, "max_distance_bp", 100000))[[1L]]
  if (!is_multiomic_object(multiomic_data)) multiomic_data <- as_multiomic_object(multiomic_data, verbose = FALSE)
  validate_multiomic_object(multiomic_data)
  if (is.character(predicted_tfbs) && length(predicted_tfbs) == 1L && file.exists(predicted_tfbs)) predicted_tfbs <- load_predicted_tfbs(predicted_tfbs)
  predicted_tfbs <- build_predicted_tfbs(predicted_tfbs)
  gene_tss <- .module2_normalize_gene_tss(gene_tss)
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
  tf_target_corr <- .module2_pair_correlations(gene_expr, gene_expr, tf_pairs, "tf", "target_gene", .module2_corr_cutoffs(cfg, "tf_target"), n_cores = n_cores)
  tf_target_pass <- tf_target_corr[tf_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 TF-target correlation: {nrow(tf_target_pass)} pair(s) passed.")
  candidates <- .module2_build_candidates(predicted_tfbs, tf_target_pass, gene_tss, regulatory_prior = regulatory_prior, max_distance_bp = max_distance_bp)
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target candidates after TF-target and TSS/prior filters: {nrow(candidates)} pair(s).")
  fp_pairs <- unique(candidates[, c("fp_id", "target_gene"), drop = FALSE])
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: testing {nrow(fp_pairs)} restricted pair(s).")
  fp_target_corr <- .module2_pair_correlations(fp_score, gene_expr, fp_pairs, "fp_id", "target_gene", .module2_corr_cutoffs(cfg, "fp_target"), n_cores = n_cores)
  fp_target_pass <- fp_target_corr[fp_target_corr$pass %in% TRUE, , drop = FALSE]
  if (isTRUE(verbose)) .log_inform("Module 2 FP-target correlation: {nrow(fp_target_pass)} pair(s) passed.")
  links <- dplyr::inner_join(predicted_tfbs[, c("fp_id", "tf"), drop = FALSE], tf_target_pass[, c("tf", "target_gene"), drop = FALSE], by = "tf")
  links <- dplyr::inner_join(links, candidates[, c("candidate_id", "fp_id", "target_gene"), drop = FALSE], by = c("fp_id", "target_gene"))
  links <- dplyr::inner_join(links, fp_target_pass[, c("fp_id", "target_gene"), drop = FALSE], by = c("fp_id", "target_gene"))
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
    manifest <- dplyr::bind_rows(.module2_write_table(predicted_tfbs, output_dir, "module1_predicted_tfbs", output_format), .module2_write_table(candidates, output_dir, "module2_fp_target_candidates", output_format), .module2_write_table(tf_target_corr, output_dir, "module2_tf_target_corr", output_format), .module2_write_table(fp_target_corr, output_dir, "module2_fp_target_corr", output_format), .module2_write_table(links, output_dir, "module2_links", output_format), .module2_write_table(activity, output_dir, "module2_condition_activity", output_format), .module2_write_table(qc_summary, output_dir, "module2_qc_summary", "csv"))
    manifest_path <- file.path(output_dir, "module2_manifest.csv")
    readr::write_csv(manifest, manifest_path)
    reports$manifest <- manifest_path
  }
  out <- list(predicted_tfbs = predicted_tfbs, candidates = candidates, tf_target_corr = tf_target_corr, fp_target_corr = fp_target_corr, links = links, condition_activity = activity, qc_summary = qc_summary, manifest = manifest, reports = reports, parameters = list(max_distance_bp = max_distance_bp, n_cores = .module2_default_cores(n_cores), output_format = output_format))
  class(out) <- c("craftgrn_module2", "list")
  out
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

#' Load compact Module 2 outputs
#'
#' @param path Module 2 output directory or module2_manifest.csv path.
#' @return A named list of tables.
#' @export
load_module2_links <- function(path) {
  if (dir.exists(path)) path <- file.path(path, "module2_manifest.csv")
  if (!file.exists(path)) .log_abort("Module 2 manifest not found: {path}")
  manifest <- data.table::fread(path, showProgress = FALSE)
  out <- list()
  for (i in seq_len(nrow(manifest))) {
    p <- as.character(manifest$path[[i]])
    nm <- as.character(manifest$table[[i]])
    if (grepl("\\.parquet$", p, ignore.case = TRUE)) {
      if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet Module 2 outputs.")
      out[[nm]] <- tibble::as_tibble(arrow::read_parquet(p))
    } else {
      out[[nm]] <- tibble::as_tibble(data.table::fread(p, showProgress = FALSE))
    }
  }
  out$manifest <- tibble::as_tibble(manifest)
  out
}

#' Query compact Module 2 links
#'
#' @param module2 Module 2 result list or loaded output list.
#' @param tf Optional TF filter.
#' @param fp_id Optional FP filter.
#' @param target_gene Optional target-gene filter.
#' @param pass_only Keep only passing links.
#' @return A tibble of matching final links.
#' @export
query_module2_links <- function(module2, tf = NULL, fp_id = NULL, target_gene = NULL, pass_only = TRUE) {
  links <- if (is.data.frame(module2$links)) module2$links else module2$module2_links
  if (!is.data.frame(links)) .log_abort("Module 2 links table not found.")
  dt <- data.table::as.data.table(links)
  if (!is.null(tf)) dt <- dt[tf %in% as.character(tf)]
  if (!is.null(fp_id)) dt <- dt[fp_id %in% as.character(fp_id)]
  if (!is.null(target_gene)) dt <- dt[target_gene %in% as.character(target_gene)]
  if (isTRUE(pass_only) && "module2_link_pass" %in% names(dt)) dt <- dt[module2_link_pass %in% TRUE]
  tibble::as_tibble(dt)
}

#' Validate compact Module 2 links
#'
#' @param module2 Module 2 result list or loaded output list.
#' @return TRUE invisibly when valid.
#' @export
validate_module2_links <- function(module2) {
  links <- if (is.data.frame(module2$links)) module2$links else module2$module2_links
  if (!is.data.frame(links)) .log_abort("Module 2 links table not found.")
  need <- c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass")
  missing <- setdiff(need, names(links))
  if (length(missing)) .log_abort("Module 2 links missing required columns: {paste(missing, collapse = \", \")}.")
  invisible(TRUE)
}

#' Export Module 2 links as BEDPE
#'
#' @param module2 Module 2 result list or loaded output list.
#' @param output_file BEDPE output file.
#' @param tf Optional TF subset.
#' @return Output path invisibly.
#' @export
export_tf_target_bedpe <- function(module2, output_file, tf = NULL) {
  links <- query_module2_links(module2, tf = tf, pass_only = TRUE)
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
