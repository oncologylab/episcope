# File: utils_step3_compact_diff_links.R
# Purpose: Compact Module 3 differential-link export from Module 2 manifests.

.module3_cfg <- function(project_config = NULL) {
  cfg <- .module2_cfg(project_config)
  if (is.list(cfg$module3)) cfg <- utils::modifyList(cfg, cfg$module3)
  cfg
}

.module3_cfg_value <- function(cfg, name, default = NULL) {
  if (is.list(cfg) && !is.null(cfg[[name]])) cfg[[name]] else default
}

.module3_safe_label <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

.module3_fmt_cut <- function(x) {
  x <- as.numeric(x)[[1L]]
  out <- formatC(x, format = "f", digits = 3)
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  gsub("\\.", "p", out)
}

.module3_default_output_dir <- function(cfg, project_dir = NULL) {
  base_dir <- .module3_cfg_value(cfg, "base_dir", project_dir)
  if (is.null(base_dir) || !nzchar(base_dir)) base_dir <- "."
  run_tag <- .module3_cfg_value(cfg, "module2_run_label", .module3_cfg_value(cfg, "run_tag", "module3"))
  diff_dir <- .module3_cfg_value(cfg, "diff_dir_name", NULL)
  if (is.null(diff_dir) || !nzchar(diff_dir)) {
    mode <- .module3_cfg_value(cfg, "fp_filter_mode", "log2fc")
    fp_cut <- if (identical(mode, "delta")) {
      .module3_cfg_value(cfg, "fp_delta_cutoff", 0.5)
    } else {
      .module3_cfg_value(cfg, "fp_log2fc_cutoff", 1)
    }
    gene_cut <- .module3_cfg_value(cfg, "gene_log2fc_cutoff", 1)
    diff_dir <- paste0("diff_grn_fp_", mode, "_", .module3_fmt_cut(fp_cut), "_gene_log2fc_", .module3_fmt_cut(gene_cut))
  }
  file.path(base_dir, "benchmark_diff_grn_and_regulatory_topics", run_tag, diff_dir, "diff_links_filtered")
}

.module3_read_table <- function(path, format = NULL, columns = NULL, allow_missing_columns = FALSE) {
  if (identical(format, "parquet") || grepl("\\.parquet$", path, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) .log_abort("Package arrow is required to read Parquet Module 2 outputs.")
    if (is.null(columns)) return(data.table::as.data.table(arrow::read_parquet(path)))
    if (isTRUE(allow_missing_columns)) {
      schema <- arrow::ParquetFileReader$create(path)$GetSchema()
      columns <- intersect(columns, names(schema))
    }
    return(data.table::as.data.table(arrow::read_parquet(path, col_select = columns)))
  }
  if (is.null(columns)) {
    return(data.table::fread(path, showProgress = FALSE))
  }
  if (isTRUE(allow_missing_columns)) {
    available <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
    columns <- intersect(columns, available)
  }
  data.table::fread(path, select = columns, showProgress = FALSE)
}

.module3_manifest_table <- function(module2, table) {
  if (is.character(module2) && length(module2) == 1L) {
    module2 <- load_module2_links(module2)
  }
  man <- module2$manifest
  if (!is.data.frame(man) || !nrow(man)) .log_abort("Module 2 manifest is missing or empty.")
  row <- man[as.character(man$table) == table, , drop = FALSE]
  if (!nrow(row)) .log_abort("Module 2 manifest is missing table `{table}`.")
  as.character(row$path[[1L]])
}

.module3_read_manifest_file <- function(path) {
  if (!file.exists(path)) .log_abort("Manifest file does not exist: {path}")
  out <- data.table::fread(path, showProgress = FALSE)
  if (!all(c("path", "format") %in% names(out))) .log_abort("Manifest must include path and format columns: {path}")
  out
}

.module3_candidate_columns <- function() {
  c(
    "candidate_id", "fp_id", "target_gene", "distance_to_tss",
    "candidate_source", "within_tss_window", "prior_supported", "prior_id",
    "prior_source", "prior_score", "prior_status"
  )
}

.module3_read_manifest_tables <- function(manifest, columns = NULL, allow_missing_columns = FALSE) {
  if (!nrow(manifest)) return(data.table::data.table())
  data.table::rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
    .module3_read_table(
      as.character(manifest$path[[i]]),
      as.character(manifest$format[[i]]),
      columns = columns,
      allow_missing_columns = allow_missing_columns
    )
  }), use.names = TRUE, fill = TRUE)
}

.module3_read_static_corr <- function(module2, table, keep_pass = TRUE) {
  path <- .module3_manifest_table(module2, table)
  man <- if (grepl("_manifest\\.csv$", path)) .module3_read_manifest_file(path) else NULL
  cols <- if (identical(table, "module2_tf_target_corr")) {
    c("tf", "target_gene", "best_r", "best_p", "best_fdr", "pass")
  } else {
    c("fp_id", "target_gene", "best_r", "best_p", "best_fdr", "pass")
  }
  if (is.null(man)) {
    fmt <- if (grepl("\\.parquet$", path, ignore.case = TRUE)) "parquet" else "csv"
    dt <- .module3_read_table(path, fmt, columns = cols)
  } else {
    dt <- data.table::rbindlist(lapply(seq_len(nrow(man)), function(i) {
      .module3_read_table(as.character(man$path[[i]]), as.character(man$format[[i]]), columns = cols)
    }), use.names = TRUE, fill = TRUE)
  }
  if (isTRUE(keep_pass) && "pass" %in% names(dt)) dt <- dt[pass %in% TRUE]
  dt
}

.module3_match_values <- function(mat, rows, col) {
  idx <- match(rows, rownames(mat))
  out <- rep(NA_real_, length(rows))
  ok <- !is.na(idx) & col %in% colnames(mat)
  if (any(ok)) out[ok] <- as.numeric(mat[idx[ok], col])
  out
}

.module3_match_flags <- function(mat, rows, col) {
  out <- .module3_match_values(mat, rows, col)
  out[!is.finite(out)] <- 0
  out > 0
}

.module3_safe_log2fc <- function(num, den, pseudocount = 1) {
  out <- suppressWarnings(log2((as.numeric(num) + pseudocount) / (as.numeric(den) + pseudocount)))
  out[!is.finite(out)] <- NA_real_
  out
}

.module3_fp_signal_mode <- function(fp_signal_mode = NULL, cfg = NULL) {
  mode <- fp_signal_mode %||% .module3_cfg_value(cfg, "fp_signal_mode", "actual")
  mode <- as.character(mode)[[1L]]
  match.arg(mode, c("actual", "link_padded"))
}

.module3_prepare_chunk_links <- function(link_dt, cand_dt, fp_corr, tf_corr) {
  if (!nrow(link_dt)) return(data.table::data.table())
  link_dt <- data.table::as.data.table(link_dt)
  link_dt <- link_dt[module2_link_pass %in% TRUE]
  if (!nrow(link_dt)) return(data.table::data.table())

  cand_dt <- data.table::as.data.table(cand_dt)
  cand_keep <- .module3_candidate_columns()
  cand_dt <- cand_dt[, intersect(cand_keep, names(cand_dt)), with = FALSE]
  data.table::setkey(cand_dt, candidate_id, fp_id, target_gene)
  link_dt <- cand_dt[link_dt, on = c("candidate_id", "fp_id", "target_gene")]

  data.table::setkey(fp_corr, fp_id, target_gene)
  data.table::setkey(tf_corr, tf, target_gene)
  link_dt <- fp_corr[link_dt, on = c("fp_id", "target_gene")]
  data.table::setnames(link_dt, c("best_r", "best_p", "best_fdr"), c("r_gene", "p_gene", "p_adj_gene"), skip_absent = TRUE)
  link_dt <- tf_corr[link_dt, on = c("tf", "target_gene")]
  data.table::setnames(link_dt, c("best_r", "best_p", "best_fdr"), c("r_rna_gene", "p_rna_gene", "p_rna_adj_gene"), skip_absent = TRUE)
  link_dt <- link_dt[is.finite(link_dt$r_gene) & is.finite(link_dt$r_rna_gene)]
  link_dt
}

.module3_build_chunk_delta_prepared <- function(link_dt,
                                                multiomic_data,
                                                case_id,
                                                ctrl_id,
                                                pseudocount = 1,
                                                fp_signal_mode = c("actual", "link_padded")) {
  fp_signal_mode <- match.arg(fp_signal_mode)
  if (!nrow(link_dt)) return(data.table::data.table())
  mats <- multiomic_data$matrices
  case_suffix <- paste0("_", case_id)
  ctrl_suffix <- paste0("_", ctrl_id)
  fp_case <- .module3_match_values(mats$fp_score, link_dt$fp_id, case_id)
  fp_ctrl <- .module3_match_values(mats$fp_score, link_dt$fp_id, ctrl_id)
  tf_case <- .module3_match_values(mats$gene_expr, link_dt$tf, case_id)
  tf_ctrl <- .module3_match_values(mats$gene_expr, link_dt$tf, ctrl_id)
  gene_case <- .module3_match_values(mats$gene_expr, link_dt$target_gene, case_id)
  gene_ctrl <- .module3_match_values(mats$gene_expr, link_dt$target_gene, ctrl_id)
  fp_bound_case <- .module3_match_flags(mats$fp_bound, link_dt$fp_id, case_id)
  fp_bound_ctrl <- .module3_match_flags(mats$fp_bound, link_dt$fp_id, ctrl_id)
  tf_flag_case <- .module3_match_flags(mats$gene_on, link_dt$tf, case_id)
  tf_flag_ctrl <- .module3_match_flags(mats$gene_on, link_dt$tf, ctrl_id)
  gene_flag_case <- .module3_match_flags(mats$gene_on, link_dt$target_gene, case_id)
  gene_flag_ctrl <- .module3_match_flags(mats$gene_on, link_dt$target_gene, ctrl_id)
  active_case <- fp_bound_case & tf_flag_case & gene_flag_case
  active_ctrl <- fp_bound_ctrl & tf_flag_ctrl & gene_flag_ctrl
  fp_case_signal <- fp_case
  fp_ctrl_signal <- fp_ctrl
  if (identical(fp_signal_mode, "link_padded")) {
    fp_case_signal <- data.table::fifelse(active_case, fp_case, 0)
    fp_ctrl_signal <- data.table::fifelse(active_ctrl, fp_ctrl, 0)
  }
  score_base <- as.numeric(link_dt$r_gene) * as.numeric(link_dt$r_rna_gene)
  score_case <- ifelse(active_case, score_base, 0)
  score_ctrl <- ifelse(active_ctrl, score_base, 0)

  out <- data.table::data.table(
    tf = as.character(link_dt$tf),
    gene_key = as.character(link_dt$target_gene),
    peak_id = as.character(link_dt$fp_id)
  )
  out[[paste0("link_score", case_suffix)]] <- score_case
  out[[paste0("link_score", ctrl_suffix)]] <- score_ctrl
  out[[paste0("link_sign", case_suffix)]] <- "+"
  out[[paste0("link_sign", ctrl_suffix)]] <- "+"
  out[[paste0("fp_score", case_suffix)]] <- fp_case_signal
  out[[paste0("fp_score", ctrl_suffix)]] <- fp_ctrl_signal
  out[[paste0("tf_expr", case_suffix)]] <- tf_case
  out[[paste0("tf_expr", ctrl_suffix)]] <- tf_ctrl
  out[[paste0("gene_expr", case_suffix)]] <- gene_case
  out[[paste0("gene_expr", ctrl_suffix)]] <- gene_ctrl
  out[[paste0("active", case_suffix)]] <- active_case
  out[[paste0("active", ctrl_suffix)]] <- active_ctrl
  out[[paste0("fp_bound", case_suffix)]] <- as.integer(fp_bound_case)
  out[[paste0("fp_bound", ctrl_suffix)]] <- as.integer(fp_bound_ctrl)
  out[[paste0("tf_expr_flag", case_suffix)]] <- as.integer(tf_flag_case)
  out[[paste0("tf_expr_flag", ctrl_suffix)]] <- as.integer(tf_flag_ctrl)
  out[[paste0("gene_expr_flag", case_suffix)]] <- as.integer(gene_flag_case)
  out[[paste0("gene_expr_flag", ctrl_suffix)]] <- as.integer(gene_flag_ctrl)
  out <- data.table::copy(out)
  out[, `:=`(
    delta_link_score = score_case - score_ctrl,
    delta_fp_score = fp_case_signal - fp_ctrl_signal,
    delta_tf_expr = tf_case - tf_ctrl,
    delta_gene_expr = gene_case - gene_ctrl,
    log2FC_fp_score = .module3_safe_log2fc(fp_case_signal, fp_ctrl_signal, pseudocount),
    log2FC_tf_expr = .module3_safe_log2fc(tf_case, tf_ctrl, pseudocount),
    log2FC_gene_expr = .module3_safe_log2fc(gene_case, gene_ctrl, pseudocount),
    fc_link_score = (score_case + pseudocount) / (score_ctrl + pseudocount),
    fc_fp_score = (fp_case_signal + pseudocount) / (fp_ctrl_signal + pseudocount),
    fc_tf_expr = (tf_case + pseudocount) / (tf_ctrl + pseudocount),
    fc_gene_expr = (gene_case + pseudocount) / (gene_ctrl + pseudocount),
    edge_change = fifelse(score_case - score_ctrl > 0, "gain", fifelse(score_case - score_ctrl < 0, "loss", "shared")),
    active_both = active_case & active_ctrl,
    active_any = active_case | active_ctrl,
    fp_signal_mode = fp_signal_mode,
    r_gene = as.numeric(link_dt$r_gene),
    p_gene = as.numeric(link_dt$p_gene),
    p_adj_gene = as.numeric(link_dt$p_adj_gene),
    r_rna_gene = as.numeric(link_dt$r_rna_gene),
    p_rna_gene = as.numeric(link_dt$p_rna_gene),
    p_rna_adj_gene = as.numeric(link_dt$p_rna_adj_gene)
  )]
  optional <- intersect(c("distance_to_tss", "candidate_source", "within_tss_window", "prior_supported", "prior_id", "prior_source", "prior_score", "prior_status"), names(link_dt))
  for (nm in optional) out[[nm]] <- link_dt[[nm]]
  out
}

.module3_build_chunk_delta <- function(link_dt,
                                       cand_dt,
                                       fp_corr,
                                       tf_corr,
                                       multiomic_data,
                                       case_id,
                                       ctrl_id,
                                       pseudocount = 1,
                                       fp_signal_mode = c("actual", "link_padded")) {
  fp_signal_mode <- match.arg(fp_signal_mode)
  prepared <- .module3_prepare_chunk_links(link_dt, cand_dt, fp_corr, tf_corr)
  .module3_build_chunk_delta_prepared(
    prepared,
    multiomic_data,
    case_id,
    ctrl_id,
    pseudocount = pseudocount,
    fp_signal_mode = fp_signal_mode
  )
}

.module3_filter_direction <- function(dt, cfg, direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (!nrow(dt)) return(dt)
  mode <- match.arg(as.character(.module3_cfg_value(cfg, "fp_filter_mode", "log2fc"))[[1L]], c("delta", "log2fc"))
  gene_cut <- as.numeric(.module3_cfg_value(cfg, "gene_log2fc_cutoff", 1))[[1L]]
  fp_cut <- if (identical(mode, "delta")) {
    as.numeric(.module3_cfg_value(cfg, "fp_delta_cutoff", 0.5))[[1L]]
  } else {
    as.numeric(.module3_cfg_value(cfg, "fp_log2fc_cutoff", 1))[[1L]]
  }
  expr_cut <- as.numeric(.module3_cfg_value(cfg, "threshold_expr", 10))[[1L]]
  gene_expr_cut <- as.numeric(.module3_cfg_value(cfg, "threshold_gene_expr", expr_cut))[[1L]]
  tf_expr_cut <- as.numeric(.module3_cfg_value(cfg, "threshold_tf_expr", expr_cut))[[1L]]
  fp_score_cut <- as.numeric(.module3_cfg_value(cfg, "threshold_fp_score", 2))[[1L]]
  tf_opposition_cut <- as.numeric(.module3_cfg_value(cfg, "tf_opposition_log2fc_cutoff", gene_cut))[[1L]]

  fp_value <- if (identical(mode, "delta")) dt$delta_fp_score else dt$log2FC_fp_score
  gene_l2 <- dt$log2FC_gene_expr
  tf_l2 <- dt$log2FC_tf_expr
  case_fp <- grep("^fp_score_", names(dt), value = TRUE)[1L]
  ctrl_fp <- grep("^fp_score_", names(dt), value = TRUE)[2L]
  case_tf <- grep("^tf_expr_", names(dt), value = TRUE)[1L]
  ctrl_tf <- grep("^tf_expr_", names(dt), value = TRUE)[2L]
  case_gene <- grep("^gene_expr_", names(dt), value = TRUE)[1L]
  ctrl_gene <- grep("^gene_expr_", names(dt), value = TRUE)[2L]
  fp_high <- ifelse(gene_l2 >= 0, dt[[case_fp]], dt[[ctrl_fp]])
  tf_high <- ifelse(gene_l2 >= 0, dt[[case_tf]], dt[[ctrl_tf]])
  gene_high <- ifelse(gene_l2 >= 0, dt[[case_gene]], dt[[ctrl_gene]])
  tf_opp <- is.finite(tf_l2) & is.finite(gene_l2) & abs(tf_l2) >= tf_opposition_cut & sign(tf_l2) == -sign(gene_l2)
  base_keep <- dt$active_any %in% TRUE &
    is.finite(gene_l2) &
    is.finite(fp_value) &
    is.finite(fp_high) & fp_high > fp_score_cut &
    is.finite(tf_high) & tf_high > tf_expr_cut &
    is.finite(gene_high) & gene_high > gene_expr_cut &
    !tf_opp
  if (identical(direction, "up")) {
    keep <- base_keep & gene_l2 >= gene_cut & fp_value >= fp_cut
  } else {
    keep <- base_keep & gene_l2 <= -gene_cut & fp_value <= -fp_cut
  }
  dt[keep]
}

.module3_comparison_specs <- function(compar, cfg) {
  if (is.null(compar)) {
    base_dir <- .module3_cfg_value(cfg, "base_dir", ".")
    compar <- file.path(base_dir, "data", "episcope_comparisons.csv")
  }
  if (is.character(compar) && length(compar) == 1L) {
    if (!file.exists(compar)) .log_abort("Comparison file does not exist: {compar}")
    compar <- data.table::fread(compar, showProgress = FALSE)
  }
  if (!is.data.frame(compar)) .log_abort("`compar` must be a data.frame or CSV path.")
  compar <- data.table::as.data.table(compar)
  if (!all(c("cond1_label", "cond2_label") %in% names(compar))) {
    .log_abort("`compar` must include cond1_label and cond2_label columns.")
  }
  unique(compar[, .(cond1_label = as.character(cond1_label), cond2_label = as.character(cond2_label))])
}

#' Prepare differential links for Module 3
#'
#' Converts Module 2 link manifests into the filtered differential-link
#' files consumed by CraftGRN topic-modeling utilities. This avoids writing full
#' per-condition GRN matrices and keeps Module 3 compatible with the existing
#' `*_filtered_links_up.csv` and `*_filtered_links_down.csv` contract.
#'
#' @param module2 Module 2 object returned by [link_tf_targets()] or a path to a
#'   Module 2 output directory containing `module2_manifest.csv`.
#' @param multiomic_data CraftGRN multiomic object returned by
#'   [load_prep_multiomic_data()].
#' @param compar Comparison table or CSV path with `cond1_label` and
#'   `cond2_label`. If `NULL`, `data/episcope_comparisons.csv` under
#'   `base_dir` is used.
#' @param project_config Project config list or YAML path.
#' @param output_dir Directory for filtered differential links. If `NULL`, a
#'   Module 3 benchmark directory is derived from the project config.
#' @param n_cores Number of data.table threads to use while reading and joining chunks. Defaults to all available cores.
#'   Comparison-level parallelism is controlled by `module3_comparison_workers` in the project config and defaults to 1 for RAM safety.
#' @param pseudocount Pseudocount for log2 fold-change calculations.
#' @param fp_signal_mode FP signal used for differential FP fold changes.
#'   actual uses the measured FP score in both conditions. link_padded
#'   sets the FP score to zero in conditions where the TF-FP-gene link is not
#'   active before calculating delta_fp_score and log2FC_fp_score.
#' @param overwrite Overwrite existing filtered link files.
#' @param verbose Emit concise progress messages.
#'
#' @return A tibble manifest with one row per comparison.
#' @export
module3_prepare_differential_links <- function(module2,
                                               multiomic_data,
                                               compar = NULL,
                                               project_config = NULL,
                                               output_dir = NULL,
                                               n_cores = NULL,
                                               pseudocount = 1,
                                               fp_signal_mode = NULL,
                                               overwrite = FALSE,
                                               verbose = TRUE) {
  cfg <- .module3_cfg(project_config)
  fp_signal_mode <- .module3_fp_signal_mode(fp_signal_mode, cfg)
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object.")
  validate_multiomic_object(multiomic_data)
  if (is.character(module2) && length(module2) == 1L) module2 <- load_module2_links(module2)
  if (!inherits(module2, "craftgrn_module2") && !is.list(module2)) .log_abort("`module2` must be a Module 2 object or output directory.")
  output_dir <- output_dir %||% .module3_default_output_dir(cfg)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(output_dir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- .module3_comparison_specs(compar, cfg)
  if (!nrow(specs)) .log_abort("No Module 3 comparisons to process.")

  link_manifest <- .module3_read_manifest_file(.module3_manifest_table(module2, "module2_links"))
  cand_manifest <- .module3_read_manifest_file(.module3_manifest_table(module2, "module2_fp_target_candidates"))
  aligned_candidate_chunks <- nrow(cand_manifest) == nrow(link_manifest)
  if (aligned_candidate_chunks && all(c("chunk_id") %in% names(cand_manifest)) && all(c("chunk_id") %in% names(link_manifest))) {
    aligned_candidate_chunks <- identical(as.character(cand_manifest$chunk_id), as.character(link_manifest$chunk_id))
  }
  if (isTRUE(verbose)) .log_inform("Module 3 bridge: loading Module 2 correlation summaries.")
  fp_corr <- .module3_read_static_corr(module2, "module2_fp_target_corr", keep_pass = TRUE)
  tf_corr <- .module3_read_static_corr(module2, "module2_tf_target_corr", keep_pass = TRUE)
  data.table::setkey(fp_corr, fp_id, target_gene)
  data.table::setkey(tf_corr, tf, target_gene)
  shared_candidates <- NULL
  if (!aligned_candidate_chunks) {
    if (isTRUE(verbose)) {
      .log_inform("Module 3 bridge: loading {nrow(cand_manifest)} shared candidate chunk(s) for {nrow(link_manifest)} link chunk(s).")
    }
    shared_candidates <- .module3_read_manifest_tables(
      cand_manifest,
      columns = .module3_candidate_columns(),
      allow_missing_columns = TRUE
    )
  }
  dt_threads <- .module2_default_cores(n_cores)
  old_dt_threads <- data.table::setDTthreads(dt_threads)
  on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
  comparison_workers <- min(as.integer(.module3_cfg_value(cfg, "module3_comparison_workers", 1))[[1L]], nrow(specs))
  if (!is.finite(comparison_workers) || comparison_workers < 1L) comparison_workers <- 1L
  if (isTRUE(verbose)) .log_inform("Module 3 bridge: {nrow(specs)} comparison(s), {nrow(link_manifest)} link chunk(s), {dt_threads} data.table thread(s), {comparison_workers} comparison worker(s), FP signal mode {fp_signal_mode}.")

  manifest_path <- file.path(output_dir, "filtered_links_manifest.csv")
  old_manifest <- if (file.exists(manifest_path)) {
    tryCatch(data.table::fread(manifest_path, showProgress = FALSE), error = function(e) data.table::data.table())
  } else {
    data.table::data.table()
  }
  comparison_info <- lapply(seq_len(nrow(specs)), function(i) {
    case_id <- specs$cond1_label[[i]]
    ctrl_id <- specs$cond2_label[[i]]
    stem <- paste0(.module3_safe_label(case_id), "_vs_", .module3_safe_label(ctrl_id))
    up_path <- file.path(output_dir, paste0(stem, "_filtered_links_up.csv"))
    down_path <- file.path(output_dir, paste0(stem, "_filtered_links_down.csv"))
    qc_path <- file.path(cache_dir, paste0(stem, "_qc.csv"))
    skipped <- !isTRUE(overwrite) && file.exists(up_path) && file.exists(down_path)
    list(case_id = case_id, ctrl_id = ctrl_id, stem = stem, up_path = up_path, down_path = down_path, qc_path = qc_path, skipped = skipped)
  })
  active_idx <- which(!vapply(comparison_info, function(x) isTRUE(x$skipped), logical(1)))
  res <- vector("list", nrow(specs))
  for (i in setdiff(seq_len(nrow(specs)), active_idx)) {
    info <- comparison_info[[i]]
    old_row <- if ("comparison_id" %in% names(old_manifest)) old_manifest[comparison_id == info$stem] else data.table::data.table()
    n_up <- if (nrow(old_row) && "n_up" %in% names(old_row)) old_row$n_up[[1L]] else NA_integer_
    n_down <- if (nrow(old_row) && "n_down" %in% names(old_row)) old_row$n_down[[1L]] else NA_integer_
    res[[i]] <- tibble::tibble(comparison_id = info$stem, case_id = info$case_id, ctrl_id = info$ctrl_id, up_path = info$up_path, down_path = info$down_path, n_up = n_up, n_down = n_down, fp_signal_mode = fp_signal_mode, skipped = TRUE)
  }

  if (length(active_idx)) {
    up_parts <- lapply(active_idx, function(i) vector("list", nrow(link_manifest)))
    down_parts <- lapply(active_idx, function(i) vector("list", nrow(link_manifest)))
    qc_parts <- lapply(active_idx, function(i) vector("list", nrow(link_manifest)))
    for (j in seq_len(nrow(link_manifest))) {
      if (isTRUE(verbose)) .log_inform("Module 3 bridge: preparing link chunk {j}/{nrow(link_manifest)} for {length(active_idx)} comparison(s).")
      link_dt <- .module3_read_table(as.character(link_manifest$path[[j]]), as.character(link_manifest$format[[j]]), columns = c("link_id", "tf", "fp_id", "target_gene", "candidate_id", "module2_link_pass"))
      cand_dt <- if (aligned_candidate_chunks) {
        .module3_read_table(
          as.character(cand_manifest$path[[j]]),
          as.character(cand_manifest$format[[j]]),
          columns = .module3_candidate_columns(),
          allow_missing_columns = TRUE
        )
      } else {
        shared_candidates
      }
      n_links <- nrow(link_dt)
      prepared <- .module3_prepare_chunk_links(link_dt, cand_dt, fp_corr, tf_corr)
      n_prepared <- nrow(prepared)
      rm(link_dt, cand_dt)
      for (k in seq_along(active_idx)) {
        i <- active_idx[[k]]
        info <- comparison_info[[i]]
        delta_dt <- .module3_build_chunk_delta_prepared(
          prepared,
          multiomic_data,
          info$case_id,
          info$ctrl_id,
          pseudocount = pseudocount,
          fp_signal_mode = fp_signal_mode
        )
        up_dt <- .module3_filter_direction(delta_dt, cfg, "up")
        down_dt <- .module3_filter_direction(delta_dt, cfg, "down")
        up_parts[[k]][[j]] <- up_dt
        down_parts[[k]][[j]] <- down_dt
        qc_parts[[k]][[j]] <- data.table::data.table(chunk_id = j, n_links = n_links, n_prepared = n_prepared, n_delta = nrow(delta_dt), n_up = nrow(up_dt), n_down = nrow(down_dt))
        rm(delta_dt, up_dt, down_dt)
      }
      rm(prepared)
      invisible(gc(FALSE))
    }
    for (k in seq_along(active_idx)) {
      i <- active_idx[[k]]
      info <- comparison_info[[i]]
      up <- data.table::rbindlist(up_parts[[k]], use.names = TRUE, fill = TRUE)
      down <- data.table::rbindlist(down_parts[[k]], use.names = TRUE, fill = TRUE)
      if (nrow(up)) up <- unique(up, by = c("tf", "gene_key", "peak_id"))
      if (nrow(down)) down <- unique(down, by = c("tf", "gene_key", "peak_id"))
      data.table::fwrite(up, info$up_path)
      data.table::fwrite(down, info$down_path)
      data.table::fwrite(data.table::rbindlist(qc_parts[[k]], use.names = TRUE), info$qc_path)
      res[[i]] <- tibble::tibble(comparison_id = info$stem, case_id = info$case_id, ctrl_id = info$ctrl_id, up_path = info$up_path, down_path = info$down_path, n_up = nrow(up), n_down = nrow(down), fp_signal_mode = fp_signal_mode, skipped = FALSE)
    }
  }
  manifest <- dplyr::bind_rows(res)
  readr::write_csv(manifest, manifest_path)
  if (isTRUE(verbose)) {
    .log_inform("Module 3 bridge wrote {sum(!manifest$skipped)} comparison(s) to {output_dir}.")
    .log_inform("Module 3 filtered links manifest: {manifest_path}")
  }
  manifest
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
