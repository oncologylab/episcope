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
  file.path(base_dir, "regulatory_topics", "differential_links")
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

.module3_manifest_row <- function(module2, table) {
  if (is.character(module2) && length(module2) == 1L) {
    module2 <- load_module2_links(module2)
  }
  man <- module2$manifest
  if (!is.data.frame(man) || !nrow(man)) .log_abort("Module 2 manifest is missing or empty.")
  row <- man[as.character(man$table) == table, , drop = FALSE]
  if (!nrow(row)) .log_abort("Module 2 manifest is missing table `{table}`.")
  data.table::as.data.table(row[1L, , drop = FALSE])
}

.module3_read_manifest_file <- function(path) {
  if (!file.exists(path)) .log_abort("Manifest file does not exist: {path}")
  out <- data.table::fread(path, showProgress = FALSE)
  if (!all(c("path", "format") %in% names(out))) .log_abort("Manifest must include path and format columns: {path}")
  out
}

.module3_table_rows <- function(module2, table) {
  row <- .module3_manifest_row(module2, table)
  fmt <- as.character(row$format[[1L]])
  path <- as.character(row$path[[1L]])
  if (identical(fmt, "manifest")) {
    return(.module3_read_manifest_file(path))
  }
  row$chunk_id <- 1L
  row[, c("chunk_id", "path", "format", "n_rows"), with = FALSE]
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
  cols <- if (identical(table, "module2_tf_target_corr")) {
    c("tf", "target_gene", "best_r", "best_p", "best_fdr", "pass")
  } else {
    c("fp_id", "target_gene", "best_r", "best_p", "best_fdr", "pass")
  }
  man <- .module3_table_rows(module2, table)
  dt <- data.table::rbindlist(lapply(seq_len(nrow(man)), function(i) {
    .module3_read_table(as.character(man$path[[i]]), as.character(man$format[[i]]), columns = cols)
  }), use.names = TRUE, fill = TRUE)
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

.module3_first_present <- function(dt, choices) {
  hit <- intersect(choices, names(dt))
  if (length(hit)) hit[[1L]] else NULL
}

.module3_read_de_results <- function(rna_de_results = NULL, cfg = NULL) {
  path <- rna_de_results %||% .module3_cfg_value(cfg, "rna_de_results", NULL)
  if (is.null(path) || !nzchar(as.character(path)[[1L]])) return(NULL)
  if (is.data.frame(path)) {
    de <- data.table::as.data.table(path)
  } else {
    path <- as.character(path)[[1L]]
    if (!file.exists(path)) .log_abort("RNA differential result file does not exist: {path}")
    de <- data.table::fread(path, showProgress = FALSE)
  }
  gene_col <- .module3_first_present(de, c("gene_key", "gene_symbol", "Symbol", "symbol", "gene", "gene_id", "ENSEMBL"))
  lfc_col <- .module3_first_present(de, c("log2fc_rna", "log2FC_rna", "log2FoldChange", "log2fc", "log2FC"))
  if (is.null(gene_col) || is.null(lfc_col)) {
    .log_abort("RNA differential results must include a gene column and a log2 fold-change column.")
  }
  p_col <- .module3_first_present(de, c("pvalue", "p_value", "p"))
  padj_col <- .module3_first_present(de, c("padj", "fdr", "qvalue", "q_value"))
  source_col <- .module3_first_present(de, c("de_source", "source"))
  test_col <- .module3_first_present(de, c("de_test_id", "test_id", "contrast", "comparison_id"))
  comp_col <- .module3_first_present(de, c("comparison_id", "contrast_id"))
  out <- data.table::data.table(
    comparison_id = if (!is.null(comp_col)) as.character(de[[comp_col]]) else NA_character_,
    gene_key = as.character(de[[gene_col]]),
    log2fc_rna = suppressWarnings(as.numeric(de[[lfc_col]])),
    pvalue_rna = if (!is.null(p_col)) suppressWarnings(as.numeric(de[[p_col]])) else NA_real_,
    padj_rna = if (!is.null(padj_col)) suppressWarnings(as.numeric(de[[padj_col]])) else NA_real_,
    de_source = if (!is.null(source_col)) as.character(de[[source_col]]) else "external",
    de_test_id = if (!is.null(test_col)) as.character(de[[test_col]]) else NA_character_
  )
  out <- out[!is.na(gene_key) & nzchar(gene_key) & is.finite(log2fc_rna)]
  if (!nrow(out)) return(NULL)
  out <- out[order(is.na(padj_rna), padj_rna)]
  unique(out, by = c("comparison_id", "gene_key"))
}

.module3_lookup_de <- function(de_dt, genes) {
  out <- list(
    log2fc = rep(NA_real_, length(genes)),
    pvalue = rep(NA_real_, length(genes)),
    padj = rep(NA_real_, length(genes)),
    source = rep("direct_fc", length(genes)),
    test_id = rep(NA_character_, length(genes))
  )
  if (is.null(de_dt) || !nrow(de_dt)) return(out)
  idx <- match(as.character(genes), de_dt$gene_key)
  ok <- !is.na(idx)
  if (any(ok)) {
    pvalue <- if ("pvalue_rna" %in% names(de_dt)) de_dt$pvalue_rna else rep(NA_real_, nrow(de_dt))
    padj <- if ("padj_rna" %in% names(de_dt)) de_dt$padj_rna else rep(NA_real_, nrow(de_dt))
    source <- if ("de_source" %in% names(de_dt)) de_dt$de_source else rep("external", nrow(de_dt))
    test_id <- if ("de_test_id" %in% names(de_dt)) de_dt$de_test_id else rep(NA_character_, nrow(de_dt))
    out$log2fc[ok] <- de_dt$log2fc_rna[idx[ok]]
    out$pvalue[ok] <- pvalue[idx[ok]]
    out$padj[ok] <- padj[idx[ok]]
    out$source[ok] <- source[idx[ok]]
    out$test_id[ok] <- test_id[idx[ok]]
  }
  out
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
                                                cond1_matrix_id,
                                                cond2_matrix_id,
                                                pseudocount = 1,
                                                fp_signal_mode = c("actual", "link_padded"),
                                                rna_de_results = NULL,
                                                comparison_id = NULL) {
  fp_signal_mode <- match.arg(fp_signal_mode)
  if (!nrow(link_dt)) return(data.table::data.table())
  mats <- multiomic_data$matrices
  fp_cond1 <- .module3_match_values(mats$fp_score, link_dt$fp_id, cond1_matrix_id)
  fp_cond2 <- .module3_match_values(mats$fp_score, link_dt$fp_id, cond2_matrix_id)
  tf_cond1 <- .module3_match_values(mats$gene_expr, link_dt$tf, cond1_matrix_id)
  tf_cond2 <- .module3_match_values(mats$gene_expr, link_dt$tf, cond2_matrix_id)
  gene_cond1 <- .module3_match_values(mats$gene_expr, link_dt$target_gene, cond1_matrix_id)
  gene_cond2 <- .module3_match_values(mats$gene_expr, link_dt$target_gene, cond2_matrix_id)
  fp_bound_cond1 <- .module3_match_flags(mats$fp_bound, link_dt$fp_id, cond1_matrix_id)
  fp_bound_cond2 <- .module3_match_flags(mats$fp_bound, link_dt$fp_id, cond2_matrix_id)
  tf_flag_cond1 <- .module3_match_flags(mats$gene_on, link_dt$tf, cond1_matrix_id)
  tf_flag_cond2 <- .module3_match_flags(mats$gene_on, link_dt$tf, cond2_matrix_id)
  gene_flag_cond1 <- .module3_match_flags(mats$gene_on, link_dt$target_gene, cond1_matrix_id)
  gene_flag_cond2 <- .module3_match_flags(mats$gene_on, link_dt$target_gene, cond2_matrix_id)
  active_cond1 <- fp_bound_cond1 & tf_flag_cond1 & gene_flag_cond1
  active_cond2 <- fp_bound_cond2 & tf_flag_cond2 & gene_flag_cond2
  fp_cond1_signal <- fp_cond1
  fp_cond2_signal <- fp_cond2
  if (identical(fp_signal_mode, "link_padded")) {
    fp_cond1_signal <- data.table::fifelse(active_cond1, fp_cond1, 0)
    fp_cond2_signal <- data.table::fifelse(active_cond2, fp_cond2, 0)
  }
  score_base <- as.numeric(link_dt$r_gene) * as.numeric(link_dt$r_rna_gene)
  score_cond1 <- ifelse(active_cond1, score_base, 0)
  score_cond2 <- ifelse(active_cond2, score_base, 0)
  de_dt <- rna_de_results
  if (!is.null(de_dt) && nrow(de_dt) && !is.null(comparison_id) && "comparison_id" %in% names(de_dt)) {
    comp_id <- as.character(comparison_id)
    scoped <- de_dt[is.na(de_dt[["comparison_id"]]) | de_dt[["comparison_id"]] == comp_id]
    de_dt <- scoped
  }
  de_gene <- .module3_lookup_de(de_dt, link_dt$target_gene)
  de_tf <- .module3_lookup_de(de_dt, link_dt$tf)
  direct_gene_lfc <- .module3_safe_log2fc(gene_cond1, gene_cond2, pseudocount)
  direct_tf_lfc <- .module3_safe_log2fc(tf_cond1, tf_cond2, pseudocount)
  gene_lfc <- de_gene$log2fc
  tf_lfc <- de_tf$log2fc
  gene_lfc[!is.finite(gene_lfc)] <- direct_gene_lfc[!is.finite(gene_lfc)]
  tf_lfc[!is.finite(tf_lfc)] <- direct_tf_lfc[!is.finite(tf_lfc)]

  out <- data.table::data.table(
    tf = as.character(link_dt$tf),
    gene_key = as.character(link_dt$target_gene),
    peak_id = as.character(link_dt$fp_id)
  )
  out$link_score_cond1 <- score_cond1
  out$link_score_cond2 <- score_cond2
  out$link_sign_cond1 <- "+"
  out$link_sign_cond2 <- "+"
  out$fp_score_cond1 <- fp_cond1_signal
  out$fp_score_cond2 <- fp_cond2_signal
  out$tf_expr_cond1 <- tf_cond1
  out$tf_expr_cond2 <- tf_cond2
  out$gene_expr_cond1 <- gene_cond1
  out$gene_expr_cond2 <- gene_cond2
  out$active_cond1 <- active_cond1
  out$active_cond2 <- active_cond2
  out$fp_bound_cond1 <- as.integer(fp_bound_cond1)
  out$fp_bound_cond2 <- as.integer(fp_bound_cond2)
  out$tf_expr_flag_cond1 <- as.integer(tf_flag_cond1)
  out$tf_expr_flag_cond2 <- as.integer(tf_flag_cond2)
  out$gene_expr_flag_cond1 <- as.integer(gene_flag_cond1)
  out$gene_expr_flag_cond2 <- as.integer(gene_flag_cond2)
  out <- data.table::copy(out)
  out[, `:=`(
    delta_link_score = score_cond1 - score_cond2,
    delta_fp_score = fp_cond1_signal - fp_cond2_signal,
    delta_tf_expr = tf_cond1 - tf_cond2,
    delta_gene_expr = gene_cond1 - gene_cond2,
    log2FC_fp_score = .module3_safe_log2fc(fp_cond1_signal, fp_cond2_signal, pseudocount),
    log2FC_tf_expr = tf_lfc,
    log2FC_gene_expr = gene_lfc,
    log2FC_tf_expr_direct = direct_tf_lfc,
    log2FC_gene_expr_direct = direct_gene_lfc,
    p_rna_de_gene = de_gene$pvalue,
    p_adj_rna_de_gene = de_gene$padj,
    p_rna_de_tf = de_tf$pvalue,
    p_adj_rna_de_tf = de_tf$padj,
    de_source_gene = de_gene$source,
    de_source_tf = de_tf$source,
    de_test_id_gene = de_gene$test_id,
    de_test_id_tf = de_tf$test_id,
    fc_link_score = (score_cond1 + pseudocount) / (score_cond2 + pseudocount),
    fc_fp_score = (fp_cond1_signal + pseudocount) / (fp_cond2_signal + pseudocount),
    fc_tf_expr = (tf_cond1 + pseudocount) / (tf_cond2 + pseudocount),
    fc_gene_expr = (gene_cond1 + pseudocount) / (gene_cond2 + pseudocount),
    edge_change = fifelse(score_cond1 - score_cond2 > 0, "gain", fifelse(score_cond1 - score_cond2 < 0, "loss", "shared")),
    active_both = active_cond1 & active_cond2,
    active_any = active_cond1 | active_cond2,
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
                                       cond1_matrix_id,
                                       cond2_matrix_id,
                                       pseudocount = 1,
                                       fp_signal_mode = c("actual", "link_padded"),
                                       rna_de_results = NULL) {
  fp_signal_mode <- match.arg(fp_signal_mode)
  prepared <- .module3_prepare_chunk_links(link_dt, cand_dt, fp_corr, tf_corr)
  .module3_build_chunk_delta_prepared(
    prepared,
    multiomic_data,
    cond1_matrix_id,
    cond2_matrix_id,
    pseudocount = pseudocount,
    fp_signal_mode = fp_signal_mode,
    rna_de_results = rna_de_results
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
  de_padj_cut <- as.numeric(.module3_cfg_value(cfg, "gene_de_padj_cutoff", NA_real_))[[1L]]
  fp_high <- ifelse(gene_l2 >= 0, dt$fp_score_cond1, dt$fp_score_cond2)
  tf_high <- ifelse(gene_l2 >= 0, dt$tf_expr_cond1, dt$tf_expr_cond2)
  gene_high <- ifelse(gene_l2 >= 0, dt$gene_expr_cond1, dt$gene_expr_cond2)
  tf_opp <- is.finite(tf_l2) & is.finite(gene_l2) & abs(tf_l2) >= tf_opposition_cut & sign(tf_l2) == -sign(gene_l2)
  base_keep <- dt$active_any %in% TRUE &
    is.finite(gene_l2) &
    is.finite(fp_value) &
    is.finite(fp_high) & fp_high > fp_score_cut &
    is.finite(tf_high) & tf_high > tf_expr_cut &
    is.finite(gene_high) & gene_high > gene_expr_cut &
    !tf_opp
  if (is.finite(de_padj_cut)) {
    base_keep <- base_keep & is.finite(dt$p_adj_rna_de_gene) & dt$p_adj_rna_de_gene <= de_padj_cut
  }
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
  if (!all(c("cond1_label", "cond2_label") %in% names(compar)) &&
      !all(c("cond1_matrix_id", "cond2_matrix_id") %in% names(compar))) {
    .log_abort("`compar` must include cond1_label and cond2_label columns.")
  }
  if (!"cond1_label" %in% names(compar)) compar[, cond1_label := as.character(cond1_matrix_id)]
  if (!"cond2_label" %in% names(compar)) compar[, cond2_label := as.character(cond2_matrix_id)]
  compar[, `:=`(cond1_label = as.character(cond1_label), cond2_label = as.character(cond2_label))]
  if ("cond1_matrix_id" %in% names(compar)) {
    data.table::set(compar, j = "cond1_matrix_id", value = as.character(compar[["cond1_matrix_id"]]))
  } else {
    data.table::set(compar, j = "cond1_matrix_id", value = compar[["cond1_label"]])
  }
  if ("cond2_matrix_id" %in% names(compar)) {
    data.table::set(compar, j = "cond2_matrix_id", value = as.character(compar[["cond2_matrix_id"]]))
  } else {
    data.table::set(compar, j = "cond2_matrix_id", value = compar[["cond2_label"]])
  }
  if ("comparison_id" %in% names(compar)) {
    compar[, comparison_id := as.character(compar[["comparison_id"]])]
  } else if ("contrast_id" %in% names(compar)) {
    compar[, comparison_id := as.character(compar[["contrast_id"]])]
  } else if (all(c("cond1_base", "cond2_base") %in% names(compar))) {
    compar[, comparison_id := paste(compar[["cond1_base"]], compar[["cond2_base"]], sep = "_vs_")]
  } else {
    compar[, comparison_id := paste(compar[["cond1_label"]], compar[["cond2_label"]], sep = "_vs_")]
  }
  if (all(c("cond1_base", "cond2_base") %in% names(compar))) {
    fallback_id <- paste(as.character(compar$cond1_base), as.character(compar$cond2_base), sep = "_vs_")
  } else {
    fallback_id <- paste(compar$cond1_label, compar$cond2_label, sep = "_vs_")
  }
  bad_id <- is.na(compar$comparison_id) | !nzchar(trimws(compar$comparison_id))
  compar[bad_id, comparison_id := fallback_id[bad_id]]
  compar[, comparison_id := .module3_safe_label(compar[["comparison_id"]])]
  if ("cond1_id" %in% names(compar)) {
    data.table::set(compar, j = "cond1_id", value = as.character(compar[["cond1_id"]]))
  } else if ("cond1_base" %in% names(compar)) {
    data.table::set(compar, j = "cond1_id", value = as.character(compar[["cond1_base"]]))
  } else if ("case_label" %in% names(compar)) {
    data.table::set(compar, j = "cond1_id", value = as.character(compar[["case_label"]]))
  } else {
    data.table::set(compar, j = "cond1_id", value = compar[["cond1_label"]])
  }
  if ("cond2_id" %in% names(compar)) {
    data.table::set(compar, j = "cond2_id", value = as.character(compar[["cond2_id"]]))
  } else if ("cond2_base" %in% names(compar)) {
    data.table::set(compar, j = "cond2_id", value = as.character(compar[["cond2_base"]]))
  } else if ("ctrl_label" %in% names(compar)) {
    data.table::set(compar, j = "cond2_id", value = as.character(compar[["ctrl_label"]]))
  } else {
    data.table::set(compar, j = "cond2_id", value = compar[["cond2_label"]])
  }
  if (!"cond1_display" %in% names(compar)) data.table::set(compar, j = "cond1_display", value = compar[["cond1_id"]])
  if (!"cond2_display" %in% names(compar)) data.table::set(compar, j = "cond2_display", value = compar[["cond2_id"]])
  if (!"comparison_display" %in% names(compar)) {
    data.table::set(
      compar,
      j = "comparison_display",
      value = paste(compar[["cond1_display"]], compar[["cond2_display"]], sep = " vs ")
    )
  }
  label_out <- as.character(compar[["comparison_display"]])
  if ("comparison_label" %in% names(compar)) {
    explicit_label <- as.character(compar[["comparison_label"]])
    has_explicit_label <- !is.na(explicit_label) & nzchar(trimws(explicit_label))
    label_out[has_explicit_label] <- explicit_label[has_explicit_label]
  }
  if (!"comparison_group" %in% names(compar)) data.table::set(compar, j = "comparison_group", value = NA_character_)
  unique(data.table::data.table(
    comparison_id = compar[["comparison_id"]],
    comparison_label = label_out,
    comparison_group = as.character(compar[["comparison_group"]]),
    cond1_id = compar[["cond1_id"]],
    cond2_id = compar[["cond2_id"]],
    cond1_label = as.character(compar[["cond1_display"]]),
    cond2_label = as.character(compar[["cond2_display"]]),
    cond1_matrix_id = compar[["cond1_matrix_id"]],
    cond2_matrix_id = compar[["cond2_matrix_id"]]
  ))
}

#' Prepare differential links for Module 3
#'
#' Converts Module 2 link manifests into the filtered differential-link
#' files consumed by CraftGRN topic-modeling utilities. This avoids writing full
#' per-condition GRN matrices and keeps Module 3 compatible with the existing
#' `*_filtered_links_up.csv` and `*_filtered_links_down.csv` contract.
#'
#' @param module2 Module 2 object returned by [predict_tf_targets()] or a path to a
#'   Module 2 output directory containing `module2_manifest.csv`.
#' @param multiomic_data CraftGRN multiomic object returned by
#'   [load_prep_multiomic_data()].
#' @param compar Comparison table or CSV path with `cond1_label` and
#'   `cond2_label`. If `NULL`, `data/episcope_comparisons.csv` under
#'   `base_dir` is used.
#' @param project_config Project config list or YAML path.
#' @param output_dir Directory for filtered differential links. If `NULL`,
#'   `regulatory_topics/differential_links` under `base_dir` is used.
#' @param n_cores Number of data.table threads to use while reading and joining chunks. Defaults to all available cores.
#'   Comparison-level parallelism is controlled by `module3_comparison_workers` in the project config and defaults to 1 for RAM safety.
#' @param pseudocount Pseudocount for log2 fold-change calculations.
#' @param rna_de_results Optional standardized RNA differential expression
#'   table or CSV. When provided, target-gene and TF log2 fold changes are read
#'   from this table and direct condition fold changes are used only for missing
#'   genes.
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
                                               rna_de_results = NULL,
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
  qc_dir <- file.path(output_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  specs <- .module3_comparison_specs(compar, cfg)
  if (!nrow(specs)) .log_abort("No Module 3 comparisons to process.")
  de_dt <- .module3_read_de_results(rna_de_results, cfg)
  if (isTRUE(verbose) && !is.null(de_dt)) .log_inform("Module 3 bridge: using RNA DE results for {nrow(de_dt)} gene(s).")

  link_manifest <- .module3_table_rows(module2, "module2_links")
  cand_manifest <- .module3_table_rows(module2, "module2_fp_target_candidates")
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
    cond1_matrix_id <- specs$cond1_matrix_id[[i]]
    cond2_matrix_id <- specs$cond2_matrix_id[[i]]
    stem <- specs$comparison_id[[i]]
    up_path <- file.path(output_dir, paste0(stem, "_filtered_links_up.csv"))
    down_path <- file.path(output_dir, paste0(stem, "_filtered_links_down.csv"))
    skipped <- !isTRUE(overwrite) && file.exists(up_path) && file.exists(down_path)
    list(
      cond1_id = specs$cond1_id[[i]],
      cond2_id = specs$cond2_id[[i]],
      cond1_label = specs$cond1_label[[i]],
      cond2_label = specs$cond2_label[[i]],
      cond1_matrix_id = cond1_matrix_id,
      cond2_matrix_id = cond2_matrix_id,
      stem = stem,
      comparison_label = specs$comparison_label[[i]],
      comparison_group = specs$comparison_group[[i]],
      up_path = up_path,
      down_path = down_path,
      skipped = skipped
    )
  })
  active_idx <- which(!vapply(comparison_info, function(x) isTRUE(x$skipped), logical(1)))
  res <- vector("list", nrow(specs))
  for (i in setdiff(seq_len(nrow(specs)), active_idx)) {
    info <- comparison_info[[i]]
    old_row <- if ("comparison_id" %in% names(old_manifest)) old_manifest[comparison_id == info$stem] else data.table::data.table()
    n_up <- if (nrow(old_row) && "n_up" %in% names(old_row)) old_row$n_up[[1L]] else NA_integer_
    n_down <- if (nrow(old_row) && "n_down" %in% names(old_row)) old_row$n_down[[1L]] else NA_integer_
    res[[i]] <- tibble::tibble(comparison_id = info$stem, comparison_label = info$comparison_label, comparison_group = info$comparison_group, cond1_id = info$cond1_id, cond2_id = info$cond2_id, cond1_label = info$cond1_label, cond2_label = info$cond2_label, cond1_matrix_id = info$cond1_matrix_id, cond2_matrix_id = info$cond2_matrix_id, up_path = info$up_path, down_path = info$down_path, n_up = n_up, n_down = n_down, fp_signal_mode = fp_signal_mode, skipped = TRUE)
  }

  qc_rows <- list()
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
          info$cond1_matrix_id,
          info$cond2_matrix_id,
          pseudocount = pseudocount,
          fp_signal_mode = fp_signal_mode,
          rna_de_results = de_dt,
          comparison_id = info$stem
        )
        up_dt <- .module3_filter_direction(delta_dt, cfg, "up")
        down_dt <- .module3_filter_direction(delta_dt, cfg, "down")
        up_parts[[k]][[j]] <- up_dt
        down_parts[[k]][[j]] <- down_dt
        qc_parts[[k]][[j]] <- data.table::data.table(comparison_id = info$stem, comparison_label = info$comparison_label, comparison_group = info$comparison_group, cond1_id = info$cond1_id, cond2_id = info$cond2_id, cond1_label = info$cond1_label, cond2_label = info$cond2_label, cond1_matrix_id = info$cond1_matrix_id, cond2_matrix_id = info$cond2_matrix_id, chunk_id = j, n_links = n_links, n_prepared = n_prepared, n_delta = nrow(delta_dt), n_up = nrow(up_dt), n_down = nrow(down_dt))
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
      up[, `:=`(comparison_id = info$stem, comparison_label = info$comparison_label, comparison_group = info$comparison_group, cond1_id = info$cond1_id, cond2_id = info$cond2_id, cond1_label = info$cond1_label, cond2_label = info$cond2_label, cond1_matrix_id = info$cond1_matrix_id, cond2_matrix_id = info$cond2_matrix_id)]
      down[, `:=`(comparison_id = info$stem, comparison_label = info$comparison_label, comparison_group = info$comparison_group, cond1_id = info$cond1_id, cond2_id = info$cond2_id, cond1_label = info$cond1_label, cond2_label = info$cond2_label, cond1_matrix_id = info$cond1_matrix_id, cond2_matrix_id = info$cond2_matrix_id)]
      data.table::fwrite(up, info$up_path)
      data.table::fwrite(down, info$down_path)
      qc_rows[[length(qc_rows) + 1L]] <- data.table::rbindlist(qc_parts[[k]], use.names = TRUE, fill = TRUE)
      res[[i]] <- tibble::tibble(comparison_id = info$stem, comparison_label = info$comparison_label, comparison_group = info$comparison_group, cond1_id = info$cond1_id, cond2_id = info$cond2_id, cond1_label = info$cond1_label, cond2_label = info$cond2_label, cond1_matrix_id = info$cond1_matrix_id, cond2_matrix_id = info$cond2_matrix_id, up_path = info$up_path, down_path = info$down_path, n_up = nrow(up), n_down = nrow(down), fp_signal_mode = fp_signal_mode, skipped = FALSE)
    }
  }
  manifest <- dplyr::bind_rows(res)
  qc_chunks <- if (length(qc_rows)) {
    data.table::rbindlist(qc_rows, use.names = TRUE, fill = TRUE)
  } else {
    data.table::data.table(
      comparison_id = character(),
      comparison_label = character(),
      comparison_group = character(),
      cond1_id = character(),
      cond2_id = character(),
      cond1_label = character(),
      cond2_label = character(),
      cond1_matrix_id = character(),
      cond2_matrix_id = character(),
      chunk_id = integer(),
      n_links = integer(),
      n_prepared = integer(),
      n_delta = integer(),
      n_up = integer(),
      n_down = integer()
    )
  }
  qc_summary <- if (nrow(qc_chunks)) {
    by_cols <- c("comparison_id", "comparison_label", "comparison_group", "cond1_id", "cond2_id", "cond1_label", "cond2_label", "cond1_matrix_id", "cond2_matrix_id")
    qc_chunks[, .(
      n_chunks = .N,
      n_links = sum(get("n_links")),
      n_prepared = sum(get("n_prepared")),
      n_delta = sum(get("n_delta")),
      n_up = sum(get("n_up")),
      n_down = sum(get("n_down"))
    ), by = by_cols]
  } else {
    manifest_dt <- data.table::as.data.table(manifest)
    data.table::data.table(
      comparison_id = manifest_dt[["comparison_id"]],
      comparison_label = manifest_dt[["comparison_label"]],
      comparison_group = manifest_dt[["comparison_group"]],
      cond1_id = manifest_dt[["cond1_id"]],
      cond2_id = manifest_dt[["cond2_id"]],
      cond1_label = manifest_dt[["cond1_label"]],
      cond2_label = manifest_dt[["cond2_label"]],
      cond1_matrix_id = manifest_dt[["cond1_matrix_id"]],
      cond2_matrix_id = manifest_dt[["cond2_matrix_id"]],
      n_chunks = 0L,
      n_links = 0,
      n_prepared = 0,
      n_delta = 0,
      n_up = manifest_dt[["n_up"]],
      n_down = manifest_dt[["n_down"]]
    )
  }
  data.table::fwrite(qc_chunks, file.path(qc_dir, "differential_link_chunks.csv"))
  data.table::fwrite(qc_summary, file.path(qc_dir, "differential_link_summary.csv"))
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
