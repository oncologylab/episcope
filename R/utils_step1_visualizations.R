# File: utils_step1_visualizations.R
# Purpose: Public Module 1 TFBS visualization helpers.

.module1_relevant_config <- function(config) {
  if (is.null(config)) return(list())
  if (is.character(config) && length(config) == 1L && file.exists(config)) {
    config <- yaml::read_yaml(config)
  }
  if (!is.list(config)) return(list())
  keep <- c(
    "project", "project_name", "project_date", "base_dir", "source_project", "multiomic_object", "db", "ref_genome",
    "label_col", "fp_root_dir", "fp_cache_dir", "fp_cache_tag", "atac_master",
    "atac_data_path", "rna_mapped", "rna_path", "sample_metadata", "metadata_path",
    "threshold_expr", "threshold_gene_expr", "threshold_fp_score",
    "threshold_fp_tf_corr_r", "threshold_fp_tf_corr_p", "threshold_fp_tf_corr_fdr",
    "filter_to_canonical_bound", "footprint_sample_scope", "mid_slop",
    "round_digits", "score_match_pct", "output_mode", "module1_default_condition",
    "fp_score_raw", "fp_score_raw_path"
  )
  out <- config[intersect(keep, names(config))]
  secret <- grepl("token|password|secret|credential|api[_-]?key", names(out), ignore.case = TRUE)
  out[secret] <- "[REDACTED]"
  out
}

.module1_tfbs_source <- function(module1) {
  if (is.data.frame(module1)) return(list(data = tibble::as_tibble(module1)))
  if (is.list(module1) && !is.data.frame(module1)) {
    if (is.data.frame(module1$predicted_tfbs) && nrow(module1$predicted_tfbs)) {
      return(list(data = tibble::as_tibble(module1$predicted_tfbs)))
    }
    manifest <- module1$reports$predicted_tfbs_manifest %||%
      module1$predicted_tfbs_paths$manifest
    if (!is.null(manifest) && file.exists(manifest)) {
      return(list(manifest = .qc_read_manifest_chunks(manifest)))
    }
  }
  if (is.character(module1) && length(module1) == 1L) {
    path <- path.expand(module1)
    if (dir.exists(path)) {
      manifest <- file.path(path, "module1_predicted_tfbs_manifest.csv")
      if (file.exists(manifest)) return(list(manifest = .qc_read_manifest_chunks(manifest)))
      .log_abort("No predicted TFBS manifest was found in `module1`.")
    }
    if (!file.exists(path)) .log_abort("Module 1 TFBS source does not exist: {path}")
    candidate <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
    if (is.data.frame(candidate) && all(c("path", "format") %in% names(candidate))) {
      return(list(manifest = candidate))
    }
    return(list(data = .qc_read_table_file(path, columns = c("tf", "fp_id"))))
  }
  .log_abort("`module1` must be a predict_tfbs() result, output directory, TFBS table, table path, or manifest path.")
}

.module1_tfbs_each <- function(source, columns = c("tf", "fp_id"), callback) {
  if (is.data.frame(source$data)) {
    callback(source$data[, intersect(columns, names(source$data)), drop = FALSE], 1L)
    return(invisible(NULL))
  }
  manifest <- source$manifest
  if (!is.data.frame(manifest) || !nrow(manifest)) return(invisible(NULL))
  for (i in seq_len(nrow(manifest))) {
    callback(
      .qc_read_table_file(as.character(manifest$path[[i]]), as.character(manifest$format[[i]]), columns = columns),
      i
    )
  }
  invisible(NULL)
}

.module1_tfbs_pairs <- function(module1, tfs = NULL, fp_ids = NULL) {
  source <- .module1_tfbs_source(module1)
  parts <- list()
  .module1_tfbs_each(source, callback = function(x, i) {
    if (!all(c("tf", "fp_id") %in% names(x))) return()
    x <- x[, c("tf", "fp_id"), drop = FALSE]
    x$tf <- as.character(x$tf)
    x$fp_id <- as.character(x$fp_id)
    keep <- !is.na(x$tf) & nzchar(x$tf) & !is.na(x$fp_id) & nzchar(x$fp_id)
    if (!is.null(tfs)) keep <- keep & toupper(x$tf) %in% toupper(tfs)
    if (!is.null(fp_ids)) keep <- keep & x$fp_id %in% fp_ids
    parts[[length(parts) + 1L]] <<- unique(x[keep, , drop = FALSE])
  })
  if (!length(parts)) return(tibble::tibble(tf = character(), fp_id = character()))
  tibble::as_tibble(unique(data.table::rbindlist(parts, use.names = TRUE)))
}

.module1_tfbs_profile_counts <- function(module1) {
  fp_id <- n <- tf <- NULL
  source <- .module1_tfbs_source(module1)
  tf_parts <- list()
  fp_parts <- list()
  .module1_tfbs_each(source, callback = function(x, i) {
    if (!all(c("tf", "fp_id") %in% names(x))) return()
    dt <- unique(data.table::as.data.table(x[, c("tf", "fp_id"), drop = FALSE]))
    dt <- dt[!is.na(tf) & nzchar(tf) & !is.na(fp_id) & nzchar(fp_id)]
    tf_parts[[length(tf_parts) + 1L]] <<- dt[, .(n = .N), by = tf]
    fp_parts[[length(fp_parts) + 1L]] <<- dt[, .(n = .N), by = fp_id]
  })
  combine <- function(parts, key) {
    if (!length(parts)) return(data.table::data.table())
    data.table::rbindlist(parts, use.names = TRUE)[, .(n = sum(n)), by = key]
  }
  list(tf = combine(tf_parts, "tf"), fp = combine(fp_parts, "fp_id"))
}

.module1_resolve_multiomic <- function(module1, multiomic_data = NULL, project_config = NULL) {
  if (is_multiomic_object(multiomic_data)) return(multiomic_data)
  if (is.list(module1) && is_multiomic_object(module1$omics_data)) return(module1$omics_data)
  if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) {
    rds <- .qc_find_module1_omics_rds(module1, project_config = project_config)
    if (!is.na(rds) && file.exists(rds)) return(readRDS(rds))
  }
  .log_abort("`multiomic_data` is required and must be a CraftGRN multiomic object.")
}

.module1_safe_log2fc <- function(num, den, pseudocount = 1) {
  num <- suppressWarnings(as.numeric(num)) + pseudocount
  den <- suppressWarnings(as.numeric(den)) + pseudocount
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & num > 0 & den > 0
  out[ok] <- log2(num[ok] / den[ok])
  out
}

#' Plot TFBS changes between two conditions
#'
#' Summarizes active predicted TFBS and footprint signal for each TF. A TFBS is
#' active when its footprint is bound and its TF is expressed in the condition.
#'
#' @param module1 A `predict_tfbs()` result, Module 1 output directory,
#'   predicted TFBS table, table path, or manifest path.
#' @param multiomic_data Optional CraftGRN multiomic object. It may be omitted
#'   when it is present in `module1`.
#' @param cond1,cond2 Condition names. Changes are computed as cond1 relative to
#'   cond2.
#' @param pseudocount Positive pseudocount for log2 fold changes.
#' @param label_top_n Number of the most divergent TFs to label.
#' @param verbose Emit concise progress messages.
#' @return A ggplot object with the TF summary attached as `tfbs_summary`.
#' @export
plot_tfbs_condition_comparison <- function(module1,
                                           multiomic_data = NULL,
                                           cond1,
                                           cond2,
                                           pseudocount = 1,
                                           label_top_n = 15L,
                                           verbose = TRUE) {
  active1 <- active2 <- expr1 <- expr2 <- fp1 <- fp2 <- fp_id <- label <- max_tf_expression <- NULL
  fp_score_cond1 <- fp_score_cond2 <- n_tfbs_cond1 <- n_tfbs_cond2 <- tf <- NULL
  tf_expr_cond1 <- tf_expr_cond2 <- NULL
  omics <- .module1_resolve_multiomic(module1, multiomic_data)
  validate_multiomic_object(omics)
  cond1 <- as.character(cond1)[[1L]]
  cond2 <- as.character(cond2)[[1L]]
  conditions <- colnames(omics$matrices$fp_score)
  missing <- setdiff(c(cond1, cond2), conditions)
  if (length(missing)) .log_abort("Unknown condition(s): {paste(missing, collapse = ', ')}.")
  pseudocount <- suppressWarnings(as.numeric(pseudocount)[[1L]])
  if (!is.finite(pseudocount) || pseudocount <= 0) .log_abort("`pseudocount` must be positive.")

  gene_map <- stats::setNames(seq_len(nrow(omics$matrices$gene_expr)), toupper(rownames(omics$matrices$gene_expr)))
  c1 <- match(cond1, conditions)
  c2 <- match(cond2, conditions)
  gene_c1 <- match(cond1, colnames(omics$matrices$gene_expr))
  gene_c2 <- match(cond2, colnames(omics$matrices$gene_expr))
  partials <- list()
  source <- .module1_tfbs_source(module1)
  .module1_tfbs_each(source, callback = function(x, i) {
    if (!all(c("tf", "fp_id") %in% names(x)) || !nrow(x)) return()
    dt <- unique(data.table::as.data.table(x[, c("tf", "fp_id"), drop = FALSE]))
    dt <- dt[!is.na(tf) & nzchar(tf) & !is.na(fp_id) & nzchar(fp_id)]
    if (!nrow(dt)) return()
    fp_idx <- match(dt$fp_id, rownames(omics$matrices$fp_score))
    tf_idx <- unname(gene_map[toupper(dt$tf)])
    valid <- !is.na(fp_idx) & !is.na(tf_idx)
    dt[, `:=`(active1 = FALSE, active2 = FALSE, fp1 = NA_real_, fp2 = NA_real_, expr1 = NA_real_, expr2 = NA_real_)]
    if (any(valid)) {
      idx <- which(valid)
      dt$active1[idx] <- omics$matrices$fp_bound[cbind(fp_idx[idx], c1)] & omics$matrices$gene_on[cbind(tf_idx[idx], gene_c1)]
      dt$active2[idx] <- omics$matrices$fp_bound[cbind(fp_idx[idx], c2)] & omics$matrices$gene_on[cbind(tf_idx[idx], gene_c2)]
      dt$fp1[idx] <- omics$matrices$fp_score[cbind(fp_idx[idx], c1)]
      dt$fp2[idx] <- omics$matrices$fp_score[cbind(fp_idx[idx], c2)]
      dt$expr1[idx] <- omics$matrices$gene_expr[cbind(tf_idx[idx], gene_c1)]
      dt$expr2[idx] <- omics$matrices$gene_expr[cbind(tf_idx[idx], gene_c2)]
    }
    dt[, union_active := active1 | active2]
    partials[[length(partials) + 1L]] <<- dt[, .(
      n_tfbs_cond1 = sum(active1, na.rm = TRUE),
      n_tfbs_cond2 = sum(active2, na.rm = TRUE),
      fp_score_cond1 = sum(data.table::fifelse(union_active, fp1, 0), na.rm = TRUE),
      fp_score_cond2 = sum(data.table::fifelse(union_active, fp2, 0), na.rm = TRUE),
      tf_expr_cond1 = suppressWarnings(max(expr1, na.rm = TRUE)),
      tf_expr_cond2 = suppressWarnings(max(expr2, na.rm = TRUE))
    ), by = tf]
  })
  if (!length(partials)) .log_abort("No predicted TFBS pairs are available.")
  summary <- data.table::rbindlist(partials, use.names = TRUE)[, .(
    n_tfbs_cond1 = sum(n_tfbs_cond1),
    n_tfbs_cond2 = sum(n_tfbs_cond2),
    fp_score_cond1 = sum(fp_score_cond1),
    fp_score_cond2 = sum(fp_score_cond2),
    tf_expr_cond1 = suppressWarnings(max(tf_expr_cond1, na.rm = TRUE)),
    tf_expr_cond2 = suppressWarnings(max(tf_expr_cond2, na.rm = TRUE))
  ), by = tf]
  summary[!is.finite(tf_expr_cond1), tf_expr_cond1 := NA_real_]
  summary[!is.finite(tf_expr_cond2), tf_expr_cond2 := NA_real_]
  summary[, `:=`(
    delta_n_predicted_tfbs = n_tfbs_cond1 - n_tfbs_cond2,
    log2fc_fp_score = .module1_safe_log2fc(fp_score_cond1, fp_score_cond2, pseudocount),
    max_tf_expression = pmax(tf_expr_cond1, tf_expr_cond2, na.rm = TRUE),
    log2fc_tf_expression = .module1_safe_log2fc(tf_expr_cond1, tf_expr_cond2, pseudocount)
  )]
  summary[!is.finite(max_tf_expression), max_tf_expression := NA_real_]
  score_x <- abs(summary$delta_n_predicted_tfbs) / max(1, max(abs(summary$delta_n_predicted_tfbs), na.rm = TRUE))
  score_y <- abs(summary$log2fc_fp_score) / max(1e-12, max(abs(summary$log2fc_fp_score), na.rm = TRUE))
  summary[, label := ""]
  label_top_n <- max(0L, as.integer(label_top_n)[[1L]])
  if (label_top_n > 0L) {
    label_idx <- head(order(score_x + score_y, decreasing = TRUE, na.last = NA), label_top_n)
    summary$label[label_idx] <- summary$tf[label_idx]
  }

  p <- ggplot2::ggplot(summary, ggplot2::aes(
    x = .data$delta_n_predicted_tfbs,
    y = .data$log2fc_fp_score,
    size = .data$max_tf_expression,
    color = .data$log2fc_tf_expression
  )) +
    ggplot2::geom_hline(yintercept = 0, color = "#94a3b8", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, color = "#94a3b8", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.82) +
    ggplot2::scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, name = "TF expression\nlog2FC") +
    ggplot2::scale_size_area(max_size = 11, name = "Max TF\nexpression") +
    ggplot2::labs(
      title = paste(cond1, "v", cond2),
      x = "delta N of predicted TFBS",
      y = "log2FC of aggregated FP scores"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  labels <- summary[nzchar(label)]
  if (nrow(labels)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = .data$label), size = 3, show.legend = FALSE, max.overlaps = Inf)
    } else {
      p <- p + ggplot2::geom_text(data = labels, ggplot2::aes(label = .data$label), size = 3, check_overlap = TRUE, show.legend = FALSE)
    }
  }
  attr(p, "tfbs_summary") <- tibble::as_tibble(summary[, setdiff(names(summary), "label"), with = FALSE])
  if (isTRUE(verbose)) .log_inform("Module 1 TFBS condition comparison summarized {nrow(summary)} TF(s).")
  p
}

.module1_jsd_matrix <- function(x) {
  counts <- as.numeric(Matrix::rowSums(x > 0))
  overlap <- as.matrix(Matrix::tcrossprod(x > 0))
  n <- nrow(overlap)
  out <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  for (i in seq_len(n)) {
    for (j in i:n) {
      a <- counts[[i]]
      b <- counts[[j]]
      shared <- overlap[i, j]
      if (a <= 0 || b <= 0) {
        value <- if (a == b) 0 else 1
      } else {
        kl1 <- (a - shared) / a + (shared / a) * log2(2 * b / (a + b))
        kl2 <- (b - shared) / b + (shared / b) * log2(2 * a / (a + b))
        value <- sqrt(max(0, 0.5 * (kl1 + kl2)))
      }
      out[i, j] <- value
      out[j, i] <- value
    }
  }
  out
}

#' Plot a TF-TF co-binding heatmap
#'
#' @param module1 A Module 1 result, output directory, predicted TFBS table,
#'   table path, or manifest path.
#' @param metric Heatmap metric: absolute shared-site counts or Jensen-Shannon
#'   distance (`jsd`).
#' @param max_tfs Maximum number of TFs, ranked by predicted TFBS count.
#' @param verbose Emit concise progress messages.
#' @return A ggplot heatmap with `cobinding_matrix` and `tf_order` attributes.
#' @export
plot_tf_tf_cobinding_heatmap <- function(module1,
                                         metric = c("absolute", "jsd"),
                                         max_tfs = 100L,
                                         verbose = TRUE) {
  n <- tf <- NULL
  metric <- match.arg(metric)
  profile_counts <- .module1_tfbs_profile_counts(module1)
  if (!nrow(profile_counts$tf)) .log_abort("No predicted TFBS pairs are available.")
  counts <- profile_counts$tf[order(-n, tf)]
  max_tfs <- suppressWarnings(as.numeric(max_tfs)[[1L]])
  if (!is.finite(max_tfs) || max_tfs < 2) .log_abort("`max_tfs` must be at least 2.")
  selected <- head(as.character(counts$tf), as.integer(max_tfs))
  pairs <- .module1_tfbs_pairs(module1, tfs = selected)
  tf_levels <- sort(unique(as.character(pairs$tf)))
  fp_levels <- sort(unique(as.character(pairs$fp_id)))
  if (length(tf_levels) < 2L) .log_abort("At least two TFs with predicted binding are required.")
  if (!requireNamespace("Matrix", quietly = TRUE)) .log_abort("Package `Matrix` is required for TF-TF co-binding heatmaps.")
  incidence <- Matrix::sparseMatrix(
    i = match(pairs$tf, tf_levels),
    j = match(pairs$fp_id, fp_levels),
    x = 1,
    dims = c(length(tf_levels), length(fp_levels)),
    dimnames = list(tf_levels, fp_levels)
  )
  incidence@x[] <- 1
  absolute <- as.matrix(Matrix::tcrossprod(incidence))
  value_matrix <- if (identical(metric, "absolute")) absolute else .module1_jsd_matrix(incidence)
  cluster_input <- if (identical(metric, "absolute")) stats::dist(log1p(value_matrix)) else stats::as.dist(value_matrix)
  order_idx <- stats::hclust(cluster_input, method = "complete")$order
  tf_order <- rownames(value_matrix)[order_idx]
  long <- as.data.frame(as.table(value_matrix), stringsAsFactors = FALSE)
  names(long) <- c("tf_row", "tf_col", "value")
  long$tf_row <- factor(long$tf_row, levels = rev(tf_order))
  long$tf_col <- factor(long$tf_col, levels = tf_order)
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$tf_col, .data$tf_row, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = if (identical(metric, "absolute")) "TF-TF co-binding: absolute shared TFBS" else "TF-TF co-binding: Jensen-Shannon distance",
      x = "TF",
      y = "TF",
      fill = if (identical(metric, "absolute")) "Shared TFBS" else "JSD"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
  if (identical(metric, "absolute")) {
    p <- p + ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08306b", trans = "log1p")
  } else {
    p <- p + ggplot2::scale_fill_gradient(low = "#fff7ec", high = "#7f0000", limits = c(0, 1))
  }
  attr(p, "cobinding_matrix") <- value_matrix
  attr(p, "tf_order") <- tf_order
  if (isTRUE(verbose)) .log_inform("Module 1 co-binding heatmap includes {length(tf_order)} TF(s).")
  p
}

.module1_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

#' Build an interactive TF UMAP report from predicted TFBS
#'
#' @param module1 A Module 1 result, output directory, predicted TFBS table,
#'   table path, or manifest path.
#' @param output_file HTML output path.
#' @param top_variable_tfbs Maximum number of variable footprint features.
#' @param cluster_range Integer cluster counts offered in the HTML selector.
#' @param default_clusters Initial cluster count.
#' @param seed Random seed for UMAP and k-means.
#' @param verbose Emit concise progress messages.
#' @return Normalized path to the HTML report.
#' @export
build_tfbs_umap_report <- function(module1,
                                   output_file,
                                   top_variable_tfbs = 2000L,
                                   cluster_range = 2:12,
                                   default_clusters = 5L,
                                   seed = 1L,
                                   verbose = TRUE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) .log_abort("Package `Matrix` is required to build the TFBS UMAP report.")
  if (!requireNamespace("uwot", quietly = TRUE)) .log_abort("Package `uwot` is required to build the TFBS UMAP report. Install it with install.packages('uwot').")
  if (!is.character(output_file) || length(output_file) != 1L || !nzchar(output_file)) .log_abort("`output_file` must be a non-empty HTML path.")
  profile_counts <- .module1_tfbs_profile_counts(module1)
  if (!nrow(profile_counts$tf) || !nrow(profile_counts$fp)) .log_abort("No predicted TFBS pairs are available.")
  tf_levels <- sort(as.character(profile_counts$tf$tf))
  if (length(tf_levels) < 3L) .log_abort("At least three TFs are required for a TFBS UMAP report.")
  site_n <- profile_counts$fp$n
  variance <- (as.numeric(site_n) / length(tf_levels)) * (1 - as.numeric(site_n) / length(tf_levels))
  names(variance) <- as.character(profile_counts$fp$fp_id)
  variance <- variance[is.finite(variance) & variance > 0]
  if (!length(variance)) .log_abort("No variable TFBS features are available for UMAP.")
  top_variable_tfbs <- max(2L, as.integer(top_variable_tfbs)[[1L]])
  selected_sites <- names(head(sort(variance, decreasing = TRUE), top_variable_tfbs))
  pairs <- .module1_tfbs_pairs(module1, fp_ids = selected_sites)
  site_levels <- sort(selected_sites)
  incidence <- Matrix::sparseMatrix(
    i = match(pairs$tf, tf_levels),
    j = match(pairs$fp_id, site_levels),
    x = 1,
    dims = c(length(tf_levels), length(site_levels)),
    dimnames = list(tf_levels, site_levels)
  )
  incidence@x[] <- 1
  tf_norm <- Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(incidence), 1))
  idf <- log1p(length(tf_levels) / (1 + Matrix::colSums(incidence > 0)))
  tfidf <- tf_norm %*% incidence %*% Matrix::Diagonal(x = as.numeric(idf))
  dense <- as.matrix(tfidf)
  l2 <- sqrt(rowSums(dense^2))
  dense <- dense / pmax(l2, 1e-12)
  rank_use <- max(2L, min(50L, nrow(dense) - 1L, ncol(dense)))
  pca <- stats::prcomp(dense, center = TRUE, scale. = FALSE, rank. = rank_use)
  scores <- pca$x[, seq_len(min(rank_use, ncol(pca$x))), drop = FALSE]
  seed <- as.integer(seed)[[1L]]
  coords <- .module1_with_seed(seed, uwot::umap(
    scores,
    n_neighbors = min(15L, nrow(scores) - 1L),
    n_components = 2L,
    metric = "cosine",
    verbose = FALSE
  ))
  rownames(coords) <- tf_levels
  distinct_n <- nrow(unique(as.data.frame(signif(scores, 12))))
  cluster_range <- sort(unique(as.integer(cluster_range)))
  cluster_range <- cluster_range[is.finite(cluster_range) & cluster_range >= 2L & cluster_range <= min(nrow(scores) - 1L, distinct_n)]
  if (!length(cluster_range)) .log_abort("`cluster_range` does not contain a valid cluster count for these TF profiles.")
  assignments <- list()
  for (k in cluster_range) {
    fit <- tryCatch(.module1_with_seed(seed + k, stats::kmeans(scores, centers = k, nstart = 20L)), error = function(e) NULL)
    if (!is.null(fit)) assignments[[as.character(k)]] <- as.integer(fit$cluster)
  }
  if (!length(assignments)) .log_abort("K-means clustering failed for every requested cluster count.")
  valid_k <- as.integer(names(assignments))
  default_clusters <- as.integer(default_clusters)[[1L]]
  if (!default_clusters %in% valid_k) default_clusters <- valid_k[[which.min(abs(valid_k - default_clusters))]]
  tf_counts <- stats::setNames(profile_counts$tf$n, profile_counts$tf$tf)
  payload <- lapply(seq_along(tf_levels), function(i) list(
    tf = tf_levels[[i]],
    x = unname(coords[i, 1L]),
    y = unname(coords[i, 2L]),
    n_tfbs = unname(as.integer(tf_counts[tf_levels[[i]]])),
    clusters = stats::setNames(lapply(assignments, `[[`, i), names(assignments))
  ))
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = NA, null = "null")
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  options_html <- paste(vapply(valid_k, function(k) sprintf("<option value=\"%d\"%s>%d</option>", k, if (k == default_clusters) " selected" else "", k), character(1L)), collapse = "")
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Module 1 TFBS UMAP</title>",
    "<style>body{font-family:Arial,Helvetica,sans-serif;margin:0;color:#17202a}header{padding:18px 24px;border-bottom:1px solid #d8dee8}h1{font-size:25px;margin:0 0 10px}label{font-weight:700}select{margin-left:8px;padding:5px}main{padding:16px 24px}.meta{color:#5d6673;margin-left:16px}svg{width:100%;height:auto;border:1px solid #d8dee8;background:#fff}.axis{stroke:#94a3b8}.dot{stroke:#fff;stroke-width:1;opacity:.88}#legend{pointer-events:none}.tip{position:fixed;display:none;background:#111827;color:#fff;padding:7px 9px;border-radius:3px;pointer-events:none;font-size:12px}</style></head><body>",
    "<header><h1>TF UMAP from highly variable predicted TFBS</h1><label>Number of clusters<select id=\"k\">", options_html, "</select></label><span id=\"meta\" class=\"meta\"></span></header>",
    "<main><svg id=\"plot\" viewBox=\"0 0 1100 760\"><g id=\"layer\"></g></svg></main><div id=\"tip\" class=\"tip\"></div>",
    "<script>const D=", json, ";const sel=document.getElementById('k'),layer=document.getElementById('layer'),tip=document.getElementById('tip'),meta=document.getElementById('meta');",
    "const colors=['#2166ac','#b2182b','#1b9e77','#762a83','#e08214','#4d9221','#c51b7d','#0571b0','#7f7f7f','#a6761d','#66a61e','#e7298a'];",
    "function svgEl(n,a={}){const e=document.createElementNS('http://www.w3.org/2000/svg',n);Object.entries(a).forEach(([k,v])=>e.setAttribute(k,v));return e}function scale(v,a,b,c,d){return c+(v-a)/(b-a||1)*(d-c)}function draw(){layer.replaceChildren();const k=sel.value,xv=D.map(d=>d.x),yv=D.map(d=>d.y),xmin=Math.min(...xv),xmax=Math.max(...xv),ymin=Math.min(...yv),ymax=Math.max(...yv);D.forEach(d=>{const c=svgEl('circle',{cx:scale(d.x,xmin,xmax,55,1045),cy:scale(d.y,ymin,ymax,705,45),r:Math.max(4,Math.min(12,4+Math.log1p(d.n_tfbs))),fill:colors[(d.clusters[k]-1)%colors.length],class:'dot'});c.addEventListener('mousemove',e=>{tip.style.display='block';tip.style.left=(e.clientX+12)+'px';tip.style.top=(e.clientY+12)+'px';tip.textContent=d.tf+' | Cluster '+d.clusters[k]+' | predicted TFBS '+d.n_tfbs});c.addEventListener('mouseout',()=>tip.style.display='none');layer.appendChild(c)});const g=svgEl('g',{id:'legend'}),bg=svgEl('rect',{x:920,y:20,width:150,height:30+Number(k)*25,fill:'#fff',stroke:'#d8dee8',rx:3});g.appendChild(bg);const title=svgEl('text',{x:936,y:43,'font-size':15,'font-weight':700});title.textContent='Cluster';g.appendChild(title);for(let i=1;i<=Number(k);i++){g.appendChild(svgEl('circle',{cx:942,cy:47+i*24,r:6,fill:colors[(i-1)%colors.length]}));const t=svgEl('text',{x:958,y:52+i*24,'font-size':13});t.textContent=String(i);g.appendChild(t)}layer.appendChild(g);meta.textContent=D.length+' TFs | K='+k}sel.addEventListener('change',draw);draw();</script></body></html>"
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  path <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  if (isTRUE(verbose)) .log_inform("Module 1 TFBS UMAP HTML written: {path}")
  path
}

.module1_pack_site_indices <- function(indices, n_sites) {
  out <- raw(ceiling(as.integer(n_sites) / 8))
  indices <- unique(suppressWarnings(as.integer(indices)))
  indices <- indices[is.finite(indices) & indices >= 1L & indices <= n_sites]
  if (!length(indices)) return(out)
  byte <- (indices - 1L) %/% 8L + 1L
  mask <- bitwShiftL(1L, (indices - 1L) %% 8L)
  packed <- data.table::data.table(byte = byte, mask = mask)[, .(mask = sum(unique(mask))), by = byte]
  out[packed$byte] <- as.raw(packed$mask)
  out
}

.module1_unpack_site_indices <- function(bits, n_sites) {
  if (!length(bits) || n_sites < 1L) return(integer())
  byte_index <- rep(seq_along(bits), each = 8L)
  bit_index <- rep(0:7, times = length(bits))
  keep <- bitwAnd(rep(as.integer(bits), each = 8L), bitwShiftL(1L, bit_index)) > 0L
  out <- (byte_index[keep] - 1L) * 8L + bit_index[keep] + 1L
  out[out <= n_sites]
}

.module1_or_packed_indices <- function(bits, indices) {
  indices <- unique(suppressWarnings(as.integer(indices)))
  indices <- indices[is.finite(indices) & indices >= 1L & indices <= length(bits) * 8L]
  if (!length(indices)) return(bits)
  byte <- (indices - 1L) %/% 8L + 1L
  mask <- bitwShiftL(1L, (indices - 1L) %% 8L)
  packed <- data.table::data.table(byte = byte, mask = mask)[, .(mask = sum(unique(mask))), by = byte]
  bits[packed$byte] <- as.raw(bitwOr(as.integer(bits[packed$byte]), packed$mask))
  bits
}

.module1_explorer_cache_path <- function(module1, output_file) {
  if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) {
    return(file.path(module1, "cache", "module1_qc_analysis.rds"))
  }
  if (is.list(module1) && !is.data.frame(module1)) {
    manifest <- module1$reports$predicted_tfbs_manifest %||% module1$predicted_tfbs_paths$manifest
    if (!is.null(manifest) && length(manifest) == 1L && nzchar(manifest)) {
      return(file.path(dirname(manifest), "cache", "module1_qc_analysis.rds"))
    }
  }
  file.path(dirname(output_file), "module1_qc_analysis.rds")
}

.module1_explorer_motif_map <- function(omics) {
  x <- omics$features$fp_motif
  if (!is.data.frame(x) || !all(c("fp_id", "motif") %in% names(x))) {
    return(data.table::data.table(fp_id = character(), motif = character(), canonical_tf = character()))
  }
  canonical <- if ("tf" %in% names(x)) as.character(x$tf) else rep("", nrow(x))
  out <- data.table::data.table(fp_id = as.character(x$fp_id), motif = as.character(x$motif), canonical_tf = canonical)
  unique(out[!is.na(fp_id) & nzchar(fp_id) & !is.na(motif) & nzchar(motif)])
}

.module1_tf_identity_key <- function(x) {
  key <- toupper(trimws(as.character(x)))
  key[key == "TBET"] <- "TBX21"
  key
}

.module1_explorer_correlation_hist <- function(module1, bins = 40L, r_cutoff = NULL) {
  parts <- list()
  cutoff <- if (!is.null(r_cutoff)) {
    suppressWarnings(as.numeric(r_cutoff))
  } else if (is.list(module1) && !is.data.frame(module1) && is.list(module1$parameters)) {
    suppressWarnings(as.numeric(module1$parameters$r_cutoff %||% 0.5))
  } else {
    0.5
  }
  if (!is.finite(cutoff)) cutoff <- 0.5
  add_table <- function(x, scope) {
    if (!is.data.frame(x) || !nrow(x) || !"tf" %in% names(x)) return()
    for (method in c("pearson_r", "spearman_r")) {
      if (!method %in% names(x)) next
      pass <- suppressWarnings(as.numeric(x[[method]])) >= cutoff
      dt <- data.table::data.table(
        tf = as.character(x$tf),
        value = suppressWarnings(as.numeric(x[[method]])),
        pass = pass
      )
      dt <- dt[is.finite(value) & value >= -1 & value <= 1]
      if (!nrow(dt)) next
      dt[, bin := pmin(bins - 1L, pmax(0L, floor((value + 1) / 2 * bins)))]
      per_tf <- dt[, .(n = .N, n_pass = sum(pass)), by = .(tf, bin)]
      overall <- dt[, .(n = .N, n_pass = sum(pass)), by = bin][, tf := "Overall"]
      out <- data.table::rbindlist(list(overall, per_tf), use.names = TRUE)
      out[, `:=`(scope = scope, method = method)]
      parts[[length(parts) + 1L]] <<- out[, .(scope, method, tf, bin, n, n_pass)]
    }
  }
  if (is.list(module1) && !is.data.frame(module1)) {
    add_table(module1$motif_supported_correlations, "Canonical motif-supported")
    add_table(module1$prediction_stats, "Full prediction")
    manifest <- module1$prediction_stats_manifest
    if (is.data.frame(manifest) && nrow(manifest)) {
      for (i in seq_len(nrow(manifest))) {
        add_table(
          .qc_read_table_file(as.character(manifest$path[[i]]), as.character(manifest$format[[i]]), columns = c("tf", "pearson_r", "spearman_r")),
          "Full prediction"
        )
      }
    }
  }
  module1_dir <- if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) module1 else NULL
  if (!is.null(module1_dir)) {
    canonical_path <- file.path(module1_dir, "module1_canonical_prediction_stats.csv.gz")
    if (file.exists(canonical_path)) add_table(.qc_read_table_file(canonical_path), "Canonical motif-supported")
    stats_manifest_path <- file.path(module1_dir, "module1_prediction_stats_chunks", "module1_prediction_stats_manifest.csv")
    if (file.exists(stats_manifest_path)) {
      manifest <- .qc_read_manifest_chunks(stats_manifest_path)
      for (i in seq_len(nrow(manifest))) {
        add_table(
          .qc_read_table_file(as.character(manifest$path[[i]]), as.character(manifest$format[[i]]), columns = c("tf", "pearson_r", "spearman_r")),
          "Full prediction"
        )
      }
    }
  }
  if (!length(parts)) return(tibble::tibble())
  data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)[
    , .(n = sum(n), n_pass = sum(n_pass, na.rm = TRUE)), by = .(scope, method, tf, bin)
  ] |>
    tibble::as_tibble()
}

.module1_explorer_r_cutoff <- function(module1, project_config = NULL, omics = NULL) {
  candidates <- list()
  if (is.list(module1) && !is.data.frame(module1) && is.list(module1$parameters)) {
    candidates <- c(candidates, list(module1$parameters$r_cutoff))
  }
  if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) {
    path <- file.path(module1, "module1_run_parameters.csv")
    if (file.exists(path)) {
      saved <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) tibble::tibble())
      if (nrow(saved) && all(c("parameter", "effective_value") %in% names(saved))) {
        candidates <- c(candidates, list(saved$effective_value[match("r_cutoff", saved$parameter)]))
      }
    }
  }
  config <- .module1_relevant_config(project_config)
  candidates <- c(
    candidates,
    list(
      config$threshold_fp_tf_corr_r,
      if (is_multiomic_object(omics)) omics$project$config$threshold_fp_tf_corr_r else NULL,
      0.5
    )
  )
  for (candidate in candidates) {
    value <- suppressWarnings(as.numeric(candidate))
    if (length(value) && is.finite(value[[1L]])) return(value[[1L]])
  }
  0.5
}

.module1_build_explorer_cache <- function(module1, omics, cache_file, project_config = NULL, rebuild = FALSE, verbose = TRUE) {
  validate_multiomic_object(omics)
  r_cutoff <- .module1_explorer_r_cutoff(module1, project_config = project_config, omics = omics)
  source <- .module1_tfbs_source(module1)
  manifest_rows <- if (is.data.frame(source$manifest) && "n_rows" %in% names(source$manifest)) {
    sum(suppressWarnings(as.numeric(source$manifest$n_rows)), na.rm = TRUE)
  } else if (is.data.frame(source$data)) nrow(source$data) else NA_real_
  fingerprint <- list(
    schema = "module1_qc_analysis_v8",
    n_sites = nrow(omics$matrices$fp_score),
    conditions = colnames(omics$matrices$fp_score),
    r_cutoff = r_cutoff,
    fp_ids_digest = digest::digest(rownames(omics$matrices$fp_score), algo = "xxhash64"),
    fp_bound_digest = digest::digest(omics$matrices$fp_bound, algo = "xxhash64"),
    gene_on_digest = digest::digest(omics$matrices$gene_on, algo = "xxhash64"),
    motif_digest = digest::digest(omics$features$fp_motif, algo = "xxhash64"),
    manifest_rows = manifest_rows,
    source_files = if (is.data.frame(source$manifest) && "path" %in% names(source$manifest)) {
      info <- file.info(as.character(source$manifest$path))
      data.frame(
        path = normalizePath(as.character(source$manifest$path), winslash = "/", mustWork = FALSE),
        size = as.numeric(info$size),
        mtime = format(info$mtime, tz = "UTC", usetz = TRUE),
        stringsAsFactors = FALSE
      )
    } else NULL
  )
  if (!isTRUE(rebuild) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.list(cached) && identical(cached$fingerprint, fingerprint)) return(cached)
  }
  if (isTRUE(verbose)) .log_inform("Building compact Module 1 TFBS explorer cache.")
  fp_ids <- rownames(omics$matrices$fp_score)
  fp_index <- stats::setNames(seq_along(fp_ids), fp_ids)
  n_sites <- length(fp_ids)
  n_bytes <- ceiling(n_sites / 8)
  tf_bits <- list()
  motif_map <- .module1_explorer_motif_map(omics)
  canonical_by_motif <- if (nrow(motif_map)) {
    motif_map[, .(canonical_tf = paste(sort(unique(canonical_tf[!is.na(canonical_tf) & nzchar(canonical_tf)])), collapse = ";")), by = motif]
  } else {
    data.table::data.table(motif = character(), canonical_tf = character())
  }
  canonical_members <- if (nrow(motif_map)) {
    members <- unique(motif_map[!is.na(canonical_tf) & nzchar(canonical_tf), .(
      motif,
      canonical_display = canonical_tf,
      tf_key = .module1_tf_identity_key(canonical_tf)
    )])
    members[, .(canonical_display = paste(sort(unique(canonical_display)), collapse = ";")), by = .(motif, tf_key)]
  } else {
    data.table::data.table(motif = character(), tf_key = character(), canonical_display = character())
  }
  motif_sites <- if (nrow(motif_map)) unique(motif_map[, .(fp_id, motif)]) else motif_map[, .(fp_id, motif)]
  motif_parts <- list()
  .module1_tfbs_each(source, columns = c("tf", "fp_id"), callback = function(x, i) {
    if (!all(c("tf", "fp_id") %in% names(x)) || !nrow(x)) return()
    dt <- unique(data.table::data.table(tf = as.character(x$tf), fp_id = as.character(x$fp_id)))
    dt <- dt[!is.na(tf) & nzchar(tf) & !is.na(fp_id) & nzchar(fp_id)]
    dt[, site_index := unname(fp_index[fp_id])]
    dt <- dt[!is.na(site_index)]
    if (!nrow(dt)) return()
    by_tf <- split(dt$site_index, dt$tf)
    for (tf_name in names(by_tf)) {
      bits <- tf_bits[[tf_name]]
      if (is.null(bits)) bits <- raw(n_bytes)
      tf_bits[[tf_name]] <<- .module1_or_packed_indices(bits, by_tf[[tf_name]])
    }
    if (nrow(motif_sites)) {
      joined <- motif_sites[dt[, .(fp_id, tf)], on = "fp_id", nomatch = 0L, allow.cartesian = TRUE]
      if (nrow(joined)) motif_parts[[length(motif_parts) + 1L]] <<- joined[, .(n = .N), by = .(motif, tf)]
    }
  })
  tf_names <- sort(names(tf_bits))
  tf_bits <- tf_bits[tf_names]
  conditions <- colnames(omics$matrices$fp_bound)
  bound_bits <- lapply(seq_along(conditions), function(i) {
    .module1_pack_site_indices(which(omics$matrices$fp_bound[, i] %in% TRUE), n_sites)
  })
  names(bound_bits) <- conditions
  gene_index <- stats::setNames(seq_len(nrow(omics$matrices$gene_on)), toupper(rownames(omics$matrices$gene_on)))
  gene_on <- matrix(FALSE, nrow = length(tf_names), ncol = length(conditions), dimnames = list(tf_names, conditions))
  tf_idx <- unname(gene_index[toupper(tf_names)])
  cond_idx <- match(conditions, colnames(omics$matrices$gene_on))
  valid_tf <- which(!is.na(tf_idx))
  valid_cond <- which(!is.na(cond_idx))
  if (length(valid_tf) && length(valid_cond)) {
    gene_on[valid_tf, valid_cond] <- omics$matrices$gene_on[tf_idx[valid_tf], cond_idx[valid_cond], drop = FALSE]
  }
  popcount <- vapply(0:255, function(value) sum(as.integer(intToBits(value))[1:8]), integer(1L))
  tf_counts <- matrix(0, nrow = length(tf_names), ncol = length(conditions) + 1L, dimnames = list(tf_names, c("Overall", conditions)))
  for (i in seq_along(tf_names)) {
    tf_counts[i, "Overall"] <- sum(popcount[as.integer(tf_bits[[i]]) + 1L])
    for (j in seq_along(conditions)) {
      if (gene_on[i, j]) {
        tf_counts[i, j + 1L] <- sum(popcount[bitwAnd(as.integer(tf_bits[[i]]), as.integer(bound_bits[[j]])) + 1L])
      }
    }
  }
  motif_counts <- if (length(motif_parts)) {
    data.table::rbindlist(motif_parts, use.names = TRUE)[, .(n = sum(n)), by = .(motif, tf)]
  } else {
    data.table::data.table(motif = character(), tf = character(), n = integer())
  }
  if (nrow(motif_counts)) {
    motif_counts[, tf_key := .module1_tf_identity_key(tf)]
    data.table::setorderv(motif_counts, c("motif", "n", "tf"), order = c(1L, -1L, 1L))
    motif_counts[, rank := seq_len(.N), by = motif]
    motif_counts[, top20 := rank <= 20L]
    observed_keys <- paste(motif_counts$motif, motif_counts$tf_key, sep = "\r")
    canonical_members[, predicted := paste(motif, tf_key, sep = "\r") %in% observed_keys]
    eligible_motifs <- canonical_members[, .(eligible = length(predicted) > 0L && all(predicted)), by = motif][eligible == TRUE, motif]
    all_motifs <- sort(unique(canonical_members$motif))
    motif_summary <- list(
      total = length(all_motifs),
      eligible = length(eligible_motifs),
      excluded = length(setdiff(all_motifs, eligible_motifs))
    )
    motif_counts <- motif_counts[motif %in% eligible_motifs]
    canonical_members <- canonical_members[motif %in% eligible_motifs]
    canonical_members[, predicted := NULL]
    canonical_by_motif <- canonical_by_motif[motif %in% eligible_motifs]
    canonical_keys <- paste(canonical_members$motif, canonical_members$tf_key, sep = "\r")
    motif_counts[, is_canonical := paste(motif, tf_key, sep = "\r") %in% canonical_keys]
    motif_counts[, canonical_status := ifelse(is_canonical & top20, "top_20", "")]
    top_rows <- motif_counts[top20 == TRUE]
    canonical_rows <- merge(
      canonical_members,
      motif_counts[, .(motif, tf_key, tf, n, rank, top20, is_canonical)],
      by = c("motif", "tf_key"), all.x = TRUE, sort = FALSE
    )
    canonical_rows <- canonical_rows[is.na(top20) | !top20]
    if (nrow(canonical_rows)) {
      canonical_rows[is.na(tf) | !nzchar(tf), tf := canonical_display]
      canonical_rows[is.na(n), n := 0]
      canonical_rows[, `:=`(
        top20 = FALSE,
        is_canonical = TRUE,
        canonical_status = "predicted_outside_top20"
      )]
      canonical_rows <- canonical_rows[, .(motif, tf, n, tf_key, rank, top20, is_canonical, canonical_status)]
      motif_counts <- data.table::rbindlist(list(top_rows, canonical_rows), use.names = TRUE, fill = TRUE)
    } else {
      motif_counts <- top_rows
    }
    motif_counts <- canonical_by_motif[motif_counts, on = "motif"]
    data.table::setorderv(motif_counts, c("motif", "top20", "rank", "tf"), order = c(1L, -1L, 1L, 1L), na.last = TRUE)
  } else {
    motif_summary <- list(total = 0L, eligible = 0L, excluded = 0L)
    motif_counts <- data.table::data.table(
      motif = character(), canonical_tf = character(), tf = character(), n = integer(),
      tf_key = character(), rank = integer(), top20 = logical(), is_canonical = logical(), canonical_status = character()
    )
  }
  used_site_bits <- if (length(tf_bits)) {
    Reduce(function(a, b) as.raw(bitwOr(as.integer(a), as.integer(b))), tf_bits)
  } else {
    raw(n_bytes)
  }
  cache <- list(
    fingerprint = fingerprint,
    tf_names = tf_names,
    tf_bits = tf_bits,
    used_site_bits = used_site_bits,
    conditions = conditions,
    bound_bits = bound_bits,
    gene_on = gene_on,
    tf_counts = tf_counts,
    motif_counts = tibble::as_tibble(motif_counts),
    motif_summary = motif_summary
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, cache_file, compress = "xz")
  if (isTRUE(verbose)) .log_inform("Finished compact Module 1 TFBS explorer cache.")
  cache
}

.module1_explorer_default_condition <- function(conditions, default_condition = NULL, project_config = NULL, omics = NULL) {
  config <- .module1_relevant_config(project_config)
  candidate <- default_condition %||% config$module1_default_condition
  candidate <- as.character(candidate %||% "")[[1L]]
  if (nzchar(candidate) && candidate %in% conditions) return(candidate)
  meta <- if (is_multiomic_object(omics) && is.data.frame(omics$samples)) omics$samples else NULL
  if (is.data.frame(meta)) {
    value_cols <- intersect(c("condition", "stress_type", "sample", "name", "condition_id"), names(meta))
    key_cols <- intersect(c("condition_id", "sample_id", "id", "name", "sample"), names(meta))
    for (value_col in value_cols) {
      ctrl <- grepl("(^|[_ -])(ctrl|control)([_ -]|$)", as.character(meta[[value_col]]), ignore.case = TRUE)
      if (!any(ctrl, na.rm = TRUE)) next
      for (key_col in key_cols) {
        hits <- intersect(conditions, as.character(meta[[key_col]][ctrl]))
        if (length(hits)) return(sort(hits)[[1L]])
      }
    }
  }
  sort(conditions)[[1L]]
}

.module1_tfbs_explorer_components <- function(cache, default_condition, standalone = FALSE) {
  if (!is.list(cache) || !length(cache$tf_names)) return(NULL)
  motif_rows <- split(as.data.frame(cache$motif_counts, stringsAsFactors = FALSE), as.character(cache$motif_counts$motif))
  payload <- list(
    tfs = cache$tf_names,
    conditions = cache$conditions,
    defaultCondition = default_condition,
    tfCounts = unname(lapply(seq_len(nrow(cache$tf_counts)), function(i) as.integer(cache$tf_counts[i, ]))),
    motifs = motif_rows,
    motifSummary = cache$motif_summary %||% list(total = length(motif_rows), eligible = length(motif_rows), excluded = 0L)
  )
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", na = "null")
  payload_json <- gsub("</", "<\\/", payload_json, fixed = TRUE)
  binary_payload <- list(
    nBytes = if (length(cache$tf_bits)) length(cache$tf_bits[[1L]]) else 0L,
    tfBits = gsub("[\r\n]", "", jsonlite::base64_enc(do.call(c, cache$tf_bits))),
    boundBits = gsub("[\r\n]", "", jsonlite::base64_enc(do.call(c, cache$bound_bits))),
    geneOn = unname(lapply(seq_len(nrow(cache$gene_on)), function(i) as.integer(cache$gene_on[i, ])))
  )
  binary_payload_b64 <- .module2_report_browser_encode_browser_json_deflate_base64(binary_payload)
  css <- paste(
    ".m1-tabs{display:flex;gap:6px;margin:20px 0}.m1-tab{border:1px solid var(--line);background:var(--soft);padding:9px 14px;border-radius:5px;cursor:pointer;font-weight:700}.m1-tab.active{background:#e8f0f8;border-color:#8aa8c8;color:#173d69}",
    ".m1-panel.standalone{display:none}.m1-panel.standalone.active{display:block}.m1-controls{display:flex;flex-wrap:wrap;gap:8px;align-items:end;padding:8px 10px;background:var(--soft);border:1px solid var(--line);border-radius:4px;margin:6px 0 8px}.m1-control{display:grid;gap:2px}.m1-control label{font-size:11px;font-weight:700;color:#475569}.m1-controls select,.m1-controls input,.m1-controls button{font:inherit;font-size:12px;border:1px solid #b9c6d6;border-radius:4px;background:#fff;padding:5px 7px}.m1-controls input{min-width:220px}.m1-controls button{cursor:pointer}.m1-controls button.primary{background:#173d69;color:#fff;border-color:#173d69}.m1-controls button:disabled{opacity:.5;cursor:not-allowed}",
    ".m1-grid{display:grid;grid-template-columns:1fr 1fr;gap:10px}.m1-plot-card{border:1px solid var(--line);border-radius:4px;padding:8px 10px;overflow:hidden}.m1-plot-card h3{margin:0 0 4px;font-size:13px}.m1-chart{width:100%;min-height:190px;max-height:360px;overflow-y:auto;overflow-x:hidden;scrollbar-gutter:stable}.m1-chart svg{display:block;width:100%;height:auto}.m1-bar{fill:#315f97;cursor:pointer}.m1-bar.co{fill:#16847a}.m1-bar.rest{fill:#dbe5ef}.m1-bar.reference{fill:#fff;stroke:#b2182b;stroke-width:2}.m1-bar.highlight,.m1-dot.highlight{stroke:#b2182b;stroke-width:3}.m1-axis-label{font-size:11px;fill:#334155}.m1-value{font-size:11px;font-weight:700;fill:#172033}.m1-canonical{font-size:11px;font-weight:800;fill:#b2182b}.m1-reference-line{stroke:#9ca3af;stroke-width:1;stroke-dasharray:4 4}.m1-stem{stroke:#315f97;stroke-width:2}.m1-dot{fill:#6486af;cursor:pointer}.m1-status{font-size:11px;color:#526071;margin:5px 0}.m1-status.error{color:#a33b20}.m1-selected{display:flex;flex-wrap:wrap;gap:5px;min-height:24px;align-items:center}.m1-chip{display:inline-flex;align-items:center;gap:4px;border:1px solid #9fb4cc;background:#eef4fa;border-radius:999px;padding:3px 7px;font-size:11px;font-weight:700}.m1-chip button{border:0;background:transparent;padding:0;color:#7f1d1d;font-weight:800}",
    "@media(max-width:800px){.m1-grid{grid-template-columns:1fr}.m1-tabs{overflow:auto}.m1-chart{max-height:340px}.m1-control{min-width:0;max-width:100%}.m1-controls input{min-width:0;width:100%}}",
    sep = "\n"
  )
  worker_source <- paste0(
    "let P=null,tfBits=null,boundBits=null;const pop=Uint8Array.from({length:256},(_,i)=>{let n=i,c=0;while(n){c+=n&1;n>>=1}return c});",
    "const b64=s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0));async function decode(s){const bytes=b64(s),stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));return JSON.parse(await new Response(stream).text())}",
    "const count=a=>{let n=0;for(const x of a)n+=pop[x];return n};const seg=(a,i)=>a.subarray(i*P.nBytes,(i+1)*P.nBytes);",
    "onmessage=async e=>{const m=e.data;if(m.type==='init'){try{P=await decode(m.payload);tfBits=b64(P.tfBits);boundBits=b64(P.boundBits);postMessage({type:'ready'})}catch(err){postMessage({type:'error',message:String(err&&err.message||err)})}return}if(m.type!=='query'||!P)return;const ci=m.condition==='Overall'?-1:m.conditionIndex,n=P.nBytes,mask=new Uint8Array(n);mask.fill(255);for(const ti of m.selected){const src=seg(tfBits,ti);if(ci>=0&&!P.geneOn[ti][ci]){mask.fill(0);break}const bound=ci>=0?seg(boundBits,ci):null;for(let j=0;j<n;j++)mask[j]&=src[j]&(bound?bound[j]:255)}const total=new Array(m.nTfs),co=new Array(m.nTfs);for(let ti=0;ti<m.nTfs;ti++){const src=seg(tfBits,ti);if(ci>=0&&!P.geneOn[ti][ci]){total[ti]=0;co[ti]=0;continue}const bound=ci>=0?seg(boundBits,ci):null;let nt=0,nc=0;for(let j=0;j<n;j++){const v=src[j]&(bound?bound[j]:255);nt+=pop[v];nc+=pop[v&mask[j]]}total[ti]=nt;co[ti]=nc}postMessage({type:'result',token:m.token,total,co})}"
  )
  worker_json <- jsonlite::toJSON(worker_source, auto_unbox = TRUE)
  js <- paste0(
    "(()=>{const D=", payload_json, ";",
    "const tfIndex=new Map(D.tfs.map((x,i)=>[x,i])),condIndex=new Map(D.conditions.map((x,i)=>[x,i]));const esc=s=>String(s).replace(/[&<>\"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;'}[c]));let selectedTf='',selectedFocals=[],workerReady=false,queryToken=0,lastCobind=null;",
    "function bars(el,rows,stack=false,lollipop=false){const w=760,rh=22,left=170,right=105,h=Math.max(150,rows.length*rh+30),mx=Math.max(1,...rows.map(x=>x.total));let s='<svg viewBox=\"0 0 '+w+' '+h+'\">',referenceSeen=false;rows.forEach((d,i)=>{const y=8+i*rh,bw=(w-left-right)*d.total/mx,cw=(w-left-right)*(d.co||0)/mx,lc=d.canonical?'m1-canonical':'m1-axis-label';if(d.reference&&!referenceSeen){s+='<line x1=\"8\" y1=\"'+(y-4)+'\" x2=\"'+(w-8)+'\" y2=\"'+(y-4)+'\" class=\"m1-reference-line\"/>';referenceSeen=true}s+='<text x=\"'+(left-8)+'\" y=\"'+(y+12)+'\" text-anchor=\"end\" class=\"'+lc+'\">'+esc(d.tf)+'</text>';if(stack)s+='<rect x=\"'+left+'\" y=\"'+y+'\" width=\"'+bw+'\" height=\"14\" class=\"m1-bar rest\"/><rect x=\"'+left+'\" y=\"'+y+'\" width=\"'+cw+'\" height=\"14\" class=\"m1-bar co\" data-tf=\"'+esc(d.tf)+'\"><title>'+esc(d.tf)+': '+d.co+' co-bound of '+d.total+'</title></rect>';else if(lollipop)s+='<line x1=\"'+left+'\" y1=\"'+(y+7)+'\" x2=\"'+(left+bw)+'\" y2=\"'+(y+7)+'\" class=\"m1-stem\"/><circle cx=\"'+(left+bw)+'\" cy=\"'+(y+7)+'\" r=\"5\" class=\"m1-dot\" data-tf=\"'+esc(d.tf)+'\"><title>'+esc(d.tf)+': '+d.total+'</title></circle>';else s+='<rect x=\"'+left+'\" y=\"'+y+'\" width=\"'+Math.max(d.reference?1:bw,0)+'\" height=\"14\" class=\"m1-bar'+(d.reference?' reference':'')+'\" data-tf=\"'+esc(d.tf)+'\"><title>'+esc(d.tf)+': '+d.total+'</title></rect>';const value=d.valueLabel||d.total.toLocaleString();s+='<text x=\"'+Math.min(w-right/2,left+bw+8)+'\" y=\"'+(y+12)+'\" class=\"m1-value\">'+esc(value)+'</text>'});el.innerHTML=s+'</svg>';el.querySelectorAll('[data-tf]').forEach(x=>x.onclick=()=>highlight(x.dataset.tf))}",
    "function highlight(tf){selectedTf=tf;document.querySelectorAll('[data-tf]').forEach(x=>x.classList.toggle('highlight',x.dataset.tf===tf))}const n=id=>+document.getElementById(id).value;",
    "const counts=cond=>{const ci=cond==='Overall'?0:condIndex.get(cond)+1;return D.tfs.map((tf,i)=>({tf,total:+D.tfCounts[i][ci]||0})).sort((a,b)=>b.total-a.total||a.tf.localeCompare(b.tf))};",
    "function drawBinding(){bars(document.getElementById('binding-overall'),counts('Overall').slice(0,n('binding-n')),false,true);bars(document.getElementById('binding-condition'),counts(document.getElementById('binding-cond').value).slice(0,n('binding-n')),false,true);if(selectedTf)highlight(selectedTf)}",
    "function renderCobinding(result){const selected=new Set(result.names),rows=D.tfs.map((tf,i)=>({tf,total:+result.total[i]||0,co:+result.co[i]||0})).filter(x=>!selected.has(x.tf)).sort((a,b)=>b.co-a.co||b.total-a.total||a.tf.localeCompare(b.tf)).slice(0,n('cobind-n'));bars(document.getElementById('cobind-chart'),rows,true);document.getElementById('cobind-status').textContent=result.names.join(' + ')+' | '+result.condition+' | '+rows.length+' partners';lastCobind=result}",
    "function applyCobinding(){if(!selectedFocals.length){document.getElementById('cobind-status').textContent='Add at least one focal TF.';return}if(!workerReady){document.getElementById('cobind-status').textContent='Preparing exact co-binding data...';return}const condition=document.getElementById('cobind-cond').value,key=condition+'::'+selectedFocals.slice().sort().join('|'),cached=cobindCache.get(key);if(cached){renderCobinding(cached);return}const token=++queryToken;document.getElementById('cobind-status').textContent='Calculating exact co-binding...';worker.postMessage({type:'query',token,condition,conditionIndex:condIndex.get(condition),selected:selectedFocals.map(x=>tfIndex.get(x)),nTfs:D.tfs.length});pending.set(token,{key,names:selectedFocals.slice(),condition})}",
    "function drawMotif(){const motif=document.getElementById('motif-search').value,raw=D.motifs[motif]||[],statusEl=document.getElementById('motif-status'),chart=document.getElementById('motif-chart');statusEl.classList.remove('error');if(!raw.length){chart.replaceChildren();statusEl.textContent='No motif with retained canonical TF binding is available.';statusEl.classList.add('error');return}const rows=raw.map(x=>{const reference=!x.top20&&x.is_canonical,status=String(x.canonical_status||''),valueLabel=(+x.n||0).toLocaleString()+(status==='predicted_outside_top20'?' | rank '+x.rank:'');return{tf:x.tf,total:+x.n||0,canonical:!!x.is_canonical,reference,valueLabel}});const canonical=raw.filter(x=>x.is_canonical).map(x=>x.tf+': '+(x.canonical_status==='top_20'?'top 20, rank '+x.rank:'rank '+x.rank)),summary=D.motifSummary||{};statusEl.textContent='Canonical motif TF reference: '+canonical.join('; ')+' | '+(+summary.eligible||Object.keys(D.motifs).length)+' eligible motifs; '+(+summary.excluded||0)+' excluded without retained canonical TF binding.';bars(chart,rows)}",
    "window.exportSvg=function(id,name){const svg=document.querySelector('#'+id+' svg');if(!svg)return;const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)],{type:'image/svg+xml'}));a.download=name+'.svg';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),500)};",
    "function renderFocals(){const el=document.getElementById('focal-selected');el.innerHTML=selectedFocals.map(tf=>'<span class=\"m1-chip\">'+esc(tf)+'<button type=\"button\" data-remove=\"'+esc(tf)+'\" aria-label=\"Remove '+esc(tf)+'\">x</button></span>').join('');el.querySelectorAll('[data-remove]').forEach(x=>x.onclick=()=>{selectedFocals=selectedFocals.filter(tf=>tf!==x.dataset.remove);renderFocals()})}function addFocal(){const input=document.getElementById('focal-search'),tf=input.value,status=document.getElementById('cobind-status');status.classList.remove('error');if(!tfIndex.has(tf)){status.textContent='Choose a TF from the available list.';status.classList.add('error')}else if(selectedFocals.includes(tf)){status.textContent=tf+' is already selected.'}else{selectedFocals.push(tf);renderFocals();status.textContent='Added '+tf+'. Click Apply to recalculate.'}input.value=''}",
    "function fill(){for(const id of ['binding-cond','cobind-cond']){const s=document.getElementById(id);['Overall',...D.conditions].forEach(x=>s.add(new Option(x,x)));s.value=id==='binding-cond'?D.defaultCondition:'Overall'}const tfList=document.getElementById('tf-options');D.tfs.forEach(x=>tfList.appendChild(new Option('',x)));selectedFocals=[D.tfs[0]];renderFocals();const motifList=document.getElementById('motif-options'),motifs=Object.keys(D.motifs).sort();motifs.forEach(x=>motifList.appendChild(new Option('',x)));document.getElementById('motif-search').value=motifs[0]||''}",
    "function tab(id){document.querySelectorAll('.m1-tab').forEach(x=>x.classList.toggle('active',x.dataset.tab===id));document.querySelectorAll('.m1-panel.standalone').forEach(x=>x.classList.toggle('active',x.dataset.panel===id));history.replaceState(null,'','#'+id)}document.querySelectorAll('.m1-tab').forEach(x=>x.onclick=()=>tab(x.dataset.tab));",
    "const cobindCache=new Map(),pending=new Map(),worker=new Worker(URL.createObjectURL(new Blob([", worker_json, "],{type:'text/javascript'})));worker.onmessage=e=>{const m=e.data;if(m.type==='ready'){workerReady=true;document.getElementById('cobind-apply').disabled=false;document.getElementById('cobind-status').textContent='Exact co-binding data ready.';applyCobinding()}else if(m.type==='error'){const s=document.getElementById('cobind-status');s.textContent='Co-binding unavailable: '+m.message;s.classList.add('error')}else if(m.type==='result'){const q=pending.get(m.token);pending.delete(m.token);if(!q)return;const result={names:q.names,condition:q.condition,total:m.total,co:m.co};cobindCache.set(q.key,result);if(m.token===queryToken)renderCobinding(result)}};",
    "fill();document.getElementById('binding-cond').addEventListener('change',drawBinding);document.getElementById('binding-n').addEventListener('change',drawBinding);document.getElementById('motif-search').addEventListener('change',drawMotif);document.getElementById('focal-add').addEventListener('click',addFocal);document.getElementById('focal-search').addEventListener('keydown',e=>{if(e.key==='Enter'){e.preventDefault();addFocal()}});document.getElementById('cobind-apply').addEventListener('click',applyCobinding);document.getElementById('cobind-n').addEventListener('change',()=>{if(lastCobind)renderCobinding(lastCobind)});drawBinding();drawMotif();if(!('DecompressionStream' in window)){const s=document.getElementById('cobind-status');s.textContent='Co-binding requires a current Chrome, Edge, Firefox, or Safari browser.';s.classList.add('error')}else worker.postMessage({type:'init',payload:'", binary_payload_b64, "'});if(document.querySelector('.m1-tab'))tab(['binding','cobinding','motif'].includes(location.hash.slice(1))?location.hash.slice(1):'binding')})()"
  )
  panel_class <- if (isTRUE(standalone)) "m1-panel standalone" else "m1-panel"
  binding <- paste0(
    "<div class=\"m1-controls\"><div class=\"m1-control\"><label for=\"binding-cond\">Condition</label><select id=\"binding-cond\"></select></div><div class=\"m1-control\"><label for=\"binding-n\">Top N</label><select id=\"binding-n\"><option selected>10</option><option>20</option><option>50</option><option>100</option></select></div><button onclick=\"exportSvg('binding-overall','tfbs_binding_overall')\">Export overall SVG</button><button onclick=\"exportSvg('binding-condition','tfbs_binding_condition')\">Export condition SVG</button></div>",
    "<div class=\"m1-grid\"><div class=\"m1-plot-card\"><h3>Overall</h3><div id=\"binding-overall\" class=\"m1-chart\"></div></div><div class=\"m1-plot-card\"><h3>Selected condition</h3><div id=\"binding-condition\" class=\"m1-chart\"></div></div></div>"
  )
  cobinding <- paste0(
    "<p class=\"subtitle\">With multiple focal TFs, co-bound sites must be shared by every selected TF.</p><div class=\"m1-controls\"><div class=\"m1-control\"><label for=\"focal-search\">Add focal TF</label><input id=\"focal-search\" list=\"tf-options\" placeholder=\"Type a TF\"><datalist id=\"tf-options\"></datalist></div><button id=\"focal-add\" type=\"button\">Add TF</button><div class=\"m1-control\"><label for=\"cobind-cond\">Scope</label><select id=\"cobind-cond\"></select></div><div class=\"m1-control\"><label for=\"cobind-n\">Top N</label><select id=\"cobind-n\"><option selected>10</option><option>20</option><option>50</option><option>100</option></select></div><button id=\"cobind-apply\" class=\"primary\" type=\"button\" disabled>Apply</button><button onclick=\"exportSvg('cobind-chart','tfbs_cobinding')\">Export SVG</button></div><div id=\"focal-selected\" class=\"m1-selected\"></div><div id=\"cobind-status\" class=\"m1-status\">Preparing exact co-binding data...</div>",
    "<div class=\"m1-plot-card\"><div id=\"cobind-chart\" class=\"m1-chart\"></div></div>"
  )
  motif <- paste0(
    "<p class=\"subtitle\">The selector includes motifs with retained binding for every annotated canonical TF. The chart shows the true top 20 predicted TFs followed by canonical references outside the top 20.</p><div class=\"m1-controls\"><div class=\"m1-control\"><label for=\"motif-search\">Motif</label><input id=\"motif-search\" list=\"motif-options\" placeholder=\"Type a motif\"><datalist id=\"motif-options\"></datalist></div><button onclick=\"exportSvg('motif-chart','motif_predicted_tf')\">Export SVG</button></div><div id=\"motif-status\" class=\"m1-status\"></div>",
    "<div class=\"m1-plot-card\"><div id=\"motif-chart\" class=\"m1-chart\"></div></div>"
  )
  list(css = css, script = paste0("<script>", js, "</script>"), binding = binding, cobinding = cobinding, motif = motif, panel_class = panel_class)
}

#' Build the interactive Module 1 TFBS Explorer
#'
#' Builds an explicitly requested standalone HTML explorer. Standard
#' `predict_tfbs()` runs embed the same controls in the Module 1 QC report.
#' Exact co-binding queries run in a browser worker, and motif rankings show
#' the true top 20 plus canonical TF references at their actual rank.
#'
#' @param module1 A `predict_tfbs()` result, Module 1 output directory,
#'   predicted TFBS table, table path, or manifest path.
#' @param multiomic_data Optional CraftGRN multiomic object.
#' @param output_file HTML output path.
#' @param default_condition Initial condition. If `NULL`, use project
#'   configuration, control metadata, then A-Z order.
#' @param project_config Optional YAML path or configuration list.
#' @param rebuild_cache Rebuild the compact analysis cache even when valid.
#' @param verbose Emit concise progress messages.
#' @return Normalized path to the HTML report.
#' @export
build_module1_tfbs_explorer <- function(module1,
                                        multiomic_data = NULL,
                                        output_file = NULL,
                                        default_condition = NULL,
                                        project_config = NULL,
                                        rebuild_cache = FALSE,
                                        verbose = TRUE) {
  omics <- .module1_resolve_multiomic(module1, multiomic_data, project_config = project_config)
  if (is.null(output_file)) {
    output_file <- if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) {
      file.path(module1, "reports", "module1_tfbs_explorer.html")
    } else {
      file.path(getwd(), "module1_tfbs_explorer.html")
    }
  }
  if (!is.character(output_file) || length(output_file) != 1L || !nzchar(output_file)) {
    .log_abort("`output_file` must be a non-empty HTML path.")
  }
  cache_file <- .module1_explorer_cache_path(module1, output_file)
  cache <- .module1_build_explorer_cache(
    module1,
    omics,
    cache_file,
    project_config = project_config,
    rebuild = rebuild_cache,
    verbose = verbose
  )
  if (!length(cache$tf_names)) .log_abort("No predicted TFs are available for the TFBS Explorer.")
  default_condition <- .module1_explorer_default_condition(cache$conditions, default_condition, project_config, omics)
  x <- .module1_tfbs_explorer_components(cache, default_condition, standalone = TRUE)
  standalone_css <- paste0(
    ":root{--ink:#172033;--muted:#64748b;--line:#dbe3ed;--soft:#f6f8fb}*{box-sizing:border-box}body{margin:0;font:14px/1.45 Arial,Helvetica,sans-serif;color:var(--ink);background:#fff}.wrap{max-width:1240px;margin:auto;padding:24px 30px}header{border-bottom:1px solid var(--line)}h1{margin:0;font-size:28px}h2{font-size:20px}.subtitle{color:var(--muted)}section{border:0;padding:0}",
    x$css
  )
  html <- paste0(
    "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>Module 1 TFBS Explorer</title><style>", standalone_css, "</style></head><body>",
    "<header><div class=\"wrap\"><h1>Module 1 TFBS Explorer</h1><div class=\"subtitle\">Explore predicted binding, exact TF co-binding, and motif-centered TF rankings.</div></div></header><main class=\"wrap\">",
    "<nav class=\"m1-tabs\"><button class=\"m1-tab\" data-tab=\"binding\">Binding</button><button class=\"m1-tab\" data-tab=\"cobinding\">Co-binding</button><button class=\"m1-tab\" data-tab=\"motif\">Motif</button></nav>",
    "<section class=\"", x$panel_class, "\" data-panel=\"binding\"><h2>Predicted binding sites per TF</h2>", x$binding, "</section>",
    "<section class=\"", x$panel_class, "\" data-panel=\"cobinding\"><h2>TF co-binding summary</h2>", x$cobinding, "</section>",
    "<section class=\"", x$panel_class, "\" data-panel=\"motif\"><h2>Top predicted TFs at motif-containing footprints</h2>", x$motif, "</section>",
    "</main>", x$script, "</body></html>"
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  path <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  if (isTRUE(verbose)) .log_inform("Module 1 TFBS Explorer written: {path}")
  path
}
