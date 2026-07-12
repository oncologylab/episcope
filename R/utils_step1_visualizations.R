# File: utils_step1_visualizations.R
# Purpose: Public Module 1 TFBS visualization helpers.

.module1_relevant_config <- function(config) {
  if (is.null(config)) return(list())
  if (is.character(config) && length(config) == 1L && file.exists(config)) {
    config <- yaml::read_yaml(config)
  }
  if (!is.list(config)) return(list())
  keep <- c(
    "project", "project_name", "project_date", "base_dir", "db", "ref_genome",
    "label_col", "fp_root_dir", "fp_cache_dir", "fp_cache_tag", "atac_master",
    "atac_data_path", "rna_mapped", "rna_path", "sample_metadata", "metadata_path",
    "threshold_expr", "threshold_gene_expr", "threshold_fp_score",
    "threshold_fp_tf_corr_r", "threshold_fp_tf_corr_p", "threshold_fp_tf_corr_fdr",
    "filter_to_canonical_bound", "footprint_sample_scope", "mid_slop",
    "round_digits", "score_match_pct", "output_mode"
  )
  out <- config[intersect(keep, names(config))]
  secret <- grepl("token|password|secret|credential|api[_-]?key", names(out), ignore.case = TRUE)
  out[secret] <- "[REDACTED]"
  out
}

.module1_tfbs_source <- function(module1) {
  if (is.data.frame(module1)) return(list(data = tibble::as_tibble(module1)))
  if (is.list(module1)) {
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

.module1_resolve_multiomic <- function(module1, multiomic_data = NULL) {
  if (is_multiomic_object(multiomic_data)) return(multiomic_data)
  if (is.list(module1) && is_multiomic_object(module1$omics_data)) return(module1$omics_data)
  if (is.character(module1) && length(module1) == 1L && dir.exists(module1)) {
    rds <- .qc_find_module1_omics_rds(module1)
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
    partials[[length(partials) + 1L]] <<- dt[, .(
      n_tfbs_cond1 = sum(active1, na.rm = TRUE),
      n_tfbs_cond2 = sum(active2, na.rm = TRUE),
      fp_score_cond1 = sum(data.table::fifelse(active1, fp1, 0), na.rm = TRUE),
      fp_score_cond2 = sum(data.table::fifelse(active2, fp2, 0), na.rm = TRUE),
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
