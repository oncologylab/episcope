# Optional batch-correction evaluation used by project preparation workflows.

.batch_require_packages <- function() {
  required <- c("Biobase", "limma", "matrixStats", "RUVSeq")
  missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing)) {
    .log_abort("RUVr batch correction requires: {paste(missing, collapse = ', ')}.")
  }
  invisible(TRUE)
}

.batch_variable_rows <- function(x, top_n) {
  variance <- matrixStats::rowVars(x)
  good <- which(is.finite(variance) & variance > 0)
  if (!length(good)) .log_abort("Batch correction found no variable features.")
  good[order(variance[good], decreasing = TRUE)[seq_len(min(length(good), as.integer(top_n)))]]
}

.batch_condition_centroids <- function(x, condition) {
  levels <- unique(as.character(condition))
  out <- vapply(levels, function(level) {
    rowMeans(x[, condition == level, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(x)))
  rownames(out) <- rownames(x)
  colnames(out) <- levels
  out
}

.batch_within_condition_dispersion <- function(x, condition) {
  values <- vapply(unique(condition), function(level) {
    subset <- x[, condition == level, drop = FALSE]
    if (ncol(subset) < 2L) return(NA_real_)
    center <- rowMeans(subset)
    mean(sqrt(colMeans((subset - center)^2)))
  }, numeric(1L))
  mean(values, na.rm = TRUE)
}

.batch_effect_metrics <- function(raw, corrected, condition, batch, comparisons = NULL) {
  if (is.null(comparisons) || !nrow(comparisons)) {
    comparisons <- data.table::rbindlist(lapply(unique(batch), function(one_batch) {
      levels <- unique(condition[batch == one_batch])
      if (length(levels) < 2L) return(NULL)
      pairs <- utils::combn(levels, 2L)
      data.table::data.table(
        condition_1 = pairs[1L, ], condition_2 = pairs[2L, ], batch = one_batch
      )
    }), fill = TRUE)
  } else {
    comparisons <- data.table::as.data.table(comparisons)
  }
  if (!nrow(comparisons)) return(c(effect_spearman = 1, direction_concordance = 1))
  metrics <- lapply(seq_len(nrow(comparisons)), function(i) {
    row <- comparisons[i]
    c1 <- as.character(row$condition_1)
    c2 <- as.character(row$condition_2)
    raw_effect <- rowMeans(raw[, condition == c1, drop = FALSE]) -
      rowMeans(raw[, condition == c2, drop = FALSE])
    corrected_effect <- rowMeans(corrected[, condition == c1, drop = FALSE]) -
      rowMeans(corrected[, condition == c2, drop = FALSE])
    finite <- is.finite(raw_effect) & is.finite(corrected_effect)
    if (sum(finite) < 3L) return(c(effect_spearman = NA_real_, direction_concordance = NA_real_))
    c(
      effect_spearman = stats::cor(raw_effect[finite], corrected_effect[finite], method = "spearman"),
      direction_concordance = mean(sign(raw_effect[finite]) == sign(corrected_effect[finite]))
    )
  })
  values <- do.call(rbind, metrics)
  colMeans(values, na.rm = TRUE)
}

.fit_ruvr_candidate <- function(log_matrix, design, residuals, k) {
  if (k == 0L) return(list(corrected = log_matrix, w = matrix(numeric(), nrow = ncol(log_matrix), ncol = 0L)))
  fit <- RUVSeq::RUVr(
    x = log_matrix,
    cIdx = rownames(log_matrix),
    k = as.integer(k),
    residuals = residuals,
    center = TRUE,
    round = FALSE,
    isLog = TRUE
  )
  w <- as.matrix(fit$W)
  corrected <- limma::removeBatchEffect(log_matrix, covariates = w, design = design)
  corrected[!is.finite(corrected)] <- 0
  list(corrected = corrected, w = w)
}

.evaluate_ruvr_batch_correction <- function(matrix,
                                             metadata,
                                             sample_col = "sample_id",
                                             condition_col = "condition_id",
                                             batch_col = "study_id",
                                             k_candidates = c(0L, 1L, 2L, 3L, 5L, 10L),
                                             comparisons = NULL,
                                             top_n = 10000L,
                                             condition_distance_spearman_min = 0.95,
                                             effect_spearman_min = 0.95,
                                             direction_concordance_min = 0.90,
                                             between_condition_variance_min = 0.80,
                                             verbose = TRUE) {
  .batch_require_packages()
  metadata <- data.table::as.data.table(metadata)
  required <- c(sample_col, condition_col, batch_col)
  missing <- setdiff(required, names(metadata))
  if (length(missing)) .log_abort("Batch metadata is missing: {paste(missing, collapse = ', ')}.")
  samples <- intersect(as.character(metadata[[sample_col]]), colnames(matrix))
  if (length(samples) < 3L) .log_abort("RUVr evaluation requires at least three matched samples.")
  metadata <- metadata[match(samples, get(sample_col))]
  matrix <- as.matrix(matrix[, samples, drop = FALSE])
  storage.mode(matrix) <- "numeric"
  .resource_preflight(
    estimated_bytes = as.numeric(utils::object.size(matrix)) * 5,
    stage = paste0("RUVr evaluation for ", nrow(matrix), " features"),
    verbose = verbose
  )
  raw <- log1p(pmax(matrix, 0))
  use_rows <- .batch_variable_rows(raw, top_n)
  raw_use <- raw[use_rows, , drop = FALSE]
  condition <- as.character(metadata[[condition_col]])
  batch <- as.character(metadata[[batch_col]])
  design <- stats::model.matrix(~ 0 + factor(condition))
  residual_df <- ncol(raw_use) - qr(design)$rank
  max_k <- max(0L, residual_df - 1L)
  k_candidates <- sort(unique(as.integer(k_candidates)))
  k_candidates <- k_candidates[is.finite(k_candidates) & k_candidates >= 0L & k_candidates <= max_k]
  if (!0L %in% k_candidates) k_candidates <- c(0L, k_candidates)
  fit <- limma::lmFit(raw_use, design)
  residuals <- raw_use - fit$coefficients %*% t(design)
  raw_centroids <- .batch_condition_centroids(raw_use, condition)
  raw_distances <- as.numeric(stats::dist(t(raw_centroids)))
  raw_between <- mean(matrixStats::rowVars(raw_centroids), na.rm = TRUE)
  raw_dispersion <- .batch_within_condition_dispersion(raw_use, condition)

  candidates <- vector("list", length(k_candidates))
  audit <- vector("list", length(k_candidates))
  for (i in seq_along(k_candidates)) {
    k <- k_candidates[[i]]
    candidate <- .fit_ruvr_candidate(raw_use, design, residuals, k)
    centroids <- .batch_condition_centroids(candidate$corrected, condition)
    distances <- as.numeric(stats::dist(t(centroids)))
    distance_r <- if (length(raw_distances) >= 3L) {
      suppressWarnings(stats::cor(raw_distances, distances, method = "spearman"))
    } else {
      1
    }
    effect <- .batch_effect_metrics(raw_use, candidate$corrected, condition, batch, comparisons)
    between_ratio <- mean(matrixStats::rowVars(centroids), na.rm = TRUE) / raw_between
    dispersion <- .batch_within_condition_dispersion(candidate$corrected, condition)
    accepted <- k == 0L || (
      is.finite(distance_r) && distance_r >= condition_distance_spearman_min &&
      is.finite(effect[["effect_spearman"]]) && effect[["effect_spearman"]] >= effect_spearman_min &&
      is.finite(effect[["direction_concordance"]]) && effect[["direction_concordance"]] >= direction_concordance_min &&
      is.finite(between_ratio) && between_ratio >= between_condition_variance_min
    )
    candidates[[i]] <- candidate
    audit[[i]] <- data.table::data.table(
      k = k,
      accepted = accepted,
      within_condition_dispersion = dispersion,
      dispersion_ratio = dispersion / raw_dispersion,
      condition_distance_spearman = distance_r,
      effect_spearman = effect[["effect_spearman"]],
      direction_concordance = effect[["direction_concordance"]],
      between_condition_variance_ratio = between_ratio
    )
  }
  audit <- data.table::rbindlist(audit)
  accepted_positive <- audit[accepted == TRUE & k > 0L]
  selected_k <- 0L
  if (nrow(accepted_positive)) {
    best <- min(accepted_positive$within_condition_dispersion, na.rm = TRUE)
    selected_k <- min(accepted_positive[within_condition_dispersion <= best * 1.05, k])
  }
  audit[, selected := k == selected_k]
  selected <- candidates[[match(selected_k, k_candidates)]]
  full_corrected <- if (selected_k == 0L) {
    raw
  } else {
    corrected <- limma::removeBatchEffect(raw, covariates = selected$w, design = design)
    corrected[!is.finite(corrected)] <- 0
    corrected
  }
  if (isTRUE(verbose)) {
    .log_inform("RUVr evaluation selected K={selected_k} from: {paste(k_candidates, collapse = ', ')}.")
  }
  list(corrected_log = full_corrected, selected_k = selected_k, audit = audit, w = selected$w)
}
