# ==============================
# Multi-condition modeling (limma / lmer)
# ==============================

#' load_multi_links_matrix
#' @keywords internal
#' Build edge × sample matrix from multiple TF–gene links tables
#'
#' This is the multi-condition analogue of `compare_links_two_conditions()`.
#' It:
#'   - Loads each input via `load_links_table()`.
#'   - Optionally collapses duplicate keys per input.
#'   - Enforces identical key sets across all inputs when `strict = TRUE`.
#'   - Returns a list containing the key data.frame and a numeric matrix for
#'     a chosen value column (default `link_score`).
#'
#' @param inputs Character vector of CSV paths or list of data.frames/tibbles.
#' @param sample_names Optional character vector of sample names; must have
#'   length equal to `length(inputs)`. If `NULL`, names are derived from
#'   inputs via `.derive_cond_name()`.
#' @param value_col Name of the numeric column to extract (e.g. `"link_score"`,
#'   `"tf_expr"`, `"gene_expr"`).
#' @param strict Logical; if TRUE, require identical key sets across all inputs.
#' @param dedupe Logical; collapse duplicate keys within each input.
#' @param dedupe_method Passed to `.collapse_links_by_key()`.
#' @param verbose Logical.
#'
#' @return A list with elements:
#'   - `keys`: data.frame with columns `tf`, `gene_key`, `peak_id`.
#'   - `mat`: numeric matrix (edges × samples) for `value_col`.
#'   - `sample_names`: character vector of column names for `mat`.
#' @keywords internal
.load_multi_links_matrix <- function(inputs,
                                     sample_names = NULL,
                                     value_col    = "link_score",
                                     strict       = TRUE,
                                     dedupe       = TRUE,
                                     dedupe_method = c("max_link_score", "first"),
                                     verbose      = TRUE) {
  dedupe_method <- match.arg(dedupe_method)

  if (is.null(inputs)) {
    cli::cli_abort("`.load_multi_links_matrix()`: `inputs` cannot be NULL.")
  }

  # Normalize inputs to list
  if (is.data.frame(inputs)) {
    inputs <- list(inputs)
  } else if (is.character(inputs) || is.list(inputs)) {
    inputs <- as.list(inputs)
  } else {
    cli::cli_abort("`.load_multi_links_matrix()`: `inputs` must be a character vector of paths, a data.frame, or a list of data.frames/paths.")
  }

  n_in <- length(inputs)
  if (n_in < 2L) {
    cli::cli_abort("`.load_multi_links_matrix()`: need at least 2 inputs for multi-condition modeling.")
  }

  # Derive sample names if needed
  if (is.null(sample_names)) {
    sample_names <- vapply(
      seq_len(n_in),
      function(i) .derive_cond_name(inputs[[i]], paste0("cond", i)),
      character(1)
    )
  }
  if (length(sample_names) != n_in) {
    cli::cli_abort("`.load_multi_links_matrix()`: `sample_names` must have same length as `inputs`.")
  }

  key_cols <- c("tf", "gene_key", "peak_id")

  # Load tables
  tbls <- vector("list", n_in)
  for (i in seq_len(n_in)) {
    tbl <- load_links_table(inputs[[i]], verbose = verbose)
    if (isTRUE(dedupe)) {
      tbl <- .collapse_links_by_key(
        tbl,
        key_cols   = key_cols,
        which_cond = sample_names[[i]],
        method     = dedupe_method,
        verbose    = verbose
      )
    }
    tbls[[i]] <- tbl
  }

  # Enforce identical key sets across all inputs (pairwise vs first)
  if (isTRUE(strict)) {
    ref_tbl <- tbls[[1L]]
    for (i in 2:n_in) {
      .assert_identical_link_keys(
        ref_tbl, tbls[[i]],
        key_cols = key_cols,
        verbose  = verbose
      )
    }
  }

  # Deterministic key ordering from the first table
  ref_tbl  <- tbls[[1L]]
  ref_keys <- ref_tbl[key_cols]
  ord_ref  <- order(ref_keys$tf, ref_keys$gene_key, ref_keys$peak_id, method = "radix")
  ref_keys <- ref_keys[ord_ref, , drop = FALSE]

  # Allocate matrix
  n_edges <- nrow(ref_keys)
  mat     <- matrix(NA_real_, nrow = n_edges, ncol = n_in)
  colnames(mat) <- sample_names

  # Fill matrix by matching keys in each table to reference
  ref_tag <- interaction(ref_keys, drop = TRUE)

  for (i in seq_len(n_in)) {
    tbl_i <- tbls[[i]]
    if (!(value_col %in% names(tbl_i))) {
      cli::cli_abort("`.load_multi_links_matrix()`: value_col '{value_col}' not found in input {i}.")
    }
    keys_i <- tbl_i[key_cols]
    tag_i  <- interaction(keys_i, drop = TRUE)
    idx    <- match(ref_tag, tag_i)

    if (any(is.na(idx))) {
      cli::cli_abort(c(
        paste0("Key mismatch when aligning inputs for sample '", sample_names[[i]], "'."),
        i = "Some reference keys were not found in this input after previous checks."
      ))
    }

    vec <- as.numeric(tbl_i[[value_col]])
    mat[, i] <- vec[idx]
  }

  # Set rownames (optional but helpful)
  rownames(mat) <- paste(ref_keys$tf, ref_keys$gene_key, ref_keys$peak_id, sep = "|")

  list(
    keys         = ref_keys,
    mat          = mat,
    sample_names = sample_names
  )
}

#' Multi-condition ANOVA-like test on TF–gene link scores (limma F-test)
#'
#' This function performs a formal multi-group test per edge using limma,
#' analogous to an ANOVA F-test across all conditions at once (no pairwise
#' comparisons). Each input is one "sample" (e.g., a condition-averaged link
#' table, or one replicate).
#'
#' @param inputs Character vector of CSVs or list of data.frames/tibbles,
#'   one per sample/condition.
#' @param sample_info Data.frame with one row per input, in the same order.
#'   Must contain a column specified by `condition_col` (default `"condition"`),
#'   which will be treated as a factor.
#' @param condition_col Column in `sample_info` defining the groups (default
#'   `"condition"`).
#' @param value_col Column in the link tables to model (default `"link_score"`).
#' @param strict, dedupe, dedupe_method Passed to `.load_multi_links_matrix()`.
#' @param verbose Logical.
#'
#' @return A tibble with:
#'   - `tf`, `gene_key`, `peak_id`.
#'   - `F_stat`, `p_value`, `fdr`.
#'   - One column per design coefficient: `coef_<design_colname>`.
#'
#' @export
grn_diff_multi_anova <- function(inputs,
                                 sample_info,
                                 condition_col = "condition",
                                 value_col     = "link_score",
                                 strict        = TRUE,
                                 dedupe        = TRUE,
                                 dedupe_method = c("max_link_score", "first"),
                                 verbose       = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    cli::cli_abort("Package 'limma' is required for grn_diff_multi_anova(). Please install it.")
  }

  dedupe_method <- match.arg(dedupe_method)

  if (!is.data.frame(sample_info)) {
    cli::cli_abort("`sample_info` must be a data.frame.")
  }

  if (!(condition_col %in% names(sample_info))) {
    cli::cli_abort("`sample_info` must contain a column '{condition_col}'.")
  }

  # One row per input, in the same order
  n_in <- if (is.list(inputs) || is.data.frame(inputs)) length(as.list(inputs)) else length(inputs)
  if (nrow(sample_info) != n_in) {
    cli::cli_abort("`sample_info` must have one row per input (nrow(sample_info) == length(inputs)).")
  }

  # Build edge × sample matrix
  multi <- .load_multi_links_matrix(
    inputs        = inputs,
    sample_names  = NULL,
    value_col     = value_col,
    strict        = strict,
    dedupe        = dedupe,
    dedupe_method = dedupe_method,
    verbose       = verbose
  )
  keys <- multi$keys
  mat  <- multi$mat

  # Design: ~ 0 + condition
  sample_info[[condition_col]] <- as.factor(sample_info[[condition_col]])
  design <- stats::model.matrix(
    stats::as.formula(paste0("~ 0 + ", condition_col)),
    data = sample_info
  )

  .log(paste0("Fitting limma model with ", ncol(design), " design column(s) across ",
              nrow(mat), " edges."), verbose = verbose)

  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)

  # Global F-test for any condition effect
  F_stat <- fit$F
  p_val  <- fit$F.p.value
  fdr    <- stats::p.adjust(p_val, method = "BH")

  # Coefficients (effect sizes in design space)
  coef_mat <- fit$coefficients
  colnames(coef_mat) <- paste0("coef_", colnames(coef_mat))

  out <- cbind(
    keys,
    data.frame(
      F_stat = as.numeric(F_stat),
      p_value = as.numeric(p_val),
      fdr = as.numeric(fdr),
      coef_mat,
      check.names = FALSE
    )
  )

  tibble::as_tibble(out)
}

#' Multi-condition trend test for time/dose using limma
#'
#' Tests for a monotonic effect of a numeric covariate (time or dose) on
#' `value_col` per edge, using limma's moderated t-statistics. Returns effect
#' size (slope), SE, t, p, and FDR for the covariate.
#'
#' @param inputs Character vector of CSVs or list of data.frames/tibbles.
#' @param sample_info Data.frame with one row per input, in same order.
#' @param covariate Name of numeric column in `sample_info` (e.g. `"time"`,
#'   `"dose"`).
#' @param add_covariates Optional character vector of additional columns in
#'   `sample_info` to include as fixed effects (e.g. `"batch"`, `"cell_line"`).
#' @param value_col Column in the link tables to model (default `"link_score"`).
#' @param strict, dedupe, dedupe_method Passed to `.load_multi_links_matrix()`.
#' @param verbose Logical.
#'
#' @return A tibble with:
#'   - `tf`, `gene_key`, `peak_id`.
#'   - `slope`, `se`, `t`, `p_value`, `fdr` for the covariate.
#'   - Optionally other design coefficients as `coef_<name>`.
#'
#' @export
grn_diff_multi_trend <- function(inputs,
                                 sample_info,
                                 covariate,
                                 add_covariates = NULL,
                                 value_col      = "link_score",
                                 strict         = TRUE,
                                 dedupe         = TRUE,
                                 dedupe_method  = c("max_link_score", "first"),
                                 verbose        = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    cli::cli_abort("Package 'limma' is required for grn_diff_multi_trend(). Please install it.")
  }

  dedupe_method <- match.arg(dedupe_method)

  if (!is.data.frame(sample_info)) {
    cli::cli_abort("`sample_info` must be a data.frame.")
  }
  if (!(covariate %in% names(sample_info))) {
    cli::cli_abort("`sample_info` must contain covariate column '{covariate}'.")
  }

  # Ensure numeric covariate
  x <- sample_info[[covariate]]
  if (!is.numeric(x)) {
    x_num <- suppressWarnings(as.numeric(x))
    if (any(is.na(x_num) & !is.na(x))) {
      cli::cli_abort("Covariate '{covariate}' is not numeric and cannot be safely coerced.")
    }
    sample_info[[covariate]] <- x_num
  }

  # One row per input
  n_in <- if (is.list(inputs) || is.data.frame(inputs)) length(as.list(inputs)) else length(inputs)
  if (nrow(sample_info) != n_in) {
    cli::cli_abort("`sample_info` must have one row per input (nrow(sample_info) == length(inputs)).")
  }

  multi <- .load_multi_links_matrix(
    inputs        = inputs,
    sample_names  = NULL,
    value_col     = value_col,
    strict        = strict,
    dedupe        = dedupe,
    dedupe_method = dedupe_method,
    verbose       = verbose
  )
  keys <- multi$keys
  mat  <- multi$mat

  # Build design: ~ 1 + covariate + add_covariates
  rhs <- c(covariate, add_covariates)
  rhs <- rhs[!is.na(rhs) & nzchar(rhs)]

  form_str <- if (length(rhs)) {
    paste("~ 1 +", paste(rhs, collapse = " + "))
  } else {
    paste("~ 1 +", covariate)
  }
  design <- stats::model.matrix(stats::as.formula(form_str), data = sample_info)

  .log(paste0("Fitting limma trend model with design: ", form_str,
              " across ", nrow(mat), " edges."), verbose = verbose)

  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)

  # Locate covariate column in design
  cov_idx <- match(covariate, colnames(design))
  if (is.na(cov_idx)) {
    cli::cli_abort("Covariate '{covariate}' not found as a column in the design matrix.")
  }

  slope <- fit$coefficients[, cov_idx]
  se    <- sqrt(fit$s2.post) * fit$stdev.unscaled[, cov_idx]
  tstat <- fit$t[, cov_idx]
  p_val <- fit$p.value[, cov_idx]
  fdr   <- stats::p.adjust(p_val, method = "BH")

  # Other coefficients (optional)
  coef_mat <- fit$coefficients
  colnames(coef_mat) <- paste0("coef_", colnames(coef_mat))

  out <- cbind(
    keys,
    data.frame(
      slope   = as.numeric(slope),
      se      = as.numeric(se),
      t       = as.numeric(tstat),
      p_value = as.numeric(p_val),
      fdr     = as.numeric(fdr),
      coef_mat,
      check.names = FALSE
    )
  )

  tibble::as_tibble(out)
}

#' Edge-wise mixed-effects modeling with lmer on TF–gene links
#'
#' This function fits a linear mixed-effects model per edge using `lmer`,
#' allowing you to:
#'   - Model time/dose and other covariates as fixed effects.
#'   - Include per-sample/subject random effects (e.g., `(1 | mouse_id)`).
#'   - Extract one or more "biological" fixed effects (e.g. `treatmentLPS`
#'     for LPS vs Ctrl) with formal tests, while adjusting for time/dose.
#'
#' **NOTE:** This is computationally heavier than the limma-based functions and
#' is best applied to a subset of edges (e.g., pre-filtered by limma).
#'
#' @param inputs Character vector of CSVs or list of data.frames/tibbles.
#' @param sample_info Data.frame with one row per input, in the same order.
#'   All variables referenced in `formula` must exist in `sample_info`.
#' @param formula Model formula for `lmer` **using response name `value`**,
#'   e.g., `value ~ treatment + time + (1 | mouse_id)`.
#' @param effect_terms Character vector of row names in the fixed-effects
#'   coefficient table to extract (e.g., `"treatmentLPS"`, `"time"`).
#' @param value_col Column in the link tables to model (default `"link_score"`).
#' @param strict, dedupe, dedupe_method Passed to `.load_multi_links_matrix()`.
#' @param verbose Logical.
#' @param ... Additional arguments passed to `lme4::lmer()` / `lmerTest::lmer()`
#'   (e.g., `REML = FALSE`).
#'
#' @return A tibble with:
#'   - `tf`, `gene_key`, `peak_id`.
#'   - For each `effect_terms[i]`, columns:
#'       - `beta_<term>`, `se_<term>`, `t_<term>`, `p_<term>`, `fdr_<term>`.
#'
#' @export
grn_diff_multi_lmm <- function(inputs,
                               sample_info,
                               formula,
                               effect_terms,
                               value_col     = "link_score",
                               strict        = TRUE,
                               dedupe        = TRUE,
                               dedupe_method = c("max_link_score", "first"),
                               verbose       = TRUE,
                               ...) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    cli::cli_abort("Package 'lme4' is required for grn_diff_multi_lmm(). Please install it.")
  }

  use_lmerTest <- FALSE
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    use_lmerTest <- TRUE
  }

  dedupe_method <- match.arg(dedupe_method)

  if (!is.data.frame(sample_info)) {
    cli::cli_abort("`sample_info` must be a data.frame.")
  }
  if (!inherits(formula, "formula")) {
    cli::cli_abort("`formula` must be a valid model formula with response 'value'.")
  }
  if (length(effect_terms) == 0L) {
    cli::cli_abort("`effect_terms` must be a non-empty character vector.")
  }

  # One row per input
  n_in <- if (is.list(inputs) || is.data.frame(inputs)) length(as.list(inputs)) else length(inputs)
  if (nrow(sample_info) != n_in) {
    cli::cli_abort("`sample_info` must have one row per input (nrow(sample_info) == length(inputs)).")
  }

  multi <- .load_multi_links_matrix(
    inputs        = inputs,
    sample_names  = NULL,
    value_col     = value_col,
    strict        = strict,
    dedupe        = dedupe,
    dedupe_method = dedupe_method,
    verbose       = verbose
  )
  keys <- multi$keys
  mat  <- multi$mat

  n_edges <- nrow(mat)
  n_terms <- length(effect_terms)

  # Allocate result containers
  beta_mat <- matrix(NA_real_, nrow = n_edges, ncol = n_terms)
  se_mat   <- matrix(NA_real_, nrow = n_edges, ncol = n_terms)
  t_mat    <- matrix(NA_real_, nrow = n_edges, ncol = n_terms)
  p_mat    <- matrix(NA_real_, nrow = n_edges, ncol = n_terms)

  colnames(beta_mat) <- effect_terms
  colnames(se_mat)   <- effect_terms
  colnames(t_mat)    <- effect_terms
  colnames(p_mat)    <- effect_terms

  .log(paste0("Fitting lmer model for ", n_edges,
              " edges; this may take time. Using lmerTest: ", use_lmerTest), verbose = verbose)

  for (i in seq_len(n_edges)) {
    y <- as.numeric(mat[i, ])

    # Skip all-NA edges
    if (all(is.na(y))) next

    dat <- sample_info
    dat[["value"]] <- y

    # Drop rows with NA in response or covariates (lmer will fail otherwise)
    complete <- stats::complete.cases(dat)
    if (!all(complete)) {
      dat <- dat[complete, , drop = FALSE]
    }

    if (nrow(dat) < 3L) next  # too few samples to fit a mixed model

    fit <- tryCatch(
      {
        if (use_lmerTest) {
          lmerTest::lmer(formula, data = dat, ...)
        } else {
          lme4::lmer(formula, data = dat, ...)
        }
      },
      error = function(e) NULL
    )

    if (is.null(fit)) next

    # Extract fixed-effects table
    if (use_lmerTest) {
      summ  <- summary(fit)
      coef_tab <- as.data.frame(summ$coefficients)
    } else {
      coef_tab <- as.data.frame(stats::coef(summary(fit)))
    }

    # Rows are terms; columns: Estimate, Std. Error, t value, (Pr(>|t|)) if lmerTest
    for (j in seq_len(n_terms)) {
      term <- effect_terms[[j]]
      if (!(term %in% rownames(coef_tab))) next

      beta_mat[i, j] <- coef_tab[term, "Estimate"]
      se_mat[i, j]   <- coef_tab[term, "Std. Error"]
      t_mat[i, j]    <- coef_tab[term, "t value"]

      if (use_lmerTest && ("Pr(>|t|)" %in% colnames(coef_tab))) {
        p_mat[i, j] <- coef_tab[term, "Pr(>|t|)"]
      } else {
        p_mat[i, j] <- NA_real_
      }
    }
  }

  # FDR correction per term
  fdr_mat <- matrix(NA_real_, nrow = n_edges, ncol = n_terms)
  colnames(fdr_mat) <- effect_terms
  if (use_lmerTest) {
    for (j in seq_len(n_terms)) {
      p_j <- p_mat[, j]
      if (all(is.na(p_j))) next
      fdr_mat[, j] <- stats::p.adjust(p_j, method = "BH")
    }
  }

  # Build output data.frame
  out <- keys

  for (j in seq_len(n_terms)) {
    term <- effect_terms[[j]]
    out[[paste0("beta_", term)]] <- beta_mat[, j]
    out[[paste0("se_",   term)]] <- se_mat[, j]
    out[[paste0("t_",    term)]] <- t_mat[, j]
    out[[paste0("p_",    term)]] <- p_mat[, j]
    out[[paste0("fdr_",  term)]] <- fdr_mat[, j]
  }

  tibble::as_tibble(out)
}

# base_dir <- "/data/homes/yl814/episcope_test/GSE87218_ATAC"
# match_mode <- "strict"  # strict lenient
# regulated_genes <- 1.5 # 1.5 2
# delta_link <- 1 # 2 1
# regulation_priors <- "300kb" # "genehancer"
# threshold_fp_tf_corr_p <- 0.3
# lighting_folder <- file.path(base_dir, paste("lighting", "fp_tf_corr_FDR", threshold_fp_tf_corr_p, regulation_priors, db, "regulated_genes", regulated_genes, "delta_link", delta_link, sep = "_"))
# print(lighting_folder)

# 1. Multi-condition ANOVA on all edges -----------------------------------
# Suppose  have 5 conditions: Ctrl, LPS_1h, LPS_4h, LPS_24h, LPS_24h_culture_5d
# For multi-condition ANOVA-like tests, you must have at least some replicates

# inputs <- c(
#   file.path(lighting_folder, "lighting_cond-Ctrl_tf_gene_links.csv"),
#   file.path(lighting_folder, "lighting_cond-LPS_1h_post_tf_gene_links.csv"),
#   file.path(lighting_folder, "lighting_cond-LPS_4h_post_tf_gene_links.csv"),
#   file.path(lighting_folder, "lighting_cond-LPS_24h_post_tf_gene_links.csv"),
#   file.path(lighting_folder, "lighting_cond-LPS_24h_culture_5d_tf_gene_links.csv")
# )

# sample_info <- tibble::tibble(
#   condition = c("Ctrl", "LPS_1h", "LPS_4h", "LPS_24h", "LPS_24h_cult5d")
# )

# res_anova <- grn_diff_multi_anova(
#   inputs       = inputs,
#   sample_info  = sample_info,
#   condition_col = "condition",
#   value_col    = "link_score"
# )

# 2. Trend analysis for time or dose --------------------------------------
# sample_info <- tibble::tibble(
#   condition = c("Ctrl", "LPS_1h", "LPS_4h", "LPS_24h", "LPS_24h_cult5d"),
#   time_h    = c(0, 1, 4, 24, 24 + 24*4)  # example encoding
# )

# res_trend <- grn_diff_multi_trend(
#   inputs        = inputs,
#   sample_info   = sample_info,
#   covariate     = "time_h",
#   add_covariates = NULL,
#   value_col     = "link_score"
# )
# res_trend |> dplyr::filter(p_value < 0.05)

# 3. LMM with biological effect adjusted for time/dose --------------------
# build_lighting_sample_info <- function(lighting_folder, verbose = TRUE) {
#   if (!dir.exists(lighting_folder)) {
#     cli::cli_abort("build_lighting_sample_info(): `lighting_folder` does not exist: {lighting_folder}")
#   }
#
#   # Only per-condition files; this automatically excludes lighting_overall_*.csv
#   pat <- "^lighting_cond-.*_tf_gene_links\\.csv$"
#   files <- list.files(lighting_folder, pattern = pat, full.names = TRUE)
#   if (!length(files)) {
#     cli::cli_abort("build_lighting_sample_info(): no lighting_cond-*_tf_gene_links.csv files found in `{lighting_folder}`.")
#   }
#
#   parse_one <- function(f) {
#     b <- basename(f)
#     # Strip prefix/suffix:
#     #  "lighting_cond-LPS_1h_post_tf_gene_links.csv" -> "LPS_1h_post"
#     stub <- sub("^lighting_cond-", "", b)
#     stub <- sub("_tf_gene_links\\.csv$", "", stub)
#
#     parts <- strsplit(stub, "_", fixed = TRUE)[[1]]
#     stimulus <- parts[1]
#     rest     <- parts[-1]
#
#     # Default: Ctrl-like samples with no explicit time -> 0h
#     time_h <- 0
#     time_label <- "0h"
#
#     if (length(rest) > 0L) {
#       # Accumulate time from tokens like "1h", "4h", "24h", "6d", "5d", etc.
#       time_h_num <- 0
#       for (tok in rest) {
#         if (grepl("^[0-9]+h$", tok)) {
#           time_h_num <- time_h_num + as.numeric(sub("h$", "", tok))
#         } else if (grepl("^[0-9]+d$", tok)) {
#           time_h_num <- time_h_num + 24 * as.numeric(sub("d$", "", tok))
#         }
#       }
#       # If no numeric tokens at all, keep 0; otherwise use accumulated hours
#       if (!is.na(time_h_num) && time_h_num > 0) {
#         time_h <- time_h_num
#       }
#       time_label <- paste(rest, collapse = "_")
#     }
#
#     data.frame(
#       file        = f,
#       sample_name = stub,
#       stimulus    = stimulus,
#       time_h      = as.numeric(time_h),
#       time_label  = time_label,
#       stringsAsFactors = FALSE
#     )
#   }
#
#   out <- do.call(rbind, lapply(files, parse_one))
#
#   # Factor-ize stimulus and set Ctrl as reference if present
#   out$stimulus <- factor(out$stimulus)
#   if (any(out$stimulus == "Ctrl")) {
#     out$stimulus <- stats::relevel(out$stimulus, ref = "Ctrl")
#   }
#
#   if (isTRUE(verbose)) {
#     .log(paste0(
#       "Detected ", nrow(out), " samples across stimuli: ",
#       paste(levels(out$stimulus), collapse = ", ")
#     ), verbose = TRUE)
#   }
#
#   tibble::as_tibble(out)
# }
# sample_info <- build_lighting_sample_info(lighting_folder)

# res_lmm <- grn_diff_multi_lmm(
#   inputs        = sample_info$file,
#   sample_info   = sample_info,
#   formula       = value ~ time_h + (1 | stimulus),
#   effect_terms  = c("time_h"),
#   value_col     = "link_score"
# )

# res_lmm |> dplyr::filter(fdr_time_h < 0.05)




