# Standalone FP-only GAM cache runner (de-duplicated by fp_peak × gene_key)
# - No reliance on global variables: all required inputs are function parameters
# - Explicit namespace usage (no library())
# - Includes all helper functions needed (including model_summary + GAM fit)
#
# Required packages (installed):
#   cli, dplyr, tidyr, tibble, readr, mgcv
# Optional:
#   parallel (for mclapply on Unix), broom, broom.mixed

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.left_join_mm <- function(x, y, by) {
  # dplyr >= 1.1.0 warns on expected many-to-many joins; silence where intended.
  tryCatch(
    dplyr::left_join(x, y, by = by, relationship = "many-to-many"),
    error = function(e) dplyr::left_join(x, y, by = by)
  )
}

# ---- model_summary() ---------------------------------------------------------
model_summary <- function(model) {
  out <- list()
  out$model_class <- class(model)[1]

  # Coefficients
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")) {
    if (requireNamespace("broom.mixed", quietly = TRUE)) {
      out$coef <- broom.mixed::tidy(model, effects = "fixed")
    } else {
      out$coef <- tibble::as_tibble(summary(model)$coefficients, rownames = "term")
    }
  } else if (inherits(model, "gam")) {
    if (requireNamespace("broom", quietly = TRUE)) {
      out$coef <- broom::tidy(model)
    } else {
      out$coef <- tibble::as_tibble(summary(model)$p.table, rownames = "term")
    }
  } else if (inherits(model, "lm")) {
    if (requireNamespace("broom", quietly = TRUE)) {
      out$coef <- broom::tidy(model)
    } else {
      out$coef <- tibble::as_tibble(summary(model)$coefficients, rownames = "term")
    }
  } else {
    out$coef <- NULL
  }

  # Model-level statistics
  if (requireNamespace("broom", quietly = TRUE)) {
    out$glance <- tryCatch(broom::glance(model), error = function(e) NULL)
  } else {
    out$glance <- NULL
  }

  # Fitted and residuals (best-effort)
  out$fitted <- tryCatch(stats::fitted(model), error = function(e) NULL)
  out$resid  <- tryCatch(stats::residuals(model), error = function(e) NULL)

  out
}

# ---- tfs_to_target_genes_fp_tobias_gam() ------------------------------------
tfs_to_target_genes_fp_tobias_gam <- function(df, k = 5L) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_fp, k = k) + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- GAM grid (method/k/bs) + gam.check ------------------------------------
.with_null_pdf <- function(expr) {
  if (!requireNamespace("grDevices", quietly = TRUE)) return(force(expr))
  fnull <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  grDevices::pdf(file = fnull)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

.gam_check_metrics <- function(fit, keep_output = FALSE) {
  tmp <- NULL
  out_txt <- utils::capture.output(
    {
      tmp <<- tryCatch(
        .with_null_pdf(mgcv::gam.check(fit, rep = 0)),
        error = function(e) e
      )
    },
    type = "output"
  )

  if (inherits(tmp, "error")) {
    return(list(
      ok = FALSE,
      error = conditionMessage(tmp),
      k_prime = NA_real_,
      edf = NA_real_,
      k_index = NA_real_,
      k_p_value = NA_real_,
      output = if (isTRUE(keep_output)) paste(out_txt, collapse = "\n") else NA_character_
    ))
  }

  # Prefer returned k.check; fall back to mgcv::k.check() if needed.
  kchk <- NULL
  if (is.list(tmp) && !is.null(tmp$k.check)) kchk <- tmp$k.check
  if (is.null(kchk)) {
    kchk2 <- tryCatch(mgcv::k.check(fit), error = function(e) NULL)
    if (is.list(kchk2) && !is.null(kchk2$k.check)) kchk <- kchk2$k.check else kchk <- kchk2
  }

  if (inherits(kchk, "data.frame")) kchk <- as.matrix(kchk)

  k_index <- NA_real_
  k_p <- NA_real_
  k_prime <- NA_real_
  edf <- NA_real_
  if (!is.null(kchk) && is.matrix(kchk) && nrow(kchk) >= 1) {
    cn <- colnames(kchk) %||% character(0)
    # mgcv uses "k'" as a column name (apostrophe); be tolerant.
    if ("k'" %in% cn) k_prime <- as.numeric(kchk[1, "k'"])
    if ("edf" %in% cn) edf <- as.numeric(kchk[1, "edf"])
    if ("k-index" %in% cn) k_index <- as.numeric(kchk[1, "k-index"])
    if ("p-value" %in% cn) k_p <- as.numeric(kchk[1, "p-value"])
  }

  list(
    ok = TRUE,
    error = NA_character_,
    k_prime = k_prime,
    edf = edf,
    k_index = k_index,
    k_p_value = k_p,
    output = if (isTRUE(keep_output)) paste(out_txt, collapse = "\n") else NA_character_
  )
}

fit_fp_gam_one <- function(df,
                           k = 5L,
                           bs = "tp",
                           method = "REML") {
  s <- mgcv::s
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) k <- 5L

  tf_fp_num <- df$tf_fp
  tf_fp_num <- tf_fp_num[is.finite(tf_fp_num)]
  n_unique_tf_fp <- length(unique(tf_fp_num))
  if (!is.finite(n_unique_tf_fp) || n_unique_tf_fp < 3L) {
    stop("tf_fp has < 3 unique finite values; GAM smooth not identifiable.")
  }

  k_use <- if (k < 3L) 3L else k
  if (k_use > n_unique_tf_fp) k_use <- n_unique_tf_fp
  if (k_use < 3L) {
    stop("k after capping is < 3; cannot fit smooth with this data.")
  }

  df$cell   <- factor(df$cell)
  df$stress <- factor(df$stress)

  rhs <- paste0("s(tf_fp, k = ", k_use, ", bs = \"", bs, "\")")
  if (nlevels(df$cell) > 1L) rhs <- paste(rhs, "cell", sep = " + ")
  if (nlevels(df$stress) > 1L) rhs <- paste(rhs, "stress", sep = " + ")
  form <- stats::as.formula(paste("expr ~", rhs))
  environment(form) <- environment()

  fit <- mgcv::gam(
    formula = form,
    data = df,
    method = method
  )

  list(
    fit = fit,
    k_requested = k,
    k_used = k_use
  )
}

run_fp_gam_grid_for_pairs <- function(fp_score,
                                      rna_tbl,
                                      sample_md,
                                      pair_tbl,
                                      max_pairs = 50L,
                                      order_by = c("p_fp", "abs_r_fp", "none"),
                                      methods = c("REML", "ML", "GCV.Cp"),
                                      k_values = 1:5,
                                      bs_values = c("tp", "cr", "ps"),
                                      min_samples = 6L,
                                      standardize = c("all", "predictors", "none"),
                                      run_gam_check = TRUE,
                                      keep_models = FALSE,
                                      n_workers = NULL,
                                      prefer_fork = TRUE,
                                      keep_gamcheck_output = FALSE) {
  standardize <- match.arg(standardize)
  order_by <- match.arg(order_by)

  # ---- package checks ----
  need_pkgs <- c("cli", "dplyr", "tidyr", "tibble", "mgcv")
  missing <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    cli::cli_abort(c(
      "Missing required packages.",
      "x" = "Install: {.val {missing}}"
    ))
  }

  # ---- input checks ----
  if (!is.data.frame(fp_score) || !"peak_ID" %in% names(fp_score)) {
    cli::cli_abort("fp_score must be a data.frame with column: peak_ID")
  }
  if (!is.data.frame(rna_tbl) || !"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("rna_tbl must be a data.frame with column: HGNC")
  }
  if (!is.data.frame(sample_md) || !all(c("id", "cell", "stress_type") %in% names(sample_md))) {
    cli::cli_abort("sample_md must have columns: id, cell, stress_type")
  }
  if (!is.data.frame(pair_tbl) || !all(c("fp_peak", "gene_key") %in% names(pair_tbl))) {
    cli::cli_abort("pair_tbl must have columns: fp_peak, gene_key")
  }

  sample_ids <- as.character(sample_md$id)
  if (!length(sample_ids)) cli::cli_abort("sample_md$id is empty; no samples to model.")
  miss_fp  <- setdiff(sample_ids, names(fp_score))
  miss_rna <- setdiff(sample_ids, names(rna_tbl))
  if (length(miss_fp)) cli::cli_abort(c("Some sample ids are missing from fp_score columns.", "x" = "Missing: {.val {miss_fp}}"))
  if (length(miss_rna)) cli::cli_abort(c("Some sample ids are missing from rna_tbl columns.", "x" = "Missing: {.val {miss_rna}}"))

  pair_tbl <- tibble::as_tibble(pair_tbl) |>
    dplyr::filter(!is.na(.data$fp_peak), .data$fp_peak != "",
                  !is.na(.data$gene_key), .data$gene_key != "") |>
    dplyr::distinct(.data$fp_peak, .data$gene_key, .keep_all = TRUE)

  if (!nrow(pair_tbl)) cli::cli_abort("pair_tbl has 0 valid pairs after filtering.")

  if (order_by == "p_fp" && "p_fp" %in% names(pair_tbl)) {
    pair_tbl <- dplyr::arrange(pair_tbl, .data$p_fp)
  } else if (order_by == "abs_r_fp" && "r_fp" %in% names(pair_tbl)) {
    pair_tbl <- dplyr::arrange(pair_tbl, dplyr::desc(abs(.data$r_fp)))
  }

  max_pairs <- as.integer(max_pairs)
  if (is.finite(max_pairs) && max_pairs > 0L && nrow(pair_tbl) > max_pairs) {
    pair_tbl <- dplyr::slice_head(pair_tbl, n = max_pairs)
  }

  pair_tbl$pair_key <- paste(pair_tbl$fp_peak, pair_tbl$gene_key, sep = "|")

  # ---- long table for these pairs ----
  peaks_use <- unique(pair_tbl$fp_peak)
  genes_use <- unique(pair_tbl$gene_key)

  fp_sub <- fp_score[fp_score$peak_ID %in% peaks_use, c("peak_ID", sample_ids), drop = FALSE]
  if (!nrow(fp_sub)) cli::cli_abort("No fp_score rows matched pair_tbl$fp_peak.")
  rna_sub <- rna_tbl[rna_tbl$HGNC %in% genes_use, c("HGNC", sample_ids), drop = FALSE]
  if (!nrow(rna_sub)) cli::cli_abort("No rna_tbl rows matched pair_tbl$gene_key.")

  fp_long <- tidyr::pivot_longer(
    data = fp_sub,
    cols = dplyr::all_of(sample_ids),
    names_to = "sample_id",
    values_to = "tf_fp"
  )
  names(fp_long)[names(fp_long) == "peak_ID"] <- "fp_peak"

  fp_long <- .left_join_mm(
    fp_long,
    pair_tbl[, c("fp_peak", "gene_key"), drop = FALSE],
    by = "fp_peak"
  )
  fp_long <- fp_long[!is.na(fp_long$gene_key) & fp_long$gene_key != "", , drop = FALSE]

  rna_long <- tidyr::pivot_longer(
    data = rna_sub,
    cols = dplyr::all_of(sample_ids),
    names_to = "sample_id",
    values_to = "expr"
  )
  names(rna_long)[names(rna_long) == "HGNC"] <- "gene_key"

  long_tbl <- dplyr::inner_join(fp_long, rna_long, by = c("gene_key", "sample_id"))
  if (!nrow(long_tbl)) cli::cli_abort("No rows after joining fp_long and rna_long.")

  md_use <- sample_md[, c("id", "cell", "stress_type"), drop = FALSE]
  names(md_use)[names(md_use) == "id"] <- "sample_id"
  long_tbl <- dplyr::left_join(long_tbl, md_use, by = "sample_id")
  long_tbl$stress <- long_tbl$stress_type
  long_tbl$stress_type <- NULL

  long_tbl$pair_key <- paste(long_tbl$fp_peak, long_tbl$gene_key, sep = "|")
  long_tbl <- tibble::as_tibble(long_tbl[, c("pair_key", "fp_peak", "gene_key", "sample_id", "expr", "tf_fp", "cell", "stress")])

  # Coerce to numeric (pivot_longer can yield character when inputs aren't numeric).
  long_tbl$expr  <- suppressWarnings(as.numeric(long_tbl$expr))
  long_tbl$tf_fp <- suppressWarnings(as.numeric(long_tbl$tf_fp))
  long_tbl <- long_tbl[is.finite(long_tbl$expr) & is.finite(long_tbl$tf_fp), , drop = FALSE]
  long_tbl <- long_tbl[!is.na(long_tbl$cell) & !is.na(long_tbl$stress), , drop = FALSE]
  if (!nrow(long_tbl)) {
    cli::cli_abort("No usable rows after coercion/filtering (check expr/tf_fp numeric and sample_md join for cell/stress).")
  }

  split_list <- split(long_tbl, long_tbl$pair_key)

  # ---- parallel selection ----
  if (is.null(n_workers) || !is.finite(n_workers) || n_workers < 1) {
    n_workers <- parallel::detectCores(logical = TRUE)
    if (is.na(n_workers) || n_workers < 1) n_workers <- 1L
  }
  n_workers <- as.integer(max(1L, n_workers))

  sysname <- Sys.info()[["sysname"]]
  can_fork <- isTRUE(prefer_fork) && !is.na(sysname) && sysname != "Windows"
  use_parallel <- n_workers > 1L && requireNamespace("parallel", quietly = TRUE)

  # ---- worker: one pair_key at a time ----
  fit_one_pair <- function(pair_key) {
    df_gene <- split_list[[pair_key]] %||% NULL
    fp_peak_i <- sub("\\|.*$", "", pair_key)
    gene_i <- sub("^.*\\|", "", pair_key)

    skip_reason <- NULL
    if (is.null(df_gene) || !nrow(df_gene)) skip_reason <- "empty df"
    else if (nrow(df_gene) < as.integer(min_samples)) skip_reason <- "nrow < min_samples"
    else if (stats::sd(df_gene$expr, na.rm = TRUE) == 0) skip_reason <- "expr has zero variance"
    else if (all(is.na(df_gene$tf_fp))) skip_reason <- "tf_fp all NA"
    else if (stats::sd(df_gene$tf_fp, na.rm = TRUE) == 0) skip_reason <- "tf_fp has zero variance"

    if (!is.null(skip_reason)) {
          return(list(
            rows = tibble::tibble(
              pair_key = pair_key,
              fp_peak = fp_peak_i,
              gene_key = gene_i,
              ok = FALSE,
              error = paste0("skipped: ", skip_reason),
              method = NA_character_,
              bs = NA_character_,
              k_requested = NA_integer_,
              k_used = NA_integer_,
              n_rows = if (is.null(df_gene)) 0L else nrow(df_gene),
              n_samples = if (is.null(df_gene)) 0L else length(unique(df_gene$sample_id)),
              edf = NA_real_,
              smooth_p = NA_real_,
              r_sq = NA_real_,
              dev_expl = NA_real_,
              AIC = NA_real_,
              BIC = NA_real_,
              gamcheck_k_prime = NA_real_,
              gamcheck_edf = NA_real_,
              gamcheck_k_index = NA_real_,
              gamcheck_p_value = NA_real_,
              gamcheck_output = NA_character_
            ),
            models = NULL
          ))
        }

    if (standardize != "none") {
      if (standardize %in% c("predictors", "all")) {
        sd_tf_fp <- stats::sd(df_gene$tf_fp, na.rm = TRUE)
        if (is.finite(sd_tf_fp) && sd_tf_fp > 0) df_gene$tf_fp <- as.numeric(scale(df_gene$tf_fp))
      }
      if (standardize == "all") {
        sd_expr <- stats::sd(df_gene$expr, na.rm = TRUE)
        if (is.finite(sd_expr) && sd_expr > 0) df_gene$expr <- as.numeric(scale(df_gene$expr))
      }
    }

    rows <- list()
    models_local <- if (isTRUE(keep_models)) list() else NULL
    idx <- 0L

    for (method in methods) {
      for (bs in bs_values) {
        for (k in k_values) {
          idx <- idx + 1L
          fit_obj <- tryCatch(
            fit_fp_gam_one(df_gene, k = k, bs = bs, method = method),
            error = function(e) e
          )

          if (inherits(fit_obj, "error")) {
            rows[[idx]] <- tibble::tibble(
              pair_key = pair_key,
              fp_peak = fp_peak_i,
              gene_key = gene_i,
              ok = FALSE,
              error = conditionMessage(fit_obj),
              method = method,
              bs = bs,
              k_requested = as.integer(k),
              k_used = NA_integer_,
              n_rows = nrow(df_gene),
              n_samples = length(unique(df_gene$sample_id)),
              edf = NA_real_,
              smooth_p = NA_real_,
              r_sq = NA_real_,
              dev_expl = NA_real_,
              AIC = NA_real_,
              BIC = NA_real_,
              gamcheck_k_prime = NA_real_,
              gamcheck_edf = NA_real_,
              gamcheck_k_index = NA_real_,
              gamcheck_p_value = NA_real_,
              gamcheck_output = NA_character_
            )
            next
          }

          fit <- fit_obj$fit
          sm <- tryCatch(summary(fit), error = function(e) NULL)
          edf <- NA_real_
          smooth_p <- NA_real_
          r_sq <- NA_real_
          dev_expl <- NA_real_
          if (!is.null(sm) && !is.null(sm$s.table) && is.matrix(sm$s.table) && nrow(sm$s.table) >= 1) {
            if ("edf" %in% colnames(sm$s.table)) edf <- as.numeric(sm$s.table[1, "edf"])
            smooth_p <- as.numeric(sm$s.table[1, ncol(sm$s.table)])
          }
          if (!is.null(sm) && !is.null(sm$r.sq)) r_sq <- as.numeric(sm$r.sq)
          if (!is.null(sm) && !is.null(sm$dev.expl)) dev_expl <- as.numeric(sm$dev.expl)

          gchk <- if (isTRUE(run_gam_check)) {
            .gam_check_metrics(fit, keep_output = keep_gamcheck_output)
          } else {
            list(ok = NA, error = NA, k_prime = NA_real_, edf = NA_real_, k_index = NA_real_, k_p_value = NA_real_, output = NA_character_)
          }

          rows[[idx]] <- tibble::tibble(
            pair_key = pair_key,
            fp_peak = fp_peak_i,
            gene_key = gene_i,
            ok = TRUE,
            error = NA_character_,
            method = method,
            bs = bs,
            k_requested = as.integer(fit_obj$k_requested),
            k_used = as.integer(fit_obj$k_used),
            n_rows = nrow(df_gene),
            n_samples = length(unique(df_gene$sample_id)),
            edf = edf,
            smooth_p = as.numeric(smooth_p),
            r_sq = r_sq,
            dev_expl = dev_expl,
            AIC = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
            BIC = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
            gamcheck_k_prime = as.numeric(gchk$k_prime),
            gamcheck_edf = as.numeric(gchk$edf),
            gamcheck_k_index = as.numeric(gchk$k_index),
            gamcheck_p_value = as.numeric(gchk$k_p_value),
            gamcheck_output = as.character(gchk$output)
          )

          if (isTRUE(keep_models)) {
            models_local[[paste(pair_key, method, bs, k, sep = "||")]] <- fit
          }
        }
      }
    }

    list(rows = dplyr::bind_rows(rows), models = models_local)
  }

  # ---- run over pairs ----
  keys <- pair_tbl$pair_key
  out_list <- NULL
  if (use_parallel && can_fork) {
    out_list <- parallel::mclapply(keys, fit_one_pair, mc.cores = n_workers, mc.silent = TRUE)
  } else if (use_parallel && n_workers > 1L) {
    cl <- parallel::makeCluster(n_workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("split_list", "min_samples", "standardize", "methods", "bs_values", "k_values",
                  "run_gam_check", "keep_models", "fit_fp_gam_one", ".gam_check_metrics", "%||%"),
      envir = environment()
    )
    out_list <- parallel::parLapply(cl, keys, fun = fit_one_pair)
  } else {
    out_list <- lapply(keys, fit_one_pair)
  }

  res_tbl <- dplyr::bind_rows(lapply(out_list, `[[`, "rows"))
  models_out <- NULL
  if (isTRUE(keep_models)) {
    models_out <- unlist(lapply(out_list, `[[`, "models"), recursive = FALSE)
  }

  # Carry link-level metadata (e.g., ko_group/log2fc) into the per-model results.
  pair_meta_cols <- c(
    "tf", "fp_peak", "gene_key",
    "ko_group", "log2fc",
    "r_fp", "p_fp", "p_adj_fp",
    "motifs", "tfs", "n_fp", "atac_peak"
  )
  pair_meta_cols <- intersect(pair_meta_cols, names(pair_tbl))
  if (all(c("fp_peak", "gene_key") %in% pair_meta_cols) && nrow(res_tbl)) {
    pair_meta <- pair_tbl |>
      dplyr::select(dplyr::all_of(pair_meta_cols)) |>
      dplyr::distinct(.data$fp_peak, .data$gene_key, .keep_all = TRUE)

    res_tbl <- dplyr::left_join(
      res_tbl,
      pair_meta,
      by = c("fp_peak", "gene_key")
    )
  }

  list(results = res_tbl, pair_tbl = pair_tbl, models = models_out)
}

# ---- helpers: make gam.check columns actionable -----------------------------
.cap_p_floor <- function(p, floor = 1e-300) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  p[p < 0] <- NA_real_
  p[p == 0] <- floor
  p[p < floor] <- floor
  p[p > 1] <- 1
  p
}

annotate_fp_gam_grid_results <- function(results_tbl,
                                        p_floor = 1e-300,
                                        k_index_good = 0.9,
                                        k_p_good = 0.05,
                                        edf_frac_warn = 0.8) {
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Requires packages: dplyr, tibble")
  }
  if (!is.data.frame(results_tbl)) stop("results_tbl must be a data.frame")

  out <- tibble::as_tibble(results_tbl)

  out <- out |>
    dplyr::mutate(
      smooth_p_cap = .cap_p_floor(.data$smooth_p, floor = p_floor),
      smooth_mlog10p = -log10(.data$smooth_p_cap),
      gamcheck_p_cap = .cap_p_floor(.data$gamcheck_p_value, floor = p_floor),
      gamcheck_mlog10p = -log10(.data$gamcheck_p_cap),
      gamcheck_edf_frac = dplyr::if_else(
        is.finite(.data$gamcheck_edf) & is.finite(.data$gamcheck_k_prime) & .data$gamcheck_k_prime > 0,
        .data$gamcheck_edf / .data$gamcheck_k_prime,
        NA_real_
      ),
      k_too_low_flag = is.finite(.data$gamcheck_k_index) &
        .data$gamcheck_k_index < 1 &
        is.finite(.data$gamcheck_p_cap) &
        .data$gamcheck_p_cap < k_p_good,
      k_borderline_flag = is.finite(.data$gamcheck_k_index) &
        .data$gamcheck_k_index < k_index_good &
        is.finite(.data$gamcheck_p_cap) &
        .data$gamcheck_p_cap < k_p_good,
      edf_near_k_flag = is.finite(.data$gamcheck_edf_frac) & .data$gamcheck_edf_frac >= edf_frac_warn
    )

  out <- out |>
    dplyr::group_by(.data$pair_key) |>
    dplyr::mutate(
      delta_AIC = .data$AIC - min(.data$AIC, na.rm = TRUE),
      delta_BIC = .data$BIC - min(.data$BIC, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  out
}

pick_best_fp_gam_grid_models <- function(results_tbl,
                                         prefer = c("AIC", "BIC"),
                                         require_ok = TRUE,
                                         require_k_ok = FALSE,
                                         k_index_min = 0.9,
                                         k_p_min = 0.05) {
  prefer <- match.arg(prefer)
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Requires packages: dplyr, tibble")
  }
  tbl <- tibble::as_tibble(results_tbl)
  if (isTRUE(require_ok) && "ok" %in% names(tbl)) tbl <- dplyr::filter(tbl, .data$ok)
  if (isTRUE(require_k_ok)) {
    tbl <- dplyr::filter(
      tbl,
      is.finite(.data$gamcheck_k_index) & .data$gamcheck_k_index >= k_index_min,
      is.finite(.data$gamcheck_p_value) & .data$gamcheck_p_value >= k_p_min
    )
  }

  score_col <- if (prefer == "AIC") "AIC" else "BIC"
  tbl |>
    dplyr::group_by(.data$pair_key) |>
    dplyr::slice_min(.data[[score_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

pick_optimal_k_by_gamcheck <- function(results_tbl,
                                       score = c("AIC", "BIC"),
                                       group_cols = c("pair_key", "method", "bs"),
                                       k_col = "k_used",
                                       k_index_min = 0.9,
                                       k_p_min = 0.05) {
  score <- match.arg(score)
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tibble", quietly = TRUE)) {
    stop("Requires packages: dplyr, tibble")
  }

  tbl <- tibble::as_tibble(results_tbl)
  miss <- setdiff(c(group_cols, k_col, "ok", "gamcheck_k_index", "gamcheck_p_value", score), names(tbl))
  if (length(miss)) stop("results_tbl missing columns: ", paste(miss, collapse = ", "))

  score_col <- score

  tbl <- tbl |>
    dplyr::filter(.data$ok) |>
    dplyr::mutate(
      .k_val = as.integer(.data[[k_col]]),
      .k_ok = is.finite(.data$gamcheck_k_index) & .data$gamcheck_k_index >= k_index_min &
        is.finite(.data$gamcheck_p_value) & .data$gamcheck_p_value >= k_p_min
    )

  pick_one_group <- function(df) {
    df <- df[is.finite(df$.k_val), , drop = FALSE]
    if (!nrow(df)) return(df[0, , drop = FALSE])

    any_ok <- any(df$.k_ok, na.rm = TRUE)

    if (any_ok) {
      df_ok <- df[df$.k_ok, , drop = FALSE]
      k_star <- min(df_ok$.k_val, na.rm = TRUE)
      df_k <- df_ok[df_ok$.k_val == k_star, , drop = FALSE]
      df_k <- df_k[order(df_k[[score_col]]), , drop = FALSE]
      out <- df_k[1, , drop = FALSE]
      out$k_ok_any <- TRUE
      out$k_ok_selected <- TRUE
      out$k_selected_reason <- "min_k_passing_kcheck"
      return(out)
    }

    # If nothing passes: take the largest k available (best chance to fix k-index issue),
    # then choose best score among those rows.
    k_max <- max(df$.k_val, na.rm = TRUE)
    df_k <- df[df$.k_val == k_max, , drop = FALSE]
    df_k <- df_k[order(df_k[[score_col]]), , drop = FALSE]
    out <- df_k[1, , drop = FALSE]
    out$k_ok_any <- FALSE
    out$k_ok_selected <- FALSE
    out$k_selected_reason <- "max_k_none_passed"
    out
  }

  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::group_modify(~ pick_one_group(.x)) |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::any_of(c(".k_val", ".k_ok")))
}

# ---- main runner -------------------------------------------------------------
run_fp_gam_for_all_pairs_cached_standalone <- function(fp_score,
                                                       rna_tbl,
                                                       sample_md,
                                                       tf_gene_links_gh,
                                                       cache_dir,
                                                       r_fp_cut        = 0.3,
                                                       p_fp_cut        = 0.05,
                                                       min_samples     = 6L,
                                                       overwrite       = FALSE,
                                                       mc_cores        = 30L,
                                                       pairs_per_block = 2000L,
                                                       standardize     = c("all", "predictors", "none"),
                                                       gam_k           = 5L) {
  standardize <- match.arg(standardize)

  # ---- package checks ----
  need_pkgs <- c("cli", "dplyr", "tidyr", "tibble", "readr", "mgcv")
  missing <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    cli::cli_abort(c(
      "Missing required packages.",
      "x" = "Install: {.val {missing}}"
    ))
  }

  # ---- input checks ----
  if (!is.data.frame(fp_score) || !is.data.frame(rna_tbl) || !is.data.frame(sample_md) || !is.data.frame(tf_gene_links_gh)) {
    cli::cli_abort("All inputs fp_score, rna_tbl, sample_md, tf_gene_links_gh must be data.frames / tibbles.")
  }
  if (!all(c("id", "cell", "stress_type") %in% names(sample_md))) {
    cli::cli_abort("sample_md must have columns: id, cell, stress_type")
  }
  if (!"peak_ID" %in% names(fp_score)) {
    cli::cli_abort("fp_score must have column: peak_ID")
  }
  if (!"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("rna_tbl must have column: HGNC")
  }
  if (!all(c("fp_peak", "gene_key", "r_fp", "p_fp") %in% names(tf_gene_links_gh))) {
    cli::cli_abort("tf_gene_links_gh must have columns: fp_peak, gene_key, r_fp, p_fp")
  }

  sample_ids <- as.character(sample_md$id)
  if (!length(sample_ids)) {
    cli::cli_abort("sample_md$id is empty; no samples to model.")
  }

  # Ensure sample columns exist in both matrices
  miss_fp  <- setdiff(sample_ids, names(fp_score))
  miss_rna <- setdiff(sample_ids, names(rna_tbl))
  if (length(miss_fp)) {
    cli::cli_abort(c("Some sample ids are missing from fp_score columns.", "x" = "Missing: {.val {miss_fp}}"))
  }
  if (length(miss_rna)) {
    cli::cli_abort(c("Some sample ids are missing from rna_tbl columns.", "x" = "Missing: {.val {miss_rna}}"))
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- helpers (local) ----
  run_fp_gam_for_unit <- function(df_gene,
                                  min_samples = 6L,
                                  standardize = c("all", "predictors", "none"),
                                  gam_k = 5L) {
    standardize <- match.arg(standardize)

    if (is.null(df_gene) || !nrow(df_gene)) {
      return(simpleError("FP:GAM skipped: empty df"))
    }
    if (nrow(df_gene) < min_samples) {
      return(simpleError("FP:GAM skipped: nrow < min_samples"))
    }
    if (stats::sd(df_gene$expr, na.rm = TRUE) == 0) {
      return(simpleError("FP:GAM skipped: expr has zero variance"))
    }
    if (all(is.na(df_gene$tf_fp))) {
      return(simpleError("FP:GAM skipped: tf_fp all NA"))
    }
    if (stats::sd(df_gene$tf_fp, na.rm = TRUE) == 0) {
      return(simpleError("FP:GAM skipped: tf_fp has zero variance"))
    }

    if (standardize != "none") {
      if (standardize %in% c("predictors", "all")) {
        sd_tf_fp <- stats::sd(df_gene$tf_fp, na.rm = TRUE)
        if (is.finite(sd_tf_fp) && sd_tf_fp > 0) {
          df_gene$tf_fp <- as.numeric(scale(df_gene$tf_fp))
        }
      }
      if (standardize == "all") {
        sd_expr <- stats::sd(df_gene$expr, na.rm = TRUE)
        if (is.finite(sd_expr) && sd_expr > 0) {
          df_gene$expr <- as.numeric(scale(df_gene$expr))
        }
      }
    }

    df_gene$cell   <- factor(df_gene$cell)
    df_gene$stress <- factor(df_gene$stress)

    tryCatch(
      tfs_to_target_genes_fp_tobias_gam(df_gene, k = gam_k),
      error = function(e) e
    )
  }

  build_fp_gene_pair_tbl <- function(tf_gene_links_gh, r_fp_cut = 0.3, p_fp_cut = 0.05) {
    pair_tbl <- tf_gene_links_gh |>
      dplyr::filter(
        is.finite(.data$r_fp), is.finite(.data$p_fp),
        .data$r_fp > r_fp_cut,
        .data$p_fp < p_fp_cut,
        !is.na(.data$fp_peak), .data$fp_peak != "",
        !is.na(.data$gene_key), .data$gene_key != ""
      ) |>
      dplyr::distinct(.data$fp_peak, .data$gene_key)

    pair_tbl$pair_key <- paste(pair_tbl$fp_peak, pair_tbl$gene_key, sep = "|")
    pair_tbl
  }

  prepare_block_long_tbl <- function(pair_block, fp_score, rna_tbl, sample_md) {
    sample_ids <- as.character(sample_md$id)

    peaks_use <- unique(pair_block$fp_peak)
    fp_sub <- fp_score[fp_score$peak_ID %in% peaks_use, c("peak_ID", sample_ids), drop = FALSE]
    if (!nrow(fp_sub)) return(tibble::tibble())

    fp_long <- tidyr::pivot_longer(
      data      = fp_sub,
      cols      = dplyr::all_of(sample_ids),
      names_to  = "sample_id",
      values_to = "tf_fp"
    )
    names(fp_long)[names(fp_long) == "peak_ID"] <- "fp_peak"

    fp_long <- dplyr::left_join(
      fp_long,
      pair_block[, c("fp_peak", "gene_key"), drop = FALSE],
      by = "fp_peak"
    )
    fp_long <- fp_long[!is.na(fp_long$gene_key) & fp_long$gene_key != "", , drop = FALSE]
    if (!nrow(fp_long)) return(tibble::tibble())

    genes_use <- unique(pair_block$gene_key)
    rna_sub <- rna_tbl[rna_tbl$HGNC %in% genes_use, c("HGNC", sample_ids), drop = FALSE]
    if (!nrow(rna_sub)) return(tibble::tibble())

    rna_long <- tidyr::pivot_longer(
      data      = rna_sub,
      cols      = dplyr::all_of(sample_ids),
      names_to  = "sample_id",
      values_to = "expr"
    )
    names(rna_long)[names(rna_long) == "HGNC"] <- "gene_key"

    long_tbl <- dplyr::inner_join(fp_long, rna_long, by = c("gene_key", "sample_id"))
    if (!nrow(long_tbl)) return(tibble::tibble())

    sample_md_use <- sample_md[, c("id", "cell", "stress_type"), drop = FALSE]
    names(sample_md_use)[names(sample_md_use) == "id"] <- "sample_id"

    long_tbl <- dplyr::left_join(long_tbl, sample_md_use, by = "sample_id")
    long_tbl$stress <- long_tbl$stress_type
    long_tbl$stress_type <- NULL

    long_tbl$pair_key <- paste(long_tbl$fp_peak, long_tbl$gene_key, sep = "|")

    tibble::as_tibble(long_tbl[, c("pair_key", "fp_peak", "gene_key", "sample_id", "expr", "tf_fp", "cell", "stress")])
  }

  # ---- pair table ----
  pair_tbl <- build_fp_gene_pair_tbl(tf_gene_links_gh, r_fp_cut = r_fp_cut, p_fp_cut = p_fp_cut)

  if (!nrow(pair_tbl)) {
    cli::cli_warn("No (fp_peak, gene_key) pairs after filtering; nothing to run.")
    return(list(pair_tbl = pair_tbl, block_files_csv = character(0), block_files_rds = character(0)))
  }

  # ---- resume bookkeeping ----
  pat_csv <- "^FPGENE_FPGAM_block_\\d+\\.csv$"
  existing_csv <- list.files(cache_dir, pattern = pat_csv, full.names = TRUE)

  if (isTRUE(overwrite) && length(existing_csv)) {
    unlink(existing_csv)
    existing_csv <- character(0)
  }

  done_keys <- character(0)
  if (!isTRUE(overwrite) && length(existing_csv)) {
    for (bf in existing_csv) {
      kvec <- tryCatch({
        tmp <- readr::read_csv(
          bf,
          col_types = readr::cols_only(pair_key = readr::col_character()),
          show_col_types = FALSE
        )
        unique(tmp$pair_key)
      }, error = function(e) character(0))
      done_keys <- unique(c(done_keys, kvec))
    }
  }

  pair_tbl_todo <- pair_tbl[!pair_tbl$pair_key %in% done_keys, , drop = FALSE]
  if (!nrow(pair_tbl_todo)) {
    pat_rds <- "^FPGENE_FPGAM_block_\\d+_models\\.rds$"
    existing_rds <- list.files(cache_dir, pattern = pat_rds, full.names = TRUE)
    return(list(pair_tbl = pair_tbl, block_files_csv = sort(existing_csv), block_files_rds = sort(existing_rds)))
  }

  # ---- parallel selection ----
  use_parallel <- !is.null(mc_cores) &&
    is.numeric(mc_cores) &&
    as.integer(mc_cores) > 1L &&
    requireNamespace("parallel", quietly = TRUE)

  mc_cores <- if (use_parallel) as.integer(mc_cores) else 1L

  block_files_csv <- existing_csv
  block_files_rds <- character(0)

  block_id_offset <- length(existing_csv)
  block_counter <- 0L

  for (start_idx in seq(1L, nrow(pair_tbl_todo), by = as.integer(pairs_per_block))) {
    end_idx <- min(start_idx + as.integer(pairs_per_block) - 1L, nrow(pair_tbl_todo))
    pair_block <- pair_tbl_todo[start_idx:end_idx, , drop = FALSE]

    cli::cli_inform("FPGENE FP:GAM: pairs {start_idx}-{end_idx} of {nrow(pair_tbl_todo)}")

    long_tbl <- prepare_block_long_tbl(pair_block, fp_score, rna_tbl, sample_md)

    split_list <- if (nrow(long_tbl)) split(long_tbl, long_tbl$pair_key) else list()
    keys_block <- pair_block$pair_key

    fit_one_key <- function(k) {
      dfk <- split_list[[k]] %||% NULL
      mr  <- run_fp_gam_for_unit(dfk, min_samples = min_samples, standardize = standardize, gam_k = gam_k)

      parts <- strsplit(k, "\\|", fixed = FALSE)[[1]]
      fp_peak_i <- parts[1]
      gene_i    <- parts[2]

      status <- if (inherits(mr, "error")) {
        tibble::tibble(
          pair_key    = k,
          fp_peak     = fp_peak_i,
          gene_key    = gene_i,
          model       = "tfs_to_target_genes_fp_tobias_gam",
          ok          = FALSE,
          model_class = NA_character_,
          error       = conditionMessage(mr),
          n_rows      = if (is.null(dfk)) 0L else nrow(dfk),
          n_samples   = if (is.null(dfk)) 0L else length(unique(dfk$sample_id))
        )
      } else {
        tibble::tibble(
          pair_key    = k,
          fp_peak     = fp_peak_i,
          gene_key    = gene_i,
          model       = "tfs_to_target_genes_fp_tobias_gam",
          ok          = TRUE,
          model_class = mr$model_class %||% NA_character_,
          error       = NA_character_,
          n_rows      = if (is.null(dfk)) 0L else nrow(dfk),
          n_samples   = if (is.null(dfk)) 0L else length(unique(dfk$sample_id))
        )
      }

      list(key = k, model = mr, status = status)
    }

    if (use_parallel) {
      out_list <- parallel::mclapply(keys_block, fit_one_key, mc.cores = mc_cores)
    } else {
      out_list <- lapply(keys_block, fit_one_key)
    }

    status_df <- dplyr::bind_rows(lapply(out_list, `[[`, "status"))
    models_block <- setNames(
      lapply(out_list, `[[`, "model"),
      vapply(out_list, `[[`, character(1), "key")
    )

    block_counter <- block_counter + 1L
    block_id <- block_id_offset + block_counter

    f_csv <- file.path(cache_dir, sprintf("FPGENE_FPGAM_block_%03d.csv", block_id))
    f_rds <- file.path(cache_dir, sprintf("FPGENE_FPGAM_block_%03d_models.rds", block_id))

    readr::write_csv(status_df, f_csv)
    saveRDS(list(models = models_block), f_rds)

    block_files_csv <- c(block_files_csv, f_csv)
    block_files_rds <- c(block_files_rds, f_rds)

    rm(out_list, status_df, models_block, long_tbl, split_list)
    gc(FALSE)
  }

  list(
    pair_tbl        = pair_tbl,
    block_files_csv = sort(unique(block_files_csv)),
    block_files_rds = sort(unique(block_files_rds))
  )
}

# Run (filtered + distinct fp_peak×gene_key
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"
cache_dir_fp_gam <- file.path(ko_dir, "TF_FP_GENE_GAM_cache")

# results_fp_gam_pairs <- run_fp_gam_for_all_pairs_cached_standalone(
#   fp_score         = grn_set$fp_score,
#   rna_tbl          = grn_set$rna,
#   sample_md        = grn_set$sample_metadata_used,
#   tf_gene_links_gh = tf_gene_links_gh,
#   cache_dir        = cache_dir_fp_gam,
#   r_fp_cut         = 0.3,
#   p_fp_cut         = 0.05,
#   min_samples      = 6L,
#   overwrite        = FALSE,
#   mc_cores         = 30L,
#   pairs_per_block  = 2000L,
#   standardize      = "all",
#   gam_k            = 5L
# )


# status_all <- dplyr::bind_rows(lapply(results_fp_gam_pairs$block_files_csv, function(f) {
#   readr::read_csv(f, show_col_types = FALSE)
# }))

# Overall completion
# status_all |>
#   dplyr::summarise(
#     n_pairs = dplyr::n(),
#     ok_rate = mean(ok),
#     n_ok    = sum(ok),
#     n_fail  = sum(!ok)
#   )

# Load and merge FPGENE_FPGAM_block_*_models.rds into ONE variable: all_models

load_fp_gam_models_cache <- function(cache_dir,
                                     pattern = "^FPGENE_FPGAM_block_\\d+_models\\.rds$",
                                     strict = TRUE) {
  if (!dir.exists(cache_dir)) {
    cli::cli_abort("cache_dir does not exist: {.path {cache_dir}}")
  }

  rds_files <- sort(list.files(cache_dir, pattern = pattern, full.names = TRUE))
  if (!length(rds_files)) {
    cli::cli_abort("No model RDS files matched in cache_dir: {.path {cache_dir}}")
  }

  all_models <- list()
  for (f in rds_files) {
    obj <- readRDS(f)

    if (is.list(obj) && !is.null(obj$models) && is.list(obj$models)) {
      block_models <- obj$models
    } else if (is.list(obj)) {
      # fallback: if someone saved the models directly
      block_models <- obj
    } else {
      if (isTRUE(strict)) {
        cli::cli_abort("Unexpected RDS structure in {.path {f}}")
      } else {
        next
      }
    }

    # later blocks overwrite earlier blocks if duplicate keys exist
    dup <- intersect(names(all_models), names(block_models))
    if (length(dup)) {
      cli::cli_warn(c(
        "Duplicate pair_key(s) encountered; later blocks will overwrite earlier ones.",
        "i" = "{length(dup)} duplicates (e.g., {.val {utils::head(dup, 5)}})"
      ))
    }

    all_models <- c(all_models, block_models)
  }

  all_models
}

# ---- usage
# cache_dir_fp_gam <- file.path(
#   "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction",
#   "TF_FP_GENE_GAM_cache"
# )
#
# all_fp_gam_models <- load_fp_gam_models_cache(cache_dir_fp_gam)

# Parallel loader: 1 row per (pair_key) across all FPGENE_FPGAM block RDS files.
# - Uses forked multicore (mclapply) on Linux/macOS when available
# - Falls back to PSOCK (parLapply) on Windows
# - Robust to per-file failures (returns ok=FALSE + error message)

.read_one_fp_gam_rds_to_tbl <- function(f,
                                        keep_vectors = FALSE,
                                        strict = TRUE) {
  .pick_row <- function(tbl, term_target = NULL) {
    if (is.null(tbl) || !inherits(tbl, "data.frame") || nrow(tbl) < 1) return(NULL)
    if (!is.null(term_target) && "term" %in% names(tbl)) {
      idx <- match(term_target, tbl[["term"]])
      if (!is.na(idx)) return(tbl[idx, , drop = FALSE])
    }
    tbl[1, , drop = FALSE]
  }

  .get_scalar <- function(tbl, col) {
    if (is.null(tbl)) return(NA_real_)
    if (!col %in% names(tbl)) return(NA_real_)
    as.numeric(tbl[[col]][1])
  }

  # Read with failure capture so one bad file doesn't kill the whole job.
  obj <- tryCatch(
    readRDS(f),
    error = function(e) e
  )
  if (inherits(obj, "error")) {
    return(tibble::tibble(
      rds_file = f,
      pair_key = NA_character_,
      peak_id  = NA_character_,
      gene_key = NA_character_,
      ok       = FALSE,
      error    = paste0("readRDS failed: ", conditionMessage(obj)),
      model_class = NA_character_,
      coef_term = NA_character_,
      coef_edf = NA_real_,
      coef_ref_df = NA_real_,
      coef_statistic = NA_real_,
      coef_p_value = NA_real_,
      glance_df = NA_real_,
      glance_logLik = NA_real_,
      glance_AIC = NA_real_,
      glance_BIC = NA_real_,
      glance_deviance = NA_real_,
      glance_df_residual = NA_real_,
      glance_nobs = NA_integer_,
      glance_adj_r_squared = NA_real_,
      glance_npar = NA_integer_
    ))
  }

  block_models <- NULL
  if (is.list(obj) && !is.null(obj$models) && is.list(obj$models)) {
    block_models <- obj$models
  } else if (is.list(obj)) {
    block_models <- obj
  }

  if (is.null(block_models) || !is.list(block_models)) {
    msg <- paste0("Unexpected RDS structure in file: ", f)
    if (isTRUE(strict)) {
      return(tibble::tibble(
        rds_file = f,
        pair_key = NA_character_,
        peak_id  = NA_character_,
        gene_key = NA_character_,
        ok       = FALSE,
        error    = msg,
        model_class = NA_character_,
        coef_term = NA_character_,
        coef_edf = NA_real_,
        coef_ref_df = NA_real_,
        coef_statistic = NA_real_,
        coef_p_value = NA_real_,
        glance_df = NA_real_,
        glance_logLik = NA_real_,
        glance_AIC = NA_real_,
        glance_BIC = NA_real_,
        glance_deviance = NA_real_,
        glance_df_residual = NA_real_,
        glance_nobs = NA_integer_,
        glance_adj_r_squared = NA_real_,
        glance_npar = NA_integer_
      ))
    }
    return(tibble::tibble())
  }

  if (!length(block_models)) return(tibble::tibble())

  rows <- lapply(names(block_models), function(pair_key) {
    m <- block_models[[pair_key]]

    parts <- strsplit(pair_key, "\\|", fixed = FALSE)[[1]]
    peak_id <- if (length(parts) >= 1) parts[[1]] else NA_character_
    gene_key <- if (length(parts) >= 2) parts[[2]] else NA_character_

    err_msg <- NA_character_
    ok <- TRUE
    if (is.null(m)) {
      ok <- FALSE
      err_msg <- "NULL model object"
    } else if (is.list(m) && !is.null(m$error)) {
      ok <- FALSE
      err_msg <- as.character(m$error)[1]
    }

    coef_row   <- if (ok) .pick_row(m$coef,   term_target = "s(tf_fp)") else NULL
    glance_row <- if (ok) .pick_row(m$glance, term_target = NULL)       else NULL

    out <- tibble::tibble(
      rds_file   = f,
      pair_key   = pair_key,
      peak_id    = peak_id,
      gene_key   = gene_key,
      ok         = ok,
      error      = err_msg,
      model_class = if (ok && !is.null(m$model_class)) as.character(m$model_class)[1] else NA_character_,

      coef_term      = if (!is.null(coef_row) && "term" %in% names(coef_row)) as.character(coef_row[["term"]][1]) else NA_character_,
      coef_edf       = .get_scalar(coef_row, "edf"),
      coef_ref_df    = .get_scalar(coef_row, "ref.df"),
      coef_statistic = .get_scalar(coef_row, "statistic"),
      coef_p_value   = .get_scalar(coef_row, "p.value"),

      glance_df            = .get_scalar(glance_row, "df"),
      glance_logLik        = .get_scalar(glance_row, "logLik"),
      glance_AIC           = .get_scalar(glance_row, "AIC"),
      glance_BIC           = .get_scalar(glance_row, "BIC"),
      glance_deviance      = .get_scalar(glance_row, "deviance"),
      glance_df_residual   = .get_scalar(glance_row, "df.residual"),
      glance_nobs          = if (is.null(glance_row) || !"nobs" %in% names(glance_row)) NA_integer_ else as.integer(glance_row[["nobs"]][1]),
      glance_adj_r_squared = .get_scalar(glance_row, "adj.r.squared"),
      glance_npar          = if (is.null(glance_row) || !"npar" %in% names(glance_row)) NA_integer_ else as.integer(glance_row[["npar"]][1])
    )

    if (isTRUE(keep_vectors)) {
      out[["fitted"]] <- list(if (ok) m$fitted else NULL)
      out[["resid"]]  <- list(if (ok) m$resid  else NULL)
    }

    out
  })

  dplyr::bind_rows(rows)
}

read_fp_gam_models_as_tibble_parallel <- function(cache_dir,
                                                  pattern = "^FPGENE_FPGAM_block_\\d+_models\\.rds$",
                                                  keep_vectors = FALSE,
                                                  strict = TRUE,
                                                  n_workers = NULL,
                                                  prefer_fork = TRUE) {
  if (!dir.exists(cache_dir)) {
    cli::cli_abort("cache_dir does not exist: {.path {cache_dir}}")
  }

  rds_files <- sort(list.files(cache_dir, pattern = pattern, full.names = TRUE))
  if (!length(rds_files)) {
    cli::cli_abort("No model RDS files matched in cache_dir: {.path {cache_dir}}")
  }

  if (is.null(n_workers) || !is.finite(n_workers) || n_workers < 1) {
    n_workers <- parallel::detectCores(logical = TRUE)
    if (is.na(n_workers) || n_workers < 1) n_workers <- 1L
  }
  n_workers <- as.integer(max(1L, n_workers))

  sysname <- Sys.info()[["sysname"]]
  can_fork <- isTRUE(prefer_fork) && !is.na(sysname) && sysname != "Windows"

  # ---- parallel over files (best bang-for-buck) ----
  tbl_list <- NULL

  if (can_fork && n_workers > 1L) {
    tbl_list <- parallel::mclapply(
      rds_files,
      FUN = .read_one_fp_gam_rds_to_tbl,
      keep_vectors = keep_vectors,
      strict = strict,
      mc.cores = n_workers
    )
  } else if (n_workers > 1L) {
    cl <- parallel::makeCluster(n_workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

    # Export only what workers need.
    parallel::clusterExport(
      cl,
      varlist = c(".read_one_fp_gam_rds_to_tbl"),
      envir = environment()
    )

    tbl_list <- parallel::parLapply(
      cl,
      rds_files,
      fun = function(f, keep_vectors, strict) {
        .read_one_fp_gam_rds_to_tbl(f, keep_vectors = keep_vectors, strict = strict)
      },
      keep_vectors = keep_vectors,
      strict = strict
    )
  } else {
    tbl_list <- lapply(
      rds_files,
      FUN = .read_one_fp_gam_rds_to_tbl,
      keep_vectors = keep_vectors,
      strict = strict
    )
  }

  out <- dplyr::bind_rows(tbl_list)

  # Optional: stable ordering (file then pair_key)
  if (nrow(out) > 0) {
    out <- out[order(out[["rds_file"]], out[["pair_key"]]), , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}

# -------------------- usage (guarded; won't run on source) -------------------
# Enable explicitly via one of:
#   options(episcope.run_dev_examples = TRUE)
#   Sys.setenv(EPISCOPE_RUN_DEV_EXAMPLES = "1")
if (isTRUE(getOption("episcope.run_dev_examples")) ||
    identical(Sys.getenv("EPISCOPE_RUN_DEV_EXAMPLES"), "1")) {

  cache_dir_fp_gam <- file.path(
    "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction",
    "TF_FP_GENE_GAM_cache"
  )
  models_tbl <- read_fp_gam_models_as_tibble_parallel(cache_dir_fp_gam, n_workers = 30)
  str(models_tbl)

  # --- Inputs ---
  # tf_gene_links_gh: tibble with columns fp_peak, gene_key, tfs, r_fp, p_fp, p_adj_fp, ...
  # models_tbl: tibble with columns pair_key, peak_id, gene_key, coef_p_value, coef_edf, ok, ...

  # --- Use data.table for scale ---
  tf_dt <- data.table::as.data.table(tf_gene_links_gh)[
    ,
    .(fp_peak, gene_key, tfs, motifs, n_fp, r_fp, p_fp, p_adj_fp)
  ]

  m_dt <- data.table::as.data.table(models_tbl)[
    ok == TRUE,
    .(peak_id, gene_key, coef_p_value, coef_edf, glance_nobs, glance_adj_r_squared)
  ]
  data.table::setnames(m_dt, "peak_id", "fp_peak")

  data.table::setkey(tf_dt, fp_peak, gene_key)
  data.table::setkey(m_dt,  fp_peak, gene_key)

  # --- Join: adds GAM p-value to each TF edge (same fp_peak|gene_key replicated across TFs) ---
  tf_gam_dt <- m_dt[tf_dt, on = .(fp_peak, gene_key)]

  # --- Define "significant" flags (edit cutoffs as desired) ---
  alpha_fdr <- 0.05
  alpha_gam <- 0.05

  tf_gam_dt[
    ,
    `:=`(
      sig_corr_fdr = !is.na(p_adj_fp) & p_adj_fp <= alpha_fdr,
      sig_gam      = !is.na(coef_p_value) & coef_p_value <= alpha_gam,
      sig_both     = (!is.na(p_adj_fp) & p_adj_fp <= alpha_fdr) &
        (!is.na(coef_p_value) & coef_p_value <= alpha_gam)
    )
  ]

  # pull all TF edges for one fp_peak->gene link and inspect
  tf_gam_dt[fp_peak == "chr10:100229144-100229170" & gene_key == "ERLIN1"][
    order(p_adj_fp),
    .(fp_peak, gene_key, tfs, r_fp, p_fp, p_adj_fp, coef_p_value, coef_edf, glance_nobs, glance_adj_r_squared,
      sig_corr_fdr, sig_gam, sig_both)
  ]
}
