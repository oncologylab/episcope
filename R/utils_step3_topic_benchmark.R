# File: utils_step3_topic_benchmark.R
# Author: Yaoxiang Li
# Purpose: Package-native Module 3 topic benchmark review helpers.

.m3tb_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

.m3tb_safe_label <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x), perl = TRUE)
}

.m3tb_html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

.m3tb_relative_path <- function(paths, from_dir) {
  paths_n <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  from_n <- normalizePath(from_dir, winslash = "/", mustWork = FALSE)
  from_n <- paste0(sub("/+$", "", from_n), "/")
  out <- paths_n
  inside <- startsWith(paths_n, from_n)
  out[inside] <- substring(paths_n[inside], nchar(from_n) + 1L)
  out
}

.m3tb_review_tables_dir <- function(review_dir) {
  file.path(review_dir, "tables")
}

.m3tb_review_html_dir <- function(review_dir) {
  review_dir
}

.m3tb_review_read_dir <- function(topic_dir) {
  candidates <- c(
    file.path(topic_dir, "review", "tables"),
    file.path(topic_dir, "review", "csv"),
    file.path(topic_dir, "review_topic_experiments", "tables"),
    file.path(topic_dir, "review_topic_experiments", "csv")
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) hit[[1L]] else candidates[[1L]]
}

.m3tb_plot_dir <- function(review_dir) {
  file.path(review_dir, "assets")
}

.m3tb_method_dictionary <- function() {
  data.table::data.table(
    method = c(
      "condition_aggr_weight_lda",
      "condition_aggr_weight_vae_mlp",
      "condition_aggr_lda",
      "condition_aggr_multivi",
      "condition_unique_lda",
      "condition_unique_multivi",
      "comparison_aggr_weight_lda",
      "comparison_aggr_weight_vae_mlp",
      "comparison_aggr_lda",
      "comparison_aggr_multivi",
      "comparison_unique_lda",
      "comparison_unique_multivi"
    ),
    method_order = seq_len(12L),
    context_type = c(
      "condition", "condition", "condition", "condition", "condition", "condition",
      "comparison", "comparison", "comparison", "comparison", "comparison", "comparison"
    ),
    fp_mode = c(
      "aggregate_weight", "aggregate_weight", "aggregate", "aggregate", "unique", "unique",
      "aggregate_weight", "aggregate_weight", "aggregate", "aggregate", "unique", "unique"
    ),
    backend = c(
      "warplda", "vae", "warplda", "vae", "warplda", "vae",
      "warplda", "vae", "warplda", "vae", "warplda", "vae"
    ),
    vae_variant = c(
      "warplda", "vae_mlp", "warplda", "multivi_encoder", "warplda", "multivi_encoder",
      "warplda", "vae_mlp", "warplda", "multivi_encoder", "warplda", "multivi_encoder"
    )
  )
}

.m3tb_model_token <- function(backend, vae_variant) {
  if (identical(as.character(backend), "warplda")) return("model_wlda")
  switch(
    as.character(vae_variant),
    multivi_encoder = "model_mve",
    vae_mlp = "model_vmlp",
    paste0("model_", .m3tb_safe_label(vae_variant))
  )
}

.m3tb_model_name <- function(backend, vae_variant) {
  if (identical(as.character(backend), "warplda")) "warplda" else as.character(vae_variant)
}

.m3tb_model_label <- function(backend, vae_variant) {
  if (identical(as.character(backend), "warplda")) return("LDA")
  switch(
    as.character(vae_variant),
    multivi_encoder = "MultiVI",
    vae_mlp = "VAE-MLP",
    as.character(vae_variant)
  )
}

.m3tb_backend_slug <- function(backend, vae_variant) {
  if (identical(as.character(backend), "warplda")) return("lda")
  switch(
    as.character(vae_variant),
    multivi_encoder = "multivi",
    vae_mlp = "vae_mlp",
    .m3tb_safe_label(vae_variant)
  )
}

.m3tb_method_slug <- function(row) {
  paste(
    as.character(row$context_type[[1L]]),
    switch(
      as.character(row$fp_mode[[1L]]),
      unique = "uniq",
      aggregate = "aggr",
      aggregate_weight = "aggr_weight",
      .m3tb_safe_label(row$fp_mode[[1L]])
    ),
    .m3tb_backend_slug(row$backend[[1L]], row$vae_variant[[1L]]),
    sep = "_"
  )
}

.m3tb_setup_info <- function(context_type, fp_mode, backend, vae_variant) {
  mode_tag <- switch(
    as.character(fp_mode),
    unique = "uniq",
    aggregate = "aggr",
    aggregate_weight = "aggr_weight",
    .log_abort("Unsupported Module 3 FP mode: {fp_mode}")
  )
  mode_label <- switch(
    as.character(fp_mode),
    unique = "fp uniq",
    aggregate = "fp aggr",
    aggregate_weight = "fp aggr weight"
  )
  combo_id <- paste(
    "doc_tf",
    paste0("fp_", mode_tag),
    "tf_on",
    "count_pcl",
    .m3tb_model_token(backend, vae_variant),
    sep = "_"
  )
  if (identical(as.character(context_type), "condition")) {
    list(
      setup = paste0("std_tf_cond_fp_", mode_tag),
      setup_label = paste("cond", mode_label),
      doc_design = "condition",
      weight_label = "peak_score_gene_expr",
      combo_id = combo_id
    )
  } else {
    list(
      setup = paste0("std_tf_diff_fp_", mode_tag),
      setup_label = paste("diff", mode_label),
      doc_design = "comparison",
      weight_label = "peak_log2fc_fp_gene_fc_expr",
      combo_id = combo_id
    )
  }
}

.module3_topic_method_plan <- function(methods = "comparison_aggr_multivi",
                                       k_grid = 10L) {
  dict <- .m3tb_method_dictionary()
  if (identical(methods, "default")) methods <- "comparison_aggr_multivi"
  if (identical(methods, "all")) methods <- dict$method
  keep <- dict[method %in% methods]
  if (!nrow(keep)) {
    .log_abort("No valid Module 3 topic benchmark methods were selected.")
  }
  rows <- lapply(seq_len(nrow(keep)), function(i) {
    info <- .m3tb_setup_info(
      context_type = keep$context_type[[i]],
      fp_mode = keep$fp_mode[[i]],
      backend = keep$backend[[i]],
      vae_variant = keep$vae_variant[[i]]
    )
    data.table::data.table(
      method = keep$method[[i]],
      method_order = keep$method_order[[i]],
      context_type = keep$context_type[[i]],
      fp_mode = keep$fp_mode[[i]],
      backend = keep$backend[[i]],
      vae_variant = keep$vae_variant[[i]],
      setup = info$setup,
      setup_label = info$setup_label,
      doc_design = info$doc_design,
      weight_label = info$weight_label,
      combo_id = info$combo_id,
      model_label = .m3tb_model_label(keep$backend[[i]], keep$vae_variant[[i]]),
      model_dir = file.path(
        info$setup,
        "01_topic_models",
        info$combo_id,
        paste0("GSE192390_vae_joint_tf_docs_", info$weight_label, "_", .m3tb_model_name(keep$backend[[i]], keep$vae_variant[[i]]), "_Kgrid")
      )
    )
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, method_setup := paste(setup_label, model_label, sep = " | ")]
  out[, k_grid := paste(as.integer(k_grid), collapse = ",")]
  data.table::setorder(out, method_order)
  out[]
}

.m3tb_legacy_model_exists <- function(output_dir, method_plan, k_grid) {
  for (i in seq_len(nrow(method_plan))) {
    row <- method_plan[i]
    models_dir <- file.path(output_dir, row$model_dir[[1L]], "vae_models")
    if (!dir.exists(models_dir)) next
    if (any(file.exists(file.path(models_dir, sprintf("theta_K%d.csv", as.integer(k_grid)))))) {
      return(TRUE)
    }
  }
  FALSE
}

.m3tb_resolve_output_layout <- function(output_layout,
                                        output_dir,
                                        method_plan,
                                        k_grid,
                                        run_training,
                                        run_extraction) {
  output_layout <- match.arg(output_layout, c("auto", "standard", "benchmark", "legacy"))
  if (!identical(output_layout, "auto")) return(output_layout)
  if (!isTRUE(run_training) && !isTRUE(run_extraction) &&
      .m3tb_legacy_model_exists(output_dir, method_plan, k_grid)) {
    return("legacy")
  }
  if (nrow(method_plan) == 1L) "standard" else "benchmark"
}

.m3tb_apply_output_layout <- function(method_plan, output_dir, output_layout) {
  out <- data.table::copy(method_plan)
  if (identical(output_layout, "legacy")) {
    run_id <- sprintf("legacy_%03d", seq_len(nrow(out)))
    run_slug <- vapply(seq_len(nrow(out)), function(i) .m3tb_method_slug(out[i]), character(1L))
    out[, `:=`(
      run_id = run_id,
      run_slug = run_slug,
      run_dir = file.path(output_dir, out[["setup"]]),
      topic_documents_dir = file.path(output_dir, out[["setup"]], "topic_documents"),
      topic_models_dir = file.path(output_dir, out[["setup"]], "01_topic_models", out[["combo_id"]]),
      topic_extraction_dir = file.path(output_dir, out[["setup"]], "02_topic_extraction")
    )]
    return(out)
  }
  if (identical(output_layout, "standard")) {
    out[, `:=`(
      run_id = "selected",
      run_slug = .m3tb_method_slug(out[1]),
      run_dir = output_dir,
      topic_documents_dir = file.path(output_dir, "topic_documents"),
      topic_models_dir = file.path(output_dir, "topic_models"),
      topic_extraction_dir = file.path(output_dir, "topic_extraction")
    )]
    return(out)
  }
  run_slug <- vapply(seq_len(nrow(out)), function(i) .m3tb_method_slug(out[i]), character(1L))
  run_id <- sprintf("run_%03d", seq_len(nrow(out)))
  run_dir <- file.path(output_dir, paste(run_id, run_slug, sep = "_"))
  out[, `:=`(
    run_id = run_id,
    run_slug = run_slug,
    run_dir = run_dir,
    topic_documents_dir = file.path(run_dir, "topic_documents"),
    topic_models_dir = file.path(run_dir, "topic_models"),
    topic_extraction_dir = file.path(run_dir, "topic_extraction")
  )]
  out[]
}

.m3tb_extraction_output_dirs <- function(row, k_grid, output_layout) {
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[is.finite(k_grid)]
  if (!length(k_grid)) return(character())
  root <- row$topic_extraction_dir[[1L]]
  if (identical(output_layout, "standard") && length(k_grid) == 1L) {
    return(root)
  }
  file.path(root, paste0("K", k_grid))
}

.m3tb_flat_output_for_layout <- function(output_layout) {
  !identical(as.character(output_layout)[[1L]], "legacy")
}

.m3tb_review_dir <- function(output_dir, output_layout) {
  if (identical(output_layout, "legacy")) {
    file.path(output_dir, "review_topic_experiments")
  } else {
    file.path(output_dir, "review")
  }
}

.m3tb_check_path_lengths <- function(root,
                                     max_chars = 240L,
                                     path_alias = NULL,
                                     out_file = NULL,
                                     verbose = TRUE,
                                     top_n = 25L) {
  .assert_pkg("data.table")
  if (is.null(root) || !nzchar(as.character(root)[[1L]]) || !dir.exists(root)) {
    return(list(n_paths = 0L, n_over_limit = 0L, max_length = 0L, report = data.table::data.table()))
  }
  max_chars <- suppressWarnings(as.integer(max_chars[[1L]]))
  if (!is.finite(max_chars) || max_chars < 1L) max_chars <- 240L
  top_n <- suppressWarnings(as.integer(top_n[[1L]]))
  if (!is.finite(top_n) || top_n < 1L) top_n <- 25L
  if (is.null(path_alias) || !nzchar(as.character(path_alias)[[1L]])) {
    env_alias <- Sys.getenv("CRAFTGRN_PATH_LENGTH_ALIAS", unset = "")
    path_alias <- if (nzchar(env_alias)) env_alias else NULL
  }

  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  paths <- c(
    root_norm,
    list.files(root_norm, all.files = TRUE, no.. = TRUE, recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
  )
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  display_path <- paths
  if (!is.null(path_alias) && nzchar(as.character(path_alias)[[1L]])) {
    alias <- sub("/+$", "", gsub("\\\\", "/", as.character(path_alias)[[1L]]))
    rel <- substr(paths, nchar(root_norm) + 1L, nchar(paths))
    rel[rel == ""] <- ""
    display_path <- paste0(alias, rel)
  }

  report <- data.table::data.table(
    path = paths,
    display_path = display_path,
    n_char = nchar(display_path, type = "chars", allowNA = FALSE, keepNA = FALSE)
  )
  report[, over_limit := n_char > max_chars]
  data.table::setorder(report, -n_char, display_path)
  over <- report[over_limit == TRUE]
  if (!is.null(out_file) && nzchar(as.character(out_file)[[1L]])) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(over, out_file)
  }
  n_over <- nrow(over)
  max_len <- if (nrow(report)) max(report$n_char, na.rm = TRUE) else 0L
  if (isTRUE(verbose) && n_over > 0L) {
    shown <- head(over$display_path, top_n)
    .log_warn(c(
      "Module 3 path length sanity check found {n_over} path(s) longer than {max_chars} characters.",
      i = "Long paths may fail to open through Windows/Samba clients.",
      i = "Longest path length: {max_len} characters.",
      i = "Examples:",
      paste0("  - ", shown)
    ))
  }
  list(
    n_paths = nrow(report),
    n_over_limit = n_over,
    max_length = max_len,
    report = over
  )
}

.m3tb_clean_stale_review_layout <- function(review_dir) {
  stale <- file.path(review_dir, c("csv", "html"))
  stale <- stale[dir.exists(stale)]
  if (length(stale)) {
    unlink(stale, recursive = TRUE, force = TRUE)
  }
  invisible(TRUE)
}

.m3tb_read_theta <- function(path) {
  dt <- data.table::fread(path, showProgress = FALSE)
  if (!"doc_id" %in% names(dt)) data.table::setnames(dt, names(dt)[[1L]], "doc_id")
  doc_id <- as.character(dt$doc_id)
  topic_cols <- setdiff(names(dt), "doc_id")
  mat <- as.matrix(dt[, topic_cols, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- doc_id
  colnames(mat) <- topic_cols
  mat
}

.m3tb_normalize_rows <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x) | x < 0] <- 0
  rs <- rowSums(x)
  rs[!is.finite(rs) | rs <= 0] <- 1
  x / rs
}

.m3tb_jsd_matrix <- function(mat) {
  p <- .m3tb_normalize_rows(mat)
  n <- nrow(p)
  out <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      m <- 0.5 * (p[i, ] + p[j, ])
      kl_i <- sum(ifelse(p[i, ] > 0, p[i, ] * log(p[i, ] / pmax(m, 1e-300)), 0))
      kl_j <- sum(ifelse(p[j, ] > 0, p[j, ] * log(p[j, ] / pmax(m, 1e-300)), 0))
      val <- sqrt(max(0, 0.5 * (kl_i + kl_j)))
      out[i, j] <- val
      out[j, i] <- val
    }
  }
  rownames(out) <- rownames(p)
  colnames(out) <- rownames(p)
  out
}

.m3tb_parse_doc_label <- function(doc_id, context_type) {
  parts <- strsplit(as.character(doc_id), "::", fixed = TRUE)
  if (identical(context_type, "condition")) {
    return(vapply(parts, function(x) x[[1L]], character(1L)))
  }
  vapply(parts, function(x) {
    if (length(x) < 2L) return(x[[1L]])
    paste(x[[1L]], x[[length(x)]], sep = "::")
  }, character(1L))
}

.m3tb_display_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_vs_", " vs ", x, fixed = TRUE)
  x <- gsub("::Target-Up", " Up", x, fixed = TRUE)
  x <- gsub("::Target-Down", " Down", x, fixed = TRUE)
  x <- gsub("::Up", " Up", x, fixed = TRUE)
  x <- gsub("::Down", " Down", x, fixed = TRUE)
  x <- gsub("_", " ", x, fixed = TRUE)
  trimws(x)
}

.m3tb_short_label <- function(x, max_chars = 18L) {
  x <- .m3tb_display_label(x)
  x <- gsub("[()]", "", x)
  x <- gsub("[[:space:]]+", " ", x)
  max_chars <- as.integer(max_chars[[1L]])
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1L, max_chars - 1L), "."), x)
}

.m3tb_color_family <- function(x) {
  out <- as.character(x)
  out[is.na(out) | !nzchar(trimws(out))] <- "Other"
  out
}

.m3tb_group_color <- function(x) {
  fam <- .m3tb_color_family(x)
  pal <- c(
    "#4E79A7", "#E15759", "#59A14F", "#F28E2B", "#B07AA1", "#76B7B2",
    "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC", "#1B9E77", "#D95F02",
    "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
    "#8DD3C7", "#FB8072", "#80B1D3", "#B3DE69", "#FDB462", "#BC80BD"
  )
  keys <- sort(unique(fam))
  map <- pal[((seq_along(keys) - 1L) %% length(pal)) + 1L]
  unname(map[match(fam, keys)])
}

.m3tb_metric_group_from_label <- function(label, direction = NULL) {
  out <- as.character(label)
  out[is.na(out) | !nzchar(trimws(out))] <- "Other"
  if (!is.null(direction)) {
    direction <- as.character(direction)
    out <- paste(out, direction, sep = "::")
  }
  out
}

.m3tb_use_inferred_metric_groups <- function(metric_group) {
  metric_group <- as.character(metric_group)
  any(!is.na(metric_group) & nzchar(metric_group) & duplicated(metric_group))
}

.m3tb_condition_base <- function(x) {
  sub("_[0-9]+$", "", as.character(x), perl = TRUE)
}

.m3tb_replicate_columns <- function(multiomic_data) {
  cols <- colnames(multiomic_data$matrices$fp_score)
  if (is.null(cols) || !length(cols)) {
    .log_abort("replicate_documents requires multiomic_data$matrices$fp_score with condition columns.")
  }
  out <- data.table::data.table(
    condition_label = as.character(cols),
    condition_base = .m3tb_condition_base(cols)
  )
  out[out[["condition_label"]] != out[["condition_base"]]]
}

.m3tb_replicate_design_table <- function(comparisons, multiomic_data) {
  dt <- data.table::as.data.table(comparisons)
  if (!all(c("cond1_label", "cond2_label") %in% names(dt))) {
    .log_abort("replicate_documents requires comparisons with cond1_label and cond2_label columns.")
  }
  comp <- data.table::copy(dt)
  comp[, `:=`(
    cond1_label = as.character(cond1_label),
    cond2_label = as.character(cond2_label)
  )]
  if (!"comparison_id" %in% names(comp)) {
    comp[, comparison_id := paste(cond1_label, cond2_label, sep = "_vs_")]
  } else {
    comp[, comparison_id := as.character(comparison_id)]
  }
  rep_map <- .m3tb_replicate_columns(multiomic_data)
  rep_map <- rep_map[rep_map[["condition_base"]] %in% unique(c(comp$cond1_label, comp$cond2_label))]
  if (!nrow(rep_map)) {
    .log_abort("No replicate condition columns matched the requested comparisons.")
  }
  condition <- unique(data.table::data.table(
    context_type = "condition",
    comparison_label = rep_map$condition_label,
    display_label = .m3tb_display_label(rep_map$condition_label),
    metric_group = rep_map$condition_base
  ))
  case_rep <- merge(
    comp,
    data.table::data.table(
      cond1_label = rep_map$condition_base,
      direction_condition_label = rep_map$condition_label
    ),
    by = "cond1_label",
    allow.cartesian = TRUE
  )
  case_rep[, direction := "Target-Up"]
  ctrl_rep <- merge(
    comp,
    data.table::data.table(
      cond2_label = rep_map$condition_base,
      direction_condition_label = rep_map$condition_label
    ),
    by = "cond2_label",
    allow.cartesian = TRUE
  )
  ctrl_rep[, direction := "Target-Down"]
  comparison <- data.table::rbindlist(list(case_rep, ctrl_rep), use.names = TRUE, fill = TRUE)
  comparison_label <- paste0(comparison$comparison_id, "__", comparison$direction_condition_label, "::", comparison$direction)
  comparison_metric_base <- if ("comparison_group" %in% names(comparison)) {
    as.character(comparison$comparison_group)
  } else {
    as.character(comparison$comparison_id)
  }
  comparison_metric_base[is.na(comparison_metric_base) | !nzchar(trimws(comparison_metric_base))] <- comparison$comparison_id[
    is.na(comparison_metric_base) | !nzchar(trimws(comparison_metric_base))
  ]
  comparison[, `:=`(
    context_type = "comparison",
    comparison_label = comparison_label,
    display_label = .m3tb_display_label(comparison_label),
    metric_group = paste(comparison_metric_base, comparison$direction, sep = "::")
  )]
  comparison <- unique(comparison[, .(context_type, comparison_label, display_label, metric_group)])
  data.table::rbindlist(list(condition, comparison), use.names = TRUE)
}

.m3tb_design_table <- function(comparisons,
                               replicate_documents = FALSE,
                               multiomic_data = NULL) {
  if (is.null(comparisons)) {
    return(data.table::data.table())
  }
  if (is.character(comparisons) && length(comparisons) == 1L) {
    if (!file.exists(comparisons)) .log_abort("Comparison file does not exist: {comparisons}")
    comparisons <- data.table::fread(comparisons, showProgress = FALSE)
  }
  dt <- data.table::as.data.table(comparisons)
  if (all(c("context_type", "comparison_label", "metric_group") %in% names(dt))) {
    out <- unique(dt[, .(
      context_type = as.character(context_type),
      comparison_label = as.character(comparison_label),
      display_label = if ("display_label" %in% names(dt)) {
        as.character(display_label)
      } else {
        .m3tb_display_label(comparison_label)
      },
      metric_group = as.character(metric_group)
    )])
    return(out)
  }
  if (isTRUE(replicate_documents)) {
    if (is.null(multiomic_data)) {
      .log_abort("replicate_documents = TRUE requires multiomic_data.")
    }
    return(.m3tb_replicate_design_table(dt, multiomic_data))
  }
  if (all(c("condition_label", "condition_group") %in% names(dt))) {
    return(unique(dt[, .(
      context_type = "condition",
      comparison_label = as.character(condition_label),
      display_label = .m3tb_display_label(condition_label),
      metric_group = as.character(condition_group)
    )]))
  }
  if (all(c("cond1_label", "cond2_label") %in% names(dt))) {
    conds <- unique(c(as.character(dt$cond1_label), as.character(dt$cond2_label)))
    condition <- data.table::data.table(
      context_type = "condition",
      comparison_label = conds,
      display_label = .m3tb_display_label(conds),
      metric_group = conds
    )
    if (!"comparison_id" %in% names(dt)) {
      dt[, comparison_id := paste(as.character(cond1_label), as.character(cond2_label), sep = "_vs_")]
    } else {
      dt[, comparison_id := as.character(comparison_id)]
    }
    if ("comparison_label" %in% names(dt)) {
      display_base <- as.character(dt$comparison_label)
    } else if ("comparison_display" %in% names(dt)) {
      display_base <- as.character(dt$comparison_display)
    } else {
      display_base <- paste(as.character(dt$cond1_label), as.character(dt$cond2_label), sep = " vs ")
    }
    display_base[is.na(display_base) | !nzchar(trimws(display_base))] <- dt$comparison_id[
      is.na(display_base) | !nzchar(trimws(display_base))
    ]
    metric_base <- if ("comparison_group" %in% names(dt)) {
      as.character(dt$comparison_group)
    } else if ("cond1_display" %in% names(dt)) {
      as.character(dt$cond1_display)
    } else {
      as.character(dt$cond1_label)
    }
    metric_base[is.na(metric_base) | !nzchar(trimws(metric_base))] <- display_base[
      is.na(metric_base) | !nzchar(trimws(metric_base))
    ]
    up_metric_group <- .m3tb_metric_group_from_label(metric_base, "Target-Up")
    down_metric_group <- .m3tb_metric_group_from_label(metric_base, "Target-Down")
    if (!"comparison_group" %in% names(dt) && !.m3tb_use_inferred_metric_groups(c(up_metric_group, down_metric_group))) {
      up_metric_group <- paste(dt$comparison_id, "Target-Up", sep = "::")
      down_metric_group <- paste(dt$comparison_id, "Target-Down", sep = "::")
    }
    comparison <- data.table::rbindlist(list(
      data.table::data.table(
        context_type = "comparison",
        comparison_label = paste(dt$comparison_id, "Target-Up", sep = "::"),
        display_label = paste(display_base, "Target-Up"),
        metric_group = up_metric_group
      ),
      data.table::data.table(
        context_type = "comparison",
        comparison_label = paste(dt$comparison_id, "Target-Down", sep = "::"),
        display_label = paste(display_base, "Target-Down"),
        metric_group = down_metric_group
      )
    ), use.names = TRUE, fill = TRUE)
    return(unique(data.table::rbindlist(list(condition, comparison), use.names = TRUE, fill = TRUE)))
  }
  .log_abort("comparisons must include either condition_label/condition_group or cond1_label/cond2_label columns.")
}

.m3tb_topic_profiles <- function(theta, context_type) {
  topic_cols <- colnames(theta)
  labels <- .m3tb_parse_doc_label(rownames(theta), context_type)
  dt <- data.table::as.data.table(theta)
  dt[, doc_id := rownames(theta)]
  dt[, comparison_label := labels]
  dt[, display_label := .m3tb_display_label(comparison_label)]
  out <- dt[, c(
    list(n_docs = data.table::uniqueN(doc_id)),
    lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE))
  ), by = .(comparison_label, display_label), .SDcols = topic_cols]
  data.table::setorder(out, comparison_label)
  out[]
}

.m3tb_score_theta_one <- function(theta, row, design, csv_dir) {
  profiles <- .m3tb_topic_profiles(theta, row$context_type[[1L]])
  topic_cols <- grep("^Topic[0-9]+$", names(profiles), value = TRUE)
  if (!length(topic_cols)) topic_cols <- setdiff(names(profiles), c("comparison_label", "display_label", "n_docs"))
  merged <- merge(
    profiles,
    design[context_type == row$context_type[[1L]], .(
      comparison_label,
      design_display_label = display_label,
      metric_group
    )],
    by = "comparison_label",
    all.x = TRUE,
    sort = FALSE
  )
  if ("design_display_label" %in% names(merged)) {
    use_design_label <- !is.na(merged$design_display_label) & nzchar(trimws(merged$design_display_label))
    merged[use_design_label, display_label := design_display_label]
    merged[, design_display_label := NULL]
  }
  missing <- merged[is.na(metric_group) | !nzchar(metric_group)]
  if (nrow(missing)) {
    merged[is.na(metric_group) | !nzchar(metric_group), metric_group := comparison_label]
  }
  x <- as.matrix(merged[, topic_cols, with = FALSE])
  rownames(x) <- merged$comparison_label
  dmat <- .m3tb_jsd_matrix(x)
  jsd_path <- file.path(
    csv_dir,
    paste0(
      "theta_",
      row$context_type[[1L]],
      "_jsd_distance_matrix_",
      .m3tb_safe_label(row$setup_label[[1L]]),
      "_K",
      as.integer(row$selected_k[[1L]]),
      "_",
      .m3tb_topic_space_value(row, "topic_space", "raw"),
      ".csv"
    )
  )
  data.table::fwrite(data.table::data.table(label = rownames(dmat), as.data.frame(dmat, check.names = FALSE)), jsd_path)
  coords <- if (nrow(dmat) <= 1L) {
    matrix(c(0, 0), ncol = 2L)
  } else {
    fit <- tryCatch(stats::cmdscale(stats::as.dist(dmat), k = 2L, eig = FALSE), error = function(e) NULL)
    if (is.null(fit)) {
      cbind(seq_len(nrow(dmat)), rep(0, nrow(dmat)))
    } else if (is.null(dim(fit))) {
      cbind(fit, rep(0, length(fit)))
    } else if (ncol(fit) < 2L) {
      cbind(fit[, 1L], rep(0, nrow(fit)))
    } else {
      fit[, 1:2, drop = FALSE]
    }
  }
  groups <- merged$metric_group
  per_label <- lapply(seq_len(nrow(dmat)), function(i) {
    same <- setdiff(which(groups == groups[[i]]), i)
    other <- which(groups != groups[[i]])
    inner <- if (length(same)) mean(dmat[i, same], na.rm = TRUE) else NA_real_
    inter <- if (length(other)) mean(dmat[i, other], na.rm = TRUE) else NA_real_
    score <- if (length(same) && length(other) && is.finite(inner) && is.finite(inter)) {
      (inter - inner) / (inter + inner + 1e-12)
    } else {
      NA_real_
    }
    data.table::data.table(
      method_order = as.integer(row$method_order[[1L]]),
      context_type = row$context_type[[1L]],
      setup = row$setup[[1L]],
      setup_label = row$setup_label[[1L]],
      model_label = row$model_label[[1L]],
      method_setup = row$method_setup[[1L]],
      k = as.integer(row$selected_k[[1L]]),
      display_label = merged$display_label[[i]],
      comparison_label = merged$comparison_label[[i]],
      metric_group = groups[[i]],
      group_size = sum(groups == groups[[i]]),
      n_docs = merged$n_docs[[i]],
      inner_distance = inner,
      inter_distance = inter,
      theta_condition_label_score = score
    )
  })
  per_label <- data.table::rbindlist(per_label, use.names = TRUE, fill = TRUE)
  score <- per_label[, .(
    n_labels = .N,
    n_scored_labels = sum(is.finite(theta_condition_label_score)),
    n_metric_groups = data.table::uniqueN(metric_group),
    mean_inner_distance = mean(inner_distance, na.rm = TRUE),
    mean_inter_distance = mean(inter_distance, na.rm = TRUE),
    theta_condition_separation_score = mean(theta_condition_label_score, na.rm = TRUE)
  ), by = .(method_order, context_type, setup, setup_label, model_label, method_setup, k)]
  mds_points <- data.table::copy(merged)
  mds_points[, `:=`(
    group_label = comparison_label,
    MDS1 = as.numeric(coords[, 1L]),
    MDS2 = as.numeric(coords[, 2L]),
    n_features = nrow(theta),
    method_order = as.integer(row$method_order[[1L]]),
    method_setup = row$method_setup[[1L]],
    k = as.integer(row$selected_k[[1L]]),
    setup = row$setup[[1L]],
    setup_label = row$setup_label[[1L]],
    doc_design = row$context_type[[1L]],
    fp_term_mode = row$fp_mode[[1L]],
    model_label = row$model_label[[1L]],
    panel_label = paste(row$setup_label[[1L]], row$model_label[[1L]], sep = "\n"),
    color_family = .m3tb_color_family(metric_group),
    color = .m3tb_group_color(metric_group)
  )]
  list(score = score, per_label = per_label, profiles = profiles, mds_points = mds_points)
}

.m3tb_candidate_model_dirs <- function(output_dir, row) {
  combo_root <- file.path(output_dir, row$setup[[1L]], "01_topic_models", row$combo_id[[1L]])
  expected <- file.path(output_dir, row$model_dir[[1L]])
  candidates <- character()
  if ("topic_models_dir" %in% names(row) && dir.exists(row$topic_models_dir[[1L]])) {
    found <- list.dirs(row$topic_models_dir[[1L]], recursive = TRUE, full.names = TRUE)
    found <- dirname(found[basename(found) == "vae_models"])
    candidates <- c(candidates, found)
  }
  if (dir.exists(file.path(expected, "vae_models"))) {
    candidates <- c(candidates, expected)
  }
  if (dir.exists(combo_root)) {
    found <- list.dirs(combo_root, recursive = TRUE, full.names = TRUE)
    found <- dirname(found[basename(found) == "vae_models"])
    candidates <- c(candidates, found)
  }
  unique(candidates[file.exists(file.path(candidates, "vae_models"))])
}

.m3tb_existing_model_rows <- function(output_dir, method_plan, k_grid) {
  rows <- lapply(seq_len(nrow(method_plan)), function(i) {
    row <- method_plan[i]
    candidates <- .m3tb_candidate_model_dirs(output_dir, row)
    if (!length(candidates)) return(NULL)
    found <- lapply(candidates, function(model_dir_abs) {
      models_dir <- file.path(model_dir_abs, "vae_models")
      present <- as.integer(k_grid[
        file.exists(file.path(models_dir, sprintf("theta_K%d.csv", as.integer(k_grid)))) &
          file.exists(file.path(models_dir, sprintf("phi_K%d.csv", as.integer(k_grid))))
      ])
      if (!length(present)) return(NULL)
      out <- data.table::copy(row)
      out <- out[rep(1L, length(present))]
      out[, `:=`(selected_k = present, model_dir = model_dir_abs)]
      out
    })
    data.table::rbindlist(found, use.names = TRUE, fill = TRUE)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) .log_abort("No existing theta/phi files were found under {output_dir}.")
  data.table::setorder(out, method_order, selected_k)
  out[]
}

.m3tb_topic_space_rows <- function(model_rows) {
  model_rows <- data.table::as.data.table(model_rows)
  if (!nrow(model_rows)) return(model_rows)
  rows <- lapply(seq_len(nrow(model_rows)), function(i) {
    row <- data.table::copy(model_rows[i])
    trained_k <- as.integer(row$selected_k[[1L]])
    extraction_dir <- .m3tb_find_extraction_subdir(row)
    map_file <- file.path(extraction_dir, "topic_optimization_map.csv")
    optimized_theta_file <- file.path(
      extraction_dir,
      "rds",
      "topic_theta_optimized.rds"
    )
    optimized_terms_file <- file.path(extraction_dir, "topic_terms.csv")
    has_combined <- file.exists(map_file) &&
      file.exists(optimized_theta_file) &&
      file.exists(optimized_terms_file)
    optimized_k <- trained_k
    if (has_combined) {
      mapping <- data.table::fread(
        map_file,
        select = "optimized_topic",
        showProgress = FALSE
      )
      optimized_k <- data.table::uniqueN(mapping$optimized_topic)
    }
    if (!is.finite(optimized_k) || optimized_k < 1L) optimized_k <- trained_k
    raw <- data.table::copy(row)
    raw[, `:=`(
      trained_k = trained_k,
      topic_space = "raw",
      topic_count = trained_k,
      report_key = sprintf("raw_k%d", trained_k),
      report_label = sprintf("Raw K%d", trained_k)
    )]
    if (!has_combined) return(raw)
    combined <- data.table::copy(row)
    combined[, `:=`(
      trained_k = trained_k,
      topic_space = "combined",
      topic_count = optimized_k,
      report_key = sprintf("combined_k%d_from_k%d", optimized_k, trained_k),
      report_label = sprintf("Combined K%d (from K%d)", optimized_k, trained_k)
    )]
    data.table::rbindlist(list(raw, combined), use.names = TRUE, fill = TRUE)
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  out[, topic_space_order := data.table::fifelse(topic_space == "raw", 1L, 2L)]
  data.table::setorder(out, method_order, trained_k, topic_space_order)
  out[, topic_space_order := NULL]
  out[]
}

.m3tb_topic_space_value <- function(row, name, default = NULL) {
  if (!(name %in% names(row))) return(default)
  value <- row[[name]][[1L]]
  if (length(value) && !is.na(value)) value else default
}

.m3tb_topic_space_theta <- function(row) {
  topic_space <- .m3tb_topic_space_value(row, "topic_space", "raw")
  trained_k <- as.integer(.m3tb_topic_space_value(
    row,
    "trained_k",
    row$selected_k[[1L]]
  ))
  if (identical(topic_space, "combined")) {
    extraction_dir <- .m3tb_find_extraction_subdir(row)
    path <- file.path(extraction_dir, "rds", "topic_theta_optimized.rds")
    if (!file.exists(path)) {
      .log_abort(
        "Missing combined theta for trained K{trained_k}: {path}"
      )
    }
    theta <- readRDS(path)
    theta <- as.matrix(theta)
    storage.mode(theta) <- "numeric"
    return(theta)
  }
  path <- file.path(
    row$model_dir[[1L]],
    "vae_models",
    sprintf("theta_K%d.csv", trained_k)
  )
  .m3tb_read_theta(path)
}

.m3tb_topic_space_file <- function(extraction_dir,
                                    topic_space = c("raw", "combined"),
                                    kind = c("terms", "overall_pathways",
                                             "condition_pathways")) {
  topic_space <- match.arg(topic_space)
  kind <- match.arg(kind)
  names_by_kind <- list(
    terms = if (identical(topic_space, "raw")) {
      c("topic_terms_raw.csv")
    } else {
      c("topic_terms.csv")
    },
    overall_pathways = if (identical(topic_space, "raw")) {
      c("topic_pathway_enrichment_topic_terms_raw.csv")
    } else {
      c(
        "topic_pathway_enrichment_topic_terms_combined.csv",
        "topic_pathway_enrichment_topic_terms.csv"
      )
    },
    condition_pathways = if (identical(topic_space, "raw")) {
      c("per_condition_topic_pathway_enrichment_raw.csv")
    } else {
      c(
        "per_condition_topic_pathway_enrichment_combined.csv",
        "per_condition_topic_pathway_enrichment.csv"
      )
    }
  )
  candidates <- file.path(extraction_dir, names_by_kind[[kind]])
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) hit[[1L]] else candidates[[1L]]
}

.m3tb_score_theta_outputs <- function(output_dir,
                                      model_rows,
                                      comparisons,
                                      score_prefix = "theta_condition_separation",
                                      replicate_documents = FALSE,
                                      multiomic_data = NULL,
                                      review_dir = NULL) {
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  csv_dir <- .m3tb_review_tables_dir(review_dir)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  design <- .m3tb_design_table(
    comparisons,
    replicate_documents = replicate_documents,
    multiomic_data = multiomic_data
  )
  data.table::fwrite(design, file.path(csv_dir, paste0(score_prefix, "_design.csv")))
  payload <- lapply(seq_len(nrow(model_rows)), function(i) {
    row <- model_rows[i]
    theta <- .m3tb_topic_space_theta(row)
    scored <- .m3tb_score_theta_one(theta, row, design, csv_dir)
    meta <- c("trained_k", "topic_space", "topic_count", "report_key", "report_label")
    meta <- meta[meta %in% names(row)]
    for (part in names(scored)) {
      if (!is.data.frame(scored[[part]])) next
      for (column in meta) scored[[part]][, (column) := row[[column]][[1L]]]
    }
    scored
  })
  scores <- data.table::rbindlist(lapply(payload, `[[`, "score"), use.names = TRUE, fill = TRUE)
  per_label <- data.table::rbindlist(lapply(payload, `[[`, "per_label"), use.names = TRUE, fill = TRUE)
  mds_points <- data.table::rbindlist(lapply(payload, `[[`, "mds_points"), use.names = TRUE, fill = TRUE)
  data.table::setorder(scores, method_order, k)
  data.table::setorder(per_label, method_order, k, comparison_label)
  if (nrow(mds_points)) data.table::setorder(mds_points, k, method_order, group_label)
  data.table::fwrite(scores, file.path(csv_dir, paste0(score_prefix, "_scores.csv")))
  data.table::fwrite(per_label, file.path(csv_dir, paste0(score_prefix, "_per_label.csv")))
  data.table::fwrite(per_label, file.path(csv_dir, paste0(score_prefix, "_score_long.csv")))
  data.table::fwrite(mds_points, file.path(csv_dir, "theta_group_mds_points.csv"))
  heatmap_columns <- intersect(
    c(
      "method_order", "context_type", "method_setup", "setup",
      "setup_label", "model_label", "k", "trained_k", "topic_space",
      "topic_count", "report_key", "report_label",
      "theta_condition_separation_score"
    ),
    names(scores)
  )
  heatmap_values <- scores[, ..heatmap_columns]
  data.table::fwrite(heatmap_values, file.path(csv_dir, paste0(score_prefix, "_score_heatmap_values.csv")))
  score_wide <- data.table::copy(scores)
  if (!"topic_space" %in% names(score_wide)) score_wide[, topic_space := "raw"]
  wide <- data.table::dcast(
    score_wide[, .(
      method_order, context_type, method_setup, topic_space,
      k = as.integer(k), theta_condition_separation_score
    )],
    method_order + context_type + method_setup + topic_space ~ k,
    value.var = "theta_condition_separation_score"
  )
  k_cols <- as.character(sort(unique(scores$k)))
  data.table::setnames(wide, intersect(k_cols, names(wide)), paste0("K", intersect(k_cols, names(wide))))
  data.table::setorder(wide, method_order)
  data.table::fwrite(wide, file.path(csv_dir, paste0(score_prefix, "_score_heatmap_values_matrix.csv")))
  list(scores = scores, per_label = per_label, matrix = wide, design = design, mds_points = mds_points, score_prefix = score_prefix)
}

.m3tb_find_topic_links <- function(output_dir, row) {
  roots <- character()
  if ("topic_extraction_dir" %in% names(row)) {
    roots <- c(roots, row$topic_extraction_dir[[1L]])
  }
  roots <- c(roots, file.path(output_dir, row$setup[[1L]], "02_topic_extraction"))
  roots <- unique(roots[dir.exists(roots)])
  if (!length(roots)) return(character())
  files <- unlist(lapply(roots, function(root) {
    list.files(root, "topic_links(_pass)?[.]csv([.]gz)?$", recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  if (!length(files) || !"selected_k" %in% names(row)) return(files)
  k <- as.integer(row$selected_k[[1L]])
  k_pattern <- paste0("(^|[^0-9])K", k, "([^0-9]|$)")
  keep <- grepl(k_pattern, files, perl = TRUE)
  if (any(keep)) return(files[keep])
  direct <- files[normalizePath(dirname(files), winslash = "/", mustWork = FALSE) %in%
    normalizePath(roots, winslash = "/", mustWork = FALSE)]
  direct
}

.m3tb_empty_pass_counts <- function() {
  data.table::data.table(
    method_order = integer(),
    method_setup = character(),
    setup = character(),
    model_label = character(),
    selected_k = integer(),
    status = character(),
    count = integer(),
    count_basis = character()
  )
}

.m3tb_empty_item_coverage_counts <- function() {
  data.table::data.table(
    method_order = integer(),
    method_setup = character(),
    setup = character(),
    model_label = character(),
    selected_k = integer(),
    unit = character(),
    status = character(),
    count = numeric(),
    total = numeric(),
    fraction = numeric(),
    percent = numeric(),
    count_basis = character()
  )
}

.m3tb_empty_shared_counts <- function() {
  data.table::data.table(
    n_topics = integer(),
    n_items = integer(),
    method_order = integer(),
    method_setup = character(),
    setup = character(),
    model_label = character(),
    selected_k = integer(),
    unit = character()
  )
}

.m3tb_as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(is.finite(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

.m3tb_count_shared_topics <- function(dt, item_col) {
  if (!nrow(dt) || !(item_col %in% names(dt)) || !("topic_num" %in% names(dt))) {
    return(data.table::data.table(n_topics = integer(), n_items = integer()))
  }
  tmp <- unique(dt[!is.na(get(item_col)) & get(item_col) != "", .(item = get(item_col), topic_num)])
  if (!nrow(tmp)) return(data.table::data.table(n_topics = integer(), n_items = integer()))
  tmp[, .(n_topics = data.table::uniqueN(topic_num)), by = item][, .(n_items = .N), by = n_topics]
}

.m3tb_read_topic_link_summary <- function(topic_dir) {
  path <- file.path(topic_dir, "topic_link_summary.csv")
  if (!file.exists(path)) return(NULL)
  dt <- data.table::fread(path, showProgress = FALSE)
  if (!nrow(dt) || !all(c("n_scored_rows", "n_pass_rows") %in% names(dt))) {
    return(NULL)
  }
  list(
    scored = sum(as.numeric(dt$n_scored_rows), na.rm = TRUE),
    pass = sum(as.numeric(dt$n_pass_rows), na.rm = TRUE),
    output_mode = paste(unique(as.character(dt$output_mode %||% NA_character_)), collapse = ",")
  )
}

.m3tb_attach_topic_row_meta <- function(dt, row) {
  if (!nrow(dt)) return(dt)
  dt[, `:=`(
    method_order = row$method_order[[1L]],
    method_setup = row$method_setup[[1L]],
    setup = row$setup[[1L]],
    model_label = row$model_label[[1L]],
    selected_k = as.integer(row$selected_k[[1L]])
  )]
  dt
}

.m3tb_read_topic_item_coverage_counts <- function(topic_dir, row) {
  path <- file.path(topic_dir, "topic_item_coverage_counts.csv")
  if (!file.exists(path)) return(data.table::data.table())
  dt <- data.table::fread(path, showProgress = FALSE)
  need <- c("unit", "status", "count")
  if (!all(need %in% names(dt))) return(data.table::data.table())
  if (!"total" %in% names(dt)) {
    totals <- dt[, .(total = sum(as.numeric(count), na.rm = TRUE)), by = unit]
    dt <- merge(dt, totals, by = "unit", all.x = TRUE, sort = FALSE)
  }
  if (!"fraction" %in% names(dt)) {
    dt[, fraction := data.table::fifelse(total > 0, as.numeric(count) / total, NA_real_)]
  }
  if (!"percent" %in% names(dt)) {
    dt[, percent := 100 * as.numeric(fraction)]
  }
  if (!"count_basis" %in% names(dt)) {
    dt[, count_basis := paste0(unit, " assigned to at least one topic")]
  }
  .m3tb_attach_topic_row_meta(dt, row)
}

.m3tb_read_topic_term_coverage_counts <- function(topic_dir, row) {
  path <- file.path(topic_dir, "topic_terms.csv")
  if (!file.exists(path)) return(data.table::data.table())
  tt <- data.table::fread(path, showProgress = FALSE)
  if (!nrow(tt) || !"term_id" %in% names(tt)) return(data.table::data.table())
  coverage <- .topic_item_coverage_from_terms(tt, score_mat = NULL)
  .m3tb_attach_topic_row_meta(coverage, row)
}

.m3tb_read_pathway_shared_counts <- function(topic_dir) {
  preferred <- file.path(
    topic_dir,
    c(
      "topic_pathway_enrichment_topic_terms.csv",
      "topic_pathway_enrichment_peak_and_gene_dotplot.csv",
      "topic_pathway_enrichment_gene_only_dotplot.csv",
      "topic_pathway_enrichment_peak_and_gene_prob_dotplot.csv"
    )
  )
  files <- preferred[file.exists(preferred)]
  if (!length(files)) {
    files <- list.files(
      topic_dir,
      pattern = "^topic_pathway_enrichment_.*_dotplot[.]csv$",
      full.names = TRUE
    )
  }
  if (!length(files)) return(data.table::data.table(n_topics = integer(), n_items = integer()))
  dt <- data.table::fread(files[[1L]], showProgress = FALSE)
  if (!nrow(dt)) return(data.table::data.table(n_topics = integer(), n_items = integer()))
  if (!"topic_num" %in% names(dt) && "topic" %in% names(dt)) {
    data.table::setnames(dt, "topic", "topic_num")
  }
  item_col <- if ("pathway_key" %in% names(dt)) {
    "pathway_key"
  } else if ("pathway_norm_key" %in% names(dt)) {
    "pathway_norm_key"
  } else if ("pathway" %in% names(dt)) {
    "pathway"
  } else {
    NA_character_
  }
  if (is.na(item_col) || !"topic_num" %in% names(dt)) {
    return(data.table::data.table(n_topics = integer(), n_items = integer()))
  }
  if ("is_sig" %in% names(dt)) {
    dt <- dt[.m3tb_as_flag(is_sig)]
  } else if ("padj" %in% names(dt)) {
    dt <- dt[is.finite(padj) & padj < 0.05]
  }
  .m3tb_count_shared_topics(dt, item_col)
}

.m3tb_existing_shared_cache <- function(review_dir) {
  if (is.null(review_dir)) return(NULL)
  path <- file.path(.m3tb_review_tables_dir(review_dir), "topic_setup_shared_topic_counts.csv")
  if (!file.exists(path)) return(NULL)
  dt <- data.table::fread(path, showProgress = FALSE)
  need <- c("method_order", "method_setup", "setup", "model_label", "selected_k", "unit", "n_topics", "n_items")
  if (!all(need %in% names(dt))) return(NULL)
  dt
}

.m3tb_cached_shared_for_row <- function(cache, row) {
  if (is.null(cache) || !nrow(cache)) return(data.table::data.table())
  out <- cache[
    as.integer(method_order) == as.integer(row$method_order[[1L]]) &
      as.character(method_setup) == as.character(row$method_setup[[1L]]) &
      as.character(setup) == as.character(row$setup[[1L]]) &
      as.character(model_label) == as.character(row$model_label[[1L]]) &
      as.integer(selected_k) == as.integer(row$selected_k[[1L]]) &
      as.character(unit) %chin% c("Links", "Genes", "TFs")
  ]
  data.table::copy(out)
}

.m3tb_standard_link_total_from_summary <- function(summary_counts, row) {
  if (is.null(summary_counts) || !is.finite(summary_counts$scored)) {
    return(NA_real_)
  }
  k <- suppressWarnings(as.numeric(row$selected_k[[1L]]))
  if (!is.finite(k) || k <= 0) return(NA_real_)
  total <- summary_counts$scored / k
  if (!is.finite(total) || total < 0) return(NA_real_)
  round(total)
}

.m3tb_unique_link_pass_from_shared <- function(shared_counts) {
  if (is.null(shared_counts) || !nrow(shared_counts) || !"unit" %in% names(shared_counts)) {
    return(NA_real_)
  }
  out <- shared_counts[as.character(unit) == "Links", sum(as.numeric(n_items), na.rm = TRUE)]
  if (!is.finite(out)) NA_real_ else out
}

.m3tb_pass_counts_from_summary <- function(summary_counts, row, pass_count) {
  total <- .m3tb_standard_link_total_from_summary(summary_counts, row)
  if (!is.finite(total) || !is.finite(pass_count)) {
    return(NULL)
  }
  pass_count <- min(as.numeric(pass_count), total)
  data.table::data.table(
    method_order = row$method_order[[1L]],
    method_setup = row$method_setup[[1L]],
    setup = row$setup[[1L]],
    model_label = row$model_label[[1L]],
    selected_k = as.integer(row$selected_k[[1L]]),
    status = c("Pass", "Fail"),
    count = c(pass_count, total - pass_count),
    count_basis = "Standard TF-FP-gene links"
  )
}

.m3tb_pass_counts_from_topic_rows <- function(summary_counts, row) {
  if (is.null(summary_counts) || !is.finite(summary_counts$scored) || !is.finite(summary_counts$pass)) {
    return(NULL)
  }
  pass_count <- min(as.numeric(summary_counts$pass), as.numeric(summary_counts$scored))
  data.table::data.table(
    method_order = row$method_order[[1L]],
    method_setup = row$method_setup[[1L]],
    setup = row$setup[[1L]],
    model_label = row$model_label[[1L]],
    selected_k = as.integer(row$selected_k[[1L]]),
    status = c("Pass", "Fail"),
    count = c(pass_count, summary_counts$scored - pass_count),
    count_basis = "Topic-link rows"
  )
}

.m3tb_summarize_topic_links <- function(output_dir, model_rows, review_dir = NULL) {
  output_dir <- normalizePath(as.character(output_dir)[[1L]], winslash = "/", mustWork = FALSE)
  model_rows <- data.table::as.data.table(model_rows)
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  review_dir <- normalizePath(as.character(review_dir)[[1L]], winslash = "/", mustWork = FALSE)
  existing_shared <- .m3tb_existing_shared_cache(review_dir)
  fast_summary <- isTRUE(getOption("craftgrn.topic_review.fast_summary", nrow(model_rows) > 24L))
  scan_link_tables <- isTRUE(getOption("craftgrn.topic_review.scan_link_tables", FALSE))
  max_link_file_bytes <- suppressWarnings(as.numeric(
    getOption("craftgrn.topic_review.max_link_file_bytes", 2 * 1024^3)
  ))
  if (is.na(max_link_file_bytes) || max_link_file_bytes <= 0) {
    max_link_file_bytes <- 2 * 1024^3
  }
  workers <- .m3tb_review_worker_count(nrow(model_rows))
  process_one <- function(i) {
    old_dt_threads <- data.table::getDTthreads()
    on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
    if (workers > 1L) data.table::setDTthreads(1L)
    link_assignment_rows <- list()
    item_coverage_rows <- list()
    shared_rows <- list()
    seen <- new.env(parent = emptyenv())
    row <- model_rows[i]
    files <- .m3tb_find_topic_links(output_dir, row)
    if (!length(files)) {
      return(list(
        link_assignment = link_assignment_rows,
        item_coverage = item_coverage_rows,
        shared = shared_rows
      ))
    }
    for (path in files) {
      key <- paste(row$method[[1L]], row$selected_k[[1L]], normalizePath(path, winslash = "/", mustWork = FALSE), sep = "\r")
      if (exists(key, envir = seen, inherits = FALSE)) next
      assign(key, TRUE, envir = seen)
      topic_dir <- dirname(path)
      summary_counts <- .m3tb_read_topic_link_summary(topic_dir)
      item_coverage <- .m3tb_read_topic_item_coverage_counts(topic_dir, row)
      if (!nrow(item_coverage)) {
        item_coverage <- .m3tb_read_topic_term_coverage_counts(topic_dir, row)
      }
      if (nrow(item_coverage)) {
        item_coverage_rows[[length(item_coverage_rows) + 1L]] <- item_coverage
      }
      cached_shared <- .m3tb_cached_shared_for_row(existing_shared, row)
      if (nrow(cached_shared)) {
        shared_rows[[length(shared_rows) + 1L]] <- cached_shared
      }
      summary_pass <- .m3tb_pass_counts_from_summary(
        summary_counts,
        row,
        pass_count = .m3tb_unique_link_pass_from_shared(cached_shared)
      )
      if (!is.null(summary_pass)) {
        link_assignment_rows[[length(link_assignment_rows) + 1L]] <- summary_pass
      }
      pathway_shared <- .m3tb_read_pathway_shared_counts(topic_dir)
      if (nrow(pathway_shared)) {
        pathway_shared[, `:=`(
          method_order = row$method_order[[1L]],
          method_setup = row$method_setup[[1L]],
          setup = row$setup[[1L]],
          model_label = row$model_label[[1L]],
          selected_k = as.integer(row$selected_k[[1L]]),
          unit = "Pathways"
        )]
        shared_rows[[length(shared_rows) + 1L]] <- pathway_shared
      }
      if (!is.null(summary_pass) && nrow(cached_shared)) {
        next
      }
      fast_pass_added <- FALSE
      if (isTRUE(fast_summary) && !nrow(cached_shared)) {
        topic_row_pass <- .m3tb_pass_counts_from_topic_rows(summary_counts, row)
        if (!is.null(topic_row_pass)) {
          link_assignment_rows[[length(link_assignment_rows) + 1L]] <- topic_row_pass
          fast_pass_added <- TRUE
        }
      }
      if (!isTRUE(scan_link_tables)) {
        next
      }
      link_file_bytes <- suppressWarnings(as.numeric(file.info(path)$size[[1L]]))
      if (is.finite(max_link_file_bytes) && is.finite(link_file_bytes) &&
          link_file_bytes > max_link_file_bytes) {
        .log_abort(c(
          "Raw topic-link table exceeds the safe review scan limit.",
          i = "File: {path}",
          i = "Size: {format(link_file_bytes, big.mark = ',')} bytes; limit: {format(max_link_file_bytes, big.mark = ',')} bytes.",
          i = "Use compact topic_link_summary.csv and topic_item_coverage_counts.csv when possible.",
          i = "For an intentional large diagnostic, set options(craftgrn.topic_review.max_link_file_bytes = Inf)."
        ))
      }
      header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
      cols <- intersect(
        c("doc_id", "topic_num", "topic", "tf", "peak_id", "gene_key", "link_pass", "peak_pass", "gene_pass"),
        header
      )
      dt <- data.table::fread(path, select = cols, showProgress = FALSE)
      if (!"topic_num" %in% names(dt) && "topic" %in% names(dt)) data.table::setnames(dt, "topic", "topic_num")
      if (!all(c("tf", "peak_id", "gene_key") %in% names(dt))) next
      if (!"link_pass" %in% names(dt)) dt[, link_pass := TRUE]
      dt[, link_pass := .m3tb_as_flag(link_pass)]
      if (!"peak_pass" %in% names(dt)) dt[, peak_pass := link_pass]
      if (!"gene_pass" %in% names(dt)) dt[, gene_pass := link_pass]
      dt[, `:=`(peak_pass = .m3tb_as_flag(peak_pass), gene_pass = .m3tb_as_flag(gene_pass))]
      link_id <- paste(dt$tf, dt$peak_id, dt$gene_key, sep = "::")
      link_pass_flag <- dt$link_pass | (dt$peak_pass & dt$gene_pass)
      link_status <- data.table::data.table(link_id = link_id, pass = link_pass_flag)
      link_status <- link_status[, .(pass = any(pass, na.rm = TRUE)), by = link_id]
      standard_total <- .m3tb_standard_link_total_from_summary(summary_counts, row)
      if (is.null(summary_pass) && is.finite(standard_total)) {
        pass_count <- sum(link_status$pass, na.rm = TRUE)
        fail_count <- pmax(0, standard_total - pass_count)
        count_basis <- "Standard TF-FP-gene links"
      } else {
        pass_count <- sum(link_status$pass, na.rm = TRUE)
        fail_count <- sum(!link_status$pass, na.rm = TRUE)
        count_basis <- if (grepl("topic_links_pass[.]csv([.]gz)?$", basename(path))) {
          "Pass-only unique links"
        } else {
          "Unique links"
        }
      }
      if (is.null(summary_pass) && !isTRUE(fast_pass_added)) {
        link_assignment_rows[[length(link_assignment_rows) + 1L]] <- data.table::data.table(
          method_order = row$method_order[[1L]],
          method_setup = row$method_setup[[1L]],
          setup = row$setup[[1L]],
          model_label = row$model_label[[1L]],
          selected_k = as.integer(row$selected_k[[1L]]),
          status = c("Pass", "Fail"),
          count = c(pass_count, fail_count),
          count_basis = count_basis
        )
      }
      keep <- dt[link_pass | (peak_pass & gene_pass)]
      if (!nrow(item_coverage) || !all(c("Genes", "TFs") %in% unique(as.character(item_coverage$unit)))) {
        link_coverage <- .topic_item_coverage_from_links(dt, keep)
        if (nrow(link_coverage)) {
          link_coverage <- .m3tb_attach_topic_row_meta(link_coverage, row)
          item_coverage_rows[[length(item_coverage_rows) + 1L]] <- link_coverage
        }
      }
      unit_map <- data.table::data.table(unit_label = c("Links", "Genes", "TFs"), unit_col = c("link_id", "gene_key", "tf"))
      for (u in seq_len(nrow(unit_map))) {
        unit_col <- unit_map$unit_col[[u]]
        if (identical(unit_col, "link_id")) keep[, link_id := paste(tf, peak_id, gene_key, sep = "::")]
        shared <- .m3tb_count_shared_topics(keep, unit_col)
        if (!nrow(shared)) next
        shared[, `:=`(
          method_order = row$method_order[[1L]],
          method_setup = row$method_setup[[1L]],
          setup = row$setup[[1L]],
          model_label = row$model_label[[1L]],
          selected_k = as.integer(row$selected_k[[1L]]),
          unit = unit_map$unit_label[[u]]
        )]
        shared_rows[[length(shared_rows) + 1L]] <- shared
      }
    }
    list(
      link_assignment = link_assignment_rows,
      item_coverage = item_coverage_rows,
      shared = shared_rows
    )
  }
  results <- .m3tb_review_lapply(seq_len(nrow(model_rows)), function(i) {
    tryCatch(process_one(i), error = function(e) .m3tb_review_row_error(i, e))
  }, workers = workers)
  failures <- vapply(results, inherits, logical(1L), "m3tb_review_row_error")
  if (any(failures)) {
    msg <- vapply(results[failures], function(x) {
      sprintf("row %s: %s", x$index, x$message)
    }, character(1L))
    .log_abort(paste(c("Failed to summarize Module 3 topic links.", msg), collapse = "\n"))
  }
  link_assignment_rows <- unlist(lapply(results, `[[`, "link_assignment"), recursive = FALSE)
  item_coverage_rows <- unlist(lapply(results, `[[`, "item_coverage"), recursive = FALSE)
  shared_rows <- unlist(lapply(results, `[[`, "shared"), recursive = FALSE)
  link_assignment <- if (length(link_assignment_rows)) {
    data.table::rbindlist(link_assignment_rows, use.names = TRUE, fill = TRUE)
  } else {
    .m3tb_empty_pass_counts()
  }
  item_coverage <- if (length(item_coverage_rows)) {
    data.table::rbindlist(item_coverage_rows, use.names = TRUE, fill = TRUE)
  } else {
    .m3tb_empty_item_coverage_counts()
  }
  if (nrow(item_coverage)) {
    item_coverage[, count := as.numeric(count)]
    item_coverage[, total := as.numeric(total)]
    item_coverage[, fraction := data.table::fifelse(total > 0, count / total, NA_real_)]
  }
  pass <- if (nrow(item_coverage)) item_coverage else link_assignment
  shared <- if (length(shared_rows)) {
    data.table::rbindlist(shared_rows, use.names = TRUE, fill = TRUE)
  } else {
    .m3tb_empty_shared_counts()
  }
  csv_dir <- .m3tb_review_tables_dir(review_dir)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(item_coverage, file.path(csv_dir, "topic_setup_item_coverage_counts.csv"))
  data.table::fwrite(link_assignment, file.path(csv_dir, "topic_setup_link_assignment_counts.csv"))
  data.table::fwrite(pass, file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  data.table::fwrite(shared, file.path(csv_dir, "topic_setup_shared_topic_counts.csv"))
  list(pass = pass, shared = shared, item_coverage = item_coverage, link_assignment = link_assignment)
}

.m3tb_plot_count_pdf <- function(dt, out_file, title, x_col = "method_setup") {
  grDevices::pdf(out_file, width = 10, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!nrow(dt)) {
    graphics::plot.new()
    graphics::title(title)
    graphics::text(0.5, 0.5, "No rows available")
    return(invisible(out_file))
  }
  plot_dt <- data.table::as.data.table(dt)
  x <- as.character(plot_dt[[x_col]])
  y_source <- if ("count" %in% names(plot_dt)) plot_dt$count else plot_dt$n_items
  label_source <- if ("status" %in% names(plot_dt)) plot_dt$status else plot_dt$unit
  y <- as.numeric(y_source)
  names(y) <- paste(x, label_source, sep = "\n")
  graphics::barplot(y, las = 2, cex.names = 0.65, col = "#4C78A8", main = title, ylab = "Count")
  invisible(out_file)
}

.m3tb_review_theme <- function(base_size = 8.5) {
  ggplot2::theme_bw(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", color = "black", linewidth = 0.7),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

.m3tb_method_label <- function(x) {
  x <- as.character(x)
  data.table::fifelse(
    grepl("vae|multivi", x, ignore.case = TRUE),
    data.table::fifelse(grepl("vae", x, ignore.case = TRUE), "VAE-MLP", "MultiVI"),
    data.table::fifelse(grepl("lda", x, ignore.case = TRUE), "LDA", x)
  )
}

.m3tb_setup_label <- function(x) {
  x <- sub("[|].*$", "", as.character(x))
  trimws(x)
}

.m3tb_json_for_html <- function(x) {
  out <- jsonlite::toJSON(x, dataframe = "rows", auto_unbox = TRUE, null = "null", digits = 8)
  gsub("</", "<\\/", as.character(out), fixed = TRUE)
}

.m3tb_polish_topic_report_html <- function(html) {
  html <- gsub(
    "return{genes:Array.from(out).sort(),byTf:rows,edgeN:edges.size}}function mergeTfSupport(parts)",
    paste0(
      "return{genes:Array.from(out).sort(),byTf:rows,edgeN:edges.size,",
      "context_id:String(ctx.subgrn_context_id||''),",
      "context_label:String(ctx.comparison_label||ctx.comparison_id||'')}}function mergeTfSupport(parts)"
    ),
    html,
    fixed = TRUE
  )
  html <- gsub(
    "return{n:genes.size,genes:detail.length?detail:Array.from(genes).sort(),target_genes:Array.from(genes).sort(),byTf:rows,tfN:rows.length,edgeN:edgeN}}function ensureTfPathwayCounts",
    paste0(
      "const ctxPart=(parts||[]).find(p=>p&&p.context_id)||{};",
      "return{n:genes.size,genes:detail.length?detail:Array.from(genes).sort(),",
      "target_genes:Array.from(genes).sort(),byTf:rows,tfN:rows.length,edgeN:edgeN,",
      "context_id:ctxPart.context_id||'',context_label:ctxPart.context_label||''}}",
      "function chooseTfSupport(parts){const out=(parts||[]).map(p=>mergeTfSupport([p]));",
      "if(!out.length)return mergeTfSupport([]);",
      "out.sort((a,b)=>(Number(b.tfN)||0)-(Number(a.tfN)||0)||",
      "(Number(b.n)||0)-(Number(a.n)||0)||",
      "String(a.context_label||'').localeCompare(String(b.context_label||'')));",
      "return out[0]}function ensureTfPathwayCounts"
    ),
    html,
    fixed = TRUE
  )
  html <- gsub(
    "tfPathwayCounts.set(key,mergeTfSupport(parts));tfPathwayPending.delete(key);drawPathways()",
    "tfPathwayCounts.set(key,chooseTfSupport(parts));tfPathwayPending.delete(key);drawPathways()",
    html,
    fixed = TRUE
  )
  html <- gsub(
    "TF-target pairs represented: '+(Number(support.edgeN)||0)+",
    "TF-target pairs represented: '+(Number(support.edgeN)||0)+(support.context_label?'\\nSubGRN context used for TF counts: '+support.context_label:'')+",
    html,
    fixed = TRUE
  )
  html <- gsub(
    "badgeText=selectedTfUppers.length>1?String(n)+'/'+String(tfN):String(n)",
    paste0(
      "badgeText=selectedTfUppers.length>1?selectedTfUppers.map(k=>{",
      "const r=(support.byTf||[]).find(x=>String(x.tf||'').toUpperCase()===k);",
      "return String(r?Number(r.n)||0:0)}).join('/'):String(n)"
    ),
    html,
    fixed = TRUE
  )
  html <- gsub(
    "selectedTfUppers.length>1?'TF genes/TFs':'TF targets'",
    "selectedTfUppers.length>1?'TF counts':'TF targets'",
    html,
    fixed = TRUE
  )
  html <- gsub(
    "support.genes.join('; ')",
    "support.genes.join(' / ')",
    html,
    fixed = TRUE
  )
  html <- gsub(
    paste0(
      "const detail=[];if(rows.length)detail.push('Summary: '+genes.size+' unique target genes; '+",
      "rows.length+' selected TFs with hits; '+edgeN+' TF-target pairs');rows.forEach"
    ),
    "const detail=[];rows.forEach",
    html,
    fixed = TRUE
  )
  html
}

.m3tb_read_probability_csv <- function(path, id_col) {
  dt <- data.table::fread(path, showProgress = FALSE)
  if (!id_col %in% names(dt)) data.table::setnames(dt, names(dt)[[1L]], id_col)
  ids <- as.character(dt[[id_col]])
  mat <- as.matrix(dt[, setdiff(names(dt), id_col), with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  mat
}

.m3tb_topic_mds_from_phi <- function(phi, theta) {
  phi <- as.matrix(phi)
  storage.mode(phi) <- "numeric"
  row_names <- rownames(phi)
  col_names <- colnames(phi)
  row_topic <- !is.null(row_names) & grepl("^Topic[0-9]+$", row_names)
  col_topic <- !is.null(col_names) & grepl("^Topic[0-9]+$", col_names)
  if (sum(row_topic, na.rm = TRUE) >= 2L && sum(row_topic, na.rm = TRUE) >= sum(col_topic, na.rm = TRUE)) {
    x <- phi
    topic_names <- row_names
  } else {
    x <- t(phi)
    topic_names <- col_names
  }
  if (is.null(topic_names)) topic_names <- paste0("Topic", seq_len(nrow(x)))
  x[!is.finite(x) | x < 0] <- 0
  coords <- if (nrow(x) <= 1L) {
    matrix(c(0, 0), ncol = 2L)
  } else {
    d <- .m3tb_jsd_matrix(x)
    fit <- tryCatch(stats::cmdscale(stats::as.dist(d), k = 2L, eig = FALSE), error = function(e) NULL)
    if (is.null(fit)) {
      cbind(seq_len(nrow(x)), rep(0, nrow(x)))
    } else if (is.null(dim(fit))) {
      cbind(fit, rep(0, length(fit)))
    } else if (ncol(fit) < 2L) {
      cbind(fit[, 1L], rep(0, nrow(fit)))
    } else {
      fit[, 1:2, drop = FALSE]
    }
  }
  if (ncol(coords) < 2L) coords <- cbind(coords[, 1L], rep(0, nrow(coords)))
  theta_mean <- colMeans(as.matrix(theta), na.rm = TRUE)
  mean_theta <- if (all(topic_names %in% names(theta_mean))) {
    as.numeric(theta_mean[topic_names])
  } else if (length(theta_mean) == length(topic_names)) {
    as.numeric(theta_mean)
  } else {
    rep(NA_real_, length(topic_names))
  }
  topic_num <- suppressWarnings(as.integer(sub("^Topic", "", topic_names)))
  topic_num[!is.finite(topic_num)] <- seq_along(topic_num)[!is.finite(topic_num)]
  data.table::data.table(
    topic = topic_names,
    topic_num = topic_num,
    MDS1 = as.numeric(coords[, 1L]),
    MDS2 = as.numeric(coords[, 2L]),
    mean_theta = mean_theta
  )
}

.m3tb_topic_waterfall <- function(theta, context_type) {
  profiles <- .m3tb_topic_profiles(theta, context_type)
  topic_cols <- grep("^Topic[0-9]+$", names(profiles), value = TRUE)
  data.table::melt(
    profiles,
    id.vars = c("comparison_label", "display_label", "n_docs"),
    measure.vars = topic_cols,
    variable.name = "topic",
    value.name = "theta_mean"
  )[, `:=`(
    topic = as.character(topic),
    topic_num = as.integer(sub("^Topic", "", as.character(topic))),
    theta_mean = as.numeric(theta_mean)
  )][order(topic_num, -theta_mean)]
}

.m3tb_tf_topic_rows <- function(theta, context_type) {
  topic_cols <- grep("^Topic[0-9]+$", colnames(theta), value = TRUE)
  if (!length(topic_cols)) return(data.table::data.table())
  doc_ids <- rownames(theta)
  if (is.null(doc_ids)) doc_ids <- paste0("doc", seq_len(nrow(theta)))
  parts <- strsplit(as.character(doc_ids), "::", fixed = TRUE)
  tf <- vapply(parts, function(x) {
    if (length(x) >= 2L) x[[2L]] else NA_character_
  }, character(1L))
  group <- .m3tb_parse_doc_label(doc_ids, context_type)
  base <- data.table::data.table(
    doc_id = as.character(doc_ids),
    comparison_label = as.character(group),
    display_label = .m3tb_display_label(group),
    tf = as.character(tf),
    tf_upper = toupper(as.character(tf))
  )
  keep <- !is.na(base$tf) & nzchar(trimws(base$tf))
  if (!any(keep)) return(data.table::data.table())
  mat <- as.matrix(theta[keep, topic_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  out <- data.table::as.data.table(mat)
  out[, row_id := seq_len(.N)]
  base <- base[keep]
  base[, row_id := seq_len(.N)]
  out <- merge(base, out, by = "row_id", sort = FALSE)
  out[, row_id := NULL]
  out <- data.table::melt(
    out,
    id.vars = c("doc_id", "comparison_label", "display_label", "tf", "tf_upper"),
    measure.vars = topic_cols,
    variable.name = "topic",
    value.name = "theta"
  )
  out[, `:=`(
    topic = as.character(topic),
    topic_num = as.integer(sub("^Topic", "", as.character(topic))),
    theta = as.numeric(theta)
  )]
  out[is.finite(theta)][order(comparison_label, topic_num, -theta, tf)]
}

.m3tb_find_extraction_subdir <- function(row) {
  roots <- unique(c(
    if ("topic_extraction_dir" %in% names(row)) row$topic_extraction_dir[[1L]] else character(),
    if ("selected_k" %in% names(row) && "topic_extraction_dir" %in% names(row)) {
      file.path(row$topic_extraction_dir[[1L]], paste0("K", as.integer(row$selected_k[[1L]])))
    } else {
      character()
    }
  ))
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(NA_character_)
  dirs <- unique(c(roots, unlist(lapply(roots, list.dirs, recursive = TRUE, full.names = TRUE), use.names = FALSE)))
  has_topic <- file.exists(file.path(dirs, "topic_terms.csv")) |
    file.exists(file.path(dirs, "topic_links.csv")) |
    file.exists(file.path(dirs, "topic_links_pass.csv"))
  if (any(has_topic)) return(dirs[which(has_topic)[[1L]]])
  dirs[[1L]]
}

.m3tb_reuse_existing_subgrn_manifest <- function(subgrn_payload_dir, row, k) {
  review_dir <- dirname(normalizePath(subgrn_payload_dir, winslash = "/", mustWork = FALSE))
  manifest_path <- file.path(review_dir, "tables", "pathway_subgrn_manifest.csv")
  if (!file.exists(manifest_path)) return(data.table::data.table())
  dt <- data.table::fread(manifest_path, showProgress = FALSE)
  if (!nrow(dt) || !"payload_file" %in% names(dt)) return(data.table::data.table())
  if ("method_setup" %in% names(dt) && "method_setup" %in% names(row)) {
    dt <- dt[method_setup == row$method_setup[[1L]]]
  }
  if ("run_id" %in% names(dt) && "run_id" %in% names(row)) {
    dt <- dt[run_id == row$run_id[[1L]]]
  }
  if ("k" %in% names(dt)) {
    selected_k <- as.integer(k[[1L]])
    dt[, k_for_filter := suppressWarnings(as.integer(k))]
    dt <- dt[k_for_filter == selected_k][, k_for_filter := NULL][]
  }
  if ("topic_space" %in% names(dt) && "topic_space" %in% names(row)) {
    dt <- dt[topic_space == row$topic_space[[1L]]]
  }
  if ("report_key" %in% names(dt) && "report_key" %in% names(row)) {
    dt <- dt[report_key == row$report_key[[1L]]]
  }
  if (!nrow(dt)) return(data.table::data.table())
  dt <- dt[file.exists(file.path(subgrn_payload_dir, payload_file))]
  dt[]
}

.m3tb_parse_overlap_hits <- function(overlap) {
  vals <- suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))
  vals[!is.finite(vals)] <- NA_integer_
  vals
}

.m3tb_parse_overlap_total <- function(overlap) {
  vals <- suppressWarnings(as.integer(sub("^.*/", "", as.character(overlap))))
  vals[!is.finite(vals)] <- NA_integer_
  vals
}

.m3tb_find_doc_term_file <- function(model_dir) {
  if (is.null(model_dir) || !nzchar(as.character(model_dir)[[1L]])) return(NA_character_)
  files <- c(
    file.path(model_dir, "rds", "doc_term.rds"),
    file.path(model_dir, "doc_term.rds"),
    file.path(model_dir, "doc_term.csv"),
    file.path(model_dir, "doc_term.arrow")
  )
  hit <- files[file.exists(files)]
  if (length(hit)) return(hit[[1L]])
  found <- list.files(model_dir, pattern = "^doc_term[.](rds|csv|arrow)$", recursive = TRUE, full.names = TRUE)
  if (length(found)) found[[1L]] else NA_character_
}

.m3tb_doc_term_gene_universe <- function(doc_term_file) {
  if (is.na(doc_term_file) || !file.exists(doc_term_file)) return(character())
  if (grepl("[.]rds$", doc_term_file, ignore.case = TRUE)) {
    dt <- data.table::as.data.table(readRDS(doc_term_file))
  } else if (grepl("[.]arrow$", doc_term_file, ignore.case = TRUE)) {
    if (!requireNamespace("arrow", quietly = TRUE)) return(character())
    dt <- data.table::as.data.table(arrow::read_ipc_file(doc_term_file))
  } else {
    dt <- data.table::fread(doc_term_file, showProgress = FALSE)
  }
  if (!"term_id" %in% names(dt)) return(character())
  genes <- as.character(dt$term_id)
  genes <- genes[grepl("^GENE:", genes)]
  genes <- sub("^GENE:", "", genes)
  sort(unique(genes[!is.na(genes) & nzchar(genes)]))
}

.m3tb_universe_pathway_cache_file <- function(model_dir, extraction_dir) {
  candidates <- unique(c(
    file.path(model_dir, "topic_pathway_enrichment_gene_universe_all.csv"),
    file.path(extraction_dir, "topic_pathway_enrichment_gene_universe_all.csv"),
    file.path(dirname(extraction_dir), "topic_pathway_enrichment_gene_universe_all.csv")
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(hit[[1L]])
  candidates[[1L]]
}

.m3tb_read_or_compute_universe_pathways <- function(model_dir,
                                                    extraction_dir,
                                                    compute_if_missing = TRUE,
                                                    pathway_species = NULL,
                                                    verbose = FALSE) {
  cache_file <- .m3tb_universe_pathway_cache_file(model_dir, extraction_dir)
  if (file.exists(cache_file)) {
    return(data.table::fread(cache_file, showProgress = FALSE))
  }
  if (!isTRUE(compute_if_missing)) return(data.table::data.table())
  genes <- .m3tb_doc_term_gene_universe(.m3tb_find_doc_term_file(model_dir))
  if (!length(genes)) return(data.table::data.table())
  pathway_backend <- .pathway_backend(NULL)
  pathway_species <- .normalize_pathway_species_mode(pathway_species)
  human_mouse_best <- identical(pathway_species, "human_mouse_best")
  if (!.pathway_backend_available(pathway_backend)) {
    backend_label <- if (identical(pathway_backend, "enrichly")) "enrichly" else "enrichR"
    .log_abort("Computing Module 3 pathway universe totals requires the selected backend: {backend_label}.")
  }
  if (identical(pathway_backend, "enrichr") && !isTRUE(.ensure_enrichr_ready(verbose = FALSE))) {
    .log_abort("Enrichr is not reachable; cannot compute Module 3 pathway universe totals.")
  }
  if (isTRUE(verbose)) {
    .log_inform("Computing Module 3 pathway universe totals for {length(genes)} document genes.")
  }
  run_one_species <- function(species) {
    enr <- tryCatch(
      .run_enrichr_cached(
        genes = genes,
        dbs = .default_pathway_databases(species),
        sleep_time = 0,
        cache_dir = file.path(extraction_dir, "cache", pathway_backend, species),
        backend = pathway_backend,
        pathway_species = species
      ),
      error = function(e) .log_abort("Module 3 pathway universe enrichment failed: {conditionMessage(e)}")
    )
    .topic_enrichr_result_to_table(enr, topic_name = 0L, genes = genes, pathway_species = species)
  }
  out <- if (human_mouse_best) {
    all_out <- data.table::rbindlist(
      lapply(c("human", "mouse"), run_one_species),
      use.names = TRUE,
      fill = TRUE
    )
    .select_best_human_mouse_pathways(all_out)
  } else {
    run_one_species(pathway_species)
  }
  if (nrow(out)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(out, cache_file)
  }
  out
}

.m3tb_apply_universe_pathway_counts <- function(pathways, universe_pathways) {
  dt <- data.table::copy(data.table::as.data.table(pathways))
  if (!nrow(dt)) return(dt)
  uni <- data.table::as.data.table(universe_pathways)
  if (!nrow(uni)) return(dt)
  if (!"pathway" %in% names(uni) && "term" %in% names(uni)) data.table::setnames(uni, "term", "pathway")
  if (!"pathway_key" %in% names(uni)) uni[, pathway_key := pathway]
  if (!"overlap_hits" %in% names(uni)) uni[, overlap_hits := .m3tb_parse_overlap_hits(overlap)]
  if (!"gene_total" %in% names(uni)) uni[, gene_total := .m3tb_parse_overlap_total(overlap)]
  uni[, `:=`(
    pathway = as.character(pathway),
    pathway_key = as.character(pathway_key),
    overlap_hits = suppressWarnings(as.integer(overlap_hits)),
    gene_total = suppressWarnings(as.integer(gene_total))
  )]
  uni[
    is.na(gene_total) | !is.finite(gene_total),
    gene_total := pmax(overlap_hits, 1L, na.rm = TRUE)
  ]
  by_pathway <- uni[
    !is.na(pathway) & nzchar(pathway) & is.finite(gene_total),
    .(gene_total_universe = max(gene_total, na.rm = TRUE)),
    by = pathway
  ]
  if ("gene_total_universe" %in% names(dt)) dt[, gene_total_universe := NULL]
  dt <- merge(dt, by_pathway, by = "pathway", all.x = TRUE, sort = FALSE)
  if (any(is.na(dt$gene_total_universe)) && "pathway_key" %in% names(dt)) {
    by_key <- uni[
      !is.na(pathway_key) & nzchar(pathway_key) & is.finite(gene_total),
      .(gene_total_universe_key = max(gene_total, na.rm = TRUE)),
      by = pathway_key
    ]
    dt <- merge(dt, by_key, by = "pathway_key", all.x = TRUE, sort = FALSE)
    dt[is.na(gene_total_universe), gene_total_universe := gene_total_universe_key]
    dt[, gene_total_universe_key := NULL]
  }
  dt[is.na(gene_total_universe) | !is.finite(gene_total_universe), gene_total_universe := pmax(gene_total, gene_in, na.rm = TRUE)]
  dt[, gene_total_universe := pmax(as.integer(gene_total_universe), as.integer(gene_total), as.integer(gene_in), na.rm = TRUE)]
  dt[, gene_out := pmax(as.integer(gene_total_universe) - as.integer(gene_in), 0L)]
  dt[]
}

.m3tb_read_pathway_tables <- function(extraction_dir,
                                      per_group = FALSE,
                                      model_dir = NULL,
                                      compute_universe = TRUE,
                                      padj_cut = 0.05,
                                      verbose = FALSE,
                                      pathway_file = NULL) {
  empty <- data.table::data.table(
    topic = integer(),
    topic_num = integer(),
    comparison_id = character(),
    direction_group = character(),
    comparison_label = character(),
    display_label = character(),
    pathway = character(),
    pathway_key = character(),
    padj = numeric(),
    combined_score = numeric(),
    overlap = character(),
    overlap_genes = character(),
    gene_in = integer(),
    gene_total = integer(),
    gene_total_universe = integer(),
    gene_out = integer(),
    genes = character()
  )
  if (is.na(extraction_dir) || !dir.exists(extraction_dir)) return(empty)
  if (!is.null(pathway_file)) {
    files <- as.character(pathway_file)
    files <- files[file.exists(files)]
  } else if (isTRUE(per_group)) {
    files <- file.path(
      extraction_dir,
      c(
        "per_comparison_topic_pathway_enrichment.csv",
        "per_condition_topic_pathway_enrichment.csv",
        "topic_term_pathway_enrichment.csv"
      )
    )
    files <- files[file.exists(files)]
    if (!length(files)) {
      files <- list.files(extraction_dir, "_dotplot[.]csv$", recursive = TRUE, full.names = TRUE)
      files <- files[grepl("per_comparison_pathway", files, fixed = TRUE)]
    }
  } else {
    files <- file.path(extraction_dir, "topic_pathway_enrichment_topic_terms.csv")
    files <- files[file.exists(files)]
    if (!length(files)) {
      files <- list.files(extraction_dir, "^topic_pathway_enrichment_.*_dotplot[.]csv$", full.names = TRUE)
    }
  }
  if (!length(files)) return(empty)
  rows <- lapply(files, function(path) {
    dt <- data.table::fread(path, showProgress = FALSE)
    if (!nrow(dt)) return(empty)
    if (!"pathway" %in% names(dt) && "term" %in% names(dt)) data.table::setnames(dt, "term", "pathway")
    if (!"padj" %in% names(dt) && "adjusted_p_value" %in% names(dt)) data.table::setnames(dt, "adjusted_p_value", "padj")
    if (!"combined_score" %in% names(dt)) dt[, combined_score := NA_real_]
    if (!"topic" %in% names(dt) && "topic_num" %in% names(dt)) data.table::setnames(dt, "topic_num", "topic")
    has_condition_id <- "condition_id" %in% names(dt)
    if (!"comparison_id" %in% names(dt) && "condition_id" %in% names(dt)) dt[, comparison_id := condition_id]
    has_comparison_id <- "comparison_id" %in% names(dt)
    has_direction_group <- "direction_group" %in% names(dt)
    if (!"overlap" %in% names(dt)) dt[, overlap := NA_character_]
    if (!"overlap_hits" %in% names(dt)) dt[, overlap_hits := .m3tb_parse_overlap_hits(overlap)]
    if (!"pathway_key" %in% names(dt) && "pathway_norm_key" %in% names(dt)) {
      data.table::setnames(dt, "pathway_norm_key", "pathway_key")
    }
    if (!"pathway_key" %in% names(dt)) dt[, pathway_key := pathway]
    if (!"genes" %in% names(dt)) dt[, genes := NA_character_]
    if (!"overlap_genes" %in% names(dt)) dt[, overlap_genes := genes]
    if (!"comparison_id" %in% names(dt)) dt[, comparison_id := NA_character_]
    if (!"direction_group" %in% names(dt)) dt[, direction_group := NA_character_]
    if (isTRUE(per_group)) {
      if (isTRUE(has_condition_id)) {
        dt[, direction_group := "All"]
        dt[, comparison_label := as.character(comparison_id)]
      } else if (isTRUE(has_comparison_id)) {
        direction_col <- if (isTRUE(has_direction_group)) "direction_group" else NA_character_
        direction_val <- if (!is.na(direction_col)) as.character(dt[[direction_col]]) else "All"
        direction_val <- ifelse(direction_val %in% c("Up", "Target-Up"), "Target-Up",
          ifelse(direction_val %in% c("Down", "Target-Down"), "Target-Down", direction_val)
        )
        dt[, direction_group := .m3tb_subgrn_direction_group(direction_val)]
        dt[, comparison_label := paste(as.character(comparison_id), direction_val, sep = "::")]
      } else {
        label <- sub("_dotplot[.]csv$", "", basename(path))
        label <- sub("_All$", "", label)
        label <- sub("_Target-Up$", "::Target-Up", label)
        label <- sub("_Target-Down$", "::Target-Down", label)
        label <- sub("_Up$", "::Target-Up", label)
        label <- sub("_Down$", "::Target-Down", label)
        dt[, comparison_label := label]
        label_chr <- as.character(dt$comparison_label)
        dt[, `:=`(
          comparison_id = sub("::.*$", "", label_chr),
          direction_group = .m3tb_subgrn_direction_group(ifelse(grepl("::", label_chr, fixed = TRUE), sub("^.*::", "", label_chr), "All"))
        )]
      }
      dt[, display_label := .m3tb_display_label(comparison_label)]
    } else {
      dt[, `:=`(comparison_id = NA_character_, direction_group = NA_character_)]
      dt[, comparison_label := NA_character_]
      dt[, display_label := NA_character_]
    }
    dt[, `:=`(
      topic = as.integer(topic),
      topic_num = as.integer(topic),
      padj = suppressWarnings(as.numeric(padj)),
      combined_score = suppressWarnings(as.numeric(combined_score)),
      gene_in = suppressWarnings(as.integer(overlap_hits)),
      gene_total = .m3tb_parse_overlap_total(overlap),
      pathway = as.character(pathway),
      pathway_key = as.character(pathway_key),
      comparison_id = as.character(comparison_id),
      direction_group = as.character(direction_group),
      overlap_genes = as.character(overlap_genes),
      genes = as.character(genes)
    )]
    dt[is.na(gene_total) | !is.finite(gene_total), gene_total := pmax(gene_in, 1L, na.rm = TRUE)]
    dt[, gene_total_universe := gene_total]
    dt[, gene_out := pmax(gene_total - gene_in, 0L)]
    keep <- intersect(names(empty), names(dt))
    dt[, ..keep]
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) return(empty)
  out <- out[
    is.finite(topic_num) &
      !is.na(pathway) & nzchar(pathway) &
      is.finite(gene_in) & gene_in > 0L &
      is.finite(padj) & padj <= as.numeric(padj_cut)
  ]
  data.table::setorder(out, topic_num, comparison_label, padj, -gene_in, pathway)
  if (!is.null(model_dir)) {
    universe <- .m3tb_read_or_compute_universe_pathways(
      model_dir = model_dir,
      extraction_dir = extraction_dir,
      compute_if_missing = compute_universe,
      verbose = verbose
    )
    out <- .m3tb_apply_universe_pathway_counts(out, universe)
  }
  out[]
}

.m3tb_read_condition_pathway_tables <- function(extraction_dir,
                                                model_dir = NULL,
                                                compute_universe = TRUE,
                                                padj_cut = 0.05,
                                                verbose = FALSE,
                                                pathway_file = NULL) {
  pathways <- .m3tb_read_pathway_tables(
    extraction_dir,
    per_group = TRUE,
    model_dir = model_dir,
    compute_universe = compute_universe,
    padj_cut = padj_cut,
    verbose = verbose,
    pathway_file = pathway_file
  )
  if (nrow(pathways) || !is.null(pathway_file)) return(pathways)
  .m3tb_read_pathway_tables(
    extraction_dir,
    per_group = FALSE,
    model_dir = model_dir,
    compute_universe = compute_universe,
    padj_cut = padj_cut,
    verbose = verbose
  )
}

.m3tb_topic_report_html <- function(title,
                                    topic_mds,
                                    waterfall,
                                    tf_topic = NULL,
                                    pathways,
                                    out_html,
                                    condition_pathways = NULL,
                                    subgrn_manifest = NULL,
                                    subgrn_payload_base = "../../pathway_subgrn_payloads") {
  topic_json <- .m3tb_json_for_html(topic_mds)
  waterfall_json <- .m3tb_json_for_html(waterfall)
  if (is.null(tf_topic)) tf_topic <- data.table::data.table()
  tf_topic_json <- .m3tb_json_for_html(tf_topic)
  pathway_json <- .m3tb_json_for_html(pathways)
  if (is.null(condition_pathways)) condition_pathways <- data.table::data.table()
  condition_pathway_json <- .m3tb_json_for_html(condition_pathways)
  if (is.null(subgrn_manifest)) subgrn_manifest <- .m3tb_subgrn_empty_manifest()
  subgrn_json <- .m3tb_json_for_html(subgrn_manifest)
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    "<meta name=\"craftgrn-module3-report-schema\" content=\"2\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    "<style>",
    "html,body{width:100%;height:100%;overflow:hidden}body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif;font-weight:700}",
    ".top{height:48px;background:#fff;border-bottom:1px solid #d6d6d0;display:flex;justify-content:space-between;gap:10px;align-items:center;padding:6px 12px;box-sizing:border-box}",
    "h1{font-size:18px;line-height:1.05;margin:0}.meta{font-size:11px;color:#475569;line-height:1.1}.controls{display:flex;gap:8px;align-items:center;white-space:nowrap}",
    "select,input,button{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #aaa;border-radius:3px;background:#fff;padding:5px 7px}button{background:#111;color:#fff}",
    ".grid{height:calc(100vh - 48px);display:grid;grid-template-columns:minmax(520px,1fr) minmax(700px,1.35fr);gap:8px;padding:8px;box-sizing:border-box}",
    ".left{display:grid;grid-template-rows:minmax(0,1fr) minmax(0,1fr);gap:8px;min-height:0}.pane{background:#fff;border:1px solid #d6d6d0;display:flex;flex-direction:column;min-height:0;overflow:hidden}",
    ".paneHead{padding:8px 10px;border-bottom:1px solid #e5e5df;display:flex;justify-content:space-between;gap:8px;align-items:center}.pane h2{font-size:20px;line-height:1.05;margin:0}.paneMeta{display:flex;flex-direction:column;gap:3px;align-items:flex-end;min-width:0}.contextBadge{max-width:560px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;background:#111;color:#fff;border-radius:3px;padding:4px 7px;font-size:12px;line-height:1.1}.barControls{display:flex;gap:8px;align-items:center;white-space:nowrap}.barControls select{max-width:180px}.barControls input{width:130px}.barControls select[multiple]{height:46px;min-width:130px}.body{position:relative;flex:1;min-height:0}.splitWaterfall{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);height:100%;min-height:0}.splitSide{position:relative;min-width:0;min-height:0}.splitSide:first-child{border-right:1px solid #e5e5df}.clearTf{display:none;margin-left:8px;padding:3px 6px;font-size:11px}.note{font-size:10px;color:#555;border-top:1px solid #e5e5df;padding:5px 8px;line-height:1.1}.tooltip{position:absolute;display:none;background:rgba(17,17,17,0.92);color:#fff;font:700 12px Arial,Helvetica,sans-serif;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:360px;line-height:1.35;z-index:5}body.embed .top{display:none}body.embed .grid{height:100vh}",
    "svg{width:100%;height:100%;display:block}.label{font-size:12px;font-weight:700}.small{font-size:10px}.barIn{fill:#cc454b}.barOut{fill:#a9cfe5}.pathAxis{stroke:#111;stroke-width:1.6;shape-rendering:crispEdges}.pathTick{stroke:#777;stroke-width:.9;shape-rendering:crispEdges}.pathLabel{font-size:16px;font-weight:700;fill:#111}.pathLabelTopicSpecific{fill:#cc2f3c}.pathLabelGroupSpecific{fill:#2563eb}.pathLabelBothSpecific{fill:#7e22ce}.pathTickText{font-size:14px;font-weight:700;fill:#111}.pathLegendText{font-size:14px;font-weight:700;fill:#334155}.mdsLeader{stroke:#334155;stroke-width:1.4;opacity:.78}.mdsLabel{font-size:22px;font-weight:700;fill:#111;paint-order:stroke;stroke:#fff;stroke-width:6px;stroke-linejoin:round}",
    .m3tb_subgrn_css(),
    "@media(max-width:1100px){.grid{grid-template-columns:1fr 1fr}.top{height:52px}.grid{height:calc(100vh - 52px)}}",
    "</style></head><body>",
    "<div class=\"top\">",
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    "<div class=\"controls\"><label>Palette <select id=\"paletteSelect\"><option value=\"default\" selected>Default</option><option value=\"aaas\">AAAS</option><option value=\"d3\">D3</option><option value=\"npg\">NPG</option></select></label><label>Topic <select id=\"topicSelect\"></select></label><label>Pathway scope <select id=\"pathwayScopeSelect\"></select></label><button id=\"exportSvgButton\">Export SVG</button></div></div>",
    "<main class=\"grid\"><div class=\"left\">",
    "<section class=\"pane\"><div class=\"paneHead\"><h2>Intertopic Distance Map</h2><span id=\"mdsStats\" class=\"meta\"></span></div><div class=\"body\"><div class=\"tooltip\" id=\"mdsTooltip\"></div><svg id=\"mdsSvg\" viewBox=\"0 0 900 620\"><g id=\"mdsLayer\"></g></svg></div><div class=\"note\">Topics are embedded by Jensen-Shannon distance.</div></section>",
    "<section class=\"pane\"><div class=\"paneHead\"><h2>Topic Bar Plots</h2><div class=\"barControls\"><label>TF <input id=\"tfSearchInput\" list=\"tfOptionList\" placeholder=\"Type TF\"/><datalist id=\"tfOptionList\"></datalist></label><label>Selected <select id=\"tfSelect\" multiple size=\"3\"></select></label><span id=\"waterfallStats\" class=\"meta\"></span><button id=\"clearTfButton\" class=\"clearTf\" type=\"button\">Clear TF</button></div></div><div class=\"body\"><div class=\"tooltip\" id=\"wfTooltip\"></div><div class=\"splitWaterfall\"><div class=\"splitSide\"><svg id=\"tfWaterfallSvg\" viewBox=\"0 0 760 760\"><g id=\"tfWaterfallLayer\"></g></svg></div><div class=\"splitSide\"><svg id=\"waterfallSvg\" viewBox=\"0 0 760 760\"><g id=\"waterfallLayer\"></g></svg></div></div></div><div class=\"note\">Left: TF-to-topic probability from theta. Right: condition/comparison mean theta for the selected topic.</div></section>",
    "</div><section class=\"pane\"><div class=\"paneHead\"><h2>Pathways</h2><div class=\"paneMeta\"><div id=\"selectedContext\" class=\"contextBadge\"></div><span id=\"pathStats\" class=\"meta\"></span></div></div><div class=\"body\"><div class=\"tooltip\" id=\"pathTooltip\"></div><svg id=\"pathSvg\" viewBox=\"0 0 1500 980\"><g id=\"pathLayer\"></g></svg>",
    .m3tb_subgrn_html_panel(),
    "</div><div class=\"note\">Filter: N_gene >= 5, adjusted p-value < 0.05, top 30 per topic or selected condition. X-axis is full gene-universe pathway hits; red is topic overlap and blue is remaining universe hits. Click a pathway bar to open the sub-GRN browser.</div></section></main>",
    "<script>",
    "if(new URLSearchParams(location.search).get('embed')==='1')document.body.classList.add('embed');",
    paste0("const TOPICS=", topic_json, ";"),
    paste0("const WATERFALL=", waterfall_json, ";"),
    paste0("const TF_TOPIC=", tf_topic_json, ";"),
    paste0("const PATHWAYS_TOPIC=", pathway_json, ";"),
    paste0("const PATHWAYS_CONDITION=", condition_pathway_json, ";"),
    "const PATHWAYS=PATHWAYS_TOPIC;",
    paste0("const SUBGRN_MANIFEST=", subgrn_json, ";"),
    paste0("const SUBGRN_PAYLOAD_BASE=", jsonlite::toJSON(subgrn_payload_base, auto_unbox = TRUE), ";"),
    "const COLOR_PRESETS={default:['#2563eb','#dc2626','#16a34a','#ca8a04','#7c3aed','#0891b2','#be123c','#4b5563','#65a30d','#c2410c'],npg:['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85'],aaas:['#EE0000','#3B4992','#008B45','#631879','#008280','#BB0021','#5F559B','#A20056','#808180','#1B1919'],d3:['#1F77B4','#FF7F0E','#2CA02C','#D62728','#9467BD','#8C564B','#E377C2','#7F7F7F','#BCBD22','#17BECF']};",
    "const topicSelect=document.getElementById('topicSelect'),paletteSelect=document.getElementById('paletteSelect'),tfSelect=document.getElementById('tfSelect'),tfSearchInput=document.getElementById('tfSearchInput'),tfOptionList=document.getElementById('tfOptionList'),pathwayScopeSelect=document.getElementById('pathwayScopeSelect'),groupSelect=pathwayScopeSelect,mdsSvg=document.getElementById('mdsSvg'),tfWaterfallSvg=document.getElementById('tfWaterfallSvg'),waterfallSvg=document.getElementById('waterfallSvg'),pathSvg=document.getElementById('pathSvg'),mdsLayer=document.getElementById('mdsLayer'),tfWaterfallLayer=document.getElementById('tfWaterfallLayer'),waterfallLayer=document.getElementById('waterfallLayer'),pathLayer=document.getElementById('pathLayer'),mdsTooltip=document.getElementById('mdsTooltip'),wfTooltip=document.getElementById('wfTooltip'),pathTooltip=document.getElementById('pathTooltip'),clearTfButton=document.getElementById('clearTfButton');let selectedTfUppers=[],selectedTfUpper='',tfLabelMap=new Map();const tfPathwayCounts=new Map(),tfPathwayPending=new Set();",
    "function esc(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;')}function el(n,a){const x=document.createElementNS('http://www.w3.org/2000/svg',n);Object.entries(a||{}).forEach(([k,v])=>x.setAttribute(k,v));return x}function sc(vals,lo,hi){const xs=vals.map(Number).filter(Number.isFinite),a=Math.min(...xs),b=Math.max(...xs);return v=>!Number.isFinite(Number(v))||a===b?(lo+hi)/2:lo+(Number(v)-a)/(b-a)*(hi-lo)}function tooltip(tt,ev,text){tt.innerHTML=esc(text).replace(/\\n/g,'<br/>');tt.style.display='block';tt.style.left=(ev.offsetX+12)+'px';tt.style.top=(ev.offsetY+12)+'px'}function svgFontSize(svg,targetPx){const m=svg.getScreenCTM();const s=m?Math.sqrt(Math.abs(m.a*m.d)):1;return Math.max(1,targetPx/Math.max(s,0.001)).toFixed(2)}",
    "function activePalette(){return COLOR_PRESETS[paletteSelect.value]||COLOR_PRESETS.default}function color(key){const colors=activePalette();let h=2166136261;const s=String(key||'NA');for(let i=0;i<s.length;i++){h^=s.charCodeAt(i);h=Math.imul(h,16777619)}return colors[(h>>>0)%colors.length]}function topicNum(t){return Number(String(t).replace(/^Topic/,''))}",
    "function exportStyle(){return 'text{font-family:Arial,Helvetica,sans-serif;font-weight:700}.label{font-size:12px;font-weight:700}.small{font-size:10px}.barIn{fill:#cc454b}.barOut{fill:#a9cfe5}.pathAxis{stroke:#111;stroke-width:1.6;shape-rendering:crispEdges}.pathTick{stroke:#777;stroke-width:.9;shape-rendering:crispEdges}.pathLabel{font-size:16px;font-weight:700;fill:#111}.pathLabelTopicSpecific{fill:#cc2f3c}.pathLabelGroupSpecific{fill:#2563eb}.pathLabelBothSpecific{fill:#7e22ce}.pathTickText{font-size:14px;font-weight:700;fill:#111}.pathLegendText{font-size:14px;font-weight:700;fill:#334155}.mdsLeader{stroke:#334155;stroke-width:1.4;opacity:.78}.mdsLabel{font-size:22px;font-weight:700;fill:#111;paint-order:stroke;stroke:#fff;stroke-width:6px;stroke-linejoin:round}'}function exportSvgTitle(root,text,x,y,size){const t=el('text',{x:x,y:y,'font-size':size,'font-weight':700,fill:'#111'});t.textContent=text;root.appendChild(t)}function addExportPanel(root,src,title,x,y,w,h){exportSvgTitle(root,title,x,y+22,24);const clone=src.cloneNode(true);clone.setAttribute('x',x);clone.setAttribute('y',y+38);clone.setAttribute('width',w);clone.setAttribute('height',h-38);clone.setAttribute('preserveAspectRatio','xMidYMid meet');clone.removeAttribute('style');root.appendChild(clone)}function exportSvg(){draw();const root=el('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 2480 1400',width:2480,height:1400});const style=el('style',{});style.textContent=exportStyle();root.appendChild(style);root.appendChild(el('rect',{x:0,y:0,width:2480,height:1400,fill:'#fff'}));exportSvgTitle(root,document.title||'Topic report',30,42,28);addExportPanel(root,mdsSvg,'Intertopic Distance Map',30,74,760,590);addExportPanel(root,waterfallSvg,'Condition Bar Plot',30,720,760,620);addExportPanel(root,pathSvg,'Pathways',830,74,1600,1266);const text=new XMLSerializer().serializeToString(root),blob=new Blob([text],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=(document.title||'topic_report').replace(/[^A-Za-z0-9_.-]+/g,'_')+'.svg';document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(a.href)}",
    "function selectedTfSet(){return new Set(selectedTfUppers)}function hasSelectedTf(){return selectedTfUppers.length>0}function selectedTfKey(){return selectedTfUppers.slice().sort().join(',')}function selectedTfLabel(){return selectedTfUppers.join(', ')}function tfKeyFromInput(x){const raw=String(x||'').trim();if(!raw)return'';const up=raw.toUpperCase();if(tfLabelMap.has(up))return up;for(const [key,label] of tfLabelMap.entries()){if(String(label).toUpperCase()===up)return key}return''}function updateTfControls(){const keep=selectedTfSet();if(tfSelect)Array.from(tfSelect.options).forEach(o=>{o.selected=keep.has(o.value)});selectedTfUpper=selectedTfUppers[0]||'';if(clearTfButton)clearTfButton.style.display=hasSelectedTf()?'inline-block':'none'}function setSelectedTfs(vals){const seen=new Set(),out=[];(vals||[]).forEach(v=>{const k=tfKeyFromInput(v)||String(v||'').toUpperCase();if(k&&tfLabelMap.has(k)&&!seen.has(k)){seen.add(k);out.push(k)}});selectedTfUppers=out;updateTfControls();draw()}function setSelectedTf(tf){const k=tfKeyFromInput(tf)||String(tf||'').toUpperCase();if(!k||!tfLabelMap.has(k))return;const s=selectedTfSet();if(s.has(k))s.delete(k);else s.add(k);setSelectedTfs(Array.from(s))}function addTfFromSearch(){const k=tfKeyFromInput(tfSearchInput?tfSearchInput.value:'');if(k){const s=selectedTfSet();s.add(k);if(tfSearchInput)tfSearchInput.value='';setSelectedTfs(Array.from(s))}}function fillTfSelect(){const old=selectedTfUppers.slice(),seen=new Map();TF_TOPIC.forEach(d=>{const key=String(d.tf_upper||d.tf||'').toUpperCase();if(key&&!seen.has(key))seen.set(key,String(d.tf||key))});const rows=Array.from(seen.entries()).map(([key,label])=>({key,label})).sort((a,b)=>a.label.localeCompare(b.label,undefined,{sensitivity:'base'})||a.key.localeCompare(b.key));tfLabelMap=new Map(rows.map(d=>[d.key,d.label]));tfSelect.replaceChildren();if(tfOptionList)tfOptionList.replaceChildren();rows.forEach(d=>{const o=document.createElement('option');o.value=d.key;o.textContent=d.label;tfSelect.appendChild(o);if(tfOptionList){const q=document.createElement('option');q.value=d.label;tfOptionList.appendChild(q)}});selectedTfUppers=old.filter(k=>tfLabelMap.has(k));updateTfControls()}",
    "function selectedScopeLabel(){const opt=pathwayScopeSelect.options[pathwayScopeSelect.selectedIndex];return pathwayScopeSelect.value==='__overall__'?'Overall topic pathways':(opt?opt.textContent:pathwayScopeSelect.value)}function updateSelectedContext(){const x=document.getElementById('selectedContext');if(x)x.textContent='Selected: Topic '+topicNum(topicSelect.value)+' | '+selectedScopeLabel()}",
    "function init(){TOPICS.sort((a,b)=>topicNum(a.topic)-topicNum(b.topic)).forEach(d=>{const o=document.createElement('option');o.value=d.topic;o.textContent='Topic '+topicNum(d.topic);topicSelect.appendChild(o)});fillPathwayScopes();fillTfSelect();topicSelect.addEventListener('change',draw);paletteSelect.addEventListener('change',draw);pathwayScopeSelect.addEventListener('change',draw);tfSelect.addEventListener('change',()=>setSelectedTfs(Array.from(tfSelect.selectedOptions).map(o=>o.value)));if(tfSearchInput){tfSearchInput.addEventListener('change',addTfFromSearch);tfSearchInput.addEventListener('keydown',ev=>{if(ev.key==='Enter'){ev.preventDefault();addTfFromSearch()}})}if(clearTfButton)clearTfButton.addEventListener('click',()=>setSelectedTfs([]));document.getElementById('exportSvgButton').addEventListener('click',exportSvg);draw()}",
    "function drawTopic(){mdsLayer.replaceChildren();const w=900,h=620,p=86,labelSize=svgFontSize(mdsSvg,34),labelNum=Number(labelSize),sx=sc(TOPICS.map(d=>+d.MDS1),p,w-p),sy=sc(TOPICS.map(d=>+d.MDS2),h-p,p),mx=Math.max(...TOPICS.map(d=>+d.mean_theta||0),1e-6),sel=topicSelect.value,cx=w/2,cy=h/2;const pts=TOPICS.map((d,i)=>{const x=sx(+d.MDS1),y=sy(+d.MDS2),r=28+52*Math.sqrt((+d.mean_theta||0)/mx);return{d,i,x,y,r,meanTheta:+d.mean_theta||0}});const selectTopic=p=>{topicSelect.value=p.d.topic;draw()};pts.slice().sort((a,b)=>((a.d.topic===sel)?1:0)-((b.d.topic===sel)?1:0)||b.r-a.r).forEach(p=>{const isSel=p.d.topic===sel,c=el('circle',{cx:p.x,cy:p.y,r:p.r,fill:color(p.d.topic),opacity:isSel?0.9:0.42,stroke:isSel?'#111':'#fff','stroke-width':isSel?6:2.6,style:'cursor:pointer'});c.addEventListener('click',()=>selectTopic(p));c.addEventListener('mousemove',ev=>tooltip(mdsTooltip,ev,p.d.topic+'\\nmean document-to-topic probability: '+p.meanTheta.toFixed(4)));c.addEventListener('mouseout',()=>mdsTooltip.style.display='none');mdsLayer.appendChild(c)});const labels=pts.map(p=>{const text=String(p.d.topic).replace('Topic','');return{p,text,x:p.x,y:p.y,ax:p.x,ay:p.y,w:Math.max(34,text.length*labelNum*.72+18),h:labelNum+18,displaced:false}});function clampLabel(a){a.x=Math.max(a.w/2+8,Math.min(w-a.w/2-8,a.x));a.y=Math.max(a.h/2+8,Math.min(h-a.h/2-8,a.y))}function setOutward(a,extra){let dx=a.p.x-cx,dy=a.p.y-cy,len=Math.hypot(dx,dy);if(len<1){const ang=(a.p.i+1)*2.399963;dx=Math.cos(ang);dy=Math.sin(ang);len=1}const gap=a.p.r+42+(extra||0);a.ax=Math.max(a.w/2+8,Math.min(w-a.w/2-8,a.p.x+dx/len*gap));a.ay=Math.max(a.h/2+8,Math.min(h-a.h/2-8,a.p.y+dy/len*gap));if(!a.displaced){a.x=a.ax;a.y=a.ay}a.displaced=true}function overlaps(a,b,pad){return Math.abs(a.x-b.x)<(a.w+b.w)/2+pad&&Math.abs(a.y-b.y)<(a.h+b.h)/2+pad}for(let i=0;i<labels.length;i++){for(let j=i+1;j<labels.length;j++){if(overlaps(labels[i],labels[j],10)){setOutward(labels[i],(i%3)*8);setOutward(labels[j],(j%3)*8)}}}for(let iter=0;iter<120;iter++){labels.forEach(a=>{a.x+=(a.ax-a.x)*(a.displaced?.055:.06);a.y+=(a.ay-a.y)*(a.displaced?.055:.06);clampLabel(a)});for(let i=0;i<labels.length;i++){for(let j=i+1;j<labels.length;j++){const a=labels[i],b=labels[j],ox=(a.w+b.w)/2+10-Math.abs(a.x-b.x),oy=(a.h+b.h)/2+10-Math.abs(a.y-b.y);if(ox>0&&oy>0){if(!a.displaced&&!b.displaced){setOutward(a,(i%3)*8);setOutward(b,(j%3)*8)}let dx=a.x-b.x,dy=a.y-b.y;if(Math.abs(dx)<.1&&Math.abs(dy)<.1){dx=(a.p.i<=b.p.i?-1:1);dy=((a.p.i+b.p.i)%2?-1:1)}const len=Math.max(.1,Math.hypot(dx,dy)),push=Math.min(ox,oy)*.22;if(a.displaced&&b.displaced){a.x+=dx/len*push;a.y+=dy/len*push;b.x-=dx/len*push;b.y-=dy/len*push}else if(a.displaced){a.x+=dx/len*push*1.5;a.y+=dy/len*push*1.5}else if(b.displaced){b.x-=dx/len*push*1.5;b.y-=dy/len*push*1.5}clampLabel(a);clampLabel(b)}}}}labels.forEach(a=>{const p=a.p,dx=a.x-p.x,dy=a.y-p.y,dist=Math.hypot(dx,dy),moved=a.displaced&&dist>Math.max(10,p.r*.32);if(moved){const len=Math.max(1,dist),x1=p.x+dx/len*p.r,y1=p.y+dy/len*p.r,line=el('line',{x1:x1,y1:y1,x2:a.x,y2:a.y,class:'mdsLeader'});line.addEventListener('click',()=>selectTopic(p));mdsLayer.appendChild(line)}const t=el('text',{x:a.x,y:a.y+labelNum*.18,'text-anchor':'middle',class:'mdsLabel','font-size':labelSize,style:'cursor:pointer'});t.textContent=a.text;t.addEventListener('click',()=>selectTopic(p));t.addEventListener('mousemove',ev=>tooltip(mdsTooltip,ev,p.d.topic+'\\nmean document-to-topic probability: '+p.meanTheta.toFixed(4)));t.addEventListener('mouseout',()=>mdsTooltip.style.display='none');mdsLayer.appendChild(t)});document.getElementById('mdsStats').textContent=TOPICS.length+' topics'}",
    "function currentTfTopicRows(){const topic=topicSelect.value,scope=pathwayScopeSelect.value,selSet=selectedTfSet();let rows=TF_TOPIC.filter(d=>d.topic===topic);function addSelected(top,sorted){selectedTfUppers.forEach(k=>{if(!top.some(d=>String(d.tf_upper||d.tf||'').toUpperCase()===k)){const sel=sorted.find(d=>String(d.tf_upper||d.tf||'').toUpperCase()===k);if(sel)top.push(sel)}});return top}if(scope&&scope!=='__overall__'){rows=rows.filter(d=>String(d.comparison_label||'')===String(scope));const sorted=rows.sort((a,b)=>(+b.theta)-(+a.theta)||String(a.tf||'').localeCompare(String(b.tf||''))),top=addSelected(sorted.slice(0,30),sorted);return top}const best=new Map();rows.forEach(d=>{const k=String(d.tf_upper||d.tf||'').toUpperCase();if(!k)return;const old=best.get(k);if(!old||(+d.theta)>(+old.theta))best.set(k,d)});const sorted=Array.from(best.values()).sort((a,b)=>(+b.theta)-(+a.theta)||String(a.tf||'').localeCompare(String(b.tf||''))),top=addSelected(sorted.slice(0,30),sorted);return top}",
    "function drawTfWaterfall(){tfWaterfallLayer.replaceChildren();const rows=currentTfTopicRows(),selSet=selectedTfSet(),w=760,h=760,left=184,right=42,top=82,labelY=42,titleY=22,labelSize=svgFontSize(tfWaterfallSvg,12),smallSize=svgFontSize(tfWaterfallSvg,12),maxVal=Math.max(...rows.map(d=>+d.theta||0),.001),ticks=niceTicks(maxVal,5),axisMax=Math.max(...ticks,maxVal),plotW=w-left-right,rowH=Math.max(18,Math.min(34,(h-58-top-10)/Math.max(1,rows.length))),x=v=>left+Math.max(0,+v||0)/axisMax*plotW,selectRow=d=>setSelectedTf(d.tf_upper||d.tf||'');ticks.forEach(t=>{const xx=x(t);tfWaterfallLayer.appendChild(el('line',{x1:xx,y1:top-7,x2:xx,y2:h-56,stroke:'#777','stroke-width':.8,opacity:t===0?1:.32}));tfWaterfallLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=Number(t).toFixed(2).replace(/\\.?0+$/,'')});tfWaterfallLayer.appendChild(el('line',{x1:left,y1:56,x2:w-right,y2:56,stroke:'#111','stroke-width':1.8}));tfWaterfallLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent='TF theta for selected topic';rows.forEach((d,i)=>{const y=top+i*rowH,val=+d.theta||0,bw=x(val)-left,key=String(d.tf_upper||d.tf||'').toUpperCase(),isSel=selSet.has(key),label=String(d.tf||'TF'),shortLabel=label.length>18?label.slice(0,16)+'..':label,msg=label+'\\ntheta: '+val.toFixed(4)+'\\n'+String(d.display_label||d.comparison_label||'');const lab=el('text',{x:left-10,y:y+rowH*.62,'font-size':labelSize,'text-anchor':'end','font-weight':700,fill:isSel?'#991b1b':'#111',style:'cursor:pointer'});lab.textContent=shortLabel;lab.addEventListener('click',()=>selectRow(d));lab.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));lab.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(lab);const r=el('rect',{x:left,y:y+4,width:Math.max(1,bw),height:Math.max(5,rowH-9),fill:isSel?'#dc2626':color(d.tf_upper||d.tf),opacity:isSel?.95:.76,stroke:isSel?'#111':'none','stroke-width':isSel?3:0,style:'cursor:pointer'});r.addEventListener('click',()=>selectRow(d));r.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));r.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(r);const vt=el('text',{x:Math.min(w-right,left+bw+6),y:y+rowH*.62,'font-size':smallSize,'font-weight':700,fill:'#334155',style:'cursor:pointer'});vt.textContent=val.toFixed(3);vt.addEventListener('click',()=>selectRow(d));vt.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));vt.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(vt)});updateTfControls()}",
    "function tfCountKey(d){return selectedTfKey()+'|'+String(d.comparison_label||d.comparison_id||'')+'|'+String(d.topic_num||d.topic||'')+'|'+String(d.pathway_key||d.pathway_norm_key||d.pathway||'')}function splitGenes(s){return String(s||'').split(/[;,]/).map(x=>x.trim()).filter(Boolean)}function rowGeneSet(d){return new Set(splitGenes(d.genes||d.overlap_genes||''))}function contextMatchEdge(e,ctx,genes){return e.subgrn_context_id?e.subgrn_context_id===ctx.subgrn_context_id:(String(e.comparison_id||'')===String(ctx.comparison_id||'')&&String(e.direction||'')===String(ctx.direction||'')&&genes.has(String(e.gene_key||'')))}function countTfGenesInPayload(payload,ctx,row){const genes=rowGeneSet(row),allowed=subgrnGeneSet(ctx),tfSet=selectedTfSet();allowed.forEach(g=>genes.add(g));const out=new Set(),byTf=new Map(),edges=new Set();((payload&&payload.tf_gene_edges)||[]).forEach(e=>{const tf=String(e.tf_upper||e.tf||'').toUpperCase(),gene=String(e.gene_key||'');if(tfSet.has(tf)&&genes.has(gene)&&contextMatchEdge(e,ctx,genes)){out.add(gene);const rec=byTf.get(tf)||{tf:tf,label:tfLabelMap.get(tf)||String(e.tf||tf),genes:new Set()};rec.genes.add(gene);byTf.set(tf,rec);edges.add(tf+'\\t'+gene)}});const rows=Array.from(byTf.values()).map(r=>({tf:r.tf,label:r.label,genes:Array.from(r.genes).sort(),n:r.genes.size})).sort((a,b)=>b.n-a.n||a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));return{genes:Array.from(out).sort(),byTf:rows,edgeN:edges.size}}function mergeTfSupport(parts){const genes=new Set(),byTf=new Map();let edgeN=0;(parts||[]).forEach(p=>{if(Array.isArray(p))p={genes:p,byTf:[],edgeN:0};(p.genes||[]).forEach(g=>genes.add(g));(p.byTf||[]).forEach(r=>{const key=String(r.tf||'').toUpperCase();if(!key)return;const rec=byTf.get(key)||{tf:key,label:r.label||tfLabelMap.get(key)||key,genes:new Set()};(r.genes||[]).forEach(g=>rec.genes.add(g));byTf.set(key,rec)});edgeN+=Number(p.edgeN)||0});const rows=Array.from(byTf.values()).map(r=>({tf:r.tf,label:r.label,genes:Array.from(r.genes).sort(),n:r.genes.size})).sort((a,b)=>b.n-a.n||a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));const detail=[];if(rows.length)detail.push('Summary: '+genes.size+' unique target genes; '+rows.length+' selected TFs with hits; '+edgeN+' TF-target pairs');rows.forEach(r=>detail.push(r.label+' ('+r.n+'): '+r.genes.join(',')));return{n:genes.size,genes:detail.length?detail:Array.from(genes).sort(),target_genes:Array.from(genes).sort(),byTf:rows,tfN:rows.length,edgeN:edgeN}}function ensureTfPathwayCounts(rows){if(!hasSelectedTf())return;rows.forEach(d=>{const key=tfCountKey(d);if(tfPathwayCounts.has(key)||tfPathwayPending.has(key))return;const ctxs=subgrnContextRows(d);if(!ctxs.length){tfPathwayCounts.set(key,mergeTfSupport([]));return}tfPathwayPending.add(key);Promise.all(ctxs.map(ctx=>loadSubgrnPayload(ctx.payload_file).then(payload=>countTfGenesInPayload(payload,ctx,d)).catch(()=>({genes:[],byTf:[],edgeN:0})))).then(parts=>{tfPathwayCounts.set(key,mergeTfSupport(parts));tfPathwayPending.delete(key);drawPathways()})})}",
    "function drawGroupWaterfall(){waterfallLayer.replaceChildren();const topic=topicSelect.value,allRows=WATERFALL.filter(d=>d.topic===topic).sort((a,b)=>Number(b.theta_mean)-Number(a.theta_mean)),rows=allRows.slice(0,25),w=760,h=760,left=345,right=30,top=96,axisY=62,labelY=52,titleY=28,labelSize=svgFontSize(waterfallSvg,10),smallSize=svgFontSize(waterfallSvg,9),maxVal=Math.max(...WATERFALL.map(d=>Number(d.theta_mean)).filter(Number.isFinite),.001),rawTick=maxVal/4,tickPow=Math.pow(10,Math.floor(Math.log10(rawTick))),tickStep=([1,2,5,10].map(v=>v*tickPow).find(v=>rawTick<=v)||10*tickPow),axisMax=Math.max(tickStep,Math.ceil(maxVal/tickStep)*tickStep,.001),ticks=[],plotW=w-left-right,rowH=Math.max(14,Math.min(26,(h-58-top-10)/Math.max(1,rows.length))),fmt=v=>{const x=Math.abs(Number(v));if(!Number.isFinite(x)||x===0)return'0';if(x>=1)return Number(v).toFixed(2).replace(/\\.?0+$/,'');if(x>=.01)return Number(v).toFixed(3).replace(/\\.?0+$/,'');return Number(v).toExponential(1)},xval=v=>left+Math.max(0,Number(v||0))/axisMax*plotW;for(let tick=0;tick<=axisMax+tickStep*.25;tick+=tickStep)ticks.push(Number(tick.toPrecision(12)));ticks.forEach(t=>{const xx=xval(t);waterfallLayer.appendChild(el('line',{x1:xx,y1:top-6,x2:xx,y2:h-58,stroke:'#777','stroke-width':.8,opacity:t===0?1:.32,'shape-rendering':'crispEdges'}));waterfallLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=fmt(t)});waterfallLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,stroke:'#111','stroke-width':1.8,'shape-rendering':'crispEdges'}));waterfallLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent='Mean document-to-topic probability';rows.forEach((d,i)=>{const y=top+i*rowH,val=Number(d.theta_mean),bw=plotW*Math.max(0,val)/axisMax,label=String(d.display_label||d.comparison_label||'').replace(/_/g,' ');waterfallLayer.appendChild(el('text',{x:left-10,y:y+rowH*.62,'font-size':labelSize,'text-anchor':'end','font-weight':700,fill:'#111'})).textContent=label.length>58?label.slice(0,56)+'..':label;const r=el('rect',{x:left,y:y+3,width:Math.max(1,bw),height:Math.max(4,rowH-7),fill:color(d.comparison_label||label),opacity:.9});r.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,label+'\\nmean document-to-topic probability: '+val.toFixed(4)+'\\ndocs: '+d.n_docs));r.addEventListener('mouseout',()=>wfTooltip.style.display='none');waterfallLayer.appendChild(r);const valueText=el('text',{x:Math.min(w-right,left+bw+6),y:y+rowH*.62,'font-size':smallSize,'font-weight':700,fill:'#334155'});if(left+bw+58>w-right)valueText.setAttribute('text-anchor','end');valueText.textContent=val.toFixed(3);waterfallLayer.appendChild(valueText)});document.getElementById('waterfallStats').textContent=rows.length+' of '+allRows.length+' groups'}",
    "function niceTicks(maxVal,n){const out=[];if(!Number.isFinite(maxVal)||maxVal<=0)return[0,1];const raw=maxVal/Math.max(1,n),pow=Math.pow(10,Math.floor(Math.log10(raw))),step=([1,2,5,10].map(v=>v*pow).find(v=>raw<=v)||10*pow);for(let x=0;x<=maxVal+step*.5;x+=step)out.push(x);return out}",
    "function fillPathwayScopes(){pathwayScopeSelect.replaceChildren();const o=document.createElement('option');o.value='__overall__';o.textContent='Overall topic pathways';pathwayScopeSelect.appendChild(o);const seen=new Map();PATHWAYS_CONDITION.forEach(d=>{const id=String(d.comparison_label||d.comparison_id||'');if(!id||seen.has(id))return;seen.set(id,String(d.display_label||d.comparison_label||d.comparison_id||id))});Array.from(seen.entries()).map(([id,label])=>({id,label})).sort((a,b)=>a.label.localeCompare(b.label,undefined,{sensitivity:'base'})||a.id.localeCompare(b.id)).forEach(d=>{const opt=document.createElement('option');opt.value=d.id;opt.textContent=d.label;pathwayScopeSelect.appendChild(opt)})}",
    "function selectedPathways(){const scope=pathwayScopeSelect.value;if(scope&&scope!=='__overall__'&&PATHWAYS_CONDITION.length)return PATHWAYS_CONDITION.filter(d=>String(d.comparison_label||d.comparison_id||'')===scope);return PATHWAYS_TOPIC.length?PATHWAYS_TOPIC:PATHWAYS_CONDITION}",
    "function wrapWords(s,maxChars,maxLines){const words=String(s||'').split(/\\s+/),lines=[];let line='';words.forEach(w=>{const next=line?line+' '+w:w;if(next.length>maxChars&&line){lines.push(line);line=w}else line=next});if(line)lines.push(line);if(lines.length>maxLines){lines.length=maxLines;lines[maxLines-1]=lines[maxLines-1].slice(0,Math.max(0,maxChars-2))+'..'}return lines}",
    "function addWrappedLabel(layer,text,x,y,maxChars,maxLines,lineHeight,fontSize,cls){const t=el('text',{x:x,y:y,'text-anchor':'end',class:cls||'pathLabel','font-size':fontSize});wrapWords(text,maxChars,maxLines).forEach((line,i)=>{const sp=el('tspan',{x:x,dy:i===0?0:lineHeight});sp.textContent=line;t.appendChild(sp)});layer.appendChild(t);return t}",
    "function drawPathways(){pathLayer.replaceChildren();const activePathways=selectedPathways(),tn=topicNum(topicSelect.value),allRows=activePathways.filter(d=>+d.topic_num===tn||+d.topic===tn),rows=allRows.slice().sort((a,b)=>(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||''))).slice(0,30),topicMap=new Map();activePathways.forEach(d=>{const k=String(d.pathway||'');if(!k)return;if(!topicMap.has(k))topicMap.set(k,new Set());topicMap.get(k).add(+d.topic_num||+d.topic)});const scope=pathwayScopeSelect.value==='__overall__'?'overall topic pathways':'selected condition';document.getElementById('pathStats').textContent=rows.length+' of '+allRows.length+' rows for Topic '+tn+' in '+scope;const w=1500,h=980,left=770,right=56,top=104,bottom=82,axisY=54,labelY=44,titleY=24,rowH=Math.min(30,(h-top-bottom-42)/Math.max(1,rows.length)),barH=Math.max(10,rowH*.5),totals=rows.map(d=>Math.max(+d.gene_total_universe||0,+d.gene_in||0));if(!allRows.length){const t=el('text',{x:40,y:80,class:'label'});t.textContent='No pathway rows are available for this topic and pathway scope.';pathLayer.appendChild(t);return}const ticks=niceTicks(Math.max(...totals,1),8),axisMax=Math.max(...ticks,...totals,1),x=v=>left+Math.max(0,+v||0)/axisMax*(w-left-right);pathLayer.appendChild(el('rect',{x:left,y:h-34,width:22,height:12,class:'barIn'}));pathLayer.appendChild(el('text',{x:left+30,y:h-22,class:'pathLegendText'})).textContent='topic genes';pathLayer.appendChild(el('rect',{x:left+190,y:h-34,width:22,height:12,class:'barOut'}));pathLayer.appendChild(el('text',{x:left+220,y:h-22,class:'pathLegendText'})).textContent='universe remainder';ticks.forEach(t=>{const xx=x(t);pathLayer.appendChild(el('line',{x1:xx,y1:top-10,x2:xx,y2:h-bottom-28,class:'pathTick',opacity:t===0?1:.35}));pathLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle',class:'pathTickText'})).textContent=String(Math.round(t))});pathLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,class:'pathAxis'}));pathLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle',class:'pathTickText'})).textContent='Genes in full document gene-universe enrichment';rows.forEach((d,i)=>{const y=top+18+i*rowH,inside=Math.max(0,+d.gene_in||0),total=Math.max(inside,+d.gene_total_universe||0),outside=Math.max(0,total-inside),cls=(topicMap.get(String(d.pathway||''))||new Set()).size===1?'pathLabel pathLabelTopicSpecific':'pathLabel';const lab=addWrappedLabel(pathLayer,d.pathway,left-12,y+barH*.72,76,rowH<27?1:2,14,16,cls);const inW=Math.max(inside>0?3:0,x(inside)-left),outW=Math.max(outside>0?2:0,x(total)-x(inside)),rin=el('rect',{x:left,y:y,width:inW,height:barH,class:'barIn',rx:1}),rout=el('rect',{x:left+inW,y:y,width:outW,height:barH,class:'barOut',rx:1}),marker=el('line',{x1:x(inside),y1:y-1,x2:x(inside),y2:y+barH+1,stroke:'#7f1d1d','stroke-width':1.2}),msg=String(d.pathway)+'\\nGene in topic: '+inside+'\\nGene out of topic: '+outside+'\\nUniverse hits: '+total+'\\nadj-p: '+(Number.isFinite(Number(d.padj))?Number(d.padj).toExponential(2):'NA')+'\\nGenes: '+String(d.genes||'');[rin,rout,marker,lab].forEach(node=>{node.style.cursor='pointer';node.addEventListener('click',()=>openSubgrn(d));node.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));node.addEventListener('mouseout',()=>pathTooltip.style.display='none')});pathLayer.appendChild(rin);pathLayer.appendChild(rout);pathLayer.appendChild(marker)})}",
    "function drawPathways(){pathLayer.replaceChildren();const activePathways=selectedPathways(),tn=topicNum(topicSelect.value),allRows=activePathways.filter(d=>+d.topic_num===tn||+d.topic===tn),baseRows=allRows.slice().sort((a,b)=>(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||''))).slice(0,30);ensureTfPathwayCounts(baseRows);let rows=baseRows;if(hasSelectedTf()){rows=baseRows.slice().sort((a,b)=>((tfPathwayCounts.get(tfCountKey(b))||{n:0}).n-(tfPathwayCounts.get(tfCountKey(a))||{n:0}).n)||((tfPathwayCounts.get(tfCountKey(b))||{tfN:0}).tfN-(tfPathwayCounts.get(tfCountKey(a))||{tfN:0}).tfN)||(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||'')))}const topicMap=new Map();activePathways.forEach(d=>{const k=String(d.pathway||'');if(!k)return;if(!topicMap.has(k))topicMap.set(k,new Set());topicMap.get(k).add(+d.topic_num||+d.topic)});const scope=pathwayScopeSelect.value==='__overall__'?'overall topic pathways':'selected condition';document.getElementById('pathStats').textContent=rows.length+' of '+allRows.length+' rows for Topic '+tn+' in '+scope+(hasSelectedTf()?' | selected TF target genes and per-TF counts shown':'');const w=1500,h=980,left=770,right=110,top=104,bottom=82,axisY=54,labelY=44,titleY=24,rowH=Math.min(30,(h-top-bottom-42)/Math.max(1,rows.length)),barH=Math.max(10,rowH*.5),totals=rows.map(d=>Math.max(+d.gene_total_universe||0,+d.gene_in||0));if(!allRows.length){const t=el('text',{x:40,y:80,class:'label'});t.textContent='No pathway rows are available for this topic and pathway scope.';pathLayer.appendChild(t);return}const ticks=niceTicks(Math.max(...totals,1),8),axisMax=Math.max(...ticks,...totals,1),x=v=>left+Math.max(0,+v||0)/axisMax*(w-left-right);pathLayer.appendChild(el('rect',{x:left,y:h-34,width:22,height:12,class:'barIn'}));pathLayer.appendChild(el('text',{x:left+30,y:h-22,class:'pathLegendText'})).textContent='topic genes';pathLayer.appendChild(el('rect',{x:left+190,y:h-34,width:22,height:12,class:'barOut'}));pathLayer.appendChild(el('text',{x:left+220,y:h-22,class:'pathLegendText'})).textContent='universe remainder';if(hasSelectedTf())pathLayer.appendChild(el('text',{x:w-8,y:h-22,'text-anchor':'end',class:'pathLegendText'})).textContent=selectedTfUppers.length>1?'TF genes/TFs':'TF targets';ticks.forEach(t=>{const xx=x(t);pathLayer.appendChild(el('line',{x1:xx,y1:top-10,x2:xx,y2:h-bottom-28,class:'pathTick',opacity:t===0?1:.35}));pathLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle',class:'pathTickText'})).textContent=String(Math.round(t))});pathLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,class:'pathAxis'}));pathLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle',class:'pathTickText'})).textContent='Genes in full document gene-universe enrichment';rows.forEach((d,i)=>{const y=top+18+i*rowH,inside=Math.max(0,+d.gene_in||0),total=Math.max(inside,+d.gene_total_universe||0),outside=Math.max(0,total-inside),cls=(topicMap.get(String(d.pathway||''))||new Set()).size===1?'pathLabel pathLabelTopicSpecific':'pathLabel',support=hasSelectedTf()?(tfPathwayCounts.get(tfCountKey(d))||{n:0,genes:[],tfN:0,edgeN:0}):null;const lab=addWrappedLabel(pathLayer,d.pathway,left-12,y+barH*.72,76,rowH<27?1:2,14,16,cls);const inW=Math.max(inside>0?3:0,x(inside)-left),outW=Math.max(outside>0?2:0,x(total)-x(inside)),rin=el('rect',{x:left,y:y,width:inW,height:barH,class:'barIn',rx:1}),rout=el('rect',{x:left+inW,y:y,width:outW,height:barH,class:'barOut',rx:1}),marker=el('line',{x1:x(inside),y1:y-1,x2:x(inside),y2:y+barH+1,stroke:'#7f1d1d','stroke-width':1.2}),extra=support?'\\nSelected TFs: '+selectedTfLabel()+'\\nUnique selected-TF target genes in pathway: '+support.n+'\\nSelected TFs with hits: '+(Number(support.tfN)||0)+' / '+selectedTfUppers.length+'\\nTF-target pairs represented: '+(Number(support.edgeN)||0)+(support.genes&&support.genes.length?'\\nPer-TF target genes: '+support.genes.join('; '):''):'',msg=String(d.pathway)+'\\nGene in topic: '+inside+'\\nGene out of topic: '+outside+'\\nUniverse hits: '+total+'\\nadj-p: '+(Number.isFinite(Number(d.padj))?Number(d.padj).toExponential(2):'NA')+extra+'\\nGenes: '+String(d.genes||'');[rin,rout,marker,lab].forEach(node=>{node.style.cursor='pointer';node.addEventListener('click',()=>openSubgrn(d));node.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));node.addEventListener('mouseout',()=>pathTooltip.style.display='none')});pathLayer.appendChild(rin);pathLayer.appendChild(rout);pathLayer.appendChild(marker);if(support){const n=Number(support.n)||0,tfN=Number(support.tfN)||0,badgeText=selectedTfUppers.length>1?String(n)+'/'+String(tfN):String(n);if(n>0){const badge=el('rect',{x:w-72,y:y-2,width:64,height:barH+4,rx:8,fill:'#dc2626',stroke:'#7f1d1d','stroke-width':1.4,opacity:.95,style:'cursor:pointer'});badge.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));badge.addEventListener('mouseout',()=>pathTooltip.style.display='none');badge.addEventListener('click',()=>openSubgrn(d));pathLayer.appendChild(badge);rin.setAttribute('stroke','#dc2626');rin.setAttribute('stroke-width',2.4);rout.setAttribute('stroke','#dc2626');rout.setAttribute('stroke-width',1.6)}const ct=el('text',{x:w-40,y:y+barH*.78,'text-anchor':'middle','font-size':15,'font-weight':700,fill:n>0?'#fff':'#94a3b8',style:'cursor:pointer'});ct.textContent=badgeText;ct.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));ct.addEventListener('mouseout',()=>pathTooltip.style.display='none');ct.addEventListener('click',()=>openSubgrn(d));pathLayer.appendChild(ct)}})}",
    .m3tb_subgrn_js(),
    "function draw(){updateSelectedContext();drawTopic();drawTfWaterfall();drawGroupWaterfall();drawPathways()}init();",
    "</script></body></html>"
  )
  html <- .m3tb_polish_topic_report_html(html)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_html, useBytes = TRUE)
  out_html
}

.m3tb_condition_mds_svg_path <- function(out_html) {
  stem <- sub("[.]html$", "", basename(out_html))
  stem <- sub("^condition_topic_report_", "condition_mds_", stem)
  file.path(dirname(out_html), "assets", paste0(stem, ".svg"))
}

.m3tb_write_condition_mds_svg <- function(group_mds, out_svg) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .log_abort("Package ggplot2 is required to render Module 3 condition MDS SVG reports.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    .log_abort("Package ggrepel is required to render Module 3 condition MDS SVG reports.")
  }
  dt <- data.table::copy(data.table::as.data.table(group_mds))
  if (!nrow(dt)) return(invisible(FALSE))
  if (!"mds_label" %in% names(dt)) dt[, mds_label := .m3tb_short_label(display_label, 18L)]
  if (!"color" %in% names(dt)) dt[, color := .m3tb_group_color(display_label)]
  if (!"shape_value" %in% names(dt)) {
    dt[, shape_value := data.table::fifelse(grepl("Down|Target-Down", paste(display_label, comparison_label), ignore.case = TRUE), 25L, 16L)]
  }
  dt[, `:=`(
    plot_label = data.table::fifelse(!is.na(mds_label) & nzchar(mds_label), mds_label, display_label),
    color_value = data.table::fifelse(grepl("^#[0-9A-Fa-f]{6}$", color), color, "#2A9D8F"),
    shape_value = data.table::fifelse(as.integer(shape_value) == 25L, 25L, 16L)
  )]
  is_comparison <- any(dt$doc_design == "comparison" | dt$shape_value == 25L, na.rm = TRUE)
  if (!"panel_label" %in% names(dt)) {
    if ("method_setup" %in% names(dt)) {
      dt[, panel_label := method_setup]
    } else {
      dt[, panel_label := "Condition"]
    }
  }
  dt[, base_panel_label := as.character(panel_label)]
  if (is_comparison) {
    dt[, direction_label := data.table::fifelse(shape_value == 25L, "Down", "Up")]
    dt[, split_panel_label := paste(base_panel_label, direction_label, sep = "\n")]
  } else {
    dt[, split_panel_label := base_panel_label]
  }
  if ("method_order" %in% names(dt)) {
    panel_order <- unique(dt[order(method_order, base_panel_label, data.table::fifelse(shape_value == 25L, 2L, 1L))]$split_panel_label)
  } else {
    panel_order <- unique(dt[order(base_panel_label, data.table::fifelse(shape_value == 25L, 2L, 1L))]$split_panel_label)
  }
  dt[, split_panel_label := factor(split_panel_label, levels = panel_order)]
  dt[, label_mode := "repel"]
  dt[, scaled_x := {
    r <- range(MDS1, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) 0.5 else (MDS1 - r[1]) / diff(r)
  }, by = split_panel_label]
  dt[, scaled_y := {
    r <- range(MDS2, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) 0.5 else (MDS2 - r[1]) / diff(r)
  }, by = split_panel_label]
  dt[, nearest_scaled_dist := {
    n <- .N
    if (n <= 1L) {
      Inf
    } else {
      vapply(seq_len(n), function(i) {
        min(sqrt((scaled_x[i] - scaled_x[-i])^2 + (scaled_y[i] - scaled_y[-i])^2), na.rm = TRUE)
      }, numeric(1L))
    }
  }, by = split_panel_label]
  max_panel_n <- max(data.table::frank(dt$split_panel_label, ties.method = "dense"), 1L)
  max_panel_n <- max(dt[, .N, by = split_panel_label]$N, na.rm = TRUE)
  dense_labels <- !is_comparison && is.finite(max_panel_n) && max_panel_n > 12L
  if (dense_labels) {
    dt[, plot_label := .m3tb_short_label(plot_label, 15L)]
  }
  point_size <- if (dense_labels) 5.4 else 6.4
  many_panels <- length(panel_order) > 4L
  repel_size <- if (is_comparison && many_panels) 2.75 else if (is_comparison) 3.2 else if (dense_labels) 3.85 else 4.75
  center_size <- if (is_comparison) 3.5 else if (dense_labels) 4.15 else 5.2
  repel_box_padding <- if (is_comparison) 0.68 else if (dense_labels) 0.82 else 1.35
  repel_point_padding <- if (is_comparison) 0.25 else if (dense_labels) 0.34 else 0.5
  repel_force <- if (is_comparison && many_panels) 8 else if (is_comparison) 44 else if (dense_labels) 38 else 24
  repel_max_iter <- if (many_panels) 8000L else 90000L
  repel_max_time <- if (many_panels) 0.8 else if (is_comparison || dense_labels) 10 else 6
  repel_dt <- dt[label_mode == "repel"]
  center_dt <- dt[label_mode == "center"]
  p <- ggplot2::ggplot(dt, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = color_value, fill = color_value, shape = shape_value),
      size = point_size,
      alpha = 0.9,
      stroke = 0.8
    ) +
    ggplot2::geom_text(
      data = center_dt,
      ggplot2::aes(label = plot_label, color = color_value),
      fontface = "bold",
      size = center_size,
      show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = repel_dt,
      ggplot2::aes(label = plot_label, color = color_value),
      fontface = "bold",
      size = repel_size,
      min.segment.length = 0,
      box.padding = repel_box_padding,
      point.padding = repel_point_padding,
      force = repel_force,
      force_pull = 0.55,
      max.iter = repel_max_iter,
      max.time = repel_max_time,
      max.overlaps = Inf,
      segment.color = "grey35",
      segment.size = 0.55,
      label.padding = grid::unit(0.16, "lines"),
      label.r = grid::unit(0.08, "lines"),
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.86),
      seed = 1
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.42)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.42)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "MDS1", y = "MDS2") +
    ggplot2::theme_classic(base_size = 14, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(face = "bold", size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 16),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey72", linewidth = 0.35),
      panel.spacing = grid::unit(if (is_comparison) 18 else 12, "pt"),
      plot.margin = ggplot2::margin(14, 36, 14, 36),
      legend.position = "none",
      aspect.ratio = if (is_comparison) 0.62 else 0.62
    )
  if (is_comparison) p <- p + ggplot2::facet_wrap(~ split_panel_label, ncol = 2, scales = "free")
  dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
  grDevices::svg(out_svg, width = if (is_comparison) 12 else 9.8, height = if (is_comparison) 7 else 5.4, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  invisible(TRUE)
}

.m3tb_condition_mds_plot <- function(group_mds, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .log_abort("Package ggplot2 is required to render Module 3 condition MDS reports.")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    .log_abort("Package ggrepel is required to render Module 3 condition MDS reports.")
  }
  dt <- data.table::copy(data.table::as.data.table(group_mds))
  if (!nrow(dt)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No MDS rows"))
  }
  if (!"group_label" %in% names(dt)) {
    if ("comparison_label" %in% names(dt)) {
      dt[, group_label := comparison_label]
    } else if ("display_label" %in% names(dt)) {
      dt[, group_label := display_label]
    } else {
      dt[, group_label := paste0("Group", seq_len(.N))]
    }
  }
  if (!"display_label" %in% names(dt)) dt[, display_label := group_label]
  if (!"comparison_label" %in% names(dt)) dt[, comparison_label := group_label]
  if (!"mds_label" %in% names(dt)) dt[, mds_label := .m3tb_short_label(display_label, 18L)]
  if (!"color" %in% names(dt)) dt[, color := .m3tb_group_color(display_label)]
  if (!"doc_design" %in% names(dt)) dt[, doc_design := "condition"]
  if (!"shape_value" %in% names(dt)) {
    dt[, shape_value := data.table::fifelse(grepl("Down|Target-Down", paste(display_label, comparison_label), ignore.case = TRUE), 25L, 16L)]
  }
  dt[, `:=`(
    plot_label = data.table::fifelse(!is.na(mds_label) & nzchar(mds_label), mds_label, display_label),
    color_value = data.table::fifelse(grepl("^#[0-9A-Fa-f]{6}$", color), color, "#2A9D8F"),
    shape_value = data.table::fifelse(as.integer(shape_value) == 25L, 25L, 16L)
  )]
  is_comparison <- any(dt$doc_design == "comparison" | dt$shape_value == 25L, na.rm = TRUE)
  if (!"panel_label" %in% names(dt)) {
    if ("method_setup" %in% names(dt)) {
      dt[, panel_label := method_setup]
    } else {
      dt[, panel_label := "Condition"]
    }
  }
  dt[, base_panel_label := as.character(panel_label)]
  if (is_comparison) {
    dt[, direction_label := data.table::fifelse(shape_value == 25L, "Down", "Up")]
    dt[, split_panel_label := paste(base_panel_label, direction_label, sep = "\n")]
  } else {
    dt[, split_panel_label := base_panel_label]
  }
  if ("method_order" %in% names(dt)) {
    panel_order <- unique(dt[order(method_order, base_panel_label, data.table::fifelse(shape_value == 25L, 2L, 1L))]$split_panel_label)
  } else {
    panel_order <- unique(dt[order(base_panel_label, data.table::fifelse(shape_value == 25L, 2L, 1L))]$split_panel_label)
  }
  dt[, split_panel_label := factor(split_panel_label, levels = panel_order)]
  dt[, label_mode := "repel"]
  dt[, scaled_x := {
    r <- range(MDS1, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) 0.5 else (MDS1 - r[1]) / diff(r)
  }, by = split_panel_label]
  dt[, scaled_y := {
    r <- range(MDS2, na.rm = TRUE)
    if (!all(is.finite(r)) || diff(r) == 0) 0.5 else (MDS2 - r[1]) / diff(r)
  }, by = split_panel_label]
  max_panel_n <- max(dt[, .N, by = split_panel_label]$N, na.rm = TRUE)
  dense_labels <- !is_comparison && is.finite(max_panel_n) && max_panel_n > 12L
  if (dense_labels) {
    dt[, plot_label := .m3tb_short_label(plot_label, 15L)]
  }
  many_panels <- length(panel_order) > 4L
  point_size <- if (many_panels) 3.4 else if (dense_labels) 5.4 else 6.4
  repel_size <- if (is_comparison && many_panels) 2.35 else if (is_comparison) 3.2 else if (dense_labels) 3.85 else 4.75
  repel_box_padding <- if (many_panels) 0.58 else if (is_comparison) 0.68 else if (dense_labels) 0.82 else 1.35
  repel_point_padding <- if (many_panels) 0.42 else if (is_comparison) 0.25 else if (dense_labels) 0.34 else 0.5
  repel_force <- if (many_panels) 5 else if (is_comparison) 44 else if (dense_labels) 38 else 24
  p <- ggplot2::ggplot(dt, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_point(
      ggplot2::aes(color = color_value, fill = color_value, shape = shape_value),
      size = point_size,
      alpha = 0.9,
      stroke = 0.8
    )
  if (many_panels) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = plot_label, color = color_value),
      fontface = "bold",
      size = repel_size,
      min.segment.length = 0,
      box.padding = repel_box_padding,
      point.padding = repel_point_padding,
      force = repel_force,
      force_pull = 0.02,
      max.iter = 20000,
      max.time = 2,
      max.overlaps = Inf,
      segment.color = "grey35",
      segment.size = 0.55,
      seed = 1
    )
  } else {
    p <- p + ggrepel::geom_label_repel(
      ggplot2::aes(label = plot_label, color = color_value),
      fontface = "bold",
      size = repel_size,
      min.segment.length = 0,
      box.padding = repel_box_padding,
      point.padding = repel_point_padding,
      force = repel_force,
      force_pull = 0.55,
      max.iter = 90000,
      max.time = if (is_comparison || dense_labels) 10 else 6,
      max.overlaps = Inf,
      segment.color = "grey35",
      segment.size = 0.55,
      label.padding = grid::unit(0.16, "lines"),
      label.r = grid::unit(0.08, "lines"),
      label.size = 0.25,
      fill = grDevices::adjustcolor("white", alpha.f = 0.86),
      seed = 1
    )
  }
  p <- p +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.42)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.42)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = "MDS1", y = "MDS2") +
    ggplot2::theme_classic(base_size = 14, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(face = "bold", size = 12),
      strip.text = ggplot2::element_text(face = "bold", size = 16),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey72", linewidth = 0.35),
      panel.spacing = grid::unit(if (is_comparison) 10 else 12, "pt"),
      plot.margin = ggplot2::margin(6, 16, 6, 16),
      legend.position = "none",
      aspect.ratio = if (is_comparison) 1.05 else 0.62
    )
  if (is_comparison || length(panel_order) > 1L) {
    p <- p + ggplot2::facet_wrap(~ split_panel_label, ncol = if (many_panels) 6L else 2L, scales = "free")
  }
  p
}

.m3tb_condition_report_html <- function(title,
                                        group_mds,
                                        group_topic,
                                        tf_topic = NULL,
                                        pathways,
                                        out_html,
                                        mds_svg_src = NULL,
                                        subgrn_manifest = NULL,
                                        subgrn_payload_base = "../../pathway_subgrn_payloads",
                                        condition_payload = NULL,
                                        report_state = list()) {
  if (!is.null(condition_payload) &&
      nzchar(as.character(condition_payload$payload_file %||% ""))) {
    return(.m3cr_condition_report_html(
      title = title,
      group_mds = group_mds,
      group_topic = group_topic,
      tf_topic = tf_topic,
      pathways = pathways,
      out_html = out_html,
      condition_payload = condition_payload,
      report_state = report_state
    ))
  }
  group_json <- .m3tb_json_for_html(group_mds)
  topic_json <- .m3tb_json_for_html(group_topic)
  if (is.null(tf_topic)) tf_topic <- data.table::data.table()
  tf_topic_json <- .m3tb_json_for_html(tf_topic)
  pathway_json <- .m3tb_json_for_html(pathways)
  if (is.null(subgrn_manifest)) subgrn_manifest <- .m3tb_subgrn_empty_manifest()
  subgrn_json <- .m3tb_json_for_html(subgrn_manifest)
  mds_img <- "<div class=\"tooltip\" id=\"mdsTooltip\"></div><svg id=\"mdsSvg\" viewBox=\"0 0 980 680\"><g id=\"mdsLayer\"></g></svg>"
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    "<meta name=\"craftgrn-module3-report-schema\" content=\"2\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    paste0("<style>html,body{width:100%;height:100%;overflow:hidden}body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif;font-weight:700}.top{height:48px;background:#fff;border-bottom:1px solid #d6d6d0;display:flex;justify-content:space-between;gap:10px;align-items:center;padding:6px 12px;box-sizing:border-box}h1{font-size:18px;line-height:1.05;margin:0}.meta{font-size:11px;color:#475569;line-height:1.1}.controls{display:flex;gap:8px;align-items:center;white-space:nowrap}select,input,button{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #aaa;border-radius:3px;background:#fff;padding:5px 7px}button{background:#111;color:#fff}.grid{height:calc(100vh - 48px);display:grid;grid-template-columns:minmax(520px,1fr) minmax(700px,1.35fr);gap:8px;padding:8px;box-sizing:border-box}.left{display:grid;grid-template-rows:minmax(0,1fr) minmax(0,1fr);gap:8px;min-height:0}.pane{background:#fff;border:1px solid #d6d6d0;display:flex;flex-direction:column;min-height:0;overflow:hidden}.paneHead{padding:8px 10px;border-bottom:1px solid #e5e5df;display:flex;justify-content:space-between;gap:8px;align-items:center}.pane h2{font-size:20px;line-height:1.05;margin:0}.barControls{display:flex;gap:8px;align-items:center;white-space:nowrap}.barControls select{max-width:180px}.barControls input{width:130px}.barControls select[multiple]{height:46px;min-width:130px}.body{position:relative;flex:1;min-height:0}.splitWaterfall{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);height:100%;min-height:0}.splitSide{position:relative;min-width:0;min-height:0}.splitSide:first-child{border-right:1px solid #e5e5df}.clearTf{display:none;margin-left:8px;padding:3px 6px;font-size:11px}.note{font-size:10px;color:#555;border-top:1px solid #e5e5df;padding:5px 8px;line-height:1.1}.tooltip{position:absolute;display:none;background:rgba(17,17,17,0.92);color:#fff;font:700 12px Arial,Helvetica,sans-serif;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:380px;line-height:1.35;z-index:5}body.embed .top{display:none}body.embed .grid{height:100vh}svg{width:100%;height:100%;display:block}.label{font-size:12px;font-weight:700}.small{font-size:10px}.mdsLeader{stroke:#555;stroke-width:1.7;opacity:.75}.mdsPointHit{cursor:pointer;fill:transparent;stroke:transparent;pointer-events:all}.mdsLabel{paint-order:stroke;stroke:#fff;stroke-width:6px;stroke-linejoin:round;cursor:pointer}.barIn{fill:#cc454b}.barOut{fill:#a9cfe5}.pathAxis{stroke:#111;stroke-width:1.6;shape-rendering:crispEdges}.pathTick{stroke:#777;stroke-width:.9;shape-rendering:crispEdges}.pathLabel{font-size:16px;font-weight:700;fill:#111}.pathLabelTopicSpecific{fill:#cc2f3c}.pathLabelGroupSpecific{fill:#2563eb}.pathLabelBothSpecific{fill:#7e22ce}.pathTickText{font-size:14px;font-weight:700;fill:#111}.pathLegendText{font-size:14px;font-weight:700;fill:#334155}", .m3tb_subgrn_css(), "@media(max-width:1100px){.grid{grid-template-columns:1fr 1fr}.top{height:52px}.grid{height:calc(100vh - 52px)}}</style>"),
    "</head><body><div class=\"top\">",
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    "<div class=\"controls\"><label>Condition/comparison <select id=\"groupSelect\"></select></label><label>Topic <select id=\"topicSelect\"></select></label><button id=\"exportSvgButton\">Export SVG</button></div></div>",
    paste0("<main class=\"grid\"><div class=\"left\"><section class=\"pane\"><div class=\"paneHead\"><h2>Condition/Comparison MDS</h2><span id=\"mdsStats\" class=\"meta\"></span></div><div class=\"body\">", mds_img, "</div><div class=\"note\">Dots are condition/comparison mean document-to-topic probability profiles (theta) using Jensen-Shannon distance.</div></section><section class=\"pane\"><div class=\"paneHead\"><h2>Topic Bar Plots</h2><div class=\"barControls\"><label>TF <input id=\"tfSearchInput\" list=\"tfOptionList\" placeholder=\"Type TF\"/><datalist id=\"tfOptionList\"></datalist></label><label>Selected <select id=\"tfSelect\" multiple size=\"3\"></select></label><span id=\"waterfallStats\" class=\"meta\"></span><button id=\"clearTfButton\" class=\"clearTf\" type=\"button\">Clear TF</button></div></div><div class=\"body\"><div class=\"tooltip\" id=\"wfTooltip\"></div><div class=\"splitWaterfall\"><div class=\"splitSide\"><svg id=\"tfWaterfallSvg\" viewBox=\"0 0 760 760\"><g id=\"tfWaterfallLayer\"></g></svg></div><div class=\"splitSide\"><svg id=\"waterfallSvg\" viewBox=\"0 0 760 760\"><g id=\"waterfallLayer\"></g></svg></div></div></div><div class=\"note\">Left: TF theta for the selected topic. Right: topics for the selected condition; selecting TFs ranks topics by their mean TF theta.</div></section></div><section class=\"pane\"><div class=\"paneHead\"><h2>Pathways</h2><span id=\"pathStats\" class=\"meta\"></span></div><div class=\"body\"><div class=\"tooltip\" id=\"pathTooltip\"></div><svg id=\"pathSvg\" viewBox=\"0 0 1500 980\"><g id=\"pathLayer\"></g></svg>", .m3tb_subgrn_html_panel(), "</div><div class=\"note\">Full per-condition or per-comparison pathway table is embedded for the selected method and K; the plot displays the top 30 rows for the selected group and topic ranked by adjusted p-value. Pathway name colors: red = topic-specific, blue = condition/comparison-specific, purple = both. Click a pathway bar to open the sub-GRN browser.</div></section></main>"),
    "<script>",
    paste0("const GROUP_MDS=", group_json, ";"),
    paste0("const GROUP_TOPIC=", topic_json, ";"),
    paste0("const TF_TOPIC=", tf_topic_json, ";"),
    paste0("const PATHWAYS=", pathway_json, ";"),
    paste0("const SUBGRN_MANIFEST=", subgrn_json, ";"),
    paste0("const SUBGRN_PAYLOAD_BASE=", jsonlite::toJSON(subgrn_payload_base, auto_unbox = TRUE), ";"),
    "if(new URLSearchParams(location.search).get('embed')==='1')document.body.classList.add('embed');",
    "const PAL=['#4E79A7','#F28E2B','#59A14F','#E15759','#B07AA1','#76B7B2','#EDC948','#FF9DA7','#9C755F','#BAB0AC','#1F77B4','#D62728'];const groupSelect=document.getElementById('groupSelect'),topicSelect=document.getElementById('topicSelect'),tfSelect=document.getElementById('tfSelect'),tfSearchInput=document.getElementById('tfSearchInput'),tfOptionList=document.getElementById('tfOptionList'),mdsSvg=document.getElementById('mdsSvg'),mdsExportSvg=null,tfWaterfallSvg=document.getElementById('tfWaterfallSvg'),waterfallSvg=document.getElementById('waterfallSvg'),pathSvg=document.getElementById('pathSvg'),mdsLayer=document.getElementById('mdsLayer'),tfWaterfallLayer=document.getElementById('tfWaterfallLayer'),waterfallLayer=document.getElementById('waterfallLayer'),pathLayer=document.getElementById('pathLayer'),mdsTooltip=document.getElementById('mdsTooltip'),wfTooltip=document.getElementById('wfTooltip'),pathTooltip=document.getElementById('pathTooltip'),clearTfButton=document.getElementById('clearTfButton');let selectedTfUppers=[],selectedTfUpper='',tfLabelMap=new Map();const tfPathwayCounts=new Map(),tfPathwayPending=new Set();function esc(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;')}function el(n,a){const x=document.createElementNS('http://www.w3.org/2000/svg',n);Object.entries(a||{}).forEach(([k,v])=>x.setAttribute(k,v));return x}function sc(vals,lo,hi){const xs=vals.map(Number).filter(Number.isFinite),a=Math.min(...xs),b=Math.max(...xs);return v=>!Number.isFinite(Number(v))||a===b?(lo+hi)/2:lo+(Number(v)-a)/(b-a)*(hi-lo)}function tooltip(tt,ev,text){tt.innerHTML=esc(text).replace(/\\n/g,'<br/>');tt.style.display='block';tt.style.left=(ev.offsetX+12)+'px';tt.style.top=(ev.offsetY+12)+'px'}function svgFontSize(svg,targetPx){if(!svg)return String(targetPx);const m=svg.getScreenCTM();const s=m?Math.sqrt(Math.abs(m.a*m.d)):1;return Math.max(1,targetPx/Math.max(s,0.001)).toFixed(2)}function color(i){return PAL[i%PAL.length]}function groupColor(d,i){const c=String(d.color||'');return /^#[0-9A-Fa-f]{6}$/.test(c)?c:color(i)}function topicNum(t){return Number(String(t).replace(/^Topic/,''))}function isDownGroup(d){const s=String(d.comparison_label||'')+' '+String(d.display_label||'')+' '+String(d.group_label||'');return Number(d.shape_value)===25||/Target-Down|\\bDown\\b/i.test(s)}",
    "function selectGroup(id){groupSelect.value=id;fillTopics();draw()}",
    "function selectedGroupRows(){return GROUP_TOPIC.filter(d=>d.comparison_label===groupSelect.value).sort((a,b)=>Number(b.theta_mean)-Number(a.theta_mean)||Number(a.topic_num)-Number(b.topic_num))}",
    "function topicHasPathwayForGroup(topic,grp){const tn=topicNum(topic),hasGroupPath=PATHWAYS.some(d=>String(d.comparison_label||'')!=='');return PATHWAYS.some(d=>(+d.topic_num===tn||+d.topic===tn)&&(!hasGroupPath||String(d.comparison_label||'')===grp))}",
    "function selectedGroupLabel(){const opt=groupSelect.options[groupSelect.selectedIndex];return opt?opt.textContent:groupSelect.value}function updateSelectedContext(){const x=document.getElementById('selectedContext');if(x)x.textContent='Selected: '+selectedGroupLabel()+' | Topic '+topicNum(topicSelect.value)}",
    "function exportStyle(){return 'text{font-family:Arial,Helvetica,sans-serif;font-weight:700}.label{font-size:12px;font-weight:700}.small{font-size:10px}.barIn{fill:#cc454b}.barOut{fill:#a9cfe5}.pathAxis{stroke:#111;stroke-width:1.6;shape-rendering:crispEdges}.pathTick{stroke:#777;stroke-width:.9;shape-rendering:crispEdges}.pathLabel{font-size:16px;font-weight:700;fill:#111}.pathLabelTopicSpecific{fill:#cc2f3c}.pathLabelGroupSpecific{fill:#2563eb}.pathLabelBothSpecific{fill:#7e22ce}.pathTickText{font-size:14px;font-weight:700;fill:#111}.pathLegendText{font-size:14px;font-weight:700;fill:#334155}.mdsLeader{stroke:#555;stroke-width:1.1;opacity:.65}'}function exportSvgTitle(root,text,x,y,size){const t=el('text',{x:x,y:y,'font-size':size,'font-weight':700,fill:'#111'});t.textContent=text;root.appendChild(t)}function addExportPanel(root,src,title,x,y,w,h){exportSvgTitle(root,title,x,y+22,24);if(src&&src.tagName&&src.tagName.toLowerCase()==='svg'){const clone=src.cloneNode(true);clone.setAttribute('x',x);clone.setAttribute('y',y+38);clone.setAttribute('width',w);clone.setAttribute('height',h-38);clone.setAttribute('preserveAspectRatio','xMidYMid meet');clone.removeAttribute('style');root.appendChild(clone);return}if(src&&src.getAttribute){const img=el('image',{x:x,y:y+38,width:w,height:h-38,preserveAspectRatio:'xMidYMid meet'});img.setAttributeNS('http://www.w3.org/1999/xlink','href',src.getAttribute('src')||'');img.setAttribute('href',src.getAttribute('src')||'');root.appendChild(img)}}function exportSvg(){draw();const root=el('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 2480 1400',width:2480,height:1400});const style=el('style',{});style.textContent=exportStyle();root.appendChild(style);root.appendChild(el('rect',{x:0,y:0,width:2480,height:1400,fill:'#fff'}));exportSvgTitle(root,document.title||'Condition topic report',30,42,28);addExportPanel(root,mdsExportSvg||mdsSvg,'Condition/Comparison MDS',30,74,760,590);addExportPanel(root,waterfallSvg,'Topic Bar Plot',30,720,760,620);addExportPanel(root,pathSvg,'Pathways',830,74,1600,1266);const text=new XMLSerializer().serializeToString(root),blob=new Blob([text],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=(document.title||'condition_topic_report').replace(/[^A-Za-z0-9_.-]+/g,'_')+'.svg';document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(a.href)}",
    "function selectedTfSet(){return new Set(selectedTfUppers)}function hasSelectedTf(){return selectedTfUppers.length>0}function selectedTfKey(){return selectedTfUppers.slice().sort().join(',')}function selectedTfLabel(){return selectedTfUppers.join(', ')}function tfKeyFromInput(x){const raw=String(x||'').trim();if(!raw)return'';const up=raw.toUpperCase();if(tfLabelMap.has(up))return up;for(const [key,label] of tfLabelMap.entries()){if(String(label).toUpperCase()===up)return key}return''}function updateTfControls(){const keep=selectedTfSet();if(tfSelect)Array.from(tfSelect.options).forEach(o=>{o.selected=keep.has(o.value)});selectedTfUpper=selectedTfUppers[0]||'';if(clearTfButton)clearTfButton.style.display=hasSelectedTf()?'inline-block':'none'}function setSelectedTfs(vals){const seen=new Set(),out=[];(vals||[]).forEach(v=>{const k=tfKeyFromInput(v)||String(v||'').toUpperCase();if(k&&tfLabelMap.has(k)&&!seen.has(k)){seen.add(k);out.push(k)}});selectedTfUppers=out;updateTfControls();draw()}function setSelectedTf(tf){const k=tfKeyFromInput(tf)||String(tf||'').toUpperCase();if(!k||!tfLabelMap.has(k))return;const s=selectedTfSet();if(s.has(k))s.delete(k);else s.add(k);setSelectedTfs(Array.from(s))}function addTfFromSearch(){const k=tfKeyFromInput(tfSearchInput?tfSearchInput.value:'');if(k){const s=selectedTfSet();s.add(k);if(tfSearchInput)tfSearchInput.value='';setSelectedTfs(Array.from(s))}}function fillTfSelect(){const old=selectedTfUppers.slice(),seen=new Map();TF_TOPIC.forEach(d=>{const key=String(d.tf_upper||d.tf||'').toUpperCase();if(key&&!seen.has(key))seen.set(key,String(d.tf||key))});const rows=Array.from(seen.entries()).map(([key,label])=>({key,label})).sort((a,b)=>a.label.localeCompare(b.label,undefined,{sensitivity:'base'})||a.key.localeCompare(b.key));tfLabelMap=new Map(rows.map(d=>[d.key,d.label]));tfSelect.replaceChildren();if(tfOptionList)tfOptionList.replaceChildren();rows.forEach(d=>{const o=document.createElement('option');o.value=d.key;o.textContent=d.label;tfSelect.appendChild(o);if(tfOptionList){const q=document.createElement('option');q.value=d.label;tfOptionList.appendChild(q)}});selectedTfUppers=old.filter(k=>tfLabelMap.has(k));updateTfControls()}",
    "function init(){GROUP_MDS.slice().sort((a,b)=>String(a.display_label||a.comparison_label||'').localeCompare(String(b.display_label||b.comparison_label||''),undefined,{sensitivity:'base'})).forEach((d,i)=>{const o=document.createElement('option');o.value=d.comparison_label;o.textContent=d.display_label||d.comparison_label;groupSelect.appendChild(o)});if(groupSelect.options.length)groupSelect.value=groupSelect.options[0].value;fillTfSelect();groupSelect.addEventListener('change',()=>selectGroup(groupSelect.value));topicSelect.addEventListener('change',draw);tfSelect.addEventListener('change',()=>setSelectedTfs(Array.from(tfSelect.selectedOptions).map(o=>o.value)));if(tfSearchInput){tfSearchInput.addEventListener('change',addTfFromSearch);tfSearchInput.addEventListener('keydown',ev=>{if(ev.key==='Enter'){ev.preventDefault();addTfFromSearch()}})}if(clearTfButton)clearTfButton.addEventListener('click',()=>setSelectedTfs([]));document.getElementById('exportSvgButton').addEventListener('click',exportSvg);fillTopics();draw()}",
    "function fillTopics(){const rows=selectedGroupRows(),grp=groupSelect.value;topicSelect.replaceChildren();rows.forEach(d=>{const o=document.createElement('option');o.value=d.topic;o.textContent='Topic '+topicNum(d.topic)+' ('+(+d.theta_mean||0).toFixed(3)+')';topicSelect.appendChild(o)});const firstWithPath=rows.find(d=>topicHasPathwayForGroup(d.topic,grp));if(firstWithPath)topicSelect.value=firstWithPath.topic}",
    "function mdsLabelText(d){return String(d.mds_label||d.display_label||d.comparison_label||'').replace(/_/g,' ')}function drawTextBox(g,text,x,y,anchor,fontSize,fill){const t=el('text',{x:x,y:y,'text-anchor':anchor||'middle','font-size':fontSize||13,'font-weight':700,fill:fill||'#111',class:'mdsLabel'});t.textContent=text;g.appendChild(t);return{t:t}}function layoutLabels(items,panel){const fs=16,labels=items.map((p,i)=>{const txt=mdsLabelText(p.d);let dx=p.x-(panel.x0+panel.x1)/2,dy=p.y-(panel.y0+panel.y1)/2,len=Math.hypot(dx,dy);if(len<1){const a=(i+1)*2.399963;dx=Math.cos(a);dy=Math.sin(a);len=1}const w=Math.max(54,Math.min(148,txt.length*fs*.62+16)),h=fs+10,gap=30+(i%3)*7;return{p,txt,x:p.x+dx/len*gap,y:p.y+dy/len*gap,w,h}});function clamp(a){a.x=Math.max(panel.x0+a.w/2+4,Math.min(panel.x1-a.w/2-4,a.x));a.y=Math.max(panel.y0+a.h/2+4,Math.min(panel.y1-a.h/2-4,a.y))}labels.forEach(clamp);for(let iter=0;iter<120;iter++){for(let i=0;i<labels.length;i++)for(let j=i+1;j<labels.length;j++){const a=labels[i],b=labels[j],ox=(a.w+b.w)/2+7-Math.abs(a.x-b.x),oy=(a.h+b.h)/2+7-Math.abs(a.y-b.y);if(ox>0&&oy>0){let dx=a.x-b.x,dy=a.y-b.y;if(Math.abs(dx)+Math.abs(dy)<.1){dx=i<j?-1:1;dy=(i+j)%2?-1:1}const len=Math.max(.1,Math.hypot(dx,dy)),push=Math.min(ox,oy)*.22;a.x+=dx/len*push;a.y+=dy/len*push;b.x-=dx/len*push;b.y-=dy/len*push;clamp(a);clamp(b)}}}return labels}function drawMarker(g,p,selected){const d=p.d,c=groupColor(d,p.i),shape=isDownGroup(d)?'down':'circle';if(shape==='down'){const s=18,pts=[[p.x,p.y+s],[p.x-s,p.y-s],[p.x+s,p.y-s]].map(v=>v.join(',')).join(' ');g.appendChild(el('polygon',{points:pts,fill:c,opacity:.86,stroke:selected?'#111':'#fff','stroke-width':selected?4:1.6}))}else{g.appendChild(el('circle',{cx:p.x,cy:p.y,r:selected?15:12,fill:c,opacity:.9,stroke:selected?'#111':'#fff','stroke-width':selected?4:1.6}))}const hit=el('circle',{cx:p.x,cy:p.y,r:36,class:'mdsPointHit'});hit.addEventListener('click',()=>selectGroup(d.comparison_label));hit.addEventListener('mousemove',ev=>tooltip(mdsTooltip,ev,String(d.display_label||d.comparison_label)+'\\ndocs: '+d.n_docs));hit.addEventListener('mouseout',()=>mdsTooltip.style.display='none');g.appendChild(hit)}function drawMds(){if(!mdsLayer)return;mdsLayer.replaceChildren();const W=980,H=680,isComparison=GROUP_MDS.some(d=>String(d.doc_design||'')==='comparison'||isDownGroup(d)),panels=isComparison?[{name:'Up',rows:GROUP_MDS.filter(d=>!isDownGroup(d)),x0:38,x1:484,y0:82,y1:582},{name:'Down',rows:GROUP_MDS.filter(d=>isDownGroup(d)),x0:496,x1:960,y0:82,y1:582}]:[{name:'Condition',rows:GROUP_MDS,x0:50,x1:940,y0:72,y1:590}];panels.forEach(panel=>{if(!panel.rows.length)return;mdsLayer.appendChild(el('rect',{x:panel.x0,y:36,width:panel.x1-panel.x0,height:28,fill:'#f0f0ee',stroke:'#ddd'}));const strip=el('text',{x:(panel.x0+panel.x1)/2,y:56,'text-anchor':'middle','font-size':18,'font-weight':700,fill:'#222'});strip.textContent=panel.name;mdsLayer.appendChild(strip);mdsLayer.appendChild(el('line',{x1:panel.x0,y1:panel.y1,x2:panel.x1,y2:panel.y1,stroke:'#111','stroke-width':1.1}));mdsLayer.appendChild(el('line',{x1:panel.x0,y1:panel.y0,x2:panel.x0,y2:panel.y1,stroke:'#111','stroke-width':1.1}));const sx=sc(panel.rows.map(d=>+d.MDS1),panel.x0+42,panel.x1-42),sy=sc(panel.rows.map(d=>+d.MDS2),panel.y1-42,panel.y0+42),pts=panel.rows.map((d,i)=>({d,i,x:sx(+d.MDS1),y:sy(+d.MDS2)}));pts.forEach(p=>drawMarker(mdsLayer,p,p.d.comparison_label===groupSelect.value));const labels=layoutLabels(pts,panel);labels.forEach(a=>{mdsLayer.appendChild(el('line',{x1:a.p.x,y1:a.p.y,x2:a.x,y2:a.y,class:'mdsLeader'}));const lab=drawTextBox(mdsLayer,a.txt.length>24?a.txt.slice(0,22)+'..':a.txt,a.x,a.y+4,'middle',16,groupColor(a.p.d,a.p.i));[lab.t].forEach(node=>{node.style.cursor='pointer';node.addEventListener('click',()=>selectGroup(a.p.d.comparison_label));node.addEventListener('mousemove',ev=>tooltip(mdsTooltip,ev,String(a.p.d.display_label||a.p.d.comparison_label)+'\\ndocs: '+a.p.d.n_docs));node.addEventListener('mouseout',()=>mdsTooltip.style.display='none')})});const xl=el('text',{x:(panel.x0+panel.x1)/2,y:H-22,'text-anchor':'middle','font-size':15,'font-weight':700,fill:'#111'});xl.textContent='MDS1';mdsLayer.appendChild(xl);const yl=el('text',{x:18,y:(panel.y0+panel.y1)/2,'text-anchor':'middle','font-size':15,'font-weight':700,fill:'#111',transform:'rotate(-90 18 '+((panel.y0+panel.y1)/2)+')'});if(!isComparison||panel.name==='Up')yl.textContent='MDS2';mdsLayer.appendChild(yl)});document.getElementById('mdsStats').textContent=GROUP_MDS.length+' groups'}",
    "function drawWaterfall(){waterfallLayer.replaceChildren();const rows=GROUP_TOPIC.filter(d=>d.comparison_label===groupSelect.value).sort((a,b)=>(+b.theta_mean)-(+a.theta_mean)),w=980,h=760,left=170,right=46,top=96,axisY=62,labelY=52,titleY=28,labelSize=svgFontSize(waterfallSvg,17),smallSize=svgFontSize(waterfallSvg,14),maxVal=Math.max(...GROUP_TOPIC.map(d=>Number(d.theta_mean)).filter(Number.isFinite),.001),ticks=niceTicks(maxVal,5),axisMax=Math.max(...ticks,maxVal),plotW=w-left-right,rowH=Math.max(18,Math.min(42,(h-58-top-10)/Math.max(1,rows.length))),x=v=>left+Math.max(0,Number(v||0))/axisMax*plotW,sel=topicSelect.value;ticks.forEach(t=>{const xx=x(t);waterfallLayer.appendChild(el('line',{x1:xx,y1:top-6,x2:xx,y2:h-58,stroke:'#777','stroke-width':.8,opacity:t===0?1:.32}));waterfallLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=Number(t).toFixed(2).replace(/\\.?0+$/,'')});waterfallLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,stroke:'#111','stroke-width':1.8}));waterfallLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent='Document-to-topic probability';rows.forEach((d,i)=>{const y=top+i*rowH,val=Number(d.theta_mean),bw=x(val)-left,isSel=d.topic===sel;waterfallLayer.appendChild(el('text',{x:left-10,y:y+rowH*.62,'font-size':labelSize,'text-anchor':'end','font-weight':700,fill:isSel?'#111':'#333'})).textContent='Topic '+String(d.topic).replace(/^Topic/,'');const r=el('rect',{x:left,y:y+4,width:Math.max(1,bw),height:Math.max(5,rowH-9),fill:color(topicNum(d.topic)),opacity:isSel?.95:.72,stroke:isSel?'#111':'none','stroke-width':isSel?4:0,style:'cursor:pointer'});r.addEventListener('click',()=>{topicSelect.value=d.topic;draw()});r.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,'Topic '+String(d.topic).replace(/^Topic/,'')+'\\ndocument-to-topic probability: '+val.toFixed(4)));r.addEventListener('mouseout',()=>wfTooltip.style.display='none');waterfallLayer.appendChild(r);waterfallLayer.appendChild(el('text',{x:Math.min(w-right,left+bw+6),y:y+rowH*.62,'font-size':smallSize,'font-weight':700,fill:'#334155'})).textContent=val.toFixed(3)});document.getElementById('waterfallStats').textContent=rows.length+' topics'}",
    "function niceTicks(maxVal,n){const out=[];if(!Number.isFinite(maxVal)||maxVal<=0)return[0,1];const raw=maxVal/Math.max(1,n),pow=Math.pow(10,Math.floor(Math.log10(raw))),step=([1,2,5,10].map(v=>v*pow).find(v=>raw<=v)||10*pow);for(let x=0;x<=maxVal+step*.5;x+=step)out.push(x);return out}",
    "function wrapWords(s,maxChars,maxLines){const words=String(s||'').split(/\\s+/),lines=[];let line='';words.forEach(w=>{const next=line?line+' '+w:w;if(next.length>maxChars&&line){lines.push(line);line=w}else line=next});if(line)lines.push(line);if(lines.length>maxLines){lines.length=maxLines;lines[maxLines-1]=lines[maxLines-1].slice(0,Math.max(0,maxChars-2))+'..'}return lines}",
    "function addWrappedLabel(layer,text,x,y,maxChars,maxLines,lineHeight,fontSize,cls){const t=el('text',{x:x,y:y,'text-anchor':'end',class:cls||'pathLabel','font-size':fontSize});wrapWords(text,maxChars,maxLines).forEach((line,i)=>{const sp=el('tspan',{x:x,dy:i===0?0:lineHeight});sp.textContent=line;t.appendChild(sp)});layer.appendChild(t);return t}",
    "function drawPathways(){pathLayer.replaceChildren();const tn=topicNum(topicSelect.value),grp=groupSelect.value,hasGroupPath=PATHWAYS.some(d=>String(d.comparison_label||'')!=='');const allRows=PATHWAYS.filter(d=>(+d.topic_num===tn||+d.topic===tn)&&(!hasGroupPath||String(d.comparison_label||'')===grp));let rows=allRows.slice().sort((a,b)=>(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||''))).slice(0,30);const topicMap=new Map(),groupMap=new Map();PATHWAYS.forEach(d=>{const k=String(d.pathway||''),g=String(d.comparison_label||'');if(!k)return;if(!topicMap.has(k))topicMap.set(k,new Set());topicMap.get(k).add(+d.topic_num||+d.topic);if(g){if(!groupMap.has(k))groupMap.set(k,new Set());groupMap.get(k).add(g)}});document.getElementById('pathStats').textContent=rows.length+' of '+allRows.length+' rows for Topic '+tn+(hasGroupPath?' in selected group':'');const w=1500,h=980,left=770,right=56,top=104,bottom=82,axisY=54,labelY=44,titleY=24,rowH=Math.min(30,(h-top-bottom-42)/Math.max(1,rows.length)),barH=Math.max(10,rowH*.5),totals=rows.map(d=>Math.max(+d.gene_total_universe||0,+d.gene_in||0));if(!allRows.length){const t=el('text',{x:40,y:80,class:'label'});t.textContent=hasGroupPath?'No pathway rows are available for this topic and selected group.':'No pathway rows are available for this topic.';pathLayer.appendChild(t);return}const ticks=niceTicks(Math.max(...totals,1),8),axisMax=Math.max(...ticks,...totals,1),x=v=>left+Math.max(0,+v||0)/axisMax*(w-left-right);pathLayer.appendChild(el('rect',{x:left,y:h-34,width:22,height:12,class:'barIn'}));pathLayer.appendChild(el('text',{x:left+30,y:h-22,class:'pathLegendText'})).textContent='topic genes';pathLayer.appendChild(el('rect',{x:left+190,y:h-34,width:22,height:12,class:'barOut'}));pathLayer.appendChild(el('text',{x:left+220,y:h-22,class:'pathLegendText'})).textContent='universe remainder';ticks.forEach(t=>{const xx=x(t);pathLayer.appendChild(el('line',{x1:xx,y1:top-10,x2:xx,y2:h-bottom-28,class:'pathTick',opacity:t===0?1:.35}));pathLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle',class:'pathTickText'})).textContent=String(Math.round(t))});pathLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,class:'pathAxis'}));pathLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle',class:'pathTickText'})).textContent='Genes in full document gene-universe enrichment';rows.forEach((d,i)=>{const y=top+18+i*rowH,inside=Math.max(0,+d.gene_in||0),total=Math.max(inside,+d.gene_total_universe||0),outside=Math.max(0,total-inside),pkey=String(d.pathway||''),topicSpecific=(topicMap.get(pkey)||new Set()).size===1,groupSpecific=(groupMap.get(pkey)||new Set()).size===1,cls=topicSpecific&&groupSpecific?'pathLabel pathLabelBothSpecific':topicSpecific?'pathLabel pathLabelTopicSpecific':groupSpecific?'pathLabel pathLabelGroupSpecific':'pathLabel',lab=addWrappedLabel(pathLayer,d.pathway,left-12,y+barH*.72,76,rowH<27?1:2,14,16,cls),inW=Math.max(inside>0?3:0,x(inside)-left),outW=Math.max(outside>0?2:0,x(total)-x(inside)),rin=el('rect',{x:left,y:y,width:inW,height:barH,class:'barIn',rx:1}),rout=el('rect',{x:left+inW,y:y,width:outW,height:barH,class:'barOut',rx:1}),msg=String(d.pathway)+'\\nGene in topic: '+inside+'\\nGene out of topic: '+outside+'\\nUniverse hits: '+total+'\\nadj-p: '+(Number.isFinite(Number(d.padj))?Number(d.padj).toExponential(2):'NA')+'\\nGenes: '+String(d.genes||'');[rin,rout,lab].forEach(node=>{node.style.cursor='pointer';node.addEventListener('click',()=>openSubgrn(d));node.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));node.addEventListener('mouseout',()=>pathTooltip.style.display='none')});pathLayer.appendChild(rin);pathLayer.appendChild(rout)})}",
    "function conditionTfRowsForTopic(){const grp=groupSelect.value,topic=topicSelect.value,sorted=TF_TOPIC.filter(d=>String(d.comparison_label||'')===String(grp)&&d.topic===topic).sort((a,b)=>(+b.theta)-(+a.theta)||String(a.tf||'').localeCompare(String(b.tf||''))),top=sorted.slice(0,30);selectedTfUppers.forEach(k=>{if(!top.some(d=>String(d.tf_upper||d.tf||'').toUpperCase()===k)){const sel=sorted.find(d=>String(d.tf_upper||d.tf||'').toUpperCase()===k);if(sel)top.push(sel)}});return top}",
    "function drawTfWaterfall(){tfWaterfallLayer.replaceChildren();const rows=conditionTfRowsForTopic(),selSet=selectedTfSet(),w=760,h=760,left=184,right=42,top=82,labelY=42,titleY=22,labelSize=svgFontSize(tfWaterfallSvg,12),smallSize=svgFontSize(tfWaterfallSvg,12),maxVal=Math.max(...rows.map(d=>+d.theta||0),.001),ticks=niceTicks(maxVal,5),axisMax=Math.max(...ticks,maxVal),plotW=w-left-right,rowH=Math.max(18,Math.min(34,(h-58-top-10)/Math.max(1,rows.length))),x=v=>left+Math.max(0,+v||0)/axisMax*plotW,selectRow=d=>setSelectedTf(d.tf_upper||d.tf||'');ticks.forEach(t=>{const xx=x(t);tfWaterfallLayer.appendChild(el('line',{x1:xx,y1:top-7,x2:xx,y2:h-56,stroke:'#777','stroke-width':.8,opacity:t===0?1:.32}));tfWaterfallLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=Number(t).toFixed(2).replace(/\\.?0+$/,'')});tfWaterfallLayer.appendChild(el('line',{x1:left,y1:56,x2:w-right,y2:56,stroke:'#111','stroke-width':1.8}));tfWaterfallLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent='TF theta for selected topic';rows.forEach((d,i)=>{const y=top+i*rowH,val=+d.theta||0,bw=x(val)-left,key=String(d.tf_upper||d.tf||'').toUpperCase(),isSel=selSet.has(key),label=String(d.tf||'TF'),shortLabel=label.length>18?label.slice(0,16)+'..':label,msg=label+'\\ntheta: '+val.toFixed(4);const lab=el('text',{x:left-10,y:y+rowH*.62,'font-size':labelSize,'text-anchor':'end','font-weight':700,fill:isSel?'#991b1b':'#111',style:'cursor:pointer'});lab.textContent=shortLabel;lab.addEventListener('click',()=>selectRow(d));lab.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));lab.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(lab);const r=el('rect',{x:left,y:y+4,width:Math.max(1,bw),height:Math.max(5,rowH-9),fill:isSel?'#dc2626':groupColor({color:''},i),opacity:isSel?.95:.76,stroke:isSel?'#111':'none','stroke-width':isSel?3:0,style:'cursor:pointer'});r.addEventListener('click',()=>selectRow(d));r.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));r.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(r);const vt=el('text',{x:Math.min(w-right,left+bw+6),y:y+rowH*.62,'font-size':smallSize,'font-weight':700,fill:'#334155',style:'cursor:pointer'});vt.textContent=val.toFixed(3);vt.addEventListener('click',()=>selectRow(d));vt.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,msg));vt.addEventListener('mouseout',()=>wfTooltip.style.display='none');tfWaterfallLayer.appendChild(vt)});updateTfControls()}",
    "function drawWaterfall(){waterfallLayer.replaceChildren();const grp=groupSelect.value;let rows;if(hasSelectedTf()){const tfSet=selectedTfSet(),byTopic=new Map();TF_TOPIC.filter(d=>String(d.comparison_label||'')===String(grp)&&tfSet.has(String(d.tf_upper||d.tf||'').toUpperCase())).forEach(d=>{const k=d.topic,v=byTopic.get(k)||{comparison_label:d.comparison_label,display_label:d.display_label,n_docs:1,topic:d.topic,topic_num:d.topic_num,total:0,n:0,tf:d.tf};v.total+=Number(d.theta)||0;v.n+=1;byTopic.set(k,v)});rows=Array.from(byTopic.values()).map(d=>({comparison_label:d.comparison_label,display_label:d.display_label,n_docs:d.n,topic:d.topic,topic_num:d.topic_num,theta_mean:d.n?d.total/d.n:0,tf:selectedTfLabel()})).sort((a,b)=>(+b.theta_mean)-(+a.theta_mean)||(+a.topic_num)-(+b.topic_num))}else{rows=GROUP_TOPIC.filter(d=>d.comparison_label===grp).sort((a,b)=>(+b.theta_mean)-(+a.theta_mean))}const w=760,h=760,left=135,right=42,top=82,labelY=42,titleY=22,labelSize=svgFontSize(waterfallSvg,15),smallSize=svgFontSize(waterfallSvg,12),maxVal=Math.max(...rows.map(d=>Number(d.theta_mean)).filter(Number.isFinite),.001),ticks=niceTicks(maxVal,5),axisMax=Math.max(...ticks,maxVal),plotW=w-left-right,rowH=Math.max(18,Math.min(38,(h-58-top-10)/Math.max(1,rows.length))),x=v=>left+Math.max(0,Number(v||0))/axisMax*plotW,sel=topicSelect.value;ticks.forEach(t=>{const xx=x(t);waterfallLayer.appendChild(el('line',{x1:xx,y1:top-6,x2:xx,y2:h-58,stroke:'#777','stroke-width':.8,opacity:t===0?1:.32}));waterfallLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=Number(t).toFixed(2).replace(/\\.?0+$/,'')});waterfallLayer.appendChild(el('line',{x1:left,y1:56,x2:w-right,y2:56,stroke:'#111','stroke-width':1.8}));waterfallLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle','font-size':smallSize,'font-weight':700,fill:'#111'})).textContent=hasSelectedTf()?'Selected TF mean theta by topic':'Document-to-topic probability';rows.forEach((d,i)=>{const y=top+i*rowH,val=Number(d.theta_mean),bw=x(val)-left,isSel=d.topic===sel;waterfallLayer.appendChild(el('text',{x:left-10,y:y+rowH*.62,'font-size':labelSize,'text-anchor':'end','font-weight':700,fill:isSel?'#111':'#333'})).textContent='Topic '+String(d.topic).replace(/^Topic/,'');const r=el('rect',{x:left,y:y+4,width:Math.max(1,bw),height:Math.max(5,rowH-9),fill:color(topicNum(d.topic)),opacity:isSel?.95:.72,stroke:isSel?'#111':'none','stroke-width':isSel?4:0,style:'cursor:pointer'});r.addEventListener('click',()=>{topicSelect.value=d.topic;draw()});r.addEventListener('mousemove',ev=>tooltip(wfTooltip,ev,'Topic '+String(d.topic).replace(/^Topic/,'')+'\\nprobability: '+val.toFixed(4)+(hasSelectedTf()?'\\nTFs: '+selectedTfLabel():'')));r.addEventListener('mouseout',()=>wfTooltip.style.display='none');waterfallLayer.appendChild(r);waterfallLayer.appendChild(el('text',{x:Math.min(w-right,left+bw+6),y:y+rowH*.62,'font-size':smallSize,'font-weight':700,fill:'#334155'})).textContent=val.toFixed(3)});document.getElementById('waterfallStats').textContent=rows.length+' topics'+(hasSelectedTf()?' for '+selectedTfLabel():'')}",
    "function tfCountKey(d){return selectedTfKey()+'|'+String(d.comparison_label||d.comparison_id||'')+'|'+String(d.topic_num||d.topic||'')+'|'+String(d.pathway_key||d.pathway_norm_key||d.pathway||'')}function splitGenes(s){return String(s||'').split(/[;,]/).map(x=>x.trim()).filter(Boolean)}function rowGeneSet(d){return new Set(splitGenes(d.genes||d.overlap_genes||''))}function contextMatchEdge(e,ctx,genes){return e.subgrn_context_id?e.subgrn_context_id===ctx.subgrn_context_id:(String(e.comparison_id||'')===String(ctx.comparison_id||'')&&String(e.direction||'')===String(ctx.direction||'')&&genes.has(String(e.gene_key||'')))}function countTfGenesInPayload(payload,ctx,row){const genes=rowGeneSet(row),allowed=subgrnGeneSet(ctx),tfSet=selectedTfSet();allowed.forEach(g=>genes.add(g));const out=new Set(),byTf=new Map(),edges=new Set();((payload&&payload.tf_gene_edges)||[]).forEach(e=>{const tf=String(e.tf_upper||e.tf||'').toUpperCase(),gene=String(e.gene_key||'');if(tfSet.has(tf)&&genes.has(gene)&&contextMatchEdge(e,ctx,genes)){out.add(gene);const rec=byTf.get(tf)||{tf:tf,label:tfLabelMap.get(tf)||String(e.tf||tf),genes:new Set()};rec.genes.add(gene);byTf.set(tf,rec);edges.add(tf+'\\t'+gene)}});const rows=Array.from(byTf.values()).map(r=>({tf:r.tf,label:r.label,genes:Array.from(r.genes).sort(),n:r.genes.size})).sort((a,b)=>b.n-a.n||a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));return{genes:Array.from(out).sort(),byTf:rows,edgeN:edges.size}}function mergeTfSupport(parts){const genes=new Set(),byTf=new Map();let edgeN=0;(parts||[]).forEach(p=>{if(Array.isArray(p))p={genes:p,byTf:[],edgeN:0};(p.genes||[]).forEach(g=>genes.add(g));(p.byTf||[]).forEach(r=>{const key=String(r.tf||'').toUpperCase();if(!key)return;const rec=byTf.get(key)||{tf:key,label:r.label||tfLabelMap.get(key)||key,genes:new Set()};(r.genes||[]).forEach(g=>rec.genes.add(g));byTf.set(key,rec)});edgeN+=Number(p.edgeN)||0});const rows=Array.from(byTf.values()).map(r=>({tf:r.tf,label:r.label,genes:Array.from(r.genes).sort(),n:r.genes.size})).sort((a,b)=>b.n-a.n||a.label.localeCompare(b.label,undefined,{sensitivity:'base'}));const detail=[];if(rows.length)detail.push('Summary: '+genes.size+' unique target genes; '+rows.length+' selected TFs with hits; '+edgeN+' TF-target pairs');rows.forEach(r=>detail.push(r.label+' ('+r.n+'): '+r.genes.join(',')));return{n:genes.size,genes:detail.length?detail:Array.from(genes).sort(),target_genes:Array.from(genes).sort(),byTf:rows,tfN:rows.length,edgeN:edgeN}}function ensureTfPathwayCounts(rows){if(!hasSelectedTf())return;rows.forEach(d=>{const key=tfCountKey(d);if(tfPathwayCounts.has(key)||tfPathwayPending.has(key))return;const ctxs=subgrnContextRows(d);if(!ctxs.length){tfPathwayCounts.set(key,mergeTfSupport([]));return}tfPathwayPending.add(key);Promise.all(ctxs.map(ctx=>loadSubgrnPayload(ctx.payload_file).then(payload=>countTfGenesInPayload(payload,ctx,d)).catch(()=>({genes:[],byTf:[],edgeN:0})))).then(parts=>{tfPathwayCounts.set(key,mergeTfSupport(parts));tfPathwayPending.delete(key);drawPathways()})})}",
    "function drawPathways(){pathLayer.replaceChildren();const tn=topicNum(topicSelect.value),grp=groupSelect.value,hasGroupPath=PATHWAYS.some(d=>String(d.comparison_label||'')!=='');const allRows=PATHWAYS.filter(d=>(+d.topic_num===tn||+d.topic===tn)&&(!hasGroupPath||String(d.comparison_label||'')===grp));const baseRows=allRows.slice().sort((a,b)=>(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||''))).slice(0,30);ensureTfPathwayCounts(baseRows);let rows=baseRows;if(hasSelectedTf()){rows=baseRows.slice().sort((a,b)=>((tfPathwayCounts.get(tfCountKey(b))||{n:0}).n-(tfPathwayCounts.get(tfCountKey(a))||{n:0}).n)||((tfPathwayCounts.get(tfCountKey(b))||{tfN:0}).tfN-(tfPathwayCounts.get(tfCountKey(a))||{tfN:0}).tfN)||(+a.padj)-(+b.padj)||(+b.gene_in)-(+a.gene_in)||String(a.pathway||'').localeCompare(String(b.pathway||'')))}const topicMap=new Map(),groupMap=new Map();PATHWAYS.forEach(d=>{const k=String(d.pathway||''),g=String(d.comparison_label||'');if(!k)return;if(!topicMap.has(k))topicMap.set(k,new Set());topicMap.get(k).add(+d.topic_num||+d.topic);if(g){if(!groupMap.has(k))groupMap.set(k,new Set());groupMap.get(k).add(g)}});document.getElementById('pathStats').textContent=rows.length+' of '+allRows.length+' rows for Topic '+tn+(hasGroupPath?' in selected group':'')+(hasSelectedTf()?' | selected TF target genes and per-TF counts shown':'');const w=1500,h=980,left=770,right=110,top=104,bottom=82,axisY=54,labelY=44,titleY=24,rowH=Math.min(30,(h-top-bottom-42)/Math.max(1,rows.length)),barH=Math.max(10,rowH*.5),totals=rows.map(d=>Math.max(+d.gene_total_universe||0,+d.gene_in||0));if(!allRows.length){const t=el('text',{x:40,y:80,class:'label'});t.textContent=hasGroupPath?'No pathway rows are available for this topic and selected group.':'No pathway rows are available for this topic.';pathLayer.appendChild(t);return}const ticks=niceTicks(Math.max(...totals,1),8),axisMax=Math.max(...ticks,...totals,1),x=v=>left+Math.max(0,+v||0)/axisMax*(w-left-right);pathLayer.appendChild(el('rect',{x:left,y:h-34,width:22,height:12,class:'barIn'}));pathLayer.appendChild(el('text',{x:left+30,y:h-22,class:'pathLegendText'})).textContent='topic genes';pathLayer.appendChild(el('rect',{x:left+190,y:h-34,width:22,height:12,class:'barOut'}));pathLayer.appendChild(el('text',{x:left+220,y:h-22,class:'pathLegendText'})).textContent='universe remainder';if(hasSelectedTf())pathLayer.appendChild(el('text',{x:w-8,y:h-22,'text-anchor':'end',class:'pathLegendText'})).textContent=selectedTfUppers.length>1?'TF genes/TFs':'TF targets';ticks.forEach(t=>{const xx=x(t);pathLayer.appendChild(el('line',{x1:xx,y1:top-10,x2:xx,y2:h-bottom-28,class:'pathTick',opacity:t===0?1:.35}));pathLayer.appendChild(el('text',{x:xx,y:labelY,'text-anchor':'middle',class:'pathTickText'})).textContent=String(Math.round(t))});pathLayer.appendChild(el('line',{x1:left,y1:axisY,x2:w-right,y2:axisY,class:'pathAxis'}));pathLayer.appendChild(el('text',{x:(left+w-right)/2,y:titleY,'text-anchor':'middle',class:'pathTickText'})).textContent='Genes in full document gene-universe enrichment';rows.forEach((d,i)=>{const y=top+18+i*rowH,inside=Math.max(0,+d.gene_in||0),total=Math.max(inside,+d.gene_total_universe||0),outside=Math.max(0,total-inside),pkey=String(d.pathway||''),topicSpecific=(topicMap.get(pkey)||new Set()).size===1,groupSpecific=(groupMap.get(pkey)||new Set()).size===1,cls=topicSpecific&&groupSpecific?'pathLabel pathLabelBothSpecific':topicSpecific?'pathLabel pathLabelTopicSpecific':groupSpecific?'pathLabel pathLabelGroupSpecific':'pathLabel',support=hasSelectedTf()?(tfPathwayCounts.get(tfCountKey(d))||{n:0,genes:[],tfN:0,edgeN:0}):null,lab=addWrappedLabel(pathLayer,d.pathway,left-12,y+barH*.72,76,rowH<27?1:2,14,16,cls),inW=Math.max(inside>0?3:0,x(inside)-left),outW=Math.max(outside>0?2:0,x(total)-x(inside)),rin=el('rect',{x:left,y:y,width:inW,height:barH,class:'barIn',rx:1}),rout=el('rect',{x:left+inW,y:y,width:outW,height:barH,class:'barOut',rx:1}),extra=support?'\\nSelected TFs: '+selectedTfLabel()+'\\nUnique selected-TF target genes in pathway: '+support.n+'\\nSelected TFs with hits: '+(Number(support.tfN)||0)+' / '+selectedTfUppers.length+'\\nTF-target pairs represented: '+(Number(support.edgeN)||0)+(support.genes&&support.genes.length?'\\nPer-TF target genes: '+support.genes.join('; '):''):'',msg=String(d.pathway)+'\\nGene in topic: '+inside+'\\nGene out of topic: '+outside+'\\nUniverse hits: '+total+'\\nadj-p: '+(Number.isFinite(Number(d.padj))?Number(d.padj).toExponential(2):'NA')+extra+'\\nGenes: '+String(d.genes||'');[rin,rout,lab].forEach(node=>{node.style.cursor='pointer';node.addEventListener('click',()=>openSubgrn(d));node.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));node.addEventListener('mouseout',()=>pathTooltip.style.display='none')});pathLayer.appendChild(rin);pathLayer.appendChild(rout);if(support){const n=Number(support.n)||0,tfN=Number(support.tfN)||0,badgeText=selectedTfUppers.length>1?String(n)+'/'+String(tfN):String(n);if(n>0){const badge=el('rect',{x:w-72,y:y-2,width:64,height:barH+4,rx:8,fill:'#dc2626',stroke:'#7f1d1d','stroke-width':1.4,opacity:.95,style:'cursor:pointer'});badge.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));badge.addEventListener('mouseout',()=>pathTooltip.style.display='none');badge.addEventListener('click',()=>openSubgrn(d));pathLayer.appendChild(badge);rin.setAttribute('stroke','#dc2626');rin.setAttribute('stroke-width',2.4);rout.setAttribute('stroke','#dc2626');rout.setAttribute('stroke-width',1.6)}const ct=el('text',{x:w-40,y:y+barH*.78,'text-anchor':'middle','font-size':15,'font-weight':700,fill:n>0?'#fff':'#94a3b8',style:'cursor:pointer'});ct.textContent=badgeText;ct.addEventListener('mousemove',ev=>tooltip(pathTooltip,ev,msg));ct.addEventListener('mouseout',()=>pathTooltip.style.display='none');ct.addEventListener('click',()=>openSubgrn(d));pathLayer.appendChild(ct)}})}",
    .m3tb_subgrn_js(),
    "function draw(){updateSelectedContext();drawMds();drawTfWaterfall();drawWaterfall();drawPathways()}init();",
    "</script></body></html>"
  )
  html <- sub(
    "</style>",
    paste0(
      ".paneMeta{display:flex;flex-direction:column;gap:3px;align-items:flex-end;min-width:0}",
      ".contextBadge{max-width:560px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;",
      "background:#111;color:#fff;border-radius:3px;padding:4px 7px;font-size:12px;line-height:1.1}",
      "</style>"
    ),
    html,
    fixed = TRUE
  )
  html <- sub(
    "<section class=\"pane\"><div class=\"paneHead\"><h2>Pathways</h2><span id=\"pathStats\" class=\"meta\"></span></div><div class=\"body\">",
    "<section class=\"pane\"><div class=\"paneHead\"><h2>Pathways</h2><div class=\"paneMeta\"><div id=\"selectedContext\" class=\"contextBadge\"></div><span id=\"pathStats\" class=\"meta\"></span></div></div><div class=\"body\">",
    html,
    fixed = TRUE
  )
  html <- .m3tb_polish_topic_report_html(html)
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, out_html, useBytes = TRUE)
  out_html
}

.m3tb_open_png <- function(filename, width, height, res = 140L) {
  opened <- tryCatch({
    grDevices::png(filename, width = width, height = height, res = res, bg = "white", type = "cairo")
    TRUE
  }, error = function(e) FALSE)
  if (!opened) {
    opened <- tryCatch({
      grDevices::png(filename, width = width, height = height, res = res, bg = "white", type = "cairo-png")
      TRUE
    }, error = function(e) FALSE)
  }
  if (!opened) {
    grDevices::png(filename, width = width, height = height, res = res, bg = "white")
  }
  invisible(TRUE)
}

.m3tb_write_review_pngs <- function(score_result, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  mds <- data.table::as.data.table(score_result$mds_points)
  scores <- data.table::as.data.table(score_result$scores)
  if (!nrow(mds)) return(data.table::data.table())
  k_values <- sort(unique(as.integer(mds$k)))
  rows <- lapply(k_values, function(k_value) {
    left <- file.path(out_dir, sprintf("theta_phi_topic_distance_correlation_k%d.png", k_value))
    right <- file.path(out_dir, sprintf("theta_group_mds_k%d.png", k_value))
    .m3tb_open_png(left, width = 1000, height = 1000, res = 140)
    graphics::plot.new()
    graphics::title(sprintf("K%d theta score summary", k_value))
    sub <- scores[as.integer(k) == as.integer(k_value)]
    vals <- if (nrow(sub)) as.numeric(sub$theta_condition_separation_score) else numeric()
    if (length(vals) && any(is.finite(vals))) {
      plot_vals <- vals
      plot_vals[!is.finite(plot_vals)] <- 0
      graphics::barplot(plot_vals, names.arg = sub$method_setup, las = 2, col = "#4E79A7", ylab = "Separation score")
      if (any(!is.finite(vals))) {
        graphics::mtext("Unscored singleton groups are plotted at 0.", side = 3, line = -1, cex = 0.75)
      }
    } else {
      graphics::text(0.5, 0.5, "No scored groups")
    }
    grDevices::dev.off()
    sub_mds <- mds[as.integer(k) == as.integer(k_value)]
    if (nrow(sub_mds)) {
      p_mds <- .m3tb_condition_mds_plot(
        sub_mds,
        title = sprintf("K%d condition/comparison MDS from theta", k_value)
      )
      n_panels <- if ("panel_label" %in% names(sub_mds)) {
        n_dirs <- if (any(sub_mds$doc_design == "comparison", na.rm = TRUE)) 2L else 1L
        data.table::uniqueN(sub_mds$panel_label) * n_dirs
      } else {
        if (any(sub_mds$doc_design == "comparison", na.rm = TRUE)) 2L else 1L
      }
      if (is.finite(n_panels) && n_panels > 4L) {
        .m3tb_open_png(right, width = 5600, height = 3600, res = 300)
      } else {
        .m3tb_open_png(right, width = 2600, height = 1500, res = 180)
      }
      print(p_mds)
    } else {
      .m3tb_open_png(right, width = 1800, height = 1100, res = 150)
      graphics::plot.new()
      graphics::text(0.5, 0.5, "No MDS rows")
    }
    grDevices::dev.off()
    data.table::data.table(
      k = k_value,
      left_src = file.path("assets", basename(left)),
      right_src = file.path("assets", basename(right)),
      has_score = length(vals) > 0L && any(is.finite(vals))
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.m3tb_write_theta_review_html <- function(score_result, out_file, image_manifest) {
  manifest_json <- .m3tb_json_for_html(image_manifest)
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/><title>Combined theta review</title>",
    "<style>body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif;overflow:hidden;font-weight:700}.top{height:48px;box-sizing:border-box;background:#fff;border-bottom:1px solid #d6d6d0;padding:5px 10px;display:flex;justify-content:space-between;gap:12px;align-items:center}h1{font-size:15px;line-height:1.1;margin:0}.controls{display:flex;gap:8px;align-items:center;font-size:13px;min-width:min(52vw,760px)}input[type=range]{width:100%;accent-color:#1f78b4}.dashboard{height:calc(100vh - 48px);display:grid;grid-template-columns:minmax(0,1.2fr) minmax(0,3.4fr);gap:6px;padding:6px;box-sizing:border-box}.dashboard.noScore{grid-template-columns:1fr}.dashboard.noScore .scorePane{display:none}.pane{min-width:0;min-height:0;background:#fff;border:1px solid #d6d6d0;overflow:hidden;display:flex;align-items:center;justify-content:center}.plotImg{width:100%;height:100%;display:block;background:#fff;object-fit:contain}.rightPlot{object-fit:contain}@media(max-width:1200px){body{overflow:auto}.dashboard{height:auto;grid-template-columns:1fr}.pane{height:760px}}</style></head><body>",
    "<div class=\"top\"><div><h1>Combined theta review</h1><div id=\"plotMeta\" style=\"font-size:10px;color:#555\"></div></div><label class=\"controls\">K <input id=\"kSlider\" type=\"range\" min=\"0\" max=\"0\" step=\"1\" value=\"0\"/></label></div>",
    "<main id=\"dashboard\" class=\"dashboard\"><section class=\"pane scorePane\"><img id=\"correlationImg\" class=\"plotImg\" alt=\"Theta score summary\"/></section><section class=\"pane\"><img id=\"groupMdsImg\" class=\"plotImg rightPlot\" alt=\"Condition and comparison MDS from theta\"/></section></main>",
    "<script>",
    paste0("const IMAGE_MANIFEST=", manifest_json, ";"),
    "const slider=document.getElementById('kSlider'),meta=document.getElementById('plotMeta'),left=document.getElementById('correlationImg'),right=document.getElementById('groupMdsImg'),dashboard=document.getElementById('dashboard');slider.max=Math.max(0,IMAGE_MANIFEST.length-1);function render(){const i=Math.max(0,Math.min(IMAGE_MANIFEST.length-1,Number(slider.value)||0));const item=IMAGE_MANIFEST[i]||{};left.src=item.left_src||'';right.src=item.right_src||'';dashboard.classList.toggle('noScore',!item.has_score);meta.textContent=item.k?(item.has_score?'Showing K'+item.k+'; score summary and MDS use the same source tables.':'Showing K'+item.k+'; no replicate-resolved separation score is available, so the MDS is shown full width.'):'No review images available.'}slider.addEventListener('input',render);render();",
    "</script></body></html>"
  )
  writeLines(html, out_file, useBytes = TRUE)
  out_file
}

.m3tb_review_row_value <- function(row, name, default = NA_character_) {
  if (name %in% names(row)) as.character(row[[name]][[1L]]) else default
}

.m3tb_review_report_slug <- function(row, k) {
  run_id <- .m3tb_review_row_value(row, "run_id", "")
  run_short <- if (grepl("^run_[0-9]+", run_id)) {
    sub("^(run_[0-9]+).*$", "\\1", run_id)
  } else {
    run_id
  }
  model <- tolower(.m3tb_review_row_value(
    row,
    "model_label",
    .m3tb_review_row_value(row, "method", "model")
  ))
  topic_space <- .m3tb_topic_space_value(row, "topic_space", "raw")
  topic_count <- suppressWarnings(as.integer(.m3tb_topic_space_value(
    row,
    "topic_count",
    k
  )))
  space_label <- if (identical(topic_space, "combined")) {
    paste0("combined_K", topic_count)
  } else {
    "raw"
  }
  parts <- c(
    run_short,
    model,
    paste0("K", as.integer(k)),
    space_label
  )
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (!length(parts)) parts <- c(.m3tb_review_row_value(row, "method_setup", "method"), paste0("K", as.integer(k)))
  .safe_filename(paste(parts, collapse = "_"))
}

.m3tb_relative_html_path <- function(path, base_dir) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  base_dir <- normalizePath(base_dir, winslash = "/", mustWork = FALSE)
  rel <- sub(paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", base_dir), "/?"), "", path)
  rel
}

.m3tb_write_k_report_index <- function(out_file, reports, title) {
  reports <- data.table::as.data.table(reports)
  if (!nrow(reports)) return(NA_character_)
  reports[, src := vapply(path, .m3tb_relative_html_path, character(1L), base_dir = dirname(out_file))]
  if (!"report_label" %in% names(reports)) reports[, report_label := paste0("K", as.integer(k))]
  if (!"trained_k" %in% names(reports)) reports[, trained_k := as.integer(k)]
  if (!"topic_space" %in% names(reports)) reports[, topic_space := "raw"]
  reports[, space_order := data.table::fifelse(topic_space == "raw", 1L, 2L)]
  data.table::setorder(reports, trained_k, space_order)
  json <- .m3tb_json_for_html(reports[, .(
    label, k, k_label = report_label, src
  )])
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    "<style>html,body{width:100%;height:100%;overflow:hidden}body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif}.top{height:46px;display:flex;align-items:center;justify-content:space-between;gap:12px;padding:6px 12px;border-bottom:1px solid #d6d6d0;background:#fff;box-sizing:border-box}h1{font-size:18px;line-height:1;margin:0;font-weight:700}.controls{display:flex;gap:8px;align-items:center;font-weight:700;white-space:nowrap}select{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #aaa;background:#fff;color:#111;border-radius:3px;padding:4px 7px}iframe{display:block;width:100vw;height:calc(100vh - 46px);border:0;background:#fff}</style></head><body>",
    "<div class=\"top\">",
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    "<div class=\"controls\"><label>K <select id=\"kSelect\"></select></label></div></div><iframe id=\"frame\"></iframe>",
    "<script>",
    paste0("const REPORTS=", json, ";"),
    "const kSelect=document.getElementById('kSelect'),frame=document.getElementById('frame');function embedSrc(src){return src+(String(src).indexOf('?')>=0?'&':'?')+'embed=1'}REPORTS.forEach((r,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=r.k_label;kSelect.appendChild(o)});function load(){const hit=REPORTS[Number(kSelect.value)||0]||REPORTS[0];if(hit){frame.removeAttribute('src'+'doc');frame.src=embedSrc(hit.src);document.title=hit.label}}kSelect.addEventListener('change',load);load();",
    "</script></body></html>"
  )
  writeLines(html, out_file, useBytes = TRUE)
  out_file
}

.m3tb_write_combined_report_index <- function(out_file, reports, title) {
  reports <- data.table::as.data.table(reports)
  if (!nrow(reports)) return(NA_character_)
  reports[, src := vapply(path, .m3tb_relative_html_path, character(1L), base_dir = dirname(out_file))]
  if (!"method_setup" %in% names(reports)) reports[, method_setup := label]
  if (!"run_id" %in% names(reports)) reports[, run_id := method_setup]
  has_k <- "k" %in% names(reports) && any(is.finite(suppressWarnings(as.numeric(reports$k))))
  if (!has_k) reports[, k := NA_integer_]
  json <- .m3tb_json_for_html(reports[, .(label, k, method_setup, run_id, src)])
  k_control <- if (has_k) "<label>K <select id=\"kSelect\"></select></label>" else ""
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    "<style>html,body{width:100%;height:100%;overflow:hidden}body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif}.top{height:46px;display:flex;align-items:center;justify-content:space-between;gap:12px;padding:6px 12px;border-bottom:1px solid #d6d6d0;background:#fff;box-sizing:border-box}h1{font-size:18px;line-height:1;margin:0;font-weight:700}.controls{display:flex;gap:8px;align-items:center;font-weight:700;white-space:nowrap}select{font:700 12px Arial,Helvetica,sans-serif;border:1px solid #aaa;background:#fff;color:#111;border-radius:3px;padding:4px 7px}iframe{display:block;width:100vw;height:calc(100vh - 46px);border:0;background:#fff}</style></head><body>",
    "<div class=\"top\">",
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    paste0("<div class=\"controls\"><label>Method <select id=\"methodSelect\"></select></label>", k_control, "<label>Report <select id=\"reportSelect\"></select></label></div></div><iframe id=\"frame\"></iframe>"),
    "<script>",
    paste0("const REPORTS=", json, ";"),
    "const methodSelect=document.getElementById('methodSelect'),kSelect=document.getElementById('kSelect'),reportSelect=document.getElementById('reportSelect'),frame=document.getElementById('frame');function embedSrc(src){return src+(String(src).indexOf('?')>=0?'&':'?')+'embed=1'}function uniq(x){return Array.from(new Set(x))}function fillSelect(sel,vals,current){if(!sel)return;sel.replaceChildren();vals.forEach(v=>{const o=document.createElement('option');o.value=String(v);o.textContent=String(v);sel.appendChild(o)});if(current&&vals.map(String).includes(String(current)))sel.value=String(current)}function filtered(){return REPORTS.filter(r=>(!methodSelect.value||r.method_setup===methodSelect.value)&&(!kSelect||!kSelect.value||String(r.k)===kSelect.value))}function refresh(changed){const curM=methodSelect.value,curK=kSelect?kSelect.value:'',curR=reportSelect.value;let rows=REPORTS;if(changed!=='method')fillSelect(methodSelect,uniq(REPORTS.map(r=>r.method_setup)),curM);rows=REPORTS.filter(r=>!methodSelect.value||r.method_setup===methodSelect.value);if(kSelect&&changed!=='k')fillSelect(kSelect,uniq(rows.map(r=>String(r.k))).sort((a,b)=>Number(a)-Number(b)),curK);rows=filtered();fillSelect(reportSelect,rows.map((r,i)=>String(i)),curR);Array.from(reportSelect.options).forEach((o,i)=>{const r=rows[i];if(r)o.textContent=r.label});load()}function load(){const rows=filtered();const hit=rows[Number(reportSelect.value)||0]||rows[0]||REPORTS[0];if(hit){methodSelect.value=hit.method_setup;if(kSelect)kSelect.value=String(hit.k);frame.removeAttribute('src'+'doc');frame.src=embedSrc(hit.src);document.title=hit.label}}methodSelect.addEventListener('change',()=>refresh('method'));if(kSelect)kSelect.addEventListener('change',()=>refresh('k'));reportSelect.addEventListener('change',load);refresh();",
    "</script></body></html>"
  )
  writeLines(html, out_file, useBytes = TRUE)
  out_file
}

.m3tb_review_worker_count <- function(n_tasks,
                                      available_bytes = NULL,
                                      detected_cores = NULL) {
  n_tasks <- suppressWarnings(as.integer(n_tasks)[[1L]])
  if (!is.finite(n_tasks) || n_tasks < 1L) return(1L)
  opt <- getOption("craftgrn.review.workers", NULL)
  env <- Sys.getenv("CRAFTGRN_REVIEW_WORKERS", unset = "")
  if (!nzchar(env)) env <- Sys.getenv("CRAFTGRN_REVIEW_THREADS", unset = "")
  raw <- if (!is.null(opt)) opt else env
  explicit_workers <- !is.null(opt) || nzchar(env)
  if (!isTRUE(explicit_workers) &&
      identical(Sys.getenv("TESTTHAT", unset = ""), "true") &&
      !isTRUE(getOption("craftgrn.review.auto_workers_in_tests", FALSE))) {
    return(1L)
  }
  requested <- suppressWarnings(as.integer(raw)[[1L]])
  max_workers <- suppressWarnings(as.integer(
    getOption("craftgrn.review.max_workers", 8L)
  )[[1L]])
  if (!is.finite(max_workers) || max_workers < 1L) max_workers <- 8L
  if (!is.finite(requested) || requested < 1L) requested <- max_workers
  if (is.null(detected_cores)) {
    detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  }
  detected_cores <- suppressWarnings(as.integer(detected_cores)[[1L]])
  if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
  if (is.null(available_bytes)) available_bytes <- .available_memory_bytes()
  available_bytes <- suppressWarnings(as.numeric(available_bytes)[[1L]])
  if (!is.finite(available_bytes) || available_bytes <= 0) return(1L)
  reserve_gb <- suppressWarnings(as.numeric(
    getOption("craftgrn.review.memory_reserve_gb", 32)
  )[[1L]])
  per_worker_gb <- suppressWarnings(as.numeric(
    getOption("craftgrn.review.memory_per_worker_gb", 6)
  )[[1L]])
  if (!is.finite(reserve_gb) || reserve_gb < 0) reserve_gb <- 32
  if (!is.finite(per_worker_gb) || per_worker_gb <= 0) per_worker_gb <- 6
  reserve_bytes <- max(reserve_gb * 1024^3, 0.25 * available_bytes)
  memory_workers <- floor(max(0, available_bytes - reserve_bytes) / (per_worker_gb * 1024^3))
  max(
    1L,
    min(
      n_tasks,
      requested,
      max_workers,
      detected_cores,
      as.integer(memory_workers)
    )
  )
}

.m3tb_review_lapply <- function(x, fun, workers) {
  workers <- max(1L, suppressWarnings(as.integer(workers)[[1L]]))
  if (workers <= 1L) {
    return(lapply(x, fun))
  }
  parent_lib_paths <- .libPaths()
  old_r_libs <- Sys.getenv("R_LIBS", unset = NA_character_)
  restore_r_libs <- function() {
    if (is.na(old_r_libs)) {
      Sys.unsetenv("R_LIBS")
    } else {
      Sys.setenv(R_LIBS = old_r_libs)
    }
  }
  on.exit(restore_r_libs(), add = TRUE)
  Sys.setenv(R_LIBS = paste(parent_lib_paths, collapse = .Platform$path.sep))
  cl <- parallel::makeCluster(workers, type = "PSOCK")
  restore_r_libs()
  cluster_running <- TRUE
  on.exit({
    if (isTRUE(cluster_running)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
  worker_bootstrap <- function(lib_paths) {
    .libPaths(unique(c(lib_paths, .libPaths())))
    library(craftgrn)
    data.table::setDTthreads(1L)
    options(craftgrn.topic_review.fast_summary = TRUE)
    NULL
  }
  environment(worker_bootstrap) <- baseenv()
  parallel::clusterCall(cl, worker_bootstrap, parent_lib_paths)
  result <- tryCatch(
    .m3tb_review_par_lapply(cl, x, fun),
    error = identity
  )
  if (!inherits(result, "error")) return(result)

  try(parallel::stopCluster(cl), silent = TRUE)
  cluster_running <- FALSE
  .log_warn(
    paste0(
      "Parallel Module 3 review workers failed; retrying ",
      length(x), " task(s) sequentially. Cause: ", conditionMessage(result)
    )
  )
  lapply(x, fun)
}

.m3tb_review_par_lapply <- function(cl, x, fun) {
  parallel::parLapplyLB(cl, x, fun)
}

.m3tb_review_row_error <- function(i, err) {
  structure(
    list(index = i, message = conditionMessage(err)),
    class = "m3tb_review_row_error"
  )
}

.m3tb_filter_pathways_to_subgrn_manifest <- function(pathways, subgrn_manifest, per_group = FALSE) {
  pathways <- data.table::as.data.table(pathways)
  subgrn_manifest <- data.table::as.data.table(subgrn_manifest)
  if (!nrow(pathways) || !nrow(subgrn_manifest)) return(pathways[0])
  if (!"pathway_key" %in% names(pathways) && "pathway_norm_key" %in% names(pathways)) {
    data.table::setnames(pathways, "pathway_norm_key", "pathway_key")
  }
  if (!"topic_num" %in% names(pathways) && "topic" %in% names(pathways)) {
    pathways[, topic_num := suppressWarnings(as.integer(topic))]
  }
  pathways[, `:=`(
    topic_num = suppressWarnings(as.integer(topic_num)),
    pathway_key = as.character(pathway_key)
  )]
  subgrn_manifest[, `:=`(
    topic_num = suppressWarnings(as.integer(topic_num)),
    pathway_key = as.character(pathway_key)
  )]
  if (isTRUE(per_group)) {
    if (!"comparison_label" %in% names(pathways)) return(pathways[0])
    keys <- unique(subgrn_manifest[, .(
      comparison_label = as.character(comparison_label),
      topic_num,
      pathway_key
    )])
    pathways[, comparison_label := as.character(comparison_label)]
    return(merge(pathways, keys, by = c("comparison_label", "topic_num", "pathway_key"), all = FALSE, sort = FALSE))
  }
  keys <- unique(subgrn_manifest[, .(topic_num, pathway_key)])
  merge(pathways, keys, by = c("topic_num", "pathway_key"), all = FALSE, sort = FALSE)
}

.m3tb_build_review_pages_for_row <- function(i,
                                             model_rows,
                                             score_result,
                                             topic_page_dir,
                                             condition_page_dir,
                                             subgrn_payload_dir,
                                             condition_pair_payload_dir) {
  old_dt_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
  data.table::setDTthreads(1L)
  report_state <- attr(score_result, "report_state") %||% list()

  row <- model_rows[i]
  k <- as.integer(row$selected_k[[1L]])
  theta <- tryCatch(
    .m3tb_topic_space_theta(row),
    error = function(e) NULL
  )
  if (is.null(theta)) {
    return(list(
      topic = data.table::data.table(),
      condition = data.table::data.table(),
      subgrn_manifest = data.table::data.table(),
      condition_payload = data.table::data.table()
    ))
  }

  extraction_dir <- .m3tb_find_extraction_subdir(row)
  topic_space <- .m3tb_topic_space_value(row, "topic_space", "combined")
  slug <- .m3tb_review_report_slug(row, k)
  run_id <- .m3tb_review_row_value(row, "run_id", row$method_setup[[1L]])

  topic_pathways <- .m3tb_read_pathway_tables(
    extraction_dir,
    per_group = FALSE,
    model_dir = row$model_dir[[1L]],
    compute_universe = TRUE,
    pathway_file = .m3tb_topic_space_file(
      extraction_dir,
      topic_space,
      "overall_pathways"
    )
  )
  condition_pathways <- .m3tb_read_condition_pathway_tables(
    extraction_dir,
    model_dir = row$model_dir[[1L]],
    compute_universe = TRUE,
    pathway_file = .m3tb_topic_space_file(
      extraction_dir,
      topic_space,
      "condition_pathways"
    )
  )
  max_context_raw <- Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_MAX_CONTEXTS", unset = "")
  payload_disabled <- nzchar(max_context_raw) &&
    is.finite(suppressWarnings(as.integer(max_context_raw[[1L]]))) &&
    suppressWarnings(as.integer(max_context_raw[[1L]])) <= 0L
  subgrn_manifest <- if (isTRUE(payload_disabled)) {
    data.table::data.table()
  } else {
    .m3tb_write_pathway_subgrn_payload(
      extraction_dir = extraction_dir,
      model_dir = row$model_dir[[1L]],
      payload_dir = subgrn_payload_dir,
      payload_name = paste0(slug, "_pathway_subgrn"),
      condition_pathways = condition_pathways,
      topic_pathways = topic_pathways
    )
  }
  if (!nrow(subgrn_manifest) && !isTRUE(payload_disabled)) {
    subgrn_manifest <- .m3tb_reuse_existing_subgrn_manifest(subgrn_payload_dir, row, k)
  }
  if (nrow(subgrn_manifest)) {
    subgrn_manifest[, `:=`(
      method_setup = row$method_setup[[1L]],
      run_id = run_id,
      k = k,
      topic_space = topic_space,
      report_key = .m3tb_topic_space_value(row, "report_key", paste0("k", k))
    )]
  }
  condition_pathways_for_report <- condition_pathways

  tf_topic <- .m3tb_tf_topic_rows(theta, row$context_type[[1L]])
  condition_report <- data.table::data.table()
  condition_payload <- .m3cr_empty_payload_spec()
  group_mds <- data.table::as.data.table(score_result$mds_points)[
    as.integer(k) == as.integer(row$selected_k[[1L]]) &
      method_setup == row$method_setup[[1L]]
  ]
  if ("report_key" %in% names(group_mds) && "report_key" %in% names(row)) {
    group_mds <- group_mds[report_key == row$report_key[[1L]]]
  }
  if (nrow(group_mds)) {
    group_mds[, mds_label := .m3tb_short_label(display_label, 18L)]
    group_topic <- .m3tb_topic_waterfall(theta, row$context_type[[1L]])
    condition_label <- sprintf(
      "%s | %s condition topic view",
      row$method_setup[[1L]],
      .m3tb_topic_space_value(row, "report_label", sprintf("K%d", k))
    )
    condition_out <- file.path(condition_page_dir, sprintf("%s_condition_topic_report.html", slug))
    if (identical(as.character(row$context_type[[1L]]), "condition")) {
      condition_payload <- .m3cr_write_condition_payload(
        extraction_dir = extraction_dir,
        model_dir = row$model_dir[[1L]],
        payload_dir = condition_pair_payload_dir,
        payload_name = paste0(slug, "_condition_pair"),
        payload_base = "../../assets/condition_pair_payloads",
        topic_space = topic_space
      )
    }
    .m3tb_condition_report_html(
      condition_label,
      group_mds,
      group_topic,
      tf_topic,
      condition_pathways_for_report,
      condition_out,
      subgrn_manifest = subgrn_manifest,
      subgrn_payload_base = "../../pathway_subgrn_payloads",
      condition_payload = condition_payload,
      report_state = report_state
    )
    .m3tb_validate_current_report_html(condition_out)
    condition_report <- data.table::data.table(
      label = condition_label,
      k = k,
      trained_k = as.integer(.m3tb_topic_space_value(row, "trained_k", k)),
      topic_count = as.integer(.m3tb_topic_space_value(row, "topic_count", ncol(theta))),
      topic_space = topic_space,
      report_key = .m3tb_topic_space_value(row, "report_key", paste0("k", k)),
      report_label = .m3tb_topic_space_value(row, "report_label", paste0("K", k)),
      method_setup = row$method_setup[[1L]],
      run_id = run_id,
      path = condition_out
    )
  }

  list(
    topic = data.table::data.table(),
    condition = condition_report,
    subgrn_manifest = subgrn_manifest,
    condition_payload = data.table::data.table(
      method_setup = row$method_setup[[1L]],
      run_id = run_id,
      k = k,
      topic_space = topic_space,
      report_key = .m3tb_topic_space_value(row, "report_key", paste0("k", k)),
      payload_file = condition_payload$payload_file,
      network_payload_file = condition_payload$network_payload_file,
      n_tf_gene = condition_payload$n_tf_gene,
      n_tf_peak_gene = condition_payload$n_tf_peak_gene
    )
  )
}

.m3tb_validate_current_report_html <- function(path) {
  if (!file.exists(path)) .log_abort("Module 3 report was not written: {path}")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  schema3 <- grepl(
    "craftgrn-module3-report-schema\" content=\"3",
    html,
    fixed = TRUE
  )
  schema4 <- grepl(
    "craftgrn-module3-report-schema\" content=\"4",
    html,
    fixed = TRUE
  )
  schema5 <- grepl(
    "craftgrn-module3-report-schema\" content=\"5",
    html,
    fixed = TRUE
  )
  schema6 <- grepl(
    "craftgrn-module3-report-schema\" content=\"6",
    html,
    fixed = TRUE
  )
  required <- if (schema5 || schema6) {
    c(
      "Condition Probability",
      "TF Probability",
      "id=\"cond1Select\"",
      "id=\"cond2Select\"",
      "id=\"tfSelect\"",
      "id=\"pathwaySelect\"",
      "id=\"tfButterflySvg\"",
      "id=\"topicPageStatus\"",
      "id=\"tfPageStatus\"",
      "id=\"pathPageStatus\"",
      "id=\"networkTabFilter\"",
      "id=\"networkTabLayout\"",
      "id=\"networkTabAppearance\"",
      "role=\"dialog\"",
      "CONDITION_PAYLOAD"
    )
  } else if (schema4) {
    c(
      "Condition Probability",
      "TF/Condition Probability",
      "id=\"cond1Select\"",
      "id=\"cond2Select\"",
      "id=\"tfSelect\"",
      "id=\"pathwaySelect\"",
      "id=\"tfButterflySvg\"",
      "CONDITION_PAYLOAD"
    )
  } else if (schema3) {
    c(
      "Condition Topic Scores",
      "id=\"cond1Select\"",
      "id=\"cond2Select\"",
      "id=\"tfSelect\"",
      "id=\"pathwaySelect\"",
      "CONDITION_PAYLOAD"
    )
  } else {
    c(
      "craftgrn-module3-report-schema\" content=\"2",
      "Topic Bar Plots",
      "id=\"tfSearchInput\"",
      "id=\"tfSelect\"",
      "multiple size=\"3\"",
      "SUBGRN_MANIFEST"
    )
  }
  missing <- required[!vapply(required, grepl, logical(1L), x = html, fixed = TRUE)]
  obsolete <- grepl("Topic Waterfall", html, fixed = TRUE) ||
    ((schema3 || schema4 || schema5 || schema6) && grepl("Topic Bar Plots", html, fixed = TRUE))
  if (length(missing) || obsolete) {
    details <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (obsolete) "contains an obsolete panel label"
    )
    .log_abort("Generated Module 3 report does not match the current schema ({paste(details, collapse = '; ')}): {path}")
  }
  invisible(TRUE)
}

.m3tb_write_review_html <- function(output_dir, score_result, link_summary, review_dir = NULL) {
  output_dir <- normalizePath(as.character(output_dir)[[1L]], winslash = "/", mustWork = FALSE)
  force(score_result)
  force(link_summary)
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  review_dir <- normalizePath(as.character(review_dir)[[1L]], winslash = "/", mustWork = FALSE)
  html_dir <- .m3tb_review_html_dir(review_dir)
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  condition_report_dir <- file.path(html_dir, "condition_topic_reports")
  condition_page_dir <- file.path(condition_report_dir, "pages")
  subgrn_payload_dir <- file.path(html_dir, "pathway_subgrn_payloads")
  asset_dir <- .m3tb_plot_dir(review_dir)
  condition_pair_payload_dir <- file.path(asset_dir, "condition_pair_payloads")
  dir.create(condition_report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(condition_page_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(subgrn_payload_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(condition_pair_payload_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(html_dir, pattern = "^(topic_report|condition_topic_report)_K[0-9]+[.]html$", full.names = TRUE))
  unlink(list.files(condition_report_dir, pattern = "_K[0-9]+_condition_topic_report[.]html$", full.names = TRUE))
  unlink(list.files(condition_page_dir, pattern = "_K[0-9]+_condition_topic_report[.]html$", full.names = TRUE))
  obsolete <- c(
    file.path(html_dir, "theta_phi_and_group_mds.html"),
    file.path(html_dir, "topic_method_k_topic_mds_report.html"),
    file.path(html_dir, "topic_reports"),
    file.path(condition_page_dir, "assets"),
    file.path(review_dir, "tf_std_six_setups_pass_state_counts.pdf"),
    file.path(review_dir, "tf_std_six_setups_shared_topic_counts.pdf")
  )
  unlink(obsolete[file.exists(obsolete)], recursive = TRUE, force = TRUE)
  obsolete_tf_std <- list.files(
    review_dir,
    pattern = "^tf_std_six_setups_",
    recursive = TRUE,
    full.names = TRUE
  )
  unlink(obsolete_tf_std, recursive = TRUE, force = TRUE)
  obsolete_theta_assets <- list.files(
    asset_dir,
    pattern = "^(theta_group_mds|theta_phi_topic_distance_correlation)_",
    full.names = TRUE
  )
  unlink(obsolete_theta_assets, force = TRUE)
  subgrn_max_context_env <- Sys.getenv("CRAFTGRN_PATHWAY_SUBGRN_MAX_CONTEXTS", unset = "")
  subgrn_skip_payloads <- nzchar(subgrn_max_context_env) &&
    suppressWarnings(as.integer(subgrn_max_context_env[[1L]])) <= 0L
  if (!isTRUE(subgrn_skip_payloads)) {
    unlink(list.files(subgrn_payload_dir, pattern = "[.]js$", full.names = TRUE))
  }
  model_rows <- attr(score_result, "model_rows")
  condition_reports <- data.table::data.table()
  condition_page_reports <- data.table::data.table()
  subgrn_manifests <- list()
  if (is.data.frame(model_rows) && nrow(model_rows)) {
    workers <- .m3tb_review_worker_count(nrow(model_rows))
    page_rows <- .m3tb_review_lapply(seq_len(nrow(model_rows)), function(i) {
      tryCatch(
        .m3tb_build_review_pages_for_row(
          i = i,
          model_rows = model_rows,
          score_result = score_result,
          topic_page_dir = "",
          condition_page_dir = condition_page_dir,
          subgrn_payload_dir = subgrn_payload_dir,
          condition_pair_payload_dir = condition_pair_payload_dir
        ),
        error = function(e) .m3tb_review_row_error(i, e)
      )
    }, workers = workers)
    failures <- vapply(page_rows, inherits, logical(1L), "m3tb_review_row_error")
    if (any(failures)) {
      msg <- vapply(page_rows[failures], function(x) {
        sprintf("row %s: %s", x$index, x$message)
      }, character(1L))
      .log_abort(paste(c("Failed to build Module 3 review HTML pages.", msg), collapse = "\n"))
    }
    condition_reports <- data.table::rbindlist(lapply(page_rows, `[[`, "condition"), use.names = TRUE, fill = TRUE)
    condition_page_reports <- data.table::copy(condition_reports)
    subgrn_manifests <- lapply(page_rows, `[[`, "subgrn_manifest")
    condition_payload_manifest <- data.table::rbindlist(
      lapply(page_rows, `[[`, "condition_payload"),
      use.names = TRUE,
      fill = TRUE
    )
    if (nrow(condition_payload_manifest)) {
      csv_dir <- .m3tb_review_tables_dir(review_dir)
      dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(
        condition_payload_manifest,
        file.path(csv_dir, "condition_pair_payload_manifest.csv")
      )
    }
    report_paths <- condition_page_reports$path
    invisible(lapply(report_paths, .m3tb_validate_current_report_html))
    existing_pages <- list.files(
      condition_page_dir,
      pattern = "[.]html$",
      full.names = TRUE
    )
    stale_pages <- setdiff(
      normalizePath(existing_pages, winslash = "/", mustWork = FALSE),
      normalizePath(report_paths, winslash = "/", mustWork = FALSE)
    )
    unlink(stale_pages[file.exists(stale_pages)], force = TRUE)
  }
  if (length(subgrn_manifests)) {
    subgrn_manifest <- data.table::rbindlist(subgrn_manifests, use.names = TRUE, fill = TRUE)
    if (nrow(subgrn_manifest)) {
      csv_dir <- .m3tb_review_tables_dir(review_dir)
      dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(subgrn_manifest, file.path(csv_dir, "pathway_subgrn_manifest.csv"))
    }
  }
  if (nrow(condition_reports)) {
    condition_reports[, method_key := paste(run_id, method_setup, sep = "\t")]
    condition_reports <- data.table::rbindlist(lapply(split(condition_reports, condition_reports$method_key), function(rows) {
      data.table::setorder(rows, k)
      slug <- .safe_filename(paste(rows$run_id[[1L]], rows$method_setup[[1L]], sep = "_"))
      out <- file.path(condition_report_dir, sprintf("%s_condition_topic_report.html", slug))
      .m3tb_write_k_report_index(out, rows, sprintf("%s condition topic view", rows$method_setup[[1L]]))
      data.table::data.table(
        label = rows$method_setup[[1L]],
        method_setup = rows$method_setup[[1L]],
        run_id = rows$run_id[[1L]],
        path = out
      )
    }), use.names = TRUE, fill = TRUE)
  }
  files <- c(condition = file.path(html_dir, "topic_method_k_condition_mds_report.html"))
  stale_global <- file.path(
    html_dir,
    c(
      "topic_method_k_topic_mds_report_global_term_group.html",
      "topic_method_k_condition_mds_report_global_term_group.html"
    )
  )
  unlink(stale_global[file.exists(stale_global)])
  if (!nrow(condition_page_reports)) condition_page_reports <- condition_reports
  .m3cr_write_combined_report_index(
    files[["condition"]],
    condition_page_reports,
    report_state = attr(score_result, "report_state") %||% list()
  )
  files
}

.m3tb_write_review_outputs <- function(output_dir, score_result, link_summary, review_dir = NULL) {
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
  html <- .m3tb_write_review_html(output_dir, score_result, link_summary, review_dir = review_dir)
  invisible(html)
}

.module3_qc_detect_differential_links_dir <- function(topic_dir, differential_links_dir = NULL) {
  if (!is.null(differential_links_dir) && nzchar(as.character(differential_links_dir)[[1L]])) {
    path <- as.character(differential_links_dir)[[1L]]
    if (dir.exists(path)) return(path)
    return(NULL)
  }
  candidates <- unique(c(
    file.path(topic_dir, "differential_links"),
    file.path(dirname(topic_dir), "differential_links"),
    topic_dir
  ))
  hit <- candidates[dir.exists(candidates) & file.exists(file.path(candidates, "filtered_links_manifest.csv"))]
  if (length(hit)) return(hit[[1L]])
  hit <- candidates[dir.exists(candidates) & length(list.files(candidates, "_filtered_links_(up|down)[.]csv$", full.names = FALSE)) > 0L]
  if (length(hit)) hit[[1L]] else NULL
}

.module3_qc_filtered_link_files <- function(differential_links_dir) {
  if (is.null(differential_links_dir) || !dir.exists(differential_links_dir)) return(character())
  manifest <- file.path(differential_links_dir, "filtered_links_manifest.csv")
  if (file.exists(manifest)) {
    man <- data.table::fread(manifest, showProgress = FALSE)
    files <- unique(c(as.character(man$up_path), as.character(man$down_path)))
    files <- ifelse(file.exists(files), files, file.path(differential_links_dir, basename(files)))
    return(files[file.exists(files)])
  }
  list.files(differential_links_dir, "_filtered_links_(up|down)[.]csv$", full.names = TRUE)
}

.module3_qc_read_filtered_link_file <- function(path) {
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  cols <- intersect(
    c(
      "tf", "gene_key", "peak_id", "delta_fp_score", "log2FC_tf_expr",
      "comparison_id", "comparison_group", "cond1_id", "cond2_id",
      "cond1_label", "cond2_label"
    ),
    header
  )
  dt <- data.table::fread(path, select = cols, showProgress = FALSE)
  if (!("tf" %in% names(dt))) return(data.table::data.table())
  if (!("gene_key" %in% names(dt))) dt[, gene_key := NA_character_]
  if (!("delta_fp_score" %in% names(dt))) dt[, delta_fp_score := NA_real_]
  if (!("log2FC_tf_expr" %in% names(dt))) dt[, log2FC_tf_expr := NA_real_]
  direction <- if (grepl("_filtered_links_up[.]csv$", basename(path))) "up" else "down"
  comparison_id <- sub("_filtered_links_(up|down)[.]csv$", "", basename(path))
  if (!("comparison_id" %in% names(dt))) dt[, comparison_id := comparison_id]
  dt[, direction := direction]
  dt
}

.module3_qc_empty_tf_summary <- function() {
  data.table::data.table(
    comparison_id = character(),
    comparison_group = character(),
    cond1_id = character(),
    cond2_id = character(),
    cond1_label = character(),
    cond2_label = character(),
    TF = character(),
    n_links_up = integer(),
    n_links_down = integer(),
    n_target_genes_up = integer(),
    n_target_genes_down = integer(),
    net_target_gene_count = integer(),
    tf_delta_sum = numeric(),
    tf_delta_sum_abs = numeric(),
    median_log2FC_tf_expr = numeric(),
    dominant_direction = character(),
    rank = integer()
  )
}

.module3_qc_summarize_differential_tfs <- function(differential_links_dir,
                                                   output_dir,
                                                   top_n = 20L,
                                                   verbose = TRUE) {
  files <- .module3_qc_filtered_link_files(differential_links_dir)
  out_dir <- file.path(output_dir, "differential_tf_summary")
  per_dir <- file.path(out_dir, "per_comparison")
  dir.create(per_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(files)) {
    empty <- .module3_qc_empty_tf_summary()
    data.table::fwrite(empty, file.path(out_dir, "module3_top_differential_tfs.csv"))
    return(list(combined = empty, files = character(), output_dir = out_dir))
  }
  dt <- data.table::rbindlist(lapply(files, .module3_qc_read_filtered_link_file), use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) {
    empty <- .module3_qc_empty_tf_summary()
    data.table::fwrite(empty, file.path(out_dir, "module3_top_differential_tfs.csv"))
    return(list(combined = empty, files = character(), output_dir = out_dir))
  }
  for (col in c("comparison_group", "cond1_id", "cond2_id", "cond1_label", "cond2_label")) {
    if (!(col %in% names(dt))) dt[, (col) := NA_character_]
  }
  dt[, `:=`(
    TF = as.character(tf),
    gene_key = as.character(gene_key),
    delta_fp_score = suppressWarnings(as.numeric(delta_fp_score)),
    log2FC_tf_expr = suppressWarnings(as.numeric(log2FC_tf_expr))
  )]
  dt <- dt[!is.na(TF) & nzchar(TF)]
  if (!nrow(dt)) {
    empty <- .module3_qc_empty_tf_summary()
    data.table::fwrite(empty, file.path(out_dir, "module3_top_differential_tfs.csv"))
    return(list(combined = empty, files = character(), output_dir = out_dir))
  }
  gene_dt <- unique(dt[, .(comparison_id, TF, gene_key, direction)])
  count_dt <- gene_dt[, .(
    n_target_genes_up = data.table::uniqueN(gene_key[direction == "up" & !is.na(gene_key)]),
    n_target_genes_down = data.table::uniqueN(gene_key[direction == "down" & !is.na(gene_key)])
  ), by = .(comparison_id, TF)]
  link_dt <- dt[, .(
    comparison_group = .SD$comparison_group[[1L]],
    cond1_id = .SD$cond1_id[[1L]],
    cond2_id = .SD$cond2_id[[1L]],
    cond1_label = .SD$cond1_label[[1L]],
    cond2_label = .SD$cond2_label[[1L]],
    n_links_up = sum(direction == "up", na.rm = TRUE),
    n_links_down = sum(direction == "down", na.rm = TRUE),
    tf_delta_sum = sum(delta_fp_score, na.rm = TRUE),
    tf_delta_sum_abs = sum(abs(delta_fp_score), na.rm = TRUE),
    median_log2FC_tf_expr = stats::median(log2FC_tf_expr, na.rm = TRUE)
  ), by = .(comparison_id, TF)]
  out <- merge(link_dt, count_dt, by = c("comparison_id", "TF"), all.x = TRUE, sort = FALSE)
  out[is.na(n_target_genes_up), n_target_genes_up := 0L]
  out[is.na(n_target_genes_down), n_target_genes_down := 0L]
  out[, net_target_gene_count := n_target_genes_up - n_target_genes_down]
  out[, dominant_direction := data.table::fifelse(net_target_gene_count > 0, "up", data.table::fifelse(net_target_gene_count < 0, "down", "balanced"))]
  out[!is.finite(median_log2FC_tf_expr), median_log2FC_tf_expr := NA_real_]
  out[, abs_net_target_gene_count := abs(net_target_gene_count)]
  data.table::setorder(out, comparison_id, -abs_net_target_gene_count, -tf_delta_sum_abs, TF)
  out[, rank := seq_len(.N), by = comparison_id]
  out[, abs_net_target_gene_count := NULL]
  top_n <- max(1L, as.integer(top_n)[[1L]])
  combined <- out[rank <= top_n]
  combined_path <- file.path(out_dir, "module3_top_differential_tfs.csv")
  data.table::fwrite(combined, combined_path)
  per_files <- vapply(split(out, out$comparison_id), function(x) {
    path <- file.path(per_dir, paste0(.m3tb_safe_label(x$comparison_id[[1L]]), "_differential_tfs.csv"))
    data.table::fwrite(x, path)
    path
  }, character(1L))
  if (isTRUE(verbose)) {
    .log_inform("Wrote Module 3 differential TF summaries for {length(per_files)} comparison(s): {out_dir}")
  }
  list(combined = combined, files = unname(per_files), output_dir = out_dir, combined_file = combined_path)
}

.module3_qc_read_differential_summary <- function(differential_links_dir) {
  path <- if (!is.null(differential_links_dir)) file.path(differential_links_dir, "qc", "differential_link_summary.csv") else NULL
  if (!is.null(path) && file.exists(path)) return(data.table::fread(path, showProgress = FALSE))
  files <- .module3_qc_filtered_link_files(differential_links_dir)
  if (!length(files)) return(data.table::data.table())
  rows <- lapply(files, function(path) {
    header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
    n <- nrow(data.table::fread(path, select = intersect("tf", header), showProgress = FALSE))
    comparison_id <- sub("_filtered_links_(up|down)[.]csv$", "", basename(path))
    direction <- if (grepl("_filtered_links_up[.]csv$", basename(path))) "up" else "down"
    data.table::data.table(comparison_id = comparison_id, direction = direction, n = n)
  })
  dt <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  dcast <- data.table::dcast(dt, comparison_id ~ direction, value.var = "n", fill = 0)
  if (!("up" %in% names(dcast))) dcast[, up := 0L]
  if (!("down" %in% names(dcast))) dcast[, down := 0L]
  data.table::setnames(dcast, c("up", "down"), c("n_up", "n_down"))
  dcast
}

.module3_diff_grn_read_browser_file <- function(path, nrows = Inf) {
  header <- names(data.table::fread(path, nrows = 0L, showProgress = FALSE))
  cols <- intersect(
    c(
      "tf", "gene_key", "peak_id", "delta_link_score", "link_score_cond1",
      "link_score_cond2", "delta_fp_score", "delta_tf_expr",
      "delta_gene_expr", "log2FC_tf_expr", "log2FC_gene_expr",
      "distance_to_tss", "candidate_source", "r_gene", "r_rna_gene",
      "comparison_id", "comparison_group", "cond1_id", "cond2_id",
      "cond1_label", "cond2_label"
    ),
    header
  )
  if (is.finite(nrows)) {
    dt <- data.table::fread(
      path,
      select = cols,
      nrows = max(1L, as.integer(nrows)[[1L]]),
      showProgress = FALSE
    )
  } else {
    dt <- data.table::fread(path, select = cols, showProgress = FALSE)
  }
  if (!nrow(dt) || !all(c("tf", "gene_key") %chin% names(dt))) return(data.table::data.table())
  direction <- if (grepl("_filtered_links_up[.]csv$", basename(path))) "up" else "down"
  comparison_id <- sub("_filtered_links_(up|down)[.]csv$", "", basename(path))
  if (!("comparison_id" %in% names(dt))) dt[, comparison_id := comparison_id]
  dt[, direction := direction]
  for (col in c("comparison_group", "cond1_id", "cond2_id", "cond1_label", "cond2_label", "peak_id", "candidate_source")) {
    if (!(col %in% names(dt))) dt[, (col) := NA_character_]
  }
  for (col in c("delta_link_score", "link_score_cond1", "link_score_cond2", "delta_fp_score", "delta_tf_expr", "delta_gene_expr", "log2FC_tf_expr", "log2FC_gene_expr", "distance_to_tss", "r_gene", "r_rna_gene")) {
    if (!(col %in% names(dt))) dt[, (col) := NA_real_]
    dt[, (col) := suppressWarnings(as.numeric(get(col)))]
  }
  dt[, `:=`(
    tf = as.character(tf),
    gene_key = as.character(gene_key),
    peak_id = as.character(peak_id)
  )]
  dt[!is.na(tf) & nzchar(tf) & !is.na(gene_key) & nzchar(gene_key)]
}

.module3_diff_grn_aggregate_network_edges <- function(differential_links_dir,
                                                       browser_edge_cap_per_group = NULL,
                                                       tf_keep_by_comparison = NULL,
                                                       browser_max_rows_per_file = 50000L) {
  files <- .module3_qc_filtered_link_files(differential_links_dir)
  if (!length(files)) return(list(edges = data.table::data.table(), browser_edges = data.table::data.table(), manifest = data.table::data.table()))
  if (is.null(browser_edge_cap_per_group)) browser_edge_cap_per_group <- 300L
  browser_edge_cap_per_group <- max(1L, as.integer(browser_edge_cap_per_group)[[1L]])
  max_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  }
  min_abs_finite <- function(x) {
    x <- suppressWarnings(abs(as.numeric(x)))
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  }
  median_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x)) stats::median(x) else NA_real_
  }
  aggregate_one <- function(path) {
    dt <- .module3_diff_grn_read_browser_file(path, nrows = browser_max_rows_per_file)
    if (!nrow(dt)) return(data.table::data.table())
    if (!is.null(tf_keep_by_comparison)) {
      comp <- as.character(dt$comparison_id[[1L]])
      keep <- tf_keep_by_comparison[[comp]]
      if (!is.null(keep) && length(keep)) {
        dt <- dt[toupper(tf) %chin% keep]
        if (!nrow(dt)) return(data.table::data.table())
      }
    }
    dt[, `:=`(
      edge_score_row = abs(data.table::fcoalesce(delta_link_score, delta_fp_score, log2FC_gene_expr, 0)),
      from = toupper(tf),
      to = as.character(gene_key)
    )]
    out <- dt[, {
      best_i <- which.max(edge_score_row)
      if (!length(best_i) || !is.finite(edge_score_row[best_i])) best_i <- 1L
      data.table::data.table(
        comparison_group = comparison_group[[1L]],
        cond1_id = cond1_id[[1L]],
        cond2_id = cond2_id[[1L]],
        cond1_label = cond1_label[[1L]],
        cond2_label = cond2_label[[1L]],
        from = from[[1L]],
        to = to[[1L]],
        edge_class = paste0(direction[[1L]], "_differential"),
        edge_score = sum(delta_link_score, na.rm = TRUE),
        abs_edge_score = max_finite(edge_score_row),
        edge_r = max_finite(edge_score_row),
        n_supporting_links = data.table::uniqueN(peak_id[!is.na(peak_id) & nzchar(peak_id)]),
        best_peak_ID = peak_id[[best_i]],
        best_distance_to_tss = min_abs_finite(distance_to_tss),
        candidate_source = candidate_source[[best_i]],
        median_log2FC_tf_expr = median_finite(log2FC_tf_expr),
        median_log2FC_gene_expr = median_finite(log2FC_gene_expr),
        max_abs_r_gene = max_finite(abs(r_gene)),
        max_abs_r_rna_gene = max_finite(abs(r_rna_gene))
      )
    }, by = .(comparison_id, direction, tf, gene_key)]
    for (col in c("abs_edge_score", "edge_r", "best_distance_to_tss", "median_log2FC_tf_expr", "median_log2FC_gene_expr", "max_abs_r_gene", "max_abs_r_rna_gene")) {
      out[!is.finite(get(col)), (col) := NA_real_]
    }
    out[is.na(n_supporting_links), n_supporting_links := 0L]
    out[, tf_target_count := data.table::uniqueN(to), by = .(comparison_id, direction, from)]
    out[, tf_edge_score_sum := sum(abs_edge_score, na.rm = TRUE), by = .(comparison_id, direction, from)]
    data.table::setorder(out, comparison_id, direction, -tf_target_count, -tf_edge_score_sum, -abs_edge_score, from, to)
    tf_rank <- unique(out[, .(comparison_id, direction, from, tf_target_count, tf_edge_score_sum)])
    data.table::setorder(tf_rank, comparison_id, direction, -tf_target_count, -tf_edge_score_sum, from)
    tf_rank[, tf_rank := seq_len(.N), by = .(comparison_id, direction)]
    out <- merge(out, tf_rank[, .(comparison_id, direction, from, tf_rank)], by = c("comparison_id", "direction", "from"), all.x = TRUE, sort = FALSE)
    data.table::setorder(out, comparison_id, direction, -tf_target_count, -tf_edge_score_sum, -abs_edge_score, from, to)
    out[, edge_rank := seq_len(.N), by = .(comparison_id, direction)]
    out
  }
  edge_dt <- data.table::rbindlist(lapply(files, aggregate_one), use.names = TRUE, fill = TRUE)
  if (!nrow(edge_dt)) return(list(edges = data.table::data.table(), browser_edges = data.table::data.table(), manifest = data.table::data.table()))
  for (col in c("abs_edge_score", "edge_r", "best_distance_to_tss", "median_log2FC_tf_expr", "median_log2FC_gene_expr", "max_abs_r_gene", "max_abs_r_rna_gene")) {
    edge_dt[!is.finite(get(col)), (col) := NA_real_]
  }
  data.table::setorder(edge_dt, comparison_id, direction, -tf_target_count, -tf_edge_score_sum, -abs_edge_score, from, to)
  browser_dt <- edge_dt[edge_rank <= browser_edge_cap_per_group]
  manifest <- unique(edge_dt[, .(
    comparison_group,
    cond1_id,
    cond2_id,
    cond1_label,
    cond2_label,
    n_edges = .N,
    n_tfs = data.table::uniqueN(from),
    n_genes = data.table::uniqueN(to)
  ), by = .(comparison_id, direction)])
  data.table::setorder(manifest, comparison_id, direction)
  list(edges = edge_dt, browser_edges = browser_dt, manifest = manifest)
}

.module3_diff_grn_write_network_browser <- function(edges,
                                                    manifest,
                                                    out_html,
                                                    top_tf_n = 10L,
                                                    top_link_n = 300L,
                                                    default_direction = "up",
                                                    title = "Differential GRN browser") {
  dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)
  if (!nrow(edges)) {
    writeLines("<!doctype html><html><head><meta charset=\"utf-8\"/></head><body><b>No differential GRN edges to plot.</b></body></html>", out_html, useBytes = TRUE)
    return(out_html)
  }
  payload <- edges[, .(
    comparison_id,
    direction,
    from,
    to,
    edge_score,
    abs_edge_score,
    n_supporting_links,
    best_peak_ID,
    best_distance_to_tss,
    tf_target_count,
    tf_rank,
    edge_rank
  )]
  payload_json <- jsonlite::toJSON(payload, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null")
  manifest_json <- jsonlite::toJSON(manifest, dataframe = "rows", auto_unbox = TRUE, null = "null", na = "null")
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"/>",
    paste0("<title>", .m3tb_html_escape(title), "</title>"),
    "<style>",
    "body{margin:0;background:#f7f7f5;color:#111;font-family:Arial,Helvetica,sans-serif;font-weight:700}",
    ".wrap{max-width:min(calc((100vh - 245px) * 1.7778),calc(100vw - 24px));margin:0 auto;padding:12px 12px 16px 12px}",
    ".top{display:block;border-bottom:1px solid #d6d6d0;padding-bottom:10px;margin-bottom:10px}",
    "h1{font-size:21px;line-height:1.18;margin:0 0 6px 0}.meta{font-size:12px;line-height:1.35;color:#555;max-width:1250px;margin-bottom:10px}",
    ".controls{display:flex;gap:8px 10px;align-items:center;flex-wrap:wrap}.control{display:flex;gap:5px;align-items:center;font-size:12px;color:#333;white-space:nowrap}",
    "select,input{font:700 13px Arial,Helvetica,sans-serif;border:1px solid #aaa;background:#fff;color:#111;border-radius:3px;padding:6px 8px}input{width:74px}input[type=range]{width:96px;padding:0}input[type=color]{width:32px;height:30px;padding:2px}input[type=checkbox]{width:auto}",
    "#comparisonSelect{width:270px}#directionSelect{width:96px}#layoutSelect{width:116px}#paletteSelect{width:220px}#nodeSelect{width:170px}",
    "button{font:700 13px Arial,Helvetica,sans-serif;border:1px solid #777;background:#222;color:#fff;border-radius:3px;padding:7px 10px;cursor:pointer}",
    ".note{font-size:11px;color:#666;margin:0 0 10px 0}.canvas{position:relative;width:100%;aspect-ratio:16/9;max-height:calc(100vh - 245px);border:1px solid #d6d6d0;background:#fff;box-shadow:0 1px 2px rgba(0,0,0,0.04);overflow:hidden;cursor:grab}.canvas.panning{cursor:grabbing}",
    "svg{width:100%;height:100%;display:block;background:#fff}.edge{fill:none}.node{stroke:transparent;stroke-width:0;cursor:grab}.node:active{cursor:grabbing}.label{font-family:Arial,Helvetica,sans-serif;font-weight:700;fill:#fff;stroke:#666;stroke-width:2.4px;paint-order:stroke;dominant-baseline:middle;text-anchor:middle;pointer-events:none}.geneLabel{fill:#111;stroke:#fff;stroke-width:2.6px;paint-order:stroke;text-anchor:start}.selected{stroke:#d7263d;stroke-width:2.6}.tooltip{position:absolute;display:none;background:rgba(17,17,17,0.92);color:#fff;font:700 12px Arial,Helvetica,sans-serif;padding:7px 8px;border-radius:3px;pointer-events:none;max-width:390px;line-height:1.35}",
    "</style></head><body><div class=\"wrap\"><div class=\"top\"><div>",
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    "<div class=\"meta\">Select a comparison and direction, then choose how many top TFs and top TF-to-gene links to draw. Edges aggregate footprint-supported Module 3 differential links; tooltips retain footprint support and differential metrics.</div></div>",
    "<div class=\"controls\">",
    "<label class=\"control\">Comparison <select id=\"comparisonSelect\"></select></label>",
    "<label class=\"control\">Direction <select id=\"directionSelect\"><option value=\"up\">Up</option><option value=\"down\">Down</option></select></label>",
    paste0("<label class=\"control\">Top TFs <input id=\"topTfN\" type=\"number\" min=\"0\" value=\"", as.integer(top_tf_n)[[1L]], "\"/></label>"),
    paste0("<label class=\"control\">Top links <input id=\"topLinkN\" type=\"number\" min=\"0\" value=\"", as.integer(top_link_n)[[1L]], "\"/></label>"),
    "<label class=\"control\">Layout <select id=\"layoutSelect\"><option value=\"force\" selected>Force</option><option value=\"radial\">Radial</option><option value=\"columns\">Columns</option><option value=\"bipartite\">Bipartite</option><option value=\"hierarchy\">Hierarchy</option><option value=\"concentric\">Concentric</option><option value=\"circle\">Circle</option><option value=\"grid\">Grid</option><option value=\"spiral\">Spiral</option><option value=\"clustered\">Clustered</option></select></label>",
    "<label class=\"control\">Spacing <input id=\"spacingRange\" type=\"range\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/><input id=\"spacingValue\" type=\"number\" min=\"0.5\" max=\"2\" step=\"0.01\" value=\"1\"/></label>",
    "<label class=\"control\">TF box <input id=\"tfBoxMin\" type=\"number\" min=\"6\" max=\"80\" step=\"1\" value=\"8\"/><input id=\"tfBoxMax\" type=\"number\" min=\"8\" max=\"110\" step=\"1\" value=\"28\"/></label>",
    "<label class=\"control\">Labels <input id=\"showGeneLabels\" type=\"checkbox\" checked/></label><label class=\"control\">Arrows <input id=\"showArrows\" type=\"checkbox\" checked/></label>",
    "<label class=\"control\">Palette <select id=\"paletteSelect\"><option value=\"default\" selected>Default by direction</option><option value=\"npg\">NPG inspired</option><option value=\"aaas\">AAAS inspired</option><option value=\"nejm\">NEJM inspired</option><option value=\"lancet\">Lancet inspired</option><option value=\"jama\">JAMA inspired</option><option value=\"bmj\">BMJ inspired</option><option value=\"jco\">JCO inspired</option><option value=\"ucscgb\">UCSCGB inspired</option><option value=\"d3\">D3 inspired</option><option value=\"gephi\">Gephi inspired</option><option value=\"observable\">Observable inspired</option><option value=\"primer\">Primer inspired</option><option value=\"atlassian\">Atlassian inspired</option><option value=\"iterm\">iTerm inspired</option><option value=\"locuszoom\">LocusZoom inspired</option><option value=\"igv\">IGV inspired</option><option value=\"cosmic\">COSMIC inspired</option><option value=\"uchicago\">UChicago inspired</option><option value=\"startrek\">Star Trek inspired</option><option value=\"tron\">Tron inspired</option><option value=\"futurama\">Futurama inspired</option><option value=\"rickandmorty\">Rick and Morty inspired</option><option value=\"simpsons\">Simpsons inspired</option><option value=\"flatui\">Flat UI inspired</option><option value=\"frontiers\">Frontiers inspired</option><option value=\"gsea\">GSEA inspired</option><option value=\"bootstrap5\">Bootstrap 5 inspired</option><option value=\"material\">Material Design inspired</option><option value=\"tailwind\">Tailwind CSS inspired</option></select></label>",
    "<label class=\"control\">Select node <select id=\"nodeSelect\"></select></label><button id=\"resetButton\" type=\"button\">Reset layout</button><button id=\"exportSvgButton\" type=\"button\">Export SVG</button>",
    "</div></div><p class=\"note\" id=\"statsLine\"></p><div class=\"canvas\" id=\"canvas\"><div class=\"tooltip\" id=\"tooltip\"></div><svg id=\"networkSvg\" viewBox=\"0 0 1600 900\" role=\"img\" aria-label=\"Differential GRN network\"><defs><marker id=\"arrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"5\" markerHeight=\"5\" orient=\"auto\"><path d=\"M 0 0 L 10 5 L 0 10 z\" fill=\"#8a8a8a\" fill-opacity=\"0.55\"/></marker></defs><g id=\"viewLayer\"><g id=\"edgeLayer\"></g><g id=\"nodeLayer\"></g><g id=\"labelLayer\"></g></g></svg></div></div>",
    "<script>",
    paste0("const FULL_EDGES=", payload_json, ";"),
    paste0("const MANIFEST=", manifest_json, ";"),
    paste0("const DEFAULT_DIRECTION='", .m3tb_html_escape(default_direction), "';"),
    "const WIDTH=1600,HEIGHT=900,svg=document.getElementById('networkSvg'),canvas=document.getElementById('canvas'),viewLayer=document.getElementById('viewLayer'),edgeLayer=document.getElementById('edgeLayer'),nodeLayer=document.getElementById('nodeLayer'),labelLayer=document.getElementById('labelLayer'),comparisonSelect=document.getElementById('comparisonSelect'),directionSelect=document.getElementById('directionSelect'),topTfN=document.getElementById('topTfN'),topLinkN=document.getElementById('topLinkN'),layoutSelect=document.getElementById('layoutSelect'),spacingRange=document.getElementById('spacingRange'),spacingValue=document.getElementById('spacingValue'),tfBoxMin=document.getElementById('tfBoxMin'),tfBoxMax=document.getElementById('tfBoxMax'),showGeneLabels=document.getElementById('showGeneLabels'),showArrows=document.getElementById('showArrows'),paletteSelect=document.getElementById('paletteSelect'),nodeSelect=document.getElementById('nodeSelect'),statsLine=document.getElementById('statsLine'),tooltip=document.getElementById('tooltip');",
    "const COLOR_PRESETS={default_up:['#B2182B','#2166AC','#A0A0A0','#B2182B'],default_down:['#2166AC','#B2182B','#A0A0A0','#2166AC'],npg:['#E64B35','#4DBBD5','#7E6148','#3C5488'],aaas:['#EE0000','#3B4992','#008B45','#631879'],nejm:['#BC3C29','#0072B5','#E18727','#20854E'],lancet:['#ED0000','#00468B','#42B540','#0099B4'],jama:['#374E55','#DF8F44','#B24745','#6A6599'],bmj:['#2B8CBE','#E34A33','#969696','#636363'],jco:['#CD534C','#0073C2','#868686','#EFC000'],ucscgb:['#FF0000','#0000FF','#999999','#666666'],d3:['#FF7F0E','#1F77B4','#7F7F7F','#2CA02C'],gephi:['#FF7F00','#377EB8','#999999','#4DAF4A'],observable:['#4269D0','#EF553B','#AAAAAA','#6CC5B0'],primer:['#CF222E','#0969DA','#8C959F','#57606A'],atlassian:['#FF5630','#0052CC','#97A0AF','#42526E'],iterm:['#CC241D','#458588','#A89984','#665C54'],locuszoom:['#D7191C','#2C7BB6','#ABD9E9','#7570B3'],igv:['#E41A1C','#377EB8','#999999','#4DAF4A'],cosmic:['#D73027','#4575B4','#969696','#542788'],uchicago:['#800000','#155F83','#8A9045','#767676'],startrek:['#CC0C00','#5C88DA','#9C9C9C','#FFCC00'],tron:['#FF410D','#0085CA','#6EE2FF','#F7C530'],futurama:['#FF6F00','#1B9E77','#A6A6A6','#7570B3'],rickandmorty:['#FAFD7C','#00B0C8','#B2DF8A','#808080'],simpsons:['#FED439','#709AE1','#8A9197','#D2AF81'],flatui:['#E74C3C','#3498DB','#95A5A6','#2C3E50'],frontiers:['#E64B35','#4DBBD5','#B09C85','#3C5488'],gsea:['#E31A1C','#1F78B4','#BDBDBD','#33A02C'],bootstrap5:['#DC3545','#0D6EFD','#6C757D','#198754'],material:['#F44336','#2196F3','#9E9E9E','#4CAF50'],tailwind:['#EF4444','#3B82F6','#9CA3AF','#10B981']};let state={nodes:[],edges:[],drag:null,pan:null,selected:'',view:{x:0,y:0,k:1}};",
    "function esc(s){return String(s==null?'':s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;')}function el(n,a){const x=document.createElementNS('http://www.w3.org/2000/svg',n);Object.entries(a||{}).forEach(([k,v])=>x.setAttribute(k,v));return x}function num(x,d){const v=Number(x);return Number.isFinite(v)?v:d}function limitCount(input,fallback){const v=Math.floor(num(input.value,fallback));return v===0?Infinity:Math.max(1,v)}function spacing(){return Math.max(.5,Math.min(2,num(spacingValue.value,1)))}function colors(){const p=paletteSelect.value;if(p==='default')return COLOR_PRESETS[directionSelect.value==='down'?'default_down':'default_up'];return COLOR_PRESETS[p]||COLOR_PRESETS.default_up}function hexRgb(hex){const x=String(hex||'#777777').replace('#','');const v=parseInt(x,16);return{r:(v>>16)&255,g:(v>>8)&255,b:v&255}}function mixRgb(a,b,t){return 'rgb('+Math.round(a.r+(b.r-a.r)*t)+','+Math.round(a.g+(b.g-a.g)*t)+','+Math.round(a.b+(b.b-a.b)*t)+')'}function edgeColor(v,lo,hi){const t=hi>lo?(v-lo)/(hi-lo):.8;return mixRgb({r:226,g:226,b:226},hexRgb(colors()[3]),Math.max(0,Math.min(1,t)))}function selectedRows(){const cmp=comparisonSelect.value,dir=directionSelect.value,tfN=limitCount(topTfN,10),linkN=limitCount(topLinkN,300),rows=FULL_EDGES.filter(e=>e.comparison_id===cmp&&e.direction===dir&&num(e.tf_rank,999999)<=tfN).sort((a,b)=>num(b.abs_edge_score,0)-num(a.abs_edge_score,0)||String(a.from).localeCompare(String(b.from))||String(a.to).localeCompare(String(b.to)));return Number.isFinite(linkN)?rows.slice(0,linkN):rows}",
    "function buildGraph(){const edges=selectedRows();const nodeMap=new Map();edges.forEach(e=>{if(!nodeMap.has(e.from))nodeMap.set(e.from,{id:e.from,node_type:'TF',shared_targets:0,score:0});if(!nodeMap.has(e.to))nodeMap.set(e.to,{id:e.to,node_type:'Gene',shared_targets:0,score:0});const tf=nodeMap.get(e.from);tf.shared_targets+=1;tf.score+=num(e.abs_edge_score,0)});return{nodes:Array.from(nodeMap.values()),edges:edges.map(e=>Object.assign({},e,{edge_r:num(e.abs_edge_score,0)})),totalEdges:FULL_EDGES.filter(e=>e.comparison_id===comparisonSelect.value&&e.direction===directionSelect.value).length}}",
    "function sortedTf(ns){return ns.filter(n=>n.node_type==='TF').sort((a,b)=>num(b.shared_targets,0)-num(a.shared_targets,0)||a.id.localeCompare(b.id))}function sortedGenes(ns){return ns.filter(n=>n.node_type==='Gene').sort((a,b)=>a.id.localeCompare(b.id))}function setColumn(items,x,top,bottom){const desired=(items.some(n=>n.node_type==='TF')?86:64)*spacing(),gap=Math.min(desired,(bottom-top)/Math.max(items.length-1,1));items.forEach((n,i)=>{n.x=x;n.y=(top+bottom)/2-gap*(items.length-1)/2+i*gap})}function placeRing(items,ring,a0){items.forEach((n,i)=>{const a=a0+2*Math.PI*i/Math.max(items.length,1);n.x=WIDTH/2+Math.cos(a)*ring*spacing();n.y=HEIGHT/2+Math.sin(a)*ring*spacing()})}function radialLayout(ns){placeRing(sortedTf(ns),210,-Math.PI/2);placeRing(sortedGenes(ns),360,-Math.PI/2)}function columnsLayout(ns){setColumn(sortedTf(ns),WIDTH/2-330*spacing(),70,830);setColumn(sortedGenes(ns),WIDTH/2+330*spacing(),70,830)}function forceLayout(ns,es){radialLayout(ns);const by=new Map(ns.map(n=>[n.id,n]));for(let iter=0;iter<180;iter++){ns.forEach(n=>{n.vx=n.vx||0;n.vy=n.vy||0});for(let i=0;i<ns.length;i++)for(let j=i+1;j<ns.length;j++){const a=ns[i],b=ns[j];let dx=a.x-b.x,dy=a.y-b.y,d2=Math.max(dx*dx+dy*dy,64),d=Math.sqrt(d2),f=Math.min(3,2500*spacing()*spacing()/d2);dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f}es.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;let dx=b.x-a.x,dy=b.y-a.y,d=Math.max(Math.sqrt(dx*dx+dy*dy),1),target=120*spacing(),f=(d-target)*.012;dx/=d;dy/=d;a.vx+=dx*f;a.vy+=dy*f;b.vx-=dx*f;b.vy-=dy*f});ns.forEach(n=>{n.vx+=(WIDTH/2-n.x)*.002;n.vy+=(HEIGHT/2-n.y)*.002;n.x=Math.max(35,Math.min(WIDTH-35,n.x+n.vx));n.y=Math.max(35,Math.min(HEIGHT-35,n.y+n.vy));n.vx*=.76;n.vy*=.76})}ns.forEach(n=>{delete n.vx;delete n.vy})}function layout(ns,es){const m=layoutSelect.value;if(m==='columns'||m==='bipartite'||m==='hierarchy')columnsLayout(ns);else if(m==='radial'||m==='concentric'||m==='circle'||m==='grid'||m==='spiral'||m==='clustered')radialLayout(ns);else forceLayout(ns,es)}",
    "function radiusScale(ns){const vals=ns.filter(n=>n.node_type==='TF').map(n=>num(n.shared_targets,1)),lo=Math.min(...vals,1),hi=Math.max(...vals,1),b0=num(tfBoxMin.value,8),b1=num(tfBoxMax.value,28);return n=>{if(n.node_type==='Gene')return 6;if(hi<=lo)return 16;const t=(num(n.shared_targets,1)-lo)/(hi-lo);return b0+Math.pow(Math.max(0,Math.min(1,t)),.85)*(b1-b0)}}function lineEnd(a,b){const dx=b.x-a.x,dy=b.y-a.y,len=Math.max(Math.sqrt(dx*dx+dy*dy),1),r=b.node_type==='TF'?b.r*1.2:b.r;return{x:b.x-dx/len*r,y:b.y-dy/len*r}}function labelPos(n){if(n.node_type==='Gene')return{x:n.x+n.r+7,y:n.y,anchor:'start'};return{x:n.x,y:n.y,anchor:'middle'}}function applyView(){viewLayer.setAttribute('transform','translate('+state.view.x+' '+state.view.y+') scale('+state.view.k+')')}function screenPoint(ev){const p=svg.createSVGPoint();p.x=ev.clientX;p.y=ev.clientY;return p.matrixTransform(svg.getScreenCTM().inverse())}function screenToWorld(ev){const p=screenPoint(ev);return{x:(p.x-state.view.x)/state.view.k,y:(p.y-state.view.y)/state.view.k}}",
    "function render(){const g=buildGraph();state.nodes=g.nodes;state.edges=g.edges;const rFor=radiusScale(state.nodes);state.nodes.forEach(n=>{n.r=rFor(n)});layout(state.nodes,state.edges);const by=new Map(state.nodes.map(n=>[n.id,n])),scores=state.edges.map(e=>num(e.abs_edge_score,0)),lo=Math.min(...scores,0),hi=Math.max(...scores,1);edgeLayer.replaceChildren();nodeLayer.replaceChildren();labelLayer.replaceChildren();nodeSelect.replaceChildren();state.edges.forEach(e=>{const a=by.get(e.from),b=by.get(e.to);if(!a||!b)return;const end=lineEnd(a,b),line=el('line',{class:'edge',x1:a.x,y1:a.y,x2:end.x,y2:end.y,stroke:edgeColor(num(e.abs_edge_score,0),lo,hi),'stroke-width':1.15,'stroke-opacity':.48,'data-from':e.from,'data-to':e.to});if(showArrows.checked)line.setAttribute('marker-end','url(#arrow)');line.dataset.title=e.from+' -> '+e.to+'\\ncomparison: '+e.comparison_id+'\\ndirection: '+e.direction+'\\nsupporting footprints: '+e.n_supporting_links+'\\nbest footprint: '+(e.best_peak_ID||'NA')+'\\nbest distance to TSS: '+(e.best_distance_to_tss==null?'NA':Math.round(e.best_distance_to_tss))+'\\nedge score: '+num(e.abs_edge_score,0).toFixed(3);edgeLayer.appendChild(line)});state.nodes.forEach(n=>{const fill=n.node_type==='TF'?colors()[1]:colors()[2];let c;if(n.node_type==='TF'){c=el('rect',{class:'node',x:n.x-n.r*1.35,y:n.y-n.r*.78,width:n.r*2.7,height:n.r*1.56,rx:3,ry:3,fill:fill})}else{c=el('circle',{class:'node',cx:n.x,cy:n.y,r:n.r,fill:fill})}c.dataset.id=n.id;c.dataset.title=n.id+'\\n'+n.node_type+'\\nselected targets: '+num(n.shared_targets,0);nodeLayer.appendChild(c);if(n.node_type!=='Gene'||showGeneLabels.checked){const p=labelPos(n),t=el('text',{class:'label '+(n.node_type==='Gene'?'geneLabel':''),x:p.x,y:p.y,'text-anchor':p.anchor,'font-size':n.node_type==='TF'?Math.max(9,Math.min(19,n.r*.48)):9});t.dataset.id=n.id;t.textContent=n.id;labelLayer.appendChild(t)}});state.nodes.slice().sort((a,b)=>a.id.localeCompare(b.id)).forEach(n=>{const o=document.createElement('option');o.value=n.id;o.textContent=n.id;nodeSelect.appendChild(o)});statsLine.textContent='Showing '+state.nodes.length+' nodes and '+state.edges.length+'/'+g.totalEdges+' embedded TF-to-gene edges for '+comparisonSelect.value+' '+directionSelect.value+'. Edges aggregate footprint-supported differential GRN links.';applyView()}function redrawPositions(){const by=new Map(state.nodes.map(n=>[n.id,n]));edgeLayer.querySelectorAll('.edge').forEach(line=>{const a=by.get(line.dataset.from),b=by.get(line.dataset.to);if(!a||!b)return;const end=lineEnd(a,b);line.setAttribute('x1',a.x);line.setAttribute('y1',a.y);line.setAttribute('x2',end.x);line.setAttribute('y2',end.y)});nodeLayer.querySelectorAll('.node').forEach(shape=>{const n=by.get(shape.dataset.id);if(!n)return;if(shape.tagName.toLowerCase()==='rect'){shape.setAttribute('x',n.x-n.r*1.35);shape.setAttribute('y',n.y-n.r*.78)}else{shape.setAttribute('cx',n.x);shape.setAttribute('cy',n.y)}});labelLayer.querySelectorAll('text').forEach(t=>{const n=by.get(t.dataset.id);if(!n)return;const p=labelPos(n);t.setAttribute('x',p.x);t.setAttribute('y',p.y);t.setAttribute('text-anchor',p.anchor)})}",
    "function init(){Array.from(new Set(MANIFEST.map(d=>d.comparison_id))).sort().forEach(id=>{const o=document.createElement('option');o.value=id;o.textContent=id;comparisonSelect.appendChild(o)});if(comparisonSelect.options.length)comparisonSelect.value=comparisonSelect.options[0].value;directionSelect.value=DEFAULT_DIRECTION;[comparisonSelect,directionSelect,topTfN,topLinkN,layoutSelect,tfBoxMin,tfBoxMax,showGeneLabels,showArrows,paletteSelect].forEach(x=>{x.addEventListener('change',render);x.addEventListener('input',render)});spacingRange.addEventListener('input',()=>{spacingValue.value=spacingRange.value;render()});spacingValue.addEventListener('change',()=>{spacingRange.value=spacingValue.value;render()});document.getElementById('resetButton').addEventListener('click',()=>{state.view={x:0,y:0,k:1};render()});document.getElementById('exportSvgButton').addEventListener('click',exportSvg);render()}function inlineSvgStyles(clone){clone.querySelectorAll('line').forEach(x=>{x.setAttribute('fill','none');x.setAttribute('stroke-linecap','round')});clone.querySelectorAll('rect,circle').forEach(x=>{x.setAttribute('stroke','none')});Array.from(clone.querySelectorAll('text')).forEach(x=>{const isGene=String(x.getAttribute('class')||'').includes('geneLabel'),fs=Number(x.getAttribute('font-size')||11);x.setAttribute('font-family','Arial, Helvetica, sans-serif');x.setAttribute('font-weight','700');x.setAttribute('dominant-baseline','alphabetic');x.setAttribute('y',Number(x.getAttribute('y')||0)+fs*0.34);x.removeAttribute('paint-order');const halo=x.cloneNode(true);halo.setAttribute('fill','none');halo.setAttribute('stroke',isGene?'#fff':'#666');halo.setAttribute('stroke-width',isGene?'3':'2.8');halo.setAttribute('stroke-linejoin','round');halo.setAttribute('stroke-linecap','round');x.parentNode.insertBefore(halo,x);x.setAttribute('fill',isGene?'#111':'#fff');x.setAttribute('stroke','none')})}function exportSvg(){const clone=svg.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('width',WIDTH);clone.setAttribute('height',HEIGHT);inlineSvgStyles(clone);const bg=document.createElementNS('http://www.w3.org/2000/svg','rect');bg.setAttribute('x',0);bg.setAttribute('y',0);bg.setAttribute('width',WIDTH);bg.setAttribute('height',HEIGHT);bg.setAttribute('fill','#fff');clone.insertBefore(bg,clone.firstChild);const text=new XMLSerializer().serializeToString(clone),blob=new Blob([text],{type:'image/svg+xml'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=(document.title||'differential_grn_network').replace(/[^A-Za-z0-9_.-]+/g,'_')+'.svg';document.body.appendChild(a);a.click();a.remove();URL.revokeObjectURL(a.href)}",
    "nodeSelect.addEventListener('change',()=>{state.selected=nodeSelect.value;nodeLayer.querySelectorAll('.node').forEach(c=>c.classList.toggle('selected',c.dataset.id===state.selected))});svg.addEventListener('wheel',ev=>{ev.preventDefault();const p=screenPoint(ev),old=state.view.k,next=Math.max(.25,Math.min(6,old*(ev.deltaY<0?1.15:1/1.15)));state.view.x=p.x-(p.x-state.view.x)*next/old;state.view.y=p.y-(p.y-state.view.y)*next/old;state.view.k=next;applyView()},{passive:false});nodeLayer.addEventListener('mousedown',ev=>{if(!(ev.target.classList&&ev.target.classList.contains('node')))return;ev.preventDefault();ev.stopPropagation();const n=state.nodes.find(x=>x.id===ev.target.dataset.id);if(n){state.drag=n;nodeSelect.value=n.id;nodeSelect.dispatchEvent(new Event('change'))}});svg.addEventListener('mousedown',ev=>{if(ev.target.classList&&ev.target.classList.contains('node'))return;const p=screenPoint(ev);state.pan={startX:p.x,startY:p.y,x:state.view.x,y:state.view.y};canvas.classList.add('panning')});svg.addEventListener('mousemove',ev=>{if(state.drag){const p=screenToWorld(ev);state.drag.x=Math.max(30,Math.min(WIDTH-30,p.x));state.drag.y=Math.max(30,Math.min(HEIGHT-30,p.y));redrawPositions()}else if(state.pan){const p=screenPoint(ev);state.view.x=state.pan.x+(p.x-state.pan.startX);state.view.y=state.pan.y+(p.y-state.pan.startY);applyView()}tooltip.style.left=(ev.offsetX+12)+'px';tooltip.style.top=(ev.offsetY+12)+'px'});svg.addEventListener('mouseup',()=>{state.drag=null;state.pan=null;canvas.classList.remove('panning')});svg.addEventListener('mouseleave',()=>{state.drag=null;state.pan=null;canvas.classList.remove('panning');tooltip.style.display='none'});svg.addEventListener('mouseover',ev=>{const title=ev.target.dataset?ev.target.dataset.title:'';if(!title)return;tooltip.innerHTML=esc(title).replace(/\\n/g,'<br/>');tooltip.style.display='block'});svg.addEventListener('mouseout',ev=>{if(ev.target.dataset&&ev.target.dataset.title)tooltip.style.display='none'});init();",
    "</script></body></html>"
  )
  writeLines(html, out_html, useBytes = TRUE)
  out_html
}

#' Run a Module 3 topic benchmark and review report
#'
#' Runs or reviews Module 3 topic models using the package-native benchmark
#' layout. The function can train/extract models, or it can score and report
#' existing theta/phi and topic-link outputs when `run_training = FALSE` and
#' `run_extraction = FALSE`.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   files, or condition-link files when `input_source = "condition_links"`.
#' @param multiomic_data Optional CraftGRN multiomic object. Required when
#'   `replicate_documents = TRUE`.
#' @param comparisons Comparison or condition grouping table, or a CSV path.
#'   For condition separation scoring, use columns `condition_label` and
#'   `condition_group`.
#' @param output_dir Topic benchmark output directory.
#' @param methods Character vector of benchmark method IDs, `"default"`, or
#'   `"all"`.
#' @param training_methods Optional subset of `methods` to train. `NULL` trains
#'   every selected method when `run_training = TRUE`.
#' @param extraction_methods Optional subset of `methods` to extract. `NULL`
#'   extracts every selected method when `run_extraction = TRUE`.
#' @param k_grid Integer topic numbers.
#' @param output_layout Output folder layout. `"standard"` is the clean
#'   single-method layout for regular use, `"benchmark"` writes shallow
#'   `run_001_*` method folders plus a top-level review folder, `"legacy"` uses
#'   the historical nested benchmark layout, and `"auto"` picks legacy when
#'   existing legacy outputs are being reviewed, otherwise standard for one
#'   method and benchmark for method grids.
#' @param replicate_documents Whether theta document labels are replicate
#'   resolved. When `TRUE`, score files use the
#'   `theta_condition_replicate_separation` prefix and condition replicates are
#'   grouped by the condition label after removing a trailing replicate suffix.
#' @param reuse_if_exists Reuse existing model outputs where possible.
#' @param local_threads Optional thread count for model training.
#' @param warplda_iterations Number of native WarpLDA iterations.
#' @param vae_python Optional Python executable for VAE training.
#' @param vae_epochs Number of VAE training epochs.
#' @param vae_batch_size VAE mini-batch size.
#' @param vae_hidden VAE hidden-layer width.
#' @param vae_lr VAE learning rate.
#' @param vae_seed VAE random seed.
#' @param vae_device VAE device, for example `"auto"`, `"cpu"`, or `"cuda"`.
#'   `"auto"` uses CUDA when PyTorch can access it and otherwise uses CPU.
#' @param input_source Input source type passed to topic training. Use
#'   `"differential_links"` for comparison links and `"condition_links"` for
#'   condition-native links.
#' @param sample_subset Optional condition/sample labels passed to the Module 3
#'   training engine. When supplied, only comparisons whose condition labels are
#'   both in this vector are used.
#' @param analysis_label Label used to name the topic-model analysis.
#' @param extraction_topic_report_args Optional named list of topic-extraction
#'   report argument overrides.
#' @param memory_safety Module 3 memory policy.
#' @param memory_max_fraction Maximum fraction of currently available memory
#'   that a Module 3 stage may use.
#' @param extraction_k_workers Number of K values to extract concurrently.
#'   Use `1` for same-process WarpLDA training and extraction.
#' @param extraction_k_max_workers Maximum concurrent K extraction workers.
#' @param memory_safety Module 3 memory policy. `"strict"` fails before an
#'   unsafe or unmeasurable allocation.
#' @param memory_max_fraction Maximum fraction of currently available memory
#'   that a Module 3 stage may use. Defaults to `0.8`.
#' @param extraction_k_workers Number of K values to extract concurrently.
#'   `NULL` selects workers adaptively from memory and CPU headroom.
#' @param extraction_k_max_workers Maximum concurrent K extraction workers.
#' @param extraction_k_memory_gb Optional conservative memory estimate per K
#'   worker in GiB. `NULL` estimates from the shared `edges_docs.rds` cache.
#' @param extraction_k_memory_reserve_gb Minimum available RAM in GiB to keep
#'   unused while scheduling K workers. At least 25 percent is also reserved.
#' @param memory_safety Module 3 memory policy passed to model training.
#' @param memory_max_fraction Maximum fraction of currently available memory
#'   that a Module 3 stage may use.
#' @param paper_benchmark_csv Optional Tsao 2022-style benchmark CSV used to
#'   score fibroblast paper-topic recovery from existing theta and phi files.
#' @param run_fibroblast_paper_recovery Whether to write the fibroblast
#'   paper-topic recovery review when `paper_benchmark_csv` is supplied.
#' @param run_training Train topic models before reporting.
#' @param run_extraction Run topic extraction before reporting.
#' @param run_reports Build score tables and review reports.
#' @param path_length_check Whether to scan output paths and warn when any are
#'   longer than `path_length_limit`.
#' @param path_length_limit Maximum displayed path length before warning. The
#'   default is conservative for Windows/Samba clients.
#' @param path_length_alias Optional replacement root used only for path-length
#'   display and counting, for example `"Z:/episcope_test/project"`.
#' @param verbose Emit concise progress messages.
#'
#' @return A list with method plan, discovered model rows, score results, and
#'   review report paths.
#' @noRd
run_module3_topic_benchmark <- function(filtered_dir,
                                        multiomic_data = NULL,
                                        comparisons,
                                        output_dir,
                                        methods = "comparison_aggr_multivi",
                                        training_methods = NULL,
                                        extraction_methods = NULL,
                                        k_grid = 10L,
                                        output_layout = c("auto", "standard", "benchmark", "legacy"),
                                        replicate_documents = FALSE,
                                        reuse_if_exists = TRUE,
                                        local_threads = NULL,
                                        warplda_iterations = 2000L,
                                        count_method = c("log", "bin"),
                                        count_input = NULL,
                                        vae_python = NULL,
                                        vae_epochs = 200L,
                                        vae_batch_size = 64L,
                                        vae_hidden = 128L,
                                        vae_lr = 1e-3,
                                        vae_seed = 123L,
                                        vae_device = "auto",
                                        input_source = c("differential_links", "condition_links"),
                                        sample_subset = NULL,
                                        analysis_label = NULL,
                                        report_state = list(),
                                        extraction_topic_report_args = list(),
                                        extraction_k_workers = NULL,
                                        extraction_k_max_workers = 4L,
                                        extraction_k_memory_gb = NULL,
                                        extraction_k_memory_reserve_gb = 32,
                                        memory_safety = c("strict", "adaptive", "off"),
                                        memory_max_fraction = 0.8,
                                        paper_benchmark_csv = NULL,
                                        run_fibroblast_paper_recovery = !is.null(paper_benchmark_csv),
                                        run_training = TRUE,
                                        run_extraction = TRUE,
                                        run_reports = TRUE,
                                        path_length_check = TRUE,
                                        path_length_limit = 240L,
                                        path_length_alias = NULL,
                                        verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[is.finite(k_grid) & k_grid > 1L]
  if (!length(k_grid)) .log_abort("k_grid must include at least one integer greater than 1.")
  count_method <- match.arg(count_method)
  input_source <- match.arg(input_source)
  memory_safety <- match.arg(memory_safety)
  count_input_effective <- .resolve_topic_count_input(count_method = count_method, count_input = count_input)
  method_plan <- .module3_topic_method_plan(methods = methods, k_grid = k_grid)
  normalize_stage_methods <- function(x, label) {
    if (is.null(x)) return(method_plan$method)
    x <- unique(as.character(x))
    bad <- setdiff(x, method_plan$method)
    if (length(bad)) .log_abort("{label} must be a subset of methods; unknown: {paste(bad, collapse = ', ')}")
    x
  }
  training_methods <- normalize_stage_methods(training_methods, "training_methods")
  extraction_methods <- normalize_stage_methods(extraction_methods, "extraction_methods")
  output_layout <- .m3tb_resolve_output_layout(
    output_layout = output_layout,
    output_dir = output_dir,
    method_plan = method_plan,
    k_grid = k_grid,
    run_training = run_training,
    run_extraction = run_extraction
  )
  method_plan <- .m3tb_apply_output_layout(method_plan, output_dir, output_layout)
  review_dir <- .m3tb_review_dir(output_dir, output_layout)
  if (!identical(output_layout, "legacy")) {
    .m3tb_clean_stale_review_layout(review_dir)
  }
  csv_dir <- .m3tb_review_tables_dir(review_dir)
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(method_plan, file.path(csv_dir, "module3_topic_method_plan.csv"))
  if (identical(output_layout, "benchmark")) {
    run_cols <- c(
      "run_id",
      "run_slug",
      "method",
      "method_order",
      "context_type",
      "fp_mode",
      "backend",
      "vae_variant",
      "method_setup",
      "run_dir",
      "topic_documents_dir",
      "topic_models_dir",
      "topic_extraction_dir"
    )
    data.table::fwrite(
      method_plan[, run_cols, with = FALSE],
      file.path(output_dir, "runs.csv")
    )
  }
  if (isTRUE(verbose)) {
    .log_inform("Module 3 topic benchmark methods: {nrow(method_plan)}.")
    .log_inform("Module 3 output layout: {output_layout}.")
  }
  if (isTRUE(run_training)) {
    training_plan <- method_plan[method %in% training_methods]
    for (i in seq_len(nrow(training_plan))) {
      row <- training_plan[i]
      if (isTRUE(verbose)) .log_inform("Training Module 3 topic method: {row$method_setup}.")
      train_topic_models(
        Kgrid = k_grid,
        input_dir = filtered_dir,
        output_dir = row$topic_models_dir[[1L]],
        input_source = input_source,
        sample_subset = sample_subset,
        analysis_label = analysis_label,
        doc_mode = "tf",
        doc_design = row$doc_design[[1L]],
        fp_term_mode = row$fp_mode[[1L]],
        gene_term_mode = "unique",
        include_tf_terms = TRUE,
        count_method = count_method,
        count_input = count_input_effective,
        backend = row$backend[[1L]],
        vae_variant = row$vae_variant[[1L]],
        reuse_if_exists = reuse_if_exists,
        local_threads = local_threads,
        memory_safety = memory_safety,
        memory_max_fraction = memory_max_fraction,
        warplda_iterations = warplda_iterations,
        vae_python = vae_python,
        vae_epochs = vae_epochs,
        vae_batch_size = vae_batch_size,
        vae_hidden = vae_hidden,
        vae_lr = vae_lr,
        vae_seed = vae_seed,
        vae_device = vae_device,
        save_full_doc_term_csv = FALSE,
        flat_output = .m3tb_flat_output_for_layout(output_layout),
        topic_report_args = list()
      )
    }
  }
  if (isTRUE(run_extraction)) {
    link_topic_requested_cores <- if (is.null(local_threads)) {
      .available_cores(logical = TRUE)
    } else {
      max(1L, as.integer(local_threads[[1L]]))
    }
    link_topic_core_cap <- getOption("craftgrn.topic_link.max_cores", NULL)
    if (is.null(link_topic_core_cap)) {
      env_cap <- Sys.getenv("CRAFTGRN_TOPIC_LINK_MAX_CORES", unset = "")
      link_topic_core_cap <- if (nzchar(env_cap)) suppressWarnings(as.integer(env_cap)) else 8L
    }
    link_topic_core_cap <- suppressWarnings(as.integer(link_topic_core_cap[[1L]]))
    if (!is.finite(link_topic_core_cap) || link_topic_core_cap < 1L) link_topic_core_cap <- 8L
    link_topic_n_cores <- min(link_topic_requested_cores, link_topic_core_cap)
    extraction_plan <- method_plan[method %in% extraction_methods]
    for (i in seq_len(nrow(extraction_plan))) {
      row <- extraction_plan[i]
      model_root <- row$topic_models_dir[[1L]]
      extract_roots <- .m3tb_extraction_output_dirs(row, k_grid, output_layout)
      assignment_method <- extraction_topic_report_args$topic_term_assignment_method %||%
        .module3_default_term_assignment_method(row$method[[1L]])
      sequential_extraction <- .module3_extraction_requires_sequential(
        extraction_topic_report_args,
        assignment_method
      )
      k_plan <- .module3_extraction_k_worker_plan(
        n_tasks = length(k_grid),
        model_dir = model_root,
        k_workers = extraction_k_workers,
        k_max_workers = extraction_k_max_workers,
        k_memory_gb = extraction_k_memory_gb,
        k_memory_reserve_gb = extraction_k_memory_reserve_gb,
        cores = link_topic_requested_cores,
        link_scoring = sequential_extraction
      )
      # Forking after the native OpenMP WarpLDA sampler has run can leave child
      # processes with an invalid OpenMP runtime state. Keep same-process
      # training plus extraction sequential; extraction-only runs can still
      # use the adaptive K worker plan.
      if (isTRUE(run_training) && identical(row$backend[[1L]], "warplda") && k_plan$workers > 1L) {
        k_plan$workers <- 1L
        k_plan$cores <- min(k_plan$cores, link_topic_requested_cores)
        k_plan$reason <- "sequential after same-process OpenMP WarpLDA training"
      }
      link_cores_per_worker <- max(1L, floor(link_topic_requested_cores / k_plan$workers))
      link_cores_per_worker <- min(link_cores_per_worker, link_topic_n_cores)
      if (isTRUE(verbose)) {
        .log_inform(
          "Extracting Module 3 topics for method: {row$method_setup}; K={paste(k_grid, collapse = ',')}; {k_plan$workers} concurrent K worker(s); estimated {(.format_bytes(k_plan$per_worker_bytes))} per worker; reserve {(.format_bytes(k_plan$reserve_bytes))}; {k_plan$reason}."
        )
      }
      tasks <- lapply(seq_along(k_grid), function(j) list(
        k = k_grid[[j]],
        output_dir = extract_roots[[j]]
      ))
      row_report_args <- modifyList(
        list(
          fp_term_mode = row$fp_mode[[1L]],
          link_topic_n_cores = link_cores_per_worker,
          topic_model_family = .module3_topic_model_family(row$method[[1L]]),
          topic_term_assignment_method = assignment_method
        ),
        extraction_topic_report_args
      )
      row_report_args$thrP <- .module3_resolve_gammafit_thrP(
        extraction_topic_report_args$thrP %||% NULL,
        method = row$method[[1L]]
      )
      worker_fun <- function(task) {
        extract_regulatory_topics(
          k = task$k,
          model_dir = model_root,
          output_dir = task$output_dir,
          backend = row$backend[[1L]],
          vae_variant = row$vae_variant[[1L]],
          doc_mode = "tf",
          weight_label = row$weight_label[[1L]],
          flatten_single_output = .m3tb_flat_output_for_layout(output_layout),
          topic_report_args = row_report_args
        )
      }
      .module3_run_extraction_k_tasks(
        tasks = tasks,
        worker_fun = worker_fun,
        workers = k_plan$workers,
        cores = k_plan$cores
      )
    }
  }
  model_rows <- data.table::data.table()
  if (isTRUE(run_reports)) {
    model_rows <- .m3tb_existing_model_rows(output_dir, method_plan, k_grid)
  } else if (isTRUE(run_training) || isTRUE(run_extraction)) {
    model_rows <- tryCatch(
      .m3tb_existing_model_rows(output_dir, method_plan, k_grid),
      error = function(e) data.table::data.table()
    )
  }
  score_result <- NULL
  link_summary <- NULL
  html <- character()
  fibroblast_paper_recovery <- NULL
  if (isTRUE(run_reports)) {
    if (isTRUE(verbose)) .log_inform("Scoring Module 3 theta separation and writing review reports.")
    report_model_rows <- .m3tb_topic_space_rows(model_rows)
    score_prefix <- if (isTRUE(replicate_documents)) {
      "theta_condition_replicate_separation"
    } else {
      "theta_condition_separation"
    }
    score_result <- .m3tb_score_theta_outputs(
      output_dir,
      report_model_rows,
      comparisons,
      score_prefix = score_prefix,
      replicate_documents = replicate_documents,
      multiomic_data = multiomic_data,
      review_dir = review_dir
    )
    attr(score_result, "model_rows") <- report_model_rows
    attr(score_result, "report_state") <- report_state
    link_summary <- .m3tb_summarize_topic_links(output_dir, model_rows, review_dir = review_dir)
    html <- .m3tb_write_review_outputs(output_dir, score_result, link_summary, review_dir = review_dir)
    if (isTRUE(run_fibroblast_paper_recovery) && !is.null(paper_benchmark_csv)) {
      fibroblast_paper_recovery <- .m3fb_score_existing_models(
        output_dir = output_dir,
        model_rows = model_rows,
        benchmark_csv = paper_benchmark_csv,
        review_dir = review_dir,
        verbose = verbose
      )
    }
  }
  path_length_summary <- NULL
  if (isTRUE(path_length_check)) {
    path_length_summary <- .m3tb_check_path_lengths(
      root = output_dir,
      max_chars = path_length_limit,
      path_alias = path_length_alias,
      out_file = file.path(review_dir, "tables", "path_length_sanity_check.csv"),
      verbose = verbose
    )
  }
  invisible(list(
    method_plan = method_plan,
    model_rows = model_rows,
    score = score_result,
    topic_link_summary = link_summary,
    fibroblast_paper_recovery = fibroblast_paper_recovery,
    html = html,
    output_layout = output_layout,
    review_dir = review_dir,
    path_length_summary = path_length_summary,
    replicate_documents = isTRUE(replicate_documents),
    multiomic_data_provided = !is.null(multiomic_data)
  ))
}

#' Prepare Module 3 topic-model inputs
#'
#' Builds and caches the document-level link table, document-term table, sparse
#' document-term matrix, and summary metadata used by Module 3 topic modeling.
#' This step is useful when you want to inspect topic inputs before training or
#' reuse the same input cache for extraction.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   CSV files.
#' @param output_dir Directory where topic input caches are written.
#' @param tf_cluster_map Optional named vector mapping TF names to motif
#'   clusters. Required only when `doc_mode = "tf_cluster"`.
#' @param doc_mode Document mode, either `"tf"` or `"tf_cluster"`.
#' @param doc_design Document design, either `"condition"` or `"comparison"`.
#' @param fp_term_mode Footprint term mode: `"aggregate_weight"`,
#'   `"aggregate"`, or `"unique"`.
#' @param gene_term_mode Gene term mode passed to comparison document-term
#'   construction.
#' @param sample_subset Optional condition/sample labels to retain.
#' @param analysis_label Label written to the input summary.
#' @param count_method Count conversion method.
#' @param count_scale Count scaling factor.
#' @param count_input Count column used for the cached sparse topic matrix. If
#'   `NULL`, CraftGRN uses `pseudo_count_log` for `count_method = "log"` and
#'   `pseudo_count_bin` otherwise.
#' @param input_source Input source type. Use `"differential_links"` for
#'   filtered comparison links and `"condition_links"` for condition-native
#'   links from `module3_prepare_condition_links()`.
#' @param threshold_gene_expr Minimum condition-level target-gene expression.
#' @param threshold_fp_score Minimum condition-level footprint score.
#' @param threshold_tf_expr Minimum condition-level TF expression.
#' @param include_tf_terms Whether to include TF self-terms.
#' @param top_terms_per_doc Optional maximum terms per document.
#' @param min_df Minimum document frequency for terms.
#' @param abs_log2fc_fp_min Minimum absolute footprint fold-change filter.
#' @param abs_delta_fp_min Minimum absolute footprint delta filter.
#' @param abs_log2fc_gene_min Minimum absolute target-gene fold-change filter.
#' @param require_fp_bound_either Require footprint binding in either condition.
#' @param require_tf_expr_either Require TF expression in either condition.
#' @param require_gene_expr_either Require target-gene expression in either
#'   condition.
#' @param direction_consistency Direction consistency filter mode.
#' @param save_full_doc_term_csv Whether to write the full document-term CSV.
#' @param check_repeated_values Warn about repeated inconsistent term values.
#'   The high-throughput default is `FALSE`; set to `TRUE` for diagnostic
#'   audits.
#' @param overwrite If FALSE, reuse an existing complete cache.
#' @param verbose Emit concise progress messages.
#'
#' @return A list with cache paths and input summary counts.
#' @noRd
module3_prepare_topic_inputs <- function(filtered_dir,
                                         output_dir,
                                         tf_cluster_map,
                                         input_source = c("differential_links", "condition_links"),
                                         doc_mode = c("tf", "tf_cluster"),
                                         doc_design = c("condition", "comparison"),
                                         fp_term_mode = c("aggregate_weight", "aggregate", "unique"),
                                         gene_term_mode = c("unique", "aggregate"),
                                         sample_subset = NULL,
                                         analysis_label = NULL,
                                         count_method = c("bin", "log"),
                                         count_scale = 50,
                                         count_input = NULL,
                                         threshold_gene_expr = 0,
                                         threshold_fp_score = 0,
                                         threshold_tf_expr = -Inf,
                                         include_tf_terms = TRUE,
                                         top_terms_per_doc = Inf,
                                         min_df = 2,
                                         abs_log2fc_fp_min = 0,
                                         abs_delta_fp_min = 1,
                                         abs_log2fc_gene_min = 1,
                                         require_fp_bound_either = TRUE,
                                         require_tf_expr_either = TRUE,
                                         require_gene_expr_either = TRUE,
                                         direction_consistency = "aligned",
                                         save_full_doc_term_csv = FALSE,
                                         check_repeated_values = FALSE,
                                         overwrite = FALSE,
                                         verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("Matrix")
  input_source <- match.arg(input_source)
  doc_mode <- match.arg(doc_mode)
  doc_design <- match.arg(doc_design)
  if (identical(input_source, "condition_links") && !identical(doc_design, "condition")) {
    .log_abort("input_source = 'condition_links' requires doc_design = 'condition'.")
  }
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  gene_term_mode <- match.arg(gene_term_mode)
  count_method <- match.arg(count_method)
  count_input_requested <- if (is.null(count_input) || !length(count_input)) NA_character_ else as.character(count_input[[1L]])
  count_input_effective <- .resolve_topic_count_input(count_method = count_method, count_input = count_input)
  sample_subset <- if (is.null(sample_subset)) NULL else unique(as.character(sample_subset))
  sample_subset <- sample_subset[!is.na(sample_subset) & nzchar(sample_subset)]
  input_signature <- .module3_topic_input_signature(
    input_dir = filtered_dir,
    input_source = input_source,
    sample_subset = sample_subset,
    settings = list(
      doc_mode = doc_mode,
      doc_design = doc_design,
      fp_term_mode = fp_term_mode,
      gene_term_mode = gene_term_mode,
      count_method = count_method,
      count_scale = as.numeric(count_scale),
      count_input_effective = count_input_effective,
      threshold_gene_expr = threshold_gene_expr,
      threshold_fp_score = threshold_fp_score,
      threshold_tf_expr = threshold_tf_expr,
      include_tf_terms = isTRUE(include_tf_terms),
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      abs_log2fc_fp_min = abs_log2fc_fp_min,
      abs_delta_fp_min = abs_delta_fp_min,
      abs_log2fc_gene_min = abs_log2fc_gene_min,
      require_fp_bound_either = isTRUE(require_fp_bound_either),
      require_tf_expr_either = isTRUE(require_tf_expr_either),
      require_gene_expr_either = isTRUE(require_gene_expr_either),
      direction_consistency = direction_consistency
    )
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rds_dir <- file.path(output_dir, "rds")
  required_cache <- file.path(rds_dir, c("edges_docs.rds", "doc_term.rds", "dtm.rds", "dtm_index.rds"))
  summary_path <- file.path(output_dir, "topic_input_summary.csv")
  if (!isTRUE(overwrite) && all(file.exists(required_cache)) && file.exists(summary_path)) {
    summary_dt <- data.table::fread(summary_path, showProgress = FALSE)
    cache_matches <- nrow(summary_dt) &&
      all(c("input_source", "doc_design", "doc_mode", "fp_term_mode", "input_signature") %in% names(summary_dt)) &&
      identical(as.character(summary_dt$input_source[[1L]]), input_source) &&
      identical(as.character(summary_dt$doc_design[[1L]]), doc_design) &&
      identical(as.character(summary_dt$doc_mode[[1L]]), doc_mode) &&
      identical(as.character(summary_dt$fp_term_mode[[1L]]), fp_term_mode) &&
      identical(as.character(summary_dt$input_signature[[1L]]), input_signature) &&
      all(c("count_method", "count_input_effective") %in% names(summary_dt)) &&
      identical(as.character(summary_dt$count_method[[1L]]), count_method) &&
      identical(as.character(summary_dt$count_input_effective[[1L]]), count_input_effective)
    if (isTRUE(cache_matches)) {
      if (isTRUE(verbose)) .log_inform("Reusing existing Module 3 topic input cache: {output_dir}")
      return(invisible(list(output_dir = output_dir, summary = summary_dt, reused = TRUE)))
    }
    if (isTRUE(verbose)) {
      .log_inform("Existing Module 3 topic input cache does not match requested document settings; rebuilding: {output_dir}")
    }
  }
  if (identical(input_source, "condition_links")) {
    if (isTRUE(verbose)) .log_inform("Loading condition-native links from {filtered_dir}.")
    edges_dt <- .module3_read_condition_links(filtered_dir, conditions = sample_subset)
    n_loaded <- nrow(edges_dt)
  } else {
    delta_files <- .module3_filtered_link_files(filtered_dir)
    if (!length(delta_files)) {
      delta_files <- list.files(filtered_dir, "_filtered_links(_(up|down))?\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) {
      delta_files <- list.files(filtered_dir, "_delta_links_filtered(_(up|down))?\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) {
      delta_files <- list.files(filtered_dir, "_delta_links\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) .log_abort("No Module 3 filtered link files found in {filtered_dir}.")
    if (isTRUE(verbose)) .log_inform("Loading {length(delta_files)} Module 3 filtered-link file(s).")
    edges_dt <- data.table::as.data.table(load_delta_links_many(delta_files, keep_original = FALSE))
    edges_dt <- .apply_module3_manifest_labels(edges_dt, filtered_dir)
    n_loaded <- nrow(edges_dt)
    if (!("comparison_id" %in% names(edges_dt))) .log_abort("Module 3 links are missing comparison_id.")
    if (length(sample_subset)) {
      if (!all(c("cond1_id", "cond2_id") %in% names(edges_dt))) {
        .log_abort("sample_subset requires cond1_id and cond2_id columns in Module 3 links.")
      }
      edges_dt <- edges_dt[cond1_id %in% sample_subset & cond2_id %in% sample_subset]
    }
  }
  if (!nrow(edges_dt)) .log_abort("No Module 3 links remain after sample subsetting.")
  edges_filt <- if (identical(input_source, "condition_links")) {
    edges_dt
  } else {
    filter_edges_for_tf_topics(
      edges_dt,
      abs_log2fc_fp_min = abs_log2fc_fp_min,
      abs_delta_fp_min = abs_delta_fp_min,
      abs_log2fc_gene_min = abs_log2fc_gene_min,
      require_fp_bound_either = require_fp_bound_either,
      require_tf_expr_either = require_tf_expr_either,
      require_gene_expr_either = require_gene_expr_either,
      direction_consistency = direction_consistency
    )
  }
  if (!nrow(edges_filt)) .log_abort("No Module 3 links passed topic-input filters.")
  if (identical(doc_design, "condition")) {
    edges_docs <- if (identical(input_source, "condition_links")) {
      data.table::copy(edges_filt)
    } else {
      add_condition_tf_docs(edges_filt, tf_cluster_map = tf_cluster_map, doc_mode = doc_mode)
    }
    doc_term <- build_doc_term_condition_union(
      edges_docs,
      count_method = count_method,
      count_scale = count_scale,
      prefix_terms = TRUE,
      threshold_gene_expr = threshold_gene_expr,
      threshold_fp_score = threshold_fp_score,
      threshold_tf_expr = threshold_tf_expr,
      include_tf_terms = isTRUE(include_tf_terms),
      require_tf_expr = identical(doc_mode, "tf"),
      fp_term_mode = fp_term_mode,
      check_repeated_values = check_repeated_values
    )
  } else {
    edges_docs <- add_tf_docs(edges_filt, doc_mode = doc_mode, direction_by = "gene", tf_cluster_map = tf_cluster_map)
    doc_term <- build_doc_term_joint(
      edges_docs,
      weight_type_peak = "log2fc_fp",
      weight_type_gene = "log2fc_gene",
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      distinct_terms = TRUE,
      gene_term_mode = gene_term_mode,
      fp_term_mode = fp_term_mode,
      include_tf_terms = isTRUE(include_tf_terms),
      tf_weight_type = "log2fc_tf",
      balance_mode = "min",
      prefix_terms = TRUE,
      threshold_gene_expr = threshold_gene_expr,
      threshold_fp_score = threshold_fp_score,
      threshold_tf_expr = threshold_tf_expr,
      require_condition_thresholds = identical(doc_mode, "tf"),
      check_repeated_values = check_repeated_values
    )
  }
  if (!nrow(doc_term)) .log_abort("Module 3 document-term table is empty.")
  if (identical(doc_design, "condition")) {
    .write_module3_document_term_qc(
      doc_term = doc_term,
      output_dir = output_dir,
      count_column = count_input_effective,
      title = paste0(
        "Module 3 document-term QC - ",
        as.character(analysis_label %||% "condition documents")[[1L]]
      ),
      verbose = verbose
    )
  }
  write_doc_term_cache(doc_term, out_dir = output_dir, save_full_doc_term_csv = isTRUE(save_full_doc_term_csv))
  .save_all(output_dir, "edges_filtered", edges_filt)
  .save_all(output_dir, "edges_docs", edges_docs)
  .save_all(output_dir, "doc_term", doc_term)
  dtm_obj <- build_sparse_dtm(doc_term, count_col = count_input_effective)
  .save_all(output_dir, "dtm", dtm_obj$dtm)
  .save_all(output_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))
  summary_dt <- data.table::data.table(
    analysis_label = if (is.null(analysis_label)) NA_character_ else as.character(analysis_label[[1L]]),
    input_signature = input_signature,
    input_source = input_source,
    doc_design = doc_design,
    doc_mode = doc_mode,
    fp_term_mode = fp_term_mode,
    count_method = count_method,
    count_scale = as.numeric(count_scale),
    count_input_requested = count_input_requested,
    count_input_effective = count_input_effective,
    n_link_rows_loaded = as.double(n_loaded),
    n_link_rows_after_subset = as.double(nrow(edges_dt)),
    n_link_rows_after_filter = as.double(nrow(edges_filt)),
    n_document_edge_rows = as.double(nrow(edges_docs)),
    n_doc_term_rows = as.double(nrow(doc_term)),
    n_documents = as.double(data.table::uniqueN(doc_term$doc_id)),
    n_terms = as.double(data.table::uniqueN(doc_term$term_id)),
    n_nonzero = as.double(Matrix::nnzero(dtm_obj$dtm)),
    n_model_tokens = as.double(sum(.safe_num(doc_term[[count_input_effective]]), na.rm = TRUE))
  )
  data.table::fwrite(summary_dt, summary_path)
  if (isTRUE(verbose)) {
    .log_inform("Prepared Module 3 topic inputs: {summary_dt$n_documents} document(s), {summary_dt$n_terms} term(s), {summary_dt$n_doc_term_rows} doc-term row(s).")
  }
  invisible(list(output_dir = output_dir, summary = summary_dt, reused = FALSE))
}

#' Construct input documents for topic modeling
#'
#' @description
#' Builds and caches the document-level link table, document-term table, sparse
#' document-term matrix, and summary metadata used by Module 3 topic modeling.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   CSV files.
#' @param output_dir Directory where topic input caches are written.
#' @param tf_cluster_map Named vector mapping TF names to motif clusters.
#' @param check_repeated_values Warn about repeated inconsistent term values.
#'   The high-throughput default is `FALSE`; set to `TRUE` for diagnostic
#'   audits.
#' @param ... Additional topic-document construction arguments passed to the
#'   internal Module 3 document builder.
#'
#' @return A list with cache paths and input summary counts.
#' @export
module3_construct_docs <- function(filtered_dir,
                                   output_dir,
                                   tf_cluster_map = NULL,
                                   check_repeated_values = FALSE,
                                   ...) {
  module3_prepare_topic_inputs(
    filtered_dir = filtered_dir,
    output_dir = output_dir,
    tf_cluster_map = tf_cluster_map,
    check_repeated_values = check_repeated_values,
    ...
  )
}

#' Build a Module 3 QC HTML report
#'
#' Writes a self-contained HTML report for Module 3 topic-model outputs. The
#' report summarizes topic-input caches, model rows, theta separation scores,
#' compact topic-link pass counts, and differential-link summaries when
#' available.
#'
#' @param topic_dir Module 3 topic output directory.
#' @param output_dir Directory where the report is written. Defaults to
#'   `topic_dir/reports`.
#' @param differential_links_dir Optional Module 3 differential-link directory.
#'   If `NULL`, CraftGRN tries to detect a sibling or nested
#'   `differential_links` directory.
#' @param title Report title.
#' @param top_n Number of top differential TFs retained per comparison in the
#'   QC summary CSV.
#' @param verbose Emit concise progress messages.
#'
#' @return Path to the HTML report.
#' @export
build_module3_qc_report <- function(topic_dir,
                                    output_dir = file.path(topic_dir, "reports"),
                                    differential_links_dir = NULL,
                                    title = "Module 3 QC report",
                                    top_n = 20L,
                                    verbose = TRUE) {
  .assert_pkg("data.table")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(output_dir, "module3_qc_report.html")
  review_csv <- .m3tb_review_read_dir(topic_dir)
  read_optional <- function(path) {
    if (file.exists(path)) data.table::fread(path, showProgress = FALSE) else data.table::data.table()
  }
  input_summary <- read_optional(file.path(topic_dir, "topic_documents", "topic_input_summary.csv"))
  if (!nrow(input_summary)) {
    candidates <- list.files(topic_dir, "topic_input_summary[.]csv$", recursive = TRUE, full.names = TRUE)
    if (length(candidates)) input_summary <- read_optional(candidates[[1L]])
  }
  method_plan <- read_optional(file.path(review_csv, "module3_topic_method_plan.csv"))
  theta_scores <- read_optional(file.path(review_csv, "theta_condition_separation_scores.csv"))
  if (!nrow(theta_scores)) theta_scores <- read_optional(file.path(review_csv, "theta_condition_replicate_separation_scores.csv"))
  pass_counts <- read_optional(file.path(review_csv, "topic_setup_pass_state_counts.csv"))
  link_summaries <- list.files(topic_dir, "topic_link_summary[.]csv$", recursive = TRUE, full.names = TRUE)
  link_summary <- if (length(link_summaries)) {
    data.table::rbindlist(lapply(link_summaries, data.table::fread, showProgress = FALSE), use.names = TRUE, fill = TRUE)
  } else {
    data.table::data.table()
  }
  differential_links_dir <- .module3_qc_detect_differential_links_dir(topic_dir, differential_links_dir)
  differential_summary <- .module3_qc_read_differential_summary(differential_links_dir)
  differential_tfs <- .module3_qc_summarize_differential_tfs(
    differential_links_dir = differential_links_dir,
    output_dir = output_dir,
    top_n = top_n,
    verbose = verbose
  )
  table_html <- function(dt, empty_label = "No rows available") {
    if (!nrow(dt)) return(paste0("<p class=\"muted\">", .m3tb_html_escape(empty_label), "</p>"))
    dt <- utils::head(dt, 20L)
    cols <- names(dt)
    rows <- apply(as.data.frame(dt), 1, function(x) {
      paste0("<tr>", paste0("<td>", .m3tb_html_escape(x), "</td>", collapse = ""), "</tr>")
    })
    paste0(
      "<table><thead><tr>",
      paste0("<th>", .m3tb_html_escape(cols), "</th>", collapse = ""),
      "</tr></thead><tbody>",
      paste(rows, collapse = ""),
      "</tbody></table>"
    )
  }
  metric_card <- function(label, value) {
    paste0("<div class=\"card\"><div class=\"label\">", .m3tb_html_escape(label), "</div><div class=\"value\">", .m3tb_html_escape(value), "</div></div>")
  }
  n_models <- if (nrow(method_plan)) nrow(method_plan) else length(list.files(topic_dir, "theta_K[0-9]+[.]csv$", recursive = TRUE))
  n_pass <- if (nrow(pass_counts) && "status" %in% names(pass_counts)) sum(pass_counts$status == "Pass" & is.finite(pass_counts$count), na.rm = TRUE) else NA_real_
  n_diff_comparisons <- if (nrow(differential_summary) && "comparison_id" %in% names(differential_summary)) data.table::uniqueN(differential_summary$comparison_id) else NA_real_
  n_diff_links <- if (nrow(differential_summary) && all(c("n_up", "n_down") %in% names(differential_summary))) {
    sum(suppressWarnings(as.numeric(differential_summary$n_up)), suppressWarnings(as.numeric(differential_summary$n_down)), na.rm = TRUE)
  } else {
    NA_real_
  }
  css <- paste0(
    "body{font-family:Arial,sans-serif;margin:30px;color:#17212b;background:#fbfbf8}",
    "h1{font-size:26px;margin-bottom:4px}h2{font-size:18px;margin-top:28px}",
    ".muted{color:#667085}.cards{display:grid;grid-template-columns:repeat(4,minmax(130px,1fr));gap:12px;margin:20px 0}",
    ".card{border:1px solid #d7d1c3;background:#fff;padding:12px;border-radius:4px}",
    ".label{font-size:12px;color:#667085}.value{font-size:22px;font-weight:700;color:#194b5f}",
    "table{border-collapse:collapse;width:100%;background:#fff}th,td{border:1px solid #d7d1c3;padding:6px 8px;font-size:12px}th{background:#eef3f1;text-align:left}"
  )
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\">",
    paste0("<title>", .m3tb_html_escape(title), "</title><style>", css, "</style></head><body>"),
    paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
    paste0("<p class=\"muted\">Topic directory: ", .m3tb_html_escape(normalizePath(topic_dir, winslash = "/", mustWork = FALSE)), "</p>"),
    "<div class=\"cards\">",
    metric_card("Model rows", n_models),
    metric_card("Input docs", if (nrow(input_summary)) input_summary$n_documents[[1L]] else "NA"),
    metric_card("Input terms", if (nrow(input_summary)) input_summary$n_terms[[1L]] else "NA"),
    metric_card("Passing links", if (is.finite(n_pass)) n_pass else "NA"),
    metric_card("Diff comparisons", if (is.finite(n_diff_comparisons)) n_diff_comparisons else "NA"),
    metric_card("Diff links", if (is.finite(n_diff_links)) n_diff_links else "NA"),
    "</div>",
    "<h2>Topic Input Summary</h2>", table_html(input_summary),
    "<h2>Differential GRN Summary</h2>", table_html(differential_summary),
    "<h2>Top Differential TFs</h2>", table_html(differential_tfs$combined),
    "<h2>Method Plan</h2>", table_html(method_plan),
    "<h2>Theta Separation Scores</h2>", table_html(theta_scores),
    "<h2>Topic-Link Pass Counts</h2>", table_html(pass_counts),
    "<h2>Topic-Link Output Summary</h2>", table_html(link_summary),
    "</body></html>"
  )
  writeLines(html, report_path, useBytes = TRUE)
  if (isTRUE(verbose)) .log_inform("Wrote Module 3 QC report: {report_path}")
  invisible(report_path)
}

#' Export interactive HTML browsers of topic modeling results
#'
#' @description
#' Builds a self-contained index browser for existing Module 3 topic-modeling
#' review outputs at the topic, condition, comparison, and pathway levels. This
#' function organizes existing outputs and does not train or extract models.
#'
#' @param topic_dir Module 3 topic output directory.
#' @param output_dir Directory where the browser HTML and manifest are written.
#' @param include Existing output families to include.
#' @param verbose Emit concise progress messages.
#'
#' @return Path to the HTML browser.
#' @export
visualize_topic_modeling_results <- function(topic_dir,
                                             output_dir = file.path(topic_dir, "reports"),
                                             include = c("topic", "condition", "comparison", "pathway"),
                                             verbose = TRUE) {
  if (!dir.exists(topic_dir)) .log_abort("`topic_dir` not found: {topic_dir}")
  include <- match.arg(include, choices = c("topic", "condition", "comparison", "pathway"), several.ok = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  html_files <- list.files(topic_dir, "[.]html$", recursive = TRUE, full.names = TRUE)
  html_files <- html_files[normalizePath(html_files, winslash = "/", mustWork = FALSE) != normalizePath(file.path(output_dir, "topic_modeling_results.html"), winslash = "/", mustWork = FALSE)]
  classify <- function(path) {
    b <- tolower(basename(path))
    p <- tolower(path)
    if (grepl("pathway", b) || grepl("pathway", p)) return("pathway")
    if (grepl("condition|group_mds", b)) return("condition")
    if (grepl("comparison|theta_phi", b)) return("comparison")
    "topic"
  }
  manifest <- data.table::data.table(
    family = vapply(html_files, classify, character(1L)),
    file = basename(html_files),
    path = html_files,
    rel_path = .m3tb_relative_path(html_files, output_dir)
  )
  manifest <- manifest[family %in% include]
  data.table::fwrite(manifest, file.path(output_dir, "topic_modeling_results_manifest.csv"))
  table_rows <- if (nrow(manifest)) {
    paste0(
      "<tr><td>", .m3tb_html_escape(manifest$family), "</td><td>",
      "<a href=\"", .m3tb_html_escape(manifest$rel_path), "\">",
      .m3tb_html_escape(manifest$file), "</a></td></tr>",
      collapse = ""
    )
  } else {
    "<tr><td colspan=\"2\">No existing topic-modeling HTML outputs found.</td></tr>"
  }
  html <- c(
    "<!doctype html><html><head><meta charset=\"utf-8\"><title>Topic modeling results</title>",
    "<style>body{font-family:Arial,sans-serif;margin:28px;color:#17212b;background:#fbfbf8}table{border-collapse:collapse;width:100%;background:#fff}th,td{border:1px solid #d7d1c3;padding:7px 9px;font-size:13px}th{background:#eef3f1;text-align:left}input{padding:7px;width:280px;margin:8px 0 14px 0}</style>",
    "</head><body><h1>Topic modeling results</h1>",
    "<input id=\"q\" placeholder=\"Filter reports\" oninput=\"filterRows()\">",
    "<table id=\"tbl\"><thead><tr><th>Family</th><th>Report</th></tr></thead><tbody>",
    table_rows,
    "</tbody></table>",
    "<script>function filterRows(){const q=document.getElementById('q').value.toLowerCase();document.querySelectorAll('#tbl tbody tr').forEach(r=>{r.style.display=r.textContent.toLowerCase().includes(q)?'':'none';});}</script>",
    "</body></html>"
  )
  out <- file.path(output_dir, "topic_modeling_results.html")
  writeLines(html, out, useBytes = TRUE)
  if (isTRUE(verbose)) .log_inform("Wrote topic-modeling results browser: {out}")
  invisible(out)
}

#' Export an interactive HTML browser of differential GRNs
#'
#' @description
#' Builds an interactive TF-to-gene network browser from Module 3 filtered
#' differential links. Users can select a comparison, choose up or down
#' differential links, adjust the number of top TFs and links to display, and
#' inspect footprint-supported edge evidence in tooltips.
#'
#' @param differential_links_dir Module 3 differential-link directory.
#' @param output_dir Directory where the browser HTML and CSV summaries are
#'   written.
#' @param top_tf_n Default number of top TFs shown in the browser.
#' @param top_link_n Default number of top TF-to-gene links shown in the
#'   browser.
#' @param default_direction Initial direction selected in the browser.
#' @param browser_max_rows_per_file Maximum filtered-link rows read from each
#'   comparison/direction file when building the browser payload. The full
#'   filtered-link CSVs remain the authoritative data source; this cap keeps
#'   the self-contained HTML browser responsive for large projects.
#' @param top_n Deprecated compatibility alias for \code{top_tf_n}.
#' @param verbose Emit concise progress messages.
#'
#' @return Path to the HTML browser.
#' @export
visualize_differential_grns <- function(differential_links_dir,
                                        output_dir = file.path(differential_links_dir, "reports"),
                                        top_tf_n = 10L,
                                        top_link_n = 300L,
                                        default_direction = "up",
                                        browser_max_rows_per_file = 50000L,
                                        top_n = NULL,
                                        verbose = TRUE) {
  if (!dir.exists(differential_links_dir)) .log_abort("`differential_links_dir` not found: {differential_links_dir}")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(top_n)) top_tf_n <- top_n
  top_tf_n <- max(1L, as.integer(top_tf_n)[[1L]])
  top_link_n <- max(1L, as.integer(top_link_n)[[1L]])
  browser_max_rows_per_file <- max(1L, as.integer(browser_max_rows_per_file)[[1L]])
  default_direction <- match.arg(default_direction, c("up", "down"))
  browser_tf_cap <- max(20L, top_tf_n * 2L)
  diff_summary <- .module3_qc_read_differential_summary(differential_links_dir)
  tf_summary <- .module3_qc_summarize_differential_tfs(
    differential_links_dir = differential_links_dir,
    output_dir = output_dir,
    top_n = browser_tf_cap,
    verbose = verbose
  )
  summary_path <- file.path(output_dir, "differential_grn_summary.csv")
  data.table::fwrite(diff_summary, summary_path)
  tf_keep <- NULL
  if (is.data.frame(tf_summary$combined) && nrow(tf_summary$combined) && all(c("comparison_id", "TF") %in% names(tf_summary$combined))) {
    tf_keep <- split(toupper(as.character(tf_summary$combined$TF)), tf_summary$combined$comparison_id)
    tf_keep <- lapply(tf_keep, unique)
  }
  network_dir <- file.path(output_dir, "differential_grn_networks")
  dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)
  browser_edge_cap <- max(300L, top_link_n)
  network <- .module3_diff_grn_aggregate_network_edges(
    differential_links_dir = differential_links_dir,
    browser_edge_cap_per_group = browser_edge_cap,
    tf_keep_by_comparison = tf_keep,
    browser_max_rows_per_file = browser_max_rows_per_file
  )
  if (nrow(network$manifest)) network$manifest[, browser_max_rows_per_file := browser_max_rows_per_file]
  edges_path <- file.path(network_dir, "differential_grn_edges.csv")
  browser_edges_path <- file.path(network_dir, "differential_grn_browser_edges.csv")
  manifest_path <- file.path(output_dir, "differential_grn_browser_manifest.csv")
  data.table::fwrite(network$edges, edges_path)
  data.table::fwrite(network$browser_edges, browser_edges_path)
  data.table::fwrite(network$manifest, manifest_path)
  out <- file.path(output_dir, "differential_grns.html")
  .module3_diff_grn_write_network_browser(
    edges = network$browser_edges,
    manifest = network$manifest,
    out_html = out,
    top_tf_n = top_tf_n,
    top_link_n = top_link_n,
    default_direction = default_direction,
    title = sprintf("Differential GRN browser (%s default)", default_direction)
  )
  if (isTRUE(verbose)) .log_inform("Wrote differential GRN browser: {out}")
  invisible(out)
}

#' Run regulatory topic modeling
#'
#' Production-oriented Module 3 wrapper for one selected topic-document method.
#' It uses the clean standard output layout, compact topic-link output by
#' default, and writes a Module 3 QC report.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   files, or condition-link files when `input_source = "condition_links"`.
#' @param multiomic_data Optional CraftGRN multiomic object. Required when
#'   `replicate_documents = TRUE`.
#' @param comparisons Comparison or condition grouping table, or a CSV path.
#' @param output_dir Topic output directory.
#' @param method Single Module 3 method ID.
#' @param k_grid Integer topic numbers. Multiple K values are allowed for
#'   standard K review within the selected method.
#' @param replicate_documents Whether theta document labels are replicate
#'   resolved for condition-separation scoring.
#' @param reuse_if_exists Reuse existing model outputs where possible.
#' @param local_threads Optional thread count for model training.
#' @param warplda_iterations Number of native WarpLDA iterations.
#' @param input_source Input source type. Use `"differential_links"` for
#'   comparison links and `"condition_links"` for condition-native links.
#' @param sample_subset Optional condition/sample labels to keep.
#' @param analysis_label Label used to name the topic-model analysis.
#' @param topic_link_output Topic-link output mode. `"pass"` writes compact
#'   passing links and summaries only; `"full"` writes exhaustive all-topic
#'   links; `"both"` writes both; `"none"` disables topic-link export.
#' @param extraction_topic_report_args Optional named list of topic-extraction
#'   report argument overrides.
#' @param run_training Train topic models before reporting.
#' @param run_extraction Run topic extraction before reporting.
#' @param run_reports Build score tables and review reports.
#' @param build_qc_report Whether to build the Module 3 QC report.
#' @param verbose Emit concise progress messages.
#'
#' @return An invisible list with topic input/model/extraction paths, review
#'   outputs, and `qc_report` when requested.
#' @noRd
run_regulatory_topics <- function(filtered_dir,
                                  multiomic_data = NULL,
                                  comparisons,
                                  output_dir,
                                  method = "comparison_aggr_multivi",
                                  k_grid = 10L,
                                  replicate_documents = FALSE,
                                  reuse_if_exists = TRUE,
                                  local_threads = NULL,
                                  warplda_iterations = 2000L,
                                  count_method = c("log", "bin"),
                                  count_input = NULL,
                                  vae_python = NULL,
                                  vae_epochs = 200L,
                                  vae_batch_size = 64L,
                                  vae_hidden = 128L,
                                  vae_lr = 1e-3,
                                  vae_seed = 123L,
                                  vae_device = "auto",
                                  input_source = c("differential_links", "condition_links"),
                                  sample_subset = NULL,
                                  analysis_label = NULL,
                                  topic_link_output = c("pass", "full", "both", "none"),
                                  pathway_backend = NULL,
                                  report_state = list(),
                                  extraction_topic_report_args = list(),
                                  extraction_k_workers = NULL,
                                  extraction_k_max_workers = 4L,
                                  memory_safety = c("strict", "adaptive", "off"),
                                  memory_max_fraction = 0.8,
                                  run_training = TRUE,
                                  run_extraction = TRUE,
                                  run_reports = TRUE,
                                  build_qc_report = TRUE,
                                  verbose = TRUE) {
  topic_link_output <- match.arg(topic_link_output)
  count_method <- match.arg(count_method)
  input_source <- match.arg(input_source)
  memory_safety <- match.arg(memory_safety)
  if (length(method) != 1L) .log_abort("run_regulatory_topics() expects one selected method.")
  extraction_args <- modifyList(
    list(
      link_topic_output = topic_link_output,
      pathway_link_scores_file = if (topic_link_output %in% c("pass", "both")) "topic_links_pass.csv" else "topic_links.csv",
      pathway_backend = .pathway_backend(pathway_backend)
    ),
    extraction_topic_report_args
  )
  res <- run_module3_topic_benchmark(
    filtered_dir = filtered_dir,
    multiomic_data = multiomic_data,
    comparisons = comparisons,
    output_dir = output_dir,
    methods = method,
    k_grid = k_grid,
    output_layout = "standard",
    replicate_documents = replicate_documents,
    reuse_if_exists = reuse_if_exists,
    local_threads = local_threads,
    warplda_iterations = warplda_iterations,
    count_method = count_method,
    count_input = count_input,
    vae_python = vae_python,
    vae_epochs = vae_epochs,
    vae_batch_size = vae_batch_size,
    vae_hidden = vae_hidden,
    vae_lr = vae_lr,
    vae_seed = vae_seed,
    vae_device = vae_device,
    input_source = input_source,
    sample_subset = sample_subset,
    analysis_label = analysis_label,
    report_state = report_state,
    extraction_topic_report_args = extraction_args,
    extraction_k_workers = extraction_k_workers,
    extraction_k_max_workers = extraction_k_max_workers,
    memory_safety = memory_safety,
    memory_max_fraction = memory_max_fraction,
    run_training = run_training,
    run_extraction = run_extraction,
    run_reports = run_reports,
    verbose = verbose
  )
  if (isTRUE(build_qc_report)) {
    res$qc_report <- build_module3_qc_report(output_dir, differential_links_dir = filtered_dir, verbose = verbose)
  }
  invisible(res)
}

.module3_read_project_config <- function(project_config) {
  if (is.null(project_config)) return(list())
  if (is.character(project_config) && length(project_config) == 1L) {
    if (!file.exists(project_config)) .log_abort("Module 3 project config not found: {project_config}")
    cfg <- yaml::read_yaml(project_config)
    .validate_project_config(cfg, source = basename(project_config))
    return(cfg)
  }
  if (is.list(project_config)) {
    .validate_project_config(project_config)
    return(project_config)
  }
  .log_abort("`project_config` must be NULL, a YAML path, or a list.")
}

.module3_cfg_value <- function(cfg, names, default = NULL) {
  for (nm in names) {
    if (!is.null(cfg[[nm]])) return(cfg[[nm]])
    if (startsWith(nm, "module3_") && is.list(cfg$module3)) {
      nested_name <- sub("^module3_", "", nm)
      if (!is.null(cfg$module3[[nested_name]])) return(cfg$module3[[nested_name]])
    }
  }
  default
}

.module3_cfg_int_vector <- function(x, default = NULL) {
  if (is.null(x)) return(default)
  if (is.character(x) && length(x) == 1L && grepl(",", x, fixed = TRUE)) {
    x <- strsplit(x, ",", fixed = TRUE)[[1L]]
  }
  vals <- suppressWarnings(as.integer(x))
  vals <- vals[is.finite(vals) & vals > 1L]
  vals <- sort(unique(vals))
  if (!length(vals)) return(default)
  vals
}

.module3_resolve_gammafit_thrP <- function(value = NULL, method = NULL) {
  family <- .module3_topic_model_family(method)
  if (is.null(value)) return(.module3_default_gammafit_thrP(family))
  if (is.numeric(value) && length(value) == 1L) {
    cutoff <- suppressWarnings(as.numeric(value[[1L]]))
  } else {
    values <- if (is.list(value)) unlist(value, use.names = TRUE) else value
    value_names <- names(values)
    if (is.null(value_names) || any(!nzchar(value_names))) {
      .log_abort(
        "topic_gammafit_thrP must be one probability or a named lda/multivi/vae_mlp/default mapping."
      )
    }
    value_names <- tolower(gsub("-", "_", value_names, fixed = TRUE))
    aliases <- c(warplda = "lda", vae = "vae_mlp")
    value_names[value_names %in% names(aliases)] <- aliases[value_names[value_names %in% names(aliases)]]
    if (anyDuplicated(value_names)) {
      .log_abort("topic_gammafit_thrP contains duplicate model-family entries after alias normalization.")
    }
    names(values) <- value_names
    selected <- if (family %in% names(values)) {
      values[[family]]
    } else if ("default" %in% names(values)) {
      values[["default"]]
    } else {
      NULL
    }
    if (is.null(selected)) return(.module3_default_gammafit_thrP(family))
    cutoff <- suppressWarnings(as.numeric(selected[[1L]]))
  }
  if (!is.finite(cutoff) || cutoff <= 0 || cutoff >= 1) {
    .log_abort("topic_gammafit_thrP must be strictly between 0 and 1.")
  }
  cutoff
}

.module3_default_term_assignment_method <- function(method) {
  value <- tolower(paste(as.character(method), collapse = " "))
  is_aggregate <- grepl("aggr|aggregate", value)
  is_weighted <- grepl("weight", value, fixed = TRUE)
  is_unique <- grepl("uniq|unique", value)
  if (is_aggregate && !is_weighted && !is_unique) "gammafit_maxprob" else "max_phi"
}

.module3_report_state <- function(cfg) {
  report <- cfg$report %||% cfg$report_state %||% list()
  if (!is.list(report)) .log_abort("`report` in the project config must be a mapping.")
  colors <- report$condition_colors %||% cfg$condition_colors %||% list()
  if (is.atomic(colors) && !is.null(names(colors))) colors <- as.list(colors)
  if (!is.list(colors)) .log_abort("`report.condition_colors` must be a named mapping of condition IDs to hex colors.")
  if (length(colors)) {
    if (is.null(names(colors)) || any(!nzchar(names(colors)))) {
      .log_abort("Every `report.condition_colors` entry must have a condition ID.")
    }
    values <- unlist(colors, use.names = FALSE)
    if (any(!grepl("^#[0-9A-Fa-f]{6}$", values))) {
      .log_abort("Every `report.condition_colors` value must be a six-digit hex color such as `#4E79A7`.")
    }
    colors <- stats::setNames(as.list(toupper(values)), names(colors))
  }
  condition_order <- unique(as.character(report$condition_order %||% character()))
  condition_order <- condition_order[!is.na(condition_order) & nzchar(condition_order)]
  defaults <- report$defaults %||% list()
  if (!is.list(defaults)) .log_abort("`report.defaults` in the project config must be a mapping.")
  list(
    condition_colors = colors,
    condition_order = condition_order,
    defaults = defaults
  )
}

.module3_resolve_topic_run_config <- function(project_config = NULL,
                                             method = NULL,
                                             k_grid = NULL,
                                             warplda_iterations = NULL,
                                             topic_link_output = NULL,
                                             count_method = NULL,
                                             count_input = NULL,
                                             topic_score_method = NULL,
                                             topic_term_assignment_method = NULL,
                                             vae_device = NULL,
                                             vae_batch_size = NULL,
                                             pathway_backend = NULL,
                                             pathway_species = NULL,
                                             pathway_databases = NULL) {
  cfg <- .module3_read_project_config(project_config)
  method <- if (is.null(method)) {
    as.character(.module3_cfg_value(cfg, c("topic_method", "module3_topic_method"), "comparison_aggr_multivi"))[[1L]]
  } else {
    as.character(method)
  }
  k_raw <- if (is.null(k_grid)) {
    .module3_cfg_value(cfg, c("topic_k", "module3_topic_k", "topic_k_grid", "module3_topic_k_grid"), 10L)
  } else {
    k_grid
  }
  k_grid <- .module3_cfg_int_vector(k_raw, default = 10L)
  iterations <- if (is.null(warplda_iterations)) {
    .module3_cfg_value(cfg, c("warplda_iterations", "topic_warplda_iterations", "module3_warplda_iterations"), 2000L)
  } else {
    warplda_iterations
  }
  iterations <- suppressWarnings(as.integer(iterations[[1L]]))
  if (!is.finite(iterations) || iterations < 1L) iterations <- 2000L
  link_output <- if (is.null(topic_link_output)) {
    as.character(.module3_cfg_value(cfg, c("topic_link_output", "module3_topic_link_output"), "pass"))[[1L]]
  } else {
    as.character(topic_link_output)[[1L]]
  }
  count_method <- if (is.null(count_method)) {
    as.character(.module3_cfg_value(cfg, c("topic_count_method", "module3_topic_count_method"), "log"))[[1L]]
  } else {
    as.character(count_method)[[1L]]
  }
  if (!count_method %in% c("log", "bin")) count_method <- "log"
  count_input <- if (is.null(count_input)) {
    .module3_cfg_value(cfg, c("topic_count_input", "module3_topic_count_input"), NULL)
  } else {
    count_input
  }
  count_input <- .resolve_topic_count_input(count_method = count_method, count_input = count_input)
  vae_device <- if (is.null(vae_device)) {
    as.character(.module3_cfg_value(cfg, c("topic_vae_device", "module3_topic_vae_device"), "auto"))[[1L]]
  } else {
    as.character(vae_device)[[1L]]
  }
  vae_batch_size <- if (is.null(vae_batch_size)) {
    suppressWarnings(as.integer(.module3_cfg_value(cfg, c("topic_vae_batch_size", "module3_topic_vae_batch_size"), 64L)[[1L]]))
  } else {
    suppressWarnings(as.integer(vae_batch_size[[1L]]))
  }
  if (!is.finite(vae_batch_size) || vae_batch_size < 1L) vae_batch_size <- 64L
  pathway_backend <- if (is.null(pathway_backend)) {
    as.character(.module3_cfg_value(cfg, c("pathway_backend", "module3_pathway_backend"), "enrichly"))[[1L]]
  } else {
    as.character(pathway_backend)[[1L]]
  }
  pathway_backend <- .pathway_backend(pathway_backend)
  pathway_species <- if (is.null(pathway_species)) {
    explicit_species <- .module3_cfg_value(cfg, c("pathway_species", "module3_pathway_species", "pathway_organism", "module3_pathway_organism"), NULL)
    if (!is.null(explicit_species)) {
      explicit_species
    } else {
      .pathway_species_from_ref_genome(.module3_cfg_value(cfg, c("ref_genome", "genome", "reference_genome"), NULL))
    }
  } else {
    pathway_species
  }
  pathway_species <- .normalize_pathway_species_mode(pathway_species)
  pathway_databases <- if (is.null(pathway_databases)) {
    .module3_cfg_value(cfg, c("pathway_databases", "module3_pathway_databases"), NULL)
  } else {
    pathway_databases
  }
  if (!is.null(pathway_databases)) {
    pathway_databases <- unique(trimws(as.character(pathway_databases)))
    pathway_databases <- pathway_databases[!is.na(pathway_databases) & nzchar(pathway_databases)]
    if (!length(pathway_databases)) pathway_databases <- NULL
  }
  gammafit_scope_raw <- .module3_cfg_value(
    cfg,
    c("topic_gammafit_scope", "module3_topic_gammafit_scope", "topic_gammafit_scopes", "module3_topic_gammafit_scopes"),
    "topic_term_group"
  )
  gammafit_scope <- as.character(gammafit_scope_raw)[[1L]]
  if (!gammafit_scope %in% c("topic_term_group", "global_term_group")) {
    gammafit_scope <- "topic_term_group"
  }
  gammafit_thrP <- .module3_resolve_gammafit_thrP(.module3_cfg_value(
    cfg,
    c("topic_gammafit_thrP", "module3_topic_gammafit_thrP"),
    NULL
  ), method = method)
  gammafit_min_terms <- suppressWarnings(as.integer(.module3_cfg_value(
    cfg,
    c("topic_gammafit_min_terms", "module3_topic_gammafit_min_terms"),
    50L
  )[[1L]]))
  if (!is.finite(gammafit_min_terms) || gammafit_min_terms < 1L) {
    gammafit_min_terms <- 50L
  }
  topic_score_method <- if (is.null(topic_score_method)) {
    as.character(.module3_cfg_value(
      cfg,
      c("topic_score_method", "module3_topic_score_method"),
      "normtop_specificity"
    ))[[1L]]
  } else {
    as.character(topic_score_method)[[1L]]
  }
  if (!topic_score_method %in% c("normtop_specificity", "rowmax_phi")) {
    topic_score_method <- "normtop_specificity"
  }
  topic_term_assignment_method <- if (is.null(topic_term_assignment_method)) {
    as.character(.module3_cfg_value(
      cfg,
      c("topic_term_assignment_method", "module3_topic_term_assignment_method"),
      .module3_default_term_assignment_method(method)
    ))[[1L]]
  } else {
    as.character(topic_term_assignment_method)[[1L]]
  }
  if (!topic_term_assignment_method %in% c("gammafit_maxprob", "max_phi", "gammafit")) {
    topic_term_assignment_method <- .module3_default_term_assignment_method(method)
  }
  topic_link_method <- as.character(.module3_cfg_value(
    cfg,
    c("topic_link_method", "module3_topic_link_method"),
    "gammafit"
  ))[[1L]]
  if (!topic_link_method %in% c("gammafit", "theta_and_terms", "gene_prob", "link_score_prob", "link_score_efdr")) {
    topic_link_method <- "gammafit"
  }
  topic_link_prob_cutoff <- .module3_cfg_value(
    cfg,
    c("topic_link_prob_cutoff", "module3_topic_link_prob_cutoff"),
    0.3
  )
  topic_tf_membership_cutoff <- suppressWarnings(as.numeric(.module3_cfg_value(
    cfg,
    c("topic_tf_membership_cutoff", "module3_topic_tf_membership_cutoff"),
    0.3
  )[[1L]]))
  if (!is.finite(topic_tf_membership_cutoff) || topic_tf_membership_cutoff < 0 || topic_tf_membership_cutoff > 1) {
    topic_tf_membership_cutoff <- 0.3
  }
  topic_tf_primary_margin_cutoff <- suppressWarnings(as.numeric(.module3_cfg_value(
    cfg,
    c("topic_tf_primary_margin_cutoff", "module3_topic_tf_primary_margin_cutoff"),
    0.1
  )[[1L]]))
  if (!is.finite(topic_tf_primary_margin_cutoff) || topic_tf_primary_margin_cutoff < 0 || topic_tf_primary_margin_cutoff > 1) {
    topic_tf_primary_margin_cutoff <- 0.1
  }
  run_raw_theta_document_heatmap <- .module3_cfg_value(
    cfg,
    c(
      "topic_raw_theta_document_heatmap",
      "module3_topic_raw_theta_document_heatmap",
      "run_raw_theta_document_heatmap"
    ),
    TRUE
  )
  run_raw_theta_document_heatmap <- .as_logical_flag(run_raw_theta_document_heatmap[[1L]])
  optimize_topics_raw <- .module3_cfg_value(
    cfg,
    c(
      "topic_optimize_topics", "module3_topic_optimize_topics",
      "module3_optimize_topics"
    ),
    NULL
  )
  optimize_topics <- if (is.null(optimize_topics_raw)) {
    NULL
  } else {
    .as_logical_flag(optimize_topics_raw[[1L]])
  }
  topic_merge_min_genes <- suppressWarnings(as.integer(.module3_cfg_value(
    cfg,
    c(
      "topic_merge_min_genes", "module3_topic_merge_min_genes",
      "module3_merge_min_genes"
    ),
    50L
  )[[1L]]))
  topic_merge_min_links <- suppressWarnings(as.integer(.module3_cfg_value(
    cfg,
    c(
      "topic_merge_min_links", "module3_topic_merge_min_links",
      "module3_merge_min_links"
    ),
    200L
  )[[1L]]))
  topic_merge_similarity_threshold <- suppressWarnings(as.numeric(.module3_cfg_value(
    cfg,
    c(
      "topic_merge_similarity_threshold",
      "module3_topic_merge_similarity_threshold",
      "module3_merge_similarity_threshold"
    ),
    0.90
  )[[1L]]))
  run_topic_assignment_qc_raw <- .module3_cfg_value(
    cfg,
    c(
      "topic_assignment_qc", "module3_topic_assignment_qc",
      "module3_run_topic_assignment_qc"
    ),
    NULL
  )
  run_topic_assignment_qc <- if (is.null(run_topic_assignment_qc_raw)) {
    NULL
  } else {
    .as_logical_flag(run_topic_assignment_qc_raw[[1L]])
  }
  topic_qc_umap_links_per_condition <- suppressWarnings(as.integer(
    .module3_cfg_value(
      cfg,
      c(
        "topic_qc_umap_links_per_condition",
        "module3_topic_qc_umap_links_per_condition",
        "module3_qc_umap_links_per_condition"
      ),
      10000L
    )[[1L]]
  ))
  topic_qc_top_tfs <- suppressWarnings(as.integer(.module3_cfg_value(
    cfg,
    c("topic_qc_top_tfs", "module3_topic_qc_top_tfs", "module3_qc_top_tfs"),
    150L
  )[[1L]]))
  topic_qc_seed <- suppressWarnings(as.integer(.module3_cfg_value(
    cfg,
    c("topic_qc_seed", "module3_topic_qc_seed", "module3_qc_seed"),
    20260716L
  )[[1L]]))
  list(
    method = method,
    k_grid = k_grid,
    warplda_iterations = iterations,
    topic_link_output = link_output,
    count_method = count_method,
    count_input = count_input,
    vae_device = vae_device,
    vae_batch_size = vae_batch_size,
    pathway_backend = pathway_backend,
    report_state = .module3_report_state(cfg),
    extraction_args = list(
      binarize_method = "gammafit",
      topic_score_method = topic_score_method,
      topic_term_assignment_method = topic_term_assignment_method,
      gammafit_scope = gammafit_scope,
      thrP = gammafit_thrP,
      in_topic_min_terms = gammafit_min_terms,
      link_topic_method = topic_link_method,
      link_topic_prob_cutoff = topic_link_prob_cutoff,
      topic_tf_membership_cutoff = topic_tf_membership_cutoff,
      topic_tf_primary_margin_cutoff = topic_tf_primary_margin_cutoff,
      pathway_species = pathway_species,
      pathway_databases = pathway_databases,
      theta_umap_condition_colors = .module3_topic_condition_colors(cfg),
      run_raw_theta_document_heatmap = run_raw_theta_document_heatmap,
      optimize_topics = optimize_topics,
      topic_merge_min_genes = topic_merge_min_genes,
      topic_merge_min_links = topic_merge_min_links,
      topic_merge_similarity_threshold = topic_merge_similarity_threshold,
      run_topic_assignment_qc = run_topic_assignment_qc,
      topic_qc_umap_links_per_condition = topic_qc_umap_links_per_condition,
      topic_qc_top_tfs = topic_qc_top_tfs,
      topic_qc_seed = topic_qc_seed
    ),
    benchmark = list(
      enabled = isTRUE(.module3_cfg_value(cfg, c("topic_benchmark_enabled", "module3_topic_benchmark_enabled"), FALSE)),
      methods = .module3_cfg_value(cfg, c("topic_benchmark_methods", "module3_topic_benchmark_methods"), character()),
      k_grid = .module3_cfg_int_vector(.module3_cfg_value(cfg, c("topic_benchmark_k_grid", "module3_topic_benchmark_k_grid"), NULL), default = integer())
    )
  )
}

#' Run topic modeling
#'
#' @description
#' Wrapper function to conduct the full regulatory topic-modeling workflow for
#' one selected topic-document construction method.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   files.
#' @param multiomic_data Optional CraftGRN multiomic object. Required when
#'   `replicate_documents = TRUE`.
#' @param comparisons Comparison or condition grouping table, or a CSV path.
#' @param output_dir Topic output directory.
#' @param project_config Optional project YAML path or config list. When
#'   supplied, `topic_method`, `topic_k` or `topic_k_grid`,
#'   `warplda_iterations`, `topic_link_output`, `topic_score_method`,
#'   `topic_term_assignment_method`,
#'   `topic_gammafit_thrP` (a scalar or an `lda`/`multivi`/`vae_mlp` mapping),
#'   `topic_tf_membership_cutoff`, `topic_tf_primary_margin_cutoff`, and
#'   `topic_raw_theta_document_heatmap` are used for arguments that are left as
#'   `NULL` or as default extraction settings. A nested `report` mapping can
#'   define `condition_colors`, `condition_order`, and initial `defaults` for
#'   generated reports without creating package or browser-persistent state.
#'   Defaults may use exact `condition_1` and `condition_2` IDs or reusable
#'   `condition_1_suffix` and `condition_2_suffix` selectors.
#' @param method Single Module 3 method ID. If `NULL`, read from
#'   `project_config` or use the package default.
#' @param k_grid Integer topic numbers. If `NULL`, read from `project_config`
#'   or use `10`.
#' @param warplda_iterations Number of native WarpLDA iterations. If `NULL`,
#'   read from `project_config` or use `2000`.
#' @param topic_link_output Topic-link output mode. If `NULL`, read from
#'   `project_config` or use `"pass"`.
#' @param count_method Topic count conversion method. If `NULL`, read from
#'   `project_config` or use `"log"`.
#' @param count_input Topic count column for model fitting. If `NULL`, inferred
#'   from `count_method`.
#' @param topic_score_method Topic-term score method for extraction. If `NULL`,
#'   read from `project_config` or use `"normtop_specificity"`.
#' @param topic_term_assignment_method Term-to-topic assignment method. If
#'   `NULL`, read from `project_config`. Aggregate Gene/Peak methods default to
#'   `"gammafit_maxprob"`, which applies GammaFit first and retains a Gene/Peak
#'   pair only when the independently selected maximum-phi passing topics agree.
#'   Unique or aggregate-weight methods retain the independent `"max_phi"`
#'   default.
#' @param optimize_topics Whether eligible condition-topic extractions merge
#'   undersized or highly similar topics before downstream reports. If `NULL`,
#'   use project config or the standard condition-mode default.
#' @param run_topic_assignment_qc Whether to write the standard per-K topic
#'   assignment QC PDF. If `NULL`, use project config.
#' @param topic_merge_min_genes Minimum assigned genes required to retain a
#'   topic without a size-based merge. If `NULL`, use project config.
#' @param topic_merge_min_links Minimum aligned TF-target links required to
#'   retain a topic without a size-based merge. If `NULL`, use project config.
#' @param topic_merge_similarity_threshold Mean Gene/Peak Hellinger similarity
#'   at or above which topics are merged. If `NULL`, use project config.
#' @param topic_qc_umap_links_per_condition Maximum deterministic UMAP sample
#'   size per condition. Full-universe counts are not sampled.
#' @param topic_qc_top_tfs Number of globally ranked TFs shown in the pooled
#'   TF-by-topic QC heatmap.
#' @param vae_device VAE device, for example `"auto"`, `"cpu"`, or `"cuda"`.
#'   If `NULL`, read from `project_config` or use `"auto"`.
#' @param vae_batch_size VAE mini-batch size. If `NULL`, read from
#'   `project_config` or use `64`.
#' @param pathway_backend Pathway enrichment backend. Use `"enrichly"` for
#'   local cached enrichment or `"enrichr"` for the Enrichr web API. If `NULL`,
#'   read from `project_config` or use `"enrichly"`.
#' @param pathway_species Species used to choose default pathway databases.
#'   Supported values include `"human"`, `"mouse"`, and
#'   `"human_mouse_best"`. The best-of mode runs human and mouse pathway
#'   databases separately and reports the better row per topic and normalized
#'   pathway name, ranked by adjusted p-value, logp, combined score, and overlap
#'   size. If `NULL`, read from `project_config` or infer from `ref_genome`.
#' @param pathway_databases Optional Enrichr database names. If `NULL`, read
#'   from `project_config` or use the species-specific package defaults.
#' @param sample_subset Optional condition/sample IDs to retain when building
#'   topic documents. Condition-link runs read only matching manifest rows.
#' @param analysis_label Stable label written to topic input metadata.
#' @param extraction_topic_report_args Optional named list of topic-extraction
#'   report argument overrides. Values here override project config values,
#'   including TF-to-topic assignment cutoffs.
#' @param extraction_k_workers Number of K values to extract concurrently.
#'   Use `1` to request sequential extraction.
#' @param extraction_k_max_workers Maximum concurrent K extraction workers.
#' @param memory_safety Module 3 memory policy. `"strict"` fails before an
#'   unsafe allocation, `"adaptive"` reduces concurrency where possible, and
#'   `"off"` disables memory preflight checks.
#' @param memory_max_fraction Maximum fraction of currently available memory
#'   that a Module 3 stage may use. Defaults to `0.8`.
#' @param input_source Input source type. Use `"differential_links"` for
#'   comparison links and `"condition_links"` for condition-native links.
#' @param ... Additional arguments passed to the internal topic-modeling
#'   wrapper.
#'
#' @return An invisible list with topic input/model/extraction paths, review
#'   outputs, and `qc_report` when requested.
#' @export
run_topic_modeling <- function(filtered_dir,
                               multiomic_data = NULL,
                               comparisons,
                               output_dir,
                               project_config = NULL,
                               method = NULL,
                               k_grid = NULL,
                               warplda_iterations = NULL,
                               topic_link_output = NULL,
                               count_method = NULL,
                               count_input = NULL,
                               topic_score_method = NULL,
                               topic_term_assignment_method = NULL,
                               optimize_topics = NULL,
                               run_topic_assignment_qc = NULL,
                               topic_merge_min_genes = NULL,
                               topic_merge_min_links = NULL,
                               topic_merge_similarity_threshold = NULL,
                               topic_qc_umap_links_per_condition = NULL,
                               topic_qc_top_tfs = NULL,
                               vae_device = NULL,
                               vae_batch_size = NULL,
                               pathway_backend = NULL,
                               pathway_species = NULL,
                               pathway_databases = NULL,
                               sample_subset = NULL,
                               analysis_label = NULL,
                               extraction_topic_report_args = list(),
                               extraction_k_workers = NULL,
                               extraction_k_max_workers = 4L,
                               memory_safety = c("strict", "adaptive", "off"),
                               memory_max_fraction = 0.8,
                               input_source = c("differential_links", "condition_links"),
                               ...) {
  input_source <- match.arg(input_source)
  memory_safety <- match.arg(memory_safety)
  resolved <- .module3_resolve_topic_run_config(
    project_config = project_config,
    method = method,
    k_grid = k_grid,
    warplda_iterations = warplda_iterations,
    topic_link_output = topic_link_output,
    count_method = count_method,
    count_input = count_input,
    topic_score_method = topic_score_method,
    topic_term_assignment_method = topic_term_assignment_method,
    vae_device = vae_device,
    vae_batch_size = vae_batch_size,
    pathway_backend = pathway_backend,
    pathway_species = pathway_species,
    pathway_databases = pathway_databases
  )
  if (!is.null(optimize_topics)) {
    extraction_topic_report_args$optimize_topics <- isTRUE(optimize_topics)
  }
  if (!is.null(run_topic_assignment_qc)) {
    extraction_topic_report_args$run_topic_assignment_qc <-
      isTRUE(run_topic_assignment_qc)
  }
  direct_extraction_controls <- list(
    topic_merge_min_genes = topic_merge_min_genes,
    topic_merge_min_links = topic_merge_min_links,
    topic_merge_similarity_threshold = topic_merge_similarity_threshold,
    topic_qc_umap_links_per_condition = topic_qc_umap_links_per_condition,
    topic_qc_top_tfs = topic_qc_top_tfs
  )
  for (control in names(direct_extraction_controls)) {
    value <- direct_extraction_controls[[control]]
    if (!is.null(value)) extraction_topic_report_args[[control]] <- value
  }
  run_regulatory_topics(
    filtered_dir = filtered_dir,
    multiomic_data = multiomic_data,
    comparisons = comparisons,
    output_dir = output_dir,
    method = resolved$method,
    k_grid = resolved$k_grid,
    warplda_iterations = resolved$warplda_iterations,
    topic_link_output = resolved$topic_link_output,
    count_method = resolved$count_method,
    count_input = resolved$count_input,
    pathway_backend = resolved$pathway_backend,
    report_state = resolved$report_state,
    extraction_topic_report_args = modifyList(resolved$extraction_args, extraction_topic_report_args),
    vae_device = resolved$vae_device,
    vae_batch_size = resolved$vae_batch_size,
    input_source = input_source,
    sample_subset = sample_subset,
    analysis_label = analysis_label,
    extraction_k_workers = extraction_k_workers,
    extraction_k_max_workers = extraction_k_max_workers,
    memory_safety = memory_safety,
    memory_max_fraction = memory_max_fraction,
    ...
  )
}
