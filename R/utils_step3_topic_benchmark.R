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

.module3_topic_method_plan <- function(methods = "condition_aggr_weight_lda",
                                       k_grid = 10L) {
  dict <- .m3tb_method_dictionary()
  if (identical(methods, "default")) methods <- "condition_aggr_weight_lda"
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
    slug <- .m3tb_backend_slug(out$backend[[1L]], out$vae_variant[[1L]])
    out[, `:=`(
      run_id = "selected",
      run_slug = .m3tb_method_slug(out[1]),
      run_dir = output_dir,
      topic_documents_dir = file.path(output_dir, "topic_documents"),
      topic_models_dir = file.path(output_dir, "topic_models", slug),
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

.m3tb_review_dir <- function(output_dir, output_layout) {
  if (identical(output_layout, "legacy")) {
    file.path(output_dir, "review_topic_experiments")
  } else {
    file.path(output_dir, "review")
  }
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
  x <- gsub("_", " ", x, fixed = TRUE)
  trimws(x)
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
  comparison[, `:=`(
    context_type = "comparison",
    comparison_label = comparison_label,
    display_label = .m3tb_display_label(comparison_label),
    metric_group = paste(comparison$comparison_id, comparison$direction, sep = "::")
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
    return(data.table::data.table(
      context_type = "condition",
      comparison_label = conds,
      display_label = .m3tb_display_label(conds),
      metric_group = conds
    ))
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
    design[context_type == row$context_type[[1L]], .(comparison_label, metric_group)],
    by = "comparison_label",
    all.x = TRUE,
    sort = FALSE
  )
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
      ".csv"
    )
  )
  data.table::fwrite(data.table::data.table(label = rownames(dmat), as.data.frame(dmat, check.names = FALSE)), jsd_path)
  groups <- merged$metric_group
  per_label <- lapply(seq_len(nrow(dmat)), function(i) {
    same <- setdiff(which(groups == groups[[i]]), i)
    other <- which(groups != groups[[i]])
    inner <- if (length(same)) mean(dmat[i, same], na.rm = TRUE) else 0
    inter <- if (length(other)) mean(dmat[i, other], na.rm = TRUE) else NA_real_
    score <- (inter - inner) / (inter + inner + 1e-12)
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
    n_metric_groups = data.table::uniqueN(metric_group),
    mean_inner_distance = mean(inner_distance, na.rm = TRUE),
    mean_inter_distance = mean(inter_distance, na.rm = TRUE),
    theta_condition_separation_score = mean(theta_condition_label_score, na.rm = TRUE)
  ), by = .(method_order, context_type, setup, setup_label, model_label, method_setup, k)]
  list(score = score, per_label = per_label, profiles = profiles)
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

.m3tb_score_theta_outputs <- function(output_dir,
                                      model_rows,
                                      comparisons,
                                      score_prefix = "theta_condition_separation",
                                      replicate_documents = FALSE,
                                      multiomic_data = NULL,
                                      review_dir = NULL) {
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  csv_dir <- file.path(review_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  design <- .m3tb_design_table(
    comparisons,
    replicate_documents = replicate_documents,
    multiomic_data = multiomic_data
  )
  data.table::fwrite(design, file.path(csv_dir, paste0(score_prefix, "_design.csv")))
  payload <- lapply(seq_len(nrow(model_rows)), function(i) {
    row <- model_rows[i]
    theta_path <- file.path(row$model_dir[[1L]], "vae_models", sprintf("theta_K%d.csv", as.integer(row$selected_k[[1L]])))
    theta <- .m3tb_read_theta(theta_path)
    .m3tb_score_theta_one(theta, row, design, csv_dir)
  })
  scores <- data.table::rbindlist(lapply(payload, `[[`, "score"), use.names = TRUE, fill = TRUE)
  per_label <- data.table::rbindlist(lapply(payload, `[[`, "per_label"), use.names = TRUE, fill = TRUE)
  data.table::setorder(scores, method_order, k)
  data.table::setorder(per_label, method_order, k, comparison_label)
  data.table::fwrite(scores, file.path(csv_dir, paste0(score_prefix, "_scores.csv")))
  data.table::fwrite(per_label, file.path(csv_dir, paste0(score_prefix, "_per_label.csv")))
  data.table::fwrite(per_label, file.path(csv_dir, paste0(score_prefix, "_score_long.csv")))
  heatmap_values <- scores[, .(
    method_order,
    context_type,
    method_setup,
    setup,
    setup_label,
    model_label,
    k,
    theta_condition_separation_score
  )]
  data.table::fwrite(heatmap_values, file.path(csv_dir, paste0(score_prefix, "_score_heatmap_values.csv")))
  wide <- data.table::dcast(
    scores[, .(method_order, context_type, method_setup, k = as.integer(k), theta_condition_separation_score)],
    method_order + context_type + method_setup ~ k,
    value.var = "theta_condition_separation_score"
  )
  k_cols <- as.character(sort(unique(scores$k)))
  data.table::setnames(wide, intersect(k_cols, names(wide)), paste0("K", intersect(k_cols, names(wide))))
  data.table::setorder(wide, method_order)
  data.table::fwrite(wide, file.path(csv_dir, paste0(score_prefix, "_score_heatmap_values_matrix.csv")))
  list(scores = scores, per_label = per_label, matrix = wide, design = design, score_prefix = score_prefix)
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
  if (any(keep)) files[keep] else character()
}

.m3tb_empty_pass_counts <- function() {
  data.table::data.table(
    method_order = integer(),
    method_setup = character(),
    setup = character(),
    model_label = character(),
    selected_k = integer(),
    status = character(),
    count = integer()
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

.m3tb_summarize_topic_links <- function(output_dir, model_rows, review_dir = NULL) {
  pass_rows <- list()
  shared_rows <- list()
  seen <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(model_rows))) {
    row <- model_rows[i]
    files <- .m3tb_find_topic_links(output_dir, row)
    if (!length(files)) next
    for (path in files) {
      key <- paste(row$method[[1L]], row$selected_k[[1L]], normalizePath(path, winslash = "/", mustWork = FALSE), sep = "\r")
      if (exists(key, envir = seen, inherits = FALSE)) next
      assign(key, TRUE, envir = seen)
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
      link_status <- data.table::data.table(link_id = link_id, pass = dt$link_pass | (dt$peak_pass & dt$gene_pass))
      link_status <- link_status[, .(pass = any(pass, na.rm = TRUE)), by = link_id]
      pass_rows[[length(pass_rows) + 1L]] <- data.table::data.table(
        method_order = row$method_order[[1L]],
        method_setup = row$method_setup[[1L]],
        setup = row$setup[[1L]],
        model_label = row$model_label[[1L]],
        selected_k = as.integer(row$selected_k[[1L]]),
        status = c("Pass", "Fail"),
        count = c(sum(link_status$pass, na.rm = TRUE), sum(!link_status$pass, na.rm = TRUE))
      )
      keep <- dt[link_pass | (peak_pass & gene_pass)]
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
  }
  pass <- if (length(pass_rows)) {
    data.table::rbindlist(pass_rows, use.names = TRUE, fill = TRUE)
  } else {
    .m3tb_empty_pass_counts()
  }
  shared <- if (length(shared_rows)) {
    data.table::rbindlist(shared_rows, use.names = TRUE, fill = TRUE)
  } else {
    .m3tb_empty_shared_counts()
  }
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  csv_dir <- file.path(review_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(pass, file.path(csv_dir, "topic_setup_pass_state_counts.csv"))
  data.table::fwrite(shared, file.path(csv_dir, "topic_setup_shared_topic_counts.csv"))
  data.table::fwrite(pass, file.path(csv_dir, "tf_std_six_setups_pass_state_counts.csv"))
  data.table::fwrite(shared, file.path(csv_dir, "tf_std_six_setups_shared_topic_counts.csv"))
  list(pass = pass, shared = shared)
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

.m3tb_write_review_html <- function(output_dir, score_result, link_summary, review_dir = NULL) {
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  html_dir <- file.path(review_dir, "html")
  dir.create(html_dir, recursive = TRUE, showWarnings = FALSE)
  matrix_html <- paste(utils::capture.output(print(score_result$matrix)), collapse = "\n")
  pass_html <- paste(utils::capture.output(print(utils::head(link_summary$pass, 20L))), collapse = "\n")
  shared_html <- paste(utils::capture.output(print(utils::head(link_summary$shared, 20L))), collapse = "\n")
  base_doc <- function(title, body) {
    c(
      "<!doctype html>",
      "<html><head><meta charset=\"utf-8\">",
      paste0("<title>", .m3tb_html_escape(title), "</title>"),
      "<style>body{font-family:Arial,sans-serif;margin:28px;color:#1f2933}h1{font-size:24px}pre{background:#f6f8fa;padding:14px;border:1px solid #d0d7de;overflow:auto}.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}</style>",
      "</head><body>",
      paste0("<h1>", .m3tb_html_escape(title), "</h1>"),
      body,
      "</body></html>"
    )
  }
  theta_body <- paste0(
    "<p>Condition separation scores are computed from Jensen-Shannon distances between mean theta profiles.</p>",
    "<div class=\"grid\"><section><h2>Score matrix</h2><pre>", .m3tb_html_escape(matrix_html), "</pre></section>",
    "<section><h2>Topic-link pass states</h2><pre>", .m3tb_html_escape(pass_html), "</pre></section></div>"
  )
  method_body <- paste0(
    "<p>Topic method and K report for completed model outputs.</p>",
    "<div class=\"grid\"><section><h2>Pass states</h2><pre>", .m3tb_html_escape(pass_html), "</pre></section>",
    "<section><h2>Shared topic counts</h2><pre>", .m3tb_html_escape(shared_html), "</pre></section></div>"
  )
  files <- c(
    theta_phi = file.path(html_dir, "theta_phi_and_group_mds.html"),
    method = file.path(html_dir, "topic_method_k_topic_mds_report.html"),
    method_global = file.path(html_dir, "topic_method_k_topic_mds_report_global_term_group.html")
  )
  writeLines(base_doc("Theta, phi, and group MDS review", theta_body), files[["theta_phi"]], useBytes = TRUE)
  writeLines(base_doc("Topic method K report", method_body), files[["method"]], useBytes = TRUE)
  writeLines(base_doc("Topic method K report - global term group", method_body), files[["method_global"]], useBytes = TRUE)
  files
}

.m3tb_write_review_outputs <- function(output_dir, score_result, link_summary, review_dir = NULL) {
  if (is.null(review_dir)) review_dir <- file.path(output_dir, "review_topic_experiments")
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
  .m3tb_plot_count_pdf(link_summary$pass, file.path(review_dir, "tf_std_six_setups_pass_state_counts.pdf"), "Topic-link pass-state counts")
  .m3tb_plot_count_pdf(link_summary$shared, file.path(review_dir, "tf_std_six_setups_shared_topic_counts.pdf"), "Shared topic counts", x_col = "method_setup")
  html <- .m3tb_write_review_html(output_dir, score_result, link_summary, review_dir = review_dir)
  invisible(html)
}

#' Run a Module 3 topic benchmark and review report
#'
#' Runs or reviews Module 3 topic models using the package-native benchmark
#' layout. The function can train/extract models, or it can score and report
#' existing theta/phi and topic-link outputs when `run_training = FALSE` and
#' `run_extraction = FALSE`.
#'
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   files.
#' @param multiomic_data Optional CraftGRN multiomic object. Required when
#'   `replicate_documents = TRUE`.
#' @param comparisons Comparison or condition grouping table, or a CSV path.
#'   For condition separation scoring, use columns `condition_label` and
#'   `condition_group`.
#' @param output_dir Topic benchmark output directory.
#' @param methods Character vector of benchmark method IDs, `"default"`, or
#'   `"all"`.
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
#' @param sample_subset Optional condition/sample labels passed to
#'   [train_topic_models()]. When supplied, only comparisons whose case and
#'   control labels are both in this vector are used.
#' @param analysis_label Label used to name the topic-model analysis.
#' @param extraction_topic_report_args Optional named list of topic-extraction
#'   report argument overrides.
#' @param run_training Train topic models before reporting.
#' @param run_extraction Run topic extraction before reporting.
#' @param run_reports Build score tables and review reports.
#' @param verbose Emit concise progress messages.
#'
#' @return A list with method plan, discovered model rows, score results, and
#'   review report paths.
#' @export
run_module3_topic_benchmark <- function(filtered_dir,
                                        multiomic_data = NULL,
                                        comparisons,
                                        output_dir,
                                        methods = "condition_aggr_weight_lda",
                                        k_grid = 10L,
                                        output_layout = c("auto", "standard", "benchmark", "legacy"),
                                        replicate_documents = FALSE,
                                        reuse_if_exists = TRUE,
                                        local_threads = NULL,
                                        sample_subset = NULL,
                                        analysis_label = NULL,
                                        extraction_topic_report_args = list(),
                                        run_training = TRUE,
                                        run_extraction = TRUE,
                                        run_reports = TRUE,
                                        verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[is.finite(k_grid) & k_grid > 1L]
  if (!length(k_grid)) .log_abort("k_grid must include at least one integer greater than 1.")
  method_plan <- .module3_topic_method_plan(methods = methods, k_grid = k_grid)
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
  csv_dir <- file.path(review_dir, "csv")
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
    for (i in seq_len(nrow(method_plan))) {
      row <- method_plan[i]
      if (isTRUE(verbose)) .log_inform("Training Module 3 topic method: {row$method_setup}.")
      train_topic_models(
        Kgrid = k_grid,
        input_dir = filtered_dir,
        output_dir = row$topic_models_dir[[1L]],
        sample_subset = sample_subset,
        analysis_label = analysis_label,
        doc_mode = "tf",
        doc_design = row$doc_design[[1L]],
        fp_term_mode = row$fp_mode[[1L]],
        gene_term_mode = "unique",
        include_tf_terms = TRUE,
        count_input = "pseudo_count_log",
        backend = row$backend[[1L]],
        vae_variant = row$vae_variant[[1L]],
        reuse_if_exists = reuse_if_exists,
        local_threads = local_threads,
        save_full_doc_term_csv = FALSE,
        topic_report_args = list()
      )
    }
  }
  if (isTRUE(run_extraction)) {
    link_topic_n_cores <- if (is.null(local_threads)) {
      .available_cores(logical = TRUE)
    } else {
      max(1L, as.integer(local_threads[[1L]]))
    }
    for (i in seq_len(nrow(method_plan))) {
      row <- method_plan[i]
      model_root <- row$topic_models_dir[[1L]]
      for (k in k_grid) {
        extract_root <- file.path(row$topic_extraction_dir[[1L]], paste0("K", k))
        if (isTRUE(verbose)) {
          .log_inform("Extracting Module 3 topics for method: {row$method_setup}; K={k}.")
        }
        extract_regulatory_topics(
          k = k,
          model_dir = model_root,
          output_dir = extract_root,
          backend = row$backend[[1L]],
          vae_variant = row$vae_variant[[1L]],
          doc_mode = "tf",
          weight_label = row$weight_label[[1L]],
          flatten_single_output = FALSE,
          topic_report_args = modifyList(
            list(
              fp_term_mode = row$fp_mode[[1L]],
              link_topic_n_cores = link_topic_n_cores
            ),
            extraction_topic_report_args
          )
        )
      }
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
  if (isTRUE(run_reports)) {
    if (isTRUE(verbose)) .log_inform("Scoring Module 3 theta separation and writing review reports.")
    score_prefix <- if (isTRUE(replicate_documents)) {
      "theta_condition_replicate_separation"
    } else {
      "theta_condition_separation"
    }
    score_result <- .m3tb_score_theta_outputs(
      output_dir,
      model_rows,
      comparisons,
      score_prefix = score_prefix,
      replicate_documents = replicate_documents,
      multiomic_data = multiomic_data,
      review_dir = review_dir
    )
    link_summary <- .m3tb_summarize_topic_links(output_dir, model_rows, review_dir = review_dir)
    html <- .m3tb_write_review_outputs(output_dir, score_result, link_summary, review_dir = review_dir)
  }
  invisible(list(
    method_plan = method_plan,
    model_rows = model_rows,
    score = score_result,
    topic_link_summary = link_summary,
    html = html,
    output_layout = output_layout,
    review_dir = review_dir,
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
#' @param tf_cluster_map Named vector mapping TF names to motif clusters.
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
#' @param overwrite If FALSE, reuse an existing complete cache.
#' @param verbose Emit concise progress messages.
#'
#' @return A list with cache paths and input summary counts.
#' @export
module3_prepare_topic_inputs <- function(filtered_dir,
                                         output_dir,
                                         tf_cluster_map,
                                         doc_mode = c("tf", "tf_cluster"),
                                         doc_design = c("condition", "comparison"),
                                         fp_term_mode = c("aggregate_weight", "aggregate", "unique"),
                                         gene_term_mode = c("unique", "aggregate"),
                                         sample_subset = NULL,
                                         analysis_label = NULL,
                                         count_method = c("bin", "log"),
                                         count_scale = 50,
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
                                         overwrite = FALSE,
                                         verbose = TRUE) {
  .assert_pkg("data.table")
  .assert_pkg("Matrix")
  doc_mode <- match.arg(doc_mode)
  doc_design <- match.arg(doc_design)
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  gene_term_mode <- match.arg(gene_term_mode)
  count_method <- match.arg(count_method)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rds_dir <- file.path(output_dir, "rds")
  required_cache <- file.path(rds_dir, c("edges_docs.rds", "doc_term.rds", "dtm.rds", "dtm_index.rds"))
  summary_path <- file.path(output_dir, "topic_input_summary.csv")
  if (!isTRUE(overwrite) && all(file.exists(required_cache)) && file.exists(summary_path)) {
    if (isTRUE(verbose)) .log_inform("Reusing existing Module 3 topic input cache: {output_dir}")
    summary_dt <- data.table::fread(summary_path, showProgress = FALSE)
    return(invisible(list(output_dir = output_dir, summary = summary_dt, reused = TRUE)))
  }
  delta_files <- list.files(filtered_dir, "_filtered_links(_(up|down))?\\.csv$", full.names = TRUE)
  if (!length(delta_files)) {
    delta_files <- list.files(filtered_dir, "_delta_links_filtered(_(up|down))?\\.csv$", full.names = TRUE)
  }
  if (!length(delta_files)) {
    delta_files <- list.files(filtered_dir, "_delta_links\\.csv$", full.names = TRUE)
  }
  if (!length(delta_files)) .log_abort("No Module 3 filtered link files found in {filtered_dir}.")
  if (isTRUE(verbose)) .log_inform("Loading {length(delta_files)} Module 3 filtered-link file(s).")
  edges_dt <- data.table::as.data.table(load_delta_links_many(delta_files, keep_original = FALSE))
  n_loaded <- nrow(edges_dt)
  if (!("comparison_id" %in% names(edges_dt))) .log_abort("Module 3 links are missing comparison_id.")
  sample_subset <- if (is.null(sample_subset)) NULL else unique(as.character(sample_subset))
  sample_subset <- sample_subset[!is.na(sample_subset) & nzchar(sample_subset)]
  if (length(sample_subset)) {
    if (!all(c("cond1_id", "cond2_id") %in% names(edges_dt))) {
      .log_abort("sample_subset requires cond1_id and cond2_id columns in Module 3 links.")
    }
    edges_dt <- edges_dt[cond1_id %in% sample_subset & cond2_id %in% sample_subset]
  }
  if (!nrow(edges_dt)) .log_abort("No Module 3 links remain after sample subsetting.")
  edges_filt <- filter_edges_for_tf_topics(
    edges_dt,
    abs_log2fc_fp_min = abs_log2fc_fp_min,
    abs_delta_fp_min = abs_delta_fp_min,
    abs_log2fc_gene_min = abs_log2fc_gene_min,
    require_fp_bound_either = require_fp_bound_either,
    require_tf_expr_either = require_tf_expr_either,
    require_gene_expr_either = require_gene_expr_either,
    direction_consistency = direction_consistency
  )
  if (!nrow(edges_filt)) .log_abort("No Module 3 links passed topic-input filters.")
  if (identical(doc_design, "condition")) {
    edges_docs <- add_condition_tf_docs(edges_filt, tf_cluster_map = tf_cluster_map, doc_mode = doc_mode)
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
      fp_term_mode = fp_term_mode
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
      require_condition_thresholds = identical(doc_mode, "tf")
    )
  }
  if (!nrow(doc_term)) .log_abort("Module 3 document-term table is empty.")
  write_doc_term_cache(doc_term, out_dir = output_dir, save_full_doc_term_csv = isTRUE(save_full_doc_term_csv))
  .save_all(output_dir, "edges_filtered", edges_filt)
  .save_all(output_dir, "edges_docs", edges_docs)
  .save_all(output_dir, "doc_term", doc_term)
  dtm_obj <- build_sparse_dtm(doc_term, count_col = "pseudo_count")
  .save_all(output_dir, "dtm", dtm_obj$dtm)
  .save_all(output_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))
  summary_dt <- data.table::data.table(
    analysis_label = if (is.null(analysis_label)) NA_character_ else as.character(analysis_label[[1L]]),
    doc_design = doc_design,
    doc_mode = doc_mode,
    fp_term_mode = fp_term_mode,
    n_link_rows_loaded = as.double(n_loaded),
    n_link_rows_after_subset = as.double(nrow(edges_dt)),
    n_link_rows_after_filter = as.double(nrow(edges_filt)),
    n_document_edge_rows = as.double(nrow(edges_docs)),
    n_doc_term_rows = as.double(nrow(doc_term)),
    n_documents = as.double(data.table::uniqueN(doc_term$doc_id)),
    n_terms = as.double(data.table::uniqueN(doc_term$term_id)),
    n_nonzero = as.double(Matrix::nnzero(dtm_obj$dtm))
  )
  data.table::fwrite(summary_dt, summary_path)
  if (isTRUE(verbose)) {
    .log_inform("Prepared Module 3 topic inputs: {summary_dt$n_documents} document(s), {summary_dt$n_terms} term(s), {summary_dt$n_doc_term_rows} doc-term row(s).")
  }
  invisible(list(output_dir = output_dir, summary = summary_dt, reused = FALSE))
}

#' Build a Module 3 QC HTML report
#'
#' Writes a lightweight self-contained HTML report for Module 3 topic-model
#' outputs. The report summarizes topic-input caches, model rows, theta
#' separation scores, and compact topic-link pass counts when available.
#'
#' @param topic_dir Module 3 topic output directory.
#' @param output_dir Directory where the report is written. Defaults to
#'   `topic_dir/reports`.
#' @param title Report title.
#' @param verbose Emit concise progress messages.
#'
#' @return Path to the HTML report.
#' @export
build_module3_qc_report <- function(topic_dir,
                                    output_dir = file.path(topic_dir, "reports"),
                                    title = "Module 3 QC report",
                                    verbose = TRUE) {
  .assert_pkg("data.table")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(output_dir, "module3_qc_report.html")
  review_csv <- file.path(topic_dir, "review", "csv")
  if (!dir.exists(review_csv)) review_csv <- file.path(topic_dir, "review_topic_experiments", "csv")
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
    "</div>",
    "<h2>Topic Input Summary</h2>", table_html(input_summary),
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

#' Run regulatory topic modeling
#'
#' Production-oriented Module 3 wrapper for one selected topic-document method.
#' It uses the clean standard output layout, compact topic-link output by
#' default, and writes a Module 3 QC report.
#'
#' @inheritParams run_module3_topic_benchmark
#' @param method Single Module 3 method ID.
#' @param topic_link_output Topic-link output mode. `"pass"` writes compact
#'   passing links and summaries only; `"full"` writes exhaustive all-topic
#'   links; `"both"` writes both; `"none"` disables topic-link export.
#' @param build_qc_report Whether to build the Module 3 QC report.
#'
#' @return The result list from `run_module3_topic_benchmark()`, with
#'   `qc_report` added when requested.
#' @export
run_regulatory_topics <- function(filtered_dir,
                                  multiomic_data = NULL,
                                  comparisons,
                                  output_dir,
                                  method = "condition_aggr_weight_lda",
                                  k_grid = 10L,
                                  replicate_documents = FALSE,
                                  reuse_if_exists = TRUE,
                                  local_threads = NULL,
                                  sample_subset = NULL,
                                  analysis_label = NULL,
                                  topic_link_output = c("pass", "full", "both", "none"),
                                  extraction_topic_report_args = list(),
                                  run_training = TRUE,
                                  run_extraction = TRUE,
                                  run_reports = TRUE,
                                  build_qc_report = TRUE,
                                  verbose = TRUE) {
  topic_link_output <- match.arg(topic_link_output)
  if (length(method) != 1L) .log_abort("run_regulatory_topics() expects one selected method.")
  extraction_args <- modifyList(
    list(
      link_topic_output = topic_link_output,
      pathway_link_scores_file = if (topic_link_output %in% c("pass", "both")) "topic_links_pass.csv" else "topic_links.csv"
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
    sample_subset = sample_subset,
    analysis_label = analysis_label,
    extraction_topic_report_args = extraction_args,
    run_training = run_training,
    run_extraction = run_extraction,
    run_reports = run_reports,
    verbose = verbose
  )
  if (isTRUE(build_qc_report)) {
    res$qc_report <- build_module3_qc_report(output_dir, verbose = verbose)
  }
  invisible(res)
}
