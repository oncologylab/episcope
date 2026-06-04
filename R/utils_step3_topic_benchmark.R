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

.m3tb_design_table <- function(comparisons) {
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

.m3tb_score_theta_outputs <- function(output_dir, model_rows, comparisons) {
  csv_dir <- file.path(output_dir, "review_topic_experiments", "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  design <- .m3tb_design_table(comparisons)
  data.table::fwrite(design, file.path(csv_dir, "theta_condition_separation_design.csv"))
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
  data.table::fwrite(scores, file.path(csv_dir, "theta_condition_separation_scores.csv"))
  data.table::fwrite(per_label, file.path(csv_dir, "theta_condition_separation_per_label.csv"))
  data.table::fwrite(per_label, file.path(csv_dir, "theta_condition_separation_score_long.csv"))
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
  data.table::fwrite(heatmap_values, file.path(csv_dir, "theta_condition_separation_score_heatmap_values.csv"))
  wide <- data.table::dcast(
    scores[, .(method_order, context_type, method_setup, k = as.integer(k), theta_condition_separation_score)],
    method_order + context_type + method_setup ~ k,
    value.var = "theta_condition_separation_score"
  )
  k_cols <- as.character(sort(unique(scores$k)))
  data.table::setnames(wide, intersect(k_cols, names(wide)), paste0("K", intersect(k_cols, names(wide))))
  data.table::setorder(wide, method_order)
  data.table::fwrite(wide, file.path(csv_dir, "theta_condition_separation_score_heatmap_values_matrix.csv"))
  list(scores = scores, per_label = per_label, matrix = wide, design = design)
}

.m3tb_find_topic_links <- function(output_dir, row) {
  root <- file.path(output_dir, row$setup[[1L]], "02_topic_extraction")
  if (!dir.exists(root)) return(character())
  files <- list.files(root, "topic_links[.]csv$", recursive = TRUE, full.names = TRUE)
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

.m3tb_summarize_topic_links <- function(output_dir, model_rows) {
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
  csv_dir <- file.path(output_dir, "review_topic_experiments", "csv")
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

.m3tb_write_review_html <- function(output_dir, score_result, link_summary) {
  html_dir <- file.path(output_dir, "review_topic_experiments", "html")
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

.m3tb_write_review_outputs <- function(output_dir, score_result, link_summary) {
  review_dir <- file.path(output_dir, "review_topic_experiments")
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
  .m3tb_plot_count_pdf(link_summary$pass, file.path(review_dir, "tf_std_six_setups_pass_state_counts.pdf"), "Topic-link pass-state counts")
  .m3tb_plot_count_pdf(link_summary$shared, file.path(review_dir, "tf_std_six_setups_shared_topic_counts.pdf"), "Shared topic counts", x_col = "method_setup")
  html <- .m3tb_write_review_html(output_dir, score_result, link_summary)
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
#' @param multiomic_data Optional CraftGRN multiomic object. Reserved for future
#'   package-level replicate document expansion.
#' @param comparisons Comparison or condition grouping table, or a CSV path.
#'   For condition separation scoring, use columns `condition_label` and
#'   `condition_group`.
#' @param output_dir Topic benchmark output directory.
#' @param methods Character vector of benchmark method IDs, `"default"`, or
#'   `"all"`.
#' @param k_grid Integer topic numbers.
#' @param replicate_documents Whether documents are replicate resolved.
#' @param reuse_if_exists Reuse existing model outputs where possible.
#' @param local_threads Optional thread count for model training.
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
                                        replicate_documents = FALSE,
                                        reuse_if_exists = TRUE,
                                        local_threads = NULL,
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
  csv_dir <- file.path(output_dir, "review_topic_experiments", "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(method_plan, file.path(csv_dir, "module3_topic_method_plan.csv"))
  if (isTRUE(verbose)) {
    .log_inform("Module 3 topic benchmark methods: {nrow(method_plan)}.")
  }
  if (isTRUE(run_training)) {
    for (i in seq_len(nrow(method_plan))) {
      row <- method_plan[i]
      if (isTRUE(verbose)) .log_inform("Training Module 3 topic method: {row$method_setup}.")
      train_topic_models(
        Kgrid = k_grid,
        input_dir = filtered_dir,
        output_dir = file.path(output_dir, row$setup[[1L]], "01_topic_models", row$combo_id[[1L]]),
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
      model_root <- file.path(output_dir, row$setup[[1L]], "01_topic_models", row$combo_id[[1L]])
      for (k in k_grid) {
        extract_root <- file.path(output_dir, row$setup[[1L]], "02_topic_extraction", paste0(row$method[[1L]], "_K", k))
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
  model_rows <- .m3tb_existing_model_rows(output_dir, method_plan, k_grid)
  score_result <- NULL
  link_summary <- NULL
  html <- character()
  if (isTRUE(run_reports)) {
    if (isTRUE(verbose)) .log_inform("Scoring Module 3 theta separation and writing review reports.")
    score_result <- .m3tb_score_theta_outputs(output_dir, model_rows, comparisons)
    link_summary <- .m3tb_summarize_topic_links(output_dir, model_rows)
    html <- .m3tb_write_review_outputs(output_dir, score_result, link_summary)
  }
  invisible(list(
    method_plan = method_plan,
    model_rows = model_rows,
    score = score_result,
    topic_link_summary = link_summary,
    html = html,
    replicate_documents = isTRUE(replicate_documents),
    multiomic_data_provided = !is.null(multiomic_data)
  ))
}
