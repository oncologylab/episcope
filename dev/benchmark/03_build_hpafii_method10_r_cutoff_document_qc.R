#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(args[[1L]]) else "all"
allowed_stages <- c("all", "footprint", "eligibility", "documents", "report")
if (!stage %in% allowed_stages) {
  stop("Stage must be one of: ", paste(allowed_stages, collapse = ", "))
}

required_packages <- c(
  "arrow", "craftgrn", "cowplot", "data.table", "digest", "ggplot2",
  "patchwork", "scales", "tidyselect"
)
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

parse_bool <- function(name, default = TRUE) {
  value <- tolower(trimws(Sys.getenv(
    name,
    unset = if (isTRUE(default)) "true" else "false"
  )))
  if (value %in% c("1", "true", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "no", "n")) return(FALSE)
  stop(name, " must be true or false.")
}

log_info <- function(...) {
  craftgrn:::.log_inform(paste0(...))
}

read_parquet_dt <- function(path, columns = NULL) {
  if (!file.exists(path)) stop("Missing Parquet input: ", path)
  value <- if (is.null(columns)) {
    arrow::read_parquet(path, as_data_frame = TRUE)
  } else {
    arrow::read_parquet(
      path,
      col_select = tidyselect::all_of(columns),
      as_data_frame = TRUE
    )
  }
  data.table::as.data.table(value)
}

write_parquet_dt <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(
    data.table::as.data.table(value),
    path,
    compression = "zstd"
  )
  invisible(path)
}

project_root <- normalizePath(
  Sys.getenv(
    "CRAFTGRN_METHOD10_PROJECT",
    unset = paste0(
      "/data/homes/yl814/episcope_test/",
      "nutrient_stress_strict_JASPAR2024_expanded"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
analysis_root <- file.path(
  project_root,
  "regulatory_topics_hpafii_condition_models_de_gene_filtered"
)
source_run <- file.path(
  project_root,
  "runs",
  "delta_fp_0p2_20260624_180950"
)
module1_root <- file.path(source_run, "predict_tf_binding_sites")
module2_root <- file.path(source_run, "connect_tf_target_genes")
condition_root <- file.path(analysis_root, "condition_links")
output_root <- Sys.getenv(
  "CRAFTGRN_METHOD10_OUTPUT",
  unset = file.path(
    analysis_root,
    "topic_runs",
    "review",
    "method10_module1_module2_r_cutoff_sensitivity"
  )
)
output_root <- path.expand(output_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

reuse <- parse_bool("CRAFTGRN_METHOD10_REUSE", default = TRUE)
minimum_r <- 0.5
module1_cutoffs <- c(0.5, 0.7)
module2_cutoffs <- c(0.5, 0.7)
score_floor <- 2

stats_manifest_path <- file.path(
  module1_root,
  "module1_prediction_stats_chunks",
  "module1_prediction_stats_manifest.csv"
)
link_manifest_path <- file.path(
  module2_root,
  "data",
  "links",
  "module2_links_manifest.csv"
)
condition_manifest_path <- file.path(
  condition_root,
  "condition_links_manifest.csv"
)
tf_corr_path <- file.path(
  module2_root,
  "data",
  "correlations",
  "tf_target_corr.parquet"
)
fp_corr_path <- file.path(
  module2_root,
  "data",
  "correlations",
  "fp_target_corr.parquet"
)
expression_path <- file.path(condition_root, "condition_gene_expression.csv")
multiomic_path <- file.path(
  project_root,
  "predict_tf_binding_sites",
  "01_multiomic_data_object_JASPAR2024.rds"
)
raw_score_path <- file.path(
  project_root,
  "cache",
  "fp_scores_JASPAR2024.csv"
)
project_config_path <- file.path(project_root, "project.yaml")

required_inputs <- c(
  stats_manifest_path, link_manifest_path, condition_manifest_path,
  tf_corr_path, fp_corr_path, expression_path, multiomic_path,
  raw_score_path, project_config_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))
}

file_identity <- function(path) {
  info <- file.info(path)
  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size = as.numeric(info$size),
    mtime = as.numeric(info$mtime)
  )
}

source_signature <- digest::digest(
  list(
    project_config_md5 = unname(tools::md5sum(project_config_path)),
    stats_manifest_md5 = unname(tools::md5sum(stats_manifest_path)),
    link_manifest_md5 = unname(tools::md5sum(link_manifest_path)),
    condition_manifest_md5 = unname(tools::md5sum(condition_manifest_path)),
    tf_corr = file_identity(tf_corr_path),
    fp_corr = file_identity(fp_corr_path),
    expression = file_identity(expression_path),
    multiomic = file_identity(multiomic_path),
    raw_score = file_identity(raw_score_path),
    minimum_r = minimum_r,
    method10 = list(
      document = "condition::TF",
      gene_weight = "condition expression specificity",
      peak_term = "unique coordinate",
      peak_weight = "TF expression",
      peak_gene_token_ratio = 0.5
    )
  ),
  algo = "xxhash64"
)
signature_path <- file.path(output_root, "input_signature.txt")
if (file.exists(signature_path)) {
  old_signature <- trimws(readLines(signature_path, warn = FALSE)[[1L]])
  if (!identical(old_signature, source_signature)) {
    stop(
      "Existing output uses a different input signature. Set a new ",
      "CRAFTGRN_METHOD10_OUTPUT path instead of mixing results."
    )
  }
} else {
  writeLines(source_signature, signature_path, useBytes = TRUE)
}

condition_manifest <- data.table::fread(
  condition_manifest_path,
  showProgress = FALSE
)
if (!all(c("condition_id", "path", "format") %in% names(condition_manifest))) {
  stop("Condition manifest is missing condition_id, path, or format.")
}
if (nrow(condition_manifest) != 17L) {
  stop("Expected 17 HPAFII condition inputs; found ", nrow(condition_manifest), ".")
}
condition_order <- as.character(condition_manifest$condition_id)

cell_line_specs <- data.table::data.table(
  cell_line = c("HPAFII", "AsPC-1", "PANC-1"),
  analysis_dir = file.path(
    project_root,
    c(
      "regulatory_topics_hpafii_condition_models_de_gene_filtered",
      "regulatory_topics_aspc1_condition_models_de_gene_filtered",
      "regulatory_topics_panc1_condition_models_de_gene_filtered"
    )
  )
)
cell_line_specs[, condition_manifest_path := file.path(
  analysis_dir,
  "condition_links",
  "condition_links_manifest.csv"
)]
if (any(!file.exists(cell_line_specs$condition_manifest_path))) {
  stop("A cell-line condition manifest is missing.")
}

report_theme <- function(base_size = 11) {
  craftgrn:::.m3_qc_theme(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        color = "#111111"
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", color = "#111111"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(face = "bold", size = base_size),
      plot.caption = ggplot2::element_text(face = "plain", size = base_size - 2),
      legend.position = "bottom"
    )
}

footprint_pdf <- file.path(output_root, "01_footprint_strength_diagnostic.pdf")
footprint_audit_path <- file.path(output_root, "footprint_strength_audit.csv")

detect_density_modes <- function(values, adjust, upper) {
  density_fit <- stats::density(
    log1p(values),
    n = 4096L,
    adjust = adjust,
    from = log1p(min(values)),
    to = log1p(upper),
    na.rm = TRUE
  )
  change <- diff(density_fit$y)
  peak_index <- which(
    utils::head(change, -1L) > 0 & utils::tail(change, -1L) <= 0
  ) + 1L
  peak_index <- peak_index[
    density_fit$y[peak_index] >= 0.03 * max(density_fit$y)
  ]
  list(
    curve = data.table::data.table(
      score_log = density_fit$x,
      score = expm1(density_fit$x),
      density = density_fit$y,
      bandwidth = paste0("adjust = ", adjust)
    ),
    modes = data.table::data.table(
      bandwidth = paste0("adjust = ", adjust),
      mode_score = expm1(density_fit$x[peak_index]),
      mode_density = density_fit$y[peak_index]
    )
  )
}

build_footprint_diagnostic <- function() {
  if (isTRUE(reuse) && file.exists(footprint_pdf) &&
      file.exists(footprint_audit_path)) {
    log_info("Reusing footprint-strength diagnostic: ", footprint_pdf)
    return(invisible(footprint_pdf))
  }
  log_info("Loading Module 1 footprint scores for the upstream strength diagnostic.")
  multiomic <- readRDS(multiomic_path)
  raw_score <- data.table::fread(raw_score_path, showProgress = FALSE)
  footprint_ids <- as.character(multiomic$features$fp$fp_id)
  row_index <- match(footprint_ids, raw_score$peak_ID)
  if (anyNA(row_index)) {
    stop("Raw footprint scores do not cover every Module 1 footprint.")
  }
  score_columns <- colnames(multiomic$matrices$fp_bound)
  absent_columns <- setdiff(score_columns, names(raw_score))
  if (length(absent_columns)) {
    stop("Raw footprint scores are missing sample(s): ", paste(absent_columns, collapse = ", "))
  }
  score_matrix <- as.matrix(raw_score[row_index, ..score_columns])
  bound_matrix <- multiomic$matrices$fp_bound
  if (!identical(dim(score_matrix), dim(bound_matrix))) {
    stop("Raw footprint score and bound matrices do not align.")
  }
  if (any(rowSums(bound_matrix > 0, na.rm = TRUE) < 1L)) {
    stop("The compact Module 1 object contains a footprint not bound in any condition.")
  }
  score_matrix[!bound_matrix] <- NA_real_
  maximum_bound_score <- apply(score_matrix, 1L, max, na.rm = TRUE)
  maximum_bound_score[!is.finite(maximum_bound_score)] <- NA_real_
  maximum_bound_score <- maximum_bound_score[is.finite(maximum_bound_score)]
  if (!length(maximum_bound_score)) stop("No finite bound footprint scores were found.")

  display_upper <- as.numeric(stats::quantile(
    maximum_bound_score,
    probs = 0.999,
    na.rm = TRUE,
    names = FALSE
  ))
  adjustments <- c(0.75, 1, 1.5, 2)
  density_parts <- lapply(
    adjustments,
    function(adjust) detect_density_modes(
      maximum_bound_score,
      adjust = adjust,
      upper = display_upper
    )
  )
  density_curves <- data.table::rbindlist(lapply(density_parts, `[[`, "curve"))
  density_modes <- data.table::rbindlist(lapply(density_parts, `[[`, "modes"))
  mode_counts <- density_modes[, .N, by = bandwidth]
  if (nrow(mode_counts) != length(adjustments) || any(mode_counts$N != 1L)) {
    stop("Footprint density did not have one stable material mode at every bandwidth.")
  }

  threshold_grid <- seq(2, 12, by = 0.1)
  retention <- data.table::data.table(
    threshold = threshold_grid,
    retained = vapply(
      threshold_grid,
      function(cutoff) sum(maximum_bound_score >= cutoff),
      numeric(1L)
    )
  )
  retention[, retained_percent := 100 * retained / length(maximum_bound_score)]
  selected_thresholds <- c(2, 3, 4, 5)
  selected_retention <- data.table::data.table(
    threshold = selected_thresholds,
    retained = vapply(
      selected_thresholds,
      function(cutoff) sum(maximum_bound_score >= cutoff),
      numeric(1L)
    )
  )
  selected_retention[, retained_percent :=
    100 * retained / length(maximum_bound_score)]

  quantile_probs <- c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)
  quantile_values <- as.numeric(stats::quantile(
    maximum_bound_score,
    probs = quantile_probs,
    na.rm = TRUE,
    names = FALSE
  ))
  audit <- data.table::rbindlist(list(
    data.table::data.table(
      metric = paste0("quantile_", quantile_probs),
      value = quantile_values,
      count = NA_real_
    ),
    selected_retention[, .(
      metric = paste0("retained_score_ge_", threshold),
      value = retained_percent,
      count = retained
    )],
    density_modes[, .(
      metric = paste0("density_mode_", gsub("[ =.]", "_", bandwidth)),
      value = mode_score,
      count = NA_real_
    )]
  ), use.names = TRUE, fill = TRUE)
  audit[, source_signature := source_signature]
  data.table::fwrite(audit, footprint_audit_path)

  score_breaks <- c(2, 3, 4, 5, 8, 12, 20)
  density_plot <- ggplot2::ggplot(
    density_curves[bandwidth == "adjust = 1"],
    ggplot2::aes(score_log, density)
  ) +
    ggplot2::geom_line(color = "#D95F02", linewidth = 1.2) +
    ggplot2::geom_vline(
      xintercept = log1p(score_floor),
      color = "#B2182B",
      linewidth = 1.1
    ) +
    ggplot2::geom_vline(
      xintercept = log1p(c(3, 4, 5)),
      color = "#777777",
      linewidth = 0.45,
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = log1p(score_breaks),
      labels = score_breaks,
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = "Footprint-score pattern",
      subtitle = NULL,
      x = "Highest footprint score in a bound condition (log scale)",
      y = "Relative frequency",
      caption = "Red line: current cutoff of 2. Gray lines: 3, 4, and 5."
    ) +
    report_theme(11) +
    ggplot2::theme(legend.position = "none")

  retention_labels <- selected_retention[, .(
    threshold,
    retained_percent,
    label = paste0(
      "Cutoff ", threshold, ": ",
      format(retained, big.mark = ",", scientific = FALSE),
      " (", sprintf("%.1f", retained_percent), "%)"
    )
  )]
  retention_plot <- ggplot2::ggplot(
    retention,
    ggplot2::aes(threshold, retained_percent)
  ) +
    ggplot2::geom_area(fill = "#67A9CF", alpha = 0.3) +
    ggplot2::geom_line(color = "#2166AC", linewidth = 1.1) +
    ggplot2::geom_point(
      data = retention_labels,
      ggplot2::aes(threshold, retained_percent),
      color = "#B2182B",
      size = 2.5
    ) +
    ggplot2::geom_text(
      data = retention_labels,
      ggplot2::aes(threshold, retained_percent, label = label),
      hjust = -0.03,
      vjust = -0.45,
      size = 3.2,
      family = "Helvetica",
      fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 2:12,
      limits = c(2, 12.7),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 105),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "What remains at each cutoff",
      subtitle = "Higher cutoffs remove many footprints.",
      x = "Footprint-score cutoff",
      y = "Footprints retained",
      caption = NULL
    ) +
    report_theme(11) +
    ggplot2::theme(legend.position = "none")

  summary_lines <- c(
    "Recommendation",
    "Keep the current cutoff of 2.",
    "",
    paste0("Footprints checked: ", format(length(maximum_bound_score), big.mark = ",")),
    paste0("Cutoff 3 keeps ", sprintf("%.0f", selected_retention[threshold == 3, retained_percent]), "%"),
    paste0("Cutoff 4 keeps ", sprintf("%.0f", selected_retention[threshold == 4, retained_percent]), "%"),
    paste0("Cutoff 5 keeps ", sprintf("%.0f", selected_retention[threshold == 5, retained_percent]), "%")
  )
  summary_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "label",
      x = 0,
      y = 1,
      label = paste(summary_lines, collapse = "\n"),
      hjust = 0,
      vjust = 1,
      size = 4.1,
      family = "Helvetica",
      fontface = "bold",
      linewidth = 0.45,
      label.padding = grid::unit(0.6, "lines"),
      fill = "#F7F4EC",
      color = "#111111"
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void()

  page <- density_plot / (retention_plot | summary_plot) +
    patchwork::plot_layout(heights = c(1.08, 1)) +
    patchwork::plot_annotation(
      title = "Should we use a higher footprint-score cutoff?",
      subtitle = "No clear split appears in the score pattern.",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = "Helvetica", face = "bold", size = 18
        ),
        plot.subtitle = ggplot2::element_text(
          family = "Helvetica", face = "bold", size = 12
        )
      )
    )
  grDevices::cairo_pdf(
    footprint_pdf,
    width = 15,
    height = 10,
    family = "Helvetica",
    bg = "white",
    onefile = TRUE
  )
  print(page)
  grDevices::dev.off()
  rm(
    multiomic, raw_score, score_matrix, bound_matrix,
    maximum_bound_score
  )
  invisible(gc())
  log_info("Wrote footprint-strength diagnostic: ", footprint_pdf)
  invisible(footprint_pdf)
}

eligibility_path <- file.path(output_root, "link_r_eligibility.parquet")
eligibility_summary_path <- file.path(output_root, "link_r_eligibility_summary.csv")
annotated_root <- file.path(output_root, "annotated_condition_links")
annotated_manifest_path <- file.path(
  annotated_root,
  "annotated_condition_links_manifest.csv"
)

build_link_eligibility <- function() {
  if (isTRUE(reuse) && file.exists(eligibility_path) &&
      file.exists(eligibility_summary_path)) {
    log_info("Reusing link R-eligibility table: ", eligibility_path)
    return(read_parquet_dt(eligibility_path))
  }
  log_info("Loading Module 2 correlation tables at R >= ", minimum_r, ".")
  tf_corr <- read_parquet_dt(
    tf_corr_path,
    c("tf", "target_gene", "best_r")
  )
  data.table::setnames(tf_corr, "best_r", "tf_target_best_r")
  tf_corr <- unique(tf_corr[
    is.finite(tf_target_best_r) & tf_target_best_r >= minimum_r
  ])
  if (anyDuplicated(tf_corr, by = c("tf", "target_gene"))) {
    stop("TF-target correlation keys are not unique.")
  }
  fp_corr <- read_parquet_dt(
    fp_corr_path,
    c("fp_id", "target_gene", "best_r")
  )
  data.table::setnames(fp_corr, "best_r", "fp_target_best_r")
  fp_corr <- unique(fp_corr[
    is.finite(fp_target_best_r) & fp_target_best_r >= minimum_r
  ])
  if (anyDuplicated(fp_corr, by = c("fp_id", "target_gene"))) {
    stop("FP-target correlation keys are not unique.")
  }
  data.table::setkey(tf_corr, tf, target_gene)
  data.table::setkey(fp_corr, fp_id, target_gene)

  stats_manifest <- data.table::fread(
    stats_manifest_path,
    showProgress = FALSE
  )
  link_manifest <- data.table::fread(
    link_manifest_path,
    showProgress = FALSE
  )
  if (!all(c("chunk_id", "path") %in% names(stats_manifest)) ||
      !all(c("chunk_id", "path") %in% names(link_manifest))) {
    stop("Module 1 or Module 2 chunk manifest is incomplete.")
  }
  chunk_plan <- merge(
    stats_manifest[, .(
      chunk_id = as.integer(chunk_id),
      stats_path = as.character(path)
    )],
    link_manifest[, .(
      chunk_id = as.integer(chunk_id),
      link_path = as.character(path)
    )],
    by = "chunk_id",
    all = TRUE,
    sort = TRUE
  )
  if (anyNA(chunk_plan$stats_path) || anyNA(chunk_plan$link_path)) {
    stop("Module 1 prediction-stat and Module 2 link chunks do not align.")
  }

  eligibility_parts <- vector("list", nrow(chunk_plan))
  chunk_audit <- vector("list", nrow(chunk_plan))
  for (i in seq_len(nrow(chunk_plan))) {
    row <- chunk_plan[i]
    log_info(
      "Filtering Module 1/2 link chunk ", i, "/", nrow(chunk_plan),
      " (chunk ", row$chunk_id, ")."
    )
    module1_stats <- read_parquet_dt(
      row$stats_path,
      c("tf", "fp_id", "best_r")
    )
    data.table::setnames(module1_stats, "best_r", "module1_best_r")
    module1_stats <- unique(module1_stats[
      is.finite(module1_best_r) & module1_best_r >= minimum_r
    ])
    if (anyDuplicated(module1_stats, by = c("tf", "fp_id"))) {
      stop("Module 1 prediction-stat keys are not unique in chunk ", row$chunk_id, ".")
    }
    links <- read_parquet_dt(
      row$link_path,
      c("link_id", "tf", "fp_id", "target_gene")
    )
    n_source <- nrow(links)
    if (anyDuplicated(links$link_id)) {
      stop("Module 2 link IDs are not unique in chunk ", row$chunk_id, ".")
    }
    links <- merge(
      links,
      module1_stats,
      by = c("tf", "fp_id"),
      all = FALSE,
      sort = FALSE,
      allow.cartesian = FALSE
    )
    links <- merge(
      links,
      tf_corr,
      by = c("tf", "target_gene"),
      all = FALSE,
      sort = FALSE,
      allow.cartesian = FALSE
    )
    links <- merge(
      links,
      fp_corr,
      by = c("fp_id", "target_gene"),
      all = FALSE,
      sort = FALSE,
      allow.cartesian = FALSE
    )
    links[, module2_min_best_r := pmin(
      tf_target_best_r,
      fp_target_best_r
    )]
    eligibility_parts[[i]] <- links[, .(
      link_id,
      module1_best_r,
      tf_target_best_r,
      fp_target_best_r,
      module2_min_best_r
    )]
    chunk_audit[[i]] <- data.table::data.table(
      chunk_id = row$chunk_id,
      source_links = n_source,
      eligible_links = nrow(links)
    )
    rm(module1_stats, links)
    invisible(gc())
  }
  eligibility <- data.table::rbindlist(
    eligibility_parts,
    use.names = TRUE
  )
  if (anyDuplicated(eligibility$link_id)) {
    stop("R-eligibility output contains duplicate link IDs.")
  }
  data.table::setkey(eligibility, link_id)
  write_parquet_dt(eligibility, eligibility_path)
  audit <- data.table::rbindlist(chunk_audit)
  audit[, source_signature := source_signature]
  data.table::fwrite(audit, eligibility_summary_path)
  log_info(
    "Retained ", format(nrow(eligibility), big.mark = ","),
    " Module 2 links at Module 1 R >= 0.5 and both Module 2 R >= 0.5."
  )
  eligibility
}

build_annotated_condition_links <- function(eligibility) {
  if (isTRUE(reuse) && file.exists(annotated_manifest_path)) {
    cached <- data.table::fread(
      annotated_manifest_path,
      showProgress = FALSE
    )
    if (nrow(cached) == nrow(condition_manifest) &&
        all(file.exists(cached$path)) &&
        all(cached$source_signature == source_signature)) {
      log_info("Reusing annotated condition-link subsets: ", annotated_root)
      return(cached)
    }
  }
  dir.create(annotated_root, recursive = TRUE, showWarnings = FALSE)
  data.table::setkey(eligibility, link_id)
  output_rows <- vector("list", nrow(condition_manifest))
  for (i in seq_len(nrow(condition_manifest))) {
    condition_id <- as.character(condition_manifest$condition_id[[i]])
    source_path <- as.character(condition_manifest$path[[i]])
    output_path <- file.path(
      annotated_root,
      paste0(condition_id, "_eligible_condition_links.parquet")
    )
    log_info(
      "Annotating condition links ", i, "/", nrow(condition_manifest),
      ": ", condition_id, "."
    )
    condition_links <- read_parquet_dt(source_path)
    source_rows <- nrow(condition_links)
    if (anyDuplicated(condition_links$link_id)) {
      stop("Condition link IDs are duplicated for ", condition_id, ".")
    }
    annotated <- eligibility[
      condition_links,
      on = "link_id",
      nomatch = 0L
    ]
    if (nrow(annotated) > source_rows) {
      stop("Eligibility join multiplied rows for ", condition_id, ".")
    }
    if (!all(annotated$condition_id == condition_id)) {
      stop("Condition identity changed while annotating ", condition_id, ".")
    }
    write_parquet_dt(annotated, output_path)
    output_rows[[i]] <- data.table::data.table(
      condition_id = condition_id,
      path = normalizePath(output_path, winslash = "/", mustWork = TRUE),
      source_rows = source_rows,
      eligible_rows = nrow(annotated),
      source_signature = source_signature
    )
    rm(condition_links, annotated)
    invisible(gc())
  }
  output <- data.table::rbindlist(output_rows)
  data.table::fwrite(output, annotated_manifest_path)
  output
}

load_annotated_condition_links <- function(manifest) {
  log_info("Loading R-annotated condition links for document construction.")
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    read_parquet_dt(manifest$path[[i]])
  })
  output <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  if (anyDuplicated(output, by = c("condition_id", "link_id"))) {
    stop("Annotated condition links contain duplicate condition/link IDs.")
  }
  output
}

combination_table <- data.table::CJ(
  module1_r = module1_cutoffs,
  module2_r = module2_cutoffs,
  sorted = TRUE
)
combination_table[, `:=`(
  combination_id = sprintf(
    "M1R%s_M2R%s",
    gsub("[.]", "p", module1_r),
    gsub("[.]", "p", module2_r)
  ),
  combination_label = sprintf(
    "M1 R >= %.1f / M2 R >= %.1f",
    module1_r,
    module2_r
  )
)]

document_summary_path <- file.path(output_root, "document_combination_summary.csv")
target_count_path <- file.path(output_root, "target_gene_counts.csv")
cell_line_target_count_path <- file.path(
  output_root,
  "cell_line_target_gene_counts.csv"
)
cell_line_target_signature_path <- file.path(
  output_root,
  "cell_line_target_gene_counts_signature.txt"
)

build_cell_line_target_counts <- function(eligibility = NULL) {
  manifests <- lapply(
    seq_len(nrow(cell_line_specs)),
    function(i) {
      value <- data.table::fread(
        cell_line_specs$condition_manifest_path[[i]],
        showProgress = FALSE
      )
      if (!all(c("condition_id", "path", "format") %in% names(value))) {
        stop("A cell-line condition manifest is incomplete.")
      }
      if (any(!file.exists(value$path))) {
        stop("A cell-line condition-link file is missing.")
      }
      value
    }
  )
  names(manifests) <- cell_line_specs$cell_line
  count_signature <- digest::digest(
    list(
      source_signature = source_signature,
      cutoffs = combination_table,
      manifests = lapply(
        seq_along(manifests),
        function(i) {
          list(
            manifest_md5 = unname(tools::md5sum(
              cell_line_specs$condition_manifest_path[[i]]
            )),
            condition_files = lapply(manifests[[i]]$path, file_identity)
          )
        }
      )
    ),
    algo = "xxhash64"
  )
  can_reuse <- isTRUE(reuse) &&
    file.exists(cell_line_target_count_path) &&
    file.exists(cell_line_target_signature_path) &&
    identical(
      trimws(readLines(cell_line_target_signature_path, warn = FALSE)[[1L]]),
      count_signature
    )
  if (can_reuse) {
    cached <- data.table::fread(
      cell_line_target_count_path,
      showProgress = FALSE
    )
    expected <- data.table::CJ(
      cell_line = cell_line_specs$cell_line,
      combination_id = combination_table$combination_id,
      unique = TRUE
    )
    observed <- unique(cached[, .(cell_line, combination_id)])
    if (!nrow(expected[!observed, on = c("cell_line", "combination_id")])) {
      log_info("Reusing three-cell-line target-gene counts.")
      return(cached)
    }
  }
  if (is.null(eligibility)) {
    eligibility <- read_parquet_dt(eligibility_path)
  }
  data.table::setkey(eligibility, link_id)

  hpafii <- data.table::fread(target_count_path, showProgress = FALSE)
  hpafii[, cell_line := "HPAFII"]
  output <- list(hpafii)
  for (cell_line in c("AsPC-1", "PANC-1")) {
    manifest <- manifests[[cell_line]]
    gene_unions <- stats::setNames(
      vector("list", nrow(combination_table)),
      combination_table$combination_id
    )
    rows <- vector("list", nrow(manifest) * nrow(combination_table))
    row_index <- 0L
    for (i in seq_len(nrow(manifest))) {
      condition_id <- as.character(manifest$condition_id[[i]])
      log_info(
        "Counting eligible genes for ", cell_line, " condition ",
        i, "/", nrow(manifest), "."
      )
      links <- read_parquet_dt(
        as.character(manifest$path[[i]]),
        c("link_id", "target_gene")
      )
      if (anyDuplicated(links$link_id)) {
        stop("Condition link IDs are duplicated for ", condition_id, ".")
      }
      annotated <- eligibility[links, on = "link_id", nomatch = 0L]
      for (j in seq_len(nrow(combination_table))) {
        combo <- combination_table[j]
        genes <- unique(annotated[
          module1_best_r >= combo$module1_r &
            tf_target_best_r >= combo$module2_r &
            fp_target_best_r >= combo$module2_r,
          as.character(target_gene)
        ])
        gene_unions[[combo$combination_id]] <- union(
          gene_unions[[combo$combination_id]],
          genes
        )
        row_index <- row_index + 1L
        rows[[row_index]] <- data.table::data.table(
          cell_line = cell_line,
          condition_id = condition_id,
          unique_target_genes = length(genes),
          module1_r = combo$module1_r,
          module2_r = combo$module2_r,
          combination_id = combo$combination_id,
          combination_label = combo$combination_label
        )
      }
      rm(links, annotated)
      invisible(gc())
    }
    cell_rows <- data.table::rbindlist(rows)
    union_rows <- combination_table[, .(
      cell_line = cell_line,
      condition_id = "All conditions (union)",
      unique_target_genes = lengths(gene_unions[combination_id]),
      module1_r,
      module2_r,
      combination_id,
      combination_label
    )]
    output[[length(output) + 1L]] <- data.table::rbindlist(
      list(cell_rows, union_rows),
      use.names = TRUE
    )
  }
  output <- data.table::rbindlist(output, use.names = TRUE, fill = TRUE)
  output[, source_signature := count_signature]
  for (current_cell_line in cell_line_specs$cell_line) {
    validate_target_count_monotonicity(
      output[cell_line == current_cell_line]
    )
  }
  data.table::setcolorder(
    output,
    c(
      "cell_line", "condition_id", "unique_target_genes",
      "module1_r", "module2_r", "combination_id", "combination_label",
      "source_signature"
    )
  )
  data.table::fwrite(output, cell_line_target_count_path)
  writeLines(count_signature, cell_line_target_signature_path, useBytes = TRUE)
  output
}

build_one_document_qc <- function(edges, row, specificity_lookup) {
  combo_dir <- file.path(output_root, "document_qc", row$combination_id)
  dir.create(combo_dir, recursive = TRUE, showWarnings = FALSE)
  qc_path <- file.path(combo_dir, "document_term_qc.pdf")
  summary_path <- file.path(combo_dir, "combination_summary.csv")
  target_path <- file.path(combo_dir, "target_gene_counts.csv")
  membership_path <- file.path(combo_dir, "term_membership.parquet")
  complete <- all(file.exists(c(
    qc_path,
    summary_path,
    target_path,
    membership_path
  )))
  if (isTRUE(reuse) && complete) {
    log_info("Reusing document QC for ", row$combination_label, ".")
    return(list(
      summary = data.table::fread(summary_path, showProgress = FALSE),
      target = data.table::fread(target_path, showProgress = FALSE),
      qc = qc_path,
      membership = membership_path
    ))
  }

  filtered <- data.table::copy(edges[
    module1_best_r >= row$module1_r &
      tf_target_best_r >= row$module2_r &
      fp_target_best_r >= row$module2_r
  ])
  if (!nrow(filtered)) {
    stop("No condition links remain for ", row$combination_label, ".")
  }
  if (data.table::uniqueN(filtered$condition_id) != length(condition_order)) {
    stop("Not all 17 conditions remain for ", row$combination_label, ".")
  }
  log_info(
    "Building method-10 documents for ", row$combination_label,
    " from ", format(nrow(filtered), big.mark = ","), " condition links."
  )
  doc_term <- craftgrn:::build_doc_term_condition_union(
    edges_condition = filtered,
    count_method = "log",
    count_scale = 50,
    prefix_terms = TRUE,
    threshold_gene_expr = 10,
    threshold_fp_score = 0,
    threshold_tf_expr = -Inf,
    include_tf_terms = FALSE,
    require_tf_expr = TRUE,
    fp_term_mode = "unique",
    condition_peak_weighting = "tf_expression",
    condition_gene_specificity = specificity_lookup,
    balance_mode = "min",
    check_repeated_values = FALSE
  )
  if (!nrow(doc_term)) {
    stop("Document construction returned no terms for ", row$combination_label, ".")
  }

  token_cap <- craftgrn:::.cap_warplda_token_counts(
    doc_term$pseudo_count_log
  )
  doc_term[, pseudo_count_log := token_cap$counts]
  doc_term <- craftgrn:::.topic_finalize_condition_tf_counts(
    doc_term = doc_term,
    count_column = "pseudo_count_log",
    final_peak_gene_token_ratio = 0.5,
    condition_term_idf = FALSE,
    condition_term_idf_floor = 0.1
  )
  doc_term[, pseudo_count := pseudo_count_log]
  token_audit <- attr(doc_term, "final_condition_token_audit")
  if (is.null(token_audit) || !nrow(token_audit)) {
    stop("Final Gene/Peak token audit is missing for ", row$combination_label, ".")
  }
  token_totals <- token_audit[, .(
    final_tokens = sum(final_tokens)
  ), by = modality]
  peak_tokens <- token_totals[modality == "peak", final_tokens]
  gene_tokens <- token_totals[modality == "gene", final_tokens]
  observed_ratio <- peak_tokens / gene_tokens
  if (!length(observed_ratio) || !is.finite(observed_ratio) ||
      abs(observed_ratio - 0.5) > 1e-4) {
    stop(
      "Final Peak/Gene token ratio is invalid for ",
      row$combination_label, ": ", observed_ratio
    )
  }

  title <- paste0(
    "HPAFII method 10 document-term QC - Module 1 R >= ",
    row$module1_r,
    "; Module 2 TF-target and FP-target R >= ",
    row$module2_r
  )
  craftgrn:::.write_module3_document_term_qc(
    doc_term = doc_term,
    output_dir = combo_dir,
    count_column = "pseudo_count_log",
    title = title,
    document_unit = "tf",
    verbose = TRUE
  )

  gene_rows <- doc_term[startsWith(term_id, "GENE:")]
  target_counts <- gene_rows[, .(
    unique_target_genes = data.table::uniqueN(sub("^GENE:", "", term_id))
  ), by = .(
    condition_id = sub("::[^:]+$", "", doc_id)
  )]
  target_counts <- merge(
    data.table::data.table(condition_id = condition_order),
    target_counts,
    by = "condition_id",
    all.x = TRUE,
    sort = FALSE
  )
  target_counts[is.na(unique_target_genes), unique_target_genes := 0L]
  all_union <- data.table::data.table(
    condition_id = "All conditions (union)",
    unique_target_genes = data.table::uniqueN(
      sub("^GENE:", "", gene_rows$term_id)
    )
  )
  target_counts <- data.table::rbindlist(list(target_counts, all_union))
  target_counts[, `:=`(
    module1_r = row$module1_r,
    module2_r = row$module2_r,
    combination_id = row$combination_id,
    combination_label = row$combination_label
  )]

  term_type <- data.table::fifelse(
    startsWith(doc_term$term_id, "GENE:"),
    "Gene",
    "Peak"
  )
  summary <- data.table::data.table(
    module1_r = row$module1_r,
    module2_r = row$module2_r,
    combination_id = row$combination_id,
    combination_label = row$combination_label,
    condition_link_rows = nrow(filtered),
    documents = data.table::uniqueN(doc_term$doc_id),
    conditions = data.table::uniqueN(sub("::[^:]+$", "", doc_term$doc_id)),
    document_term_rows = nrow(doc_term),
    unique_gene_terms = data.table::uniqueN(doc_term$term_id[term_type == "Gene"]),
    unique_peak_terms = data.table::uniqueN(doc_term$term_id[term_type == "Peak"]),
    raw_model_tokens = token_cap$raw_tokens,
    capped_model_tokens = token_cap$tokens,
    global_token_scale_factor = token_cap$scale_factor,
    final_gene_tokens = gene_tokens,
    final_peak_tokens = peak_tokens,
    observed_peak_gene_ratio = observed_ratio,
    source_signature = source_signature
  )
  data.table::fwrite(summary, summary_path)
  data.table::fwrite(target_counts, target_path)
  membership <- unique(doc_term[, .(
    doc_id = as.character(doc_id),
    term_id = as.character(term_id)
  )])
  write_parquet_dt(membership, membership_path)
  rm(filtered, doc_term, gene_rows, membership, token_audit)
  invisible(gc())
  list(
    summary = summary,
    target = target_counts,
    qc = qc_path,
    membership = membership_path
  )
}

validate_nested_membership <- function(results) {
  relationship <- data.table::data.table(
    strict = c(
      "M1R0p5_M2R0p7",
      "M1R0p7_M2R0p5",
      "M1R0p7_M2R0p7",
      "M1R0p7_M2R0p7"
    ),
    loose = c(
      "M1R0p5_M2R0p5",
      "M1R0p5_M2R0p5",
      "M1R0p5_M2R0p7",
      "M1R0p7_M2R0p5"
    )
  )
  membership_cache <- new.env(parent = emptyenv())
  fetch <- function(combination_id) {
    if (!exists(combination_id, envir = membership_cache, inherits = FALSE)) {
      path <- results[[combination_id]]$membership
      value <- read_parquet_dt(path, c("doc_id", "term_id"))
      data.table::setkey(value, doc_id, term_id)
      assign(combination_id, value, envir = membership_cache)
    }
    get(combination_id, envir = membership_cache, inherits = FALSE)
  }
  audit <- lapply(seq_len(nrow(relationship)), function(i) {
    strict <- fetch(relationship$strict[[i]])
    loose <- fetch(relationship$loose[[i]])
    unexpected <- strict[!loose, on = c("doc_id", "term_id")]
    data.table::data.table(
      strict = relationship$strict[[i]],
      loose = relationship$loose[[i]],
      strict_rows = nrow(strict),
      loose_rows = nrow(loose),
      unexpected_rows = nrow(unexpected),
      passed = nrow(unexpected) == 0L
    )
  })
  audit <- data.table::rbindlist(audit)
  if (any(!audit$passed)) {
    stop("A stricter document-term result is not nested in its looser result.")
  }
  audit[, source_signature := source_signature]
  data.table::fwrite(
    audit,
    file.path(output_root, "document_nesting_audit.csv")
  )
  invisible(audit)
}

validate_target_count_monotonicity <- function(target_counts) {
  wide <- data.table::dcast(
    target_counts,
    condition_id ~ combination_id,
    value.var = "unique_target_genes"
  )
  required <- combination_table$combination_id
  if (!all(required %in% names(wide))) {
    stop("Target-gene comparison is missing a cutoff combination.")
  }
  valid <- with(wide,
    M1R0p5_M2R0p7 <= M1R0p5_M2R0p5 &
      M1R0p7_M2R0p5 <= M1R0p5_M2R0p5 &
      M1R0p7_M2R0p7 <= M1R0p5_M2R0p7 &
      M1R0p7_M2R0p7 <= M1R0p7_M2R0p5
  )
  if (any(!valid)) {
    stop("Target-gene counts increased under a stricter cutoff.")
  }
  invisible(TRUE)
}

build_documents <- function() {
  eligibility <- build_link_eligibility()
  annotated_manifest <- build_annotated_condition_links(eligibility)
  rm(eligibility)
  invisible(gc())
  annotated <- load_annotated_condition_links(annotated_manifest)
  specificity_lookup <- craftgrn:::.module3_condition_gene_specificity_lookup(
    expression_file = expression_path,
    expression_min = 10,
    temperature = 0.5,
    uniform_floor = 0.1
  )
  results <- vector("list", nrow(combination_table))
  names(results) <- combination_table$combination_id
  for (i in seq_len(nrow(combination_table))) {
    row <- combination_table[i]
    results[[row$combination_id]] <- build_one_document_qc(
      edges = annotated,
      row = row,
      specificity_lookup = specificity_lookup
    )
  }
  document_summary <- data.table::rbindlist(
    lapply(results, `[[`, "summary"),
    use.names = TRUE,
    fill = TRUE
  )
  target_counts <- data.table::rbindlist(
    lapply(results, `[[`, "target"),
    use.names = TRUE,
    fill = TRUE
  )
  if (any(document_summary$conditions != length(condition_order))) {
    stop("A document combination does not cover all 17 conditions.")
  }
  validate_nested_membership(results)
  validate_target_count_monotonicity(target_counts)
  data.table::fwrite(document_summary, document_summary_path)
  data.table::fwrite(target_counts, target_count_path)
  rm(annotated, specificity_lookup)
  invisible(gc())
  results
}

comparison_pdf <- file.path(output_root, "02_target_gene_comparison.pdf")
final_pdf <- file.path(
  output_root,
  "15_Method10_K30_RCutoff_and_FPStrength_DocumentTermQC.pdf"
)

compact_qc_summary_path <- file.path(
  output_root,
  "compact_document_qc_summary.csv"
)
compact_qc_paths <- file.path(
  output_root,
  paste0("compact_document_qc_", combination_table$combination_id, ".pdf")
)

simple_condition_label <- function(value) {
  value <- sub("^(HPAFII|Aspc1|Panc1)_", "", value)
  value <- sub("_Ctrl$", "", value)
  value <- gsub("Gln[.]Arg", "Gln + Arg", value)
  value <- gsub("Met[.]Cys", "Met + Cys", value)
  gsub("_", " ", value, fixed = TRUE)
}

cutoff_labels <- stats::setNames(
  sprintf(
    "%.1f / %.1f",
    combination_table$module1_r,
    combination_table$module2_r
  ),
  combination_table$combination_id
)

build_compact_qc_summary <- function() {
  if (isTRUE(reuse) && file.exists(compact_qc_summary_path)) {
    cached <- data.table::fread(
      compact_qc_summary_path,
      showProgress = FALSE
    )
    if (nrow(cached) && all(cached$source_signature == source_signature)) {
      log_info("Reusing compact document summaries.")
      return(cached)
    }
  }
  output <- vector("list", nrow(combination_table))
  for (i in seq_len(nrow(combination_table))) {
    combo <- combination_table[i]
    membership_path <- file.path(
      output_root,
      "document_qc",
      combo$combination_id,
      "term_membership.parquet"
    )
    membership <- read_parquet_dt(
      membership_path,
      c("doc_id", "term_id")
    )
    membership[, `:=`(
      condition_id = sub("::[^:]+$", "", doc_id),
      term_type = data.table::fifelse(
        startsWith(term_id, "GENE:"),
        "Gene",
        "Peak"
      )
    )]
    per_document <- membership[, .(
      terms_in_document = .N
    ), by = .(condition_id, doc_id, term_type)]
    typical <- per_document[, .(
      median_terms = as.numeric(stats::median(terms_in_document)),
      lower_quartile = as.numeric(stats::quantile(
        terms_in_document,
        probs = 0.25,
        names = FALSE
      )),
      upper_quartile = as.numeric(stats::quantile(
        terms_in_document,
        probs = 0.75,
        names = FALSE
      )),
      tf_documents = data.table::uniqueN(doc_id)
    ), by = .(condition_id, term_type)]
    totals <- membership[, .(
      unique_terms = data.table::uniqueN(term_id),
      document_term_rows = .N
    ), by = .(condition_id, term_type)]
    result <- merge(
      typical,
      totals,
      by = c("condition_id", "term_type"),
      all = TRUE,
      sort = FALSE
    )
    result[, `:=`(
      module1_r = combo$module1_r,
      module2_r = combo$module2_r,
      combination_id = combo$combination_id,
      combination_label = combo$combination_label,
      source_signature = source_signature
    )]
    output[[i]] <- result
    rm(membership, per_document, typical, totals, result)
    invisible(gc())
  }
  output <- data.table::rbindlist(output, use.names = TRUE)
  data.table::fwrite(output, compact_qc_summary_path)
  output
}

build_comparison_page <- function(target_counts) {
  plot_data <- target_counts[
    cell_line == "HPAFII" & condition_id != "All conditions (union)"
  ]
  plot_data[, condition_display := simple_condition_label(condition_id)]
  plot_data[, condition_display := factor(
    condition_display,
    levels = rev(simple_condition_label(condition_order))
  )]
  plot_data[, cutoff := factor(
    cutoff_labels[combination_id],
    levels = cutoff_labels
  )]
  colors <- stats::setNames(
    c("#2166AC", "#67A9CF", "#B2182B", "#EF8A62"),
    cutoff_labels
  )
  per_condition_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      condition_display,
      unique_target_genes,
      fill = cutoff
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.82, preserve = "single"),
      width = 0.72,
      color = "#222222",
      linewidth = 0.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      title = "HPAFII conditions",
      subtitle = "Unique target genes kept in each condition",
      x = NULL,
      y = "Unique target genes",
      fill = "Module 1 / Module 2"
    ) +
    report_theme(10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 8)
    )

  union_data <- target_counts[condition_id == "All conditions (union)"]
  union_data[, `:=`(
    cell_line = factor(cell_line, levels = cell_line_specs$cell_line),
    cutoff = factor(cutoff_labels[combination_id], levels = rev(cutoff_labels))
  )]
  union_plot <- ggplot2::ggplot(
    union_data,
    ggplot2::aes(cutoff, unique_target_genes, fill = cutoff)
  ) +
    ggplot2::geom_col(width = 0.68, color = "#222222", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = format(unique_target_genes, big.mark = ",", scientific = FALSE)
      ),
      hjust = -0.08,
      family = "Helvetica",
      fontface = "bold",
      size = 3.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_grid(cell_line ~ .) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::labs(
      title = "All conditions combined",
      subtitle = "Rows show Module 1 / Module 2 cutoffs",
      x = NULL,
      y = "Unique target genes"
    ) +
    report_theme(10) +
    ggplot2::theme(legend.position = "none")

  page <- (per_condition_plot | union_plot) +
    patchwork::plot_layout(widths = c(1.45, 1))
  page <- page + patchwork::plot_annotation(
    title = "How stricter cutoffs reduce target-gene coverage",
    subtitle = "A higher cutoff keeps only stronger links and fewer genes.",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Helvetica", face = "bold", size = 18
      ),
      plot.subtitle = ggplot2::element_text(
        family = "Helvetica", face = "bold", size = 11
      )
    )
  )
  grDevices::cairo_pdf(
    comparison_pdf,
    width = 15,
    height = 10,
    family = "Helvetica",
    bg = "white",
    onefile = TRUE
  )
  print(page)
  grDevices::dev.off()
  invisible(comparison_pdf)
}

build_compact_qc_pages <- function(document_summary, compact_summary) {
  term_colors <- c(Gene = "#4C78A8", Peak = "#F58518")
  for (i in seq_len(nrow(combination_table))) {
    combo <- combination_table[i]
    plot_data <- data.table::copy(
      compact_summary[combination_id == combo$combination_id]
    )
    plot_data[, `:=`(
      condition_display = factor(
        simple_condition_label(condition_id),
        levels = rev(simple_condition_label(condition_order))
      ),
      term_type = factor(term_type, levels = c("Gene", "Peak"))
    )]
    typical_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(condition_display, median_terms, fill = term_type)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge2(width = 0.8, preserve = "single"),
        width = 0.7,
        color = "#222222",
        linewidth = 0.2
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = term_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(big.mark = ","),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::labs(
        title = "Terms in a typical document",
        subtitle = "Middle value across transcription factors",
        x = NULL,
        y = "Terms",
        fill = NULL
      ) +
      report_theme(10) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "bottom"
      )
    unique_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(condition_display, unique_terms, fill = term_type)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge2(width = 0.8, preserve = "single"),
        width = 0.7,
        color = "#222222",
        linewidth = 0.2
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = term_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_number(big.mark = ","),
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::labs(
        title = "Unique terms in each condition",
        subtitle = "Each term is counted once",
        x = NULL,
        y = "Unique terms",
        fill = NULL
      ) +
      report_theme(10) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "bottom"
      )

    summary_row <- document_summary[
      combination_id == combo$combination_id
    ][1L]
    subtitle <- paste0(
      format(summary_row$documents, big.mark = ","), " documents | ",
      format(summary_row$unique_gene_terms, big.mark = ","), " genes | ",
      format(summary_row$unique_peak_terms, big.mark = ","), " peaks"
    )
    page <- (typical_plot | unique_plot) +
      patchwork::plot_layout(widths = c(1, 1)) +
      patchwork::plot_annotation(
        title = paste0(
          "HPAFII documents: Module 1 cutoff ", combo$module1_r,
          ", Module 2 cutoff ", combo$module2_r
        ),
        subtitle = subtitle,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            family = "Helvetica", face = "bold", size = 18
          ),
          plot.subtitle = ggplot2::element_text(
            family = "Helvetica", face = "bold", size = 11
          )
        )
      )
    grDevices::cairo_pdf(
      compact_qc_paths[[i]],
      width = 15,
      height = 10,
      family = "Helvetica",
      bg = "white",
      onefile = TRUE
    )
    print(page)
    grDevices::dev.off()
  }
  invisible(compact_qc_paths)
}

assemble_report <- function() {
  required <- c(
    footprint_pdf,
    document_summary_path,
    target_count_path
  )
  if (any(!file.exists(required))) {
    stop("Report assembly inputs are missing: ", paste(required[!file.exists(required)], collapse = ", "))
  }
  document_summary <- data.table::fread(
    document_summary_path,
    showProgress = FALSE
  )
  target_counts <- build_cell_line_target_counts()
  build_comparison_page(target_counts)
  compact_summary <- build_compact_qc_summary()
  build_compact_qc_pages(document_summary, compact_summary)
  qpdf <- Sys.which("qpdf")
  if (!nzchar(qpdf)) stop("qpdf is required to assemble the final report.")
  source_pdfs <- c(footprint_pdf, comparison_pdf, compact_qc_paths)
  qpdf_args <- c(
    "--empty",
    "--pages",
    unlist(lapply(source_pdfs, function(path) c(path, "1-z"))),
    "--",
    final_pdf
  )
  status <- system2(qpdf, qpdf_args)
  if (!identical(status, 0L) || !file.exists(final_pdf)) {
    stop("qpdf failed to assemble the final report.")
  }
  page_count <- suppressWarnings(as.integer(system2(
    qpdf,
    c("--show-npages", final_pdf),
    stdout = TRUE
  )[[1L]]))
  if (!identical(page_count, 6L)) {
    stop("Final report must contain 6 pages; found ", page_count, ".")
  }
  check_status <- system2(qpdf, c("--check", final_pdf))
  if (!identical(check_status, 0L)) stop("qpdf integrity validation failed.")
  report_manifest <- data.table::data.table(
    file = basename(final_pdf),
    path = normalizePath(final_pdf, winslash = "/", mustWork = TRUE),
    pages = page_count,
    bytes = as.numeric(file.info(final_pdf)$size),
    md5 = unname(tools::md5sum(final_pdf)),
    source_signature = source_signature
  )
  data.table::fwrite(
    report_manifest,
    file.path(output_root, "final_report_manifest.csv")
  )
  log_info("Wrote final 6-page report: ", final_pdf)
  invisible(final_pdf)
}

if (stage %in% c("all", "footprint")) {
  build_footprint_diagnostic()
}
if (stage %in% c("all", "eligibility")) {
  eligibility <- build_link_eligibility()
  build_annotated_condition_links(eligibility)
  rm(eligibility)
  invisible(gc())
}
if (stage %in% c("all", "documents")) {
  invisible(build_documents())
}
if (stage %in% c("all", "report")) {
  invisible(assemble_report())
}

invisible(TRUE)
