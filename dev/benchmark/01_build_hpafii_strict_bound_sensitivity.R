#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(args[[1L]]) else "all"
allowed_stages <- c("all", "plots", "prepare", "module1")
if (!stage %in% allowed_stages) {
  stop("Stage must be one of: ", paste(allowed_stages, collapse = ", "))
}

required_packages <- c(
  "craftgrn", "data.table", "ggplot2", "scales"
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

log_info <- function(...) {
  craftgrn:::.log_inform(paste0(...))
}

parse_workers <- function() {
  workers <- suppressWarnings(as.integer(Sys.getenv(
    "CRAFTGRN_STRICT_BOUND_WORKERS",
    unset = "8"
  )))
  if (!is.finite(workers) || workers < 1L) stop("Invalid worker count.")
  min(workers, 12L)
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
output_root <- file.path(
  analysis_root,
  "topic_runs",
  "review",
  "method10_module1_module2_r_cutoff_sensitivity"
)
work_root <- file.path(output_root, "hpa_strict_bound_sensitivity")
module1_root <- file.path(work_root, "module1_r0p5_p1e10")
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

project_config_path <- file.path(project_root, "project.yaml")
metadata_path <- file.path(project_root, "data", "sample_metadata_expanded.csv")
atac_path <- file.path(project_root, "data", "atac_master.csv")
cache_root <- file.path(project_root, "cache")
adaptive_path <- file.path(
  work_root,
  "adaptive_bound_cutoffs.csv"
)
front_pdf <- file.path(work_root, "footprint_bound_sensitivity_front.pdf")
histogram_path <- file.path(work_root, "footprint_score_histogram.csv")
strict_cutoff_path <- file.path(work_root, "strict_sample_cutoffs.csv")
multiomic_path <- file.path(work_root, "multiomic_p1e10.rds")
module1_marker <- file.path(module1_root, "strict_bound_module1_complete.txt")

required_inputs <- c(
  project_config_path,
  metadata_path,
  atac_path,
  adaptive_path,
  file.path(cache_root, "fp_scores_JASPAR2024.csv"),
  file.path(cache_root, "fp_bounds_JASPAR2024.csv"),
  file.path(cache_root, "fp_sites_JASPAR2024.csv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "))
}

p_grid <- data.table::data.table(
  p_id = c("p1e10", "p1e12", "p1e15", "p1e20"),
  p_label = c("p = 1e-10", "p = 1e-12", "p = 1e-15", "p = 1e-20"),
  p_value = c(1e-10, 1e-12, 1e-15, 1e-20)
)
p_colors <- stats::setNames(
  c("#2C7BB6", "#00A6CA", "#F4A261", "#B2182B"),
  p_grid$p_label
)

load_metadata <- function() {
  metadata <- data.table::fread(metadata_path, showProgress = FALSE)
  metadata <- metadata[cell == "HPAFII"]
  if (nrow(metadata) != 33L) {
    stop("Expected 33 HPAFII samples; found ", nrow(metadata), ".")
  }
  if (anyDuplicated(metadata$id)) stop("HPAFII sample ids are duplicated.")
  for (column in c("id_fp", "id_atac")) {
    if (column %in% names(metadata) && any(metadata[[column]] != metadata$id)) {
      stop("This benchmark requires HPAFII ", column, " to match metadata id.")
    }
  }
  metadata
}

build_strict_cutoffs <- function(metadata) {
  adaptive <- data.table::fread(adaptive_path, showProgress = FALSE)
  required <- c("sample_id", "null_mode", "null_sd", "cutoff_p1e3")
  missing <- setdiff(required, names(adaptive))
  if (length(missing)) {
    stop("Adaptive cutoff table is missing: ", paste(missing, collapse = ", "))
  }
  if (!setequal(adaptive$sample_id, metadata$id)) {
    stop("Adaptive cutoff samples do not match the 33 HPAFII samples.")
  }
  adaptive <- adaptive[match(metadata$id, sample_id)]
  adaptive[, join_key__ := 1L]
  p_grid_join <- data.table::copy(p_grid)
  p_grid_join[, join_key__ := 1L]
  strict <- merge(
    adaptive[, .(sample_id, null_mode, null_sd, cutoff_p1e3, join_key__)],
    p_grid_join,
    by = "join_key__",
    allow.cartesian = TRUE
  )
  strict[, join_key__ := NULL]
  strict[, cutoff := null_mode + stats::qnorm(
    p_value,
    lower.tail = FALSE
  ) * null_sd]
  if (any(!is.finite(strict$cutoff)) || any(strict$cutoff <= 0)) {
    stop("Strict sample cutoffs are not positive finite values.")
  }
  strict[, p_order__ := match(p_id, p_grid$p_id)]
  data.table::setorder(strict, sample_id, p_order__)
  strict[, p_order__ := NULL]
  data.table::fwrite(strict, strict_cutoff_path)
  strict
}

report_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        face = "bold",
        color = "#111111"
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", color = "#111111"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 5),
      plot.subtitle = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.caption = ggplot2::element_text(face = "plain", size = base_size - 2),
      legend.position = "bottom"
    )
}

simple_condition_label <- function(value) {
  value <- sub("^HPAFII_", "", value)
  value <- gsub("Gln[.]Arg", "Gln + Arg", value)
  value <- gsub("Met[.]Cys", "Met + Cys", value)
  gsub("_", " ", value, fixed = TRUE)
}

histogram_one <- function(scores, bound, breaks, sample_id, condition_label) {
  score_log <- log1p(pmax(0, as.numeric(scores)))
  bound <- as.logical(bound)
  bound[is.na(bound)] <- FALSE
  status_values <- c("Bound", "Unbound")
  parts <- lapply(status_values, function(status) {
    keep <- if (identical(status, "Bound")) bound else !bound
    counts <- graphics::hist(
      score_log[keep],
      breaks = breaks,
      plot = FALSE,
      include.lowest = TRUE,
      right = TRUE
    )$counts
    data.table::data.table(
      sample_id = sample_id,
      condition_label = condition_label,
      status = status,
      score_log = (breaks[-1L] + breaks[-length(breaks)]) / 2,
      count = as.numeric(counts)
    )
  })
  data.table::rbindlist(parts)
}

build_front_pages <- function() {
  metadata <- load_metadata()
  strict <- build_strict_cutoffs(metadata)
  if (file.exists(histogram_path)) {
    histogram <- data.table::fread(histogram_path, showProgress = FALSE)
  } else {
    log_info("Loading full HPAFII footprint scores for the distribution pages.")
    score <- data.table::fread(
      file.path(cache_root, "fp_scores_JASPAR2024.csv"),
      select = c("peak_ID", metadata$id),
      showProgress = FALSE
    )
    bound <- data.table::fread(
      file.path(cache_root, "fp_bounds_JASPAR2024.csv"),
      select = c("peak_ID", metadata$id),
      showProgress = FALSE
    )
    if (!identical(score$peak_ID, bound$peak_ID)) {
      stop("Raw footprint score and bound rows do not align.")
    }
    score_max <- max(
      as.matrix(score[, metadata$id, with = FALSE]),
      na.rm = TRUE
    )
    breaks <- seq(0, log1p(score_max) + 1e-8, length.out = 181L)
    histogram <- vector("list", nrow(metadata))
    for (i in seq_len(nrow(metadata))) {
      sample_id <- metadata$id[[i]]
      histogram[[i]] <- histogram_one(
        scores = score[[sample_id]],
        bound = bound[[sample_id]],
        breaks = breaks,
        sample_id = sample_id,
        condition_label = metadata$name[[i]]
      )
    }
    histogram <- data.table::rbindlist(histogram)
    histogram[, condition_short := simple_condition_label(condition_label)]
    histogram[, total_sample := sum(count), by = sample_id]
    histogram[, share := count / total_sample]
    data.table::fwrite(histogram, histogram_path)
    rm(score, bound)
    invisible(gc())
  }

  combined <- histogram[, .(count = sum(count)), by = .(status, score_log)]
  combined[, share := count / sum(count)]
  score_breaks <- c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256)
  status_colors <- c(Bound = "#D95F02", Unbound = "#4C78A8")
  combined_plot <- ggplot2::ggplot(
    combined,
    ggplot2::aes(score_log, share, color = status, fill = status)
  ) +
    ggplot2::geom_area(alpha = 0.2, position = "identity") +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::scale_fill_manual(values = status_colors) +
    ggplot2::scale_x_continuous(
      breaks = log1p(score_breaks),
      labels = score_breaks,
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = "HPAFII footprint-score distribution",
      subtitle = "All 33 samples and all footprint sites",
      x = "Footprint score (log scale)",
      y = "Share of all sites",
      color = NULL,
      fill = NULL,
      caption = "Bound and unbound are the original fp-tools calls."
    ) +
    report_theme(13)

  metadata[, facet_label := paste(id, simple_condition_label(name), sep = " | ")]
  facet_levels <- rev(metadata$facet_label)
  histogram[, facet_label := factor(
    paste(sample_id, simple_condition_label(condition_label), sep = " | "),
    levels = facet_levels
  )]
  strict_plot_data <- merge(
    strict,
    metadata[, .(sample_id = id, condition_label = name)],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  strict_plot_data[, facet_label := factor(
    paste(sample_id, simple_condition_label(condition_label), sep = " | "),
    levels = facet_levels
  )]
  strict_plot_data[, cutoff_log := log1p(cutoff)]
  strict_plot_data[, p_label := factor(p_label, levels = p_grid$p_label)]
  sample_plot <- ggplot2::ggplot(
    histogram,
    ggplot2::aes(score_log, share, color = status)
  ) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.9) +
    ggplot2::geom_vline(
      data = strict_plot_data,
      ggplot2::aes(xintercept = cutoff_log, linetype = p_label),
      color = "#111111",
      linewidth = 0.34
    ) +
    ggplot2::geom_vline(
      xintercept = log1p(2),
      color = "#B2182B",
      linewidth = 0.34,
      linetype = "solid"
    ) +
    ggplot2::facet_wrap(~facet_label, ncol = 6, scales = "free_y") +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "22", "42", "13"),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = log1p(c(0, 1, 2, 4, 16, 64, 256)),
      labels = c(0, 1, 2, 4, 16, 64, 256),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::labs(
      title = "Stricter learned cutoffs in each sample",
      subtitle = "Black lines show p = 1e-10 to 1e-20. Red line shows score 2.",
      x = "Footprint score (log scale)",
      y = "Share of sites in the sample",
      color = NULL,
      linetype = "Learned cutoff"
    ) +
    report_theme(8) +
    ggplot2::theme(
      panel.spacing = grid::unit(0.35, "lines"),
      axis.text = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  multiplier <- merge(
    strict[, .(cutoff_p1e3 = unique(cutoff_p1e3)), by = sample_id],
    metadata[, .(sample_id = id, condition_label = name)],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  multiplier <- multiplier[, .(
    multiplier = c("1.5x", "2x", "3x"),
    cutoff = cutoff_p1e3 * c(1.5, 2, 3)
  ), by = .(sample_id, condition_label)]
  multiplier[, `:=`(
    facet_label = factor(
      paste(sample_id, simple_condition_label(condition_label), sep = " | "),
      levels = facet_levels
    ),
    cutoff_log = log1p(cutoff),
    multiplier = factor(multiplier, levels = c("1.5x", "2x", "3x"))
  )]
  multiplier_plot <- ggplot2::ggplot(
    histogram,
    ggplot2::aes(score_log, share, color = status)
  ) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.9) +
    ggplot2::geom_vline(
      data = multiplier,
      ggplot2::aes(xintercept = cutoff_log, linetype = multiplier),
      color = "#111111",
      linewidth = 0.34
    ) +
    ggplot2::geom_vline(
      xintercept = log1p(2),
      color = "#B2182B",
      linewidth = 0.34
    ) +
    ggplot2::facet_wrap(~facet_label, ncol = 6, scales = "free_y") +
    ggplot2::scale_color_manual(values = status_colors) +
    ggplot2::scale_linetype_manual(values = c("solid", "22", "42")) +
    ggplot2::scale_x_continuous(
      breaks = log1p(c(0, 1, 2, 4, 16, 64, 256)),
      labels = c(0, 1, 2, 4, 16, 64, 256),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::labs(
      title = "Multiples of the original learned cutoff",
      subtitle = "Black lines show 1.5x, 2x, and 3x. Red line shows score 2.",
      x = "Footprint score (log scale)",
      y = "Share of sites in the sample",
      color = NULL,
      linetype = "Multiplier"
    ) +
    report_theme(8) +
    ggplot2::theme(
      panel.spacing = grid::unit(0.35, "lines"),
      axis.text = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  grDevices::cairo_pdf(
    front_pdf,
    width = 15,
    height = 10,
    family = "Helvetica",
    bg = "white",
    onefile = TRUE
  )
  print(combined_plot)
  print(sample_plot)
  print(multiplier_plot)
  grDevices::dev.off()
  log_info("Wrote three footprint-bound sensitivity pages: ", front_pdf)
  invisible(front_pdf)
}

bound_path <- function(p_id) {
  file.path(work_root, paste0("fp_bound_", p_id, ".rds"))
}

prepare_multiomic <- function() {
  metadata <- load_metadata()
  strict <- build_strict_cutoffs(metadata)
  output_files <- c(multiomic_path, vapply(p_grid$p_id, bound_path, character(1L)))
  if (all(file.exists(output_files))) {
    log_info("Reusing strict-bound HPAFII multiomic inputs.")
    return(invisible(output_files))
  }
  log_info("Loading aligned footprint cache for the HPAFII p-only rebuild.")
  fp_aligned <- craftgrn:::load_fp_aligned_from_cache(
    cache_dir = cache_root,
    cache_tag = "JASPAR2024",
    output_mode = "full",
    load_id_map = FALSE,
    cache_format = "csv",
    expand_compact_annotation = TRUE,
    verbose = TRUE
  )
  p10 <- strict[p_id == "p1e10"]
  p10_threshold <- stats::setNames(p10$cutoff, p10$sample_id)
  if (file.exists(multiomic_path)) {
    log_info("Reusing the common HPAFII multiomic universe at p = 1e-10.")
    multiomic <- readRDS(multiomic_path)
  } else {
    log_info("Preparing the common HPAFII multiomic universe at p = 1e-10.")
    multiomic <- craftgrn::load_prep_multiomic_data(
      config = project_config_path,
      fp_aligned = fp_aligned,
      metadata = as.data.frame(metadata),
      label_col = "name",
      expected_n = 33L,
      threshold_gene_expr = 10,
      threshold_fp_score = p10_threshold,
      write_outputs = FALSE,
      use_parallel = TRUE,
      verbose = TRUE
    )
    saveRDS(multiomic, multiomic_path, compress = FALSE)
  }

  fp_ids <- rownames(multiomic$matrices$fp_score)
  fp_index <- match(fp_ids, fp_aligned$fp_score$peak_ID)
  if (anyNA(fp_index)) stop("Prepared footprint ids are missing from the raw cache.")
  sample_columns <- c("peak_ID", metadata$id)
  score_sub <- fp_aligned$fp_score[fp_index, sample_columns, drop = FALSE]
  bound_sub <- fp_aligned$fp_bound[fp_index, sample_columns, drop = FALSE]
  fp_annotation <- data.frame(
    fp_peak = as.character(multiomic$features$fp$fp_id),
    atac_peak = as.character(multiomic$features$fp$atac_peak),
    stringsAsFactors = FALSE
  )
  atac <- data.table::fread(atac_path, showProgress = FALSE)
  atac_loaded <- craftgrn:::load_atac(as.data.frame(atac), sort_peaks = TRUE)
  conditions <- colnames(multiomic$matrices$fp_bound)
  for (i in seq_len(nrow(p_grid))) {
    p_id_use <- p_grid$p_id[[i]]
    one <- strict[p_id == p_id_use]
    thresholds <- stats::setNames(one$cutoff, one$sample_id)
    log_info("Building strict condition-bound matrix for ", p_grid$p_label[[i]], ".")
    condition_bound <- craftgrn:::make_fp_bound_condition(
      fp_bound_tbl = bound_sub,
      fp_score_tbl = score_sub,
      atac_overlap_tbl = atac_loaded$overlap,
      fp_annotation_tbl = fp_annotation,
      metadata = as.data.frame(metadata),
      label_col = "name",
      threshold_fp_score = thresholds,
      min_samples = 1L
    )
    row_index <- match(fp_ids, condition_bound$peak_ID)
    if (anyNA(row_index)) stop("Strict condition-bound rows do not align.")
    missing_conditions <- setdiff(conditions, names(condition_bound))
    if (length(missing_conditions)) {
      stop("Strict bound output is missing conditions: ", paste(missing_conditions, collapse = ", "))
    }
    bound_matrix <- as.matrix(condition_bound[row_index, conditions, drop = FALSE]) > 0
    rownames(bound_matrix) <- fp_ids
    colnames(bound_matrix) <- conditions
    if (identical(p_id_use, "p1e10") &&
        !identical(bound_matrix, multiomic$matrices$fp_bound)) {
      stop("Recomputed p = 1e-10 bound matrix does not match the prepared object.")
    }
    saveRDS(bound_matrix, bound_path(p_id_use), compress = FALSE)
  }
  strict_counts <- data.table::rbindlist(lapply(p_grid$p_id, function(p_id) {
    matrix <- readRDS(bound_path(p_id))
    data.table::data.table(
      p_id = p_id,
      condition_id = colnames(matrix),
      bound_footprints = colSums(matrix)
    )
  }))
  strict_counts[, p_order := match(p_id, p_grid$p_id)]
  data.table::setorder(strict_counts, condition_id, p_order)
  nested <- strict_counts[, all(diff(bound_footprints) <= 0), by = condition_id]
  if (any(!nested$V1)) stop("Stricter p cutoffs increased a condition bound count.")
  strict_counts[, p_order := NULL]
  data.table::fwrite(
    strict_counts,
    file.path(work_root, "strict_condition_bound_counts.csv")
  )
  rm(fp_aligned, score_sub, bound_sub, atac, atac_loaded, multiomic)
  invisible(gc())
  log_info("Saved strict-bound HPAFII multiomic inputs.")
  invisible(output_files)
}

run_module1 <- function() {
  prepare_multiomic()
  predicted_manifest <- file.path(
    module1_root,
    "module1_predicted_tfbs_manifest.csv"
  )
  if (file.exists(module1_marker) && file.exists(predicted_manifest)) {
    log_info("Reusing completed Module 1 strict-bound run: ", module1_root)
    return(invisible(predicted_manifest))
  }
  multiomic <- readRDS(multiomic_path)
  dir.create(module1_root, recursive = TRUE, showWarnings = FALSE)
  log_info("Running Module 1 at R >= 0.5 on the p = 1e-10 universe.")
  result <- craftgrn::predict_tfbs(
    omics_data = multiomic,
    out_dir = module1_root,
    db = "JASPAR2024",
    label_col = "name",
    project_config = project_config_path,
    r_cutoff = 0.5,
    p_cutoff = NULL,
    fdr_cutoff = NULL,
    filter_to_canonical_bound = TRUE,
    write_outputs = TRUE,
    write_stats = FALSE,
    write_bed = FALSE,
    write_qc_report = FALSE,
    output_format = "parquet",
    return_prediction_stats = FALSE,
    cores = parse_workers(),
    verbose = TRUE
  )
  manifest <- result$reports$predicted_tfbs_manifest
  if (!is.character(manifest) || !file.exists(manifest)) {
    stop("Module 1 did not write the predicted TFBS manifest.")
  }
  writeLines(
    c(
      paste0("completed=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "module1_r=0.5",
      "bound_cutoff=p1e10",
      paste0("predicted_tfbs_manifest=", normalizePath(manifest, winslash = "/"))
    ),
    module1_marker,
    useBytes = TRUE
  )
  rm(result, multiomic)
  invisible(gc())
  log_info("Completed Module 1 strict-bound run.")
  invisible(manifest)
}

if (stage %in% c("all", "plots")) build_front_pages()
if (stage %in% c("all", "prepare")) prepare_multiomic()
if (stage %in% c("all", "module1")) run_module1()

invisible(TRUE)
