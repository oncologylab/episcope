#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(args[[1L]]) else "all"
allowed_stages <- c("all", "module2", "documents", "report")
if (!stage %in% allowed_stages) {
  stop("Stage must be one of: ", paste(allowed_stages, collapse = ", "))
}

required_packages <- c(
  "arrow", "craftgrn", "data.table", "ggplot2", "patchwork", "scales",
  "tidyselect"
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
module2_root <- file.path(work_root, "module2_tf0p5_fp0p3_p1e10")
raw_doc_root <- file.path(work_root, "document_terms_raw")
final_doc_root <- file.path(work_root, "document_terms_final")
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

multiomic_path <- file.path(work_root, "multiomic_p1e10.rds")
predicted_manifest_path <- file.path(
  module1_root,
  "module1_predicted_tfbs_manifest.csv"
)
module2_link_manifest_path <- file.path(
  module2_root,
  "data",
  "links",
  "module2_links_manifest.csv"
)
module2_marker <- file.path(module2_root, "strict_bound_module2_complete.txt")
condition_root <- file.path(analysis_root, "condition_links")
condition_manifest_path <- file.path(condition_root, "condition_links_manifest.csv")
expression_path <- file.path(condition_root, "condition_gene_expression.csv")
gene_union_path <- file.path(condition_root, "global_differential_gene_union.csv")
project_config_path <- file.path(project_root, "project.yaml")
gene_tss_path <- file.path(
  project_root,
  "data",
  "gene_tss_hg38_GENCODEv50_project.csv"
)
document_summary_path <- file.path(work_root, "p_only_document_summary.csv")
compact_summary_path <- file.path(work_root, "p_only_compact_qc_summary.csv")
union_summary_path <- file.path(work_root, "p_only_union_summary.csv")
document_audit_path <- file.path(work_root, "p_only_document_build_audit.csv")
page_pdf <- file.path(work_root, "p_only_document_comparison.pdf")

p_grid <- data.table::data.table(
  p_id = c("p1e10", "p1e12", "p1e15", "p1e20"),
  p_label = c("p = 1e-10", "p = 1e-12", "p = 1e-15", "p = 1e-20"),
  p_order = 1:4
)
p_colors <- stats::setNames(
  c("#2C7BB6", "#00A6CA", "#F4A261", "#B2182B"),
  p_grid$p_label
)

bound_path <- function(p_id) {
  file.path(work_root, paste0("fp_bound_", p_id, ".rds"))
}

required_inputs <- c(
  multiomic_path,
  predicted_manifest_path,
  vapply(p_grid$p_id, bound_path, character(1L)),
  condition_manifest_path,
  expression_path,
  gene_union_path,
  project_config_path,
  gene_tss_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Missing strict-bound input(s). Run 01_build_hpafii_strict_bound_sensitivity.R first: ",
    paste(missing_inputs, collapse = ", ")
  )
}

run_module2 <- function() {
  if (file.exists(module2_marker) && file.exists(module2_link_manifest_path)) {
    log_info("Reusing completed Module 2 strict-bound run: ", module2_root)
    return(invisible(module2_link_manifest_path))
  }
  multiomic <- readRDS(multiomic_path)
  config <- list(
    ref_genome = "hg38",
    gene_tss = gene_tss_path,
    link_window_bp = 100000,
    module2 = list(
      threshold_tf_target_corr_r = 0.5,
      threshold_fp_target_corr_r = 0.3
    )
  )
  dir.create(module2_root, recursive = TRUE, showWarnings = FALSE)
  log_info("Running Module 2 at TF-target R >= 0.5 and FP-target R >= 0.3.")
  result <- craftgrn::predict_tf_targets(
    multiomic_data = multiomic,
    predicted_tfbs = predicted_manifest_path,
    gene_tss = gene_tss_path,
    project_config = config,
    output_dir = module2_root,
    max_distance_bp = 100000,
    n_cores = parse_workers(),
    output_format = "parquet",
    verbose = TRUE,
    write_qc_report = FALSE,
    qc_report_validate = FALSE
  )
  manifest <- result$reports$links_manifest
  if (!is.character(manifest) || !file.exists(manifest)) {
    stop("Module 2 did not write its link manifest.")
  }
  writeLines(
    c(
      paste0("completed=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "module1_r=0.5",
      "tf_target_r=0.5",
      "fp_target_r=0.3",
      "bound_cutoff=p1e10",
      paste0("links_manifest=", normalizePath(manifest, winslash = "/"))
    ),
    module2_marker,
    useBytes = TRUE
  )
  rm(result, multiomic)
  invisible(gc())
  log_info("Completed Module 2 strict-bound run.")
  invisible(manifest)
}

read_link_chunk <- function(path, format) {
  columns <- c("tf", "fp_id", "target_gene", "module2_link_pass")
  if (identical(format, "parquet") || grepl("[.]parquet$", path, ignore.case = TRUE)) {
    return(data.table::as.data.table(arrow::read_parquet(
      path,
      col_select = tidyselect::all_of(columns),
      as_data_frame = TRUE
    )))
  }
  data.table::fread(path, select = columns, showProgress = FALSE)
}

raw_doc_path <- function(p_id, condition_id) {
  file.path(
    raw_doc_root,
    p_id,
    paste0(gsub("[^A-Za-z0-9_.-]+", "_", condition_id), ".parquet")
  )
}

condition_audit_path <- function(condition_id) {
  file.path(
    raw_doc_root,
    "audit",
    paste0(gsub("[^A-Za-z0-9_.-]+", "_", condition_id), ".csv")
  )
}

build_condition_doc_terms <- function(
    condition_id,
    multiomic,
    bounds,
    link_manifest,
    gene_union,
    specificity) {
  cached_paths <- vapply(
    p_grid$p_id,
    raw_doc_path,
    condition_id = condition_id,
    FUN.VALUE = character(1L)
  )
  cached_audit_path <- condition_audit_path(condition_id)
  if (all(file.exists(cached_paths)) && file.exists(cached_audit_path)) {
    log_info("Reusing raw document terms for ", condition_id, ".")
    return(data.table::fread(cached_audit_path, showProgress = FALSE))
  }
  fp_ids <- rownames(multiomic$matrices$fp_score)
  gene_ids <- rownames(multiomic$matrices$gene_expr)
  fp_score <- multiomic$matrices$fp_score[, condition_id]
  gene_expr <- multiomic$matrices$gene_expr[, condition_id]
  edge_parts <- vector("list", nrow(link_manifest))
  source_links <- 0
  for (i in seq_len(nrow(link_manifest))) {
    links <- read_link_chunk(
      as.character(link_manifest$path[[i]]),
      as.character(link_manifest$format[[i]])
    )
    source_links <- source_links + nrow(links)
    links <- links[module2_link_pass %in% TRUE]
    if (!nrow(links)) next
    fp_index <- match(links$fp_id, fp_ids)
    target_index <- match(links$target_gene, gene_ids)
    tf_index <- match(links$tf, gene_ids)
    valid <- !is.na(fp_index) & !is.na(target_index) & !is.na(tf_index) &
      links$target_gene %in% gene_union &
      bounds[[1L]][cbind(fp_index, match(condition_id, colnames(bounds[[1L]])))] &
      is.finite(gene_expr[target_index]) & gene_expr[target_index] >= 10 &
      is.finite(gene_expr[tf_index])
    if (!any(valid)) next
    links <- links[valid]
    fp_index <- fp_index[valid]
    target_index <- target_index[valid]
    tf_index <- tf_index[valid]
    strict_level <- rep(1L, nrow(links))
    for (p_i in 2:4) {
      strict_level <- strict_level + as.integer(
        bounds[[p_i]][cbind(fp_index, match(condition_id, colnames(bounds[[p_i]])))]
      )
    }
    edge_parts[[i]] <- data.table::data.table(
      condition_label = condition_id,
      tf_doc = as.character(links$tf),
      tf = as.character(links$tf),
      gene_key = as.character(links$target_gene),
      peak_id = as.character(links$fp_id),
      fp_score_condition = as.numeric(fp_score[fp_index]),
      gene_expr_condition = as.numeric(gene_expr[target_index]),
      tf_expr_condition = as.numeric(gene_expr[tf_index]),
      strict_level = strict_level
    )
    rm(links)
  }
  edges <- data.table::rbindlist(edge_parts, use.names = TRUE, fill = TRUE)
  if (!nrow(edges)) stop("No eligible links remain for ", condition_id, ".")
  audit <- vector("list", nrow(p_grid))
  for (p_i in seq_len(nrow(p_grid))) {
    p_id <- p_grid$p_id[[p_i]]
    selected <- edges[strict_level >= p_i]
    if (!nrow(selected)) {
      stop("No links remain for ", condition_id, " at ", p_grid$p_label[[p_i]], ".")
    }
    log_info(
      "Building ", p_grid$p_label[[p_i]], " terms for ", condition_id,
      " from ", format(nrow(selected), big.mark = ","), " links."
    )
    doc_term <- craftgrn:::build_doc_term_condition_union(
      edges_condition = selected,
      count_method = "log",
      count_scale = 50,
      prefix_terms = TRUE,
      threshold_gene_expr = 10,
      threshold_fp_score = -Inf,
      threshold_tf_expr = -Inf,
      include_tf_terms = FALSE,
      require_tf_expr = TRUE,
      fp_term_mode = "unique",
      condition_peak_weighting = "none",
      condition_gene_specificity = specificity,
      balance_mode = "min",
      check_repeated_values = FALSE
    )
    if (!nrow(doc_term)) {
      stop("Document construction returned no terms for ", condition_id, ".")
    }
    path <- raw_doc_path(p_id, condition_id)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(doc_term, path, compression = "zstd")
    audit[[p_i]] <- data.table::data.table(
      p_id = p_id,
      condition_id = condition_id,
      raw_doc_path = normalizePath(path, winslash = "/", mustWork = TRUE),
      source_links = source_links,
      eligible_links = nrow(selected)
    )
    rm(selected, doc_term)
    invisible(gc())
  }
  rm(edges, edge_parts)
  invisible(gc())
  audit <- data.table::rbindlist(audit)
  dir.create(dirname(cached_audit_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(audit, cached_audit_path)
  audit
}

apply_tf_expression_peak_weight <- function(doc_term, multiomic) {
  documents <- unique(doc_term[, .(doc_id)])
  documents[, `:=`(
    condition_id = sub("::[^:]+$", "", doc_id),
    tf = sub("^.*::", "", doc_id)
  )]
  gene_ids <- rownames(multiomic$matrices$gene_expr)
  documents[, tf_index := match(tf, gene_ids)]
  documents[, condition_index := match(
    condition_id,
    colnames(multiomic$matrices$gene_expr)
  )]
  if (anyNA(documents$tf_index) || anyNA(documents$condition_index)) {
    stop("Unable to align a document TF or condition for Peak weighting.")
  }
  documents[, tf_expr_condition := as.numeric(
    multiomic$matrices$gene_expr[cbind(tf_index, condition_index)]
  )]
  invalid_documents <- documents[
    !is.finite(tf_expr_condition) | tf_expr_condition <= 0,
    doc_id
  ]
  if (length(invalid_documents)) {
    log_info(
      "Removing ", length(invalid_documents),
      " TF document(s) with zero condition expression before Peak weighting."
    )
    doc_term <- doc_term[!doc_id %in% invalid_documents]
    documents <- documents[!doc_id %in% invalid_documents]
  }
  factor_result <- craftgrn:::.topic_condition_tf_expression_peak_factors(
    documents[, .(doc_id, tf, tf_expr_condition)]
  )
  doc_term[
    factor_result$factors,
    peak_multiplier := i.peak_multiplier,
    on = "doc_id"
  ]
  peak_rows <- startsWith(doc_term$term_id, "PEAK:")
  if (any(peak_rows & (!is.finite(doc_term$peak_multiplier) |
      doc_term$peak_multiplier <= 0))) {
    stop("A Peak term is missing its TF-expression multiplier.")
  }
  doc_term[peak_rows, weight := weight * peak_multiplier]
  doc_term[, peak_multiplier := NULL]
  craftgrn:::.topic_apply_counts(
    doc_term,
    count_method = "log",
    count_scale = 50
  )
}

finalize_one_p <- function(p_id_use, audit, multiomic) {
  final_path <- file.path(final_doc_root, paste0("document_terms_", p_id_use, ".parquet"))
  if (file.exists(final_path)) {
    log_info("Reusing finalized document terms for ", p_id_use, ".")
    return(data.table::as.data.table(arrow::read_parquet(final_path)))
  }
  paths <- audit[p_id == p_id_use, raw_doc_path]
  doc_term <- data.table::rbindlist(lapply(paths, function(path) {
    data.table::as.data.table(arrow::read_parquet(path))
  }), use.names = TRUE, fill = TRUE)
  if (anyDuplicated(doc_term, by = c("doc_id", "term_id"))) {
    stop("Raw document terms are duplicated across conditions for ", p_id_use, ".")
  }
  doc_term <- apply_tf_expression_peak_weight(doc_term, multiomic)
  token_cap <- craftgrn:::.cap_warplda_token_counts(doc_term$pseudo_count_log)
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
  token_totals <- token_audit[, .(final_tokens = sum(final_tokens)), by = modality]
  ratio <- token_totals[modality == "peak", final_tokens] /
    token_totals[modality == "gene", final_tokens]
  if (!length(ratio) || !is.finite(ratio) || abs(ratio - 0.5) > 1e-4) {
    stop("Final Peak/Gene token ratio is invalid for ", p_id_use, ".")
  }
  dir.create(final_doc_root, recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(doc_term, final_path, compression = "zstd")
  doc_term
}

summarize_documents <- function(doc_terms, audit) {
  compact_parts <- vector("list", nrow(p_grid))
  union_parts <- vector("list", nrow(p_grid))
  summary_parts <- vector("list", nrow(p_grid))
  for (i in seq_len(nrow(p_grid))) {
    p_id_use <- p_grid$p_id[[i]]
    p_label <- p_grid$p_label[[i]]
    doc_term <- doc_terms[[p_id_use]]
    doc_term[, `:=`(
      condition_id = sub("::[^:]+$", "", doc_id),
      term_type = data.table::fifelse(
        startsWith(term_id, "GENE:"),
        "Gene",
        "Peak"
      )
    )]
    per_document <- doc_term[, .(terms_in_document = .N),
      by = .(condition_id, doc_id, term_type)]
    typical <- per_document[, .(
      median_terms = as.numeric(stats::median(terms_in_document)),
      lower_quartile = as.numeric(stats::quantile(
        terms_in_document,
        0.25,
        names = FALSE
      )),
      upper_quartile = as.numeric(stats::quantile(
        terms_in_document,
        0.75,
        names = FALSE
      )),
      tf_documents = data.table::uniqueN(doc_id)
    ), by = .(condition_id, term_type)]
    unique_terms <- doc_term[, .(
      unique_terms = data.table::uniqueN(term_id)
    ), by = .(condition_id, term_type)]
    compact <- merge(
      typical,
      unique_terms,
      by = c("condition_id", "term_type"),
      all = TRUE,
      sort = FALSE
    )
    compact[, `:=`(p_id = p_id_use, p_label = p_label)]
    compact_parts[[i]] <- compact
    union_parts[[i]] <- doc_term[, .(
      condition_id = "All conditions",
      unique_terms = data.table::uniqueN(term_id)
    ), by = term_type][, `:=`(p_id = p_id_use, p_label = p_label)]
    token_type <- data.table::fifelse(
      startsWith(doc_term$term_id, "GENE:"),
      "Gene",
      "Peak"
    )
    summary_parts[[i]] <- data.table::data.table(
      p_id = p_id_use,
      p_label = p_label,
      conditions = data.table::uniqueN(doc_term$condition_id),
      documents = data.table::uniqueN(doc_term$doc_id),
      document_term_rows = nrow(doc_term),
      unique_gene_terms = data.table::uniqueN(doc_term$term_id[token_type == "Gene"]),
      unique_peak_terms = data.table::uniqueN(doc_term$term_id[token_type == "Peak"]),
      eligible_links = sum(audit[p_id == p_id_use, eligible_links], na.rm = TRUE)
    )
  }
  compact <- data.table::rbindlist(compact_parts, use.names = TRUE, fill = TRUE)
  union <- data.table::rbindlist(union_parts, use.names = TRUE, fill = TRUE)
  summary <- data.table::rbindlist(summary_parts, use.names = TRUE, fill = TRUE)
  if (any(summary$conditions != 17L)) stop("A p-only setup is missing a condition.")
  numeric_fields <- c(
    "documents", "document_term_rows", "unique_gene_terms",
    "unique_peak_terms", "eligible_links"
  )
  for (field in numeric_fields) {
    if (any(diff(summary[[field]]) > 0)) {
      stop("A stricter p cutoff increased ", field, ".")
    }
  }
  data.table::fwrite(compact, compact_summary_path)
  data.table::fwrite(union, union_summary_path)
  data.table::fwrite(summary, document_summary_path)
  invisible(list(compact = compact, union = union, summary = summary))
}

build_documents <- function() {
  run_module2()
  complete <- all(file.exists(c(
    document_summary_path,
    compact_summary_path,
    union_summary_path,
    document_audit_path
  )))
  if (complete) {
    log_info("Reusing completed p-only document summaries.")
    return(invisible(TRUE))
  }
  multiomic <- readRDS(multiomic_path)
  bounds <- lapply(p_grid$p_id, function(p_id) readRDS(bound_path(p_id)))
  names(bounds) <- p_grid$p_id
  link_manifest <- data.table::fread(
    module2_link_manifest_path,
    showProgress = FALSE
  )
  if (!all(c("path", "format") %in% names(link_manifest)) ||
      any(!file.exists(link_manifest$path))) {
    stop("Module 2 link manifest is incomplete.")
  }
  condition_manifest <- data.table::fread(
    condition_manifest_path,
    showProgress = FALSE
  )
  conditions <- as.character(condition_manifest$condition_id)
  if (length(conditions) != 17L || anyDuplicated(conditions)) {
    stop("Expected 17 unique HPAFII conditions.")
  }
  missing_conditions <- setdiff(conditions, colnames(multiomic$matrices$fp_score))
  if (length(missing_conditions)) {
    stop("Multiomic object is missing: ", paste(missing_conditions, collapse = ", "))
  }
  gene_union <- data.table::fread(gene_union_path, showProgress = FALSE)
  gene_union <- unique(as.character(gene_union[pass_abs_log2fc %in% TRUE, gene_key]))
  if (!length(gene_union)) stop("Global differential-gene union is empty.")
  specificity <- craftgrn:::.module3_condition_gene_specificity_lookup(
    expression_file = expression_path,
    expression_min = 10,
    temperature = 0.5,
    uniform_floor = 0.1
  )
  audits <- vector("list", length(conditions))
  for (i in seq_along(conditions)) {
    log_info("Preparing strict-bound documents ", i, "/", length(conditions), ": ", conditions[[i]], ".")
    audits[[i]] <- build_condition_doc_terms(
      condition_id = conditions[[i]],
      multiomic = multiomic,
      bounds = bounds,
      link_manifest = link_manifest,
      gene_union = gene_union,
      specificity = specificity
    )
  }
  audit <- data.table::rbindlist(audits, use.names = TRUE, fill = TRUE)
  if (anyNA(audit$eligible_links)) {
    stop("Cached raw document terms lack a complete link-count audit.")
  }
  data.table::fwrite(audit, document_audit_path)
  doc_terms <- lapply(p_grid$p_id, function(p_id) {
    finalize_one_p(p_id, audit = audit, multiomic = multiomic)
  })
  names(doc_terms) <- p_grid$p_id
  summarize_documents(doc_terms, audit)
  rm(doc_terms, audit, audits, bounds, multiomic)
  invisible(gc())
  log_info("Completed p-only strict-bound document summaries.")
  invisible(TRUE)
}

simple_condition_label <- function(value) {
  value <- sub("^HPAFII_", "", value)
  value <- sub("_Ctrl$", "", value)
  value <- gsub("Gln[.]Arg", "Gln + Arg", value)
  value <- gsub("Met[.]Cys", "Met + Cys", value)
  gsub("_", " ", value, fixed = TRUE)
}

report_theme <- function(base_size = 9) {
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
      legend.position = "bottom"
    )
}

build_report_page <- function() {
  build_documents()
  compact <- data.table::fread(compact_summary_path, showProgress = FALSE)
  union <- data.table::fread(union_summary_path, showProgress = FALSE)
  condition_manifest <- data.table::fread(
    condition_manifest_path,
    showProgress = FALSE
  )
  condition_levels <- rev(simple_condition_label(condition_manifest$condition_id))
  compact[, `:=`(
    condition_display = factor(
      simple_condition_label(condition_id),
      levels = condition_levels
    ),
    term_type = factor(term_type, levels = c("Gene", "Peak")),
    p_label = factor(p_label, levels = p_grid$p_label)
  )]
  typical_plot <- ggplot2::ggplot(
    compact,
    ggplot2::aes(condition_display, median_terms, fill = p_label)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.82, preserve = "single"),
      width = 0.72,
      color = "#222222",
      linewidth = 0.15
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~term_type, scales = "free_x", ncol = 2) +
    ggplot2::scale_fill_manual(values = p_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      title = "Terms in a typical TF document",
      x = NULL,
      y = "Terms",
      fill = NULL
    ) +
    report_theme(8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7)
    )

  unique_data <- data.table::rbindlist(list(
    compact[, .(condition_id, term_type, p_id, p_label, unique_terms)],
    union[, .(condition_id, term_type, p_id, p_label, unique_terms)]
  ))
  unique_levels <- c("All conditions", simple_condition_label(condition_manifest$condition_id))
  unique_data[, `:=`(
    condition_display = factor(
      data.table::fifelse(
        condition_id == "All conditions",
        condition_id,
        simple_condition_label(condition_id)
      ),
      levels = rev(unique_levels)
    ),
    term_type = factor(term_type, levels = c("Gene", "Peak")),
    p_label = factor(p_label, levels = p_grid$p_label)
  )]
  unique_plot <- ggplot2::ggplot(
    unique_data,
    ggplot2::aes(condition_display, unique_terms, fill = p_label)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.82, preserve = "single"),
      width = 0.72,
      color = "#222222",
      linewidth = 0.15
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~term_type, scales = "free_x", ncol = 2) +
    ggplot2::scale_fill_manual(values = p_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(big.mark = ","),
      expand = ggplot2::expansion(mult = c(0, 0.04))
    ) +
    ggplot2::labs(
      title = "Unique terms in each condition",
      x = NULL,
      y = "Unique terms",
      fill = NULL
    ) +
    report_theme(8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7)
    )

  page <- (typical_plot | unique_plot) +
    patchwork::plot_layout(widths = c(1, 1.05), guides = "collect") +
    patchwork::plot_annotation(
      title = "How stricter bound calls change the documents",
      subtitle = paste0(
        "Fixed cutoffs: Module 1 / TF-target / FP-target = 0.5 / 0.5 / 0.3. ",
        "No extra footprint-score cutoff."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          family = "Helvetica",
          face = "bold",
          size = 18
        ),
        plot.subtitle = ggplot2::element_text(
          family = "Helvetica",
          face = "bold",
          size = 11
        )
      )
    ) &
    ggplot2::theme(legend.position = "bottom")
  grDevices::cairo_pdf(
    page_pdf,
    width = 15,
    height = 10,
    family = "Helvetica",
    bg = "white",
    onefile = TRUE
  )
  print(page)
  grDevices::dev.off()
  log_info("Wrote p-only document comparison page: ", page_pdf)
  invisible(page_pdf)
}

if (stage %in% c("all", "module2")) run_module2()
if (stage %in% c("all", "documents")) build_documents()
if (stage %in% c("all", "report")) build_report_page()

invisible(TRUE)
