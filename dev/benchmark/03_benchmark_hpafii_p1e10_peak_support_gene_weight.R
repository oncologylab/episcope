#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) tolower(args[[1L]]) else "all"
allowed_stages <- c("all", "documents", "train", "qc", "compare", "deliver")
if (!stage %in% allowed_stages) {
  stop("Stage must be one of: ", paste(allowed_stages, collapse = ", "))
}

required_packages <- c(
  "arrow", "craftgrn", "data.table", "digest", "ggplot2", "Matrix",
  "patchwork", "scales", "yaml"
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

log_warn <- function(...) {
  craftgrn:::.log_warn(paste0(...))
}

detect_physical_cores <- function() {
  topology <- suppressWarnings(system2(
    "lscpu",
    "-p=CORE,SOCKET",
    stdout = TRUE,
    stderr = FALSE
  ))
  topology <- topology[!startsWith(topology, "#") & nzchar(topology)]
  physical <- if (length(topology)) {
    length(unique(topology))
  } else {
    parallel::detectCores(logical = FALSE)
  }
  if (!is.finite(physical) || physical < 1L) physical <- 1L
  as.integer(physical)
}

parse_workers <- function() {
  physical <- detect_physical_cores()
  concurrent_models <- suppressWarnings(as.integer(Sys.getenv(
    "CRAFTGRN_P1E10_CONCURRENT_MODELS",
    unset = "2"
  )))
  if (!is.finite(concurrent_models) || concurrent_models < 1L) {
    stop("Invalid concurrent-model count.")
  }
  workers <- suppressWarnings(as.integer(Sys.getenv(
    "CRAFTGRN_P1E10_BENCHMARK_WORKERS",
    unset = as.character(max(1L, floor(physical / concurrent_models)))
  )))
  if (!is.finite(workers) || workers < 1L) stop("Invalid worker count.")
  min(workers, physical, 18L)
}

parse_concurrent_models <- function() {
  jobs <- suppressWarnings(as.integer(Sys.getenv(
    "CRAFTGRN_P1E10_CONCURRENT_MODELS",
    unset = "2"
  )))
  if (!is.finite(jobs) || jobs < 1L) stop("Invalid concurrent-model count.")
  min(jobs, 2L)
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
repository_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)
script_path <- normalizePath(
  file.path(
    repository_root,
    "dev",
    "benchmark",
    "03_benchmark_hpafii_p1e10_peak_support_gene_weight.R"
  ),
  winslash = "/",
  mustWork = TRUE
)
analysis_root <- file.path(
  project_root,
  "regulatory_topics_hpafii_condition_models_de_gene_filtered"
)
topic_runs <- file.path(analysis_root, "topic_runs")
strict_root <- file.path(
  topic_runs,
  "review",
  "method10_module1_module2_r_cutoff_sensitivity",
  "hpa_strict_bound_sensitivity"
)
benchmark_root <- file.path(
  topic_runs,
  "review",
  "method10_p1e10_peak_support_gene_weight_benchmark"
)
cache_root <- file.path(
  topic_runs,
  "cache",
  "method10_p1e10_peak_support_gene_weight_benchmark"
)
delivery_root <- file.path(benchmark_root, "delivery")
box_parent <- paste0(
  "yilab:Yaoxiang/thesis_project/manuscript/module3_topic_model/",
  "nutrient_hpafii_input_design_qc"
)
box_root <- paste0(
  box_parent,
  "/16_Method10_P1e10_PeakSupport_GeneWeight_Benchmark"
)
box_plans_root <- "yilab:Yaoxiang/thesis_project/plans"
plan_date <- "2026-09-02"
box_plan_names <- c(
  condition_specific = paste0(
    plan_date,
    "_HPAFII_Method10_P1e10_ConditionSpecificPeaks.yaml"
  ),
  tf_union = paste0(
    plan_date,
    "_HPAFII_Method10_P1e10_TFUnionPeaks.yaml"
  )
)

raw_doc_root <- file.path(strict_root, "document_terms_raw", "p1e10")
baseline_doc_path <- file.path(
  strict_root,
  "document_terms_final",
  "document_terms_p1e10.parquet"
)
multiomic_path <- file.path(strict_root, "multiomic_p1e10.rds")
bound_path <- file.path(strict_root, "fp_bound_p1e10.rds")
module2_manifest_path <- file.path(
  strict_root,
  "module2_tf0p5_fp0p3_p1e10",
  "data",
  "links",
  "module2_links_manifest.csv"
)
gene_union_path <- file.path(
  analysis_root,
  "condition_links",
  "global_differential_gene_union.csv"
)
condition_manifest_path <- file.path(
  analysis_root,
  "condition_links",
  "condition_links_manifest.csv"
)
project_config_path <- file.path(project_root, "project.yaml")
reference_qc_root <- file.path(
  topic_runs,
  "review",
  "condition_unique_fp_gene_separation_benchmark",
  "topic_assignment_qc"
)
reference_tf_sets_path <- file.path(
  reference_qc_root,
  "tf_target_gene_eligible_pairs.rds"
)
reference_tf_panel_path <- file.path(
  reference_qc_root,
  "tf_target_gene_umap_tf_panel.csv"
)

support_paths <- c(
  condition_specific = file.path(cache_root, "condition_specific_support.parquet"),
  tf_union = file.path(cache_root, "tf_union_support.parquet")
)
budget_path <- file.path(cache_root, "common_document_token_budget.csv")
support_summary_path <- file.path(cache_root, "document_support_summary.csv")
peak_gene_map_path <- file.path(cache_root, "peak_gene_map.csv")
support_manifest_path <- file.path(cache_root, "support_manifest.csv")
source_signature_path <- file.path(cache_root, "source_signature.txt")
source_identity_cache_path <- file.path(cache_root, "source_file_identities.rds")

ratio_levels <- c(2L, 3L, 4L, 6L, 8L, 10L)
run_table <- data.table::rbindlist(list(
  data.table::data.table(
    run_number = c(27:29, 33:35),
    support = "condition_specific",
    gene_to_peak = ratio_levels
  ),
  data.table::data.table(
    run_number = c(30:32, 36:38),
    support = "tf_union",
    gene_to_peak = ratio_levels
  )
))
data.table::setorder(run_table, run_number)
run_table[, support_label := data.table::fifelse(
  support == "condition_specific",
  "Condition-specific Peaks",
  "TF-union Peaks"
)]
run_table[, run_id := sprintf(
  "run_%03d_lda_condition_tf_unique_fp_p1e10_%s_K30_g%dp1",
  run_number,
  support,
  gene_to_peak
)]
run_table[, run_root := file.path(topic_runs, run_id)]
run_table[, model_root := file.path(run_root, "topic_models")]
run_table[, extraction_root := file.path(run_root, "topic_extraction", "K30")]
run_table[, qc_name := sprintf(
  "LDA_K30_PeakSupport-%s_Ratio-Gene%d-Peak1_QC.pdf",
  data.table::fifelse(
    support == "condition_specific",
    "ConditionSpecific",
    "TFUnion"
  ),
  gene_to_peak
)]

required_inputs <- c(
  baseline_doc_path,
  multiomic_path,
  bound_path,
  module2_manifest_path,
  gene_union_path,
  condition_manifest_path,
  project_config_path
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing benchmark input(s): ", paste(missing_inputs, collapse = ", "))
}
raw_doc_paths <- sort(list.files(
  raw_doc_root,
  pattern = "[.]parquet$",
  full.names = TRUE
))
if (length(raw_doc_paths) != 17L) {
  stop("Expected 17 p = 1e-10 raw condition document files.")
}

source_identity_cache <- if (file.exists(source_identity_cache_path)) {
  readRDS(source_identity_cache_path)
} else {
  list()
}

file_identity <- function(path) {
  info <- file.info(path)
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  bytes <- as.numeric(info$size)
  modified <- format(info$mtime, tz = "UTC", usetz = TRUE)
  cached <- source_identity_cache[[normalized_path]]
  if (!is.null(cached) && identical(cached$bytes, bytes) &&
      identical(cached$modified, modified)) {
    return(cached)
  }
  value <- list(
    path = normalized_path,
    bytes = bytes,
    modified = modified,
    md5 = unname(tools::md5sum(path))
  )
  source_identity_cache[[normalized_path]] <<- value
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  saveRDS(source_identity_cache, source_identity_cache_path)
  value
}

source_signature <- digest::digest(
  list(
    baseline = file_identity(baseline_doc_path),
    raw_documents = lapply(raw_doc_paths, file_identity),
    multiomic = file_identity(multiomic_path),
    bound = file_identity(bound_path),
    module2_manifest = file_identity(module2_manifest_path),
    gene_union = file_identity(gene_union_path),
    filters = list(
      module1_r = 0.5,
      tf_target_r = 0.5,
      fp_target_r = 0.3,
      footprint_p = 1e-10,
      target_expression_min = 10
    ),
    ratios = c(2, 3, 4),
    fixed_token_budget = TRUE
  ),
  algo = "xxhash64"
)

check_source_signature <- function() {
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(source_signature_path)) {
    old <- trimws(readLines(source_signature_path, warn = FALSE)[[1L]])
    if (!identical(old, source_signature)) {
      stop(
        "Existing p = 1e-10 Peak-support cache has a different source ",
        "signature. Use a new benchmark root instead of mixing results."
      )
    }
  } else {
    writeLines(source_signature, source_signature_path, useBytes = TRUE)
  }
  invisible(source_signature)
}

read_probability <- function(path) {
  value <- data.table::fread(path, showProgress = FALSE)
  ids <- as.character(value[[1L]])
  value[[1L]] <- NULL
  matrix <- as.matrix(value)
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- ids
  matrix
}

write_probability <- function(matrix, id_name, path) {
  value <- data.table::as.data.table(matrix, keep.rownames = id_name)
  data.table::fwrite(value, path)
  invisible(path)
}

load_baseline_documents <- function() {
  data.table::as.data.table(arrow::read_parquet(
    baseline_doc_path,
    col_select = c("doc_id", "term_id", "pseudo_count_log"),
    as_data_frame = TRUE
  ))
}

load_raw_documents <- function(valid_documents) {
  value <- data.table::rbindlist(lapply(raw_doc_paths, function(path) {
    data.table::as.data.table(arrow::read_parquet(
      path,
      col_select = c("doc_id", "term_id", "weight"),
      as_data_frame = TRUE
    ))
  }), use.names = TRUE)
  value <- value[doc_id %in% valid_documents]
  if (!nrow(value) || anyDuplicated(value, by = c("doc_id", "term_id"))) {
    stop("Sparse source documents are empty or contain duplicate terms.")
  }
  value[, modality := data.table::fifelse(
    startsWith(term_id, "GENE:"),
    "gene",
    data.table::fifelse(startsWith(term_id, "PEAK:"), "peak", NA_character_)
  )]
  if (anyNA(value$modality) || any(!is.finite(value$weight) | value$weight <= 0)) {
    stop("Sparse source documents contain invalid terms or weights.")
  }
  value[]
}

document_index <- function(baseline) {
  docs <- unique(baseline[, .(doc_id)])
  docs[, `:=`(
    condition_id = sub("::[^:]+$", "", doc_id),
    tf = sub("^.*::", "", doc_id)
  )]
  if (nrow(docs) != 7381L || data.table::uniqueN(docs$tf) != 479L) {
    stop("Unexpected p = 1e-10 document universe.")
  }
  docs[]
}

tf_expression_factors <- function(docs, multiomic) {
  gene_expr <- multiomic$matrices$gene_expr
  tf_index <- match(docs$tf, rownames(gene_expr))
  condition_index <- match(docs$condition_id, colnames(gene_expr))
  if (anyNA(tf_index) || anyNA(condition_index)) {
    stop("Could not align a document TF or condition to the RNA matrix.")
  }
  factors <- data.table::copy(docs)
  factors[, tf_expr_condition := as.numeric(
    gene_expr[cbind(tf_index, condition_index)]
  )]
  if (any(!is.finite(factors$tf_expr_condition) |
      factors$tf_expr_condition <= 0)) {
    stop("A retained document has non-positive TF expression.")
  }
  result <- craftgrn:::.topic_condition_tf_expression_peak_factors(
    factors[, .(doc_id, tf, tf_expr_condition)]
  )
  data.table::fwrite(
    result$audit,
    file.path(cache_root, "tf_expression_peak_weighting_audit.csv")
  )
  result$factors
}

apply_peak_factor_and_counts <- function(value, factors) {
  value <- data.table::copy(value)
  value[factors, peak_multiplier := i.peak_multiplier, on = "doc_id"]
  invalid <- value[
    modality == "peak" &
      (!is.finite(peak_multiplier) | peak_multiplier <= 0)
  ]
  if (nrow(invalid)) stop("A Peak term is missing its TF-expression factor.")
  value[modality == "peak", weight := weight * peak_multiplier]
  value[, peak_multiplier := NULL]
  value[, source_count := craftgrn:::weight_to_count(
    weight,
    method = "log",
    scale = 50
  )]
  if (any(!is.finite(value$source_count) | value$source_count <= 0)) {
    stop("Log count conversion produced an invalid source count.")
  }
  value[]
}

balance_union_modalities <- function(value) {
  totals <- value[, .(total = sum(weight_raw)), by = .(doc_id, modality)]
  wide <- data.table::dcast(
    totals,
    doc_id ~ modality,
    value.var = "total",
    fill = 0
  )
  if (!all(c("gene", "peak") %in% names(wide)) ||
      any(wide$gene <= 0 | wide$peak <= 0)) {
    stop("TF-union balancing requires both modalities in every document.")
  }
  wide[, common_mass := pmin(gene, peak)]
  scale <- data.table::rbindlist(list(
    wide[, .(doc_id, modality = "gene", multiplier = common_mass / gene)],
    wide[, .(doc_id, modality = "peak", multiplier = common_mass / peak)]
  ))
  value[scale, multiplier := i.multiplier, on = c("doc_id", "modality")]
  if (any(!is.finite(value$multiplier) | value$multiplier <= 0)) {
    stop("TF-union modality balancing produced an invalid multiplier.")
  }
  value[, weight := weight_raw * multiplier]
  value[, multiplier := NULL]
  value[]
}

build_union_peak_rows <- function(sparse, docs, multiomic) {
  union_terms <- unique(sparse[modality == "peak", .(
    tf = sub("^.*::", "", doc_id),
    term_id
  )])
  by_tf_terms <- split(union_terms$term_id, union_terms$tf)
  by_tf_docs <- split(docs, docs$tf)
  fp_score <- multiomic$matrices$fp_score
  pieces <- vector("list", length(by_tf_terms))
  names_tf <- names(by_tf_terms)
  for (i in seq_along(names_tf)) {
    tf <- names_tf[[i]]
    terms <- sort(unique(by_tf_terms[[tf]]))
    one_docs <- by_tf_docs[[tf]]
    peak_ids <- sub("^PEAK:", "", terms)
    row_index <- match(peak_ids, rownames(fp_score))
    column_index <- match(one_docs$condition_id, colnames(fp_score))
    if (anyNA(row_index) || anyNA(column_index)) {
      stop("Could not align TF-union Peak scores for TF ", tf, ".")
    }
    score <- fp_score[row_index, column_index, drop = FALSE]
    pieces[[i]] <- data.table::data.table(
      doc_id = rep(one_docs$doc_id, each = length(terms)),
      term_id = rep(terms, times = nrow(one_docs)),
      modality = "peak",
      weight_raw = as.numeric(score)
    )
  }
  output <- data.table::rbindlist(pieces, use.names = TRUE)
  if (any(!is.finite(output$weight_raw) | output$weight_raw <= 0)) {
    stop("TF-union Peak scores must all be positive and finite.")
  }
  output[]
}

build_peak_gene_map <- function(vocabulary, multiomic, docs) {
  if (file.exists(peak_gene_map_path)) {
    return(data.table::fread(peak_gene_map_path, showProgress = FALSE))
  }
  manifest <- data.table::fread(module2_manifest_path, showProgress = FALSE)
  gene_union <- data.table::fread(gene_union_path, showProgress = FALSE)$gene_key
  peak_ids <- sub("^PEAK:", "", vocabulary[startsWith(vocabulary, "PEAK:")])
  gene_ids <- sub("^GENE:", "", vocabulary[startsWith(vocabulary, "GENE:")])
  fp_bound <- readRDS(bound_path)
  gene_expr <- multiomic$matrices$gene_expr
  fp_index_names <- rownames(fp_bound)
  if (is.null(fp_index_names)) fp_index_names <- rownames(multiomic$matrices$fp_score)
  doc_keys <- paste(docs$condition_id, docs$tf, sep = "\r")
  pieces <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    links <- data.table::as.data.table(arrow::read_parquet(
      as.character(manifest$path[[i]]),
      col_select = c("tf", "fp_id", "target_gene", "module2_link_pass"),
      as_data_frame = TRUE
    ))
    links <- links[
      module2_link_pass %in% TRUE & fp_id %in% peak_ids &
        target_gene %in% gene_ids & target_gene %in% gene_union
    ]
    if (!nrow(links)) next
    fp_index <- match(links$fp_id, fp_index_names)
    gene_index <- match(links$target_gene, rownames(gene_expr))
    if (anyNA(fp_index) || anyNA(gene_index)) {
      stop("Could not align a retained Module 2 Peak or target gene.")
    }
    eligible <- logical(nrow(links))
    for (condition in intersect(colnames(fp_bound), colnames(gene_expr))) {
      valid_doc <- paste(condition, links$tf, sep = "\r") %in% doc_keys
      if (!any(valid_doc & !eligible)) next
      condition_index_bound <- match(condition, colnames(fp_bound))
      condition_index_gene <- match(condition, colnames(gene_expr))
      eligible <- eligible |
        (valid_doc & !is.na(fp_index) & !is.na(gene_index) &
          fp_bound[cbind(fp_index, condition_index_bound)] > 0 &
          is.finite(gene_expr[cbind(gene_index, condition_index_gene)]) &
          gene_expr[cbind(gene_index, condition_index_gene)] >= 10)
    }
    pieces[[i]] <- unique(links[eligible, .(
      tf = as.character(tf),
      peak_id = as.character(fp_id),
      gene_key = as.character(target_gene)
    )])
  }
  map <- unique(data.table::rbindlist(pieces, use.names = TRUE, fill = TRUE))
  if (!nrow(map)) stop("No strict p = 1e-10 Peak-to-gene pairs were retained.")
  data.table::setorder(map, peak_id, gene_key, tf)
  data.table::fwrite(map, peak_gene_map_path)
  map[]
}

build_document_supports <- function() {
  check_source_signature()
  outputs <- c(
    support_paths,
    budget_path,
    support_summary_path,
    peak_gene_map_path,
    support_manifest_path
  )
  if (all(file.exists(outputs))) {
    log_info("Reusing validated p = 1e-10 document supports.")
    return(invisible(outputs))
  }
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  log_info("Loading the p = 1e-10 Method 10 documents.")
  baseline <- load_baseline_documents()
  docs <- document_index(baseline)
  sparse <- load_raw_documents(docs$doc_id)
  vocabulary <- sort(unique(sparse$term_id))
  gene_count <- sum(startsWith(vocabulary, "GENE:"))
  peak_count <- sum(startsWith(vocabulary, "PEAK:"))
  if (gene_count != 8085L || peak_count != 21298L) {
    stop(
      "Unexpected strict vocabulary: ", gene_count, " Genes and ",
      peak_count, " Peaks."
    )
  }
  multiomic <- readRDS(multiomic_path)
  factors <- tf_expression_factors(docs, multiomic)

  condition_specific <- sparse[, .(
    doc_id,
    term_id,
    modality,
    weight
  )]
  condition_specific <- apply_peak_factor_and_counts(
    condition_specific,
    factors
  )
  arrow::write_parquet(
    condition_specific,
    support_paths[["condition_specific"]],
    compression = "zstd"
  )

  log_info("Expanding TF-specific Peak unions across existing documents.")
  union_peaks <- build_union_peak_rows(sparse, docs, multiomic)
  union <- data.table::rbindlist(list(
    sparse[modality == "gene", .(
      doc_id,
      term_id,
      modality,
      weight_raw = weight
    )],
    union_peaks
  ), use.names = TRUE)
  if (anyDuplicated(union, by = c("doc_id", "term_id"))) {
    stop("TF-union construction produced duplicate document terms.")
  }
  union <- balance_union_modalities(union)
  union <- apply_peak_factor_and_counts(
    union[, .(doc_id, term_id, modality, weight)],
    factors
  )
  arrow::write_parquet(
    union,
    support_paths[["tf_union"]],
    compression = "zstd"
  )

  base_budget <- baseline[, .(
    baseline_tokens = sum(as.numeric(pseudo_count_log))
  ), by = doc_id]
  gene_terms <- sparse[modality == "gene", .(gene_terms = .N), by = doc_id]
  union_terms_by_tf <- unique(sparse[modality == "peak", .(
    tf = sub("^.*::", "", doc_id),
    term_id
  )])[, .(union_peak_terms = .N), by = tf]
  budget <- merge(docs, base_budget, by = "doc_id")
  budget <- merge(budget, gene_terms, by = "doc_id")
  budget <- merge(budget, union_terms_by_tf, by = "tf")
  budget[, minimum_tokens := pmax(
    5 * union_peak_terms,
    ceiling(1.5 * gene_terms)
  )]
  budget[, common_tokens := pmax(baseline_tokens, minimum_tokens)]
  budget[, tokens_added := common_tokens - baseline_tokens]
  data.table::setorder(budget, doc_id)
  data.table::fwrite(budget, budget_path)

  support_summary <- data.table::rbindlist(lapply(names(support_paths), function(id) {
    value <- if (identical(id, "condition_specific")) condition_specific else union
    value[, .(
      documents = data.table::uniqueN(doc_id),
      document_term_rows = .N,
      gene_document_term_rows = sum(modality == "gene"),
      peak_document_term_rows = sum(modality == "peak"),
      unique_gene_terms = data.table::uniqueN(term_id[modality == "gene"]),
      unique_peak_terms = data.table::uniqueN(term_id[modality == "peak"]),
      source_tokens = sum(source_count)
    )][, support := id]
  }), use.names = TRUE)
  data.table::setcolorder(support_summary, c("support", setdiff(
    names(support_summary),
    "support"
  )))
  data.table::fwrite(support_summary, support_summary_path)

  map <- build_peak_gene_map(vocabulary, multiomic, docs)
  support_manifest <- data.table::data.table(
    artifact = c(
      "condition_specific_support", "tf_union_support", "common_token_budget",
      "peak_gene_map"
    ),
    path = normalizePath(
      c(support_paths, budget_path, peak_gene_map_path),
      winslash = "/",
      mustWork = TRUE
    )
  )
  support_manifest[, `:=`(
    bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path)),
    source_signature = source_signature
  )]
  data.table::fwrite(support_manifest, support_manifest_path)

  peak_by_tf_doc <- union[modality == "peak", .(
    peak_signature = digest::digest(sort(term_id), algo = "xxhash64")
  ), by = .(tf = sub("^.*::", "", doc_id), doc_id)]
  union_check <- peak_by_tf_doc[, .(
    signatures = data.table::uniqueN(peak_signature)
  ), by = tf]
  if (any(union_check$signatures != 1L)) {
    stop("A TF has different Peak support across its union documents.")
  }
  log_info(
    "Built p = 1e-10 supports: ",
    format(nrow(condition_specific), big.mark = ","), " sparse rows and ",
    format(nrow(union), big.mark = ","), " TF-union rows."
  )
  rm(baseline, sparse, condition_specific, union, union_peaks, multiomic)
  invisible(gc())
  invisible(outputs)
}

software_versions <- function() {
  packages <- c("craftgrn", "arrow", "data.table", "Matrix", "yaml")
  package_versions <- stats::setNames(lapply(packages, function(package) {
    as.character(utils::packageVersion(package))
  }), packages)
  commit <- suppressWarnings(system2(
    "git",
    c("-C", shQuote(repository_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ))
  list(
    R = R.version.string,
    platform = R.version$platform,
    packages = package_versions,
    repository_commit = if (length(commit)) trimws(commit[[1L]]) else NA_character_
  )
}

artifact_record <- function(path) {
  if (!file.exists(path)) {
    return(list(path = path, exists = FALSE))
  }
  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    exists = TRUE,
    bytes = as.numeric(file.info(path)$size),
    md5 = unname(tools::md5sum(path))
  )
}

run_configuration <- function(row, status = "configured") {
  support <- as.character(row$support[[1L]])
  gene_to_peak <- as.integer(row$gene_to_peak[[1L]])
  run_root <- as.character(row$run_root[[1L]])
  model_root <- as.character(row$model_root[[1L]])
  extraction_root <- as.character(row$extraction_root[[1L]])
  qc_path <- file.path(benchmark_root, "qc", as.character(row$qc_name[[1L]]))
  metrics_path <- file.path(model_root, "model_metrics.csv")
  recorded_threads <- if (file.exists(metrics_path)) {
    as.integer(data.table::fread(
      metrics_path,
      nrows = 1L,
      showProgress = FALSE
    )$threads[[1L]])
  } else {
    parse_workers()
  }
  list(
    schema_version = 1L,
    plan_date = plan_date,
    run_id = as.character(row$run_id[[1L]]),
    status = status,
    purpose = paste0(
      "HPAFII K30 LDA comparison of condition-specific and TF-union Peak ",
      "support under a sample-specific p = 1e-10 footprint cutoff."
    ),
    project = list(
      name = "nutrient_stress_strict_JASPAR2024_expanded",
      cell_line = "HPAFII",
      local_root = project_root,
      project_yaml = project_config_path,
      analysis_root = analysis_root,
      conditions = 17L
    ),
    source_artifacts = list(
      source_signature = source_signature,
      strict_bound_work_root = strict_root,
      baseline_p1e10_document_terms = baseline_doc_path,
      raw_p1e10_document_root = raw_doc_root,
      multiomic_object = multiomic_path,
      p1e10_condition_bound_matrix = bound_path,
      module2_links_manifest = module2_manifest_path,
      differential_gene_union = gene_union_path,
      condition_manifest = condition_manifest_path,
      shared_support_table = unname(support_paths[[support]]),
      common_token_budget = budget_path,
      peak_gene_map = peak_gene_map_path
    ),
    upstream_filters = list(
      module1_tf_peak_min_r = 0.5,
      module2_tf_target_min_r = 0.5,
      module2_fp_target_min_r = 0.3,
      correlation_statistic = "positive maximum of Pearson and Spearman",
      correlation_p_or_fdr_gate = FALSE,
      footprint_bound_call_required_for_union_entry = TRUE,
      atac_support_required_for_union_entry = TRUE,
      footprint_cutoff = list(
        method = "sample-specific learned null distribution",
        p_value = 1e-10,
        condition_minimum_passing_samples = 1L
      ),
      target_gene_set = "global union with abs(log2FC) >= 1 in at least one comparison",
      target_condition_expression_min = 10,
      tf_condition_expression = "positive finite expression required"
    ),
    document_construction = list(
      document_id = "condition::TF",
      document_count = 7381L,
      gene_terms = list(
        identifier = "GENE:<target_gene>",
        support = "condition-specific eligible target genes",
        weight = "condition expression with condition-gene specificity weighting",
        specificity_temperature = 0.5,
        specificity_floor = 0.1,
        expression_min = 10
      ),
      peak_terms = list(
        identifier = "PEAK:<genomic_coordinate>",
        mode = "unique coordinate",
        support = if (support == "condition_specific") {
          "Only Peaks passing p = 1e-10 in the current condition"
        } else {
          paste0(
            "For each TF, union Peaks passing p = 1e-10 in at least one ",
            "condition; include that union in every existing document for the TF"
          )
        },
        propagated_peak_weight = if (support == "condition_specific") {
          "not applicable"
        } else {
          "positive raw footprint score in the current condition, even below p = 1e-10"
        },
        creates_new_gene_terms = FALSE,
        creates_missing_tf_documents = FALSE,
        weight = "condition footprint score, followed by condition TF-expression multiplier",
        tf_expression_scaling = "per-TF log2(expression + 1) relative to its condition median"
      ),
      modality_balance_before_counting = "equal Gene and Peak weight mass within each document",
      count_conversion = "ceiling(log1p(weight) * 50)",
      minimum_document_frequency = 1L,
      vocabulary = list(
        genes = 8085L,
        peaks = 21298L,
        shared_across_all_models = TRUE
      )
    ),
    token_allocation = list(
      fixed_between_support_designs_at_each_ratio = TRUE,
      budget_source = "p = 1e-10 condition-specific Gene2:Peak1 documents",
      budget_rule = paste0(
        "max(original common budget, (Gene:Peak ratio + 1) * ",
        "TF-union Peak terms)"
      ),
      higher_ratio_expansion = gene_to_peak > 4L,
      gene_to_peak_ratio = paste0(gene_to_peak, ":1"),
      peak_to_gene_ratio_argument = 1 / gene_to_peak,
      minimum_tokens_per_present_term = 1L,
      allocation = "deterministic largest remainder"
    ),
    model = list(
      family = "LDA",
      backend = "WarpLDA OpenMP",
      topics = 30L,
      alpha_sum = 20,
      alpha_per_topic = 20 / 30,
      beta = 1 / 30,
      iterations = 200L,
      convergence_tolerance = 0.001,
      seed = 123L,
      worker_models = 1L,
      threads_per_model = recorded_threads,
      memory_safety = "strict",
      maximum_memory_fraction = 0.8
    ),
    topic_assignment = list(
      score = "normtop_specificity",
      candidate_filter = "GammaFit by topic and term group",
      gammafit_thrP = 0.7,
      gammafit_min_terms = 50L,
      assignment = "gammafit_maxprob",
      unique_peak_assignment = "independent coordinate Peaks",
      tf_retention_filter = paste0(
        "A TF passes for a target topic when at least one existing ",
        "condition::TF document has raw theta >= 0.3 for that topic"
      )
    ),
    outputs = list(
      local_run_root = run_root,
      topic_models = model_root,
      topic_extraction = extraction_root,
      qc_pdf = qc_path,
      comparison_root = benchmark_root,
      box_delivery_root = box_root,
      box_plan = paste0(box_plans_root, "/", box_plan_names[[support]]),
      files = list(
        dtm = artifact_record(file.path(model_root, "rds", "dtm.rds")),
        metrics = artifact_record(file.path(model_root, "model_metrics.csv")),
        theta = artifact_record(file.path(model_root, "vae_models", "theta_K30.csv")),
        phi = artifact_record(file.path(model_root, "vae_models", "phi_K30.csv")),
        qc = artifact_record(qc_path)
      )
    ),
    execution = list(
      script = script_path,
      command = paste(
        "Rscript dev/benchmark/03_benchmark_hpafii_p1e10_peak_support_gene_weight.R",
        stage
      ),
      software = software_versions()
    )
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || all(is.na(x))) y else x
}

write_run_config <- function(row, status = "configured") {
  run_root <- as.character(row$run_root[[1L]])
  dir.create(run_root, recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(
    run_configuration(row, status = status),
    file.path(run_root, "run_config.yaml"),
    indent.mapping.sequence = TRUE
  )
  invisible(file.path(run_root, "run_config.yaml"))
}

write_design_configs <- function(status = "configured") {
  dir.create(benchmark_root, recursive = TRUE, showWarnings = FALSE)
  for (support_id in unique(run_table$support)) {
    rows <- run_table[support == support_id]
    representative <- run_configuration(rows[1L], status)
    shared_model <- representative$model
    shared_model$threads_per_model <- NULL
    shared_model$threads_recorded_per_run <- TRUE
    shared_allocation <- representative$token_allocation
    shared_allocation$higher_ratio_expansion <- NULL
    shared_allocation$gene_to_peak_ratio <- NULL
    shared_allocation$peak_to_gene_ratio_argument <- NULL
    shared_allocation$ratios_above_4_expand_budget_to_keep_all_peaks <- TRUE
    config <- list(
      schema_version = 1L,
      plan_date = plan_date,
      design_id = support_id,
      status = status,
      source_signature = source_signature,
      description = if (support_id == "condition_specific") {
        "Each condition::TF document keeps its own p = 1e-10 passing Peaks."
      } else {
        paste0(
          "Each TF uses its cross-condition union of p = 1e-10 passing Peaks; ",
          "each existing document uses its own raw footprint score."
        )
      },
      source_artifacts = representative$source_artifacts,
      upstream_filters = representative$upstream_filters,
      document_construction = representative$document_construction,
      token_allocation = shared_allocation,
      ratios = lapply(seq_len(nrow(rows)), function(i) list(
        gene_to_peak = paste0(rows$gene_to_peak[[i]], ":1"),
        expanded_budget = rows$gene_to_peak[[i]] > 4L,
        run_id = rows$run_id[[i]],
        threads_per_model = run_configuration(
          rows[i],
          status
        )$model$threads_per_model,
        local_run_root = rows$run_root[[i]],
        qc_pdf = file.path(benchmark_root, "qc", rows$qc_name[[i]])
      )),
      model = shared_model,
      topic_assignment = representative$topic_assignment,
      execution = representative$execution,
      box_delivery_root = box_root,
      box_plan = paste0(box_plans_root, "/", box_plan_names[[support_id]])
    )
    yaml::write_yaml(
      config,
      file.path(benchmark_root, paste0(support_id, "_config.yaml")),
      indent.mapping.sequence = TRUE
    )
  }
  invisible(TRUE)
}

allocate_common_document_budget <- function(dtm, target_budget) {
  output <- numeric(length(dtm@x))
  for (row_index in seq_len(nrow(dtm))) {
    begin <- dtm@p[[row_index]] + 1L
    end <- dtm@p[[row_index + 1L]]
    positions <- seq.int(begin, end)
    term_count <- length(positions)
    target <- as.integer(target_budget[[row_index]])
    if (target < term_count) {
      stop("The common token budget cannot retain every document term.")
    }
    score <- dtm@x[positions]
    remaining <- target - term_count
    quota <- score / sum(score) * remaining
    fractional <- quota - floor(quota)
    cumulative <- cumsum(fractional)
    residual <- floor(cumulative + 1e-8) - c(
      0,
      head(floor(cumulative + 1e-8), -1L)
    )
    counts <- 1 + floor(quota) + residual
    correction <- target - sum(counts)
    if (correction != 0) counts[[which.max(counts)]] <- counts[[which.max(counts)]] + correction
    if (any(counts < 1) || sum(counts) != target) {
      stop("Common document-budget allocation failed.")
    }
    output[positions] <- counts
  }
  dtm@x <- output
  dtm
}

make_source_dtm <- function(support_id, documents, vocabulary, budget) {
  support <- data.table::as.data.table(arrow::read_parquet(
    support_paths[[support_id]],
    col_select = c("doc_id", "term_id", "source_count"),
    as_data_frame = TRUE
  ))
  document_index <- match(support$doc_id, documents)
  term_index <- match(support$term_id, vocabulary)
  if (anyNA(document_index) || anyNA(term_index)) {
    stop("Support rows do not align to the common DTM index.")
  }
  dtm <- Matrix::sparseMatrix(
    i = document_index,
    j = term_index,
    x = as.numeric(support$source_count),
    dims = c(length(documents), length(vocabulary)),
    dimnames = list(documents, vocabulary),
    repr = "R"
  )
  dtm <- methods::as(dtm, "RsparseMatrix")
  source_totals <- as.numeric(Matrix::rowSums(dtm))
  budget_index <- match(documents, budget$doc_id)
  common_budget <- budget$common_tokens[budget_index]
  union_peak_terms <- budget$union_peak_terms[budget_index]
  if (anyNA(common_budget) || anyNA(union_peak_terms) ||
      any(source_totals <= 0)) {
    stop("A document is missing its source count or common token budget.")
  }
  list(
    dtm = dtm,
    support_rows = nrow(support),
    common_budget = common_budget,
    union_peak_terms = union_peak_terms
  )
}

write_document_qc <- function(dtm, row, overwrite = FALSE) {
  model_root <- as.character(row$model_root[[1L]])
  qc_path <- file.path(model_root, "document_term_qc.pdf")
  summary_path <- file.path(model_root, "document_term_qc_summary.csv")
  if (!overwrite && file.exists(qc_path) && file.exists(summary_path)) {
    return(invisible(qc_path))
  }
  entries <- Matrix::summary(dtm)
  doc_term <- data.table::data.table(
    doc_id = rownames(dtm)[entries$i],
    term_id = colnames(dtm)[entries$j],
    pseudo_count_log = as.numeric(entries$x)
  )
  craftgrn:::.write_module3_document_term_qc(
    doc_term = doc_term,
    output_dir = model_root,
    count_column = "pseudo_count_log",
    title = paste0(
      "HPAFII p = 1e-10 | ", row$support_label[[1L]],
      " | Gene:Peak ", row$gene_to_peak[[1L]], ":1"
    ),
    document_unit = "tf",
    verbose = TRUE
  )
  rm(entries, doc_term)
  invisible(gc())
  invisible(qc_path)
}

prepare_design_matrices <- function(support_id) {
  build_document_supports()
  rows <- run_table[support == support_id]
  dtm_is_current <- function(row) {
    dtm_path <- file.path(row$model_root[[1L]], "rds", "dtm.rds")
    if (!file.exists(dtm_path)) return(FALSE)
    if (row$gene_to_peak[[1L]] <= 4L) return(TRUE)
    summary_path <- file.path(row$model_root[[1L]], "topic_input_summary.csv")
    if (!file.exists(summary_path)) return(FALSE)
    summary <- data.table::fread(summary_path, nrows = 1L, showProgress = FALSE)
    qc_paths <- file.path(
      row$model_root[[1L]],
      c("document_term_qc.pdf", "document_term_qc_summary.csv")
    )
    "ratio_budget_version" %in% names(summary) &&
      identical(as.integer(summary$ratio_budget_version[[1L]]), 2L) &&
      all(file.exists(qc_paths)) &&
      all(file.info(qc_paths)$mtime >= file.info(dtm_path)$mtime)
  }
  current <- vapply(
    seq_len(nrow(rows)),
    function(i) dtm_is_current(rows[i]),
    logical(1L)
  )
  if (all(current)) {
    log_info("Reusing all DTM files for ", support_id, ".")
    return(invisible(TRUE))
  }
  condition_support <- data.table::as.data.table(arrow::read_parquet(
    support_paths[["condition_specific"]],
    col_select = c("doc_id", "term_id"),
    as_data_frame = TRUE
  ))
  documents <- sort(unique(condition_support$doc_id))
  vocabulary <- c(
    sort(unique(condition_support$term_id[startsWith(
      condition_support$term_id,
      "GENE:"
    )])),
    sort(unique(condition_support$term_id[startsWith(
      condition_support$term_id,
      "PEAK:"
    )]))
  )
  budget <- data.table::fread(budget_path, showProgress = FALSE)
  source <- make_source_dtm(support_id, documents, vocabulary, budget)
  gene_term <- startsWith(vocabulary, "GENE:")
  for (i in seq_len(nrow(rows))) {
    row <- rows[i]
    model_root <- as.character(row$model_root[[1L]])
    dtm_path <- file.path(model_root, "rds", "dtm.rds")
    if (current[[i]]) next
    dir.create(file.path(model_root, "rds"), recursive = TRUE, showWarnings = FALSE)
    log_info(
      "Allocating ", support_id, " DTM at Gene:Peak ",
      row$gene_to_peak[[1L]], ":1."
    )
    target_budget <- pmax(
      source$common_budget,
      (as.integer(row$gene_to_peak[[1L]]) + 1L) * source$union_peak_terms
    )
    dtm <- allocate_common_document_budget(source$dtm, target_budget)
    finalized <- craftgrn:::.topic_finalize_sparse_counts_cpp(
      row_pointer = dtm@p,
      column_index = dtm@j,
      source_count = dtm@x,
      gene_term = gene_term,
      idf_multiplier = rep(1, length(vocabulary)),
      peak_gene_ratio = 1 / row$gene_to_peak[[1L]]
    )
    dtm@x <- as.numeric(finalized$counts)
    dtm <- methods::as(dtm, "CsparseMatrix")
    observed_budget <- as.numeric(Matrix::rowSums(dtm))
    if (!identical(observed_budget, as.numeric(target_budget))) {
      stop("Final DTM did not preserve the ratio-specific document budget.")
    }
    peak_tokens <- as.numeric(Matrix::rowSums(dtm[, !gene_term, drop = FALSE]))
    peak_gene_ratio <- 1 / row$gene_to_peak[[1L]]
    expected_peak <- floor(
      observed_budget * peak_gene_ratio / (1 + peak_gene_ratio) + 0.5
    )
    if (any(peak_tokens != expected_peak)) {
      stop("Final DTM did not realize the requested Gene:Peak ratio.")
    }
    saveRDS(dtm, dtm_path, compress = TRUE)
    input_summary <- data.table::data.table(
      analysis_label = paste0(
        "HPAFII_p1e10_", support_id, "_K30_Gene",
        row$gene_to_peak[[1L]], "_Peak1"
      ),
      source_signature = source_signature,
      support = support_id,
      gene_peak_ratio = paste0(row$gene_to_peak[[1L]], ":1"),
      documents = nrow(dtm),
      gene_terms = sum(gene_term),
      peak_terms = sum(!gene_term),
      nonzero_document_terms = length(dtm@x),
      model_tokens = sum(dtm@x),
      common_budget_added_tokens = sum(
        budget$tokens_added[match(rownames(dtm), budget$doc_id)]
      ),
      ratio_budget_added_tokens = sum(target_budget - source$common_budget),
      ratio_budget_version = 2L,
      count_method = "log",
      count_scale = 50,
      min_df = 1L
    )
    data.table::fwrite(
      input_summary,
      file.path(model_root, "topic_input_summary.csv")
    )
    file.copy(
      peak_gene_map_path,
      file.path(model_root, "peak_gene_map.csv"),
      overwrite = TRUE,
      copy.mode = TRUE
    )
    write_run_config(row, status = "documents_complete")
    write_document_qc(dtm, row, overwrite = TRUE)
    rm(dtm, finalized)
    invisible(gc())
  }
  rm(source, condition_support)
  invisible(gc())
  invisible(TRUE)
}

prepare_documents <- function() {
  build_document_supports()
  prepare_design_matrices("condition_specific")
  prepare_design_matrices("tf_union")
  write_design_configs(status = "documents_complete")
  invisible(TRUE)
}

train_one_model <- function(row) {
  model_root <- as.character(row$model_root[[1L]])
  phi_path <- file.path(model_root, "vae_models", "phi_K30.csv")
  theta_path <- file.path(model_root, "vae_models", "theta_K30.csv")
  if (file.exists(phi_path) && file.exists(theta_path)) {
    log_info("Reusing trained model: ", row$run_id[[1L]])
    return(invisible(TRUE))
  }
  dtm_path <- file.path(model_root, "rds", "dtm.rds")
  if (!file.exists(dtm_path)) stop("Missing DTM: ", dtm_path)
  dtm <- readRDS(dtm_path)
  log_info("Training K30 WarpLDA model: ", row$run_id[[1L]])
  fits <- craftgrn:::run_warplda_models(
    dtm,
    K_grid = 30L,
    iterations = 200L,
    alpha_by_topic = TRUE,
    alpha_sum = 20,
    beta = 1 / 30,
    seed = 123L,
    save_tmp_dir = file.path(model_root, "tmp_models"),
    workers = 1L,
    threads_per_model = parse_workers(),
    sampler = "warp_omp",
    metrics_file = file.path(model_root, "model_metrics.csv"),
    memory_safety = "strict",
    memory_max_fraction = 0.8
  )
  fit <- readRDS(fits$fit_files[[1L]])
  dir.create(file.path(model_root, "vae_models"), recursive = TRUE, showWarnings = FALSE)
  write_probability(fit$theta, "doc_id", theta_path)
  write_probability(fit$phi, "term_id", phi_path)
  metrics <- data.table::as.data.table(fits$metrics)
  metrics[, `:=`(
    source_signature = source_signature,
    support = row$support[[1L]],
    gene_peak_ratio = paste0(row$gene_to_peak[[1L]], ":1"),
    run_id = row$run_id[[1L]]
  )]
  data.table::fwrite(metrics, file.path(model_root, "model_metrics.csv"))
  selection <- craftgrn:::plot_model_selection_cistopic(
    metrics,
    file.path(model_root, "model_selection.pdf"),
    title_prefix = paste0(
      "HPAFII p = 1e-10 | ", row$support_label[[1L]],
      " | Gene:Peak ", row$gene_to_peak[[1L]], ":1"
    )
  )
  saveRDS(selection, file.path(model_root, "rds", "model_selection.rds"))
  write_run_config(row, status = "model_complete")
  rm(dtm, fits, fit, selection)
  invisible(gc())
  invisible(TRUE)
}

train_models <- function() {
  prepare_documents()
  pending <- which(!file.exists(file.path(
    run_table$model_root,
    "vae_models",
    "phi_K30.csv"
  )))
  jobs <- min(parse_concurrent_models(), length(pending))
  if (jobs > 1L && .Platform$OS.type != "windows") {
    log_info(
      "Training ", length(pending), " remaining models with ", jobs,
      " bounded jobs and ", parse_workers(), " threads per model."
    )
    result <- parallel::mclapply(
      pending,
      function(i) train_one_model(run_table[i]),
      mc.cores = jobs,
      mc.preschedule = FALSE
    )
    failed <- vapply(result, inherits, logical(1L), what = "try-error")
    if (any(failed)) {
      stop("At least one bounded K30 training job failed: ", paste(
        vapply(result[failed], as.character, character(1L)),
        collapse = " | "
      ))
    }
  } else {
    for (i in pending) train_one_model(run_table[i])
  }
  write_design_configs(status = "models_complete")
  invisible(TRUE)
}

make_retention_assignments <- function(theta, pair_assignment, peak_map) {
  pairs <- data.table::copy(data.table::as.data.table(pair_assignment))
  pairs[, pair_index := seq_len(.N)]
  pair_keys <- pairs[, .(
    peak_id = sub("^PEAK:", "", peak_term_id),
    gene_key = sub("^GENE:", "", gene_term_id),
    pair_index
  )]
  links <- data.table::copy(peak_map)
  links[pair_keys, pair_index := i.pair_index, on = c("peak_id", "gene_key")]
  if (anyNA(links$pair_index)) {
    stop("A TF-Peak-Gene record is missing its topic-assignment pair.")
  }

  doc_tf <- sub("^.*::", "", rownames(theta))
  tf_levels <- sort(unique(links$tf))
  topic_ids <- craftgrn:::.m3_opt_topic_ids(theta, "column")
  tf_topic_max <- matrix(
    0,
    nrow = length(tf_levels),
    ncol = length(topic_ids),
    dimnames = list(tf_levels, as.character(topic_ids))
  )
  for (tf in intersect(tf_levels, unique(doc_tf))) {
    tf_topic_max[tf, ] <- apply(theta[doc_tf == tf, , drop = FALSE], 2L, max)
  }
  assigned_topic <- suppressWarnings(as.integer(
    pairs$assigned_topic[links$pair_index]
  ))
  topic_position <- match(assigned_topic, topic_ids)
  valid <- is.finite(topic_position)
  raw_aligned <- rep(FALSE, nrow(links))
  raw_aligned[valid] <- tf_topic_max[cbind(
    match(links$tf[valid], tf_levels),
    topic_position[valid]
  )] >= 0.3
  data.table::data.table(
    tf_index = match(links$tf, tf_levels),
    target_index = match(links$gene_key, sort(unique(links$gene_key))),
    pair_index = as.integer(links$pair_index),
    raw_aligned = raw_aligned
  )
}

make_identity_qc <- function(theta,
                             phi,
                             topic_terms,
                             pair_assignment,
                             peak_map) {
  gene_ids <- which(startsWith(colnames(phi), "GENE:"))
  peak_ids <- which(startsWith(colnames(phi), "PEAK:"))
  similarity <- craftgrn:::.m3_opt_hellinger_similarity(
    phi,
    gene_ids,
    peak_ids
  )
  diag(similarity) <- 1
  counts <- pair_assignment[
    craftgrn:::.as_logical_flag(assigned) & is.finite(assigned_topic),
    .(
      peaks = data.table::uniqueN(peak_term_id),
      genes = data.table::uniqueN(gene_term_id)
    ),
    by = .(raw_topic = as.integer(assigned_topic))
  ]
  topics <- craftgrn:::.m3_opt_topic_ids(phi, "row")
  counts <- counts[data.table::data.table(raw_topic = topics), on = "raw_topic"]
  counts[is.na(peaks), peaks := 0L]
  counts[is.na(genes), genes := 0L]
  list(
    document_design = "condition_tf",
    assignment_mode = "gene_peak",
    theta = theta,
    raw_phi = phi,
    raw_topic_terms = topic_terms,
    raw_pair_assignment = pair_assignment,
    raw_to_optimized = stats::setNames(topics, topics),
    qc = list(
      raw_topic_ids = topics,
      optimized_topic_ids = topics,
      raw_topic_similarity = similarity,
      raw_counts = counts,
      assignments = make_retention_assignments(
        theta,
        pair_assignment,
        peak_map
      )
    )
  )
}

mean_off_diagonal <- function(value) {
  if (!is.matrix(value) || nrow(value) < 2L || ncol(value) < 2L) {
    return(NA_real_)
  }
  entries <- value[row(value) != col(value)]
  entries <- entries[is.finite(entries)]
  if (length(entries)) mean(entries) else NA_real_
}

write_qc_one <- function(row) {
  model_root <- as.character(row$model_root[[1L]])
  extraction_root <- as.character(row$extraction_root[[1L]])
  qc_dir <- file.path(benchmark_root, "qc")
  qc_path <- file.path(qc_dir, as.character(row$qc_name[[1L]]))
  metrics_path <- file.path(extraction_root, "benchmark_metrics.csv")
  existing_pages <- if (file.exists(qc_path)) {
    suppressWarnings(as.integer(system2(
      "qpdf",
      c("--show-npages", shQuote(qc_path)),
      stdout = TRUE,
      stderr = TRUE
    )[[1L]]))
  } else {
    NA_integer_
  }
  if (identical(existing_pages, 7L) && file.exists(metrics_path)) {
    log_info("Reusing QC report: ", row$qc_name[[1L]])
    return(data.table::fread(metrics_path, showProgress = FALSE))
  }
  theta <- read_probability(file.path(
    model_root,
    "vae_models",
    "theta_K30.csv"
  ))
  phi <- read_probability(file.path(
    model_root,
    "vae_models",
    "phi_K30.csv"
  ))
  if (ncol(theta) != 30L || nrow(phi) != 30L) {
    stop("Expected a K30 theta and phi matrix for ", row$run_id[[1L]], ".")
  }
  score <- craftgrn:::score_terms_normtop(
    phi,
    method = "normtop_specificity"
  )
  candidates <- craftgrn:::binarize_topics(
    score,
    method = "gammafit",
    thrP = 0.7,
    min_terms = 50L,
    gammafit_scope = "topic_term_group"
  )
  terms <- craftgrn:::.assign_topic_terms(
    raw_phi = phi,
    score_mat = score,
    candidate_terms = candidates,
    method = "gammafit_maxprob",
    independent_unique_peaks = TRUE
  )
  peak_map <- data.table::fread(peak_gene_map_path, showProgress = FALSE)
  pairs <- craftgrn:::.topic_unique_peak_gene_assignment_table(
    terms,
    peak_map[, .(peak_id, gene_key)]
  )
  genes <- craftgrn:::.topic_gene_assignment_table(terms)
  diagnostics <- craftgrn:::.gammafit_diagnostics_by_termclass(
    score,
    topic_terms = candidates,
    topic_score_method = "normtop_specificity",
    thrP = 0.7,
    min_terms = 50L,
    gammafit_scope = "topic_term_group"
  )
  assignment_summary <- craftgrn:::.topic_gene_peak_assignment_summary(
    pairs,
    thrP = 0.7,
    assignment_method = "gammafit_maxprob_coordinate_peak",
    model_family = "lda",
    expected_topic_count = 30L
  )
  optimization <- make_identity_qc(theta, phi, terms, pairs, peak_map)

  dir.create(extraction_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(extraction_root, "rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(terms, file.path(extraction_root, "topic_terms.csv"))
  data.table::fwrite(
    pairs,
    file.path(extraction_root, "topic_coordinate_peak_gene_assignment.csv")
  )
  data.table::fwrite(
    genes,
    file.path(extraction_root, "topic_gene_assignment.csv")
  )
  data.table::fwrite(
    diagnostics,
    file.path(extraction_root, "topic_gammafit_diagnostics.csv")
  )
  data.table::fwrite(
    assignment_summary,
    file.path(extraction_root, "topic_term_assignment_summary.csv")
  )
  saveRDS(score, file.path(
    extraction_root,
    "rds",
    "topic_term_scores_normtop_specificity.rds"
  ))

  config <- yaml::read_yaml(project_config_path)
  tf_sets <- if (file.exists(reference_tf_sets_path)) {
    readRDS(reference_tf_sets_path)
  } else {
    NULL
  }
  tf_panel <- if (file.exists(reference_tf_panel_path)) {
    data.table::fread(reference_tf_panel_path, showProgress = FALSE)
  } else {
    NULL
  }
  title <- paste0(
    "HPAFII | K30 | ", row$support_label[[1L]],
    " | Gene:Peak ", row$gene_to_peak[[1L]], ":1 | p = 1e-10"
  )
  log_info("Writing QC report: ", row$qc_name[[1L]])
  craftgrn:::.write_module3_topic_assignment_qc(
    optimization,
    out_file = qc_path,
    title_prefix = title,
    condition_colors = craftgrn:::.module3_topic_condition_colors(config),
    gene_term_assignment = genes,
    gene_umap_genes = genes$target_gene,
    peak_term_assignment = terms,
    tf_target_gene_sets = tf_sets,
    tf_target_panel = tf_panel,
    top_n_tfs = 25L,
    seed = 20260902L + as.integer(row$run_number[[1L]]),
    sections = "standard_hellinger",
    peak_umap_top_n = Inf,
    report_scope = "condition_correlation",
    sidebar_mode = "terms"
  )
  extraction_qc <- file.path(extraction_root, "topic_assignment_qc_K30.pdf")
  if (!file.copy(qc_path, extraction_qc, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not copy the QC report into the extraction folder.")
  }
  pages <- suppressWarnings(as.integer(system2(
    "qpdf",
    c("--show-npages", shQuote(qc_path)),
    stdout = TRUE,
    stderr = TRUE
  )[[1L]]))
  if (!identical(pages, 7L)) {
    stop("Expected seven QC pages but found ", pages, " in ", qc_path, ".")
  }

  condition_theta <- craftgrn:::.m3_qc_condition_tf_topic_theta(
    theta,
    matched_tfs = TRUE
  )
  condition_correlation <- craftgrn:::.m3_qc_theta_condition_correlation(
    condition_theta
  )
  assigned_pairs <- pairs[
    craftgrn:::.as_logical_flag(assigned) & is.finite(assigned_topic)
  ]
  topic_gene_counts <- assigned_pairs[, .(
    genes = data.table::uniqueN(target_gene)
  ), by = .(topic = as.integer(assigned_topic))]
  assigned_term_rows <- terms[craftgrn:::.as_logical_flag(in_topic)]
  model_metrics <- data.table::fread(
    file.path(model_root, "model_metrics.csv"),
    showProgress = FALSE
  )[K == 30L]
  metrics <- data.table::data.table(
    run_id = as.character(row$run_id[[1L]]),
    support = as.character(row$support[[1L]]),
    support_label = as.character(row$support_label[[1L]]),
    gene_to_peak = as.integer(row$gene_to_peak[[1L]]),
    K = 30L,
    documents = nrow(theta),
    gene_terms = sum(startsWith(colnames(phi), "GENE:")),
    peak_terms = sum(startsWith(colnames(phi), "PEAK:")),
    assigned_target_genes = data.table::uniqueN(assigned_pairs$target_gene),
    assigned_peak_terms = data.table::uniqueN(
      assigned_term_rows[startsWith(term_id, "PEAK:"), term_id]
    ),
    assigned_topics = data.table::uniqueN(topic_gene_counts$topic),
    largest_topic_gene_fraction = if (nrow(topic_gene_counts)) {
      max(topic_gene_counts$genes) / sum(topic_gene_counts$genes)
    } else {
      NA_real_
    },
    mean_topic_similarity = mean_off_diagonal(
      optimization$qc$raw_topic_similarity
    ),
    mean_condition_correlation = mean_off_diagonal(condition_correlation),
    perplexity = as.numeric(model_metrics$perplexity[[1L]]),
    log_likelihood = as.numeric(model_metrics$loglik[[1L]]),
    tokens = as.numeric(model_metrics$n_tokens[[1L]]),
    iterations_completed = as.integer(model_metrics$iterations_completed[[1L]]),
    converged = isTRUE(model_metrics$converged[[1L]]),
    qc_pages = pages,
    qc_pdf = normalizePath(qc_path, winslash = "/", mustWork = TRUE)
  )
  data.table::fwrite(metrics, metrics_path)
  write_run_config(row, status = "qc_complete")
  rm(
    theta, phi, score, candidates, terms, peak_map, pairs, genes,
    diagnostics, optimization, condition_theta, condition_correlation,
    assigned_pairs, topic_gene_counts, assigned_term_rows
  )
  invisible(gc())
  metrics
}

write_qc_reports <- function() {
  train_models()
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Package devtools is required to load the current QC report code.")
  }
  devtools::load_all(
    repository_root,
    quiet = TRUE,
    export_all = FALSE,
    helpers = FALSE
  )
  report_arguments <- names(formals(craftgrn:::.write_module3_topic_assignment_qc))
  if (!all(c("peak_term_assignment", "peak_umap_top_n") %in% report_arguments)) {
    stop("The current CraftGRN source does not provide the required Peak QC pages.")
  }
  metrics <- data.table::rbindlist(lapply(
    seq_len(nrow(run_table)),
    function(i) write_qc_one(run_table[i])
  ), use.names = TRUE, fill = TRUE)
  data.table::fwrite(metrics, file.path(benchmark_root, "benchmark_metrics.csv"))
  write_design_configs(status = "qc_complete")
  invisible(metrics)
}

normalize_topic_rows <- function(phi, term_index) {
  value <- phi[, term_index, drop = FALSE]
  totals <- rowSums(value)
  if (any(!is.finite(totals) | totals <= 0)) {
    stop("A topic has no finite Gene probability mass.")
  }
  value / totals
}

topic_alignment_one <- function(gene_to_peak) {
  ratio_value <- as.integer(gene_to_peak)
  rows <- run_table[run_table$gene_to_peak == ratio_value]
  if (nrow(rows) != 2L) stop("Expected two support designs per ratio.")
  sparse_row <- rows[support == "condition_specific"]
  union_row <- rows[support == "tf_union"]
  sparse_phi <- read_probability(file.path(
    sparse_row$model_root[[1L]],
    "vae_models",
    "phi_K30.csv"
  ))
  union_phi <- read_probability(file.path(
    union_row$model_root[[1L]],
    "vae_models",
    "phi_K30.csv"
  ))
  if (!identical(colnames(sparse_phi), colnames(union_phi))) {
    stop("The two designs do not share an identical vocabulary and order.")
  }
  gene_index <- which(startsWith(colnames(sparse_phi), "GENE:"))
  sparse_gene <- normalize_topic_rows(sparse_phi, gene_index)
  union_gene <- normalize_topic_rows(union_phi, gene_index)
  similarity <- tcrossprod(sqrt(sparse_gene), sqrt(union_gene))
  sparse_best <- max.col(similarity, ties.method = "first")
  union_best <- max.col(t(similarity), ties.method = "first")
  data.table::rbindlist(list(
    data.table::data.table(
      gene_to_peak = ratio_value,
      source_support = "condition_specific",
      source_topic = seq_len(nrow(similarity)),
      best_matching_support = "tf_union",
      best_matching_topic = sparse_best,
      hellinger_similarity = similarity[cbind(
        seq_len(nrow(similarity)),
        sparse_best
      )]
    ),
    data.table::data.table(
      gene_to_peak = ratio_value,
      source_support = "tf_union",
      source_topic = seq_len(ncol(similarity)),
      best_matching_support = "condition_specific",
      best_matching_topic = union_best,
      hellinger_similarity = t(similarity)[cbind(
        seq_len(ncol(similarity)),
        union_best
      )]
    )
  ))
}

comparison_theme <- function(base_size = 10) {
  craftgrn:::.m3_qc_theme(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(
        family = "Helvetica",
        color = "#111111"
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(size = base_size),
      legend.position = "bottom"
    )
}

write_comparison_report <- function() {
  metrics_path <- file.path(benchmark_root, "benchmark_metrics.csv")
  if (!file.exists(metrics_path)) write_qc_reports()
  metrics <- data.table::fread(metrics_path, showProgress = FALSE)
  if (nrow(metrics) != nrow(run_table) || any(!metrics$converged) ||
      any(metrics$qc_pages != 7L)) {
    stop("Benchmark metrics do not describe all converged seven-page runs.")
  }
  support <- data.table::fread(support_summary_path, showProgress = FALSE)
  alignment <- data.table::rbindlist(lapply(
    sort(unique(run_table$gene_to_peak)),
    topic_alignment_one
  ))
  data.table::fwrite(
    alignment,
    file.path(benchmark_root, "topic_alignment.csv")
  )
  alignment_summary <- alignment[, .(
    mean_best_topic_similarity = mean(hellinger_similarity),
    median_best_topic_similarity = stats::median(hellinger_similarity),
    minimum_best_topic_similarity = min(hellinger_similarity)
  ), by = gene_to_peak]
  comparison <- merge(
    metrics,
    alignment_summary,
    by = "gene_to_peak",
    all.x = TRUE,
    sort = FALSE
  )
  comparison <- merge(
    comparison,
    support[, .(
      support,
      document_term_rows,
      gene_document_term_rows,
      peak_document_term_rows
    )],
    by = "support",
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setorder(comparison, support, gene_to_peak)
  comparison_path <- file.path(benchmark_root, "comparison_metrics.csv")
  data.table::fwrite(comparison, comparison_path)

  colors <- c(
    "Condition-specific Peaks" = "#0072B2",
    "TF-union Peaks" = "#D55E00"
  )
  support_plot <- ggplot2::ggplot(
    unique(comparison[, .(
      support_label,
      peak_document_term_rows
    )]),
    ggplot2::aes(support_label, peak_document_term_rows, fill = support_label)
  ) +
    ggplot2::geom_col(width = 0.62, color = "#222222", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(peak_document_term_rows)),
      vjust = -0.4,
      size = 4
    ) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(0, 0.16))
    ) +
    ggplot2::labs(
      title = "Peak support in the documents",
      subtitle = "The TF-union design keeps the same Peak set across conditions for each TF.",
      x = NULL,
      y = "Peak entries"
    ) +
    comparison_theme()

  fit_plot <- ggplot2::ggplot(
    comparison,
    ggplot2::aes(
      factor(gene_to_peak),
      perplexity,
      color = support_label,
      group = support_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(
      title = "Model fit",
      subtitle = "Lower perplexity means the model explains its own documents better.",
      x = "Gene:Peak ratio",
      y = "Perplexity",
      color = NULL
    ) +
    comparison_theme()

  page_one <- support_plot + fit_plot +
    patchwork::plot_layout(widths = c(1, 1.15)) +
    patchwork::plot_annotation(
      title = "HPAFII p = 1e-10 Peak-support benchmark",
      subtitle = "Both Peak-support designs use the same token budget at each ratio."
    )

  metric_data <- data.table::melt(
    comparison[, c(
      "assigned_target_genes",
      "largest_topic_gene_fraction",
      "mean_condition_correlation",
      "mean_best_topic_similarity"
    ) := lapply(.SD, as.numeric), .SDcols = c(
      "assigned_target_genes",
      "largest_topic_gene_fraction",
      "mean_condition_correlation",
      "mean_best_topic_similarity"
    )],
    id.vars = c("support_label", "gene_to_peak"),
    measure.vars = c(
      "assigned_target_genes",
      "largest_topic_gene_fraction",
      "mean_condition_correlation",
      "mean_best_topic_similarity"
    ),
    variable.name = "metric",
    value.name = "value"
  )
  labels <- c(
    assigned_target_genes = "Assigned Genes",
    largest_topic_gene_fraction = "Largest Topic share",
    mean_condition_correlation = "Condition similarity",
    mean_best_topic_similarity = "Cross-design Topic match"
  )
  metric_data[, metric := factor(metric, levels = names(labels), labels = labels)]
  metric_plot <- ggplot2::ggplot(
    metric_data,
    ggplot2::aes(
      factor(gene_to_peak),
      value,
      fill = support_label
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge2(width = 0.78, preserve = "single"),
      width = 0.68,
      color = "#222222",
      linewidth = 0.18
    ) +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", ncol = 2L) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(
      title = "Topic results",
      subtitle = "Cross-design Topic match is the best Gene-profile match between the two support rules.",
      x = "Gene:Peak ratio",
      y = NULL,
      fill = NULL
    ) +
    comparison_theme()

  comparison_pdf <- file.path(
    benchmark_root,
    "HPAFII_K30_PeakSupport_GeneWeight_Comparison.pdf"
  )
  grDevices::cairo_pdf(comparison_pdf, width = 12, height = 7.5)
  print(page_one)
  print(metric_plot)
  grDevices::dev.off()
  pages <- as.integer(system2(
    "qpdf",
    c("--show-npages", shQuote(comparison_pdf)),
    stdout = TRUE
  )[[1L]])
  if (!identical(pages, 2L)) stop("Comparison report must have two pages.")
  write_design_configs(status = "comparison_complete")
  log_info("Wrote the two-page benchmark comparison report.")
  invisible(comparison)
}

write_delivery_readme <- function() {
  text <- c(
    "# HPAFII K30 Peak-support benchmark",
    "",
    "This folder compares two ways to keep Peak terms in condition::TF documents.",
    "",
    "- Condition-specific: each document keeps only Peaks that pass p = 1e-10 in that condition.",
    "- TF-union: each TF keeps its union of passing Peaks across conditions. Each condition uses its own raw footprint score.",
    "- Gene:Peak ratios tested: 2:1, 3:1, 4:1, 6:1, 8:1, and 10:1.",
    "- Both designs use the same document budget at each ratio and K = 30.",
    "- Ratios above 4:1 use a larger budget so every Peak stays in the document.",
    "",
    "Start with `HPAFII_K30_PeakSupport_GeneWeight_Comparison.pdf`.",
    "Each QC PDF adds matched-TF Hellinger separation as page 7.",
    "The two YAML files contain the full input, filter, document, model, and output records.",
    "",
    paste0("Source signature: `", source_signature, "`")
  )
  path <- file.path(benchmark_root, "README.md")
  writeLines(text, path, useBytes = TRUE)
  path
}

delivery_file_table <- function() {
  qc_paths <- file.path(benchmark_root, "qc", run_table$qc_name)
  paths <- c(
    file.path(benchmark_root, "HPAFII_K30_PeakSupport_GeneWeight_Comparison.pdf"),
    file.path(benchmark_root, "comparison_metrics.csv"),
    file.path(benchmark_root, "topic_alignment.csv"),
    file.path(benchmark_root, "condition_specific_config.yaml"),
    file.path(benchmark_root, "tf_union_config.yaml"),
    file.path(benchmark_root, "README.md"),
    qc_paths
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing delivery file(s): ", paste(missing, collapse = ", "))
  }
  data.table::data.table(
    file = basename(paths),
    source = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    md5 = unname(tools::md5sum(paths))
  )
}

update_box_parent_index <- function() {
  update_root <- file.path(benchmark_root, "parent_index_update")
  dir.create(update_root, recursive = TRUE, showWarnings = FALSE)
  remote_index <- paste0(box_parent, "/INPUT_DESIGN_INDEX.csv")
  remote_readme <- paste0(box_parent, "/README.md")
  local_index <- file.path(update_root, "INPUT_DESIGN_INDEX.csv")
  local_readme <- file.path(update_root, "README.md")
  status_index <- system2(
    "rclone",
    c("copyto", shQuote(remote_index), shQuote(local_index)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!identical(attr(status_index, "status") %||% 0L, 0L) ||
      !file.exists(local_index)) {
    stop("Could not download the Box input-design index.")
  }
  status_readme <- system2(
    "rclone",
    c("copyto", shQuote(remote_readme), shQuote(local_readme)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!identical(attr(status_readme, "status") %||% 0L, 0L) ||
      !file.exists(local_readme)) {
    stop("Could not download the Box input-design README.")
  }

  index <- data.table::fread(local_index, showProgress = FALSE)
  new_files <- file.path(
    "16_Method10_P1e10_PeakSupport_GeneWeight_Benchmark",
    run_table$qc_name
  )
  index <- index[!file %in% new_files]
  new_rows <- data.table::data.table(
    order = max(as.integer(index$order), na.rm = TRUE) + seq_len(nrow(run_table)),
    file = new_files,
    source_run = run_table$run_id,
    model = "LDA",
    K = 30L,
    document_unit = "condition::TF (7381 documents)",
    input_terms = data.table::fifelse(
      run_table$support == "condition_specific",
      "Gene expression + condition-specific coordinate Peaks",
      "Gene expression + per-TF union coordinate Peaks"
    ),
    gene_weight = "Condition expression + specificity",
    peak_weight = "Condition footprint score + TF expression",
    gene_peak_ratio = paste0(run_table$gene_to_peak, ":1"),
    status = "active_experimental",
    pages = 7L,
    role = "p = 1e-10 Peak-support and Gene-weight benchmark"
  )
  new_rows <- new_rows[, intersect(names(index), names(new_rows)), with = FALSE]
  index <- data.table::rbindlist(list(index, new_rows), use.names = TRUE, fill = TRUE)
  data.table::setorder(index, order)
  data.table::fwrite(index, local_index)

  readme <- readLines(local_readme, warn = FALSE)
  marker <- "## Method 10 p = 1e-10 Peak-support benchmark"
  benchmark_note <- paste0(
    "See `16_Method10_P1e10_PeakSupport_GeneWeight_Benchmark/` for ",
    "the condition-specific and TF-union Peak comparison at Gene:Peak ",
    "ratios 2:1, 3:1, 4:1, 6:1, 8:1, and 10:1."
  )
  if (!marker %in% readme) {
    readme <- c(
      readme,
      "",
      marker,
      "",
      benchmark_note
    )
  } else {
    marker_index <- match(marker, readme)
    note_index <- which(
      seq_along(readme) > marker_index &
        startsWith(readme, "See `16_Method10_P1e10_PeakSupport_GeneWeight_Benchmark/")
    )
    if (length(note_index)) {
      readme[note_index[[1L]]] <- benchmark_note
    } else {
      readme <- append(readme, c("", benchmark_note), after = marker_index)
    }
  }
  writeLines(readme, local_readme, useBytes = TRUE)

  for (path in c(local_index, local_readme)) {
    remote <- paste0(box_parent, "/", basename(path))
    status <- system2(
      "rclone",
      c("copyto", shQuote(path), shQuote(remote)),
      stdout = TRUE,
      stderr = TRUE
    )
    if (!identical(attr(status, "status") %||% 0L, 0L)) {
      stop("Could not update Box file: ", basename(path))
    }
  }
  invisible(TRUE)
}

update_box_plans <- function() {
  plans <- data.table::data.table(
    source = file.path(
      benchmark_root,
      c("condition_specific_config.yaml", "tf_union_config.yaml")
    ),
    remote_name = unname(box_plan_names)
  )
  if (any(!file.exists(plans$source))) {
    stop("A Box plans YAML is missing locally.")
  }
  for (i in seq_len(nrow(plans))) {
    remote <- paste0(box_plans_root, "/", plans$remote_name[[i]])
    status <- system2(
      "rclone",
      c("copyto", shQuote(plans$source[[i]]), shQuote(remote), "--checksum"),
      stdout = TRUE,
      stderr = TRUE
    )
    if (!identical(attr(status, "status") %||% 0L, 0L)) {
      stop("Could not update Box plan: ", plans$remote_name[[i]])
    }
  }
  invisible(TRUE)
}

deliver_results <- function() {
  write_comparison_report()
  write_delivery_readme()
  write_design_configs(status = "delivered")
  files <- delivery_file_table()
  dir.create(delivery_root, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(
    files$source,
    file.path(delivery_root, files$file),
    overwrite = TRUE,
    copy.mode = TRUE
  )
  if (!all(copied)) stop("Could not assemble the compact local delivery.")
  manifest <- data.table::copy(files)
  manifest[, source := NULL]
  manifest[, `:=`(
    source_signature = source_signature,
    box_folder = box_root
  )]
  manifest_path <- file.path(delivery_root, "delivery_manifest.csv")
  data.table::fwrite(manifest, manifest_path)

  log_info("Uploading the compact benchmark delivery to Box.")
  status <- system2(
    "rclone",
    c("copy", shQuote(delivery_root), shQuote(box_root), "--checksum"),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!identical(attr(status, "status") %||% 0L, 0L)) {
    stop("Box upload failed: ", paste(status, collapse = " "))
  }
  check <- system2(
    "rclone",
    c("check", shQuote(delivery_root), shQuote(box_root), "--checksum"),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!identical(attr(check, "status") %||% 0L, 0L)) {
    stop("Box checksum validation failed: ", paste(check, collapse = " "))
  }
  expected_remote <- sort(c(files$file, basename(manifest_path)))
  remote_files <- system2(
    "rclone",
    c("lsf", shQuote(box_root), "--files-only"),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!identical(attr(remote_files, "status") %||% 0L, 0L) ||
      !identical(sort(remote_files), expected_remote)) {
    stop("The Box folder contains an unexpected or missing delivery file.")
  }
  update_box_parent_index()
  update_box_plans()
  for (i in seq_len(nrow(run_table))) {
    write_run_config(run_table[i], status = "delivered")
  }
  log_info("Validated the compact Box delivery: ", box_root)
  invisible(manifest)
}

if (stage == "documents") {
  prepare_documents()
} else if (stage == "train") {
  train_models()
} else if (stage == "qc") {
  write_qc_reports()
} else if (stage == "compare") {
  write_comparison_report()
} else if (stage == "deliver") {
  deliver_results()
} else {
  deliver_results()
}
