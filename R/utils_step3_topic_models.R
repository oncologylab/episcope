# utils_grn_topic_tf_docs_warplda.R
# TF-centered topic modeling for differential GRN edges (WarpLDA-only; cisTopic-like model selection)
# Author: Yaoxiang Li
# Updated: 2026-01-13

# =============================================================================
# 0) Small helpers + assertions
# =============================================================================

if (!exists(".log_inform", mode = "function")) {
  .log_inform <- function(msg, ..., .envir = parent.frame()) {
    cli::cli_inform(msg, ..., .envir = .envir)
  }
}
if (!exists(".log_warn", mode = "function")) {
  .log_warn <- function(msg, ..., .envir = parent.frame()) {
    cli::cli_warn(msg, ..., .envir = .envir)
  }
}
if (!exists(".log_abort", mode = "function")) {
  .log_abort <- function(msg, ..., .envir = parent.frame()) {
    cli::cli_abort(msg, ..., .envir = .envir)
  }
}

.assert_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .log_abort("Package {.pkg {pkg}} is required but not installed.")
  }
  invisible(TRUE)
}

.topic_extraction_step_names <- function() {
  c(
    "topic_terms",
    "topic_links",
    "gammafit_summary",
    "link_efdr_summary",
    "tf_topic_assignment",
    "raw_theta_documents",
    "document_theta_umap",
    "topic_term_heatmap",
    "topic_by_comparison",
    "pathway",
    "intertopic_distance"
  )
}

.normalize_topic_extraction_steps <- function(steps = NULL) {
  if (is.null(steps)) return(NULL)
  if (is.list(steps)) steps <- unlist(steps, use.names = FALSE)
  steps <- unlist(strsplit(as.character(steps), ",", fixed = TRUE), use.names = FALSE)
  steps <- unique(trimws(steps))
  steps <- steps[nzchar(steps)]
  if (!length(steps)) return(NULL)
  aliases <- list(
    all = .topic_extraction_step_names(),
    default = c(
      "topic_terms",
      "topic_links",
      "gammafit_summary",
      "tf_topic_assignment",
      "document_theta_umap",
      "topic_term_heatmap",
      "topic_by_comparison",
      "pathway"
    ),
    core = c("topic_terms", "topic_links", "gammafit_summary"),
    plots = c(
      "gammafit_summary",
      "tf_topic_assignment",
      "document_theta_umap",
      "topic_term_heatmap",
      "topic_by_comparison",
      "intertopic_distance"
    ),
    heatmaps = c("tf_topic_assignment", "topic_term_heatmap", "topic_by_comparison"),
    reports = c("gammafit_summary", "tf_topic_assignment", "document_theta_umap", "topic_term_heatmap", "topic_by_comparison", "pathway")
  )
  expanded <- unlist(lapply(steps, function(step) {
    if (step %in% names(aliases)) aliases[[step]] else step
  }), use.names = FALSE)
  expanded <- unique(expanded)
  bad <- setdiff(expanded, .topic_extraction_step_names())
  if (length(bad)) {
    .log_abort("Unsupported extraction step(s): {paste(bad, collapse = ', ')}.")
  }
  expanded
}

.topic_step_enabled <- function(steps, step, default = TRUE) {
  if (is.null(steps)) return(isTRUE(default))
  step %in% steps
}

.assert_has_cols <- function(df, cols, context = NULL) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    msg <- if (!is.null(context)) sprintf("{.strong %s}: missing columns.", context) else "Missing columns."
    .log_abort(c(msg, i = paste(miss, collapse = ", ")))
  }
  invisible(TRUE)
}

.safe_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- NA_real_
  x
}

.safe_log2fc <- function(num, den) {
  num <- .safe_num(num)
  den <- .safe_num(den)
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & num > 0 & den > 0
  out[ok] <- log2(num[ok] / den[ok])
  out
}

.safe_sign <- function(x) {
  x <- .safe_num(x)
  x[!is.finite(x)] <- 0
  ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
}

# fold-change magnitude from log2FC (symmetric for up/down)
.fc_mag_from_log2fc <- function(log2fc) {
  log2fc <- .safe_num(log2fc)
  out <- 2 ^ abs(log2fc)
  out[!is.finite(out)] <- NA_real_
  out
}

# =============================================================================
# 1) Load + standardize your Step2 *_delta_links.csv files
# =============================================================================

parse_delta_links_filename <- function(file) {
  b <- basename(file)
  direction <- NULL
  if (grepl("_filtered_links_up\\.csv$", b, ignore.case = TRUE)) {
    direction <- "up"
  } else if (grepl("_filtered_links_down\\.csv$", b, ignore.case = TRUE)) {
    direction <- "down"
  }
  if (grepl("_delta_links_filtered_up\\.csv$", b, ignore.case = TRUE)) {
    direction <- "up"
  } else if (grepl("_delta_links_filtered_down\\.csv$", b, ignore.case = TRUE)) {
    direction <- "down"
  }
  b <- sub("_filtered_links_(up|down)\\.csv$", "", b, ignore.case = TRUE)
  b <- sub("_delta_links_filtered_(up|down)\\.csv$", "", b, ignore.case = TRUE)
  b <- sub("_filtered_links\\.csv$", "", b, ignore.case = TRUE)
  b <- sub("_delta_links_filtered\\.csv$", "", b, ignore.case = TRUE)
  b <- sub("_delta_links\\.csv$", "", b, ignore.case = TRUE)
  parts <- strsplit(b, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(list(comparison_id = b, cond1_id = NA_character_, cond2_id = NA_character_, direction = direction))
  }
  list(comparison_id = b, cond1_id = parts[[1]], cond2_id = parts[[2]], direction = direction)
}

standardize_delta_links_one <- function(file, keep_original = TRUE) {
  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("readr")

  ids <- parse_delta_links_filename(file)
  dt <- readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
  dt <- data.table::as.data.table(dt)

  .assert_has_cols(dt, c("tf", "gene_key", "peak_id"), context = "standardize_delta_links_one")
  if ("comparison_id" %in% names(dt) && any(nzchar(as.character(dt$comparison_id)))) {
    ids$comparison_id <- as.character(dt$comparison_id[which(nzchar(as.character(dt$comparison_id)))[[1L]]])
  }
  if ("cond1_matrix_id" %in% names(dt) && any(nzchar(as.character(dt$cond1_matrix_id)))) {
    ids$cond1_id <- as.character(dt$cond1_matrix_id[which(nzchar(as.character(dt$cond1_matrix_id)))[[1L]]])
  } else if ("cond1_id" %in% names(dt) && any(nzchar(as.character(dt$cond1_id)))) {
    ids$cond1_id <- as.character(dt$cond1_id[which(nzchar(as.character(dt$cond1_id)))[[1L]]])
  } else if ("case_id" %in% names(dt) && any(nzchar(as.character(dt$case_id)))) {
    ids$cond1_id <- as.character(dt$case_id[which(nzchar(as.character(dt$case_id)))[[1L]]])
  }
  if ("cond2_matrix_id" %in% names(dt) && any(nzchar(as.character(dt$cond2_matrix_id)))) {
    ids$cond2_id <- as.character(dt$cond2_matrix_id[which(nzchar(as.character(dt$cond2_matrix_id)))[[1L]]])
  } else if ("cond2_id" %in% names(dt) && any(nzchar(as.character(dt$cond2_id)))) {
    ids$cond2_id <- as.character(dt$cond2_id[which(nzchar(as.character(dt$cond2_id)))[[1L]]])
  } else if ("ctrl_id" %in% names(dt) && any(nzchar(as.character(dt$ctrl_id)))) {
    ids$cond2_id <- as.character(dt$ctrl_id[which(nzchar(as.character(dt$ctrl_id)))[[1L]]])
  }
  if (is.na(ids$cond1_id) || is.na(ids$cond2_id) || !nzchar(ids$cond1_id) || !nzchar(ids$cond2_id)) {
    if (!nrow(dt)) {
      if (!"comparison_id" %in% names(dt)) dt[, comparison_id := ids$comparison_id]
      if (!"comparison_label" %in% names(dt)) dt[, comparison_label := ids$comparison_id]
      if (!"cond1_id" %in% names(dt)) dt[, cond1_id := ids$cond1_id]
      if (!"cond2_id" %in% names(dt)) dt[, cond2_id := ids$cond2_id]
      if (!is.null(ids$direction) && !"direction_group" %in% names(dt)) dt[, direction_group := ids$direction]
      if (!isTRUE(keep_original)) {
        keep <- c(
          "comparison_id", "comparison_label", "cond1_id", "cond2_id", "cond1_label", "cond2_label", "direction_group",
          "tf", "gene_key", "peak_id",
          "fp_bound_cond1", "fp_bound_cond2",
          "tf_expr_flag_cond1", "tf_expr_flag_cond2",
          "gene_expr_flag_cond1", "gene_expr_flag_cond2",
          "tf_expr_cond1", "tf_expr_cond2", "gene_expr_cond1", "gene_expr_cond2",
          "fp_score_cond1", "fp_score_cond2",
          "delta_fp", "delta_gene", "log2fc_fp", "log2fc_gene", "fc_mag_fp", "fc_mag_gene",
          "log2fc_tf", "fc_mag_tf"
        )
        missing <- setdiff(keep, names(dt))
        for (nm_missing in missing) dt[[nm_missing]] <- if (nm_missing %in% c("comparison_id", "comparison_label", "cond1_id", "cond2_id", "cond1_label", "cond2_label", "direction_group", "tf", "gene_key", "peak_id")) character() else numeric()
        return(dt[, keep, with = FALSE])
      }
      return(dt)
    }
    .log_abort(c(
      "Cannot determine condition labels for Module 3 differential links.",
      i = "Use <COND1>_vs_<COND2> filenames or include cond1_id and cond2_id columns.",
      i = paste0("Got: ", basename(file))
    ))
  }

  nm <- function(prefix, cond) paste0(prefix, "_", cond)
  col_first <- function(choices) {
    hit <- intersect(choices, names(dt))
    if (length(hit)) hit[[1L]] else NULL
  }

  # expected wide names
  fp_bound_cond1_nm <- col_first(c("fp_bound_cond1", "fp_bound_case", nm("fp_bound", ids$cond1_id)))
  fp_bound_cond2_nm <- col_first(c("fp_bound_cond2", "fp_bound_ctrl", nm("fp_bound", ids$cond2_id)))
  tf_expr_flag_cond1_nm <- col_first(c("tf_expr_flag_cond1", "tf_expr_flag_case", nm("tf_expr_flag", ids$cond1_id)))
  tf_expr_flag_cond2_nm <- col_first(c("tf_expr_flag_cond2", "tf_expr_flag_ctrl", nm("tf_expr_flag", ids$cond2_id)))
  gene_expr_flag_cond1_nm <- col_first(c("gene_expr_flag_cond1", "gene_expr_flag_case", nm("gene_expr_flag", ids$cond1_id)))
  gene_expr_flag_cond2_nm <- col_first(c("gene_expr_flag_cond2", "gene_expr_flag_ctrl", nm("gene_expr_flag", ids$cond2_id)))

  tf_expr_cond1_nm <- col_first(c("tf_expr_cond1", "tf_expr_case", nm("tf_expr", ids$cond1_id)))
  tf_expr_cond2_nm <- col_first(c("tf_expr_cond2", "tf_expr_ctrl", nm("tf_expr", ids$cond2_id)))
  gene_expr_cond1_nm <- col_first(c("gene_expr_cond1", "gene_expr_case", nm("gene_expr", ids$cond1_id)))
  gene_expr_cond2_nm <- col_first(c("gene_expr_cond2", "gene_expr_ctrl", nm("gene_expr", ids$cond2_id)))
  fp_score_cond1_nm <- col_first(c("fp_score_cond1", "fp_score_case", nm("fp_score", ids$cond1_id), nm("fp_bed_score", ids$cond1_id)))
  fp_score_cond2_nm <- col_first(c("fp_score_cond2", "fp_score_ctrl", nm("fp_score", ids$cond2_id), nm("fp_bed_score", ids$cond2_id)))

  has <- function(x) length(x) == 1L && !is.na(x) && x %in% names(dt)

  dt[, comparison_id := ids$comparison_id]
  dt[, cond1_id := if ("cond1_id" %in% names(dt)) as.character(cond1_id) else ids$cond1_id]
  dt[, cond2_id := if ("cond2_id" %in% names(dt)) as.character(cond2_id) else ids$cond2_id]
  if (!"cond1_label" %in% names(dt)) dt[, cond1_label := cond1_id]
  if (!"cond2_label" %in% names(dt)) dt[, cond2_label := cond2_id]
  if (!"comparison_label" %in% names(dt)) {
    dt[, comparison_label := comparison_id]
  }
  dt[, comparison_label := as.character(comparison_label)]
  dt[is.na(comparison_label) | !nzchar(trimws(comparison_label)), comparison_label := comparison_id]
  if (!is.null(ids$direction)) {
    dt[, direction_group := ids$direction]
  }

  # flags/bound (missing -> 0)
  dt[, fp_bound_cond1 := if (has(fp_bound_cond1_nm)) as.integer(get(fp_bound_cond1_nm)) else 0L]
  dt[, fp_bound_cond2 := if (has(fp_bound_cond2_nm)) as.integer(get(fp_bound_cond2_nm)) else 0L]
  dt[, tf_expr_flag_cond1 := if (has(tf_expr_flag_cond1_nm)) as.integer(get(tf_expr_flag_cond1_nm)) else 0L]
  dt[, tf_expr_flag_cond2 := if (has(tf_expr_flag_cond2_nm)) as.integer(get(tf_expr_flag_cond2_nm)) else 0L]
  dt[, gene_expr_flag_cond1 := if (has(gene_expr_flag_cond1_nm)) as.integer(get(gene_expr_flag_cond1_nm)) else 0L]
  dt[, gene_expr_flag_cond2 := if (has(gene_expr_flag_cond2_nm)) as.integer(get(gene_expr_flag_cond2_nm)) else 0L]

  # raw expr + raw fp (optional)
  dt[, tf_expr_cond1 := if (has(tf_expr_cond1_nm)) .safe_num(get(tf_expr_cond1_nm)) else NA_real_]
  dt[, tf_expr_cond2 := if (has(tf_expr_cond2_nm)) .safe_num(get(tf_expr_cond2_nm)) else NA_real_]
  dt[, gene_expr_cond1 := if (has(gene_expr_cond1_nm)) .safe_num(get(gene_expr_cond1_nm)) else NA_real_]
  dt[, gene_expr_cond2 := if (has(gene_expr_cond2_nm)) .safe_num(get(gene_expr_cond2_nm)) else NA_real_]
  dt[, fp_score_cond1 := if (has(fp_score_cond1_nm)) .safe_num(get(fp_score_cond1_nm)) else NA_real_]
  dt[, fp_score_cond2 := if (has(fp_score_cond2_nm)) .safe_num(get(fp_score_cond2_nm)) else NA_real_]

  # standardized deltas/log2FC and FC-magnitude (gene/fp)
  dt[, delta_fp := if (has("delta_fp_score")) .safe_num(delta_fp_score) else if (has("delta_fp_bed_score")) .safe_num(delta_fp_bed_score) else NA_real_]
  dt[, delta_gene := if (has("delta_gene_expr")) .safe_num(delta_gene_expr) else NA_real_]

  dt[, log2fc_fp := {
    if (has("log2FC_fp_score")) .safe_num(log2FC_fp_score)
    else if (has("log2FC_fp_bed_score")) .safe_num(log2FC_fp_bed_score)
    else if (has("fc_fp_score")) .safe_num(log2(fc_fp_score))
    else if (has("fc_fp_bed_score")) .safe_num(log2(fc_fp_bed_score))
    else if (has("delta_fp_score")) .safe_num(delta_fp_score)
    else if (has("delta_fp_bed_score")) .safe_num(delta_fp_bed_score)
    else NA_real_
  }]

  dt[, log2fc_gene := {
    if (has("log2FC_gene_expr")) .safe_num(log2FC_gene_expr)
    else if (has("fc_gene_expr")) .safe_num(log2(fc_gene_expr))
    else if (has("delta_gene_expr")) .safe_num(delta_gene_expr)
    else NA_real_
  }]

  dt[, fc_mag_fp := {
    if (has("fc_fp_score") || has("fc_fp_bed_score")) {
      # symmetric magnitude
      v <- if (has("fc_fp_score")) .safe_num(fc_fp_score) else .safe_num(fc_fp_bed_score)
      v[v <= 0] <- NA_real_
      w <- pmax(v, 1 / v)
      w[!is.finite(w)] <- NA_real_
      w
    } else {
      .fc_mag_from_log2fc(log2fc_fp)
    }
  }]

  dt[, fc_mag_gene := {
    if (has("fc_gene_expr")) {
      v <- .safe_num(fc_gene_expr)
      v[v <= 0] <- NA_real_
      w <- pmax(v, 1 / v)
      w[!is.finite(w)] <- NA_real_
      w
    } else {
      .fc_mag_from_log2fc(log2fc_gene)
    }
  }]

  dt[, log2fc_tf := if (has("log2FC_tf_expr")) .safe_num(log2FC_tf_expr) else {
    if (has("tf_expr_cond1") && has("tf_expr_cond2")) {
      .safe_log2fc(tf_expr_cond1, tf_expr_cond2)
    } else {
      NA_real_
    }
  }]
  dt[, fc_mag_tf := .fc_mag_from_log2fc(log2fc_tf)]

  if (!isTRUE(keep_original)) {
    keep <- c(
      "comparison_id","comparison_label","cond1_id","cond2_id","cond1_label","cond2_label","direction_group",
      "tf","gene_key","peak_id",
      "fp_bound_cond1","fp_bound_cond2",
      "tf_expr_flag_cond1","tf_expr_flag_cond2",
      "gene_expr_flag_cond1","gene_expr_flag_cond2",
      "tf_expr_cond1","tf_expr_cond2","gene_expr_cond1","gene_expr_cond2",
      "fp_score_cond1","fp_score_cond2",
      "delta_fp","delta_gene","log2fc_fp","log2fc_gene","fc_mag_fp","fc_mag_gene",
      "log2fc_tf","fc_mag_tf"
    )
    keep <- intersect(keep, names(dt))
    dt <- dt[, ..keep]
  }

  dt[]
}

#' Load filtered differential links for Module 3 topics
#'
#' Read one or more filtered differential-link CSV files and standardize the
#' columns used by Module 3 topic-model input construction.
#'
#' @param files Vector of delta-link CSV paths.
#' @param keep_original Keep original columns when loading delta links.
#' @param n_max_files Optional maximum number of files to load.
#'
#' @return A data.table of standardized differential links.
#' @noRd
load_delta_links_many <- function(files, keep_original = TRUE, n_max_files = Inf) {
  .assert_pkg("cli")
  .assert_pkg("data.table")

  files <- as.character(files)
  files <- files[file.exists(files)]
  if (!length(files)) .log_abort("No existing files provided to load_delta_links_many().")
  if (is.finite(n_max_files)) files <- head(files, as.integer(n_max_files))

  lst <- lapply(files, function(f) standardize_delta_links_one(f, keep_original = keep_original))
  data.table::rbindlist(lst, use.names = TRUE, fill = TRUE)
}

.apply_module3_manifest_labels <- function(edges_dt, filtered_dir) {
  .assert_pkg("data.table")
  dt <- data.table::as.data.table(edges_dt)
  manifest <- file.path(as.character(filtered_dir)[[1L]], "filtered_links_manifest.csv")
  if (!file.exists(manifest) || !"comparison_id" %in% names(dt)) {
    return(dt[])
  }
  man <- data.table::fread(manifest, showProgress = FALSE)
  if (!all(c("comparison_id", "comparison_label") %in% names(man))) {
    return(dt[])
  }
  label_map <- unique(man[
    !is.na(comparison_id) & nzchar(as.character(comparison_id)) &
      !is.na(comparison_label) & nzchar(trimws(as.character(comparison_label))),
    .(
      comparison_id = as.character(comparison_id),
      comparison_label_manifest = as.character(comparison_label)
    )
  ])
  if (!nrow(label_map)) {
    return(dt[])
  }
  dt[, comparison_id := as.character(comparison_id)]
  dt <- merge(dt, label_map, by = "comparison_id", all.x = TRUE, sort = FALSE)
  if (!"comparison_label" %in% names(dt)) {
    dt[, comparison_label := comparison_id]
  }
  use_manifest <- !is.na(dt$comparison_label_manifest) & nzchar(trimws(dt$comparison_label_manifest))
  dt[use_manifest, comparison_label := comparison_label_manifest]
  dt[, comparison_label_manifest := NULL]
  dt[]
}

.module3_filtered_link_files <- function(path) {
  .assert_pkg("data.table")
  path <- as.character(path)
  if (!length(path)) .log_abort("`path` must contain at least one file or directory.")
  if (length(path) == 1L && dir.exists(path)) {
    manifest <- file.path(path, "filtered_links_manifest.csv")
    if (file.exists(manifest)) {
      man <- data.table::fread(manifest, showProgress = FALSE)
      path_cols <- intersect(c("up_path", "down_path"), names(man))
      files <- if (length(path_cols)) {
        unique(as.character(unlist(man[, path_cols, with = FALSE], use.names = FALSE)))
      } else {
        character()
      }
      files <- files[!is.na(files) & nzchar(files)]
      files <- ifelse(file.exists(files), files, file.path(path, basename(files)))
    } else {
      files <- list.files(path, "_filtered_links_(up|down)[.]csv$", full.names = TRUE)
    }
  } else if (length(path) == 1L && file.exists(path) && grepl("manifest[.]csv$", basename(path))) {
    man <- data.table::fread(path, showProgress = FALSE)
    path_cols <- intersect(c("up_path", "down_path"), names(man))
    files <- if (length(path_cols)) {
      unique(as.character(unlist(man[, path_cols, with = FALSE], use.names = FALSE)))
    } else {
      character()
    }
    files <- files[!is.na(files) & nzchar(files)]
    manifest_dir <- dirname(path)
    files <- ifelse(file.exists(files), files, file.path(manifest_dir, basename(files)))
  } else {
    files <- path
  }
  files[file.exists(files)]
}

.module3_topic_input_signature <- function(input_dir,
                                           input_source,
                                           sample_subset = NULL,
                                           settings = list()) {
  .assert_pkg("data.table")
  .assert_pkg("digest")
  input_source <- match.arg(input_source, c("differential_links", "condition_links"))
  sample_subset <- sort(unique(as.character(sample_subset)))
  sample_subset <- sample_subset[!is.na(sample_subset) & nzchar(sample_subset)]
  files <- if (identical(input_source, "condition_links")) {
    manifest_path <- .module3_condition_link_manifest_path(input_dir)
    manifest <- data.table::fread(manifest_path, showProgress = FALSE)
    if (!all(c("condition_id", "path") %in% names(manifest))) {
      .log_abort("Condition-link manifest must include condition_id and path for cache validation.")
    }
    if (length(sample_subset)) manifest <- manifest[condition_id %in% sample_subset]
    paths <- as.character(manifest$path)
    paths[!grepl("^/", paths)] <- file.path(dirname(manifest_path), paths[!grepl("^/", paths)])
    c(manifest_path, paths)
  } else {
    manifest_path <- file.path(input_dir, "filtered_links_manifest.csv")
    c(if (file.exists(manifest_path)) manifest_path else character(), .module3_filtered_link_files(input_dir))
  }
  files <- sort(unique(normalizePath(files[file.exists(files)], winslash = "/", mustWork = TRUE)))
  info <- file.info(files)
  file_rows <- data.table::data.table(
    path = files,
    size = as.numeric(info$size),
    mtime = as.numeric(info$mtime)
  )
  digest::digest(
    list(
      signature_version = 1L,
      input_source = input_source,
      sample_subset = sample_subset,
      files = file_rows,
      settings = settings
    ),
    algo = "xxhash64",
    serialize = TRUE
  )
}

.module3_apply_condition_gene_specificity <- function(doc_term,
                                                      expression_file,
                                                      count_column,
                                                      expression_min,
                                                      temperature = 0.5,
                                                      uniform_floor = 0.1) {
  .assert_pkg("data.table")
  dt <- data.table::as.data.table(doc_term)
  required <- c("doc_id", "term_id", count_column)
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    .log_abort(
      "Condition gene-specificity weighting requires doc-term columns: {paste(missing, collapse = ', ')}."
    )
  }
  if (!file.exists(expression_file)) {
    .log_abort("Condition gene expression file not found: {expression_file}")
  }
  temperature <- suppressWarnings(as.numeric(temperature[[1L]]))
  uniform_floor <- suppressWarnings(as.numeric(uniform_floor[[1L]]))
  expression_min <- suppressWarnings(as.numeric(expression_min[[1L]]))
  if (!is.finite(temperature) || temperature <= 0) {
    .log_abort("Condition specificity temperature must be positive.")
  }
  if (!is.finite(uniform_floor) || uniform_floor <= 0 || uniform_floor >= 1) {
    .log_abort("Condition specificity floor must be strictly between 0 and 1.")
  }
  if (!is.finite(expression_min) || expression_min < 0) {
    .log_abort("Condition specificity expression minimum must be non-negative.")
  }

  expression <- data.table::fread(expression_file, showProgress = FALSE)
  expression_required <- c("condition_id", "gene_key", "expression")
  expression_missing <- setdiff(expression_required, names(expression))
  if (length(expression_missing)) {
    .log_abort(
      "Condition gene expression file is missing columns: {paste(expression_missing, collapse = ', ')}."
    )
  }
  expression <- expression[, .(
    condition_id = as.character(condition_id),
    gene_key = as.character(gene_key),
    expression = suppressWarnings(as.numeric(expression))
  )]
  expression <- expression[
    !is.na(condition_id) & nzchar(condition_id) &
      !is.na(gene_key) & nzchar(gene_key) & is.finite(expression)
  ]
  if (!nrow(expression)) {
    .log_abort("Condition gene expression file has no usable rows: {expression_file}")
  }
  duplicate_rows <- expression[, .N, by = .(condition_id, gene_key)][N > 1L]
  if (nrow(duplicate_rows)) {
    .log_abort(
      "Condition gene expression must contain one row per condition and gene; found {nrow(duplicate_rows)} duplicate pair(s)."
    )
  }

  expression[, expressed := expression >= expression_min]
  expression[, scaled := log2(expression + 1) / temperature]
  expression[expressed == FALSE, scaled := -Inf]
  expression[, centered := scaled - if (any(is.finite(scaled))) {
    max(scaled[is.finite(scaled)])
  } else {
    0
  }, by = gene_key]
  expression[!is.finite(centered), centered := -Inf]
  expression[, specificity_mass := exp(centered)]
  expression[!is.finite(specificity_mass), specificity_mass := 0]
  expression[, `:=`(
    specificity_total = sum(specificity_mass),
    expressed_conditions = sum(expressed)
  ), by = gene_key]
  expression[, multiplier := uniform_floor]
  expression[expressed & specificity_total > 0, multiplier :=
    uniform_floor +
      (1 - uniform_floor) *
      specificity_mass / specificity_total * expressed_conditions]

  doc_id <- as.character(dt$doc_id)
  term_id <- as.character(dt$term_id)
  doc_condition <- sub("::[^:]+$", "", doc_id)
  target_term <- startsWith(term_id, "GENE:") | startsWith(term_id, "PEAK:")
  if (!any(target_term)) {
    .log_abort("Condition gene-specificity weighting found no GENE: or PEAK: terms.")
  }
  term_gene <- sub("^(GENE|PEAK):", "", term_id)
  lookup <- expression[, .(condition_id, gene_key, multiplier)]
  data.table::setkey(lookup, condition_id, gene_key)
  matched_multiplier <- lookup[
    data.table::data.table(condition_id = doc_condition, gene_key = term_gene),
    on = .(condition_id, gene_key),
    multiplier
  ]
  multiplier <- rep(1, nrow(dt))
  multiplier[target_term] <- uniform_floor
  matched <- target_term & is.finite(matched_multiplier)
  multiplier[matched] <- matched_multiplier[matched]

  original <- suppressWarnings(as.numeric(dt[[count_column]]))
  if (any(!is.finite(original) | original <= 0)) {
    .log_abort("Condition gene-specificity weighting requires positive finite topic counts.")
  }
  document_levels <- unique(doc_id)
  document_index <- match(doc_id, document_levels)
  original_mass <- as.numeric(rowsum(original, document_index, reorder = FALSE))
  weighted <- original * multiplier
  weighted_mass <- as.numeric(rowsum(weighted, document_index, reorder = FALSE))
  weighted <- weighted * (original_mass / weighted_mass)[document_index]
  weighted <- pmax(1, round(weighted))

  rounded_mass <- as.numeric(rowsum(weighted, document_index, reorder = FALSE))
  mass_delta <- round(original_mass) - rounded_mass
  if (any(mass_delta != 0)) {
    row_index <- data.table::data.table(
      document_index = document_index,
      row_index = seq_along(document_index),
      count = weighted
    )[, row_index[which.max(count)], by = document_index]$V1
    candidate <- weighted[row_index] + mass_delta
    direct <- candidate >= 1
    weighted[row_index[direct]] <- candidate[direct]
    difficult <- which(!direct)
    for (document in difficult) {
      remaining <- -mass_delta[[document]]
      rows <- which(document_index == document & weighted > 1)
      rows <- rows[order(weighted[rows], decreasing = TRUE)]
      for (row in rows) {
        remove <- min(weighted[[row]] - 1, remaining)
        weighted[[row]] <- weighted[[row]] - remove
        remaining <- remaining - remove
        if (remaining <= 0) break
      }
      if (remaining > 0) {
        .log_abort("Unable to preserve document token mass after specificity weighting.")
      }
    }
  }
  final_mass <- as.numeric(rowsum(weighted, document_index, reorder = FALSE))
  if (any(final_mass != round(original_mass))) {
    .log_abort("Condition gene-specificity weighting did not preserve document token totals.")
  }
  dt[, (count_column) := weighted]

  audit <- data.table::data.table(
    weighting = "gene_specificity",
    expression_file = normalizePath(expression_file, winslash = "/", mustWork = TRUE),
    expression_min = expression_min,
    temperature = temperature,
    uniform_floor = uniform_floor,
    n_documents = length(document_levels),
    n_doc_term_rows = nrow(dt),
    n_target_term_rows = sum(target_term),
    n_expression_matched_rows = sum(matched),
    expression_matched_fraction = sum(matched) / max(1, sum(target_term)),
    multiplier_min = min(multiplier[target_term]),
    multiplier_p01 = unname(stats::quantile(multiplier[target_term], 0.01)),
    multiplier_median = stats::median(multiplier[target_term]),
    multiplier_p99 = unname(stats::quantile(multiplier[target_term], 0.99)),
    multiplier_max = max(multiplier[target_term]),
    original_tokens = sum(round(original_mass)),
    weighted_tokens = sum(final_mass),
    max_document_token_difference = max(abs(final_mass - round(original_mass)))
  )
  list(doc_term = dt, audit = audit)
}

#' Load Module 3 differential links
#'
#' Load filtered differential-link files produced by
#' [module3_prepare_differential_links()]. The input can be a Module 3
#' differential-link directory, a manifest path, or one or more filtered-link
#' CSV files.
#'
#' @param path Module 3 differential-link directory, manifest CSV, or filtered
#'   link CSV path(s).
#' @param keep_original Keep original columns in addition to standardized
#'   Module 3 columns.
#' @param n_max_files Optional maximum number of files to load.
#'
#' @return A data.table of standardized differential links.
#' @noRd
load_differential_links <- function(path, keep_original = FALSE, n_max_files = Inf) {
  files <- .module3_filtered_link_files(path)
  if (!length(files)) .log_abort("No Module 3 filtered-link CSV files were found.")
  load_delta_links_many(files, keep_original = keep_original, n_max_files = n_max_files)
}

#' Query Module 3 differential links
#'
#' Filter Module 3 differential links by comparison, direction, TF, target gene,
#' footprint, or distance to TSS.
#'
#' @param links Differential-link data, a Module 3 differential-link directory,
#'   a manifest path, or filtered-link CSV path(s).
#' @param comparison_id Optional comparison IDs to keep.
#' @param direction Optional direction labels, usually `"up"` and/or `"down"`.
#' @param tf Optional TF names to keep.
#' @param gene Optional target-gene names to keep.
#' @param fp_id Optional footprint IDs to keep.
#' @param max_distance_to_tss Optional maximum absolute distance to TSS.
#' @param keep_original Keep original columns when `links` is a path.
#'
#' @return A data.table of matching Module 3 differential links.
#' @noRd
query_differential_links <- function(links,
                                     comparison_id = NULL,
                                     direction = NULL,
                                     tf = NULL,
                                     gene = NULL,
                                     fp_id = NULL,
                                     max_distance_to_tss = NULL,
                                     keep_original = FALSE) {
  .assert_pkg("data.table")
  dt <- if (is.data.frame(links)) {
    data.table::as.data.table(links)
  } else {
    load_differential_links(links, keep_original = keep_original)
  }
  if (!is.null(comparison_id) && "comparison_id" %in% names(dt)) {
    keep <- as.character(comparison_id)
    dt <- dt[dt[["comparison_id"]] %in% keep]
  }
  if (!is.null(direction) && "direction_group" %in% names(dt)) {
    keep <- as.character(direction)
    dt <- dt[dt[["direction_group"]] %in% keep]
  }
  if (!is.null(tf) && "tf" %in% names(dt)) {
    keep <- as.character(tf)
    dt <- dt[dt[["tf"]] %in% keep]
  }
  if (!is.null(gene) && "gene_key" %in% names(dt)) {
    keep <- as.character(gene)
    dt <- dt[dt[["gene_key"]] %in% keep]
  }
  if (!is.null(fp_id) && "peak_id" %in% names(dt)) {
    keep <- as.character(fp_id)
    dt <- dt[dt[["peak_id"]] %in% keep]
  }
  if (!is.null(max_distance_to_tss) && "distance_to_tss" %in% names(dt)) {
    max_dist <- as.numeric(max_distance_to_tss)[[1L]]
    dt <- dt[is.finite(dt[["distance_to_tss"]]) & abs(dt[["distance_to_tss"]]) <= max_dist]
  }
  dt[]
}

# =============================================================================
# 2) Step 1: Per-comparison GRN link filtering (edge-level QC)
# =============================================================================

filter_edges_for_tf_topics <- function(edges,
                                       abs_log2fc_fp_min = 1,
                                       abs_delta_fp_min = NA_real_,
                                       abs_log2fc_gene_min = 1,
                                       require_fp_bound_either = TRUE,
                                       require_tf_expr_either = TRUE,
                                       require_gene_expr_either = TRUE,
                                       direction_consistency = c("aligned", "none")) {
  direction_consistency <- match.arg(direction_consistency)
  .assert_pkg("data.table")

  dt <- data.table::as.data.table(edges)
  if (!"fp_bound_cond1" %in% names(dt) && "fp_bound_case" %in% names(dt)) dt[, fp_bound_cond1 := fp_bound_case]
  if (!"fp_bound_cond2" %in% names(dt) && "fp_bound_ctrl" %in% names(dt)) dt[, fp_bound_cond2 := fp_bound_ctrl]
  if (!"tf_expr_flag_cond1" %in% names(dt) && "tf_expr_flag_case" %in% names(dt)) dt[, tf_expr_flag_cond1 := tf_expr_flag_case]
  if (!"tf_expr_flag_cond2" %in% names(dt) && "tf_expr_flag_ctrl" %in% names(dt)) dt[, tf_expr_flag_cond2 := tf_expr_flag_ctrl]
  if (!"gene_expr_flag_cond1" %in% names(dt) && "gene_expr_flag_case" %in% names(dt)) dt[, gene_expr_flag_cond1 := gene_expr_flag_case]
  if (!"gene_expr_flag_cond2" %in% names(dt) && "gene_expr_flag_ctrl" %in% names(dt)) dt[, gene_expr_flag_cond2 := gene_expr_flag_ctrl]
  .assert_has_cols(
    dt,
    c("comparison_id","tf","gene_key","peak_id","log2fc_fp","log2fc_gene"),
    context = "filter_edges_for_tf_topics"
  )

  # flags/bound: missing -> 0
  get_flag <- function(nm) {
    if (nm %in% names(dt)) {
      x <- suppressWarnings(as.integer(dt[[nm]]))
      x[!is.finite(x)] <- 0L
      x
    } else {
      rep.int(0L, nrow(dt))
    }
  }

  tf_either <- (get_flag("tf_expr_flag_cond1") >= 1L) | (get_flag("tf_expr_flag_cond2") >= 1L)
  gene_either <- (get_flag("gene_expr_flag_cond1") >= 1L) | (get_flag("gene_expr_flag_cond2") >= 1L)
  fp_either <- (get_flag("fp_bound_cond1") >= 1L) | (get_flag("fp_bound_cond2") >= 1L)

  keep <- rep(TRUE, nrow(dt))
  if (isTRUE(require_tf_expr_either)) keep <- keep & tf_either
  if (isTRUE(require_gene_expr_either)) keep <- keep & gene_either
  if (isTRUE(require_fp_bound_either)) keep <- keep & fp_either

  l2fp <- .safe_num(dt[["log2fc_fp"]])
  l2gn <- .safe_num(dt[["log2fc_gene"]])

  if (is.finite(abs_log2fc_fp_min) && abs_log2fc_fp_min > 0) {
    keep <- keep & is.finite(l2fp) & (abs(l2fp) >= abs_log2fc_fp_min)
  }
  if (is.finite(abs_delta_fp_min) && abs_delta_fp_min > 0) {
    if ("delta_fp" %in% names(dt)) {
      dlt <- .safe_num(dt[["delta_fp"]])
    } else if ("delta_fp_score" %in% names(dt)) {
      dlt <- .safe_num(dt[["delta_fp_score"]])
    } else if ("delta_fp_bed_score" %in% names(dt)) {
      dlt <- .safe_num(dt[["delta_fp_bed_score"]])
    } else {
      .log_abort("abs_delta_fp_min requires delta_fp or delta_fp_score in edges.")
    }
    keep <- keep & is.finite(dlt) & (abs(dlt) >= abs_delta_fp_min)
  }
  keep <- keep & is.finite(l2gn) & (abs(l2gn) >= abs_log2fc_gene_min)

  if (direction_consistency == "aligned") {
    fp_dir <- .safe_sign(l2fp)
    gene_dir <- .safe_sign(l2gn)
    if ("delta_fp" %in% names(dt)) {
      alt <- .safe_sign(dt[["delta_fp"]])
      fp_dir[fp_dir == 0L] <- alt[fp_dir == 0L]
    }
    if ("delta_gene" %in% names(dt)) {
      alt <- .safe_sign(dt[["delta_gene"]])
      gene_dir[gene_dir == 0L] <- alt[gene_dir == 0L]
    }
    keep <- keep & (fp_dir == gene_dir) & (fp_dir != 0L)
  }

  dt[keep]
}

# =============================================================================
# 3) Sparse DTM + WarpLDA-only model fitting across K
# =============================================================================

build_sparse_dtm <- function(doc_term, count_col = "pseudo_count") {
  .assert_pkg("cli")
  .assert_pkg("Matrix")
  .assert_pkg("data.table")

  dt <- data.table::as.data.table(doc_term)
  .assert_has_cols(dt, c("doc_id","term_id", count_col), context = "build_sparse_dtm")
  dt <- dt[is.finite(get(count_col)) & get(count_col) > 0]
  if (!nrow(dt)) .log_abort("No non-zero entries in doc_term; check filters/weights.")

  docs <- unique(dt$doc_id)
  terms <- unique(dt$term_id)

  doc_index <- stats::setNames(seq_along(docs), docs)
  term_index <- stats::setNames(seq_along(terms), terms)

  i <- unname(doc_index[dt$doc_id])
  j <- unname(term_index[dt$term_id])
  x <- as.numeric(dt[[count_col]])

  dtm <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(docs), length(terms)),
    dimnames = list(docs, terms)
  )

  list(dtm = dtm, doc_index = doc_index, term_index = term_index)
}

.warplda_available_threads <- function() {
  cores <- if (requireNamespace("parallelly", quietly = TRUE)) {
    suppressWarnings(parallelly::availableCores())
  } else {
    suppressWarnings(parallel::detectCores(logical = TRUE))
  }
  cores <- suppressWarnings(as.integer(cores[[1L]]))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  cores
}

.warplda_thread_cap <- function(cores) {
  opt_cap <- getOption("craftgrn.warplda.max_threads", NULL)
  if (!is.null(opt_cap) && length(opt_cap) > 0L && !is.na(opt_cap[[1L]])) {
    out <- suppressWarnings(as.integer(opt_cap[[1L]]))
    if (is.finite(out) && out > 0L) return(min(cores, out))
  }

  env_cap <- Sys.getenv("CRAFTGRN_WARPLDA_MAX_THREADS", unset = "")
  if (nzchar(env_cap)) {
    out <- suppressWarnings(as.integer(env_cap))
    if (is.finite(out) && out > 0L) return(min(cores, out))
  }

  cores
}

.warplda_default_threads <- function(n_threads = NULL) {
  if (!is.null(n_threads) && length(n_threads) > 0L && !is.na(n_threads[[1L]])) {
    out <- suppressWarnings(as.integer(n_threads[[1L]]))
    if (is.finite(out) && out > 0L) return(out)
  }
  .warplda_thread_cap(.warplda_available_threads())
}

.warplda_resolve_threads_per_model <- function(threads_per_model = NULL,
                                               workers = 1L) {
  workers <- suppressWarnings(as.integer(workers[[1L]]))
  if (!is.finite(workers) || workers < 1L) workers <- 1L
  is_null <- is.null(threads_per_model) ||
    length(threads_per_model) == 0L ||
    is.na(threads_per_model[[1L]])
  if (isTRUE(is_null) && workers > 1L) {
    return(as.integer(max(1L, floor(.warplda_default_threads() / workers))))
  }
  .warplda_default_threads(threads_per_model)
}

.module3_memory_policy <- function(memory_safety = getOption("craftgrn.module3.memory_safety", "strict"),
                                   memory_max_fraction = getOption("craftgrn.module3.memory_max_fraction", 0.8)) {
  memory_safety <- match.arg(as.character(memory_safety)[[1L]], c("strict", "adaptive", "off"))
  memory_max_fraction <- suppressWarnings(as.numeric(memory_max_fraction)[[1L]])
  if (!is.finite(memory_max_fraction) || memory_max_fraction < 0.1 || memory_max_fraction > 0.9) {
    .log_abort("memory_max_fraction must be between 0.1 and 0.9.")
  }
  list(mode = memory_safety, max_fraction = memory_max_fraction)
}

.module3_memory_preflight <- function(estimated_bytes,
                                      stage,
                                      memory_safety = getOption("craftgrn.module3.memory_safety", "strict"),
                                      memory_max_fraction = getOption("craftgrn.module3.memory_max_fraction", 0.8),
                                      available_bytes = .available_memory_bytes(),
                                      verbose = TRUE) {
  policy <- .module3_memory_policy(memory_safety, memory_max_fraction)
  estimated_bytes <- suppressWarnings(as.numeric(estimated_bytes)[[1L]])
  available_bytes <- suppressWarnings(as.numeric(available_bytes)[[1L]])
  budget_bytes <- if (is.finite(available_bytes) && available_bytes > 0) {
    available_bytes * policy$max_fraction
  } else {
    NA_real_
  }
  result <- list(
    stage = as.character(stage)[[1L]],
    mode = policy$mode,
    max_fraction = policy$max_fraction,
    estimated_bytes = estimated_bytes,
    available_bytes = available_bytes,
    budget_bytes = budget_bytes,
    allowed = TRUE
  )
  if (identical(policy$mode, "off")) return(result)
  if (!is.finite(available_bytes) || available_bytes <= 0) {
    result$allowed <- !identical(policy$mode, "strict")
    if (!result$allowed) {
      .log_abort("Module 3 strict memory preflight could not measure available RAM for stage: {stage}.")
    }
    .log_warn("Module 3 could not measure available RAM for stage {stage}; continuing in adaptive mode with one worker.")
    return(result)
  }
  if (!is.finite(estimated_bytes) || estimated_bytes <= 0) {
    result$allowed <- !identical(policy$mode, "strict")
    if (!result$allowed) {
      .log_abort("Module 3 strict memory preflight has no valid peak estimate for stage: {stage}.")
    }
    return(result)
  }
  if (isTRUE(verbose)) {
    .log_inform(
      "Module 3 memory preflight [{stage}]: estimate {(.format_bytes(estimated_bytes))}; available {(.format_bytes(available_bytes))}; budget {(.format_bytes(budget_bytes))} ({round(100 * policy$max_fraction)}%)."
    )
  }
  if (estimated_bytes > budget_bytes) {
    result$allowed <- FALSE
    .log_abort(c(
      "Module 3 memory preflight refused an unsafe allocation.",
      i = paste0("Stage: ", stage),
      i = paste0("Estimated peak: ", .format_bytes(estimated_bytes)),
      i = paste0("Current safe budget: ", .format_bytes(budget_bytes))
    ))
  }
  result
}

.cap_warplda_token_counts <- function(counts, max_tokens = 1.9e9) {
  counts <- suppressWarnings(as.numeric(counts))
  counts[!is.finite(counts) | counts < 0] <- 0
  raw_tokens <- sum(counts)
  if (raw_tokens <= max_tokens) {
    return(list(counts = as.integer(counts), raw_tokens = raw_tokens, tokens = raw_tokens, scale_factor = 1))
  }
  positive <- counts > 0
  lo <- 0
  hi <- 1
  for (i in seq_len(40L)) {
    mid <- (lo + hi) / 2
    candidate <- ifelse(positive, pmax(1, floor(counts * mid)), 0)
    if (sum(candidate) > max_tokens) hi <- mid else lo <- mid
  }
  adjusted <- ifelse(positive, pmax(1, floor(counts * lo)), 0)
  tokens <- sum(adjusted)
  if (tokens > max_tokens) .log_abort("Unable to cap WarpLDA expanded token counts below the supported limit.")
  list(
    counts = as.integer(adjusted),
    raw_tokens = raw_tokens,
    tokens = tokens,
    scale_factor = lo
  )
}

.warplda_memory_estimate <- function(dtm,
                                     K,
                                     n_threads = NULL,
                                     safety_factor = 2.5) {
  .assert_pkg("Matrix")
  if (!inherits(dtm, "dgCMatrix")) dtm <- methods::as(dtm, "dgCMatrix")

  K <- suppressWarnings(as.integer(K[[1L]]))
  if (!is.finite(K) || K < 2L) .log_abort("K must be an integer >= 2.")

  n_threads <- .warplda_default_threads(n_threads)
  safety_factor <- suppressWarnings(as.numeric(safety_factor[[1L]]))
  if (!is.finite(safety_factor) || safety_factor < 1) safety_factor <- 2.5

  n_doc <- as.numeric(nrow(dtm))
  n_word <- as.numeric(ncol(dtm))
  n_nz <- as.numeric(Matrix::nnzero(dtm))
  n_tokens <- as.numeric(sum(dtm))

  bytes_int <- 4
  bytes_double <- 8

  # Native WarpLDA simultaneously retains three token-index vectors and three
  # token-topic work vectors. The peak factor covers allocator and conversion
  # overhead observed during production fits.
  token_bytes <- n_tokens * bytes_int * 6
  sparse_index_bytes <- ((n_doc + 1) + (n_word + 1) + (n_nz * 3)) * bytes_int
  count_bytes <- ((n_doc * K) + (n_word * K) + K) * bytes_int
  output_bytes <- ((n_doc * K) + (n_word * K) + (n_doc * K)) * bytes_double
  thread_bytes <- as.numeric(n_threads) * K * bytes_double * 2
  dtm_bytes <- as.numeric(utils::object.size(dtm))
  base_bytes <- token_bytes + sparse_index_bytes + count_bytes + output_bytes + thread_bytes + dtm_bytes
  estimated_peak_bytes <- ceiling(base_bytes * safety_factor)

  list(
    n_doc = n_doc,
    n_terms = n_word,
    n_nonzero = n_nz,
    n_tokens = n_tokens,
    K = as.integer(K),
    n_threads = as.integer(n_threads),
    base_bytes = as.numeric(base_bytes),
    estimated_peak_bytes = as.numeric(estimated_peak_bytes),
    estimated_peak_gb = as.numeric(estimated_peak_bytes) / 1024^3,
    label = .format_bytes(estimated_peak_bytes)
  )
}

.warplda_reference_init <- function(dtm, K, seed) {
  .assert_pkg("Matrix")
  if (!inherits(dtm, "dgCMatrix")) dtm <- methods::as(dtm, "dgCMatrix")
  n_tokens <- as.numeric(sum(dtm))
  if (!is.finite(n_tokens) || n_tokens <= 0) {
    .log_abort("warp_ref requires a DTM with positive integer token counts.")
  }
  if (abs(n_tokens - round(n_tokens)) > 1e-7) {
    .log_abort("warp_ref requires positive integer-like DTM counts.")
  }
  if (n_tokens > .Machine$integer.max) {
    .log_abort("warp_ref expanded token count exceeds the supported integer limit.")
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed[[1L]]))
  ref_seeds <- as.integer(stats::runif(2L, 1, 2^31 - 1))
  set.seed(ref_seeds[[1L]])
  n_tokens_i <- as.integer(round(n_tokens))
  list(
    topic = sample.int(n = as.integer(K), size = n_tokens_i, replace = TRUE) - 1L,
    proposal = sample.int(n = as.integer(K), size = n_tokens_i, replace = TRUE) - 1L,
    seeds = ref_seeds
  )
}

# Fit ONE native WarpLDA-compatible model; returns theta/phi + metrics
fit_warplda_one <- function(dtm,
                            K,
                            iterations = 2000L,
                            alpha = NULL,
                            beta = NULL,
                            seed = 1L,
                            convergence_tol = 1e-3,
                            n_check_convergence = 10L,
                            n_iter_inference = 10L,
                            n_threads = NULL,
                            sampler = c("warp_omp", "warp_ref", "warp_mh", "gibbs_sync"),
                            progressbar = interactive()) {
  .assert_pkg("cli")
  .assert_pkg("Matrix")

  sampler <- match.arg(sampler)
  if (!inherits(dtm, "dgCMatrix")) dtm <- methods::as(dtm, "dgCMatrix")
  if (is.null(alpha)) alpha <- 50 / as.numeric(K)
  if (is.null(beta)) beta <- 1 / as.numeric(K)
  beta <- as.numeric(beta)
  n_threads <- .warplda_default_threads(n_threads)
  ref_init <- NULL
  if (identical(sampler, "warp_ref")) {
    ref_init <- .warplda_reference_init(dtm, K = K, seed = seed)
  }
  memory_estimate <- .warplda_memory_estimate(dtm, K = K, n_threads = n_threads)

  if (isTRUE(progressbar)) {
    .log_inform("Native WarpLDA does not use an interactive progress bar; use verbose Module 3 logs for progress.")
  }

  fit <- .craftgrn_warplda_fit_cpp(
    dtm = dtm,
    K = as.integer(K),
    iterations = as.integer(iterations),
    alpha = as.numeric(alpha),
    beta = as.numeric(beta),
    seed = as.integer(seed),
    convergence_tol = as.numeric(convergence_tol),
    n_check_convergence = as.integer(n_check_convergence),
    n_iter_inference = as.integer(n_iter_inference),
    n_threads = as.integer(n_threads),
    sampler = sampler,
    ref_topic = if (is.null(ref_init)) NULL else ref_init$topic,
    ref_proposal = if (is.null(ref_init)) NULL else ref_init$proposal,
    ref_seeds = if (is.null(ref_init)) NULL else ref_init$seeds
  )
  trace <- data.table::as.data.table(fit$loglik_trace)
  relative_change <- NA_real_
  if (nrow(trace) >= 2L) {
    previous <- as.numeric(trace$loglikelihood[[nrow(trace) - 1L]])
    current <- as.numeric(trace$loglikelihood[[nrow(trace)]])
    relative_change <- abs(current - previous) / max(1, abs(previous))
  }
  converged <- is.finite(relative_change) && convergence_tol >= 0 &&
    relative_change < convergence_tol
  theta <- as.matrix(fit$theta)
  phi <- as.matrix(fit$phi)

  if (!is.null(rownames(dtm)) && nrow(theta) == nrow(dtm)) rownames(theta) <- rownames(dtm)
  if (!is.null(colnames(dtm)) && ncol(phi) == ncol(dtm)) colnames(phi) <- colnames(dtm)
  if (is.null(colnames(theta))) colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  if (is.null(rownames(phi))) rownames(phi) <- colnames(theta)

  list(
    model = list(
      backend = "craftgrn_native_warplda",
      loglik_trace = fit$loglik_trace,
      threads = as.integer(fit$threads),
      sampler = sampler
    ),
    K = as.integer(K),
    iterations = as.integer(fit$iterations),
    iterations_requested = as.integer(iterations),
    alpha = alpha,
    beta = beta,
    seed = as.integer(seed),
    theta = theta,
    phi = phi,
    metrics = list(
      n_tokens = as.numeric(fit$n_tokens),
      perplexity = as.numeric(fit$perplexity),
      loglik_approx = as.numeric(fit$loglik),
      backend = "craftgrn_native_warplda",
      threads = as.integer(fit$threads),
      sampler = sampler,
      iterations_requested = as.integer(iterations),
      iterations_completed = as.integer(fit$iterations),
      convergence_tolerance = as.numeric(convergence_tol),
      convergence_check_interval = as.integer(n_check_convergence),
      final_check_relative_change = as.numeric(relative_change),
      converged = isTRUE(converged),
      estimated_peak_memory_gb = as.numeric(memory_estimate$estimated_peak_gb),
      estimated_peak_memory = memory_estimate$label
    )
  )
}

.existing_file <- function(path) {
  is.character(path) & !is.na(path) & nzchar(path) & file.exists(path)
}

.warplda_completed_from_cache <- function(K_grid,
                                          existing_metrics,
                                          fit_files_all,
                                          iterations = NULL,
                                          alpha_by_topic = NULL,
                                          alpha = NULL,
                                          beta = NULL,
                                          seed = NULL,
                                          sampler = NULL) {
  sampler_value <- if (is.null(sampler)) NULL else as.character(sampler[[1L]])
  fit_matches <- function(path, K) {
    if (!.existing_file(path)) return(FALSE)
    fit <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(fit)) return(FALSE)
    fit_sampler <- as.character(fit$model$sampler %||% fit$metrics$sampler %||% "gibbs_sync")
    settings_match <- TRUE
    if (!is.null(iterations)) {
      settings_match <- settings_match && isTRUE(all.equal(
        as.integer(fit$iterations_requested %||% fit$iterations),
        as.integer(iterations)
      ))
    }
    if (!is.null(alpha_by_topic)) {
      expected_alpha <- if (isTRUE(alpha_by_topic)) 50 / as.numeric(K) else as.numeric(alpha)
      settings_match <- settings_match && isTRUE(all.equal(as.numeric(fit$alpha), expected_alpha))
      expected_beta <- if (is.null(beta)) 1 / as.numeric(K) else as.numeric(beta)
      settings_match <- settings_match && isTRUE(all.equal(as.numeric(fit$beta), expected_beta))
    }
    if (!is.null(seed)) {
      settings_match <- settings_match && identical(as.integer(fit$seed), as.integer(seed))
    }
    settings_match && (is.null(sampler_value) || identical(fit_sampler, sampler_value))
  }
  done_from_file <- vapply(
    seq_along(K_grid),
    function(i) fit_matches(fit_files_all[[i]], K_grid[[i]]),
    logical(1L)
  )
  done_from_metrics <- rep(FALSE, length(K_grid))
  if (nrow(existing_metrics) && "K" %in% names(existing_metrics)) {
    metrics_dt <- data.table::as.data.table(existing_metrics)
    for (idx in seq_along(K_grid)) {
      k_value <- as.integer(K_grid[[idx]])
      rows <- metrics_dt[as.integer(metrics_dt[["K"]]) == k_value]
      if (!nrow(rows)) next
      if (!is.null(sampler_value)) {
        if ("sampler" %in% names(rows)) {
          rows <- rows[as.character(rows[["sampler"]]) == sampler_value]
        } else if (!identical(sampler_value, "gibbs_sync")) {
          rows <- rows[0L]
        }
      }
      if (!nrow(rows)) next
      metric_paths <- character()
      if ("fit_file" %in% names(rows)) {
        metric_paths <- as.character(rows[["fit_file"]])
      }
      expected_path <- fit_files_all[[idx]]
      metric_paths <- unique(c(metric_paths, expected_path))
      done_from_metrics[[idx]] <- any(vapply(
        metric_paths,
        fit_matches,
        logical(1L),
        K = k_value
      ))
    }
  }
  list(done_from_file = done_from_file, done_from_metrics = done_from_metrics)
}

# Fit WarpLDA models across K grid (cisTopic-like runWarpLDAModels)
run_warplda_models <- function(dtm,
                               K_grid,
                               iterations = 2000L,
                               alpha_by_topic = TRUE,
                               alpha = NULL,
                               beta = 0.1,
                               seed = 123,
                               save_tmp_dir = NULL,
                               workers = 1L,
                               threads_per_model = NULL,
                               sampler = c("warp_omp", "warp_ref", "warp_mh", "gibbs_sync"),
                               metrics_file = NULL,
                               memory_safety = getOption("craftgrn.module3.memory_safety", "strict"),
                               memory_max_fraction = getOption("craftgrn.module3.memory_max_fraction", 0.8),
                               verbose = TRUE) {
  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("Matrix")

  sampler <- match.arg(sampler)
  K_grid <- as.integer(K_grid)
  K_grid <- sort(unique(K_grid[is.finite(K_grid) & K_grid > 1L]))
  if (!length(K_grid)) .log_abort("K_grid must include integers > 1.")

  if (!is.null(save_tmp_dir)) dir.create(save_tmp_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(metrics_file)) {
    dir.create(dirname(metrics_file), recursive = TRUE, showWarnings = FALSE)
  }

  workers <- suppressWarnings(as.integer(workers))
  if (!is.finite(workers) || workers < 1L) {
    workers <- 1L
  }
  workers <- min(workers, length(K_grid))
  threads_per_model <- .warplda_resolve_threads_per_model(threads_per_model, workers = workers)

  existing_metrics <- data.table::data.table()
  if (!is.null(metrics_file) && file.exists(metrics_file)) {
    existing_metrics <- data.table::fread(metrics_file)
  }

  fit_file_for_k <- function(K) {
    if (is.null(save_tmp_dir)) {
      return(NA_character_)
    }
    file.path(save_tmp_dir, sprintf("fit_K%d_%s.rds", as.integer(K), sampler))
  }

  metric_from_fit <- function(K, fit_file) {
    fit <- readRDS(fit_file)
    fit_sampler <- if (!is.null(fit$metrics$sampler)) as.character(fit$metrics$sampler) else "gibbs_sync"
    metric_value <- function(name, default = NA) {
      value <- fit$metrics[[name]]
      if (is.null(value) || !length(value)) default else value[[1L]]
    }
    data.table::data.table(
      K = as.integer(fit$K),
      perplexity = as.numeric(fit$metrics$perplexity),
      loglik = as.numeric(fit$metrics$loglik_approx),
      n_tokens = as.numeric(fit$metrics$n_tokens),
      sampler = fit_sampler,
      alpha = as.numeric(fit$alpha),
      beta = as.numeric(fit$beta),
      seed = as.integer(fit$seed),
      threads = as.integer(metric_value("threads", fit$model$threads)),
      iterations_requested = as.integer(metric_value(
        "iterations_requested",
        if (!is.null(fit$iterations_requested)) fit$iterations_requested else fit$iterations
      )),
      iterations_completed = as.integer(metric_value("iterations_completed", fit$iterations)),
      convergence_tolerance = as.numeric(metric_value("convergence_tolerance")),
      final_check_relative_change = as.numeric(metric_value("final_check_relative_change")),
      converged = as.logical(metric_value("converged")),
      fit_file = fit_file
    )
  }

  fit_one_k <- function(K) {
    a <- if (isTRUE(alpha_by_topic)) 50 / as.numeric(K) else alpha
    beta_display <- if (is.null(beta)) "1/K" else signif(as.numeric(beta), 4)
    fit_file <- fit_file_for_k(K)
    if (!is.na(fit_file) && file.exists(fit_file)) {
      return(metric_from_fit(K, fit_file))
    }
    if (isTRUE(verbose)) {
      .log_inform("WarpLDA: fitting K={K}, iterations={iterations}, alpha={signif(a, 4)}, beta={beta_display}, seed={seed}")
    }
    one_memory <- .warplda_memory_estimate(dtm, K = K, n_threads = threads_per_model)
    .module3_memory_preflight(
      one_memory$estimated_peak_bytes,
      stage = paste0("WarpLDA K=", K),
      memory_safety = memory_safety,
      memory_max_fraction = memory_max_fraction,
      verbose = verbose
    )
    fit <- tryCatch(
      fit_warplda_one(
        dtm,
        K = K,
        iterations = iterations,
        alpha = a,
        beta = beta,
        seed = seed,
        n_threads = threads_per_model,
        sampler = sampler
      ),
      error = function(e) {
        .log_abort(c(
          "WarpLDA fit failed.",
          i = paste0("K=", K, ", iterations=", iterations, ", alpha=", signif(a, 4), ", beta=", beta_display, ", seed=", seed),
          i = paste0("dtm dims=", nrow(dtm), "x", ncol(dtm), ", nnzero=", Matrix::nnzero(dtm), ", tokens=", sum(dtm))
        ), parent = e)
      }
    )
    if (isTRUE(verbose) &&
        is.finite(fit$metrics$final_check_relative_change) &&
        !isTRUE(fit$metrics$converged)) {
      .log_warn(
        "WarpLDA K={K} exhausted {fit$metrics$iterations_completed} iteration(s) without meeting its convergence tolerance; final checked relative log-likelihood change was {signif(fit$metrics$final_check_relative_change, 4)}."
      )
    }
    if (!is.na(fit_file)) {
      saveRDS(fit, fit_file)
    }
    out <- data.table::data.table(
      K = as.integer(fit$K),
      perplexity = as.numeric(fit$metrics$perplexity),
      loglik = as.numeric(fit$metrics$loglik_approx),
      n_tokens = as.numeric(fit$metrics$n_tokens),
      sampler = sampler,
      alpha = as.numeric(fit$alpha),
      beta = as.numeric(fit$beta),
      seed = as.integer(fit$seed),
      threads = as.integer(fit$metrics$threads),
      iterations_requested = as.integer(fit$metrics$iterations_requested),
      iterations_completed = as.integer(fit$metrics$iterations_completed),
      convergence_tolerance = as.numeric(fit$metrics$convergence_tolerance),
      final_check_relative_change = as.numeric(fit$metrics$final_check_relative_change),
      converged = as.logical(fit$metrics$converged),
      fit_file = fit_file
    )
    rm(fit)
    invisible(gc())
    out
  }

  fit_files_all <- vapply(K_grid, fit_file_for_k, character(1))
  cache_status <- .warplda_completed_from_cache(
    K_grid,
    existing_metrics,
    fit_files_all,
    iterations = iterations,
    alpha_by_topic = alpha_by_topic,
    alpha = alpha,
    beta = beta,
    seed = seed,
    sampler = sampler
  )
  done_from_file <- cache_status$done_from_file
  done_from_metrics <- cache_status$done_from_metrics
  stale_fit_files <- fit_files_all[.existing_file(fit_files_all) & !done_from_file]
  if (length(stale_fit_files)) unlink(stale_fit_files, force = TRUE)
  todo <- K_grid[!(done_from_file | done_from_metrics)]
  max_k_memory <- .warplda_memory_estimate(dtm, K = max(K_grid), n_threads = threads_per_model)
  concurrent_peak_bytes <- max_k_memory$estimated_peak_bytes * max(1L, workers)
  available_bytes <- .available_memory_bytes()
  per_model_memory_label <- .format_bytes(max_k_memory$estimated_peak_bytes)
  concurrent_memory_label <- .format_bytes(concurrent_peak_bytes)
  available_memory_label <- .format_bytes(available_bytes)

  if (isTRUE(verbose)) {
    .log_inform(
      "WarpLDA: dtm dims = {nrow(dtm)} x {ncol(dtm)}, nnzero = {Matrix::nnzero(dtm)}, tokens = {sum(dtm)}; completed K = {sum(done_from_file | done_from_metrics)}/{length(K_grid)}; workers = {workers}; threads/model = {threads_per_model}"
    )
    .log_inform("WarpLDA sampler = {sampler}.")
    .log_inform(
      "WarpLDA memory estimate: peak per model <= {per_model_memory_label} at max K={max(K_grid)}; concurrent peak <= {concurrent_memory_label}."
    )
    if (is.finite(available_bytes)) {
      .log_inform("WarpLDA memory available now: {available_memory_label}.")
      .module3_memory_preflight(
        concurrent_peak_bytes,
        stage = "WarpLDA model grid",
        memory_safety = memory_safety,
        memory_max_fraction = memory_max_fraction,
        available_bytes = available_bytes,
        verbose = FALSE
      )
    }
  }

  metrics_list <- list()
  if (any(done_from_file)) {
    metrics_list <- c(metrics_list, lapply(K_grid[done_from_file], function(K) metric_from_fit(K, fit_file_for_k(K))))
  }
  if (nrow(existing_metrics)) {
    keep <- data.table::as.data.table(existing_metrics)[K %in% K_grid]
    if (nrow(keep)) {
      if (!"fit_file" %in% names(keep)) {
        keep[, fit_file := fit_file_for_k(K)]
      }
      keep <- keep[.existing_file(fit_file)]
      if (nrow(keep)) {
        if (!"sampler" %in% names(keep)) keep[, sampler := "gibbs_sync"]
        metric_columns <- intersect(
          c(
            "K", "perplexity", "loglik", "n_tokens", "sampler",
            "alpha", "beta", "seed", "threads",
            "iterations_requested", "iterations_completed",
            "convergence_tolerance", "final_check_relative_change",
            "converged", "fit_file"
          ),
          names(keep)
        )
        metrics_list <- c(metrics_list, list(keep[, metric_columns, with = FALSE]))
      }
    }
  }

  if (length(todo)) {
    if (workers > 1L) {
      old_dt_threads <- data.table::getDTthreads()
      data.table::setDTthreads(1L)
      on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
      fit_new <- parallel::mclapply(
        todo,
        fit_one_k,
        mc.cores = workers,
        mc.preschedule = FALSE
      )
    } else {
      fit_new <- lapply(todo, fit_one_k)
    }
    metrics_list <- c(metrics_list, fit_new)
  }

  metrics_tbl <- data.table::rbindlist(metrics_list, use.names = TRUE, fill = TRUE)
  if (!nrow(metrics_tbl)) {
    .log_abort("WarpLDA did not produce any fitted models.")
  }
  metrics_tbl <- unique(metrics_tbl, by = "K")
  data.table::setorder(metrics_tbl, K)
  missing_k <- setdiff(K_grid, metrics_tbl$K)
  if (length(missing_k)) {
    .log_abort("WarpLDA missing fitted K value(s): {paste(missing_k, collapse = ', ')}")
  }
  if (!is.null(metrics_file)) {
    data.table::fwrite(metrics_tbl, metrics_file)
  }
  fit_files <- stats::setNames(as.character(metrics_tbl$fit_file), as.character(metrics_tbl$K))
  public_metric_columns <- intersect(
    c(
      "K", "perplexity", "loglik", "n_tokens", "sampler",
      "alpha", "beta", "seed", "threads",
      "iterations_requested", "iterations_completed",
      "convergence_tolerance", "final_check_relative_change", "converged"
    ),
    names(metrics_tbl)
  )
  list(
    fit_files = fit_files,
    metrics = metrics_tbl[, public_metric_columns, with = FALSE]
  )
}

# =============================================================================
# 6) Step 6: Model selection + cisTopic-like plots (loglik, 2nd-derivative, perplexity)
# =============================================================================

.second_derivative <- function(x) {
  # discrete 2nd derivative: d2[i] aligns to x[i] (endpoints NA)
  x <- .safe_num(x)
  n <- length(x)
  d2 <- rep(NA_real_, n)
  if (n < 3) return(d2)
  for (i in 2:(n - 1)) d2[i] <- x[i + 1] - 2 * x[i] + x[i - 1]
  d2
}

select_model_indices <- function(metrics_tbl) {
  .assert_pkg("data.table")
  m <- data.table::as.data.table(metrics_tbl)
  .assert_has_cols(m, c("K","loglik","perplexity"), context = "select_model_indices")

  data.table::setorder(m, K)

  m[, d2_loglik := .second_derivative(loglik)]

  idx_maxlik <- if (any(is.finite(m$loglik))) which.max(ifelse(is.finite(m$loglik), m$loglik, -Inf)) else NA_integer_
  idx_perp <- if (any(is.finite(m$perplexity))) which.min(ifelse(is.finite(m$perplexity), m$perplexity, Inf)) else NA_integer_
  idx_knee <- if (any(is.finite(m$d2_loglik))) which.max(ifelse(is.finite(m$d2_loglik), m$d2_loglik, -Inf)) else NA_integer_

  list(
    table = m,
    idx = list(maximum = idx_maxlik, perplexity = idx_perp, derivative = idx_knee)
  )
}

plot_model_selection_cistopic <- function(metrics_tbl, out_pdf, title_prefix = NULL) {
  .assert_pkg("data.table")

  sel <- select_model_indices(metrics_tbl)
  m <- sel$table
  idx <- sel$idx
  data.table::setorder(m, K)

  safe_panel <- function(x, y, xlab, ylab, main) {
    if (!length(y) || !any(is.finite(y))) {
      plot.new()
      title(main = main)
      text(0.5, 0.5, labels = "No finite values to plot")
      return(invisible(NULL))
    }
    plot(x, y, type = "b", pch = 16, xlab = xlab, ylab = ylab, main = main, col = "grey20")
  }

  top_n_idx <- function(values, n = 3L, decreasing = TRUE) {
    idx <- which(is.finite(values))
    if (!length(idx)) return(integer(0))
    ord <- order(values[idx], decreasing = decreasing)
    idx[ord[seq_len(min(n, length(ord)))]]
  }

  top_loglik <- top_n_idx(m$loglik, n = 3L, decreasing = TRUE)
  top_d2 <- top_n_idx(m$d2_loglik, n = 3L, decreasing = TRUE)
  top_perp <- top_n_idx(m$perplexity, n = 3L, decreasing = FALSE)

  styles <- list(
    loglik = list(col = "#1f77b4", pch = 16),
    d2 = list(col = "#ff7f0e", pch = 17),
    perp = list(col = "#2ca02c", pch = 15)
  )

  .assign_label_pos <- function(x, y, idx) {
    pos <- rep(3, length(idx))
    if (length(idx) < 2) return(pos)
    ord <- order(x[idx], y[idx])
    pos_alt <- c(1, 3, 2, 4)
    pos[ord] <- pos_alt[(seq_along(ord) - 1L) %% length(pos_alt) + 1L]
    pos
  }

  overlay_marks <- function(x, y) {
    add_one <- function(idx, style, label_offset = 0) {
      if (!length(idx)) return(invisible(NULL))
      idx <- idx[is.finite(y[idx])]
      if (!length(idx)) return(invisible(NULL))
      points(x[idx], y[idx], pch = style$pch, cex = 1.3, col = style$col)
      pos <- .assign_label_pos(x, y, idx)
      text(x[idx], y[idx], labels = m$K[idx], pos = pos, offset = 0.6 + label_offset, cex = 0.7, col = style$col)
    }
    add_one(top_loglik, styles$loglik, label_offset = 0.2)
    add_one(top_d2, styles$d2, label_offset = 0.5)
    add_one(top_perp, styles$perp, label_offset = 0.8)
  }

  grDevices::pdf(out_pdf, width = 11, height = 4)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); grDevices::dev.off() }, add = TRUE)

  par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.2, 1))

  prefix <- if (!is.null(title_prefix) && nzchar(title_prefix)) paste0(title_prefix, " | ") else ""

  # (1) loglik
  safe_panel(m$K, m$loglik, "Number of topics", "log P(D|M,T)", paste0(prefix, "Model selection: loglik"))
  if (any(is.finite(m$loglik))) overlay_marks(m$K, m$loglik)

  # (2) 2nd derivative on loglik
  safe_panel(m$K, m$d2_loglik, "Number of topics",
             "2nd derivative on the log-likelihood", paste0(prefix, "Model selection: 2nd derivative"))
  if (any(is.finite(m$d2_loglik))) overlay_marks(m$K, m$d2_loglik)

  # (3) perplexity (lower is better)
  safe_panel(m$K, m$perplexity, "Number of topics", "Perplexity", paste0(prefix, "Model selection: perplexity"))
  if (any(is.finite(m$perplexity))) overlay_marks(m$K, m$perplexity)

  legend(
    "topright",
    legend = c("Top3 loglik", "Top3 d2 loglik", "Top3 perplexity"),
    col = c(styles$loglik$col, styles$d2$col, styles$perp$col),
    pch = c(styles$loglik$pch, styles$d2$pch, styles$perp$pch),
    bty = "n",
    cex = 0.8
  )

  invisible(sel)
}

# =============================================================================
# 7) Topicterm scoring (cisTopic-like NormTop) + binarization (Gamma-fit or TopN)
# =============================================================================

.module3_topic_model_family <- function(x) {
  value <- tolower(paste(as.character(x), collapse = " "))
  if (grepl("multivi", value, fixed = TRUE)) return("multivi")
  if (grepl("vae_mlp", value, fixed = TRUE) || grepl("vae-mlp", value, fixed = TRUE)) {
    return("vae_mlp")
  }
  if (grepl("warplda", value, fixed = TRUE) || grepl("lda", value, fixed = TRUE)) {
    return("lda")
  }
  "default"
}

.module3_default_gammafit_thrP <- function(model_family) {
  family <- .module3_topic_model_family(model_family)
  switch(
    family,
    lda = 0.70,
    multivi = 0.50,
    vae_mlp = 0.80,
    0.80
  )
}

.gamma_moments <- function(x) {
  x <- .safe_num(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) return(NULL)
  m <- mean(x)
  v <- stats::var(x)
  if (!is.finite(m) || !is.finite(v) || v <= 0) return(NULL)
  shape <- m * m / v
  rate <- m / v
  if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) return(NULL)
  list(shape = shape, rate = rate)
}

score_terms_normtop <- function(phi,
                                method = c("normtop_specificity", "rowmax_phi")) {
  method <- match.arg(method)
  phi <- as.matrix(phi)
  rn <- rownames(phi)
  cn <- colnames(phi)
  phi[!is.finite(phi) | phi < 0] <- 0
  rs <- rowSums(phi)
  rs[!is.finite(rs) | rs <= 0] <- 1
  phi_prob <- phi / rs
  phi_prob[!is.finite(phi_prob) | phi_prob < 0] <- 0

  if (identical(method, "normtop_specificity")) {
    eps <- 1e-12
    log_phi <- log(phi_prob + eps)
    mean_log_phi <- colMeans(log_phi, na.rm = TRUE)
    sc <- phi_prob * sweep(log_phi, 2L, mean_log_phi, "-")
    sc[!is.finite(sc) | sc < 0] <- 0
    mx <- apply(sc, 1L, max, na.rm = TRUE)
    mx[!is.finite(mx) | mx <= 0] <- 1
    sc <- sc / mx
    sc[!is.finite(sc)] <- 0
    sc[sc < 0] <- 0
    sc[sc > 1] <- 1
    rownames(sc) <- rn
    colnames(sc) <- cn
    return(sc)
  }

  # Legacy behavior: normalize each topic by its largest term probability.
  mx <- apply(phi, 1, max, na.rm = TRUE)
  mx[!is.finite(mx) | mx <= 0] <- 1
  sc <- phi / mx
  sc[!is.finite(sc)] <- 0
  sc[sc < 0] <- 0
  sc[sc > 1] <- 1
  rownames(sc) <- rn
  colnames(sc) <- cn
  sc
}

.validate_topic_probability_matrix <- function(x, name, source = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  finite_pos <- is.finite(x) & x > 0
  if (!any(finite_pos)) {
    where <- if (!is.null(source) && nzchar(source)) paste0(" in ", source) else ""
    .log_abort("Invalid topic model {name}{where}: no finite positive values.")
  }
  bad_frac <- mean(!is.finite(x))
  if (is.finite(bad_frac) && bad_frac > 0.5) {
    where <- if (!is.null(source) && nzchar(source)) paste0(" in ", source) else ""
    .log_abort("Invalid topic model {name}{where}: more than half of values are non-finite.")
  }
  x
}

.term_group <- function(term_id) {
  term_id <- as.character(term_id)
  ifelse(grepl("^PEAK:", term_id), "PEAK",
         ifelse(grepl("^GENE:", term_id), "GENE", "OTHER"))
}

.gamma_cutoff_with_fallback <- function(sc, thrP = 0.975, min_terms = 50L) {
  sc <- as.numeric(sc)
  sc[!is.finite(sc)] <- 0
  ord <- order(sc, decreasing = TRUE)
  fit <- .gamma_moments(sc)
  thr <- NA_real_
  if (!is.null(fit)) {
    thr <- stats::qgamma(as.numeric(thrP), shape = fit$shape, rate = fit$rate)
  }
  if (!is.finite(thr) && length(ord)) {
    take_idx <- min(length(ord), as.integer(min_terms))
    thr <- sc[ord[take_idx]]
  }
  if (!is.finite(thr)) thr <- NA_real_
  thr
}

.gammafit_cutoffs_by_termclass <- function(score_mat,
                                           thrP = 0.975,
                                           min_terms = 50L,
                                           gammafit_scope = c("topic_term_group", "global_term_group")) {
  gammafit_scope <- match.arg(gammafit_scope)
  score_mat <- as.matrix(score_mat)
  K <- nrow(score_mat)
  if (!K) return(data.table::data.table(
    topic_num = integer(0),
    peaks_gamma_cutoff = numeric(0),
    gene_gamma_cutoff = numeric(0),
    other_gamma_cutoff = numeric(0),
    gammafit_scope = character(0)
  ))
  terms <- colnames(score_mat)
  if (is.null(terms)) {
    terms <- paste0("term_", seq_len(ncol(score_mat)))
    colnames(score_mat) <- terms
  }
  grp <- .term_group(terms)
  idx_peak <- which(grp == "PEAK")
  idx_gene <- which(grp == "GENE")
  idx_other <- which(grp == "OTHER")
  out <- data.table::data.table(
    topic_num = seq_len(K),
    peaks_gamma_cutoff = NA_real_,
    gene_gamma_cutoff = NA_real_,
    other_gamma_cutoff = NA_real_,
    gammafit_scope = gammafit_scope
  )
  if (identical(gammafit_scope, "global_term_group")) {
    if (length(idx_peak)) {
      out[, peaks_gamma_cutoff := .gamma_cutoff_with_fallback(as.numeric(score_mat[, idx_peak, drop = FALSE]), thrP = thrP, min_terms = min_terms)]
    }
    if (length(idx_gene)) {
      out[, gene_gamma_cutoff := .gamma_cutoff_with_fallback(as.numeric(score_mat[, idx_gene, drop = FALSE]), thrP = thrP, min_terms = min_terms)]
    }
    if (length(idx_other)) {
      out[, other_gamma_cutoff := .gamma_cutoff_with_fallback(as.numeric(score_mat[, idx_other, drop = FALSE]), thrP = thrP, min_terms = min_terms)]
    }
    return(out[])
  }
  for (k in seq_len(K)) {
    sc <- as.numeric(score_mat[k, ])
    sc[!is.finite(sc)] <- 0
    if (length(idx_peak)) {
      out$peaks_gamma_cutoff[k] <- .gamma_cutoff_with_fallback(sc[idx_peak], thrP = thrP, min_terms = min_terms)
    }
    if (length(idx_gene)) {
      out$gene_gamma_cutoff[k] <- .gamma_cutoff_with_fallback(sc[idx_gene], thrP = thrP, min_terms = min_terms)
    }
    if (length(idx_other)) {
      out$other_gamma_cutoff[k] <- .gamma_cutoff_with_fallback(sc[idx_other], thrP = thrP, min_terms = min_terms)
    }
  }
  out
}

.gammafit_stats_for_scores <- function(sc) {
  sc <- .safe_num(sc)
  sc[!is.finite(sc)] <- 0
  pos <- sc[is.finite(sc) & sc > 0]
  fit <- .gamma_moments(sc)
  list(
    term_count = length(sc),
    positive_count = length(pos),
    zero_fraction = if (length(sc)) mean(sc <= 0 | !is.finite(sc)) else NA_real_,
    gamma_shape = if (!is.null(fit)) fit$shape else NA_real_,
    gamma_rate = if (!is.null(fit)) fit$rate else NA_real_
  )
}

.gammafit_diagnostics_by_termclass <- function(score_mat,
                                               topic_terms = NULL,
                                               topic_score_method = c("normtop_specificity", "rowmax_phi"),
                                               thrP = 0.975,
                                               min_terms = 50L,
                                               gammafit_scope = c("topic_term_group", "global_term_group")) {
  topic_score_method <- match.arg(topic_score_method)
  gammafit_scope <- match.arg(gammafit_scope)
  score_mat <- as.matrix(score_mat)
  K <- nrow(score_mat)
  if (!K) return(data.table::data.table())
  terms <- colnames(score_mat)
  if (is.null(terms)) {
    terms <- paste0("term_", seq_len(ncol(score_mat)))
    colnames(score_mat) <- terms
  }
  min_terms <- as.integer(min_terms)
  if (!is.finite(min_terms) || min_terms < 1L) min_terms <- 1L
  term_grp <- .term_group(terms)
  group_idx <- list(
    PEAK = which(term_grp == "PEAK"),
    GENE = which(term_grp == "GENE"),
    OTHER = which(term_grp == "OTHER")
  )
  cut_tbl <- .gammafit_cutoffs_by_termclass(
    score_mat,
    thrP = thrP,
    min_terms = min_terms,
    gammafit_scope = gammafit_scope
  )
  selected_after <- data.table::data.table()
  if (!is.null(topic_terms) && is.data.frame(topic_terms) && nrow(topic_terms)) {
    selected_after <- data.table::as.data.table(topic_terms)
    if (!"topic_num" %in% names(selected_after)) {
      if ("topic" %in% names(selected_after)) {
        selected_after[, topic_num := as.integer(gsub("^Topic", "", as.character(topic)))]
      } else {
        selected_after[, topic_num := NA_integer_]
      }
    }
    if ("in_topic" %in% names(selected_after)) {
      selected_after <- selected_after[.as_logical_flag(in_topic)]
    }
    selected_after[, term_group := .term_group(term_id)]
    selected_after <- selected_after[
      is.finite(topic_num),
      .(selected_after_min_terms = .N),
      by = .(topic_num, term_group)
    ]
  }

  global_stats <- NULL
  if (identical(gammafit_scope, "global_term_group")) {
    global_stats <- lapply(group_idx, function(idx) {
      if (!length(idx)) return(.gammafit_stats_for_scores(numeric()))
      .gammafit_stats_for_scores(as.numeric(score_mat[, idx, drop = FALSE]))
    })
  }

  out <- vector("list", K * length(group_idx))
  pos <- 0L
  for (k in seq_len(K)) {
    sc <- as.numeric(score_mat[k, ])
    sc[!is.finite(sc)] <- 0
    for (grp_name in names(group_idx)) {
      idx <- group_idx[[grp_name]]
      pos <- pos + 1L
      cutoff_col <- switch(
        grp_name,
        PEAK = "peaks_gamma_cutoff",
        GENE = "gene_gamma_cutoff",
        OTHER = "other_gamma_cutoff"
      )
      cutoff <- cut_tbl[[cutoff_col]][k]
      local_stats <- if (identical(gammafit_scope, "global_term_group")) {
        global_stats[[grp_name]]
      } else if (length(idx)) {
        .gammafit_stats_for_scores(sc[idx])
      } else {
        .gammafit_stats_for_scores(numeric())
      }
      selected_by_gamma <- if (length(idx) && is.finite(cutoff)) {
        sum(sc[idx] >= cutoff & sc[idx] > 0, na.rm = TRUE)
      } else {
        0L
      }
      out[[pos]] <- data.table::data.table(
        topic_score_method = topic_score_method,
        gammafit_scope = gammafit_scope,
        thrP = as.numeric(thrP),
        min_terms = as.integer(min_terms),
        topic_num = as.integer(k),
        term_group = grp_name,
        term_count = as.integer(local_stats$term_count),
        positive_count = as.integer(local_stats$positive_count),
        zero_fraction = as.numeric(local_stats$zero_fraction),
        gamma_shape = as.numeric(local_stats$gamma_shape),
        gamma_rate = as.numeric(local_stats$gamma_rate),
        gamma_cutoff = as.numeric(cutoff),
        selected_by_gamma = as.integer(selected_by_gamma),
        selected_after_min_terms = NA_integer_,
        forced_min_terms = FALSE
      )
    }
  }
  out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  if (nrow(selected_after)) {
    out[, selected_after_min_terms := NULL]
    out <- merge(
      out,
      selected_after,
      by = c("topic_num", "term_group"),
      all.x = TRUE,
      sort = FALSE
    )
    out[is.na(selected_after_min_terms), selected_after_min_terms := 0L]
  } else {
    out[, selected_after_min_terms := 0L]
  }
  scope_is_topic <- identical(gammafit_scope, "topic_term_group")
  out[, forced_min_terms := scope_is_topic &
        selected_after_min_terms > selected_by_gamma]
  data.table::setcolorder(out, c(
    "topic_score_method", "gammafit_scope", "thrP", "min_terms",
    "topic_num", "term_group", "term_count", "positive_count",
    "zero_fraction", "gamma_shape", "gamma_rate", "gamma_cutoff",
    "selected_by_gamma", "selected_after_min_terms", "forced_min_terms"
  ))
  out[]
}

binarize_topics <- function(score_mat,
                            method = c("gammafit", "topn"),
                            thrP = 0.975,
                            top_n_terms = 500L,
                            min_terms = 50L,
                            gammafit_scope = c("topic_term_group", "global_term_group")) {
  method <- match.arg(method)
  gammafit_scope <- match.arg(gammafit_scope)
  min_terms <- as.integer(min_terms)
  if (!is.finite(min_terms) || min_terms < 1L) min_terms <- 1L
  score_mat <- as.matrix(score_mat)
  K <- nrow(score_mat)
  terms <- colnames(score_mat)
  if (is.null(terms)) {
    terms <- paste0("term_", seq_len(ncol(score_mat)))
    colnames(score_mat) <- terms
  }
  term_grp <- .term_group(terms)
  idx_peak <- which(term_grp == "PEAK")
  idx_gene <- which(term_grp == "GENE")
  idx_other <- which(term_grp == "OTHER")
  cut_tbl <- NULL
  if (method == "gammafit") {
    cut_tbl <- .gammafit_cutoffs_by_termclass(
      score_mat,
      thrP = thrP,
      min_terms = min_terms,
      gammafit_scope = gammafit_scope
    )
  }

  out_list <- vector("list", K)
  for (k in seq_len(K)) {
    sc <- as.numeric(score_mat[k, ])
    sc[!is.finite(sc)] <- 0
    ord <- order(sc, decreasing = TRUE)

    if (method == "topn") {
      n_take <- min(as.integer(top_n_terms), length(ord))
      in_set <- rep(FALSE, length(sc))
      if (n_take > 0) in_set[ord[seq_len(n_take)]] <- TRUE
    } else {
      in_set <- rep(FALSE, length(sc))
      .apply_group <- function(idx, thr) {
        if (!length(idx)) return(invisible(NULL))
        if (is.finite(thr)) {
          keep <- (sc[idx] >= thr) & (sc[idx] > 0)
        } else {
          keep <- rep(FALSE, length(idx))
        }
        if (identical(gammafit_scope, "topic_term_group") && sum(keep, na.rm = TRUE) < as.integer(min_terms)) {
          ord_local <- idx[order(sc[idx], decreasing = TRUE)]
          take <- min(length(ord_local), as.integer(min_terms))
          if (take > 0) keep <- idx %in% ord_local[seq_len(take)]
        }
        in_set[idx] <<- keep
      }
      .apply_group(idx_peak, cut_tbl$peaks_gamma_cutoff[k])
      .apply_group(idx_gene, cut_tbl$gene_gamma_cutoff[k])
      .apply_group(idx_other, cut_tbl$other_gamma_cutoff[k])
    }

    out_one <- data.frame(
      topic = k,
      term_id = terms,
      score = sc,
      in_topic = as.logical(in_set),
      stringsAsFactors = FALSE
    )
    if (identical(method, "gammafit") && identical(gammafit_scope, "global_term_group")) {
      cutoff <- rep(NA_real_, length(terms))
      cutoff[idx_peak] <- cut_tbl$peaks_gamma_cutoff[k]
      cutoff[idx_gene] <- cut_tbl$gene_gamma_cutoff[k]
      cutoff[idx_other] <- cut_tbl$other_gamma_cutoff[k]
      out_one$term_group <- term_grp
      out_one$gammafit_scope <- gammafit_scope
      out_one$gamma_cutoff <- cutoff
    }
    out_list[[k]] <- out_one
  }

  do.call(rbind, out_list)
}

.assign_topic_terms <- function(raw_phi,
                                score_mat,
                                candidate_terms,
                                method = c("gammafit_maxprob", "max_phi", "gammafit")) {
  method <- match.arg(method)
  raw_phi <- as.matrix(raw_phi)
  score_mat <- as.matrix(score_mat)
  if (!identical(dim(raw_phi), dim(score_mat))) {
    .log_abort("raw_phi and score_mat must have identical dimensions for topic-term assignment.")
  }
  if (is.null(rownames(raw_phi))) rownames(raw_phi) <- paste0("Topic", seq_len(nrow(raw_phi)))
  if (is.null(colnames(raw_phi))) colnames(raw_phi) <- paste0("term_", seq_len(ncol(raw_phi)))
  if (is.null(rownames(score_mat))) rownames(score_mat) <- rownames(raw_phi)
  if (is.null(colnames(score_mat))) colnames(score_mat) <- colnames(raw_phi)
  if (!identical(rownames(raw_phi), rownames(score_mat)) ||
      !identical(colnames(raw_phi), colnames(score_mat))) {
    .log_abort("raw_phi and score_mat dimnames must match for topic-term assignment.")
  }

  out <- data.table::as.data.table(candidate_terms)
  required <- c("topic", "term_id", "score", "in_topic")
  missing <- setdiff(required, names(out))
  if (length(missing)) {
    .log_abort("candidate_terms is missing required columns: {paste(missing, collapse = ', ')}.")
  }
  out[, topic_num := suppressWarnings(as.integer(gsub("^Topic", "", as.character(topic))))]
  out[, term_group := .term_group(term_id)]
  out[, gammafit_candidate := .as_logical_flag(in_topic)]
  out[, `:=`(
    assignment_method = data.table::fcase(
      term_group %in% c("GENE", "PEAK") & identical(method, "gammafit_maxprob"), "gammafit_maxprob",
      term_group %in% c("GENE", "PEAK") & identical(method, "max_phi"), "max_phi",
      default = "gammafit"
    ),
    term_pair_key = NA_character_,
    paired_term_id = NA_character_,
    shared_gammafit_candidate = NA,
    common_candidate_count = NA_integer_,
    term_maxprob_selected = NA,
    term_maxprob_topic = NA_integer_,
    term_maxprob_probability = NA_real_,
    term_maxprob_margin = NA_real_,
    assignment_status = data.table::fcase(
      identical(method, "max_phi"), "assigned_by_raw_phi",
      identical(method, "gammafit"), "gammafit_membership",
      default = "not_applicable"
    )
  )]
  topic_index <- match(out$topic_num, seq_len(nrow(raw_phi)))
  term_index <- match(out$term_id, colnames(raw_phi))
  valid <- is.finite(topic_index) & is.finite(term_index)
  out[, phi := 0]
  out$phi[valid] <- raw_phi[cbind(topic_index[valid], term_index[valid])]
  out[!is.finite(phi) | phi < 0, phi := 0]

  if (identical(method, "max_phi")) {
    phi_clean <- raw_phi
    phi_clean[!is.finite(phi_clean) | phi_clean < 0] <- 0
    max_index <- max.col(t(phi_clean), ties.method = "first")
    max_topic <- stats::setNames(seq_len(nrow(phi_clean))[max_index], colnames(phi_clean))
    eligible <- out$term_group %in% c("GENE", "PEAK") & out$term_id %in% names(max_topic)
    out[eligible, in_topic := topic_num == unname(max_topic[term_id])]
  } else if (identical(method, "gammafit_maxprob")) {
    phi_clean <- raw_phi
    phi_clean[!is.finite(phi_clean) | phi_clean < 0] <- 0
    phi_sum <- colSums(phi_clean)
    phi_sum[!is.finite(phi_sum) | phi_sum <= 0] <- 1
    phi_by_term <- sweep(phi_clean, 2L, phi_sum, "/")
    phi_by_term[!is.finite(phi_by_term) | phi_by_term < 0] <- 0

    term_ids <- colnames(phi_clean)
    gene_keys <- sub("^GENE:", "", term_ids[grepl("^GENE:", term_ids)])
    peak_keys <- sub("^PEAK:", "", term_ids[grepl("^PEAK:", term_ids)])
    pair_keys <- sort(intersect(gene_keys, peak_keys))
    gene_idx <- match(paste0("GENE:", pair_keys), term_ids)
    peak_idx <- match(paste0("PEAK:", pair_keys), term_ids)

    candidate_mat <- matrix(FALSE, nrow = nrow(phi_clean), ncol = ncol(phi_clean))
    valid_candidate <- valid & out$gammafit_candidate
    candidate_mat[cbind(topic_index[valid_candidate], term_index[valid_candidate])] <- TRUE

    select_candidate_max <- function(indices) {
      candidate <- candidate_mat[, indices, drop = FALSE]
      probability <- phi_by_term[, indices, drop = FALSE]
      probability[!candidate | !is.finite(probability)] <- 0
      probability_sum <- colSums(probability)
      has_candidate <- colSums(candidate) > 0L &
        is.finite(probability_sum) & probability_sum > 0
      probability <- sweep(
        probability,
        2L,
        ifelse(has_candidate, probability_sum, 1),
        "/"
      )
      probability[!is.finite(probability) | probability < 0] <- 0
      selected_topic <- max.col(t(probability), ties.method = "first")
      selected_topic[!has_candidate] <- NA_integer_
      selected_probability <- rep(NA_real_, length(indices))
      selected_probability[has_candidate] <- probability[cbind(
        selected_topic[has_candidate],
        which(has_candidate)
      )]
      runner_up_probability <- probability
      if (any(has_candidate)) {
        runner_up_probability[cbind(
          selected_topic[has_candidate],
          which(has_candidate)
        )] <- -Inf
      }
      runner_up <- apply(runner_up_probability, 2L, max, na.rm = TRUE)
      runner_up[!is.finite(runner_up)] <- 0
      list(
        candidate = candidate,
        has_candidate = has_candidate,
        topic = selected_topic,
        probability = selected_probability,
        margin = selected_probability - runner_up
      )
    }

    gene_only_idx <- which(grepl("^GENE:", term_ids))
    gene_only_mode <- length(gene_only_idx) && !any(grepl("^PEAK:", term_ids))

    term_pair_index <- rep(NA_integer_, length(term_ids))
    if (length(pair_keys)) {
      term_pair_index[gene_idx] <- seq_along(pair_keys)
      term_pair_index[peak_idx] <- seq_along(pair_keys)
    }
    row_pair_index <- term_pair_index[term_index]
    paired_rows <- valid & !is.na(row_pair_index) & out$term_group %in% c("GENE", "PEAK")
    unmatched_rows <- out$term_group %in% c("GENE", "PEAK") & !paired_rows

    out[out$term_group %in% c("GENE", "PEAK"), in_topic := FALSE]
    if (gene_only_mode) {
      gene_max <- select_candidate_max(gene_only_idx)
      gene_rows <- valid & out$term_group == "GENE"
      gene_position <- match(term_index[gene_rows], gene_only_idx)
      row_topic <- topic_index[gene_rows]
      selected_topic <- gene_max$topic[gene_position]
      out[gene_rows, `:=`(
        in_topic = !is.na(selected_topic) & row_topic == selected_topic,
        term_pair_key = sub("^GENE:", "", term_id),
        common_candidate_count = colSums(gene_max$candidate)[gene_position],
        term_maxprob_selected = !is.na(selected_topic) & row_topic == selected_topic,
        term_maxprob_topic = selected_topic,
        term_maxprob_probability = gene_max$probability[gene_position],
        term_maxprob_margin = gene_max$margin[gene_position],
        assignment_status = ifelse(
          gene_max$has_candidate[gene_position],
          "assigned_gammafit_maxprob_gene",
          "no_gammafit_candidate"
        )
      )]
    } else {
      out$assignment_status[unmatched_rows] <- "missing_matching_gene_or_peak_term"
    }

    if (length(pair_keys)) {
      shared_candidate <- candidate_mat[, gene_idx, drop = FALSE] &
        candidate_mat[, peak_idx, drop = FALSE]
      common_count <- colSums(shared_candidate)
      # Select each term independently after GammaFit; retain only exact agreement.
      gene_max <- select_candidate_max(gene_idx)
      peak_max <- select_candidate_max(peak_idx)
      has_assignment <- gene_max$has_candidate & peak_max$has_candidate &
        gene_max$topic == peak_max$topic
      selected_topic <- gene_max$topic
      selected_topic[!has_assignment] <- NA_integer_
      pair_status <- data.table::fcase(
        !gene_max$has_candidate & !peak_max$has_candidate, "no_gammafit_candidate",
        !gene_max$has_candidate, "gene_no_gammafit_candidate",
        !peak_max$has_candidate, "peak_no_gammafit_candidate",
        has_assignment, "assigned_gammafit_maxprob_agreement",
        common_count == 0L, "no_shared_gammafit_topic",
        default = "maxprob_topic_disagreement"
      )
      pair_row <- row_pair_index[paired_rows]
      pair_topic <- topic_index[paired_rows]
      is_gene_row <- out$term_group[paired_rows] == "GENE"
      row_max_topic <- ifelse(
        is_gene_row,
        gene_max$topic[pair_row],
        peak_max$topic[pair_row]
      )
      row_max_probability <- ifelse(
        is_gene_row,
        gene_max$probability[pair_row],
        peak_max$probability[pair_row]
      )
      row_max_margin <- ifelse(
        is_gene_row,
        gene_max$margin[pair_row],
        peak_max$margin[pair_row]
      )
      out$term_pair_key[paired_rows] <- pair_keys[pair_row]
      out$paired_term_id[paired_rows] <- ifelse(
        is_gene_row,
        paste0("PEAK:", pair_keys[pair_row]),
        paste0("GENE:", pair_keys[pair_row])
      )
      out$shared_gammafit_candidate[paired_rows] <- shared_candidate[cbind(pair_topic, pair_row)]
      out$common_candidate_count[paired_rows] <- common_count[pair_row]
      out$term_maxprob_selected[paired_rows] <- !is.na(row_max_topic) &
        pair_topic == row_max_topic
      out$term_maxprob_topic[paired_rows] <- row_max_topic
      out$term_maxprob_probability[paired_rows] <- row_max_probability
      out$term_maxprob_margin[paired_rows] <- row_max_margin
      out$assignment_status[paired_rows] <- pair_status[pair_row]
      out$in_topic[paired_rows] <- has_assignment[pair_row] &
        pair_topic == selected_topic[pair_row]
    }
  }
  out[, in_topic := .as_logical_flag(in_topic)]
  data.table::setcolorder(out, c(
    "topic", "topic_num", "term_id", "term_group", "score", "phi",
    "in_topic", "assignment_method", "gammafit_candidate",
    "term_pair_key", "paired_term_id", "shared_gammafit_candidate",
    "common_candidate_count", "term_maxprob_selected",
    "term_maxprob_topic", "term_maxprob_probability",
    "term_maxprob_margin", "assignment_status",
    setdiff(names(out), c(
      "topic", "topic_num", "term_id", "term_group", "score", "phi",
      "in_topic", "assignment_method", "gammafit_candidate",
      "term_pair_key", "paired_term_id", "shared_gammafit_candidate",
      "common_candidate_count", "term_maxprob_selected",
      "term_maxprob_topic", "term_maxprob_probability",
      "term_maxprob_margin", "assignment_status"
    ))
  ))
  out[]
}

.topic_gene_peak_assignment_table <- function(topic_terms) {
  tt <- data.table::as.data.table(topic_terms)
  required <- c(
    "topic_num", "term_id", "term_group", "phi", "in_topic",
    "gammafit_candidate", "term_maxprob_topic",
    "term_maxprob_probability", "term_maxprob_margin"
  )
  if (!nrow(tt) || length(setdiff(required, names(tt)))) {
    return(data.table::data.table())
  }
  if (!"term_pair_key" %in% names(tt)) {
    tt[, term_pair_key := data.table::fcase(
      term_group == "GENE", sub("^GENE:", "", term_id),
      term_group == "PEAK", sub("^PEAK:", "", term_id),
      default = NA_character_
    )]
  }
  tt <- tt[term_group %in% c("GENE", "PEAK") & !is.na(term_pair_key) & nzchar(term_pair_key)]
  if (!nrow(tt)) return(data.table::data.table())

  collapse_topics <- function(x) {
    x <- sort(unique(as.integer(x[is.finite(x)])))
    if (!length(x)) "" else paste(x, collapse = ";")
  }
  summarize_group <- function(group_name, prefix) {
    x <- tt[term_group == group_name]
    if (!nrow(x)) {
      out <- data.table::data.table(
        target_gene = character(),
        term_id = character(),
        gammafit_topic_count = integer(),
        gammafit_topics = character(),
        maxprob_topic = character(),
        maxprob_probability = numeric(),
        maxprob_margin = numeric(),
        assigned_topic = character(),
        assigned_phi = numeric()
      )
      data.table::setnames(
        out,
        setdiff(names(out), "target_gene"),
        paste0(prefix, "_", setdiff(names(out), "target_gene"))
      )
      return(out)
    }
    out <- x[, .(
      term_id = as.character(term_id[[1L]]),
      gammafit_topic_count = sum(.as_logical_flag(gammafit_candidate)),
      gammafit_topics = collapse_topics(topic_num[.as_logical_flag(gammafit_candidate)]),
      maxprob_topic = collapse_topics(term_maxprob_topic),
      maxprob_probability = {
        value <- unique(term_maxprob_probability[is.finite(term_maxprob_probability)])
        if (length(value)) as.numeric(value[[1L]]) else NA_real_
      },
      maxprob_margin = {
        value <- unique(term_maxprob_margin[is.finite(term_maxprob_margin)])
        if (length(value)) as.numeric(value[[1L]]) else NA_real_
      },
      assigned_topic = collapse_topics(topic_num[.as_logical_flag(in_topic)]),
      assigned_phi = {
        value <- phi[.as_logical_flag(in_topic)]
        if (length(value)) as.numeric(value[[1L]]) else NA_real_
      }
    ), by = .(target_gene = term_pair_key)]
    data.table::setnames(
      out,
      setdiff(names(out), "target_gene"),
      paste0(prefix, "_", setdiff(names(out), "target_gene"))
    )
    out
  }

  gene <- summarize_group("GENE", "gene")
  peak <- summarize_group("PEAK", "peak")
  out <- merge(gene, peak, by = "target_gene", all = TRUE, sort = TRUE)
  paired <- tt[term_group == "GENE" & !is.na(paired_term_id), .(
    common_candidate_count = suppressWarnings(max(common_candidate_count, na.rm = TRUE)),
    shared_gammafit_topics = collapse_topics(topic_num[.as_logical_flag(shared_gammafit_candidate)]),
    assignment_status = as.character(assignment_status[[1L]])
  ), by = .(target_gene = term_pair_key)]
  paired[!is.finite(common_candidate_count), common_candidate_count := 0L]
  out <- merge(out, paired, by = "target_gene", all.x = TRUE, sort = TRUE)
  out[, assigned_topic := data.table::fcase(
    !is.na(gene_assigned_topic) & nzchar(gene_assigned_topic) &
      gene_assigned_topic == peak_assigned_topic,
    gene_assigned_topic,
    default = NA_character_
  )]
  out[, assignment_status := data.table::fcase(
    is.na(gene_term_id) | is.na(peak_term_id), "missing_matching_gene_or_peak_term",
    !is.na(assigned_topic) & nzchar(assigned_topic), "assigned_gammafit_maxprob_agreement",
    !is.na(assignment_status) & nzchar(assignment_status), assignment_status,
    default = "no_shared_gammafit_topic"
  )]
  out[, assigned := !is.na(assigned_topic) & nzchar(assigned_topic)]
  out[, min_maxprob_probability := pmin(
    gene_maxprob_probability,
    peak_maxprob_probability,
    na.rm = FALSE
  )]
  out[, min_maxprob_margin := pmin(
    gene_maxprob_margin,
    peak_maxprob_margin,
    na.rm = FALSE
  )]
  data.table::setcolorder(out, c(
    "target_gene", "assigned", "assigned_topic", "assignment_status",
    "gene_term_id", "peak_term_id", "gene_gammafit_topic_count",
    "peak_gammafit_topic_count", "common_candidate_count",
    "gene_gammafit_topics", "peak_gammafit_topics",
    "shared_gammafit_topics", "gene_maxprob_topic", "peak_maxprob_topic",
    "gene_maxprob_probability", "peak_maxprob_probability",
    "min_maxprob_probability", "gene_maxprob_margin",
    "peak_maxprob_margin", "min_maxprob_margin",
    "gene_assigned_phi", "peak_assigned_phi",
    "gene_assigned_topic", "peak_assigned_topic"
  ))
  data.table::setorder(out, -assigned, target_gene)
  out[]
}

.topic_gene_assignment_table <- function(topic_terms) {
  tt <- data.table::as.data.table(topic_terms)
  required <- c(
    "topic_num", "term_id", "term_group", "phi", "score", "in_topic",
    "gammafit_candidate", "term_maxprob_probability", "term_maxprob_margin",
    "assignment_status"
  )
  if (!nrow(tt) || length(setdiff(required, names(tt)))) {
    return(data.table::data.table())
  }
  tt <- tt[term_group == "GENE"]
  if (!nrow(tt)) return(data.table::data.table())
  tt[, .(
    target_gene = sub("^GENE:", "", term_id[[1L]]),
    assigned = any(.as_logical_flag(in_topic)),
    assigned_topic = {
      value <- topic_num[.as_logical_flag(in_topic)]
      if (length(value)) as.integer(value[[1L]]) else NA_integer_
    },
    assigned_phi = {
      value <- phi[.as_logical_flag(in_topic)]
      if (length(value)) as.numeric(value[[1L]]) else NA_real_
    },
    assigned_score = {
      value <- score[.as_logical_flag(in_topic)]
      if (length(value)) as.numeric(value[[1L]]) else NA_real_
    },
    gammafit_topic_count = sum(.as_logical_flag(gammafit_candidate)),
    maxprob_probability = {
      value <- unique(term_maxprob_probability[is.finite(term_maxprob_probability)])
      if (length(value)) as.numeric(value[[1L]]) else NA_real_
    },
    maxprob_margin = {
      value <- unique(term_maxprob_margin[is.finite(term_maxprob_margin)])
      if (length(value)) as.numeric(value[[1L]]) else NA_real_
    },
    assignment_status = {
      value <- unique(as.character(assignment_status))
      value <- value[!is.na(value) & nzchar(value)]
      if (length(value)) value[[1L]] else NA_character_
    }
  ), by = term_id][, term_id := NULL][]
}

.topic_gene_assignment_for_optimization <- function(gene_assignment,
                                                    topic_terms) {
  assignments <- data.table::copy(data.table::as.data.table(gene_assignment))
  terms <- data.table::as.data.table(topic_terms)
  if (!nrow(assignments) || !nrow(terms)) return(assignments)
  candidates <- terms[
    term_group == "GENE" & .as_logical_flag(gammafit_candidate),
    .(gammafit_topics = paste(sort(unique(as.integer(topic_num))), collapse = ";")),
    by = .(target_gene = sub("^GENE:", "", term_id))
  ]
  assignments <- merge(assignments, candidates, by = "target_gene", all.x = TRUE, sort = TRUE)
  candidate_topics <- assignments[["gammafit_topics"]]
  candidate_topics[is.na(candidate_topics)] <- ""
  assignments[, `:=`(
    gene_term_id = paste0("GENE:", target_gene),
    peak_term_id = paste0("GENE:", target_gene),
    gene_gammafit_topics = candidate_topics,
    peak_gammafit_topics = candidate_topics
  )]
  assignments[, gammafit_topics := NULL]
  assignments[]
}

.topic_gene_assignment_summary <- function(gene_assignment,
                                           thrP,
                                           assignment_method,
                                           model_family = NA_character_) {
  x <- data.table::as.data.table(gene_assignment)
  assigned <- x[.as_logical_flag(assigned)]
  data.table::data.table(
    assignment_method = as.character(assignment_method),
    model_family = as.character(model_family),
    fp_term_mode = "gene_expression",
    gammafit_thrP = as.numeric(thrP),
    total_genes = nrow(x),
    assigned_genes = nrow(assigned),
    assigned_percent = if (nrow(x)) 100 * nrow(assigned) / nrow(x) else NA_real_,
    no_gammafit_candidate_count = sum(x$gammafit_topic_count < 1L, na.rm = TRUE),
    recovered_after_topic_merge_count = if ("recovered_after_merge" %in% names(x)) {
      sum(.as_logical_flag(x$recovered_after_merge), na.rm = TRUE)
    } else {
      0L
    },
    lost_after_topic_merge_count = if ("raw_assigned_topic" %in% names(x)) {
      sum(
        is.finite(suppressWarnings(as.integer(x$raw_assigned_topic))) &
          !.as_logical_flag(x$assigned),
        na.rm = TRUE
      )
    } else {
      0L
    }
  )
}

.topic_gene_peak_assignment_summary <- function(pair_assignment,
                                                thrP,
                                                assignment_method,
                                                model_family = NA_character_) {
  x <- data.table::as.data.table(pair_assignment)
  matched <- x[!is.na(gene_term_id) & !is.na(peak_term_id)]
  assigned <- matched[.as_logical_flag(assigned)]
  reason_status <- if ("raw_assignment_status" %in% names(matched)) {
    as.character(matched$raw_assignment_status)
  } else {
    as.character(matched$assignment_status)
  }
  finite_quantile <- function(value, probability) {
    value <- as.numeric(value[is.finite(value)])
    if (!length(value)) return(NA_real_)
    as.numeric(stats::quantile(value, probability, names = FALSE, type = 7L))
  }
  data.table::data.table(
    assignment_method = as.character(assignment_method),
    model_family = as.character(model_family),
    gammafit_thrP = as.numeric(thrP),
    gene_term_count = sum(!is.na(x$gene_term_id)),
    peak_term_count = sum(!is.na(x$peak_term_id)),
    matched_target_count = nrow(matched),
    assigned_target_count = nrow(assigned),
    assigned_target_percent = if (nrow(matched)) 100 * nrow(assigned) / nrow(matched) else NA_real_,
    no_gene_gammafit_candidate_count = sum(
      reason_status %in% c("no_gammafit_candidate", "gene_no_gammafit_candidate"),
      na.rm = TRUE
    ),
    no_peak_gammafit_candidate_count = sum(
      reason_status %in% c("no_gammafit_candidate", "peak_no_gammafit_candidate"),
      na.rm = TRUE
    ),
    no_shared_topic_count = sum(reason_status == "no_shared_gammafit_topic", na.rm = TRUE),
    maxprob_topic_disagreement_count = sum(
      reason_status == "maxprob_topic_disagreement",
      na.rm = TRUE
    ),
    recovered_after_topic_merge_count = if ("recovered_after_merge" %in% names(matched)) {
      sum(.as_logical_flag(matched$recovered_after_merge), na.rm = TRUE)
    } else {
      0L
    },
    lost_after_topic_merge_count = if ("raw_assigned_topic" %in% names(matched)) {
      sum(
        is.finite(suppressWarnings(as.integer(matched$raw_assigned_topic))) &
          !.as_logical_flag(matched$assigned),
        na.rm = TRUE
      )
    } else {
      0L
    },
    multiple_common_topic_count = sum(matched$common_candidate_count > 1L, na.rm = TRUE),
    multiple_common_topic_percent = if (nrow(assigned)) {
      100 * sum(assigned$common_candidate_count > 1L, na.rm = TRUE) / nrow(assigned)
    } else {
      NA_real_
    },
    min_maxprob_probability_p10 = finite_quantile(assigned$min_maxprob_probability, 0.1),
    min_maxprob_probability_median = finite_quantile(assigned$min_maxprob_probability, 0.5),
    min_maxprob_margin_p10 = finite_quantile(assigned$min_maxprob_margin, 0.1),
    min_maxprob_margin_median = finite_quantile(assigned$min_maxprob_margin, 0.5)
  )
}

# Gamma-fit diagnostic plots (cisTopic-like)
plot_gammafit_binarize <- function(score_mat,
                                   out_file,
                                   thrP = 0.975,
                                   min_terms = 50L,
                                   gammafit_scope = c("topic_term_group", "global_term_group"),
                                   title_prefix = NULL,
                                   tf_list = NULL,
                                   edges_docs = NULL,
                                   topic_terms = NULL,
                                   topic_links = NULL,
                                   item_coverage = NULL,
                                   show_peak_expanded_link_coverage = TRUE,
                                   show_gammafit_pages = TRUE,
                                   panels_per_row = 5L,
                                   panels_per_col = 2L,
                                   breaks = 100L) {
  score_mat <- as.matrix(score_mat)
  gammafit_scope <- match.arg(gammafit_scope)
  K <- nrow(score_mat)
  if (!K) return(invisible(NULL))

  min_terms <- as.integer(min_terms)
  if (!is.finite(min_terms) || min_terms < 1L) min_terms <- 1L
  terms <- colnames(score_mat)
  if (is.null(terms)) {
    terms <- paste0("term_", seq_len(ncol(score_mat)))
    colnames(score_mat) <- terms
  }
  show_gammafit_pages <- isTRUE(show_gammafit_pages)
  # Keep gammafit small multiples visually consistent across runs.
  gamma_cols <- 5L
  gamma_rows <- 5L
  per_page <- gamma_cols * gamma_rows

  term_grp <- .term_group(terms)
  idx_peak <- which(term_grp == "PEAK")
  idx_gene <- which(term_grp == "GENE")
  cut_tbl <- .gammafit_cutoffs_by_termclass(
    score_mat,
    thrP = thrP,
    min_terms = min_terms,
    gammafit_scope = gammafit_scope
  )
  peak_cutoffs <- cut_tbl$peaks_gamma_cutoff
  gene_cutoffs <- cut_tbl$gene_gamma_cutoff
  n_selected_peak <- integer(K)
  n_selected_gene <- integer(K)
  for (k in seq_len(K)) {
    sc <- as.numeric(score_mat[k, ])
    sc[!is.finite(sc)] <- 0
    if (length(idx_peak)) {
      keep_peak <- is.finite(peak_cutoffs[k]) & sc[idx_peak] >= peak_cutoffs[k] & sc[idx_peak] > 0
      if (sum(keep_peak, na.rm = TRUE) < as.integer(min_terms)) {
        ord_peak <- idx_peak[order(sc[idx_peak], decreasing = TRUE)]
        take <- min(length(ord_peak), as.integer(min_terms))
        keep_peak <- idx_peak %in% ord_peak[seq_len(take)]
      }
      n_selected_peak[k] <- sum(keep_peak, na.rm = TRUE)
    }
    if (length(idx_gene)) {
      keep_gene <- is.finite(gene_cutoffs[k]) & sc[idx_gene] >= gene_cutoffs[k] & sc[idx_gene] > 0
      if (sum(keep_gene, na.rm = TRUE) < as.integer(min_terms)) {
        ord_gene <- idx_gene[order(sc[idx_gene], decreasing = TRUE)]
        take <- min(length(ord_gene), as.integer(min_terms))
        keep_gene <- idx_gene %in% ord_gene[seq_len(take)]
      }
      n_selected_gene[k] <- sum(keep_gene, na.rm = TRUE)
    }
  }
  coverage_tbl <- .topic_assignment_coverage_summary_table(
    topic_terms,
    score_mat,
    item_coverage,
    show_peak_expanded_link_coverage = show_peak_expanded_link_coverage
  )
  .safe_hist_breaks <- function(x, n = 100L) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(0, 1))
    xr <- range(x, finite = TRUE)
    if (!all(is.finite(xr))) return(c(0, 1))
    if (xr[1] == xr[2]) {
      eps <- if (xr[1] == 0) 1e-6 else abs(xr[1]) * 1e-6
      return(c(xr[1] - eps, xr[2] + eps))
    }
    br <- pretty(xr, n = max(10L, as.integer(n)))
    br <- unique(as.numeric(br))
    br <- br[is.finite(br)]
    if (length(br) < 2L) return(c(xr[1], xr[2]))
    br
  }

  grDevices::pdf(out_file, width = 15, height = 15)
  device_id <- grDevices::dev.cur()
  par_opts <- graphics::par(no.readonly = TRUE)
  device_open <- TRUE
  on.exit({
    open_devices <- grDevices::dev.list()
    if (isTRUE(device_open) && !is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.set(device_id)
      suppressWarnings(graphics::par(par_opts))
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)

  if (nrow(coverage_tbl)) {
    graphics::par(mfrow = c(1, 1), mar = c(5, 10, 4, 2), oma = c(0, 0, 0, 0))
    coverage_plot <- data.table::copy(rev(coverage_tbl))
    coverage_plot[, plot_label := data.table::fcase(
      label == "GENE terms", "Gene terms",
      label == "PEAK terms", "Peak terms",
      default = label
    )]
    pct <- pmax(0, pmin(100, coverage_plot$percent))
    bar_cols <- c(
      Terms = "#4C78A8",
      "GENE terms" = "#54A24B",
      "PEAK terms" = "#72B7B2",
      Genes = "#59A14F",
      Links = "#E15759",
      TFs = "#F58518"
    )
    fill_cols <- unname(bar_cols[coverage_plot$label])
    fill_cols[is.na(fill_cols)] <- "#4C78A8"
    bp <- graphics::barplot(
      pct,
      names.arg = coverage_plot$plot_label,
      horiz = TRUE,
      las = 1,
      xlim = c(0, 110),
      col = fill_cols,
      border = NA,
      xlab = "Assigned to at least one topic (%)",
      main = paste(c(title_prefix, "Topic Assignment Coverage"), collapse = " - "),
      cex.main = 1.15,
      cex.names = 0.95,
      font.main = 2,
      font.lab = 2,
      font.axis = 2
    )
    graphics::text(
      x = pmin(104, pct + 2),
      y = bp,
      labels = coverage_plot$label_text,
      adj = c(0, 0.5),
      cex = 0.9,
      font = 2
    )
  }

  if (show_gammafit_pages) {
    total_panels <- K * 2L
    n_pages <- max(1L, ceiling(total_panels / per_page))
    for (page_idx in seq_len(n_pages)) {
      start <- (page_idx - 1L) * per_page + 1L
      graphics::par(
        mfrow = c(gamma_rows, gamma_cols),
        mar = c(4.2, 4.6, 2.4, 1.2),
        mgp = c(2.4, 0.8, 0),
        oma = c(0, 0, 2, 0)
      )
      for (slot_idx in seq_len(per_page)) {
        panel_idx <- start + slot_idx - 1L
        if (panel_idx > total_panels) {
          graphics::plot.new()
          next
        }
        k <- ((panel_idx - 1L) %/% 2L) + 1L
        panel_type <- if ((panel_idx %% 2L) == 1L) "GENE" else "PEAK"
        sc <- as.numeric(score_mat[k, ])
        sc[!is.finite(sc)] <- 0
        peak_sc <- if (length(idx_peak)) sc[idx_peak] else numeric(0)
        gene_sc <- if (length(idx_gene)) sc[idx_gene] else numeric(0)
        sc_panel <- if (identical(panel_type, "GENE")) gene_sc else peak_sc
        if (!length(sc_panel) || !any(sc_panel > 0)) {
          graphics::plot.new()
          graphics::title(main = paste0("Topic", k, " ", panel_type), font.main = 2, cex.main = 0.9)
          graphics::text(0.5, 0.5, paste0("No positive ", tolower(panel_type), " scores"), cex = 0.8)
          next
        }
        col_fill <- if (identical(panel_type, "GENE")) "#54A24B" else "#4C78A8"
        col_line <- if (identical(panel_type, "GENE")) "#1B7F2A" else "#1F4E8C"
        cutoff_k <- if (identical(panel_type, "GENE")) gene_cutoffs[k] else peak_cutoffs[k]
        n_sel_k <- if (identical(panel_type, "GENE")) n_selected_gene[k] else n_selected_peak[k]
        graphics::hist(
          sc_panel,
          breaks = .safe_hist_breaks(sc_panel, n = as.integer(breaks)),
          prob = TRUE,
          main = paste0("Topic", k, " ", panel_type),
          col = grDevices::adjustcolor(col_fill, alpha.f = 0.55),
          xlab = "Score",
          ylab = "Density",
          cex.main = 0.9,
          font.main = 2
        )
        fit_panel <- .gamma_moments(sc_panel)
        if (!is.null(fit_panel)) {
          graphics::curve(
            stats::dgamma(x, rate = fit_panel$rate, shape = fit_panel$shape),
            add = TRUE,
            col = col_line,
            lwd = 2
          )
        }
        graphics::abline(v = cutoff_k, lty = 2, lwd = 2, col = col_line)
        ann_txt <- paste0(
          panel_type, " > ", signif(cutoff_k, 3), " (n=", n_sel_k, ")"
        )
        ann_txt <- paste(strwrap(ann_txt, width = 28), collapse = "\n")
        usr <- graphics::par("usr")
        graphics::text(
          x = usr[2],
          y = usr[4],
          labels = ann_txt,
          adj = c(1, 1),
          cex = 0.7,
          font = 2
        )
      }
      if (!is.null(title_prefix) && nzchar(title_prefix)) {
        graphics::mtext(
          paste0(title_prefix, " gammafit thresholds"),
          outer = TRUE,
          cex = 1.1,
          font = 2,
          line = 0.5
        )
      }
    }
  }

  grDevices::dev.set(device_id)
  suppressWarnings(graphics::par(par_opts))
  grDevices::dev.off(device_id)
  device_open <- FALSE

  invisible(list(peaks_gamma_cutoff = peak_cutoffs, gene_gamma_cutoff = gene_cutoffs))
}

plot_link_efdr_summary <- function(topic_links,
                                   out_file,
                                   title_prefix = NULL,
                                   fdr_q = 0.2,
                                   panels_per_page = 25L,
                                   ncol = 5L) {
  dt <- data.table::as.data.table(topic_links)
  req <- c("topic_num", "link_score", "link_efdr_q", "link_pass")
  if (!nrow(dt) || !all(req %in% names(dt))) return(invisible(NULL))
  dt[, topic_num := suppressWarnings(as.integer(topic_num))]
  dt[, link_score := .safe_num(link_score)]
  dt[, link_efdr_q := .safe_num(link_efdr_q)]
  dt[, link_pass := .as_logical_flag(link_pass)]
  dt <- dt[is.finite(topic_num) & topic_num >= 1L]
  if (!nrow(dt)) return(invisible(NULL))

  topics <- sort(unique(dt$topic_num))
  panels_per_page <- as.integer(panels_per_page)
  if (!is.finite(panels_per_page) || panels_per_page < 1L) panels_per_page <- 25L
  ncol <- as.integer(ncol)
  if (!is.finite(ncol) || ncol < 1L) ncol <- 5L
  nrow <- as.integer(ceiling(panels_per_page / ncol))
  if (!is.finite(nrow) || nrow < 1L) nrow <- 5L
  .safe_hist_breaks <- function(x, n = 80L) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(0, 1))
    xr <- range(x, finite = TRUE)
    if (!all(is.finite(xr))) return(c(0, 1))
    if (xr[1] == xr[2]) {
      eps <- if (xr[1] == 0) 1e-6 else abs(xr[1]) * 1e-6
      return(c(xr[1] - eps, xr[2] + eps))
    }
    br <- pretty(xr, n = max(10L, as.integer(n)))
    br <- unique(as.numeric(br))
    br <- br[is.finite(br)]
    if (length(br) < 2L) return(c(xr[1], xr[2]))
    br
  }

  grDevices::pdf(out_file, width = 15, height = 15)
  op <- graphics::par(no.readonly = TRUE)
  on.exit({
    suppressWarnings(graphics::par(op))
    grDevices::dev.off()
  }, add = TRUE)

  # Page A: link_score distributions by topic, highlighting eFDR pass links.
  n_pages <- max(1L, ceiling(length(topics) / panels_per_page))
  for (page_idx in seq_len(n_pages)) {
    start <- (page_idx - 1L) * panels_per_page + 1L
    graphics::par(mfrow = c(nrow, ncol), mar = c(4.2, 4.6, 2.4, 1.2), mgp = c(2.4, 0.8, 0), oma = c(0, 0, 2, 0))
    for (slot_idx in seq_len(panels_per_page)) {
      k_idx <- start + slot_idx - 1L
      if (k_idx > length(topics)) {
        graphics::plot.new()
        next
      }
      k <- topics[k_idx]
      dtk <- dt[topic_num == k]
      sc_all <- dtk$link_score[is.finite(dtk$link_score) & dtk$link_score >= 0]
      sc_pass <- dtk$link_score[dtk$link_pass & is.finite(dtk$link_score) & dtk$link_score >= 0]
      if (!length(sc_all)) {
        graphics::plot.new()
        graphics::title(main = paste0("Topic", k), font.main = 2, cex.main = 0.9)
        graphics::text(0.5, 0.5, "No link_score", cex = 0.8)
        next
      }
      graphics::hist(
        sc_all,
        breaks = .safe_hist_breaks(sc_all, n = 80L),
        main = paste0("Topic", k),
        xlab = "link_score",
        ylab = "Count",
        col = grDevices::adjustcolor("#8da0cb", alpha.f = 0.55),
        border = NA
      )
      if (length(sc_pass)) {
        graphics::hist(
          sc_pass,
          breaks = .safe_hist_breaks(sc_pass, n = 80L),
          add = TRUE,
          col = grDevices::adjustcolor("#fc8d62", alpha.f = 0.55),
          border = NA
        )
        sc_cut <- min(sc_pass, na.rm = TRUE)
        if (is.finite(sc_cut)) graphics::abline(v = sc_cut, lty = 2, lwd = 2, col = "#d7301f")
      }
      graphics::legend(
        "topright",
        legend = c("All links", "eFDR pass", "min pass score"),
        fill = c(grDevices::adjustcolor("#8da0cb", alpha.f = 0.55), grDevices::adjustcolor("#fc8d62", alpha.f = 0.55), NA),
        border = c(NA, NA, NA),
        lty = c(NA, NA, 2),
        col = c(NA, NA, "#d7301f"),
        bty = "n",
        cex = 0.65
      )
      ann <- sprintf("pass=%d/%d", length(sc_pass), length(sc_all))
      usr <- graphics::par("usr")
      graphics::text(usr[2], usr[4], labels = ann, adj = c(1, 1), cex = 0.7, font = 2)
    }
    main_txt <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
      paste0(title_prefix, " link_score eFDR diagnostics")
    } else {
      "Link_score eFDR diagnostics"
    }
    graphics::mtext(main_txt, outer = TRUE, cex = 1.1, font = 2, line = 0.5)
  }

  # Page B: link_efdr_q distributions by topic with threshold line.
  n_pages_q <- max(1L, ceiling(length(topics) / panels_per_page))
  for (page_idx in seq_len(n_pages_q)) {
    start <- (page_idx - 1L) * panels_per_page + 1L
    graphics::par(mfrow = c(nrow, ncol), mar = c(4.2, 4.6, 2.4, 1.2), mgp = c(2.4, 0.8, 0), oma = c(0, 0, 2, 0))
    for (slot_idx in seq_len(panels_per_page)) {
      k_idx <- start + slot_idx - 1L
      if (k_idx > length(topics)) {
        graphics::plot.new()
        next
      }
      k <- topics[k_idx]
      dtk <- dt[topic_num == k]
      q_all <- dtk$link_efdr_q[is.finite(dtk$link_efdr_q) & dtk$link_efdr_q >= 0]
      if (!length(q_all)) {
        graphics::plot.new()
        graphics::title(main = paste0("Topic", k), font.main = 2, cex.main = 0.9)
        graphics::text(0.5, 0.5, "No link_efdr_q", cex = 0.8)
        next
      }
      graphics::hist(
        q_all,
        breaks = .safe_hist_breaks(q_all, n = 80L),
        main = paste0("Topic", k),
        xlab = "link_efdr_q",
        ylab = "Count",
        col = grDevices::adjustcolor("#66c2a5", alpha.f = 0.6),
        border = NA
      )
      if (is.finite(fdr_q)) graphics::abline(v = fdr_q, lty = 2, lwd = 2, col = "#1f78b4")
      ann <- sprintf("q<=%.3g: %d/%d", fdr_q, sum(dtk$link_pass, na.rm = TRUE), nrow(dtk))
      usr <- graphics::par("usr")
      graphics::text(usr[2], usr[4], labels = ann, adj = c(1, 1), cex = 0.7, font = 2)
    }
    main_txt <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
      paste0(title_prefix, " link_efdr_q diagnostics")
    } else {
      "Link_efdr_q diagnostics"
    }
    graphics::mtext(main_txt, outer = TRUE, cex = 1.1, font = 2, line = 0.5)
  }

  # Page C: per-topic pass/fail stacked bars.
  cnt <- dt[, .(
    pass = sum(link_pass, na.rm = TRUE),
    total = .N
  ), by = topic_num]
  cnt[, fail := pmax(total - pass, 0L)]
  cnt <- cnt[order(topic_num)]
  mat <- rbind(Pass = cnt$pass, Fail = cnt$fail)
  colnames(mat) <- paste0("Topic", cnt$topic_num)
  graphics::par(mfrow = c(1, 1), mar = c(5, 8, 4, 1), oma = c(0, 0, 1.5, 0))
  graphics::barplot(
    mat,
    horiz = TRUE,
    beside = FALSE,
    col = c(Pass = "#1b9e77", Fail = "#d95f02"),
    border = NA,
    xlab = "Links",
    ylab = "",
    main = "Link pass/fail by topic",
    las = 1
  )
  graphics::legend("topright", legend = c("Pass", "Fail"), fill = c("#1b9e77", "#d95f02"), border = NA, bty = "n")
  if (!is.null(title_prefix) && nzchar(title_prefix)) {
    graphics::mtext(title_prefix, outer = TRUE, cex = 1.05, font = 2, line = 0.4)
  }

  invisible(TRUE)
}

# =============================================================================
# 8) Saving + reloadable artifacts (cisTopic tmp-like behavior)
# =============================================================================

.save_all <- function(out_dir, name, obj) {
  rds_dir <- file.path(out_dir, "rds")
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(rds_dir, sprintf("%s.rds", name)))
  invisible(TRUE)
}

# =============================================================================
# 9) Plot helpers (heatmaps + intertopic distance)
# =============================================================================

.safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

.topic_model_selection_title <- function(out_dir, backend_label = NULL, vae_variant = NULL) {
  run_base <- basename(out_dir)
  cell <- sub("_vae_.*$", "", run_base)
  if (!nzchar(cell) || identical(cell, run_base)) {
    cell <- NA_character_
  }
  combo <- basename(dirname(out_dir))
  combo <- sub("^doc_", "doc_", combo)
  combo <- sub("_model_.*$", "", combo)
  backend <- backend_label
  if (is.null(backend) || !nzchar(backend)) {
    backend <- if (grepl("warplda|wlda", run_base, ignore.case = TRUE) ||
      grepl("model_wlda", basename(dirname(out_dir)), ignore.case = TRUE)) {
      "WarpLDA"
    } else if (!is.null(vae_variant) && nzchar(vae_variant)) {
      paste("VAE", vae_variant)
    } else {
      "Topic model"
    }
  }
  pieces <- c(backend, cell, combo)
  pieces <- pieces[!is.na(pieces) & nzchar(pieces)]
  paste(pieces, collapse = " | ")
}

.as_snake_token <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)
  x
}

.as_safe_id_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)
  x
}

.fmt_topic_token_num <- function(x) {
  x_chr <- format(as.numeric(x), trim = TRUE, scientific = FALSE)
  gsub("\\.", "p", x_chr)
}

.map_topic_count_token <- function(x) {
  switch(
    as.character(x),
    pseudo_count_log = "count_pcl",
    pseudo_count_bin = "count_pcb",
    weight = "count_w",
    as.character(x)
  )
}

.resolve_topic_count_input <- function(count_method = c("bin", "log"), count_input = NULL) {
  count_method <- match.arg(count_method)
  choices <- c("pseudo_count_bin", "pseudo_count_log", "weight")
  if (is.null(count_input) || !length(count_input) || is.na(count_input[[1L]]) || !nzchar(as.character(count_input[[1L]]))) {
    return(if (identical(count_method, "log")) "pseudo_count_log" else "pseudo_count_bin")
  }
  match.arg(as.character(count_input[[1L]]), choices)
}

.map_topic_gene_token <- function(x) {
  switch(
    as.character(x),
    aggregate = "gene_agg",
    unique = "gene_uniq",
    as.character(x)
  )
}

.map_topic_model_token <- function(backend, vae_variant) {
  if (identical(backend, "vae")) {
    if (identical(vae_variant, "multivi_encoder")) "model_mve" else paste0("model_", .as_snake_token(vae_variant))
  } else {
    "model_wlda"
  }
}

.short_link_method_tag <- function(method) {
  switch(
    as.character(method),
    peak_and_gene = "pg",
    peak_and_gene_prob = "pgm",
    gene_only = "gene",
    link_score_prob = "lsp",
    link_score_efdr = "lsefdr",
    as.character(method)
  )
}

.pathway_method_suffix <- function(method) {
  switch(
    as.character(method),
    peak_and_gene = "peak_pass_gene_pass",
    peak_and_gene_prob = "peak_pass_gene_prob_pass",
    gene_only = "gene_pass",
    link_score_prob = "link_score_prob_pass",
    link_score_efdr = "link_score_efdr_pass",
    as.character(method)
  )
}

build_topic_compact_run_dirname <- function(thrP_use,
                                            cut_mode,
                                            gene_term_mode,
                                            include_tf_terms,
                                            count_input,
                                            dataset_tag,
                                            doc_mode,
                                            backend,
                                            vae_variant,
                                            k_use) {
  lnk_token <- if (identical(cut_mode, "max")) "link_pg_vs_pgm_max" else paste0("link_pg_vs_pgm_", .fmt_topic_token_num(cut_mode))
  d_token <- if (identical(doc_mode, "tf")) "doc_tf" else "doc_ctf"
  tf_token <- if (isTRUE(include_tf_terms)) "tf_on" else "tf_off"
  dataset_token <- .as_safe_id_token(dataset_tag)
  paste0(
    dataset_token,
    "_", d_token,
    "_w_peak_delta_gene_fc",
    "_", .map_topic_gene_token(gene_term_mode),
    "_", tf_token,
    "_", .map_topic_count_token(count_input),
    "_", .map_topic_model_token(backend, vae_variant),
    "_k", as.integer(k_use),
    "_ft_tp", .fmt_topic_token_num(thrP_use),
    "_", lnk_token
  )
}

write_topic_directory_name_readme <- function(out_dir, fields) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readme_path <- file.path(out_dir, "DIRECTORY_NAME_README.txt")
  thrP_val <- if ("thrP" %in% names(fields)) fields$thrP else NA
  cut_val <- if ("link_cutoff_mode" %in% names(fields)) fields$link_cutoff_mode else NA
  lnk_val <- if (is.na(cut_val)) NA_character_ else {
    if (identical(as.character(cut_val), "max")) "link_pg_vs_pgm_max" else paste0("link_pg_vs_pgm_", .fmt_topic_token_num(cut_val))
  }
  g_val <- if ("gene_term_mode" %in% names(fields)) .map_topic_gene_token(fields$gene_term_mode) else NA_character_
  tf_val <- if ("include_tf_terms" %in% names(fields)) if (isTRUE(fields$include_tf_terms)) "tf_on" else "tf_off" else NA_character_
  c_val <- if ("count_input" %in% names(fields)) .map_topic_count_token(fields$count_input) else NA_character_
  id_val <- if ("dataset_tag" %in% names(fields)) .as_safe_id_token(fields$dataset_tag) else NA_character_
  d_val <- if ("doc_mode" %in% names(fields)) if (identical(fields$doc_mode, "tf")) "doc_tf" else "doc_ctf" else NA_character_
  m_val <- if (all(c("backend", "vae_variant") %in% names(fields))) .map_topic_model_token(fields$backend, fields$vae_variant) else NA_character_
  k_val <- if ("selected_k" %in% names(fields)) as.integer(fields$selected_k) else NA_integer_
  lines <- c(
    "Directory Naming Guide",
    "======================",
    sprintf("directory_basename: %s", basename(out_dir)),
    sprintf("directory_path: %s", out_dir),
    "",
    "Compact naming template",
    "{dataset}_{doc_ctf|doc_tf}_w_peak_delta_gene_fc_{gene_agg|gene_uniq}_{tf_off|tf_on}_{count_pcl|count_pcb|count_w}_{model_mve|model_wlda}_k{K}_ft_tp{thrP}_link_pg_vs_pgm_{0p3|0p5|max}",
    "",
    "Token meanings and possible options",
    "- {dataset}: user-defined ID (placed first)",
    "- document construction: doc_* + w_peak_delta_gene_fc + gene_* + tf_* + count_*",
    "- model/training: model_* + k*",
    "- extraction settings: ft_tp* + link_pg_vs_pgm_*",
    "- gene_*: gene term mode (gene_agg, gene_uniq)",
    "- tf_*: include non-target TF terms in docs (tf_off=no, tf_on=yes)",
    "- count_*: term count input (count_pcl=pseudo_count_log, count_pcb=pseudo_count_bin, count_w=weight)",
    "- doc_*: document mode (doc_ctf=tf-cluster docs, doc_tf=tf docs)",
    "- w_peak_delta_gene_fc: term-weight recipe (peak delta footprint + gene fold change)",
    "- model_*: model backend (model_mve=multivi_encoder VAE, model_wlda=WarpLDA)",
    "- k*: selected number of topics (positive integer)",
    "",
    "Values used for this directory",
    sprintf("- dataset: %s", if (is.na(id_val)) "NA" else id_val),
    sprintf("- tp: %s", if (is.na(thrP_val)) "NA" else .fmt_topic_token_num(thrP_val)),
    sprintf("- link: %s", if (is.na(lnk_val)) "NA" else lnk_val),
    sprintf("- gene: %s", if (is.na(g_val)) "NA" else g_val),
    sprintf("- tf: %s", if (is.na(tf_val)) "NA" else tf_val),
    sprintf("- count: %s", if (is.na(c_val)) "NA" else c_val),
    sprintf("- doc: %s", if (is.na(d_val)) "NA" else d_val),
    "- weight: w_peak_delta_gene_fc",
    sprintf("- model: %s", if (is.na(m_val)) "NA" else m_val),
    sprintf("- K: %s", if (is.na(k_val)) "NA" else as.character(k_val))
  )
  writeLines(lines, con = readme_path, useBytes = TRUE)
  invisible(readme_path)
}

summarize_topic_combo_failures <- function(combo_failures, combo_error_log = NULL) {
  if (!length(combo_failures)) {
    .log_inform("All topic-model combinations completed successfully.")
    return(invisible(NULL))
  }
  fail_df <- do.call(rbind, combo_failures)
  if (!is.data.frame(fail_df) || !nrow(fail_df)) {
    .log_warn("Completed with failures, but failed-combo table is empty.")
    return(invisible(NULL))
  }
  if (!is.null(combo_error_log) && nzchar(combo_error_log)) {
    .log_warn("Completed with {nrow(fail_df)} failed combination(s). Error log: {combo_error_log}")
  } else {
    .log_warn("Completed with {nrow(fail_df)} failed combination(s).")
  }
  for (j in seq_len(nrow(fail_df))) {
    .log_warn(sprintf(
      "FAILED row=%d | combo=%s | backend=%s | gene_term_mode=%s | include_tf_terms=%s | count_input=%s | error=%s",
      fail_df$row[j],
      fail_df$combo_tag[j],
      fail_df$backend[j],
      fail_df$gene_term_mode[j],
      fail_df$include_tf_terms[j],
      fail_df$count_input[j],
      fail_df$error[j]
    ))
  }
  invisible(fail_df)
}

make_topic_report_args_simple <- function(thrP,
                                          link_prob_cutoff,
                                          link_fdr_p,
                                          link_topic_method = "gammafit",
                                          topic_score_method = c("normtop_specificity", "rowmax_phi"),
                                          gammafit_scope = c("topic_term_group", "global_term_group"),
                                          modules = list(
                                            pathway = TRUE,
                                            topic_by_comparison = TRUE,
                                            intertopic_distance = TRUE
                                          ),
                                          overwrite = list(
                                            link_topic = TRUE,
                                            pathway = TRUE
                                          ),
                                          extraction_steps = NULL) {
  topic_score_method <- match.arg(topic_score_method)
  gammafit_scope <- match.arg(gammafit_scope)
  extraction_steps <- .normalize_topic_extraction_steps(extraction_steps)
  list(
    extraction_steps = extraction_steps,
    pathway_source = "topic_terms",
    topic_score_method = topic_score_method,
    thrP = thrP,
    gammafit_scope = gammafit_scope,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
    pathway_overwrite = isTRUE(overwrite$pathway),
    pathway_per_comparison = TRUE,
    pathway_per_comparison_dir = ".",
    pathway_per_comparison_flat = TRUE,
    pathway_split_direction = TRUE,
    run_pathway_gsea = FALSE,
    run_link_topic_scores = .topic_step_enabled(extraction_steps, "topic_links", FALSE),
    run_gammafit_summary = .topic_step_enabled(extraction_steps, "gammafit_summary", TRUE),
    run_link_efdr_summary = .topic_step_enabled(extraction_steps, "link_efdr_summary", TRUE),
    link_topic_gate_mode = "none",
    link_topic_overwrite = isTRUE(overwrite$link_topic),
    link_topic_method = link_topic_method,
    link_topic_prob_cutoff = link_prob_cutoff,
    link_topic_fdr_q = 0.5,
    link_topic_fdr_p = link_fdr_p,
    pathway_link_scores_file = "topic_links.csv",
    pathway_link_scores_file_tf = "topic_links.csv",
    pathway_link_gene_terms_file = NULL,
    pathway_link_min_prob = 0,
    pathway_link_include_tf = FALSE,
    pathway_link_include_gene = TRUE,
    pathway_link_gene_min_prob = 0,
    pathway_link_tf_min_prob = 0,
    pathway_link_tf_max_topics = Inf,
    pathway_link_tf_top_n_per_topic = NA_integer_,
    top_n_per_topic = Inf,
    dot_top_n_per_topic = 25L,
    max_pathways = Inf,
    run_pathway_enrichment = .topic_step_enabled(extraction_steps, "pathway", isTRUE(modules$pathway)),
    run_raw_theta_document_heatmap = .topic_step_enabled(extraction_steps, "raw_theta_documents", FALSE),
    run_document_theta_umap = .topic_step_enabled(extraction_steps, "document_theta_umap", TRUE),
    run_topic_term_heatmap = .topic_step_enabled(extraction_steps, "topic_term_heatmap", TRUE),
    run_topic_by_comparison_heatmaps = .topic_step_enabled(extraction_steps, "topic_by_comparison", isTRUE(modules$topic_by_comparison)),
    run_intertopic_distance_map = .topic_step_enabled(extraction_steps, "intertopic_distance", isTRUE(modules$intertopic_distance))
  )
}

.parse_doc_id <- function(doc_id, doc_design = c("comparison", "condition")) {
  doc_design <- match.arg(doc_design)
  doc_id <- as.character(doc_id)
  parts <- data.table::tstrsplit(doc_id, "::", fixed = TRUE)
  mat <- do.call(cbind, lapply(parts, as.character))
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol = 1)
  }

  n_parts <- rowSums(!is.na(mat) & nzchar(mat))
  n <- length(doc_id)
  comparison_id <- if (ncol(mat) >= 1L) mat[, 1] else rep(NA_character_, n)
  direction <- rep(NA_character_, n)
  tf_doc <- rep(NA_character_, n)

  if (identical(doc_design, "condition")) {
    has_tf <- n_parts >= 2L
    if (any(has_tf)) {
      tf_doc[has_tf] <- vapply(
        which(has_tf),
        function(i) paste(mat[i, 2:n_parts[i]], collapse = "::"),
        character(1)
      )
    }
    return(data.table::data.table(comparison_id = comparison_id, tf_doc = tf_doc, direction = direction))
  }

  has_dir <- n_parts >= 2L
  if (any(has_dir)) {
    direction[has_dir] <- mat[cbind(which(has_dir), n_parts[has_dir])]
  }
  has_tf <- n_parts >= 3L
  if (any(has_tf)) {
    tf_doc[has_tf] <- vapply(
      which(has_tf),
      function(i) {
        if (n_parts[i] <= 2L) return(NA_character_)
        paste(mat[i, 2:(n_parts[i] - 1L)], collapse = "::")
      },
      character(1)
    )
  }

  data.table::data.table(comparison_id = comparison_id, tf_doc = tf_doc, direction = direction)
}

.module3_tf_display_label <- function(tf) {
  tf <- as.character(tf)
  data.table::fifelse(tf == "Tbx21", "Tbet", tf)
}

.empty_tf_direction_topic_summary <- function() {
  data.table::data.table(
    comparison_id = character(),
    tf = character(),
    tf_display = character(),
    up_primary_topic = character(),
    up_primary_theta = numeric(),
    up_margin = numeric(),
    up_ambiguous = logical(),
    down_primary_topic = character(),
    down_primary_theta = numeric(),
    down_margin = numeric(),
    down_ambiguous = logical(),
    direction_topic_status = character()
  )
}

.build_tf_topic_assignment_tables <- function(theta,
                                              doc_design = c("comparison", "condition"),
                                              membership_cutoff = 0.3,
                                              primary_margin_cutoff = 0.1) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  membership_cutoff <- suppressWarnings(as.numeric(membership_cutoff))
  primary_margin_cutoff <- suppressWarnings(as.numeric(primary_margin_cutoff))
  if (!is.finite(membership_cutoff) || membership_cutoff < 0 || membership_cutoff > 1) {
    .log_abort("membership_cutoff must be a number between 0 and 1.")
  }
  if (!is.finite(primary_margin_cutoff) || primary_margin_cutoff < 0 || primary_margin_cutoff > 1) {
    .log_abort("primary_margin_cutoff must be a number between 0 and 1.")
  }

  theta <- as.matrix(theta)
  if (!nrow(theta) || !ncol(theta)) {
    empty <- data.table::data.table()
    return(list(membership = empty, pass = empty, primary = empty, direction_summary = empty))
  }
  doc_id <- rownames(theta)
  if (is.null(doc_id) || anyNA(doc_id) || any(!nzchar(doc_id))) {
    .log_abort("theta must have non-empty rownames containing document IDs.")
  }
  topic <- colnames(theta)
  if (is.null(topic) || anyNA(topic) || any(!nzchar(topic))) {
    topic <- paste0("Topic", seq_len(ncol(theta)))
    colnames(theta) <- topic
  }
  topic_num <- suppressWarnings(as.integer(sub("^Topic", "", topic)))
  bad_topic_num <- !is.finite(topic_num)
  if (any(bad_topic_num)) topic_num[bad_topic_num] <- seq_along(topic)[bad_topic_num]

  doc_info <- .parse_doc_id(doc_id, doc_design = doc_design)
  doc_tbl <- data.table::data.table(doc_id = doc_id)
  doc_tbl <- cbind(doc_tbl, doc_info)
  doc_tbl[, tf := as.character(tf_doc)]
  doc_tbl[, tf_display := .module3_tf_display_label(tf)]
  doc_tbl <- doc_tbl[!is.na(tf) & nzchar(tf)]
  if (!nrow(doc_tbl)) {
    empty <- data.table::data.table()
    return(list(membership = empty, pass = empty, primary = empty, direction_summary = empty))
  }

  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, doc_id := doc_id]
  dt <- merge(doc_tbl, theta_dt, by = "doc_id", all.x = TRUE, sort = FALSE)
  topic_cols <- intersect(topic, names(dt))
  membership <- data.table::melt(
    dt,
    id.vars = c("doc_id", "comparison_id", "tf", "tf_display", "direction"),
    measure.vars = topic_cols,
    variable.name = "topic",
    value.name = "theta"
  )
  membership[, theta := as.numeric(theta)]
  membership <- membership[is.finite(theta)]
  topic_lookup <- data.table::data.table(topic = topic, topic_num = topic_num)
  membership <- merge(membership, topic_lookup, by = "topic", all.x = TRUE, sort = FALSE)
  membership[, membership_pass := theta >= membership_cutoff]
  data.table::setcolorder(
    membership,
    c("doc_id", "comparison_id", "tf", "tf_display", "direction", "topic", "topic_num", "theta", "membership_pass")
  )
  data.table::setorder(membership, comparison_id, tf_display, direction, topic_num)
  pass <- membership[membership_pass == TRUE]

  primary <- membership[, {
    ord <- order(-theta, topic_num)
    theta_ord <- theta[ord]
    topic_ord <- topic[ord]
    topic_num_ord <- topic_num[ord]
    total <- sum(theta, na.rm = TRUE)
    p <- if (is.finite(total) && total > 0) theta / total else rep(0, .N)
    entropy <- -sum(p[p > 0] * log(p[p > 0]), na.rm = TRUE)
    primary_theta <- theta_ord[[1L]]
    second_theta <- if (length(theta_ord) >= 2L) theta_ord[[2L]] else NA_real_
    margin <- primary_theta - ifelse(is.finite(second_theta), second_theta, 0)
    .(
      primary_topic = topic_ord[[1L]],
      primary_topic_num = topic_num_ord[[1L]],
      primary_theta = primary_theta,
      second_topic = if (length(topic_ord) >= 2L) topic_ord[[2L]] else NA_character_,
      second_topic_num = if (length(topic_num_ord) >= 2L) topic_num_ord[[2L]] else NA_integer_,
      second_theta = second_theta,
      margin = margin,
      entropy = entropy,
      effective_topics = exp(entropy),
      n_pass_topics = sum(theta >= membership_cutoff, na.rm = TRUE),
      ambiguous = !is.finite(primary_theta) ||
        primary_theta < membership_cutoff ||
        !is.finite(margin) ||
        margin < primary_margin_cutoff
    )
  }, by = .(doc_id, comparison_id, tf, tf_display, direction)]
  data.table::setorder(primary, comparison_id, tf_display, direction)

  direction_summary <- .empty_tf_direction_topic_summary()
  if (identical(doc_design, "comparison")) {
    up <- primary[direction == "Target-Up", .(
      comparison_id, tf, tf_display,
      up_primary_topic = primary_topic,
      up_primary_theta = primary_theta,
      up_margin = margin,
      up_ambiguous = ambiguous
    )]
    down <- primary[direction == "Target-Down", .(
      comparison_id, tf, tf_display,
      down_primary_topic = primary_topic,
      down_primary_theta = primary_theta,
      down_margin = margin,
      down_ambiguous = ambiguous
    )]
    direction_summary <- merge(up, down, by = c("comparison_id", "tf", "tf_display"), all = TRUE, sort = FALSE)
    if (nrow(direction_summary)) {
      direction_summary[, direction_topic_status := data.table::fcase(
        !is.na(up_primary_topic) & !is.na(down_primary_topic) &
          (.as_logical_flag(up_ambiguous) | .as_logical_flag(down_ambiguous)),
        "ambiguous",
        !is.na(up_primary_topic) & !is.na(down_primary_topic) & up_primary_topic == down_primary_topic,
        "same_topic",
        !is.na(up_primary_topic) & !is.na(down_primary_topic) & up_primary_topic != down_primary_topic,
        "different_topic",
        !is.na(up_primary_topic) & is.na(down_primary_topic),
        "up_only",
        is.na(up_primary_topic) & !is.na(down_primary_topic),
        "down_only",
        default = "ambiguous"
      )]
      data.table::setorder(direction_summary, comparison_id, tf_display)
    }
  }

  list(
    membership = membership[],
    pass = pass[],
    primary = primary[],
    direction_summary = direction_summary[]
  )
}

.tf_topic_doc_label <- function(comparison_id, tf_display, direction = NA_character_) {
  out <- paste(comparison_id, tf_display, sep = " | ")
  has_dir <- !is.na(direction) & nzchar(direction)
  out[has_dir] <- paste(out[has_dir], direction[has_dir], sep = " | ")
  out
}

.raw_theta_document_label <- function(comparison_id, tf_display, direction = NA_character_) {
  group <- as.character(comparison_id)
  has_dir <- !is.na(direction) & nzchar(direction)
  group[has_dir] <- paste(group[has_dir], direction[has_dir], sep = " | ")
  paste(tf_display, group, sep = " | ")
}

.topic_factor_palette <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(x)])))
  if (!length(x)) return(character())
  base_cols <- c(
    "#3b5b92", "#b04a58", "#3f7f5f", "#8a6f2a", "#6a4c93",
    "#2f7f9f", "#b66a2c", "#5f6f52", "#9a4f7a", "#4f6f8f",
    "#7c5c3e", "#4b7f7a"
  )
  cols <- grDevices::colorRampPalette(base_cols)(length(x))
  stats::setNames(cols, x)
}

.topic_term_assignment_table <- function(score_mat, topic_terms = NULL) {
  .assert_pkg("data.table")
  score_mat <- as.matrix(score_mat)
  if (!nrow(score_mat) || !ncol(score_mat)) return(data.table::data.table())
  if (is.null(rownames(score_mat))) rownames(score_mat) <- paste0("Topic", seq_len(nrow(score_mat)))
  if (is.null(colnames(score_mat))) colnames(score_mat) <- paste0("term_", seq_len(ncol(score_mat)))
  storage.mode(score_mat) <- "double"
  score_mat[!is.finite(score_mat) | score_mat < 0] <- 0

  terms <- colnames(score_mat)
  topic_names <- rownames(score_mat)
  max_idx <- max.col(t(score_mat), ties.method = "first")
  max_score <- apply(score_mat, 2L, max, na.rm = TRUE)
  max_score[!is.finite(max_score)] <- 0
  out <- data.table::data.table(
    term_id = terms,
    term_group = .term_group(terms),
    max_score = as.numeric(max_score),
    max_score_topic = topic_names[max_idx],
    in_any_topic = FALSE,
    primary_topic = NA_character_,
    primary_phi = NA_real_,
    assignment_method = NA_character_,
    assigned_topics = NA_character_
  )

  if (!is.null(topic_terms) && is.data.frame(topic_terms) && nrow(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (!"term_id" %in% names(tt)) return(out[])
    if ("in_topic" %in% names(tt)) {
      tt <- tt[.as_logical_flag(in_topic)]
    }
    if (!nrow(tt)) return(out[])
    if (!"topic_num" %in% names(tt)) {
      if ("topic" %in% names(tt)) {
        topic_chr <- as.character(tt$topic)
        tt[, topic_num := suppressWarnings(as.integer(gsub("^Topic", "", topic_chr)))]
      } else {
        tt[, topic_num := NA_integer_]
      }
    }
    tt <- tt[term_id %in% terms & is.finite(topic_num) & topic_num >= 1L & topic_num <= length(topic_names)]
    if (!nrow(tt)) return(out[])
    tt[, topic_label := topic_names[as.integer(topic_num)]]
    tt[, term_score := score_mat[cbind(match(topic_label, topic_names), match(term_id, terms))]]
    tt[, term_score := .safe_num(term_score)]
    if (!"phi" %in% names(tt)) tt[, phi := NA_real_]
    if (!"assignment_method" %in% names(tt)) tt[, assignment_method := NA_character_]
    data.table::setorder(tt, term_id, -term_score, topic_num)
    primary <- tt[, .(
      primary_topic = topic_label[[1L]],
      primary_phi = .safe_num(phi[[1L]]),
      assignment_method = as.character(assignment_method[[1L]]),
      assigned_topics = paste(unique(topic_label), collapse = ";")
    ), by = term_id]
    out[primary, `:=`(
      in_any_topic = TRUE,
      primary_topic = i.primary_topic,
      primary_phi = i.primary_phi,
      assignment_method = i.assignment_method,
      assigned_topics = i.assigned_topics
    ), on = "term_id"]
  }
  out[]
}

.write_topic_term_primary_assignment <- function(score_mat,
                                                 topic_terms = NULL,
                                                 assignment_file) {
  .assert_pkg("data.table")
  score_mat <- as.matrix(score_mat)
  if (!nrow(score_mat) || !ncol(score_mat)) return(invisible(NULL))
  if (is.null(rownames(score_mat))) rownames(score_mat) <- paste0("Topic", seq_len(nrow(score_mat)))
  if (is.null(colnames(score_mat))) colnames(score_mat) <- paste0("term_", seq_len(ncol(score_mat)))
  storage.mode(score_mat) <- "double"
  score_mat[!is.finite(score_mat) | score_mat < 0] <- 0
  score_mat[score_mat > 1] <- 1

  assignment <- .topic_term_assignment_table(score_mat, topic_terms = topic_terms)
  if (!nrow(assignment)) return(invisible(NULL))
  assignment[, topic_assignment := data.table::fifelse(in_any_topic, primary_topic, "Not assigned")]
  data.table::setorder(assignment, term_group, topic_assignment, -max_score, term_id)
  dir.create(dirname(assignment_file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(assignment, assignment_file)
  invisible(assignment_file)
}

.plot_topic_term_phi_score_comparison_heatmap <- function(phi,
                                                          score_mat,
                                                          topic_terms = NULL,
                                                          out_file,
                                                          topic_score_method = "normtop_specificity",
                                                          title_prefix = NULL,
                                                          cluster_topics = TRUE,
                                                          show_rownames = FALSE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    return(invisible(NULL))
  }
  .assert_pkg("data.table")
  phi <- as.matrix(phi)
  score_mat <- as.matrix(score_mat)
  if (!nrow(phi) || !ncol(phi) || !nrow(score_mat) || !ncol(score_mat)) return(invisible(NULL))
  if (is.null(rownames(phi))) rownames(phi) <- paste0("Topic", seq_len(nrow(phi)))
  if (is.null(colnames(phi))) colnames(phi) <- paste0("term_", seq_len(ncol(phi)))
  if (is.null(rownames(score_mat))) rownames(score_mat) <- rownames(phi)
  if (is.null(colnames(score_mat))) colnames(score_mat) <- colnames(phi)

  common_topics <- intersect(rownames(score_mat), rownames(phi))
  common_terms <- intersect(colnames(score_mat), colnames(phi))
  if (!length(common_topics) || !length(common_terms)) return(invisible(NULL))
  phi <- phi[common_topics, common_terms, drop = FALSE]
  score_mat <- score_mat[common_topics, common_terms, drop = FALSE]
  storage.mode(phi) <- "double"
  storage.mode(score_mat) <- "double"
  phi[!is.finite(phi) | phi < 0] <- 0
  score_mat[!is.finite(score_mat) | score_mat < 0] <- 0
  phi[phi > 1] <- 1
  score_mat[score_mat > 1] <- 1

  assignment <- .topic_term_assignment_table(score_mat, topic_terms = topic_terms)
  assignment <- assignment[term_id %in% common_terms]
  if (!nrow(assignment)) return(invisible(NULL))
  assignment[, topic_assignment := data.table::fifelse(in_any_topic, primary_topic, "Not assigned")]
  data.table::setorder(assignment, term_group, topic_assignment, -max_score, term_id)
  term_order <- assignment$term_id
  topic_order <- common_topics
  if (isTRUE(cluster_topics) && length(topic_order) > 1L) {
    topic_dist <- try(stats::dist(score_mat[topic_order, term_order, drop = FALSE]), silent = TRUE)
    if (!inherits(topic_dist, "try-error")) {
      topic_hc <- try(stats::hclust(topic_dist), silent = TRUE)
      if (!inherits(topic_hc, "try-error")) {
        topic_order <- topic_hc$labels[topic_hc$order]
      }
    }
  }

  score_display <- t(score_mat[topic_order, term_order, drop = FALSE])
  score_display[!is.finite(score_display)] <- 0

  ann <- data.frame(
    Term_group = factor(assignment$term_group, levels = sort(unique(assignment$term_group))),
    Topic_assignment = factor(assignment$topic_assignment, levels = sort(unique(assignment$topic_assignment))),
    Max_score = assignment$max_score,
    check.names = FALSE
  )
  rownames(ann) <- term_order
  annotation_colors <- list(
    Term_group = .topic_factor_palette(assignment$term_group),
    Topic_assignment = .topic_factor_palette(assignment$topic_assignment)
  )

  k <- ncol(score_display)
  title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, sprintf("%s term score K%d", topic_score_method, k), sep = " | ")
  } else {
    sprintf("%s term score K%d", topic_score_method, k)
  }
  width <- max(8, min(16, 5.8 + 0.46 * k))
  height <- max(6.2, min(60, 3.8 + 0.03 * nrow(score_display)))
  fontsize_col <- if (k > 30L) 5.8 else if (k > 18L) 6.8 else 8
  heatmap_colors <- grDevices::colorRampPalette(c("#f7fbff", "#d6e5f3", "#8fb9d9", "#3478b6", "#08306b"))(101)
  heatmap_breaks <- seq(0, 1, length.out = 102)

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(out_file, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  device_id <- grDevices::dev.cur()
  device_open <- TRUE
  on.exit({
    open_devices <- grDevices::dev.list()
    if (isTRUE(device_open) && !is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)

  ph_score <- pheatmap::pheatmap(
    score_display,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    border_color = NA,
    annotation_row = ann,
    annotation_colors = annotation_colors,
    show_rownames = isTRUE(show_rownames),
    show_colnames = TRUE,
    fontsize = 9,
    fontsize_col = fontsize_col,
    main = title,
    angle_col = 45,
    silent = TRUE
  )

  grDevices::dev.set(device_id)
  grid::grid.newpage()
  grid::grid.draw(ph_score$gtable)
  grDevices::dev.off(device_id)
  device_open <- FALSE
  invisible(out_file)
}

.plot_raw_theta_document_heatmap <- function(theta,
                                             out_file,
                                             doc_design = c("comparison", "condition"),
                                             title_prefix = NULL,
                                             cluster_topics = TRUE,
                                             show_rownames = FALSE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    return(invisible(NULL))
  }
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  theta <- as.matrix(theta)
  if (!nrow(theta) || !ncol(theta)) return(invisible(NULL))
  if (is.null(rownames(theta))) {
    .log_abort("theta must have rownames containing document IDs.")
  }
  if (is.null(colnames(theta))) {
    colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  }
  storage.mode(theta) <- "double"
  theta[!is.finite(theta)] <- 0

  doc_info <- .parse_doc_id(rownames(theta), doc_design = doc_design)
  meta <- data.table::data.table(doc_id = rownames(theta))
  meta <- cbind(meta, doc_info)
  meta[, tf := as.character(tf_doc)]
  meta[, tf_display := .module3_tf_display_label(tf)]
  meta <- meta[!is.na(tf_display) & nzchar(tf_display)]
  if (!nrow(meta)) return(invisible(NULL))
  theta <- theta[meta$doc_id, , drop = FALSE]
  row_max <- apply(theta, 1L, max, na.rm = TRUE)
  row_max[!is.finite(row_max)] <- 0
  meta[, row_max_theta := as.numeric(row_max)]
  meta[, group_label := as.character(comparison_id)]
  if (identical(doc_design, "comparison")) {
    has_dir <- !is.na(meta[["direction"]]) & nzchar(meta[["direction"]])
    meta[has_dir, group_label := paste(comparison_id, direction, sep = " | ")]
  }
  meta[, row_label := .raw_theta_document_label(
    comparison_id = meta[["comparison_id"]],
    tf_display = meta[["tf_display"]],
    direction = meta[["direction"]]
  )]
  meta[, row_label := make.unique(row_label, sep = " #")]
  data.table::setorder(meta, tf_display, -row_max_theta, group_label, doc_id)

  mat <- theta[meta$doc_id, , drop = FALSE]
  rownames(mat) <- meta$row_label
  gaps_row <- integer()
  if (nrow(meta) > 1L) {
    tf_runs <- rle(meta$tf_display)
    gaps_row <- cumsum(tf_runs$lengths)
    gaps_row <- gaps_row[gaps_row < nrow(meta)]
  }
  annotation_row <- data.frame(
    Group = factor(meta$group_label, levels = sort(unique(meta$group_label))),
    row_max_theta = meta$row_max_theta,
    check.names = FALSE
  )
  rownames(annotation_row) <- meta$row_label
  annotation_colors <- list(Group = .topic_factor_palette(meta$group_label))

  k <- ncol(mat)
  title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, sprintf("raw theta documents K%d", k), sep = " | ")
  } else {
    sprintf("Raw theta documents K%d", k)
  }
  width <- max(7.2, min(12.5, 4.8 + 0.38 * k))
  height <- max(6.2, min(60, 3.8 + 0.105 * nrow(mat)))
  fontsize_row <- if (nrow(mat) > 500L) {
    2.1
  } else if (nrow(mat) > 250L) {
    2.8
  } else if (nrow(mat) > 120L) {
    4.0
  } else {
    6.2
  }
  fontsize_col <- if (k > 30L) 5.8 else if (k > 18L) 6.8 else 8
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(out_file, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  device_id <- grDevices::dev.cur()
  device_open <- TRUE
  on.exit({
    open_devices <- grDevices::dev.list()
    if (isTRUE(device_open) && !is.null(open_devices) && device_id %in% open_devices) {
      grDevices::dev.off(device_id)
    }
  }, add = TRUE)
  ph <- pheatmap::pheatmap(
    mat,
    color = grDevices::colorRampPalette(c("#f7fbff", "#d6e5f3", "#8fb9d9", "#3478b6", "#08306b"))(101),
    breaks = seq(0, 1, length.out = 102),
    cluster_rows = FALSE,
    cluster_cols = isTRUE(cluster_topics) && ncol(mat) > 1L,
    scale = "none",
    border_color = NA,
    gaps_row = gaps_row,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    show_rownames = isTRUE(show_rownames),
    show_colnames = TRUE,
    fontsize = 9,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    main = title,
    angle_col = 45,
    silent = TRUE
  )
  grDevices::dev.set(device_id)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off(device_id)
  device_open <- FALSE
  invisible(out_file)
}

.module3_topic_condition_colors <- function(cfg) {
  report <- cfg$report %||% cfg$report_state %||% list()
  if (!is.list(report)) .log_abort("report in the project config must be a mapping.")
  colors <- report$condition_colors %||% cfg$condition_colors %||% list()
  if (is.list(colors)) colors <- unlist(colors, use.names = TRUE)
  if (!length(colors)) return(character())
  color_names <- names(colors)
  colors <- stats::setNames(as.character(unname(colors)), color_names)
  if (is.null(color_names) || any(!nzchar(color_names))) {
    .log_abort("condition_colors must be a named mapping of condition IDs to hex colors.")
  }
  if (any(!grepl("^#[0-9A-Fa-f]{6}$", colors))) {
    .log_abort("condition_colors values must be six-digit hex colors such as #4E79A7.")
  }
  stats::setNames(toupper(colors), color_names)
}

.module3_bright_topic_palette <- function(topics) {
  topics <- unique(as.character(topics[!is.na(topics) & nzchar(topics)]))
  if (!length(topics)) return(character())
  colors <- c(
    "#E15759", "#4E79A7", "#59A14F", "#F28E2B", "#B07AA1",
    "#00A6D6", "#EDC948", "#9C755F", "#FF5DA2", "#00A878",
    "#D95F02", "#7655C5", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#1F78B4", "#33A02C", "#E31A1C", "#6A3D9A",
    "#B15928", "#17BECF", "#BCBD22", "#F05A9D"
  )
  stats::setNames(rep(colors, length.out = length(topics)), topics)
}

.plot_document_theta_umap <- function(theta,
                                      out_file,
                                      doc_design = c("comparison", "condition"),
                                      selected_tfs = NULL,
                                      top_n_tfs = 12L,
                                      seed = 123L,
                                      n_neighbors = 30L,
                                      condition_colors = NULL,
                                      title_prefix = NULL) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    .log_warn("Skipping document theta UMAP because package {.pkg uwot} is not installed.")
    return(invisible(NULL))
  }
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  doc_design <- match.arg(doc_design)
  theta <- as.matrix(theta)
  if (nrow(theta) < 4L || ncol(theta) < 2L) {
    return(invisible(NULL))
  }
  if (is.null(rownames(theta)) || any(!nzchar(rownames(theta)))) {
    .log_abort("theta must have non-empty rownames containing document IDs.")
  }
  storage.mode(theta) <- "double"
  theta[!is.finite(theta) | theta < 0] <- 0
  row_total <- rowSums(theta)
  keep <- is.finite(row_total) & row_total > 0
  theta <- theta[keep, , drop = FALSE]
  row_total <- row_total[keep]
  if (nrow(theta) < 4L) return(invisible(NULL))
  theta <- theta / row_total

  doc_info <- .parse_doc_id(rownames(theta), doc_design = doc_design)
  meta <- data.table::data.table(doc_id = rownames(theta))
  meta <- cbind(meta, doc_info)
  meta[, tf := as.character(tf_doc)]
  meta[, tf_display := .module3_tf_display_label(tf)]
  meta[, group_label := as.character(comparison_id)]
  if (identical(doc_design, "comparison")) {
    has_direction <- !is.na(direction) & nzchar(direction)
    meta[has_direction, group_label := paste(comparison_id, direction, sep = " | ")]
  }

  top_n_tfs <- suppressWarnings(as.integer(top_n_tfs[[1L]]))
  if (!is.finite(top_n_tfs) || top_n_tfs < 1L) top_n_tfs <- 12L
  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, tf_display := meta$tf_display]
  topic_cols <- setdiff(names(theta_dt), "tf_display")
  tf_scores <- theta_dt[!is.na(tf_display) & nzchar(tf_display), {
    values <- as.matrix(.SD)
    score <- if (nrow(values) > 1L) sum(apply(values, 2L, stats::var), na.rm = TRUE) else 0
    .(theta_variation = as.numeric(score), n_documents = .N)
  }, by = tf_display, .SDcols = topic_cols]
  data.table::setorder(tf_scores, -theta_variation, -n_documents, tf_display)

  available_tfs <- sort(unique(meta$tf_display[!is.na(meta$tf_display) & nzchar(meta$tf_display)]))
  if (is.null(selected_tfs) || !length(selected_tfs)) {
    selected_tfs <- head(tf_scores$tf_display, top_n_tfs)
  } else {
    requested <- unique(as.character(selected_tfs))
    matched <- match(toupper(requested), toupper(available_tfs))
    selected_tfs <- available_tfs[stats::na.omit(matched)]
    missing_tfs <- requested[is.na(matched)]
    if (length(missing_tfs)) {
      .log_warn("Ignoring selected TFs absent from theta documents: {paste(missing_tfs, collapse = ', ')}.")
    }
  }
  selected_tfs <- unique(selected_tfs)

  n_neighbors <- suppressWarnings(as.integer(n_neighbors[[1L]]))
  if (!is.finite(n_neighbors)) n_neighbors <- 30L
  n_neighbors <- max(2L, min(n_neighbors, nrow(theta) - 1L))
  seed <- suppressWarnings(as.integer(seed[[1L]]))
  if (!is.finite(seed)) seed <- 123L
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  coords <- uwot::umap(
    sqrt(theta),
    n_neighbors = n_neighbors,
    min_dist = 0.15,
    metric = "euclidean",
    n_components = 2L,
    n_threads = 1L,
    ret_model = FALSE,
    verbose = FALSE
  )
  coords <- data.table::data.table(
    doc_id = rownames(theta),
    UMAP1 = as.numeric(coords[, 1L]),
    UMAP2 = as.numeric(coords[, 2L])
  )
  coords <- cbind(coords, meta[, .(group_label, comparison_id, direction, tf, tf_display)])
  coords[, primary_topic := colnames(theta)[max.col(theta, ties.method = "first")]]
  coords[, primary_theta := apply(theta, 1L, max)]
  coords[, selected_tf := tf_display %in% selected_tfs]
  coords <- merge(coords, tf_scores, by = "tf_display", all.x = TRUE, sort = FALSE)

  group_colors <- .topic_factor_palette(coords$group_label)
  if (!is.null(condition_colors) && length(condition_colors)) {
    if (is.list(condition_colors)) condition_colors <- unlist(condition_colors, use.names = TRUE)
    color_names <- names(condition_colors)
    condition_colors <- stats::setNames(as.character(unname(condition_colors)), color_names)
    if (is.null(color_names) || any(!nzchar(color_names))) {
      .log_abort("condition_colors must be a named mapping of condition IDs to hex colors.")
    }
    if (any(!grepl("^#[0-9A-Fa-f]{6}$", condition_colors))) {
      .log_abort("condition_colors values must be six-digit hex colors such as #4E79A7.")
    }
    matched_groups <- intersect(names(group_colors), names(condition_colors))
    group_colors[matched_groups] <- toupper(condition_colors[matched_groups])
  }
  coords[, condition_color := unname(group_colors[as.character(group_label)])]
  topic_levels <- colnames(theta)
  topic_colors <- .module3_bright_topic_palette(topic_levels)
  coords[, topic_color := unname(topic_colors[as.character(primary_topic)])]

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  csv_file <- sub("[.]pdf$", ".csv", out_file, ignore.case = TRUE)
  selected_file <- sub("[.]pdf$", "_selected_tfs.csv", out_file, ignore.case = TRUE)
  data.table::fwrite(coords, csv_file)
  data.table::fwrite(
    tf_scores[, selected := tf_display %in% selected_tfs][order(-selected, -theta_variation, tf_display)],
    selected_file
  )

  title_base <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, sprintf("document theta UMAP K%d", ncol(theta)), sep = " | ")
  } else {
    sprintf("Document theta UMAP K%d", ncol(theta))
  }
  title_base <- paste(strwrap(title_base, width = 95L), collapse = "\n")
  common_theme <- ggplot2::theme_minimal(base_family = "Helvetica", base_size = 9) +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "black", size = 9),
      axis.title = ggplot2::element_text(face = "bold", size = 9),
      axis.text = ggplot2::element_text(face = "bold", size = 8, color = "black"),
      plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 8, hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    )
  x_range <- range(coords$UMAP1, na.rm = TRUE)
  y_range <- range(coords$UMAP2, na.rm = TRUE)
  x_pad <- max(diff(x_range) * 0.04, 0.1)
  y_pad <- max(diff(y_range) * 0.04, 0.1)
  x_limits <- x_range + c(-x_pad, x_pad)
  y_limits <- y_range + c(-y_pad, y_pad)
  centroids <- coords[, .(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), by = group_label]
  p_condition <- ggplot2::ggplot(coords, ggplot2::aes(UMAP1, UMAP2, color = group_label)) +
    ggplot2::geom_point(size = 0.32, alpha = 0.62) +
    ggplot2::scale_color_manual(values = group_colors, guide = "none") +
    ggplot2::labs(
      title = "Colored by condition",
      subtitle = "Labels mark condition centroids",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::coord_fixed(xlim = x_limits, ylim = y_limits, ratio = 1, clip = "off") +
    common_theme
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_condition <- p_condition + ggrepel::geom_label_repel(
      data = centroids,
      ggplot2::aes(label = group_label),
      color = "black",
      fill = "white",
      fontface = "bold",
      family = "Helvetica",
      size = 2.25,
      box.padding = 0.28,
      point.padding = 0.2,
      segment.color = "white",
      max.overlaps = Inf,
      seed = seed,
      show.legend = FALSE
    )
  } else {
    .log_warn("Package ggrepel is unavailable; condition UMAP labels cannot be repelled.")
    p_condition <- p_condition + ggplot2::geom_label(
      data = centroids,
      ggplot2::aes(label = group_label),
      color = "black",
      fill = "white",
      fontface = "bold",
      size = 2.25,
      linewidth = 0.25,
      show.legend = FALSE
    )
  }

  coords[, primary_topic := factor(primary_topic, levels = topic_levels)]
  topic_centroids <- coords[, .(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2)), by = primary_topic]
  p_topic <- ggplot2::ggplot(coords[order(primary_theta)], ggplot2::aes(UMAP1, UMAP2, color = primary_topic)) +
    ggplot2::geom_point(size = 0.32, alpha = 0.62) +
    ggplot2::scale_color_manual(values = topic_colors, drop = FALSE, guide = "none") +
    ggplot2::labs(
      title = "Colored by primary topic",
      subtitle = "Labels mark primary-topic centroids",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::coord_fixed(xlim = x_limits, ylim = y_limits, ratio = 1, clip = "off") +
    common_theme
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_topic <- p_topic + ggrepel::geom_label_repel(
      data = topic_centroids,
      ggplot2::aes(label = primary_topic),
      color = "black",
      fill = "white",
      fontface = "bold",
      family = "Helvetica",
      size = 2.25,
      box.padding = 0.28,
      point.padding = 0.2,
      segment.color = "white",
      max.overlaps = Inf,
      seed = seed,
      show.legend = FALSE
    )
  } else {
    p_topic <- p_topic + ggplot2::geom_label(
      data = topic_centroids,
      ggplot2::aes(label = primary_topic),
      color = "black",
      fill = "white",
      fontface = "bold",
      size = 2.25,
      linewidth = 0.25,
      show.legend = FALSE
    )
  }

  combined <- gridExtra::arrangeGrob(
    p_condition,
    p_topic,
    ncol = 2L,
    widths = c(1, 1),
    top = grid::textGrob(title_base, gp = grid::gpar(fontfamily = "Helvetica", fontface = "bold", fontsize = 14))
  )
  grDevices::pdf(out_file, width = 16, height = 9, family = "Helvetica", useDingbats = FALSE, bg = "white")
  device_open <- TRUE
  on.exit(if (isTRUE(device_open)) grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(combined)
  grDevices::dev.off()
  device_open <- FALSE
  invisible(out_file)
}

.select_tf_topic_primary_docs <- function(primary, max_docs = 120L) {
  dt <- data.table::as.data.table(primary)
  if (!nrow(dt)) return(character())
  max_docs <- as.integer(max_docs)
  if (!is.finite(max_docs) || max_docs <= 0L) max_docs <- 120L
  topics <- sort(unique(dt$primary_topic_num[is.finite(dt$primary_topic_num)]))
  if (!length(topics)) return(head(dt$doc_id, max_docs))
  per_topic <- max(1L, ceiling(max_docs / length(topics)))
  picked <- dt[order(primary_topic_num, -primary_theta, comparison_id, tf_display), head(.SD, per_topic), by = primary_topic_num]
  if (nrow(picked) > max_docs) {
    picked <- picked[order(primary_topic_num, -primary_theta, comparison_id, tf_display)][seq_len(max_docs)]
  }
  picked$doc_id
}

.plot_tf_topic_assignment_heatmap <- function(assign,
                                              out_file,
                                              title_prefix = NULL,
                                              max_docs = 120L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
    return(invisible(NULL))
  }
  membership <- data.table::as.data.table(assign$membership)
  primary <- data.table::as.data.table(assign$primary)
  if (!nrow(membership) || !nrow(primary)) return(invisible(NULL))
  keep_docs <- .select_tf_topic_primary_docs(primary, max_docs = max_docs)
  plot_dt <- membership[doc_id %in% keep_docs]
  plot_dt[, doc_label := .tf_topic_doc_label(comparison_id, tf_display, direction)]
  doc_levels <- unique(primary[doc_id %in% keep_docs][order(primary_topic_num, -primary_theta, comparison_id, tf_display), .tf_topic_doc_label(comparison_id, tf_display, direction)])
  plot_dt[, doc_label := factor(doc_label, levels = rev(doc_levels))]
  plot_dt[, topic := factor(as.character(topic), levels = unique(as.character(topic[order(topic_num)])))]
  title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "TF-topic assignment heatmap", sep = " | ")
  } else {
    "TF-topic assignment heatmap"
  }
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = topic, y = doc_label, fill = theta)) +
    ggplot2::geom_tile(color = "grey88", linewidth = 0.08) +
    ggplot2::geom_point(
      data = plot_dt[membership_pass == TRUE],
      ggplot2::aes(x = topic, y = doc_label),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.4,
      stroke = 0.35,
      fill = "black",
      color = "black"
    ) +
    ggplot2::scale_fill_gradientn(colors = c("#f7fbff", "#9ecae1", "#08519c"), limits = c(0, 1), oob = scales::squish, name = "theta") +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Topic-balanced sample of ", length(keep_docs), " TF-direction documents. Dots mark theta >= cutoff."),
      x = "Topic",
      y = "TF direction document"
    ) +
    ggplot2::theme_minimal(base_size = 8.5, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 5.5, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = 9.8, height = max(6, min(16, length(keep_docs) * 0.12 + 2)), limitsize = FALSE)
  invisible(out_file)
}

.plot_tf_primary_topic_dotplot <- function(primary,
                                           out_file,
                                           title_prefix = NULL,
                                           max_docs = 120L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
    return(invisible(NULL))
  }
  dt <- data.table::as.data.table(primary)
  if (!nrow(dt)) return(invisible(NULL))
  max_docs <- as.integer(max_docs)
  if (!is.finite(max_docs) || max_docs <= 0L) max_docs <- 120L
  keep_docs <- .select_tf_topic_primary_docs(dt, max_docs = max_docs)
  dt <- dt[doc_id %in% keep_docs]
  dt <- dt[order(primary_topic_num, -primary_theta, comparison_id, tf_display)]
  dt[, doc_label := .tf_topic_doc_label(comparison_id, tf_display, direction)]
  dt[, doc_label := factor(doc_label, levels = rev(unique(doc_label)))]
  dt[, primary_topic := factor(as.character(primary_topic), levels = unique(as.character(primary_topic[order(primary_topic_num)])))]
  title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "primary TF-topic labels", sep = " | ")
  } else {
    "Primary TF-topic labels"
  }
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = primary_topic, y = doc_label)) +
    ggplot2::geom_point(ggplot2::aes(size = primary_theta, fill = direction, shape = ambiguous), color = "black", stroke = 0.28, alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(1.2, 5.2), limits = c(0, 1), name = "Primary theta") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), name = "Ambiguous") +
    ggplot2::scale_fill_manual(values = c("Target-Up" = "#b2182b", "Target-Down" = "#2166ac", "FP-Up" = "#ef8a62", "FP-Down" = "#67a9cf"), na.value = "#636363", name = "Direction") +
    ggplot2::labs(
      title = title,
      subtitle = "One primary topic is assigned per TF-direction document.",
      x = "Primary topic",
      y = "TF direction document"
    ) +
    ggplot2::theme_minimal(base_size = 8.5, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 5.5, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 7, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = 9.8, height = max(6, min(16, nrow(dt) * 0.12 + 2)), limitsize = FALSE)
  invisible(out_file)
}

.plot_tf_direction_topic_map <- function(direction_summary,
                                         out_file,
                                         title_prefix = NULL,
                                         max_pairs = 120L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("scales", quietly = TRUE)) {
    return(invisible(NULL))
  }
  dt <- data.table::as.data.table(direction_summary)
  if (!nrow(dt)) return(invisible(NULL))
  max_pairs <- as.integer(max_pairs)
  if (!is.finite(max_pairs) || max_pairs <= 0L) max_pairs <- 120L
  dt[, max_theta := pmax(up_primary_theta, down_primary_theta, na.rm = TRUE)]
  dt[!is.finite(max_theta), max_theta := 0]
  dt <- dt[order(direction_topic_status != "different_topic", comparison_id, -max_theta, tf_display)]
  per_status_comp <- 4L
  dt <- dt[, head(.SD, per_status_comp), by = .(direction_topic_status, comparison_id)]
  dt <- dt[order(direction_topic_status != "different_topic", comparison_id, -max_theta, tf_display)]
  if (nrow(dt) > max_pairs) dt <- dt[seq_len(max_pairs)]
  dt[, row_label := paste(comparison_id, tf_display, sep = " | ")]
  dt[, row_label := factor(row_label, levels = rev(unique(row_label)))]
  long <- data.table::rbindlist(list(
    dt[, .(row_label, comparison_id, tf_display, side = "Target-Up", topic = up_primary_topic, theta = up_primary_theta, status = direction_topic_status)],
    dt[, .(row_label, comparison_id, tf_display, side = "Target-Down", topic = down_primary_topic, theta = down_primary_theta, status = direction_topic_status)]
  ), use.names = TRUE, fill = TRUE)
  long <- long[!is.na(topic) & nzchar(topic)]
  long[, side_x := ifelse(side == "Target-Up", 1, 2)]
  title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "Up/Down TF-topic map", sep = " | ")
  } else {
    "Up/Down TF-topic map"
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = dt[!is.na(up_primary_topic) & !is.na(down_primary_topic)],
      ggplot2::aes(x = 1, xend = 2, y = row_label, yend = row_label, color = direction_topic_status),
      linewidth = 0.38,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = side_x, y = row_label, size = theta, fill = topic),
      shape = 21,
      color = "black",
      stroke = 0.28,
      alpha = 0.92
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), limits = c(0, 1), name = "Primary theta") +
    ggplot2::scale_color_manual(values = c(same_topic = "#4d9221", different_topic = "#762a83", up_only = "#b2182b", down_only = "#2166ac", ambiguous = "#636363"), name = "Status") +
    ggplot2::scale_x_continuous(breaks = c(1, 2), labels = c("Target-Up", "Target-Down"), limits = c(0.75, 2.25)) +
    ggplot2::labs(
      title = title,
      subtitle = "Each row keeps Target-Up and Target-Down as separate topic assignments.",
      x = NULL,
      y = "Comparison and TF"
    ) +
    ggplot2::theme_minimal(base_size = 8.5, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "black"),
      axis.text.y = ggplot2::element_text(size = 5.5, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(size = 8, face = "bold", color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, p, width = 10.2, height = max(6, min(16, nrow(dt) * 0.12 + 2)), limitsize = FALSE)
  invisible(out_file)
}

.write_tf_topic_assignment_outputs <- function(theta,
                                               out_dir,
                                               doc_design = c("comparison", "condition"),
                                               membership_cutoff = 0.3,
                                               primary_margin_cutoff = 0.1) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  assign <- .build_tf_topic_assignment_tables(
    theta = theta,
    doc_design = doc_design,
    membership_cutoff = membership_cutoff,
    primary_margin_cutoff = primary_margin_cutoff
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(assign$pass, file.path(out_dir, "tf_topic_membership_pass.csv"))
  data.table::fwrite(assign$primary, file.path(out_dir, "tf_topic_primary.csv"))
  data.table::fwrite(assign$direction_summary, file.path(out_dir, "tf_direction_topic_summary.csv"))
  tf_list <- unique(assign$primary[, .(tf, tf_display)])
  data.table::setorder(tf_list, tf_display, tf)
  data.table::fwrite(tf_list, file.path(out_dir, "tf_topic_tf_list.csv"))
  assign_flat <- merge(
    assign$pass,
    assign$primary[, .(doc_id, tf, tf_display, primary_topic, primary_theta, margin, ambiguous)],
    by = c("doc_id", "tf", "tf_display"),
    all.x = TRUE,
    sort = FALSE
  )
  data.table::setorder(assign_flat, tf_display, doc_id, topic_num)
  data.table::fwrite(assign_flat, file.path(out_dir, "tf_topic_assignments.csv"))
  summary <- data.table::data.table(
    metric = c(
      "n_documents",
      "n_tfs",
      "n_membership_rows",
      "n_pass_membership_rows",
      "n_ambiguous_documents",
      "membership_cutoff",
      "primary_margin_cutoff"
    ),
    value = c(
      data.table::uniqueN(assign$primary$doc_id),
      data.table::uniqueN(assign$primary$tf_display),
      nrow(assign$membership),
      nrow(assign$pass),
      sum(assign$primary$ambiguous, na.rm = TRUE),
      membership_cutoff,
      primary_margin_cutoff
    )
  )
  data.table::fwrite(summary, file.path(out_dir, "tf_topic_assignment_summary.csv"))
  invisible(assign)
}

.prepare_tf_doc_topic_assignment_heatmap_data <- function(theta,
                                                          phi,
                                                          score_mat = NULL,
                                                          topic_terms = NULL,
                                                          doc_design = c("comparison", "condition"),
                                                          aggregate_fun = c("max", "mean")) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  aggregate_fun <- match.arg(aggregate_fun)
  theta <- as.matrix(theta)
  phi <- as.matrix(phi)
  if (is.null(rownames(theta))) .log_abort("theta must have rownames containing doc_id.")
  if (is.null(colnames(theta))) colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  if (is.null(rownames(phi))) rownames(phi) <- colnames(theta)
  if (is.null(colnames(phi))) .log_abort("phi must have term_id column names.")

  doc_info <- .parse_doc_id(rownames(theta), doc_design = doc_design)
  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, doc_id := rownames(theta)]
  theta_dt <- cbind(data.table::data.table(doc_id = rownames(theta)), doc_info, theta_dt[, setdiff(names(theta_dt), "doc_id"), with = FALSE])
  theta_dt <- theta_dt[!is.na(tf_doc) & nzchar(tf_doc) & !is.na(comparison_id) & nzchar(comparison_id)]
  if (!nrow(theta_dt)) return(data.table::data.table())
  topic_cols <- setdiff(names(theta_dt), c("doc_id", "comparison_id", "tf_doc", "direction"))
  topic_cols <- topic_cols[grepl("^Topic", topic_cols)]
  if (!length(topic_cols)) return(data.table::data.table())
  theta_long <- data.table::melt(
    theta_dt,
    id.vars = c("doc_id", "comparison_id", "tf_doc", "direction"),
    measure.vars = topic_cols,
    variable.name = "topic",
    value.name = "value"
  )
  theta_long[, value := as.numeric(value)]
  theta_long <- theta_long[is.finite(value)]
  theta_long[, topic_num := as.integer(gsub("^Topic", "", as.character(topic)))]
  theta_long[, value_raw := value]
  theta_long[, topic_max := max(value, na.rm = TRUE), by = topic]
  theta_long[!is.finite(topic_max) | topic_max <= 0, topic_max := 1]
  theta_long[, value := value / topic_max]
  theta_long[!is.finite(value) | value < 0, value := 0]
  theta_long[value > 1, value := 1]
  theta_long[, topic_max := NULL]
  theta_long[, tf := as.character(tf_doc)]
  if (identical(doc_design, "condition")) {
    theta_long[, page_label := comparison_id]
  } else {
    theta_long[, page_label := ifelse(!is.na(direction) & nzchar(direction), paste(comparison_id, direction, sep = "::"), comparison_id)]
  }
  theta_long[, panel := "TF doc score"]
  theta_long[, term_pass := FALSE]

  agg_fun <- if (identical(aggregate_fun, "max")) max else mean
  theta_agg <- theta_long[, .(value = agg_fun(value, na.rm = TRUE)), by = .(tf, topic, topic_num)]
  theta_agg[, `:=`(
    doc_id = NA_character_,
    comparison_id = "Aggregate",
    direction = NA_character_,
    page_label = "Aggregate",
    panel = "TF doc score",
    term_pass = FALSE
  )]

  if (is.null(score_mat)) {
    score_mat <- score_terms_normtop(phi)
  } else {
    score_mat <- as.matrix(score_mat)
  }
  if (!all(topic_cols %in% rownames(score_mat))) {
    keep_topics <- intersect(topic_cols, rownames(score_mat))
  } else {
    keep_topics <- topic_cols
  }
  tf_vals <- sort(unique(theta_long$tf))
  term_ids <- paste0("GENE:", tf_vals)
  keep_terms <- intersect(term_ids, colnames(score_mat))
  term_dt <- data.table::data.table()
  if (length(keep_topics) && length(keep_terms)) {
    sub_score <- score_mat[keep_topics, keep_terms, drop = FALSE]
    term_dt <- data.table::as.data.table(as.table(sub_score))
    data.table::setnames(term_dt, c("topic", "term_id", "value"))
    term_dt[, `:=`(
      topic = as.character(topic),
      term_id = as.character(term_id),
      value = as.numeric(value),
      tf = sub("^GENE:", "", as.character(term_id)),
      topic_num = as.integer(gsub("^Topic", "", as.character(topic)))
    )]
  }
  pass_dt <- data.table::data.table(topic = character(), term_id = character(), term_pass = logical())
  if (!is.null(topic_terms) && is.data.frame(topic_terms) && nrow(topic_terms)) {
    pass_dt <- data.table::as.data.table(topic_terms)
    if (!"topic" %in% names(pass_dt) && "topic_num" %in% names(pass_dt)) {
      pass_dt[, topic := paste0("Topic", as.integer(topic_num))]
    }
    if ("topic_num" %in% names(pass_dt)) {
      pass_dt[, topic := paste0("Topic", as.integer(topic_num))]
    } else if ("topic" %in% names(pass_dt)) {
      pass_dt[, topic := as.character(topic)]
      pass_dt[grepl("^[0-9]+$", topic), topic := paste0("Topic", topic)]
    }
    if ("in_topic" %in% names(pass_dt)) {
      pass_dt[, term_pass := .as_logical_flag(in_topic)]
    } else {
      pass_dt[, term_pass := TRUE]
    }
    pass_dt <- pass_dt[term_pass == TRUE & term_id %in% term_ids, .(topic = as.character(topic), term_id = as.character(term_id), term_pass = TRUE)]
    pass_dt <- unique(pass_dt)
  }
  if (nrow(term_dt)) {
    term_dt <- merge(term_dt, pass_dt, by = c("topic", "term_id"), all.x = TRUE, sort = FALSE)
    term_dt[is.na(term_pass), term_pass := FALSE]
    pages <- unique(c("Aggregate", theta_long$page_label))
    term_base <- data.table::copy(term_dt)
    term_dt <- data.table::rbindlist(lapply(pages, function(pg) {
      x <- data.table::copy(term_base)
      x[, page_label := pg]
      x
    }), use.names = TRUE, fill = TRUE)
    term_dt[, `:=`(
      doc_id = NA_character_,
      comparison_id = page_label,
      direction = NA_character_,
      panel = "TF term score"
    )]
  }

  out <- data.table::rbindlist(
    list(theta_agg, theta_long, term_dt),
    use.names = TRUE,
    fill = TRUE
  )
  out[, `:=`(tf = as.character(tf), topic = as.character(topic))]
  out <- out[!is.na(tf) & nzchar(tf) & !is.na(topic) & nzchar(topic)]
  if (!nrow(out)) return(out)

  topic_order_dt <- theta_agg[, .(max_doc_score = max(value, na.rm = TRUE)), by = .(topic, topic_num)]
  data.table::setorder(topic_order_dt, -max_doc_score, topic_num)
  topic_order <- topic_order_dt$topic
  tf_order_dt <- theta_agg[, .(max_doc_score = max(value, na.rm = TRUE)), by = .(tf, topic_num)]
  data.table::setorder(tf_order_dt, tf, -max_doc_score, topic_num)
  tf_order_dt <- tf_order_dt[, .SD[1L], by = tf]
  data.table::setorder(tf_order_dt, topic_num, -max_doc_score, tf)
  tf_order <- tf_order_dt$tf
  page_order <- c("Aggregate", setdiff(unique(theta_long$page_label), "Aggregate"))

  out[, topic := factor(as.character(topic), levels = topic_order)]
  out[, tf := factor(as.character(tf), levels = rev(tf_order))]
  out[, page_label := factor(as.character(page_label), levels = page_order)]
  out[, panel := factor(as.character(panel), levels = c("TF doc score", "TF term score"))]
  out[]
}

plot_tf_doc_topic_assignment_heatmaps <- function(theta,
                                                  phi,
                                                  out_file,
                                                  score_mat = NULL,
                                                  topic_terms = NULL,
                                                  title_prefix = NULL,
                                                  doc_design = c("comparison", "condition"),
                                                  aggregate_fun = c("max", "mean")) {
  .assert_pkg("data.table")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .log_inform("Skipping TF topic assignment heatmap: {.pkg ggplot2} not installed.")
    return(invisible(NULL))
  }
  doc_design <- match.arg(doc_design)
  aggregate_fun <- match.arg(aggregate_fun)
  plot_dt <- .prepare_tf_doc_topic_assignment_heatmap_data(
    theta = theta,
    phi = phi,
    score_mat = score_mat,
    topic_terms = topic_terms,
    doc_design = doc_design,
    aggregate_fun = aggregate_fun
  )
  if (!nrow(plot_dt)) return(invisible(NULL))
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  n_tf <- length(levels(plot_dt$tf))
  n_topic <- length(levels(plot_dt$topic))
  width <- max(8, min(22, n_topic * 0.28 + 4.5))
  height <- max(6, min(36, n_tf * 0.11 + 2.6))
  row_size <- if (n_tf > 180) 3.2 else if (n_tf > 100) 4.2 else if (n_tf > 60) 5.5 else 7
  col_size <- if (n_topic > 40) 5 else if (n_topic > 25) 6 else 8
  title_main <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "TF topic assignment", sep = " | ")
  } else {
    "TF topic assignment"
  }
  caption_txt <- paste(
    "Left panel: TF document theta normalized within each topic.",
    "Right panel: normalized GENE:<TF> topic-term score; black dots mark gammafit in-topic calls.",
    "Aggregate page uses max TF doc score across pages."
  )

  grDevices::pdf(out_file, width = width, height = height, family = "Helvetica", onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (pg in levels(plot_dt$page_label)) {
    sub <- plot_dt[page_label == pg]
    if (!nrow(sub)) next
    pass_sub <- sub[panel == "TF term score" & term_pass == TRUE & value > 0]
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = topic, y = tf, fill = value)) +
      ggplot2::geom_tile(color = "grey86", linewidth = 0.08) +
      ggplot2::facet_grid(. ~ panel) +
      ggplot2::scale_fill_gradient(low = "white", high = "#3a78af", name = "Value", limits = c(0, 1), oob = scales::squish) +
      ggplot2::labs(
        title = title_main,
        subtitle = as.character(pg),
        x = "Topic",
        y = "TF",
        caption = caption_txt
      ) +
      ggplot2::theme_minimal(base_family = "Helvetica", base_size = 9) +
      ggplot2::theme(
        text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold", size = 9),
        axis.text.x = ggplot2::element_text(face = "bold", size = col_size, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = ggplot2::element_text(face = "bold", size = row_size),
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 9),
        plot.title = ggplot2::element_text(face = "bold", size = 10),
        plot.subtitle = ggplot2::element_text(face = "bold", size = 9),
        plot.caption = ggplot2::element_text(face = "bold", size = 7),
        legend.title = ggplot2::element_text(face = "bold", size = 8),
        legend.text = ggplot2::element_text(face = "bold", size = 7)
      )
    if (nrow(pass_sub)) {
      p <- p + ggplot2::geom_point(
        data = pass_sub,
        ggplot2::aes(x = topic, y = tf),
        inherit.aes = FALSE,
        shape = 21,
        size = 0.85,
        stroke = 0.25,
        color = "black",
        fill = "black"
      )
    }
    print(p)
  }
  invisible(out_file)
}

.direction_group <- function(direction) {
  direction <- as.character(direction)
  out <- rep(NA_character_, length(direction))
  up <- grepl("up", direction, ignore.case = TRUE)
  down <- grepl("down", direction, ignore.case = TRUE)
  out[up] <- "Up"
  out[down] <- "Down"
  out[!up & !down] <- direction[!up & !down]
  out[is.na(out) | !nzchar(out)] <- "Unknown"
  out
}

.as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    lx <- tolower(trimws(x))
    return(lx %in% c("true", "t", "1", "yes", "y"))
  }
  rep(FALSE, length(x))
}

.get_cellline_from_comparison <- function(x) sub("_.*$", "", x)

.set_grob_fontface <- function(grob, fontface = "bold") {
  if (!is.null(grob$gp)) {
    if (!is.null(grob$gp$font)) grob$gp$font <- NULL
    grob$gp$fontface <- fontface
  }
  if (!is.null(grob$children) && length(grob$children)) {
    grob$children <- lapply(grob$children, .set_grob_fontface, fontface = fontface)
  }
  grob
}

.zscale_rowcol <- function(mat) {
  mat <- as.matrix(mat)
  if (!nrow(mat) || !ncol(mat)) return(mat)
  mat <- scale(mat)
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  mat
}

plot_doc_topic_heatmaps_by_comparison <- function(theta,
                                                  out_dir,
                                                  edges_docs = NULL,
                                                  title_prefix = NULL,
                                                  option_label = NULL,
                                                  doc_design = c("comparison", "condition")) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_inform("Skipping doc-topic heatmaps: {.pkg pheatmap} not installed.")
    return(invisible(NULL))
  }
  if (is.null(rownames(theta))) {
    .log_inform("Skipping doc-topic heatmaps: theta has no rownames.")
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  doc_info <- .parse_doc_id(rownames(theta), doc_design = doc_design)
  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, doc_id := rownames(theta)]
  dt <- data.table::data.table(doc_id = rownames(theta))
  dt <- cbind(dt, doc_info, theta_dt[, setdiff(names(theta_dt), "doc_id"), with = FALSE])

  topic_cols <- setdiff(names(dt), c("doc_id", "comparison_id", "tf_doc", "direction"))
  if (!length(topic_cols)) return(invisible(NULL))
  is_cmpr_mode <- all(is.na(dt$tf_doc)) || !any(nzchar(dt$tf_doc))

  ed <- NULL
  if (!is.null(edges_docs)) {
    ed <- data.table::as.data.table(edges_docs)
  }

  if (is_cmpr_mode) {
    if (is.null(ed)) {
      .log_inform("Skipping doc-topic heatmaps: edges_docs required for comparison-doc mode.")
      return(invisible(NULL))
    }
    weight_col <- if (option_label == "opt3_gene_fc_expr") "fc_mag_gene" else if (option_label == "opt2_peak_fc_fp") "fc_mag_fp" else "delta_fp"
    term_col <- if (option_label == "opt3_gene_fc_expr") "gene_key" else "peak_id"
    if (!(weight_col %in% names(ed)) || !(term_col %in% names(ed)) || !"doc_id" %in% names(ed)) {
      .log_inform("Skipping doc-topic heatmaps: missing columns for comparison-doc mode.")
      return(invisible(NULL))
    }
    ed <- ed[!is.na(get(term_col)) & nzchar(get(term_col))]
    ed[, term_id := as.character(get(term_col))]
    ed[, w := abs(.safe_num(get(weight_col)))]
    ed <- ed[is.finite(w) & w > 0]
    if (!nrow(ed)) return(invisible(NULL))
    ed_unique <- ed[, .(weight = max(w, na.rm = TRUE)), by = .(doc_id, comparison_id, tf, term_id)]
    tf_w <- ed_unique[, .(tf_weight = sum(weight, na.rm = TRUE)), by = .(doc_id, comparison_id, tf)]
    tf_w[, direction := .parse_doc_id(doc_id, doc_design = doc_design)$direction]

    theta_dt <- data.table::as.data.table(theta)
    theta_dt[, doc_id := rownames(theta)]

    for (comp in unique(na.omit(tf_w$comparison_id))) {
      tf_comp <- tf_w[comparison_id == comp]
      if (!nrow(tf_comp)) next
      dirs <- unique(na.omit(tf_comp$direction))
      if (!length(dirs)) dirs <- NA_character_
      for (dir in dirs) {
        tf_sub <- if (is.na(dir)) tf_comp else tf_comp[direction == dir]
        if (!nrow(tf_sub)) next
        tf_sub[, doc_weight_total := sum(tf_weight, na.rm = TRUE), by = doc_id]
        tf_sub[!is.finite(doc_weight_total) | doc_weight_total <= 0, doc_weight_total := 1]
        tf_sub[, tf_weight_norm := tf_weight / doc_weight_total]

        mat_doc <- merge(tf_sub, theta_dt, by = "doc_id", all.x = TRUE)
        if (!nrow(mat_doc)) next
        for (col in topic_cols) {
          mat_doc[, (col) := get(col) * tf_weight_norm]
        }
        agg <- mat_doc[, lapply(.SD, sum, na.rm = TRUE), by = tf, .SDcols = topic_cols]
        mat <- as.matrix(agg[, ..topic_cols])
        rownames(mat) <- agg$tf
        if (!nrow(mat) || !ncol(mat)) next

        has_dir <- length(dir) == 1L && !is.na(dir) && nzchar(dir)
        name_part <- if (has_dir) paste(comp, dir, sep = "__") else comp
        out_file <- file.path(out_dir, sprintf("%s_tf_topic_heatmap.pdf", .safe_filename(name_part)))
        font_row <- if (nrow(mat) > 80) 6 else if (nrow(mat) > 40) 8 else 11
        font_col <- if (ncol(mat) > 40) 8 else 12
        width <- max(7, ncol(mat) * 0.12)
        height <- max(6, nrow(mat) * 0.12)
        main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
          if (has_dir) paste(title_prefix, comp, dir, sep = " | ") else paste(title_prefix, comp, sep = " | ")
        } else {
          if (has_dir) paste(comp, dir, sep = " | ") else comp
        }

        grDevices::pdf(out_file, width = width, height = height)
        tryCatch(
          {
            ph <- pheatmap::pheatmap(
              mat,
              cluster_rows = nrow(mat) > 1L,
              cluster_cols = ncol(mat) > 1L,
              show_rownames = TRUE,
              show_colnames = TRUE,
              fontsize_row = font_row,
              fontsize_col = font_col,
              main = main_title,
              border_color = NA,
              silent = TRUE
            )
            if (!is.null(ph$gtable)) {
              idx <- which(ph$gtable$layout$name %in% c("row_names", "col_names"))
              for (i in idx) {
                ph$gtable$grobs[[i]] <- .set_grob_fontface(ph$gtable$grobs[[i]], "bold")
              }
              grid::grid.newpage()
              grid::grid.draw(ph$gtable)
            } else {
              pheatmap::pheatmap(mat, border_color = NA)
            }
          },
          finally = grDevices::dev.off()
        )
      }
    }
    return(invisible(TRUE))
  }

  for (comp in unique(na.omit(dt$comparison_id))) {
    comp_sub <- dt[comparison_id == comp]
    if (!nrow(comp_sub)) next
    dirs <- unique(na.omit(comp_sub$direction))
    if (!length(dirs)) dirs <- NA_character_
    for (dir in dirs) {
      sub <- if (is.na(dir)) comp_sub else comp_sub[direction == dir]
      if (!nrow(sub)) next
      agg <- sub[, lapply(.SD, mean, na.rm = TRUE), by = tf_doc, .SDcols = topic_cols]
      if (!is.null(ed) && "tf_doc" %in% names(ed)) {
        if (identical(doc_design, "condition") && "condition_label" %in% names(ed)) {
          tf_keep <- unique(ed[condition_label == comp, tf_doc])
        } else if ("comparison_id" %in% names(ed)) {
          if ("direction" %in% names(ed) && !is.na(dir)) {
            tf_keep <- unique(ed[comparison_id == comp & direction == dir, tf_doc])
          } else {
            tf_keep <- unique(ed[comparison_id == comp, tf_doc])
          }
        } else {
          tf_keep <- character()
        }
        agg <- agg[tf_doc %in% tf_keep]
      }
      mat <- as.matrix(agg[, ..topic_cols])
      rownames(mat) <- agg$tf_doc
      if (!nrow(mat) || !ncol(mat)) next

      has_dir <- length(dir) == 1L && !is.na(dir) && nzchar(dir)
      name_part <- if (has_dir) paste(comp, dir, sep = "__") else comp
      out_file <- file.path(out_dir, sprintf("%s_tf_topic_heatmap.pdf", .safe_filename(name_part)))
      font_row <- if (nrow(mat) > 80) 6 else if (nrow(mat) > 40) 8 else 11
      font_col <- if (ncol(mat) > 40) 8 else 12
      width <- max(7, ncol(mat) * 0.1)
      height <- max(6, nrow(mat) * 0.12)
      main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
        if (has_dir) paste(title_prefix, comp, dir, sep = " | ") else paste(title_prefix, comp, sep = " | ")
      } else {
        if (has_dir) paste(comp, dir, sep = " | ") else comp
      }

      grDevices::pdf(out_file, width = width, height = height)
      tryCatch(
        {
          cluster_rows <- nrow(mat) > 1L
          cluster_cols <- ncol(mat) > 1L
          ph <- pheatmap::pheatmap(
            mat,
            cluster_rows = cluster_rows,
            cluster_cols = cluster_cols,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize_row = font_row,
            fontsize_col = font_col,
            main = main_title,
            border_color = NA,
            silent = TRUE
          )
          if (!is.null(ph$gtable)) {
            idx <- which(ph$gtable$layout$name %in% c("row_names", "col_names"))
            for (i in idx) {
              ph$gtable$grobs[[i]] <- .set_grob_fontface(ph$gtable$grobs[[i]], "bold")
            }
            grid::grid.newpage()
            grid::grid.draw(ph$gtable)
          } else {
            pheatmap::pheatmap(mat, border_color = NA)
          }
        },
        finally = grDevices::dev.off()
      )
    }
  }

  invisible(TRUE)
}

plot_tf_topic_heatmaps_from_link_scores <- function(link_scores,
                                                    out_dir,
                                                    title_prefix = NULL,
                                                    value_col = c("prob", "score"),
                                                    min_value = 0,
                                                    per_comparison = TRUE,
                                                    split_direction = TRUE) {
  .assert_pkg("data.table")
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_inform("Skipping link-score TF-topic heatmaps: {.pkg pheatmap} not installed.")
    return(invisible(NULL))
  }

  dt <- data.table::as.data.table(link_scores)
  if (!nrow(dt)) return(invisible(NULL))
  if (!("topic_num" %in% names(dt))) {
    if (!("topic" %in% names(dt))) {
      .log_abort("link_scores must have topic_num or topic.")
    }
    dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
  }
  dt <- dt[is.finite(topic_num)]
  if (!nrow(dt)) return(invisible(NULL))
  dt <- dt[!is.na(tf) & nzchar(tf)]
  if ("gene_key" %in% names(dt)) {
    dt <- dt[!is.na(gene_key) & nzchar(gene_key)]
  }
  if (!nrow(dt)) return(invisible(NULL))
  if (!("gene_key" %in% names(dt))) {
    .log_abort("link_scores missing gene_key for TF->gene link heatmaps.")
  }
  dt[, link_id := paste(tf, gene_key, sep = "::")]
  if (!nrow(dt)) return(invisible(NULL))

  value_col <- match.arg(value_col)
  if (!value_col %in% names(dt)) {
    value_col <- if ("prob" %in% names(dt)) "prob" else if ("score" %in% names(dt)) "score" else value_col
  }
  if (!value_col %in% names(dt)) {
    .log_abort("link_scores missing value column: {value_col}.")
  }
  dt <- dt[is.finite(get(value_col))]
  if (is.finite(min_value) && min_value > 0) {
    dt <- dt[get(value_col) >= min_value]
  }
  if (!nrow(dt)) return(invisible(NULL))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(out_dir, "link_scores_doc_topic_heatmaps_debug.txt")
  log_msg <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
  }

  doc_info <- .parse_doc_id(dt$doc_id)
  dt <- cbind(dt, doc_info)
  dt[, direction_group := if (isTRUE(split_direction)) .direction_group(direction) else "All"]
  dt <- dt[!is.na(comparison_id) & nzchar(comparison_id)]
  if (!nrow(dt)) {
    log_msg("No comparison_id parsed from doc_id; skipping TF-topic heatmaps.")
    return(invisible(NULL))
  }

  if (isTRUE(per_comparison)) {
    comps <- unique(dt$comparison_id)
  } else {
    dt[, comparison_id := "All"]
    dt[, direction_group := "All"]
    comps <- "All"
  }
  for (cmp in comps) {
    cmp_dt <- dt[comparison_id == cmp]
    if (!nrow(cmp_dt)) next
    dirs <- unique(cmp_dt$direction_group)
    for (dir_lab in dirs) {
      sub_dt <- cmp_dt[direction_group == dir_lab]
      if (!nrow(sub_dt)) next

      # Aggregate TF->gene link -> topic strength as max membership across links.
      agg <- sub_dt[, .(value = max(get(value_col), na.rm = TRUE)), by = .(link_id, topic_num)]
      if (!nrow(agg)) next
      mat_dt <- data.table::dcast(agg, link_id ~ topic_num, value.var = "value", fill = 0)
      topic_cols <- setdiff(names(mat_dt), "link_id")
      topic_nums <- suppressWarnings(as.integer(topic_cols))
      ord <- order(topic_nums)
      topic_cols <- topic_cols[ord]
      mat <- as.matrix(mat_dt[, ..topic_cols])
      rownames(mat) <- mat_dt$link_id
      colnames(mat) <- paste0("Topic", topic_nums[ord])
      if (!nrow(mat) || !ncol(mat)) next

      width <- max(7, ncol(mat) * 0.12)
      height <- max(6, nrow(mat) * 0.12)
      font_row <- if (nrow(mat) > 80) 6 else if (nrow(mat) > 40) 8 else 11
      font_col <- if (ncol(mat) > 40) 8 else 12
      main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
        paste(title_prefix, cmp, dir_lab, sep = " | ")
      } else {
        paste(cmp, dir_lab, sep = " | ")
      }

      out_file <- file.path(
        out_dir,
        sprintf("%s_tf_topic_heatmap.pdf", .safe_filename(paste(cmp, dir_lab, sep = "__")))
      )
      log_msg(sprintf("Writing %s (%s, %s)", out_file, cmp, dir_lab))
      grDevices::pdf(out_file, width = width, height = height)
      tryCatch(
        {
          ph <- pheatmap::pheatmap(
            mat,
            cluster_rows = nrow(mat) > 1L,
            cluster_cols = ncol(mat) > 1L,
            show_rownames = TRUE,
            show_colnames = TRUE,
            fontsize_row = font_row,
            fontsize_col = font_col,
            main = main_title,
            border_color = NA,
            silent = TRUE
          )
          if (!is.null(ph$gtable)) {
            idx <- which(ph$gtable$layout$name %in% c("row_names", "col_names"))
            for (i in idx) {
              ph$gtable$grobs[[i]] <- .set_grob_fontface(ph$gtable$grobs[[i]], "bold")
            }
            grid::grid.newpage()
            grid::grid.draw(ph$gtable)
          } else {
            pheatmap::pheatmap(mat, border_color = NA)
          }
        },
        finally = grDevices::dev.off()
      )
    }
  }

  invisible(TRUE)
}

.resolve_topic_comparison_labels <- function(dt,
                                             edges_docs = NULL,
                                             label_cleaner = NULL) {
  .assert_pkg("data.table")
  dt <- data.table::as.data.table(dt)
  if (!"direction_label" %in% names(dt)) dt[, direction_label := NA_character_]
  dt[, comparison_display := as.character(comparison_id)]
  dt[, doc_display_label := NA_character_]

  if (!is.null(edges_docs) && is.data.frame(edges_docs) && nrow(edges_docs)) {
    ed <- data.table::as.data.table(edges_docs)
    if ("comparison_label" %in% names(ed) && "doc_id" %in% names(ed)) {
      ed[, comparison_label := as.character(comparison_label)]
      doc_map <- ed[
        !is.na(doc_id) & nzchar(as.character(doc_id)),
        .(
          comparison_display = {
            x <- comparison_label[!is.na(comparison_label) & nzchar(trimws(comparison_label))]
            if (length(x)) x[[1L]] else NA_character_
          },
          doc_display_label = {
            if ("doc_display_label" %in% names(ed)) {
              x <- as.character(doc_display_label)
              x <- x[!is.na(x) & nzchar(trimws(x))]
              if (length(x)) x[[1L]] else NA_character_
            } else {
              NA_character_
            }
          }
        ),
        by = doc_id
      ]
      dt <- merge(dt, doc_map, by = "doc_id", all.x = TRUE, suffixes = c("", "_from_edges"))
      if ("comparison_display_from_edges" %in% names(dt)) {
        use <- !is.na(dt$comparison_display_from_edges) & nzchar(trimws(dt$comparison_display_from_edges))
        dt[use, comparison_display := comparison_display_from_edges]
        dt[, comparison_display_from_edges := NULL]
      }
      if ("doc_display_label_from_edges" %in% names(dt)) {
        use <- !is.na(dt$doc_display_label_from_edges) & nzchar(trimws(dt$doc_display_label_from_edges))
        dt[use, doc_display_label := doc_display_label_from_edges]
        dt[, doc_display_label_from_edges := NULL]
      }
    } else if ("comparison_label" %in% names(ed) && "comparison_id" %in% names(ed)) {
      ed[, comparison_id := as.character(comparison_id)]
      ed[, comparison_label := as.character(comparison_label)]
      comp_map <- ed[
        !is.na(comparison_id) & nzchar(comparison_id),
        .(
          comparison_display = {
            x <- comparison_label[!is.na(comparison_label) & nzchar(trimws(comparison_label))]
            if (length(x)) x[[1L]] else NA_character_
          }
        ),
        by = comparison_id
      ]
      dt <- merge(dt, comp_map, by = "comparison_id", all.x = TRUE, suffixes = c("", "_from_edges"))
      if ("comparison_display_from_edges" %in% names(dt)) {
        use <- !is.na(dt$comparison_display_from_edges) & nzchar(trimws(dt$comparison_display_from_edges))
        dt[use, comparison_display := comparison_display_from_edges]
        dt[, comparison_display_from_edges := NULL]
      }
    }
  }

  dt[, comparison_label := ifelse(
    !is.na(direction_label) & nzchar(direction_label),
    paste(comparison_display, direction_label, sep = "::"),
    comparison_display
  )]
  has_doc_display <- !is.na(dt$doc_display_label) & nzchar(trimws(dt$doc_display_label))
  same_direction <- is.na(dt$direction_label) | !nzchar(dt$direction_label) |
    endsWith(dt$doc_display_label, paste0("::", dt$direction_label))
  dt[has_doc_display & same_direction, comparison_label := doc_display_label]
  if (!is.null(label_cleaner)) {
    dt[, comparison_label := label_cleaner(comparison_label)]
  }
  dt
}

plot_topic_by_comparison_heatmaps <- function(theta,
                                              out_dir,
                                              edges_docs = NULL,
                                              direction_mode = c("gene", "fp"),
                                              title_prefix = NULL,
                                              label_cleaner = NULL,
                                              topic_links = NULL,
                                              annotate_unique_genes = TRUE,
                                              doc_design = c("comparison", "condition")) {
  .assert_pkg("data.table")
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_inform("Skipping topic-by-comparison heatmaps: {.pkg pheatmap} not installed.")
    return(invisible(NULL))
  }
  if (isTRUE(annotate_unique_genes) && !is.null(topic_links) && !requireNamespace("gtable", quietly = TRUE)) {
    .log_inform("Skipping topic-by-comparison gene-count annotation: {.pkg gtable} not installed.")
    annotate_unique_genes <- FALSE
  }
  if (is.null(rownames(theta))) {
    .log_inform("Skipping topic-by-comparison heatmaps: theta has no rownames.")
    return(invisible(NULL))
  }

  direction_mode <- match.arg(direction_mode)
  doc_design <- match.arg(doc_design)
  doc_info <- .parse_doc_id(rownames(theta), doc_design = doc_design)
  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, doc_id := rownames(theta)]
  dt <- data.table::data.table(doc_id = rownames(theta))
  dt <- cbind(dt, doc_info, theta_dt[, setdiff(names(theta_dt), "doc_id"), with = FALSE])

  topic_cols <- setdiff(names(dt), c("doc_id", "comparison_id", "tf_doc", "direction"))
  if (!length(topic_cols)) return(invisible(NULL))

  .topic_gene_annotation <- function(topic_cols_use) {
    if (!isTRUE(annotate_unique_genes) || is.null(topic_links) || !is.data.frame(topic_links) || !nrow(topic_links)) {
      return(NULL)
    }
    tl <- data.table::as.data.table(topic_links)
    req <- c("topic_num", "gene_key", "peak_pass", "gene_pass")
    if (!all(req %in% names(tl))) return(NULL)
    tl[, topic_num := as.integer(topic_num)]
    cnt <- tl[
      .as_logical_flag(peak_pass) & .as_logical_flag(gene_pass) & !is.na(gene_key) & nzchar(gene_key),
      .(unique_genes = data.table::uniqueN(gene_key)),
      by = topic_num
    ]
    ann <- data.frame("Unique genes" = rep(0L, length(topic_cols_use)), check.names = FALSE)
    rownames(ann) <- topic_cols_use
    if (nrow(cnt)) {
      idx <- match(paste0("Topic", cnt$topic_num), rownames(ann))
      keep <- !is.na(idx)
      ann[idx[keep], "Unique genes"] <- cnt$unique_genes[keep]
    }
    ann
  }

  .draw_topic_pheatmap <- function(mat, out_file, main_title, width, height) {
    ann_row <- .topic_gene_annotation(rownames(mat))
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(out_file, width = width, height = height, family = "Helvetica")
    device_id <- grDevices::dev.cur()
    device_open <- TRUE
    on.exit({
      open_devices <- grDevices::dev.list()
      if (isTRUE(device_open) && !is.null(open_devices) && device_id %in% open_devices) {
        grDevices::dev.off(device_id)
      }
    }, add = TRUE)
    ph <- pheatmap::pheatmap(
      mat,
      cluster_rows = nrow(mat) > 1L,
      cluster_cols = ncol(mat) > 1L,
      show_rownames = TRUE,
      show_colnames = ncol(mat) <= 50,
      fontsize = 9,
      fontsize_row = 9,
      fontsize_col = 9,
      main = main_title,
      border_color = NA,
      annotation_row = ann_row,
      annotation_colors = list("Unique genes" = grDevices::colorRampPalette(c("#f7fcf5", "#74c476", "#006d2c"))(100)),
      silent = TRUE
    )
    if (!is.null(ann_row) && nrow(ann_row)) {
      row_order <- if (!is.null(ph$tree_row)) rownames(mat)[ph$tree_row$order] else rownames(mat)
      values <- ann_row[row_order, "Unique genes"]
      labels <- as.character(values)
      ann_pos <- ph$gtable$layout[ph$gtable$layout$name == "row_annotation", ]
      if (nrow(ann_pos) == 1L) {
        ph$gtable$widths[ann_pos$l] <- grid::unit(0.65, "in")
        ann_pal <- grDevices::colorRampPalette(c("#f7fcf5", "#74c476", "#006d2c"))(100)
        value_rng <- range(values, na.rm = TRUE)
        if (!all(is.finite(value_rng)) || diff(value_rng) == 0) {
          color_idx <- rep(50L, length(values))
        } else {
          color_idx <- as.integer(round((values - value_rng[1]) / diff(value_rng) * 99)) + 1L
        }
        color_idx <- pmin(100L, pmax(1L, color_idx))
        y_pos <- (rev(seq_along(labels)) - 0.5) / length(labels)
        tile_grob <- grid::rectGrob(
          x = grid::unit(0.5, "npc"),
          y = grid::unit(y_pos, "npc"),
          width = grid::unit(1, "npc"),
          height = grid::unit(1 / length(labels), "npc"),
          gp = grid::gpar(fill = ann_pal[color_idx], col = NA)
        )
        label_col <- ifelse(values >= max(values, na.rm = TRUE) * 0.55, "white", "black")
        label_grob <- grid::textGrob(
          labels,
          x = grid::unit(0.5, "npc"),
          y = grid::unit(y_pos, "npc"),
          gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica", col = label_col)
        )
        ph$gtable <- gtable::gtable_add_grob(
          ph$gtable,
          tile_grob,
          t = ann_pos$t,
          l = ann_pos$l,
          b = ann_pos$b,
          r = ann_pos$r,
          name = "row_annotation_gene_count_tiles",
          clip = "off"
        )
        ph$gtable <- gtable::gtable_add_grob(
          ph$gtable,
          label_grob,
          t = ann_pos$t,
          l = ann_pos$l,
          b = ann_pos$b,
          r = ann_pos$r,
          name = "row_annotation_gene_count_labels",
          clip = "off"
        )
      }
    }
    .bold_text_grobs <- function(g) {
      if (inherits(g, "text")) {
        gp_list <- as.list(g$gp)
        gp_list$font <- NULL
        g$gp <- do.call(grid::gpar, modifyList(gp_list, list(fontface = "bold", fontfamily = "Helvetica", fontsize = 9)))
      }
      if (!is.null(g$grobs)) {
        g$grobs <- lapply(g$grobs, .bold_text_grobs)
      }
      if (!is.null(g$children)) {
        g$children <- do.call(grid::gList, lapply(g$children, .bold_text_grobs))
      }
      g
    }
    ph$gtable <- .bold_text_grobs(ph$gtable)
    grDevices::dev.set(device_id)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    grDevices::dev.off(device_id)
    device_open <- FALSE
    invisible(out_file)
  }

  dt[, direction_label := direction]
  if (identical(doc_design, "condition")) {
    dt[, direction_label := NA_character_]
  }
  if (identical(doc_design, "comparison") && direction_mode == "fp" && !is.null(edges_docs)) {
    ed <- data.table::as.data.table(edges_docs)
    if ("doc_id" %in% names(ed)) {
      fp_sign <- .safe_sign(ed[["log2fc_fp"]])
      if ("delta_fp" %in% names(ed)) {
        alt <- .safe_sign(ed[["delta_fp"]])
        fp_sign[fp_sign == 0L] <- alt[fp_sign == 0L]
      }
      ed[, fp_dir_sum := fp_sign]
      fp_dir <- ed[, .(fp_dir = sum(fp_dir_sum, na.rm = TRUE)), by = doc_id]
      fp_dir[, fp_dir := ifelse(fp_dir > 0, "FP-Up", ifelse(fp_dir < 0, "FP-Down", NA_character_))]
      dt <- merge(dt, fp_dir, by = "doc_id", all.x = TRUE)
      dt[!is.na(fp_dir), direction_label := fp_dir]
    }
  }
  dt <- .resolve_topic_comparison_labels(dt, edges_docs = edges_docs, label_cleaner = label_cleaner)

  comp_avg <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(comparison_id, comparison_label), .SDcols = topic_cols]
  comp_avg[, cellline := .get_cellline_from_comparison(comparison_id)]

  cells <- unique(na.omit(comp_avg$cellline))
  if (length(cells) > 1L) {
    mat <- t(as.matrix(comp_avg[, ..topic_cols]))
    colnames(mat) <- comp_avg$comparison_label
    rownames(mat) <- topic_cols
    out_file <- file.path(out_dir, "All_topic_by_comparison.pdf")
    main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
      paste(title_prefix, "All comparisons", sep = " | ")
    } else {
      "All comparisons"
    }
    .draw_topic_pheatmap(mat, out_file, main_title, width = max(7, ncol(mat) * 0.35), height = max(7, nrow(mat) * 0.25))
  } else {
    for (cell in cells) {
      sub <- comp_avg[cellline == cell]
      if (!nrow(sub)) next
      mat <- t(as.matrix(sub[, ..topic_cols]))
      colnames(mat) <- sub$comparison_label
      rownames(mat) <- topic_cols
      out_file <- file.path(out_dir, sprintf("%s_topic_by_comparison.pdf", .safe_filename(cell)))
      main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
        paste(title_prefix, cell, sep = " | ")
      } else {
        cell
      }
      .draw_topic_pheatmap(mat, out_file, main_title, width = max(7, ncol(mat) * 0.35), height = max(7, nrow(mat) * 0.25))
    }
  }
  invisible(TRUE)
}

.tf_from_theta_by_topic <- function(theta,
                                    edges_docs,
                                    topics,
                                    top_n = 50L,
                                    min_theta = NA_real_) {
  if (is.null(theta) || is.null(edges_docs)) return(data.table::data.table())
  ed <- data.table::as.data.table(edges_docs)
  if (!all(c("doc_id", "tf") %in% names(ed))) return(data.table::data.table())
  ed <- ed[!is.na(doc_id) & nzchar(doc_id) & !is.na(tf) & nzchar(tf)]
  if (!nrow(ed)) return(data.table::data.table())

  theta_mat <- as.matrix(theta)
  doc_ids <- rownames(theta_mat)
  if (is.null(doc_ids) || anyNA(doc_ids) || any(doc_ids == "")) {
    doc_ids <- unique(ed$doc_id)
    if (length(doc_ids) != nrow(theta_mat)) return(data.table::data.table())
  }
  rownames(theta_mat) <- doc_ids

  topics <- sort(unique(as.integer(topics)))
  topics <- topics[is.finite(topics)]
  if (!length(topics)) return(data.table::data.table())

  topic_cols <- colnames(theta_mat)
  if (!is.null(topic_cols) && all(grepl("^Topic\\d+$", topic_cols))) {
    col_ids <- as.integer(sub("^Topic", "", topic_cols))
    col_map <- stats::setNames(seq_along(col_ids), col_ids)
    topic_idx <- unname(col_map[as.character(topics)])
  } else {
    topic_idx <- topics
  }

  out <- vector("list", length(topics))
  for (i in seq_along(topics)) {
    idx <- topic_idx[[i]]
    if (is.na(idx) || idx < 1L || idx > ncol(theta_mat)) next
    vec <- theta_mat[, idx]
    doc_dt <- data.table::data.table(doc_id = doc_ids, theta = as.numeric(vec))
    doc_dt <- doc_dt[is.finite(theta)]
    if (is.finite(min_theta)) doc_dt <- doc_dt[theta >= min_theta]
    if (!nrow(doc_dt)) next
    data.table::setorder(doc_dt, -theta)
    if (is.finite(top_n) && top_n > 0L && nrow(doc_dt) > top_n) {
      doc_dt <- doc_dt[seq_len(as.integer(top_n))]
    }
    tf_dt <- merge(doc_dt, ed[, .(doc_id, tf)], by = "doc_id", allow.cartesian = TRUE)
    if (!nrow(tf_dt)) next
    tf_dt <- tf_dt[!is.na(tf) & nzchar(tf)]
    if (!nrow(tf_dt)) next
    tf_dt[, topic := topics[[i]]]
    out[[i]] <- tf_dt[, .(score = max(theta, na.rm = TRUE)), by = .(topic, gene = tf)]
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

topic_gene_sets <- function(topic_terms,
                            edges_docs,
                            option_label,
                            use_all_terms = FALSE,
                            theta = NULL,
                            tf_link_mode = c("theta", "none"),
                            tf_theta_top_n = 50L,
                            tf_theta_min = NA_real_) {
  .assert_pkg("data.table")
  tf_link_mode <- match.arg(tf_link_mode)
  tt <- data.table::as.data.table(topic_terms)
  if (!nrow(tt)) return(list())
  if (isTRUE(use_all_terms)) {
    if (!("term_id" %in% names(tt))) return(list())
    tt <- tt[!is.na(term_id) & nzchar(term_id)]
  } else {
    if (!("in_topic" %in% names(tt))) return(list())
    in_set_vec <- .as_logical_flag(tt$in_topic)
    tt <- tt[in_set_vec]
  }
  if (!nrow(tt)) return(list())

  if (option_label == "opt3_gene_fc_expr") {
    tt[, gene_key := gsub("^GENE:", "", term_id)]
    genes_dt <- tt[!is.na(gene_key) & nzchar(gene_key), .(gene = unique(gene_key)), by = topic]
  } else if (option_label == "joint") {
    if (is.null(edges_docs)) return(list())
    ed <- data.table::as.data.table(edges_docs)
    if (!("peak_id" %in% names(ed)) || !("gene_key" %in% names(ed))) return(list())
    ed <- ed[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key)]
    tt_gene <- tt[grepl("^GENE:", term_id)]
    tt_gene[, gene_key := gsub("^GENE:", "", term_id)]
    genes_a <- tt_gene[!is.na(gene_key) & nzchar(gene_key), .(gene = unique(gene_key)), by = topic]
    tt_peak <- tt[grepl("^PEAK:", term_id)]
    tt_peak[, peak_id := gsub("^PEAK:", "", term_id)]
    peak_dt <- tt_peak[!is.na(peak_id) & nzchar(peak_id), .(topic, peak_id)]
    genes_b <- merge(peak_dt, ed[, .(peak_id, gene = gene_key)], by = "peak_id", allow.cartesian = TRUE)
    genes_b <- genes_b[, .(gene = unique(gene)), by = topic]
    genes_dt <- data.table::rbindlist(list(genes_a, genes_b), use.names = TRUE, fill = TRUE)
    genes_dt <- genes_dt[, .(gene = unique(gene)), by = topic]
  } else {
    if (is.null(edges_docs)) return(list())
    ed <- data.table::as.data.table(edges_docs)
    if (!("peak_id" %in% names(ed)) || !("gene_key" %in% names(ed))) return(list())
    tt[, peak_id := gsub("^PEAK:", "", term_id)]
    peak_dt <- tt[!is.na(peak_id) & nzchar(peak_id), .(topic, peak_id)]
    ed <- ed[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key)]
    genes_dt <- merge(peak_dt, ed[, .(peak_id, gene = gene_key)], by = "peak_id", allow.cartesian = TRUE)
    genes_dt <- genes_dt[, .(gene = unique(gene)), by = topic]
  }

  if (tf_link_mode == "theta") {
    tf_dt <- .tf_from_theta_by_topic(
      theta = theta,
      edges_docs = edges_docs,
      topics = unique(genes_dt$topic),
      top_n = tf_theta_top_n,
      min_theta = tf_theta_min
    )
    if (nrow(tf_dt)) {
      tf_dt <- tf_dt[, .(gene = unique(gene)), by = topic]
      genes_dt <- data.table::rbindlist(list(genes_dt, tf_dt), use.names = TRUE, fill = TRUE)
      genes_dt <- genes_dt[, .(gene = unique(gene)), by = topic]
    }
  }

  if (!nrow(genes_dt)) return(list())
  out <- split(genes_dt$gene, genes_dt$topic)
  out <- lapply(out, function(x) unique(x[!is.na(x) & nzchar(x)]))
  out[lengths(out) > 0]
}

topic_gene_sets_from_terms <- function(topic_terms,
                                       edges_docs = NULL,
                                       option_label = c("opt1_peak_delta_fp", "opt2_peak_fc_fp",
                                                        "opt3_gene_fc_expr", "joint"),
                                       use_all_terms = FALSE,
                                       include_peak_terms = TRUE) {
  .assert_pkg("data.table")
  option_label <- match.arg(option_label)
  tt <- data.table::as.data.table(topic_terms)
  if (!nrow(tt) || !("term_id" %in% names(tt))) return(list())
  if (!("topic_num" %in% names(tt))) {
    if ("topic" %in% names(tt)) {
      tt[, topic_num := as.integer(gsub("^Topic", "", as.character(topic)))]
    } else {
      .log_abort("topic_terms must contain topic_num or topic.")
    }
  }
  if (!isTRUE(use_all_terms)) {
    if (!("in_topic" %in% names(tt))) return(list())
    tt <- tt[.as_logical_flag(in_topic)]
  }
  tt <- tt[is.finite(topic_num) & !is.na(term_id) & nzchar(term_id)]
  if (!nrow(tt)) return(list())

  gene_dt <- tt[grepl("^GENE:", term_id), .(
    topic_num = as.integer(topic_num),
    gene = sub("^GENE:", "", as.character(term_id))
  )]

  if (isTRUE(include_peak_terms) && !identical(option_label, "opt3_gene_fc_expr")) {
    peak_dt <- tt[grepl("^PEAK:", term_id), .(
      topic_num = as.integer(topic_num),
      peak_or_gene = sub("^PEAK:", "", as.character(term_id))
    )]
    if (nrow(peak_dt)) {
      peak_gene <- data.table::data.table()
      if (!is.null(edges_docs)) {
        ed <- data.table::as.data.table(edges_docs)
        if (all(c("peak_id", "gene_key") %in% names(ed))) {
          ed <- unique(ed[
            !is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key),
            .(peak_or_gene = as.character(peak_id), gene = as.character(gene_key))
          ])
          peak_gene <- merge(peak_dt, ed, by = "peak_or_gene", allow.cartesian = TRUE)
          peak_gene <- peak_gene[, .(topic_num, gene)]
        }
      }
      mapped <- if (nrow(peak_gene)) unique(ed$peak_or_gene) else character()
      direct_peak_gene <- peak_dt[
        !(peak_or_gene %in% mapped) & !is.na(peak_or_gene) & nzchar(peak_or_gene),
        .(topic_num, gene = peak_or_gene)
      ]
      gene_dt <- data.table::rbindlist(
        list(gene_dt, peak_gene, direct_peak_gene),
        use.names = TRUE,
        fill = TRUE
      )
    }
  }

  gene_dt <- gene_dt[!is.na(gene) & nzchar(gene)]
  if (!nrow(gene_dt)) return(list())
  gene_dt <- unique(gene_dt[, .(topic_num, gene)])
  out <- split(gene_dt$gene, gene_dt$topic_num)
  out <- lapply(out, function(x) unique(as.character(x)))
  out[lengths(out) > 0]
}

topic_gene_sets_by_comparison_terms <- function(topic_terms,
                                                edges_docs,
                                                theta,
                                                theta_min = 0.3,
                                                include_peak_terms = TRUE,
                                                use_all_terms = FALSE,
                                                doc_design = c("comparison", "condition"),
                                                split_direction = TRUE) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  tt <- data.table::as.data.table(topic_terms)
  ed <- data.table::as.data.table(edges_docs)
  theta <- as.matrix(theta)
  if (!nrow(tt) || !nrow(ed) || !nrow(theta) || !ncol(theta)) {
    return(data.table::data.table())
  }
  if (is.null(rownames(theta))) {
    .log_abort("theta must have document row names for stratified topic-term pathway gene sets.")
  }
  if (!all(c("doc_id", "gene_key") %in% names(ed))) {
    .log_abort("edges_docs must contain doc_id and gene_key.")
  }
  if (!("term_id" %in% names(tt))) {
    .log_abort("topic_terms must contain term_id.")
  }
  if (!("topic_num" %in% names(tt))) {
    if ("topic" %in% names(tt)) {
      tt[, topic_num := as.integer(gsub("^Topic", "", as.character(topic)))]
    } else {
      .log_abort("topic_terms must contain topic_num or topic.")
    }
  }
  if (!isTRUE(use_all_terms)) {
    if (!("in_topic" %in% names(tt))) {
      .log_abort("topic_terms must contain in_topic unless use_all_terms = TRUE.")
    }
    tt <- tt[.as_logical_flag(in_topic)]
  }
  tt <- tt[is.finite(topic_num) & !is.na(term_id) & nzchar(term_id)]
  if (!nrow(tt)) return(data.table::data.table())

  gene_terms <- tt[grepl("^GENE:", term_id), .(
    topic = as.integer(topic_num),
    gene = sub("^GENE:", "", as.character(term_id)),
    term_source = "GENE"
  )]
  if (isTRUE(include_peak_terms)) {
    peak_terms <- tt[grepl("^PEAK:", term_id), .(
      topic = as.integer(topic_num),
      gene = sub("^PEAK:", "", as.character(term_id)),
      term_source = "PEAK"
    )]
    term_genes <- data.table::rbindlist(list(gene_terms, peak_terms), use.names = TRUE, fill = TRUE)
  } else {
    term_genes <- gene_terms
  }
  term_genes <- term_genes[!is.na(gene) & nzchar(gene)]
  term_genes <- unique(term_genes, by = c("topic", "gene", "term_source"))
  if (!nrow(term_genes)) return(data.table::data.table())

  theta_min <- .safe_num(theta_min)
  if (!is.finite(theta_min)) theta_min <- 0.3
  theta_colnames <- colnames(theta)
  topic_ids <- if (!is.null(theta_colnames) && all(grepl("^Topic\\d+$", theta_colnames))) {
    as.integer(sub("^Topic", "", theta_colnames))
  } else {
    seq_len(ncol(theta))
  }
  theta_dt <- data.table::as.data.table(theta)
  theta_dt[, doc_id := rownames(theta)]
  theta_long <- data.table::melt(
    theta_dt,
    id.vars = "doc_id",
    variable.name = "topic_col",
    value.name = "theta"
  )
  theta_long[, topic_idx__ := match(as.character(topic_col), theta_colnames)]
  theta_long[, topic := topic_ids[topic_idx__]]
  theta_long[, theta := .safe_num(theta)]
  theta_long <- theta_long[is.finite(topic) & is.finite(theta) & theta >= theta_min]
  theta_long[, topic_idx__ := NULL]
  if (!nrow(theta_long)) return(data.table::data.table())

  info <- .parse_doc_id(theta_long$doc_id, doc_design = doc_design)
  theta_long <- cbind(theta_long, info)
  if (identical(doc_design, "comparison")) {
    theta_long[, direction_group := if (isTRUE(split_direction)) .direction_group(direction) else "All"]
  } else {
    theta_long[, direction_group := "All"]
  }
  if (!"tf" %in% names(theta_long)) {
    if ("tf_doc" %in% names(theta_long)) {
      theta_long[, tf := tf_doc]
    } else {
      theta_long[, tf := NA_character_]
    }
  }
  theta_long <- theta_long[!is.na(comparison_id) & nzchar(comparison_id)]
  if (!nrow(theta_long)) return(data.table::data.table())

  doc_genes <- unique(ed[
    !is.na(doc_id) & nzchar(doc_id) & !is.na(gene_key) & nzchar(gene_key),
    .(doc_id, gene = as.character(gene_key))
  ])
  doc_topic_genes <- merge(
    theta_long[, .(doc_id, topic, theta, comparison_id, direction, direction_group, tf)],
    doc_genes,
    by = "doc_id",
    allow.cartesian = TRUE
  )
  if (!nrow(doc_topic_genes)) return(data.table::data.table())

  out <- merge(
    doc_topic_genes,
    term_genes,
    by = c("topic", "gene"),
    allow.cartesian = TRUE
  )
  if (!nrow(out)) return(data.table::data.table())
  out <- out[, .(
    n_docs = data.table::uniqueN(doc_id),
    max_theta = max(theta, na.rm = TRUE),
    term_source = paste(sort(unique(term_source)), collapse = ";")
  ), by = .(comparison_id, direction_group, topic, gene)]
  data.table::setorder(out, comparison_id, direction_group, topic, gene)
  out[]
}

.split_pathway_gene_string <- function(x) {
  x <- as.character(x %||% "")
  if (!length(x) || is.na(x[[1L]]) || !nzchar(x[[1L]])) return(character())
  out <- unlist(strsplit(x[[1L]], "[;,]", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  unique(out[!is.na(out) & nzchar(out)])
}

.topic_pathway_retest_from_overall <- function(topic_terms,
                                               edges_docs,
                                               out_dir,
                                               option_label = c("opt1_peak_delta_fp", "opt2_peak_fc_fp",
                                                                "opt3_gene_fc_expr", "joint"),
                                               include_peak_terms = TRUE,
                                               use_all_terms = FALSE,
                                               per_comparison_dir = ".",
                                               split_direction = TRUE,
                                               background_size = 20000L,
                                               overall_pathway_file = NULL,
                                               pathway_species = NULL,
                                               doc_design = c("comparison", "condition")) {
  .assert_pkg("data.table")
  option_label <- match.arg(option_label)
  doc_design <- match.arg(doc_design)
  group_id_col <- if (identical(doc_design, "condition")) "condition_id" else "comparison_id"
  group_label <- if (identical(doc_design, "condition")) "condition" else "comparison"
  group_prefix <- if (identical(doc_design, "condition")) "per_condition" else "per_comparison"
  group_topic_genes_col <- if (identical(doc_design, "condition")) "condition_topic_genes" else "comparison_topic_genes"
  group_topic_symbols_col <- if (identical(doc_design, "condition")) "condition_topic_gene_symbols" else "comparison_topic_gene_symbols"
  pathway_species_mode <- .normalize_pathway_species_mode(pathway_species)
  human_mouse_best <- identical(pathway_species_mode, "human_mouse_best")
  key_species <- if (human_mouse_best) c("human", "mouse") else pathway_species_mode
  add_gene_keys <- function(dt) {
    dt <- data.table::as.data.table(dt)
    out <- lapply(key_species, function(sp) {
      tmp <- data.table::copy(dt)
      keys <- .gene_symbol_key_table(tmp$gene, species = sp)
      tmp[, `:=`(
        pathway_species = sp,
        gene_key__ = keys$gene_key__,
        gene_canonical = keys$gene_canonical,
        gene_match_type = keys$gene_match_type,
        gene_matched = keys$gene_matched,
        gene_ambiguous = keys$gene_ambiguous
      )]
      tmp
    })
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }
  out_dir_pc <- if (is.null(per_comparison_dir) || !nzchar(as.character(per_comparison_dir)[[1L]]) || identical(as.character(per_comparison_dir)[[1L]], ".")) {
    out_dir
  } else {
    file.path(out_dir, per_comparison_dir)
  }
  dir.create(out_dir_pc, recursive = TRUE, showWarnings = FALSE)
  if (is.null(overall_pathway_file)) {
    overall_pathway_file <- file.path(out_dir, "topic_pathway_enrichment_topic_terms.csv")
  }
  if (is.null(overall_pathway_file) || !file.exists(overall_pathway_file)) {
    .log_abort("Missing overall topic pathway table: {overall_pathway_file}. Run overall_topic_pathway first.")
  }

  topic_sets <- topic_gene_sets_from_terms(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    option_label = option_label,
    use_all_terms = use_all_terms,
    include_peak_terms = include_peak_terms
  )
  if (!length(topic_sets)) return(data.table::data.table())
  topic_gene_dt <- data.table::rbindlist(
    lapply(names(topic_sets), function(tp) {
      data.table::data.table(topic = as.integer(tp), gene = unique(as.character(topic_sets[[tp]])))
    }),
    use.names = TRUE,
    fill = TRUE
  )
  topic_gene_dt <- topic_gene_dt[is.finite(topic) & !is.na(gene) & nzchar(gene)]
  topic_gene_dt <- add_gene_keys(topic_gene_dt)
  topic_gene_dt <- topic_gene_dt[!is.na(gene_key__) & nzchar(gene_key__)]
  topic_gene_dt <- unique(topic_gene_dt)
  if (!nrow(topic_gene_dt)) return(data.table::data.table())

  ed <- data.table::as.data.table(edges_docs)
  .assert_has_cols(ed, c("doc_id", "gene_key"), context = sprintf("per-%s topic pathway retest", group_label))
  ed <- ed[!is.na(doc_id) & nzchar(doc_id) & !is.na(gene_key) & nzchar(gene_key)]
  if (!nrow(ed)) return(data.table::data.table())
  doc_info <- .parse_doc_id(ed$doc_id, doc_design = doc_design)
  ed <- cbind(ed, doc_info)
  if (identical(doc_design, "comparison")) {
    ed[, direction_group := if (isTRUE(split_direction)) .direction_group(direction) else "All"]
  } else {
    ed[, direction_group := "All"]
  }
  if (identical(doc_design, "condition")) {
    ed[, condition_id := comparison_id]
  }
  if (identical(doc_design, "condition") && "condition_label" %in% names(ed)) {
    comp_genes <- unique(ed[
      !is.na(get(group_id_col)) & nzchar(get(group_id_col)),
      .(
        group_id = get(group_id_col),
        condition_label = as.character(condition_label),
        direction_group,
        gene = as.character(gene_key)
      )
    ])
  } else {
    comp_genes <- unique(ed[
      !is.na(get(group_id_col)) & nzchar(get(group_id_col)),
      .(group_id = get(group_id_col), direction_group, gene = as.character(gene_key))
    ])
  }
  comp_genes <- comp_genes[!is.na(gene) & nzchar(gene)]
  comp_genes <- add_gene_keys(comp_genes)
  comp_genes <- comp_genes[!is.na(gene_key__) & nzchar(gene_key__)]
  if (!nrow(comp_genes)) return(data.table::data.table())

  query_gene_dt <- merge(
    comp_genes,
    unique(topic_gene_dt[, .(topic, pathway_species, gene_key__)]),
    by = c("pathway_species", "gene_key__"),
    allow.cartesian = TRUE,
    sort = FALSE
  )
  if (!nrow(query_gene_dt)) return(data.table::data.table())
  query_by_cols <- c("group_id", if ("condition_label" %in% names(query_gene_dt)) "condition_label", "direction_group", "topic", "pathway_species")
  query_dt <- query_gene_dt[, .(
    query_genes = list(sort(unique(gene))),
    query_gene_symbols = list(sort(unique(gene_canonical[!is.na(gene_canonical) & nzchar(gene_canonical)]))),
    query_gene_keys = list(sort(unique(gene_key__))),
    query_size = data.table::uniqueN(gene_key__),
    topic_gene_count = data.table::uniqueN(topic_gene_dt[topic == .BY$topic & pathway_species == .BY$pathway_species, gene_key__]),
    gene_match_summary = .gene_symbol_match_summary(gene_match_type)
  ), by = query_by_cols]
  query_dt <- query_dt[query_size > 0]
  if (!nrow(query_dt)) return(data.table::data.table())

  overall <- data.table::fread(overall_pathway_file, showProgress = FALSE)
  if (!nrow(overall)) return(data.table::data.table())
  .assert_has_cols(overall, c("topic", "pathway", "genes", "overlap"), context = "overall topic pathway table")
  overall[, topic := as.integer(topic)]
  overall <- overall[is.finite(topic) & !is.na(pathway) & nzchar(pathway)]
  if (!nrow(overall)) return(data.table::data.table())
  if (!"database" %in% names(overall)) {
    overall[, database := sub(":.*$", "", as.character(pathway))]
  }
  if (human_mouse_best) {
    if ("selected_pathway_species" %in% names(overall)) {
      overall[, pathway_species := as.character(selected_pathway_species)]
    } else if ("pathway_species" %in% names(overall)) {
      overall[, pathway_species := as.character(pathway_species)]
    } else {
      overall[, pathway_species := vapply(database, .pathway_database_species, character(1))]
    }
    overall <- overall[pathway_species %in% c("human", "mouse")]
  } else {
    overall[, pathway_species := pathway_species_mode]
  }
  if (!nrow(overall)) return(data.table::data.table())
  if (!"pathway_norm_key" %in% names(overall)) {
    db_label_for_key <- if ("database_label" %in% names(overall)) overall$database_label else sub(":.*$", "", as.character(overall$pathway))
    term_for_key <- if ("pathway_term" %in% names(overall)) overall$pathway_term else sub("^[^:]+:\\s*", "", as.character(overall$pathway))
    overall[, pathway_norm_key := mapply(
      .normalize_topic_pathway_term,
      pathway = pathway,
      database_label = db_label_for_key,
      pathway_term = term_for_key,
      USE.NAMES = FALSE
    )]
  }
  if (!"selected_pathway_species" %in% names(overall)) overall[, selected_pathway_species := pathway_species]
  if (!"selected_database" %in% names(overall)) overall[, selected_database := database]
  if (!"pval" %in% names(overall)) overall[, pval := NA_real_]
  if (!"padj" %in% names(overall)) overall[, padj := NA_real_]
  if (!"overlap_hits" %in% names(overall)) {
    overall[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
  }
  overlap_term_size <- suppressWarnings(as.integer(sub("^.*/", "", as.character(overall$overlap))))
  if (!"term_size" %in% names(overall)) {
    overall[, term_size := overlap_term_size]
  } else {
    overall[, term_size := suppressWarnings(as.integer(term_size))]
    overall[!is.finite(term_size) | is.na(term_size), term_size := overlap_term_size[.I]]
  }
  if (!"background_size" %in% names(overall)) {
    overall[, background_size := NA_integer_]
  }
  bg_default <- suppressWarnings(as.integer(background_size[[1L]]))
  if (!is.finite(bg_default) || bg_default < 1L) bg_default <- 20000L
  overall[, background_size := suppressWarnings(as.integer(background_size))]
  overall[!is.finite(background_size) | is.na(background_size), background_size := bg_default]
  overall[, term_size := as.integer(term_size)]
  overall[, background_size := as.integer(background_size)]
  if (any(is.finite(overall$term_size) & overall$term_size > overall$background_size, na.rm = TRUE)) {
    .log_abort("background_size must be at least as large as every pathway term.")
  }
  overall[, overall_overlap_genes_list := lapply(genes, .split_pathway_gene_string)]
  overall_keys <- Map(function(x, sp) {
    .gene_symbol_key_table(x, species = sp)
  }, overall$overall_overlap_genes_list, overall$pathway_species)
  overall_gene_key_dt <- data.table::rbindlist(overall_keys, use.names = TRUE, fill = TRUE)
  if (!all(c("gene", "gene_canonical", "gene_match_type", "gene_matched", "gene_ambiguous") %in% names(overall_gene_key_dt))) {
    overall_gene_key_dt <- data.table::data.table(
      gene = character(),
      gene_canonical = character(),
      gene_match_type = character(),
      gene_matched = logical(),
      gene_ambiguous = logical()
    )
  }
  .write_gene_symbol_conversion_audit(
    out_dir,
    list(
      topic_term_genes = unique(topic_gene_dt[, .(gene, gene_canonical, gene_match_type, gene_matched, gene_ambiguous)]),
      document_genes = unique(comp_genes[, .(gene, gene_canonical, gene_match_type, gene_matched, gene_ambiguous)]),
      overall_pathway_overlap_genes = unique(overall_gene_key_dt[, .(gene, gene_canonical, gene_match_type, gene_matched, gene_ambiguous)])
    )
  )
  overall[, overall_overlap_gene_keys_list := lapply(overall_keys, function(x) unique(x[!is.na(gene_key__) & nzchar(gene_key__), gene_key__]))]
  overall[, overall_overlap_gene_symbols_list := lapply(overall_keys, function(x) unique(x[!is.na(gene_canonical) & nzchar(gene_canonical), gene_canonical]))]
  overall[, overall_gene_match_summary := vapply(overall_keys, function(x) .gene_symbol_match_summary(x$gene_match_type), character(1))]
  overall[, overall_overlap_genes := vapply(overall_overlap_genes_list, paste, character(1), collapse = ";")]
  overall[, overall_overlap_gene_symbols := vapply(overall_overlap_gene_symbols_list, paste, character(1), collapse = ";")]
  overall[, overall_overlap_hits := lengths(overall_overlap_gene_keys_list)]
  overall_small <- overall[, .(
    topic,
    pathway_species,
    pathway_norm_key,
    pathway,
    database,
    selected_pathway_species,
    selected_database,
    overall_pval = as.numeric(pval),
    overall_padj = as.numeric(padj),
    overall_overlap = as.character(overlap),
    overall_overlap_hits = as.integer(overall_overlap_hits),
    overall_overlap_genes,
    overall_overlap_gene_symbols,
    overall_gene_match_summary,
    overall_overlap_genes_list,
    overall_overlap_gene_keys_list,
    overall_overlap_gene_symbols_list,
    term_size = as.integer(term_size),
    background_size = as.integer(background_size),
    overall_cluster_size = if ("cluster_size" %in% names(overall)) as.integer(cluster_size) else NA_integer_,
    overall_combined_score = if ("combined_score" %in% names(overall)) as.numeric(combined_score) else NA_real_,
    overall_odds_ratio = if ("odds_ratio" %in% names(overall)) as.numeric(odds_ratio) else NA_real_,
    human_padj = if ("human_padj" %in% names(overall)) as.numeric(human_padj) else NA_real_,
    mouse_padj = if ("mouse_padj" %in% names(overall)) as.numeric(mouse_padj) else NA_real_,
    human_logp = if ("human_logp" %in% names(overall)) as.numeric(human_logp) else NA_real_,
    mouse_logp = if ("mouse_logp" %in% names(overall)) as.numeric(mouse_logp) else NA_real_,
    human_overlap_hits = if ("human_overlap_hits" %in% names(overall)) as.integer(human_overlap_hits) else NA_integer_,
    mouse_overlap_hits = if ("mouse_overlap_hits" %in% names(overall)) as.integer(mouse_overlap_hits) else NA_integer_
  )]
  overall_small <- overall_small[is.finite(term_size) & term_size > 0 & is.finite(background_size) & background_size > 0]
  if (!nrow(overall_small)) return(data.table::data.table())

  res <- merge(
    query_dt,
    overall_small,
    by = c("topic", "pathway_species"),
    allow.cartesian = TRUE,
    sort = FALSE
  )
  if (!nrow(res)) return(data.table::data.table())
  res[, overlap_gene_key_list := Map(intersect, query_gene_keys, overall_overlap_gene_keys_list)]
  res[, overlap_gene_list := Map(function(genes, keys, keep_keys) {
    out <- genes[keys %in% keep_keys]
    sort(unique(out[!is.na(out) & nzchar(out)]))
  }, query_genes, query_gene_keys, overlap_gene_key_list)]
  res[, overlap_gene_symbol_list := Map(function(genes, keys, keep_keys) {
    out <- genes[keys %in% keep_keys]
    sort(unique(out[!is.na(out) & nzchar(out)]))
  }, query_gene_symbols, query_gene_keys, overlap_gene_key_list)]
  res[, overlap_hits := lengths(overlap_gene_key_list)]
  res[, overlap_genes := vapply(overlap_gene_list, paste, character(1), collapse = ";")]
  res[, overlap_gene_symbols := vapply(overlap_gene_symbol_list, paste, character(1), collapse = ";")]
  res[, pval := stats::phyper(
    q = overlap_hits - 1,
    m = term_size,
    n = background_size - term_size,
    k = query_size,
    lower.tail = FALSE
  )]
  res[!is.finite(pval), pval := NA_real_]
  res[, padj := stats::p.adjust(pval, method = "BH"), by = .(group_id, direction_group, topic, database)]
  res[, odds_ratio := {
    a <- as.numeric(overlap_hits)
    b <- as.numeric(query_size) - a
    c <- as.numeric(term_size) - a
    d <- as.numeric(background_size) - as.numeric(term_size) - b
    out <- (a * d) / (b * c)
    out[!is.finite(out)] <- Inf
    out
  }]
  res[, combined_score := odds_ratio * -log(pmax(pval, 1e-300))]
  res[, logp := -log10(pmax(padj, .Machine$double.xmin))]
  res[, overlap := paste0(overlap_hits, "/", term_size)]
  res[, (group_topic_genes_col) := vapply(query_genes, paste, character(1), collapse = ";")]
  res[, (group_topic_symbols_col) := vapply(query_gene_symbols, paste, character(1), collapse = ";")]
  res[, c(
    "query_genes", "query_gene_symbols", "query_gene_keys", "overall_overlap_genes_list",
    "overall_overlap_gene_keys_list", "overall_overlap_gene_symbols_list",
    "overlap_gene_key_list", "overlap_gene_list", "overlap_gene_symbol_list"
  ) := NULL]
  data.table::setnames(res, "group_id", group_id_col)
  cols <- c(
    group_id_col, if (identical(doc_design, "condition")) "condition_label", "direction_group", "topic", "pathway", "database",
    "pathway_species", "pathway_norm_key", "selected_pathway_species", "selected_database",
    "pval", "padj", "logp", "overlap", "overlap_hits", "overlap_genes",
    "overlap_gene_symbols",
    "query_size", "term_size", "background_size", "odds_ratio", "combined_score",
    group_topic_genes_col, group_topic_symbols_col, "topic_gene_count",
    "gene_match_summary", "overall_pval", "overall_padj",
    "overall_overlap", "overall_overlap_hits", "overall_overlap_genes",
    "overall_overlap_gene_symbols", "overall_gene_match_summary",
    "overall_cluster_size", "overall_combined_score", "overall_odds_ratio",
    "human_padj", "mouse_padj", "human_logp", "mouse_logp",
    "human_overlap_hits", "mouse_overlap_hits"
  )
  cols <- cols[cols %in% names(res)]
  data.table::setcolorder(res, cols)
  data.table::setorderv(res, c(group_id_col, "direction_group", "topic", "padj", "pval", "overlap_hits", "pathway"), order = c(1L, 1L, 1L, 1L, 1L, -1L, 1L))
  out_file <- file.path(out_dir_pc, paste0(group_prefix, "_topic_pathway_enrichment.csv"))
  data.table::fwrite(res, out_file)
  res[]
}

topic_gene_ranks <- function(topic_terms,
                             edges_docs,
                             option_label,
                             score_col = "score",
                             peak_gene_agg = c("max", "sum"),
                             theta = NULL,
                             tf_link_mode = c("theta", "none"),
                             tf_theta_top_n = 50L,
                             tf_theta_min = NA_real_) {
  .assert_pkg("data.table")
  peak_gene_agg <- match.arg(peak_gene_agg)
  tf_link_mode <- match.arg(tf_link_mode)
  tt <- data.table::as.data.table(topic_terms)
  if (!nrow(tt) || !(score_col %in% names(tt))) return(list())
  tt <- tt[!is.na(term_id) & nzchar(term_id)]
  tt[, score := .safe_num(get(score_col))]
  tt <- tt[is.finite(score)]
  if (!nrow(tt)) return(list())

  if (option_label == "opt3_gene_fc_expr") {
    tt[, gene := gsub("^GENE:", "", term_id)]
    genes_dt <- tt[!is.na(gene) & nzchar(gene), .(score = max(score, na.rm = TRUE)), by = .(topic, gene)]
  } else if (option_label == "joint") {
    if (is.null(edges_docs)) return(list())
    ed <- data.table::as.data.table(edges_docs)
    if (!all(c("peak_id", "gene_key") %in% names(ed))) return(list())
    ed <- ed[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key)]

    tt_gene <- tt[grepl("^GENE:", term_id)]
    tt_gene[, gene := gsub("^GENE:", "", term_id)]
    genes_a <- tt_gene[!is.na(gene) & nzchar(gene), .(score = max(score, na.rm = TRUE)), by = .(topic, gene)]

    tt_peak <- tt[grepl("^PEAK:", term_id)]
    tt_peak[, peak_id := gsub("^PEAK:", "", term_id)]
    peak_dt <- tt_peak[!is.na(peak_id) & nzchar(peak_id), .(topic, peak_id, score)]
    genes_b <- merge(peak_dt, ed[, .(peak_id, gene = gene_key)], by = "peak_id", allow.cartesian = TRUE)
    if (nrow(genes_b)) {
      genes_b <- genes_b[, .(
        score = if (peak_gene_agg == "sum") sum(score, na.rm = TRUE) else max(score, na.rm = TRUE)
      ), by = .(topic, gene)]
    }

    genes_dt <- data.table::rbindlist(list(genes_a, genes_b), use.names = TRUE, fill = TRUE)
    genes_dt <- genes_dt[!is.na(gene) & nzchar(gene)]
    genes_dt <- genes_dt[, .(
      score = if (peak_gene_agg == "sum") sum(score, na.rm = TRUE) else max(score, na.rm = TRUE)
    ), by = .(topic, gene)]
  } else {
    if (is.null(edges_docs)) return(list())
    ed <- data.table::as.data.table(edges_docs)
    if (!all(c("peak_id", "gene_key") %in% names(ed))) return(list())
    ed <- ed[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key)]

    tt[, peak_id := gsub("^PEAK:", "", term_id)]
    peak_dt <- tt[!is.na(peak_id) & nzchar(peak_id), .(topic, peak_id, score)]
    genes_dt <- merge(peak_dt, ed[, .(peak_id, gene = gene_key)], by = "peak_id", allow.cartesian = TRUE)
    genes_dt <- genes_dt[!is.na(gene) & nzchar(gene)]
    genes_dt <- genes_dt[, .(
      score = if (peak_gene_agg == "sum") sum(score, na.rm = TRUE) else max(score, na.rm = TRUE)
    ), by = .(topic, gene)]
  }

  if (tf_link_mode == "theta") {
    tf_dt <- .tf_from_theta_by_topic(
      theta = theta,
      edges_docs = edges_docs,
      topics = unique(tt$topic),
      top_n = tf_theta_top_n,
      min_theta = tf_theta_min
    )
    if (nrow(tf_dt)) {
      genes_dt <- data.table::rbindlist(list(genes_dt, tf_dt), use.names = TRUE, fill = TRUE)
      genes_dt <- genes_dt[!is.na(gene) & nzchar(gene)]
      genes_dt <- genes_dt[, .(score = max(score, na.rm = TRUE)), by = .(topic, gene)]
    }
  }

  if (!nrow(genes_dt)) return(list())
  out <- split(genes_dt, genes_dt$topic)
  out <- lapply(out, function(df) {
    vec <- df$score
    names(vec) <- df$gene
    vec <- vec[!is.na(names(vec)) & nzchar(names(vec))]
    vec <- vec[is.finite(vec)]
    vec <- sort(vec, decreasing = TRUE)
    vec
  })
  out[lengths(out) > 0]
}

# =============================================================================
# Link -> topic membership (posterior-style scoring)
# =============================================================================

compute_link_topic_scores <- function(edges_docs,
                                      theta,
                                      phi,
                                      topic_terms = NULL,
                                      out_file = NULL,
                                      gate_mode = c("none", "peak_in_set", "gene_in_set", "peak_and_gene_in_set"),
                                      top_k = 3L,
                                      min_prob = 0,
                                      include_tf = FALSE,
                                      overwrite = FALSE,
                                      chunk_size = 5000L,
                                      n_cores = 1L) {
  # Principle (baseline):
  # For a link (doc_id, PEAK, GENE, optional TF), we score topics by
  #   score_k ?theta[doc,k] * phi[k,PEAK] * phi[k,GENE] * (phi[k,TF] if include_tf)
  # and normalize over topics to get link -> topic posterior probabilities.
  #
  # Option C (hard confirmation gate):
  # If gate_mode != "none", a topic only contributes if the PEAK/GENE term(s)
  # are marked in_topic for that topic in topic_terms (e.g., peak_and_gene_in_set).
  .assert_pkg("data.table")
  .assert_pkg("cli")

  gate_mode <- match.arg(gate_mode)
  top_k <- as.integer(top_k)
  min_prob <- as.numeric(min_prob)
  chunk_size <- as.integer(chunk_size)
  n_cores <- as.integer(n_cores)

  if (!is.null(out_file) && file.exists(out_file) && !isTRUE(overwrite)) {
    .log_inform("Skipping link->topic scoring; file exists: {out_file}")
    return(invisible(FALSE))
  }

  dt <- data.table::as.data.table(edges_docs)
  req <- c("doc_id", "peak_id", "gene_key")
  if (include_tf) req <- c(req, "tf")
  .assert_has_cols(dt, req, context = "compute_link_topic_scores")

  theta_mat <- as.matrix(theta)
  phi_mat <- as.matrix(phi)
  if (is.null(colnames(theta_mat))) colnames(theta_mat) <- paste0("Topic", seq_len(ncol(theta_mat)))
  if (is.null(rownames(phi_mat))) rownames(phi_mat) <- colnames(theta_mat)
  if (!all(rownames(phi_mat) == colnames(theta_mat))) {
    phi_mat <- phi_mat[colnames(theta_mat), , drop = FALSE]
  }
  topic_labels <- colnames(theta_mat)
  topic_nums <- .m3_opt_topic_ids(theta_mat, "column")

  doc_ids <- rownames(theta_mat)
  if (is.null(doc_ids) || anyNA(doc_ids) || any(doc_ids == "")) {
    .log_abort("theta must have rownames (doc_id) for link->topic scoring.")
  }
  doc_index <- stats::setNames(seq_along(doc_ids), doc_ids)

  term_ids <- colnames(phi_mat)
  if (is.null(term_ids)) .log_abort("phi must have colnames (term_id).")
  term_index <- stats::setNames(seq_along(term_ids), term_ids)

  dt <- dt[!is.na(doc_id) & nzchar(doc_id)]
  dt[, peak_term := paste0("PEAK:", peak_id)]
  dt[, gene_term := paste0("GENE:", gene_key)]
  if (include_tf) dt[, tf_term := paste0("GENE:", tf)]
  dt[, doc_idx := doc_index[doc_id]]
  dt[, peak_idx := term_index[peak_term]]
  dt[, gene_idx := term_index[gene_term]]
  if (include_tf) dt[, tf_idx := term_index[tf_term]]

  dt <- dt[!is.na(doc_idx) & !is.na(peak_idx) & !is.na(gene_idx)]
  if (!nrow(dt)) {
    .log_inform("No valid links for link->topic scoring after term/doc matching.")
    return(invisible(FALSE))
  }

  inset_map <- NULL
  if (gate_mode != "none") {
    tt <- data.table::as.data.table(topic_terms)
    if (!nrow(tt)) {
      .log_abort("gate_mode requires non-empty topic_terms.")
    }
    in_set_vec <- .as_logical_flag(tt$in_topic)
    tt <- tt[in_set_vec]
    if (!nrow(tt)) {
      .log_abort("gate_mode requires in_topic terms in topic_terms.")
    }
    tt[, topic := as.integer(topic)]
    inset_map <- split(tt$topic, tt$term_id)
  }

  score_chunk <- function(chunk_dt) {
    out <- vector("list", nrow(chunk_dt))
    for (i in seq_len(nrow(chunk_dt))) {
      row <- chunk_dt[i]
      theta_vec <- theta_mat[row$doc_idx, ]
      score <- theta_vec *
        phi_mat[, row$peak_idx] *
        phi_mat[, row$gene_idx]
      if (include_tf && !is.na(row$tf_idx)) {
        score <- score * phi_mat[, row$tf_idx]
      }

      if (gate_mode != "none") {
        mask <- rep(TRUE, length(score))
        if (gate_mode %in% c("peak_in_set", "peak_and_gene_in_set")) {
          allowed <- inset_map[[row$peak_term]]
          mask <- mask & topic_nums %in% allowed
        }
        if (gate_mode %in% c("gene_in_set", "peak_and_gene_in_set")) {
          allowed <- inset_map[[row$gene_term]]
          mask <- mask & topic_nums %in% allowed
        }
        score[!mask] <- 0
      }

      total <- sum(score)
      if (!is.finite(total) || total <= 0) next
      prob <- score / total
      ord <- order(prob, decreasing = TRUE)
      if (is.finite(top_k) && top_k > 0L) ord <- ord[seq_len(min(top_k, length(ord)))]
      if (is.finite(min_prob) && min_prob > 0) {
        ord <- ord[prob[ord] >= min_prob]
      }
      if (!length(ord)) next
      tf_val <- if ("tf" %in% names(row)) row$tf else NA_character_
      out[[i]] <- data.table::data.table(
        doc_id = row$doc_id,
        tf = tf_val,
        peak_id = row$peak_id,
        gene_key = row$gene_key,
        topic = topic_labels[ord],
        topic_num = topic_nums[ord],
        score = score[ord],
        prob = prob[ord]
      )
    }
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  }

  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size <= 0L) {
    chunk_size <- nrow(dt)
  }
  idx <- split(seq_len(nrow(dt)), ceiling(seq_len(nrow(dt)) / chunk_size))

  if (.Platform$OS.type != "windows" && n_cores > 1L && length(idx) > 1L) {
    res <- parallel::mclapply(
      idx,
      function(ii) score_chunk(dt[ii]),
      mc.cores = min(n_cores, length(idx))
    )
  } else {
    res <- lapply(idx, function(ii) score_chunk(dt[ii]))
  }
  out <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  if (is.null(out) || !nrow(out)) {
    .log_inform("No link->topic scores passed filters.")
    return(invisible(FALSE))
  }

  if (!is.null(out_file)) {
    data.table::fwrite(out, out_file)
  }
  out[]
}

.gammafit_cutoffs <- function(score_mat,
                              thrP = 0.975,
                              min_terms = 50L,
                              gammafit_scope = c("topic_term_group", "global_term_group")) {
  # Combined cutoff using the maximum class-specific cutoff.
  gammafit_scope <- match.arg(gammafit_scope)
  cut_tbl <- .gammafit_cutoffs_by_termclass(
    score_mat,
    thrP = thrP,
    min_terms = min_terms,
    gammafit_scope = gammafit_scope
  )
  if (!nrow(cut_tbl)) return(numeric(0))
  apply(as.matrix(cut_tbl[, .(peaks_gamma_cutoff, gene_gamma_cutoff, other_gamma_cutoff)]), 1, function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else max(x)
  })
}

.resolve_fp_term_mode <- function(x) {
  if (exists(".topic_term_mode", mode = "function")) {
    return(.topic_term_mode(x))
  }
  x <- as.character(x)[1]
  x <- tolower(gsub("-", "_", x))
  if (x %in% c("fp_uniq", "fp_unique", "uniq", "unique")) return("unique")
  if (x %in% c("fp_aggr", "fp_agg", "agg", "aggregate")) return("aggregate")
  if (x %in% c("fp_aggr_weight", "fp_agg_weight", "aggregate_weight")) return("aggregate_weight")
  if (x %in% c("gene_expression", "expression", "gene_only_expression")) return("gene_expression")
  .log_abort("Unsupported fp_term_mode: {x}")
}

.topic_links_to_link_scores <- function(topic_links,
                                        method = c("peak_and_gene", "peak_and_gene_prob", "gene_only", "link_score_prob", "link_score_efdr"),
                                        min_prob = 0) {
  .assert_pkg("data.table")
  method <- match.arg(method)
  dt <- data.table::as.data.table(topic_links)
  if (!nrow(dt)) return(data.table::data.table())
  req <- c("doc_id", "peak_id", "gene_key", "topic_num", "peak_score", "gene_score")
  if (method == "peak_and_gene") {
    req <- c(req, "peak_pass", "gene_pass")
  } else if (method == "peak_and_gene_prob") {
    req <- c(req, "peak_pass", "gene_prob_pass")
  } else if (method == "gene_only") {
    req <- c(req, "gene_pass")
  } else {
    req <- c(req, "link_pass")
  }
  .assert_has_cols(dt, req, context = "topic_links_to_link_scores")
  if (method == "peak_and_gene") {
    dt <- dt[.as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
  } else if (method == "peak_and_gene_prob") {
    dt <- dt[.as_logical_flag(peak_pass) & .as_logical_flag(gene_prob_pass)]
  } else if (method == "gene_only") {
    dt <- dt[.as_logical_flag(gene_pass)]
  } else {
    dt <- dt[.as_logical_flag(link_pass)]
  }
  if (!nrow(dt)) return(data.table::data.table())
  if ("link_score" %in% names(dt)) {
    dt[, score := .safe_num(link_score)]
    dt[!is.finite(score), score := 0]
  } else if (method == "gene_only") {
    dt[, score := .safe_num(gene_score)]
    dt[!is.finite(score), score := 0]
  } else {
    dt[, score := .safe_num(peak_score) * .safe_num(gene_score)]
    dt[!is.finite(score), score := 0]
  }
  if (method == "link_score_prob" && "link_score_prob" %in% names(dt)) {
    dt[, prob := .safe_num(link_score_prob)]
    dt[!is.finite(prob), prob := 0]
  } else if (method == "peak_and_gene_prob" && "gene_prob" %in% names(dt)) {
    dt[, prob := .safe_num(gene_prob)]
    dt[!is.finite(prob), prob := 0]
  } else if (method == "gene_only" && "gene_score" %in% names(dt)) {
    dt[, prob := .safe_num(gene_score)]
    dt[!is.finite(prob), prob := 0]
  } else {
    dt[, prob := 1]
  }
  if (is.finite(min_prob) && min_prob > 1) {
    dt <- dt[FALSE]
  }
  keep <- c("doc_id", "tf", "peak_id", "gene_key", "topic_num", "score", "prob")
  if (!"tf" %in% names(dt)) dt[, tf := NA_character_]
  dt[, ..keep]
}

.topic_links_path <- function(out_dir, prefer = c("pass", "full")) {
  prefer <- match.arg(prefer)
  pass_path <- file.path(out_dir, "topic_links_pass.csv")
  pass_gz_path <- file.path(out_dir, "topic_links_pass.csv.gz")
  full_path <- file.path(out_dir, "topic_links.csv")
  if (identical(prefer, "pass") && file.exists(pass_path)) return(pass_path)
  if (identical(prefer, "pass") && file.exists(pass_gz_path)) return(pass_gz_path)
  if (file.exists(full_path)) return(full_path)
  if (file.exists(pass_path)) return(pass_path)
  if (file.exists(pass_gz_path)) return(pass_gz_path)
  full_path
}

.empirical_link_fdr <- function(peak_score,
                                gene_score,
                                B = 100L,
                                seed = 1L) {
  peak_score <- .safe_num(peak_score)
  gene_score <- .safe_num(gene_score)
  ok <- is.finite(peak_score) & is.finite(gene_score)
  n <- sum(ok)
  p_emp <- rep(NA_real_, length(peak_score))
  q_emp <- rep(NA_real_, length(peak_score))
  if (n <= 0L) {
    return(list(p_emp = p_emp, q_emp = q_emp))
  }
  if (n == 1L || !is.finite(B) || B < 1L) {
    p_emp[ok] <- 1
    q_emp[ok] <- 1
    return(list(p_emp = p_emp, q_emp = q_emp))
  }

  B <- as.integer(B)
  ps <- peak_score[ok]
  gs <- gene_score[ok]
  obs <- ps * gs
  ge_count <- numeric(n)

  if (is.finite(seed)) {
    seed <- as.integer(seed)
    withr::with_seed(seed, {
      for (b in seq_len(B)) {
        perm_gs <- sample(gs, length(gs), replace = FALSE)
        null_scores <- ps * perm_gs
        null_scores <- sort(null_scores)
        n_lt <- findInterval(obs, null_scores, left.open = TRUE, rightmost.closed = TRUE)
        ge_count <- ge_count + (n - n_lt)
      }
    })
  } else {
    for (b in seq_len(B)) {
      perm_gs <- sample(gs, length(gs), replace = FALSE)
      null_scores <- ps * perm_gs
      null_scores <- sort(null_scores)
      n_lt <- findInterval(obs, null_scores, left.open = TRUE, rightmost.closed = TRUE)
      ge_count <- ge_count + (n - n_lt)
    }
  }

  p <- (1 + ge_count) / (1 + as.double(B) * n)
  q <- stats::p.adjust(p, method = "BH")
  p_emp[ok] <- p
  q_emp[ok] <- q
  list(p_emp = p_emp, q_emp = q_emp)
}

.topic_item_coverage_row <- function(unit, total, pass, count_basis) {
  total <- as.numeric(total)
  pass <- as.numeric(pass)
  if (!is.finite(total) || total < 0) total <- 0
  if (!is.finite(pass) || pass < 0) pass <- 0
  pass <- min(pass, total)
  fraction <- if (total > 0) c(pass / total, (total - pass) / total) else c(NA_real_, NA_real_)
  data.table::data.table(
    unit = unit,
    status = c("Pass", "Fail"),
    count = c(pass, total - pass),
    total = total,
    fraction = fraction,
    percent = 100 * fraction,
    count_basis = count_basis
  )
}

.topic_item_coverage_from_terms <- function(topic_terms, score_mat = NULL) {
  total_terms <- character()
  pass_terms <- character()
  if (!is.null(score_mat) && !is.null(colnames(score_mat))) {
    total_terms <- unique(as.character(colnames(score_mat)))
  }
  if (!is.null(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (nrow(tt) && "term_id" %in% names(tt)) {
      total_terms <- union(total_terms, unique(as.character(tt$term_id)))
      if ("in_topic" %in% names(tt)) {
        pass_terms <- unique(as.character(tt[.as_logical_flag(in_topic), term_id]))
      }
    }
  }
  total_terms <- total_terms[!is.na(total_terms) & nzchar(total_terms)]
  pass_terms <- pass_terms[!is.na(pass_terms) & nzchar(pass_terms)]
  .topic_item_coverage_row(
    unit = "Terms",
    total = length(total_terms),
    pass = length(pass_terms),
    count_basis = "Model terms assigned to at least one topic"
  )
}

.topic_term_coverage_summary_lines <- function(topic_terms, score_mat = NULL) {
  total_terms <- character()
  pass_terms <- character()
  if (!is.null(score_mat) && !is.null(colnames(score_mat))) {
    total_terms <- unique(as.character(colnames(score_mat)))
  }
  if (!is.null(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (nrow(tt) && "term_id" %in% names(tt)) {
      total_terms <- union(total_terms, unique(as.character(tt$term_id)))
      if ("in_topic" %in% names(tt)) {
        pass_terms <- unique(as.character(tt[.as_logical_flag(in_topic), term_id]))
      }
    }
  }
  total_terms <- total_terms[!is.na(total_terms) & nzchar(total_terms)]
  pass_terms <- pass_terms[!is.na(pass_terms) & nzchar(pass_terms)]

  .fmt_line <- function(label, total, pass) {
    pct <- if (total > 0) 100 * pass / total else NA_real_
    pct_txt <- if (is.finite(pct)) sprintf("%.2f", pct) else "NA"
    sprintf("%s: %d / %d = %s%%", label, as.integer(pass), as.integer(total), pct_txt)
  }

  total_gene <- total_terms[grepl("^GENE:", total_terms)]
  pass_gene <- intersect(pass_terms, total_gene)
  total_peak <- total_terms[grepl("^PEAK:", total_terms)]
  pass_peak <- intersect(pass_terms, total_peak)
  c(
    .fmt_line("Terms", length(total_terms), length(intersect(pass_terms, total_terms))),
    .fmt_line("GENE terms", length(total_gene), length(pass_gene)),
    .fmt_line("PEAK terms", length(total_peak), length(pass_peak))
  )
}

.topic_assignment_coverage_summary_table <- function(topic_terms,
                                                     score_mat = NULL,
                                                     item_coverage = NULL,
                                                     show_peak_expanded_link_coverage = TRUE) {
  total_terms <- character()
  pass_terms <- character()
  if (!is.null(score_mat) && !is.null(colnames(score_mat))) {
    total_terms <- unique(as.character(colnames(score_mat)))
  }
  if (!is.null(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (nrow(tt) && "term_id" %in% names(tt)) {
      total_terms <- union(total_terms, unique(as.character(tt$term_id)))
      if ("in_topic" %in% names(tt)) {
        pass_terms <- unique(as.character(tt[.as_logical_flag(in_topic), term_id]))
      }
    }
  }
  total_terms <- total_terms[!is.na(total_terms) & nzchar(total_terms)]
  pass_terms <- pass_terms[!is.na(pass_terms) & nzchar(pass_terms)]

  make_row <- function(label, pass, total) {
    pass <- as.numeric(pass)
    total <- as.numeric(total)
    if (!is.finite(pass) || pass < 0) pass <- 0
    if (!is.finite(total) || total < 0) total <- 0
    pass <- min(pass, total)
    percent <- if (total > 0) 100 * pass / total else NA_real_
    data.table::data.table(
      label = label,
      pass = pass,
      total = total,
      percent = percent,
      label_text = sprintf(
        "%d / %d = %s%%",
        as.integer(pass),
        as.integer(total),
        if (is.finite(percent)) sprintf("%.2f", percent) else "NA"
      )
    )
  }

  total_gene_terms <- total_terms[grepl("^GENE:", total_terms)]
  total_peak_terms <- total_terms[grepl("^PEAK:", total_terms)]
  out <- list(
    make_row("Terms", length(intersect(pass_terms, total_terms)), length(total_terms)),
    make_row("GENE terms", length(intersect(pass_terms, total_gene_terms)), length(total_gene_terms)),
    make_row("PEAK terms", length(intersect(pass_terms, total_peak_terms)), length(total_peak_terms))
  )
  if (!is.null(item_coverage)) {
    cov <- data.table::as.data.table(item_coverage)
    req <- c("unit", "status", "count", "total")
    if (nrow(cov) && all(req %in% names(cov))) {
      cov <- cov[status == "Pass" & unit %in% c("Genes", "TF-gene-doc links", "Links", "TFs")]
      if ("count_basis" %in% names(cov)) {
        cov <- cov[unit != "TFs" | grepl("raw theta", count_basis, ignore.case = TRUE)]
      }
      if (nrow(cov)) {
        cov[, unit_order__ := match(unit, c("Genes", "TF-gene-doc links", "Links", "TFs"))]
        data.table::setorder(cov, unit_order__)
        cov[, unit_order__ := NULL]
        for (i in seq_len(nrow(cov))) {
          if (identical(cov$unit[[i]], "Links") && !isTRUE(show_peak_expanded_link_coverage)) {
            next
          }
          label_i <- if (identical(cov$unit[[i]], "Links")) "TF-peak-gene links" else cov$unit[[i]]
          out[[length(out) + 1L]] <- make_row(label_i, cov$count[[i]], cov$total[[i]])
        }
      }
    }
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

.topic_item_coverage_from_links <- function(scored_links, pass_links) {
  scored_links <- data.table::as.data.table(scored_links)
  pass_links <- data.table::as.data.table(pass_links)
  out <- list()
  if ("gene_key" %in% names(scored_links)) {
    total_genes <- unique(as.character(scored_links[!is.na(gene_key) & nzchar(gene_key), gene_key]))
    pass_genes <- if ("gene_key" %in% names(pass_links)) {
      unique(as.character(pass_links[!is.na(gene_key) & nzchar(gene_key), gene_key]))
    } else {
      character()
    }
    out[[length(out) + 1L]] <- .topic_item_coverage_row(
      unit = "Genes",
      total = length(total_genes),
      pass = length(pass_genes),
      count_basis = "Genes assigned to at least one topic"
    )
  }
  link_key_cols <- intersect(c("doc_id", "tf", "peak_id", "gene_key"), names(scored_links))
  if (length(link_key_cols) >= 3L) {
    make_link_key <- function(x) {
      x <- data.table::as.data.table(x)
      if (!nrow(x) || !all(link_key_cols %in% names(x))) return(character())
      key_dt <- unique(x[, ..link_key_cols])
      for (cc in link_key_cols) {
        key_dt[[cc]] <- as.character(key_dt[[cc]])
        key_dt[[cc]][is.na(key_dt[[cc]])] <- ""
      }
      key <- do.call(paste, c(key_dt[, ..link_key_cols], sep = "\r"))
      unique(key[nzchar(key)])
    }
    total_links <- make_link_key(scored_links)
    pass_links_key <- make_link_key(pass_links)
    out[[length(out) + 1L]] <- .topic_item_coverage_row(
      unit = "Links",
      total = length(total_links),
      pass = length(intersect(pass_links_key, total_links)),
      count_basis = "TF-peak-gene links assigned to at least one topic"
    )
  }
  tf_gene_doc_cols <- intersect(c("doc_id", "tf", "gene_key"), names(scored_links))
  if (length(tf_gene_doc_cols) >= 2L) {
    make_tf_gene_doc_key <- function(x) {
      x <- data.table::as.data.table(x)
      if (!nrow(x) || !all(tf_gene_doc_cols %in% names(x))) return(character())
      key_dt <- unique(x[, ..tf_gene_doc_cols])
      for (cc in tf_gene_doc_cols) {
        key_dt[[cc]] <- as.character(key_dt[[cc]])
        key_dt[[cc]][is.na(key_dt[[cc]])] <- ""
      }
      key <- do.call(paste, c(key_dt[, ..tf_gene_doc_cols], sep = "\r"))
      unique(key[nzchar(key)])
    }
    total_tf_gene_doc <- make_tf_gene_doc_key(scored_links)
    pass_tf_gene_doc <- make_tf_gene_doc_key(pass_links)
    out[[length(out) + 1L]] <- .topic_item_coverage_row(
      unit = "TF-gene-doc links",
      total = length(total_tf_gene_doc),
      pass = length(intersect(pass_tf_gene_doc, total_tf_gene_doc)),
      count_basis = "TF-gene-document links assigned to at least one topic"
    )
  }
  if ("tf" %in% names(scored_links)) {
    total_tfs <- unique(as.character(scored_links[!is.na(tf) & nzchar(tf), tf]))
    pass_tfs <- if ("tf" %in% names(pass_links)) {
      unique(as.character(pass_links[!is.na(tf) & nzchar(tf), tf]))
    } else {
      character()
    }
    out[[length(out) + 1L]] <- .topic_item_coverage_row(
      unit = "TFs",
      total = length(total_tfs),
      pass = length(pass_tfs),
      count_basis = "TFs represented in at least one passing topic link"
    )
  }
  if (!length(out)) return(data.table::data.table())
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

.topic_item_coverage_from_tf_assignment <- function(assign) {
  if (is.null(assign) || !is.list(assign)) return(data.table::data.table())
  primary <- data.table::as.data.table(assign$primary %||% data.table::data.table())
  pass <- data.table::as.data.table(assign$pass %||% data.table::data.table())
  if (!nrow(primary) || !"tf" %in% names(primary)) {
    return(data.table::data.table())
  }
  total_tfs <- unique(as.character(primary[!is.na(tf) & nzchar(tf), tf]))
  pass_tfs <- if (nrow(pass) && "tf" %in% names(pass)) {
    unique(as.character(pass[!is.na(tf) & nzchar(tf), tf]))
  } else {
    character()
  }
  .topic_item_coverage_row(
    unit = "TFs",
    total = length(total_tfs),
    pass = length(pass_tfs),
    count_basis = "TFs assigned to at least one topic from raw theta documents"
  )
}

.replace_tf_item_coverage_with_assignment <- function(out_dir, assign) {
  .assert_pkg("data.table")
  tf_coverage <- .topic_item_coverage_from_tf_assignment(assign)
  if (!nrow(tf_coverage)) return(invisible(FALSE))
  coverage_path <- file.path(out_dir, "topic_item_coverage_counts.csv")
  if (file.exists(coverage_path)) {
    coverage <- data.table::fread(coverage_path, showProgress = FALSE)
    coverage <- coverage[unit != "TFs"]
  } else {
    coverage <- data.table::data.table()
  }
  coverage <- data.table::rbindlist(
    list(coverage, tf_coverage),
    use.names = TRUE,
    fill = TRUE
  )
  unit_levels <- c("Terms", "Genes", "TF-gene-doc links", "Links", "TFs")
  coverage[, unit_order__ := match(unit, unit_levels)]
  coverage[is.na(unit_order__), unit_order__ := length(unit_levels) + seq_len(sum(is.na(unit_order__)))]
  status_levels <- c("Pass", "Fail")
  coverage[, status_order__ := match(status, status_levels)]
  coverage[is.na(status_order__), status_order__ := length(status_levels) + seq_len(sum(is.na(status_order__)))]
  data.table::setorder(coverage, unit_order__, status_order__)
  coverage[, c("unit_order__", "status_order__") := NULL]
  data.table::fwrite(coverage, coverage_path)
  invisible(TRUE)
}

.topic_item_coverage_from_scored_objects <- function(scored_links,
                                                     pass_links,
                                                     topic_terms = NULL,
                                                     score_mat = NULL) {
  data.table::rbindlist(
    list(
      .topic_item_coverage_from_terms(topic_terms, score_mat),
      .topic_item_coverage_from_links(scored_links, pass_links)
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

compute_topic_links <- function(edges_docs,
                                score_mat,
                                raw_score_mat = NULL,
                                topic_terms = NULL,
                                use_final_term_assignment = FALSE,
                                theta = NULL,
                                topic_tf_membership_cutoff = 0.3,
                                fp_term_mode = c("unique", "aggregate", "aggregate_weight", "gene_expression"),
                                binarize_method = c("gammafit", "topn"),
                                gammafit_scope = c("topic_term_group", "global_term_group"),
                                link_method = c("gammafit", "theta_and_terms", "link_score_efdr", "link_score_prob", "gene_prob"),
                                link_prob_cutoff = 0.3,
                                thrP = 0.975,
                                min_terms = 50L,
                                fdr_q = 0.2,
                                fdr_p = NA_real_,
                                efdr_scope = c("per_topic", "global"),
                                efdr_B = 100L,
                                efdr_seed = 1L,
                                out_file = NULL,
                                output_mode = c("pass", "full", "both", "none"),
                                pass_file = NULL,
                                chunk_size = 5000L,
                                n_cores = 1L,
                                overwrite = FALSE) {
  .assert_pkg("data.table")
  .assert_pkg("cli")
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  binarize_method <- match.arg(binarize_method)
  gammafit_scope <- match.arg(gammafit_scope)
  link_method <- match.arg(link_method)
  output_mode <- match.arg(output_mode)
  efdr_scope <- match.arg(efdr_scope)
  chunk_size <- as.integer(chunk_size)
  n_cores <- as.integer(n_cores)
  efdr_B <- as.integer(efdr_B)
  link_prob_cutoff_raw <- link_prob_cutoff
  fdr_q <- .safe_num(fdr_q)
  fdr_p <- .safe_num(fdr_p)
  use_prob_max <- is.character(link_prob_cutoff_raw) && length(link_prob_cutoff_raw) &&
    identical(tolower(trimws(link_prob_cutoff_raw[1])), "max")
  if (identical(link_method, "link_score_prob") || identical(link_method, "gene_prob")) {
    if (!use_prob_max) {
      link_prob_cutoff <- .safe_num(link_prob_cutoff_raw)
      if (!is.finite(link_prob_cutoff) || link_prob_cutoff <= 0 || link_prob_cutoff >= 1) {
        .log_abort("link_prob_cutoff must be in (0,1) or 'max' for {link_method}.")
      }
    } else {
      link_prob_cutoff <- "max"
    }
  } else if (!identical(link_method, "gammafit") && !identical(link_method, "theta_and_terms")) {
    use_p_cut <- is.finite(fdr_p) && fdr_p > 0 && fdr_p < 1
    use_q_cut <- is.finite(fdr_q) && fdr_q > 0 && fdr_q < 1
    if (!use_p_cut && !use_q_cut) {
      .log_abort("Provide at least one valid cutoff: fdr_p in (0,1) or fdr_q in (0,1).")
    }
  }

  if (!is.null(out_file) && output_mode %in% c("full", "both") && file.exists(out_file) && !isTRUE(overwrite)) {
    .log_inform("Skipping topic_links; file exists: {out_file}")
    return(invisible(FALSE))
  }
  if (!is.null(pass_file) && output_mode %in% c("pass", "both") && file.exists(pass_file) && !isTRUE(overwrite)) {
    .log_inform("Skipping compact topic_links; file exists: {pass_file}")
    return(invisible(FALSE))
  }

  dt <- data.table::as.data.table(edges_docs)
  req <- c("doc_id", "peak_id", "gene_key")
  .assert_has_cols(dt, req, context = "compute_topic_links")
  if (!"tf" %in% names(dt)) dt[, tf := NA_character_]

  score_mat <- as.matrix(score_mat)
  K <- nrow(score_mat)
  if (!K) {
    .log_abort("score_mat has zero topics.")
  }
  # Raw matrix is used for exported peak_score/gene_score; normalized score_mat is
  # used for gammafit cutoffs and peak_pass/gene_pass decisions.
  if (is.null(raw_score_mat)) {
    raw_score_mat <- score_mat
  } else {
    raw_score_mat <- as.matrix(raw_score_mat)
    raw_score_mat[!is.finite(raw_score_mat) | raw_score_mat < 0] <- 0
    if (nrow(raw_score_mat) == K) {
      if (!is.null(rownames(raw_score_mat)) && !is.null(rownames(score_mat))) {
        rr <- match(rownames(score_mat), rownames(raw_score_mat))
        if (all(is.finite(rr))) raw_score_mat <- raw_score_mat[rr, , drop = FALSE]
      }
    } else {
      .log_warn("raw_score_mat row count does not match score_mat; falling back to score_mat.")
      raw_score_mat <- score_mat
    }
  }
  term_ids <- colnames(score_mat)
  if (is.null(term_ids)) .log_abort("score_mat must have colnames (term_id).")
  term_index <- stats::setNames(seq_along(term_ids), term_ids)
  raw_terms <- colnames(raw_score_mat)
  if (is.null(raw_terms) || !all(term_ids %in% raw_terms)) {
    .log_warn("raw_score_mat colnames do not fully match score_mat; falling back to score_mat.")
    raw_score_mat <- score_mat
  } else {
    raw_score_mat <- raw_score_mat[, term_ids, drop = FALSE]
  }
  theta_gate <- NULL
  theta_cutoff <- .safe_num(topic_tf_membership_cutoff)
  if (identical(link_method, "theta_and_terms")) {
    if (is.null(theta)) {
      .log_abort("theta is required when link_method = 'theta_and_terms'.")
    }
    theta_gate <- as.matrix(theta)
    if (is.null(rownames(theta_gate))) {
      .log_abort("theta must have document row names when link_method = 'theta_and_terms'.")
    }
    if (ncol(theta_gate) != K) {
      .log_abort("theta must have the same number of topic columns as score_mat rows.")
    }
    if (!is.null(colnames(theta_gate)) && !is.null(rownames(score_mat))) {
      theta_match <- match(rownames(score_mat), colnames(theta_gate))
      if (all(is.finite(theta_match))) {
        theta_gate <- theta_gate[, theta_match, drop = FALSE]
      }
    }
    theta_gate[!is.finite(theta_gate) | theta_gate < 0] <- 0
    if (!is.finite(theta_cutoff) || theta_cutoff < 0 || theta_cutoff > 1) {
      .log_abort("topic_tf_membership_cutoff must be between 0 and 1 for theta_and_terms.")
    }
  }

  dt <- dt[!is.na(doc_id) & nzchar(doc_id)]
  if (identical(fp_term_mode, "aggregate")) {
    dt[, peak_term := paste0("PEAK:", gene_key)]
  } else {
    dt[, peak_term := paste0("PEAK:", peak_id)]
  }
  dt[, gene_term := paste0("GENE:", gene_key)]
  gene_only_terms <- fp_term_mode %in% c("aggregate_weight", "gene_expression")
  if (gene_only_terms) {
    dt[, peak_idx := NA_integer_]
  } else {
    dt[, peak_idx := term_index[peak_term]]
  }
  dt[, gene_idx := term_index[gene_term]]
  if (gene_only_terms) {
    dt <- dt[!is.na(gene_idx)]
  } else {
    dt <- dt[!is.na(peak_idx) & !is.na(gene_idx)]
  }
  if (!nrow(dt)) {
    .log_inform("No valid links for topic_links after term matching.")
    return(invisible(FALSE))
  }
  n_candidate_rows <- as.double(nrow(dt)) * as.double(K)
  max_rows_env <- Sys.getenv("CRAFTGRN_TOPIC_LINK_MAX_SCORED_ROWS", unset = "")
  max_scored_rows <- if (nzchar(max_rows_env)) {
    suppressWarnings(as.numeric(max_rows_env[[1L]]))
  } else {
    suppressWarnings(as.numeric(getOption("craftgrn.topic_link.max_scored_rows", 1e8)))
  }
  if (is.na(max_scored_rows) || max_scored_rows <= 0) max_scored_rows <- 1e8
  .log_inform(
    "Topic-link preflight: {format(nrow(dt), big.mark = ',')} valid edges x {K} topics = {format(n_candidate_rows, scientific = FALSE, big.mark = ',')} candidate rows."
  )
  if (is.finite(max_scored_rows) && n_candidate_rows > max_scored_rows) {
    .log_abort(c(
      "Topic-link scoring would materialize too many edge-topic rows.",
      i = "Candidate rows: {format(n_candidate_rows, scientific = FALSE, big.mark = ',')}; safety limit: {format(max_scored_rows, scientific = FALSE, big.mark = ',')}.",
      i = "Standard extraction does not require topic-link enumeration.",
      i = "For an intentional large diagnostic, set CRAFTGRN_TOPIC_LINK_MAX_SCORED_ROWS or options(craftgrn.topic_link.max_scored_rows = Inf)."
    ))
  }

  gamma_cutoffs_peak <- rep(NA_real_, K)
  gamma_cutoffs_gene <- rep(NA_real_, K)
  in_set_map <- NULL
  final_topic_map <- NULL
  if (isTRUE(use_final_term_assignment)) {
    if (is.null(topic_terms) || !is.data.frame(topic_terms) || !nrow(topic_terms)) {
      .log_abort("topic_terms is required when use_final_term_assignment = TRUE.")
    }
    tt <- data.table::as.data.table(topic_terms)
    .assert_has_cols(tt, c("term_id", "in_topic"), context = "compute_topic_links topic_terms")
    if (!"topic_num" %in% names(tt)) {
      if (!"topic" %in% names(tt)) {
        .log_abort("topic_terms must contain topic_num or topic.")
      }
      tt[, topic_num := suppressWarnings(as.integer(gsub("^Topic", "", as.character(topic))))]
    }
    tt <- tt[
      .as_logical_flag(in_topic) & is.finite(topic_num) &
        grepl("^(GENE|PEAK):", term_id)
    ]
    final_assignment <- unique(tt[, .(term_id, topic_num = as.integer(topic_num))])
    duplicate_terms <- final_assignment[, .N, by = term_id][N > 1L]
    if (nrow(duplicate_terms)) {
      .log_abort(
        "Final topic-term assignment contains {nrow(duplicate_terms)} term(s) assigned to multiple topics."
      )
    }
    final_topic_map <- stats::setNames(
      final_assignment$topic_num,
      final_assignment$term_id
    )
  } else if (binarize_method == "gammafit") {
    cut_tbl <- .gammafit_cutoffs_by_termclass(
      score_mat,
      thrP = thrP,
      min_terms = min_terms,
      gammafit_scope = gammafit_scope
    )
    gamma_cutoffs_peak <- cut_tbl$peaks_gamma_cutoff
    gamma_cutoffs_gene <- cut_tbl$gene_gamma_cutoff
  } else if (!is.null(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (nrow(tt)) {
      in_set_vec <- .as_logical_flag(tt$in_topic)
      tt <- tt[in_set_vec]
      tt[, topic := as.integer(topic)]
      in_set_map <- split(tt$topic, tt$term_id)
    }
  }

  score_chunk <- function(chunk_dt) {
    if (!nrow(chunk_dt)) return(data.table::data.table())
    peak_idx <- chunk_dt$peak_idx
    gene_idx <- chunk_dt$gene_idx
    has_peak_terms <- !gene_only_terms
    if (has_peak_terms) {
      peak_scores_gate <- score_mat[, peak_idx, drop = FALSE]
      peak_scores_raw <- raw_score_mat[, peak_idx, drop = FALSE]
    } else {
      peak_scores_gate <- matrix(NA_real_, nrow = K, ncol = nrow(chunk_dt))
      peak_scores_raw <- peak_scores_gate
    }
    gene_scores_gate <- score_mat[, gene_idx, drop = FALSE]
    gene_scores_raw <- raw_score_mat[, gene_idx, drop = FALSE]
    n <- nrow(chunk_dt)
    rep_doc <- rep(chunk_dt$doc_id, each = K)
    rep_tf <- rep(chunk_dt$tf, each = K)
    rep_peak <- rep(chunk_dt$peak_id, each = K)
    rep_gene <- rep(chunk_dt$gene_key, each = K)
    topic_ids <- .m3_opt_topic_ids(score_mat, "row")
    rep_topic <- rep(topic_ids, times = n)
    peak_score <- as.vector(peak_scores_raw)
    gene_score <- as.vector(gene_scores_raw)
    peak_score_gate <- as.vector(peak_scores_gate)
    gene_score_gate <- as.vector(gene_scores_gate)
    peaks_gamma_cutoff <- rep(gamma_cutoffs_peak, times = n)
    gene_gamma_cutoff <- rep(gamma_cutoffs_gene, times = n)
    theta_value <- rep(NA_real_, length(rep_topic))
    theta_pass <- rep(TRUE, length(rep_topic))
    if (identical(link_method, "theta_and_terms")) {
      theta_rows <- match(chunk_dt$doc_id, rownames(theta_gate))
      theta_chunk <- matrix(NA_real_, nrow = n, ncol = K)
      ok_theta <- !is.na(theta_rows)
      if (any(ok_theta)) {
        theta_chunk[ok_theta, ] <- theta_gate[theta_rows[ok_theta], , drop = FALSE]
      }
      theta_value <- as.vector(t(theta_chunk))
      theta_pass <- is.finite(theta_value) & theta_value >= theta_cutoff
    }
    if (isTRUE(use_final_term_assignment)) {
      gene_topic <- unname(final_topic_map[chunk_dt$gene_term])
      gene_pass <- rep_topic == rep(gene_topic, each = K)
      gene_pass[is.na(gene_pass)] <- FALSE
      if (has_peak_terms) {
        peak_topic <- unname(final_topic_map[chunk_dt$peak_term])
        peak_pass <- rep_topic == rep(peak_topic, each = K)
        peak_pass[is.na(peak_pass)] <- FALSE
      } else {
        peak_pass <- rep(TRUE, length(rep_topic))
      }
    } else if (binarize_method == "gammafit") {
      if (has_peak_terms) {
        peak_pass <- is.finite(peaks_gamma_cutoff) & peak_score_gate >= peaks_gamma_cutoff & peak_score_gate > 0
      } else {
        peak_pass <- rep(TRUE, length(rep_topic))
      }
      gene_pass <- is.finite(gene_gamma_cutoff) & gene_score_gate >= gene_gamma_cutoff & gene_score_gate > 0
    } else if (!is.null(in_set_map)) {
      peak_allowed <- in_set_map[chunk_dt$peak_term]
      gene_allowed <- in_set_map[chunk_dt$gene_term]
      peak_pass <- rep(gene_only_terms, length(rep_topic))
      gene_pass <- rep(FALSE, length(rep_topic))
      if (has_peak_terms) {
        for (i in seq_len(n)) {
          if (length(peak_allowed[[i]])) {
            idx <- ((i - 1L) * K + 1L):(i * K)
            peak_pass[idx] <- rep_topic[idx] %in% peak_allowed[[i]]
          }
        }
      }
      for (i in seq_len(n)) {
        if (length(gene_allowed[[i]])) {
          idx <- ((i - 1L) * K + 1L):(i * K)
          gene_pass[idx] <- rep_topic[idx] %in% gene_allowed[[i]]
        }
      }
    } else {
      peak_pass <- rep(FALSE, length(rep_topic))
      gene_pass <- rep(FALSE, length(rep_topic))
    }
    if (identical(link_method, "gammafit")) {
      keep <- peak_pass & gene_pass
      if (!any(keep)) return(data.table::data.table())
      return(data.table::data.table(
        doc_id = rep_doc[keep],
        tf = rep_tf[keep],
        peak_id = rep_peak[keep],
        gene_key = rep_gene[keep],
        topic_num = rep_topic[keep],
        peak_score = peak_score[keep],
        gene_score = gene_score[keep],
        peaks_gamma_cutoff = peaks_gamma_cutoff[keep],
        gene_gamma_cutoff = gene_gamma_cutoff[keep],
        peak_pass = TRUE,
        gene_pass = TRUE,
        link_pass = TRUE
      ))
    }
    if (identical(link_method, "theta_and_terms")) {
      keep <- peak_pass & gene_pass & theta_pass
      if (!any(keep)) return(data.table::data.table())
      return(data.table::data.table(
        doc_id = rep_doc[keep],
        tf = rep_tf[keep],
        peak_id = rep_peak[keep],
        gene_key = rep_gene[keep],
        topic_num = rep_topic[keep],
        peak_score = peak_score[keep],
        gene_score = gene_score[keep],
        peaks_gamma_cutoff = peaks_gamma_cutoff[keep],
        gene_gamma_cutoff = gene_gamma_cutoff[keep],
        peak_pass = peak_pass[keep],
        gene_pass = gene_pass[keep],
        theta = theta_value[keep],
        theta_cutoff = theta_cutoff,
        theta_pass = TRUE,
        link_pass = TRUE
      ))
    }
    data.table::data.table(
      doc_id = rep_doc,
      tf = rep_tf,
      peak_id = rep_peak,
      gene_key = rep_gene,
      topic_num = rep_topic,
      peak_score = peak_score,
      gene_score = gene_score,
      peaks_gamma_cutoff = peaks_gamma_cutoff,
      gene_gamma_cutoff = gene_gamma_cutoff,
      peak_pass = peak_pass,
      gene_pass = gene_pass
    )
  }

  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size <= 0L) {
    chunk_size <- nrow(dt)
  }
  idx <- split(seq_len(nrow(dt)), ceiling(seq_len(nrow(dt)) / chunk_size))

  if (.Platform$OS.type != "windows" && n_cores > 1L && length(idx) > 1L) {
    res <- parallel::mclapply(
      idx,
      function(ii) score_chunk(dt[ii]),
      mc.cores = min(n_cores, length(idx))
    )
  } else {
    res <- lapply(idx, function(ii) score_chunk(dt[ii]))
  }
  out <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  if (!nrow(out)) {
    .log_inform("No topic_links generated.")
    return(invisible(FALSE))
  }

  if (identical(link_method, "gammafit")) {
    n_before_gammafit <- n_candidate_rows
    n_pass_gammafit <- nrow(out)
    gammafit_label <- if (gene_only_terms) "gene_only" else "peak_and_gene"
    .log_inform(
      "topic_links gammafit: {gammafit_label} pass {n_pass_gammafit}/{n_before_gammafit} rows."
    )
  } else if (identical(link_method, "theta_and_terms")) {
    n_before_theta <- n_candidate_rows
    n_pass_theta <- nrow(out)
    theta_label <- if (gene_only_terms) "theta_and_gene" else "theta_peak_and_gene"
    .log_inform(
      "topic_links theta_and_terms: {theta_label} pass {n_pass_theta}/{n_before_theta} rows (theta>={theta_cutoff})."
    )
  } else if (identical(link_method, "link_score_prob")) {
    out[, link_score := peak_score * gene_score]
    grp <- c("doc_id", "tf", "peak_id", "gene_key")
    out[, link_score_prob := {
      s <- sum(link_score, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) rep(1 / .N, .N) else link_score / s
    }, by = grp]
    if (isTRUE(use_prob_max)) {
      out[, link_pass := FALSE]
      out[, rid___ := .I]
      max_idx <- out[, {
        cand <- which(is.finite(link_score_prob) & .as_logical_flag(peak_pass) & .as_logical_flag(gene_pass))
        if (!length(cand)) .(rid___ = integer(0)) else .(rid___ = rid___[cand[which.max(link_score_prob[cand])]])
      }, by = grp]
      if (nrow(max_idx)) {
        out[rid___ %in% max_idx$rid___, link_pass := TRUE]
      }
      out[, rid___ := NULL]
    } else {
      out[, link_pass := is.finite(link_score_prob) & link_score_prob >= link_prob_cutoff &
                       .as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
    }
    n_pass <- sum(out$link_pass, na.rm = TRUE)
    .log_inform(
      "topic_links prob: pass {n_pass}/{nrow(out)} rows (link_score_prob {if (isTRUE(use_prob_max)) '= max' else paste0('>=', link_prob_cutoff)})."
    )
  } else if (identical(link_method, "gene_prob")) {
    gene_topic <- unique(out[, .(doc_id, gene_key, topic_num, gene_score, gene_pass)])
    gene_topic[, gene_score := .safe_num(gene_score)]
    gene_topic[, gene_prob := 0]
    gene_topic[, gene_prob_pass := FALSE]
    grp <- c("doc_id", "gene_key")
    gene_topic[, `:=`(gene_prob = {
      keep <- .as_logical_flag(gene_pass)
      s <- sum(gene_score[keep], na.rm = TRUE)
      out_prob <- rep(0, .N)
      if (is.finite(s) && s > 0) out_prob[keep] <- gene_score[keep] / s
      out_prob
    }), by = grp]
    if (isTRUE(use_prob_max)) {
      gene_topic[, rid___ := .I]
      keep_idx <- gene_topic[, {
        cand <- which(.as_logical_flag(gene_pass) & is.finite(gene_prob))
        if (!length(cand)) .(rid___ = integer(0)) else .(rid___ = rid___[cand[which.max(gene_prob[cand])]])
      }, by = grp]
      if (nrow(keep_idx)) {
        gene_topic[rid___ %in% keep_idx$rid___, gene_prob_pass := TRUE]
      }
      gene_topic[, rid___ := NULL]
    } else {
      gene_topic[, gene_prob_pass := .as_logical_flag(gene_pass) & is.finite(gene_prob) & gene_prob >= link_prob_cutoff]
    }
    out <- merge(
      out,
      gene_topic[, .(doc_id, gene_key, topic_num, gene_prob, gene_prob_pass)],
      by = c("doc_id", "gene_key", "topic_num"),
      all.x = TRUE,
      sort = FALSE
    )
    out[is.na(gene_prob), gene_prob := 0]
    out[is.na(gene_prob_pass), gene_prob_pass := FALSE]
    n_pass <- sum(out$gene_prob_pass & .as_logical_flag(out$peak_pass), na.rm = TRUE)
    .log_inform(
      "topic_links gene_prob: peak_and_gene_prob pass {n_pass}/{nrow(out)} rows (gene_prob {if (isTRUE(use_prob_max)) '= max' else paste0('>=', link_prob_cutoff)})."
    )
  } else {
    out[, link_score := peak_score * gene_score]
    use_p_cut <- is.finite(fdr_p) && fdr_p > 0 && fdr_p < 1
    out[, `:=`(link_efdr_p = NA_real_, link_efdr_q = NA_real_)]
    if (efdr_scope == "per_topic") {
      topic_ids <- sort(unique(out$topic_num))
      for (k in topic_ids) {
        idxk <- which(out$topic_num == k)
        if (!length(idxk)) next
        efdr <- .empirical_link_fdr(
          peak_score = out$peak_score[idxk],
          gene_score = out$gene_score[idxk],
          B = efdr_B,
          seed = as.integer(efdr_seed) + as.integer(k)
        )
        out$link_efdr_p[idxk] <- efdr$p_emp
        out$link_efdr_q[idxk] <- efdr$q_emp
      }
    } else {
      efdr <- .empirical_link_fdr(
        peak_score = out$peak_score,
        gene_score = out$gene_score,
        B = efdr_B,
        seed = as.integer(efdr_seed)
      )
      out[, `:=`(link_efdr_p = efdr$p_emp, link_efdr_q = efdr$q_emp)]
    }

    n_before_efdr <- nrow(out)
    if (use_p_cut) {
      out[, link_pass := is.finite(link_efdr_p) & link_efdr_p <= fdr_p &
                       .as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
    } else {
      out[, link_pass := is.finite(link_efdr_q) & link_efdr_q <= fdr_q &
                       .as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
    }
    n_pass_efdr <- sum(out$link_pass, na.rm = TRUE)
    if (use_p_cut) {
      .log_inform(
        "topic_links eFDR ({efdr_scope}): pass {n_pass_efdr}/{n_before_efdr} rows (p<={fdr_p}, B={efdr_B}, seed={efdr_seed})."
      )
    } else {
      .log_inform(
        "topic_links eFDR ({efdr_scope}): pass {n_pass_efdr}/{n_before_efdr} rows (q<={fdr_q}, B={efdr_B}, seed={efdr_seed})."
      )
    }
  }

  pass_dt <- if (identical(link_method, "gammafit") || identical(link_method, "theta_and_terms")) {
    out[.as_logical_flag(link_pass)]
  } else if (identical(link_method, "gene_prob")) {
    out[.as_logical_flag(peak_pass) & .as_logical_flag(gene_prob_pass)]
  } else {
    out[.as_logical_flag(link_pass)]
  }
  if (is.null(pass_file) && !is.null(out_file) && output_mode %in% c("pass", "both")) {
    root <- sub("[.]csv([.]gz)?$", "", out_file)
    pass_file <- paste0(root, "_pass.csv")
  }
  if (!is.null(out_file) && output_mode %in% c("full", "both")) {
    data.table::fwrite(out, out_file)
  }
  if (!is.null(pass_file) && output_mode %in% c("pass", "both")) {
    data.table::fwrite(pass_dt, pass_file)
  }
  if (!is.null(out_file) || !is.null(pass_file)) {
    summary_dir <- dirname(if (!is.null(pass_file)) pass_file else out_file)
    summary_file <- file.path(summary_dir, "topic_link_summary.csv")
    n_scored_rows <- if (identical(link_method, "gammafit") || identical(link_method, "theta_and_terms")) {
      n_candidate_rows
    } else {
      as.double(nrow(out))
    }
    summary_dt <- data.table::data.table(
      link_method = link_method,
      output_mode = output_mode,
      n_scored_rows = n_scored_rows,
      n_pass_rows = as.double(nrow(pass_dt)),
      full_file = if (output_mode %in% c("full", "both") && !is.null(out_file)) basename(out_file) else NA_character_,
      pass_file = if (output_mode %in% c("pass", "both") && !is.null(pass_file)) basename(pass_file) else NA_character_
    )
    data.table::fwrite(summary_dt, summary_file)
    coverage_dt <- .topic_item_coverage_from_scored_objects(
      scored_links = dt,
      pass_links = pass_dt,
      topic_terms = topic_terms,
      score_mat = score_mat
    )
    data.table::fwrite(coverage_dt, file.path(summary_dir, "topic_item_coverage_counts.csv"))
  }
  if (identical(link_method, "gammafit") || identical(link_method, "theta_and_terms")) {
    out[.as_logical_flag(link_pass)]
  } else {
    out[]
  }
}

link_scores_to_gene_sets <- function(link_scores,
                                     include_tf = TRUE,
                                     include_gene = TRUE,
                                     min_prob = 0,
                                     gene_min_prob = NULL,
                                     tf_min_prob = NULL,
                                     tf_max_topics = Inf,
                                     tf_top_n_per_topic = NA_integer_,
                                     tf_link_scores = NULL,
                                     gene_terms = NULL) {
  .assert_pkg("data.table")
  to_dt <- function(x) {
    dt <- data.table::as.data.table(x)
    if (!nrow(dt)) return(dt)
    if (!("topic_num" %in% names(dt))) {
      if (!("topic" %in% names(dt))) {
        .log_abort("link_scores must have topic_num or topic.")
      }
      dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
    }
    if (!("prob" %in% names(dt))) dt[, prob := 1]
    dt <- dt[is.finite(topic_num)]
    dt
  }

  dt <- to_dt(link_scores)
  if (!nrow(dt)) return(list())
  if (is.null(gene_min_prob) || !is.finite(gene_min_prob)) gene_min_prob <- min_prob
  if (is.null(tf_min_prob) || !is.finite(tf_min_prob)) tf_min_prob <- min_prob

  if (!nrow(dt)) return(list())

  gene_terms_dt <- NULL
  if (!is.null(gene_terms)) {
    gt <- data.table::as.data.table(gene_terms)
    if (nrow(gt)) {
      if (!("topic_num" %in% names(gt))) {
        if (!("topic" %in% names(gt))) {
          .log_abort("gene_terms must have topic_num or topic.")
        }
        gt[, topic_num := as.integer(gsub("^Topic", "", topic))]
      }
      if ("term_id" %in% names(gt)) {
        gt <- gt[grepl("^GENE:", term_id)]
        gt[, gene := sub("^GENE:", "", term_id)]
      } else if ("gene" %in% names(gt)) {
        gt[, gene := as.character(gene)]
      } else {
        .log_abort("gene_terms must have term_id or gene.")
      }
      if ("in_topic" %in% names(gt)) {
        gt <- gt[isTRUE(in_topic)]
      }
      gt <- gt[is.finite(topic_num)]
      gt <- gt[!is.na(gene) & nzchar(gene)]
      if (nrow(gt)) {
        gene_terms_dt <- unique(gt[, .(topic_num, gene)])
      }
    }
  }

  gene_dt <- list()
  if (isTRUE(include_gene) && "gene_key" %in% names(dt)) {
    gene_tbl <- dt[prob >= gene_min_prob, .(topic_num, gene = gene_key)]
    if (!is.null(gene_terms_dt) && nrow(gene_terms_dt)) {
      gene_tbl <- merge(gene_tbl, gene_terms_dt, by = c("topic_num", "gene"))
    }
    gene_dt[[length(gene_dt) + 1L]] <- gene_tbl
  }
  if (isTRUE(include_tf) && "tf" %in% names(dt)) {
    tf_source <- if (is.null(tf_link_scores)) dt else to_dt(tf_link_scores)
    if (!nrow(tf_source)) {
      tf_dt <- data.table::data.table()
    } else {
      tf_dt <- tf_source[prob >= tf_min_prob, .(topic_num, gene = tf, prob)]
    }
    tf_dt <- tf_dt[!is.na(gene) & nzchar(gene)]
    if (nrow(tf_dt)) {
      if (is.finite(tf_top_n_per_topic) && tf_top_n_per_topic > 0L) {
        data.table::setorder(tf_dt, topic_num, -prob)
        tf_dt <- tf_dt[, head(.SD, as.integer(tf_top_n_per_topic)), by = topic_num]
      }
      if (is.finite(tf_max_topics)) {
        tf_counts <- tf_dt[, .(n_topics = data.table::uniqueN(topic_num)), by = gene]
        tf_keep <- tf_counts[n_topics <= as.integer(tf_max_topics), gene]
        tf_dt <- tf_dt[gene %in% tf_keep]
      }
    }
    if (nrow(tf_dt)) {
      gene_dt[[length(gene_dt) + 1L]] <- tf_dt[, .(topic_num, gene)]
    }
  }
  if (!length(gene_dt)) return(list())
  gene_dt <- data.table::rbindlist(gene_dt, use.names = TRUE, fill = TRUE)
  gene_dt <- gene_dt[!is.na(gene) & nzchar(gene)]
  if (!nrow(gene_dt)) return(list())
  gene_dt <- unique(gene_dt)
  sets <- split(gene_dt$gene, gene_dt$topic_num)
  sets <- lapply(sets, function(x) unique(as.character(x)))
  sets[lengths(sets) > 0]
}

.topic_pathway_db_label <- function(db) {
  db_short <- c(
    GO_Biological_Process_2023 = "GO:BP",
    GO_Cellular_Component_2023 = "GO:CC",
    GO_Molecular_Function_2023 = "GO:MF",
    Reactome_2022 = "Reactome",
    ImmuneSigDB = "ImmuneSigDB",
    WikiPathways_2024_Human = "WikiPathways",
    WikiPathways_2024_Mouse = "WikiPathways",
    MSigDB_Hallmark_2020 = "Hallmark",
    KEGG_2021_Human = "KEGG",
    KEGG_2019_Mouse = "KEGG"
  )
  if (db %in% names(db_short)) db_short[[db]] else db
}

.normalize_topic_pathway_term <- function(pathway, database_label = NULL, pathway_term = NULL) {
  label <- as.character(database_label %||% NA_character_)
  term <- as.character(pathway_term %||% NA_character_)
  if (!length(label) || is.na(label[[1L]]) || !nzchar(label[[1L]])) {
    label <- sub(":.*$", "", as.character(pathway %||% ""))
  }
  if (!length(term) || is.na(term[[1L]]) || !nzchar(term[[1L]])) {
    term <- sub("^[^:]+:\\s*", "", as.character(pathway %||% ""))
  }
  label <- tolower(label[[1L]])
  term <- tolower(term[[1L]])
  term <- gsub("\\b(homo sapiens|mus musculus|human|mouse)\\b", " ", term, perl = TRUE)
  term <- gsub("\\b(hsa|mmu)[0-9]{4,6}\\b", " ", term, perl = TRUE)
  term <- gsub("\\bwp[0-9]+\\b", " ", term, perl = TRUE)
  term <- gsub("[^a-z0-9]+", " ", term, perl = TRUE)
  term <- trimws(gsub("\\s+", " ", term, perl = TRUE))
  label <- trimws(gsub("[^a-z0-9]+", " ", label, perl = TRUE))
  paste(label, term, sep = ":")
}

.select_best_human_mouse_pathways <- function(res_dt) {
  .assert_pkg("data.table")
  res_dt <- data.table::as.data.table(res_dt)
  if (!nrow(res_dt)) return(res_dt)
  if (!"pathway_species" %in% names(res_dt)) res_dt[, pathway_species := NA_character_]
  if (!"database" %in% names(res_dt)) res_dt[, database := sub(":.*$", "", as.character(pathway))]
  if (!"database_label" %in% names(res_dt)) res_dt[, database_label := sub(":.*$", "", as.character(pathway))]
  if (!"pathway_term" %in% names(res_dt)) res_dt[, pathway_term := sub("^[^:]+:\\s*", "", as.character(pathway))]
  if (!"overlap_hits" %in% names(res_dt)) {
    res_dt[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
  }
  if (!"combined_score" %in% names(res_dt)) res_dt[, combined_score := NA_real_]
  if (!"logp" %in% names(res_dt)) res_dt[, logp := -log10(pmax(as.numeric(padj), .Machine$double.xmin))]

  res_dt[, pathway_norm_key := mapply(
    .normalize_topic_pathway_term,
    pathway = pathway,
    database_label = database_label,
    pathway_term = pathway_term,
    USE.NAMES = FALSE
  )]
  res_dt[, selected_pathway_species := as.character(pathway_species)]
  res_dt[, selected_database := as.character(database)]

  species_summary <- res_dt[
    pathway_species %in% c("human", "mouse"),
    {
      ranked <- .SD[order(
        as.numeric(padj),
        -as.numeric(logp),
        -as.numeric(combined_score),
        -as.integer(overlap_hits),
        as.character(database),
        as.character(pathway)
      )]
      ranked[1L, .(
        species_padj = as.numeric(padj),
        species_pval = as.numeric(pval),
        species_logp = as.numeric(logp),
        species_combined_score = as.numeric(combined_score),
        species_overlap_hits = as.integer(overlap_hits),
        species_database = as.character(database),
        species_pathway = as.character(pathway)
      )]
    },
    by = .(topic, pathway_norm_key, pathway_species)
  ]
  if (nrow(species_summary)) {
    species_summary <- data.table::dcast(
      species_summary,
      topic + pathway_norm_key ~ pathway_species,
      value.var = c(
        "species_padj", "species_pval", "species_logp", "species_combined_score",
        "species_overlap_hits", "species_database", "species_pathway"
      )
    )
    data.table::setnames(
      species_summary,
      old = grep("^species_", names(species_summary), value = TRUE),
      new = sub("^species_([^_]+(?:_[^_]+)*)_(human|mouse)$", "\\2_\\1", grep("^species_", names(species_summary), value = TRUE), perl = TRUE)
    )
  }

  ranked <- res_dt[order(
    topic,
    pathway_norm_key,
    as.numeric(padj),
    -as.numeric(logp),
    -as.numeric(combined_score),
    -as.integer(overlap_hits),
    as.character(database),
    as.character(pathway)
  )]
  out <- ranked[, .SD[1L], by = .(topic, pathway_norm_key)]
  if (nrow(species_summary)) {
    out <- merge(out, species_summary, by = c("topic", "pathway_norm_key"), all.x = TRUE, sort = FALSE)
  }
  data.table::setorder(out, topic, padj, pval, -overlap_hits, pathway)
  out[]
}

.topic_enrichr_result_to_table <- function(enr, topic_name, genes, pathway_species = NULL) {
  if (is.null(enr) || !length(enr)) return(data.table::data.table())
  rows <- lapply(names(enr), function(db) {
    df <- enr[[db]]
    if (is.null(df) || !nrow(df)) return(NULL)
    if (!("Adjusted.P.value" %in% names(df)) || !("Term" %in% names(df))) return(NULL)
    df <- df[is.finite(df$Adjusted.P.value), , drop = FALSE]
    if (!nrow(df)) return(NULL)
    db_label <- .topic_pathway_db_label(db)
    term_clean <- gsub("\\s*\\([^)]*\\)$", "", df$Term)
    combined_score <- if ("Combined.Score" %in% names(df)) df$Combined.Score else NA_real_
    odds_ratio <- if ("Odds.Ratio" %in% names(df)) df$Odds.Ratio else NA_real_
    overlap_hits <- if ("Overlap" %in% names(df)) {
      suppressWarnings(as.integer(sub("/.*$", "", as.character(df$Overlap))))
    } else {
      rep(NA_integer_, nrow(df))
    }
    data.table::data.table(
      topic = as.integer(topic_name),
      pathway = paste(db_label, term_clean, sep = ": "),
      database = as.character(db),
      database_label = as.character(db_label),
      pathway_term = as.character(term_clean),
      pathway_species = as.character(pathway_species %||% NA_character_),
      padj = as.numeric(df$Adjusted.P.value),
      pval = if ("P.value" %in% names(df)) as.numeric(df$P.value) else NA_real_,
      overlap = if ("Overlap" %in% names(df)) as.character(df$Overlap) else NA_character_,
      overlap_hits = overlap_hits,
      genes = if ("Genes" %in% names(df)) as.character(df$Genes) else NA_character_,
      logp = -log10(df$Adjusted.P.value),
      combined_score = as.numeric(combined_score),
      odds_ratio = as.numeric(odds_ratio),
      cluster_size = length(genes),
      query_size = if ("query_size" %in% names(df)) as.integer(df$query_size) else length(genes),
      term_size = if ("term_size" %in% names(df)) as.integer(df$term_size) else NA_integer_,
      background_size = if ("background_size" %in% names(df)) as.integer(df$background_size) else NA_integer_
    )
  })
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

plot_topic_pathway_enrichment_gsea <- function(topic_terms,
                                               edges_docs,
                                               option_label,
                                               out_dir,
                                               theta = NULL,
                                               dbs = c(
                                                 "GO_Biological_Process_2023",
                                                 "GO_Cellular_Component_2023",
                                                 "GO_Molecular_Function_2023",
                                                 "Reactome_2022",
                                                 "WikiPathways_2024_Human",
                                                 "MSigDB_Hallmark_2020",
                                                 "KEGG_2021_Human"
                                               ),
                                               species = "Homo sapiens",
                                               padj_cut = 0.05,
                                               min_size = 10L,
                                               max_size = 500L,
                                               nperm = 1000L,
                                               top_n_per_topic = 20L,
                                               max_pathways = 200L,
                                               peak_gene_agg = c("max", "sum"),
                                               tf_link_mode = c("theta", "none"),
                                               tf_theta_top_n = 50L,
                                               tf_theta_min = NA_real_,
                                               title_prefix = NULL) {
  .assert_pkg("data.table")
  log_path <- file.path(out_dir, "topic_pathway_enrichment_gsea_debug.txt")
  log_msg <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
  }

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    msg <- "Skipping GSEA pathway analysis: fgsea not installed."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    msg <- "Skipping GSEA pathway analysis: msigdbr not installed."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    msg <- "Skipping GSEA pathway analysis: pheatmap not installed."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  db_map <- list(
    GO_Biological_Process_2023 = list(category = "C5", subcategory = "GO:BP", label = "GO:BP"),
    GO_Cellular_Component_2023 = list(category = "C5", subcategory = "GO:CC", label = "GO:CC"),
    GO_Molecular_Function_2023 = list(category = "C5", subcategory = "GO:MF", label = "GO:MF"),
    Reactome_2022 = list(category = "C2", subcategory = "CP:REACTOME", label = "Reactome"),
    WikiPathways_2024_Human = list(category = "C2", subcategory = "CP:WIKIPATHWAYS", label = "WikiPathways"),
    MSigDB_Hallmark_2020 = list(category = "H", subcategory = NULL, label = "Hallmark"),
    KEGG_2021_Human = list(category = "C2", subcategory = "CP:KEGG", label = "KEGG")
  )
  dbs_use <- intersect(dbs, names(db_map))
  if (!length(dbs_use)) {
    msg <- "Skipping GSEA pathway analysis: no supported dbs in db_map."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  msig_formals <- names(formals(msigdbr::msigdbr))
  use_collection <- all(c("collection", "subcollection") %in% msig_formals)
  msig_list <- lapply(dbs_use, function(db) {
    cfg <- db_map[[db]]
    msig_args <- if (use_collection) {
      list(species = species, collection = cfg$category)
    } else {
      list(species = species, category = cfg$category)
    }
    if (!is.null(cfg$subcategory) && nzchar(cfg$subcategory)) {
      msig_args[[if (use_collection) "subcollection" else "subcategory"]] <- cfg$subcategory
    }
    msig <- do.call(msigdbr::msigdbr, msig_args)
    if (!nrow(msig)) return(NULL)
    msig <- data.table::as.data.table(msig)
    msig[, pathway := paste(cfg$label, gs_name, sep = ": ")]
    msig[, .(pathway, gene_symbol)]
  })
  msig_dt <- data.table::rbindlist(msig_list, use.names = TRUE, fill = TRUE)
  if (!nrow(msig_dt)) {
    msg <- "Skipping GSEA pathway analysis: msigdbr returned no gene sets."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  pathways <- split(msig_dt$gene_symbol, msig_dt$pathway)

  gene_ranks <- topic_gene_ranks(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    option_label = option_label,
    score_col = "score",
    peak_gene_agg = peak_gene_agg,
    theta = theta,
    tf_link_mode = tf_link_mode,
    tf_theta_top_n = tf_theta_top_n,
    tf_theta_min = tf_theta_min
  )
  if (!length(gene_ranks)) {
    msg <- "Skipping GSEA pathway analysis: no ranked gene lists."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  res_list <- vector("list", length(gene_ranks))
  names(res_list) <- names(gene_ranks)
  log_msg(sprintf("GSEA species: %s", species))
  log_msg(sprintf("GSEA dbs: %s", paste(dbs_use, collapse = ",")))
  for (nm in names(gene_ranks)) {
    stats <- gene_ranks[[nm]]
    if (length(stats) < as.integer(min_size)) {
      log_msg(sprintf("Topic %s skipped: gene rank size < %d", nm, min_size))
      next
    }
    fg <- tryCatch(
      fgsea::fgsea(
        pathways = pathways,
        stats = stats,
        minSize = as.integer(min_size),
        maxSize = as.integer(max_size),
        nperm = as.integer(nperm)
      ),
      error = function(e) {
        log_msg(sprintf("Topic %s fgsea error: %s", nm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fg) || !nrow(fg)) {
      log_msg(sprintf("Topic %s fgsea returned NULL/empty.", nm))
      next
    }
    fg <- data.table::as.data.table(fg)
    fg <- fg[is.finite(padj) & padj <= padj_cut]
    if (!nrow(fg)) next
    fg[, logp := -log10(pmax(padj, .Machine$double.xmin))]
    fg[, topic := as.integer(nm)]
    lead <- if ("leadingEdge" %in% names(fg)) {
      vapply(fg$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
    } else {
      rep(NA_character_, nrow(fg))
    }
    res_list[[nm]] <- fg[, .(topic, pathway = pathway, NES = NES, padj = padj, logp = logp)]
    res_list[[nm]][, leading_edge := lead]
  }

  res_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (!nrow(res_dt)) {
    msg <- "Skipping GSEA pathway analysis: no enriched terms at padj_cut."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) {
    res_dt <- res_dt[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
  } else {
    res_dt <- res_dt[order(topic, -logp)]
  }
  if (nrow(res_dt) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
    path_rank <- res_dt[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
    if (nrow(path_rank) > as.integer(max_pathways)) {
      keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
      res_dt <- res_dt[pathway %in% keep]
      log_msg(sprintf("Filtered pathways to top %d by max logp.", as.integer(max_pathways)))
    }
  }

  out_csv <- file.path(out_dir, "topic_pathway_enrichment_gsea.csv")
  data.table::fwrite(res_dt, out_csv)

  mat_dt <- data.table::dcast(res_dt, pathway ~ topic, value.var = "logp", fill = 0)
  mat <- as.matrix(mat_dt[, -1, with = FALSE])
  rownames(mat) <- mat_dt$pathway
  colnames(mat) <- paste0("Topic", colnames(mat))

  main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "GSEA enrichment", sep = " | ")
  } else {
    "GSEA enrichment"
  }

  n_topics <- ncol(mat)
  width <- max(10, n_topics * 0.6)
  height <- max(8, min(160, nrow(mat) * 0.25))
  font_row <- if (nrow(mat) > 200) 4 else if (nrow(mat) > 120) 5 else 7
  font_col <- if (n_topics > 60) 5 else if (n_topics > 40) 6 else if (n_topics > 25) 8 else 10

  gsea_heatmap <- file.path(out_dir, "topic_pathway_enrichment_gsea_heatmap.pdf")
  grDevices::pdf(gsea_heatmap, width = width, height = height)
  tryCatch(
    {
      ph <- suppressWarnings(pheatmap::pheatmap(
        mat,
        cluster_rows = nrow(mat) > 1L,
        cluster_cols = ncol(mat) > 1L,
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize_row = font_row,
        fontsize_col = font_col,
        angle_col = 90,
        main = main_title,
        legend_labels = expression(-log[10]~"(adj.P)"),
        border_color = NA,
        silent = TRUE
      ))
      if (!is.null(ph$gtable)) {
        idx <- which(ph$gtable$layout$name %in% c("row_names", "col_names", "main"))
        for (i in idx) {
          ph$gtable$grobs[[i]] <- .set_grob_fontface(ph$gtable$grobs[[i]], "bold")
        }
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
      } else {
        pheatmap::pheatmap(mat, border_color = NA)
      }
    },
    finally = grDevices::dev.off()
  )

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    dot_path <- file.path(out_dir, "topic_pathway_enrichment_gsea_dotplot.pdf")
    plot_dt <- data.table::copy(res_dt)
    plot_dt[, topic_num := as.integer(topic)]
    plot_dt[, topic := paste0("Topic", topic_num)]
    plot_dt[, size_val := abs(NES)]
    path_order <- plot_dt[, .(
      min_topic = min(topic_num, na.rm = TRUE),
      max_score = max(size_val, na.rm = TRUE)
    ), by = pathway]
    path_order <- path_order[order(min_topic, -max_score)]
    plot_dt[, pathway := factor(pathway, levels = rev(path_order$pathway))]
    plot_dt[, topic := factor(topic, levels = paste0("Topic", sort(unique(topic_num))))]

    n_topics_plot <- length(unique(plot_dt$topic))
    n_paths_plot <- length(unique(plot_dt$pathway))
    wrap_label <- function(x, width = 70) {
      vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
    }
    size_range <- if (n_paths_plot > 80) c(0.6, 5) else c(1, 8)

    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = topic, y = pathway, color = logp, size = size_val)
    ) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::scale_color_gradient(
        low = "#2c7bb6",
        high = "#d7191c",
        name = expression(-log[10]~"(adj.P)")
      ) +
      ggplot2::scale_size_continuous(name = "abs(NES)", range = size_range) +
      ggplot2::scale_y_discrete(labels = function(x) wrap_label(x, width = 70)) +
      ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
      ggplot2::labs(
        x = "Topic",
        y = NULL,
        title = if (!is.null(title_prefix) && nzchar(title_prefix)) {
          paste(title_prefix, "GSEA dot plot", sep = " | ")
        } else {
          "GSEA dot plot"
        }
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        axis.text.y = ggplot2::element_text(size = if (n_paths_plot > 80) 6 else 8, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 9, face = "bold", angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major = ggplot2::element_line(color = "#e6e6e6"),
        plot.margin = ggplot2::margin(10, 30, 10, 10)
      )

    ggplot2::ggsave(
      dot_path,
      p,
      width = max(8, n_topics_plot * 0.6),
      height = min(50, max(8, n_paths_plot * 0.35)),
      limitsize = FALSE
    )
    dot_csv <- sub("\\.pdf$", ".csv", dot_path)
    ord_dt <- plot_dt[order(
      match(pathway, levels(pathway)),
      match(topic, levels(topic))
    )]
    data.table::fwrite(
      ord_dt[, .(topic, pathway, NES, logp, padj, leading_edge)],
      dot_csv
    )
    log_msg(sprintf("GSEA dot plot table saved to: %s", dot_csv))
    log_msg(sprintf("GSEA dot plot saved to: %s", dot_path))
  } else {
    log_msg("Skipping GSEA dot plot: ggplot2 not installed.")
  }

  invisible(TRUE)
}

plot_topic_pathway_enrichment_heatmap <- function(topic_terms,
                                                  edges_docs,
                                                  option_label,
                                                  out_file,
                                                  theta = NULL,
                                                  dbs = NULL,
                                                  pathway_species = NULL,
                                                  padj_cut = 0.05,
                                                  min_genes = 5L,
                                                  top_n_per_topic = 20L,
                                                  dot_top_n_per_topic = 25L,
                                                  max_pathways = 200L,
                                                  title_prefix = NULL,
                                                  use_all_terms = FALSE,
                                                  make_heatmap = FALSE,
                                                  make_dotplot = TRUE,
                                                  tf_link_mode = c("theta", "none"),
                                                  tf_theta_top_n = 50L,
                                                  tf_theta_min = NA_real_,
                                                  enrichr_sleep_time = 0,
                                                  enrichr_cache_dir = NULL,
                                                  enrichr_n_cores = 1L,
                                                  pathway_backend = NULL) {
  .assert_pkg("data.table")
  tf_link_mode <- match.arg(tf_link_mode)
  enrichr_sleep_time <- .normalize_enrichr_sleep_time(enrichr_sleep_time)
  enrichr_n_cores <- .normalize_enrichr_n_cores(enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  pathway_species_mode <- .normalize_pathway_species_mode(pathway_species)
  human_mouse_best <- identical(pathway_species_mode, "human_mouse_best")
  if (is.null(dbs)) dbs <- .default_pathway_databases(pathway_species_mode)
  dbs_by_species <- if (human_mouse_best) .split_pathway_databases_by_species(dbs) else NULL
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- .module3_default_enrichr_cache_dir(dirname(out_file), backend = pathway_backend)
  }
  log_path <- file.path(dirname(out_file), "topic_pathway_enrichment_debug.txt")
  log_msg <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
  }
  .quiet_enrichr_call <- function(expr) {
    val <- NULL
    utils::capture.output(
      val <- suppressMessages(eval.parent(substitute(expr))),
      type = "output"
    )
    val
  }

  if (!.pathway_backend_available(pathway_backend)) {
    backend_label <- if (identical(pathway_backend, "enrichly")) "enrichly" else "enrichR"
    msg <- sprintf("Skipping pathway enrichment heatmap: selected backend is not installed: %s.", backend_label)
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (identical(pathway_backend, "enrichr")) {
    .ensure_enrichr_ready(
      site = "Enrichr",
      verbose = TRUE,
      log_fun = log_msg
    )
  }

  if (is.null(topic_terms) || !nrow(topic_terms)) {
    msg <- "Skipping pathway enrichment heatmap: topic_terms empty."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (!isTRUE(use_all_terms)) {
    if (!("in_topic" %in% names(topic_terms))) {
      msg <- "Skipping pathway enrichment heatmap: topic_terms missing in_topic."
      .log_inform(msg)
      log_msg(msg)
      return(invisible(NULL))
    }
    if (!any(.as_logical_flag(topic_terms$in_topic))) {
      msg <- "Skipping pathway enrichment heatmap: topic_terms has no in_topic == TRUE."
      .log_inform(msg)
      log_msg(msg)
      return(invisible(NULL))
    }
  }
  if (option_label != "opt3_gene_fc_expr") {
    if (is.null(edges_docs)) {
      msg <- "Skipping pathway enrichment heatmap: edges_docs is NULL for peak/joint option."
      .log_inform(msg)
      log_msg(msg)
      return(invisible(NULL))
    }
    if (!all(c("peak_id", "gene_key") %in% names(edges_docs))) {
      msg <- "Skipping pathway enrichment heatmap: edges_docs missing peak_id or gene_key."
      .log_inform(msg)
      log_msg(msg)
      return(invisible(NULL))
    }
  }

  gene_sets <- topic_gene_sets_from_terms(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    option_label = option_label,
    use_all_terms = use_all_terms,
    include_peak_terms = !identical(option_label, "opt3_gene_fc_expr")
  )
  if (!length(gene_sets)) {
    msg <- "Skipping pathway enrichment heatmap: no topic gene sets after mapping."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  log_msg(sprintf("Option: %s", option_label))
  log_msg("Gene set source: assigned topic terms.")
  log_msg(sprintf("Topics with gene sets: %s", paste(names(gene_sets), collapse = ",")))
  if (human_mouse_best) {
    log_msg("Pathway species mode: human_mouse_best.")
    log_msg(sprintf("Human DBs: %s", paste(dbs_by_species$human, collapse = ",")))
    log_msg(sprintf("Mouse DBs: %s", paste(dbs_by_species$mouse, collapse = ",")))
  } else {
    log_msg(sprintf("Pathway species: %s", pathway_species_mode))
    log_msg(sprintf("DBs: %s", paste(dbs, collapse = ",")))
  }
  gene_set_dt <- data.table::rbindlist(
    lapply(names(gene_sets), function(nm) {
      data.table::data.table(topic = as.integer(nm), gene = gene_sets[[nm]])
    }),
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(gene_set_dt)) {
    data.table::fwrite(
      gene_set_dt,
      file.path(dirname(out_file), "topic_pathway_gene_sets.csv")
    )
  }

  res_list <- vector("list", length(gene_sets))
  names(res_list) <- names(gene_sets)
  run_one_enrichr_species <- function(nm, species, dbs_use, cache_dir_use) {
    genes <- gene_sets[[nm]]
    log_msg(sprintf("Topic %s %s gene count: %d", nm, species, length(genes)))
    if (length(genes) < as.integer(min_genes)) {
      log_msg(sprintf("Topic %s %s skipped: gene count < %d", nm, species, min_genes))
      return(NULL)
    }
    if (!length(dbs_use)) {
      log_msg(sprintf("Topic %s %s skipped: no databases.", nm, species))
      return(NULL)
    }
    enr <- tryCatch(
      .quiet_enrichr_call(.run_enrichr_cached(
        genes = genes,
        dbs = dbs_use,
        sleep_time = enrichr_sleep_time,
        cache_dir = cache_dir_use,
        backend = pathway_backend,
        pathway_species = species
      )),
      error = function(e) {
        log_msg(sprintf("Topic %s %s enrichr error: %s", nm, species, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(enr)) {
      log_msg(sprintf("Topic %s %s enrichr returned NULL.", nm, species))
      return(NULL)
    }
    out <- .topic_enrichr_result_to_table(enr, topic_name = nm, genes = genes, pathway_species = species)
    out <- out[is.finite(padj) & padj <= as.numeric(padj_cut)]
    n_hits <- nrow(out)
    log_msg(sprintf("Topic %s %s enriched pathways: %d (padj<=%s)", nm, species, n_hits, padj_cut))
    out
  }
  run_one_enrichr <- function(nm) {
    if (human_mouse_best) {
      out <- lapply(c("human", "mouse"), function(species) {
        run_one_enrichr_species(
          nm = nm,
          species = species,
          dbs_use = dbs_by_species[[species]],
          cache_dir_use = file.path(enrichr_cache_dir, species)
        )
      })
      return(data.table::rbindlist(out, use.names = TRUE, fill = TRUE))
    }
    run_one_enrichr_species(
      nm = nm,
      species = pathway_species_mode,
      dbs_use = dbs,
      cache_dir_use = enrichr_cache_dir
    )
  }
  set_names <- names(gene_sets)
  if (.Platform$OS.type != "windows" && enrichr_n_cores > 1L && length(set_names) > 1L) {
    res_list <- parallel::mclapply(
      set_names,
      run_one_enrichr,
      mc.cores = min(enrichr_n_cores, length(set_names))
    )
    names(res_list) <- set_names
  } else {
    for (nm in set_names) {
      res_list[[nm]] <- run_one_enrichr(nm)
    }
  }

  res_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  if (!nrow(res_dt)) {
    msg <- "Skipping pathway enrichment heatmap: no enriched terms at padj_cut."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  res_dt <- res_dt[is.finite(logp) & logp > 0]
  if (!nrow(res_dt)) {
    msg <- "Skipping pathway enrichment heatmap: all logp values non-finite or zero."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (human_mouse_best) {
    data.table::fwrite(
      res_dt,
      file.path(dirname(out_file), "topic_pathway_enrichment_topic_terms_human_mouse_all.csv")
    )
    res_dt <- .select_best_human_mouse_pathways(res_dt)
    log_msg(sprintf("Human/mouse best rows selected: %d.", nrow(res_dt)))
  }
  if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) {
    res_dt <- res_dt[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
  } else {
    res_dt <- res_dt[order(topic, -logp)]
  }
  if (nrow(res_dt) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
    path_rank <- res_dt[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
    if (nrow(path_rank) > as.integer(max_pathways)) {
      keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
      res_dt <- res_dt[pathway %in% keep]
      log_msg(sprintf("Filtered pathways to top %d by max logp.", as.integer(max_pathways)))
    }
  }
  res_dt[, topic := as.integer(topic)]
  data.table::setorder(res_dt, topic)
  data.table::fwrite(
    res_dt,
    file.path(dirname(out_file), "topic_pathway_enrichment_topic_terms.csv")
  )
  log_msg(sprintf("Total enriched pathways (unique): %d", length(unique(res_dt$pathway))))
  log_msg(sprintf("Debug log written to: %s", log_path))

  if (isTRUE(make_heatmap)) {
    log_msg("Skipping pathway heatmap: heatmap output has been removed from standard Module 3 reports.")
  }

  if (isTRUE(make_dotplot) && requireNamespace("ggplot2", quietly = TRUE)) {
    dot_path <- file.path(dirname(out_file), "topic_pathway_enrichment_dotplot.pdf")
    plot_dt <- data.table::copy(res_dt)
    plot_dt <- plot_dt[is.finite(padj)]
    plot_dt[, topic_num := as.integer(topic)]
    plot_dt <- plot_dt[is.finite(topic_num)]
    if (!("overlap_hits" %in% names(plot_dt))) {
      plot_dt[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
    }
    plot_dt[is.na(overlap_hits), overlap_hits := 0L]
    plot_dt[, pathway_key := pathway]
    dot_n <- suppressWarnings(as.integer(dot_top_n_per_topic[[1L]]))
    if (!is.finite(dot_n) || dot_n <= 0L) {
      dot_n <- suppressWarnings(as.integer(top_n_per_topic[[1L]]))
    }
    if (!is.finite(dot_n) || dot_n <= 0L) dot_n <- 12L
    top_dt <- if (nrow(plot_dt)) {
      plot_dt[, {
        ranked <- .SD[order(padj, -overlap_hits, pathway)]
        sig <- ranked[padj <= as.numeric(padj_cut)]
        if (nrow(sig)) head(sig, dot_n) else head(ranked, dot_n)
      }, by = topic_num]
    } else {
      plot_dt
    }
    if (nrow(top_dt)) {
      plot_dt <- plot_dt[pathway_key %in% unique(top_dt$pathway_key)]
      sig_keys <- unique(plot_dt[padj <= as.numeric(padj_cut), pathway_key])
      plot_dt <- plot_dt[pathway_key %in% sig_keys]
    }
    if (nrow(plot_dt) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
      path_rank <- plot_dt[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway_key]
      if (nrow(path_rank) > as.integer(max_pathways)) {
        keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway_key]
        plot_dt <- plot_dt[pathway_key %in% keep]
      }
    }
    if (!nrow(plot_dt)) {
      log_msg("Skipping pathway dot plot: no significant pathway terms.")
      return(invisible(TRUE))
    }
    top_key <- unique(top_dt[, .(topic_num, pathway_key)])
    top_key[, is_topic_top_pathway := TRUE]
    plot_dt <- merge(plot_dt, top_key, by = c("topic_num", "pathway_key"), all.x = TRUE, sort = FALSE)
    plot_dt[is.na(is_topic_top_pathway), is_topic_top_pathway := FALSE]
    topic_levels <- sort(unique(plot_dt$topic_num))
    row_levels <- character(0)
    for (tp in topic_levels) {
      owned <- plot_dt[topic_num == tp & is_topic_top_pathway == TRUE]
      if (!nrow(owned)) next
      data.table::setorder(owned, -logp, -overlap_hits, pathway)
      for (key in owned$pathway_key) {
        if (!(key %in% row_levels)) row_levels <- c(row_levels, key)
      }
    }
    remaining <- setdiff(unique(plot_dt$pathway_key), row_levels)
    if (length(remaining)) {
      remain_rank <- plot_dt[pathway_key %in% remaining, .(
        max_logp = max(logp, na.rm = TRUE),
        max_overlap = max(overlap_hits, na.rm = TRUE)
      ), by = pathway_key]
      data.table::setorder(remain_rank, -max_logp, -max_overlap, pathway_key)
      row_levels <- c(row_levels, remain_rank$pathway_key)
    }
    if (is.finite(max_pathways) && as.numeric(max_pathways) > 0 &&
        length(row_levels) > as.integer(max_pathways)) {
      row_levels <- row_levels[seq_len(as.integer(max_pathways))]
      plot_dt <- plot_dt[pathway_key %in% row_levels]
    }
    plot_dt[, pathway_key := factor(pathway_key, levels = rev(row_levels))]
    plot_dt[, topic_label := factor(paste0("Topic", topic_num), levels = paste0("Topic", topic_levels))]
    plot_dt[, score_val := data.table::fifelse(is.finite(combined_score), combined_score, NA_real_)]
    plot_dt[!is.finite(score_val) & is.finite(odds_ratio), score_val := odds_ratio]
    plot_dt[!is.finite(score_val) & is.finite(logp), score_val := logp]
    plot_dt[, score_val := pmax(score_val, 0)]
    size_cap <- suppressWarnings(stats::quantile(plot_dt$score_val[is.finite(plot_dt$score_val)], probs = 0.95, na.rm = TRUE, names = FALSE))
    if (is.finite(size_cap)) plot_dt[, score_val := pmin(score_val, size_cap)]
    plot_dt[, is_sig := is.finite(logp) & logp >= 1.3]
    dot_csv <- sub("\\.pdf$", ".csv", dot_path)
    ord_dt <- plot_dt[order(
      match(pathway_key, levels(pathway_key)),
      match(topic_label, levels(topic_label))
    )]
    data.table::fwrite(ord_dt, dot_csv)

    family <- if (exists(".diff_grn_pathway_font_family", mode = "function")) {
      .diff_grn_pathway_font_family()
    } else {
      "sans"
    }
    n_topics_plot <- data.table::uniqueN(plot_dt$topic_label)
    n_paths_plot <- data.table::uniqueN(plot_dt$pathway_key)
    max_val <- max(plot_dt$logp[is.finite(plot_dt$logp)], 1.300001, na.rm = TRUE)
    color_values <- if (exists(".diff_grn_pathway_rescale", mode = "function")) {
      .diff_grn_pathway_rescale(c(0, 1.3, seq(1.3 + 1e-6, max_val, length.out = 6)))
    } else {
      color_breaks <- c(0, 1.3, seq(1.3 + 1e-6, max_val, length.out = 6))
      color_rng <- range(color_breaks, finite = TRUE)
      if (!all(is.finite(color_rng)) || diff(color_rng) == 0) {
        rep(0, length(color_breaks))
      } else {
        (color_breaks - color_rng[1]) / diff(color_rng)
      }
    }

    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = topic_label, y = pathway_key, size = score_val)
    ) +
      ggplot2::geom_point(ggplot2::aes(fill = logp), shape = 21, color = "transparent", alpha = 0.85) +
      ggplot2::geom_point(
        data = plot_dt[is_sig == TRUE],
        ggplot2::aes(fill = logp, color = logp),
        shape = 21,
        stroke = 0.7,
        alpha = 0.95
      ) +
      ggplot2::scale_fill_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        name = "-log10 adjusted p-value"
      ) +
      ggplot2::scale_color_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        guide = "none"
      ) +
      ggplot2::scale_size_continuous(name = "Combined score", range = c(0.8, 3.2)) +
      ggplot2::scale_y_discrete(labels = function(x) {
        y <- as.character(x)
        y <- gsub("\\s+", " ", y)
        too_long <- nchar(y) > 58L
        y[too_long] <- paste0(substr(y[too_long], 1L, 55L), "...")
        y
      }) +
      ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
      ggplot2::labs(
        x = "Topic",
        y = "Pathway term",
        title = if (!is.null(title_prefix) && nzchar(title_prefix)) {
          paste(title_prefix, "Pathway dot plot", sep = " | ")
        } else {
          "Pathway dot plot"
        },
        caption = if (is.finite(dot_n) && as.numeric(dot_n) > 0) {
          sprintf("Union of top %d pathways per topic; all occurrences shown.", as.integer(dot_n))
        } else {
          sprintf("All significant pathways per topic; adjusted p-value <= %s.", as.character(padj_cut))
        }
      ) +
      ggplot2::theme_bw(base_size = 9, base_family = family) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", hjust = 0.5),
        plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
      )

    if (exists(".diff_grn_pathway_save_fixed_panel_pdf", mode = "function")) {
      .diff_grn_pathway_save_fixed_panel_pdf(
        p,
        dot_path,
        panel_width = max(1.2, 0.11 * n_topics_plot),
        panel_height = max(1.1, 0.17 * n_paths_plot),
        family = family
      )
    } else {
      ggplot2::ggsave(
        dot_path,
        p,
        width = max(8, n_topics_plot * 0.25),
        height = max(6, n_paths_plot * 0.2),
        limitsize = FALSE
      )
    }
    data.table::fwrite(
      ord_dt[, .(
        topic = topic_num,
        pathway_key,
        pathway,
        padj,
        pval,
        overlap,
        overlap_hits,
        genes,
        logp,
        combined_score,
        odds_ratio,
        cluster_size,
        is_topic_top_pathway,
        topic_num,
        topic_label,
        score_val,
        is_sig
      )],
      dot_csv
    )
    log_msg(sprintf("Dot plot table saved to: %s", dot_csv))
    log_msg(sprintf("Dot plot saved to: %s", dot_path))
  } else if (isTRUE(make_dotplot)) {
    log_msg("Skipping pathway dot plot: ggplot2 not installed.")
  } else {
    log_msg("Skipping pathway dot plot: make_dotplot = FALSE.")
  }

  invisible(TRUE)
}

plot_topic_pathway_enrichment_by_comparison_terms <- function(topic_terms,
                                                              edges_docs,
                                                              theta,
                                                              out_dir,
                                                              title_prefix = NULL,
                                                              dbs = NULL,
                                                              pathway_species = NULL,
                                                              padj_cut = 0.05,
                                                              theta_min = 0.3,
                                                              include_peak_terms = TRUE,
                                                              use_all_terms = FALSE,
                                                              per_comparison_dir = ".",
                                                              per_comparison_flat = TRUE,
                                                              split_direction = TRUE,
                                                              min_genes = 5L,
                                                              top_n_per_topic = 20L,
                                                              dot_top_n_per_topic = Inf,
                                                              max_pathways = 200L,
                                                              make_dotplot = TRUE,
                                                              enrichr_sleep_time = 0,
                                                              enrichr_cache_dir = NULL,
                                                              enrichr_n_cores = 1L,
                                                              pathway_backend = NULL,
                                                              background_size = 20000L,
                                                              overall_pathway_file = NULL,
                                                              overwrite = FALSE,
                                                              doc_design = c("comparison", "condition")) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  group_label <- if (identical(doc_design, "condition")) "condition" else "comparison"
  group_prefix <- if (identical(doc_design, "condition")) "per_condition" else "per_comparison"
  enrichr_sleep_time <- .normalize_enrichr_sleep_time(enrichr_sleep_time)
  enrichr_n_cores <- .normalize_enrichr_n_cores(enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  pathway_species_mode <- .normalize_pathway_species_mode(pathway_species)
  human_mouse_best <- identical(pathway_species_mode, "human_mouse_best")
  if (is.null(dbs)) dbs <- .default_pathway_databases(pathway_species_mode)
  dbs_by_species <- if (human_mouse_best) .split_pathway_databases_by_species(dbs) else NULL
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir, backend = pathway_backend)
  }
  out_dir_pc <- if (is.null(per_comparison_dir) || !nzchar(as.character(per_comparison_dir)[[1L]]) || identical(as.character(per_comparison_dir)[[1L]], ".")) {
    out_dir
  } else {
    file.path(out_dir, per_comparison_dir)
  }
  dir.create(out_dir_pc, recursive = TRUE, showWarnings = FALSE)
  if (is.null(title_prefix)) title_prefix <- basename(out_dir)
  log_path <- file.path(out_dir_pc, paste0(group_prefix, "_topic_pathway_debug.txt"))
  log_msg <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
  }
  retest <- .topic_pathway_retest_from_overall(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    out_dir = out_dir,
    option_label = "joint",
    include_peak_terms = include_peak_terms,
    use_all_terms = use_all_terms,
    per_comparison_dir = per_comparison_dir,
    split_direction = split_direction,
    background_size = background_size,
    overall_pathway_file = overall_pathway_file,
    pathway_species = pathway_species,
    doc_design = doc_design
  )
  log_msg(sprintf(
    "Per-%s pathway retest from overall table wrote %d row(s).",
    group_label,
    nrow(retest)
  ))
  log_msg(sprintf("No per-%s Enrichr/enrichly calls were run in this step.", group_label))
  return(invisible(retest))

  .quiet_enrichr_call <- function(expr) {
    val <- NULL
    utils::capture.output(
      val <- suppressMessages(eval.parent(substitute(expr))),
      type = "output"
    )
    val
  }
  .select_topic_dotplot_terms <- function(res_dt) {
    dt_sel <- data.table::copy(res_dt)
    dt_sel <- dt_sel[is.finite(padj)]
    if (!nrow(dt_sel)) return(dt_sel)
    dt_sel[, topic := as.integer(topic)]
    dt_sel <- dt_sel[is.finite(topic)]
    if (!nrow(dt_sel)) return(dt_sel)
    if (!("overlap_hits" %in% names(dt_sel))) {
      dt_sel[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
    }
    dt_sel[is.na(overlap_hits), overlap_hits := 0L]
    dt_sel[, logp := -log10(pmax(padj, 1e-300))]
    dt_sel[, pathway_key := pathway]
    if (is.finite(dot_top_n_per_topic) && as.numeric(dot_top_n_per_topic) > 0) {
      top_n <- max(1L, as.integer(dot_top_n_per_topic))
      top_dt <- dt_sel[, {
        ranked <- .SD[order(padj, -overlap_hits, pathway)]
        sig <- ranked[padj <= as.numeric(padj_cut)]
        if (nrow(sig)) head(sig, top_n) else head(ranked, top_n)
      }, by = topic]
      if (!nrow(top_dt)) return(top_dt)
      out <- dt_sel[pathway_key %in% unique(top_dt$pathway_key)]
      sig_keys <- unique(out[padj <= as.numeric(padj_cut), pathway_key])
      out <- out[pathway_key %in% sig_keys]
    } else {
      top_dt <- dt_sel[padj <= as.numeric(padj_cut)]
      if (!nrow(top_dt)) return(top_dt)
      out <- data.table::copy(top_dt)
    }
    if (!nrow(out)) return(out)
    top_key <- unique(top_dt[, .(topic, pathway_key)])
    top_key[, is_topic_top_pathway := TRUE]
    out <- merge(out, top_key, by = c("topic", "pathway_key"), all.x = TRUE, sort = FALSE)
    out[is.na(is_topic_top_pathway), is_topic_top_pathway := FALSE]
    topic_levels <- sort(unique(out$topic))
    row_levels <- character(0)
    for (tp in topic_levels) {
      owned <- out[topic == tp & is_topic_top_pathway == TRUE]
      if (!nrow(owned)) next
      data.table::setorder(owned, -logp, -overlap_hits, pathway)
      for (key in owned$pathway_key) {
        if (!(key %in% row_levels)) row_levels <- c(row_levels, key)
      }
    }
    remaining <- setdiff(unique(out$pathway_key), row_levels)
    if (length(remaining)) {
      remain_rank <- out[pathway_key %in% remaining, .(
        max_logp = max(logp, na.rm = TRUE),
        max_overlap = max(overlap_hits, na.rm = TRUE)
      ), by = pathway_key]
      data.table::setorder(remain_rank, -max_logp, -max_overlap, pathway_key)
      row_levels <- c(row_levels, remain_rank$pathway_key)
    }
    if (nrow(out) && is.finite(max_pathways) && as.numeric(max_pathways) > 0 && length(row_levels) > as.integer(max_pathways)) {
      row_levels <- row_levels[seq_len(as.integer(max_pathways))]
      out <- out[pathway_key %in% row_levels]
    }
    out[, pathway_key := factor(pathway_key, levels = rev(row_levels))]
    data.table::setorder(out, topic, -is_topic_top_pathway, padj, -overlap_hits, pathway)
    out[]
  }
  .write_dotplot <- function(res_dt, dot_prefix, plot_title) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      log_msg("Skipping dot plot: ggplot2 not installed.")
      return(invisible(NULL))
    }
    plot_dt <- .select_topic_dotplot_terms(res_dt)
    if (!nrow(plot_dt)) {
      log_msg("Skipping dot plot: no significant pathway terms.")
      return(invisible(NULL))
    }
    plot_dt[, topic_num := as.integer(topic)]
    topic_levels <- sort(unique(plot_dt$topic_num))
    plot_dt[, topic_label := factor(paste0("Topic", topic_num), levels = paste0("Topic", topic_levels))]
    plot_dt[, score_val := data.table::fifelse(is.finite(combined_score), combined_score, NA_real_)]
    plot_dt[!is.finite(score_val) & is.finite(odds_ratio), score_val := odds_ratio]
    plot_dt[!is.finite(score_val) & is.finite(logp), score_val := logp]
    plot_dt[, score_val := pmax(score_val, 0)]
    size_cap <- suppressWarnings(stats::quantile(plot_dt$score_val[is.finite(plot_dt$score_val)], probs = 0.95, na.rm = TRUE, names = FALSE))
    if (is.finite(size_cap)) plot_dt[, score_val := pmin(score_val, size_cap)]
    plot_dt[, is_sig := is.finite(logp) & logp >= 1.3]
    dot_path <- paste0(dot_prefix, ".pdf")
    dot_csv <- paste0(dot_prefix, ".csv")
    ord_dt <- plot_dt[order(
      match(pathway_key, levels(pathway_key)),
      match(topic_label, levels(topic_label))
    )]
    data.table::fwrite(ord_dt, dot_csv)

    family <- if (exists(".diff_grn_pathway_font_family", mode = "function")) {
      .diff_grn_pathway_font_family()
    } else {
      "sans"
    }
    n_topics_plot <- data.table::uniqueN(plot_dt$topic_label)
    n_paths_plot <- data.table::uniqueN(plot_dt$pathway_key)
    max_val <- max(plot_dt$logp[is.finite(plot_dt$logp)], 1.300001, na.rm = TRUE)
    color_breaks <- c(0, 1.3, seq(1.3 + 1e-6, max_val, length.out = 6))
    color_rng <- range(color_breaks, finite = TRUE)
    color_values <- if (!all(is.finite(color_rng)) || diff(color_rng) == 0) {
      rep(0, length(color_breaks))
    } else {
      (color_breaks - color_rng[1]) / diff(color_rng)
    }
    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = topic_label, y = pathway_key, size = score_val)
    ) +
      ggplot2::geom_point(ggplot2::aes(fill = logp), shape = 21, color = "transparent", alpha = 0.85) +
      ggplot2::geom_point(
        data = plot_dt[is_sig == TRUE],
        ggplot2::aes(fill = logp, color = logp),
        shape = 21,
        stroke = 0.7,
        alpha = 0.95
      ) +
      ggplot2::scale_fill_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        name = "-log10 adjusted p-value"
      ) +
      ggplot2::scale_color_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        guide = "none"
      ) +
      ggplot2::scale_size_continuous(name = "Combined score", range = c(0.8, 3.2)) +
      ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
      ggplot2::scale_y_discrete(labels = function(x) {
        y <- as.character(x)
        y <- gsub("\\s+", " ", y)
        too_long <- nchar(y) > 58L
        y[too_long] <- paste0(substr(y[too_long], 1L, 55L), "...")
        y
      }) +
      ggplot2::labs(
        x = "Topic",
        y = "Pathway term",
        title = plot_title,
        caption = if (is.finite(dot_top_n_per_topic) && as.numeric(dot_top_n_per_topic) > 0) {
          sprintf("Union of top %d pathways per topic; all occurrences shown.", as.integer(dot_top_n_per_topic))
        } else {
          sprintf("All significant pathways per topic; adjusted p-value <= %s.", as.character(padj_cut))
        }
      ) +
      ggplot2::theme_bw(base_size = 9, base_family = family) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", hjust = 0.5),
        plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
      )
    tryCatch(
      {
        if (exists(".diff_grn_pathway_save_fixed_panel_pdf", mode = "function")) {
          .diff_grn_pathway_save_fixed_panel_pdf(
            p,
            dot_path,
            panel_width = max(1.2, 0.11 * n_topics_plot),
            panel_height = max(1.1, 0.17 * n_paths_plot),
            family = family
          )
        } else {
          ggplot2::ggsave(dot_path, p, width = max(8, n_topics_plot * 0.25), height = max(6, n_paths_plot * 0.2), limitsize = FALSE)
        }
      },
      error = function(e) log_msg(sprintf("Dot plot failed: %s", conditionMessage(e)))
    )
    log_msg(sprintf("Dot plot CSV saved: %s", dot_csv))
    if (file.exists(dot_path)) log_msg(sprintf("Dot plot PDF saved: %s", dot_path))
    invisible(TRUE)
  }

  gene_sets_dt <- topic_gene_sets_by_comparison_terms(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    theta = theta,
    theta_min = theta_min,
    include_peak_terms = include_peak_terms,
    use_all_terms = use_all_terms,
    doc_design = doc_design,
    split_direction = split_direction
  )
  gene_set_file <- file.path(out_dir_pc, "topic_term_pathway_gene_sets.csv")
  data.table::fwrite(gene_sets_dt, gene_set_file)
  if (!nrow(gene_sets_dt)) {
    log_msg("Skipping per-comparison topic-term pathway enrichment: no theta-term gene sets.")
    return(invisible(NULL))
  }
  gene_set_summary <- gene_sets_dt[, .(
    n_genes = data.table::uniqueN(gene),
    n_gene_doc_pairs = sum(n_docs, na.rm = TRUE),
    max_gene_docs = max(n_docs, na.rm = TRUE),
    max_theta = max(max_theta, na.rm = TRUE)
  ), by = .(comparison_id, direction_group, topic)]
  data.table::setorder(gene_set_summary, comparison_id, direction_group, topic)
  data.table::fwrite(gene_set_summary, file.path(out_dir_pc, "topic_term_pathway_gene_set_summary.csv"))

  eligible <- gene_set_summary[n_genes >= as.integer(min_genes)]
  if (!nrow(eligible)) {
    log_msg(sprintf("Skipping enrichment: no comparison/direction/topic gene set has at least %d genes.", as.integer(min_genes)))
    return(invisible(TRUE))
  }
  if (!.pathway_backend_available(pathway_backend)) {
    backend_label <- if (identical(pathway_backend, "enrichly")) "enrichly" else "enrichR"
    msg <- sprintf("Skipping per-comparison topic-term pathway enrichment: selected backend is not installed: %s.", backend_label)
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (identical(pathway_backend, "enrichr")) {
    .ensure_enrichr_ready(site = "Enrichr", verbose = TRUE, log_fun = log_msg)
  }

  all_res <- list()
  key_dt <- unique(eligible[, .(comparison_id, direction_group)])
  for (i in seq_len(nrow(key_dt))) {
    cmp <- key_dt$comparison_id[[i]]
    dir_lab <- key_dt$direction_group[[i]]
    cmp_dir <- if (isTRUE(per_comparison_flat)) {
      out_dir_pc
    } else {
      file.path(out_dir_pc, .safe_filename(cmp))
    }
    dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)
    prefix <- file.path(cmp_dir, paste0(.safe_filename(cmp), "_", .safe_filename(dir_lab), "_topic_term_dotplot"))
    enrich_csv <- paste0(prefix, "_enrichment.csv")
    if (!isTRUE(overwrite) && file.exists(enrich_csv) && (!isTRUE(make_dotplot) || file.exists(paste0(prefix, ".pdf")))) {
      log_msg(sprintf("%s | %s: existing outputs found; skipping.", cmp, dir_lab))
      all_res[[length(all_res) + 1L]] <- data.table::fread(enrich_csv, showProgress = FALSE)
      next
    }
    sub_summary <- eligible[comparison_id == cmp & direction_group == dir_lab]
    sub_genes <- gene_sets_dt[comparison_id == cmp & direction_group == dir_lab]
    gene_sets <- lapply(sub_summary$topic, function(tp) {
      unique(sub_genes[topic == tp, gene])
    })
    names(gene_sets) <- as.character(sub_summary$topic)
    run_one <- function(nm) {
      genes <- gene_sets[[nm]]
      log_msg(sprintf("%s | %s | Topic %s gene count: %d", cmp, dir_lab, nm, length(genes)))
      enr <- tryCatch(
        .quiet_enrichr_call(.run_enrichr_cached(
          genes = genes,
          dbs = dbs,
          sleep_time = enrichr_sleep_time,
          cache_dir = enrichr_cache_dir,
          backend = pathway_backend,
          pathway_species = pathway_species
        )),
        error = function(e) {
          log_msg(sprintf("%s | %s | Topic %s enrichr error: %s", cmp, dir_lab, nm, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(enr)) return(NULL)
      .topic_enrichr_result_to_table(enr, topic_name = nm, genes = genes, pathway_species = pathway_species)
    }
    set_names <- names(gene_sets)
    if (.Platform$OS.type != "windows" && enrichr_n_cores > 1L && length(set_names) > 1L) {
      res_list <- parallel::mclapply(
        set_names,
        run_one,
        mc.cores = min(enrichr_n_cores, length(set_names))
      )
      names(res_list) <- set_names
    } else {
      res_list <- vector("list", length(gene_sets))
      names(res_list) <- set_names
      for (nm in set_names) res_list[[nm]] <- run_one(nm)
    }
    res_sub <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
    if (!nrow(res_sub)) {
      log_msg(sprintf("%s | %s: no enriched pathways.", cmp, dir_lab))
      next
    }
    res_sub <- res_sub[is.finite(logp) & logp > 0]
    if (!nrow(res_sub)) next
    if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) {
      res_sub <- res_sub[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
    } else {
      res_sub <- res_sub[order(topic, -logp)]
    }
    if (nrow(res_sub) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
      path_rank <- res_sub[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
      if (nrow(path_rank) > as.integer(max_pathways)) {
        keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
        res_sub <- res_sub[pathway %in% keep]
      }
    }
    res_sub[, `:=`(
      comparison_id = cmp,
      direction_group = dir_lab,
      pathway_source = "theta_topic_terms"
    )]
    data.table::setorder(res_sub, topic, padj)
    data.table::fwrite(res_sub, enrich_csv)
    all_res[[length(all_res) + 1L]] <- res_sub
    if (isTRUE(make_dotplot)) {
      plot_title <- paste(cmp, dir_lab, "Pathway enrichment (theta topic terms)", sep = " | ")
      .write_dotplot(res_sub, dot_prefix = prefix, plot_title = plot_title)
    }
  }
  res_dt <- data.table::rbindlist(all_res, use.names = TRUE, fill = TRUE)
  if (nrow(res_dt)) {
    data.table::setorder(res_dt, comparison_id, direction_group, topic, padj)
    data.table::fwrite(res_dt, file.path(out_dir_pc, "per_comparison_topic_pathway_enrichment.csv"))
  }
  invisible(TRUE)
}

plot_topic_pathway_enrichment_by_condition_terms <- function(topic_terms,
                                                             edges_docs,
                                                             theta,
                                                             out_dir,
                                                             ...) {
  plot_topic_pathway_enrichment_by_comparison_terms(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    theta = theta,
    out_dir = out_dir,
    doc_design = "condition",
    split_direction = FALSE,
    ...
  )
}

.module3_pathway_settings_from_audit <- function(extraction_dir) {
  path <- file.path(extraction_dir, "topic_pathway_enrichment_debug.txt")
  lines <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  value_after <- function(prefix) {
    hit <- grep(prefix, lines, fixed = TRUE, value = TRUE)
    if (!length(hit)) return("")
    sub(paste0("^.*", prefix), "", hit[[1L]])
  }
  species <- if (any(grepl("Pathway species mode: human_mouse_best", lines, fixed = TRUE))) {
    "human_mouse_best"
  } else {
    value_after("Pathway species: ")
  }
  if (!nzchar(species)) species <- "human"
  db_lines <- if (identical(species, "human_mouse_best")) {
    c(value_after("Human DBs: "), value_after("Mouse DBs: "))
  } else {
    value_after("DBs: ")
  }
  dbs <- unique(trimws(unlist(strsplit(db_lines, ",", fixed = TRUE))))
  dbs <- dbs[!is.na(dbs) & nzchar(dbs)]
  if (!length(dbs)) dbs <- .default_pathway_databases(species)
  list(species = species, databases = dbs)
}

.module3_copy_topic_space_pathway_outputs <- function(from_dir,
                                                       extraction_dir,
                                                       topic_space) {
  suffix <- paste0("_", topic_space)
  files <- c(
    topic_pathway_enrichment_topic_terms = "csv",
    topic_pathway_enrichment_dotplot = "csv",
    topic_pathway_enrichment_dotplot = "pdf",
    topic_pathway_gene_sets = "csv"
  )
  source_names <- paste0(names(files), ".", unname(files))
  target_names <- paste0(names(files), suffix, ".", unname(files))
  copied <- logical(length(source_names))
  for (i in seq_along(source_names)) {
    source <- file.path(from_dir, source_names[[i]])
    target <- file.path(extraction_dir, target_names[[i]])
    if (!file.exists(source)) next
    copied[[i]] <- file.copy(source, target, overwrite = TRUE)
  }
  if (identical(topic_space, "combined")) {
    standard <- c(
      "topic_pathway_enrichment_topic_terms.csv",
      "topic_pathway_enrichment_dotplot.csv",
      "topic_pathway_enrichment_dotplot.pdf",
      "topic_pathway_gene_sets.csv"
    )
    for (name in standard) {
      source <- file.path(from_dir, name)
      if (file.exists(source)) {
        file.copy(source, file.path(extraction_dir, name), overwrite = TRUE)
      }
    }
  }
  invisible(copied)
}

.module3_extraction_trained_k <- function(extraction_dir, model_dir = NULL) {
  directory_k <- suppressWarnings(as.integer(sub(
    "^K",
    "",
    basename(extraction_dir)
  )))
  if (is.finite(directory_k) && directory_k >= 2L) return(directory_k)

  settings_file <- file.path(
    extraction_dir,
    "topic_assignment_qc_settings.csv"
  )
  if (file.exists(settings_file)) {
    settings <- data.table::fread(settings_file, showProgress = FALSE)
    settings_k <- if ("raw_k" %in% names(settings)) {
      unique(as.integer(settings$raw_k))
    } else {
      integer()
    }
    settings_k <- settings_k[is.finite(settings_k) & settings_k >= 2L]
    if (length(settings_k) == 1L) return(settings_k)
  }

  if (!is.null(model_dir) && dir.exists(model_dir)) {
    theta_files <- list.files(
      file.path(model_dir, "vae_models"),
      pattern = "^theta_K[0-9]+[.]csv$",
      full.names = FALSE
    )
    model_k <- unique(suppressWarnings(as.integer(sub(
      "^theta_K([0-9]+)[.]csv$",
      "\\1",
      theta_files
    ))))
    model_k <- model_k[is.finite(model_k) & model_k >= 2L]
    if (length(model_k) == 1L) return(model_k)
  }

  .log_abort(
    paste0(
      "Could not infer one trained K for extraction directory: ",
      extraction_dir,
      ". Use a K<integer> directory or retain topic_assignment_qc_settings.csv."
    )
  )
}

.module3_rebuild_topic_space_pathways <- function(extraction_dir,
                                                   model_dir,
                                                   topic_space = c("raw", "combined"),
                                                   overwrite = FALSE,
                                                   enrichr_n_cores = 1L,
                                                   pathway_backend = NULL,
                                                   pathway_species = NULL,
                                                   pathway_databases = NULL,
                                                   verbose = TRUE) {
  topic_space <- match.arg(topic_space)
  extraction_dir <- normalizePath(extraction_dir, winslash = "/", mustWork = TRUE)
  model_dir <- normalizePath(model_dir, winslash = "/", mustWork = TRUE)
  settings <- .module3_pathway_settings_from_audit(extraction_dir)
  if (!is.null(pathway_species)) {
    settings$species <- .normalize_pathway_species_mode(pathway_species)
  }
  if (!is.null(pathway_databases)) {
    pathway_databases <- unique(as.character(pathway_databases))
    settings$databases <- pathway_databases[
      !is.na(pathway_databases) & nzchar(pathway_databases)
    ]
  }
  if (!length(settings$databases)) {
    settings$databases <- .default_pathway_databases(settings$species)
  }
  terms_file <- if (identical(topic_space, "raw")) {
    file.path(extraction_dir, "topic_terms_raw.csv")
  } else {
    file.path(extraction_dir, "topic_terms.csv")
  }
  if (!file.exists(terms_file)) {
    .log_abort("Missing {topic_space} topic terms: {terms_file}")
  }
  target_file <- file.path(
    extraction_dir,
    paste0("topic_pathway_enrichment_topic_terms_", topic_space, ".csv")
  )
  if (!isTRUE(overwrite) && file.exists(target_file)) {
    if (isTRUE(verbose)) {
      .log_inform(
        "Keeping existing {topic_space} pathways: {extraction_dir}"
      )
    }
    return(invisible(list(overall = target_file)))
  }
  topic_terms <- data.table::fread(terms_file, showProgress = FALSE)
  topic_term_edges <- data.table::data.table(
    peak_id = character(),
    gene_key = character()
  )
  trained_k <- .module3_extraction_trained_k(extraction_dir, model_dir)
  if (identical(topic_space, "raw")) {
    theta_file <- file.path(
      model_dir,
      "vae_models",
      sprintf("theta_K%d.csv", trained_k)
    )
    theta <- .m3tb_read_probability_csv(theta_file, "doc_id")
  } else {
    theta_file <- file.path(
      extraction_dir,
      "rds",
      "topic_theta_optimized.rds"
    )
    if (!file.exists(theta_file)) {
      .log_abort("Missing combined theta: {theta_file}")
    }
    theta <- readRDS(theta_file)
  }
  work_dir <- tempfile(
    paste0("module3_pathways_", topic_space, "_"),
    tmpdir = extraction_dir
  )
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  cache_dir <- file.path(
    extraction_dir,
    "cache",
    .pathway_backend(pathway_backend)
  )
  plot_topic_pathway_enrichment_heatmap(
    topic_terms = topic_terms,
    edges_docs = topic_term_edges,
    option_label = "joint",
    out_file = file.path(work_dir, "topic_pathway_enrichment_heatmap.pdf"),
    theta = theta,
    dbs = settings$databases,
    pathway_species = settings$species,
    padj_cut = 0.05,
    min_genes = 5L,
    top_n_per_topic = Inf,
    dot_top_n_per_topic = 30L,
    max_pathways = 1000L,
    title_prefix = sprintf("K%d | %s", trained_k, topic_space),
    use_all_terms = FALSE,
    make_heatmap = FALSE,
    make_dotplot = TRUE,
    tf_link_mode = "none",
    enrichr_cache_dir = cache_dir,
    enrichr_n_cores = enrichr_n_cores,
    pathway_backend = pathway_backend
  )
  overall_file <- file.path(
    work_dir,
    "topic_pathway_enrichment_topic_terms.csv"
  )
  if (!file.exists(overall_file)) {
    .log_abort(
      "No overall {topic_space} pathway table was produced for {extraction_dir}."
    )
  }
  .module3_copy_topic_space_pathway_outputs(
    work_dir,
    extraction_dir,
    topic_space
  )
  manifest <- data.table::data.table(
    topic_space = topic_space,
    trained_k = trained_k,
    topic_count = ncol(theta),
    terms_file = basename(terms_file),
    terms_hash = digest::digest(
      file = terms_file,
      algo = "xxhash64",
      serialize = FALSE
    ),
    theta_file = basename(theta_file),
    theta_hash = digest::digest(
      file = theta_file,
      algo = "xxhash64",
      serialize = FALSE
    ),
    pathway_species = settings$species,
    pathway_databases = paste(settings$databases, collapse = ";"),
    pathway_backend = .pathway_backend(pathway_backend),
    built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  manifest_file <- file.path(extraction_dir, "topic_space_pathway_manifest.csv")
  previous <- if (file.exists(manifest_file)) {
    data.table::fread(manifest_file, showProgress = FALSE)
  } else {
    data.table::data.table()
  }
  previous <- previous[topic_space != manifest$topic_space]
  data.table::fwrite(
    data.table::rbindlist(list(previous, manifest), use.names = TRUE, fill = TRUE),
    manifest_file
  )
  invisible(list(
    overall = target_file,
    manifest = manifest_file
  ))
}

plot_topic_pathway_enrichment_from_link_scores <- function(link_scores,
                                                           out_dir,
                                                           title_prefix = NULL,
                                                           file_tag = NULL,
                                                           dbs = NULL,
                                                           pathway_species = NULL,
                                                           padj_cut = 0.05,
                                                           min_genes = 5L,
                                                           top_n_per_topic = 20L,
                                                           dot_top_n_per_topic = Inf,
                                                           max_pathways = 200L,
                                                           include_tf = TRUE,
                                                           include_gene = TRUE,
                                                           min_prob = 0,
                                                           gene_min_prob = NULL,
                                                           tf_min_prob = NULL,
                                                           tf_max_topics = Inf,
                                              tf_top_n_per_topic = NA_integer_,
                                              tf_link_scores = NULL,
                                                           gene_terms = NULL,
                                                           per_comparison = FALSE,
                                                           per_comparison_dir = "per_comparison_pathway",
                                                           per_comparison_flat = FALSE,
                                                           split_direction = TRUE,
                                                           make_heatmap = FALSE,
                                                           make_dotplot = TRUE,
                                                           enrichr_sleep_time = 0,
                                                           enrichr_cache_dir = NULL,
                                                           enrichr_n_cores = 1L,
                                                           pathway_backend = NULL,
                                                           overwrite = FALSE,
                                                           doc_design = c("comparison", "condition")) {
  .assert_pkg("data.table")
  doc_design <- match.arg(doc_design)
  enrichr_sleep_time <- .normalize_enrichr_sleep_time(enrichr_sleep_time)
  enrichr_n_cores <- .normalize_enrichr_n_cores(enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  pathway_species_mode <- .normalize_pathway_species_mode(pathway_species)
  human_mouse_best <- identical(pathway_species_mode, "human_mouse_best")
  if (is.null(dbs)) dbs <- .default_pathway_databases(pathway_species_mode)
  dbs_by_species <- if (human_mouse_best) .split_pathway_databases_by_species(dbs) else NULL
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir, backend = pathway_backend)
  }
  overwrite <- isTRUE(overwrite)
  log_tag <- if (!is.null(file_tag) && nzchar(file_tag)) paste0("_", file_tag) else ""
  log_path <- file.path(out_dir, paste0("topic_pathway_enrichment_links_debug", log_tag, ".txt"))
  log_msg <- function(msg) {
    stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] %s\n", stamp, msg), file = log_path, append = TRUE)
  }
  .quiet_enrichr_call <- function(expr) {
    val <- NULL
    utils::capture.output(
      val <- suppressMessages(eval.parent(substitute(expr))),
      type = "output"
    )
    val
  }

  if (!.pathway_backend_available(pathway_backend)) {
    backend_label <- if (identical(pathway_backend, "enrichly")) "enrichly" else "enrichR"
    msg <- sprintf("Skipping link-score pathway enrichment: selected backend is not installed: %s.", backend_label)
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (identical(pathway_backend, "enrichr")) {
    .ensure_enrichr_ready(
      site = "Enrichr",
      verbose = TRUE,
      log_fun = log_msg
    )
  }

  dt <- data.table::as.data.table(link_scores)
  if (!nrow(dt)) {
    msg <- "Skipping link-score pathway enrichment: link_scores empty."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (!("topic_num" %in% names(dt))) {
    if (!("topic" %in% names(dt))) {
      .log_abort("link_scores must have topic_num or topic.")
    }
    dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
  }
  dt <- dt[is.finite(topic_num)]
  if (!nrow(dt)) {
    msg <- "Skipping link-score pathway enrichment: no valid topic_num."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  tf_dt_all <- NULL
  if (!is.null(tf_link_scores)) {
    tf_dt_all <- data.table::as.data.table(tf_link_scores)
    if (nrow(tf_dt_all)) {
      if (!("topic_num" %in% names(tf_dt_all))) {
        if (!("topic" %in% names(tf_dt_all))) {
          .log_abort("tf_link_scores must have topic_num or topic.")
        }
        tf_dt_all[, topic_num := as.integer(gsub("^Topic", "", topic))]
      }
      if (!("prob" %in% names(tf_dt_all))) tf_dt_all[, prob := 1]
      tf_dt_all <- tf_dt_all[is.finite(topic_num)]
    }
    if (!nrow(tf_dt_all)) tf_dt_all <- NULL
  }

  .collect_enrichr <- function(gene_sets, log_fun = log_msg) {
    if (!length(gene_sets)) return(data.table::data.table())
    run_one_enrichr_species <- function(nm, species, dbs_use, cache_dir_use) {
      genes <- gene_sets[[nm]]
      if (length(genes) < as.integer(min_genes)) {
        log_fun(sprintf("Topic %s %s skipped: gene count < %d", nm, species, min_genes))
        return(NULL)
      }
      if (!length(dbs_use)) {
        log_fun(sprintf("Topic %s %s skipped: no databases.", nm, species))
        return(NULL)
      }
      log_fun(sprintf("Topic %s %s gene count: %d", nm, species, length(genes)))
      enr <- tryCatch(
        .quiet_enrichr_call(.run_enrichr_cached(
          genes = genes,
          dbs = dbs_use,
          sleep_time = enrichr_sleep_time,
          cache_dir = cache_dir_use,
          backend = pathway_backend,
          pathway_species = species
        )),
        error = function(e) {
          log_fun(sprintf("Topic %s %s enrichr error: %s", nm, species, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(enr)) return(NULL)
      out <- .topic_enrichr_result_to_table(enr, topic_name = nm, genes = genes, pathway_species = species)
      log_fun(sprintf("Topic %s %s enriched pathways: %d finite adjusted p-values", nm, species, nrow(out)))
      out
    }
    run_one_enrichr <- function(nm) {
      if (human_mouse_best) {
        out <- lapply(c("human", "mouse"), function(species) {
          run_one_enrichr_species(
            nm = nm,
            species = species,
            dbs_use = dbs_by_species[[species]],
            cache_dir_use = file.path(enrichr_cache_dir, species)
          )
        })
        out <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
        if (nrow(out)) out <- .select_best_human_mouse_pathways(out)
        return(out)
      }
      run_one_enrichr_species(
        nm = nm,
        species = pathway_species_mode,
        dbs_use = dbs,
        cache_dir_use = enrichr_cache_dir
      )
    }
    set_names <- names(gene_sets)
    if (.Platform$OS.type != "windows" && enrichr_n_cores > 1L && length(set_names) > 1L) {
      res_list <- parallel::mclapply(
        set_names,
        run_one_enrichr,
        mc.cores = min(enrichr_n_cores, length(set_names))
      )
      names(res_list) <- set_names
    } else {
      res_list <- vector("list", length(gene_sets))
      names(res_list) <- set_names
      for (nm in set_names) {
        res_list[[nm]] <- run_one_enrichr(nm)
      }
    }
    data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
  }

  .select_topic_dotplot_terms <- function(res_dt, top_n = dot_top_n_per_topic) {
    dt_sel <- data.table::copy(res_dt)
    dt_sel <- dt_sel[is.finite(padj)]
    if (!nrow(dt_sel)) return(dt_sel)
    dt_sel[, topic := as.integer(topic)]
    dt_sel <- dt_sel[is.finite(topic)]
    if (!nrow(dt_sel)) return(dt_sel)
    if (!("overlap_hits" %in% names(dt_sel))) {
      dt_sel[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
    }
    dt_sel[is.na(overlap_hits), overlap_hits := 0L]
    dt_sel[, logp := -log10(pmax(padj, 1e-300))]
    dt_sel[, pathway_key := pathway]
    if (is.finite(top_n) && as.numeric(top_n) > 0) {
      top_n <- max(1L, as.integer(top_n))
      top_dt <- dt_sel[, {
        ranked <- .SD[order(padj, -overlap_hits, pathway)]
        sig <- ranked[padj <= as.numeric(padj_cut)]
        if (nrow(sig)) head(sig, top_n) else head(ranked, top_n)
      }, by = topic]
      if (!nrow(top_dt)) return(top_dt)
      selected_keys <- unique(top_dt$pathway_key)
      out <- dt_sel[pathway_key %in% selected_keys]
      sig_keys <- unique(out[padj <= as.numeric(padj_cut), pathway_key])
      out <- out[pathway_key %in% sig_keys]
    } else {
      top_dt <- dt_sel[padj <= as.numeric(padj_cut)]
      if (!nrow(top_dt)) return(top_dt)
      out <- data.table::copy(top_dt)
    }
    if (!nrow(out)) return(out)
    top_key <- unique(top_dt[, .(topic, pathway_key)])
    top_key[, is_topic_top_pathway := TRUE]
    out <- merge(out, top_key, by = c("topic", "pathway_key"), all.x = TRUE, sort = FALSE)
    out[is.na(is_topic_top_pathway), is_topic_top_pathway := FALSE]
    topic_levels <- sort(unique(out$topic))
    row_levels <- character(0)
    for (tp in topic_levels) {
      owned <- out[topic == tp & is_topic_top_pathway == TRUE]
      if (!nrow(owned)) next
      data.table::setorder(owned, -logp, -overlap_hits, pathway)
      for (key in owned$pathway_key) {
        if (!(key %in% row_levels)) row_levels <- c(row_levels, key)
      }
    }
    remaining <- setdiff(unique(out$pathway_key), row_levels)
    if (length(remaining)) {
      remain_rank <- out[pathway_key %in% remaining, .(
        max_logp = max(logp, na.rm = TRUE),
        max_overlap = max(overlap_hits, na.rm = TRUE)
      ), by = pathway_key]
      data.table::setorder(remain_rank, -max_logp, -max_overlap, pathway_key)
      row_levels <- c(row_levels, remain_rank$pathway_key)
    }
    if (nrow(out) && is.finite(max_pathways) && as.numeric(max_pathways) > 0 && length(row_levels) > as.integer(max_pathways)) {
      row_levels <- row_levels[seq_len(as.integer(max_pathways))]
      out <- out[pathway_key %in% row_levels]
    }
    out[, pathway_key := factor(pathway_key, levels = rev(row_levels))]
    data.table::setorder(out, topic, -is_topic_top_pathway, padj, -overlap_hits, pathway)
    out[]
  }

  .write_dotplot <- function(res_dt, dot_prefix, plot_title, log_fun = log_msg) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      log_fun("Skipping dot plot: ggplot2 not installed.")
      return(invisible(NULL))
    }
    plot_dt <- .select_topic_dotplot_terms(res_dt)
    if (!nrow(plot_dt)) {
      log_fun("Skipping dot plot: no finite pathway terms.")
      return(invisible(NULL))
    }
    plot_dt[, topic_num := as.integer(topic)]
    topic_levels <- sort(unique(plot_dt$topic_num))
    plot_dt[, topic_label := factor(paste0("Topic", topic_num), levels = paste0("Topic", topic_levels))]
    plot_dt[, score_val := data.table::fifelse(is.finite(combined_score), combined_score, NA_real_)]
    plot_dt[!is.finite(score_val) & is.finite(odds_ratio), score_val := odds_ratio]
    plot_dt[!is.finite(score_val) & is.finite(logp), score_val := logp]
    plot_dt[, score_val := pmax(score_val, 0)]
    size_cap <- suppressWarnings(stats::quantile(plot_dt$score_val[is.finite(plot_dt$score_val)], probs = 0.95, na.rm = TRUE, names = FALSE))
    if (is.finite(size_cap)) plot_dt[, score_val := pmin(score_val, size_cap)]
    plot_dt[, is_sig := is.finite(logp) & logp >= 1.3]

    dot_path <- paste0(dot_prefix, ".pdf")
    dot_csv <- paste0(dot_prefix, ".csv")
    ord_dt <- plot_dt[order(
      match(pathway_key, levels(pathway_key)),
      match(topic_label, levels(topic_label))
    )]
    data.table::fwrite(ord_dt, dot_csv)

    family <- if (exists(".diff_grn_pathway_font_family", mode = "function")) {
      .diff_grn_pathway_font_family()
    } else {
      "sans"
    }
    n_topics_plot <- data.table::uniqueN(plot_dt$topic_label)
    n_paths_plot <- data.table::uniqueN(plot_dt$pathway_key)
    max_val <- max(plot_dt$logp[is.finite(plot_dt$logp)], 1.300001, na.rm = TRUE)
    color_values <- if (exists(".diff_grn_pathway_rescale", mode = "function")) {
      .diff_grn_pathway_rescale(c(0, 1.3, seq(1.3 + 1e-6, max_val, length.out = 6)))
    } else {
      color_breaks <- c(0, 1.3, seq(1.3 + 1e-6, max_val, length.out = 6))
      color_rng <- range(color_breaks, finite = TRUE)
      if (!all(is.finite(color_rng)) || diff(color_rng) == 0) {
        rep(0, length(color_breaks))
      } else {
        (color_breaks - color_rng[1]) / diff(color_rng)
      }
    }
    p <- ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(x = topic_label, y = pathway_key, size = score_val)
    ) +
      ggplot2::geom_point(ggplot2::aes(fill = logp), shape = 21, color = "transparent", alpha = 0.85) +
      ggplot2::geom_point(
        data = plot_dt[is_sig == TRUE],
        ggplot2::aes(fill = logp, color = logp),
        shape = 21,
        stroke = 0.7,
        alpha = 0.95
      ) +
      ggplot2::scale_fill_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        name = "-log10 adjusted p-value"
      ) +
      ggplot2::scale_color_gradientn(
        colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
        values = color_values,
        limits = c(0, max_val),
        oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
        guide = "none"
      ) +
      ggplot2::scale_size_continuous(name = "Combined score", range = c(0.8, 3.2)) +
      ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
      ggplot2::scale_y_discrete(labels = function(x) {
        y <- as.character(x)
        y <- gsub("\\s+", " ", y)
        too_long <- nchar(y) > 58L
        y[too_long] <- paste0(substr(y[too_long], 1L, 55L), "...")
        y
      }) +
      ggplot2::labs(
        x = "Topic",
        y = "Pathway term",
        title = plot_title,
        caption = if (is.finite(dot_top_n_per_topic) && as.numeric(dot_top_n_per_topic) > 0) {
          sprintf("Union of top %d pathways per topic; all occurrences shown.", as.integer(dot_top_n_per_topic))
        } else {
          sprintf("All significant pathways per topic; adjusted p-value <= %s.", as.character(padj_cut))
        }
      ) +
      ggplot2::theme_bw(base_size = 9, base_family = family) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", hjust = 0.5),
        plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
        legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
      )
    ok_plot <- TRUE
    tryCatch(
      {
        if (exists(".diff_grn_pathway_save_fixed_panel_pdf", mode = "function")) {
          .diff_grn_pathway_save_fixed_panel_pdf(
            p,
            dot_path,
            panel_width = max(1.2, 0.11 * n_topics_plot),
            panel_height = max(1.1, 0.17 * n_paths_plot),
            family = family
          )
        } else {
          ggplot2::ggsave(dot_path, p, width = max(8, n_topics_plot * 0.25), height = max(6, n_paths_plot * 0.2), limitsize = FALSE)
        }
      },
      error = function(e) {
        ok_plot <<- FALSE
        log_fun(sprintf("Dot plot failed: %s", conditionMessage(e)))
      }
    )
    log_fun(sprintf("Dot plot CSV saved: %s", dot_csv))
    if (ok_plot) {
      log_fun(sprintf("Dot plot PDF saved: %s", dot_path))
    }
    invisible(TRUE)
  }

  gene_sets <- link_scores_to_gene_sets(
    dt,
    include_tf = include_tf,
    include_gene = include_gene,
    min_prob = min_prob,
    gene_min_prob = gene_min_prob,
    tf_min_prob = tf_min_prob,
    tf_max_topics = tf_max_topics,
    tf_top_n_per_topic = tf_top_n_per_topic,
    tf_link_scores = tf_dt_all,
    gene_terms = gene_terms
  )
  if (!length(gene_sets)) {
    msg <- "Skipping link-score pathway enrichment: no gene sets after filtering."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  res_dt <- .collect_enrichr(gene_sets, log_fun = log_msg)
  if (!nrow(res_dt)) {
    msg <- "Skipping link-score pathway enrichment: no enriched terms at padj_cut."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }

  res_dt <- res_dt[is.finite(logp) & logp > 0]
  if (!nrow(res_dt)) {
    msg <- "Skipping link-score pathway enrichment: all logp values non-finite or zero."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) {
    res_dt <- res_dt[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
  } else {
    res_dt <- res_dt[order(topic, -logp)]
  }
  if (nrow(res_dt) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
    path_rank <- res_dt[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
    if (nrow(path_rank) > as.integer(max_pathways)) {
      keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
      res_dt <- res_dt[pathway %in% keep]
      log_msg(sprintf("Filtered pathways to top %d by max logp.", as.integer(max_pathways)))
    }
  }
  res_dt[, topic := as.integer(topic)]
  data.table::setorder(res_dt, topic)

  if (isTRUE(make_heatmap)) {
    log_msg("Skipping heatmap: pathway heatmap output has been removed from standard Module 3 reports.")
  }

  if (isTRUE(make_dotplot)) {
    dot_base <- if (!is.null(file_tag) && nzchar(file_tag)) {
      paste0("topic_pathway_enrichment_", file_tag, "_dotplot")
    } else {
      "topic_pathway_enrichment_dotplot"
    }
    .write_dotplot(
      res_dt,
      dot_prefix = file.path(out_dir, dot_base),
      plot_title = main_title
    )
  } else {
    log_msg("Skipping dot plot: make_dotplot = FALSE.")
  }

  if (isTRUE(per_comparison)) {
    doc_info <- .parse_doc_id(dt$doc_id, doc_design = doc_design)
    dt <- cbind(dt, doc_info)
    dt[, direction_group := if (isTRUE(split_direction) && identical(doc_design, "comparison")) .direction_group(direction) else "All"]
    dt <- dt[!is.na(comparison_id) & nzchar(comparison_id)]
    if (!is.null(tf_dt_all) && nrow(tf_dt_all) && "doc_id" %in% names(tf_dt_all)) {
      tf_info <- .parse_doc_id(tf_dt_all$doc_id, doc_design = doc_design)
      tf_dt_all <- cbind(tf_dt_all, tf_info)
      tf_dt_all[, direction_group := if (isTRUE(split_direction) && identical(doc_design, "comparison")) .direction_group(direction) else "All"]
      tf_dt_all <- tf_dt_all[!is.na(comparison_id) & nzchar(comparison_id)]
      if (!nrow(tf_dt_all)) tf_dt_all <- NULL
    }
    if (nrow(dt)) {
      out_dir_pc <- file.path(out_dir, per_comparison_dir)
      dir.create(out_dir_pc, recursive = TRUE, showWarnings = FALSE)
      log_msg(sprintf("Per-comparison pathway: %d comparisons", length(unique(dt$comparison_id))))
      for (cmp in unique(dt$comparison_id)) {
        cmp_dt <- dt[comparison_id == cmp]
        if (!nrow(cmp_dt)) next
        cmp_dir <- if (isTRUE(per_comparison_flat)) {
          out_dir_pc
        } else {
          file.path(out_dir_pc, .safe_filename(cmp))
        }
        if (!isTRUE(per_comparison_flat)) {
          dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)
        }
        for (dir_lab in unique(cmp_dt$direction_group)) {
          sub_dt <- cmp_dt[direction_group == dir_lab]
          if (!nrow(sub_dt)) next
          log_msg(sprintf("Per-comparison %s | %s: links=%d", cmp, dir_lab, nrow(sub_dt)))
          log_local <- function(msg) log_msg(sprintf("%s | %s: %s", cmp, dir_lab, msg))
          prefix <- file.path(cmp_dir, paste0(.safe_filename(cmp), "_", .safe_filename(dir_lab), "_dotplot"))
          if (!overwrite && isTRUE(make_dotplot)) {
            have_csv <- file.exists(paste0(prefix, ".csv"))
            have_pdf <- file.exists(paste0(prefix, ".pdf"))
            if (isTRUE(have_csv) && isTRUE(have_pdf)) {
              log_local("Existing dot plot outputs found; skipping.")
              next
            }
          }
          tf_sub <- NULL
          if (!is.null(tf_dt_all) && nrow(tf_dt_all)) {
            if (all(c("comparison_id", "direction_group") %in% names(tf_dt_all))) {
              tf_sub <- tf_dt_all[comparison_id == cmp & direction_group == dir_lab]
            } else {
              tf_sub <- tf_dt_all
            }
          }
          gs <- link_scores_to_gene_sets(
            sub_dt,
            include_tf = include_tf,
            include_gene = include_gene,
            min_prob = min_prob,
            gene_min_prob = gene_min_prob,
            tf_min_prob = tf_min_prob,
            tf_max_topics = tf_max_topics,
            tf_top_n_per_topic = tf_top_n_per_topic,
            tf_link_scores = tf_sub,
            gene_terms = gene_terms
          )
          if (!length(gs)) {
            log_local("No gene sets after filtering; skipping.")
            next
          }
          res_sub <- .collect_enrichr(gs, log_fun = log_local)
          if (!nrow(res_sub)) {
            log_local("No enriched terms at padj_cut; skipping.")
            next
          }
          res_sub <- res_sub[is.finite(logp) & logp > 0]
          if (!nrow(res_sub)) next
          if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) {
            res_sub <- res_sub[order(-logp), .SD[seq_len(min(.N, as.integer(top_n_per_topic)))], by = topic]
          } else {
            res_sub <- res_sub[order(topic, -logp)]
          }
          if (nrow(res_sub) && is.finite(max_pathways) && as.numeric(max_pathways) > 0) {
            path_rank <- res_sub[, .(max_logp = max(logp, na.rm = TRUE)), by = pathway]
            if (nrow(path_rank) > as.integer(max_pathways)) {
              keep <- path_rank[order(-max_logp)][seq_len(as.integer(max_pathways)), pathway]
              res_sub <- res_sub[pathway %in% keep]
            }
          }
          res_sub[, topic := as.integer(topic)]
          data.table::setorder(res_sub, topic)
          plot_title <- paste(cmp, dir_lab, "Pathway enrichment (link scores)", sep = " | ")
          if (isTRUE(make_dotplot)) {
            .write_dotplot(res_sub, dot_prefix = prefix, plot_title = plot_title, log_fun = log_local)
          } else {
            log_local("Skipping dot plot: make_dotplot = FALSE.")
          }
        }
      }
    }
  }

  invisible(TRUE)
}

rerun_pathway_from_link_scores <- function(out_dir,
                                           link_scores_file = NULL,
                                           tf_link_scores_file = NULL,
                                           gene_terms_file = NULL,
                                           allow_missing = FALSE,
                                           ...) {
  resolve_rel <- function(path, base_dir) {
    if (is.null(path) || !nzchar(path)) return(path)
    if (grepl("^/", path)) return(path)
    cand <- file.path(base_dir, path)
    if (file.exists(cand)) return(cand)
    path
  }
  if (is.null(link_scores_file)) {
    cand <- file.path(out_dir, "link_topic_scores_baseline.csv")
    if (file.exists(cand)) {
      link_scores_file <- cand
    } else {
    cand <- file.path(out_dir, "link_topic_scores_gate_peak_and_gene_in_set.csv")
      if (file.exists(cand)) link_scores_file <- cand
    }
  }
  link_scores_file <- resolve_rel(link_scores_file, out_dir)
  if (is.null(link_scores_file) || !file.exists(link_scores_file)) {
    if (isTRUE(allow_missing)) {
      .log_inform("Skipping link-score pathway enrichment: link_scores_file not found.")
      return(invisible(NULL))
    }
    .log_abort("link_scores_file not found for rerun_pathway_from_link_scores.")
  }
  link_scores <- data.table::fread(link_scores_file)
  tf_link_scores <- NULL
  tf_link_scores_file <- resolve_rel(tf_link_scores_file, out_dir)
  if (!is.null(tf_link_scores_file) && file.exists(tf_link_scores_file)) {
    tf_link_scores <- data.table::fread(tf_link_scores_file)
  }
  gene_terms <- NULL
  gene_terms_file <- resolve_rel(gene_terms_file, out_dir)
  if (!is.null(gene_terms_file) && file.exists(gene_terms_file)) {
    gene_terms <- data.table::fread(gene_terms_file)
  }
  dots <- list(...)
  if (!"title_prefix" %in% names(dots)) {
    dots$title_prefix <- basename(out_dir)
  }
  do.call(
    plot_topic_pathway_enrichment_from_link_scores,
    c(list(link_scores = link_scores, tf_link_scores = tf_link_scores, gene_terms = gene_terms, out_dir = out_dir), dots)
  )
}

rerun_pathway_from_topic_links <- function(out_dir,
                                           topic_links_file = NULL,
                                           method = c("peak_and_gene", "peak_and_gene_prob", "gene_only", "link_score_prob", "link_score_efdr"),
                                           allow_missing = FALSE,
                                           ...) {
  resolve_rel <- function(path, base_dir) {
    if (is.null(path) || !nzchar(path)) return(path)
    if (grepl("^/", path)) return(path)
    cand <- file.path(base_dir, path)
    if (file.exists(cand)) return(cand)
    path
  }
  method <- match.arg(method)
  if (is.null(topic_links_file)) {
    cand <- .topic_links_path(out_dir, prefer = "pass")
    if (file.exists(cand)) topic_links_file <- cand
  }
  topic_links_file <- resolve_rel(topic_links_file, out_dir)
  if (is.null(topic_links_file) || !file.exists(topic_links_file)) {
    fallback <- .topic_links_path(out_dir, prefer = "pass")
    if (!is.null(fallback) && file.exists(fallback)) {
      .log_inform("Topic-link pathway enrichment: using pass-only topic links at {fallback}.")
      topic_links_file <- fallback
    }
  }
  if (is.null(topic_links_file) || !file.exists(topic_links_file)) {
    if (isTRUE(allow_missing)) {
      .log_inform("Skipping topic-link pathway enrichment: topic_links_file not found.")
      return(invisible(NULL))
    }
    .log_abort("topic_links_file not found for rerun_pathway_from_topic_links.")
  }
  topic_links <- data.table::fread(topic_links_file)
  link_scores <- .topic_links_to_link_scores(topic_links, method = method)
  dots <- list(...)
  if (!"title_prefix" %in% names(dots) || is.null(dots$title_prefix)) {
    dots$title_prefix <- basename(out_dir)
  }
  dots$title_prefix <- paste(dots$title_prefix, method, sep = " | ")
  if (!"file_tag" %in% names(dots) || is.null(dots$file_tag) || !nzchar(dots$file_tag)) {
    dots$file_tag <- method
  }
  do.call(
    plot_topic_pathway_enrichment_from_link_scores,
    c(list(link_scores = link_scores, out_dir = out_dir), dots)
  )
}

filter_topic_pathway_dotplot_csv <- function(dotplot_csv,
                                             out_prefix,
                                             padj_cut = 0.001,
                                             title = NULL,
                                             max_pathways = Inf) {
  .assert_pkg("data.table")
  .assert_pkg("ggplot2")
  if (!file.exists(dotplot_csv)) {
    .log_abort("Missing pathway dotplot CSV: {dotplot_csv}")
  }
  dt <- data.table::fread(dotplot_csv, showProgress = FALSE)
  req <- c("topic", "pathway_key", "pathway", "padj", "logp", "score_val")
  .assert_has_cols(dt, req, context = "filter_topic_pathway_dotplot_csv")
  dt <- dt[is.finite(padj) & padj <= as.numeric(padj_cut)]
  if (!nrow(dt)) {
    data.table::fwrite(dt, paste0(out_prefix, ".csv"))
    grDevices::pdf(paste0(out_prefix, ".pdf"), width = 8, height = 5, family = "Helvetica")
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::plot.new()
    graphics::text(0.5, 0.55, if (is.null(title)) "Pathway enrichment" else title, font = 2)
    graphics::text(0.5, 0.42, sprintf("No pathway terms with adjusted p-value <= %s.", as.character(padj_cut)))
    return(invisible(FALSE))
  }
  dt[, topic := as.integer(topic)]
  dt <- dt[is.finite(topic)]
  if (!("overlap_hits" %in% names(dt))) {
    dt[, overlap_hits := suppressWarnings(as.integer(sub("/.*$", "", as.character(overlap))))]
  }
  dt[is.na(overlap_hits), overlap_hits := 0L]
  dt[, logp := -log10(pmax(padj, 1e-300))]
  dt[, pathway_key := as.character(pathway_key)]
  topic_levels <- sort(unique(dt$topic))
  row_levels <- character(0)
  for (tp in topic_levels) {
    owned <- dt[topic == tp]
    if (!nrow(owned)) next
    data.table::setorder(owned, -logp, -overlap_hits, pathway)
    for (key in owned$pathway_key) {
      if (!(key %in% row_levels)) row_levels <- c(row_levels, key)
    }
  }
  remaining <- setdiff(unique(dt$pathway_key), row_levels)
  if (length(remaining)) {
    remain_rank <- dt[pathway_key %in% remaining, .(
      max_logp = max(logp, na.rm = TRUE),
      max_overlap = max(overlap_hits, na.rm = TRUE)
    ), by = pathway_key]
    data.table::setorder(remain_rank, -max_logp, -max_overlap, pathway_key)
    row_levels <- c(row_levels, remain_rank$pathway_key)
  }
  if (is.finite(max_pathways) && as.numeric(max_pathways) > 0 && length(row_levels) > as.integer(max_pathways)) {
    row_levels <- row_levels[seq_len(as.integer(max_pathways))]
    dt <- dt[pathway_key %in% row_levels]
  }
  dt[, pathway_key := factor(pathway_key, levels = rev(row_levels))]
  dt[, topic_label := factor(paste0("Topic", topic), levels = paste0("Topic", topic_levels))]
  if (!("score_val" %in% names(dt))) {
    dt[, score_val := data.table::fifelse(is.finite(combined_score), combined_score, logp)]
  }
  dt[, score_val := pmax(as.numeric(score_val), 0)]
  dt[, pathway_order := match(as.character(pathway_key), row_levels)]
  data.table::setorder(dt, pathway_order, topic, padj, -overlap_hits, pathway)
  dt[, pathway_order := NULL]
  data.table::fwrite(dt, paste0(out_prefix, ".csv"))

  family <- if (exists(".diff_grn_pathway_font_family", mode = "function")) {
    .diff_grn_pathway_font_family()
  } else {
    "sans"
  }
  max_val <- max(dt$logp[is.finite(dt$logp)], -log10(as.numeric(padj_cut)), na.rm = TRUE)
  color_values <- if (exists(".diff_grn_pathway_rescale", mode = "function")) {
    .diff_grn_pathway_rescale(c(0, -log10(as.numeric(padj_cut)), seq(-log10(as.numeric(padj_cut)) + 1e-6, max_val, length.out = 6)))
  } else {
    color_breaks <- c(0, -log10(as.numeric(padj_cut)), seq(-log10(as.numeric(padj_cut)) + 1e-6, max_val, length.out = 6))
    color_rng <- range(color_breaks, finite = TRUE)
    if (!all(is.finite(color_rng)) || diff(color_rng) == 0) {
      rep(0, length(color_breaks))
    } else {
      (color_breaks - color_rng[1]) / diff(color_rng)
    }
  }
  n_topics_plot <- data.table::uniqueN(dt$topic_label)
  n_paths_plot <- data.table::uniqueN(dt$pathway_key)
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = topic_label, y = pathway_key, size = score_val)) +
    ggplot2::geom_point(ggplot2::aes(fill = logp, color = logp), shape = 21, stroke = 0.7, alpha = 0.95) +
    ggplot2::scale_fill_gradientn(
      colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
      values = color_values,
      limits = c(0, max_val),
      oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
      name = "-log10 adjusted p-value"
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("grey85", "grey85", grDevices::hcl.colors(6, "Plasma")),
      values = color_values,
      limits = c(0, max_val),
      oob = function(x, range = c(0, 1)) pmin(pmax(x, range[1]), range[2]),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(name = "Combined score", range = c(0.8, 3.2)) +
    ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
    ggplot2::scale_y_discrete(labels = function(x) {
      y <- as.character(x)
      y <- gsub("\\s+", " ", y)
      too_long <- nchar(y) > 58L
      y[too_long] <- paste0(substr(y[too_long], 1L, 55L), "...")
      y
    }) +
    ggplot2::labs(
      x = "Topic",
      y = "Pathway term",
      title = if (is.null(title)) "Pathway enrichment" else title,
      caption = sprintf("All significant pathways per topic; adjusted p-value <= %s.", as.character(padj_cut))
    ) +
    ggplot2::theme_bw(base_size = 9, base_family = family) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      plot.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", hjust = 0.5),
      plot.caption = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      axis.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      axis.text.x = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black", angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      legend.title = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black"),
      legend.text = ggplot2::element_text(size = 9, face = "bold", family = family, color = "black")
    )
  if (exists(".diff_grn_pathway_save_fixed_panel_pdf", mode = "function")) {
    .diff_grn_pathway_save_fixed_panel_pdf(
      p,
      paste0(out_prefix, ".pdf"),
      panel_width = max(1.2, 0.11 * n_topics_plot),
      panel_height = max(1.1, 0.17 * n_paths_plot),
      family = family
    )
  } else {
    ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = max(8, n_topics_plot * 0.25), height = max(6, n_paths_plot * 0.2), limitsize = FALSE)
  }
  invisible(TRUE)
}

plot_intertopic_distance_map <- function(phi, topic_terms, out_file, option_label, title_prefix = NULL) {
  phi <- as.matrix(phi)
  if (nrow(phi) < 2L) return(invisible(NULL))

  topic_names <- rownames(phi)
  if (is.null(topic_names)) topic_names <- paste0("Topic", seq_len(nrow(phi)))
  rownames(phi) <- topic_names

  phi[!is.finite(phi)] <- 0
  rs <- rowSums(phi)
  rs[!is.finite(rs) | rs == 0] <- 1
  phi_prob <- phi / rs

  kl_div <- function(p, q) {
    idx <- (p > 0) & (q > 0)
    sum(p[idx] * log(p[idx] / q[idx]))
  }

  K <- nrow(phi_prob)
  dmat <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i >= j) next
      p <- phi_prob[i, ]
      q <- phi_prob[j, ]
      m <- 0.5 * (p + q)
      jsd <- 0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
      if (!is.finite(jsd)) jsd <- 0
      dmat[i, j] <- jsd
      dmat[j, i] <- jsd
    }
  }
  diag(dmat) <- 0
  if (!any(is.finite(dmat))) {
    .log_inform("Skipping intertopic distance map: no finite distances.")
    return(invisible(NULL))
  }
  if (any(!is.finite(dmat))) {
    max_finite <- max(dmat[is.finite(dmat)], na.rm = TRUE)
    dmat[!is.finite(dmat)] <- max_finite
  }
  if (max(dmat, na.rm = TRUE) == 0) {
    .log_inform("Skipping intertopic distance map: all distances are zero.")
    return(invisible(NULL))
  }

  coords <- stats::cmdscale(stats::as.dist(dmat), k = 2)
  if (is.null(coords)) return(invisible(NULL))

  sizes <- rep(1, K)
  if (!is.null(topic_terms) && nrow(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    tt <- tt[.as_logical_flag(in_topic)]
    if (option_label == "joint") {
      counts <- tt[, .N, by = topic]
    } else {
      term_prefix <- if (grepl("opt3", option_label)) "GENE:" else "PEAK:"
      tt <- tt[grepl(paste0("^", term_prefix), term_id)]
      counts <- tt[, .N, by = topic]
    }
    topic_ids <- .m3_opt_topic_ids(phi, "row")
    count_ids <- suppressWarnings(as.integer(sub("^Topic", "", counts$topic)))
    sizes <- counts$N[match(topic_ids, count_ids)]
    sizes[!is.finite(sizes)] <- 1
  }
  sizes <- sqrt(pmax(sizes, 1))
  sizes <- sizes / max(sizes, na.rm = TRUE)
  sizes <- 0.6 + 2.2 * sizes

  grDevices::pdf(out_file, width = 7.5, height = 6.2)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); grDevices::dev.off() }, add = TRUE)
  x <- coords[, 1]
  y <- coords[, 2]
  x_pad <- diff(range(x)) * 0.1
  y_pad <- diff(range(y)) * 0.1
  main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste("Intertopic Distance Map (MDS)", title_prefix, sep = " | ")
  } else {
    "Intertopic Distance Map (MDS)"
  }
  plot(x, y,
       type = "n",
       xlab = "Dim 1",
       ylab = "Dim 2",
       main = main_title,
       xlim = range(x) + c(-x_pad, x_pad),
       ylim = range(y) + c(-y_pad, y_pad),
       asp = 1)
  abline(h = 0, v = 0, col = "grey85", lty = "dashed")
  symbols(x, y,
          circles = sizes,
          inches = 0.3,
          bg = grDevices::adjustcolor("#5ab4e6", alpha.f = 0.6),
          fg = "grey30",
          add = TRUE)
  topic_nums <- gsub("^Topic", "", topic_names)
  text(x, y, labels = topic_nums, cex = 0.8, font = 2)
  label_text <- if (option_label == "joint") {
    "#Features (topic)"
  } else if (grepl("opt3", option_label)) {
    "#Target genes (topic)"
  } else {
    "#FP peaks (topic)"
  }
  mtext(label_text, side = 3, line = 0.2, cex = 0.85)
  invisible(TRUE)
}

# =============================================================================
# 9) Report-only runner (use precomputed theta/phi)
# =============================================================================

run_tfdocs_report_from_topic_base <- function(topic_base,
                                              dtm,
                                              edges_docs,
                                              out_dir,
                                              option_label = c("opt1_peak_delta_fp", "opt2_peak_fc_fp",
                                                               "opt3_gene_fc_expr", "joint"),
                                              direction_by = c("gene", "fp", "none"),
                                              topic_score_method = c("normtop_specificity", "rowmax_phi"),
                                              topic_term_assignment_method = c("max_phi", "gammafit_maxprob", "gammafit"),
                                              topic_model_family = NULL,
                                              binarize_method = c("gammafit", "topn"),
                                              gammafit_scope = c("topic_term_group", "global_term_group"),
                                              thrP = 0.975,
                                              top_n_terms = 500L,
                                              in_topic_min_terms = 50L,
                                              pathway_use_all_terms = FALSE,
                                              pathway_make_heatmap = FALSE,
                                              pathway_make_dotplot = TRUE,
                                              pathway_overwrite = FALSE,
                                              top_n_per_topic = 20L,
                                              dot_top_n_per_topic = 25L,
                                              max_pathways = 200L,
                                              pathway_tf_link_mode = c("theta", "none"),
                                              pathway_tf_top_n_docs = 50L,
                                              pathway_tf_min_theta = NA_real_,
                                              run_pathway_gsea = FALSE,
                                              gsea_species = "Homo sapiens",
                                              gsea_nperm = 1000L,
                                              gsea_peak_gene_agg = c("max", "sum"),
                                              pathway_source = c("topic_terms", "link_scores"),
                                              pathway_link_scores_file = NULL,
                                              pathway_link_scores_file_tf = NULL,
                                              pathway_link_gene_terms_file = NULL,
                                              pathway_link_min_prob = 0,
                                              pathway_link_include_tf = TRUE,
                                              pathway_link_include_gene = TRUE,
                                              pathway_link_gene_min_prob = NULL,
                                              pathway_link_tf_min_prob = NULL,
                                              pathway_link_tf_max_topics = Inf,
                                              pathway_link_tf_top_n_per_topic = NA_integer_,
                                              pathway_per_comparison = FALSE,
                                              pathway_per_comparison_dir = ".",
                                              pathway_per_comparison_flat = TRUE,
                                              pathway_split_direction = TRUE,
                                              pathway_enrichr_sleep_time = 0,
                                              pathway_enrichr_cache_dir = NULL,
                                              pathway_enrichr_n_cores = NULL,
                                              pathway_backend = NULL,
                                              pathway_species = NULL,
                                              pathway_databases = NULL,
                                              run_link_topic_scores = FALSE,
                                              fp_term_mode = c("unique", "aggregate", "aggregate_weight", "gene_expression"),
                                              link_topic_gate_mode = "none",
                                              link_topic_top_k = 3L,
                                              link_topic_min_prob = 0,
                                              link_topic_include_tf = FALSE,
                                              link_topic_chunk_size = 5000L,
                                              link_topic_n_cores = 1L,
                                              link_topic_overwrite = FALSE,
                                              link_topic_method = c("gammafit", "theta_and_terms", "gene_prob", "link_score_prob", "link_score_efdr"),
                                              link_topic_prob_cutoff = 0.3,
                                              link_topic_fdr_q = 0.2,
                                              link_topic_fdr_p = NA_real_,
                                              link_topic_efdr_scope = c("per_topic", "global"),
                                              link_topic_efdr_B = 100L,
                                              link_topic_efdr_seed = 1L,
                                              link_topic_report_methods = NULL,
                                              link_topic_output = c("pass", "full", "both", "none"),
                                              run_gammafit_summary = TRUE,
                                              run_link_efdr_summary = TRUE,
                                              run_pathway_enrichment = TRUE,
                                              run_tf_topic_assignment = TRUE,
                                              topic_tf_membership_cutoff = 0.3,
                                              topic_tf_primary_margin_cutoff = 0.1,
                                              run_raw_theta_document_heatmap = FALSE,
                                              run_document_theta_umap = TRUE,
                                              theta_umap_selected_tfs = NULL,
                                              theta_umap_top_n_tfs = 12L,
                                              theta_umap_seed = 123L,
                                              theta_umap_n_neighbors = 30L,
                                              theta_umap_condition_colors = NULL,
                                              run_topic_term_heatmap = TRUE,
                                              run_topic_by_comparison_heatmaps = TRUE,
                                              run_intertopic_distance_map = TRUE,
                                              optimize_topics = NULL,
                                              topic_merge_min_genes = 50L,
                                              topic_merge_min_links = 200L,
                                              topic_merge_similarity_threshold = 0.90,
                                              run_topic_assignment_qc = NULL,
                                              topic_qc_umap_links_per_condition = 3000L,
                                              topic_qc_top_tfs = 150L,
                                              topic_qc_condition_expression_file = NULL,
                                              topic_qc_reference_condition = NULL,
                                              topic_qc_upregulated_log2fc_min = 1,
                                              topic_qc_upregulated_pseudocount = 1,
                                              topic_qc_seed = 20260716L,
                                              extraction_steps = NULL,
                                              topic_by_comparison_label_cleaner = NULL,
                                              doc_design = c("comparison", "condition"),
                                              title_prefix = NULL) {
  option_label <- match.arg(option_label)
  doc_design <- match.arg(doc_design)
  direction_by <- match.arg(direction_by)
  topic_score_method <- match.arg(topic_score_method)
  topic_term_assignment_method <- match.arg(topic_term_assignment_method)
  binarize_method <- match.arg(binarize_method)
  gammafit_scope <- match.arg(gammafit_scope)
  pathway_tf_link_mode <- match.arg(pathway_tf_link_mode)
  gsea_peak_gene_agg <- match.arg(gsea_peak_gene_agg)
  pathway_source <- match.arg(pathway_source)
  pathway_enrichr_sleep_time <- .normalize_enrichr_sleep_time(pathway_enrichr_sleep_time)
  pathway_enrichr_n_cores <- .normalize_enrichr_n_cores(pathway_enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  if (is.null(pathway_enrichr_cache_dir)) {
    pathway_enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir, backend = pathway_backend)
  }
  pathway_species <- .normalize_pathway_species_mode(pathway_species)
  pathway_dbs <- if (is.null(pathway_databases)) {
    .default_pathway_databases(pathway_species)
  } else {
    unique(as.character(pathway_databases))
  }
  pathway_dbs <- pathway_dbs[!is.na(pathway_dbs) & nzchar(pathway_dbs)]
  if (!length(pathway_dbs)) .log_abort("pathway_databases must contain at least one database.")
  pathway_topic_term_theta_min <- .safe_num(pathway_tf_min_theta)
  if (!is.finite(pathway_topic_term_theta_min)) {
    pathway_topic_term_theta_min <- .safe_num(topic_tf_membership_cutoff)
  }
  if (!is.finite(pathway_topic_term_theta_min)) pathway_topic_term_theta_min <- 0.3
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  if (identical(topic_term_assignment_method, "gammafit_maxprob") &&
      !fp_term_mode %in% c("aggregate", "gene_expression")) {
    .log_abort(
      "gammafit_maxprob requires paired aggregate terms or gene_expression terms."
    )
  }
  allowed_gate_modes <- c("none", "peak_in_set", "gene_in_set", "peak_and_gene_in_set")
  link_topic_gate_mode <- unique(as.character(link_topic_gate_mode))
  link_topic_method <- match.arg(link_topic_method)
  link_topic_output <- match.arg(link_topic_output)
  link_topic_efdr_scope <- match.arg(link_topic_efdr_scope)
  extraction_steps <- .normalize_topic_extraction_steps(extraction_steps)
  if (!is.null(extraction_steps)) {
    run_link_topic_scores <- .topic_step_enabled(extraction_steps, "topic_links", TRUE)
    run_gammafit_summary <- .topic_step_enabled(extraction_steps, "gammafit_summary", TRUE)
    run_link_efdr_summary <- .topic_step_enabled(extraction_steps, "link_efdr_summary", TRUE)
    run_pathway_enrichment <- .topic_step_enabled(extraction_steps, "pathway", TRUE)
    run_tf_topic_assignment <- .topic_step_enabled(extraction_steps, "tf_topic_assignment", TRUE)
    run_raw_theta_document_heatmap <- .topic_step_enabled(extraction_steps, "raw_theta_documents", FALSE)
    run_document_theta_umap <- .topic_step_enabled(extraction_steps, "document_theta_umap", TRUE)
    run_topic_term_heatmap <- .topic_step_enabled(extraction_steps, "topic_term_heatmap", TRUE)
    run_topic_by_comparison_heatmaps <- .topic_step_enabled(extraction_steps, "topic_by_comparison", TRUE)
    run_intertopic_distance_map <- .topic_step_enabled(extraction_steps, "intertopic_distance", TRUE)
  }
  if (!all(link_topic_gate_mode %in% allowed_gate_modes)) {
    .log_abort("link_topic_gate_mode must be one of: {paste(allowed_gate_modes, collapse = ', ')}.")
  }
  optimization_eligible <- identical(doc_design, "condition") &&
    fp_term_mode %in% c("aggregate", "gene_expression") &&
    identical(topic_term_assignment_method, "gammafit_maxprob")
  optimize_topics <- if (is.null(optimize_topics)) {
    optimization_eligible
  } else {
    isTRUE(optimize_topics)
  }
  run_topic_assignment_qc <- if (is.null(run_topic_assignment_qc)) {
    optimize_topics
  } else {
    isTRUE(run_topic_assignment_qc)
  }
  if (optimize_topics && !optimization_eligible) {
    .log_abort(
      paste(
        "Topic optimization requires condition documents, aggregate or gene_expression terms,",
        "and topic_term_assignment_method = 'gammafit_maxprob'."
      )
    )
  }

  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("Matrix")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(title_prefix)) title_prefix <- basename(out_dir)
  topic_model_family <- .module3_topic_model_family(c(topic_model_family, title_prefix))
  step_timing <- list()
  .record_step <- function(step, code) {
    start_time <- Sys.time()
    start_elapsed <- proc.time()[["elapsed"]]
    status <- "ok"
    error_message <- NA_character_
    value <- tryCatch(
      force(code),
      error = function(e) {
        status <<- "error"
        error_message <<- conditionMessage(e)
        stop(e)
      },
      finally = {
        elapsed <- proc.time()[["elapsed"]] - start_elapsed
        step_timing[[length(step_timing) + 1L]] <<- data.table::data.table(
          step = step,
          start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
          elapsed_seconds = as.numeric(elapsed),
          status = status,
          error_message = error_message
        )
      }
    )
    value
  }
  .write_step_timing <- function() {
    if (!length(step_timing)) return(invisible(NULL))
    timing_dt <- data.table::rbindlist(step_timing, use.names = TRUE, fill = TRUE)
    data.table::fwrite(timing_dt, file.path(out_dir, "topic_extraction_step_timing.csv"))
    invisible(timing_dt)
  }
  on.exit(.write_step_timing(), add = TRUE)

  theta <- .validate_topic_probability_matrix(topic_base$theta, "theta")
  phi <- .validate_topic_probability_matrix(topic_base$phi, "phi")
  model_k <- nrow(phi)

  if (is.null(rownames(theta)) && !is.null(rownames(dtm))) {
    rownames(theta) <- rownames(dtm)
  }
  if (is.null(colnames(theta))) colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  if (is.null(rownames(phi))) rownames(phi) <- colnames(theta)
  if (is.null(colnames(phi)) && !is.null(colnames(dtm))) colnames(phi) <- colnames(dtm)

  topic_base$theta <- theta
  topic_base$phi <- phi

  score_mat <- .record_step(
    "score_terms",
    score_terms_normtop(phi, method = topic_score_method)
  )
  gamma_score_mat <- score_mat
  candidate_topic_terms <- .record_step(
    "binarize_topics",
    binarize_topics(
      score_mat,
      method = binarize_method,
      thrP = thrP,
      top_n_terms = top_n_terms,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope
    )
  )
  topic_terms <- .record_step(
    "assign_topic_terms",
    .assign_topic_terms(
      raw_phi = phi,
      score_mat = score_mat,
      candidate_terms = candidate_topic_terms,
      method = topic_term_assignment_method
    )
  )
  data.table::fwrite(topic_terms, file.path(out_dir, "topic_terms.csv"))
  .save_all(out_dir, "topic_terms", topic_terms)
  .save_all(out_dir, "topic_term_scores_normtop", score_mat)
  .save_all(out_dir, paste0("topic_term_scores_", topic_score_method), score_mat)
  writeLines(topic_score_method, file.path(out_dir, "topic_score_method.txt"))
  writeLines(topic_term_assignment_method, file.path(out_dir, "topic_term_assignment_method.txt"))
  if (identical(topic_term_assignment_method, "gammafit_maxprob")) {
    gene_only_assignment <- identical(fp_term_mode, "gene_expression")
    pair_assignment <- if (gene_only_assignment) {
      .topic_gene_assignment_table(topic_terms)
    } else {
      .topic_gene_peak_assignment_table(topic_terms)
    }
    if (gene_only_assignment) {
      if (!nrow(pair_assignment) || !any(pair_assignment$assigned)) {
        .log_abort("gene_expression GammaFit-MaxProb did not assign any target genes.")
      }
      data.table::fwrite(
        pair_assignment,
        file.path(out_dir, "topic_gene_assignment.csv")
      )
      data.table::fwrite(
        .topic_gene_assignment_summary(
          pair_assignment,
          thrP = thrP,
          assignment_method = topic_term_assignment_method,
          model_family = topic_model_family
        ),
        file.path(out_dir, "topic_term_assignment_summary.csv")
      )
    } else {
    pair_columns <- c("gene_term_id", "peak_term_id")
    matched_pair_count <- if (nrow(pair_assignment) && all(pair_columns %in% names(pair_assignment))) {
      pair_assignment[!is.na(gene_term_id) & !is.na(peak_term_id), .N]
    } else {
      0L
    }
    if (!is.finite(matched_pair_count) || matched_pair_count < 1L) {
      .log_abort(
        "gammafit_maxprob requires matched aggregate GENE:<gene> and PEAK:<gene> terms."
      )
    }
    data.table::fwrite(
      pair_assignment,
      file.path(out_dir, "topic_gene_peak_assignment.csv")
    )
    data.table::fwrite(
      .topic_gene_peak_assignment_summary(
        pair_assignment,
        thrP = thrP,
        assignment_method = topic_term_assignment_method,
        model_family = topic_model_family
      ),
      file.path(out_dir, "topic_term_assignment_summary.csv")
    )
    }
    if (isTRUE(optimize_topics)) {
      raw_k <- nrow(phi)
      data.table::fwrite(
        topic_terms,
        file.path(out_dir, "topic_terms_raw.csv")
      )
      raw_assignment_file <- if (gene_only_assignment) {
        "topic_gene_assignment_raw.csv"
      } else {
        "topic_gene_peak_assignment_raw.csv"
      }
      data.table::fwrite(pair_assignment, file.path(out_dir, raw_assignment_file))
      optimizer_assignment <- if (gene_only_assignment) {
        .topic_gene_assignment_for_optimization(pair_assignment, topic_terms)
      } else {
        pair_assignment
      }
      qc_comparison_audit <- NULL
      qc_condition_gene_expression <- if (!is.null(topic_qc_condition_expression_file)) {
        expression_file <- path.expand(as.character(topic_qc_condition_expression_file[[1L]]))
        if (!file.exists(expression_file)) {
          .log_abort("topic_qc_condition_expression_file does not exist: {expression_file}")
        }
        qc_comparison_audit <- data.table::fread(expression_file, showProgress = FALSE)
        .m3_qc_condition_expression_from_comparisons(
          qc_comparison_audit,
          genes = pair_assignment$target_gene
        )
      } else {
        .m3_qc_condition_gene_expression(
          edges_docs,
          genes = pair_assignment$target_gene
        )
      }
      optimization <- .record_step(
        "optimize_topics",
        .module3_optimize_condition_topics(
          theta = theta,
          phi = phi,
          dtm = dtm,
          topic_terms = topic_terms,
          pair_assignment = optimizer_assignment,
          assignment_mode = if (gene_only_assignment) "gene_only" else "gene_peak",
          condition_gene_expression = qc_condition_gene_expression,
          condition_upregulated_genes = if (!is.null(qc_comparison_audit)) {
            .m3_qc_upregulated_genes_from_comparisons(
              qc_comparison_audit,
              genes = pair_assignment$target_gene,
              reference_condition = topic_qc_reference_condition,
              log2fc_min = topic_qc_upregulated_log2fc_min
            )
          } else {
            .m3_qc_upregulated_genes_vs_reference(
              qc_condition_gene_expression,
              reference_condition = topic_qc_reference_condition,
              log2fc_min = topic_qc_upregulated_log2fc_min,
              pseudocount = topic_qc_upregulated_pseudocount
            )
          },
          upregulation_reference_condition = topic_qc_reference_condition,
          upregulated_log2fc_min = topic_qc_upregulated_log2fc_min,
          min_genes = topic_merge_min_genes,
          min_links = topic_merge_min_links,
          similarity_threshold = topic_merge_similarity_threshold,
          tf_topic_cutoff = topic_tf_membership_cutoff,
          umap_max_links_per_condition = topic_qc_umap_links_per_condition,
          seed = topic_qc_seed
        )
      )
      theta <- optimization$theta
      phi <- optimization$phi
      score_mat <- optimization$score
      topic_terms <- optimization$topic_terms
      pair_assignment <- optimization$pair_assignment
      topic_base$theta <- theta
      topic_base$phi <- phi
      data.table::fwrite(topic_terms, file.path(out_dir, "topic_terms.csv"))
      final_assignment_file <- if (gene_only_assignment) {
        "topic_gene_assignment.csv"
      } else {
        "topic_gene_peak_assignment.csv"
      }
      data.table::fwrite(pair_assignment, file.path(out_dir, final_assignment_file))
      .save_all(out_dir, "topic_terms", topic_terms)
      .save_all(out_dir, "topic_term_scores_normtop", score_mat)
      .save_all(out_dir, paste0("topic_term_scores_", topic_score_method), score_mat)
      assignment_summary <- if (gene_only_assignment) {
        .topic_gene_assignment_summary(
          pair_assignment,
          thrP = thrP,
          assignment_method = "gammafit_maxprob_optimized",
          model_family = topic_model_family
        )
      } else {
        .topic_gene_peak_assignment_summary(
          pair_assignment,
          thrP = thrP,
          assignment_method = "gammafit_maxprob_optimized",
          model_family = topic_model_family
        )
      }
      data.table::fwrite(
        assignment_summary,
        file.path(out_dir, "topic_term_assignment_summary.csv")
      )
      if (isTRUE(run_topic_assignment_qc)) {
        .record_step(
          "topic_assignment_qc",
          .write_module3_topic_optimization_outputs(
            optimization = optimization,
            out_dir = out_dir,
            raw_k = raw_k,
            title_prefix = title_prefix,
            condition_colors = theta_umap_condition_colors,
            top_n_tfs = topic_qc_top_tfs,
            seed = topic_qc_seed
          )
        )
      }
    }
  }
  if (identical(binarize_method, "gammafit")) {
    gamma_diagnostics_tbl <- .gammafit_diagnostics_by_termclass(
      gamma_score_mat,
      topic_terms = candidate_topic_terms,
      topic_score_method = topic_score_method,
      thrP = thrP,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope
    )
    data.table::fwrite(
      gamma_diagnostics_tbl,
      file.path(out_dir, "topic_gammafit_diagnostics.csv")
    )
    gamma_cutoffs_tbl <- .gammafit_cutoffs_by_termclass(
      gamma_score_mat,
      thrP = thrP,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope
    )
    data.table::fwrite(
      gamma_cutoffs_tbl,
      file.path(out_dir, "topic_gamma_cutoffs.csv")
    )
  }

  topic_links_tbl <- NULL
  if (isTRUE(run_link_topic_scores)) {
    topic_links_tbl <- .record_step(
      "topic_links",
      {
        link_tbl <- compute_topic_links(
          edges_docs = edges_docs,
          score_mat = score_mat,
          raw_score_mat = topic_base$phi,
          topic_terms = topic_terms,
          use_final_term_assignment = identical(
            topic_term_assignment_method,
            "gammafit_maxprob"
          ),
          theta = topic_base$theta,
          topic_tf_membership_cutoff = topic_tf_membership_cutoff,
          fp_term_mode = fp_term_mode,
          binarize_method = binarize_method,
          gammafit_scope = gammafit_scope,
          link_method = link_topic_method,
          link_prob_cutoff = link_topic_prob_cutoff,
          thrP = thrP,
          min_terms = in_topic_min_terms,
          fdr_q = link_topic_fdr_q,
          fdr_p = link_topic_fdr_p,
          efdr_scope = link_topic_efdr_scope,
          efdr_B = link_topic_efdr_B,
          efdr_seed = link_topic_efdr_seed,
          out_file = file.path(out_dir, "topic_links.csv"),
          output_mode = link_topic_output,
          pass_file = file.path(out_dir, "topic_links_pass.csv"),
          chunk_size = link_topic_chunk_size,
          n_cores = link_topic_n_cores,
          overwrite = link_topic_overwrite
        )
        if (!is.data.frame(link_tbl)) {
          topic_links_path <- .topic_links_path(out_dir, prefer = "pass")
          if (file.exists(topic_links_path)) {
            link_tbl <- data.table::fread(topic_links_path)
          }
        }
        link_tbl
      }
    )
  }

  if (binarize_method == "gammafit" && isTRUE(run_gammafit_summary)) {
    .record_step(
      "gammafit_summary",
      {
        tf_terms <- if (is.data.frame(edges_docs) && "tf" %in% names(edges_docs)) unique(edges_docs$tf) else NULL
        topic_terms_tbl <- data.table::fread(file.path(out_dir, "topic_terms.csv"))
        item_coverage_tbl <- if (file.exists(file.path(out_dir, "topic_item_coverage_counts.csv"))) {
          data.table::fread(file.path(out_dir, "topic_item_coverage_counts.csv"))
        } else {
          NULL
        }
        plot_gammafit_binarize(
          score_mat,
          out_file = file.path(out_dir, "topic_terms_and_cutoffs_summary.pdf"),
          thrP = thrP,
          min_terms = in_topic_min_terms,
          gammafit_scope = gammafit_scope,
          title_prefix = title_prefix,
          tf_list = tf_terms,
          edges_docs = edges_docs,
          topic_terms = topic_terms_tbl,
          topic_links = topic_links_tbl,
          item_coverage = item_coverage_tbl,
          show_peak_expanded_link_coverage = identical(fp_term_mode, "unique")
        )
      }
    )
  }
  if (isTRUE(run_link_efdr_summary) && identical(link_topic_method, "link_score_efdr") && is.data.frame(topic_links_tbl) && nrow(topic_links_tbl)) {
    plot_link_efdr_summary(
      topic_links = topic_links_tbl,
      out_file = file.path(out_dir, "topic_link_efdr_summary.pdf"),
      title_prefix = title_prefix,
      fdr_q = link_topic_fdr_q
    )
  }

  if (isTRUE(run_tf_topic_assignment)) {
    tf_topic_assign <- .record_step(
      "tf_topic_assignment",
      .write_tf_topic_assignment_outputs(
        theta = topic_base$theta,
        out_dir = out_dir,
        doc_design = doc_design,
        membership_cutoff = topic_tf_membership_cutoff,
        primary_margin_cutoff = topic_tf_primary_margin_cutoff
      )
    )
    .replace_tf_item_coverage_with_assignment(out_dir, tf_topic_assign)
    if (binarize_method == "gammafit" && isTRUE(run_gammafit_summary)) {
      .record_step(
        "gammafit_summary_after_tf_assignment",
        {
          tf_terms <- if (is.data.frame(edges_docs) && "tf" %in% names(edges_docs)) unique(edges_docs$tf) else NULL
          topic_terms_tbl <- data.table::fread(file.path(out_dir, "topic_terms.csv"))
          item_coverage_tbl <- if (file.exists(file.path(out_dir, "topic_item_coverage_counts.csv"))) {
            data.table::fread(file.path(out_dir, "topic_item_coverage_counts.csv"))
          } else {
            NULL
          }
          plot_gammafit_binarize(
            score_mat,
            out_file = file.path(out_dir, "topic_terms_and_cutoffs_summary.pdf"),
            thrP = thrP,
            min_terms = in_topic_min_terms,
            gammafit_scope = gammafit_scope,
            title_prefix = title_prefix,
            tf_list = tf_terms,
            edges_docs = edges_docs,
            topic_terms = topic_terms_tbl,
            topic_links = topic_links_tbl,
            item_coverage = item_coverage_tbl,
            show_gammafit_pages = identical(topic_term_assignment_method, "gammafit"),
            show_peak_expanded_link_coverage = identical(fp_term_mode, "unique")
          )
        }
      )
    }
  }
  if (isTRUE(run_raw_theta_document_heatmap)) {
    .record_step(
      "raw_theta_document_heatmap",
      {
        raw_theta_path <- file.path(out_dir, sprintf("raw_theta_documents_K%d.pdf", model_k))
        .plot_raw_theta_document_heatmap(
          theta = topic_base$theta,
          out_file = raw_theta_path,
          doc_design = doc_design,
          title_prefix = title_prefix
        )
      }
    )
  }
  if (isTRUE(run_document_theta_umap)) {
    .record_step(
      "document_theta_umap",
      {
        theta_umap_path <- file.path(out_dir, sprintf("document_theta_umap_K%d.pdf", model_k))
        .plot_document_theta_umap(
          theta = topic_base$theta,
          out_file = theta_umap_path,
          doc_design = doc_design,
          selected_tfs = theta_umap_selected_tfs,
          top_n_tfs = theta_umap_top_n_tfs,
          seed = theta_umap_seed,
          n_neighbors = theta_umap_n_neighbors,
          condition_colors = theta_umap_condition_colors,
          title_prefix = title_prefix
        )
      }
    )
  }
  if (isTRUE(run_topic_term_heatmap)) {
    .record_step(
      "topic_term_heatmaps",
      {
        term_compare_path <- file.path(out_dir, sprintf("topic_term_phi_score_heatmap_K%d.pdf", model_k))
        .write_topic_term_primary_assignment(
          score_mat = score_mat,
          topic_terms = topic_terms,
          assignment_file = file.path(out_dir, "topic_term_primary_assignment.csv")
        )
        .plot_topic_term_phi_score_comparison_heatmap(
          phi = phi,
          score_mat = score_mat,
          topic_terms = topic_terms,
          out_file = term_compare_path,
          topic_score_method = topic_score_method,
          title_prefix = title_prefix
        )
      }
    )
  }
  direction_mode <- if (option_label == "opt3_gene_fc_expr") {
    "gene"
  } else if (option_label == "joint") {
    if (direction_by == "fp") "fp" else "gene"
  } else {
    "fp"
  }

  if (isTRUE(run_topic_by_comparison_heatmaps)) {
    .record_step(
      "topic_by_comparison_heatmaps",
      plot_topic_by_comparison_heatmaps(
        theta = topic_base$theta,
        out_dir = out_dir,
        edges_docs = edges_docs,
        direction_mode = direction_mode,
        title_prefix = title_prefix,
        label_cleaner = topic_by_comparison_label_cleaner,
        topic_links = topic_links_tbl,
        annotate_unique_genes = TRUE,
        doc_design = doc_design
      )
    )
  }

  if (isTRUE(run_pathway_enrichment) && isTRUE(pathway_overwrite)) {
    .record_step(
      "pathway_cleanup",
      {
        unlink(list.files(out_dir, pattern = "^topic_pathway_enrichment", full.names = TRUE), recursive = TRUE, force = TRUE)
        if (!is.null(pathway_per_comparison_dir) &&
            nzchar(as.character(pathway_per_comparison_dir)[[1L]]) &&
            !identical(as.character(pathway_per_comparison_dir)[[1L]], ".")) {
          unlink(file.path(out_dir, pathway_per_comparison_dir), recursive = TRUE, force = TRUE)
        } else {
          unlink(
            file.path(
              out_dir,
              c(
                "per_comparison_topic_pathway_enrichment.csv",
                "per_comparison_topic_pathway_debug.txt",
                "per_condition_topic_pathway_enrichment.csv",
                "per_condition_topic_pathway_debug.txt",
                "topic_term_pathway_enrichment.csv",
                "topic_term_pathway_debug.txt"
              )
            ),
            force = TRUE
          )
        }
        unlink(list.files(out_dir, pattern = paste0("^", .safe_filename(pathway_per_comparison_dir), "_"), full.names = TRUE), recursive = TRUE, force = TRUE)
      }
    )
  }

  if (isTRUE(run_pathway_enrichment)) {
    .record_step(
      "pathway_enrichment",
      {
    if (pathway_source == "link_scores") {
      method_secondary <- if (identical(link_topic_method, "gene_prob")) "peak_and_gene_prob" else link_topic_method
      report_methods <- unique(c("peak_and_gene", method_secondary))
      if (!is.null(link_topic_report_methods)) {
        report_methods <- unique(as.character(link_topic_report_methods))
      }
      report_methods <- intersect(
        report_methods,
        c("peak_and_gene", "peak_and_gene_prob", "gene_only", "link_score_prob", "link_score_efdr")
      )
      for (method in report_methods) {
        rerun_pathway_from_topic_links(
          out_dir = out_dir,
          topic_links_file = pathway_link_scores_file,
          method = method,
          allow_missing = TRUE,
          title_prefix = title_prefix,
          include_tf = pathway_link_include_tf,
          include_gene = pathway_link_include_gene,
          min_prob = pathway_link_min_prob,
          gene_min_prob = pathway_link_gene_min_prob,
          tf_min_prob = pathway_link_tf_min_prob,
          tf_max_topics = pathway_link_tf_max_topics,
          tf_top_n_per_topic = pathway_link_tf_top_n_per_topic,
          top_n_per_topic = top_n_per_topic,
          dot_top_n_per_topic = dot_top_n_per_topic,
          max_pathways = max_pathways,
          per_comparison = pathway_per_comparison,
          per_comparison_dir = paste0(pathway_per_comparison_dir, "_", .pathway_method_suffix(method)),
          per_comparison_flat = pathway_per_comparison_flat,
          split_direction = pathway_split_direction,
          make_heatmap = pathway_make_heatmap,
          make_dotplot = pathway_make_dotplot,
          enrichr_sleep_time = pathway_enrichr_sleep_time,
          enrichr_cache_dir = pathway_enrichr_cache_dir,
          enrichr_n_cores = pathway_enrichr_n_cores,
          dbs = pathway_dbs,
          pathway_species = pathway_species,
          pathway_backend = pathway_backend,
          overwrite = pathway_overwrite,
          doc_design = doc_design
        )
      }
    } else {
      plot_topic_pathway_enrichment_heatmap(
        topic_terms = topic_terms,
        edges_docs = edges_docs,
        option_label = option_label,
        out_file = file.path(out_dir, "topic_pathway_enrichment_heatmap.pdf"),
        title_prefix = title_prefix,
        use_all_terms = pathway_use_all_terms,
        make_heatmap = pathway_make_heatmap,
        make_dotplot = pathway_make_dotplot,
        top_n_per_topic = top_n_per_topic,
        dot_top_n_per_topic = dot_top_n_per_topic,
        max_pathways = max_pathways,
        theta = topic_base$theta,
        tf_link_mode = pathway_tf_link_mode,
        tf_theta_top_n = pathway_tf_top_n_docs,
        tf_theta_min = pathway_tf_min_theta,
        enrichr_sleep_time = pathway_enrichr_sleep_time,
        enrichr_cache_dir = pathway_enrichr_cache_dir,
        enrichr_n_cores = pathway_enrichr_n_cores,
        dbs = pathway_dbs,
        pathway_species = pathway_species,
        pathway_backend = pathway_backend
      )

      if (isTRUE(pathway_per_comparison)) {
        retest_doc_design <- if (identical(doc_design, "condition")) "condition" else "comparison"
        plot_topic_pathway_enrichment_by_comparison_terms(
          topic_terms = topic_terms,
          edges_docs = edges_docs,
          theta = topic_base$theta,
          out_dir = out_dir,
          title_prefix = title_prefix,
          theta_min = pathway_topic_term_theta_min,
          include_peak_terms = !identical(option_label, "opt3_gene_fc_expr"),
          use_all_terms = pathway_use_all_terms,
          per_comparison_dir = pathway_per_comparison_dir,
          per_comparison_flat = pathway_per_comparison_flat,
          split_direction = pathway_split_direction,
          min_genes = 5L,
          top_n_per_topic = top_n_per_topic,
          dot_top_n_per_topic = dot_top_n_per_topic,
          max_pathways = max_pathways,
          make_dotplot = pathway_make_dotplot,
          enrichr_sleep_time = pathway_enrichr_sleep_time,
          enrichr_cache_dir = pathway_enrichr_cache_dir,
          enrichr_n_cores = pathway_enrichr_n_cores,
          dbs = pathway_dbs,
          pathway_species = pathway_species,
          pathway_backend = pathway_backend,
          overwrite = pathway_overwrite,
          doc_design = retest_doc_design
        )
      }

      if (isTRUE(run_pathway_gsea)) {
        plot_topic_pathway_enrichment_gsea(
          topic_terms = topic_terms,
          edges_docs = edges_docs,
          option_label = option_label,
          out_dir = out_dir,
          theta = topic_base$theta,
          species = gsea_species,
          padj_cut = 0.05,
          min_size = 10L,
          max_size = 500L,
          nperm = gsea_nperm,
          top_n_per_topic = 20L,
          max_pathways = 200L,
          peak_gene_agg = gsea_peak_gene_agg,
          tf_link_mode = pathway_tf_link_mode,
          tf_theta_top_n = pathway_tf_top_n_docs,
          tf_theta_min = pathway_tf_min_theta,
          title_prefix = title_prefix
        )
      }
    }
      }
    )
  }

  if (isTRUE(run_intertopic_distance_map)) {
    .record_step(
      "intertopic_distance_map",
      plot_intertopic_distance_map(
        phi = topic_base$phi,
        topic_terms = topic_terms,
        out_file = file.path(out_dir, "intertopic_distance_map.pdf"),
        option_label = option_label,
        title_prefix = title_prefix
      )
    )
  }
  invisible(list(
    out_dir = out_dir,
    edges_docs = edges_docs,
    dtm = dtm,
    topic_terms = topic_terms
  ))
}

# =============================================================================
# 9) Main runner for ONE term option (WarpLDA-only; saves everything)
# =============================================================================

run_tfdocs_warplda_one_option <- function(edges_all,
                                          out_dir,
                                          option_label = c("opt1_peak_delta_fp", "opt2_peak_fc_fp", "opt3_gene_fc_expr", "joint"),
                                          # QC + docs
                                          doc_mode = c("tf", "tf_cluster", "comparison"),
                                          tf_cluster_map = NULL,
                                          direction_by = c("gene", "fp", "none"),
                                          abs_log2fc_fp_min = 1,
                                          abs_delta_fp_min = NA_real_,
                                          abs_log2fc_gene_min = 1,
                                          require_fp_bound_either = TRUE,
                                          require_tf_expr_either = TRUE,
                                          require_gene_expr_either = TRUE,
                                          direction_consistency = c("aligned", "none"),
                                          # doc-term
                                          top_terms_per_doc = 500L,
                                          min_df = 2L,
                                          count_method = c("bin", "log"),
                                          count_scale = 50,
                                          distinct_terms = FALSE,
                                          joint_weight_fp = c("delta_fp", "fc_mag_fp"),
                                          joint_weight_gene = c("fc_mag_gene"),
                                          joint_balance = c("min", "none"),
                                          # WarpLDA
                                          K_grid = c(2:15, 20, 25, 35, 40, 45, 50),
                                          iterations = 2000L,
                                          alpha_by_topic = TRUE,
                                          beta = NULL,
                                          sampler = c("warp_omp", "warp_ref", "warp_mh", "gibbs_sync"),
                                          seed = 123,
                                          # topic definition
                                          topic_score_method = c("normtop_specificity", "rowmax_phi"),
                                          topic_term_assignment_method = c("max_phi", "gammafit"),
                                          binarize_method = c("gammafit", "topn"),
                                          gammafit_scope = c("topic_term_group", "global_term_group"),
                                          thrP = 0.975,
                                          top_n_terms = 500L,
                                          in_topic_min_terms = 50L,
                                          pathway_use_all_terms = FALSE,
                                          pathway_make_heatmap = FALSE,
                                          pathway_make_dotplot = TRUE,
                                          top_n_per_topic = 20L,
                                          max_pathways = 200L,
                                          pathway_tf_link_mode = c("theta", "none"),
                                          pathway_tf_top_n_docs = 50L,
                                          pathway_tf_min_theta = NA_real_,
                                          topic_tf_membership_cutoff = 0.3,
                                          run_pathway_gsea = FALSE,
                                          gsea_species = "Homo sapiens",
                                          gsea_nperm = 1000L,
                                          gsea_peak_gene_agg = c("max", "sum"),
                                          pathway_source = c("topic_terms", "link_scores"),
                                          pathway_link_scores_file = NULL,
                                          pathway_link_scores_file_tf = NULL,
                                          pathway_link_gene_terms_file = NULL,
                                          pathway_link_min_prob = 0,
                                          pathway_link_include_tf = TRUE,
                                          pathway_link_include_gene = TRUE,
                                          pathway_link_gene_min_prob = NULL,
                                          pathway_link_tf_min_prob = NULL,
                                          pathway_link_tf_max_topics = Inf,
                                          pathway_link_tf_top_n_per_topic = NA_integer_,
                                          pathway_per_comparison = FALSE,
                                          pathway_per_comparison_dir = ".",
                                          pathway_per_comparison_flat = TRUE,
                                          pathway_split_direction = TRUE,
                                          pathway_overwrite = FALSE,
                                          pathway_enrichr_sleep_time = 0,
                                          pathway_enrichr_cache_dir = NULL,
                                          pathway_enrichr_n_cores = NULL,
                                          pathway_backend = NULL,
                                          pathway_species = NULL,
                                          pathway_databases = NULL,
                                          run_link_topic_scores = FALSE,
                                          fp_term_mode = c("unique", "aggregate", "aggregate_weight", "gene_expression"),
                                          link_topic_gate_mode = "none",
                                          link_topic_top_k = 3L,
                                          link_topic_min_prob = 0,
                                          link_topic_include_tf = FALSE,
                                          link_topic_chunk_size = 5000L,
                                          link_topic_n_cores = 1L,
                                          link_topic_overwrite = FALSE,
                                          link_topic_method = c("gammafit", "theta_and_terms", "gene_prob", "link_score_prob", "link_score_efdr"),
                                          link_topic_prob_cutoff = 0.3,
                                          link_topic_fdr_q = 0.2,
                                          link_topic_fdr_p = NA_real_,
                                          link_topic_efdr_scope = c("per_topic", "global"),
                                          link_topic_efdr_B = 100L,
                                          link_topic_efdr_seed = 1L,
                                          link_topic_report_methods = NULL,
                                          topic_by_comparison_label_cleaner = NULL) {
  option_label <- match.arg(option_label)
  doc_mode <- match.arg(doc_mode)
  direction_by <- match.arg(direction_by)
  direction_consistency <- match.arg(direction_consistency)
  count_method <- match.arg(count_method)
  joint_weight_fp <- match.arg(joint_weight_fp)
  joint_weight_gene <- match.arg(joint_weight_gene)
  joint_balance <- match.arg(joint_balance)
  sampler <- match.arg(sampler)
  topic_score_method <- match.arg(topic_score_method)
  topic_term_assignment_method <- match.arg(topic_term_assignment_method)
  binarize_method <- match.arg(binarize_method)
  gammafit_scope <- match.arg(gammafit_scope)
  gsea_peak_gene_agg <- match.arg(gsea_peak_gene_agg)
  pathway_tf_link_mode <- match.arg(pathway_tf_link_mode)
  pathway_source <- match.arg(pathway_source)
  pathway_enrichr_sleep_time <- .normalize_enrichr_sleep_time(pathway_enrichr_sleep_time)
  pathway_enrichr_n_cores <- .normalize_enrichr_n_cores(pathway_enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  if (is.null(pathway_enrichr_cache_dir)) {
    pathway_enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir, backend = pathway_backend)
  }
  pathway_species <- .normalize_pathway_species_mode(pathway_species)
  pathway_dbs <- if (is.null(pathway_databases)) {
    .default_pathway_databases(pathway_species)
  } else {
    unique(as.character(pathway_databases))
  }
  pathway_dbs <- pathway_dbs[!is.na(pathway_dbs) & nzchar(pathway_dbs)]
  if (!length(pathway_dbs)) .log_abort("pathway_databases must contain at least one database.")
  pathway_topic_term_theta_min <- .safe_num(pathway_tf_min_theta)
  if (!is.finite(pathway_topic_term_theta_min)) {
    pathway_topic_term_theta_min <- .safe_num(topic_tf_membership_cutoff)
  }
  if (!is.finite(pathway_topic_term_theta_min)) pathway_topic_term_theta_min <- 0.3
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  link_topic_efdr_scope <- match.arg(link_topic_efdr_scope)
  link_topic_method <- match.arg(link_topic_method)
  allowed_gate_modes <- c("none", "peak_in_set", "gene_in_set", "peak_and_gene_in_set")
  link_topic_gate_mode <- unique(as.character(link_topic_gate_mode))
  if (!all(link_topic_gate_mode %in% allowed_gate_modes)) {
    .log_abort("link_topic_gate_mode must be one of: {paste(allowed_gate_modes, collapse = ', ')}.")
  }

  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("Matrix")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "tmp_models"), recursive = TRUE, showWarnings = FALSE)

  # record calc params (cisTopic-like)
  calc_params <- list(
    option_label = option_label,
    doc_mode = doc_mode,
    direction_by = direction_by,
    abs_log2fc_fp_min = abs_log2fc_fp_min,
    abs_delta_fp_min = abs_delta_fp_min,
    abs_log2fc_gene_min = abs_log2fc_gene_min,
    require_fp_bound_either = require_fp_bound_either,
    require_tf_expr_either = require_tf_expr_either,
    require_gene_expr_either = require_gene_expr_either,
    direction_consistency = direction_consistency,
    doc_term = list(
      top_terms_per_doc = as.integer(top_terms_per_doc),
      min_df = as.integer(min_df),
      count_method = count_method,
      count_scale = as.numeric(count_scale),
      distinct_terms = isTRUE(distinct_terms)
    ),
    runWarpLDAModels = list(
      seed = as.integer(seed),
      iterations = as.integer(iterations),
      alpha = if (isTRUE(alpha_by_topic)) "50/K" else NA_character_,
      alphaByTopic = isTRUE(alpha_by_topic),
      beta = if (is.null(beta)) "1/K" else as.numeric(beta),
      sampler = sampler,
      K_grid = as.integer(K_grid)
    ),
    topic_definition = list(
      NormTop = TRUE,
      topic_score_method = topic_score_method,
      topic_term_assignment_method = topic_term_assignment_method,
      binarize_method = binarize_method,
      gammafit_scope = gammafit_scope,
      thrP = as.numeric(thrP),
      top_n_terms = as.integer(top_n_terms)
    ),
    pathway = list(
      top_n_per_topic = if (is.finite(top_n_per_topic) && as.numeric(top_n_per_topic) > 0) as.integer(top_n_per_topic) else Inf,
      max_pathways = if (is.finite(max_pathways) && as.numeric(max_pathways) > 0) as.integer(max_pathways) else Inf,
      use_all_terms = isTRUE(pathway_use_all_terms),
      make_heatmap = isTRUE(pathway_make_heatmap),
      make_dotplot = isTRUE(pathway_make_dotplot),
      source = pathway_source,
      link_min_prob = as.numeric(pathway_link_min_prob),
      link_include_tf = isTRUE(pathway_link_include_tf),
      link_include_gene = isTRUE(pathway_link_include_gene),
      link_gene_min_prob = as.numeric(pathway_link_gene_min_prob),
      link_tf_min_prob = as.numeric(pathway_link_tf_min_prob),
      link_tf_max_topics = as.numeric(pathway_link_tf_max_topics),
      link_tf_top_n_per_topic = as.integer(pathway_link_tf_top_n_per_topic),
      per_comparison = isTRUE(pathway_per_comparison),
      per_comparison_dir = pathway_per_comparison_dir,
      split_direction = isTRUE(pathway_split_direction),
      tf_link_mode = pathway_tf_link_mode,
      tf_theta_top_n = as.integer(pathway_tf_top_n_docs),
      tf_theta_min = as.numeric(pathway_tf_min_theta)
    )
  )
  if (option_label == "joint") {
    calc_params$doc_term$joint <- list(
      weight_fp = joint_weight_fp,
      weight_gene = joint_weight_gene,
      balance = joint_balance
    )
  }
  .save_all(out_dir, "calc_params", calc_params)

  edges_in <- data.table::as.data.table(edges_all)

  # 1) QC filter
  edges_filt <- filter_edges_for_tf_topics(
    edges_in,
    abs_log2fc_fp_min = abs_log2fc_fp_min,
    abs_delta_fp_min = abs_delta_fp_min,
    abs_log2fc_gene_min = abs_log2fc_gene_min,
    require_fp_bound_either = require_fp_bound_either,
    require_tf_expr_either = require_tf_expr_either,
    require_gene_expr_either = require_gene_expr_either,
    direction_consistency = direction_consistency
  )
  if (!nrow(edges_filt)) .log_abort("No edges passed filtering; relax thresholds or inspect standardized columns.")
  .save_all(out_dir, "edges_filtered", edges_filt)

  # 2) docs
  edges_docs <- add_tf_docs(
    edges_filt,
    doc_mode = doc_mode,
    direction_by = direction_by,
    tf_cluster_map = tf_cluster_map
  )
  if (!nrow(edges_docs)) .log_abort("No edges remained after doc assignment; check direction_by/log2fc_gene.")
  .save_all(out_dir, "edges_docs", edges_docs)

  # 3) doc-term (requested options + joint)
  doc_term <- switch(
    option_label,
    opt1_peak_delta_fp = build_doc_term_opt1_peak_delta_fp(
      edges_docs,
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      prefix_terms = TRUE,
      distinct_terms = distinct_terms
    ),
    opt2_peak_fc_fp = build_doc_term_opt2_peak_fc_fp(
      edges_docs,
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      prefix_terms = TRUE,
      distinct_terms = distinct_terms
    ),
    opt3_gene_fc_expr = build_doc_term_opt3_gene_fc_expr(
      edges_docs,
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      prefix_terms = TRUE,
      distinct_terms = distinct_terms
    ),
    joint = build_doc_term_joint(
      edges_docs,
      weight_type_peak = joint_weight_fp,
      weight_type_gene = joint_weight_gene,
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      distinct_terms = distinct_terms,
      balance_mode = joint_balance,
      prefix_terms = TRUE
    )
  )

  if (!nrow(doc_term)) .log_abort("doc_term is empty; adjust thresholds/top_terms_per_doc/min_df.")

  data.table::fwrite(doc_term, file.path(out_dir, "doc_term.csv"))
  .save_all(out_dir, "doc_term", doc_term)

  # 4) DTM
  dtm_obj <- build_sparse_dtm(doc_term)
  dtm <- dtm_obj$dtm
  .save_all(out_dir, "dtm", dtm)
  .save_all(out_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))

  # 4b) runWarpLDAModels (multi-K) with tmp saves
  fits_out <- run_warplda_models(
    dtm,
    K_grid = K_grid,
    iterations = iterations,
    alpha_by_topic = alpha_by_topic,
    alpha = NULL,
    beta = beta,
    seed = seed,
    save_tmp_dir = file.path(out_dir, "tmp_models"),
    workers = 1L,
    threads_per_model = NULL,
    sampler = sampler,
    metrics_file = file.path(out_dir, "model_metrics.csv")
  )
  fit_files <- fits_out$fit_files
  metrics_tbl <- fits_out$metrics
  data.table::fwrite(metrics_tbl, file.path(out_dir, "model_metrics.csv"))
  .save_all(out_dir, "model_metrics", metrics_tbl)

  # 6) selection plots (cisTopic-like 3 panels)
  title_prefix <- .topic_model_selection_title(out_dir, backend_label = "WarpLDA")
  sel <- plot_model_selection_cistopic(metrics_tbl, file.path(out_dir, "model_selection.pdf"), title_prefix = title_prefix)
  .save_all(out_dir, "model_selection", sel)

  # Save "selected models" for each criterion (maximum, derivative knee, perplexity)
  m <- sel$table
  idx <- sel$idx

  pick_fit_by_K <- function(K) {
    key <- as.character(as.integer(K))
    fit_file <- fit_files[[key]]
    if (is.null(fit_file) || is.na(fit_file) || !nzchar(fit_file) || !file.exists(fit_file)) {
      .log_abort("Missing WarpLDA fit artifact for K={K}; rerun run_warplda_models to recreate tmp_models/fit_K*.rds.")
    }
    readRDS(fit_file)
  }

  selected <- list()

  if (is.finite(idx$maximum)) {
    bestK <- m$K[idx$maximum]
    fit <- pick_fit_by_K(bestK)
    selected$maximum <- list(K = bestK, theta = fit$theta, phi = fit$phi, metrics = fit$metrics)
    .save_all(out_dir, sprintf("theta_maximum_K%d", bestK), fit$theta)
    .save_all(out_dir, sprintf("phi_maximum_K%d", bestK), fit$phi)
  }

  if (is.finite(idx$derivative)) {
    bestK <- m$K[idx$derivative]
    fit <- pick_fit_by_K(bestK)
    selected$derivative <- list(K = bestK, theta = fit$theta, phi = fit$phi, metrics = fit$metrics)
    .save_all(out_dir, sprintf("theta_derivative_K%d", bestK), fit$theta)
    .save_all(out_dir, sprintf("phi_derivative_K%d", bestK), fit$phi)
  }

  if (is.finite(idx$perplexity)) {
    bestK <- m$K[idx$perplexity]
    fit <- pick_fit_by_K(bestK)
    selected$perplexity <- list(K = bestK, theta = fit$theta, phi = fit$phi, metrics = fit$metrics)
    .save_all(out_dir, sprintf("theta_perplexity_K%d", bestK), fit$theta)
    .save_all(out_dir, sprintf("phi_perplexity_K%d", bestK), fit$phi)
  }

  if (!length(selected)) {
    .log_inform("No valid model-selection indices; using first fitted model as fallback.")
    fit <- pick_fit_by_K(metrics_tbl$K[[1L]])
    selected$fallback <- list(K = fit$K, theta = fit$theta, phi = fit$phi, metrics = fit$metrics)
    .save_all(out_dir, sprintf("theta_fallback_K%d", fit$K), fit$theta)
    .save_all(out_dir, sprintf("phi_fallback_K%d", fit$K), fit$phi)
  }

  .save_all(out_dir, "selected_models", selected)

  # 7) topic definition (NormTop + binarize)
  # Do it for the perplexity-selected model by default (common WarpLDA practice),
  # but keep all selected outputs saved above.
  topic_base <- if (!is.null(selected$perplexity)) selected$perplexity else selected[[1]]
  if (is.null(topic_base)) .log_abort("No selected model available to define topics; check metrics/selection.")
  phi <- topic_base$phi

  score_mat <- score_terms_normtop(phi, method = topic_score_method)
  candidate_topic_terms <- binarize_topics(
    score_mat,
    method = binarize_method,
    thrP = thrP,
    top_n_terms = top_n_terms,
    min_terms = in_topic_min_terms,
    gammafit_scope = gammafit_scope
  )
  topic_terms <- .assign_topic_terms(
    raw_phi = phi,
    score_mat = score_mat,
    candidate_terms = candidate_topic_terms,
    method = topic_term_assignment_method
  )
  data.table::fwrite(topic_terms, file.path(out_dir, "topic_terms.csv"))
  .save_all(out_dir, "topic_terms", topic_terms)
  .save_all(out_dir, paste0("topic_term_scores_", topic_score_method), score_mat)
  writeLines(topic_score_method, file.path(out_dir, "topic_score_method.txt"))
  writeLines(topic_term_assignment_method, file.path(out_dir, "topic_term_assignment_method.txt"))
  if (identical(binarize_method, "gammafit")) {
    gamma_diagnostics_tbl <- .gammafit_diagnostics_by_termclass(
      score_mat,
      topic_terms = candidate_topic_terms,
      topic_score_method = topic_score_method,
      thrP = thrP,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope
    )
    data.table::fwrite(
      gamma_diagnostics_tbl,
      file.path(out_dir, "topic_gammafit_diagnostics.csv")
    )
    gamma_cutoffs_tbl <- .gammafit_cutoffs_by_termclass(
      score_mat,
      thrP = thrP,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope
    )
    data.table::fwrite(
      gamma_cutoffs_tbl,
      file.path(out_dir, "topic_gamma_cutoffs.csv")
    )
  }

  topic_links_tbl <- NULL
  if (isTRUE(run_link_topic_scores)) {
    topic_links_tbl <- compute_topic_links(
      edges_docs = edges_docs,
      score_mat = score_mat,
      raw_score_mat = phi,
      topic_terms = topic_terms,
      theta = topic_base$theta,
      topic_tf_membership_cutoff = topic_tf_membership_cutoff,
      fp_term_mode = fp_term_mode,
      binarize_method = binarize_method,
      gammafit_scope = gammafit_scope,
      link_method = link_topic_method,
      link_prob_cutoff = link_topic_prob_cutoff,
      thrP = thrP,
      min_terms = in_topic_min_terms,
      fdr_q = link_topic_fdr_q,
      fdr_p = link_topic_fdr_p,
      efdr_scope = link_topic_efdr_scope,
      efdr_B = link_topic_efdr_B,
      efdr_seed = link_topic_efdr_seed,
      out_file = file.path(out_dir, "topic_links.csv"),
      output_mode = "full",
      chunk_size = link_topic_chunk_size,
      n_cores = link_topic_n_cores,
      overwrite = link_topic_overwrite
    )
    if (!is.data.frame(topic_links_tbl)) {
      topic_links_path <- file.path(out_dir, "topic_links.csv")
      if (file.exists(topic_links_path)) {
        topic_links_tbl <- data.table::fread(topic_links_path)
      }
    }
  }

  if (binarize_method == "gammafit") {
    tf_terms <- if (is.data.frame(edges_docs) && "tf" %in% names(edges_docs)) unique(edges_docs$tf) else NULL
    topic_terms_tbl <- data.table::fread(file.path(out_dir, "topic_terms.csv"))
    item_coverage_tbl <- if (file.exists(file.path(out_dir, "topic_item_coverage_counts.csv"))) {
      data.table::fread(file.path(out_dir, "topic_item_coverage_counts.csv"))
    } else {
      NULL
    }
    plot_gammafit_binarize(
      score_mat,
      out_file = file.path(out_dir, "topic_terms_and_cutoffs_summary.pdf"),
      thrP = thrP,
      min_terms = in_topic_min_terms,
      gammafit_scope = gammafit_scope,
      title_prefix = title_prefix,
      tf_list = tf_terms,
      edges_docs = edges_docs,
      topic_terms = topic_terms_tbl,
      topic_links = topic_links_tbl,
      item_coverage = item_coverage_tbl,
      show_peak_expanded_link_coverage = identical(fp_term_mode, "unique")
    )
  }
  if (identical(link_topic_method, "link_score_efdr") && is.data.frame(topic_links_tbl) && nrow(topic_links_tbl)) {
    plot_link_efdr_summary(
      topic_links = topic_links_tbl,
      out_file = file.path(out_dir, "topic_link_efdr_summary.pdf"),
      title_prefix = title_prefix,
      fdr_q = link_topic_fdr_q
    )
  }

  # Also save full score matrix (can be large)
  .save_all(out_dir, "topic_term_scores_normtop", score_mat)

  direction_mode <- if (option_label == "opt3_gene_fc_expr") {
    "gene"
  } else if (option_label == "joint") {
    if (direction_by == "fp") "fp" else "gene"
  } else {
    "fp"
  }
  plot_topic_by_comparison_heatmaps(
    theta = topic_base$theta,
    out_dir = out_dir,
    edges_docs = edges_docs,
    direction_mode = direction_mode,
    title_prefix = title_prefix,
    label_cleaner = topic_by_comparison_label_cleaner,
    topic_links = topic_links_tbl,
    annotate_unique_genes = TRUE
  )
  if (pathway_source == "link_scores") {
    method_secondary <- if (identical(link_topic_method, "gene_prob")) "peak_and_gene_prob" else link_topic_method
    report_methods <- unique(c("peak_and_gene", method_secondary))
    if (!is.null(link_topic_report_methods)) {
      report_methods <- unique(as.character(link_topic_report_methods))
    }
    report_methods <- intersect(
      report_methods,
      c("peak_and_gene", "peak_and_gene_prob", "gene_only", "link_score_prob", "link_score_efdr")
    )
    for (method in report_methods) {
      rerun_pathway_from_topic_links(
        out_dir = out_dir,
        topic_links_file = pathway_link_scores_file,
        method = method,
        allow_missing = TRUE,
        title_prefix = title_prefix,
        include_tf = pathway_link_include_tf,
        include_gene = pathway_link_include_gene,
        min_prob = pathway_link_min_prob,
        gene_min_prob = pathway_link_gene_min_prob,
        tf_min_prob = pathway_link_tf_min_prob,
        tf_max_topics = pathway_link_tf_max_topics,
        tf_top_n_per_topic = pathway_link_tf_top_n_per_topic,
        top_n_per_topic = top_n_per_topic,
        dot_top_n_per_topic = dot_top_n_per_topic,
        max_pathways = max_pathways,
        per_comparison = pathway_per_comparison,
        per_comparison_dir = paste0(pathway_per_comparison_dir, "_", .pathway_method_suffix(method)),
        per_comparison_flat = pathway_per_comparison_flat,
        split_direction = pathway_split_direction,
        make_heatmap = pathway_make_heatmap,
        enrichr_sleep_time = pathway_enrichr_sleep_time,
        enrichr_cache_dir = pathway_enrichr_cache_dir,
        enrichr_n_cores = pathway_enrichr_n_cores,
        dbs = pathway_dbs,
        pathway_species = pathway_species,
        pathway_backend = pathway_backend,
        overwrite = pathway_overwrite
      )
    }
  } else {
    plot_topic_pathway_enrichment_heatmap(
      topic_terms = topic_terms,
      edges_docs = edges_docs,
      option_label = option_label,
      out_file = file.path(out_dir, "topic_pathway_enrichment_heatmap.pdf"),
      title_prefix = title_prefix,
      use_all_terms = pathway_use_all_terms,
      make_heatmap = pathway_make_heatmap,
      top_n_per_topic = top_n_per_topic,
      dot_top_n_per_topic = dot_top_n_per_topic,
      max_pathways = max_pathways,
      theta = topic_base$theta,
      tf_link_mode = pathway_tf_link_mode,
      tf_theta_top_n = pathway_tf_top_n_docs,
      tf_theta_min = pathway_tf_min_theta,
      enrichr_sleep_time = pathway_enrichr_sleep_time,
      enrichr_cache_dir = pathway_enrichr_cache_dir,
      enrichr_n_cores = pathway_enrichr_n_cores,
      dbs = pathway_dbs,
      pathway_species = pathway_species,
      pathway_backend = pathway_backend
    )

    if (isTRUE(pathway_per_comparison)) {
      plot_topic_pathway_enrichment_by_comparison_terms(
        topic_terms = topic_terms,
        edges_docs = edges_docs,
        theta = topic_base$theta,
        out_dir = out_dir,
        title_prefix = title_prefix,
        theta_min = pathway_topic_term_theta_min,
        include_peak_terms = !identical(option_label, "opt3_gene_fc_expr"),
        use_all_terms = pathway_use_all_terms,
        per_comparison_dir = pathway_per_comparison_dir,
        per_comparison_flat = pathway_per_comparison_flat,
        split_direction = pathway_split_direction,
        min_genes = 5L,
        top_n_per_topic = top_n_per_topic,
        dot_top_n_per_topic = dot_top_n_per_topic,
        max_pathways = max_pathways,
        make_dotplot = pathway_make_dotplot,
        enrichr_sleep_time = pathway_enrichr_sleep_time,
        enrichr_cache_dir = pathway_enrichr_cache_dir,
        enrichr_n_cores = pathway_enrichr_n_cores,
        dbs = pathway_dbs,
        pathway_species = pathway_species,
        pathway_backend = pathway_backend,
        overwrite = pathway_overwrite,
        doc_design = "comparison"
      )
    }

    if (isTRUE(run_pathway_gsea)) {
      plot_topic_pathway_enrichment_gsea(
        topic_terms = topic_terms,
        edges_docs = edges_docs,
        option_label = option_label,
        out_dir = out_dir,
        theta = topic_base$theta,
        species = gsea_species,
        padj_cut = 0.05,
        min_size = 10L,
        max_size = 500L,
        nperm = gsea_nperm,
        top_n_per_topic = 20L,
        max_pathways = 200L,
        peak_gene_agg = gsea_peak_gene_agg,
        tf_link_mode = pathway_tf_link_mode,
        tf_theta_top_n = pathway_tf_top_n_docs,
        tf_theta_min = pathway_tf_min_theta,
        title_prefix = title_prefix
      )
    }
  }
  plot_intertopic_distance_map(
    phi = topic_base$phi,
    topic_terms = topic_terms,
    out_file = file.path(out_dir, "intertopic_distance_map.pdf"),
    option_label = option_label,
    title_prefix = title_prefix
  )
  list(
    out_dir = out_dir,
    calc_params = calc_params,
    edges_filtered = edges_filt,
    edges_docs = edges_docs,
    doc_term = doc_term,
    dtm = dtm,
    fits_tmp_dir = file.path(out_dir, "tmp_models"),
    metrics = metrics_tbl,
    selection = sel,
    selected_models = selected,
    topic_terms = topic_terms
  )
}

# =============================================================================
# 12) VAE helpers for regulatory topic modeling
# =============================================================================

#' Build a TF cluster map from motif metadata
#'
#' Convert motif metadata with gene-symbol and sub-cluster columns into a named
#' TF-to-cluster vector for Module 3 topic documents.
#'
#' @param motif_path Path to the motif database TSV file, a motif metadata
#'   data.frame, or a list containing `motif_db`.
#'
#' @return A named character vector mapping TF symbols to cluster labels.
#' @noRd
build_tf_cluster_map_from_motif <- function(motif_path) {
  .assert_pkg("cli")
  .assert_pkg("readr")
  .assert_pkg("data.table")

  if (is.list(motif_path) && !is.data.frame(motif_path) && is.data.frame(motif_path$motif_db)) {
    motif_db <- motif_path$motif_db
  } else if (is.data.frame(motif_path)) {
    motif_db <- motif_path
  } else {
    if (!is.character(motif_path) || length(motif_path) != 1L) {
      .log_abort("`motif_path` must be a file path or a motif_db data.frame.")
    }
    if (!file.exists(motif_path)) .log_abort("Missing motif db: {motif_path}")
    motif_db <- readr::read_tsv(motif_path, show_col_types = FALSE)
  }
  if ("HGNC" %in% names(motif_db) && !("gene_symbol" %in% names(motif_db))) {
    motif_db <- dplyr::rename(motif_db, gene_symbol = HGNC)
  }
  motif_dt <- data.table::as.data.table(motif_db)
  if (!("gene_symbol" %in% names(motif_dt)) || !("sub_cluster_name" %in% names(motif_dt))) {
    .log_abort("Motif db missing gene_symbol or sub_cluster_name columns.")
  }

  tf_exclude <- unique(unlist(strsplit(motif_dt[sub_cluster_name == "#N/A", gene_symbol], "::", fixed = TRUE)))
  tf_exclude <- toupper(tf_exclude[!is.na(tf_exclude) & nzchar(tf_exclude)])
  motif_dt <- motif_dt[!is.na(gene_symbol) & nzchar(gene_symbol) & !is.na(sub_cluster_name) & nzchar(sub_cluster_name)]
  motif_dt <- motif_dt[sub_cluster_name != "#N/A"]
  tf_map <- motif_dt[, .(tf = unlist(strsplit(gene_symbol, "::", fixed = TRUE))), by = sub_cluster_name]
  tf_map <- tf_map[!is.na(tf) & nzchar(tf)]
  tf_map[, tf := toupper(tf)]
  tf_best <- tf_map[, .N, by = .(tf, sub_cluster_name)]
  data.table::setorder(tf_best, tf, -N, sub_cluster_name)
  tf_best <- tf_best[, .SD[1], by = tf]
  tf_cluster_map <- stats::setNames(tf_best$sub_cluster_name, tf_best$tf)
  list(tf_cluster_map = tf_cluster_map, tf_exclude = tf_exclude)
}

.resolve_vae_python <- function() {
  py <- Sys.getenv("VAE_PYTHON", unset = "")
  if (!nzchar(py)) {
    default_py <- "/data/homes/yl814/miniconda3/bin/python"
    if (file.exists(default_py)) py <- default_py
  }
  if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) py <- Sys.which("python3")
  if (!nzchar(py)) py <- Sys.which("py")
  if (!nzchar(py)) .log_abort("Python not found on PATH. Set VAE_PYTHON to full path.")
  py
}

.python_has_pyarrow <- function(python) {
  if (!is.character(python) || length(python) != 1L || !nzchar(python)) {
    return(FALSE)
  }
  status <- tryCatch(
    suppressWarnings(system2(python, c("-c", shQuote("import pyarrow")), stdout = FALSE, stderr = FALSE)),
    error = function(e) 1L
  )
  identical(status, 0L)
}

.vae_doc_term_cache_plan <- function(has_r_arrow, python_has_pyarrow, save_full_doc_term_csv = FALSE) {
  write_arrow <- isTRUE(has_r_arrow) && isTRUE(python_has_pyarrow)
  list(
    write_arrow = write_arrow,
    save_full_doc_term_csv = isTRUE(save_full_doc_term_csv) || !write_arrow
  )
}

.normalize_vae_python_variant <- function(vae_variant) {
  variant <- as.character(vae_variant %||% "")
  if (!length(variant) || !nzchar(variant[[1]])) {
    .log_abort("vae_variant must be a non-empty string.")
  }
  variant <- variant[[1]]
  alias_map <- c(
    vmlp = "vae_mlp",
    model_vmlp = "vae_mlp",
    mve = "multivi_encoder",
    model_mve = "multivi_encoder",
    moetm = "moetm_encoder_decoder",
    model_moetm = "moetm_encoder_decoder"
  )
  if (variant %in% names(alias_map)) {
    return(unname(alias_map[[variant]]))
  }
  valid <- c("vae_mlp", "moetm_encoder_decoder", "multivi_encoder")
  if (!variant %in% valid) {
    .log_abort("Unsupported VAE variant: {variant}. Supported variants: {paste(valid, collapse = ', ')}.")
  }
  variant
}

.reset_topic_model_artifacts <- function(out_dir, backend, reuse_if_exists) {
  if (isTRUE(reuse_if_exists)) {
    return(invisible(FALSE))
  }
  backend <- as.character(backend %||% "")
  paths <- c(
    file.path(out_dir, "model_metrics.csv"),
    file.path(out_dir, "model_selection.pdf"),
    file.path(out_dir, "vae_progress.tsv"),
    file.path(out_dir, "vae_train.log"),
    file.path(out_dir, "tmp_models"),
    file.path(out_dir, "vae_models")
  )
  if (identical(backend, "vae")) {
    paths <- c(paths, file.path(out_dir, "vae_models"))
  }
  rds_dir <- file.path(out_dir, "rds")
  if (dir.exists(rds_dir)) {
    paths <- c(
      paths,
      list.files(
        rds_dir,
        pattern = "^(model_metrics|model_selection|theta_|phi_|selected)",
        full.names = TRUE
      )
    )
  }
  existing <- unique(paths[file.exists(paths) | dir.exists(paths)])
  if (length(existing)) {
    unlink(existing, recursive = TRUE, force = TRUE)
  }
  invisible(length(existing) > 0L)
}

run_vae_topic_report_py <- function(doc_term,
                                    edges_docs,
                                    out_dir,
                                    option_label,
                                    direction_by,
                                    vae_script,
                                    k_grid,
                                    vae_variant = "multivi_encoder",
                                    vae_python = NULL,
                                    vae_epochs = 200L,
                                    vae_batch_size = 64L,
                                    vae_hidden = 128L,
                                    vae_lr = 1e-3,
                                    vae_seed = 123L,
                                    vae_device = "auto",
                                    reuse_if_exists = TRUE,
                                    do_report = TRUE,
                                    chosen_K = NULL,
                                    count_input = c("pseudo_count_bin", "pseudo_count_log", "weight"),
                                    save_full_doc_term_csv = FALSE,
                                    topic_report_args = list()) {
  count_input <- match.arg(count_input)
  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("readr")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!nrow(doc_term)) .log_abort("doc_term is empty for VAE; relax filters.")

  count_col <- switch(
    count_input,
    pseudo_count_bin = "pseudo_count_bin",
    pseudo_count_log = "pseudo_count_log",
    weight = "weight"
  )
  if (!count_col %in% names(doc_term)) {
    .log_abort("doc_term missing required column: {count_col}")
  }
  t_cache <- proc.time()[["elapsed"]]
  doc_term$count <- doc_term[[count_col]]
  if (is.null(vae_python) || !nzchar(vae_python)) vae_python <- .resolve_vae_python()
  has_r_arrow <- requireNamespace("arrow", quietly = TRUE)
  cache_plan <- .vae_doc_term_cache_plan(
    has_r_arrow = has_r_arrow,
    python_has_pyarrow = .python_has_pyarrow(vae_python),
    save_full_doc_term_csv = save_full_doc_term_csv
  )
  cache_paths <- if (exists("write_doc_term_cache", mode = "function")) {
    write_doc_term_cache(
      doc_term,
      out_dir = out_dir,
      save_full_doc_term_csv = cache_plan$save_full_doc_term_csv,
      write_arrow = cache_plan$write_arrow
    )
  } else {
    data.table::fwrite(utils::head(doc_term, 100L), file.path(out_dir, "doc_term_first100.csv"))
    csv_path <- file.path(out_dir, "doc_term.csv")
    data.table::fwrite(doc_term, csv_path)
    list(arrow = NA_character_, preview = file.path(out_dir, "doc_term_first100.csv"), csv = csv_path)
  }
  .log_inform("VAE doc-term cache writing finished in {round(proc.time()[['elapsed']] - t_cache, 1)} sec for {out_dir}.")
  t_rds <- proc.time()[["elapsed"]]
  .save_all(out_dir, "doc_term", doc_term)
  if (!is.null(edges_docs) && nrow(edges_docs)) {
    .save_all(out_dir, "edges_docs", edges_docs)
  }
  .log_inform("VAE RDS cache writing finished in {round(proc.time()[['elapsed']] - t_rds, 1)} sec for {out_dir}.")
  t_dtm <- proc.time()[["elapsed"]]
  doc_term$pseudo_count <- doc_term[[count_col]]
  dtm_obj <- build_sparse_dtm(doc_term, count_col = "pseudo_count")
  dtm <- dtm_obj$dtm
  .save_all(out_dir, "dtm", dtm)
  .save_all(out_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))
  .log_inform("VAE sparse DTM cache writing finished in {round(proc.time()[['elapsed']] - t_dtm, 1)} sec for {out_dir}.")

  if (is.null(k_grid) || !length(k_grid)) .log_abort("k_grid required for VAE.")
  k_grid <- sort(unique(as.integer(k_grid)))
  k_grid <- k_grid[is.finite(k_grid)]
  if (!length(k_grid)) .log_abort("k_grid must contain at least one finite K value.")
  vae_python_variant <- .normalize_vae_python_variant(vae_variant)

  metrics_path <- file.path(out_dir, "model_metrics.csv")
  manifest_path <- file.path(out_dir, "vae_model_manifest.csv")
  models_dir <- file.path(out_dir, "vae_models")
  complete_k <- integer()
  old_metrics <- data.table::data.table()
  old_manifest <- data.table::data.table()
  if (file.exists(metrics_path)) {
    old_metrics <- tryCatch(data.table::fread(metrics_path, showProgress = FALSE), error = function(e) data.table::data.table())
  }
  if (file.exists(manifest_path)) {
    old_manifest <- tryCatch(data.table::fread(manifest_path, showProgress = FALSE), error = function(e) data.table::data.table())
  }
  if (nrow(old_metrics) && "K" %in% names(old_metrics) && dir.exists(models_dir)) {
    metrics_k <- unique(as.integer(old_metrics$K))
    file_k <- k_grid[
      file.exists(file.path(models_dir, sprintf("theta_K%d.csv", k_grid))) &
        file.exists(file.path(models_dir, sprintf("phi_K%d.csv", k_grid)))
    ]
    complete_k <- intersect(metrics_k, file_k)
  }
  run_k_grid <- if (isTRUE(reuse_if_exists)) setdiff(k_grid, complete_k) else k_grid
  reuse_ok <- isTRUE(reuse_if_exists) && !length(run_k_grid)

  if (!reuse_ok) {
    doc_term_input <- if (isTRUE(cache_plan$write_arrow)) cache_paths$arrow else cache_paths$csv
    if (!is.character(doc_term_input) || !length(doc_term_input) || is.na(doc_term_input) || !file.exists(doc_term_input)) {
      .log_abort("Missing doc-term cache for VAE: {doc_term_input}")
    }
    k_grid_txt <- paste(run_k_grid, collapse = ",")
    if (isTRUE(reuse_if_exists) && length(complete_k)) {
      .log_inform("Reusing VAE outputs for K={paste(sort(complete_k), collapse = ',')}; training missing K={k_grid_txt}.")
    }
    py_args <- c(
      vae_script,
      "--doc-term", doc_term_input,
      "--out-dir", out_dir,
      "--k-grid", k_grid_txt,
      "--epochs", as.character(vae_epochs),
      "--batch-size", as.character(vae_batch_size),
      "--hidden", as.character(vae_hidden),
      "--lr", as.character(vae_lr),
      "--seed", as.character(vae_seed),
      "--device", vae_device,
      "--variant", vae_python_variant,
      "--progress-log", file.path(out_dir, "vae_progress.tsv")
    )
    .log_inform("Running VAE Python with doc-term input {basename(doc_term_input)}; K grid: {k_grid_txt}.")
    t_py <- proc.time()[["elapsed"]]
    py_out <- tryCatch(
      system2(vae_python, py_args, stdout = TRUE, stderr = TRUE),
      error = function(e) .log_abort(c("Failed to run python for VAE.", i = conditionMessage(e)))
    )
    .log_inform("VAE Python process finished in {round(proc.time()[['elapsed']] - t_py, 1)} sec for {out_dir}.")
    writeLines(py_out, file.path(out_dir, "vae_train.log"))
    status <- attr(py_out, "status")
    if (!is.null(status) && is.finite(status) && status != 0) {
      .log_abort(c("VAE python exited with non-zero status.", i = paste0("status=", status)))
    }
    if (isTRUE(reuse_if_exists) && length(complete_k)) {
      new_metrics <- if (file.exists(metrics_path)) {
        data.table::fread(metrics_path, showProgress = FALSE)
      } else {
        data.table::data.table()
      }
      if (nrow(old_metrics) && nrow(new_metrics) && "K" %in% names(old_metrics) && "K" %in% names(new_metrics)) {
        old_metrics <- old_metrics[!as.integer(K) %in% as.integer(new_metrics$K)]
        merged_metrics <- data.table::rbindlist(list(old_metrics, new_metrics), use.names = TRUE, fill = TRUE)
        data.table::setorder(merged_metrics, K)
        data.table::fwrite(merged_metrics, metrics_path)
      }
      new_manifest <- if (file.exists(manifest_path)) {
        data.table::fread(manifest_path, showProgress = FALSE)
      } else {
        data.table::data.table()
      }
      if (nrow(old_manifest) && nrow(new_manifest) && "K" %in% names(old_manifest) && "K" %in% names(new_manifest)) {
        old_manifest <- old_manifest[!as.integer(K) %in% as.integer(new_manifest$K)]
        merged_manifest <- data.table::rbindlist(list(old_manifest, new_manifest), use.names = TRUE, fill = TRUE)
        data.table::setorder(merged_manifest, K)
        data.table::fwrite(merged_manifest, manifest_path)
      }
    }
  } else {
    .log_inform("Reusing existing VAE outputs for K={paste(k_grid, collapse = ',')}; skipping training for {out_dir}")
  }

  if (!file.exists(metrics_path)) .log_abort("VAE did not produce model_metrics.csv in {out_dir}")
  metrics_tbl <- data.table::as.data.table(readr::read_csv(metrics_path, show_col_types = FALSE))
  metrics_tbl[, `:=`(
    count_input_effective = count_input,
    n_model_tokens = as.double(sum(.safe_num(doc_term[[count_col]]), na.rm = TRUE))
  )]
  data.table::fwrite(metrics_tbl, metrics_path)
  .save_all(out_dir, "model_metrics", metrics_tbl)

  title_prefix <- .topic_model_selection_title(out_dir, backend_label = paste("VAE", vae_variant), vae_variant = vae_variant)
  sel <- plot_model_selection_cistopic(metrics_tbl, file.path(out_dir, "model_selection.pdf"), title_prefix = title_prefix)
  .save_all(out_dir, "model_selection", sel)

  if (!isTRUE(do_report)) {
    return(invisible(TRUE))
  }

  m <- sel$table
  idx <- sel$idx
  pick_K <- function(ix) {
    if (!is.finite(ix)) return(NULL)
    m$K[ix]
  }
  if (is.null(chosen_K)) {
    chosen_K <- pick_K(idx$perplexity)
  }
  if (is.null(chosen_K)) chosen_K <- pick_K(idx$maximum)
  if (is.null(chosen_K)) chosen_K <- m$K[1]

  theta_path <- file.path(out_dir, "vae_models", sprintf("theta_K%d.csv", chosen_K))
  phi_path <- file.path(out_dir, "vae_models", sprintf("phi_K%d.csv", chosen_K))
  if (!file.exists(theta_path) || !file.exists(phi_path)) {
    .log_abort("Missing theta/phi outputs for K={chosen_K} in {out_dir}")
  }
  theta_df <- readr::read_csv(theta_path, show_col_types = FALSE)
  phi_df <- readr::read_csv(phi_path, show_col_types = FALSE)

  theta <- as.matrix(theta_df[, -1, drop = FALSE])
  rownames(theta) <- theta_df[[1]]
  phi <- as.matrix(phi_df[, -1, drop = FALSE])
  rownames(phi) <- phi_df[[1]]

  topic_base <- list(K = chosen_K, theta = theta, phi = phi, metrics = metrics_tbl[metrics_tbl$K == chosen_K, ])

  defaults <- list(
    binarize_method = "gammafit",
    topic_score_method = "normtop_specificity",
    thrP = 0.9,
    top_n_terms = 500L,
    in_topic_min_terms = 1L,
    pathway_use_all_terms = FALSE,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
    top_n_per_topic = 100L,
    max_pathways = 1000L,
    pathway_tf_link_mode = "theta",
    pathway_tf_top_n_docs = 50L,
    pathway_tf_min_theta = NA_real_,
    run_pathway_gsea = FALSE,
    gsea_species = "Homo sapiens",
    gsea_nperm = 1000L,
    gsea_peak_gene_agg = "max",
    pathway_source = "topic_terms",
    pathway_link_scores_file = NULL,
    pathway_link_scores_file_tf = NULL,
    pathway_link_gene_terms_file = NULL,
    pathway_link_min_prob = 0,
    pathway_link_include_tf = TRUE,
    pathway_link_include_gene = TRUE,
    pathway_link_gene_min_prob = 0,
    pathway_link_tf_min_prob = 0.5,
    pathway_link_tf_max_topics = 5L,
    pathway_link_tf_top_n_per_topic = 30L,
    pathway_per_comparison = TRUE,
    pathway_per_comparison_dir = ".",
    pathway_per_comparison_flat = TRUE,
    pathway_split_direction = TRUE,
    run_link_topic_scores = FALSE,
    link_topic_gate_mode = c("none", "peak_and_gene_in_set"),
    link_topic_top_k = 3L,
    link_topic_min_prob = 0,
    link_topic_include_tf = FALSE,
    link_topic_chunk_size = 5000L,
    link_topic_n_cores = 1L,
    link_topic_overwrite = FALSE,
    link_topic_method = "gene_prob",
    link_topic_prob_cutoff = 0.3,
    link_topic_fdr_q = 0.2,
    link_topic_fdr_p = NA_real_,
    link_topic_efdr_scope = "per_topic",
    link_topic_efdr_B = 100L,
    link_topic_efdr_seed = 1L
  )
  args <- modifyList(defaults, topic_report_args)
  do.call(
    run_tfdocs_report_from_topic_base,
    c(list(
      topic_base = topic_base,
      dtm = dtm,
      edges_docs = edges_docs,
      out_dir = out_dir,
      option_label = option_label,
      direction_by = direction_by,
      title_prefix = title_prefix
    ), args)
  )

  invisible(TRUE)
}

run_vae_ctf_multivi <- function(edges_all,
                                out_root,
                                celllines = NULL,
                                tf_cluster_map,
                                tf_exclude = NULL,
                                k_grid_default,
                                k_single_map = list(),
                                abs_log2fc_fp_min = 0,
                                abs_delta_fp_min = 1,
                                abs_log2fc_gene_min = 1,
                                require_fp_bound_either = TRUE,
                                require_tf_expr_either = TRUE,
                                require_gene_expr_either = TRUE,
                                direction_consistency = "aligned",
                                top_terms_per_doc = Inf,
                                min_df = 2,
                                count_method = "bin",
                                count_scale = 50,
                                binarize_method = "gammafit",
                                thrP = 0.9,
                                top_n_terms = 500L,
                                in_topic_min_terms = 1,
                                topic_report_args = list()) {
  .assert_pkg("data.table")
  .assert_pkg("cli")

  edges_dt <- data.table::as.data.table(edges_all)
  if (!("comparison_id" %in% names(edges_dt))) .log_abort("edges_all missing comparison_id.")
  if (!("cellline" %in% names(edges_dt))) {
    edges_dt[, cellline := sub("_.*$", "", comparison_id)]
  }
  if (is.null(celllines)) {
    celllines <- sort(unique(as.character(edges_dt$cellline)))
  }
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  vae_script <- Sys.getenv("VAE_SCRIPT", unset = "")
  if (!nzchar(vae_script)) {
    vae_script <- system.file("python", "logistic_normal_vae_topics.py", package = "craftgrn")
  }
  if (!nzchar(vae_script) || !file.exists(vae_script)) {
    cand <- file.path("dev", "logistic_normal_vae_topics.py")
    if (file.exists(cand)) vae_script <- cand
  }
  if (!file.exists(vae_script)) .log_abort("Missing VAE script: {vae_script}")

  for (cell in celllines) {
    edges_sub <- edges_dt[cellline == cell]
    if (!nrow(edges_sub)) next
    if (!is.null(tf_exclude) && length(tf_exclude)) {
      edges_sub <- edges_sub[!toupper(tf) %in% tf_exclude]
    }
    if (!nrow(edges_sub)) next

    edges_filt <- filter_edges_for_tf_topics(
      edges_sub,
      abs_log2fc_fp_min = abs_log2fc_fp_min,
      abs_delta_fp_min = abs_delta_fp_min,
      abs_log2fc_gene_min = abs_log2fc_gene_min,
      require_fp_bound_either = require_fp_bound_either,
      require_tf_expr_either = require_tf_expr_either,
      require_gene_expr_either = require_gene_expr_either,
      direction_consistency = direction_consistency
    )
    if (!nrow(edges_filt)) next

    edges_docs <- add_tf_docs(
      edges_filt,
      doc_mode = "tf_cluster",
      direction_by = "gene",
      tf_cluster_map = tf_cluster_map
    )
    doc_term <- build_doc_term_joint(
      edges_docs,
      weight_type_peak = "delta_fp",
      weight_type_gene = "fc_mag_gene",
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      count_method = count_method,
      count_scale = count_scale,
      distinct_terms = TRUE,
      balance_mode = "min",
      prefix_terms = TRUE
    )
    if (!nrow(doc_term)) {
      .log_inform("Skipping VAE joint ctf_docs: no doc_term for {cell}")
      next
    }

    local_topic_args <- modifyList(list(
      binarize_method = binarize_method,
      thrP = thrP,
      top_n_terms = top_n_terms,
      in_topic_min_terms = in_topic_min_terms
    ), topic_report_args)

    out_grid <- file.path(
      out_root,
      paste0(cell, "_vae_joint_ctf_docs_peak_delta_fp_gene_fc_expr_multivi_encoder_Kgrid_default")
    )
    run_vae_topic_report_py(
      doc_term = doc_term,
      edges_docs = edges_docs,
      out_dir = out_grid,
      option_label = "joint",
      direction_by = "gene",
      vae_script = vae_script,
      k_grid = k_grid_default,
      vae_variant = "multivi_encoder",
      do_report = TRUE,
      topic_report_args = local_topic_args
    )

    if (cell %in% names(k_single_map)) {
      k_vals <- as.integer(k_single_map[[cell]])
      k_vals <- k_vals[is.finite(k_vals)]
      k_vals <- unique(k_vals)
      for (k_val in k_vals) {
        out_single <- file.path(
          out_root,
          paste0(cell, "_vae_joint_ctf_docs_peak_delta_fp_gene_fc_expr_multivi_encoder_K", k_val)
        )
        run_vae_topic_report_py(
          doc_term = doc_term,
          edges_docs = edges_docs,
          out_dir = out_single,
          option_label = "joint",
          direction_by = "gene",
          vae_script = vae_script,
          k_grid = c(k_val),
          vae_variant = "multivi_encoder",
          do_report = TRUE,
          topic_report_args = local_topic_args
        )
      }
    }
  }

  invisible(TRUE)
}

run_vae_doc_topic_heatmaps <- function(topic_root,
                                       backend = c("vae", "warplda"),
                                       vae_variant = "multivi_encoder",
                                       doc_mode = c("tf_cluster", "tf")) {
  .assert_pkg("data.table")
  if (!exists("plot_tf_topic_heatmaps_from_link_scores")) return(invisible(NULL))
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  out_dirs <- unique(c(topic_root, list.dirs(topic_root, recursive = FALSE, full.names = TRUE)))
  if (backend == "vae") {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_"), basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_"), basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(topic_root, "topic_links.csv"))) {
    out_dirs <- topic_root
  }
  for (d in out_dirs) {
    topic_links_path <- file.path(d, "topic_links.csv")
    if (!file.exists(topic_links_path)) next
    topic_links <- data.table::fread(topic_links_path)
    for (method in c("peak_and_gene", "peak_and_gene_prob")) {
      link_scores <- .topic_links_to_link_scores(topic_links, method = method)
      if (!nrow(link_scores)) next
      out_base <- file.path(d, paste0("dthm_ls_", .short_link_method_tag(method)))
      plot_tf_topic_heatmaps_from_link_scores(
        link_scores = link_scores,
        out_dir = out_base,
        title_prefix = paste("Link scores", method, basename(d), sep = " | "),
        value_col = "prob",
        min_value = 0,
        per_comparison = TRUE,
        split_direction = TRUE
      )
    }
  }
  invisible(TRUE)
}

plot_topic_delta_networks_from_link_scores <- function(link_scores,
                                                       step2_out_dir,
                                                       out_root,
                                                       min_prob = 0.5,
                                                       filter_same_direction = TRUE,
                                                       gene_fc_thresh = 1.5,
                                                       size_by = "expr_max",
                                                       motif_db = "jaspar2024") {
  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("readr")
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    .log_inform("Skipping topic delta network plots: htmlwidgets not installed.")
    return(invisible(NULL))
  }
  plot_tf_network_delta_fn <- get0("plot_tf_network_delta", mode = "function")
  if (is.null(plot_tf_network_delta_fn)) {
    cand <- file.path("R", "utils_plot_tf_network_delta.R")
    if (!file.exists(cand)) .log_abort("Missing utils_plot_tf_network_delta.R")
    source(cand)
    plot_tf_network_delta_fn <- get0("plot_tf_network_delta", mode = "function")
  }
  if (is.null(plot_tf_network_delta_fn)) {
    .log_abort("Missing plot_tf_network_delta() after sourcing network helper.")
  }

  link_dt <- data.table::as.data.table(link_scores)
  if (!nrow(link_dt)) return(invisible(NULL))
  if (!("topic_num" %in% names(link_dt))) {
    if ("topic" %in% names(link_dt)) {
      link_dt[, topic_num := as.integer(gsub("^Topic", "", topic))]
    } else {
      .log_abort("link_scores missing topic_num/topic.")
    }
  }
  if (!("prob" %in% names(link_dt))) link_dt[, prob := 1]
  doc_info <- .parse_doc_id(link_dt$doc_id)
  link_dt <- cbind(link_dt, doc_info)
  link_dt[, direction_group := .direction_group(direction)]
  link_dt <- link_dt[!is.na(comparison_id) & nzchar(comparison_id)]
  if (!nrow(link_dt)) return(invisible(NULL))

  gene_col_link <- if ("gene_key" %in% names(link_dt)) "gene_key" else if ("gene" %in% names(link_dt)) "gene" else "target_gene"
  tf_col_link <- if ("tf" %in% names(link_dt)) "tf" else "TF"
  if (!gene_col_link %in% names(link_dt) || !tf_col_link %in% names(link_dt)) {
    .log_abort("link_scores missing tf/gene columns.")
  }

  detect_mapping <- function(df_names) {
    sc <- df_names[grepl("^link_score_", df_names)]
    if (length(sc) < 2L) return(NULL)
    sc <- sc[seq_len(min(2L, length(sc)))]
    suf <- sub("^link_score_", "", sc)
    ctrl_idx <- if (grepl("10_fbs|control|baseline|ctrl", tolower(suf[1])) &&
      !grepl("10_fbs|control|baseline|ctrl", tolower(suf[2]))) 1L else
        if (!grepl("10_fbs|control|baseline|ctrl", tolower(suf[1])) &&
          grepl("10_fbs|control|baseline|ctrl", tolower(suf[2]))) 2L else 2L
    str_idx <- if (ctrl_idx == 1L) 2L else 1L
    list(
      score_ctrl_col = sc[ctrl_idx],
      score_str_col  = sc[str_idx],
      sign_ctrl_col  = paste0("link_sign_", suf[ctrl_idx]),
      sign_str_col   = paste0("link_sign_", suf[str_idx]),
      tf_expr_ctrl_col   = paste0("tf_expr_",    suf[ctrl_idx]),
      tf_expr_str_col    = paste0("tf_expr_",    suf[str_idx]),
      gene_expr_ctrl_col = paste0("gene_expr_",  suf[ctrl_idx]),
      gene_expr_str_col  = paste0("gene_expr_",  suf[str_idx])
    )
  }

  delta_cache <- new.env(parent = emptyenv())
  fetch_delta_links <- function(comp) {
    if (exists(comp, envir = delta_cache, inherits = FALSE)) {
      return(get(comp, envir = delta_cache))
    }
    source_type <- "raw"
    cand <- file.path(step2_out_dir, paste0(comp, "_filtered_links.csv"))
    if (!file.exists(cand)) {
      cand <- file.path(step2_out_dir, paste0(comp, "_delta_links_filtered.csv"))
    }
    if (file.exists(cand)) {
      source_type <- "filtered"
      df <- readr::read_csv(cand, show_col_types = FALSE)
    } else {
      cand <- file.path(step2_out_dir, paste0(comp, "_delta_links.csv"))
      if (!file.exists(cand)) {
        cand_up <- file.path(step2_out_dir, paste0(comp, "_filtered_links_up.csv"))
        cand_down <- file.path(step2_out_dir, paste0(comp, "_filtered_links_down.csv"))
        if (!file.exists(cand_up) && !file.exists(cand_down)) {
          cand_up <- file.path(step2_out_dir, paste0(comp, "_delta_links_filtered_up.csv"))
          cand_down <- file.path(step2_out_dir, paste0(comp, "_delta_links_filtered_down.csv"))
        }
        have_up <- file.exists(cand_up)
        have_down <- file.exists(cand_down)
        if (!have_up && !have_down) {
          .log_inform("Missing delta links for comparison: {comp}")
          assign(comp, NULL, envir = delta_cache)
          return(NULL)
        }
        source_type <- "filtered_split"
        dfs <- list()
        if (have_up) {
          df_up <- readr::read_csv(cand_up, show_col_types = FALSE)
          df_up$direction_group <- "up"
          dfs <- c(dfs, list(df_up))
        }
        if (have_down) {
          df_down <- readr::read_csv(cand_down, show_col_types = FALSE)
          df_down$direction_group <- "down"
          dfs <- c(dfs, list(df_down))
        }
        df <- dplyr::bind_rows(dfs)
      } else {
        source_type <- "raw"
        df <- readr::read_csv(cand, show_col_types = FALSE)
      }
    }

    # If filtered tables are missing plotting columns, enrich from raw while
    # preserving the filtered edge set via key intersection.
    if (sum(grepl("^link_score_", names(df))) < 2L &&
      !identical(source_type, "raw")) {
      raw_path <- file.path(step2_out_dir, paste0(comp, "_delta_links.csv"))
      if (file.exists(raw_path)) {
        raw_df <- readr::read_csv(raw_path, show_col_types = FALSE)
        fdt <- data.table::as.data.table(df)
        rdt <- data.table::as.data.table(raw_df)

        pick_col <- function(nm, choices) {
          hit <- intersect(choices, nm)
          if (!length(hit)) return(NA_character_)
          hit[[1]]
        }
        tf_f <- pick_col(names(fdt), c("tf", "TF"))
        gn_f <- pick_col(names(fdt), c("gene_key", "gene", "target_gene"))
        pk_f <- pick_col(names(fdt), c("peak_id", "peak_ID", "fp_peak"))
        tf_r <- pick_col(names(rdt), c("tf", "TF"))
        gn_r <- pick_col(names(rdt), c("gene_key", "gene", "target_gene"))
        pk_r <- pick_col(names(rdt), c("peak_id", "peak_ID", "fp_peak"))

        if (all(is.character(c(tf_f, gn_f, tf_r, gn_r))) &&
          all(nzchar(c(tf_f, gn_f, tf_r, gn_r)))) {
          fkey <- data.table::data.table(
            tf_key = as.character(fdt[[tf_f]]),
            gene_key = as.character(fdt[[gn_f]])
          )
          rdt[, tf_key := as.character(get(tf_r))]
          rdt[, gene_key := as.character(get(gn_r))]
          use_peak <- is.character(pk_f) && nzchar(pk_f) && is.character(pk_r) && nzchar(pk_r)
          by_cols <- c("tf_key", "gene_key")
          if (isTRUE(use_peak)) {
            fkey[, peak_key := as.character(fdt[[pk_f]])]
            rdt[, peak_key := as.character(get(pk_r))]
            by_cols <- c(by_cols, "peak_key")
          }
          fkey <- unique(fkey[!is.na(tf_key) & nzchar(tf_key) & !is.na(gene_key) & nzchar(gene_key)])
          keep_cols <- unique(c(by_cols, names(df), names(raw_df)))
          keep_cols <- intersect(keep_cols, names(rdt))
          if (nrow(fkey) && length(keep_cols)) {
            rsub <- unique(rdt[, ..keep_cols])
            merged <- data.table::merge.data.table(
              fkey,
              rsub,
              by = by_cols,
              all.x = FALSE,
              all.y = FALSE
            )
            if (nrow(merged)) {
              if ("tf_key" %in% names(merged)) merged[, tf_key := NULL]
              if ("gene_key" %in% names(merged)) merged[, gene_key := NULL]
              if ("peak_key" %in% names(merged)) merged[, peak_key := NULL]
              df <- as.data.frame(merged)
              source_type <- "filtered_enriched_raw"
            }
          }
        }
      }
    }

    if ("active_any" %in% names(df)) {
      df <- df[which(isTRUE(df$active_any) | df$active_any %in% c(1, "1", "TRUE", "true")), , drop = FALSE]
    }
    already_filtered <- identical(source_type, "filtered") ||
      identical(source_type, "filtered_split") ||
      identical(source_type, "filtered_enriched_raw")
    if (isTRUE(filter_same_direction) && !isTRUE(already_filtered) &&
      all(c("log2FC_tf_expr", "log2FC_gene_expr") %in% names(df))) {
      same_dir <- is.finite(df$log2FC_tf_expr) & is.finite(df$log2FC_gene_expr) &
        (df$log2FC_tf_expr * df$log2FC_gene_expr > 0)
      df <- df[which(same_dir), , drop = FALSE]
    }
    if (isTRUE(filter_same_direction) && isTRUE(already_filtered)) {
      .log_inform("Skipping same-direction filter for {comp}: using {source_type} delta links.")
    }

    assign(comp, df, envir = delta_cache)
    df
  }

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  comps <- unique(link_dt$comparison_id)
  topics_all <- sort(unique(link_dt$topic_num))
  summary_rows <- list()
  for (cmp in comps) {
    cmp_delta <- fetch_delta_links(cmp)
    if (is.null(cmp_delta) || !nrow(cmp_delta)) next
    ns <- names(cmp_delta)
    mp <- detect_mapping(ns)
    if (is.null(mp)) {
      .log_warn("Skipping comparison {cmp}: missing required link_score_* columns in delta links.")
      next
    }

    gene_col <- if ("gene_key" %in% ns) "gene_key" else if ("gene" %in% ns) "gene" else "target_gene"
    if (!gene_col %in% ns) next
    tf_col <- if ("tf" %in% ns) "tf" else "TF"
    if (!tf_col %in% ns) next
    peak_col <- if ("peak_id" %in% ns) {
      "peak_id"
    } else if ("peak_ID" %in% ns) {
      "peak_ID"
    } else if ("fp_peak" %in% ns) {
      "fp_peak"
    } else {
      NA_character_
    }

    for (dir_lab in unique(link_dt[comparison_id == cmp, direction_group])) {
      if (!nzchar(dir_lab)) next
      for (topic_id in topics_all) {
        link_sub <- link_dt[
          comparison_id == cmp & direction_group == dir_lab & topic_num == topic_id & prob >= min_prob
        ]
        link_sub <- link_sub[!is.na(get(tf_col_link)) & nzchar(get(tf_col_link)) &
          !is.na(get(gene_col_link)) & nzchar(get(gene_col_link))]
        if (!nrow(link_sub)) next
        link_pairs <- unique(link_sub[, .(
          tf = get(tf_col_link),
          gene_key = get(gene_col_link),
          peak_id = if ("peak_id" %in% names(link_sub)) as.character(peak_id) else NA_character_,
          prob,
          score
        )])
        if (!nrow(link_pairs)) next

        cmp_dt <- data.table::as.data.table(cmp_delta)
        pair_dt <- data.table::data.table(
          tf = link_pairs$tf,
          gene = link_pairs$gene_key,
          peak = link_pairs$peak_id
        )
        join_cols <- c(tf_col, gene_col)
        from_cols <- c("tf", "gene")
        if (is.character(peak_col) && nzchar(peak_col)) {
          pair_dt <- pair_dt[!is.na(peak) & nzchar(peak)]
          if (nrow(pair_dt)) {
            join_cols <- c(join_cols, peak_col)
            from_cols <- c(from_cols, "peak")
          }
        }
        pair_dt <- unique(pair_dt[, ..from_cols])
        data.table::setnames(pair_dt, from_cols, join_cols)
        sub_links <- data.table::merge.data.table(
          cmp_dt,
          pair_dt,
          by = join_cols,
          all = FALSE
        )
        if (!nrow(sub_links)) next

        out_dir <- file.path(out_root, .safe_filename(cmp), .safe_filename(dir_lab))
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        plot_title <- paste(cmp, dir_lab, paste0("Topic", topic_id), sep = " | ")
        args <- list(
          data = sub_links,
          plot_title = plot_title,
          layout_algo = "fr",
          physics = TRUE,
          add_direct = TRUE,
          edge_filter_min = 0,
          min_delta_abs = 0,
          keep_top_edges_per_tf = 6000,
          peak_mode = "show_all",
          show_peaks = FALSE,
          gene_fc_thresh = gene_fc_thresh,
          de_reference = "str_over_ctrl",
          motif_db = motif_db,
          score_ctrl_col = if (mp$score_ctrl_col %in% ns) mp$score_ctrl_col else NULL,
          score_str_col = if (mp$score_str_col %in% ns) mp$score_str_col else NULL,
          sign_ctrl_col = if (mp$sign_ctrl_col %in% ns) mp$sign_ctrl_col else NULL,
          sign_str_col = if (mp$sign_str_col %in% ns) mp$sign_str_col else NULL,
          tf_expr_ctrl_col = if (mp$tf_expr_ctrl_col %in% ns) mp$tf_expr_ctrl_col else NULL,
          tf_expr_str_col = if (mp$tf_expr_str_col %in% ns) mp$tf_expr_str_col else NULL,
          gene_expr_ctrl_col = if (mp$gene_expr_ctrl_col %in% ns) mp$gene_expr_ctrl_col else NULL,
          gene_expr_str_col = if (mp$gene_expr_str_col %in% ns) mp$gene_expr_str_col else NULL,
          size_by = size_by
        )
        w <- try(do.call(plot_tf_network_delta_fn, args), silent = TRUE)
        if (inherits(w, "try-error")) next
        out_html <- file.path(out_dir, paste0("Topic", topic_id, ".html"))
        htmlwidgets::saveWidget(w, out_html, selfcontained = TRUE)
        set_html_title_fn <- get0(".set_html_title", mode = "function")
        if (!is.null(set_html_title_fn)) {
          set_html_title_fn(out_html, plot_title)
        }

        tf_vals <- unique(as.character(sub_links[[tf_col]]))
        gene_vals <- unique(as.character(sub_links[[gene_col]]))
        tf_vals <- tf_vals[!is.na(tf_vals) & nzchar(tf_vals)]
        gene_vals <- gene_vals[!is.na(gene_vals) & nzchar(gene_vals)]
        link_pairs_unique <- unique(data.table::data.table(
          tf = sub_links[[tf_col]],
          gene = sub_links[[gene_col]]
        ))
        link_pairs_unique <- unique(data.table::data.table(
          tf = sub_links[[tf_col]],
          gene = sub_links[[gene_col]]
        ))
        link_pairs_unique <- link_pairs_unique[!is.na(tf) & !is.na(gene) & nzchar(tf) & nzchar(gene)]
        link_list <- if (nrow(link_pairs_unique)) {
          paste(paste0(link_pairs_unique$tf, "::", link_pairs_unique$gene), collapse = ",")
        } else {
          ""
        }
        summary_rows[[length(summary_rows) + 1L]] <- data.table::data.table(
          comparison = cmp,
          direction = dir_lab,
          topic_id = as.integer(topic_id),
          link_count = nrow(link_pairs_unique),
          peak_count = if ("peak_id" %in% names(sub_links)) length(unique(sub_links$peak_id)) else NA_integer_,
          total_target_gene_count = length(unique(gene_vals)),
          tf_count = length(unique(tf_vals)),
          gene_list = paste(sort(unique(gene_vals)), collapse = ","),
          tf_list = paste(sort(unique(tf_vals)), collapse = ","),
          link_list = link_list
        )
      }
    }
  }

  if (length(summary_rows)) {
    summary_tbl <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
    data.table::fwrite(summary_tbl, file.path(out_root, "topic_sub_network_summary.csv"))
  }

  invisible(TRUE)
}

run_vae_topic_delta_network_plots <- function(topic_root,
                                              step2_out_dir,
                                              min_prob = 0.5,
                                              filter_same_direction = TRUE,
                                              methods = c("peak_and_gene", "peak_and_gene_prob", "gene_only"),
                                              backend = c("vae", "warplda"),
                                              vae_variant = "multivi_encoder",
                                              doc_mode = c("tf_cluster", "tf")) {
  out_dirs <- unique(c(topic_root, list.dirs(topic_root, recursive = FALSE, full.names = TRUE)))
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  methods <- unique(as.character(methods))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  if (!length(methods)) .log_abort("methods must include at least one plotting mode.")
  allowed_methods <- c("peak_and_gene", "peak_and_gene_prob", "gene_only")
  bad_methods <- setdiff(methods, allowed_methods)
  if (length(bad_methods)) {
    .log_abort("Unknown methods in run_vae_topic_delta_network_plots: {paste(bad_methods, collapse = ', ')}")
  }
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  if (backend == "vae") {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_"), basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_"), basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(topic_root, "topic_links.csv"))) {
    out_dirs <- topic_root
  }
  if (!length(out_dirs)) return(invisible(NULL))
  for (d in out_dirs) {
    topic_links_path <- file.path(d, "topic_links.csv")
    if (!file.exists(topic_links_path)) next
    topic_links <- data.table::fread(topic_links_path)
    for (method in methods) {
      link_dt <- .topic_links_to_link_scores(topic_links, method = method)
      if (!nrow(link_dt)) next
      subnetwork_dir <- switch(
        as.character(method),
        peak_and_gene = "subnet_peak_gene",
        peak_and_gene_prob = "subnet_peak_gene_prob",
        paste0("subnet_", .short_link_method_tag(method))
      )
      out_root <- file.path(d, subnetwork_dir)
      plot_topic_delta_networks_from_link_scores(
        link_scores = link_dt,
        step2_out_dir = step2_out_dir,
        out_root = out_root,
        min_prob = min_prob,
        filter_same_direction = filter_same_direction
      )
    }
  }
  invisible(TRUE)
}

run_vae_topic_delta_network_pathway <- function(topic_root,
                                                backend = c("vae", "warplda"),
                                                vae_variant = "multivi_encoder",
                                                top_n_per_topic = Inf,
                                                dot_top_n_per_topic = Inf,
                                                max_pathways = Inf,
                                                per_comparison = TRUE,
                                                split_direction = TRUE,
                                                enrichr_sleep_time = 0,
                                                enrichr_cache_dir = NULL,
                                                enrichr_n_cores = NULL,
                                                pathway_backend = NULL,
                                                doc_mode = c("tf_cluster", "tf")) {
  out_dirs <- unique(c(topic_root, list.dirs(topic_root, recursive = FALSE, full.names = TRUE)))
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  enrichr_sleep_time <- .normalize_enrichr_sleep_time(enrichr_sleep_time)
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- file.path(topic_root, "cache", "enrichr")
  }
  enrichr_n_cores <- .normalize_enrichr_n_cores(enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  if (backend == "vae") {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_", vae_variant, "_"), basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[grepl(paste0("_vae_joint_", doc_tag, "_docs_peak_delta_fp_gene_fc_expr_warplda_"), basename(out_dirs))]
  }
  if (!length(out_dirs) && file.exists(file.path(topic_root, "topic_links.csv"))) {
    out_dirs <- topic_root
  }
  if (!length(out_dirs)) return(invisible(NULL))
  for (d in out_dirs) {
    for (method in c("peak_and_gene", "peak_and_gene_prob")) {
      rerun_pathway_from_topic_links(
        out_dir = d,
        topic_links_file = file.path(d, "topic_links.csv"),
        method = method,
        allow_missing = TRUE,
        include_tf = FALSE,
        include_gene = TRUE,
        min_prob = 0,
        per_comparison = isTRUE(per_comparison),
        per_comparison_dir = paste0("per_comparison_pathway_", .pathway_method_suffix(method)),
        per_comparison_flat = TRUE,
        split_direction = isTRUE(split_direction),
        make_heatmap = FALSE,
        top_n_per_topic = top_n_per_topic,
        dot_top_n_per_topic = dot_top_n_per_topic,
        max_pathways = max_pathways,
        enrichr_sleep_time = enrichr_sleep_time,
        enrichr_cache_dir = enrichr_cache_dir,
        enrichr_n_cores = enrichr_n_cores,
        pathway_backend = pathway_backend,
        overwrite = FALSE
      )
    }
  }
  invisible(TRUE)
}

# =============================================================================
# Module 3 public APIs
# =============================================================================

#' Train joint RNA+footprint topic models (Module 3)
#'
#' Builds documents from differential links and trains topic models across a
#' user-supplied K grid. The default backend is the native WarpLDA sampler
#' (`backend = "warplda"` with `warplda_sampler = "warp_omp"`). VAE backends
#' remain available when explicitly requested. The function writes model
#' outputs and model-selection summaries without running downstream topic
#' extraction.
#'
#' @param Kgrid Integer vector of K values for training.
#' @param input_dir Directory containing differential links (filtered up/down).
#' @param output_dir Directory to write topic model outputs.
#' @param input_source Input source type. Use `"differential_links"` for the
#'   historical comparison/differential-link path and `"condition_links"` for
#'   condition-native Module 3 links prepared by
#'   `module3_prepare_condition_links()`.
#' @param sample_subset Optional condition/sample labels to keep. When supplied,
#'   only comparisons whose case and control labels are both in this vector are
#'   used for topic training.
#' @param analysis_label Label used to name this topic-model analysis. If
#'   `NULL`, the label is inferred from `sample_subset` or set to
#'   `"all_samples"`.
#' @param tf_cluster_map Named vector mapping TFs to motif clusters.
#' @param doc_mode Document mode, either `"tf_cluster"` or `"tf"`.
#' @param doc_design Document design, either `"comparison"` or `"condition"`.
#' @param tf_exclude Optional TFs to exclude.
#' @param abs_log2fc_fp_min Minimum |log2FC| footprint change for QC filtering.
#' @param abs_delta_fp_min Minimum |delta_fp| for QC filtering.
#' @param abs_log2fc_gene_min Minimum |log2FC| gene expression change for QC filtering.
#' @param require_fp_bound_either Require fp_bound in either condition.
#' @param require_tf_expr_either Require TF expression in either condition.
#' @param require_gene_expr_either Require gene expression in either condition.
#' @param direction_consistency Direction consistency mode for edge filtering.
#' @param top_terms_per_doc Max terms per document.
#' @param min_df Minimum document frequency for terms.
#' @param count_method Count method ("bin" or "log").
#' @param count_scale Count scaling factor.
#' @param threshold_gene_expr Minimum condition-level target-gene expression.
#' @param threshold_fp_score Minimum condition-level footprint score.
#' @param threshold_tf_expr Minimum condition-level TF expression.
#' @param condition_gene_weighting Optional condition-mode target-gene
#'   weighting. `"none"` retains the constructed document counts;
#'   `"specificity"` emphasizes genes that are relatively more expressed in a
#'   condition while preserving every document's token total.
#' @param condition_peak_weighting Optional condition-mode aggregated-Peak
#'   weighting. `"tf_expression"` multiplies balanced Peak weights by
#'   `log2(TF expression + 1)` relative to the same TF's median across included
#'   conditions before pseudo-count conversion.
#' @param condition_gene_expression_file Long condition-gene expression CSV
#'   with `condition_id`, `gene_key`, and `expression`. When specificity
#'   weighting is selected, the default is `condition_gene_expression.csv`
#'   under `input_dir`.
#' @param condition_specificity_temperature Positive softmax temperature.
#' @param condition_specificity_floor Fraction of uniform expressed-condition
#'   weight retained for every target gene.
#' @param condition_specificity_expression_min Minimum expression used to
#'   define the conditions in which a target gene is expressed. NULL reuses
#'   `threshold_gene_expr`.
#' @param gene_term_mode Gene term mode.
#' @param fp_term_mode Footprint term mode.
#' @param include_tf_terms Whether to include TF self-terms.
#' @param count_input Count column passed to the topic backend.
#' @param vae_variant VAE variant name.
#' @param backend Topic model backend. `"warplda"` is the default native
#'   WarpLDA backend; `"vae"` runs the optional VAE backend.
#' @param vae_python Optional Python executable for VAE training.
#' @param vae_epochs Number of VAE training epochs.
#' @param vae_batch_size VAE mini-batch size.
#' @param vae_hidden VAE hidden-layer width.
#' @param vae_lr VAE learning rate.
#' @param vae_seed VAE random seed.
#' @param vae_device VAE device, for example `"auto"`, `"cpu"`, or `"cuda"`.
#'   `"auto"` uses CUDA when PyTorch can access it and otherwise uses CPU.
#' @param warplda_iterations Number of native WarpLDA iterations.
#' @param warplda_sampler Native WarpLDA sampler. `"warp_omp"` is the default
#'   OpenMP-accelerated doc/word Metropolis-Hastings sampler, `"warp_ref"` is
#'   the slower sequential fixed-seed reference sampler,
#'   `"warp_mh"` is the earlier CraftGRN approximation, and `"gibbs_sync"` is
#'   the faster deterministic collapsed sampler.
#' @param warplda_beta Topic-word prior for native WarpLDA. NULL uses `1/K`,
#'   matching the legacy WarpLDA default.
#' @param warplda_seed Integer random seed for native WarpLDA.
#' @param reuse_if_exists Reuse existing model outputs when all requested K
#'   values are present.
#' @param save_full_doc_term_csv Whether to save the full document-term CSV.
#' @param flat_output Whether to write one setup directly under `output_dir`
#'   instead of adding a long analysis-specific subdirectory. This is intended
#'   for standard package runs; benchmarks keep the nested layout.
#' @param local_threads Optional local thread count for data.table, BLAS, and native WarpLDA. NULL uses all available cores. Set options(craftgrn.warplda.max_threads = n) or CRAFTGRN_WARPLDA_MAX_THREADS to cap automatic thread use.
#' @param memory_safety Module 3 memory policy: `"strict"` fails closed,
#'   `"adaptive"` reduces concurrency where possible, and `"off"` disables
#'   memory preflight checks.
#' @param memory_max_fraction Maximum fraction of currently available memory
#'   that a Module 3 stage may use. Defaults to `0.8`.
#' @param check_repeated_values Warn about repeated inconsistent term values.
#'   The high-throughput default is `FALSE`; set to `TRUE` for diagnostic
#'   audits.
#' @param binarize_method Topic binarization method.
#' @param thrP Topic term probability threshold.
#' @param top_n_terms Number of terms per topic.
#' @param in_topic_min_terms Minimum terms per topic set.
#' @param topic_report_args Additional args for downstream extraction (saved in calc_params).
#' @noRd
train_topic_models <- function(Kgrid,
                               input_dir,
                               output_dir,
                               input_source = c("differential_links", "condition_links"),
                               sample_subset = NULL,
                               analysis_label = NULL,
                               tf_cluster_map,
                               doc_mode = c("tf_cluster", "tf"),
                               doc_design = c("comparison", "condition"),
                               tf_exclude = NULL,
                               abs_log2fc_fp_min = 0,
                               abs_delta_fp_min = 1,
                               abs_log2fc_gene_min = 1,
                               require_fp_bound_either = TRUE,
                               require_tf_expr_either = TRUE,
                               require_gene_expr_either = TRUE,
                               direction_consistency = "aligned",
                               top_terms_per_doc = Inf,
                               min_df = 2,
                               count_method = "bin",
                               count_scale = 50,
                               threshold_gene_expr = 0,
                               threshold_fp_score = 0,
                               threshold_tf_expr = -Inf,
                               condition_gene_weighting = c("none", "specificity"),
                               condition_peak_weighting = c("none", "tf_expression"),
                               condition_gene_expression_file = NULL,
                               condition_specificity_temperature = 0.5,
                               condition_specificity_floor = 0.1,
                               condition_specificity_expression_min = NULL,
                               gene_term_mode = c("aggregate", "unique"),
                               fp_term_mode = c("aggregate", "unique", "aggregate_weight", "gene_expression"),
                               include_tf_terms = FALSE,
                               count_input = NULL,
                               vae_variant = "multivi_encoder",
                               backend = c("warplda", "vae"),
                               vae_python = NULL,
                               vae_epochs = 200L,
                               vae_batch_size = 64L,
                               vae_hidden = 128L,
                               vae_lr = 1e-3,
                               vae_seed = 123L,
                               vae_device = "auto",
                               warplda_iterations = 2000L,
                               warplda_sampler = c("warp_omp", "warp_ref", "warp_mh", "gibbs_sync"),
                               warplda_beta = NULL,
                               warplda_seed = 123L,
                               reuse_if_exists = TRUE,
                               save_full_doc_term_csv = FALSE,
                               flat_output = FALSE,
                               local_threads = NULL,
                               memory_safety = c("strict", "adaptive", "off"),
                               memory_max_fraction = 0.8,
                               check_repeated_values = FALSE,
                               binarize_method = "gammafit",
                               thrP = 0.9,
                               top_n_terms = 500L,
                               in_topic_min_terms = 1,
                               topic_report_args = list()) {
  .assert_pkg("data.table")
  .assert_pkg("cli")
  memory_safety <- match.arg(memory_safety)
  .module3_memory_policy(memory_safety, memory_max_fraction)
  local_threads <- .warplda_default_threads(local_threads)
  data.table::setDTthreads(local_threads)
  Sys.setenv(
    OMP_NUM_THREADS = as.character(local_threads),
    OPENBLAS_NUM_THREADS = as.character(local_threads),
    MKL_NUM_THREADS = as.character(local_threads),
    VECLIB_MAXIMUM_THREADS = as.character(local_threads),
    NUMEXPR_NUM_THREADS = as.character(local_threads)
  )
  .log_inform("train_topic_models using {data.table::getDTthreads()} data.table thread(s).")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  gene_term_mode <- match.arg(gene_term_mode)
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
  backend <- match.arg(backend)
  warplda_iterations <- as.integer(warplda_iterations[[1L]])
  if (!is.finite(warplda_iterations) || warplda_iterations < 1L) .log_abort("`warplda_iterations` must be a positive integer.")
  warplda_sampler <- match.arg(warplda_sampler)
  warplda_seed <- suppressWarnings(as.integer(warplda_seed[[1L]]))
  if (!is.finite(warplda_seed)) .log_abort("`warplda_seed` must be one finite integer.")
  doc_mode <- match.arg(doc_mode)
  doc_design <- match.arg(doc_design)
  input_source <- match.arg(input_source)
  condition_gene_weighting <- match.arg(condition_gene_weighting)
  condition_peak_weighting <- match.arg(condition_peak_weighting)
  if (identical(input_source, "condition_links") && !identical(doc_design, "condition")) {
    .log_abort("input_source = 'condition_links' requires doc_design = 'condition'.")
  }
  if (!identical(condition_gene_weighting, "none") &&
      (!identical(input_source, "condition_links") || !identical(doc_design, "condition"))) {
    .log_abort("Condition gene weighting is available only for condition-link topic models.")
  }
  if (!identical(condition_peak_weighting, "none") &&
      (!identical(doc_design, "condition") ||
       !identical(doc_mode, "tf") ||
       !identical(fp_term_mode, "aggregate"))) {
    .log_abort(
      "TF-expression Peak weighting requires condition::TF documents with fp_term_mode = 'aggregate'."
    )
  }
  condition_specificity_expression_min <- if (is.null(condition_specificity_expression_min)) {
    threshold_gene_expr
  } else {
    suppressWarnings(as.numeric(condition_specificity_expression_min[[1L]]))
  }
  if (!is.finite(condition_specificity_expression_min) ||
      condition_specificity_expression_min < 0) {
    .log_abort("`condition_specificity_expression_min` must be non-negative.")
  }
  if (!identical(condition_gene_weighting, "none")) {
    if (is.null(condition_gene_expression_file) ||
        !length(condition_gene_expression_file) ||
        !nzchar(as.character(condition_gene_expression_file[[1L]]))) {
      condition_gene_expression_file <- file.path(input_dir, "condition_gene_expression.csv")
    } else {
      condition_gene_expression_file <- path.expand(as.character(condition_gene_expression_file[[1L]]))
      if (!.path_is_absolute(condition_gene_expression_file)) {
        condition_gene_expression_file <- file.path(input_dir, condition_gene_expression_file)
      }
    }
    if (!file.exists(condition_gene_expression_file)) {
      .log_abort("Condition gene expression file not found: {condition_gene_expression_file}")
    }
    condition_gene_expression_file <- normalizePath(
      condition_gene_expression_file,
      winslash = "/",
      mustWork = TRUE
    )
  }
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  weight_label <- if (identical(doc_design, "condition")) "peak_score_gene_expr" else "peak_log2fc_fp_gene_fc_expr"
  count_input_requested <- if (is.null(count_input) || !length(count_input)) NA_character_ else as.character(count_input[[1L]])
  count_input_effective <- .resolve_topic_count_input(count_method = count_method, count_input = count_input)
  sample_subset <- if (is.null(sample_subset)) NULL else unique(as.character(sample_subset))
  sample_subset <- sample_subset[!is.na(sample_subset) & nzchar(sample_subset)]
  input_signature <- .module3_topic_input_signature(
    input_dir = input_dir,
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
      condition_gene_weighting = condition_gene_weighting,
      condition_peak_weighting = condition_peak_weighting,
      condition_gene_expression_file = if (!identical(condition_gene_weighting, "none")) {
        info <- file.info(condition_gene_expression_file)
        c(
          path = condition_gene_expression_file,
          size = as.numeric(info$size),
          mtime = as.numeric(info$mtime)
        )
      } else {
        NULL
      },
      condition_specificity_temperature = as.numeric(condition_specificity_temperature),
      condition_specificity_floor = as.numeric(condition_specificity_floor),
      condition_specificity_expression_min = condition_specificity_expression_min,
      include_tf_terms = isTRUE(include_tf_terms),
      top_terms_per_doc = top_terms_per_doc,
      min_df = min_df,
      abs_log2fc_fp_min = abs_log2fc_fp_min,
      abs_delta_fp_min = abs_delta_fp_min,
      abs_log2fc_gene_min = abs_log2fc_gene_min,
      require_fp_bound_either = isTRUE(require_fp_bound_either),
      require_tf_expr_either = isTRUE(require_tf_expr_either),
      require_gene_expr_either = isTRUE(require_gene_expr_either),
      direction_consistency = direction_consistency,
      tf_exclude = sort(unique(toupper(as.character(tf_exclude)))),
      analysis_label = as.character(analysis_label),
      backend = backend,
      warplda_sampler = warplda_sampler,
      warplda_iterations = warplda_iterations,
      warplda_beta = warplda_beta,
      warplda_seed = warplda_seed,
      vae_variant = vae_variant,
      vae_epochs = vae_epochs,
      vae_batch_size = vae_batch_size,
      vae_hidden = vae_hidden,
      vae_lr = vae_lr,
      vae_seed = vae_seed
    )
  )
  if (identical(input_source, "condition_links")) {
    .log_inform("Loading condition-native links from {input_dir}.")
    edges_dt <- .module3_read_condition_links(input_dir, conditions = sample_subset)
    n_loaded <- nrow(edges_dt)
    .log_inform("Loaded {n_loaded} condition-link row(s).")
  } else {
    delta_files <- .module3_filtered_link_files(input_dir)
    if (!length(delta_files)) {
      delta_files <- list.files(input_dir, "_filtered_links(_(up|down))?\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) {
      delta_files <- list.files(input_dir, "_delta_links_filtered(_(up|down))?\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) {
      delta_files <- list.files(input_dir, "_delta_links\\.csv$", full.names = TRUE)
    }
    if (!length(delta_files)) .log_abort("No delta link files found in {input_dir}")

    .log_inform("Loading {length(delta_files)} delta-link file(s) from {input_dir}.")
    edges_all <- load_delta_links_many(delta_files, keep_original = FALSE)
    edges_dt <- data.table::as.data.table(edges_all)
    edges_dt <- .apply_module3_manifest_labels(edges_dt, input_dir)
    n_loaded <- nrow(edges_dt)
    .log_inform("Loaded {n_loaded} delta-link row(s).")
    if (!("comparison_id" %in% names(edges_dt))) .log_abort("edges_all missing comparison_id.")
    if (length(sample_subset)) {
      if (all(c("cond1_id", "cond2_id") %in% names(edges_dt))) {
        n_before_subset <- nrow(edges_dt)
        edges_dt <- edges_dt[cond1_id %in% sample_subset & cond2_id %in% sample_subset]
        .log_inform("Sample subset retained {nrow(edges_dt)}/{n_before_subset} delta-link row(s).")
      } else {
        .log_abort("sample_subset requires delta links with cond1_id and cond2_id columns.")
      }
    }
  }
  if (!nrow(edges_dt)) {
    .log_abort("No topic-input rows remain after applying input filters.")
  }
  if (is.null(analysis_label) || !nzchar(as.character(analysis_label[[1L]]))) {
    analysis_label <- if (length(sample_subset)) "sample_subset" else "all_samples"
  }
  analysis_label <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(analysis_label[[1L]]))
  if (!nzchar(analysis_label)) analysis_label <- "analysis"

  vae_script <- Sys.getenv("VAE_SCRIPT", unset = "")
  if (!nzchar(vae_script)) {
    vae_script <- system.file("python", "logistic_normal_vae_topics.py", package = "craftgrn")
  }
  if (!nzchar(vae_script) || !file.exists(vae_script)) {
    cand <- file.path("dev", "logistic_normal_vae_topics.py")
    if (file.exists(cand)) vae_script <- cand
  }
  if (!file.exists(vae_script)) .log_abort("Missing VAE script: {vae_script}")

  for (analysis_id in analysis_label) {
    t_cell <- proc.time()[["elapsed"]]
    .log_inform("Preparing topic input for {analysis_id}: doc_design={doc_design}, doc_mode={doc_mode}, fp_term_mode={fp_term_mode}, backend={backend}.")
    edges_sub <- data.table::copy(edges_dt)
    if (!is.null(tf_exclude) && length(tf_exclude)) {
      edges_sub <- edges_sub[!toupper(tf) %in% tf_exclude]
    }
    if (!nrow(edges_sub)) next

    .log_inform("{analysis_id}: filtering {nrow(edges_sub)} edge row(s).")
    t_filter <- proc.time()[["elapsed"]]
    edges_filt <- if (identical(input_source, "condition_links")) {
      edges_sub
    } else {
      filter_edges_for_tf_topics(
        edges_sub,
        abs_log2fc_fp_min = abs_log2fc_fp_min,
        abs_delta_fp_min = abs_delta_fp_min,
        abs_log2fc_gene_min = abs_log2fc_gene_min,
        require_fp_bound_either = require_fp_bound_either,
        require_tf_expr_either = require_tf_expr_either,
        require_gene_expr_either = require_gene_expr_either,
        direction_consistency = direction_consistency
      )
    }
    if (!nrow(edges_filt)) next
    .log_inform("{analysis_id}: retained {nrow(edges_filt)} edge row(s) after topic filters in {round(proc.time()[['elapsed']] - t_filter, 1)} sec.")

    if (identical(doc_design, "condition")) {
      .log_inform("{analysis_id}: building condition-level TF documents.")
      t_docs <- proc.time()[["elapsed"]]
      edges_docs <- if (identical(input_source, "condition_links")) {
        data.table::copy(edges_filt)
      } else {
        add_condition_tf_docs(
          edges_filt,
          tf_cluster_map = tf_cluster_map,
          doc_mode = doc_mode
        )
      }
      .log_inform("{analysis_id}: built {nrow(edges_docs)} condition-level document edge row(s) in {round(proc.time()[['elapsed']] - t_docs, 1)} sec.")
      .log_inform("{analysis_id}: building condition-level doc_term.")
      t_doc_term <- proc.time()[["elapsed"]]
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
        condition_peak_weighting = condition_peak_weighting,
        check_repeated_values = check_repeated_values
      )
      .log_inform("{analysis_id}: condition-level doc_term build finished in {round(proc.time()[['elapsed']] - t_doc_term, 1)} sec.")
    } else {
      .log_inform("{analysis_id}: building comparison-level TF documents.")
      t_docs <- proc.time()[["elapsed"]]
      edges_docs <- add_tf_docs(
        edges_filt,
        doc_mode = doc_mode,
        direction_by = "gene",
        tf_cluster_map = tf_cluster_map
      )
      .log_inform("{analysis_id}: built {nrow(edges_docs)} comparison-level document edge row(s) in {round(proc.time()[['elapsed']] - t_docs, 1)} sec.")
      .log_inform("{analysis_id}: building comparison-level doc_term.")
      t_doc_term <- proc.time()[["elapsed"]]
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
        include_tf_terms = include_tf_terms,
        tf_weight_type = "log2fc_tf",
        balance_mode = "min",
        prefix_terms = TRUE,
        threshold_gene_expr = threshold_gene_expr,
        threshold_fp_score = threshold_fp_score,
        threshold_tf_expr = threshold_tf_expr,
        require_condition_thresholds = identical(doc_mode, "tf"),
        check_repeated_values = check_repeated_values
      )
      .log_inform("{analysis_id}: comparison-level doc_term build finished in {round(proc.time()[['elapsed']] - t_doc_term, 1)} sec.")
    }
    if (!nrow(doc_term)) {
      .log_inform("Skipping topic training: no doc_term for {analysis_id}")
      next
    }
    peak_weighting_audit <- attr(doc_term, "condition_peak_weighting_audit")
    specificity_audit <- NULL
    if (identical(condition_gene_weighting, "specificity")) {
      .log_inform("{analysis_id}: applying condition target-gene specificity weighting.")
      weighted <- .module3_apply_condition_gene_specificity(
        doc_term = doc_term,
        expression_file = condition_gene_expression_file,
        count_column = count_input_effective,
        expression_min = condition_specificity_expression_min,
        temperature = condition_specificity_temperature,
        uniform_floor = condition_specificity_floor
      )
      doc_term <- weighted$doc_term
      specificity_audit <- weighted$audit
      rm(weighted)
      invisible(gc())
    }
    .log_inform("{analysis_id}: doc_term has {nrow(doc_term)} row(s), {data.table::uniqueN(doc_term$doc_id)} document(s), and {data.table::uniqueN(doc_term$term_id)} term(s).")

    token_cap <- list(
      raw_tokens = as.double(sum(.safe_num(doc_term[[count_input_effective]]), na.rm = TRUE)),
      tokens = as.double(sum(.safe_num(doc_term[[count_input_effective]]), na.rm = TRUE)),
      scale_factor = 1
    )
    if (identical(backend, "warplda")) {
      token_cap <- .cap_warplda_token_counts(doc_term[[count_input_effective]])
      if (token_cap$scale_factor < 1) {
        doc_term[, (count_input_effective) := token_cap$counts]
        .log_warn(
          "{analysis_id}: reduced WarpLDA expanded tokens from {format(token_cap$raw_tokens, scientific = FALSE)} to {format(token_cap$tokens, scientific = FALSE)} (scale factor {signif(token_cap$scale_factor, 4)}) to stay below the native integer index limit."
        )
      }
    }
    if (!is.null(peak_weighting_audit)) {
      peak_weighting_audit[, `:=`(
        model_tokens_before_global_cap = as.double(token_cap$raw_tokens),
        global_token_scale_factor = as.double(token_cap$scale_factor),
        model_tokens_after_global_cap = as.double(token_cap$tokens)
      )]
    }

    local_topic_args <- modifyList(list(
      binarize_method = binarize_method,
      thrP = thrP,
      top_n_terms = top_n_terms,
      in_topic_min_terms = in_topic_min_terms
    ), topic_report_args)

    model_name <- if (backend == "vae") vae_variant else "warplda"
    out_dir <- if (isTRUE(flat_output)) {
      output_dir
    } else {
      file.path(
        output_dir,
        paste0(analysis_id, "_vae_joint_", doc_tag, "_docs_", weight_label, "_", model_name, "_Kgrid")
      )
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!is.null(specificity_audit)) {
      data.table::fwrite(
        specificity_audit,
        file.path(out_dir, "condition_gene_specificity_weighting_audit.csv")
      )
    }
    if (!is.null(peak_weighting_audit)) {
      data.table::fwrite(
        peak_weighting_audit,
        file.path(out_dir, "condition_peak_tf_expression_weighting_audit.csv")
      )
    }
    if (identical(doc_design, "condition")) {
      .write_module3_document_term_qc(
        doc_term = doc_term,
        output_dir = out_dir,
        count_column = count_input_effective,
        title = paste0("Module 3 document-term QC - ", analysis_id),
        verbose = TRUE
      )
    }
    topic_input_summary <- data.table::data.table(
      analysis_label = analysis_id,
      input_signature = input_signature,
      input_source = input_source,
      doc_design = doc_design,
      doc_mode = doc_mode,
      fp_term_mode = fp_term_mode,
      condition_gene_weighting = condition_gene_weighting,
      condition_peak_weighting = condition_peak_weighting,
      condition_peak_scaling = if (identical(condition_peak_weighting, "tf_expression")) {
        "per_tf_log2p1_median"
      } else {
        "none"
      },
      count_method = count_method,
      count_scale = as.numeric(count_scale),
      count_input_requested = count_input_requested,
      count_input_effective = count_input_effective,
      n_model_tokens_raw = as.double(token_cap$raw_tokens),
      token_scale_factor = as.double(token_cap$scale_factor),
      n_link_rows_after_filter = as.double(nrow(edges_filt)),
      n_document_edge_rows = as.double(nrow(edges_docs)),
      n_doc_term_rows = as.double(nrow(doc_term)),
      n_documents = as.double(data.table::uniqueN(doc_term$doc_id)),
      n_terms = as.double(data.table::uniqueN(doc_term$term_id)),
      n_nonzero = as.double(sum(.safe_num(doc_term[[count_input_effective]]) > 0, na.rm = TRUE)),
      n_model_tokens = as.double(token_cap$tokens)
    )
    data.table::fwrite(topic_input_summary, file.path(out_dir, "topic_input_summary.csv"))
    metrics_path <- file.path(out_dir, "model_metrics.csv")
    models_dir <- file.path(out_dir, "vae_models")
    required_k <- sort(unique(as.integer(Kgrid)))
    required_k <- required_k[is.finite(required_k)]
    reset_done <- .reset_topic_model_artifacts(out_dir, backend, reuse_if_exists)
    if (isTRUE(reset_done)) {
      .log_inform("{analysis_id}: cleared existing topic model artifacts because reuse_if_exists = FALSE.")
    }
    metrics_have_required_k <- FALSE
    model_metadata_ok <- !file.exists(metrics_path)
    if (file.exists(metrics_path)) {
      old_metrics <- tryCatch(data.table::fread(metrics_path), error = function(e) data.table::data.table())
      if (nrow(old_metrics) && "K" %in% names(old_metrics)) {
        model_metadata_ok <- all(c("count_method", "count_input_effective", "input_signature") %in% names(old_metrics)) &&
          all(as.character(old_metrics$count_method) == count_method) &&
          all(as.character(old_metrics$count_input_effective) == count_input_effective) &&
          all(as.character(old_metrics$input_signature) == input_signature)
        if (identical(backend, "warplda")) {
          model_metadata_ok <- model_metadata_ok &&
            "sampler" %in% names(old_metrics) &&
            all(as.character(old_metrics$sampler) == warplda_sampler)
        }
        metrics_have_required_k <- all(required_k %in% as.integer(old_metrics$K))
        old_required <- old_metrics[as.integer(old_metrics$K) %in% required_k]
        count_metadata_ok <- all(c("count_method", "count_input_effective", "input_signature") %in% names(old_required)) &&
          all(as.character(old_required$count_method) == count_method) &&
          all(as.character(old_required$count_input_effective) == count_input_effective) &&
          all(as.character(old_required$input_signature) == input_signature)
        metrics_have_required_k <- metrics_have_required_k && isTRUE(count_metadata_ok)
        if (identical(backend, "warplda")) {
          if ("sampler" %in% names(old_metrics)) {
            metrics_have_required_k <- metrics_have_required_k &&
              all(as.character(old_metrics$sampler[as.integer(old_metrics$K) %in% required_k]) == warplda_sampler)
          } else {
            metrics_have_required_k <- metrics_have_required_k && identical(warplda_sampler, "gibbs_sync")
          }
        }
      }
    }
    if (isTRUE(reuse_if_exists) && !isTRUE(model_metadata_ok)) {
      reset_invalid <- .reset_topic_model_artifacts(out_dir, backend, reuse_if_exists = FALSE)
      if (isTRUE(reset_invalid)) {
        .log_warn("{analysis_id}: cleared topic model artifacts because the input or model settings changed.")
      }
      metrics_have_required_k <- FALSE
    }
    model_files_have_required_k <- if (dir.exists(models_dir)) {
      all(file.exists(file.path(models_dir, sprintf("theta_K%d.csv", required_k)))) &&
        all(file.exists(file.path(models_dir, sprintf("phi_K%d.csv", required_k))))
    } else {
      FALSE
    }
    reuse_ok <- isTRUE(reuse_if_exists) &&
      file.exists(metrics_path) &&
      metrics_have_required_k &&
      model_files_have_required_k

    if (reuse_ok) {
      .log_inform("Reusing existing topic model outputs for K={paste(required_k, collapse = ',')} in {out_dir}")
      next
    }

    if (backend == "vae") {
      .log_inform("{analysis_id}: starting VAE training in {out_dir}.")
      t_model <- proc.time()[["elapsed"]]
      run_vae_topic_report_py(
        doc_term = doc_term,
        edges_docs = edges_docs,
        out_dir = out_dir,
        option_label = "joint",
        direction_by = "gene",
        vae_script = vae_script,
        k_grid = Kgrid,
        vae_variant = vae_variant,
        vae_python = vae_python,
        vae_epochs = vae_epochs,
        vae_batch_size = vae_batch_size,
        vae_hidden = vae_hidden,
        vae_lr = vae_lr,
        vae_seed = vae_seed,
        vae_device = vae_device,
        do_report = FALSE,
        reuse_if_exists = reuse_if_exists,
        count_input = count_input_effective,
        save_full_doc_term_csv = save_full_doc_term_csv,
        topic_report_args = local_topic_args
      )
      if (file.exists(metrics_path)) {
        metrics_tbl <- data.table::fread(metrics_path, showProgress = FALSE)
        metrics_tbl[, `:=`(
          input_signature = input_signature,
          count_method = count_method,
          count_scale = as.numeric(count_scale),
          count_input_requested = count_input_requested,
          count_input_effective = count_input_effective
        )]
        data.table::fwrite(metrics_tbl, metrics_path)
        .save_all(out_dir, "model_metrics", metrics_tbl)
      }
      .log_inform("{analysis_id}: finished VAE training in {round(proc.time()[['elapsed']] - t_model, 1)} sec: {out_dir}.")
    } else {
      .log_inform("{analysis_id}: writing WarpLDA input caches in {out_dir}.")
      t_cache <- proc.time()[["elapsed"]]
      if (exists("write_doc_term_cache", mode = "function")) {
        write_doc_term_cache(
          doc_term,
          out_dir = out_dir,
          save_full_doc_term_csv = isTRUE(save_full_doc_term_csv)
        )
      } else {
        data.table::fwrite(utils::head(doc_term, 100L), file.path(out_dir, "doc_term_first100.csv"))
        if (isTRUE(save_full_doc_term_csv)) {
          data.table::fwrite(doc_term, file.path(out_dir, "doc_term.csv"))
        }
      }
      .save_all(out_dir, "doc_term", doc_term)
      .save_all(out_dir, "edges_docs", edges_docs)
      .log_inform("{analysis_id}: WarpLDA input cache writing finished in {round(proc.time()[['elapsed']] - t_cache, 1)} sec.")
      .log_inform("{analysis_id}: building sparse DTM for WarpLDA.")
      t_dtm <- proc.time()[["elapsed"]]
      dtm_obj <- build_sparse_dtm(doc_term, count_col = count_input_effective)
      dtm <- dtm_obj$dtm
      .save_all(out_dir, "dtm", dtm)
      .save_all(out_dir, "dtm_index", list(doc_index = dtm_obj$doc_index, term_index = dtm_obj$term_index))
      .log_inform("{analysis_id}: WarpLDA sparse DTM cache finished in {round(proc.time()[['elapsed']] - t_dtm, 1)} sec.")

      .log_inform("{analysis_id}: starting WarpLDA model fits.")
      t_model <- proc.time()[["elapsed"]]
      fits_out <- run_warplda_models(
        dtm,
        K_grid = Kgrid,
        iterations = warplda_iterations,
        alpha_by_topic = TRUE,
        alpha = NULL,
        beta = warplda_beta,
        seed = warplda_seed,
        save_tmp_dir = file.path(out_dir, "tmp_models"),
        workers = 1L,
        threads_per_model = local_threads,
        sampler = warplda_sampler,
        metrics_file = file.path(out_dir, "model_metrics.csv"),
        memory_safety = memory_safety,
        memory_max_fraction = memory_max_fraction
      )
      .log_inform("{analysis_id}: finished WarpLDA model fits in {round(proc.time()[['elapsed']] - t_model, 1)} sec.")
      metrics_tbl <- data.table::as.data.table(fits_out$metrics)
      metrics_tbl[, `:=`(
        input_signature = input_signature,
        count_method = count_method,
        count_scale = as.numeric(count_scale),
        count_input_requested = count_input_requested,
        count_input_effective = count_input_effective
      )]
      data.table::fwrite(metrics_tbl, file.path(out_dir, "model_metrics.csv"))
      .save_all(out_dir, "model_metrics", metrics_tbl)

      dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
      for (fit_file in fits_out$fit_files) {
        if (is.na(fit_file) || !file.exists(fit_file)) next
        fit <- readRDS(fit_file)
        K <- as.integer(fit$K)
        theta <- fit$theta
        phi <- fit$phi
        if (!is.null(theta)) {
          theta_df <- data.frame(doc_id = rownames(theta), theta, check.names = FALSE)
          readr::write_csv(theta_df, file.path(models_dir, sprintf("theta_K%d.csv", K)))
        }
        if (!is.null(phi)) {
          phi_df <- data.frame(term_id = rownames(phi), phi, check.names = FALSE)
          readr::write_csv(phi_df, file.path(models_dir, sprintf("phi_K%d.csv", K)))
        }
        rm(fit, theta, phi)
        invisible(gc())
      }

      title_prefix <- .topic_model_selection_title(out_dir, backend_label = "WarpLDA")
      sel <- plot_model_selection_cistopic(metrics_tbl, file.path(out_dir, "model_selection.pdf"), title_prefix = title_prefix)
      .save_all(out_dir, "model_selection", sel)
    }
    .log_inform("{analysis_id}: finished train_topic_models analysis stage in {round(proc.time()[['elapsed']] - t_cell, 1)} sec.")
  }

  invisible(TRUE)
}

#' Train Module 3 topic models
#'
#' Public step function for training one Module 3 topic-model setup after
#' [module3_prepare_differential_links()] has produced filtered differential
#' links. This is a thin Module 3-named wrapper around the internal training
#' engine.
#'
#' @param k_grid Integer vector of K values for training.
#' @param filtered_dir Directory containing Module 3 filtered differential-link
#'   files.
#' @param output_dir Directory to write topic model outputs.
#' @param flat_output Whether to write this selected setup directly under
#'   `output_dir`. Defaults to `TRUE` for the public step API.
#' @param ... Additional arguments passed to the internal training engine, such
#'   as `doc_design`, `fp_term_mode`, `backend`, and `local_threads`.
#'
#' @return Invisibly returns TRUE when training completes.
#' @export
module3_train_topic_models <- function(k_grid,
                                       filtered_dir,
                                       output_dir,
                                       flat_output = TRUE,
                                       ...) {
  train_topic_models(
    Kgrid = k_grid,
    input_dir = filtered_dir,
    output_dir = output_dir,
    flat_output = flat_output,
    ...
  )
}

#' Extract regulatory topics from trained models (Module 3)
#'
#' Uses precomputed topic models to compute link-topic scores, topic
#' assignments, and downstream reports for a user-selected K. The default
#' backend is native WarpLDA to match `train_topic_models()`.
#'
#' @param k Integer K selected by the user.
#' @param model_dir Directory containing trained topic model outputs.
#' @param output_dir Directory to write final topic reports.
#' @param topic_input_dir Optional directory containing shared topic input
#'   caches. When supplied, `dtm.rds` and `edges_docs.rds` can be loaded from
#'   this directory if they are not present in the model directory.
#' @param backend Topic model backend. `"warplda"` is the default native
#'   WarpLDA backend; `"vae"` reads optional VAE outputs.
#' @param vae_variant VAE variant name used in trained output directories.
#' @param doc_mode Document mode used during training.
#' @param weight_label Weight label used in trained output directories.
#' @param flatten_single_output If `TRUE` and only one trained output directory
#'   is matched, write reports directly under `output_dir`.
#' @param topic_report_args Optional list of overrides for report settings.
#'   Standard extraction keeps pathway enrichment disabled by default for speed.
#'   Set `run_pathway_enrichment = TRUE` for the overall topic pathway review.
#'   Comparison-topic runs can also set `pathway_per_comparison = TRUE`.
#' @noRd
extract_regulatory_topics <- function(k,
                                      model_dir,
                                      output_dir,
                                      topic_input_dir = NULL,
                                      backend = c("warplda", "vae"),
                                      vae_variant = "multivi_encoder",
                                      doc_mode = c("tf_cluster", "tf"),
                                      weight_label = c("peak_log2fc_fp_gene_fc_expr", "peak_delta_fp_gene_fc_expr", "peak_score_gene_expr", "gene_expression"),
                                      flatten_single_output = FALSE,
                                      topic_score_method = NULL,
                                      topic_term_assignment_method = NULL,
                                      topic_report_args = list()) {
  .assert_pkg("data.table")
  .assert_pkg("readr")

  k <- as.integer(k)
  if (!is.finite(k) || k <= 0L) .log_abort("`k` must be a positive integer.")
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  weight_label <- match.arg(weight_label)
  if (!is.null(topic_score_method)) {
    topic_report_args$topic_score_method <- match.arg(
      topic_score_method,
      choices = c("normtop_specificity", "rowmax_phi")
    )
  }
  if (!is.null(topic_term_assignment_method)) {
    topic_report_args$topic_term_assignment_method <- match.arg(
      topic_term_assignment_method,
      choices = c("gammafit_maxprob", "max_phi", "gammafit")
    )
  }
  topic_term_assignment_explicit <- !is.null(topic_report_args$topic_term_assignment_method)
  report_doc_design <- if (weight_label %in% c("peak_score_gene_expr", "gene_expression")) "condition" else "comparison"
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_dirs <- list.dirs(model_dir, recursive = FALSE, full.names = TRUE)
  if (dir.exists(file.path(model_dir, "vae_models"))) {
    out_dirs <- c(model_dir, out_dirs)
  }
  if (backend == "vae") {
    patt <- paste0("_vae_joint_", doc_tag, "_docs_", weight_label, "_", vae_variant, "_")
    out_dirs <- out_dirs[out_dirs == model_dir | grepl(patt, basename(out_dirs))]
  } else {
    out_dirs <- out_dirs[out_dirs == model_dir | grepl(paste0("_vae_joint_", doc_tag, "_docs_", weight_label, "_warplda_"), basename(out_dirs))]
  }
  out_dirs <- unique(out_dirs)
  if (!length(out_dirs)) .log_abort("No trained topic model directories found in {model_dir}")
  flatten_single_output <- isTRUE(flatten_single_output) && length(out_dirs) == 1L

  for (d in out_dirs) {
    rds_dir <- file.path(d, "rds")
    dtm_path <- file.path(rds_dir, "dtm.rds")
    edges_docs_path <- file.path(rds_dir, "edges_docs.rds")
    if ((!file.exists(dtm_path) || !file.exists(edges_docs_path)) && !is.null(topic_input_dir)) {
      shared_rds_dir <- file.path(topic_input_dir, "rds")
      shared_dtm_path <- file.path(shared_rds_dir, "dtm.rds")
      shared_edges_docs_path <- file.path(shared_rds_dir, "edges_docs.rds")
      if (!file.exists(dtm_path) && file.exists(shared_dtm_path)) dtm_path <- shared_dtm_path
      if (!file.exists(edges_docs_path) && file.exists(shared_edges_docs_path)) edges_docs_path <- shared_edges_docs_path
    }
    metrics_path <- file.path(d, "model_metrics.csv")
    topic_input_summary_path <- file.path(d, "topic_input_summary.csv")
    theta_path <- file.path(d, "vae_models", sprintf("theta_K%d.csv", k))
    phi_path <- file.path(d, "vae_models", sprintf("phi_K%d.csv", k))
    if (!file.exists(dtm_path) || !file.exists(edges_docs_path)) {
      .log_abort("Missing dtm or edges_docs for {d}")
    }
    if (!file.exists(theta_path) || !file.exists(phi_path)) {
      .log_abort("Missing theta/phi for K={k} in {d}")
    }

    dtm <- readRDS(dtm_path)
    edges_docs <- readRDS(edges_docs_path)
    metrics_tbl <- if (file.exists(metrics_path)) readr::read_csv(metrics_path, show_col_types = FALSE) else NULL
    theta_df <- readr::read_csv(theta_path, show_col_types = FALSE)
    phi_df <- readr::read_csv(phi_path, show_col_types = FALSE)
    theta <- as.matrix(theta_df[, -1, drop = FALSE]); rownames(theta) <- theta_df[[1]]
    phi <- as.matrix(phi_df[, -1, drop = FALSE]); rownames(phi) <- phi_df[[1]]
    theta <- .validate_topic_probability_matrix(theta, "theta", theta_path)
    phi <- .validate_topic_probability_matrix(phi, "phi", phi_path)

    topic_base <- list(K = k, theta = theta, phi = phi,
                       metrics = if (is.data.frame(metrics_tbl)) metrics_tbl[metrics_tbl$K == k, ] else NULL)

    if (isTRUE(flatten_single_output)) {
      out_dir <- output_dir
    } else {
      base_name <- basename(d)
      base_name <- sub("_Kgrid$", paste0("_K", k), base_name)
      out_dir <- file.path(output_dir, base_name)
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  wrap_comp_label <- function(x) {
    vapply(x, function(s) {
      gsub("::", " :: ", s, fixed = TRUE)
    }, character(1))
  }

  defaults <- list(
      binarize_method = "gammafit",
      topic_score_method = "normtop_specificity",
      topic_term_assignment_method = "gammafit_maxprob",
      topic_model_family = if (identical(backend, "warplda")) "lda" else vae_variant,
      thrP = .module3_default_gammafit_thrP(
        if (identical(backend, "warplda")) "lda" else vae_variant
      ),
      top_n_terms = 500L,
      in_topic_min_terms = 1L,
      pathway_use_all_terms = FALSE,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
      top_n_per_topic = 100L,
      max_pathways = 1000L,
      pathway_tf_link_mode = "theta",
      pathway_tf_top_n_docs = 50L,
      pathway_tf_min_theta = NA_real_,
      run_pathway_gsea = FALSE,
      gsea_species = "Homo sapiens",
      gsea_nperm = 1000L,
      gsea_peak_gene_agg = "max",
      pathway_source = "topic_terms",
      pathway_link_scores_file = NULL,
      pathway_link_scores_file_tf = NULL,
      pathway_link_gene_terms_file = NULL,
      pathway_link_min_prob = 0,
      pathway_link_include_tf = TRUE,
      pathway_link_include_gene = TRUE,
      pathway_link_gene_min_prob = 0,
      pathway_link_tf_min_prob = 0.5,
      pathway_link_tf_max_topics = 5L,
      pathway_link_tf_top_n_per_topic = 30L,
      pathway_per_comparison = FALSE,
      pathway_per_comparison_dir = ".",
      pathway_per_comparison_flat = TRUE,
      pathway_split_direction = TRUE,
      run_link_topic_scores = FALSE,
      link_topic_gate_mode = "none",
      link_topic_top_k = 3L,
      link_topic_min_prob = 0,
      link_topic_include_tf = FALSE,
      link_topic_chunk_size = 5000L,
      link_topic_n_cores = 1L,
      link_topic_overwrite = FALSE,
      link_topic_method = "gene_prob",
      link_topic_prob_cutoff = 0.3,
      link_topic_fdr_q = 0.2,
      link_topic_fdr_p = NA_real_,
      link_topic_efdr_scope = "per_topic",
      link_topic_efdr_B = 100L,
      link_topic_efdr_seed = 1L,
      link_topic_output = "pass",
      run_gammafit_summary = TRUE,
      run_link_efdr_summary = TRUE,
      run_pathway_enrichment = FALSE,
      run_raw_theta_document_heatmap = FALSE,
      run_document_theta_umap = TRUE,
      run_topic_by_comparison_heatmaps = TRUE,
      run_intertopic_distance_map = TRUE,
      fp_term_mode = "aggregate",
      topic_by_comparison_label_cleaner = wrap_comp_label
    )
    if (is.null(topic_report_args$fp_term_mode) && file.exists(topic_input_summary_path)) {
      topic_input_summary <- tryCatch(
        data.table::fread(topic_input_summary_path, nrows = 1L),
        error = function(e) data.table::data.table()
      )
      if (nrow(topic_input_summary) && "fp_term_mode" %in% names(topic_input_summary)) {
        inferred_fp_term_mode <- as.character(topic_input_summary$fp_term_mode[[1L]])
        if (nzchar(inferred_fp_term_mode)) {
          topic_report_args$fp_term_mode <- inferred_fp_term_mode
          .log_inform(
            "Using fp_term_mode={inferred_fp_term_mode} from topic input summary for topic extraction."
          )
        }
      }
    }
    if (!isTRUE(topic_term_assignment_explicit)) {
      resolved_fp_term_mode <- .resolve_fp_term_mode(
        topic_report_args$fp_term_mode %||% defaults$fp_term_mode
      )
      defaults$topic_term_assignment_method <- if (resolved_fp_term_mode %in% c("aggregate", "gene_expression")) {
        "gammafit_maxprob"
      } else {
        "max_phi"
      }
    }
    args <- modifyList(defaults, topic_report_args)
    combo_tag <- basename(dirname(d))
    model_label <- if (backend == "vae") paste("VAE", vae_variant) else "WarpLDA"
    do.call(
      run_tfdocs_report_from_topic_base,
      c(list(
        topic_base = topic_base,
        dtm = dtm,
        edges_docs = edges_docs,
        out_dir = out_dir,
        option_label = "joint",
        direction_by = "gene",
        doc_design = report_doc_design,
        title_prefix = if (backend == "vae") {
          paste0(model_label, "\n", combo_tag)
        } else {
          paste0(model_label, "\n", combo_tag)
        }
      ), args)
    )
    writeLines(report_doc_design, file.path(out_dir, "report_doc_design.txt"))
  }

  invisible(TRUE)
}

.module3_extraction_memory_bytes <- function(model_dir, memory_gb = NULL) {
  gib <- 1024^3
  memory_gb <- if (is.null(memory_gb) || !length(memory_gb)) {
    NA_real_
  } else {
    suppressWarnings(as.numeric(memory_gb)[[1L]])
  }
  model_dir <- normalizePath(as.character(model_dir)[[1L]], winslash = "/", mustWork = FALSE)
  candidates <- unique(c(
    file.path(model_dir, "rds", "edges_docs.rds"),
    list.files(model_dir, pattern = "^edges_docs[.]rds$", recursive = TRUE, full.names = TRUE)
  ))
  candidates <- candidates[file.exists(candidates)]
  compressed_bytes <- if (length(candidates)) {
    max(as.numeric(file.info(candidates)$size), na.rm = TRUE)
  } else {
    0
  }
  auto_bytes <- max(16 * gib, compressed_bytes * 24 + 4 * gib)
  if (is.finite(memory_gb) && memory_gb >= 1) {
    return(max(auto_bytes, memory_gb * gib))
  }
  auto_bytes
}

.module3_extraction_k_worker_plan <- function(n_tasks,
                                              model_dir,
                                              k_workers = NULL,
                                              k_max_workers = 4L,
                                              k_memory_gb = NULL,
                                              k_memory_reserve_gb = 32,
                                              cores = NULL,
                                              link_scoring = FALSE,
                                              available_bytes = .available_memory_bytes(),
                                              os_type = .Platform$OS.type) {
  gib <- 1024^3
  n_tasks <- suppressWarnings(as.integer(n_tasks)[[1L]])
  if (!is.finite(n_tasks) || n_tasks < 1L) n_tasks <- 1L
  k_max_workers <- suppressWarnings(as.integer(k_max_workers)[[1L]])
  if (!is.finite(k_max_workers) || k_max_workers < 1L) k_max_workers <- 4L
  cores <- if (is.null(cores)) .available_cores(logical = TRUE) else suppressWarnings(as.integer(cores)[[1L]])
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  requested <- if (is.null(k_workers) || !length(k_workers)) {
    n_tasks
  } else {
    suppressWarnings(as.integer(k_workers)[[1L]])
  }
  if (!is.finite(requested) || requested < 1L) requested <- 1L
  per_worker_bytes <- .module3_extraction_memory_bytes(model_dir, memory_gb = k_memory_gb)
  reserve_gb <- suppressWarnings(as.numeric(k_memory_reserve_gb)[[1L]])
  if (!is.finite(reserve_gb) || reserve_gb < 0) reserve_gb <- 32
  reserve_bytes <- if (is.finite(available_bytes) && available_bytes > 0) {
    max(reserve_gb * gib, available_bytes * 0.25)
  } else {
    NA_real_
  }
  memory_workers <- if (is.finite(available_bytes) && available_bytes > reserve_bytes) {
    max(1L, floor((available_bytes - reserve_bytes) / per_worker_bytes))
  } else {
    1L
  }
  reason <- "adaptive memory and CPU limits"
  workers <- min(n_tasks, requested, k_max_workers, cores, memory_workers)
  if (!identical(os_type, "unix")) {
    workers <- 1L
    reason <- "parallel K extraction currently requires fork support"
  } else if (!is.finite(available_bytes) || available_bytes <= 0) {
    workers <- 1L
    reason <- "available memory could not be measured"
  } else if (isTRUE(link_scoring)) {
    workers <- 1L
    reason <- "topic-link scoring or topic optimization requires sequential K extraction"
  } else if (memory_workers <= 1L) {
    workers <- 1L
    reason <- "available memory headroom supports one K worker"
  }
  list(
    workers = as.integer(max(1L, workers)),
    n_tasks = n_tasks,
    requested_workers = as.integer(requested),
    max_workers = as.integer(k_max_workers),
    cores = as.integer(cores),
    memory_workers = as.integer(memory_workers),
    per_worker_bytes = as.numeric(per_worker_bytes),
    available_bytes = as.numeric(available_bytes),
    reserve_bytes = as.numeric(reserve_bytes),
    reason = reason
  )
}

.module3_run_extraction_k_tasks <- function(tasks,
                                            worker_fun,
                                            workers,
                                            cores = NULL) {
  workers <- max(1L, suppressWarnings(as.integer(workers)[[1L]]))
  cores <- if (is.null(cores)) workers else max(1L, suppressWarnings(as.integer(cores)[[1L]]))
  threads_per_worker <- max(1L, floor(cores / workers))
  run_one <- function(i) {
    warnings <- character()
    old_threads <- if (requireNamespace("data.table", quietly = TRUE)) data.table::getDTthreads() else 1L
    on.exit({
      if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(old_threads)
    }, add = TRUE)
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(threads_per_worker)
    value <- tryCatch(
      withCallingHandlers(
        suppressMessages(worker_fun(tasks[[i]])),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    list(index = i, value = value, warnings = unique(warnings))
  }
  results <- if (workers > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(
      seq_along(tasks),
      run_one,
      mc.cores = min(workers, length(tasks)),
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_along(tasks), run_one)
  }
  .module3_validate_extraction_k_results(results, tasks)
  worker_warnings <- unique(unlist(lapply(results, `[[`, "warnings"), use.names = FALSE))
  if (length(worker_warnings)) {
    .log_warn("Module 3 K extraction warning(s): {paste(worker_warnings, collapse = '; ')}")
  }
  invisible(results)
}

.module3_extraction_requires_sequential <- function(topic_report_args,
                                                    assignment_method) {
  isTRUE(topic_report_args$run_link_topic_scores) ||
    isTRUE(topic_report_args$run_topic_assignment_qc) ||
    (
      identical(assignment_method, "gammafit_maxprob") &&
        !identical(topic_report_args$optimize_topics, FALSE)
    )
}

.module3_validate_extraction_k_results <- function(results, tasks) {
  delivered <- vapply(results, function(x) {
    is.list(x) &&
      all(c("index", "value", "warnings") %in% names(x)) &&
      length(x$index) == 1L &&
      is.finite(suppressWarnings(as.integer(x$index)))
  }, logical(1L))
  if (any(!delivered)) {
    missing_k <- vapply(
      which(!delivered),
      function(i) as.character(tasks[[i]]$k),
      character(1L)
    )
    .log_abort(
      paste0(
        "Module 3 K extraction worker(s) did not return a result for K=",
        paste(missing_k, collapse = ", "),
        ". Rerun extraction sequentially."
      )
    )
  }
  failures <- vapply(
    results,
    function(x) inherits(x$value, "error"),
    logical(1L)
  )
  if (any(failures)) {
    messages <- vapply(results[failures], function(x) {
      sprintf("K=%s: %s", tasks[[x$index]]$k, conditionMessage(x$value))
    }, character(1L))
    .log_abort(paste(c("Module 3 K extraction failed.", messages), collapse = "\n"))
  }
  invisible(TRUE)
}

#' Extract Module 3 regulatory topics
#'
#' Public step function for extracting regulatory topics, pathway summaries,
#' topic-link tables, direction-specific TF-to-topic assignment tables, and
#' review outputs from trained Module 3 topic models.
#'
#' Topic assignment uses unit-specific evidence. By default, GammaFit first
#' identifies candidate topics for each aggregate `GENE:<gene>` and
#' `PEAK:<gene>` pair. Each term independently selects its maximum-phi passing
#' topic, and the pair is retained only when those topics agree. The explicit
#' `"max_phi"` and legacy `"gammafit"` methods remain available for comparison.
#' TFs are
#' assigned from raw document-topic `theta` with the TF membership and primary
#' margin cutoffs. Per-comparison pathway gene sets are built from model
#' outputs only: documents with theta above the TF membership cutoff intersect
#' genes represented by topic-assigned `GENE:<gene>` and aggregate
#' `PEAK:<gene>` terms. Physical TF-peak-gene links are not used to define
#' pathway topic membership at extraction time; they can be projected later
#' onto selected comparison/topic/pathway genes for subnetworks.
#'
#' @param k Integer K value or vector of K values selected for extraction.
#' @param model_dir Directory containing trained topic model outputs.
#' @param output_dir Directory to write extracted topic outputs.
#' @param flatten_single_output Whether to write a single selected model
#'   directly under `output_dir`. Defaults to `TRUE` for the public step API.
#' @param topic_score_method Topic-term score method. `"normtop_specificity"`
#'   is the default for new extractions; `"rowmax_phi"` preserves the legacy
#'   row-maximum-scaled phi score.
#' @param topic_term_assignment_method Term-to-topic assignment.
#'   `"gammafit_maxprob"` applies GammaFit first, independently selects the
#'   maximum-phi passing topic for each aggregate Gene and Peak term, and keeps
#'   the target only when those topics agree. `"max_phi"` assigns terms
#'   independently from raw `phi`; `"gammafit"` retains multi-topic cutoff
#'   membership.
#' @param optimize_topics Whether eligible condition-topic extractions merge
#'   undersized or highly similar topics before downstream reports.
#' @param topic_merge_min_genes Minimum assigned genes required to retain a
#'   topic without a size-based merge.
#' @param topic_merge_min_links Minimum aligned TF-target links required to
#'   retain a topic without a size-based merge.
#' @param topic_merge_similarity_threshold Mean Gene/Peak Hellinger similarity
#'   at or above which two topics are merged.
#' @param run_topic_assignment_qc Whether to write the standard per-K topic
#'   assignment QC PDF and optimization audit tables.
#' @param topic_qc_umap_links_per_condition Maximum deterministic UMAP sample
#'   size per condition. Full-universe counts are never sampled.
#' @param topic_qc_top_tfs Number of globally ranked TFs shown in the pooled
#'   TF-by-topic QC heatmap.
#' @param topic_qc_condition_expression_file Optional complete condition-pair
#'   gene-expression audit CSV. When supplied, QC expression summaries use
#'   this complete matrix instead of the condition-filtered link rows.
#' @param topic_qc_reference_condition Optional reference condition used to
#'   count assigned target genes with increased expression in the QC report.
#' @param topic_qc_upregulated_log2fc_min Minimum log2 fold change above the
#'   reference condition for the upregulated-target QC panel.
#' @param topic_qc_upregulated_pseudocount Positive expression pseudocount used
#'   for the reference-condition log2 fold change.
#' @param k_workers Number of K values to extract concurrently. `NULL` selects
#'   workers adaptively from current memory and CPU headroom.
#' @param k_max_workers Maximum concurrent K workers in adaptive mode.
#' @param k_memory_gb Optional conservative memory estimate per K worker in GiB.
#'   When `NULL`, estimate from `edges_docs.rds` with a 16 GiB minimum.
#' @param k_memory_reserve_gb Minimum RAM in GiB to leave unused. Adaptive
#'   scheduling also reserves at least 25 percent of currently available RAM.
#' @param cores CPU cores available to the extraction scheduler.
#' @param verbose Emit concise extraction scheduler messages.
#' @param ... Additional arguments passed to the internal extraction engine,
#'   such as `backend`, `doc_mode`, `weight_label`, and `topic_report_args`.
#'
#' @return Invisibly returns TRUE when extraction completes.
#' @export
module3_extract_topics <- function(k,
                                   model_dir,
                                   output_dir,
                                   flatten_single_output = TRUE,
                                   topic_score_method = c("normtop_specificity", "rowmax_phi"),
                                   topic_term_assignment_method = c("gammafit_maxprob", "max_phi", "gammafit"),
                                   optimize_topics = NULL,
                                   topic_merge_min_genes = 50L,
                                   topic_merge_min_links = 200L,
                                   topic_merge_similarity_threshold = 0.90,
                                   run_topic_assignment_qc = NULL,
                                   topic_qc_umap_links_per_condition = 3000L,
                                   topic_qc_top_tfs = 150L,
                                   topic_qc_condition_expression_file = NULL,
                                   topic_qc_reference_condition = NULL,
                                   topic_qc_upregulated_log2fc_min = 1,
                                   topic_qc_upregulated_pseudocount = 1,
                                   k_workers = NULL,
                                   k_max_workers = 4L,
                                   k_memory_gb = NULL,
                                   k_memory_reserve_gb = 32,
                                   cores = NULL,
                                   verbose = TRUE,
                                   ...) {
  k <- sort(unique(as.integer(k)))
  k <- k[is.finite(k) & k > 1L]
  if (!length(k)) .log_abort("k must contain at least one integer greater than 1.")
  topic_score_method <- match.arg(topic_score_method)
  topic_term_assignment_method <- match.arg(topic_term_assignment_method)
  dots <- list(...)
  topic_report_args <- dots$topic_report_args %||% list()
  dots$topic_report_args <- modifyList(
    list(
      optimize_topics = optimize_topics,
      topic_merge_min_genes = topic_merge_min_genes,
      topic_merge_min_links = topic_merge_min_links,
      topic_merge_similarity_threshold = topic_merge_similarity_threshold,
      run_topic_assignment_qc = run_topic_assignment_qc,
      topic_qc_umap_links_per_condition = topic_qc_umap_links_per_condition,
      topic_qc_top_tfs = topic_qc_top_tfs,
      topic_qc_condition_expression_file = topic_qc_condition_expression_file,
      topic_qc_reference_condition = topic_qc_reference_condition,
      topic_qc_upregulated_log2fc_min = topic_qc_upregulated_log2fc_min,
      topic_qc_upregulated_pseudocount = topic_qc_upregulated_pseudocount
    ),
    topic_report_args
  )
  if (length(k) == 1L) {
    return(do.call(
      extract_regulatory_topics,
      c(list(
        k = k,
        model_dir = model_dir,
        output_dir = output_dir,
        flatten_single_output = flatten_single_output,
        topic_score_method = topic_score_method,
        topic_term_assignment_method = topic_term_assignment_method
      ), dots)
    ))
  }
  link_scoring <- isTRUE(dots$topic_report_args$run_link_topic_scores) ||
    !identical(dots$topic_report_args$optimize_topics, FALSE)
  plan <- .module3_extraction_k_worker_plan(
    n_tasks = length(k),
    model_dir = model_dir,
    k_workers = k_workers,
    k_max_workers = k_max_workers,
    k_memory_gb = k_memory_gb,
    k_memory_reserve_gb = k_memory_reserve_gb,
    cores = cores,
    link_scoring = link_scoring
  )
  if (isTRUE(verbose)) {
    .log_inform(
      "Module 3 K extraction: {plan$workers}/{length(k)} concurrent worker(s); estimated {(.format_bytes(plan$per_worker_bytes))} per worker; reserve {(.format_bytes(plan$reserve_bytes))}; {plan$reason}."
    )
  }
  tasks <- lapply(k, function(kk) list(
    k = kk,
    output_dir = file.path(output_dir, paste0("K", kk))
  ))
  worker_fun <- function(task) {
    do.call(
      extract_regulatory_topics,
      c(list(
        k = task$k,
        model_dir = model_dir,
        output_dir = task$output_dir,
        flatten_single_output = TRUE,
        topic_score_method = topic_score_method,
        topic_term_assignment_method = topic_term_assignment_method
      ), dots)
    )
  }
  .module3_run_extraction_k_tasks(
    tasks = tasks,
    worker_fun = worker_fun,
    workers = plan$workers,
    cores = plan$cores
  )
  invisible(TRUE)
}


# =============================================================================
# Package-side topic document construction overrides
# =============================================================================

# File: utils_step3_topic_models.R
# Author: Yaoxiang Li
# Created: 2026-05-13
#
# Purpose:
# Package-level utilities for Module 3 topic-model document construction.
# This file intentionally starts small. Add functions here gradually as they are
# needed by benchmark and package workflows instead of copying the full
# topic-model helper file.
#
# Design notes:
# - `doc_mode = "tf_cluster"` builds condition documents at the TF-cluster
#   level and preserves the previous condition-document behavior by default.
# - `doc_mode = "tf"` builds condition documents at the individual TF level:
#   `condition_label::tf`.
# - Condition-level thresholds are explicit arguments. Defaults are permissive
#   so existing workflows keep their previous behavior unless the caller opts in
#   to stricter config-derived thresholds.

#' Module 3 Topic-Model Helpers
#'
#' Internal helpers used to construct document-term inputs for topic models.
#'
#' @name topic_model_helpers
#' @noRd
NULL

.topic_safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

.topic_assert_has_cols <- function(df, cols, context = NULL) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    msg <- if (is.null(context)) {
      "Missing required columns."
    } else {
      paste0(context, " missing required columns.")
    }
    .log_abort(c(msg, i = paste(miss, collapse = ", ")))
  }
  invisible(TRUE)
}

#' Convert continuous topic weights to integer counts
#'
#' Converts non-negative continuous weights to integer counts for topic-model
#' training. The `"bin"` method rescales by the maximum observed weight, while
#' the `"log"` method applies `ceiling(log1p(weight) * scale)`.
#'
#' @param w Numeric vector of non-negative weights.
#' @param method Count conversion method, either `"bin"` or `"log"`.
#' @param scale Numeric scale factor.
#' @param min_count Minimum positive count to emit.
#'
#' @return Integer vector of counts with zeros for non-positive weights.
#' @noRd
weight_to_count <- function(w,
                            method = c("bin", "log"),
                            scale = 50,
                            min_count = 1L) {
  method <- match.arg(method)
  w <- .topic_safe_num(w)
  w[!is.finite(w) | w < 0] <- 0

  if (identical(method, "bin")) {
    mx <- suppressWarnings(max(w, na.rm = TRUE))
    if (!is.finite(mx) || mx <= 0) {
      return(rep.int(0L, length(w)))
    }
    cts <- ceiling((w / mx) * as.numeric(scale))
  } else {
    cts <- ceiling(base::log1p(w) * as.numeric(scale))
  }

  cts[w <= 0] <- 0L
  cts[cts > 0 & cts < min_count] <- min_count
  as.integer(cts)
}

.topic_condition_tf_doc <- function(tf, tf_cluster_map, doc_mode) {
  tf <- as.character(tf)
  if (identical(doc_mode, "tf")) {
    return(tf)
  }
  if (is.null(tf_cluster_map)) {
    .log_abort("doc_mode='tf_cluster' requires tf_cluster_map.")
  }
  map <- tf_cluster_map
  if (is.list(map) && !is.null(map$tf_cluster_map)) {
    map <- map$tf_cluster_map
  }
  if (!is.null(names(map))) {
    names(map) <- toupper(names(map))
    out <- unname(map[toupper(tf)])
  } else {
    out <- unname(map[tf])
  }
  out[is.na(out) | !nzchar(out)] <- tf[is.na(out) | !nzchar(out)]
  as.character(out)
}

#' Add comparison-level topic-model document IDs to links
#'
#' Assigns comparison-mode topic documents. For `doc_mode = "tf"`, documents
#' are individual TF documents (`comparison_id::tf::Target-Up` or
#' `comparison_id::tf::Target-Down`). For `doc_mode = "tf_cluster"`, documents
#' are TF-cluster documents using the same direction suffix.
#'
#' When `direction_by = "gene"`, `Target-Up` means condition 1 is the
#' direction-relevant condition and `Target-Down` means condition 2 is
#' the direction-relevant condition. The returned table includes
#' `tf_expr_condition`, `gene_expr_condition`, and `fp_score_condition` columns
#' for threshold-aware term construction.
#'
#' @param edges Link table.
#' @param doc_mode Document mode, either `"tf"` or `"tf_cluster"`.
#' @param direction_by Direction source, currently `"gene"`, `"fp"`, or
#'   `"none"` for compatibility.
#' @param tf_cluster_map Named TF-to-cluster vector required for
#'   `doc_mode = "tf_cluster"`.
#'
#' @return A data.table with `doc_id` and condition-aware score columns.
#' @noRd
add_tf_docs <- function(edges,
                        doc_mode = c("tf", "tf_cluster", "comparison"),
                        direction_by = c("gene", "fp", "none"),
                        tf_cluster_map = NULL) {
  .assert_pkg("data.table")
  doc_mode <- match.arg(doc_mode)
  direction_by <- match.arg(direction_by)

  dt <- data.table::as.data.table(edges)
  if (!"tf_expr_cond1" %in% names(dt) && "tf_expr_case" %in% names(dt)) dt[, tf_expr_cond1 := tf_expr_case]
  if (!"tf_expr_cond2" %in% names(dt) && "tf_expr_ctrl" %in% names(dt)) dt[, tf_expr_cond2 := tf_expr_ctrl]
  if (!"gene_expr_cond1" %in% names(dt) && "gene_expr_case" %in% names(dt)) dt[, gene_expr_cond1 := gene_expr_case]
  if (!"gene_expr_cond2" %in% names(dt) && "gene_expr_ctrl" %in% names(dt)) dt[, gene_expr_cond2 := gene_expr_ctrl]
  if (!"fp_score_cond1" %in% names(dt) && "fp_score_case" %in% names(dt)) dt[, fp_score_cond1 := fp_score_case]
  if (!"fp_score_cond2" %in% names(dt) && "fp_score_ctrl" %in% names(dt)) dt[, fp_score_cond2 := fp_score_ctrl]
  req <- c("comparison_id", "tf")
  if (identical(direction_by, "gene")) {
    req <- c(req, "log2fc_gene")
  }
  if (identical(direction_by, "fp")) {
    req <- c(req, "log2fc_fp")
  }
  .topic_assert_has_cols(dt, req, context = "add_tf_docs")
  dt[, comparison_id := as.character(comparison_id)]
  if (!"comparison_label" %in% names(dt)) {
    dt[, comparison_label := comparison_id]
  }
  dt[, comparison_label := as.character(comparison_label)]
  dt[is.na(comparison_label) | !nzchar(trimws(comparison_label)), comparison_label := comparison_id]

  if (identical(direction_by, "none")) {
    dt[, direction := NA_character_]
    direction_sign <- rep.int(0L, nrow(dt))
  } else if (identical(direction_by, "gene")) {
    direction_sign <- sign(.topic_safe_num(dt[["log2fc_gene"]]))
    dt[, direction := ifelse(direction_sign > 0L, "Target-Up", ifelse(direction_sign < 0L, "Target-Down", NA_character_))]
  } else {
    direction_sign <- sign(.topic_safe_num(dt[["log2fc_fp"]]))
    if ("delta_fp" %in% names(dt)) {
      alt <- sign(.topic_safe_num(dt[["delta_fp"]]))
      direction_sign[direction_sign == 0L] <- alt[direction_sign == 0L]
    }
    dt[, direction := ifelse(direction_sign > 0L, "FP-Up", ifelse(direction_sign < 0L, "FP-Down", NA_character_))]
  }

  if (identical(doc_mode, "tf")) {
    dt[, tf_doc := as.character(tf)]
  } else if (identical(doc_mode, "tf_cluster")) {
    dt[, tf_doc := .topic_condition_tf_doc(tf, tf_cluster_map, "tf_cluster")]
  } else {
    dt[, tf_doc := NA_character_]
  }

  if ("doc_modality" %in% names(dt)) {
    mod <- as.character(dt$doc_modality)
    has_mod <- !is.na(mod) & nzchar(mod)
    if (any(has_mod)) {
      dt[has_mod, tf_doc := paste(mod[has_mod], tf_doc[has_mod], sep = "::")]
    }
  }

  if (identical(doc_mode, "comparison")) {
    dt[, doc_id := if (identical(direction_by, "none")) comparison_id else paste(comparison_id, direction, sep = "::")]
    dt[, doc_display_label := if (identical(direction_by, "none")) comparison_label else paste(comparison_label, direction, sep = "::")]
  } else {
    dt[, doc_id := paste(comparison_id, tf_doc, direction, sep = "::")]
    dt[, doc_display_label := paste(comparison_label, direction, sep = "::")]
  }

  dt[, direction_sign := direction_sign]
  if (!identical(direction_by, "none")) {
    dt <- dt[!is.na(direction) & nzchar(direction)]
  }

  cond1_is_relevant <- dt$direction_sign > 0L
  cond2_is_relevant <- dt$direction_sign < 0L
  pick_direction_value <- function(cond1_col, cond2_col) {
    cond1_val <- if (cond1_col %in% names(dt)) .topic_safe_num(dt[[cond1_col]]) else rep(NA_real_, nrow(dt))
    cond2_val <- if (cond2_col %in% names(dt)) .topic_safe_num(dt[[cond2_col]]) else rep(NA_real_, nrow(dt))
    out <- rep(NA_real_, nrow(dt))
    out[cond1_is_relevant] <- cond1_val[cond1_is_relevant]
    out[cond2_is_relevant] <- cond2_val[cond2_is_relevant]
    out
  }
  dt[, `:=`(
    tf_expr_condition = pick_direction_value("tf_expr_cond1", "tf_expr_cond2"),
    gene_expr_condition = pick_direction_value("gene_expr_cond1", "gene_expr_cond2"),
    fp_score_condition = pick_direction_value("fp_score_cond1", "fp_score_cond2")
  )]
  dt[, direction_sign := NULL]

  dt[]
}

#' Add condition-level topic-model document IDs to links
#'
#' Expands per-comparison links into condition-specific rows for condition 1 and
#' condition 2. For `doc_mode = "tf_cluster"`, the document key is
#' `condition_label::tf_cluster`. For `doc_mode = "tf"`, the document key is
#' `condition_label::tf`.
#'
#' This helper does not itself filter target genes, peaks, or TF documents by
#' thresholds. It records condition-specific TF expression, target-gene
#' expression, and footprint score so `build_doc_term_condition_union()` can
#' apply threshold rules consistently during term construction.
#'
#' Required input columns are `comparison_id`, `cond1_id`, `cond2_id`, `tf`,
#' `gene_key`, `peak_id`, `fp_score_cond1`, `fp_score_cond2`,
#' `gene_expr_cond1`, and `gene_expr_cond2`. TF expression columns
#' `tf_expr_cond1` and `tf_expr_cond2` are optional for cluster-level
#' workflows but required when condition-level TF thresholds are applied.
#'
#' @param edges Link table.
#' @param tf_cluster_map Named vector mapping TF names to clusters. Required
#'   for `doc_mode = "tf_cluster"`.
#' @param doc_mode Document mode, either `"tf_cluster"` or `"tf"`.
#'
#' @return A data.table with condition-specific rows and `doc_id`.
#' @noRd
add_condition_tf_docs <- function(edges,
                                  tf_cluster_map = NULL,
                                  doc_mode = c("tf_cluster", "tf")) {
  .assert_pkg("data.table")
  doc_mode <- match.arg(doc_mode)

  dt <- data.table::as.data.table(edges)
  if (!"cond1_id" %in% names(dt) && "case_id" %in% names(dt)) dt[, cond1_id := case_id]
  if (!"cond2_id" %in% names(dt) && "ctrl_id" %in% names(dt)) dt[, cond2_id := ctrl_id]
  if (!"tf_expr_cond1" %in% names(dt) && "tf_expr_case" %in% names(dt)) dt[, tf_expr_cond1 := tf_expr_case]
  if (!"tf_expr_cond2" %in% names(dt) && "tf_expr_ctrl" %in% names(dt)) dt[, tf_expr_cond2 := tf_expr_ctrl]
  if (!"gene_expr_cond1" %in% names(dt) && "gene_expr_case" %in% names(dt)) dt[, gene_expr_cond1 := gene_expr_case]
  if (!"gene_expr_cond2" %in% names(dt) && "gene_expr_ctrl" %in% names(dt)) dt[, gene_expr_cond2 := gene_expr_ctrl]
  if (!"fp_score_cond1" %in% names(dt) && "fp_score_case" %in% names(dt)) dt[, fp_score_cond1 := fp_score_case]
  if (!"fp_score_cond2" %in% names(dt) && "fp_score_ctrl" %in% names(dt)) dt[, fp_score_cond2 := fp_score_ctrl]
  .topic_assert_has_cols(
    dt,
    c(
      "comparison_id", "cond1_id", "cond2_id", "tf", "gene_key", "peak_id",
      "fp_score_cond1", "fp_score_cond2", "gene_expr_cond1", "gene_expr_cond2"
    ),
    context = "add_condition_tf_docs"
  )

  dt[, tf_doc := .topic_condition_tf_doc(tf, tf_cluster_map, doc_mode)]

  cond1_dt <- data.table::copy(dt)
  cond1_dt[, `:=`(
    condition_label = as.character(cond1_id),
    fp_score_condition = .topic_safe_num(fp_score_cond1),
    gene_expr_condition = .topic_safe_num(gene_expr_cond1),
    tf_expr_condition = if ("tf_expr_cond1" %in% names(cond1_dt)) .topic_safe_num(tf_expr_cond1) else NA_real_
  )]

  cond2_dt <- data.table::copy(dt)
  cond2_dt[, `:=`(
    condition_label = as.character(cond2_id),
    fp_score_condition = .topic_safe_num(fp_score_cond2),
    gene_expr_condition = .topic_safe_num(gene_expr_cond2),
    tf_expr_condition = if ("tf_expr_cond2" %in% names(cond2_dt)) .topic_safe_num(tf_expr_cond2) else NA_real_
  )]

  out <- data.table::rbindlist(list(cond1_dt, cond2_dt), use.names = TRUE, fill = TRUE)
  out <- out[!is.na(condition_label) & nzchar(condition_label)]
  out[, doc_id := paste(condition_label, tf_doc, sep = "::")]
  out[]
}

#' Add condition-level TF-cluster document IDs to links
#'
#' Convenience wrapper for condition-level TF-cluster document construction used by
#' condition-level TF-cluster documents.
#'
#' @param edges Link table.
#' @param tf_cluster_map Named vector mapping TF names to clusters.
#'
#' @return A data.table with condition-specific TF-cluster document IDs.
#' @noRd
add_condition_tf_cluster_docs <- function(edges, tf_cluster_map) {
  add_condition_tf_docs(edges, tf_cluster_map = tf_cluster_map, doc_mode = "tf_cluster")
}

#' Build document-term rows from comparison-mode edges
#'
#' Converts comparison-mode topic documents into gene or peak term rows. The
#' function preserves the previous behavior by default. When condition-aware
#' thresholds are supplied, target genes and peaks are additionally filtered by
#' the direction-relevant condition values produced by `add_tf_docs()`.
#'
#' @param edges_docs Link table with `doc_id`, `tf`, `gene_key`, and `peak_id`.
#' @param term_type Term type, either `"peak"` or `"gene"`.
#' @param weight_type Edge weight column.
#' @param top_terms_per_doc Optional maximum terms per document.
#' @param min_df Minimum document frequency after document-term aggregation.
#' @param count_method Count conversion method.
#' @param count_scale Count scale.
#' @param prefix_terms Whether to prefix terms with `GENE:` or `PEAK:`.
#' @param distinct_terms Whether to keep distinct terms before aggregation.
#' @param gene_term_mode Gene aggregation mode.
#' @param include_tf_terms Whether to include current TFs as `GENE:<tf>` terms.
#' @param tf_weight_type Legacy TF weight column to use when condition-specific
#'   TF expression is absent.
#' @param threshold_gene_expr Minimum direction-relevant target expression.
#' @param threshold_fp_score Minimum direction-relevant footprint score.
#' @param threshold_tf_expr Minimum direction-relevant TF expression.
#' @param require_condition_thresholds Whether to enforce condition-specific
#'   gene/peak thresholds when the relevant columns are present.
#'
#' @return A document-term data.table.
#' @noRd
build_doc_term_from_edges <- function(edges_docs,
                                      term_type = c("peak", "gene"),
                                      weight_type = c("delta_fp", "fc_mag_fp", "fc_mag_gene", "log2fc_fp", "log2fc_gene"),
                                      top_terms_per_doc = 500L,
                                      min_df = 2L,
                                      count_method = c("bin", "log"),
                                      count_scale = 50,
                                      prefix_terms = TRUE,
                                      distinct_terms = FALSE,
                                      gene_term_mode = c("aggregate", "unique"),
                                      include_tf_terms = FALSE,
                                      tf_weight_type = c("fc_mag_tf", "log2fc_tf"),
                                      threshold_gene_expr = -Inf,
                                      threshold_fp_score = -Inf,
                                      threshold_tf_expr = -Inf,
                                      require_condition_thresholds = FALSE) {
  .assert_pkg("data.table")
  term_type <- match.arg(term_type)
  weight_type <- match.arg(weight_type)
  count_method <- match.arg(count_method)
  gene_term_mode <- match.arg(gene_term_mode)
  tf_weight_type <- match.arg(tf_weight_type)

  dt <- data.table::as.data.table(edges_docs)
  .topic_assert_has_cols(dt, c("doc_id", "tf", "gene_key", "peak_id"), context = "build_doc_term_from_edges")

  if (identical(term_type, "peak")) {
    .topic_assert_has_cols(dt, c("peak_id", weight_type), context = "build_doc_term_from_edges/peak")
    tmp <- dt[!is.na(peak_id) & nzchar(peak_id)]
    if (isTRUE(require_condition_thresholds)) {
      tmp <- tmp[
        is.finite(.topic_safe_num(fp_score_condition)) &
          .topic_safe_num(fp_score_condition) >= threshold_fp_score &
          is.finite(.topic_safe_num(gene_expr_condition)) &
          .topic_safe_num(gene_expr_condition) >= threshold_gene_expr &
          is.finite(.topic_safe_num(tf_expr_condition)) &
          .topic_safe_num(tf_expr_condition) >= threshold_tf_expr
      ]
    }
    tmp[, term_id := as.character(peak_id)]
  } else {
    .topic_assert_has_cols(dt, c("gene_key", weight_type), context = "build_doc_term_from_edges/gene")
    tmp <- dt[!is.na(gene_key) & nzchar(gene_key)]
    if (isTRUE(require_condition_thresholds)) {
      tmp <- tmp[
        is.finite(.topic_safe_num(fp_score_condition)) &
          .topic_safe_num(fp_score_condition) >= threshold_fp_score &
          is.finite(.topic_safe_num(gene_expr_condition)) &
          .topic_safe_num(gene_expr_condition) >= threshold_gene_expr &
          is.finite(.topic_safe_num(tf_expr_condition)) &
          .topic_safe_num(tf_expr_condition) >= threshold_tf_expr
      ]
    }
    tmp[, term_id := as.character(gene_key)]
  }

  if (!nrow(tmp)) {
    return(data.table::data.table())
  }
  tmp[, weight := abs(.topic_safe_num(get(weight_type)))]
  tmp <- tmp[is.finite(weight) & weight > 0]
  if (!nrow(tmp)) {
    return(data.table::data.table())
  }

  if (isTRUE(prefix_terms)) {
    tmp[, term_id := paste0(ifelse(identical(term_type, "peak"), "PEAK:", "GENE:"), term_id)]
  }

  if (identical(term_type, "gene")) {
    check_gene_values <- function(x, by_cols, label) {
      check <- x[, .(
        value_min = min(weight, na.rm = TRUE),
        value_max = max(weight, na.rm = TRUE)
      ), by = by_cols]
      check[, value_diff := abs(value_max - value_min)]
      check[, value_tol := 1e-8 * pmax(1, abs(value_min), abs(value_max))]
      bad <- check[is.finite(value_diff) & value_diff > value_tol]
      if (nrow(bad)) {
        ex <- bad[1]
        ex_txt <- paste(
          vapply(
            by_cols,
            function(col) paste0(col, "=", as.character(ex[[col]][1])),
            character(1)
          ),
          collapse = ", "
        )
        n_bad <- nrow(bad)
        .log_warn(
          "{n_bad} repeated comparison-mode gene {label} group(s) have inconsistent values; using the maximum only as a deterministic collapse reducer after this warning. Example: {ex_txt}"
        )
      }
      invisible(TRUE)
    }

    if (identical(gene_term_mode, "aggregate") && "peak_id" %in% names(tmp)) {
      check_gene_values(tmp, c("doc_id", "term_id", "peak_id"), "peak")
      check_gene_values(tmp, c("doc_id", "term_id"), "term")
      # Repeated rows here are multiple link observations projected to the same
      # comparison-specific gene term. A gene has one differential value per
      # comparison; max is only a deterministic reducer after the checks above.
      tmp <- tmp[, .(weight = max(weight, na.rm = TRUE)), by = .(doc_id, term_id, peak_id)]
      out <- tmp[, .(weight = sum(weight, na.rm = TRUE)), by = .(doc_id, term_id)]
    } else {
      check_gene_values(tmp, c("doc_id", "term_id"), "term")
      # A gene has one differential value per comparison document. Max is only
      # a deterministic reducer after the consistency check above.
      out <- tmp[, .(weight = max(weight, na.rm = TRUE)), by = .(doc_id, term_id)]
    }
  } else {
    if (isTRUE(distinct_terms)) {
      tmp <- tmp[, .(weight = max(weight, na.rm = TRUE)), by = .(doc_id, term_id)]
    }
    out <- tmp[, .(weight = sum(weight, na.rm = TRUE)), by = .(doc_id, term_id)]
  }

  if (isTRUE(include_tf_terms) && identical(term_type, "gene")) {
    tf_tmp <- dt[!is.na(tf) & nzchar(tf)]
    if (isTRUE(require_condition_thresholds) && "tf_expr_condition" %in% names(tf_tmp)) {
      tf_tmp <- tf_tmp[
        is.finite(.topic_safe_num(fp_score_condition)) &
          .topic_safe_num(fp_score_condition) >= threshold_fp_score &
          is.finite(.topic_safe_num(gene_expr_condition)) &
          .topic_safe_num(gene_expr_condition) >= threshold_gene_expr &
          is.finite(.topic_safe_num(tf_expr_condition)) &
          .topic_safe_num(tf_expr_condition) >= threshold_tf_expr
      ]
      tf_tmp[, weight := .topic_safe_num(tf_expr_condition)]
    } else {
      tf_col <- if (identical(tf_weight_type, "fc_mag_tf")) "fc_mag_tf" else "log2fc_tf"
      if (!tf_col %in% names(tf_tmp)) {
        .log_warn("TF weight column {tf_col} not found; skipping include_tf_terms.")
        tf_tmp <- tf_tmp[0]
      } else {
        tf_tmp[, weight := abs(.topic_safe_num(get(tf_col)))]
      }
    }
    if (nrow(tf_tmp)) {
      tf_tmp[, term_id := as.character(tf)]
      if (isTRUE(prefix_terms)) {
        tf_tmp[, term_id := paste0("GENE:", term_id)]
      }
      tf_tmp <- tf_tmp[is.finite(weight) & weight > 0]
      tf_out <- tf_tmp[, .(weight = max(weight, na.rm = TRUE)), by = .(doc_id, term_id)]
      out <- data.table::rbindlist(list(out, tf_out), use.names = TRUE, fill = TRUE)
    }
  }

  if (is.finite(top_terms_per_doc) && top_terms_per_doc > 0) {
    data.table::setorder(out, doc_id, -weight)
    out <- out[, head(.SD, as.integer(top_terms_per_doc)), by = doc_id]
  }

  df_tbl <- unique(out[, .(doc_id, term_id)])
  term_df <- df_tbl[, .N, by = term_id]
  keep_terms <- term_df[N >= as.integer(min_df), term_id]
  out <- out[term_id %in% keep_terms]
  if (!nrow(out)) {
    return(data.table::data.table())
  }

  out[, pseudo_count_bin := weight_to_count(weight, method = "bin", scale = count_scale)]
  out[, pseudo_count_log := weight_to_count(weight, method = "log", scale = count_scale)]
  out[, pseudo_count := if (identical(count_method, "bin")) pseudo_count_bin else pseudo_count_log]
  out[pseudo_count > 0]
}

#' Build peak terms from delta footprint scores
#'
#' Convenience wrapper for topic-report term construction.
#'
#' @inheritParams build_doc_term_from_edges
#' @return A document-term data.table.
#' @noRd
build_doc_term_opt1_peak_delta_fp <- function(edges_docs, ...) {
  build_doc_term_from_edges(edges_docs, term_type = "peak", weight_type = "delta_fp", ...)
}

#' Build peak terms from footprint fold-change magnitude
#'
#' Convenience wrapper for topic-report term construction.
#'
#' @inheritParams build_doc_term_from_edges
#' @return A document-term data.table.
#' @noRd
build_doc_term_opt2_peak_fc_fp <- function(edges_docs, ...) {
  build_doc_term_from_edges(edges_docs, term_type = "peak", weight_type = "fc_mag_fp", ...)
}

#' Build gene terms from expression fold-change magnitude
#'
#' Convenience wrapper for topic-report term construction.
#'
#' @inheritParams build_doc_term_from_edges
#' @return A document-term data.table.
#' @noRd
build_doc_term_opt3_gene_fc_expr <- function(edges_docs, ...) {
  build_doc_term_from_edges(edges_docs, term_type = "gene", weight_type = "fc_mag_gene", ...)
}

.topic_term_mode <- function(x) {
  x <- as.character(x)[1]
  switch(
    x,
    fp_uniq = "unique",
    fp_unique = "unique",
    uniq = "unique",
    unique = "unique",
    fp_aggr = "aggregate",
    fp_agg = "aggregate",
    agg = "aggregate",
    aggregate = "aggregate",
    fp_aggr_weight = "aggregate_weight",
    fp_agg_weight = "aggregate_weight",
    aggregate_weight = "aggregate_weight",
    gene_expression = "gene_expression",
    expression = "gene_expression",
    gene_only_expression = "gene_expression",
    .log_abort("Unsupported fp_term_mode: {x}")
  )
}

.topic_collapse_value <- function(x, value_col, by_cols, label, check_repeated_values = TRUE) {
  if (!nrow(x)) {
    return(data.table::data.table())
  }
  if (isTRUE(check_repeated_values)) {
    check <- x[, .(
      value_min = min(get(value_col), na.rm = TRUE),
      value_max = max(get(value_col), na.rm = TRUE)
    ), by = by_cols]
    check[, value_diff := abs(value_max - value_min)]
    check[, value_tol := 1e-8 * pmax(1, abs(value_min), abs(value_max))]
    bad <- check[is.finite(value_diff) & value_diff > value_tol]
    if (nrow(bad)) {
      ex <- bad[1]
      ex_txt <- paste(
        vapply(
          by_cols,
          function(col) paste0(col, "=", as.character(ex[[col]][1])),
          character(1)
        ),
        collapse = ", "
      )
      n_bad <- nrow(bad)
      .log_warn(
        "{n_bad} repeated {label} group(s) have inconsistent values; using the maximum only as a deterministic collapse reducer after this warning. Example: {ex_txt}"
      )
    }
  }

  # Repeated rows here are multiple link observations projected to the same
  # document-term value. They should carry the same biological value.
  # The max is only a deterministic reducer after the consistency check.
  x[, .(weight_raw = max(get(value_col), na.rm = TRUE)), by = by_cols]
}

.topic_apply_counts <- function(dt, count_method, count_scale) {
  dt[, pseudo_count_bin := weight_to_count(weight, method = "bin", scale = count_scale)]
  dt[, pseudo_count_log := weight_to_count(weight, method = "log", scale = count_scale)]
  dt[, pseudo_count := if (identical(count_method, "bin")) pseudo_count_bin else pseudo_count_log]
  dt[pseudo_count > 0]
}

.topic_balance_modalities <- function(dt, balance_mode) {
  if (!identical(balance_mode, "min") || !nrow(dt)) {
    dt[, weight := weight_raw]
    return(dt)
  }
  mods <- unique(dt$modality)
  if (!all(c("gene", "peak") %in% mods)) {
    dt[, weight := weight_raw]
    return(dt)
  }

  totals <- dt[, .(total = sum(weight_raw, na.rm = TRUE)), by = .(doc_id, modality)]
  totals_w <- data.table::dcast(totals, doc_id ~ modality, value.var = "total", fill = 0)
  if (!("gene" %in% names(totals_w))) {
    totals_w[, gene := 0]
  }
  if (!("peak" %in% names(totals_w))) {
    totals_w[, peak := 0]
  }
  totals_w[, target_total := pmin(gene, peak)]
  totals_w[, gene_scale := ifelse(gene > 0, target_total / gene, 0)]
  totals_w[, peak_scale := ifelse(peak > 0, target_total / peak, 0)]
  dt <- merge(dt, totals_w[, .(doc_id, gene_scale, peak_scale)], by = "doc_id", all.x = TRUE)
  dt[, scale := ifelse(modality == "gene", gene_scale, peak_scale)]
  dt[!is.finite(scale) | scale < 0, scale := 0]
  dt[, weight := weight_raw * scale]
  dt[, c("gene_scale", "peak_scale", "scale") := NULL]
  dt
}

.topic_condition_tf_expression_peak_factors <- function(dt) {
  .topic_assert_has_cols(
    dt,
    c("doc_id", "tf", "tf_expr_condition"),
    context = "condition TF-expression Peak weighting"
  )
  values <- unique(dt[, .(
    doc_id = as.character(doc_id),
    tf = as.character(tf),
    tf_expression = .topic_safe_num(tf_expr_condition)
  )])
  invalid <- values[
    is.na(doc_id) | !nzchar(doc_id) |
      is.na(tf) | !nzchar(tf) |
      !is.finite(tf_expression) | tf_expression <= 0
  ]
  if (nrow(invalid)) {
    .log_abort(
      "TF-expression Peak weighting requires one positive finite TF expression for every condition::TF document; found {nrow(invalid)} invalid row(s)."
    )
  }

  consistency <- values[, .(
    n_tfs = data.table::uniqueN(tf),
    expression_min = min(tf_expression),
    expression_max = max(tf_expression)
  ), by = doc_id]
  consistency[, expression_tolerance :=
    1e-8 * pmax(1, abs(expression_min), abs(expression_max))]
  inconsistent <- consistency[
    n_tfs != 1L |
      abs(expression_max - expression_min) > expression_tolerance
  ]
  if (nrow(inconsistent)) {
    .log_abort(
      "TF-expression Peak weighting requires one TF and one TF-expression value per document; found {nrow(inconsistent)} inconsistent document(s)."
    )
  }

  factors <- values[, .(
    tf = tf[[1L]],
    tf_expression = tf_expression[[1L]]
  ), by = doc_id]
  factors[, transformed_expression := log2(tf_expression + 1)]
  factors[, tf_expression_median := stats::median(transformed_expression), by = tf]
  if (any(!is.finite(factors$tf_expression_median) |
          factors$tf_expression_median <= 0)) {
    .log_abort("Unable to calculate positive per-TF expression medians for Peak weighting.")
  }
  factors[, peak_multiplier := transformed_expression / tf_expression_median]
  if (any(!is.finite(factors$peak_multiplier) | factors$peak_multiplier <= 0)) {
    .log_abort("TF-expression Peak weighting produced a non-positive or non-finite multiplier.")
  }

  tf_counts <- factors[, .(n_conditions = .N), by = tf]
  quantile_value <- function(x, probability) {
    unname(stats::quantile(x, probability, names = FALSE, na.rm = TRUE))
  }
  audit <- data.table::data.table(
    weighting = "tf_expression",
    scaling = "per_tf_log2p1_median",
    n_documents = nrow(factors),
    n_tfs = data.table::uniqueN(factors$tf),
    n_single_condition_tfs = tf_counts[n_conditions == 1L, .N],
    tf_expression_min = min(factors$tf_expression),
    tf_expression_p01 = quantile_value(factors$tf_expression, 0.01),
    tf_expression_median = stats::median(factors$tf_expression),
    tf_expression_p99 = quantile_value(factors$tf_expression, 0.99),
    tf_expression_max = max(factors$tf_expression),
    multiplier_min = min(factors$peak_multiplier),
    multiplier_p01 = quantile_value(factors$peak_multiplier, 0.01),
    multiplier_median = stats::median(factors$peak_multiplier),
    multiplier_p99 = quantile_value(factors$peak_multiplier, 0.99),
    multiplier_max = max(factors$peak_multiplier)
  )
  list(
    factors = factors[, .(doc_id, tf, tf_expression, peak_multiplier)],
    audit = audit
  )
}

.topic_modality_mass_audit <- function(dt, value_col, prefix) {
  mass <- dt[, .(mass = sum(get(value_col), na.rm = TRUE)), by = .(doc_id, modality)]
  mass <- data.table::dcast(mass, doc_id ~ modality, value.var = "mass", fill = 0)
  if (!"gene" %in% names(mass)) mass[, gene := 0]
  if (!"peak" %in% names(mass)) mass[, peak := 0]
  ratios <- mass[gene > 0, peak / gene]
  values <- list(
    sum(mass$gene),
    sum(mass$peak),
    if (length(ratios)) stats::median(ratios) else NA_real_,
    if (sum(mass$gene) > 0) sum(mass$peak) / sum(mass$gene) else NA_real_
  )
  names(values) <- paste0(
    prefix,
    c("_gene_mass", "_peak_mass", "_median_peak_gene_ratio", "_total_peak_gene_ratio")
  )
  data.table::as.data.table(values)
}

.topic_filter_term_df <- function(dt, min_df) {
  if (!nrow(dt) || !is.finite(min_df) || min_df <= 1) {
    return(dt)
  }
  df_tbl <- unique(dt[, .(doc_id, term_id)])
  term_df <- df_tbl[, .N, by = term_id]
  keep_terms <- term_df[N >= as.integer(min_df), term_id]
  dt[term_id %in% keep_terms]
}

.topic_trim_terms <- function(dt, top_terms_per_doc) {
  if (!nrow(dt) || !is.finite(top_terms_per_doc) || top_terms_per_doc <= 0) {
    return(dt)
  }
  data.table::setorder(dt, doc_id, -weight_raw, term_id)
  dt[, head(.SD, as.integer(top_terms_per_doc)), by = doc_id]
}

# Shared implementation for comparison and condition document-term builders.
# Signed differential columns are converted to positive model weights with
# abs(); direction remains encoded by the document ID.
.topic_build_doc_term_unified <- function(edges_docs,
                                          weight_type_peak,
                                          weight_type_gene,
                                          weight_type_tf = NULL,
                                          fp_term_mode = c("aggregate", "unique", "aggregate_weight", "gene_expression"),
                                          condition_peak_weighting = c("none", "tf_expression"),
                                          include_tf_terms = FALSE,
                                          count_method = c("bin", "log"),
                                          count_scale = 50,
                                          top_terms_per_doc = Inf,
                                          min_df = 1L,
                                          balance_mode = c("min", "none"),
                                          prefix_terms = TRUE,
                                          check_repeated_values = TRUE) {
  .assert_pkg("data.table")
  count_method <- match.arg(count_method)
  balance_mode <- match.arg(balance_mode)
  fp_term_mode <- .topic_term_mode(fp_term_mode)
  condition_peak_weighting <- match.arg(condition_peak_weighting)
  if (!identical(condition_peak_weighting, "none") &&
      !identical(fp_term_mode, "aggregate")) {
    .log_abort("TF-expression Peak weighting requires fp_term_mode = 'aggregate'.")
  }

  dt <- data.table::as.data.table(edges_docs)
  .topic_assert_has_cols(dt, c("doc_id", "tf", "gene_key", "peak_id", weight_type_peak, weight_type_gene), context = "build_doc_term_unified")
  peak_factor_result <- if (identical(condition_peak_weighting, "tf_expression")) {
    .topic_condition_tf_expression_peak_factors(dt)
  } else {
    NULL
  }
  if (identical(fp_term_mode, "gene_expression")) {
    keep <- unique(c(
      "doc_id", "tf", "gene_key", weight_type_gene,
      if (isTRUE(include_tf_terms)) weight_type_tf else NULL
    ))
    keep <- keep[!is.na(keep) & keep %in% names(dt)]
    dt <- dt[, keep, with = FALSE]
    dt[, `:=`(peak_id = NA_character_, fp_weight_source__ = 0)]
    weight_type_peak <- "fp_weight_source__"
  }
  dt[, `:=`(
    doc_id = as.character(doc_id),
    tf = as.character(tf),
    gene_key = as.character(gene_key),
    peak_id = as.character(peak_id),
    fp_weight = abs(.topic_safe_num(get(weight_type_peak))),
    gene_weight = abs(.topic_safe_num(get(weight_type_gene)))
  )]

  gene_dt <- dt[!is.na(gene_key) & nzchar(gene_key) & is.finite(gene_weight) & gene_weight > 0]
  if (nrow(gene_dt)) {
    gene_obs <- .topic_collapse_value(
      gene_dt,
      value_col = "gene_weight",
      by_cols = c("doc_id", "gene_key"),
      label = "gene document-term",
      check_repeated_values = check_repeated_values
    )
    data.table::setnames(gene_obs, "gene_key", "term_id")
    if (isTRUE(prefix_terms)) {
      gene_obs[, term_id := paste0("GENE:", term_id)]
    }
    gene_obs[, modality := "gene"]
  } else {
    gene_obs <- data.table::data.table(doc_id = character(), term_id = character(), weight_raw = numeric(), modality = character())
  }

  if (identical(fp_term_mode, "gene_expression")) {
    peak_obs <- data.table::data.table(
      doc_id = character(), term_id = character(),
      weight_raw = numeric(), modality = character()
    )
  } else {
    fp_dt <- dt[
      !is.na(peak_id) & nzchar(peak_id) &
        !is.na(gene_key) & nzchar(gene_key) &
        is.finite(fp_weight) & fp_weight > 0
    ]
  if (identical(fp_term_mode, "unique")) {
    if (nrow(fp_dt)) {
      peak_obs <- .topic_collapse_value(
        fp_dt,
        value_col = "fp_weight",
        by_cols = c("doc_id", "peak_id"),
        label = "peak document-term",
        check_repeated_values = check_repeated_values
      )
      data.table::setnames(peak_obs, "peak_id", "term_id")
      if (isTRUE(prefix_terms)) {
        peak_obs[, term_id := paste0("PEAK:", term_id)]
      }
      peak_obs[, modality := "peak"]
    } else {
      peak_obs <- data.table::data.table(doc_id = character(), term_id = character(), weight_raw = numeric(), modality = character())
    }
  } else {
    if (nrow(fp_dt)) {
      peak_by_gene_peak <- .topic_collapse_value(
        fp_dt,
        value_col = "fp_weight",
        by_cols = c("doc_id", "gene_key", "peak_id"),
        label = "peak-gene document-term",
        check_repeated_values = check_repeated_values
      )
      peak_obs <- peak_by_gene_peak[, .(weight_raw = sum(weight_raw, na.rm = TRUE)), by = .(doc_id, gene_key)]
      data.table::setnames(peak_obs, "gene_key", "term_id")
      if (isTRUE(prefix_terms)) {
        peak_obs[, term_id := paste0("PEAK:", term_id)]
      }
      peak_obs[, modality := "peak"]
    } else {
      peak_obs <- data.table::data.table(doc_id = character(), term_id = character(), weight_raw = numeric(), modality = character())
    }
  }
  }

  if (identical(fp_term_mode, "gene_expression")) {
    balance_mode <- "none"
  } else if (identical(fp_term_mode, "aggregate_weight")) {
    if (nrow(gene_obs) && nrow(peak_obs)) {
      peak_weight <- data.table::copy(peak_obs)
      peak_weight[, gene_term_id := sub("^PEAK:", "GENE:", term_id)]
      gene_obs <- merge(
        gene_obs,
        peak_weight[, .(doc_id, term_id = gene_term_id, fp_weight_raw = weight_raw)],
        by = c("doc_id", "term_id"),
        all = FALSE
      )
      gene_obs[, weight_raw := weight_raw * fp_weight_raw]
      gene_obs[, fp_weight_raw := NULL]
    } else {
      gene_obs <- gene_obs[0]
    }
    peak_obs <- peak_obs[0]
    balance_mode <- "none"
  }

  tf_obs <- data.table::data.table(doc_id = character(), term_id = character(), weight_raw = numeric(), modality = character())
  if (isTRUE(include_tf_terms)) {
    tf_dt <- dt[!is.na(tf) & nzchar(tf)]
    if (nrow(tf_dt)) {
      if (!is.null(weight_type_tf) && weight_type_tf %in% names(tf_dt)) {
        tf_dt[, tf_weight := abs(.topic_safe_num(get(weight_type_tf)))]
      } else if ("tf_expr_condition" %in% names(tf_dt)) {
        tf_dt[, tf_weight := .topic_safe_num(tf_expr_condition)]
      } else if ("fc_mag_tf" %in% names(tf_dt)) {
        tf_dt[, tf_weight := abs(.topic_safe_num(fc_mag_tf))]
      } else if ("log2fc_tf" %in% names(tf_dt)) {
        tf_dt[, tf_weight := abs(.topic_safe_num(log2fc_tf))]
      } else {
        tf_dt[, tf_weight := NA_real_]
      }
      tf_dt <- tf_dt[is.finite(tf_weight) & tf_weight > 0]
      if (nrow(tf_dt)) {
        tf_obs <- .topic_collapse_value(
          tf_dt,
          value_col = "tf_weight",
          by_cols = c("doc_id", "tf"),
          label = "TF self document-term",
          check_repeated_values = check_repeated_values
        )
        data.table::setnames(tf_obs, "tf", "term_id")
        if (isTRUE(prefix_terms)) {
          tf_obs[, term_id := paste0("GENE:", term_id)]
        }
        tf_obs[, modality := "gene"]
      }
    }
  }

  out <- data.table::rbindlist(list(gene_obs, tf_obs, peak_obs), use.names = TRUE, fill = TRUE)
  if (!nrow(out)) {
    return(data.table::data.table())
  }
  out <- out[is.finite(weight_raw) & weight_raw > 0]
  out <- .topic_trim_terms(out, top_terms_per_doc)
  out <- .topic_filter_term_df(out, min_df)
  if (!nrow(out)) {
    return(data.table::data.table())
  }
  out <- .topic_balance_modalities(out, balance_mode)
  out <- out[is.finite(weight) & weight > 0]
  if (!nrow(out)) {
    return(data.table::data.table())
  }
  peak_weighting_audit <- NULL
  if (!is.null(peak_factor_result)) {
    before_mass <- .topic_modality_mass_audit(
      out,
      value_col = "weight",
      prefix = "balanced_before_tf_weighting"
    )
    out[
      peak_factor_result$factors,
      peak_multiplier := i.peak_multiplier,
      on = "doc_id"
    ]
    unmatched <- out[modality == "peak" &
      (!is.finite(peak_multiplier) | peak_multiplier <= 0)]
    if (nrow(unmatched)) {
      .log_abort(
        "TF-expression Peak weighting could not match {data.table::uniqueN(unmatched$doc_id)} Peak document(s) to a multiplier."
      )
    }
    out[modality == "peak", weight := weight * peak_multiplier]
    out[, peak_multiplier := NULL]
    after_mass <- .topic_modality_mass_audit(
      out,
      value_col = "weight",
      prefix = "after_tf_weighting"
    )
    peak_weighting_audit <- cbind(
      peak_factor_result$audit,
      before_mass,
      after_mass
    )
  }
  out <- .topic_apply_counts(out, count_method = count_method, count_scale = count_scale)
  if (!is.null(peak_weighting_audit)) {
    count_mass <- .topic_modality_mass_audit(
      out,
      value_col = "pseudo_count",
      prefix = "after_count_conversion"
    )
    peak_weighting_audit <- cbind(peak_weighting_audit, count_mass)
  }
  result <- out[, .(
    weight = sum(weight, na.rm = TRUE),
    pseudo_count = sum(pseudo_count, na.rm = TRUE),
    pseudo_count_bin = sum(pseudo_count_bin, na.rm = TRUE),
    pseudo_count_log = sum(pseudo_count_log, na.rm = TRUE)
  ), by = .(doc_id, term_id)]
  if (!is.null(peak_weighting_audit)) {
    attr(result, "condition_peak_weighting_audit") <- peak_weighting_audit
  }
  result
}

#' Build joint gene-plus-peak document-term rows
#'
#' Builds comparison-level document-term rows from links that already have a
#' `doc_id`. Gene, peak, and optional TF self-term weights are converted to
#' positive values before count conversion. `fp_term_mode` controls whether FP
#' terms remain unique peaks, are aggregated by target gene, are multiplied
#' into gene weights, or are excluded for a gene-expression-only experiment.
#'
#' @inheritParams build_doc_term_from_edges
#' @param weight_type_peak Peak weight column, usually `log2fc_fp` for
#'   comparison-mode runs.
#' @param weight_type_gene Gene weight column, usually `log2fc_gene` for
#'   comparison-mode runs.
#' @param fp_term_mode Term mode: unique peaks, target-gene aggregated peaks,
#'   gene-only aggregate weights, or gene-expression-only terms.
#' @param condition_peak_weighting Optional condition-mode aggregated-Peak
#'   weighting. `"tf_expression"` multiplies balanced Peak weights by the
#'   condition-specific log TF expression relative to that TF's median across
#'   conditions before count conversion.
#' @param balance_mode Modality balancing mode, either `"min"` or `"none"`.
#' @param check_repeated_values Warn if repeated document-term rows have
#'   inconsistent source values before deterministic collapse.
#'
#' @return A data.table with `doc_id`, `term_id`, `weight`, and count columns.
#' @noRd
build_doc_term_joint <- function(edges_docs,
                                 weight_type_peak = c("delta_fp", "fc_mag_fp", "log2fc_fp"),
                                 weight_type_gene = c("fc_mag_gene", "log2fc_gene"),
                                 top_terms_per_doc = 500L,
                                 min_df = 2L,
                                 count_method = c("bin", "log"),
                                 count_scale = 50,
                                 distinct_terms = FALSE,
                                 gene_term_mode = c("aggregate", "unique"),
                                 fp_term_mode = NULL,
                                 include_tf_terms = FALSE,
                                 tf_weight_type = c("fc_mag_tf", "log2fc_tf"),
                                 balance_mode = c("min", "none"),
                                 prefix_terms = TRUE,
                                 threshold_gene_expr = -Inf,
                                 threshold_fp_score = -Inf,
                                 threshold_tf_expr = -Inf,
                                 require_condition_thresholds = FALSE,
                                 check_repeated_values = TRUE) {
  weight_type_peak <- match.arg(weight_type_peak)
  weight_type_gene <- match.arg(weight_type_gene)
  count_method <- match.arg(count_method)
  tf_weight_type <- match.arg(tf_weight_type)
  balance_mode <- match.arg(balance_mode)
  if (is.null(fp_term_mode)) {
    fp_term_mode <- if (isTRUE(distinct_terms)) "unique" else "aggregate"
  }

  dt <- data.table::as.data.table(edges_docs)
  if (isTRUE(require_condition_thresholds)) {
    dt <- dt[
      is.finite(.topic_safe_num(fp_score_condition)) &
        .topic_safe_num(fp_score_condition) >= threshold_fp_score &
        is.finite(.topic_safe_num(gene_expr_condition)) &
        .topic_safe_num(gene_expr_condition) >= threshold_gene_expr &
        is.finite(.topic_safe_num(tf_expr_condition)) &
        .topic_safe_num(tf_expr_condition) >= threshold_tf_expr
    ]
  }
  if (!nrow(dt)) {
    return(data.table::data.table())
  }

  .topic_build_doc_term_unified(
    dt,
    weight_type_peak = weight_type_peak,
    weight_type_gene = weight_type_gene,
    weight_type_tf = if (identical(tf_weight_type, "fc_mag_tf")) "fc_mag_tf" else "log2fc_tf",
    fp_term_mode = fp_term_mode,
    include_tf_terms = include_tf_terms,
    count_method = count_method,
    count_scale = count_scale,
    top_terms_per_doc = top_terms_per_doc,
    min_df = min_df,
    balance_mode = balance_mode,
    prefix_terms = prefix_terms,
    check_repeated_values = check_repeated_values
  )
}

#' Build condition-level union document-term rows
#'
#' Builds condition-level document terms after `add_condition_tf_docs()`.
#' Target genes require `gene_expr_condition >= threshold_gene_expr`, FP terms
#' require `fp_score_condition >= threshold_fp_score`, and TF documents can be
#' filtered by `tf_expr_condition >= threshold_tf_expr`. The same document
#' construction logic is used for TF and TF-cluster condition documents.
#'
#' @param edges_condition Condition-specific link table from
#'   `add_condition_tf_docs()`.
#' @param count_method Count conversion method, either `"bin"` or `"log"`.
#' @param count_scale Numeric count scale.
#' @param prefix_terms Whether to prefix terms with `GENE:` or `PEAK:`.
#' @param threshold_gene_expr Minimum target-gene expression.
#' @param threshold_fp_score Minimum footprint score.
#' @param threshold_tf_expr Minimum TF expression.
#' @param include_tf_terms Whether to add the current TF as a gene-like
#'   self-term in TF documents.
#' @param require_tf_expr Whether the document row requires TF expression.
#' @param fp_term_mode Term mode: unique peaks, target-gene aggregated peaks,
#'   gene-only aggregate weights, or gene-expression-only terms.
#' @param balance_mode Modality balancing mode, either `"min"` or `"none"`.
#' @param check_repeated_values Warn if repeated document-term rows have
#'   inconsistent source values before deterministic collapse.
#'
#' @return A data.table with `doc_id`, `term_id`, `weight`, and count columns.
#' @noRd
build_doc_term_condition_union <- function(edges_condition,
                                           count_method = c("bin", "log"),
                                           count_scale = 50,
                                           prefix_terms = TRUE,
                                           threshold_gene_expr = 0,
                                           threshold_fp_score = 0,
                                           threshold_tf_expr = -Inf,
                                           include_tf_terms = FALSE,
                                           require_tf_expr = FALSE,
                                           fp_term_mode = c("unique", "aggregate", "aggregate_weight", "gene_expression"),
                                           condition_peak_weighting = c("none", "tf_expression"),
                                           balance_mode = c("min", "none"),
                                           check_repeated_values = TRUE) {
  .assert_pkg("data.table")
  count_method <- match.arg(count_method)
  balance_mode <- match.arg(balance_mode)
  fp_term_mode <- .topic_term_mode(fp_term_mode)
  condition_peak_weighting <- match.arg(condition_peak_weighting)

  dt <- data.table::as.data.table(edges_condition)
  .topic_assert_has_cols(
    dt,
    c("condition_label", "tf_doc", "gene_key", "peak_id", "fp_score_condition", "gene_expr_condition"),
    context = "build_doc_term_condition_union"
  )

  chr <- function(x) {
    if (is.list(x)) {
      return(vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1)))
    }
    as.character(x)
  }

  dt[, `:=`(
    condition_label = chr(condition_label),
    tf_doc = chr(tf_doc),
    tf = if ("tf" %in% names(dt)) chr(tf) else NA_character_,
    gene_key = chr(gene_key),
    peak_id = chr(peak_id),
    fp_score_condition = .topic_safe_num(fp_score_condition),
    gene_expr_condition = .topic_safe_num(gene_expr_condition),
    tf_expr_condition = if ("tf_expr_condition" %in% names(dt)) .topic_safe_num(tf_expr_condition) else NA_real_
  )]

  if (isTRUE(require_tf_expr)) {
    dt <- dt[is.finite(tf_expr_condition) & tf_expr_condition >= threshold_tf_expr]
  }
  dt <- dt[
    is.finite(gene_expr_condition) & gene_expr_condition >= threshold_gene_expr &
      is.finite(fp_score_condition) & fp_score_condition >= threshold_fp_score
  ]
  if (!nrow(dt)) {
    return(data.table::data.table())
  }
  dt[, doc_id := paste(condition_label, tf_doc, sep = "::")]

  out <- .topic_build_doc_term_unified(
    dt,
    weight_type_peak = "fp_score_condition",
    weight_type_gene = "gene_expr_condition",
    weight_type_tf = "tf_expr_condition",
    fp_term_mode = fp_term_mode,
    condition_peak_weighting = condition_peak_weighting,
    include_tf_terms = include_tf_terms,
    count_method = count_method,
    count_scale = count_scale,
    top_terms_per_doc = Inf,
    min_df = 1L,
    balance_mode = balance_mode,
    prefix_terms = prefix_terms,
    check_repeated_values = check_repeated_values
  )

  if (identical(fp_term_mode, "aggregate") && isTRUE(prefix_terms)) {
    doc_tf <- unique(dt[, .(doc_id, tf_doc)])
    peak_terms <- unique(out[grepl("^PEAK:", term_id), .(
      doc_id,
      key = sub("^PEAK:", "", term_id)
    )])
    target_terms <- unique(out[grepl("^GENE:", term_id), .(
      doc_id,
      key = sub("^GENE:", "", term_id)
    )])
    target_terms <- doc_tf[target_terms, on = "doc_id"]
    # TF self terms are expected gene-like terms and do not have matching
    # aggregated peak terms. Only remove the same-document TF self term when
    # that TF is not also present as a real aggregated-peak target.
    target_terms[, has_matching_peak := paste(doc_id, key, sep = "\r") %in%
      paste(peak_terms$doc_id, peak_terms$key, sep = "\r")]
    target_terms <- target_terms[is.na(tf_doc) | key != tf_doc | has_matching_peak]
    target_terms[, has_matching_peak := NULL]
    target_counts <- target_terms[, .(N_gene = .N), by = doc_id]
    peak_counts <- peak_terms[, .(N_peak = .N), by = doc_id]
    check <- merge(target_counts, peak_counts, by = "doc_id", all = TRUE)
    check[is.na(N_gene), N_gene := 0L]
    check[is.na(N_peak), N_peak := 0L]
    check[, diff := N_gene - N_peak]
    bad <- check[diff != 0L]
    if (nrow(bad)) {
      max_abs_diff <- max(abs(bad$diff))
      .log_warn(
        "{nrow(bad)} fp_aggr document(s) have unmatched target-gene and aggregated-peak term counts after excluding expected TF self terms; max abs difference = {max_abs_diff}."
      )
    }
  }

  out
}

.module3_document_term_qc_summary <- function(doc_term, count_column = NULL) {
  dt <- data.table::as.data.table(doc_term)
  .topic_assert_has_cols(dt, c("doc_id", "term_id"), context = "document-term QC")
  count_column <- as.character(count_column %||% "")
  if (!length(count_column) || !nzchar(count_column) || !count_column %in% names(dt)) {
    count_column <- if ("count" %in% names(dt)) "count" else if ("weight" %in% names(dt)) "weight" else NULL
  }
  dt <- dt[startsWith(term_id, "GENE:") | startsWith(term_id, "PEAK:")]
  if (!nrow(dt)) return(data.table::data.table())
  dt[, `:=`(
    condition_id = sub("::[^:]+$", "", as.character(doc_id)),
    term_type = ifelse(startsWith(term_id, "GENE:"), "Gene", "Peak"),
    term_value = if (is.null(count_column)) 1 else .topic_safe_num(get(count_column))
  )]
  dt[!is.finite(term_value) | term_value < 0, term_value := 0]
  dt[, .(
    n_terms = data.table::uniqueN(term_id),
    term_mass = sum(term_value, na.rm = TRUE)
  ), by = .(condition_id, doc_id, term_type)]
}

.write_module3_document_term_qc <- function(doc_term,
                                            output_dir,
                                            count_column = NULL,
                                            title = "Module 3 document-term QC",
                                            verbose = TRUE) {
  .assert_pkg("ggplot2")
  summary <- .module3_document_term_qc_summary(doc_term, count_column = count_column)
  if (!nrow(summary)) {
    if (isTRUE(verbose)) .log_warn("Document-term QC skipped because no gene or peak terms were available.")
    return(invisible(NULL))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary_path <- file.path(output_dir, "document_term_qc_summary.csv")
  pdf_path <- file.path(output_dir, "document_term_qc.pdf")
  data.table::fwrite(summary, summary_path)
  condition_levels <- unique(summary$condition_id)
  title <- paste(strwrap(as.character(title[[1L]]), width = 82L), collapse = "\n")
  summary[, condition_id := factor(condition_id, levels = condition_levels)]
  summary[, term_type := factor(term_type, levels = c("Gene", "Peak"))]
  colors <- c(Gene = "#4477AA", Peak = "#EE7733")
  ncol_facets <- ceiling(sqrt(length(condition_levels)))
  base_theme <- ggplot2::theme_bw(base_size = 9, base_family = "Helvetica") +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "#111111"),
      axis.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.text = ggplot2::element_text(size = 9, face = "bold", color = "#111111"),
      strip.text = ggplot2::element_text(size = 9, face = "bold"),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(8, "pt"),
      aspect.ratio = 1,
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  distribution_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = term_type, y = n_terms, fill = term_type)
  )
  smallest_group <- min(summary[, .N, by = .(condition_id, term_type)]$N)
  if (smallest_group >= 2L) {
    distribution_plot <- distribution_plot + ggplot2::geom_violin(
      width = 0.82,
      trim = TRUE,
      color = "#222222",
      linewidth = 0.35,
      alpha = 0.82
    )
  } else {
    distribution_plot <- distribution_plot + ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      shape = 21,
      size = 2,
      color = "#222222",
      alpha = 0.82
    )
  }
  distribution_plot <- distribution_plot +
    ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", color = "#111111", linewidth = 0.45) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 21, size = 2.2, fill = "#111111", color = "white") +
    ggplot2::facet_wrap(~condition_id, ncol = ncol_facets, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = "Term counts per TF document; each condition uses its own y-axis",
      x = "Term type",
      y = "Distinct terms per document",
      fill = "Term type"
    ) +
    base_theme
  balance <- summary[, .(
    documents = as.double(data.table::uniqueN(doc_id)),
    document_term_rows = sum(as.double(n_terms)),
    median_terms_per_document = as.double(stats::median(n_terms)),
    total_term_mass = sum(as.double(term_mass))
  ), by = .(condition_id, term_type)]
  balance_plot <- ggplot2::ggplot(
    balance,
    ggplot2::aes(x = term_type, y = document_term_rows, fill = term_type)
  ) +
    ggplot2::geom_col(width = 0.68, color = "#222222", linewidth = 0.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(document_term_rows, big.mark = ",", scientific = FALSE)),
      vjust = -0.25,
      size = 3.1,
      family = "Helvetica",
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(~condition_id, ncol = ncol_facets, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::labs(
      title = paste0(title, " - condition and modality balance"),
      subtitle = "Total distinct document-term rows after all document filters",
      x = "Term type",
      y = "Document-term rows",
      fill = "Term type"
    ) +
    base_theme
  grDevices::cairo_pdf(
    pdf_path,
    width = max(8, 4.8 * ncol_facets),
    height = max(5.5, 4.8 * ceiling(length(condition_levels) / ncol_facets)),
    family = "Helvetica",
    bg = "white",
    onefile = TRUE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  print(distribution_plot)
  print(balance_plot)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  if (isTRUE(verbose)) .log_inform("Wrote document-term QC: {pdf_path}")
  invisible(list(pdf = pdf_path, summary = summary_path, balance = balance))
}

write_doc_term_cache <- function(doc_term,
                                 out_dir,
                                 save_full_doc_term_csv = FALSE,
                                 write_arrow = TRUE,
                                 arrow_name = "doc_term.arrow",
                                 preview_name = "doc_term_first100.csv") {
  .assert_pkg("data.table")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dt <- data.table::as.data.table(doc_term)
  arrow_path <- file.path(out_dir, arrow_name)
  preview_path <- file.path(out_dir, preview_name)

  if (isTRUE(write_arrow) && requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_feather(dt, arrow_path)
  } else if (isTRUE(write_arrow)) {
    .log_warn("Package {.pkg arrow} is not installed; skipping {arrow_path}.")
    arrow_path <- NA_character_
  } else {
    arrow_path <- NA_character_
  }

  data.table::fwrite(utils::head(dt, 100L), preview_path)
  csv_path <- NA_character_
  if (isTRUE(save_full_doc_term_csv)) {
    csv_path <- file.path(out_dir, "doc_term.csv")
    data.table::fwrite(dt, csv_path)
  }
  invisible(list(arrow = arrow_path, preview = preview_path, csv = csv_path))
}
