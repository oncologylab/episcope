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
    "doc_topic_heatmaps",
    "tf_topic_heatmaps",
    "topic_by_comparison",
    "topic_marker_heatmap",
    "pathway",
    "intertopic_distance",
    "ldavis"
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
      "doc_topic_heatmaps",
      "tf_topic_heatmaps",
      "topic_by_comparison",
      "topic_marker_heatmap",
      "pathway"
    ),
    core = c("topic_terms", "topic_links", "gammafit_summary"),
    plots = c(
      "gammafit_summary",
      "doc_topic_heatmaps",
      "tf_topic_heatmaps",
      "topic_by_comparison",
      "topic_marker_heatmap",
      "intertopic_distance",
      "ldavis"
    ),
    heatmaps = c("doc_topic_heatmaps", "tf_topic_heatmaps", "topic_by_comparison", "topic_marker_heatmap"),
    reports = c("gammafit_summary", "doc_topic_heatmaps", "tf_topic_heatmaps", "topic_by_comparison", "topic_marker_heatmap", "pathway")
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

.format_bytes <- function(bytes) {
  bytes <- suppressWarnings(as.numeric(bytes[[1L]]))
  if (!is.finite(bytes) || bytes < 0) return("unknown")
  units <- c("B", "KB", "MB", "GB", "TB")
  unit_i <- 1L
  while (bytes >= 1024 && unit_i < length(units)) {
    bytes <- bytes / 1024
    unit_i <- unit_i + 1L
  }
  sprintf("%.2f %s", bytes, units[[unit_i]])
}

.available_memory_bytes <- function() {
  meminfo <- "/proc/meminfo"
  if (!file.exists(meminfo)) return(NA_real_)
  lines <- readLines(meminfo, warn = FALSE)
  mem_avail <- grep("^MemAvailable:", lines, value = TRUE)
  if (!length(mem_avail)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", mem_avail[[1L]])))
  if (!is.finite(kb) || kb <= 0) return(NA_real_)
  kb * 1024
}

.warplda_memory_estimate <- function(dtm,
                                     K,
                                     n_threads = NULL,
                                     safety_factor = 1.6) {
  .assert_pkg("Matrix")
  if (!inherits(dtm, "dgCMatrix")) dtm <- methods::as(dtm, "dgCMatrix")

  K <- suppressWarnings(as.integer(K[[1L]]))
  if (!is.finite(K) || K < 2L) .log_abort("K must be an integer >= 2.")

  n_threads <- .warplda_default_threads(n_threads)
  safety_factor <- suppressWarnings(as.numeric(safety_factor[[1L]]))
  if (!is.finite(safety_factor) || safety_factor < 1) safety_factor <- 1.6

  n_doc <- as.numeric(nrow(dtm))
  n_word <- as.numeric(ncol(dtm))
  n_nz <- as.numeric(Matrix::nnzero(dtm))
  n_tokens <- as.numeric(sum(dtm))

  bytes_int <- 4
  bytes_double <- 8

  token_bytes <- n_tokens * bytes_int * 5
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
      estimated_peak_memory_gb = as.numeric(memory_estimate$estimated_peak_gb),
      estimated_peak_memory = memory_estimate$label
    )
  )
}

.existing_file <- function(path) {
  is.character(path) & !is.na(path) & nzchar(path) & file.exists(path)
}

.warplda_completed_from_cache <- function(K_grid, existing_metrics, fit_files_all, sampler = NULL) {
  sampler_value <- if (is.null(sampler)) NULL else as.character(sampler[[1L]])
  done_from_file <- .existing_file(fit_files_all)
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
      done_from_metrics[[idx]] <- any(.existing_file(metric_paths))
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
    data.table::data.table(
      K = as.integer(fit$K),
      perplexity = as.numeric(fit$metrics$perplexity),
      loglik = as.numeric(fit$metrics$loglik_approx),
      n_tokens = as.numeric(fit$metrics$n_tokens),
      sampler = fit_sampler,
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
    if (!is.na(fit_file)) {
      saveRDS(fit, fit_file)
    }
    out <- data.table::data.table(
      K = as.integer(fit$K),
      perplexity = as.numeric(fit$metrics$perplexity),
      loglik = as.numeric(fit$metrics$loglik_approx),
      n_tokens = as.numeric(fit$metrics$n_tokens),
      sampler = sampler,
      fit_file = fit_file
    )
    rm(fit)
    invisible(gc())
    out
  }

  fit_files_all <- vapply(K_grid, fit_file_for_k, character(1))
  cache_status <- .warplda_completed_from_cache(K_grid, existing_metrics, fit_files_all, sampler = sampler)
  done_from_file <- cache_status$done_from_file
  done_from_metrics <- cache_status$done_from_metrics
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
      if (concurrent_peak_bytes > 0.8 * available_bytes) {
        .log_warn(
          "WarpLDA estimated concurrent peak uses >80% of currently available memory; reduce workers or threads_per_model if this is not intentional."
        )
      }
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
        metrics_list <- c(metrics_list, list(keep[, .(K, perplexity, loglik, n_tokens, sampler, fit_file)]))
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
  list(fit_files = fit_files, metrics = metrics_tbl[, .(K, perplexity, loglik, n_tokens, sampler)])
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

score_terms_normtop <- function(phi) {
  # cisTopic NormTop (practical proxy): normalize each topic by its max term prob -> [0,1]
  phi <- as.matrix(phi)
  phi[!is.finite(phi) | phi < 0] <- 0
  mx <- apply(phi, 1, max, na.rm = TRUE)
  mx[!is.finite(mx) | mx <= 0] <- 1
  sc <- phi / mx
  sc[!is.finite(sc)] <- 0
  sc[sc < 0] <- 0
  sc[sc > 1] <- 1
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
  tf_list <- unique(as.character(tf_list))
  tf_list <- tf_list[!is.na(tf_list) & nzchar(tf_list)]
  tf_upper <- toupper(tf_list)
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
  cats <- c("PEAK", "Target", "TF")
  topic_labels <- paste0("Topic", seq_len(K))
  topic_labels_plot <- rev(topic_labels)
  col_map <- c(PEAK = "#4C78A8", Target = "#54A24B", TF = "#F58518")

  .empty_mat <- function() {
    m <- matrix(0L, nrow = length(cats), ncol = K, dimnames = list(cats, topic_labels))
    m
  }

  .count_terms <- function(tt) {
    m <- .empty_mat()
    if (is.null(tt) || !is.data.frame(tt) || !nrow(tt)) return(m)
    dt <- data.table::as.data.table(tt)
    if (!"in_topic" %in% names(dt) || !"term_id" %in% names(dt)) return(m)
    if (!"topic_num" %in% names(dt)) {
      if ("topic" %in% names(dt)) {
        dt[, topic_num := as.integer(gsub("^Topic", "", as.character(topic)))]
      } else {
        return(m)
      }
    }
    dt <- dt[.as_logical_flag(in_topic)]
    dt <- dt[is.finite(topic_num) & topic_num >= 1L & topic_num <= K]
    if (!nrow(dt)) return(m)
    dt[, term_str := as.character(term_id)]
    dt[, gene_name := sub("^GENE:", "", term_str)]
    dt[, category := data.table::fifelse(
      grepl("^PEAK:", term_str), "PEAK",
      data.table::fifelse(grepl("^GENE:", term_str) & (gene_name %in% tf_list | toupper(gene_name) %in% tf_upper), "TF",
              data.table::fifelse(grepl("^GENE:", term_str), "Target", NA_character_))
    )]
    dt <- dt[!is.na(category), .(topic_num, term_str, category)]
    if (!nrow(dt)) return(m)
    dt <- unique(dt, by = c("topic_num", "term_str", "category"))
    cnt <- dt[, .N, by = .(topic_num, category)]
    for (i in seq_len(nrow(cnt))) {
      m[cnt$category[i], cnt$topic_num[i]] <- cnt$N[i]
    }
    m
  }

  .count_links <- function(tl, use_and = FALSE) {
    m <- .empty_mat()
    if (is.null(tl) || !is.data.frame(tl) || !nrow(tl)) return(m)
    dt <- data.table::as.data.table(tl)
    req <- c("topic_num", "peak_pass", "gene_pass", "peak_id", "gene_key")
    if (!all(req %in% names(dt))) return(m)
    dt <- dt[is.finite(topic_num) & topic_num >= 1L & topic_num <= K]
    if (!nrow(dt)) return(m)
    if (isTRUE(use_and)) {
      dt <- dt[.as_logical_flag(peak_pass) & .as_logical_flag(gene_pass)]
    } else {
      dt <- dt[.as_logical_flag(peak_pass) | .as_logical_flag(gene_pass)]
    }
    if (!nrow(dt)) return(m)
    peak_cnt <- unique(dt[!is.na(peak_id) & nzchar(peak_id), .(topic_num, id = as.character(peak_id))])[, .N, by = topic_num]
    gene_cnt <- unique(dt[!is.na(gene_key) & nzchar(gene_key), .(topic_num, id = as.character(gene_key))])[, .N, by = topic_num]
    tf_cnt <- data.table::data.table()
    if ("tf" %in% names(dt)) {
      tf_cnt <- unique(dt[!is.na(tf) & nzchar(tf), .(topic_num, id = as.character(tf))])[, .N, by = topic_num]
    }
    if (nrow(peak_cnt)) m["PEAK", peak_cnt$topic_num] <- peak_cnt$N
    if (nrow(gene_cnt)) m["Target", gene_cnt$topic_num] <- gene_cnt$N
    if (nrow(tf_cnt)) m["TF", tf_cnt$topic_num] <- tf_cnt$N
    m
  }

  term_mat <- .count_terms(topic_terms)
  link_or_mat <- .count_links(topic_links, use_and = FALSE)
  link_and_mat <- .count_links(topic_links, use_and = TRUE)
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

  global_max <- max(
    c(colSums(term_mat), colSums(link_or_mat), colSums(link_and_mat)),
    na.rm = TRUE
  )
  if (!is.finite(global_max) || global_max <= 0) global_max <- 1
  x_limit <- c(0, global_max * 1.05)

  .plot_stack <- function(mat, main_title, add_legend = FALSE) {
    mat_plot <- mat[, rev(seq_len(ncol(mat))), drop = FALSE]
    graphics::barplot(
      mat_plot,
      horiz = TRUE,
      beside = FALSE,
      col = col_map[rownames(mat_plot)],
      border = NA,
      main = main_title,
      xlab = "Count",
      ylab = "",
      names.arg = topic_labels_plot,
      las = 1,
      xlim = x_limit,
      cex.main = 1.1,
      font.main = 2,
      font.lab = 2,
      font.axis = 2
    )
    if (isTRUE(add_legend)) {
      graphics::legend(
        "topright",
        legend = names(col_map),
        fill = unname(col_map),
        border = NA,
        bty = "n",
        cex = 0.9
      )
    }
  }

  par_opts <- graphics::par(no.readonly = TRUE)
  grDevices::pdf(out_file, width = 15, height = 15)
  on.exit({
    grDevices::dev.off()
    suppressWarnings(graphics::par(par_opts))
  }, add = TRUE)

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

  graphics::par(mfrow = c(3, 1), mar = c(4.5, 8, 4, 1), oma = c(0, 0, 2.5, 0))
  .plot_stack(term_mat, "Terms In Topic", add_legend = TRUE)
  .plot_stack(link_or_mat, "Links In Topic: Peak or Gene Pass")
  .plot_stack(link_and_mat, "Links In Topic: Peak and Gene Pass")
  if (!is.null(title_prefix) && nzchar(title_prefix)) {
    graphics::mtext(title_prefix, outer = TRUE, cex = 1.1, font = 2, line = 1)
  }

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
# 9) Plot helpers (heatmaps + intertopic distance + LDAvis)
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
                                          gammafit_scope = c("topic_term_group", "global_term_group"),
                                          modules = list(
                                            pathway = TRUE,
                                            doc_topic_heatmaps = TRUE,
                                            topic_by_comparison = TRUE,
                                            topic_marker_heatmap = TRUE,
                                            intertopic_distance = TRUE,
                                            ldavis = TRUE
                                          ),
                                          overwrite = list(
                                            link_topic = TRUE,
                                            pathway = TRUE
                                          ),
                                          extraction_steps = NULL) {
  gammafit_scope <- match.arg(gammafit_scope)
  extraction_steps <- .normalize_topic_extraction_steps(extraction_steps)
  list(
    extraction_steps = extraction_steps,
    pathway_source = "link_scores",
    thrP = thrP,
    gammafit_scope = gammafit_scope,
    pathway_make_heatmap = FALSE,
    pathway_make_dotplot = TRUE,
    pathway_overwrite = isTRUE(overwrite$pathway),
    pathway_per_comparison = TRUE,
    pathway_per_comparison_dir = "per_comparison_pathway",
    pathway_per_comparison_flat = TRUE,
    pathway_split_direction = TRUE,
    run_pathway_gsea = FALSE,
    run_link_topic_scores = .topic_step_enabled(extraction_steps, "topic_links", TRUE),
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
    dot_top_n_per_topic = Inf,
    max_pathways = Inf,
    run_pathway_enrichment = .topic_step_enabled(extraction_steps, "pathway", isTRUE(modules$pathway)),
    run_doc_topic_heatmaps = .topic_step_enabled(extraction_steps, "doc_topic_heatmaps", isTRUE(modules$doc_topic_heatmaps)) ||
      .topic_step_enabled(extraction_steps, "tf_topic_heatmaps", isTRUE(modules$doc_topic_heatmaps)),
    run_tf_topic_heatmaps = .topic_step_enabled(extraction_steps, "tf_topic_heatmaps", isTRUE(modules$doc_topic_heatmaps)),
    run_topic_by_comparison_heatmaps = .topic_step_enabled(extraction_steps, "topic_by_comparison", isTRUE(modules$topic_by_comparison)),
    run_topic_marker_heatmap = .topic_step_enabled(extraction_steps, "topic_marker_heatmap", isTRUE(modules$topic_marker_heatmap)),
    run_intertopic_distance_map = .topic_step_enabled(extraction_steps, "intertopic_distance", isTRUE(modules$intertopic_distance)),
    run_ldavis = .topic_step_enabled(extraction_steps, "ldavis", isTRUE(modules$ldavis))
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

.prepare_tf_doc_topic_assignment_heatmap_data <- function(theta,
                                                          phi,
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

  score_mat <- score_terms_normtop(phi)
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

.ldavis_mds_safe <- function(phi, theta) {
  phi <- as.matrix(phi)
  if (nrow(phi) < 2L) {
    return(data.frame(x = 0, y = 0))
  }
  phi[!is.finite(phi)] <- 0
  rs <- rowSums(phi)
  rs[!is.finite(rs) | rs == 0] <- 1
  phi <- phi / rs
  eps <- 1e-12
  phi <- pmax(phi, eps)
  phi <- phi / rowSums(phi)

  kl_div <- function(p, q) {
    idx <- (p > 0) & (q > 0)
    sum(p[idx] * log(p[idx] / q[idx]))
  }

  K <- nrow(phi)
  dmat <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i >= j) next
      p <- phi[i, ]
      q <- phi[j, ]
      m <- 0.5 * (p + q)
      jsd <- 0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
      if (!is.finite(jsd)) jsd <- 0
      dmat[i, j] <- jsd
      dmat[j, i] <- jsd
    }
  }
  diag(dmat) <- 0
  if (max(dmat, na.rm = TRUE) == 0) {
    coords <- cbind(seq_len(K), rep(0, K))
  } else {
    coords <- stats::cmdscale(stats::as.dist(dmat), k = 2)
  }
  if (is.null(dim(coords))) coords <- cbind(coords, rep(0, length(coords)))
  data.frame(x = coords[, 1], y = coords[, 2])
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
    grDevices::pdf(out_file, width = width, height = height, family = "Helvetica")
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
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

plot_topic_marker_heatmap <- function(phi,
                                      out_file,
                                      top_n = 20L,
                                      title_prefix = NULL,
                                      edges_docs = NULL,
                                      option_label = NULL,
                                      topic_terms = NULL) {
  .assert_pkg("data.table")
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    .log_inform("Skipping topic marker heatmap: {.pkg pheatmap} not installed.")
    return(invisible(NULL))
  }
  phi <- as.matrix(phi)
  if (!nrow(phi) || !ncol(phi)) return(invisible(NULL))
  if (is.null(colnames(phi))) {
    .log_inform("Skipping topic marker heatmap: phi has no column names.")
    return(invisible(NULL))
  }
  if (is.null(rownames(phi))) {
    rownames(phi) <- paste0("Topic", seq_len(nrow(phi)))
  }
  K <- nrow(phi)
  if (K < 2L) {
    .log_inform("Skipping topic marker heatmap: need >= 2 topics.")
    return(invisible(NULL))
  }

  phi[!is.finite(phi)] <- 0
  rs <- rowSums(phi)
  rs[!is.finite(rs) | rs == 0] <- 1
  phi <- phi / rs

  sum_all <- colSums(phi)
  denom <- pmax(K - 1L, 1L)
  mean_other <- (matrix(sum_all, nrow = K, ncol = ncol(phi), byrow = TRUE) - phi) / denom
  eps <- 1e-12
  log2fc <- log2(phi + eps) - log2(mean_other + eps)

  term_names <- colnames(phi)
  term_index <- stats::setNames(seq_along(term_names), term_names)
  top_idx <- lapply(seq_len(K), function(k) {
    vals <- log2fc[k, ]
    o <- order(vals, decreasing = TRUE, na.last = NA)
    head(o, min(as.integer(top_n), length(o)))
  })

  rows <- character()
  row_owner <- integer()
  for (k in seq_len(K)) {
    terms <- term_names[top_idx[[k]]]
    if (length(terms)) {
      terms <- terms[!terms %in% rows]
      if (length(terms)) {
        rows <- c(rows, terms)
        row_owner <- c(row_owner, rep.int(k, length(terms)))
      }
    }
  }
  rows <- rows[!is.na(rows) & nzchar(rows)]
  if (!length(rows)) return(invisible(NULL))
  row_owner <- row_owner[seq_along(rows)]

  mat <- t(phi[, rows, drop = FALSE])
  colnames(mat) <- rownames(phi)
  # Use per-term row z-score (not global scaling) to avoid near-uniform dark maps.
  mat_plot <- t(scale(t(mat)))
  mat_plot[!is.finite(mat_plot)] <- 0
  z_cap <- 2.5
  mat_plot <- pmax(pmin(mat_plot, z_cap), -z_cap)
  main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "Topic markers (see topic_marker_features.csv)", sep = " | ")
  } else {
    "Topic markers (see topic_marker_features.csv)"
  }

  width <- max(7, ncol(mat) * 0.25)
  height <- max(6, min(100, nrow(mat) * 0.02))
  marker_colors <- grDevices::colorRampPalette(
    c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b")
  )(100)
  breaks <- seq(-z_cap, z_cap, length.out = 101)

  feature_path <- file.path(dirname(out_file), "topic_marker_features.csv")
  marker_dt <- data.table::data.table(
    owner_topic = as.integer(row_owner),
    term_id = rows,
    marker_row_order = seq_along(rows)
  )
  marker_dt <- unique(marker_dt, by = c("owner_topic", "term_id"))

  feature_tbl <- NULL
  if (!is.null(topic_terms) && nrow(topic_terms)) {
    tt <- data.table::as.data.table(topic_terms)
    if (all(c("topic", "term_id") %in% names(tt))) {
      if ("in_topic" %in% names(tt)) {
        tt <- tt[.as_logical_flag(in_topic)]
      }
      tt <- tt[!is.na(topic) & !is.na(term_id) & nzchar(term_id)]
      tt[, topic := as.integer(topic)]
      tt <- tt[is.finite(topic)]
      if (nrow(tt)) {
        if ("score" %in% names(tt)) {
          tt[, score_num := .safe_num(score)]
          data.table::setorder(tt, topic, -score_num, term_id)
        } else {
          data.table::setorder(tt, topic, term_id)
        }
        tt <- unique(tt[, .(owner_topic = topic, term_id)], by = c("owner_topic", "term_id"))
        feature_tbl <- tt
      }
    }
  }
  if (is.null(feature_tbl) || !nrow(feature_tbl)) {
    feature_tbl <- marker_dt[, .(owner_topic, term_id)]
  }

  owner_idx <- term_index[feature_tbl$term_id]
  owner_idx[!is.finite(owner_idx)] <- NA_integer_
  owner_log2fc <- rep(NA_real_, nrow(feature_tbl))
  owner_prob <- rep(NA_real_, nrow(feature_tbl))
  keep_eval <- !is.na(owner_idx) &
    !is.na(feature_tbl$owner_topic) &
    feature_tbl$owner_topic >= 1L & feature_tbl$owner_topic <= nrow(log2fc)
  if (any(keep_eval)) {
    owner_log2fc[keep_eval] <- log2fc[cbind(feature_tbl$owner_topic[keep_eval], owner_idx[keep_eval])]
    owner_prob[keep_eval] <- phi[cbind(feature_tbl$owner_topic[keep_eval], owner_idx[keep_eval])]
  }
  feature_tbl[, term_base := gsub("^(GENE:|PEAK:)", "", term_id)]
  feature_tbl[, term_prefix := ifelse(grepl("^GENE:", term_id), "GENE",
                                      ifelse(grepl("^PEAK:", term_id), "PEAK", NA_character_))]
  feature_tbl[, owner_log2fc := as.numeric(owner_log2fc)]
  feature_tbl[, owner_prob := as.numeric(owner_prob)]
  feature_tbl <- merge(
    feature_tbl,
    marker_dt,
    by = c("owner_topic", "term_id"),
    all.x = TRUE
  )
  feature_tbl[, in_marker_heatmap := is.finite(marker_row_order)]
  data.table::setorder(feature_tbl, owner_topic, -in_marker_heatmap, marker_row_order, term_id)
  feature_tbl[, row_order := seq_len(.N)]
  data.table::setcolorder(
    feature_tbl,
    c("row_order", "term_id", "term_base", "term_prefix", "owner_topic",
      "owner_log2fc", "owner_prob", "in_marker_heatmap", "marker_row_order")
  )
  if (!is.null(edges_docs) && !is.null(option_label)) {
    ed <- data.table::as.data.table(edges_docs)
    ed <- ed[!is.na(tf) & nzchar(tf)]
    if (option_label == "opt3_gene_fc_expr") {
      ed <- ed[!is.na(gene_key) & nzchar(gene_key) & !is.na(peak_id) & nzchar(peak_id)]
      map_dt <- ed[, .(
        tf_list = paste(sort(unique(tf)), collapse = ";"),
        peak_list = paste(sort(unique(peak_id)), collapse = ";")
      ), by = gene_key]
      feature_tbl <- merge(
        feature_tbl,
        map_dt,
        by.x = "term_base",
        by.y = "gene_key",
        all.x = TRUE
      )
    } else if (option_label == "joint") {
      ed <- ed[!is.na(gene_key) & nzchar(gene_key) & !is.na(peak_id) & nzchar(peak_id)]
      map_gene <- ed[, .(
        tf_list = paste(sort(unique(tf)), collapse = ";"),
        peak_list = paste(sort(unique(peak_id)), collapse = ";")
      ), by = gene_key]
      map_peak <- ed[, .(
        tf_list = paste(sort(unique(tf)), collapse = ";"),
        gene_list = paste(sort(unique(gene_key)), collapse = ";")
      ), by = peak_id]
      gene_tbl <- merge(
        feature_tbl[term_prefix == "GENE"],
        map_gene,
        by.x = "term_base",
        by.y = "gene_key",
        all.x = TRUE
      )
      peak_tbl <- merge(
        feature_tbl[term_prefix == "PEAK"],
        map_peak,
        by.x = "term_base",
        by.y = "peak_id",
        all.x = TRUE
      )
      other_tbl <- feature_tbl[is.na(term_prefix) | !term_prefix %in% c("GENE", "PEAK")]
      feature_tbl <- data.table::rbindlist(list(gene_tbl, peak_tbl, other_tbl), use.names = TRUE, fill = TRUE)
      data.table::setorder(feature_tbl, row_order)
    } else {
      ed <- ed[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key)]
      map_dt <- ed[, .(
        tf_list = paste(sort(unique(tf)), collapse = ";"),
        gene_list = paste(sort(unique(gene_key)), collapse = ";")
      ), by = peak_id]
      feature_tbl <- merge(
        feature_tbl,
        map_dt,
        by.x = "term_base",
        by.y = "peak_id",
        all.x = TRUE
      )
    }
  }
  data.table::fwrite(feature_tbl, feature_path)

  grDevices::pdf(out_file, width = width, height = height)
  tryCatch(
    {
      pheatmap::pheatmap(
        mat_plot,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        color = marker_colors,
        breaks = breaks,
        border_color = NA,
        main = main_title
      )
    },
    finally = grDevices::dev.off()
  )
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
  topic_nums <- seq_len(ncol(theta_mat))

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

compute_topic_links <- function(edges_docs,
                                score_mat,
                                raw_score_mat = NULL,
                                topic_terms = NULL,
                                fp_term_mode = c("unique", "aggregate", "aggregate_weight"),
                                binarize_method = c("gammafit", "topn"),
                                gammafit_scope = c("topic_term_group", "global_term_group"),
                                link_method = c("gammafit", "link_score_efdr", "link_score_prob", "gene_prob"),
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
  } else {
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

  dt <- dt[!is.na(doc_id) & nzchar(doc_id)]
  if (identical(fp_term_mode, "aggregate")) {
    dt[, peak_term := paste0("PEAK:", gene_key)]
  } else {
    dt[, peak_term := paste0("PEAK:", peak_id)]
  }
  dt[, gene_term := paste0("GENE:", gene_key)]
  if (identical(fp_term_mode, "aggregate_weight")) {
    dt[, peak_idx := NA_integer_]
  } else {
    dt[, peak_idx := term_index[peak_term]]
  }
  dt[, gene_idx := term_index[gene_term]]
  if (identical(fp_term_mode, "aggregate_weight")) {
    dt <- dt[!is.na(gene_idx)]
  } else {
    dt <- dt[!is.na(peak_idx) & !is.na(gene_idx)]
  }
  if (!nrow(dt)) {
    .log_inform("No valid links for topic_links after term matching.")
    return(invisible(FALSE))
  }

  gamma_cutoffs_peak <- rep(NA_real_, K)
  gamma_cutoffs_gene <- rep(NA_real_, K)
  in_set_map <- NULL
  if (binarize_method == "gammafit") {
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
    has_peak_terms <- !identical(fp_term_mode, "aggregate_weight")
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
    rep_topic <- rep(seq_len(K), times = n)
    peak_score <- as.vector(peak_scores_raw)
    gene_score <- as.vector(gene_scores_raw)
    peak_score_gate <- as.vector(peak_scores_gate)
    gene_score_gate <- as.vector(gene_scores_gate)
    peaks_gamma_cutoff <- rep(gamma_cutoffs_peak, times = n)
    gene_gamma_cutoff <- rep(gamma_cutoffs_gene, times = n)
    if (binarize_method == "gammafit") {
      if (has_peak_terms) {
        peak_pass <- is.finite(peaks_gamma_cutoff) & peak_score_gate >= peaks_gamma_cutoff & peak_score_gate > 0
      } else {
        peak_pass <- rep(TRUE, length(rep_topic))
      }
      gene_pass <- is.finite(gene_gamma_cutoff) & gene_score_gate >= gene_gamma_cutoff & gene_score_gate > 0
    } else if (!is.null(in_set_map)) {
      peak_allowed <- in_set_map[chunk_dt$peak_term]
      gene_allowed <- in_set_map[chunk_dt$gene_term]
      peak_pass <- rep(identical(fp_term_mode, "aggregate_weight"), length(rep_topic))
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
    n_before_gammafit <- nrow(dt) * K
    n_pass_gammafit <- nrow(out)
    gammafit_label <- if (identical(fp_term_mode, "aggregate_weight")) "gene_only" else "peak_and_gene"
    .log_inform(
      "topic_links gammafit: {gammafit_label} pass {n_pass_gammafit}/{n_before_gammafit} rows."
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

  pass_dt <- if (identical(link_method, "gammafit")) {
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
    summary_file <- file.path(dirname(if (!is.null(pass_file)) pass_file else out_file), "topic_link_summary.csv")
    n_scored_rows <- if (identical(link_method, "gammafit")) {
      as.double(nrow(dt) * K)
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
  }
  if (identical(link_method, "gammafit")) {
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
    WikiPathways_2024_Human = "WikiPathways",
    MSigDB_Hallmark_2020 = "Hallmark",
    KEGG_2021_Human = "KEGG"
  )
  if (db %in% names(db_short)) db_short[[db]] else db
}

.topic_enrichr_result_to_table <- function(enr, topic_name, genes) {
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
      padj = as.numeric(df$Adjusted.P.value),
      pval = if ("P.value" %in% names(df)) as.numeric(df$P.value) else NA_real_,
      overlap = if ("Overlap" %in% names(df)) as.character(df$Overlap) else NA_character_,
      overlap_hits = overlap_hits,
      genes = if ("Genes" %in% names(df)) as.character(df$Genes) else NA_character_,
      logp = -log10(df$Adjusted.P.value),
      combined_score = as.numeric(combined_score),
      odds_ratio = as.numeric(odds_ratio),
      cluster_size = length(genes)
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
                                                  dbs = c(
                                                    "GO_Biological_Process_2023",
                                                    "GO_Cellular_Component_2023",
                                                    "GO_Molecular_Function_2023",
                                                    "Reactome_2022",
                                                    "WikiPathways_2024_Human",
                                                    "MSigDB_Hallmark_2020",
                                                    "KEGG_2021_Human"
                                                  ),
                                                  padj_cut = 0.05,
                                                  min_genes = 5L,
                                                  top_n_per_topic = 20L,
                                                  max_pathways = 200L,
                                                  title_prefix = NULL,
                                                  use_all_terms = FALSE,
                                                  make_heatmap = TRUE,
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
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- file.path(dirname(out_file), "cache", "enrichr")
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
    msg <- "Skipping pathway enrichment heatmap: neither enrichly nor enrichR is installed."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  if (isTRUE(make_heatmap) && !requireNamespace("pheatmap", quietly = TRUE)) {
    msg <- "Skipping pathway enrichment heatmap: pheatmap not installed."
    .log_inform(msg)
    log_msg(msg)
    make_heatmap <- FALSE
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

  gene_sets <- topic_gene_sets(
    topic_terms = topic_terms,
    edges_docs = edges_docs,
    option_label = option_label,
    use_all_terms = use_all_terms,
    theta = theta,
    tf_link_mode = tf_link_mode,
    tf_theta_top_n = tf_theta_top_n,
    tf_theta_min = tf_theta_min
  )
  if (!length(gene_sets)) {
    msg <- "Skipping pathway enrichment heatmap: no topic gene sets after mapping."
    .log_inform(msg)
    log_msg(msg)
    return(invisible(NULL))
  }
  log_msg(sprintf("Option: %s", option_label))
  log_msg(sprintf("Topics with gene sets: %s", paste(names(gene_sets), collapse = ",")))
  log_msg(sprintf("DBs: %s", paste(dbs, collapse = ",")))

  res_list <- vector("list", length(gene_sets))
  names(res_list) <- names(gene_sets)
  run_one_enrichr <- function(nm) {
    genes <- gene_sets[[nm]]
    log_msg(sprintf("Topic %s gene count: %d", nm, length(genes)))
    if (length(genes) < as.integer(min_genes)) {
      log_msg(sprintf("Topic %s skipped: gene count < %d", nm, min_genes))
      return(NULL)
    }
    enr <- tryCatch(
      .quiet_enrichr_call(.run_enrichr_cached(
        genes = genes,
        dbs = dbs,
        sleep_time = enrichr_sleep_time,
        cache_dir = enrichr_cache_dir,
        backend = pathway_backend
      )),
      error = function(e) {
        log_msg(sprintf("Topic %s enrichr error: %s", nm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(enr)) {
      log_msg(sprintf("Topic %s enrichr returned NULL.", nm))
      return(NULL)
    }
    rows <- lapply(names(enr), function(db) {
      df <- enr[[db]]
      if (is.null(df) || !nrow(df)) return(NULL)
      if (!("Adjusted.P.value" %in% names(df)) || !("Term" %in% names(df))) return(NULL)
      df <- df[is.finite(df$Adjusted.P.value) & df$Adjusted.P.value <= padj_cut, , drop = FALSE]
      if (!nrow(df)) return(NULL)
      db_short <- c(
        GO_Biological_Process_2023 = "GO:BP",
        GO_Cellular_Component_2023 = "GO:CC",
        GO_Molecular_Function_2023 = "GO:MF",
        Reactome_2022 = "Reactome",
        WikiPathways_2024_Human = "WikiPathways",
        MSigDB_Hallmark_2020 = "Hallmark",
        KEGG_2021_Human = "KEGG"
      )
      db_label <- if (db %in% names(db_short)) db_short[[db]] else db
      term_clean <- gsub("\\s*\\([^)]*\\)$", "", df$Term)
      combined_score <- if ("Combined.Score" %in% names(df)) df$Combined.Score else NA_real_
      odds_ratio <- if ("Odds.Ratio" %in% names(df)) df$Odds.Ratio else NA_real_
      data.table::data.table(
        topic = as.integer(nm),
        pathway = paste(db_label, term_clean, sep = ": "),
        padj = as.numeric(df$Adjusted.P.value),
        pval = if ("P.value" %in% names(df)) as.numeric(df$P.value) else NA_real_,
        overlap = if ("Overlap" %in% names(df)) as.character(df$Overlap) else NA_character_,
        genes = if ("Genes" %in% names(df)) as.character(df$Genes) else NA_character_,
        logp = -log10(df$Adjusted.P.value),
        combined_score = as.numeric(combined_score),
        odds_ratio = as.numeric(odds_ratio)
      )
    })
    out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    n_hits <- nrow(out)
    log_msg(sprintf("Topic %s enriched pathways: %d (padj<=%s)", nm, n_hits, padj_cut))
    out
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
  log_msg(sprintf("Total enriched pathways (unique): %d", length(unique(res_dt$pathway))))
  log_msg(sprintf("Debug log written to: %s", log_path))

  main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "Pathway enrichment", sep = " | ")
  } else {
    "Pathway enrichment"
  }

  if (isTRUE(make_heatmap)) {
    mat_dt <- data.table::dcast(res_dt, pathway ~ topic, value.var = "logp", fill = 0)
    mat <- as.matrix(mat_dt[, -1, with = FALSE])
    rownames(mat) <- mat_dt$pathway
    colnames(mat) <- paste0("Topic", colnames(mat))

    n_topics <- ncol(mat)
    width <- max(10, n_topics * 0.6)
    height <- max(8, min(160, nrow(mat) * 0.25))
    font_row <- if (nrow(mat) > 200) 4 else if (nrow(mat) > 120) 5 else 7
    font_col <- if (n_topics > 60) 5 else if (n_topics > 40) 6 else if (n_topics > 25) 8 else 10

    grDevices::pdf(out_file, width = width, height = height)
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
  } else {
    log_msg("Skipping pathway heatmap: make_heatmap = FALSE.")
  }

  if (isTRUE(make_dotplot) && requireNamespace("ggplot2", quietly = TRUE)) {
    dot_path <- file.path(dirname(out_file), "topic_pathway_enrichment_dotplot.pdf")
    plot_dt <- data.table::copy(res_dt)
    plot_dt[, topic_num := as.integer(topic)]
    plot_dt[, topic := paste0("Topic", topic_num)]
    plot_dt[, score_val := ifelse(is.finite(combined_score), combined_score, ifelse(is.finite(odds_ratio), odds_ratio, logp))]
    plot_dt[, score_val := pmax(score_val, 0)]
    plot_dt[, size_val := sqrt(score_val)]
    path_order <- plot_dt[, .(
      min_topic = min(topic_num, na.rm = TRUE),
      max_score = max(score_val, na.rm = TRUE)
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
      ggplot2::scale_size_continuous(name = "Score", range = size_range) +
      ggplot2::scale_y_discrete(labels = function(x) wrap_label(x, width = 70)) +
      ggplot2::scale_x_discrete(labels = function(x) gsub("^Topic", "", x)) +
      ggplot2::labs(
        x = "Topic",
        y = NULL,
        title = if (!is.null(title_prefix) && nzchar(title_prefix)) {
          paste(title_prefix, "Pathway dot plot", sep = " | ")
        } else {
          "Pathway dot plot"
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
      ord_dt[, .(topic, pathway, logp, padj, pval, combined_score, odds_ratio, overlap, genes)],
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

plot_topic_pathway_enrichment_from_link_scores <- function(link_scores,
                                                           out_dir,
                                                           title_prefix = NULL,
                                                           file_tag = NULL,
                                                           dbs = c(
                                                             "GO_Biological_Process_2023",
                                                             "GO_Cellular_Component_2023",
                                                             "GO_Molecular_Function_2023",
                                                             "Reactome_2022",
                                                             "WikiPathways_2024_Human",
                                                             "MSigDB_Hallmark_2020",
                                                             "KEGG_2021_Human"
                                                           ),
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
                                                           make_heatmap = TRUE,
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
  if (is.null(enrichr_cache_dir)) {
    enrichr_cache_dir <- file.path(out_dir, "cache", "enrichr")
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
    msg <- "Skipping link-score pathway enrichment: neither enrichly nor enrichR is installed."
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
    run_one_enrichr <- function(nm) {
      genes <- gene_sets[[nm]]
      if (length(genes) < as.integer(min_genes)) {
        log_fun(sprintf("Topic %s skipped: gene count < %d", nm, min_genes))
        return(NULL)
      }
      log_fun(sprintf("Topic %s gene count: %d", nm, length(genes)))
      enr <- tryCatch(
        .quiet_enrichr_call(.run_enrichr_cached(
          genes = genes,
          dbs = dbs,
          sleep_time = enrichr_sleep_time,
          cache_dir = enrichr_cache_dir,
          backend = pathway_backend
        )),
        error = function(e) {
          log_fun(sprintf("Topic %s enrichr error: %s", nm, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(enr)) return(NULL)
      out <- .topic_enrichr_result_to_table(enr, topic_name = nm, genes = genes)
      log_fun(sprintf("Topic %s enriched pathways: %d finite adjusted p-values", nm, nrow(out)))
      out
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

  main_title <- if (!is.null(title_prefix) && nzchar(title_prefix)) {
    paste(title_prefix, "Pathway enrichment (link scores)", sep = " | ")
  } else {
    "Pathway enrichment (link scores)"
  }

  if (isTRUE(make_heatmap) && requireNamespace("pheatmap", quietly = TRUE)) {
    mat_dt <- data.table::dcast(res_dt, pathway ~ topic, value.var = "logp", fill = 0)
    mat <- as.matrix(mat_dt[, -1, with = FALSE])
    rownames(mat) <- mat_dt$pathway
    colnames(mat) <- paste0("Topic", colnames(mat))

    n_topics <- ncol(mat)
    width <- max(10, n_topics * 0.6)
    height <- max(8, min(160, nrow(mat) * 0.25))
    font_row <- if (nrow(mat) > 200) 4 else if (nrow(mat) > 120) 5 else 7
    font_col <- if (n_topics > 60) 5 else if (n_topics > 40) 6 else if (n_topics > 25) 8 else 10

    out_file <- file.path(out_dir, "topic_pathway_enrichment_heatmap.pdf")
    grDevices::pdf(out_file, width = width, height = height)
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
  } else if (isTRUE(make_heatmap)) {
    log_msg("Skipping heatmap: pheatmap not installed.")
  } else {
    log_msg("Skipping heatmap: make_heatmap = FALSE.")
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
    sizes <- counts$N[match(seq_len(K), counts$topic)]
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

save_ldavis_html <- function(theta, phi, dtm, out_dir) {
  if (!requireNamespace("LDAvis", quietly = TRUE)) {
    .log_inform("Skipping LDAvis HTML: {.pkg LDAvis} not installed.")
    return(invisible(NULL))
  }
  if (is.null(rownames(theta)) || is.null(colnames(phi))) return(invisible(NULL))

  theta <- as.matrix(theta)
  phi <- as.matrix(phi)
  theta[!is.finite(theta)] <- 0
  phi[!is.finite(phi)] <- 0
  eps <- 1e-12
  theta <- theta + eps
  phi <- phi + eps
  rs_theta <- rowSums(theta)
  rs_theta[!is.finite(rs_theta) | rs_theta == 0] <- 1
  theta_prob <- theta / rs_theta
  rs_phi <- rowSums(phi)
  rs_phi[!is.finite(rs_phi) | rs_phi == 0] <- 1
  phi_prob <- phi / rs_phi
  theta_prob[!is.finite(theta_prob)] <- 0
  phi_prob[!is.finite(phi_prob)] <- 0

  doc_len <- as.numeric(Matrix::rowSums(dtm))
  term_freq <- as.numeric(Matrix::colSums(dtm))
  vocab <- colnames(dtm)

  json <- suppressWarnings(LDAvis::createJSON(
    phi = phi_prob,
    theta = theta_prob,
    doc.length = doc_len,
    vocab = vocab,
    term.frequency = term_freq,
    mds.method = .ldavis_mds_safe
  ))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  src_dir <- system.file("htmljs", package = "LDAvis")
  if (dir.exists(src_dir)) {
    assets <- list.files(src_dir, full.names = TRUE)
    file.copy(assets, out_dir, overwrite = TRUE, recursive = TRUE)
  }
  json_path <- file.path(out_dir, "lda.json")
  con <- file(json_path, encoding = getOption("encoding"))
  on.exit(close(con), add = TRUE)
  cat(json, file = con)
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
                                              binarize_method = c("gammafit", "topn"),
                                              gammafit_scope = c("topic_term_group", "global_term_group"),
                                              thrP = 0.975,
                                              top_n_terms = 500L,
                                              in_topic_min_terms = 50L,
                                              pathway_use_all_terms = FALSE,
                                              pathway_make_heatmap = TRUE,
                                              pathway_make_dotplot = TRUE,
                                              pathway_overwrite = FALSE,
                                              top_n_per_topic = 20L,
                                              dot_top_n_per_topic = Inf,
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
                                              pathway_per_comparison_dir = "per_comparison_pathway",
                                              pathway_per_comparison_flat = FALSE,
                                              pathway_split_direction = TRUE,
                                              pathway_enrichr_sleep_time = 0,
                                              pathway_enrichr_cache_dir = NULL,
                                              pathway_enrichr_n_cores = NULL,
                                              pathway_backend = NULL,
                                              run_link_topic_scores = FALSE,
                                              fp_term_mode = c("unique", "aggregate", "aggregate_weight"),
                                              link_topic_gate_mode = "none",
                                              link_topic_top_k = 3L,
                                              link_topic_min_prob = 0,
                                              link_topic_include_tf = FALSE,
                                              link_topic_chunk_size = 5000L,
                                              link_topic_n_cores = 1L,
                                              link_topic_overwrite = FALSE,
                                              link_topic_method = c("gammafit", "gene_prob", "link_score_prob", "link_score_efdr"),
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
                                              run_doc_topic_heatmaps = TRUE,
                                              run_tf_topic_heatmaps = TRUE,
                                              run_topic_by_comparison_heatmaps = TRUE,
                                              run_topic_marker_heatmap = TRUE,
                                              run_intertopic_distance_map = TRUE,
                                              run_ldavis = TRUE,
                                              extraction_steps = NULL,
                                              topic_by_comparison_label_cleaner = NULL,
                                              doc_design = c("comparison", "condition"),
                                              title_prefix = NULL) {
  option_label <- match.arg(option_label)
  doc_design <- match.arg(doc_design)
  direction_by <- match.arg(direction_by)
  binarize_method <- match.arg(binarize_method)
  gammafit_scope <- match.arg(gammafit_scope)
  pathway_tf_link_mode <- match.arg(pathway_tf_link_mode)
  gsea_peak_gene_agg <- match.arg(gsea_peak_gene_agg)
  pathway_source <- match.arg(pathway_source)
  pathway_enrichr_sleep_time <- .normalize_enrichr_sleep_time(pathway_enrichr_sleep_time)
  if (is.null(pathway_enrichr_cache_dir)) {
    pathway_enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir)
  }
  pathway_enrichr_n_cores <- .normalize_enrichr_n_cores(pathway_enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
  fp_term_mode <- .resolve_fp_term_mode(fp_term_mode)
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
    run_doc_topic_heatmaps <- .topic_step_enabled(extraction_steps, "doc_topic_heatmaps", TRUE)
    run_tf_topic_heatmaps <- .topic_step_enabled(extraction_steps, "tf_topic_heatmaps", TRUE)
    run_topic_by_comparison_heatmaps <- .topic_step_enabled(extraction_steps, "topic_by_comparison", TRUE)
    run_topic_marker_heatmap <- .topic_step_enabled(extraction_steps, "topic_marker_heatmap", TRUE)
    run_intertopic_distance_map <- .topic_step_enabled(extraction_steps, "intertopic_distance", TRUE)
    run_ldavis <- .topic_step_enabled(extraction_steps, "ldavis", TRUE)
  }
  if (!all(link_topic_gate_mode %in% allowed_gate_modes)) {
    .log_abort("link_topic_gate_mode must be one of: {paste(allowed_gate_modes, collapse = ', ')}.")
  }

  .assert_pkg("cli")
  .assert_pkg("data.table")
  .assert_pkg("Matrix")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(title_prefix)) title_prefix <- basename(out_dir)

  theta <- .validate_topic_probability_matrix(topic_base$theta, "theta")
  phi <- .validate_topic_probability_matrix(topic_base$phi, "phi")

  if (is.null(rownames(theta)) && !is.null(rownames(dtm))) {
    rownames(theta) <- rownames(dtm)
  }
  if (is.null(colnames(theta))) colnames(theta) <- paste0("Topic", seq_len(ncol(theta)))
  if (is.null(rownames(phi))) rownames(phi) <- colnames(theta)
  if (is.null(colnames(phi)) && !is.null(colnames(dtm))) colnames(phi) <- colnames(dtm)

  topic_base$theta <- theta
  topic_base$phi <- phi

  score_mat <- score_terms_normtop(phi)
  topic_terms <- binarize_topics(
    score_mat,
    method = binarize_method,
    thrP = thrP,
    top_n_terms = top_n_terms,
    min_terms = in_topic_min_terms,
    gammafit_scope = gammafit_scope
  )
  data.table::fwrite(topic_terms, file.path(out_dir, "topic_terms.csv"))
  .save_all(out_dir, "topic_terms", topic_terms)
  .save_all(out_dir, "topic_term_scores_normtop", score_mat)

  topic_links_tbl <- NULL
  if (isTRUE(run_link_topic_scores)) {
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
    topic_links_tbl <- compute_topic_links(
      edges_docs = edges_docs,
      score_mat = score_mat,
      raw_score_mat = topic_base$phi,
      topic_terms = topic_terms,
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
    if (!is.data.frame(topic_links_tbl)) {
      topic_links_path <- .topic_links_path(out_dir, prefer = "pass")
      if (file.exists(topic_links_path)) {
        topic_links_tbl <- data.table::fread(topic_links_path)
      }
    }
  } else {
    topic_links_path <- .topic_links_path(out_dir, prefer = "pass")
    if (file.exists(topic_links_path)) {
      topic_links_tbl <- data.table::fread(topic_links_path)
    }
  }

  if (binarize_method == "gammafit" && isTRUE(run_gammafit_summary)) {
    tf_terms <- if (is.data.frame(edges_docs) && "tf" %in% names(edges_docs)) unique(edges_docs$tf) else NULL
    topic_terms_tbl <- data.table::fread(file.path(out_dir, "topic_terms.csv"))
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
      topic_links = topic_links_tbl
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

  if (isTRUE(run_doc_topic_heatmaps)) {
    unlink(file.path(out_dir, "doc_topic_heatmaps"), recursive = TRUE, force = TRUE)
    plot_doc_topic_heatmaps_by_comparison(
      theta = topic_base$theta,
      out_dir = file.path(out_dir, "doc_topic_heatmaps"),
      edges_docs = edges_docs,
      title_prefix = title_prefix,
      option_label = option_label,
      doc_design = doc_design
    )
  }
  if (isTRUE(run_tf_topic_heatmaps)) {
    plot_tf_doc_topic_assignment_heatmaps(
      theta = topic_base$theta,
      phi = topic_base$phi,
      out_file = file.path(out_dir, "tf_topic_assignment_heatmaps.pdf"),
      topic_terms = topic_terms,
      title_prefix = title_prefix,
      doc_design = doc_design,
      aggregate_fun = "max"
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
  }

  if (isTRUE(run_topic_marker_heatmap)) {
    plot_topic_marker_heatmap(
      phi = topic_base$phi,
      out_file = file.path(out_dir, "topic_marker_term_heatmap.pdf"),
      top_n = 20L,
      title_prefix = title_prefix,
      edges_docs = edges_docs,
      option_label = option_label,
      topic_terms = topic_terms
    )
  }

  if (isTRUE(run_pathway_enrichment) && isTRUE(pathway_overwrite)) {
    unlink(list.files(out_dir, pattern = "^topic_pathway_enrichment", full.names = TRUE), recursive = TRUE, force = TRUE)
    unlink(file.path(out_dir, pathway_per_comparison_dir), recursive = TRUE, force = TRUE)
    unlink(list.files(out_dir, pattern = paste0("^", .safe_filename(pathway_per_comparison_dir), "_"), full.names = TRUE), recursive = TRUE, force = TRUE)
  }

  if (isTRUE(run_pathway_enrichment)) {
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
        max_pathways = max_pathways,
        theta = topic_base$theta,
        tf_link_mode = pathway_tf_link_mode,
        tf_theta_top_n = pathway_tf_top_n_docs,
        tf_theta_min = pathway_tf_min_theta,
        enrichr_sleep_time = pathway_enrichr_sleep_time,
        enrichr_cache_dir = pathway_enrichr_cache_dir,
        enrichr_n_cores = pathway_enrichr_n_cores,
        pathway_backend = pathway_backend
      )

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

  if (isTRUE(run_intertopic_distance_map)) {
    plot_intertopic_distance_map(
      phi = topic_base$phi,
      topic_terms = topic_terms,
      out_file = file.path(out_dir, "intertopic_distance_map.pdf"),
      option_label = option_label,
      title_prefix = title_prefix
    )
  }
  if (isTRUE(run_ldavis)) {
    save_ldavis_html(
      theta = topic_base$theta,
      phi = topic_base$phi,
      dtm = dtm,
      out_dir = file.path(out_dir, "ldavis")
    )
  }

  list(
    out_dir = out_dir,
    edges_docs = edges_docs,
    dtm = dtm,
    topic_terms = topic_terms
  )
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
                                          binarize_method = c("gammafit", "topn"),
                                          gammafit_scope = c("topic_term_group", "global_term_group"),
                                          thrP = 0.975,
                                          top_n_terms = 500L,
                                          in_topic_min_terms = 50L,
                                          pathway_use_all_terms = FALSE,
                                          pathway_make_heatmap = TRUE,
                                          pathway_make_dotplot = TRUE,
                                          top_n_per_topic = 20L,
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
                                          pathway_per_comparison_dir = "per_comparison_pathway",
                                          pathway_per_comparison_flat = FALSE,
                                          pathway_split_direction = TRUE,
                                          pathway_overwrite = FALSE,
                                          pathway_enrichr_sleep_time = 0,
                                          pathway_enrichr_cache_dir = NULL,
                                          pathway_enrichr_n_cores = NULL,
                                          pathway_backend = NULL,
                                          run_link_topic_scores = FALSE,
                                          fp_term_mode = c("unique", "aggregate", "aggregate_weight"),
                                          link_topic_gate_mode = "none",
                                          link_topic_top_k = 3L,
                                          link_topic_min_prob = 0,
                                          link_topic_include_tf = FALSE,
                                          link_topic_chunk_size = 5000L,
                                          link_topic_n_cores = 1L,
                                          link_topic_overwrite = FALSE,
                                          link_topic_method = c("gammafit", "gene_prob", "link_score_prob", "link_score_efdr"),
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
  binarize_method <- match.arg(binarize_method)
  gammafit_scope <- match.arg(gammafit_scope)
  gsea_peak_gene_agg <- match.arg(gsea_peak_gene_agg)
  pathway_tf_link_mode <- match.arg(pathway_tf_link_mode)
  pathway_source <- match.arg(pathway_source)
  pathway_enrichr_sleep_time <- .normalize_enrichr_sleep_time(pathway_enrichr_sleep_time)
  if (is.null(pathway_enrichr_cache_dir)) {
    pathway_enrichr_cache_dir <- .module3_default_enrichr_cache_dir(out_dir)
  }
  pathway_enrichr_n_cores <- .normalize_enrichr_n_cores(pathway_enrichr_n_cores)
  pathway_backend <- .pathway_backend(pathway_backend)
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

  score_mat <- score_terms_normtop(phi)
  topic_terms <- binarize_topics(
    score_mat,
    method = binarize_method,
    thrP = thrP,
    top_n_terms = top_n_terms,
    min_terms = in_topic_min_terms,
    gammafit_scope = gammafit_scope
  )
  data.table::fwrite(topic_terms, file.path(out_dir, "topic_terms.csv"))
  .save_all(out_dir, "topic_terms", topic_terms)

  topic_links_tbl <- NULL
  if (isTRUE(run_link_topic_scores)) {
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
    topic_links_tbl <- compute_topic_links(
      edges_docs = edges_docs,
      score_mat = score_mat,
      raw_score_mat = phi,
      topic_terms = topic_terms,
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
  } else {
    topic_links_path <- file.path(out_dir, "topic_links.csv")
    if (file.exists(topic_links_path)) {
      topic_links_tbl <- data.table::fread(topic_links_path)
    }
  }

  if (binarize_method == "gammafit") {
    tf_terms <- if (is.data.frame(edges_docs) && "tf" %in% names(edges_docs)) unique(edges_docs$tf) else NULL
    topic_terms_tbl <- data.table::fread(file.path(out_dir, "topic_terms.csv"))
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
      topic_links = topic_links_tbl
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

  plot_doc_topic_heatmaps_by_comparison(
    theta = topic_base$theta,
    out_dir = file.path(out_dir, "doc_topic_heatmaps"),
    edges_docs = edges_docs,
    title_prefix = title_prefix,
    option_label = option_label
  )
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
  plot_topic_marker_heatmap(
    phi = topic_base$phi,
    out_file = file.path(out_dir, "topic_marker_term_heatmap.pdf"),
    top_n = 20L,
    title_prefix = title_prefix,
    edges_docs = edges_docs,
    option_label = option_label,
    topic_terms = topic_terms
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
      max_pathways = max_pathways,
      theta = topic_base$theta,
      tf_link_mode = pathway_tf_link_mode,
      tf_theta_top_n = pathway_tf_top_n_docs,
      tf_theta_min = pathway_tf_min_theta,
      enrichr_sleep_time = pathway_enrichr_sleep_time,
      enrichr_cache_dir = pathway_enrichr_cache_dir,
      enrichr_n_cores = pathway_enrichr_n_cores,
      pathway_backend = pathway_backend
    )

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
  save_ldavis_html(
    theta = topic_base$theta,
    phi = topic_base$phi,
    dtm = dtm,
    out_dir = file.path(out_dir, "ldavis")
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
# 12) VAE helpers (nutrient pipeline)
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
  k_grid <- unique(as.integer(k_grid))
  k_grid <- k_grid[is.finite(k_grid)]
  k_grid_txt <- paste(k_grid, collapse = ",")
  vae_python_variant <- .normalize_vae_python_variant(vae_variant)

  metrics_path <- file.path(out_dir, "model_metrics.csv")
  reuse_ok <- isTRUE(reuse_if_exists) &&
    file.exists(metrics_path) &&
    dir.exists(file.path(out_dir, "vae_models"))

  if (!reuse_ok) {
    doc_term_input <- if (isTRUE(cache_plan$write_arrow)) cache_paths$arrow else cache_paths$csv
    if (!is.character(doc_term_input) || !length(doc_term_input) || is.na(doc_term_input) || !file.exists(doc_term_input)) {
      .log_abort("Missing doc-term cache for VAE: {doc_term_input}")
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
  } else {
    .log_inform("Reusing existing VAE outputs; skipping training for {out_dir}")
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
    pathway_source = "link_scores",
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
    pathway_per_comparison_dir = "per_comparison_pathway",
    pathway_per_comparison_flat = TRUE,
    pathway_split_direction = TRUE,
    run_link_topic_scores = TRUE,
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
#' @param reuse_if_exists Reuse existing model outputs when all requested K
#'   values are present.
#' @param save_full_doc_term_csv Whether to save the full document-term CSV.
#' @param flat_output Whether to write one setup directly under `output_dir`
#'   instead of adding a long analysis-specific subdirectory. This is intended
#'   for standard package runs; benchmarks keep the nested layout.
#' @param local_threads Optional local thread count for data.table, BLAS, and native WarpLDA. NULL uses all available cores. Set options(craftgrn.warplda.max_threads = n) or CRAFTGRN_WARPLDA_MAX_THREADS to cap automatic thread use.
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
                               gene_term_mode = c("aggregate", "unique"),
                               fp_term_mode = c("aggregate", "unique", "aggregate_weight"),
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
                               reuse_if_exists = TRUE,
                               save_full_doc_term_csv = FALSE,
                               flat_output = FALSE,
                               local_threads = NULL,
                               check_repeated_values = FALSE,
                               binarize_method = "gammafit",
                               thrP = 0.9,
                               top_n_terms = 500L,
                               in_topic_min_terms = 1,
                               topic_report_args = list()) {
  .assert_pkg("data.table")
  .assert_pkg("cli")
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
  doc_mode <- match.arg(doc_mode)
  doc_design <- match.arg(doc_design)
  doc_tag <- if (identical(doc_mode, "tf")) "tf" else "ctf"
  weight_label <- if (identical(doc_design, "condition")) "peak_score_gene_expr" else "peak_log2fc_fp_gene_fc_expr"
  count_input_requested <- if (is.null(count_input) || !length(count_input)) NA_character_ else as.character(count_input[[1L]])
  count_input_effective <- .resolve_topic_count_input(count_method = count_method, count_input = count_input)
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
  .log_inform("Loaded {nrow(edges_dt)} delta-link row(s).")
  if (!("comparison_id" %in% names(edges_dt))) .log_abort("edges_all missing comparison_id.")

  sample_subset <- if (is.null(sample_subset)) NULL else unique(as.character(sample_subset))
  sample_subset <- sample_subset[!is.na(sample_subset) & nzchar(sample_subset)]
  if (length(sample_subset)) {
    if (all(c("cond1_id", "cond2_id") %in% names(edges_dt))) {
      n_before_subset <- nrow(edges_dt)
      edges_dt <- edges_dt[cond1_id %in% sample_subset & cond2_id %in% sample_subset]
      .log_inform("Sample subset retained {nrow(edges_dt)}/{n_before_subset} delta-link row(s).")
    } else {
      .log_abort("sample_subset requires delta links with cond1_id and cond2_id columns.")
    }
  }
  if (!nrow(edges_dt)) {
    .log_abort("No delta-link rows remain after applying sample_subset.")
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
    .log_inform("{analysis_id}: retained {nrow(edges_filt)} edge row(s) after topic filters in {round(proc.time()[['elapsed']] - t_filter, 1)} sec.")

    if (identical(doc_design, "condition")) {
      .log_inform("{analysis_id}: building condition-level TF documents.")
      t_docs <- proc.time()[["elapsed"]]
      edges_docs <- add_condition_tf_docs(
        edges_filt,
        tf_cluster_map = tf_cluster_map,
        doc_mode = doc_mode
      )
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
    .log_inform("{analysis_id}: doc_term has {nrow(doc_term)} row(s), {data.table::uniqueN(doc_term$doc_id)} document(s), and {data.table::uniqueN(doc_term$term_id)} term(s).")

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
    topic_input_summary <- data.table::data.table(
      analysis_label = analysis_id,
      doc_design = doc_design,
      doc_mode = doc_mode,
      fp_term_mode = fp_term_mode,
      count_method = count_method,
      count_scale = as.numeric(count_scale),
      count_input_requested = count_input_requested,
      count_input_effective = count_input_effective,
      n_link_rows_after_filter = as.double(nrow(edges_filt)),
      n_document_edge_rows = as.double(nrow(edges_docs)),
      n_doc_term_rows = as.double(nrow(doc_term)),
      n_documents = as.double(data.table::uniqueN(doc_term$doc_id)),
      n_terms = as.double(data.table::uniqueN(doc_term$term_id)),
      n_nonzero = as.double(sum(.safe_num(doc_term[[count_input_effective]]) > 0, na.rm = TRUE)),
      n_model_tokens = as.double(sum(.safe_num(doc_term[[count_input_effective]]), na.rm = TRUE))
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
    if (file.exists(metrics_path)) {
      old_metrics <- tryCatch(data.table::fread(metrics_path), error = function(e) data.table::data.table())
      if (nrow(old_metrics) && "K" %in% names(old_metrics)) {
        metrics_have_required_k <- all(required_k %in% as.integer(old_metrics$K))
        old_required <- old_metrics[as.integer(old_metrics$K) %in% required_k]
        count_metadata_ok <- all(c("count_method", "count_input_effective") %in% names(old_required)) &&
          all(as.character(old_required$count_method) == count_method) &&
          all(as.character(old_required$count_input_effective) == count_input_effective)
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
        seed = 123,
        save_tmp_dir = file.path(out_dir, "tmp_models"),
        workers = 1L,
        threads_per_model = local_threads,
        sampler = warplda_sampler,
        metrics_file = file.path(out_dir, "model_metrics.csv")
      )
      .log_inform("{analysis_id}: finished WarpLDA model fits in {round(proc.time()[['elapsed']] - t_model, 1)} sec.")
      metrics_tbl <- data.table::as.data.table(fits_out$metrics)
      metrics_tbl[, `:=`(
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
#'   Standard extraction keeps pathway enrichment and per-comparison pathway
#'   reports disabled by default for speed; set `run_pathway_enrichment = TRUE`
#'   and `pathway_per_comparison = TRUE` inside this list for the full pathway
#'   review.
#' @noRd
extract_regulatory_topics <- function(k,
                                      model_dir,
                                      output_dir,
                                      topic_input_dir = NULL,
                                      backend = c("warplda", "vae"),
                                      vae_variant = "multivi_encoder",
                                      doc_mode = c("tf_cluster", "tf"),
                                      weight_label = c("peak_log2fc_fp_gene_fc_expr", "peak_delta_fp_gene_fc_expr", "peak_score_gene_expr"),
                                      flatten_single_output = FALSE,
                                      topic_report_args = list()) {
  .assert_pkg("data.table")
  .assert_pkg("readr")

  k <- as.integer(k)
  if (!is.finite(k) || k <= 0L) .log_abort("`k` must be a positive integer.")
  backend <- match.arg(backend)
  doc_mode <- match.arg(doc_mode)
  weight_label <- match.arg(weight_label)
  report_doc_design <- if (identical(weight_label, "peak_score_gene_expr")) "condition" else "comparison"
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
      pathway_source = "link_scores",
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
      pathway_per_comparison_dir = "per_comparison_pathway",
      pathway_per_comparison_flat = FALSE,
      pathway_split_direction = TRUE,
      run_link_topic_scores = TRUE,
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
      run_doc_topic_heatmaps = TRUE,
      run_tf_topic_heatmaps = TRUE,
      run_topic_by_comparison_heatmaps = TRUE,
      run_topic_marker_heatmap = TRUE,
      run_intertopic_distance_map = TRUE,
      run_ldavis = TRUE,
      topic_by_comparison_label_cleaner = wrap_comp_label
    )
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

#' Extract Module 3 regulatory topics
#'
#' Public step function for extracting regulatory topics, pathway summaries,
#' topic-link tables, and review outputs from trained Module 3 topic models.
#'
#' @param k Integer K selected for extraction.
#' @param model_dir Directory containing trained topic model outputs.
#' @param output_dir Directory to write extracted topic outputs.
#' @param flatten_single_output Whether to write a single selected model
#'   directly under `output_dir`. Defaults to `TRUE` for the public step API.
#' @param ... Additional arguments passed to the internal extraction engine,
#'   such as `backend`, `doc_mode`, `weight_label`, and `topic_report_args`.
#'
#' @return Invisibly returns TRUE when extraction completes.
#' @export
module3_extract_topics <- function(k,
                                   model_dir,
                                   output_dir,
                                   flatten_single_output = TRUE,
                                   ...) {
  extract_regulatory_topics(
    k = k,
    model_dir = model_dir,
    output_dir = output_dir,
    flatten_single_output = flatten_single_output,
    ...
  )
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
                                          fp_term_mode = c("aggregate", "unique", "aggregate_weight"),
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

  dt <- data.table::as.data.table(edges_docs)
  .topic_assert_has_cols(dt, c("doc_id", "tf", "gene_key", "peak_id", weight_type_peak, weight_type_gene), context = "build_doc_term_unified")
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

  fp_dt <- dt[!is.na(peak_id) & nzchar(peak_id) & !is.na(gene_key) & nzchar(gene_key) & is.finite(fp_weight) & fp_weight > 0]
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

  if (identical(fp_term_mode, "aggregate_weight")) {
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
  out <- .topic_apply_counts(out, count_method = count_method, count_scale = count_scale)
  out[, .(
    weight = sum(weight, na.rm = TRUE),
    pseudo_count = sum(pseudo_count, na.rm = TRUE),
    pseudo_count_bin = sum(pseudo_count_bin, na.rm = TRUE),
    pseudo_count_log = sum(pseudo_count_log, na.rm = TRUE)
  ), by = .(doc_id, term_id)]
}

#' Build joint gene-plus-peak document-term rows
#'
#' Builds comparison-level document-term rows from links that already have a
#' `doc_id`. Gene, peak, and optional TF self-term weights are converted to
#' positive values before count conversion. `fp_term_mode` controls whether FP
#' terms remain unique peaks, are aggregated by target gene, or are multiplied
#' into gene weights.
#'
#' @inheritParams build_doc_term_from_edges
#' @param weight_type_peak Peak weight column, usually `log2fc_fp` for
#'   comparison-mode runs.
#' @param weight_type_gene Gene weight column, usually `log2fc_gene` for
#'   comparison-mode runs.
#' @param fp_term_mode FP term mode: unique peaks, target-gene aggregated
#'   peaks, or gene-only aggregate weights.
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
#' @param fp_term_mode FP term mode: unique peaks, target-gene aggregated
#'   peaks, or gene-only aggregate weights.
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
                                           fp_term_mode = c("unique", "aggregate", "aggregate_weight"),
                                           balance_mode = c("min", "none"),
                                           check_repeated_values = TRUE) {
  .assert_pkg("data.table")
  count_method <- match.arg(count_method)
  balance_mode <- match.arg(balance_mode)
  fp_term_mode <- .topic_term_mode(fp_term_mode)

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
