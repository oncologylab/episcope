# File: utils_step3_condition_links.R
# Purpose: Condition-native Module 3 link preparation from Module 2 outputs.

.module3_module2_link_manifests <- function(module2) {
  if (is.character(module2) && length(module2) == 1L) {
    module2 <- load_module2_links(module2)
  }
  links <- NULL
  if (is.list(module2)) {
    if (is.data.frame(module2$links) && ncol(module2$links)) {
      links <- data.table::as.data.table(module2$links)
    } else if (is.data.frame(module2$module2_links) && ncol(module2$module2_links)) {
      links <- data.table::as.data.table(module2$module2_links)
    }
  }
  if (is.data.frame(links) && all(c("tf", "fp_id", "target_gene") %in% names(links))) {
    return(list(in_memory = links, manifest = data.table::data.table()))
  }
  if (is.list(module2) && is.data.frame(module2$module2_links_manifest)) {
    return(list(in_memory = NULL, manifest = data.table::as.data.table(module2$module2_links_manifest)))
  }
  list(in_memory = NULL, manifest = .module3_table_rows(module2, "module2_links"))
}

.module3_read_module2_link_chunk <- function(source, i, columns = NULL) {
  if (!is.null(source$in_memory)) {
    dt <- data.table::copy(source$in_memory)
    if (!is.null(columns)) {
      keep <- intersect(columns, names(dt))
      dt <- dt[, keep, with = FALSE]
    }
    return(dt)
  }
  man <- source$manifest
  .module3_read_table(
    as.character(man$path[[i]]),
    as.character(man$format[[i]]),
    columns = columns,
    allow_missing_columns = TRUE
  )
}

.module3_condition_link_manifest_path <- function(input_dir) {
  candidates <- c(
    file.path(input_dir, "condition_links_manifest.csv"),
    file.path(input_dir, "condition_links", "condition_links_manifest.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) .log_abort("No condition_links_manifest.csv found in {input_dir}.")
  hit[[1L]]
}

.module3_read_condition_links <- function(input_dir, conditions = NULL) {
  manifest_path <- .module3_condition_link_manifest_path(input_dir)
  manifest <- data.table::fread(manifest_path, showProgress = FALSE)
  if (!all(c("path", "format") %in% names(manifest))) {
    .log_abort("condition_links_manifest.csv must include path and format columns.")
  }
  if (!"condition_id" %in% names(manifest)) manifest[, condition_id := NA_character_]
  conditions <- if (is.null(conditions)) NULL else unique(as.character(conditions))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  if (length(conditions)) {
    manifest <- manifest[is.na(condition_id) | condition_id %in% conditions]
  }
  if (!nrow(manifest)) .log_abort("No condition-link manifest rows remain after condition filtering.")

  unique_files <- unique(manifest[, .(path = as.character(path), format = as.character(format))])
  out <- data.table::rbindlist(lapply(seq_len(nrow(unique_files)), function(i) {
    .module3_read_table(
      unique_files$path[[i]],
      unique_files$format[[i]],
      allow_missing_columns = TRUE
    )
  }), use.names = TRUE, fill = TRUE)
  if (!nrow(out)) .log_abort("Condition-link input is empty.")
  if (length(conditions)) out <- out[condition_id %in% conditions]
  if (!nrow(out)) .log_abort("No condition-link rows remain after condition filtering.")
  out
}

.module3_condition_matrix_map <- function(multiomic_data, conditions) {
  mats <- multiomic_data$matrices
  fp_conditions <- colnames(mats$fp_score)
  gene_conditions <- colnames(mats$gene_expr)
  matrix_conditions <- intersect(fp_conditions, gene_conditions)
  if (is.null(conditions)) {
    return(data.table::data.table(
      condition_id = matrix_conditions,
      matrix_condition_id = matrix_conditions
    ))
  }
  conditions <- unique(as.character(conditions))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  project_id <- as.character(multiomic_data$project$project_id %||% "")
  stripped <- matrix_conditions
  if (nzchar(project_id)) {
    stripped <- sub(paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", project_id), "_"), "", stripped)
  }
  gse_stripped <- sub("^[A-Za-z]+[0-9]+_", "", matrix_conditions)
  resolve_one <- function(condition) {
    exact <- which(matrix_conditions == condition)
    if (length(exact)) return(matrix_conditions[[exact[[1L]]]])
    if (nzchar(project_id)) {
      prefixed <- paste0(project_id, "_", condition)
      hit <- which(matrix_conditions == prefixed)
      if (length(hit)) return(matrix_conditions[[hit[[1L]]]])
    }
    stripped_hit <- which(stripped == condition)
    if (length(stripped_hit)) return(matrix_conditions[[stripped_hit[[1L]]]])
    gse_hit <- which(gse_stripped == condition)
    if (length(gse_hit)) return(matrix_conditions[[gse_hit[[1L]]]])
    NA_character_
  }
  mapped <- vapply(conditions, resolve_one, character(1))
  missing <- conditions[is.na(mapped)]
  if (length(missing)) {
    .log_abort("Condition(s) not found in multiomic matrices: {paste(missing, collapse = ', ')}")
  }
  data.table::data.table(condition_id = conditions, matrix_condition_id = unname(mapped))
}

#' Prepare condition-native Module 3 links
#'
#' Builds per-condition TF-peak-gene links from Module 2 final links and the
#' CraftGRN multiomic condition matrices. This path is used for condition topic
#' models that should not depend on differential-link comparisons.
#'
#' @param module2 Module 2 output object or output directory.
#' @param multiomic_data CraftGRN multiomic object.
#' @param conditions Optional condition IDs to keep.
#' @param output_dir Directory for per-condition link files and manifest.
#' @param threshold_fp_score Minimum condition footprint score.
#' @param threshold_gene_expr Minimum target-gene expression.
#' @param threshold_tf_expr Minimum TF expression.
#' @param output_format Output format: auto, parquet, or csv.
#' @param overwrite Overwrite existing condition-link files.
#' @param verbose Emit concise progress messages.
#'
#' @return A list containing the manifest and summary tables.
#' @export
module3_prepare_condition_links <- function(module2,
                                            multiomic_data,
                                            conditions = NULL,
                                            output_dir,
                                            threshold_fp_score = 0,
                                            threshold_gene_expr = 0,
                                            threshold_tf_expr = -Inf,
                                            output_format = c("auto", "parquet", "csv"),
                                            overwrite = FALSE,
                                            verbose = TRUE) {
  .assert_pkg("data.table")
  output_format <- match.arg(output_format)
  if (!is_multiomic_object(multiomic_data)) .log_abort("`multiomic_data` must be a CraftGRN multiomic object.")
  validate_multiomic_object(multiomic_data)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(output_dir, "condition_links_manifest.csv")
  summary_path <- file.path(output_dir, "condition_links_summary.csv")
  if (!isTRUE(overwrite) && file.exists(manifest_path) && file.exists(summary_path)) {
    manifest <- data.table::fread(manifest_path, showProgress = FALSE)
    summary <- data.table::fread(summary_path, showProgress = FALSE)
    if (isTRUE(verbose)) .log_inform("Reusing condition-link manifest: {manifest_path}")
    return(invisible(list(manifest = manifest, summary = summary, manifest_path = manifest_path)))
  }

  fmt <- if (identical(output_format, "auto")) {
    if (requireNamespace("arrow", quietly = TRUE)) "parquet" else "csv"
  } else {
    output_format
  }
  if (identical(fmt, "parquet") && !requireNamespace("arrow", quietly = TRUE)) {
    .log_abort("Package arrow is required to write Parquet condition links.")
  }

  mats <- multiomic_data$matrices
  condition_map <- .module3_condition_matrix_map(multiomic_data, conditions)
  if (!nrow(condition_map)) .log_abort("No conditions supplied for condition-link preparation.")
  conditions <- condition_map$condition_id

  source <- .module3_module2_link_manifests(module2)
  n_chunks <- if (!is.null(source$in_memory)) 1L else nrow(source$manifest)
  if (!n_chunks) .log_abort("Module 2 links are missing or empty.")
  link_cols <- c(
    "link_id", "tf", "fp_id", "target_gene", "module2_link_pass",
    "tf_expression_target_r", "fp_target_rna_r", "distance_to_tss",
    "candidate_source", "within_tss_window", "prior_supported",
    "prior_id", "prior_source", "prior_score", "prior_status"
  )
  manifest_rows <- vector("list", length(conditions))
  summary_rows <- vector("list", length(conditions))
  if (isTRUE(verbose)) {
    .log_inform("Preparing condition-native links for {length(conditions)} condition(s) from {n_chunks} Module 2 link chunk(s).")
  }
  for (ci in seq_along(conditions)) {
    condition <- conditions[[ci]]
    matrix_condition <- condition_map$matrix_condition_id[[ci]]
    cond_safe <- .module3_safe_label(condition)
    out_path <- file.path(output_dir, paste0(cond_safe, "_condition_links.", if (identical(fmt, "parquet")) "parquet" else "csv"))
    if (file.exists(out_path) && isTRUE(overwrite)) unlink(out_path)
    cond_chunks <- vector("list", n_chunks)
    total_scanned <- 0
    total_kept <- 0
    for (i in seq_len(n_chunks)) {
      dt <- .module3_read_module2_link_chunk(source, i, columns = link_cols)
      dt <- data.table::as.data.table(dt)
      if (!nrow(dt)) next
      .topic_assert_has_cols(dt, c("tf", "fp_id", "target_gene"), context = "module3_prepare_condition_links")
      total_scanned <- total_scanned + nrow(dt)
      if ("module2_link_pass" %in% names(dt)) dt <- dt[module2_link_pass %in% TRUE]
      if (!nrow(dt)) next
      fp_idx <- match(as.character(dt$fp_id), rownames(mats$fp_score))
      tf_idx <- match(as.character(dt$tf), rownames(mats$gene_expr))
      gene_idx <- match(as.character(dt$target_gene), rownames(mats$gene_expr))
      fp_score <- as.numeric(mats$fp_score[fp_idx, matrix_condition])
      fp_bound <- as.logical(mats$fp_bound[fp_idx, matrix_condition])
      tf_expr <- as.numeric(mats$gene_expr[tf_idx, matrix_condition])
      gene_expr <- as.numeric(mats$gene_expr[gene_idx, matrix_condition])
      keep <- is.finite(fp_score) & fp_score >= threshold_fp_score &
        fp_bound %in% TRUE &
        is.finite(gene_expr) & gene_expr >= threshold_gene_expr &
        is.finite(tf_expr) & tf_expr >= threshold_tf_expr
      if (!any(keep)) next
      one <- data.table::copy(dt[keep])
      one[, `:=`(
        condition_id = condition,
        condition_label = condition,
        matrix_condition_id = matrix_condition,
        tf = as.character(tf),
        tf_doc = as.character(tf),
        gene_key = as.character(target_gene),
        peak_id = as.character(fp_id),
        fp_score_condition = fp_score[keep],
        fp_bound_condition = fp_bound[keep],
        gene_expr_condition = gene_expr[keep],
        tf_expr_condition = tf_expr[keep],
        active_link = TRUE
      )]
      one[, doc_id := paste(condition_id, tf_doc, sep = "::")]
      cond_chunks[[i]] <- one
      total_kept <- total_kept + nrow(one)
    }
    cond_dt <- data.table::rbindlist(cond_chunks, use.names = TRUE, fill = TRUE)
    if (!nrow(cond_dt)) {
      cond_dt <- data.table::data.table(
        condition_id = character(),
        condition_label = character(),
        matrix_condition_id = character(),
        doc_id = character(),
        tf_doc = character(),
        tf = character(),
        gene_key = character(),
        peak_id = character(),
        fp_score_condition = numeric(),
        fp_bound_condition = logical(),
        gene_expr_condition = numeric(),
        tf_expr_condition = numeric(),
        active_link = logical()
      )
    }
    if (identical(fmt, "parquet")) {
      arrow::write_parquet(cond_dt, out_path)
    } else {
      data.table::fwrite(cond_dt, out_path)
    }
    manifest_rows[[ci]] <- data.table::data.table(
      condition_id = condition,
      matrix_condition_id = matrix_condition,
      path = out_path,
      format = fmt,
      n_rows_scanned = as.double(total_scanned),
      n_links = as.double(total_kept)
    )
    summary_rows[[ci]] <- data.table::data.table(
      condition_id = condition,
      n_links = as.double(nrow(cond_dt)),
      n_tfs = as.double(data.table::uniqueN(cond_dt$tf)),
      n_target_genes = as.double(data.table::uniqueN(cond_dt$gene_key)),
      n_peaks = as.double(data.table::uniqueN(cond_dt$peak_id)),
      threshold_fp_score = as.numeric(threshold_fp_score),
      threshold_gene_expr = as.numeric(threshold_gene_expr),
      threshold_tf_expr = as.numeric(threshold_tf_expr)
    )
    rm(cond_chunks, cond_dt)
    invisible(gc())
  }
  manifest <- data.table::rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
  summary <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  data.table::fwrite(manifest, manifest_path)
  data.table::fwrite(summary, summary_path)
  if (isTRUE(verbose)) {
    .log_inform("Wrote condition-link manifest: {manifest_path}")
  }
  invisible(list(manifest = manifest, summary = summary, manifest_path = manifest_path, summary_path = summary_path))
}
