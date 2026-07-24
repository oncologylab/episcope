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

.module3_condition_gene_log2fc <- function(gene_expression,
                                           conditions,
                                           pseudocount = 1) {
  gene_expression <- as.matrix(gene_expression)
  conditions <- as.character(conditions)
  if (length(conditions) != 2L || anyNA(conditions) || any(!nzchar(conditions))) {
    .log_abort("conditions must contain exactly two non-empty condition IDs.")
  }
  if (is.null(rownames(gene_expression)) || is.null(colnames(gene_expression))) {
    .log_abort("gene_expression must have gene row names and condition column names.")
  }
  missing_conditions <- setdiff(conditions, colnames(gene_expression))
  if (length(missing_conditions)) {
    .log_abort("Condition(s) missing from gene_expression: {paste(missing_conditions, collapse = ', ')}")
  }
  pseudocount <- as.numeric(pseudocount[[1L]])
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    .log_abort("pseudocount must be a positive number.")
  }
  first <- as.numeric(gene_expression[, conditions[[1L]]])
  second <- as.numeric(gene_expression[, conditions[[2L]]])
  first[!is.finite(first) | first < 0] <- 0
  second[!is.finite(second) | second < 0] <- 0
  data.table::data.table(
    gene_key = rownames(gene_expression),
    condition_1 = conditions[[1L]],
    condition_2 = conditions[[2L]],
    expression_condition_1 = first,
    expression_condition_2 = second,
    log2fc_condition_1_vs_2 = log2((first + pseudocount) / (second + pseudocount))
  )
}

.module3_write_condition_table <- function(x, path, format) {
  if (identical(format, "parquet")) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      .log_abort("Package arrow is required to write Parquet condition links.")
    }
    arrow::write_parquet(x, path)
  } else if (identical(format, "csv")) {
    data.table::fwrite(x, path)
  } else {
    .log_abort("Unsupported condition-link format: {format}")
  }
  invisible(path)
}

# Filter condition-native links to genes with a pairwise expression change.
#
# This helper stages a new manifest rather than modifying source links. The
# fold change is calculated from the condition-level RNA matrix so a missing
# regulatory link is not mistaken for zero gene expression.
#
# @noRd
.module3_filter_condition_links_by_gene_log2fc <- function(input_dir,
                                                            output_dir,
                                                            multiomic_data,
                                                            conditions,
                                                            abs_log2fc_min = 1,
                                                            pseudocount = 1,
                                                            overwrite = FALSE,
                                                            verbose = TRUE) {
  .assert_pkg("data.table")
  if (!is_multiomic_object(multiomic_data)) {
    .log_abort("multiomic_data must be a CraftGRN multiomic object.")
  }
  conditions <- as.character(conditions)
  if (length(conditions) != 2L) {
    .log_abort("Pairwise condition-link filtering requires exactly two conditions.")
  }
  abs_log2fc_min <- as.numeric(abs_log2fc_min[[1L]])
  if (!is.finite(abs_log2fc_min) || abs_log2fc_min < 0) {
    .log_abort("abs_log2fc_min must be a non-negative number.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(output_dir, "condition_links_manifest.csv")
  summary_path <- file.path(output_dir, "condition_links_summary.csv")
  de_path <- file.path(output_dir, "condition_gene_log2fc.csv")
  if (!isTRUE(overwrite) && all(file.exists(c(manifest_path, summary_path, de_path)))) {
    if (isTRUE(verbose)) .log_inform("Reusing differential-gene-filtered condition links in {output_dir}.")
    return(invisible(list(
      manifest = data.table::fread(manifest_path, showProgress = FALSE),
      summary = data.table::fread(summary_path, showProgress = FALSE),
      differential_genes = data.table::fread(de_path, showProgress = FALSE)
    )))
  }

  source_manifest_path <- .module3_condition_link_manifest_path(input_dir)
  source_manifest <- data.table::fread(source_manifest_path, showProgress = FALSE)
  .topic_assert_has_cols(source_manifest, c("condition_id", "path", "format"), context = "condition-link manifest")
  source_manifest <- source_manifest[condition_id %in% conditions]
  source_manifest <- source_manifest[match(conditions, condition_id)]
  if (nrow(source_manifest) != 2L || anyNA(source_manifest$condition_id)) {
    .log_abort("The source condition-link manifest does not contain both requested conditions.")
  }

  condition_map <- .module3_condition_matrix_map(multiomic_data, conditions)
  matrix_conditions <- condition_map$matrix_condition_id[
    match(conditions, condition_map$condition_id)
  ]
  gene_fc <- .module3_condition_gene_log2fc(
    multiomic_data$matrices$gene_expr,
    conditions = matrix_conditions,
    pseudocount = pseudocount
  )
  gene_fc[, `:=`(
    condition_1 = conditions[[1L]],
    condition_2 = conditions[[2L]],
    abs_log2fc = abs(log2fc_condition_1_vs_2),
    pass_abs_log2fc = abs(log2fc_condition_1_vs_2) >= abs_log2fc_min,
    abs_log2fc_min = abs_log2fc_min,
    pseudocount = as.numeric(pseudocount)
  )]
  data.table::setorder(gene_fc, -abs_log2fc, gene_key)
  data.table::fwrite(gene_fc, de_path)
  passing_genes <- gene_fc[pass_abs_log2fc %in% TRUE, gene_key]
  if (!length(passing_genes)) {
    .log_abort("No genes pass abs(log2FC) >= {abs_log2fc_min} for the requested condition pair.")
  }

  manifest_rows <- vector("list", nrow(source_manifest))
  summary_rows <- vector("list", nrow(source_manifest))
  for (i in seq_len(nrow(source_manifest))) {
    condition <- as.character(source_manifest$condition_id[[i]])
    format <- tolower(as.character(source_manifest$format[[i]]))
    source <- as.character(source_manifest$path[[i]])
    links <- .module3_read_table(source, format, allow_missing_columns = TRUE)
    links <- data.table::as.data.table(links)
    .topic_assert_has_cols(links, c("condition_id", "gene_key"), context = "condition links")
    n_source <- nrow(links)
    links <- links[condition_id == condition & gene_key %in% passing_genes]
    extension <- if (identical(format, "parquet")) "parquet" else "csv"
    out_path <- file.path(
      output_dir,
      paste0(.module3_safe_label(condition), "_condition_links.", extension)
    )
    .module3_write_condition_table(links, out_path, format)
    manifest_rows[[i]] <- data.table::data.table(
      condition_id = condition,
      matrix_condition_id = if ("matrix_condition_id" %in% names(source_manifest)) {
        as.character(source_manifest$matrix_condition_id[[i]])
      } else {
        NA_character_
      },
      path = normalizePath(out_path, winslash = "/", mustWork = TRUE),
      format = format,
      n_rows_scanned = as.double(n_source),
      n_links = as.double(nrow(links)),
      filter = "target_gene_abs_log2fc",
      abs_log2fc_min = abs_log2fc_min,
      condition_1 = conditions[[1L]],
      condition_2 = conditions[[2L]]
    )
    summary_rows[[i]] <- data.table::data.table(
      condition_id = condition,
      n_links_before_gene_filter = as.double(n_source),
      n_links = as.double(nrow(links)),
      n_tfs = as.double(data.table::uniqueN(links$tf)),
      n_target_genes = as.double(data.table::uniqueN(links$gene_key)),
      n_peaks = as.double(data.table::uniqueN(links$peak_id)),
      n_pairwise_de_genes = as.double(length(passing_genes)),
      abs_log2fc_min = abs_log2fc_min,
      pseudocount = as.numeric(pseudocount)
    )
    rm(links)
    invisible(gc())
  }
  manifest <- data.table::rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
  summary <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  data.table::fwrite(manifest, manifest_path)
  data.table::fwrite(summary, summary_path)
  if (isTRUE(verbose)) {
    .log_inform("Retained {length(passing_genes)} genes with abs(log2FC) >= {abs_log2fc_min}; wrote filtered links to {output_dir}.")
  }
  invisible(list(
    manifest = manifest,
    summary = summary,
    differential_genes = gene_fc,
    manifest_path = manifest_path,
    summary_path = summary_path,
    differential_genes_path = de_path
  ))
}

.module3_normalize_condition_comparisons <- function(comparisons) {
  comparisons <- data.table::as.data.table(comparisons)
  first_col <- intersect(c("cond1_id", "condition1_id", "condition_1"), names(comparisons))
  second_col <- intersect(c("cond2_id", "condition2_id", "condition_2"), names(comparisons))
  if (!length(first_col) || !length(second_col)) {
    .log_abort(
      "comparisons must include cond1_id/cond2_id or condition1_id/condition2_id columns."
    )
  }
  comparison_col <- intersect(c("comparison_id", "comparison", "comparison_label"), names(comparisons))
  out <- comparisons[, .(
    comparison_id = if (length(comparison_col)) as.character(get(comparison_col[[1L]])) else NA_character_,
    condition_1 = as.character(get(first_col[[1L]])),
    condition_2 = as.character(get(second_col[[1L]]))
  )]
  out <- out[
    !is.na(condition_1) & nzchar(condition_1) &
      !is.na(condition_2) & nzchar(condition_2) &
      condition_1 != condition_2
  ]
  if (!nrow(out)) .log_abort("comparisons does not contain any usable condition pairs.")
  out[is.na(comparison_id) | !nzchar(comparison_id),
      comparison_id := paste0(condition_1, "_vs_", condition_2)]
  out[, pair_key := vapply(
    seq_len(.N),
    function(i) paste(sort(c(condition_1[[i]], condition_2[[i]])), collapse = "\r"),
    character(1)
  )]
  out <- unique(out, by = "pair_key")
  out[, pair_key := NULL]
  out[]
}

# Build each condition's gene union from comparisons involving that condition,
# then retain union genes expressed in that condition. Link files are read and
# written one condition at a time to keep peak memory bounded.
#
# @noRd
.module3_filter_condition_links_by_comparison_union <- function(input_dir,
                                                                 output_dir,
                                                                 gene_expression,
                                                                 comparisons,
                                                                 conditions = NULL,
                                                                 abs_log2fc_min = 1,
                                                                 expression_min = 0,
                                                                 pseudocount = 1,
                                                                 overwrite = FALSE,
                                                                 verbose = TRUE) {
  .assert_pkg("data.table")
  gene_expression <- as.matrix(gene_expression)
  if (is.null(rownames(gene_expression)) || is.null(colnames(gene_expression))) {
    .log_abort("gene_expression must have gene row names and condition column names.")
  }
  if (anyDuplicated(rownames(gene_expression))) {
    .log_abort("gene_expression gene row names must be unique.")
  }
  abs_log2fc_min <- as.numeric(abs_log2fc_min[[1L]])
  expression_min <- as.numeric(expression_min[[1L]])
  pseudocount <- as.numeric(pseudocount[[1L]])
  if (!is.finite(abs_log2fc_min) || abs_log2fc_min < 0) {
    .log_abort("abs_log2fc_min must be a non-negative number.")
  }
  if (!is.finite(expression_min) || expression_min < 0) {
    .log_abort("expression_min must be a non-negative number.")
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    .log_abort("pseudocount must be a positive number.")
  }

  comparison_map <- .module3_normalize_condition_comparisons(comparisons)
  requested <- if (is.null(conditions)) {
    unique(c(comparison_map$condition_1, comparison_map$condition_2))
  } else {
    unique(as.character(conditions))
  }
  requested <- requested[!is.na(requested) & nzchar(requested)]
  comparison_map <- comparison_map[
    condition_1 %in% requested & condition_2 %in% requested
  ]
  involved <- unique(c(comparison_map$condition_1, comparison_map$condition_2))
  missing_comparisons <- setdiff(requested, involved)
  if (length(missing_comparisons)) {
    .log_abort(
      "Condition(s) have no configured comparison: {paste(missing_comparisons, collapse = ', ')}"
    )
  }
  missing_expression <- setdiff(involved, colnames(gene_expression))
  if (length(missing_expression)) {
    .log_abort(
      "Condition(s) missing from gene_expression: {paste(missing_expression, collapse = ', ')}"
    )
  }
  signature_rows <- comparison_map[
    order(comparison_id, condition_1, condition_2),
    .(comparison_id, condition_1, condition_2)
  ]
  comparison_signature <- digest::digest(
    paste(
      signature_rows$comparison_id,
      signature_rows$condition_1,
      signature_rows$condition_2,
      sep = "\t",
      collapse = "\n"
    ),
    algo = "xxhash64",
    serialize = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(output_dir, "condition_links_manifest.csv")
  summary_path <- file.path(output_dir, "condition_links_summary.csv")
  comparison_path <- file.path(output_dir, "condition_comparison_gene_log2fc.csv")
  global_union_path <- file.path(output_dir, "global_differential_gene_union.csv")
  condition_union_path <- file.path(output_dir, "condition_differential_gene_union.csv")
  expression_path <- file.path(output_dir, "condition_gene_expression.csv")
  comparison_summary_path <- file.path(output_dir, "condition_comparison_gene_filter_summary.csv")
  required_outputs <- c(
    manifest_path, summary_path, comparison_path, global_union_path,
    condition_union_path, expression_path, comparison_summary_path
  )
  if (!isTRUE(overwrite) && all(file.exists(required_outputs))) {
    existing_manifest <- data.table::fread(manifest_path, showProgress = FALSE)
    reusable <- all(c(
      "filter_scope", "comparison_signature", "abs_log2fc_min",
      "expression_min", "pseudocount"
    ) %in% names(existing_manifest)) &&
      all(existing_manifest$filter_scope == "condition_comparison_union_then_condition_expression") &&
      all(existing_manifest$comparison_signature == comparison_signature) &&
      all(as.numeric(existing_manifest$abs_log2fc_min) == abs_log2fc_min) &&
      all(as.numeric(existing_manifest$expression_min) == expression_min) &&
      all(as.numeric(existing_manifest$pseudocount) == pseudocount)
    if (isTRUE(reusable)) {
      if (isTRUE(verbose)) {
        .log_inform("Reusing condition-comparison-union-filtered links in {output_dir}.")
      }
      return(invisible(list(
        manifest = existing_manifest,
        summary = data.table::fread(summary_path, showProgress = FALSE),
        comparison_genes = data.table::fread(comparison_path, showProgress = FALSE),
        global_gene_union = data.table::fread(global_union_path, showProgress = FALSE),
        condition_gene_candidates = data.table::fread(
          condition_union_path,
          showProgress = FALSE
        ),
        condition_gene_expression = data.table::fread(expression_path, showProgress = FALSE)
      )))
    }
    if (isTRUE(verbose)) {
      .log_warn("Existing condition links use an obsolete filter scope or expression threshold; rebuilding them.")
    }
  }

  comparison_rows <- lapply(seq_len(nrow(comparison_map)), function(i) {
    row <- comparison_map[i]
    values <- .module3_condition_gene_log2fc(
      gene_expression,
      conditions = c(row$condition_1[[1L]], row$condition_2[[1L]]),
      pseudocount = pseudocount
    )
    data.table::setnames(
      values,
      c(
        "expression_condition_1", "expression_condition_2",
        "log2fc_condition_1_vs_2"
      ),
      c("expression_1", "expression_2", "log2fc")
    )
    values[, `:=`(
      comparison_id = row$comparison_id[[1L]],
      condition_1 = row$condition_1[[1L]],
      condition_2 = row$condition_2[[1L]],
      abs_log2fc = abs(log2fc),
      pass_abs_log2fc = abs(log2fc) >= abs_log2fc_min,
      abs_log2fc_min = abs_log2fc_min,
      pseudocount = pseudocount
    )]
    values[, .(
      comparison_id, condition_1, condition_2, gene_key,
      expression_1, expression_2, log2fc, abs_log2fc,
      pass_abs_log2fc, abs_log2fc_min, pseudocount
    )]
  })
  comparison_genes <- data.table::rbindlist(comparison_rows, use.names = TRUE)
  data.table::setorder(comparison_genes, comparison_id, -abs_log2fc, gene_key)
  data.table::fwrite(comparison_genes, comparison_path)

  passing <- comparison_genes[pass_abs_log2fc %in% TRUE]
  if (!nrow(passing)) {
    .log_abort("No genes pass abs(log2FC) >= {abs_log2fc_min} in the configured comparisons.")
  }
  global_union <- passing[, .(
    max_abs_log2fc = max(abs_log2fc),
    n_comparisons_passed = data.table::uniqueN(comparison_id),
    supporting_comparisons = paste(sort(unique(comparison_id)), collapse = ";")
  ), by = gene_key]
  global_union[, `:=`(
    pass_abs_log2fc = TRUE,
    abs_log2fc_min = abs_log2fc_min,
    pseudocount = pseudocount
  )]
  data.table::setorder(global_union, -max_abs_log2fc, gene_key)
  data.table::fwrite(global_union, global_union_path)

  condition_union <- data.table::rbindlist(list(
    passing[, .(
      condition_id = condition_1,
      comparison_id,
      gene_key,
      abs_log2fc
    )],
    passing[, .(
      condition_id = condition_2,
      comparison_id,
      gene_key,
      abs_log2fc
    )]
  ))
  condition_union <- condition_union[, .(
    max_abs_log2fc = max(abs_log2fc),
    n_comparisons_passed = data.table::uniqueN(comparison_id),
    supporting_comparisons = paste(sort(unique(comparison_id)), collapse = ";")
  ), by = .(condition_id, gene_key)]
  condition_union[, `:=`(
    pass_abs_log2fc = TRUE,
    abs_log2fc_min = abs_log2fc_min,
    pseudocount = pseudocount
  )]
  data.table::setorder(condition_union, condition_id, -max_abs_log2fc, gene_key)
  data.table::fwrite(condition_union, condition_union_path)

  expression_rows <- lapply(involved, function(condition) {
    condition_genes <- condition_union[condition_id == condition, gene_key]
    data.table::data.table(
      condition_id = condition,
      gene_key = condition_genes,
      expression = as.numeric(gene_expression[condition_genes, condition]),
      expression_min = expression_min
    )
  })
  condition_expression <- data.table::rbindlist(expression_rows)
  condition_expression[, expressed := is.finite(expression) & expression >= expression_min]
  data.table::setorder(condition_expression, condition_id, gene_key)
  data.table::fwrite(condition_expression, expression_path)

  comparison_summary <- comparison_genes[, .(
    n_genes_tested = .N,
    n_genes_passing = sum(pass_abs_log2fc %in% TRUE),
    passing_fraction = mean(pass_abs_log2fc %in% TRUE)
  ), by = .(comparison_id, condition_1, condition_2)]
  comparison_summary[, `:=`(
    abs_log2fc_min = abs_log2fc_min,
    pseudocount = pseudocount
  )]
  data.table::fwrite(comparison_summary, comparison_summary_path)

  source_manifest_path <- .module3_condition_link_manifest_path(input_dir)
  source_manifest <- data.table::fread(source_manifest_path, showProgress = FALSE)
  .topic_assert_has_cols(
    source_manifest,
    c("condition_id", "path", "format"),
    context = "condition-link manifest"
  )
  source_manifest <- source_manifest[condition_id %in% involved]
  source_manifest <- source_manifest[match(involved, condition_id)]
  if (nrow(source_manifest) != length(involved) || anyNA(source_manifest$condition_id)) {
    .log_abort("The source condition-link manifest does not contain every compared condition.")
  }

  manifest_rows <- vector("list", nrow(source_manifest))
  summary_rows <- vector("list", nrow(source_manifest))
  for (i in seq_len(nrow(source_manifest))) {
    condition <- as.character(source_manifest$condition_id[[i]])
    format <- tolower(as.character(source_manifest$format[[i]]))
    source <- as.character(source_manifest$path[[i]])
    genes <- condition_expression[
      condition_id == condition & expressed %in% TRUE,
      gene_key
    ]
    links <- data.table::as.data.table(.module3_read_table(
      source,
      format,
      allow_missing_columns = TRUE
    ))
    .topic_assert_has_cols(links, c("condition_id", "gene_key"), context = "condition links")
    n_source <- nrow(links)
    links <- links[condition_id == condition & gene_key %in% genes]
    extension <- if (identical(format, "parquet")) "parquet" else "csv"
    out_path <- file.path(
      output_dir,
      paste0(.module3_safe_label(condition), "_condition_links.", extension)
    )
    .module3_write_condition_table(links, out_path, format)
    manifest_rows[[i]] <- data.table::data.table(
      condition_id = condition,
      matrix_condition_id = if ("matrix_condition_id" %in% names(source_manifest)) {
        as.character(source_manifest$matrix_condition_id[[i]])
      } else {
        NA_character_
      },
      path = normalizePath(out_path, winslash = "/", mustWork = TRUE),
      format = format,
      n_rows_scanned = as.double(n_source),
      n_links = as.double(nrow(links)),
      filter = "target_gene_condition_comparison_union_abs_log2fc_and_expression",
      filter_scope = "condition_comparison_union_then_condition_expression",
      comparison_signature = comparison_signature,
      abs_log2fc_min = abs_log2fc_min,
      expression_min = expression_min,
      pseudocount = pseudocount,
      n_global_differential_genes = as.double(nrow(global_union)),
      n_condition_differential_genes = as.double(
        condition_union[condition_id == condition, .N]
      ),
      n_expressed_condition_genes = as.double(length(genes))
    )
    summary_rows[[i]] <- data.table::data.table(
      condition_id = condition,
      n_links_before_gene_filter = as.double(n_source),
      n_links = as.double(nrow(links)),
      retained_link_fraction = if (n_source) nrow(links) / n_source else NA_real_,
      n_tfs = as.double(data.table::uniqueN(links$tf)),
      n_target_genes = as.double(data.table::uniqueN(links$gene_key)),
      n_peaks = if ("peak_id" %in% names(links)) {
        as.double(data.table::uniqueN(links$peak_id))
      } else {
        NA_real_
      },
      n_global_differential_genes = as.double(nrow(global_union)),
      n_condition_differential_genes = as.double(
        condition_union[condition_id == condition, .N]
      ),
      n_expressed_condition_genes = as.double(length(genes)),
      abs_log2fc_min = abs_log2fc_min,
      expression_min = expression_min,
      pseudocount = pseudocount
    )
    rm(links)
    invisible(gc())
  }
  manifest <- data.table::rbindlist(manifest_rows, use.names = TRUE, fill = TRUE)
  summary <- data.table::rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
  data.table::fwrite(manifest, manifest_path)
  data.table::fwrite(summary, summary_path)
  if (isTRUE(verbose)) {
    .log_inform(
      "Filtered {nrow(manifest)} condition-link files using condition-specific comparison-gene unions followed by condition expression in {output_dir}."
    )
  }
  invisible(list(
    manifest = manifest,
    summary = summary,
    comparison_genes = comparison_genes,
    global_gene_union = global_union,
    condition_gene_candidates = condition_union,
    condition_gene_expression = condition_expression,
    condition_gene_union = condition_expression[expressed %in% TRUE, .(
      condition_id, gene_key
    )],
    manifest_path = manifest_path,
    summary_path = summary_path,
    comparison_genes_path = comparison_path,
    global_gene_union_path = global_union_path,
    condition_gene_union_path = condition_union_path,
    condition_gene_expression_path = expression_path,
    comparison_summary_path = comparison_summary_path
  ))
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
