#' Prepare a GRN set for condition-aware linking
#'
#' Ensures the `grn_set` carries the right slots for `light_by_condition()`
#' (normalized signal tables and labeled metadata) and disambiguates any
#' duplicate sample labels by appending the sample id.
#'
#' @param grn_set List returned by `build_grn_set()` plus downstream slots.
#' @param label_col Character scalar naming the metadata column whose values
#'   should label conditions (default `"strict_match_rna"`).
#' @return `grn_set` clone ready for `light_by_condition()`.
#' @export
prepare_grn_set_for_light_by_condition <- function(grn_set, label_col = "strict_match_rna") {
  stopifnot(is.list(grn_set), !is.null(grn_set$sample_metadata_used))
  ds <- grn_set
  if (is.data.frame(grn_set$fp_score_condition_qn)) {
    ds$fp_score <- grn_set$fp_score_condition_qn
  }
  if (is.data.frame(grn_set$rna_condition)) {
    ds$rna <- grn_set$rna_condition
  }
  ds$sample_metadata_used <- ensure_grn_sample_metadata_label(grn_set, label_col)
  ds
}

ensure_grn_sample_metadata_label <- function(grn_set, label_col) {
  meta <- grn_set$sample_metadata_used
  stopifnot(is.data.frame(meta), "id" %in% names(meta))
  pick <- label_col
  if (!(pick %in% names(meta))) {
    pick <- "id"
  }
  label <- meta[[pick]]
  missing <- which(is.na(label) | label == "")
  if (length(missing)) {
    label[missing] <- meta$id[missing]
  }
  if (anyDuplicated(label)) {
    dup_idx <- which(duplicated(label) | duplicated(label, fromLast = TRUE))
    for (i in dup_idx) {
      label[i] <- paste0(label[i], "_", meta$id[i])
    }
  }
  # Keep sample id stable; only update the label column.
  meta[[pick]] <- label
  meta
}

#' Extract per-condition link tables from lighting output (or run wrapper)
#'
#' This function supports two modes:
#' 1. Reader mode (existing behavior): provide `out_dir` and it discovers and
#'    summarizes per-condition tables written by [light_by_condition()].
#' 2. Wrapper mode (Module 2 quick-start): provide `omics_data` and `links`
#'    (or call positionally as `extract_link_info_by_condition(omics_data, links, ...)`).
#'    It will run:
#'    [make_basal_links()] -> [prepare_grn_set_for_light_by_condition()] ->
#'    [light_by_condition()] before extracting per-condition outputs.
#'
#' @param out_dir Directory containing per-condition link CSVs (reader mode), or
#'   output directory to write to (wrapper mode). If omitted in wrapper mode,
#'   inferred from `base_dir` in `.GlobalEnv`.
#' @param prefix Optional prefix used by [light_by_condition()] output files.
#'   Default is empty string (no prefix).
#' @param read_tables If TRUE, read and return a combined table with
#'   a `condition` column.
#' @param verbose Emit status messages.
#' @param omics_data Optional GRN object used in wrapper mode.
#' @param links Optional FP->gene link table for wrapper mode. Can be either a
#'   data.frame (e.g., `fp_gene_corr_use`) or a list containing
#'   `fp_gene_corr_kept`.
#' @param qc If TRUE, write a simple per-condition row-count CSV/PDF summary.
#' @param label_col Condition label column used by
#'   [prepare_grn_set_for_light_by_condition()] and [light_by_condition()].
#' @param link_score_threshold Passed to [light_by_condition()].
#' @param fp_score_threshold Passed to [light_by_condition()].
#' @param tf_expr_threshold Passed to [light_by_condition()].
#' @param fp_bound_tbl Optional; defaults to `omics_data$fp_bound_condition`.
#' @param rna_expressed_tbl Optional; defaults to `omics_data$rna_expressed`.
#' @param atac_score_tbl Optional ATAC table for [light_by_condition()].
#' @param atac_score_threshold Passed to [light_by_condition()].
#' @param require_atac_score Passed to [light_by_condition()].
#' @param fp_annotation_tbl Optional; defaults to `omics_data$fp_annotation`.
#' @param require_fp_bound Passed to [light_by_condition()].
#' @param require_gene_expr Passed to [light_by_condition()].
#' @param gene_expr_threshold Passed to [light_by_condition()].
#' @param filter_active Passed to [light_by_condition()].
#' @param use_parallel Passed to [light_by_condition()].
#' @param workers Passed to [light_by_condition()].
#' @param fp_variance_tbl Optional; defaults to `omics_data$fp_variance`.
#' @param rna_variance_tbl Optional; defaults to `omics_data$rna_variance`.
#' @param rna_tbl_for_basal Optional RNA table passed to [make_basal_links()].
#'   Defaults to `omics_data$rna_condition` if present, else `omics_data$rna`.
#' @param rna_method_for_basal Passed to [make_basal_links()].
#' @param rna_cores_for_basal Passed to [make_basal_links()].
#' @param tf_gene_cache_dir Optional cache directory for wrapper-mode TF-gene links.
#'   Defaults to `file.path(out_dir, "cache", "tf_gene_corr")`.
#' @param tf_gene_cache_tag Optional cache tag for wrapper-mode TF-gene links.
#' @param tf_gene_cache_resume If TRUE, reuse cached TF-gene links when available.
#' @param tf_gene_cache_verbose Emit cache status messages for TF-gene links.
#' @param basal_cache_dir Deprecated alias for `tf_gene_cache_dir`.
#' @param basal_cache_tag Deprecated alias for `tf_gene_cache_tag`.
#' @param basal_cache_resume Deprecated alias for `tf_gene_cache_resume`.
#' @param basal_cache_verbose Deprecated alias for `tf_gene_cache_verbose`.
#'
#' @return In reader mode, returns a manifest tibble (or `list(manifest, links)`
#'   when `read_tables = TRUE`). In wrapper mode, returns the same outputs after
#'   writing per-condition link matrices.
#' @export
extract_link_info_by_condition <- function(
  out_dir = NULL,
  prefix = "",
  read_tables = FALSE,
  verbose = TRUE,
  omics_data = NULL,
  links = NULL,
  qc = FALSE,
  label_col = "strict_match_rna",
  link_score_threshold = 0,
  fp_score_threshold = 1,
  tf_expr_threshold = 10,
  fp_bound_tbl = NULL,
  rna_expressed_tbl = NULL,
  atac_score_tbl = NULL,
  atac_score_threshold = 0,
  require_atac_score = FALSE,
  fp_annotation_tbl = NULL,
  require_fp_bound = TRUE,
  require_gene_expr = TRUE,
  gene_expr_threshold = 1L,
  filter_active = FALSE,
  use_parallel = TRUE,
  workers = 2L,
  fp_variance_tbl = NULL,
  rna_variance_tbl = NULL,
  rna_tbl_for_basal = NULL,
  rna_method_for_basal = "pearson",
  rna_cores_for_basal = workers,
  tf_gene_cache_dir = NULL,
  tf_gene_cache_tag = NULL,
  tf_gene_cache_resume = TRUE,
  tf_gene_cache_verbose = TRUE,
  basal_cache_dir = NULL,
  basal_cache_tag = NULL,
  basal_cache_resume = NULL,
  basal_cache_verbose = NULL
) {
  # Backward-compatible quick-start positional call:
  # extract_link_info_by_condition(omics_data, links, qc = TRUE)
  if (is.null(omics_data) && is.list(out_dir)) {
    omics_data <- out_dir
    out_dir <- NULL
    if (!is.character(prefix)) {
      links <- prefix
      prefix <- ""
    }
  }

  # Wrapper mode: generate per-condition tables first.
  if (!is.null(omics_data)) {
    if (!is.list(omics_data)) .log_abort("`omics_data` must be a list in wrapper mode.")

    if (is.null(out_dir) || !nzchar(out_dir)) {
      base_dir <- if (exists("base_dir", envir = .GlobalEnv, inherits = FALSE)) {
        get("base_dir", envir = .GlobalEnv)
      } else {
        NULL
      }
      if (!is.character(base_dir) || !nzchar(base_dir)) {
        .log_abort("`out_dir` is required in wrapper mode when `base_dir` is not set.")
      }
      tf_mode <- if (exists("tf_binding_mode", envir = .GlobalEnv, inherits = FALSE)) {
        as.character(get("tf_binding_mode", envir = .GlobalEnv))
      } else {
        "canonical"
      }
      tf_mode <- match.arg(tf_mode, c("canonical", "all"))
      out_dir <- if (identical(tf_mode, "all")) {
        file.path(base_dir, "connect_tf_target_genes_all")
      } else {
        file.path(base_dir, "connect_tf_target_genes")
      }
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    fp_gene_corr_kept <- NULL
    if (is.data.frame(links)) {
      fp_gene_corr_kept <- links
    } else if (is.list(links) && is.data.frame(links$fp_gene_corr_kept)) {
      fp_gene_corr_kept <- links$fp_gene_corr_kept
    }
    if (!is.data.frame(fp_gene_corr_kept)) {
      .log_abort("Wrapper mode requires `links` as a data.frame or list containing `fp_gene_corr_kept`.")
    }

    if (is.null(fp_annotation_tbl)) fp_annotation_tbl <- omics_data$fp_annotation
    if (!is.data.frame(fp_annotation_tbl)) {
      .log_abort("Wrapper mode requires `fp_annotation_tbl` or `omics_data$fp_annotation`.")
    }
    if (is.null(fp_bound_tbl)) fp_bound_tbl <- omics_data$fp_bound_condition
    if (is.null(rna_expressed_tbl)) rna_expressed_tbl <- omics_data$rna_expressed
    if (is.null(fp_variance_tbl)) fp_variance_tbl <- omics_data$fp_variance
    if (is.null(rna_variance_tbl)) rna_variance_tbl <- omics_data$rna_variance
    if (is.null(rna_tbl_for_basal)) {
      if (is.data.frame(omics_data$rna_condition)) {
        rna_tbl_for_basal <- omics_data$rna_condition
      } else {
        rna_tbl_for_basal <- omics_data$rna
      }
    }
    if (!is.data.frame(rna_tbl_for_basal)) {
      .log_abort("Wrapper mode requires `rna_tbl_for_basal` (or `omics_data$rna_condition` / `omics_data$rna`).")
    }

    basal_prefix <- if (is.character(prefix) && nzchar(prefix)) prefix else "step2"
    tf_mode_for_cache <- if (exists("tf_binding_mode", envir = .GlobalEnv, inherits = FALSE)) {
      as.character(get("tf_binding_mode", envir = .GlobalEnv))
    } else {
      "canonical"
    }
    tf_mode_for_cache <- match.arg(tf_mode_for_cache, c("canonical", "all"))
    link_mode_for_cache <- if (exists("tf_target_link_mode", envir = .GlobalEnv, inherits = FALSE)) {
      as.character(get("tf_target_link_mode", envir = .GlobalEnv))
    } else {
      "genehancer"
    }
    link_mode_for_cache <- match.arg(link_mode_for_cache, c("genehancer", "window", "loops"))

    # Backward-compatible aliases.
    if (is.null(tf_gene_cache_dir) && !is.null(basal_cache_dir)) tf_gene_cache_dir <- basal_cache_dir
    if (is.null(tf_gene_cache_tag) && !is.null(basal_cache_tag)) tf_gene_cache_tag <- basal_cache_tag
    if (is.null(basal_cache_resume)) basal_cache_resume <- tf_gene_cache_resume
    if (is.null(basal_cache_verbose)) basal_cache_verbose <- tf_gene_cache_verbose
    tf_gene_cache_resume <- isTRUE(basal_cache_resume)
    tf_gene_cache_verbose <- isTRUE(basal_cache_verbose)

    if (is.null(tf_gene_cache_dir) || !nzchar(tf_gene_cache_dir)) {
      tf_gene_cache_dir <- file.path(out_dir, "cache", "tf_gene_corr")
    }
    dir.create(tf_gene_cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(tf_gene_cache_tag) || !nzchar(tf_gene_cache_tag)) {
      tf_gene_cache_tag <- paste0("nutrient_", tf_mode_for_cache, "_", link_mode_for_cache, "_", rna_method_for_basal)
    }
    tf_gene_cache_root <- file.path(tf_gene_cache_dir, tf_gene_cache_tag)
    dir.create(tf_gene_cache_root, recursive = TRUE, showWarnings = FALSE)
    tf_gene_cache_file <- file.path(tf_gene_cache_root, "tf_gene_corr.rds")

    use_cached_basal <- isTRUE(tf_gene_cache_resume) && file.exists(tf_gene_cache_file)
    if (isTRUE(use_cached_basal)) {
      basal_links <- tryCatch(readRDS(tf_gene_cache_file), error = function(e) NULL)
      if (!is.data.frame(basal_links) || !all(c("TF", "gene_key", "peak_ID") %in% names(basal_links))) {
        if (isTRUE(tf_gene_cache_verbose)) {
          .log_warn("TF-gene cache invalid at {tf_gene_cache_file}; recomputing.")
        }
        use_cached_basal <- FALSE
      } else if (isTRUE(tf_gene_cache_verbose)) {
        .log_inform("Using cached TF-gene links: {tf_gene_cache_file}")
      }
    }
    if (!isTRUE(use_cached_basal)) {
      basal_links <- make_basal_links(
        fp_gene_corr_kept = fp_gene_corr_kept,
        fp_annotation = fp_annotation_tbl,
        out_dir = file.path(out_dir, "cache"),
        prefix = basal_prefix,
        rna_tbl = rna_tbl_for_basal,
        rna_method = rna_method_for_basal,
        rna_cores = rna_cores_for_basal,
        fp_variance = fp_variance_tbl,
        rna_variance = rna_variance_tbl
      )
      saveRDS(basal_links, tf_gene_cache_file)
      if (isTRUE(tf_gene_cache_verbose)) {
        .log_inform("Saved TF-gene cache: {tf_gene_cache_file}")
      }
    }

    ds <- prepare_grn_set_for_light_by_condition(
      omics_data,
      label_col = label_col
    )

    light_by_condition(
      ds = ds,
      basal_links = basal_links,
      out_dir = out_dir,
      prefix = prefix,
      label_col = label_col,
      link_score_threshold = link_score_threshold,
      fp_score_threshold = fp_score_threshold,
      tf_expr_threshold = tf_expr_threshold,
      fp_bound_tbl = fp_bound_tbl,
      rna_expressed_tbl = rna_expressed_tbl,
      atac_score_tbl = atac_score_tbl,
      atac_score_threshold = atac_score_threshold,
      require_atac_score = require_atac_score,
      fp_annotation_tbl = fp_annotation_tbl,
      require_fp_bound = require_fp_bound,
      require_gene_expr = require_gene_expr,
      gene_expr_threshold = gene_expr_threshold,
      filter_active = filter_active,
      use_parallel = use_parallel,
      workers = workers,
      fp_variance_tbl = fp_variance_tbl,
      rna_variance_tbl = rna_variance_tbl
    )
  }

  if (!is.character(out_dir) || !nzchar(out_dir)) .log_abort("`out_dir` must be a non-empty path.")
  if (!dir.exists(out_dir)) .log_abort("`out_dir` does not exist: {out_dir}")
  matrices_dir <- file.path(out_dir, "per_condition_link_matrices")
  search_dir <- if (dir.exists(matrices_dir)) matrices_dir else out_dir

  pattern_new <- if (is.character(prefix) && nzchar(prefix)) {
    sprintf("^%s_.*_tf_gene_links\\.csv$", prefix)
  } else {
    "^.*_tf_gene_links\\.csv$"
  }
  pattern_old <- if (is.character(prefix) && nzchar(prefix)) {
    sprintf("^%s_cond-.*_tf_gene_links\\.csv$", prefix)
  } else {
    "^cond-.*_tf_gene_links\\.csv$"
  }
  csvs <- unique(c(
    list.files(search_dir, pattern = pattern_new, full.names = TRUE),
    list.files(search_dir, pattern = pattern_old, full.names = TRUE)
  ))
  if (!length(csvs)) .log_abort("No per-condition link CSVs found in {out_dir} with prefix {prefix}.")

  parse_label <- function(path) {
    base <- basename(path)
    x <- sub("_tf_gene_links\\.csv$", "", base)
    if (is.character(prefix) && nzchar(prefix)) {
      x <- sub(sprintf("^%s_cond-", prefix), "", x)
      x <- sub(sprintf("^%s_", prefix), "", x)
    } else {
      x <- sub("^cond-", "", x)
    }
    x
  }

  manifest <- tibble::tibble(
    condition = vapply(csvs, parse_label, character(1)),
    path = csvs
  )

  manifest$n_links <- vapply(csvs, function(p) {
    tryCatch(nrow(readr::read_csv(p, show_col_types = FALSE)), error = function(e) NA_integer_)
  }, integer(1))

  if (isTRUE(qc)) {
    qc_stub <- if (is.character(prefix) && nzchar(prefix)) paste0(prefix, "_") else ""
    qc_csv <- file.path(out_dir, paste0(qc_stub, "per_condition_link_counts.csv"))
    readr::write_csv(manifest, qc_csv)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(manifest, ggplot2::aes(x = condition, y = n_links)) +
        ggplot2::geom_col(fill = "#3182bd") +
        ggplot2::labs(
          title = "Per-condition TF->TFBS->target link rows",
          x = "Condition",
          y = "Link rows"
        ) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
        )
      qc_pdf <- file.path(out_dir, paste0(qc_stub, "per_condition_link_counts.pdf"))
      ggplot2::ggsave(qc_pdf, p, width = 10, height = 6, units = "in")
      if (isTRUE(verbose)) .log_inform("Per-condition link-count QC saved: {qc_pdf}")
    } else if (isTRUE(verbose)) {
      .log_warn("Skipping per-condition link-count PDF because ggplot2 is not installed.")
    }
    if (isTRUE(verbose)) .log_inform("Per-condition link-count summary written: {qc_csv}")
  }

  if (isTRUE(read_tables)) {
    combined <- purrr::map2_dfr(csvs, manifest$condition, function(p, cond) {
      readr::read_csv(p, show_col_types = FALSE) |>
        dplyr::mutate(condition = cond)
    })
    if (isTRUE(verbose)) .log_inform("Loaded {nrow(combined)} per-condition links.")
    return(invisible(list(manifest = manifest, links = combined)))
  }

  if (isTRUE(verbose)) {
    .log_inform("Found {nrow(manifest)} per-condition link tables in {out_dir}.")
  }
  invisible(manifest)
}
