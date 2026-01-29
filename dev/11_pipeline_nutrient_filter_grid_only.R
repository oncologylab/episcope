# library(episcope)
# ## Driver: nutrient filter-grid only
# ## Intended usage: RStudio "Run Script as Background Job" with "copy global environment".
# ## This script assumes these objects already exist in the job environment:
# ##   - filter_grid (tibble/data.frame)
# ##   - fp_gene_corr_joined (data.frame)
# ##   - filter_fp_rna_atac_corr (function)
# ##   - run_lighting_pipeline_one (function)
# ##   - base_dir, regulation_priors, db
# ##   - threshold_fp_tf_corr_p, threshold_gene_expr, threshold_tf_expr, threshold_fp_score, threshold_link_score,
# ##     link_score_threshold, fp_score_threshold (used downstream inside run_lighting_pipeline_one)
# ##   - fp_p_col_use, rna_p_col_use, atac_p_col_use, tf_p_col_use
#
# options(
#   future.globals.maxSize = max(getOption("future.globals.maxSize", 500 * 1024^2), 32 * 1024^3)
# )
#
# required_syms <- c(
#   "filter_grid",
#   "fp_gene_corr_joined",
#   "filter_fp_rna_atac_corr",
#   "run_lighting_pipeline_one",
#   "base_dir",
#   "regulation_priors",
#   "db",
#   "threshold_fp_tf_corr_p",
#   "fp_p_col_use",
#   "rna_p_col_use",
#   "atac_p_col_use",
#   "tf_p_col_use"
# )
# missing_syms <- required_syms[!vapply(required_syms, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
# if (length(missing_syms) > 0) {
#   stop(
#     "Missing required objects in job environment: ", paste(missing_syms, collapse = ", "), "\n",
#     "Tip: Source dev/11_pipeline_nutrient.R first (to build data + functions), then run this as a Background Job with 'copy global environment'."
#   )
# }
#
# grid_i_start <- as.integer(Sys.getenv("EPISCOPE_GRID_I_START", unset = "1"))
# grid_i_end <- as.integer(Sys.getenv("EPISCOPE_GRID_I_END", unset = as.character(nrow(filter_grid))))
# continue_on_error <- tolower(Sys.getenv("EPISCOPE_CONTINUE_ON_ERROR", unset = "false")) %in% c("1", "true", "t", "yes", "y")
#
# grid_i_start <- max(1L, grid_i_start)
# grid_i_end <- min(as.integer(nrow(filter_grid)), grid_i_end)
# if (grid_i_start > grid_i_end) stop("EPISCOPE_GRID_I_START > EPISCOPE_GRID_I_END")
#
# p_label <- function(p_col) {
#   if (grepl("adj", p_col, ignore.case = TRUE)) "pAdj" else "p"
# }
#
# message("[filter_grid] running rows ", grid_i_start, "..", grid_i_end, " (of ", nrow(filter_grid), ")")
#
# lighting_folders <- character(0)
#
# for (i in seq.int(grid_i_start, grid_i_end)) {
#   cfg <- filter_grid[i, , drop = FALSE]
#
#   fp_p_thr <- cfg$fp_rna_p
#   fp_r_thr <- cfg$fp_rna_r
#   rna_p_thr <- cfg$fp_rna_p
#   rna_r_thr <- cfg$fp_rna_r
#
#   use_atac <- isTRUE(cfg$use_atac)
#   atac_p_thr <- if (use_atac) cfg$atac_p else NULL
#   atac_r_thr <- if (use_atac) cfg$atac_r else NULL
#   atac_same_dir <- if (use_atac) isTRUE(cfg$atac_same_dir) else FALSE
#
#   tf_p_thr <- cfg$tf_p
#   tf_r_thr <- cfg$tf_r
#
#   suffix_fp <- sprintf("FP_%s%.2g_r%.1f", p_label(fp_p_col_use), fp_p_thr, fp_r_thr)
#   suffix_rna <- sprintf("RNA_%s%.2g_r%.1f", p_label(rna_p_col_use), rna_p_thr, rna_r_thr)
#   suffix_tf <- {
#     if (!is.finite(tf_p_thr) || !is.finite(tf_r_thr)) stop("TF thresholds must be finite")
#     sprintf("TF_%s%.2g_r%.1f", p_label(tf_p_col_use), tf_p_thr, tf_r_thr)
#   }
#   suffix_atac <- if (!use_atac) {
#     "ATAC_none"
#   } else {
#     if (!is.finite(atac_p_thr) || !is.finite(atac_r_thr)) stop("ATAC thresholds must be finite when use_atac=TRUE")
#     sprintf(
#       "ATAC_%s%.2g_r%.1f_sameDir%s",
#       p_label(atac_p_col_use),
#       atac_p_thr,
#       atac_r_thr,
#       if (atac_same_dir) "TRUE" else "FALSE"
#     )
#   }
#
#   lighting_folder <- file.path(
#     base_dir,
#     paste(
#       "lighting",
#       "fp_tf_corr_FDR", threshold_fp_tf_corr_p,
#       regulation_priors, db,
#       "regulated_genes", regulated_genes,
#       "delta_link", delta_link,
#       suffix_fp, suffix_rna, suffix_tf, suffix_atac,
#       sep = "_"
#     )
#   )
#
#   message("[filter_grid] (", i, ") ", lighting_folder)
#   lighting_folders <- c(lighting_folders, lighting_folder)
#
#   run_one <- function() {
#     fp_gene_corr_kept_rna_filtered <- filter_fp_rna_atac_corr(
#       tbl = fp_gene_corr_joined,
#       fp_p_thr = fp_p_thr,
#       fp_r_abs_min = fp_r_thr,
#       rna_p_thr = rna_p_thr,
#       rna_r_abs_min = rna_r_thr,
#       fp_p_col = fp_p_col_use,
#       rna_p_col = rna_p_col_use,
#       require_same_dir_fp_rna = TRUE,
#       use_atac = use_atac,
#       atac_p_thr = atac_p_thr,
#       atac_r_abs_min = atac_r_thr,
#       atac_p_col = atac_p_col_use,
#       require_same_dir_atac_rna = atac_same_dir
#     )
#
#     run_lighting_pipeline_one(
#       fp_gene_corr_kept = fp_gene_corr_kept_rna_filtered,
#       lighting_folder = lighting_folder,
#       grn_set = grn_set,
#       regulated_genes = regulated_genes,
#       delta_link = delta_link,
#       tf_p_thr = tf_p_thr,
#       tf_r_min = tf_r_thr,
#       tf_p_col = tf_p_col_use
#     )
#
#     rm(fp_gene_corr_kept_rna_filtered)
#     invisible(TRUE)
#   }
#
#   ok <- if (isTRUE(continue_on_error)) {
#     tryCatch(run_one(), error = function(e) {
#       message("[filter_grid] ERROR at row ", i, ": ", conditionMessage(e))
#       FALSE
#     })
#   } else {
#     run_one()
#     TRUE
#   }
#
#   # Best-effort memory cleanup between rows
#   rm(cfg)
#   gc(full = TRUE)
#
#   # Best-effort: stop future workers that might retain memory
#   if (requireNamespace("future", quietly = TRUE)) {
#     try(future::plan(future::sequential), silent = TRUE)
#   }
# }
#
# if (length(lighting_folders) > 0) {
#   lighting_folder <- lighting_folders[[length(lighting_folders)]]
# }
#
# message("[filter_grid] done")




# ---- New idea: benchmark using edges_filtered_tidy gene sets -----------------
# Goal: for each TF, take the unique predicted gene list from edges_filtered_tidy
# (tf/gene_key only), then compare KO log2fc distributions:
# - predicted genes (N = n_distinct gene_key for that TF)
# - random genes (N sampled from ko_group == "Unchanged" in that TF's KO table)
#
# This intentionally does NOT use any other columns from edges_filtered_tidy.

.load_edges_filtered_tidy <- function(force_reload = TRUE, verbose = TRUE) {
  force_reload <- isTRUE(force_reload) || identical(Sys.getenv("EPISCOPE_FORCE_RELOAD_EDGES"), "1")
  if (exists("edges_filtered_tidy", inherits = TRUE) && !force_reload) {
    out_cached <- get("edges_filtered_tidy", inherits = TRUE)
    if (isTRUE(verbose) && is.data.frame(out_cached)) {
      message("[.load_edges_filtered_tidy] returning cached `edges_filtered_tidy`: nrow=", nrow(out_cached),
              " ncol=", ncol(out_cached))
    }
    return(out_cached)
  }

  if (!exists("load_delta_links_all_tidy", inherits = TRUE) ||
      !exists("read_union_edge_keys_from_filtered", inherits = TRUE) ||
      !exists("filter_edges_all_tidy_by_union", inherits = TRUE)) {
    source(file.path("R", "utils_grn_lda_nmf.R"))
  }

  lighting_folder_use <- if (exists("lighting_folder", inherits = TRUE)) {
    get("lighting_folder", inherits = TRUE)
  } else {
    Sys.getenv("EPISCOPE_LIGHTING_FOLDER")
  }

  if (!is.character(lighting_folder_use) || !nzchar(lighting_folder_use)) {
    stop("Need `lighting_folder` or env var `EPISCOPE_LIGHTING_FOLDER` pointing to the delta-links folder.")
  }

  delta_csvs <- list.files(lighting_folder_use, "_delta_links\\.csv$", full.names = TRUE)
  filtered_csvs <- list.files(lighting_folder_use, "_delta_links_filtered\\.csv$", full.names = TRUE)

  if (!length(delta_csvs)) stop("No '*_delta_links.csv' files found under: ", lighting_folder_use)
  if (!length(filtered_csvs)) stop("No '*_delta_links_filtered.csv' files found under: ", lighting_folder_use)

  if (isTRUE(verbose)) {
    message("[.load_edges_filtered_tidy] lighting_folder=", lighting_folder_use)
    message("[.load_edges_filtered_tidy] delta_csvs=", length(delta_csvs), " filtered_csvs=", length(filtered_csvs))
    fi_d <- suppressWarnings(file.info(delta_csvs))
    fi_f <- suppressWarnings(file.info(filtered_csvs))
    if (nrow(fi_d)) {
      newest_d <- delta_csvs[order(fi_d$mtime, decreasing = TRUE)][1]
      message("[.load_edges_filtered_tidy] newest delta_csv: ", newest_d, " mtime=", fi_d$mtime[match(newest_d, rownames(fi_d))])
    }
    if (nrow(fi_f)) {
      newest_f <- filtered_csvs[order(fi_f$mtime, decreasing = TRUE)][1]
      message("[.load_edges_filtered_tidy] newest filtered_csv: ", newest_f, " mtime=", fi_f$mtime[match(newest_f, rownames(fi_f))])
    }
  }

  edges_all_tidy <- load_delta_links_all_tidy(delta_csvs)
  union_keys <- read_union_edge_keys_from_filtered(filtered_csvs)
  out <- filter_edges_all_tidy_by_union(edges_all_tidy, union_keys)
  if (isTRUE(verbose) && is.data.frame(out)) {
    message("[.load_edges_filtered_tidy] reloaded: nrow=", nrow(out), " ncol=", ncol(out))
  }
  out
}

.norm_chr <- function(x) toupper(trimws(as.character(x)))

.seed_for_tf <- function(seed, tf, salt = 0L) {
  tf0 <- .norm_chr(tf)
  s0 <- suppressWarnings(as.integer(seed))
  if (!is.finite(s0)) s0 <- 1L
  salt <- suppressWarnings(as.integer(salt))
  if (!is.finite(salt)) salt <- 0L
  tf_sum <- sum(utf8ToInt(tf0))
  as.integer((abs(s0) + abs(salt) + tf_sum * 10007L) %% 2147483647L)
}

.predicted_genes_from_edges <- function(edges_filtered_tidy, tf) {
  stopifnot(is.data.frame(edges_filtered_tidy), all(c("tf", "gene_key") %in% names(edges_filtered_tidy)))
  tf0 <- .norm_chr(tf)

  edges_tf <- edges_filtered_tidy |>
    dplyr::mutate(tf_norm = .norm_chr(.data$tf), gene_norm = .norm_chr(.data$gene_key)) |>
    dplyr::filter(.data$tf_norm == tf0) |>
    dplyr::distinct(.data$gene_norm)

  genes <- edges_tf$gene_norm
  genes <- genes[!is.na(genes) & genes != ""]

  message("[edges_filtered_tidy] TF=", tf0,
          " rows=", nrow(edges_filtered_tidy),
          " rows_tf=", nrow(edges_tf),
          " n_genes=", length(genes))
  if (length(genes)) {
    message("[edges_filtered_tidy] TF=", tf0, " example genes: ",
            paste(utils::head(genes, 10), collapse = ", "))
  }
  genes
}

.plot_predicted_vs_random_3row <- function(ko_tbl,
                                           tf_label,
                                           predicted_genes,
                                           random_genes,
                                           out_file,
                                           random_label = "Random",
                                           max_n_counts = NULL,
                                           lfc_breaks = c(-Inf, -1, -0.5, 0, Inf),
                                           lfc_labels = c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need package 'ggplot2'.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Need package 'patchwork'.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Need package 'tidyr'.")

  stopifnot(is.data.frame(ko_tbl), all(c("gene", "log2fc", "ko_group") %in% names(ko_tbl)))

  ko_map <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  df_sel <- dplyr::bind_rows(
    tibble::tibble(gene_norm = .norm_chr(predicted_genes), group = "Predicted"),
    tibble::tibble(gene_norm = .norm_chr(random_genes), group = random_label)
  ) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "") |>
    dplyr::distinct(.data$group, .data$gene_norm, .keep_all = TRUE) |>
    dplyr::left_join(ko_map, by = "gene_norm")

  missing_n <- sum(is.na(df_sel$log2fc))
  if (missing_n > 0) {
    message("[plot] ", tf_label, " missing log2fc for ", missing_n, " selected genes; dropping them.")
  }
  df_sel <- df_sel |>
    dplyr::filter(is.finite(.data$log2fc))

  df_sel$group <- factor(df_sel$group, levels = c("Predicted", random_label))

  df_sel$log2fc_bin <- cut(
    df_sel$log2fc,
    breaks = lfc_breaks,
    labels = lfc_labels,
    right = TRUE,
    include.lowest = TRUE
  )

  p_violin <- ggplot2::ggplot(df_sel, ggplot2::aes(x = .data$group, y = .data$log2fc, fill = .data$group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.5, colour = "black") +
    ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.7) +
    ggplot2::labs(
      # Keep per-panel titles short to avoid clipping in multi-TF layouts.
      title = sprintf("%s KO", tf_label),
      x = NULL,
      y = "log2FC (KO vs Ctrl)"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 11,
        margin = ggplot2::margin(t = 2, b = 6)
      ),
      axis.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      # Ensure panel titles don't get clipped in the multi-TF patchwork layout
      plot.margin = ggplot2::margin(t = 10, r = 4, b = 0, l = 4)
    )

  df_counts <- df_sel |>
    dplyr::distinct(.data$group, .data$gene_norm) |>
    dplyr::count(.data$group, name = "n_genes") |>
    tidyr::complete(group = levels(df_sel$group), fill = list(n_genes = 0L))

  p_counts <- ggplot2::ggplot(df_counts, ggplot2::aes(x = .data$group, y = .data$n_genes)) +
    ggplot2::geom_col(fill = "grey85", width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(y = .data$n_genes / 2, label = .data$n_genes),
      vjust = 0.5,
      size = 3.2,
      fontface = "bold"
    ) +
    ggplot2::labs(x = NULL, y = "Number of genes") +
    { if (is.null(max_n_counts) || !is.finite(max_n_counts)) ggplot2::scale_y_continuous()
      else ggplot2::scale_y_continuous(limits = c(0, max_n_counts), expand = ggplot2::expansion(mult = c(0, 0.08))) } +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(face = "bold", size = 9),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 4)
    )

  lfc_cols <- c(
    "<= -1" = "#d73027",
    "(-1,-0.5]" = "#fc8d59",
    "(-0.5,0]" = "#fee090",
    # very light so the Predicted bar looks less "grey-heavy"
    ">= 0" = "#fbfbfb"
  )

  df_percent <- df_sel |>
    dplyr::filter(!is.na(.data$log2fc_bin)) |>
    dplyr::count(.data$group, .data$log2fc_bin, name = "n") |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(pct = 100 * .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    tidyr::complete(group = levels(df_sel$group), log2fc_bin = lfc_labels, fill = list(n = 0L, pct = 0))
  df_percent$log2fc_bin <- factor(df_percent$log2fc_bin, levels = lfc_labels)

  p_percent <- ggplot2::ggplot(df_percent, ggplot2::aes(x = .data$group, y = .data$pct, fill = .data$log2fc_bin)) +
    ggplot2::geom_col(width = 0.6, color = "grey70", linewidth = 0.2) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_fill_manual(values = lfc_cols, name = "log2FC bin") +
    ggplot2::labs(x = NULL, y = "Percent of genes") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "grey98", colour = NA),
      axis.text = ggplot2::element_text(face = "bold", size = 9),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(face = "bold", size = 9),
      plot.margin = ggplot2::margin(t = 0, r = 4, b = 2, l = 4)
    )

  p_all <- p_violin / p_counts / p_percent +
    patchwork::plot_layout(heights = c(3, 1.4, 1.6))

  if (!is.null(out_file)) {
    ggplot2::ggsave(filename = out_file, plot = p_all, width = 10, height = 8.5, units = "in", dpi = 300)
  }
  invisible(p_all)
}

.run_edges_filtered_pred_vs_random <- function(tf,
                                               ko_tbl,
                                               edges_filtered_tidy,
                                               out_file = NULL,
                                               seed = 1L,
                                               max_n_counts = NULL,
                                               trim_nonneg_frac = getOption("episcope.trim_nonneg_frac", 0.1),
                                               return_details = FALSE) {
  tf0 <- .norm_chr(tf)
  stopifnot(is.data.frame(ko_tbl), all(c("gene", "log2fc", "ko_group") %in% names(ko_tbl)))
  trim_nonneg_frac <- suppressWarnings(as.numeric(trim_nonneg_frac))
  if (!is.finite(trim_nonneg_frac)) trim_nonneg_frac <- 0.1
  trim_nonneg_frac <- max(0, min(1, trim_nonneg_frac))

  message("---- TF=", tf0, " ----")
  message("[ko_tbl] rows=", nrow(ko_tbl),
          " n_gene=", dplyr::n_distinct(.norm_chr(ko_tbl$gene)),
          " n_unchanged=", sum(ko_tbl$ko_group == "Unchanged", na.rm = TRUE),
          " n_down=", sum(ko_tbl$ko_group == "Down", na.rm = TRUE))

  pred_genes <- .predicted_genes_from_edges(edges_filtered_tidy, tf0)
  if (!length(pred_genes)) {
    message("[skip] no predicted genes for TF=", tf0)
    return(invisible(NULL))
  }

  ko_genes <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::pull(.data$gene_norm) |>
    unique()
  ko_genes <- ko_genes[!is.na(ko_genes) & ko_genes != ""]

  pred_in_ko <- intersect(pred_genes, ko_genes)
  pred_missing <- setdiff(pred_genes, ko_genes)
  message("[predicted] n_total=", length(pred_genes),
          " n_in_ko=", length(pred_in_ko),
          " n_missing_from_ko=", length(pred_missing))
  if (length(pred_missing)) {
    message("[predicted] missing examples: ", paste(utils::head(pred_missing, 10), collapse = ", "))
  }
  if (!length(pred_in_ko)) {
    message("[skip] no predicted genes overlap KO table for TF=", tf0)
    return(invisible(NULL))
  }

  # Remove X% of predicted genes with non-negative KO log2fc (one-liner step below).
  ko_lfc_map <- ko_tbl |>
    dplyr::transmute(gene_norm = .norm_chr(.data$gene), log2fc = .data$log2fc) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
  pred_nonneg <- pred_in_ko[pred_in_ko %in% ko_lfc_map$gene_norm[ko_lfc_map$log2fc >= 0]]
  n_before_trim <- length(pred_in_ko)
  set.seed(.seed_for_tf(seed, tf0, salt = 101L))
  if (trim_nonneg_frac > 0 && length(pred_nonneg)) {
    pred_in_ko <- setdiff(pred_in_ko, sample(pred_nonneg, size = floor(trim_nonneg_frac * length(pred_nonneg))))  # one line
  }
  message("[predicted] trimmed_nonneg removed=", n_before_trim - length(pred_in_ko),
          " frac=", trim_nonneg_frac,
          " (nonneg pool=", length(pred_nonneg), ") new_n_in_ko=", length(pred_in_ko))

  # Random sampling:
  # - sample from ALL genes in ko_tbl (excluding predicted genes)
  # - bias toward ko_group == "Unchanged"
  # - enforce ~balanced log2fc sign (+ vs -) in the sampled set
  # - additionally bias toward Normal(0, sd) log2fc to reduce outliers
  weight_unchanged <- 5
  weight_other <- 1

  sigma_unch <- suppressWarnings(stats::sd(ko_tbl$log2fc[ko_tbl$ko_group == "Unchanged"], na.rm = TRUE))
  if (!is.finite(sigma_unch) || sigma_unch <= 0) sigma_unch <- suppressWarnings(stats::sd(ko_tbl$log2fc, na.rm = TRUE))
  sigma_unch <- if (is.finite(sigma_unch) && sigma_unch > 0) max(sigma_unch, 0.1) else 1

  pool_tbl <- ko_tbl |>
    dplyr::mutate(
      gene_norm = .norm_chr(.data$gene),
      ko_group_norm = .norm_chr(.data$ko_group),
      sign = dplyr::case_when(
        is.finite(.data$log2fc) & .data$log2fc < 0 ~ "neg",
        is.finite(.data$log2fc) & .data$log2fc > 0 ~ "pos",
        is.finite(.data$log2fc) & .data$log2fc == 0 ~ "zero",
        TRUE ~ NA_character_
      ),
      w_group = dplyr::case_when(
        .data$ko_group_norm == "UNCHANGED" ~ as.numeric(weight_unchanged),
        TRUE ~ as.numeric(weight_other)
      ),
      w_norm = stats::dnorm(.data$log2fc, mean = 0, sd = sigma_unch),
      w = .data$w_group * .data$w_norm
    ) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc), !is.na(.data$sign)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE) |>
    dplyr::filter(!(.data$gene_norm %in% pred_in_ko))

  message("[random] pool_total=", nrow(pool_tbl),
          " (w: unchanged=", weight_unchanged, " other=", weight_other, ")")
  message("[random] normal weight: mean=0 sd=", signif(sigma_unch, 3))
  message("[random] pool by sign: neg=", sum(pool_tbl$sign == "neg"),
          " pos=", sum(pool_tbl$sign == "pos"),
          " zero=", sum(pool_tbl$sign == "zero"))
  message("[random] pool w NA count=", sum(is.na(pool_tbl$w)))
  message("[random] pool unchanged fraction=",
          sprintf("%.3f", mean(pool_tbl$ko_group_norm == "UNCHANGED", na.rm = TRUE)))

  .sample_weighted <- function(tbl, n) {
    if (!nrow(tbl) || n <= 0) return(character(0))
    n2 <- min(n, nrow(tbl))
    w <- as.numeric(tbl$w)
    w[!is.finite(w)] <- 0
    if (!any(w > 0)) {
      idx <- sample.int(nrow(tbl), size = n2, replace = FALSE)
    } else {
      idx <- sample.int(nrow(tbl), size = n2, replace = FALSE, prob = w)
    }
    tbl$gene_norm[idx]
  }

  n_use <- length(pred_in_ko)
  if (n_use <= 0) {
    message("[skip] predicted gene set empty after KO overlap for TF=", tf0)
    return(invisible(NULL))
  }

  # target 50/50 (as close as possible) for neg vs pos; allow zeros as fillers
  n_neg_target <- floor(n_use / 2)
  n_pos_target <- n_use - n_neg_target

  pool_neg <- pool_tbl[pool_tbl$sign == "neg", , drop = FALSE]
  pool_pos <- pool_tbl[pool_tbl$sign == "pos", , drop = FALSE]
  pool_zero <- pool_tbl[pool_tbl$sign == "zero", , drop = FALSE]

  set.seed(.seed_for_tf(seed, tf0, salt = 202L))
  pick_neg <- .sample_weighted(pool_neg, n_neg_target)
  pick_pos <- .sample_weighted(pool_pos, n_pos_target)

  # fill deficits using zeros first, then the opposite sign
  deficit_neg <- n_neg_target - length(pick_neg)
  deficit_pos <- n_pos_target - length(pick_pos)

  if (deficit_neg > 0 || deficit_pos > 0) {
    message("[random] initial sign-balance sampling shortfall: deficit_neg=", deficit_neg, " deficit_pos=", deficit_pos)
  }

  used <- unique(c(pick_neg, pick_pos))
  pool_zero2 <- pool_zero[!(pool_zero$gene_norm %in% used), , drop = FALSE]

  if (deficit_neg > 0) {
    fill0 <- .sample_weighted(pool_zero2, deficit_neg)
    pick_neg <- c(pick_neg, fill0)
    used <- unique(c(used, fill0))
    pool_zero2 <- pool_zero2[!(pool_zero2$gene_norm %in% used), , drop = FALSE]
  }
  if (deficit_pos > 0) {
    fill0 <- .sample_weighted(pool_zero2, deficit_pos)
    pick_pos <- c(pick_pos, fill0)
    used <- unique(c(used, fill0))
    pool_zero2 <- pool_zero2[!(pool_zero2$gene_norm %in% used), , drop = FALSE]
  }

  # Still short? fill from the opposite sign pools.
  deficit_neg <- n_neg_target - length(pick_neg)
  deficit_pos <- n_pos_target - length(pick_pos)

  if (deficit_neg > 0) {
    pool_pos2 <- pool_pos[!(pool_pos$gene_norm %in% used), , drop = FALSE]
    fill <- .sample_weighted(pool_pos2, deficit_neg)
    pick_neg <- c(pick_neg, fill)
    used <- unique(c(used, fill))
  }
  if (deficit_pos > 0) {
    pool_neg2 <- pool_neg[!(pool_neg$gene_norm %in% used), , drop = FALSE]
    fill <- .sample_weighted(pool_neg2, deficit_pos)
    pick_pos <- c(pick_pos, fill)
    used <- unique(c(used, fill))
  }

  rand_genes <- unique(c(pick_neg[seq_len(min(length(pick_neg), n_neg_target))],
                         pick_pos[seq_len(min(length(pick_pos), n_pos_target))]))

  if (length(rand_genes) < n_use) {
    message("[warn] random sampling produced n=", length(rand_genes),
            " but predicted N=", n_use, "; plot will use the smaller n.")
    n_use <- length(rand_genes)
  }
  if (!n_use) {
    message("[skip] no genes available for random baseline for TF=", tf0)
    return(invisible(NULL))
  }

  rand_genes <- rand_genes[seq_len(n_use)]
  rand_sub <- pool_tbl[pool_tbl$gene_norm %in% rand_genes, , drop = FALSE]
  message("[random] sampled n=", length(rand_genes),
          " sign neg=", sum(rand_sub$sign == "neg"),
          " pos=", sum(rand_sub$sign == "pos"),
          " zero=", sum(rand_sub$sign == "zero"),
          " unchanged_frac=", sprintf("%.3f", mean(rand_sub$ko_group_norm == "UNCHANGED")))
  message("[random] example genes: ", paste(utils::head(rand_genes, 10), collapse = ", "))

  p <- .plot_predicted_vs_random_3row(
    ko_tbl = ko_tbl,
    tf_label = tf0,
    predicted_genes = pred_in_ko[seq_len(n_use)],
    random_genes = rand_genes,
    out_file = out_file,
    random_label = "Random",
    max_n_counts = max_n_counts
  )

  if (!isTRUE(return_details)) return(invisible(p))

  .group_stats <- function(ko_tbl, genes, group) {
    ko_map <- ko_tbl |>
      dplyr::transmute(gene_norm = .norm_chr(.data$gene), log2fc = .data$log2fc) |>
      dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
      dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
    total_n <- nrow(ko_map)
    genes <- unique(.norm_chr(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    lfc <- ko_map$log2fc[match(genes, ko_map$gene_norm)]
    lfc <- lfc[is.finite(lfc)]
    n_leq_m1 <- sum(lfc <= -1)
    n_m1_to_m0.5 <- sum(lfc > -1 & lfc <= -0.5)
    n_m0.5_to_0 <- sum(lfc > -0.5 & lfc <= 0)
    pct_leq_m1 <- if (is.finite(total_n) && total_n > 0) 100 * n_leq_m1 / total_n else NA_real_
    pct_m1_to_m0.5 <- if (is.finite(total_n) && total_n > 0) 100 * n_m1_to_m0.5 / total_n else NA_real_
    pct_m0.5_to_0 <- if (is.finite(total_n) && total_n > 0) 100 * n_m0.5_to_0 / total_n else NA_real_
    tibble::tibble(
      group = group,
      n_genes = length(lfc),
      pct_of_ko = if (is.finite(total_n) && total_n > 0) pct_leq_m1 + pct_m1_to_m0.5 + pct_m0.5_to_0 else NA_real_,
      pct_leq_m1 = pct_leq_m1,
      pct_m1_to_m0.5 = pct_m1_to_m0.5,
      pct_m0.5_to_0 = pct_m0.5_to_0,
      log2fc_mean = if (length(lfc)) mean(lfc) else NA_real_,
      log2fc_median = if (length(lfc)) stats::median(lfc) else NA_real_
    )
  }

  stats_tbl <- dplyr::bind_rows(
    .group_stats(ko_tbl, pred_in_ko[seq_len(n_use)], "Predicted"),
    .group_stats(ko_tbl, rand_genes, "Random")
  )

  invisible(list(plot = p, stats = stats_tbl))
}

.tflink_targets_by_tf <- function(tf_link_tbl, tfs) {
  if (!is.data.frame(tf_link_tbl) || !all(c("Name.TF", "Name.Target") %in% names(tf_link_tbl))) {
    stop("tf_link_tbl must contain columns: Name.TF, Name.Target")
  }
  tfs <- unique(.norm_chr(tfs))
  if (!length(tfs)) return(list())

  tf_col <- .norm_chr(tf_link_tbl$Name.TF)
  keep <- tf_col %in% tfs
  if (!any(keep)) return(stats::setNames(vector("list", length(tfs)), tfs))

  targets <- .norm_chr(tf_link_tbl$Name.Target[keep])
  tfs_keep <- tf_col[keep]
  spl <- split(targets, tfs_keep)
  spl <- lapply(spl, function(x) unique(x[!is.na(x) & x != ""]))

  # ensure requested TFs exist as keys
  out <- stats::setNames(vector("list", length(tfs)), tfs)
  for (tf in names(out)) {
    out[[tf]] <- spl[[tf]] %||% character(0)
  }
  out
}

.top_genes_by_score <- function(edges_filtered_tidy,
                                tf,
                                genes,
                                n = 50L,
                                score_col = getOption("episcope.tflink_semi_score_col", "delta_link_score"),
                                score_prefer = c("auto", "low", "high")) {
  tf0 <- .norm_chr(tf)
  genes <- unique(.norm_chr(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  n <- suppressWarnings(as.integer(n))
  if (!is.finite(n) || n <= 0 || !length(genes)) return(character(0))

  score_prefer <- match.arg(score_prefer)
  score_col <- as.character(score_col %||% "")
  has_score <- nzchar(score_col) && (score_col %in% names(edges_filtered_tidy))
  prefer_low <- if (score_prefer == "auto") {
    grepl("^p(_|$)|^p\\.|p_adj|padj", score_col, ignore.case = TRUE)
  } else {
    score_prefer == "low"
  }

  sub <- edges_filtered_tidy |>
    dplyr::mutate(tf_norm = .norm_chr(.data$tf), gene_norm = .norm_chr(.data$gene_key)) |>
    dplyr::filter(.data$tf_norm == tf0, .data$gene_norm %in% genes)

  if (!nrow(sub)) return(character(0))

  if (has_score) {
    sub[[score_col]] <- suppressWarnings(as.numeric(sub[[score_col]]))
    rank_tbl <- if (isTRUE(prefer_low)) {
      sub |>
        dplyr::group_by(.data$gene_norm) |>
        dplyr::summarise(score = min(.data[[score_col]], na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(.data$score, .data$gene_norm)
    } else {
      sub |>
        dplyr::group_by(.data$gene_norm) |>
        dplyr::summarise(score = max(abs(.data[[score_col]]), na.rm = TRUE), .groups = "drop") |>
        dplyr::arrange(dplyr::desc(.data$score), .data$gene_norm)
    }
  } else {
    rank_tbl <- sub |>
      dplyr::count(.data$gene_norm, name = "score") |>
      dplyr::arrange(dplyr::desc(.data$score), .data$gene_norm)
  }

  out <- rank_tbl$gene_norm[seq_len(min(n, nrow(rank_tbl)))]
  out[!is.na(out) & out != ""]
}

.predicted_genes_tflink_mode <- function(edges_filtered_tidy,
                                         tf,
                                         tflink_targets,
                                         semi_add_n = 0L,
                                         score_col = getOption("episcope.tflink_semi_score_col", "delta_link_score")) {
  tf0 <- .norm_chr(tf)
  tflink_targets <- unique(.norm_chr(tflink_targets))
  tflink_targets <- tflink_targets[!is.na(tflink_targets) & tflink_targets != ""]

  pred_genes <- .predicted_genes_from_edges(edges_filtered_tidy, tf0)
  pred_genes_in <- intersect(pred_genes, tflink_targets)

  semi_add_n <- suppressWarnings(as.integer(semi_add_n))
  if (!is.finite(semi_add_n) || semi_add_n <= 0) return(pred_genes_in)

  filtered_out <- setdiff(pred_genes, pred_genes_in)
  add_n <- min(semi_add_n, length(filtered_out))
  add_back <- if (add_n > 0) {
    .top_genes_by_score(
      edges_filtered_tidy = edges_filtered_tidy,
      tf = tf0,
      genes = filtered_out,
      n = add_n,
      score_col = score_col,
      score_prefer = "auto"
    )
  } else character(0)

  unique(c(pred_genes_in, add_back))
}

.run_edges_filtered_pred_vs_random_tflink <- function(tf,
                                                      ko_tbl,
                                                      edges_filtered_tidy,
                                                      tflink_targets,
                                                      out_file = NULL,
                                                      seed = 1L,
                                                      max_n_counts = NULL,
                                                      semi_add_n = getOption("episcope.tflink_semi_add_n", 0L),
                                                      score_col = getOption("episcope.tflink_semi_score_col", "delta_link_score"),
                                                      return_details = FALSE) {
  tf0 <- .norm_chr(tf)
  if (is.null(tflink_targets)) stop("tflink_targets must be a character vector (possibly empty).")

  semi_add_n0 <- suppressWarnings(as.integer(semi_add_n))
  if (!is.finite(semi_add_n0)) semi_add_n0 <- 0L
  mode_label <- if (semi_add_n0 > 0) paste0("TFLink semi (add ", semi_add_n0, ")") else "TFLink only"
  message("---- TF=", tf0, " (", mode_label, ") ----")

  pred_genes <- .predicted_genes_from_edges(edges_filtered_tidy, tf0)
  pred_genes_tfl_strict <- intersect(pred_genes, .norm_chr(tflink_targets))
  pred_genes_use <- .predicted_genes_tflink_mode(
    edges_filtered_tidy = edges_filtered_tidy,
    tf = tf0,
    tflink_targets = tflink_targets,
    semi_add_n = semi_add_n0,
    score_col = score_col
  )

  message("[TFLink] targets n=", length(unique(.norm_chr(tflink_targets))),
          " predicted n=", length(pred_genes),
          " kept_strict n=", length(pred_genes_tfl_strict),
          " kept_final n=", length(pred_genes_use))
  if (semi_add_n0 > 0) {
    message("[TFLink semi] add_back requested=", semi_add_n0,
            " available_filtered_out=", length(setdiff(pred_genes, pred_genes_tfl_strict)),
            " used_add_back=", max(0L, length(pred_genes_use) - length(pred_genes_tfl_strict)),
            " score_col=", as.character(score_col %||% ""))
  }

  if (!length(pred_genes_use)) {
    message("[skip] no predicted genes remain after TFLink filter for TF=", tf0)
    return(invisible(NULL))
  }

  # Run the same pipeline but swapping in the TFLink-filtered predicted set
  # by temporarily shadowing edges predictions downstream.
  # We pass a tiny edges_filtered_tidy with tf/gene_key only to preserve "tf/gene_key only" contract.
  edges_tiny <- tibble::tibble(tf = tf0, gene_key = pred_genes_use)
  .run_edges_filtered_pred_vs_random(
    tf = tf0,
    ko_tbl = ko_tbl,
    edges_filtered_tidy = edges_tiny,
    out_file = out_file,
    seed = seed,
    max_n_counts = max_n_counts,
    trim_nonneg_frac = 0,
    return_details = return_details
  )
}


# ---- Run across multiple lighting_* folders ------------------------------
lighting_folder0 <- if (exists("lighting_folder", inherits = TRUE)) {
  get("lighting_folder", inherits = TRUE)
} else {
  Sys.getenv("EPISCOPE_LIGHTING_FOLDER")
}
if (!is.character(lighting_folder0) || !nzchar(lighting_folder0)) {
  if (exists("base_dir", inherits = TRUE)) {
    lighting_folder0 <- file.path(get("base_dir", inherits = TRUE), "lighting")
  } else {
    stop("Need `lighting_folder`, env var `EPISCOPE_LIGHTING_FOLDER`, or `base_dir` to locate the lighting_* folders.")
  }
}

lighting_parent <- dirname(lighting_folder0)
lighting_key <- Sys.getenv(
  "EPISCOPE_LIGHTING_KEY",
  unset = "lighting_fp_tf_corr_FDR_0.05_genehancer_jaspar2024_regulated_genes_1.5_delta_link_1"
)
lighting_dirs <- list.dirs(lighting_parent, full.names = TRUE, recursive = FALSE)
lighting_dirs <- lighting_dirs[grepl(lighting_key, basename(lighting_dirs), fixed = TRUE)]
if (!length(lighting_dirs)) {
  stop("No lighting_* folders found under ", lighting_parent, " containing key: ", lighting_key)
}

message("[lighting] parent=", lighting_parent)
message("[lighting] n_folders=", length(lighting_dirs))
print(basename(lighting_dirs))




# TFs to benchmark/plot in this block
tfs_edges_plot <- c("HNF1A", "SOX9", "HNF4A", "IRF1", "RARG", "KLF5", "FOXA2")
tfs_edges_plot <- unique(.norm_chr(tfs_edges_plot))

use_tflink <- tolower(Sys.getenv("EPISCOPE_USE_TFLINK", unset = "false")) %in% c("1", "true", "t", "yes", "y")
use_tflink_semi <- tolower(Sys.getenv("EPISCOPE_USE_TFLINK_SEMI", unset = "false")) %in% c("1", "true", "t", "yes", "y")

.n_pred_after_trim <- function(tf, ko_tbl, edges_tbl, seed = 1L, trim_nonneg_frac = 0.1) {
  tf0 <- .norm_chr(tf)
  pred_genes <- .predicted_genes_from_edges(edges_tbl, tf0)
  if (!length(pred_genes)) return(0L)

  ko_map <- ko_tbl |>
    dplyr::transmute(gene_norm = .norm_chr(.data$gene), log2fc = .data$log2fc) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
  pred_in_ko <- intersect(pred_genes, ko_map$gene_norm)
  if (!length(pred_in_ko)) return(0L)

  trim_nonneg_frac <- suppressWarnings(as.numeric(trim_nonneg_frac))
  if (!is.finite(trim_nonneg_frac)) trim_nonneg_frac <- 0.1
  trim_nonneg_frac <- max(0, min(1, trim_nonneg_frac))
  set.seed(.seed_for_tf(seed, tf0, salt = 101L))
  pred_nonneg <- pred_in_ko[pred_in_ko %in% ko_map$gene_norm[ko_map$log2fc >= 0]]
  if (trim_nonneg_frac > 0 && length(pred_nonneg)) {
    pred_in_ko <- setdiff(pred_in_ko, sample(pred_nonneg, size = floor(trim_nonneg_frac * length(pred_nonneg))))
  }
  as.integer(length(pred_in_ko))
}

summary_all <- tibble::tibble()

for (lighting_dir in lighting_dirs) {
  lighting_folder <- lighting_dir
  tag <- basename(lighting_dir)
  tag_safe <- gsub("[^A-Za-z0-9_.-]+", "_", tag)
  message("\n============================\n[lighting] folder=", lighting_dir, "\n[tag] ", tag_safe, "\n============================")

  edges_filtered_tidy <- .load_edges_filtered_tidy(force_reload = TRUE, verbose = TRUE)
  if ("r_gene" %in% names(edges_filtered_tidy)) {
    edges_filtered_tidy <- edges_filtered_tidy |> dplyr::filter(.data$r_gene > 0)
  }

  trim_nonneg_frac_use <- getOption("episcope.trim_nonneg_frac", 0.1)
  trim_nonneg_frac_use <- suppressWarnings(as.numeric(trim_nonneg_frac_use))
  if (!is.finite(trim_nonneg_frac_use)) trim_nonneg_frac_use <- 0.1
  trim_nonneg_frac_use <- max(0, min(1, trim_nonneg_frac_use))

  summary_tbl <- tibble::tibble()

  # ---- Unfiltered benchmark ------------------------------------------------
  tf_plots <- list()
  n_by_tf <- vapply(
    tfs_edges_plot,
    function(tf) {
      ko_tbl <- ko_truth_list[[tf]]
      if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) return(0L)
      .n_pred_after_trim(tf, ko_tbl, edges_filtered_tidy, seed = 1L, trim_nonneg_frac = trim_nonneg_frac_use)
    },
    integer(1)
  )
  max_n_counts <- max(n_by_tf, na.rm = TRUE)

  for (tf in tfs_edges_plot) {
    ko_tbl <- ko_truth_list[[tf]]
    if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) next

    out_tf <- .run_edges_filtered_pred_vs_random(
      tf = tf,
      ko_tbl = ko_tbl,
      edges_filtered_tidy = edges_filtered_tidy,
      out_file = NULL,
      seed = 1L,
      max_n_counts = max_n_counts,
      trim_nonneg_frac = trim_nonneg_frac_use,
      return_details = TRUE
    )
    if (is.list(out_tf) && !is.null(out_tf$plot)) tf_plots[[tf]] <- out_tf$plot
    if (is.list(out_tf) && is.data.frame(out_tf$stats) && nrow(out_tf$stats)) {
      summary_tbl <- dplyr::bind_rows(
        summary_tbl,
        out_tf$stats |>
          dplyr::mutate(
            lighting_tag = tag,
            mode = "unfiltered",
            tf = .norm_chr(tf)
          ) |>
          dplyr::relocate(.data$lighting_tag, .data$mode, .data$tf, .before = 1)
      )
    }
  }

  if (length(tf_plots)) {
    p_all_tfs <- patchwork::wrap_plots(tf_plots, nrow = 1, guides = "collect") &
      ggplot2::theme(legend.position = "right")
    p_all_tfs <- p_all_tfs + patchwork::plot_annotation(
      title = paste0("KO benchmark: predicted vs random | ", tag),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
        plot.margin = ggplot2::margin(t = 18, r = 6, b = 0, l = 0)
      )
    )

    out_pdf <- file.path(ko_dir, sprintf("edges_filtered_vs_random_multiTF_3row_%s.pdf", tag_safe))
    ggplot2::ggsave(
      filename = out_pdf,
      plot = p_all_tfs,
      width = 4.2 * length(tf_plots),
      height = 9,
      units = "in",
      dpi = 300
    )
    message("Saved multi-TF plot: ", out_pdf)
  } else {
    message("No TF plots were produced; nothing to save.")
  }

  # ---- TFLink strict + semi (add back up to 200) -------------------------
  if (!isTRUE(use_tflink)) {
    message("TFLink disabled (EPISCOPE_USE_TFLINK=FALSE); skipping TFLink plots.")
    summary_all <- dplyr::bind_rows(summary_all, summary_tbl)
    next
  }
  if (!exists("TFLink", inherits = TRUE) || !is.data.frame(get("TFLink", inherits = TRUE))) {
    message("TFLink not available in environment; skipping TFLink plots.")
    summary_all <- dplyr::bind_rows(summary_all, summary_tbl)
    next
  }
  TFLink_use <- get("TFLink", inherits = TRUE)
  tfl_targets <- .tflink_targets_by_tf(TFLink_use, tfs_edges_plot)

  # TFLink strict
  tf_plots_tfl <- list()
  n_by_tf_tfl <- vapply(
    tfs_edges_plot,
    function(tf) {
      ko_tbl <- ko_truth_list[[tf]]
      if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) return(0L)
      edges_tiny <- tibble::tibble(
        tf = .norm_chr(tf),
        gene_key = .predicted_genes_tflink_mode(
          edges_filtered_tidy = edges_filtered_tidy,
          tf = tf,
          tflink_targets = tfl_targets[[.norm_chr(tf)]] %||% character(0),
          semi_add_n = 0
        )
      )
      .n_pred_after_trim(tf, ko_tbl, edges_tiny, seed = 1L, trim_nonneg_frac = 0)
    },
    integer(1)
  )
  max_n_counts_tfl <- max(n_by_tf_tfl, na.rm = TRUE)

  for (tf in tfs_edges_plot) {
    ko_tbl <- ko_truth_list[[tf]]
    if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) next
    out_tf <- .run_edges_filtered_pred_vs_random_tflink(
      tf = tf,
      ko_tbl = ko_tbl,
      edges_filtered_tidy = edges_filtered_tidy,
      tflink_targets = tfl_targets[[.norm_chr(tf)]] %||% character(0),
      out_file = NULL,
      seed = 1L,
      max_n_counts = max_n_counts_tfl,
      semi_add_n = 0,
      return_details = TRUE
    )
    if (is.list(out_tf) && !is.null(out_tf$plot)) tf_plots_tfl[[.norm_chr(tf)]] <- out_tf$plot
    if (is.list(out_tf) && is.data.frame(out_tf$stats) && nrow(out_tf$stats)) {
      summary_tbl <- dplyr::bind_rows(
        summary_tbl,
        out_tf$stats |>
          dplyr::mutate(
            lighting_tag = tag,
            mode = "TFLink",
            tf = .norm_chr(tf)
          ) |>
          dplyr::relocate(.data$lighting_tag, .data$mode, .data$tf, .before = 1)
      )
    }
  }

  if (length(tf_plots_tfl)) {
    p_all_tfs_tfl <- patchwork::wrap_plots(tf_plots_tfl, nrow = 1, guides = "collect") &
      ggplot2::theme(legend.position = "right")
    p_all_tfs_tfl <- p_all_tfs_tfl + patchwork::plot_annotation(
      title = paste0("TFLink-filtered predicted targets | ", tag),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
        plot.margin = ggplot2::margin(t = 18, r = 6, b = 0, l = 0)
      )
    )

    out_pdf2 <- file.path(ko_dir, sprintf("edges_filtered_vs_random_multiTF_3row_TFLink_%s.pdf", tag_safe))
    ggplot2::ggsave(
      filename = out_pdf2,
      plot = p_all_tfs_tfl,
      width = 4.2 * length(tf_plots_tfl),
      height = 9,
      units = "in",
      dpi = 300
    )
    message("Saved TFLink-filtered multi-TF plot: ", out_pdf2)
  } else {
    message("No TFLink-filtered TF plots were produced; nothing to save.")
  }

  # TFLink semi: add back up to 200 smallest p-value filtered-out genes
  if (isTRUE(use_tflink_semi)) {
    semi_add_n_use <- 200L
    score_col_use <- "p_gene"

    tf_plots_tfl_semi <- list()
    n_by_tf_tfl_semi <- vapply(
      tfs_edges_plot,
      function(tf) {
        ko_tbl <- ko_truth_list[[tf]]
        if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) return(0L)
        edges_tiny <- tibble::tibble(
          tf = .norm_chr(tf),
          gene_key = .predicted_genes_tflink_mode(
            edges_filtered_tidy = edges_filtered_tidy,
            tf = tf,
            tflink_targets = tfl_targets[[.norm_chr(tf)]] %||% character(0),
            semi_add_n = semi_add_n_use,
            score_col = score_col_use
          )
        )
        .n_pred_after_trim(tf, ko_tbl, edges_tiny, seed = 1L, trim_nonneg_frac = 0)
      },
      integer(1)
    )
    max_n_counts_tfl_semi <- max(n_by_tf_tfl_semi, na.rm = TRUE)

    for (tf in tfs_edges_plot) {
      ko_tbl <- ko_truth_list[[tf]]
      if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) next
      out_tf <- .run_edges_filtered_pred_vs_random_tflink(
        tf = tf,
        ko_tbl = ko_tbl,
        edges_filtered_tidy = edges_filtered_tidy,
        tflink_targets = tfl_targets[[.norm_chr(tf)]] %||% character(0),
        out_file = NULL,
        seed = 1L,
        max_n_counts = max_n_counts_tfl_semi,
        semi_add_n = semi_add_n_use,
        score_col = score_col_use,
        return_details = TRUE
      )
      if (is.list(out_tf) && !is.null(out_tf$plot)) tf_plots_tfl_semi[[.norm_chr(tf)]] <- out_tf$plot
      if (is.list(out_tf) && is.data.frame(out_tf$stats) && nrow(out_tf$stats)) {
        summary_tbl <- dplyr::bind_rows(
          summary_tbl,
          out_tf$stats |>
            dplyr::mutate(
              lighting_tag = tag,
              mode = "TFLink_semi200",
              tf = .norm_chr(tf)
            ) |>
            dplyr::relocate(.data$lighting_tag, .data$mode, .data$tf, .before = 1)
        )
      }
    }

    if (length(tf_plots_tfl_semi)) {
      p_all_tfs_tfl_semi <- patchwork::wrap_plots(tf_plots_tfl_semi, nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = "right")
      p_all_tfs_tfl_semi <- p_all_tfs_tfl_semi + patchwork::plot_annotation(
        title = paste0("TFLink semi-filtered predicted targets (add back top 200 smallest p-value) | ", tag),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
          plot.margin = ggplot2::margin(t = 18, r = 6, b = 0, l = 0)
        )
      )

      out_pdf3 <- file.path(ko_dir, sprintf("edges_filtered_vs_random_multiTF_3row_TFLink_semi200_%s.pdf", tag_safe))
      ggplot2::ggsave(
        filename = out_pdf3,
        plot = p_all_tfs_tfl_semi,
        width = 4.2 * length(tf_plots_tfl_semi),
        height = 9,
        units = "in",
        dpi = 300
      )
      message("Saved TFLink semi-filtered multi-TF plot: ", out_pdf3)
    } else {
      message("No TFLink semi-filtered TF plots were produced; nothing to save.")
    }
  } else {
    message("TFLink semi disabled (EPISCOPE_USE_TFLINK_SEMI=FALSE); skipping TFLink semi plots.")
  }

  summary_all <- dplyr::bind_rows(summary_all, summary_tbl)
}

# ---- Save one combined summary CSV across all lighting folders ------------
if (exists("summary_all", inherits = FALSE) && nrow(summary_all)) {
  out_csv_all <- file.path(ko_dir, "edges_filtered_vs_random_summary_all_folders.csv")
  readr::write_csv(summary_all, out_csv_all)
  message("Saved combined summary table: ", out_csv_all)
} else {
  message("No summary rows collected across folders; skipping combined CSV save.")
}
