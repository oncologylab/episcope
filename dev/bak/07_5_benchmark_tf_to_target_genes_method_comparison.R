library(ggplot2)
library(patchwork)
ko_dir  <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

## ------------------------------------------------------------------
## Robust TF–RNA correlation helper (all genes)
## ------------------------------------------------------------------
compute_tf_rna_corr <- function(tf,
                                strict_rna,
                                method = c("pearson", "spearman")) {
  method <- match.arg(method)

  if (!all(c("ensembl_gene_id", "HGNC") %in% names(strict_rna))) {
    stop("strict_rna must contain columns 'ensembl_gene_id' and 'HGNC'.")
  }

  sample_cols <- setdiff(names(strict_rna), c("ensembl_gene_id", "HGNC"))
  if (!length(sample_cols)) {
    stop("No sample columns found in strict_rna.")
  }

  rna_mat <- as.matrix(strict_rna[, sample_cols, drop = FALSE])

  tf_rows <- which(strict_rna$HGNC == tf)
  if (!length(tf_rows)) {
    stop("TF ", tf, " not found in strict_rna$HGNC.")
  }

  ## if multiple rows for TF, average them
  tf_expr <- colMeans(rna_mat[tf_rows, , drop = FALSE])
  tf_expr <- suppressWarnings(as.numeric(tf_expr))

  ## basic QC on TF expression
  if (sum(is.finite(tf_expr)) < 3L || stats::sd(tf_expr, na.rm = TRUE) == 0) {
    return(tibble::tibble(
      gene       = character(0),
      r_rna      = numeric(0),
      p_rna      = numeric(0),
      p_adj_rna  = numeric(0)
    ))
  }

  cor_vec <- apply(
    rna_mat,
    1L,
    function(x) {
      x_num <- suppressWarnings(as.numeric(x))
      ok    <- is.finite(x_num) & is.finite(tf_expr)
      n_ok  <- sum(ok)
      if (n_ok < 3L) return(NA_real_)
      suppressWarnings(stats::cor(x_num[ok], tf_expr[ok],
                                  use = "complete.obs",
                                  method = method))
    }
  )

  ## per-gene n for p-values
  n_vec <- apply(
    rna_mat,
    1L,
    function(x) {
      x_num <- suppressWarnings(as.numeric(x))
      sum(is.finite(x_num) & is.finite(tf_expr))
    }
  )

  p_vec <- vapply(
    seq_along(cor_vec),
    function(i) {
      r <- cor_vec[i]
      n <- n_vec[i]
      if (!is.finite(r) || n < 3L) return(NA_real_)
      if (abs(r) >= 1) return(0)
      tval <- r * sqrt((n - 2) / (1 - r^2))
      2 * stats::pt(-abs(tval), df = n - 2)
    },
    numeric(1L)
  )

  padj_vec <- stats::p.adjust(p_vec, method = "BH")

  tibble::tibble(
    gene       = strict_rna$HGNC,
    r_rna      = as.numeric(cor_vec),
    p_rna      = as.numeric(p_vec),
    p_adj_rna  = as.numeric(padj_vec)
  ) |>
    dplyr::filter(!is.na(gene), gene != "")
}

## ------------------------------------------------------------------
## KO (perturbation) method-summary plot: FP / FP+ATAC / ATAC / RNA + distance to TSS
## ------------------------------------------------------------------
plot_ko_method_summary <- function(tf,
                                   ko_truth_tbl,
                                   tf_tfbs,
                                   tf_gene_links_gh,
                                   atac_gene_links_gh,
                                   tf_gene_links_30kb,
                                   atac_gene_links_30kb,
                                   strict_rna,
                                   out_dir,
                                   rna_cor_method = c("pearson", "spearman"),
                                   perturb_label = NULL,
                                   ## NEW: optional list of additional window sources
                                   ## named like "Window10kb", "Window20kb", ...
                                   window_sources = NULL) {

  rna_cor_method <- match.arg(rna_cor_method)

  ## if user doesn't specify, default to "<TF> KO" (keeps old behaviour)
  if (is.null(perturb_label)) {
    perturb_label <- sprintf("%s KO", tf)
  }

  ## ----------------------------
  ## Gene annotation / TSS helper
  ## ----------------------------
  if (!exists("gene_annot_ref_hg38", inherits = TRUE)) {
    cli::cli_abort("Object {.val gene_annot_ref_hg38} not found in environment.")
  }

  gene_tss_ref <- gene_annot_ref_hg38 |>
    dplyr::mutate(is_ensg = startsWith(ensembl_gene_id, "ENSG")) |>
    dplyr::group_by(HGNC) |>
    dplyr::arrange(dplyr::desc(is_ensg), ensembl_gene_id, .by_group = TRUE) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(HGNC, chrom, tss)

  tf_row   <- gene_tss_ref[gene_tss_ref$HGNC == tf, , drop = FALSE]
  tf_chrom <- if (nrow(tf_row)) tf_row$chrom[1] else NA_character_
  tf_tss   <- if (nrow(tf_row)) tf_row$tss[1]   else NA_integer_

  ## peak center parser: "chr:start-end" -> center
  compute_peak_center_tbl <- function(peak_vec) {
    peak_vec <- unique(peak_vec)
    peak_vec <- peak_vec[!is.na(peak_vec)]
    if (!length(peak_vec)) {
      return(tibble::tibble(
        peak   = character(0),
        chrom  = character(0),
        center = integer(0)
      ))
    }
    parts  <- strsplit(peak_vec, ":", fixed = TRUE)
    chrom  <- vapply(parts, `[`, character(1), 1L)
    range  <- vapply(parts, `[`, character(1), 2L)
    se     <- strsplit(range, "-", fixed = TRUE)
    start  <- as.integer(vapply(se, `[`, character(1), 1L))
    end    <- as.integer(vapply(se, `[`, character(1), 2L))
    center <- floor((start + end) / 2)
    tibble::tibble(
      peak   = peak_vec,
      chrom  = chrom,
      center = center
    )
  }

  ## ----------------------------
  ## Pre-compute TF–RNA correlations ONCE (full universe)
  ## ----------------------------
  rna_all <- compute_tf_rna_corr(
    tf         = tf,
    strict_rna = strict_rna,
    method     = rna_cor_method
  )

  ## per-gene RNA corr from precomputed table + distance TF↔target TSS
  compute_tf_rna_corr_for_genes <- function(gene_universe) {
    if (!nrow(rna_all)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }

    genes_use <- intersect(gene_universe, rna_all$gene)
    if (!length(genes_use)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }

    sub <- rna_all[rna_all$gene %in% genes_use, , drop = FALSE]

    ## distance: TF TSS to target gene TSS (same chrom only)
    sub2 <- dplyr::left_join(
      sub,
      gene_tss_ref,
      by = c("gene" = "HGNC")
    )
    names(sub2)[names(sub2) == "chrom"] <- "gene_chrom"
    names(sub2)[names(sub2) == "tss"]   <- "gene_tss"

    dist_vec <- rep(NA_real_, nrow(sub2))
    if (!is.na(tf_chrom) && !is.na(tf_tss)) {
      same_chr <- !is.na(sub2$gene_chrom) & sub2$gene_chrom == tf_chrom &
        !is.na(sub2$gene_tss)
      dist_vec[same_chr] <- abs(sub2$gene_tss[same_chr] - tf_tss)
    }

    tibble::tibble(
      gene     = sub2$gene,
      p_adj    = sub2$p_adj_rna,
      r        = sub2$r_rna,
      dist_tss = dist_vec
    )
  }

  ## ----------------------------
  ## Helper: dedup per (peak,gene)
  ## ----------------------------
  dedup_fp_tbl <- function(tf_gene_links_tbl, tf_symbol) {
    tbl <- tf_gene_links_tbl
    if ("tfs" %in% names(tbl)) {
      tbl <- tbl[tbl$tfs == tf_symbol, , drop = FALSE]
    }
    tbl <- tbl[!is.na(tbl$p_adj_fp), , drop = FALSE]

    if (!nrow(tbl)) {
      return(tbl[0, , drop = FALSE])
    }

    dplyr::group_by(tbl, fp_peak, gene_key) |>
      dplyr::summarise(
        atac_peak = atac_peak[which.min(p_adj_fp)],
        r_fp      = r_fp[which.min(p_adj_fp)],
        p_adj_fp  = min(p_adj_fp, na.rm = TRUE),
        .groups   = "drop"
      )
  }

  dedup_atac_tbl <- function(atac_gene_links_tbl) {
    tbl <- atac_gene_links_tbl
    tbl <- tbl[!is.na(tbl$p_adj_atac), , drop = FALSE]

    if (!nrow(tbl)) {
      return(tbl[0, , drop = FALSE])
    }

    dplyr::group_by(tbl, atac_peak, gene_key) |>
      dplyr::summarise(
        r_atac     = r_atac[which.min(p_adj_atac)],
        p_adj_atac = min(p_adj_atac, na.rm = TRUE),
        .groups    = "drop"
      )
  }

  ## ----------------------------
  ## Helper: per-gene best p/r + distance from FP / ATAC
  ## ----------------------------
  per_gene_from_fp <- function(fp_tbl, genes_univ) {
    if (!nrow(fp_tbl)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }
    fp_tbl <- fp_tbl[fp_tbl$gene_key %in% genes_univ, , drop = FALSE]
    if (!nrow(fp_tbl)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }

    dplyr::group_by(fp_tbl, gene = gene_key) |>
      dplyr::summarise(
        p_adj    = min(p_adj_fp, na.rm = TRUE),
        r        = r_fp[which.min(p_adj_fp)],
        dist_tss = dist_tss_fp[which.min(p_adj_fp)],
        .groups  = "drop"
      )
  }

  per_gene_from_atac <- function(atac_tbl, genes_univ) {
    if (!nrow(atac_tbl)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }
    atac_tbl <- atac_tbl[atac_tbl$gene_key %in% genes_univ, , drop = FALSE]
    if (!nrow(atac_tbl)) {
      return(tibble::tibble(
        gene     = character(0),
        p_adj    = numeric(0),
        r        = numeric(0),
        dist_tss = numeric(0)
      ))
    }

    dplyr::group_by(atac_tbl, gene = gene_key) |>
      dplyr::summarise(
        p_adj    = min(p_adj_atac, na.rm = TRUE),
        r        = r_atac[which.min(p_adj_atac)],
        dist_tss = dist_tss_atac[which.min(p_adj_atac)],
        .groups  = "drop"
      )
  }

  ## ----------------------------
  ## Build per-gene data for each source
  ## ----------------------------
  build_per_gene_for_source <- function(source_name,
                                        tf_gene_links_tbl,
                                        atac_gene_links_tbl) {

    fp_dedup   <- dedup_fp_tbl(tf_gene_links_tbl, tf)
    atac_dedup <- dedup_atac_tbl(atac_gene_links_tbl)

    ## add FP peak centers + gene TSS
    if (nrow(fp_dedup)) {
      fp_centers <- compute_peak_center_tbl(fp_dedup$fp_peak)
      fp_dedup <- dplyr::left_join(
        fp_dedup,
        fp_centers,
        by = c("fp_peak" = "peak")
      )
      names(fp_dedup)[names(fp_dedup) == "chrom"]  <- "fp_chrom"
      names(fp_dedup)[names(fp_dedup) == "center"] <- "fp_center"

      fp_dedup <- dplyr::left_join(
        fp_dedup,
        gene_tss_ref,
        by = c("gene_key" = "HGNC")
      )
      names(fp_dedup)[names(fp_dedup) == "chrom"] <- "gene_chrom"
      names(fp_dedup)[names(fp_dedup) == "tss"]   <- "gene_tss"

      fp_dedup$dist_tss_fp <- NA_real_
      ok_fp <- !is.na(fp_dedup$fp_center) & !is.na(fp_dedup$gene_tss) &
        !is.na(fp_dedup$fp_chrom) & !is.na(fp_dedup$gene_chrom) &
        fp_dedup$fp_chrom == fp_dedup$gene_chrom
      fp_dedup$dist_tss_fp[ok_fp] <-
        abs(fp_dedup$fp_center[ok_fp] - fp_dedup$gene_tss[ok_fp])
    } else {
      fp_dedup$dist_tss_fp <- numeric(0)
    }

    ## add ATAC peak centers + gene TSS
    if (nrow(atac_dedup)) {
      atac_centers <- compute_peak_center_tbl(atac_dedup$atac_peak)
      atac_dedup <- dplyr::left_join(
        atac_dedup,
        atac_centers,
        by = c("atac_peak" = "peak")
      )
      names(atac_dedup)[names(atac_dedup) == "chrom"]  <- "atac_chrom"
      names(atac_dedup)[names(atac_dedup) == "center"] <- "atac_center"

      atac_dedup <- dplyr::left_join(
        atac_dedup,
        gene_tss_ref,
        by = c("gene_key" = "HGNC")
      )
      names(atac_dedup)[names(atac_dedup) == "chrom"] <- "gene_chrom"
      names(atac_dedup)[names(atac_dedup) == "tss"]   <- "gene_tss"

      atac_dedup$dist_tss_atac <- NA_real_
      ok_atac <- !is.na(atac_dedup$atac_center) & !is.na(atac_dedup$gene_tss) &
        !is.na(atac_dedup$atac_chrom) & !is.na(atac_dedup$gene_chrom) &
        atac_dedup$atac_chrom == atac_dedup$gene_chrom
      atac_dedup$dist_tss_atac[ok_atac] <-
        abs(atac_dedup$atac_center[ok_atac] - atac_dedup$gene_tss[ok_atac])
    } else {
      atac_dedup$dist_tss_atac <- numeric(0)
    }

    genes_univ <- unique(c(fp_dedup$gene_key, atac_dedup$gene_key))
    genes_univ <- intersect(genes_univ, ko_truth_tbl$gene)
    genes_univ <- genes_univ[!is.na(genes_univ) & nzchar(genes_univ)]

    ## FP only
    per_fp <- per_gene_from_fp(fp_dedup, genes_univ)
    per_fp$method_label  <- "FP"

    ## FP + ATAC with filter on ATAC stats:
    ## keep only (atac_peak, gene_key) where |r_atac| > 0.3 & p_adj_atac < 0.05
    if (nrow(atac_dedup)) {
      atac_keep <- atac_dedup[
        is.finite(atac_dedup$r_atac) &
          is.finite(atac_dedup$p_adj_atac) &
          abs(atac_dedup$r_atac) > 0.3 &
          atac_dedup$p_adj_atac < 0.05,
        ,
        drop = FALSE
      ]

      if (nrow(atac_keep)) {
        fp_atac <- dplyr::inner_join(
          fp_dedup,
          atac_keep[, c("atac_peak", "gene_key"), drop = FALSE],
          by = c("atac_peak", "gene_key")
        )
      } else {
        fp_atac <- fp_dedup[0, , drop = FALSE]
      }
    } else {
      fp_atac <- fp_dedup[0, , drop = FALSE]
    }

    per_fp_atac <- per_gene_from_fp(fp_atac, genes_univ)
    per_fp_atac$method_label <- "FP+ATAC"

    ## ATAC alone (no extra filter)
    per_atac <- per_gene_from_atac(atac_dedup, genes_univ)
    per_atac$method_label <- "ATAC"

    ## RNA-only on same gene universe (from precomputed rna_all)
    per_rna <- compute_tf_rna_corr_for_genes(genes_univ)
    per_rna$method_label <- "RNA"

    ## combine
    res <- dplyr::bind_rows(per_fp, per_fp_atac, per_atac, per_rna)
    if (!nrow(res)) {
      return(res)
    }

    res <- dplyr::inner_join(
      res,
      ko_truth_tbl[, c("gene", "ko_group", "log2fc"), drop = FALSE],
      by = "gene"
    )

    res$data_source <- source_name
    res$neg_log10_p <- -log10(res$p_adj + 1e-300)

    res
  }

  ## GeneHancer
  per_gene_gh <- build_per_gene_for_source(
    source_name         = "GeneHancer",
    tf_gene_links_tbl   = tf_gene_links_gh,
    atac_gene_links_tbl = atac_gene_links_gh
  )

  ## canonical 30kb
  per_gene_30kb <- build_per_gene_for_source(
    source_name         = "Window30kb",
    tf_gene_links_tbl   = tf_gene_links_30kb,
    atac_gene_links_tbl = atac_gene_links_30kb
  )

  per_gene_list <- list(per_gene_gh, per_gene_30kb)

  ## additional window-based sources (e.g. 10kb,20kb,...,100kb)
  if (!is.null(window_sources) && length(window_sources)) {
    for (nm in names(window_sources)) {
      ## skip if we already have this source name
      if (nm %in% c("GeneHancer", "Window30kb")) {
        next
      }
      src <- window_sources[[nm]]
      per_gene_list[[length(per_gene_list) + 1L]] <- build_per_gene_for_source(
        source_name         = nm,
        tf_gene_links_tbl   = src$tf_gene_links,
        atac_gene_links_tbl = src$atac_gene_links
      )
    }
  }

  per_gene_df <- dplyr::bind_rows(per_gene_list)

  if (!nrow(per_gene_df)) {
    cli::cli_warn("No per-gene data available for plotting.")
    return(invisible(list(per_gene = per_gene_df)))
  }

  ## keep ONLY KO groups of interest; drop NA/Other/etc
  keep_groups <- c("Down", "Unchanged", "Up")
  per_gene_df <- per_gene_df[per_gene_df$ko_group %in% keep_groups, , drop = FALSE]

  per_gene_df$ko_group <- factor(per_gene_df$ko_group,
                                 levels = keep_groups)

  ## dynamic data_source order: GeneHancer first, then Window10kb,20kb,..., then others
  ds_raw    <- unique(per_gene_df$data_source)
  win_mask  <- grepl("^Window", ds_raw)
  win_names <- ds_raw[win_mask]
  win_nums  <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", win_names)))
  win_ord   <- order(win_nums, na.last = TRUE)
  ordered_windows <- win_names[win_ord]

  ds_levels <- character(0L)
  if ("GeneHancer" %in% ds_raw) {
    ds_levels <- c(ds_levels, "GeneHancer")
  }
  if (length(ordered_windows)) {
    ds_levels <- c(ds_levels, ordered_windows)
  }
  other_names <- setdiff(ds_raw, c("GeneHancer", ordered_windows))
  if (length(other_names)) {
    ds_levels <- c(ds_levels, sort(other_names))
  }

  per_gene_df$data_source <- factor(per_gene_df$data_source,
                                    levels = ds_levels)

  per_gene_df$method_label <- factor(
    per_gene_df$method_label,
    levels = c("FP", "FP+ATAC", "ATAC", "RNA")
  )

  ## ----- distance to TSS: log10 transform for plotting -----
  per_gene_df$log10_dist_tss <- ifelse(
    is.finite(per_gene_df$dist_tss) & per_gene_df$dist_tss >= 0,
    log10(per_gene_df$dist_tss + 1),
    NA_real_
  )

  ## preview distance tibble
  dist_tbl <- per_gene_df[, c("gene", "data_source", "method_label", "ko_group", "dist_tss"), drop = FALSE]
  message("Preview of distance-to-TSS tibble (first 20 rows):")
  print(utils::head(dist_tbl, 20))

  ## ----------------------------
  ## For p-value plot: trim extreme RNA outliers ONLY
  ## ----------------------------
  pval_df  <- per_gene_df
  rna_mask <- !is.na(pval_df$neg_log10_p) & pval_df$method_label == "RNA"

  if (any(rna_mask)) {
    cap <- stats::quantile(pval_df$neg_log10_p[rna_mask],
                           probs = 0.995, na.rm = TRUE)
    pval_df <- pval_df[!(rna_mask & pval_df$neg_log10_p > cap), , drop = FALSE]
  }

  ## Small helper: Wilcoxon p for Down vs Unchanged in a subset
  wilcox_du <- function(df, value_col) {
    dsub <- df[df$ko_group %in% c("Down", "Unchanged") &
                 is.finite(df[[value_col]]), , drop = FALSE]
    if (!nrow(dsub)) return(NA_real_)
    groups <- unique(dsub$ko_group)
    if (!all(c("Down", "Unchanged") %in% groups)) return(NA_real_)
    x <- dsub[dsub$ko_group == "Down", value_col, drop = TRUE]
    y <- dsub[dsub$ko_group == "Unchanged", value_col, drop = TRUE]
    if (length(x) < 3L || length(y) < 3L) return(NA_real_)
    suppressWarnings(stats::wilcox.test(x, y, alternative = "two.sided")$p.value)
  }

  ## --- stats for p-value panel ---
  stats_pval <- pval_df |>
    dplyr::ungroup() |>
    dplyr::group_by(data_source, method_label) |>
    dplyr::summarise(
      p_value = wilcox_du(dplyr::cur_data_all(), "neg_log10_p"),
      y_pos   = max(neg_log10_p, na.rm = TRUE) * 1.05,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = ifelse(is.na(p_value), "p = NA", paste0("p = ", signif(p_value, 2)))
    )

  ## --- stats for correlation panel ---
  stats_cor <- per_gene_df |>
    dplyr::ungroup() |>
    dplyr::group_by(data_source, method_label) |>
    dplyr::summarise(
      p_value = wilcox_du(dplyr::cur_data_all(), "r"),
      y_pos   = max(r, na.rm = TRUE) * 1.05,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = ifelse(is.na(p_value), "p = NA", paste0("p = ", signif(p_value, 2)))
    )

  ## ----------------------------
  ## PLOT 1: p-values (boxplots)
  ## ----------------------------
  p_pval <- ggplot2::ggplot(
    pval_df,
    ggplot2::aes(x = method_label, y = neg_log10_p, fill = ko_group)
  ) +
    ggplot2::geom_boxplot(
      position     = ggplot2::position_dodge(width = 0.7),
      width        = 0.6,
      outlier.size = 0.3,
      alpha        = 0.8
    ) +
    ggplot2::geom_text(
      data        = stats_pval,
      ggplot2::aes(x = method_label, y = y_pos, label = label),
      inherit.aes = FALSE,
      size        = 3
    ) +
    ggplot2::facet_wrap(~ data_source, nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c(Down = "#4daf4a", Unchanged = "grey60", Up = "#e41a1c"),
      name   = "KO group"
    ) +
    ggplot2::labs(
      title = sprintf("%s: per–gene p–value distributions", perturb_label),
      x     = "Method",
      y     = "-log10(p_adj)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )

  ## ----------------------------
  ## PLOT 2: correlations (boxplots)
  ## ----------------------------
  p_cor <- ggplot2::ggplot(
    per_gene_df,
    ggplot2::aes(x = method_label, y = r, fill = ko_group)
  ) +
    ggplot2::geom_boxplot(
      position     = ggplot2::position_dodge(width = 0.7),
      width        = 0.6,
      outlier.size = 0.3,
      alpha        = 0.8
    ) +
    ggplot2::geom_text(
      data        = stats_cor,
      ggplot2::aes(x = method_label, y = y_pos, label = label),
      inherit.aes = FALSE,
      size        = 3
    ) +
    ggplot2::facet_wrap(~ data_source, nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c(Down = "#4daf4a", Unchanged = "grey60", Up = "#e41a1c"),
      name   = "KO group"
    ) +
    ggplot2::labs(
      title = sprintf("%s: per–gene correlation distributions", perturb_label),
      x     = "Method",
      y     = "Correlation (r)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )

  ## ----------------------------
  ## PLOT 3: distance to TSS (boxplots of log10 distance)
  ## ----------------------------
  p_dist <- ggplot2::ggplot(
    per_gene_df,
    ggplot2::aes(x = method_label, y = log10_dist_tss, fill = ko_group)
  ) +
    ggplot2::geom_boxplot(
      position     = ggplot2::position_dodge(width = 0.7),
      width        = 0.6,
      outlier.size = 0.3,
      alpha        = 0.8
    ) +
    ggplot2::facet_wrap(~ data_source, nrow = 1) +
    ggplot2::scale_fill_manual(
      values = c(Down = "#4daf4a", Unchanged = "grey60", Up = "#e41a1c"),
      name   = "KO group"
    ) +
    ggplot2::labs(
      title = sprintf("%s: distance to TSS distributions", perturb_label),
      x     = "Method",
      y     = "log10(distance to TSS + 1 bp)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )

  combined <- p_pval / p_cor / p_dist +
    patchwork::plot_layout(heights = c(1, 1, 1)) +
    patchwork::plot_annotation(
      caption = "P-values in panels: two-sided Wilcoxon rank–sum test comparing Down vs Unchanged genes for each method and data source (for p-values and correlations)."
    )

  outfile <- file.path(
    out_dir,
    sprintf("%s_KO_per_gene_pval_corr_distTSS_boxplot.pdf", tf)
  )
  ## ---- NEW: scale width by number of facet columns -----------------
  n_facets    <- length(levels(per_gene_df$data_source))
  base_facets <- 2L  # originally: GeneHancer + Window30kb
  base_width  <- 8   # inches
  plot_width  <- base_width * (n_facets / base_facets)
  # optional: avoid too small width if n_facets < 2
  plot_width  <- max(base_width, plot_width)
  ## ------------------------------------------------------------------
  ggplot2::ggsave(
    filename = outfile,
    plot     = combined,
    width    = plot_width,
    height   = 11,
    units    = "in",
    dpi      = 300
  )

  message("Per-gene distribution plot saved to: ", outfile)

  invisible(list(
    per_gene = per_gene_df,
    plot     = combined,
    outfile  = outfile
  ))
}

## ------------------------------------------------------------------
## Load link tables (GeneHancer + multiple windows)
## ------------------------------------------------------------------

## helper: load arbitrary window-based link tables (10kb–100kb etc.)
base_dir_window <- "/data/homes/yl814/episcope_test/nutrient_stress"

## GeneHancer-based
tf_gene_links_gh   <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_jaspar2024.csv"))
atac_gene_links_gh <- readr::read_csv(file.path(base_dir_window, "atac_gene_corr_full_jaspar2024.csv"))


load_window_sources <- function(base_dir,
                                windows = c("10kb","20kb","30kb","40kb","50kb",
                                            "60kb","70kb","80kb","90kb","100kb")) {
  res <- list()
  for (w in windows) {
    tf_file   <- file.path(base_dir, sprintf("fp_gene_corr_full_%s_jaspar2024.csv", w))
    atac_file <- file.path(base_dir, sprintf("atac_gene_corr_full_%s_jaspar2024.csv", w))

    if (!file.exists(tf_file) || !file.exists(atac_file)) {
      cli::cli_warn(sprintf("Missing window files for %s; skipping.", w))
      next
    }

    res[[paste0("Window", w)]] <- list(
      tf_gene_links   = readr::read_csv(tf_file,   show_col_types = FALSE),
      atac_gene_links = readr::read_csv(atac_file, show_col_types = FALSE)
    )
  }
  res
}

## choose which windows you actually want to plot
window_sizes   <- c("10kb","20kb","30kb","40kb","50kb",
                    "60kb","70kb","80kb","90kb","100kb")
window_sources <- load_window_sources(base_dir_window, windows = window_sizes)

## keep 30kb single objects so old arguments still work
if (!"Window30kb" %in% names(window_sources)) {
  stop("Window30kb not found in window_sources – check window_sizes / files.")
}
tf_gene_links_30kb   <- window_sources[["Window30kb"]]$tf_gene_links
atac_gene_links_30kb <- window_sources[["Window30kb"]]$atac_gene_links

## Examples
res_hnf1a_summary <- plot_ko_method_summary(
  tf                   = "HNF1A",
  ko_truth_tbl         = ko_truth_HNF1A,
  tf_tfbs              = tf_tfbs_HNF1A,
  tf_gene_links_gh     = tf_gene_links_gh,
  atac_gene_links_gh   = atac_gene_links_gh,
  tf_gene_links_30kb   = tf_gene_links_30kb,
  atac_gene_links_30kb = atac_gene_links_30kb,
  strict_rna           = strict_rna,
  out_dir              = ko_dir,
  rna_cor_method       = "pearson",
  window_sources       = window_sources
)

res_hnf4a_summary <- plot_ko_method_summary(
  tf                   = "HNF4A",
  ko_truth_tbl         = ko_truth_HNF4A,
  tf_tfbs              = tf_tfbs_HNF4A,
  tf_gene_links_gh     = tf_gene_links_gh,
  atac_gene_links_gh   = atac_gene_links_gh,
  tf_gene_links_30kb   = tf_gene_links_30kb,
  atac_gene_links_30kb = atac_gene_links_30kb,
  strict_rna           = strict_rna,
  out_dir              = ko_dir,
  rna_cor_method       = "pearson",
  window_sources       = window_sources
)

## thresholds matching your HNF1A/HNF4A definition
lfc_strong <- 1        # |log2FC| > 1 for regulated
padj_sig   <- 0.05     # significant
lfc_unch   <- 0.25     # |log2FC| < 0.25 for unchanged
padj_unch  <- 0.5      # padj > 0.5 for unchanged

make_ko_truth_from_db <- function(tf_symbol,
                                  tf_perturb_db,
                                  subtype_filter = NULL,
                                  lfc_strong = 1,
                                  padj_sig   = 0.05,
                                  lfc_unch   = 0.25,
                                  padj_unch  = 0.5) {

  tbl <- get_tf_perturbation_tbl(tf_perturb_db, tf_symbol)

  if (!is.null(subtype_filter)) {
    tbl <- tbl[tbl$sub_type %in% subtype_filter, , drop = FALSE]
  }

  tbl |>
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") |>
    dplyr::transmute(
      gene    = gene_symbol,
      log2fc  = log2fc,
      ko_group = dplyr::case_when(
        !is.na(p_adj) & p_adj <= padj_sig & log2fc <= -lfc_strong ~ "Down",
        !is.na(p_adj) & p_adj <= padj_sig & log2fc >=  lfc_strong ~ "Up",
        !is.na(p_adj) & p_adj >= padj_unch & abs(log2fc) < lfc_unch ~ "Unchanged",
        TRUE ~ NA_character_  # intermediate / ignore
      )
    )
}
tf_tfbs_HNF1A  <- get_tf_tfbs("HNF1A",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_HNF4A  <- get_tf_tfbs("HNF4A",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")

tf_tfbs_IRF1  <- get_tf_tfbs("IRF1",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_RARG  <- get_tf_tfbs("RARG",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_SOX9  <- get_tf_tfbs("SOX9",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_KLF5  <- get_tf_tfbs("KLF5",  tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_FOXA2 <- get_tf_tfbs("FOXA2", tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")

## Build KO truth tables with those thresholds
ko_truth_HNF1A  <- make_ko_truth_from_db("HNF1A",  tf_perturb_db)

ko_truth_IRF1  <- make_ko_truth_from_db("IRF1",  tf_perturb_db)
ko_truth_RARG  <- make_ko_truth_from_db("RARG",  tf_perturb_db)
ko_truth_SOX9  <- make_ko_truth_from_db("SOX9",  tf_perturb_db)
ko_truth_KLF5  <- make_ko_truth_from_db("KLF5",  tf_perturb_db)

## FOXA2: only PANC1 subtype
ko_truth_FOXA2 <- make_ko_truth_from_db(
  tf_symbol      = "FOXA2",
  tf_perturb_db  = tf_perturb_db,
  subtype_filter = "PANC1"
)

tfs <- c("IRF1", "RARG", "SOX9", "KLF5", "FOXA2")

ko_truth_list <- list(
  IRF1  = ko_truth_IRF1,
  RARG  = ko_truth_RARG,
  SOX9  = ko_truth_SOX9,
  KLF5  = ko_truth_KLF5,
  FOXA2 = ko_truth_FOXA2
)

tf_tfbs_list <- list(
  IRF1  = tf_tfbs_IRF1,
  RARG  = tf_tfbs_RARG,
  SOX9  = tf_tfbs_SOX9,
  KLF5  = tf_tfbs_KLF5,
  FOXA2 = tf_tfbs_FOXA2
)

res_list <- lapply(tfs, function(tf_symbol) {
  plot_ko_method_summary(
    tf                   = tf_symbol,
    ko_truth_tbl         = ko_truth_list[[tf_symbol]],
    tf_tfbs              = tf_tfbs_list[[tf_symbol]],
    tf_gene_links_gh     = tf_gene_links_gh,
    atac_gene_links_gh   = atac_gene_links_gh,
    tf_gene_links_30kb   = tf_gene_links_30kb,
    atac_gene_links_30kb = atac_gene_links_30kb,
    strict_rna           = strict_rna,
    out_dir              = ko_dir,
    rna_cor_method       = "pearson",
    window_sources       = window_sources
  )
})
names(res_list) <- tfs



## ------------------------------------------------------------------
## HiC-based sources: load and re-use plot_ko_method_summary
## ------------------------------------------------------------------

## Load HiC-based link tables we previously saved:
##   fp_gene_corr_full_<gh_std_name>.csv
##   atac_gene_corr_full_<gh_std_name>.csv
load_hic_sources <- function(base_dir, gh_std_names) {
  res <- list()
  for (nm in gh_std_names) {
    tf_file   <- file.path(base_dir, sprintf("fp_gene_corr_full_%s.csv", nm))
    atac_file <- file.path(base_dir, sprintf("atac_gene_corr_full_%s.csv", nm))

    if (!file.exists(tf_file) || !file.exists(atac_file)) {
      cli::cli_warn(sprintf("Missing HiC files for %s; skipping.", nm))
      next
    }

    res[[nm]] <- list(
      tf_gene_links   = readr::read_csv(tf_file,   show_col_types = FALSE),
      atac_gene_links = readr::read_csv(atac_file, show_col_types = FALSE)
    )
  }
  res
}

## Discover HiC datasets: prefer names(gh_std_hic_list) if available,
## otherwise infer from filenames fp_gene_corr_full_*.csv
hic_names <- if (exists("gh_std_hic_list", inherits = TRUE)) {
  names(gh_std_hic_list)
} else {
  hic_fp_files <- list.files(
    base_dir_window,
    pattern = "^fp_gene_corr_full_.*\\.csv$",
    full.names = FALSE
  )
  sub("^fp_gene_corr_full_(.*)\\.csv$", "\\1", hic_fp_files)
}

hic_sources <- load_hic_sources(base_dir_window, hic_names)

## Empty 30kb tables so Window30kb contributes no facet in HiC plots
tf_gene_links_30kb_empty   <- tf_gene_links_gh[0, , drop = FALSE]
atac_gene_links_30kb_empty <- atac_gene_links_gh[0, , drop = FALSE]

## HiC-specific KO summaries:
## facets = GeneHancer (column 1) + each HiC cell line; no 30–100 kb windows
res_list_hic <- lapply(tfs, function(tf_symbol) {
  plot_ko_method_summary(
    tf                   = tf_symbol,
    ko_truth_tbl         = ko_truth_list[[tf_symbol]],
    tf_tfbs              = tf_tfbs_list[[tf_symbol]],
    tf_gene_links_gh     = tf_gene_links_gh,
    atac_gene_links_gh   = atac_gene_links_gh,
    tf_gene_links_30kb   = tf_gene_links_30kb_empty,
    atac_gene_links_30kb = atac_gene_links_30kb_empty,
    strict_rna           = strict_rna,
    out_dir              = ko_dir,
    rna_cor_method       = "pearson",
    window_sources       = hic_sources
  )
})
names(res_list_hic) <- paste0(tfs, "_HiC")


## ------------------------------------------------------------------
## Add HiC versions for HNF1A and HNF4A
## ------------------------------------------------------------------

res_hnf1a_hic_summary <- plot_ko_method_summary(
  tf                   = "HNF1A",
  ko_truth_tbl         = ko_truth_HNF1A,
  tf_tfbs              = tf_tfbs_HNF1A,
  tf_gene_links_gh     = tf_gene_links_gh,
  atac_gene_links_gh   = atac_gene_links_gh,
  tf_gene_links_30kb   = tf_gene_links_30kb_empty,
  atac_gene_links_30kb = atac_gene_links_30kb_empty,
  strict_rna           = strict_rna,
  out_dir              = ko_dir,
  rna_cor_method       = "pearson",
  window_sources       = hic_sources
)

res_hnf4a_hic_summary <- plot_ko_method_summary(
  tf                   = "HNF4A",
  ko_truth_tbl         = ko_truth_HNF4A,
  tf_tfbs              = tf_tfbs_HNF4A,
  tf_gene_links_gh     = tf_gene_links_gh,
  atac_gene_links_gh   = atac_gene_links_gh,
  tf_gene_links_30kb   = tf_gene_links_30kb_empty,
  atac_gene_links_30kb = atac_gene_links_30kb_empty,
  strict_rna           = strict_rna,
  out_dir              = ko_dir,
  rna_cor_method       = "pearson",
  window_sources       = hic_sources
)
