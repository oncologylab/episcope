suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(plyr)
})

base_dir_window <- "/data/homes/yl814/episcope_test/nutrient_stress"
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

## GeneHancer-based
tf_gene_links_gh   <- readr::read_csv(file.path(base_dir_window, "connect_tfs_to_target_genes", "fp_gene_corr_full_jaspar2024.csv"))
atac_gene_links_gh <- readr::read_csv(file.path(base_dir_window, "connect_tfs_to_target_genes", "atac_gene_corr_full_jaspar2024.csv"))
## Cell-line specific GeneHancer-based GRN sets
tf_gene_links_gh_Panc1  <- readr::read_csv(file.path(base_dir_window, "grn_per_cellline", "fp_gene_corr_full_grn_set_Panc1.csv"))
tf_gene_links_gh_AsPC1  <- readr::read_csv(file.path(base_dir_window, "grn_per_cellline", "fp_gene_corr_full_grn_set_AsPC1.csv"))
tf_gene_links_gh_HPAFII <- readr::read_csv(file.path(base_dir_window, "grn_per_cellline", "fp_gene_corr_full_grn_set_HPAFII.csv"))

## ==================================================================
## 1) TF → cell / lineage mapping
## ==================================================================

tf_cell_lineage <- tibble::tribble(
  ~tf,     ~cell,     ~lineage,
  "HNF1A", "AsPC1",   "inter",
  "HNF4A", "Mayo5289","classical",
  "IRF1",  "CFPAC1",  "classical",
  "SOX9",  "PANC1",   "basal",
  "KLF5",  "AsPC1",   "inter",
  "FOXA2", "PANC1",   "basal",
  "RARG",  "PANC1",   "basal"
)

get_tf_cell_lineage <- function(tf) {
  hit <- tf_cell_lineage[tf_cell_lineage$tf == tf, , drop = FALSE]
  if (nrow(hit) == 0L) {
    warning("No tf_cell_lineage entry for TF = ", tf,
            "; using NAcell / NAlineage in filenames.")
    return(list(cell = "NAcell", lineage = "NAlineage"))
  }
  list(cell = hit$cell[1], lineage = hit$lineage[1])
}

## ==================================================================
## 2) Plot helpers (split violin, p-value formatting)
## ==================================================================

geom_split_violin <- function (mapping = NULL,
                               data = NULL,
                               stat = "ydensity",
                               position = "identity", ...,
                               draw_quantiles = NULL,
                               trim = TRUE,
                               scale = "area",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE) {

  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
      data <- transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )
      grp <- data[1, "group"]

      newdata <- plyr::arrange(
        transform(
          data,
          x = if (grp %% 2 == 1) xminv else xmaxv
        ),
        if (grp %% 2 == 1) y else -y
      )

      newdata <- rbind(
        newdata[1, ],
        newdata,
        newdata[nrow(newdata), ],
        newdata[1, ]
      )
      newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

      if (length(draw_quantiles) > 0 && !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)),
                           setdiff(names(data), c("x", "y")),
                           drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)

        ggplot2:::ggname(
          "geom_split_violin",
          grid::grobTree(
            ggplot2::GeomPolygon$draw_panel(newdata, ...),
            quantile_grob
          )
        )
      } else {
        ggplot2:::ggname(
          "geom_split_violin",
          ggplot2::GeomPolygon$draw_panel(newdata, ...)
        )
      }
    }
  )

  ggplot2::layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomSplitViolin,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      trim           = trim,
      scale          = scale,
      draw_quantiles = draw_quantiles,
      na.rm          = na.rm,
      ...
    )
  )
}

format_p_short <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("N.S.")
    if (x >= 0.05) return("N.S.")
    if (x >= 0.001) {
      as.character(signif(x, 2))
    } else if (x == 0) {
      "0"
    } else {
      e <- floor(log10(x))
      sprintf("1e%d", e)
    }
  }, character(1))
}

## ==================================================================
## 3) RNA helpers: map gene keys and compute TF–RNA correlation
## ==================================================================

map_gene_keys_to_rna <- function(keys, strict_rna) {
  idx_hgnc <- match(keys, strict_rna$HGNC)
  idx_ens  <- match(keys, strict_rna$ensembl_gene_id)
  row_idx  <- ifelse(!is.na(idx_hgnc), idx_hgnc, idx_ens)
  keep     <- !is.na(row_idx)

  tibble(
    gene_key = keys[keep],
    row_idx  = row_idx[keep]
  )
}

compute_rna_corr_for_tf <- function(tf,
                                    genes,
                                    strict_rna,
                                    cor_method = c("pearson", "spearman")) {
  cor_method <- match.arg(cor_method)

  sample_cols <- setdiff(names(strict_rna),
                         c("ensembl_gene_id", "HGNC"))
  if (!length(sample_cols)) {
    stop("strict_rna has no sample columns.")
  }

  rna_expr_mat <- as.matrix(strict_rna[, sample_cols, drop = FALSE])

  idx_tf <- which(strict_rna$HGNC == tf)
  if (!length(idx_tf)) {
    idx_tf <- which(strict_rna$ensembl_gene_id == tf)
  }
  if (!length(idx_tf)) {
    stop(sprintf("TF '%s' not found in strict_rna$HGNC or ensembl_gene_id.", tf))
  }

  idx_tf  <- idx_tf[1L]
  tf_expr <- as.numeric(rna_expr_mat[idx_tf, ])

  map_univ <- map_gene_keys_to_rna(genes, strict_rna)
  if (!nrow(map_univ)) {
    stop("None of the requested genes are present in strict_rna.")
  }

  corr_list <- lapply(seq_len(nrow(map_univ)), function(i) {
    g      <- map_univ$gene_key[i]
    row_i  <- map_univ$row_idx[i]
    g_expr <- as.numeric(rna_expr_mat[row_i, ])

    ok <- is.finite(tf_expr) & is.finite(g_expr)
    n  <- sum(ok)
    if (n < 3L) {
      return(list(
        gene_key = g,
        n_rna    = n,
        r_rna    = NA_real_,
        p_rna    = NA_real_
      ))
    }

    r_val <- suppressWarnings(
      stats::cor(tf_expr, g_expr, use = "complete.obs", method = cor_method)
    )

    if (is.na(r_val)) {
      p_val <- NA_real_
    } else if (abs(r_val) >= 1) {
      p_val <- 0
    } else {
      tt <- r_val * sqrt((n - 2) / (1 - r_val^2))
      p_val <- 2 * stats::pt(-abs(tt), df = n - 2)
    }

    list(
      gene_key = g,
      n_rna    = n,
      r_rna    = as.numeric(r_val),
      p_rna    = as.numeric(p_val)
    )
  })

  dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )
}

## ==================================================================
## 4) Build KO summary table for combined filters
## ==================================================================

build_ko_summary_table_for_tf <- function(tf,
                                          ko_tbl,
                                          tf_gene_links,
                                          strict_rna,
                                          tf_tfbs,
                                          TFLink,
                                          mode = c("all", "canonical"),
                                          r_cut_main = 0.3,
                                          p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                                          cor_method = c("pearson", "spearman"),
                                          fp_p_col = c("p_adj_fp", "p_fp"),
                                          verbose = TRUE) {

  mode       <- match.arg(mode)
  cor_method <- match.arg(cor_method)
  fp_p_col   <- match.arg(fp_p_col)

  if (verbose) {
    message("Building KO summary table for TF = ", tf, " (mode = ", mode, ")")
  }

  tf_meta <- get_tf_cell_lineage(tf)
  tf_cell <- tf_meta$cell
  tf_line <- tf_meta$lineage

  if (!all(c("Name.TF", "Name.Target") %in% names(TFLink))) {
    stop("TFLink object must contain columns 'Name.TF' and 'Name.Target'.")
  }

  tflink_targets <- unique(TFLink$Name.Target[TFLink$Name.TF == tf])
  tflink_targets <- tflink_targets[!is.na(tflink_targets) & nzchar(tflink_targets)]

  tf_links_all <- tf_gene_links
  if (mode == "canonical") {
    if (!"tfs" %in% names(tf_gene_links)) {
      stop("mode = 'canonical' requires a 'tfs' column in tf_gene_links.")
    }
    tf_links_all <- tf_links_all[tf_links_all$tfs == tf, , drop = FALSE]
  }

  tf_links_mode <- tf_links_all[tf_links_all$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links_mode) == 0L) {
    stop("No tf_gene_links rows for TF ", tf, " and supplied tf_tfbs.")
  }

  if (!"gene_key" %in% names(tf_links_mode)) {
    stop("tf_gene_links must have a 'gene_key' column.")
  }
  if (!"gene" %in% names(ko_tbl) || !"log2FC" %in% names(ko_tbl)) {
    stop("ko_tbl must have columns 'gene' and 'log2FC'.")
  }
  if (!"r_fp" %in% names(tf_links_mode) || !fp_p_col %in% names(tf_links_mode)) {
    stop("tf_gene_links must contain columns 'r_fp' and '", fp_p_col, "'.")
  }

  genes_gh0 <- sort(unique(tf_links_mode$gene_key))
  genes_gh  <- intersect(genes_gh0, ko_tbl$gene)
  if (!length(genes_gh)) {
    stop("No overlapping genes between KO table and tf_gene_links for TF ", tf, ".")
  }

  rna_corr_full <- compute_rna_corr_for_tf(
    tf          = tf,
    genes       = genes_gh,
    strict_rna  = strict_rna,
    cor_method  = cor_method
  )

  genes_rna_ok <- rna_corr_full$gene_key[!is.na(rna_corr_full$r_rna)]
  genes_base   <- intersect(genes_gh, genes_rna_ok)
  if (!length(genes_base)) {
    stop("After RNA correlation, no genes remain in GeneHancer × KO for TF ", tf, ".")
  }

  ko_sub        <- ko_tbl[ko_tbl$gene %in% genes_base, , drop = FALSE]
  tf_links_mode <- tf_links_mode[tf_links_mode$gene_key %in% genes_base, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_base, , drop = FALSE]

  ## FP “significant” genes (used for RNA + FP sig universe) ----------
  links_fp  <- tf_links_mode
  pvals_fp  <- links_fp[[fp_p_col]]
  keep_sig  <- !is.na(links_fp$r_fp) &
    !is.na(pvals_fp) &
    links_fp$r_fp > r_cut_main &
    pvals_fp < 0.05

  fp_sig_genes <- unique(links_fp$gene_key[keep_sig])

  rna_sig_genes <- rna_corr_full %>%
    dplyr::filter(
      !is.na(r_rna),
      !is.na(p_adj_rna),
      r_rna > r_cut_main,
      p_adj_rna < 0.05
    ) %>%
    dplyr::pull(gene_key) %>%
    unique()

  ## Base universes incl. "both" filters ------------------------------

  # FP side
  base_FP_only    <- genes_base
  base_FP_TFLink  <- intersect(genes_base, tflink_targets)
  base_FP_RNAfil  <- intersect(genes_base, rna_sig_genes)
  base_FP_both    <- intersect(base_FP_TFLink, base_FP_RNAfil)

  # RNA side
  base_RNA_GH     <- genes_base                      # RNA + GeneHancer
  base_RNA_TFLink <- intersect(genes_base, tflink_targets)
  base_RNA_FPfil  <- intersect(genes_base, fp_sig_genes)
  base_RNA_both   <- intersect(base_RNA_TFLink, base_RNA_FPfil)

  if (verbose) {
    message("Universe sizes for TF ", tf, " (mode ", mode, "):")
    message("  GH base (GeneHancer)            : ", length(genes_base))
    message("  FP only (GH base)               : ", length(base_FP_only))
    message("  FP + TFLink                      : ", length(base_FP_TFLink))
    message("  FP + RNA sig                     : ", length(base_FP_RNAfil))
    message("  RNA + GH (GeneHancer)           : ", length(base_RNA_GH))
    message("  RNA + TFLink                     : ", length(base_RNA_TFLink))
    message("  RNA + FP sig                     : ", length(base_RNA_FPfil))
  }

  get_fp_sets <- function(base_genes, p_cut) {
    if (!length(base_genes)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    links <- tf_links_mode[tf_links_mode$gene_key %in% base_genes, , drop = FALSE]
    if (!nrow(links)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    pvals <- links[[fp_p_col]]

    pass_rows <- !is.na(links$r_fp) &
      !is.na(pvals) &
      links$r_fp > r_cut_main &
      pvals < p_cut

    genes_pass <- sort(unique(links$gene_key[pass_rows]))
    genes_all  <- sort(unique(links$gene_key))
    genes_fail <- setdiff(genes_all, genes_pass)

    list(passing = genes_pass, non_passing = genes_fail)
  }

  get_rna_sets <- function(base_genes, p_cut) {
    if (!length(base_genes)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    g_sub <- rna_corr_full[rna_corr_full$gene_key %in% base_genes, , drop = FALSE]
    if (!nrow(g_sub)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    pass_rows <- !is.na(g_sub$r_rna) &
      !is.na(g_sub$p_adj_rna) &
      g_sub$r_rna > r_cut_main &
      g_sub$p_adj_rna < p_cut

    genes_pass <- sort(g_sub$gene_key[pass_rows])
    genes_all  <- sort(g_sub$gene_key)
    genes_fail <- setdiff(genes_all, genes_pass)

    list(passing = genes_pass, non_passing = genes_fail)
  }

  ## Methods exactly as requested -------------------------------------

  method_defs <- tibble::tibble(
    method_id   = c("FP_TFLink",
                    "FP_RNAfiltered",
                    "FP_both",
                    "RNA_TFLink",
                    "RNA_FPfiltered",
                    "RNA_both"),
    method_type = c("FP", "FP", "FP", "RNA", "RNA", "RNA"),
    method_lab  = factor(
      c(
        "FP (filtered by TFLink)",
        "FP (filtered by RNA r>0.3, p<0.05)",
        "FP (filtered by TFLink & RNA)",
        "RNA (filtered by TFLink)",
        "RNA (filtered by FP r>0.3, p<0.05)",
        "RNA (filtered by TFLink & FP)"
      ),
      levels = c(
        "FP (filtered by TFLink)",
        "FP (filtered by RNA r>0.3, p<0.05)",
        "FP (filtered by TFLink & RNA)",
        "RNA (filtered by TFLink)",
        "RNA (filtered by FP r>0.3, p<0.05)",
        "RNA (filtered by TFLink & FP)"
      )
    )
  )

  get_base_universe <- function(method_id) {
    switch(
      method_id,
      FP_TFLink      = base_FP_TFLink,
      FP_RNAfiltered = base_FP_RNAfil,
      FP_both        = base_FP_both,
      RNA_TFLink     = base_RNA_TFLink,
      RNA_FPfiltered = base_RNA_FPfil,
      RNA_both       = base_RNA_both,
      character(0)
    )
  }

  ## Build long-format table ------------------------------------------

  res_list <- list()
  idx      <- 1L

  for (i in seq_len(nrow(method_defs))) {
    mid   <- method_defs$method_id[i]
    mtype <- method_defs$method_type[i]
    mlab  <- as.character(method_defs$method_lab[i])

    base_genes <- get_base_universe(mid)
    if (!length(base_genes)) {
      if (verbose) message("Method ", mlab, ": base universe empty (skipping).")
      next
    }

    for (pc in p_cuts) {
      sets <- if (mtype == "FP") {
        get_fp_sets(base_genes, p_cut = pc)
      } else {
        get_rna_sets(base_genes, p_cut = pc)
      }

      genes_pred <- intersect(sets$passing, ko_sub$gene)
      genes_non  <- intersect(sets$non_passing, ko_sub$gene)

      if (!length(genes_pred) && !length(genes_non)) next

      ko_pred <- ko_sub[ko_sub$gene %in% genes_pred, , drop = FALSE]
      ko_non  <- ko_sub[ko_sub$gene %in% genes_non,  , drop = FALSE]

      if (!nrow(ko_pred) && !nrow(ko_non)) next

      res_list[[idx]] <- dplyr::bind_rows(
        tibble::tibble(
          tf         = tf,
          cell       = tf_cell,
          lineage    = tf_line,
          mode       = mode,
          method_id  = mid,
          method_lab = mlab,
          p_cut      = pc,
          group      = "Predicted",
          gene       = ko_pred$gene,
          log2FC     = ko_pred$log2FC
        ),
        tibble::tibble(
          tf         = tf,
          cell       = tf_cell,
          lineage    = tf_line,
          mode       = mode,
          method_id  = mid,
          method_lab = mlab,
          p_cut      = pc,
          group      = "Non-predicted",
          gene       = ko_non$gene,
          log2FC     = ko_non$log2FC
        )
      )
      idx <- idx + 1L
    }
  }

  if (!length(res_list)) {
    stop("No data generated for TF ", tf, " (mode ", mode, ").")
  }

  out <- dplyr::bind_rows(res_list)

  p_levels <- sort(unique(out$p_cut))
  out$p_lab <- factor(
    out$p_cut,
    levels = p_levels,
    labels = scales::label_scientific(digits = 1)(p_levels)
  )

  out$group <- factor(out$group, levels = c("Predicted", "Non-predicted"))

  out
}


## ==================================================================
## 5) Plotter: 4-row layout (violin + three bar summaries)
## ==================================================================

plot_ko_summary_violin_for_tf <- function(tf,
                                          ko_tbl,
                                          tf_gene_links,
                                          strict_rna,
                                          tf_tfbs,
                                          TFLink,
                                          mode = c("all", "canonical"),
                                          r_cut_main = 0.3,
                                          p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                                          cor_method = c("pearson", "spearman"),
                                          ko_type = "KO",
                                          out_dir = ".",
                                          verbose = TRUE,
                                          filename_suffix = "",
                                          fp_p_col = c("p_adj_fp", "p_fp")) {

  mode       <- match.arg(mode)
  cor_method <- match.arg(cor_method)
  fp_p_col   <- match.arg(fp_p_col)

  df <- build_ko_summary_table_for_tf(
    tf            = tf,
    ko_tbl        = ko_tbl,
    tf_gene_links = tf_gene_links,
    strict_rna    = strict_rna,
    tf_tfbs       = tf_tfbs,
    TFLink        = TFLink,
    mode          = mode,
    r_cut_main    = r_cut_main,
    p_cuts        = p_cuts,
    cor_method    = cor_method,
    fp_p_col      = fp_p_col,
    verbose       = verbose
  )

  tf_meta <- get_tf_cell_lineage(tf)
  tf_cell <- tf_meta$cell
  tf_line <- tf_meta$lineage

  ## 5A) Violin + p-values --------------------------------------------

  p_test_tbl <- df %>%
    dplyr::group_by(method_lab, p_lab) %>%
    dplyr::summarise(
      p_val = {
        sub  <- dplyr::cur_data_all()
        sub  <- sub[!is.na(sub$log2FC), , drop = FALSE]
        gtab <- table(sub$group)
        if (length(gtab) != 2L || any(gtab == 0L)) {
          NA_real_
        } else {
          stats::wilcox.test(log2FC ~ group, data = sub)$p.value
        }
      },
      y = max(log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = format_p_short(p_val),
      y     = y + 0.05 * diff(range(df$log2FC, na.rm = TRUE))
    )

  fill_vals <- c(
    Predicted       = "#1b9e77",
    `Non-predicted` = "#d95f02"
  )

  p_violin <- ggplot(
    df,
    aes(x = p_lab, y = log2FC, fill = group)
  ) +
    geom_split_violin(alpha = 0.4, trim = FALSE, colour = "black") +
    geom_boxplot(
      width          = 0.2,
      alpha          = 0.6,
      show.legend    = FALSE,
      outlier.size   = 0.1,
      outlier.stroke = 0,
      outlier.alpha  = 0.3
    ) +
    geom_text(
      data        = p_test_tbl,
      aes(x = p_lab, y = y, label = label),
      inherit.aes = FALSE,
      size        = 2.6,
      vjust       = 0,
      fontface    = "bold"
    ) +
    facet_grid(cols = vars(method_lab)) +
    scale_fill_manual(values = fill_vals, name = "Group") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = sprintf("log2FC (%s %s vs Ctrl)", tf, ko_type)) +
    ggtitle(sprintf(
      "%s %s in %s (%s) - %s mode\nR > %.1f across p-value cutoffs",
      tf, ko_type, tf_cell, tf_line, mode, r_cut_main
    )) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      legend.text      = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.title.y     = element_text(face = "bold"),
      axis.text.y      = element_text(face = "bold")
    )

  ## 5B) Bar: median |log2FC| -----------------------------------------

  bar_median_df <- df %>%
    dplyr::group_by(method_lab, p_lab, group) %>%
    dplyr::summarise(
      median_abs_log2FC = median(abs(log2FC), na.rm = TRUE),
      .groups           = "drop"
    )

  p_bar_median <- ggplot(
    bar_median_df,
    aes(x = p_lab, y = median_abs_log2FC, fill = group)
  ) +
    geom_col(
      position = position_dodge(width = 0.7),
      width    = 0.6,
      colour   = "black",
      linewidth = 0.2
    ) +
    facet_grid(cols = vars(method_lab)) +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(
      name   = "Median |log2FC| (Predicted vs Non-predicted)",
      expand = expansion(mult = c(0, 0.1))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(face = "bold", size = 8)
    )

  ## 5C) Bar: number of predicted genes -------------------------------

  n_pred_df <- df %>%
    dplyr::filter(group == "Predicted") %>%
    dplyr::group_by(method_lab, p_lab) %>%
    dplyr::summarise(
      n_pred = dplyr::n_distinct(gene),
      .groups = "drop"
    )

  p_bar_n <- ggplot(
    n_pred_df,
    aes(x = p_lab, y = n_pred)
  ) +
    geom_col(fill = "grey60") +
    geom_text(
      aes(label = n_pred),
      vjust    = -0.2,
      size     = 2.4,
      fontface = "bold"
    ) +
    facet_grid(cols = vars(method_lab)) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(
      name   = "Number of predicted target genes",
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(face = "bold", size = 7)
    )

  ## 5D) Bar: percent of log2FC bins per group ------------------------

  df_bins <- df %>%
    dplyr::mutate(
      log2FC_bin = dplyr::case_when(
        is.na(log2FC)                ~ NA_character_,
        log2FC <= -1                 ~ "<= -1",
        log2FC > -1 & log2FC <= -0.5 ~ "(-1,-0.5]",
        log2FC > -0.5 & log2FC < 0   ~ "(-0.5,0)",
        log2FC >= 0                  ~ ">= 0"
      )
    ) %>%
    dplyr::filter(!is.na(log2FC_bin))

  bin_df <- df_bins %>%
    dplyr::group_by(method_lab, p_lab, group, log2FC_bin) %>%
    dplyr::summarise(
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::group_by(method_lab, p_lab, group) %>%
    dplyr::mutate(
      pct = 100 * n / sum(n)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      x_lab = factor(
        paste0(as.character(p_lab), "\n", as.character(group)),
        levels = unique(paste0(as.character(p_lab), "\n", as.character(group)))
      )
    )

  log2fc_levels <- c(">= 0", "(-0.5,0)", "(-1,-0.5]", "<= -1")
  bin_df$log2FC_bin <- factor(bin_df$log2FC_bin, levels = log2fc_levels)

  log2fc_fill <- c(
    ">= 0"       = "#80b1d3",
    "(-0.5,0)"   = "#fdb462",
    "(-1,-0.5]"  = "#fb8072",
    "<= -1"      = "#e41a1c"
  )

  p_bar_pct <- ggplot(
    bin_df,
    aes(x = x_lab, y = pct, fill = log2FC_bin)
  ) +
    geom_col(colour = "black", linewidth = 0.2) +
    facet_grid(cols = vars(method_lab)) +
    scale_fill_manual(values = log2fc_fill, name = "log2FC bin") +
    scale_x_discrete(name = "Adjusted p-value cutoff × Group") +
    scale_y_continuous(
      name   = "Percent of genes within group",
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_text(face = "bold", size = 6, angle = 45, hjust = 1),
      axis.text.y      = element_text(face = "bold", size = 7),
      legend.title     = element_text(face = "bold"),
      legend.text      = element_text(face = "bold")
    )

  ## 5E) Combine rows and save ----------------------------------------

  combined <- p_violin / p_bar_median / p_bar_n / p_bar_pct +
    patchwork::plot_layout(heights = c(3, 1.4, 1.1, 1.8))

  outfile <- file.path(
    out_dir,
    sprintf(
      "%s_%s_%s_%s_%s_FP_RNA_corr_combined_R%.1f_summary_violin_plots%s.pdf",
      tf,
      ko_type,
      tf_cell,
      tf_line,
      mode,
      r_cut_main,
      filename_suffix
    )
  )

  if (verbose) {
    message("Saving summary violin + bar plot to: ", outfile)
  }

  ggplot2::ggsave(
    filename = outfile,
    plot     = combined,
    width    = 16,
    height   = 9,
    dpi      = 600
  )

  invisible(outfile)
}

## ==================================================================
## 6) KO-table helper
## ==================================================================

make_ko_tbl <- function(ko_raw, tf_label = "TF") {
  if (!is.data.frame(ko_raw)) {
    stop("KO table for ", tf_label, " is not a data.frame.")
  }

  gene_vec <- NULL

  if ("gene_symbol" %in% names(ko_raw) && "target_hgnc" %in% names(ko_raw)) {
    gene_vec <- dplyr::coalesce(ko_raw$gene_symbol, ko_raw$target_hgnc)
  } else {
    gene_cols <- c("symbol", "gene_symbol", "gene", "Gene", "HGNC",
                   "hgnc_symbol", "target_hgnc")
    gcol <- gene_cols[gene_cols %in% names(ko_raw)][1]

    if (is.na(gcol) || !length(gcol)) {
      stop("Could not find a gene column for ", tf_label,
           " (tried: ", paste(gene_cols, collapse = ", "), ").")
    }

    gene_vec <- ko_raw[[gcol]]
  }

  lfc_cols <- c("log2fc", "log2FoldChange", "logFC", "LFC")
  lcol <- lfc_cols[lfc_cols %in% names(ko_raw)][1]

  if (is.na(lcol) || !length(lcol)) {
    stop("Could not find a log2FC column for ", tf_label,
         " (tried: ", paste(lfc_cols, collapse = ", "), ").")
  }

  tibble::tibble(
    gene   = gene_vec,
    log2FC = ko_raw[[lcol]]
  ) %>%
    dplyr::filter(!is.na(gene), !is.na(log2FC))
}

## ==================================================================
## 7) Batch runner for KO violin plots (DB-based KO truth)
## ==================================================================

run_all_ko_violin <- function(tfs,
                              ko_truth_list,
                              tf_tfbs_list,
                              tf_gene_links_gh,
                              strict_rna,
                              TFLink,
                              ko_type = "KO",
                              out_dir = ".",
                              r_cut_main = 0.3,
                              p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                              cor_method = "pearson",
                              verbose = TRUE,
                              filename_suffix = "",
                              fp_p_col = c("p_adj_fp", "p_fp")) {

  fp_p_col <- match.arg(fp_p_col)

  # keep only TFs that have both KO truth and TFBS
  tfs_use <- tfs[tfs %in% intersect(names(ko_truth_list), names(tf_tfbs_list))]
  if (!length(tfs_use)) {
    stop("run_all_ko_violin(): no TFs present in both ko_truth_list and tf_tfbs_list.")
  }

  if (verbose) {
    message("Running KO violin summaries for TFs: ",
            paste(tfs_use, collapse = ", "))
  }

  res <- lapply(tfs_use, function(tf_symbol) {
    if (verbose) {
      message("  -> TF = ", tf_symbol)
    }

    ko_raw  <- ko_truth_list[[tf_symbol]]
    tf_tfbs <- tf_tfbs_list[[tf_symbol]]

    # make compact KO table (gene, log2FC)
    ko_tbl <- make_ko_tbl(ko_raw, tf_label = tf_symbol)

    # 1) All GeneHancer links
    out_all <- plot_ko_summary_violin_for_tf(
      tf              = tf_symbol,
      ko_tbl          = ko_tbl,
      tf_gene_links   = tf_gene_links_gh,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs,
      TFLink          = TFLink,
      mode            = "all",
      r_cut_main      = r_cut_main,
      p_cuts          = p_cuts,
      cor_method      = cor_method,
      ko_type         = ko_type,
      out_dir         = out_dir,
      verbose         = verbose,
      filename_suffix = filename_suffix,
      fp_p_col        = fp_p_col
    )

    # 2) Canonical (TF-centric) subset
    out_canonical <- plot_ko_summary_violin_for_tf(
      tf              = tf_symbol,
      ko_tbl          = ko_tbl,
      tf_gene_links   = tf_gene_links_gh,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs,
      TFLink          = TFLink,
      mode            = "canonical",
      r_cut_main      = r_cut_main,
      p_cuts          = p_cuts,
      cor_method      = cor_method,
      ko_type         = ko_type,
      out_dir         = out_dir,
      verbose         = verbose,
      filename_suffix = filename_suffix,
      fp_p_col        = fp_p_col
    )

    list(all = out_all, canonical = out_canonical)
  })

  names(res) <- tfs_use
  invisible(res)
}



ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

## ==================================================================
## 9) TF perturbation DB helpers and TFBS loading
## ==================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Parse TF from dataset_id like "..._<TF>_(KD|KO|OE|CRISPRa|CRISPRi)"
.extract_tf_from_dataset_id <- function(dataset_id) {
  if (length(dataset_id) == 0L) return(character(0))

  out <- vapply(dataset_id, function(id) {
    toks <- strsplit(as.character(id), "_", fixed = TRUE)[[1]]
    if (length(toks) < 3) return(NA_character_)

    last <- toupper(toks[length(toks)])
    if (!grepl("^(KO|KD|OE|CRISPR[AI]?)$", last)) return(NA_character_)

    toupper(toks[length(toks) - 1])
  }, character(1))

  out
}

# given TF name, return tibble of DE rows for that TF
get_tf_perturbation_tbl <- function(db_path, tf) {
  tf <- toupper(as.character(tf))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # read tables
  md <- DBI::dbReadTable(con, "metadata")
  de <- DBI::dbReadTable(con, "de_result")

  md <- tibble::as_tibble(md)
  de <- tibble::as_tibble(de)

  # decide how to get TF per dataset
  if ("target" %in% names(md)) {
    tf_col <- toupper(as.character(md$target))
  } else if ("TF" %in% names(md)) {
    tf_col <- toupper(as.character(md$TF))
  } else {
    # fall back to parsing from dataset_id
    if (!"dataset_id" %in% names(md)) {
      stop("metadata table has no 'dataset_id' column; can't infer TF.")
    }
    tf_col <- .extract_tf_from_dataset_id(md$dataset_id)
  }

  md$tf_parsed <- tf_col

  md_tf <- md[!is.na(md$tf_parsed) & md$tf_parsed == tf, , drop = FALSE]

  if (!nrow(md_tf)) {
    message("No metadata rows found for TF = ", tf,
            " (check TF name or dataset_id naming convention).")
    return(tibble::tibble())
  }

  # keep only DE rows for those dataset_ids
  if (!"dataset_id" %in% names(de)) {
    stop("de_result table has no 'dataset_id' column; cannot link to metadata.")
  }

  de_tf <- de[de$dataset_id %in% md_tf$dataset_id, , drop = FALSE]

  if (!nrow(de_tf)) {
    message("No DE rows found in de_result for TF = ", tf, ".")
    return(tibble::tibble())
  }

  # join metadata onto DE
  res <- dplyr::left_join(
    de_tf,
    md_tf,
    by = "dataset_id"
  )

  # Move a few useful columns to the front if they exist
  front_cols <- c(
    "dataset_id", "tf_parsed",
    "hgnc", "ensembl_id", "log2fc", "p_value", "p_adj", "base_mean",
    "sub_type", "accession", "perturbation_type", "comparison"
  )
  front_cols <- intersect(front_cols, names(res))

  res <- dplyr::relocate(res, dplyr::all_of(front_cols), .before = 1)

  res
}

TFLink <- readr::read_tsv(
  "/data/homes/cy232/genome/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv"
)

tf_perturb_db <- "/data/homes/yl814/episcope/tf_perturb.db"

get_tf_tfbs <- function(tf,
                        tfbs_dir      = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs",
                        tfbs_r_cut    = 0.3,
                        tfbs_padj_cut = 0.05) {

  overview_file  <- file.path(tfbs_dir, sprintf("%s_overview.txt", tf))
  predicted_tfbs <- readr::read_tsv(overview_file, show_col_types = FALSE)

  predicted_tfbs_sig <- predicted_tfbs[
    !is.na(predicted_tfbs$corr_fp_tf_r) &
      !is.na(predicted_tfbs$corr_fp_tf_p_adj) &
      predicted_tfbs$corr_fp_tf_r     > tfbs_r_cut &
      predicted_tfbs$corr_fp_tf_p_adj < tfbs_padj_cut,
    ,
    drop = FALSE
  ]

  unique(
    paste0(
      predicted_tfbs_sig$TFBS_chr, ":",
      predicted_tfbs_sig$TFBS_start, "-",
      predicted_tfbs_sig$TFBS_end
    )
  )
}

## ==================================================================
## 10) Build KO truth tables from DB
## ==================================================================

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

## TFBS per TF ---------------------------------------------------------

tf_tfbs_HNF1A  <- get_tf_tfbs("HNF1A",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_HNF4A  <- get_tf_tfbs("HNF4A",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")

tf_tfbs_IRF1   <- get_tf_tfbs("IRF1",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_RARG   <- get_tf_tfbs("RARG",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_SOX9   <- get_tf_tfbs("SOX9",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_KLF5   <- get_tf_tfbs("KLF5",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_FOXA2  <- get_tf_tfbs("FOXA2",
                              tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")

## KO truth tables -----------------------------------------------------

ko_truth_HNF1A <- make_ko_truth_from_db("HNF1A", tf_perturb_db)

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

## Lists used by run_all_ko_violin ------------------------------------

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

## TFs to plot ---------------------------------------------------------

tfs_all <- c("IRF1", "RARG", "SOX9", "KLF5", "FOXA2", "HNF1A", "HNF4A")

## ==================================================================
## 11) Main call: run KO violins for all TFs
## ==================================================================

# run_all_ko_violin(
#   tfs             = tfs_all,
#   ko_truth_list   = ko_truth_list,   # named list: IRF1, RARG, SOX9, ...
#   tf_tfbs_list    = tf_tfbs_list,    # named list: IRF1, RARG, SOX9, ...
#   tf_gene_links_gh = tf_gene_links_gh,
#   strict_rna      = strict_rna,
#   TFLink          = TFLink,
#   ko_type         = "KO",
#   out_dir         = ko_dir,
#   r_cut_main      = 0.3,
#   p_cuts          = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
#   cor_method      = "pearson",
#   verbose         = TRUE
# )
## Cell-line–specific GRN sets
run_all_ko_violin(
  tfs              = tfs_all,
  ko_truth_list    = ko_truth_list,
  tf_tfbs_list     = tf_tfbs_list,
  tf_gene_links_gh = tf_gene_links_gh_Panc1,
  strict_rna       = strict_rna,
  TFLink           = TFLink,
  ko_type          = "KO",
  out_dir          = ko_dir,
  r_cut_main       = 0.3,
  p_cuts           = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
  cor_method       = "pearson",
  verbose          = TRUE,
  filename_suffix  = "_grn_set_Panc1",
  fp_p_col = "p_fp"
)

run_all_ko_violin(
  tfs              = tfs_all,
  ko_truth_list    = ko_truth_list,
  tf_tfbs_list     = tf_tfbs_list,
  tf_gene_links_gh = tf_gene_links_gh_AsPC1,
  strict_rna       = strict_rna,
  TFLink           = TFLink,
  ko_type          = "KO",
  out_dir          = ko_dir,
  r_cut_main       = 0.3,
  p_cuts           = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
  cor_method       = "pearson",
  verbose          = TRUE,
  filename_suffix  = "_grn_set_AsPC1",
  fp_p_col = "p_fp"
)

run_all_ko_violin(
  tfs              = tfs_all,
  ko_truth_list    = ko_truth_list,
  tf_tfbs_list     = tf_tfbs_list,
  tf_gene_links_gh = tf_gene_links_gh_HPAFII,
  strict_rna       = strict_rna,
  TFLink           = TFLink,
  ko_type          = "KO",
  out_dir          = ko_dir,
  r_cut_main       = 0.3,
  p_cuts           = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
  cor_method       = "pearson",
  verbose          = TRUE,
  filename_suffix  = "_grn_set_HPAFII",
  fp_p_col = "p_fp"
)

run_all_ko_violin(
  tfs              = tfs_all,
  ko_truth_list    = ko_truth_list,
  tf_tfbs_list     = tf_tfbs_list,
  tf_gene_links_gh = tf_gene_links_gh,
  strict_rna       = strict_rna,
  TFLink           = TFLink,
  ko_type          = "KO",
  out_dir          = ko_dir,
  r_cut_main       = 0.3,
  p_cuts           = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
  cor_method       = "pearson",
  verbose          = TRUE,
  filename_suffix  = "_grn_set",
  fp_p_col = "p_fp"
)

HNF1A_KO <- get_tf_perturbation_tbl(tf_perturb_db, "HNF1A")
HNF4A_KO <- readr::read_csv(
  "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/Mayo 5289 siHNF4A RNA-seq.csv"
)

## ==================================================================
## 8) Legacy: per-TF manual plotting (requires *_KO in workspace)
##      – kept for backward compatibility
## ==================================================================

## ------------------------------------------------------------------
## Legacy GRN sets for per-TF manual plotting
## ------------------------------------------------------------------

legacy_grn_sets <- setNames(
  list(
    tf_gene_links_gh,        # full 23-sample GRN (original)
    tf_gene_links_gh_Panc1,  # Panc1-specific GRN
    tf_gene_links_gh_AsPC1,  # AsPC1-specific GRN
    tf_gene_links_gh_HPAFII  # HPAFII-specific GRN
  ),
  c(
    "_grn_set",                      # suffix for full-cohort GRN
    "_grn_set_Panc1",
    "_grn_set_AsPC1",
    "_grn_set_HPAFII"
  )
)

## HNF1A ----------------------------------------------------------------
if (exists("HNF1A_KO")) {
  ko_tbl_HNF1A <- make_ko_tbl(HNF1A_KO, "HNF1A")

  for (suffix in names(legacy_grn_sets)) {
    tf_links_curr <- legacy_grn_sets[[suffix]]

    # skip invalid / empty GRN sets
    if (!is.data.frame(tf_links_curr) || nrow(tf_links_curr) == 0L) {
      message("Skipping GRN set '", suffix, "' for HNF1A: not a data.frame or has zero rows.")
      next
    }

    # mode = "all"
    plot_ko_summary_violin_for_tf(
      tf              = "HNF1A",
      ko_tbl          = ko_tbl_HNF1A,
      tf_gene_links   = tf_links_curr,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs_HNF1A,
      TFLink          = TFLink,
      mode            = "all",
      r_cut_main      = 0.3,
      p_cuts          = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
      cor_method      = "pearson",
      ko_type         = "KO",
      out_dir         = ko_dir,
      filename_suffix = suffix,
      fp_p_col = "p_fp"
    )

    # mode = "canonical"
    plot_ko_summary_violin_for_tf(
      tf              = "HNF1A",
      ko_tbl          = ko_tbl_HNF1A,
      tf_gene_links   = tf_links_curr,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs_HNF1A,
      TFLink          = TFLink,
      mode            = "canonical",
      r_cut_main      = 0.3,
      p_cuts          = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
      cor_method      = "pearson",
      ko_type         = "KO",
      out_dir         = ko_dir,
      filename_suffix = suffix,
      fp_p_col = "p_fp"
    )


  }
}

## HNF4A ----------------------------------------------------------------
if (exists("HNF4A_KO")) {
  ko_tbl_HNF4A <- make_ko_tbl(HNF4A_KO, "HNF4A")

  for (suffix in names(legacy_grn_sets)) {
    tf_links_curr <- legacy_grn_sets[[suffix]]

    # skip invalid / empty GRN sets
    if (!is.data.frame(tf_links_curr) || nrow(tf_links_curr) == 0L) {
      message("Skipping GRN set '", suffix,
              "' for HNF4A: not a data.frame or has zero rows.")
      next
    }

    # mode = "all"
    plot_ko_summary_violin_for_tf(
      tf              = "HNF4A",
      ko_tbl          = ko_tbl_HNF4A,
      tf_gene_links   = tf_links_curr,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs_HNF4A,
      TFLink          = TFLink,
      mode            = "all",
      r_cut_main      = 0.3,
      p_cuts          = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
      cor_method      = "pearson",
      ko_type         = "KO",
      out_dir         = ko_dir,
      filename_suffix = suffix,
      fp_p_col = "p_fp"
    )

    # mode = "canonical"
    plot_ko_summary_violin_for_tf(
      tf              = "HNF4A",
      ko_tbl          = ko_tbl_HNF4A,
      tf_gene_links   = tf_links_curr,
      strict_rna      = strict_rna,
      tf_tfbs         = tf_tfbs_HNF4A,
      TFLink          = TFLink,
      mode            = "canonical",
      r_cut_main      = 0.3,
      p_cuts          = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
      cor_method      = "pearson",
      ko_type         = "KO",
      out_dir         = ko_dir,
      filename_suffix = suffix,
      fp_p_col = "p_fp"
    )
  }
}



## ==================================================================
## Minimal HNF1A / SOX9 plots for FP + TFLink, p = 5e-5, 5e-4
## ==================================================================

## 1) Build KO tables for HNF1A and SOX9 -----------------------------

# compact KO tables (gene, log2FC)
ko_tbl_HNF1A <- make_ko_tbl(HNF1A_KO, tf_label = "HNF1A")
ko_tbl_SOX9  <- make_ko_tbl(ko_truth_SOX9, tf_label = "SOX9")

## 2) Build summary tables restricted to FP + TFLink & p-cuts -------

p_use <- c(5e-5, 5e-4)

df_HNF1A <- build_ko_summary_table_for_tf(
  tf            = "HNF1A",
  ko_tbl        = ko_tbl_HNF1A,
  tf_gene_links = tf_gene_links_gh,
  strict_rna    = strict_rna,
  tf_tfbs       = tf_tfbs_HNF1A,
  TFLink        = TFLink,
  mode          = "canonical",
  r_cut_main    = 0.3,
  p_cuts        = p_use,
  cor_method    = "pearson",
  fp_p_col      = "p_fp",
  verbose       = FALSE
)

df_SOX9 <- build_ko_summary_table_for_tf(
  tf            = "SOX9",
  ko_tbl        = ko_tbl_SOX9,
  tf_gene_links = tf_gene_links_gh,
  strict_rna    = strict_rna,
  tf_tfbs       = tf_tfbs_SOX9,
  TFLink        = TFLink,
  mode          = "canonical",
  r_cut_main    = 0.3,
  p_cuts        = p_use,
  cor_method    = "pearson",
  fp_p_col      = "p_fp",
  verbose       = FALSE
)

df_fp_tflink <- dplyr::bind_rows(df_HNF1A, df_SOX9) %>%
  dplyr::filter(
    method_id == "FP_TFLink",
    p_cut %in% p_use
  )

## 3) Simple violins: Predicted vs Non-predicted ---------------------

fill_vals <- c(
  "Predicted"       = "#1b9e77",
  "Non-predicted"   = "#d95f02"
)

p_violin_simple <- ggplot(
  df_fp_tflink,
  aes(x = group, y = log2FC, fill = group)
) +
  geom_violin(trim = FALSE, alpha = 0.4, colour = "black") +
  stat_summary(
    fun        = median,
    geom       = "point",
    shape      = 95,      # horizontal line
    size       = 8,       # thickness/length of the line
    colour     = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = fill_vals, name = "Group") +
  facet_grid(rows = vars(tf), cols = vars(p_lab)) +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(
    name = "log2FC (KO vs Ctrl)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = NA),
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    axis.text.y      = element_text(face = "bold")
  )

## 4) Row-4 style stacked bar (log2FC bins, ≥0 in grey) --------------

df_bins_small <- df_fp_tflink %>%
  dplyr::mutate(
    log2FC_bin = dplyr::case_when(
      is.na(log2FC)                ~ NA_character_,
      log2FC <= -1                 ~ "<= -1",
      log2FC > -1 & log2FC <= -0.5 ~ "(-1,-0.5]",
      log2FC > -0.5 & log2FC < 0   ~ "(-0.5,0)",
      log2FC >= 0                  ~ ">= 0"
    )
  ) %>%
  dplyr::filter(!is.na(log2FC_bin))

bin_df_small <- df_bins_small %>%
  dplyr::group_by(tf, method_lab, p_lab, group, log2FC_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::group_by(tf, method_lab, p_lab, group) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  ) %>%
  dplyr::ungroup()

log2fc_levels <- c( "<= -1", "(-1,-0.5]", "(-0.5,0)",">= 0" )
bin_df_small$log2FC_bin <- factor(bin_df_small$log2FC_bin, levels = log2fc_levels)

log2fc_fill_grey <- c(
  ">= 0"       = "grey70",   # was blue, now grey
  "(-0.5,0)"   = "#fdb462",
  "(-1,-0.5]"  = "#fb8072",
  "<= -1"      = "#e41a1c"
)

p_bar_pct_small <- ggplot(
  bin_df_small,
  aes(x = group, y = pct, fill = log2FC_bin)
) +
  geom_col(colour = "black", linewidth = 0.2) +
  facet_grid(rows = vars(tf), cols = vars(p_lab)) +
  scale_fill_manual(values = log2fc_fill_grey, name = "log2FC bin") +
  scale_x_discrete(name = "Group") +
  scale_y_continuous(
    name   = "Percent of genes within group",
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 100)
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "grey", colour = NA),
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold", angle = 45, hjust = 1),
    axis.text.y      = element_text(face = "bold", size = 7),
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(face = "bold")
  )

## 5) (Optional) save to PDF -----------------------------------------
ggplot2::ggsave(file.path(ko_dir, "HNF1A_SOX9_FP_TFLink_simple_violins.pdf"),
                p_violin_simple, width = 6, height = 4, dpi = 600)
ggplot2::ggsave(file.path(ko_dir, "HNF1A_SOX9_FP_TFLink_simple_row4.pdf"),
                p_bar_pct_small, width = 3, height = 4, dpi = 600)
