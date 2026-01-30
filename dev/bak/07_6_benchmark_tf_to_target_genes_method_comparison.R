suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

## TF / KO cell line / lineage mapping -------------------------------
tf_cell_lineage <- tibble::tribble(
  ~tf,     ~cell,    ~lineage, ~type,
  "HNF1A", "AsPC1",  "inter", "KD",
  "HNF4A", "Mayo5289","classical", "KD",
  "IRF1",  "CFPAC1", "classical", "KD",
  "SOX9",  "PANC1",  "basal", "KD",
  "KLF5",  "AsPC1",  "inter", "KD",
  "FOXA2", "PANC1",  "basal", "KO",
  "RARG",  "PANC1",  "basal", "KO"
)

## Nutrient-stress cell → lineage mapping (for fp/rna sig sets) --------
cell_lineage_tbl <- tibble::tribble(
  ~cell,   ~lineage,
  "HPAFII","classical",
  "AsPC1", "inter",
  "Panc1", "basal"
)

## Small helper: get meta for one cell (only samples used in GRN) ------
get_cell_meta <- function(cell_name, strict_metadata) {
  strict_metadata[strict_metadata$cell == cell_name &
                    strict_metadata$run_grn, , drop = FALSE]
}

## FP: significant peaks per cell, then per lineage ------------------
##    - log2FC(stress / ctrl) >  1 & bound(stress) == 1
##    - log2FC(stress / ctrl) < -1 & bound(ctrl)   == 1
##    - significant if this holds in ≥1 stress vs Ctrl comparison
get_sig_fp_peaks_for_cell <- function(cell_name,
                                      fp_score,
                                      fp_bound,
                                      strict_metadata,
                                      logfc_cut = 1,
                                      pseudo = 1e-3) {
  meta_cell <- get_cell_meta(cell_name, strict_metadata)

  ctrl_id <- meta_cell$id[meta_cell$stress_type == "Ctrl"]
  if (length(ctrl_id) != 1L) {
    stop("Expect exactly one Ctrl sample for cell = ", cell_name)
  }
  stress_ids <- meta_cell$id[meta_cell$stress_type != "Ctrl"]

  sig <- rep(FALSE, nrow(fp_score))

  for (sid in stress_ids) {
    s_val <- fp_score[[sid]]
    c_val <- fp_score[[ctrl_id]]

    logfc <- log2((s_val + pseudo) / (c_val + pseudo))

    up   <- logfc >  logfc_cut & fp_bound[[sid]]      == 1
    down <- logfc < -logfc_cut & fp_bound[[ctrl_id]] == 1

    sig <- sig | up | down
  }

  fp_score$peak_ID[sig]
}



## RNA: significant genes per cell, then per lineage -----------------
##    - abs(log2FC(stress / ctrl)) > 1 in ≥1 comparison
##    - no bound/unbound requirement
get_sig_rna_genes_for_cell <- function(cell_name,
                                       strict_rna,
                                       strict_metadata,
                                       logfc_cut = 1,
                                       pseudo = 1e-3) {
  meta_cell <- get_cell_meta(cell_name, strict_metadata)

  ctrl_id <- meta_cell$id[meta_cell$stress_type == "Ctrl"]
  if (length(ctrl_id) != 1L) {
    stop("Expect exactly one Ctrl RNA sample for cell = ", cell_name)
  }
  stress_ids <- meta_cell$id[meta_cell$stress_type != "Ctrl"]

  expr <- strict_rna
  sig  <- rep(FALSE, nrow(expr))

  for (sid in stress_ids) {
    s_val <- expr[[sid]]
    c_val <- expr[[ctrl_id]]

    logfc <- log2((s_val + pseudo) / (c_val + pseudo))
    sig   <- sig | (abs(logfc) > logfc_cut)
  }

  gene_id <- dplyr::coalesce(expr$HGNC, expr$ensembl_gene_id)
  gene_id[sig]
}




# get list of fp_peak changed, i.e., abs(log2FC) > 1 in at least one condition vs ctrl for that cell line. per cell line list
# per-cell significant FP peak_IDs
fp_sig_AsPC1  <- get_sig_fp_peaks_for_cell("AsPC1",  fp_score, fp_bound, strict_metadata)
fp_sig_HPAFII <- get_sig_fp_peaks_for_cell("HPAFII", fp_score, fp_bound, strict_metadata)
fp_sig_PANC1  <- get_sig_fp_peaks_for_cell("Panc1",  fp_score, fp_bound, strict_metadata)

# collapse to per-lineage vectors
fp_sig_classical <- unique(fp_sig_HPAFII)
fp_sig_inter     <- unique(fp_sig_AsPC1)
fp_sig_basal     <- unique(fp_sig_PANC1)

# get list of gene change in RNA dataset for that cell line i.e., abs(log2FC) > 1
# per-cell significant genes
rna_sig_AsPC1  <- get_sig_rna_genes_for_cell("AsPC1",  strict_rna, strict_metadata)
rna_sig_HPAFII <- get_sig_rna_genes_for_cell("HPAFII", strict_rna, strict_metadata)
rna_sig_PANC1  <- get_sig_rna_genes_for_cell("Panc1",  strict_rna, strict_metadata)

# per-lineage gene vectors
rna_sig_classical <- unique(rna_sig_HPAFII)
rna_sig_inter     <- unique(rna_sig_AsPC1)
rna_sig_basal     <- unique(rna_sig_PANC1)


## ================================================================
## Shared helpers for KO split-violin grids
##   - use these to replace your existing geom_split_violin +
##     plot_fp_ko_grid / plot_atac_ko_grid / plot_rna_ko_grid
##   - plus a new FP+RNA combined grid: plot_fp_rna_ko_grid()
## ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(plyr)
  library(patchwork)
})

# Helper: format p-values as "1e-5" (no "p=" and no extra zero)
format_pval_short <- function(p) {
  s <- sprintf("%.0e", p)          # "1e-05"
  s <- sub("e-0", "e-", s)         # "1e-5"
  s <- sub("e\\+0", "e+", s)
  s
}

# Helper: add log2FC bins for KO tables
add_logfc_bin <- function(df) {
  df$logfc_bin <- cut(
    df$log2FC,
    breaks = c(-Inf, -1, -0.5, 0, Inf),
    labels = c("<= -1", "(-1,-0.5]", "(-0.5,0)", ">= 0"),
    right = TRUE,
    include.lowest = TRUE
  )
  df
}

# -------------------------------------------------------------------
# geom_split_violin helper (same behaviour as before)
# -------------------------------------------------------------------
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
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
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

# -------------------------------------------------------------------
# Helper: compact p-value labels, no "p = ", 1-digit exponent
#   e.g. 1e-05 -> "1e-5"; 3.2e-4 -> "1e-4"; 0.012 -> "0.012"
# -------------------------------------------------------------------
format_p_short <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
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

# -------------------------------------------------------------------
# Core plotting function used by FP / ATAC / RNA / FP+RNA
#   plot_df must have columns:
#     cut_r   (numeric), cut_p (numeric),
#     r_lab   (factor),  p_lab (factor),
#     group   (factor with 3 levels),
#     gene, log2FC
# -------------------------------------------------------------------
plot_ko_grid_core <- function(tf,
                              mode,
                              method_label,
                              plot_df,
                              p_xlab,
                              file_stub,
                              out_dir = ".",
                              test = c("wilcox", "t"),
                              alternative = c("two.sided", "greater", "less"),
                              verbose = TRUE) {

  test        <- match.arg(test)
  alternative <- match.arg(alternative)

  # Make sure factor levels are consistent
  plot_df$group <- factor(
    plot_df$group,
    levels = c("Predicted_regulated",
               "Random_background",
               "Predicted_nonregulated")
  )

  # Universe size = unique genes across groups (all in KO table)
  n_univ <- length(unique(plot_df$gene))

  # Small helper: run chosen test on one facet
  run_test <- function(df_sub) {
    df_sub <- df_sub[!is.na(df_sub$log2FC), , drop = FALSE]
    if (nrow(df_sub) == 0L) return(NA_real_)

    gtab <- table(droplevels(df_sub$group))
    if (length(gtab) != 2L || any(gtab == 0L)) return(NA_real_)

    if (test == "wilcox") {
      stats::wilcox.test(
        log2FC ~ group,
        data        = df_sub,
        alternative = alternative
      )$p.value
    } else {
      stats::t.test(
        log2FC ~ group,
        data        = df_sub,
        alternative = alternative
      )$p.value
    }
  }

  # -------------------------------------------------------------------
  # p-values (row 1 & 2)
  # -------------------------------------------------------------------
  global_range <- range(plot_df$log2FC, na.rm = TRUE)
  offset       <- 0.05 * diff(global_range)

  y_pos <- plot_df %>%
    dplyr::group_by(r_lab, p_lab) %>%
    dplyr::summarise(
      y = max(log2FC, na.rm = TRUE),
      .groups = "drop"
    )

  # Row1: predicted vs random
  df_rand_all <- plot_df[plot_df$group %in% c("Predicted_regulated",
                                              "Random_background"),
                         , drop = FALSE]

  pval_rand <- df_rand_all %>%
    dplyr::group_by(r_lab, p_lab) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all()),
      .groups = "drop"
    )

  annot_rand <- dplyr::left_join(pval_rand, y_pos,
                                 by = c("r_lab", "p_lab")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(p_val)  ~ "N.S.",
        p_val >= 0.05 ~ "N.S.",
        TRUE          ~ format_p_short(p_val)
      ),
      y = y + offset
    )

  # Row2: predicted vs non-regulated
  df_non_all <- plot_df[plot_df$group %in% c("Predicted_regulated",
                                             "Predicted_nonregulated"),
                        , drop = FALSE]

  pval_non <- df_non_all %>%
    dplyr::group_by(r_lab, p_lab) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all()),
      .groups = "drop"
    )

  annot_non <- dplyr::left_join(pval_non, y_pos,
                                by = c("r_lab", "p_lab")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(p_val)  ~ "N.S.",
        p_val >= 0.05 ~ "N.S.",
        TRUE          ~ format_p_short(p_val)
      ),
      y = y + offset
    )

  # -------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------
  if (verbose) {
    message("=== Diagnostics for TF = ", tf,
            " (mode = ", mode, ", method = ", method_label, ") ===")
    message("Universe size (n_univ): ", n_univ)
    message("Test: ", test, ", alternative = ", alternative)

    message("Random vs predicted: group counts per facet (head):")
    counts_rand <- df_rand_all %>%
      dplyr::group_by(r_lab, p_lab, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_rand, 10))

    message("Predicted vs non-regulated: group counts per facet (head):")
    counts_non <- df_non_all %>%
      dplyr::group_by(r_lab, p_lab, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_non, 10))

    if (all(is.na(pval_rand$p_val))) {
      message("Random vs predicted: all p-values are NA.")
    } else {
      message("Random vs predicted: ",
              "n facets with data = ", sum(!is.na(pval_rand$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_rand$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_rand$p_val, na.rm = TRUE), 3))
    }

    if (all(is.na(pval_non$p_val))) {
      message("Predicted vs non-regulated: all p-values are NA.")
    } else {
      message("Predicted vs non-regulated: ",
              "n facets with data = ", sum(!is.na(pval_non$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_non$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_non$p_val, na.rm = TRUE), 3))
    }
  }

  # -------------------------------------------------------------------
  # Shared palette
  # -------------------------------------------------------------------
  fill_vals <- c(
    Predicted_regulated    = "#1b9e77",
    Random_background      = "#7570b3",
    Predicted_nonregulated = "#d95f02"
  )

  # -------------------------------------------------------------------
  # Row 1: predicted vs random (split violins)
  # -------------------------------------------------------------------
  df_rand <- df_rand_all

  p_rand <- ggplot(
    df_rand,
    aes(x = p_lab, y = log2FC, fill = group)
  ) +
    geom_split_violin(alpha = 0.4, trim = FALSE, colour = "black") +
    geom_boxplot(
      width        = 0.2,
      alpha        = 0.6,
      show.legend  = FALSE,
      outlier.size = 0.1,
      outlier.stroke = 0,
      outlier.alpha  = 0.3
    ) +
    geom_text(
      data = annot_rand,
      aes(x = p_lab, y = y, label = label),
      inherit.aes = FALSE,
      size        = 2.8,
      vjust       = 0,
      fontface    = "bold"
    ) +
    facet_wrap(~ r_lab, nrow = 1) +
    scale_fill_manual(
      values = fill_vals[c("Predicted_regulated", "Random_background")],
      name   = "Group",
      labels = c("Predicted regulated", "Random background")
    ) +
    scale_x_discrete(name = p_xlab) +
    scale_y_continuous(name = sprintf("log2FC (%s KO vs Ctrl)", tf)) +
    ggtitle(
      sprintf(
        "%s knockout: %s predicted regulated vs random background genes\n(mode = %s; n = %d genes in universe)",
        tf,
        method_label,
        mode,
        n_univ
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      legend.text      = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_text(face = "bold"),
      axis.text.y      = element_text(face = "bold")
    )

  # -------------------------------------------------------------------
  # Row 2: predicted vs non-regulated (split violins)
  # -------------------------------------------------------------------
  df_non <- df_non_all

  p_non <- ggplot(
    df_non,
    aes(x = p_lab, y = log2FC, fill = group)
  ) +
    geom_split_violin(alpha = 0.4, trim = FALSE, colour = "black") +
    geom_boxplot(
      width        = 0.2,
      alpha        = 0.6,
      show.legend  = FALSE,
      outlier.size = 0.1,
      outlier.stroke = 0,
      outlier.alpha  = 0.3
    ) +
    geom_text(
      data = annot_non,
      aes(x = p_lab, y = y, label = label),
      inherit.aes = FALSE,
      size        = 2.8,
      vjust       = 0,
      fontface    = "bold"
    ) +
    facet_wrap(~ r_lab, nrow = 1) +
    scale_fill_manual(
      values = fill_vals[c("Predicted_regulated", "Predicted_nonregulated")],
      name   = "Group",
      labels = c("Predicted regulated", "Predicted non-regulated")
    ) +
    scale_x_discrete(name = p_xlab) +
    scale_y_continuous(name = sprintf("log2FC (%s KO vs Ctrl)", tf)) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      legend.text      = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_text(face = "bold"),
      axis.text.y      = element_text(face = "bold")
    )

  # -------------------------------------------------------------------
  # Row 3: counts (Predicted_regulated ONLY)
  # -------------------------------------------------------------------
  counts_df <- df_non_all %>%
    dplyr::filter(group == "Predicted_regulated") %>%
    dplyr::group_by(r_lab, p_lab) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    )

  p_counts <- ggplot(
    counts_df,
    aes(x = p_lab, y = n_genes)
  ) +
    geom_col(
      width     = 0.6,
      colour    = "black",
      linewidth = 0.2
    ) +
    geom_text(
      aes(label = n_genes),
      vjust = -0.2,
      size  = 2.5
    ) +
    facet_wrap(~ r_lab, nrow = 1) +
    scale_x_discrete(name = p_xlab) +
    scale_y_continuous(
      name   = "Number of unique predicted regulated genes",
      expand = expansion(mult = c(0, 0.15))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "none",
      strip.background = element_rect(fill = "grey90", colour = NA),
      strip.text       = element_text(face = "bold"),
      axis.title.x     = element_text(face = "bold"),
      axis.title.y     = element_text(face = "bold"),
      axis.text.x      = element_text(face = "bold"),
      axis.text.y      = element_text(face = "bold")
    )

  # -------------------------------------------------------------------
  # Row 4: grouped stacked percent bars of log2FC bins
  #   Predicted_regulated vs Predicted_nonregulated
  #   - SAME format for FP / ATAC / RNA / FP+RNA
  #   - groups are now grouped on x instead of facet rows
  # -------------------------------------------------------------------
  df_bins <- df_non_all %>%
    dplyr::filter(group %in% c("Predicted_regulated",
                               "Predicted_nonregulated")) %>%
    add_logfc_bin() %>%
    dplyr::filter(!is.na(logfc_bin))

  if (nrow(df_bins) > 0L) {
    counts_pct <- df_bins %>%
      dplyr::group_by(r_lab, p_lab, group, logfc_bin) %>%
      dplyr::summarise(
        n_genes = dplyr::n_distinct(gene),
        .groups = "drop"
      ) %>%
      dplyr::group_by(r_lab, p_lab, group) %>%
      dplyr::mutate(
        pct = if (sum(n_genes) > 0) 100 * n_genes / sum(n_genes) else 0
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        logfc_bin = factor(
          logfc_bin,
          levels = c("<= -1", "(-1,-0.5]", "(-0.5,0)", ">= 0")
        ),
        # grouped x-axis: each p_cut × group as a bar
        x_gp = interaction(p_lab, group, sep = "\n")
      )

    logfc_cols <- c(
      "<= -1"      = "#b2182b",
      "(-1,-0.5]"  = "#ef8a62",
      "(-0.5,0)"   = "#fddbc7",
      ">= 0"       = "#d1e5f0"
    )

    p_stack <- ggplot(
      counts_pct,
      aes(x = x_gp, y = pct, fill = logfc_bin)
    ) +
      geom_col(
        width     = 0.6,
        colour    = "black",
        linewidth = 0.2
      ) +
      facet_wrap(~ r_lab, nrow = 1) +
      scale_fill_manual(
        values = logfc_cols,
        name   = "log2FC bin"
      ) +
      scale_x_discrete(name = paste0(p_xlab, " × group")) +
      scale_y_continuous(
        name   = "Percent of genes",
        limits = c(0, 100),
        expand = expansion(mult = c(0, 0.02))
      ) +
      theme_minimal(base_size = 10) +
      theme(
        legend.position  = "right",
        strip.background = element_rect(fill = "grey90", colour = NA),
        strip.text       = element_text(face = "bold"),
        axis.title.x     = element_text(face = "bold"),
        axis.title.y     = element_text(face = "bold"),
        axis.text.x      = element_text(face = "bold", size = 7, angle = 45, hjust = 1),
        axis.text.y      = element_text(face = "bold", size = 7)
      )
  } else {
    p_stack <- ggplot() + theme_void()
  }

  # -------------------------------------------------------------------
  # Combine rows and save
  # -------------------------------------------------------------------
  combined <- p_rand / p_non / p_counts / p_stack +
    patchwork::plot_layout(heights = c(3, 3, 1.5, 1.5))

  outfile <- file.path(out_dir, paste0(file_stub, ".pdf"))

  ggplot2::ggsave(
    filename = outfile,
    plot     = combined,
    width    = 18,
    height   = 12,
    dpi      = 600
  )

  invisible(outfile)
}



# ===================================================================
# 1) FP-correlation KO grid (slightly trimmed version of your old one)
#     - filter_atac retained, but filenames no longer contain
#       the literal string "_atac_filtered"
# ===================================================================
plot_fp_ko_grid <- function(tf,
                            ko_tbl,
                            tf_gene_links,
                            atac_gene_links,
                            tf_tfbs,
                            mode = c("canonical", "all"),
                            filter_atac = FALSE,
                            r_fp_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                            p_fp_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                            threshold_r_atac = 0.3,
                            threshold_p_adj_atac = 0.05,
                            seed = 1L,
                            out_dir = ".",
                            test = c("wilcox", "t"),
                            alternative = c("two.sided", "greater", "less"),
                            r_sign = c("abs", "pos", "neg"),
                            verbose = TRUE) {

  mode        <- match.arg(mode)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)
  r_sign      <- match.arg(r_sign)

  # ---------------- TF links + optional ATAC filtering ----------------
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs.")
  }

  if (filter_atac) {
    atac_keep <- atac_gene_links[
      !is.na(atac_gene_links$r_atac) &
        !is.na(atac_gene_links$p_adj_atac) &
        abs(atac_gene_links$r_atac) >= threshold_r_atac &
        atac_gene_links$p_adj_atac <= threshold_p_adj_atac,
      c("atac_peak", "gene_key"),
      drop = FALSE
    ]

    if (nrow(atac_keep) == 0L) {
      stop("ATAC filter removed all rows; relax r_atac / p_adj_atac cutoffs.")
    }

    tf_links_use <- dplyr::inner_join(
      tf_links0,
      atac_keep,
      by = c("atac_peak", "gene_key")
    )
  } else {
    tf_links_use <- tf_links0
  }

  if (nrow(tf_links_use) == 0L) {
    stop("No TF links left after ATAC filtering.")
  }

  # Universe of genes for this filter setting
  genes_univ <- intersect(unique(tf_links_use$gene_key), ko_tbl$gene)
  ko_univ <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  if (!length(genes_univ)) {
    stop("No overlapping genes between KO table and tf_gene_links.")
  }

  # Background gene pool for random (always unfiltered, excluding TFBS peaks)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) set.seed(seed)

  # Helper: for one (r_fp, p_fp) pair, get passing / non-passing / random genes
  get_sets_one_cut <- function(r_cut, p_cut) {
    tf_links_curr <- tf_links_use

    if (mode == "canonical") {
      motifs_vec <- toupper(as.character(tf_links_curr$motifs))
      tf_upper   <- toupper(tf)
      keep_motif <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)
      tf_links_curr <- tf_links_curr[keep_motif, , drop = FALSE]
    }

    if (nrow(tf_links_curr) == 0L) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    pass_rows <- !is.na(tf_links_curr$r_fp) &
      !is.na(tf_links_curr$p_adj_fp) &
      tf_links_curr$p_adj_fp < p_cut

    if (r_sign == "abs") {
      pass_rows <- pass_rows & abs(tf_links_curr$r_fp) > r_cut
    } else if (r_sign == "pos") {
      pass_rows <- pass_rows & tf_links_curr$r_fp > r_cut
    } else {
      pass_rows <- pass_rows & tf_links_curr$r_fp < -r_cut
    }

    genes_all  <- sort(unique(tf_links_curr$gene_key))
    genes_pass <- sort(unique(tf_links_curr$gene_key[pass_rows]))
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else character(0)

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  # ---------------- Build grid over (r_fp, p_fp) ----------------
  grid <- expand.grid(
    r_fp_cut = r_fp_cuts,
    p_fp_cut = p_fp_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_fp_cut[i]
    p_cut <- grid$p_fp_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Predicted_regulated",
        gene     = ko_pass$gene,
        log2FC   = ko_pass$log2FC
      ),
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Predicted_nonregulated",
        gene     = ko_non$gene,
        log2FC   = ko_non$log2FC
      ),
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Random_background",
        gene     = ko_rand$gene,
        log2FC   = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  # Factor labels for facets and x-axis
  r_fp_levels <- r_fp_cuts
  p_fp_levels <- sort(unique(p_fp_cuts))

  # Factor labels for facets and x-axis
  r_fp_levels <- r_fp_cuts
  p_fp_levels <- sort(unique(p_fp_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_fp_cut,
        levels = r_fp_levels,
        labels = vapply(r_fp_levels, function(rc) {
          if (r_sign == "abs") {
            paste0("|r_fp| > ", rc)
          } else if (r_sign == "pos") {
            paste0("r_fp > ", rc)
          } else {
            paste0("r_fp < -", rc)
          }
        }, character(1))
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
      )
    )


  # Filenames: label sign mode
  atac_lab  <- if (filter_atac) "filtered" else "unfiltered"
  sign_tag  <- switch(r_sign,
                      abs = "",
                      pos = "_rposOnly",
                      neg = "_rnegOnly")
  file_stub <- sprintf(
    "%s_KO_%s_%s_fpCorr%s_grid_split_violin_vsRand_vsNon",
    tf,
    atac_lab,
    mode,
    sign_tag
  )

  sign_label <- switch(
    r_sign,
    abs = "",
    pos = " (r_fp>0 only)",
    neg = " (r_fp<0 only)"
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = paste0(if (filter_atac) "FP (ATAC-filtered)" else "FP",
                          sign_label),
    plot_df      = plot_df,
    p_xlab       = "p_adj_fp cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}


# ===================================================================
# 2) ATAC-correlation KO grid
# ===================================================================
plot_atac_ko_grid <- function(tf,
                              ko_tbl,
                              tf_gene_links,
                              atac_gene_links,
                              tf_tfbs,
                              mode = c("canonical", "all"),
                              r_atac_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                              p_atac_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                              seed = 1L,
                              out_dir = ".",
                              test = c("wilcox", "t"),
                              alternative = c("two.sided", "greater", "less"),
                              r_sign = c("abs", "pos", "neg"),
                              verbose = TRUE) {

  mode        <- match.arg(mode)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)
  r_sign      <- match.arg(r_sign)

  # TF links for this TF + its TFBS
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (ATAC grid).")
  }

  # Universe of genes (as in FP grid): TF-linked genes × KO table
  genes_univ <- intersect(unique(tf_links0$gene_key), ko_tbl$gene)
  ko_univ <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  if (!length(genes_univ)) {
    stop("No overlapping genes between KO table and tf_gene_links (ATAC grid).")
  }

  # Background gene pool for random (TF-unbound peaks)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) set.seed(seed)

  # Left join ATAC correlations
  tf_links_atac <- dplyr::left_join(
    tf_links0,
    atac_gene_links,
    by = c("atac_peak", "gene_key")
  )

  # Helper: for one (r_atac, p_atac) pair
  get_sets_one_cut <- function(r_cut, p_cut) {
    tf_links_curr <- tf_links_atac

    if (mode == "canonical") {
      motifs_vec <- toupper(as.character(tf_links_curr$motifs))
      tf_upper   <- toupper(tf)
      keep_motif <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)
      tf_links_curr <- tf_links_curr[keep_motif, , drop = FALSE]
    }

    if (nrow(tf_links_curr) == 0L) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    pass_rows <- !is.na(tf_links_curr$r_atac) &
      !is.na(tf_links_curr$p_adj_atac) &
      tf_links_curr$p_adj_atac < p_cut

    if (r_sign == "abs") {
      pass_rows <- pass_rows & abs(tf_links_curr$r_atac) > r_cut
    } else if (r_sign == "pos") {
      pass_rows <- pass_rows & tf_links_curr$r_atac > r_cut
    } else {
      pass_rows <- pass_rows & tf_links_curr$r_atac < -r_cut
    }

    genes_all  <- sort(unique(tf_links_curr$gene_key))
    genes_pass <- sort(unique(tf_links_curr$gene_key[pass_rows]))
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else character(0)

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  grid <- expand.grid(
    r_atac_cut = r_atac_cuts,
    p_atac_cut = p_atac_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_atac_cut[i]
    p_cut <- grid$p_atac_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_atac_cut = r_cut,
        p_atac_cut = p_cut,
        group      = "Predicted_regulated",
        gene       = ko_pass$gene,
        log2FC     = ko_pass$log2FC
      ),
      tibble::tibble(
        r_atac_cut = r_cut,
        p_atac_cut = p_cut,
        group      = "Predicted_nonregulated",
        gene       = ko_non$gene,
        log2FC     = ko_non$log2FC
      ),
      tibble::tibble(
        r_atac_cut = r_cut,
        p_atac_cut = p_cut,
        group      = "Random_background",
        gene       = ko_rand$gene,
        log2FC     = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  r_atac_levels <- r_atac_cuts
  p_atac_levels <- sort(unique(p_atac_cuts))

  r_atac_levels <- r_atac_cuts
  p_atac_levels <- sort(unique(p_atac_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_atac_cut,
        levels = r_atac_levels,
        labels = vapply(r_atac_levels, function(rc) {
          if (r_sign == "abs") {
            paste0("|r_atac| > ", rc)
          } else if (r_sign == "pos") {
            paste0("r_atac > ", rc)
          } else {
            paste0("r_atac < -", rc)
          }
        }, character(1))
      ),
      p_lab = factor(
        p_atac_cut,
        levels = p_atac_levels,
        labels = scales::label_scientific(digits = 1)(p_atac_levels)
      )
    )


  sign_tag <- switch(r_sign,
                     abs = "",
                     pos = "_rposOnly",
                     neg = "_rnegOnly")

  file_stub <- sprintf(
    "%s_KO_ATACcorr_tfFiltered_%s%s_grid_split_violin_vsRand_vsNon",
    tf,
    mode,
    sign_tag
  )

  sign_label <- switch(
    r_sign,
    abs = "",
    pos = " (r_atac>0 only)",
    neg = " (r_atac<0 only)"
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = paste0("ATAC-correlated", sign_label),
    plot_df      = plot_df,
    p_xlab       = "p_adj_atac cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}

# ===================================================================
# 3) RNA-correlation KO grid
# ===================================================================
plot_rna_ko_grid <- function(tf,
                             ko_tbl,
                             tf_gene_links,
                             strict_rna,
                             tf_tfbs,
                             mode = c("canonical", "all"),
                             r_rna_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                             p_rna_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                             cor_method = c("pearson", "spearman"),
                             seed = 1L,
                             out_dir = ".",
                             test = c("wilcox", "t"),
                             alternative = c("two.sided", "greater", "less"),
                             r_sign = c("abs", "pos", "neg"),
                             verbose = TRUE) {

  mode        <- match.arg(mode)
  cor_method  <- match.arg(cor_method)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)
  r_sign      <- match.arg(r_sign)

  # Map gene keys to RNA rows
  map_gene_keys_to_rna <- function(keys, strict_rna) {
    idx_hgnc <- match(keys, strict_rna$HGNC)
    idx_ens  <- match(keys, strict_rna$ensembl_gene_id)
    row_idx  <- ifelse(!is.na(idx_hgnc), idx_hgnc, idx_ens)
    keep     <- !is.na(row_idx)
    tibble::tibble(
      gene_key = keys[keep],
      row_idx  = row_idx[keep]
    )
  }

  # Extract TF-linked genes & background
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (RNA grid).")
  }

  genes_tf    <- sort(unique(tf_links0$gene_key))
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes    <- sort(unique(tf_links_bg$gene_key))
  bg_genes    <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  genes_univ0 <- intersect(genes_tf, ko_tbl$gene)
  if (!length(genes_univ0)) {
    stop("No overlapping genes between KO table and tf_gene_links (RNA grid).")
  }

  if (!is.null(seed)) set.seed(seed)

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

  map_univ <- map_gene_keys_to_rna(genes_univ0, strict_rna)
  if (!nrow(map_univ)) {
    stop("No universe genes found in strict_rna (RNA grid).")
  }

  genes_univ <- map_univ$gene_key
  ko_univ    <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]

  # Compute TF–gene RNA correlations
  corr_list <- lapply(seq_len(nrow(map_univ)), function(i) {
    g      <- map_univ$gene_key[i]
    row_i  <- map_univ$row_idx[i]
    g_expr <- as.numeric(rna_expr_mat[row_i, ])

    ok <- is.finite(tf_expr) & is.finite(g_expr)
    n  <- sum(ok)
    if (n < 3L) {
      return(list(gene_key = g,
                  n_rna    = n,
                  r_rna    = NA_real_,
                  p_rna    = NA_real_))
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

  rna_corr_full <- dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )

  genes_univ <- intersect(genes_univ, rna_corr_full$gene_key)
  ko_univ    <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                                 , drop = FALSE]
  if (!length(genes_univ)) {
    stop("After RNA correlation filtering, no genes remain in universe (RNA grid).")
  }

  # Precompute mode-specific TF-linked genes
  if (mode == "canonical") {
    motifs_vec   <- toupper(as.character(tf_links0$motifs))
    tf_upper     <- toupper(tf)
    keep_motif   <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)
    tf_links_mode <- tf_links0[keep_motif, , drop = FALSE]
  } else {
    tf_links_mode <- tf_links0
  }

  genes_mode_all <- intersect(
    sort(unique(tf_links_mode$gene_key)),
    genes_univ
  )

  get_sets_one_cut <- function(r_cut, p_cut) {
    if (!length(genes_mode_all)) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    g_sub <- rna_corr_full[rna_corr_full$gene_key %in% genes_mode_all,
                           , drop = FALSE]

    pass_rows <- !is.na(g_sub$r_rna) &
      !is.na(g_sub$p_adj_rna) &
      g_sub$p_adj_rna < p_cut

    if (r_sign == "abs") {
      pass_rows <- pass_rows & abs(g_sub$r_rna) > r_cut
    } else if (r_sign == "pos") {
      pass_rows <- pass_rows & g_sub$r_rna > r_cut
    } else {
      pass_rows <- pass_rows & g_sub$r_rna < -r_cut
    }

    genes_pass <- sort(g_sub$gene_key[pass_rows])
    genes_all  <- genes_mode_all
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else character(0)

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  grid <- expand.grid(
    r_rna_cut = r_rna_cuts,
    p_rna_cut = p_rna_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_rna_cut[i]
    p_cut <- grid$p_rna_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Predicted_regulated",
        gene      = ko_pass$gene,
        log2FC    = ko_pass$log2FC
      ),
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Predicted_nonregulated",
        gene      = ko_non$gene,
        log2FC    = ko_non$log2FC
      ),
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Random_background",
        gene      = ko_rand$gene,
        log2FC    = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  r_rna_levels <- r_rna_cuts
  p_rna_levels <- sort(unique(p_rna_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = vapply(r_rna_levels, function(rc) {
          if (r_sign == "abs") {
            paste0("|r_rna| > ", rc)
          } else if (r_sign == "pos") {
            paste0("r_rna > ", rc)
          } else {
            paste0("r_rna < -", rc)
          }
        }, character(1))
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
      )
    )


  sign_tag <- switch(r_sign,
                     abs = "",
                     pos = "_rposOnly",
                     neg = "_rnegOnly")

  file_stub <- sprintf(
    "%s_KO_rna_corr_tfFiltered_%s%s_grid_split_violin_vsRand_vsNon",
    tf,
    mode,
    sign_tag
  )

  sign_label <- switch(
    r_sign,
    abs = "",
    pos = " (r_rna>0 only)",
    neg = " (r_rna<0 only)"
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = paste0("RNA TF–gene correlated", sign_label),
    plot_df      = plot_df,
    p_xlab       = "p_adj_rna cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}


# ===================================================================
# 4) NEW: FP + RNA combined KO grid
#     - predicted regulated genes must pass BOTH:
#         FP correlation (r_fp, p_adj_fp)
#         RNA TF–gene correlation (r_rna, p_adj_rna)
# ===================================================================
# ===================================================================
# 4) NEW: FP + RNA combined KO grid
#     - predicted regulated genes must pass BOTH:
#         FP correlation (r_fp, p_adj_fp)
#         RNA TF–gene correlation (r_rna, p_adj_rna)
#     - plotting is delegated to plot_ko_grid_core() so that
#       FP / ATAC / RNA / FP+RNA share the same layout
# ===================================================================
plot_fp_rna_ko_grid <- function(tf,
                                ko_tbl,
                                tf_gene_links,
                                strict_rna,
                                tf_tfbs,
                                mode = c("canonical", "all"),
                                r_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                                p_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                                cor_method = c("pearson", "spearman"),
                                seed = 1L,
                                out_dir = ".",
                                test = c("wilcox", "t"),
                                alternative = c("two.sided", "greater", "less"),
                                r_sign = c("abs", "pos", "neg"),
                                verbose = TRUE) {

  mode        <- match.arg(mode)
  cor_method  <- match.arg(cor_method)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)
  r_sign      <- match.arg(r_sign)

  # ---------------- Helper: map gene_keys -> strict_rna row indices ---------
  map_gene_keys_to_rna <- function(keys, strict_rna) {
    idx_hgnc <- match(keys, strict_rna$HGNC)
    idx_ens  <- match(keys, strict_rna$ensembl_gene_id)
    row_idx  <- ifelse(!is.na(idx_hgnc), idx_hgnc, idx_ens)
    keep     <- !is.na(row_idx)
    tibble::tibble(
      gene_key = keys[keep],
      row_idx  = row_idx[keep]
    )
  }

  # ---------------- Extract TF-linked genes & background --------------------
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (FP+RNA grid).")
  }

  genes_tf <- sort(unique(tf_links0$gene_key))

  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  genes_univ0 <- intersect(genes_tf, ko_tbl$gene)
  if (!length(genes_univ0)) {
    stop("No overlapping genes between KO table and tf_gene_links (FP+RNA grid).")
  }

  if (!is.null(seed)) set.seed(seed)

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

  # ---------------- Map universe genes to RNA rows -------------------------
  map_univ <- map_gene_keys_to_rna(genes_univ0, strict_rna)
  if (!nrow(map_univ)) {
    stop("No universe genes found in strict_rna (FP+RNA grid).")
  }

  genes_univ <- map_univ$gene_key
  ko_univ    <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]

  # ---------------- Compute TF–gene RNA correlations -----------------------
  corr_list <- lapply(seq_len(nrow(map_univ)), function(i) {
    g      <- map_univ$gene_key[i]
    row_i  <- map_univ$row_idx[i]
    g_expr <- as.numeric(rna_expr_mat[row_i, ])

    ok <- is.finite(tf_expr) & is.finite(g_expr)
    n  <- sum(ok)
    if (n < 3L) {
      return(list(gene_key = g,
                  n_rna    = n,
                  r_rna    = NA_real_,
                  p_rna    = NA_real_))
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

  rna_corr_full <- dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )

  genes_univ <- intersect(genes_univ, rna_corr_full$gene_key)
  ko_univ    <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                                 , drop = FALSE]

  if (!length(genes_univ)) {
    stop("After RNA correlation filtering, no genes remain in universe (FP+RNA grid).")
  }

  # ---------------- Precompute mode-specific FP universe -------------------
  if (mode == "canonical") {
    motifs_vec   <- toupper(as.character(tf_links0$motifs))
    tf_upper     <- toupper(tf)
    keep_motif   <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)
    tf_links_mode <- tf_links0[keep_motif, , drop = FALSE]
  } else {
    tf_links_mode <- tf_links0
  }

  tf_links_mode <- tf_links_mode[tf_links_mode$gene_key %in% genes_univ,
                                 , drop = FALSE]

  genes_mode_all <- intersect(
    sort(unique(tf_links_mode$gene_key)),
    genes_univ
  )

  # Helper: for one (r, p) pair, get FP+RNA passing / non-passing / random genes
  get_sets_one_cut <- function(r_cut, p_cut) {
    if (!length(genes_mode_all)) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    # RNA filter
    g_sub <- rna_corr_full[rna_corr_full$gene_key %in% genes_mode_all,
                           , drop = FALSE]
    pass_rna <- !is.na(g_sub$r_rna) &
      !is.na(g_sub$p_adj_rna) &
      g_sub$p_adj_rna < p_cut

    # FP filter (on the same mode-specific links)
    if (nrow(tf_links_mode) == 0L) {
      genes_pass_fp <- character(0)
    } else {
      pass_fp <- !is.na(tf_links_mode$r_fp) &
        !is.na(tf_links_mode$p_adj_fp) &
        tf_links_mode$p_adj_fp < p_cut
    }

    if (r_sign == "abs") {
      pass_rna <- pass_rna & abs(g_sub$r_rna) > r_cut
      if (nrow(tf_links_mode) > 0L) {
        pass_fp <- pass_fp & abs(tf_links_mode$r_fp) > r_cut
      }
    } else if (r_sign == "pos") {
      pass_rna <- pass_rna & g_sub$r_rna > r_cut
      if (nrow(tf_links_mode) > 0L) {
        pass_fp <- pass_fp & tf_links_mode$r_fp > r_cut
      }
    } else { # "neg"
      pass_rna <- pass_rna & g_sub$r_rna < -r_cut
      if (nrow(tf_links_mode) > 0L) {
        pass_fp <- pass_fp & tf_links_mode$r_fp < -r_cut
      }
    }

    genes_pass_rna <- sort(g_sub$gene_key[pass_rna])
    genes_pass_fp  <- if (nrow(tf_links_mode) > 0L) {
      sort(unique(tf_links_mode$gene_key[pass_fp]))
    } else {
      character(0)
    }

    genes_pass <- intersect(genes_pass_rna, genes_pass_fp)
    genes_all  <- genes_mode_all
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else {
      character(0)
    }

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  # ---------------- Build grid over shared (r, p) --------------------------
  grid <- expand.grid(
    r_cut  = r_cuts,
    p_cut  = p_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_cut[i]
    p_cut <- grid$p_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Predicted_regulated",
        gene    = ko_pass$gene,
        log2FC  = ko_pass$log2FC
      ),
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Predicted_nonregulated",
        gene    = ko_non$gene,
        log2FC  = ko_non$log2FC
      ),
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Random_background",
        gene    = ko_rand$gene,
        log2FC  = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  # Factor labels for facets and x-axis (shared r / p)
  r_levels <- r_cuts
  p_levels <- sort(unique(p_cuts))


  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_cut,
        levels = r_levels,
        labels = vapply(r_levels, function(rc) {
          if (r_sign == "abs") {
            paste0("|r_fp| & |r_rna| > ", rc)
          } else if (r_sign == "pos") {
            paste0("r_fp > ", rc, " & r_rna > ", rc)
          } else {
            paste0("r_fp < -", rc, " & r_rna < -", rc)
          }
        }, character(1))
      ),
      p_lab = factor(
        p_cut,
        levels = p_levels,
        labels = scales::label_scientific(digits = 1)(p_levels)
      ),
      group = factor(
        group,
        levels = c(
          "Predicted_regulated",
          "Random_background",
          "Predicted_nonregulated"
        )
      )
    )


  sign_tag <- switch(r_sign,
                     abs = "",
                     pos = "_rposOnly",
                     neg = "_rnegOnly")

  file_stub <- sprintf(
    "%s_KO_FPplusRNA_%s%s_grid_split_violin_vsRand_vsNon",
    tf,
    mode,
    sign_tag
  )

  sign_label <- switch(
    r_sign,
    abs = "",
    pos = " (r_fp>0 & r_rna>0 only)",
    neg = " (r_fp<0 & r_rna<0 only)"
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = paste0("FP & RNA correlated", sign_label),
    plot_df      = plot_df,
    p_xlab       = "shared p_adj cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}




## ================================================================
## Example usage of the 4 KO grid functions
##   - FP only
##   - ATAC corr only
##   - RNA corr only
##   - FP + RNA (both filters must pass)
## ================================================================

## ------------------------------------------------
## 0) Quick reminder of objects you already have
## ------------------------------------------------
# HNF1A KO table
# ko_tbl <- HNF1A_KO %>%
#   dplyr::transmute(
#     gene   = dplyr::coalesce(gene_symbol, symbol),
#     log2FC = log2fc
#   ) %>%
#   dplyr::filter(!is.na(gene))

# HNF4A KO table
# ko_tbl_HNF4A <- HNF4A_KO %>%
#   dplyr::transmute(
#     gene   = symbol,
#     log2FC = log2FoldChange
#   ) %>%
#   dplyr::filter(!is.na(gene))

# tf_gene_links      <- readr::read_csv("fp_gene_corr_full_jaspar2024.csv")
# atac_gene_links    <- readr::read_csv("atac_gene_corr_full_jaspar2024.csv")
# strict_rna         <- strict_rna   # already built in your script
# tf_tfbs_HNF1A      <- get_tf_tfbs("HNF1A", tfbs_dir = ...)
# tf_tfbs_HNF4A      <- get_tf_tfbs("HNF4A", tfbs_dir = ...)
# ko_dir             <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"


## ------------------------------------------------
## 1) FP-correlation alone
##    - unfiltered (all FP links)
##    - ATAC-filtered version (same FP cut, but restricted to ATAC-corr edges)
##    Filenames:
##      HNF1A_KO_unfiltered_all_fpCorr_grid_split_violin_vsRand_vsNon.pdf
##      HNF1A_KO_filtered_all_fpCorr_grid_split_violin_vsRand_vsNon.pdf
## ------------------------------------------------

sign_modes <- c("abs", "pos", "neg")

# HNF1A
for (sm in sign_modes) {
  plot_fp_ko_grid(
    tf              = "HNF1A",
    ko_tbl          = ko_tbl,
    tf_gene_links   = tf_gene_links,
    atac_gene_links = atac_gene_links,
    tf_tfbs         = tf_tfbs_HNF1A,
    mode            = "all",
    filter_atac     = FALSE,
    r_fp_cuts       = c(0, 0.1, 0.3, 0.5, 0.7),
    p_fp_cuts       = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    test            = "wilcox",
    alternative     = "two.sided",
    r_sign          = sm,
    verbose         = TRUE,
    out_dir         = ko_dir
  )

  plot_atac_ko_grid(
    tf              = "HNF1A",
    ko_tbl          = ko_tbl,
    tf_gene_links   = tf_gene_links,
    atac_gene_links = atac_gene_links,
    tf_tfbs         = tf_tfbs_HNF1A,
    mode            = "all",
    r_atac_cuts     = c(0, 0.1, 0.3, 0.5, 0.7),
    p_atac_cuts     = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    test            = "wilcox",
    alternative     = "two.sided",
    r_sign          = sm,
    verbose         = TRUE,
    out_dir         = ko_dir
  )

  plot_rna_ko_grid(
    tf            = "HNF1A",
    ko_tbl        = ko_tbl,
    tf_gene_links = tf_gene_links,
    strict_rna    = strict_rna,
    tf_tfbs       = tf_tfbs_HNF1A,
    mode          = "all",
    r_rna_cuts    = c(0, 0.1, 0.3, 0.5, 0.7),
    p_rna_cuts    = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    cor_method    = "pearson",
    test          = "wilcox",
    alternative   = "two.sided",
    r_sign        = sm,
    verbose       = TRUE,
    out_dir       = ko_dir
  )

  plot_fp_rna_ko_grid(
    tf            = "HNF1A",
    ko_tbl        = ko_tbl,
    tf_gene_links = tf_gene_links,
    strict_rna    = strict_rna,
    tf_tfbs       = tf_tfbs_HNF1A,
    mode          = "all",
    r_cuts        = c(0, 0.1, 0.3, 0.5, 0.7),
    p_cuts        = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
    cor_method    = "pearson",
    test          = "wilcox",
    alternative   = "two.sided",
    r_sign        = sm,
    verbose       = TRUE,
    out_dir       = ko_dir
  )
}



## ------------------------------------------------
## Extra TFs from ko_truth_list / tf_tfbs_list
##   IRF1, RARG, SOX9, KLF5, FOXA2
##   - assumes ko_truth_* have columns: gene, log2FC
##   - uses GeneHancer-based links: tf_gene_links_gh / atac_gene_links_gh
## ------------------------------------------------

run_ko_grids_for_tf <- function(tf_symbol) {
  message("=== Running KO grids for ", tf_symbol, " ===")

  sign_modes <- c("abs", "pos", "neg")

  # 1) Build KO table in the format expected by the grid functions
  ko_tbl_tf <- ko_truth_list[[tf_symbol]] %>%
    dplyr::transmute(
      gene   = gene,
      log2FC = log2fc
    ) %>%
    dplyr::filter(!is.na(gene), !is.na(log2FC))

  # 2) TFBS peaks for this TF
  tf_tfbs_tf <- tf_tfbs_list[[tf_symbol]]

  for (sm in sign_modes) {
    # FP-only grid (no ATAC filter)
    plot_fp_ko_grid(
      tf              = tf_symbol,
      ko_tbl          = ko_tbl_tf,
      tf_gene_links   = tf_gene_links_gh,
      atac_gene_links = atac_gene_links_gh,
      tf_tfbs         = tf_tfbs_tf,
      mode            = "all",
      filter_atac     = FALSE,
      r_fp_cuts       = c(0, 0.1, 0.3, 0.5, 0.7),
      p_fp_cuts       = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      test            = "wilcox",
      alternative     = "two.sided",
      r_sign          = sm,
      verbose         = TRUE,
      out_dir         = ko_dir
    )

    # ATAC-corr-only grid
    plot_atac_ko_grid(
      tf              = tf_symbol,
      ko_tbl          = ko_tbl_tf,
      tf_gene_links   = tf_gene_links_gh,
      atac_gene_links = atac_gene_links_gh,
      tf_tfbs         = tf_tfbs_tf,
      mode            = "all",
      r_atac_cuts     = c(0, 0.1, 0.3, 0.5, 0.7),
      p_atac_cuts     = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      test            = "wilcox",
      alternative     = "two.sided",
      r_sign          = sm,
      verbose         = TRUE,
      out_dir         = ko_dir
    )

    # RNA-corr-only grid
    plot_rna_ko_grid(
      tf            = tf_symbol,
      ko_tbl        = ko_tbl_tf,
      tf_gene_links = tf_gene_links_gh,
      strict_rna    = strict_rna,
      tf_tfbs       = tf_tfbs_tf,
      mode          = "all",
      r_rna_cuts    = c(0, 0.1, 0.3, 0.5, 0.7),
      p_rna_cuts    = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      cor_method    = "pearson",
      test          = "wilcox",
      alternative   = "two.sided",
      r_sign        = sm,
      verbose       = TRUE,
      out_dir       = ko_dir
    )

    # FP + RNA combined grid
    plot_fp_rna_ko_grid(
      tf            = tf_symbol,
      ko_tbl        = ko_tbl_tf,
      tf_gene_links = tf_gene_links_gh,
      strict_rna    = strict_rna,
      tf_tfbs       = tf_tfbs_tf,
      mode          = "all",
      r_cuts        = c(0, 0.1, 0.3, 0.5, 0.7),
      p_cuts        = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      cor_method    = "pearson",
      test          = "wilcox",
      alternative   = "two.sided",
      r_sign        = sm,
      verbose       = TRUE,
      out_dir       = ko_dir
    )
  }

  invisible(NULL)
}


## Run for all TFs defined:
invisible(lapply(tfs, run_ko_grids_for_tf))












base_dir_window <- "/data/homes/yl814/episcope_test/nutrient_stress"

## GeneHancer-based
tf_gene_links_gh   <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_jaspar2024.csv"))
atac_gene_links_gh <- readr::read_csv(file.path(base_dir_window, "atac_gene_corr_full_jaspar2024.csv"))





suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

## TF / KO cell line / lineage mapping -------------------------------
tf_cell_lineage <- tibble::tribble(
  ~tf,     ~cell,    ~lineage,
  "HNF1A", "AsPC1",   "inter",
  "HNF4A", "Mayo5289","classical",
  "IRF1",  "CFPAC1",  "classical",
  "SOX9",  "PANC1",   "basal",
  "KLF5",  "AsPC1",   "inter",
  "FOXA2", "PANC1",   "basal",
  "RARG",  "PANC1",   "basal"
)

## Nutrient-stress cell → lineage (if needed elsewhere) --------------
cell_lineage_tbl <- tibble::tribble(
  ~cell,   ~lineage,
  "HPAFII","classical",
  "AsPC1", "inter",
  "Panc1", "basal"
)

## Small helper: get meta for one cell (only samples used in GRN) ----
get_cell_meta <- function(cell_name, strict_metadata) {
  strict_metadata[strict_metadata$cell == cell_name &
                    strict_metadata$run_grn, , drop = FALSE]
}

## TF → (cell, lineage) ----------------------------------------------
get_tf_cell_lineage <- function(tf) {
  hit <- tf_cell_lineage[match(tf, tf_cell_lineage$tf), , drop = FALSE]
  if (nrow(hit) == 0L || all(is.na(hit$cell))) {
    warning("No tf_cell_lineage entry for TF = ", tf,
            "; using 'NAcell' / 'NAlineage' in filenames.")
    return(list(cell = "NAcell", lineage = "NAlineage"))
  }
  list(cell = hit$cell[1], lineage = hit$lineage[1])
}

## Lineage-specific FP / RNA sig sets --------------------------------
get_fp_sig_for_lineage <- function(lineage) {
  if (!exists("fp_sig_classical", inherits = TRUE) ||
      !exists("fp_sig_inter", inherits = TRUE) ||
      !exists("fp_sig_basal", inherits = TRUE)) {
    return(NULL)
  }
  switch(
    lineage,
    classical = fp_sig_classical,
    inter     = fp_sig_inter,
    basal     = fp_sig_basal,
    NULL
  )
}

get_rna_sig_for_lineage <- function(lineage) {
  if (!exists("rna_sig_classical", inherits = TRUE) ||
      !exists("rna_sig_inter", inherits = TRUE) ||
      !exists("rna_sig_basal", inherits = TRUE)) {
    return(NULL)
  }
  switch(
    lineage,
    classical = rna_sig_classical,
    inter     = rna_sig_inter,
    basal     = rna_sig_basal,
    NULL
  )
}

# Tags for filenames ------------------------------------------------
make_fp_filter_tag <- function(lineage, applied) {
  if (!applied) {
    "fp_unfiltered"
  } else {
    paste0("fp_sig_", lineage, "_filtered")
  }
}

make_rna_filter_tag <- function(lineage, applied) {
  if (!applied) {
    "rna_unfiltered"
  } else {
    paste0("rna_sig_", lineage, "_filtered")
  }
}

make_r_sign_tag <- function(r_sign) {
  # we only support "pos" for now
  if (!identical(r_sign, "pos")) {
    warning("Only r_sign = 'pos' is currently supported; overriding to 'pos'.")
  }
  "pos"
}

plot_fp_ko_grid <- function(tf,
                            ko_tbl,
                            tf_gene_links,
                            atac_gene_links,
                            tf_tfbs,
                            mode = c("canonical", "all"),
                            r_fp_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                            p_fp_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                            seed = 1L,
                            out_dir = ".",
                            test = c("wilcox", "t"),
                            alternative = c("two.sided", "greater", "less"),
                            r_sign = "pos",       # we force 'pos'
                            ko_type = "KO",
                            tf_links_scope = c("mode", "all"),
                            filter_by_fp_sig = FALSE,
                            filter_by_rna_sig = FALSE,  # now actually used
                            verbose = TRUE) {

  mode           <- match.arg(mode)
  test           <- match.arg(test)
  alternative    <- match.arg(alternative)
  tf_links_scope <- match.arg(tf_links_scope)
  r_sign_tag     <- make_r_sign_tag(r_sign)  # always "pos"
  r_sign         <- "pos"

  # -------------------------------------------------
  # TF-specific metadata (cell, lineage, filters)
  # -------------------------------------------------
  tf_meta  <- get_tf_cell_lineage(tf)
  tf_cell  <- tf_meta$cell
  tf_line  <- tf_meta$lineage

  fp_filter_tag  <- make_fp_filter_tag(tf_line, filter_by_fp_sig)
  rna_filter_tag <- make_rna_filter_tag(tf_line, filter_by_rna_sig)

  # -------------------------------------------------
  # TF links: canonical vs all
  #   - 'canonical' = tf_gene_links with tfs == tf
  #   - 'all'       = all rows
  # -------------------------------------------------
  tf_links_all <- tf_gene_links
  if (mode == "canonical") {
    if (!"tfs" %in% names(tf_gene_links)) {
      stop("mode = 'canonical' requires a 'tfs' column in tf_gene_links (FP grid).")
    }
    tf_links_all <- tf_links_all[tf_links_all$tfs == tf, , drop = FALSE]
  }
  if (tf_links_scope == "all") {
    tf_links_all <- tf_gene_links
  }

  # Restrict to TFBS peaks for this TF
  tf_links0 <- tf_links_all[tf_links_all$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (FP grid).")
  }

  # Universe of genes: TF-linked genes ∩ KO table
  genes_univ <- intersect(unique(tf_links0$gene_key), ko_tbl$gene)
  ko_univ    <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  if (!length(genes_univ)) {
    stop("No overlapping genes between KO table and tf_gene_links (FP grid).")
  }

  # FP sig filter (lineage-specific peaks) ------------------------------
  tf_links_use <- tf_links0
  if (filter_by_fp_sig) {
    fp_sig <- get_fp_sig_for_lineage(tf_line)
    if (is.null(fp_sig)) {
      if (verbose) {
        warning("filter_by_fp_sig = TRUE but no fp_sig_* object found for lineage = ", tf_line)
      }
    } else {
      n_before_fp  <- nrow(tf_links_use)
      tf_links_use <- tf_links_use[tf_links_use$fp_peak %in% fp_sig, , drop = FALSE]
      if (verbose) {
        message("FP grid (", tf, "): FP sig filter kept ",
                nrow(tf_links_use), " / ", n_before_fp, " TF links.")
      }
      if (nrow(tf_links_use) == 0L) {
        stop("FP sig filter removed all TF links for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
      genes_univ <- intersect(unique(tf_links_use$gene_key), ko_univ$gene)
      ko_univ    <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
    }
  }

  # --- NEW: optional RNA sig filter on target gene universe -----------
  genes_univ0_raw <- genes_univ
  if (filter_by_rna_sig) {
    rna_sig <- get_rna_sig_for_lineage(tf_line)
    if (is.null(rna_sig)) {
      if (verbose) {
        warning("filter_by_rna_sig = TRUE but no rna_sig_* object found for lineage = ", tf_line)
      }
    } else {
      genes_univ <- intersect(genes_univ0_raw, rna_sig)
      if (verbose) {
        message("FP grid (", tf, "): RNA sig filter kept ",
                length(genes_univ), " / ", length(genes_univ0_raw),
                " KO×TF-linked genes (lineage = ", tf_line, ").")
      }
      if (!length(genes_univ)) {
        stop("RNA sig filter removed all genes for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
      ko_univ <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
    }
  } else if (verbose) {
    message("FP grid (", tf, "): no RNA sig filter; universe size = ",
            length(genes_univ), ".")
  }

  # Background genes: TF-unbound peaks (still from tf_gene_links)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes    <- sort(unique(tf_links_bg$gene_key))
  bg_genes    <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) set.seed(seed)

  # Helper: for one (r_fp, p_fp) pair, get gene sets
  get_sets_one_cut <- function(r_cut, p_cut) {
    tf_links_curr <- tf_links_use
    if (nrow(tf_links_curr) == 0L || !length(genes_univ)) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    pass_rows <- !is.na(tf_links_curr$r_fp) &
      !is.na(tf_links_curr$p_adj_fp) &
      tf_links_curr$p_adj_fp < p_cut

    # r_sign == "pos"
    pass_rows <- pass_rows & tf_links_curr$r_fp > r_cut

    genes_all  <- sort(unique(tf_links_curr$gene_key))
    genes_pass <- sort(unique(tf_links_curr$gene_key[pass_rows]))
    genes_fail <- setdiff(genes_all, genes_pass)

    # restrict to RNA-filtered universe
    genes_pass <- intersect(genes_pass, genes_univ)
    genes_fail <- intersect(genes_fail, genes_univ)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else character(0)

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  # Grid over (r_fp, p_fp)
  grid <- expand.grid(
    r_fp_cut = r_fp_cuts,
    p_fp_cut = p_fp_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_fp_cut[i]
    p_cut <- grid$p_fp_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- sets_i$passing
    non_i  <- sets_i$non_passing
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Predicted_regulated",
        gene     = ko_pass$gene,
        log2FC   = ko_pass$log2FC
      ),
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Predicted_nonregulated",
        gene     = ko_non$gene,
        log2FC   = ko_non$log2FC
      ),
      tibble::tibble(
        r_fp_cut = r_cut,
        p_fp_cut = p_cut,
        group    = "Random_background",
        gene     = ko_rand$gene,
        log2FC   = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  r_fp_levels <- r_fp_cuts
  p_fp_levels <- sort(unique(p_fp_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_fp_cut,
        levels = r_fp_levels,
        labels = vapply(r_fp_levels, function(rc) {
          paste0("r_fp > ", rc)
        }, character(1))
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
      )
    )

  # File naming: AsPC1_inter_HNF1A_KO_canonical_FP_corr_pos_fp_unfiltered_rna_unfiltered_...
  grid_type_tag <- "FP_corr"

  file_stub <- sprintf(
    "%s_%s_%s_%s_%s_%s_%s_%s_%s_summary_violin_plots",
    tf_cell,        # e.g. AsPC1
    tf_line,        # e.g. inter
    tf,             # e.g. HNF1A
    ko_type,        # e.g. KO
    mode,           # e.g. canonical
    grid_type_tag,  # FP_corr
    r_sign_tag,     # pos
    fp_filter_tag,  # fp_unfiltered or fp_sig_lineage_filtered
    rna_filter_tag  # rna_unfiltered or rna_sig_lineage_filtered
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = "FP-correlated (r_fp>0)",
    plot_df      = plot_df,
    p_xlab       = "p_adj_fp cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}


plot_rna_ko_grid <- function(tf,
                             ko_tbl,
                             tf_gene_links,
                             strict_rna,
                             tf_tfbs,
                             mode = c("canonical", "all"),
                             r_rna_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                             p_rna_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                             cor_method = c("pearson", "spearman"),
                             seed = 1L,
                             out_dir = ".",
                             test = c("wilcox", "t"),
                             alternative = c("two.sided", "greater", "less"),
                             r_sign = "pos",      # we force 'pos'
                             ko_type = "KO",
                             tf_links_scope = c("mode", "all"),
                             filter_by_fp_sig = FALSE,
                             filter_by_rna_sig = FALSE,
                             verbose = TRUE) {

  mode           <- match.arg(mode)
  cor_method     <- match.arg(cor_method)
  test           <- match.arg(test)
  alternative    <- match.arg(alternative)
  tf_links_scope <- match.arg(tf_links_scope)
  r_sign_tag     <- make_r_sign_tag(r_sign)
  r_sign         <- "pos"

  tf_meta  <- get_tf_cell_lineage(tf)
  tf_cell  <- tf_meta$cell
  tf_line  <- tf_meta$lineage

  fp_filter_tag  <- make_fp_filter_tag(tf_line, filter_by_fp_sig)
  rna_filter_tag <- make_rna_filter_tag(tf_line, filter_by_rna_sig)

  # Map gene keys to RNA rows
  map_gene_keys_to_rna <- function(keys, strict_rna) {
    idx_hgnc <- match(keys, strict_rna$HGNC)
    idx_ens  <- match(keys, strict_rna$ensembl_gene_id)
    row_idx  <- ifelse(!is.na(idx_hgnc), idx_hgnc, idx_ens)
    keep     <- !is.na(row_idx)
    tibble::tibble(
      gene_key = keys[keep],
      row_idx  = row_idx[keep]
    )
  }

  # TF links: canonical vs all
  tf_links_all <- tf_gene_links
  if (mode == "canonical") {
    if (!"tfs" %in% names(tf_gene_links)) {
      stop("mode = 'canonical' requires a 'tfs' column in tf_gene_links (RNA grid).")
    }
    tf_links_all <- tf_links_all[tf_links_all$tfs == tf, , drop = FALSE]
  }
  if (tf_links_scope == "all") {
    tf_links_all <- tf_gene_links
  }

  # Restrict to TFBS peaks
  tf_links0 <- tf_links_all[tf_links_all$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (RNA grid).")
  }

  # Optional FP-significant peak filter (lineage-specific)
  if (filter_by_fp_sig) {
    fp_sig <- get_fp_sig_for_lineage(tf_line)
    if (is.null(fp_sig)) {
      warning("filter_by_fp_sig = TRUE but no fp_sig_* object found for lineage = ", tf_line)
    } else {
      n_before_fp <- nrow(tf_links0)
      tf_links0 <- tf_links0[tf_links0$fp_peak %in% fp_sig, , drop = FALSE]
      if (verbose) {
        message("RNA grid (", tf, "): FP sig filter kept ",
                nrow(tf_links0), " / ", n_before_fp, " TF links.")
      }
      if (nrow(tf_links0) == 0L) {
        stop("FP sig filter removed all TF links for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
    }
  }

  genes_tf <- sort(unique(tf_links0$gene_key))

  # Background gene pool (TF-unbound peaks)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes    <- sort(unique(tf_links_bg$gene_key))
  bg_genes    <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  # Universe before RNA sig filter
  genes_univ0_raw <- intersect(genes_tf, ko_tbl$gene)
  if (!length(genes_univ0_raw)) {
    stop("No overlapping genes between KO table and tf_gene_links (RNA grid).")
  }

  genes_univ0 <- genes_univ0_raw

  # RNA sig filter (lineage-specific)
  if (filter_by_rna_sig) {
    rna_sig <- get_rna_sig_for_lineage(tf_line)
    if (is.null(rna_sig)) {
      warning("filter_by_rna_sig = TRUE but no rna_sig_* object found for lineage = ", tf_line)
    } else {
      genes_univ0 <- intersect(genes_univ0, rna_sig)
      if (verbose) {
        message("RNA grid (", tf, "): RNA sig filter kept ",
                length(genes_univ0), " / ", length(genes_univ0_raw),
                " KO×TF-linked genes (lineage = ", tf_line, ").")
      }
      if (!length(genes_univ0)) {
        stop("RNA sig filter removed all genes for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
    }
  } else if (verbose) {
    message("RNA grid (", tf, "): no RNA sig filter; universe size = ",
            length(genes_univ0_raw), ".")
  }

  if (!is.null(seed)) set.seed(seed)

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

  map_univ <- map_gene_keys_to_rna(genes_univ0, strict_rna)
  if (!nrow(map_univ)) {
    stop("No universe genes found in strict_rna (RNA grid).")
  }

  genes_univ <- map_univ$gene_key
  ko_univ    <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]

  # Compute TF–gene RNA correlations
  corr_list <- lapply(seq_len(nrow(map_univ)), function(i) {
    g      <- map_univ$gene_key[i]
    row_i  <- map_univ$row_idx[i]
    g_expr <- as.numeric(rna_expr_mat[row_i, ])

    ok <- is.finite(tf_expr) & is.finite(g_expr)
    n  <- sum(ok)
    if (n < 3L) {
      return(list(gene_key = g,
                  n_rna    = n,
                  r_rna    = NA_real_,
                  p_rna    = NA_real_))
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

  rna_corr_full <- dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )

  genes_univ <- intersect(genes_univ, rna_corr_full$gene_key)
  ko_univ    <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                                 , drop = FALSE]

  if (!length(genes_univ)) {
    stop("After RNA correlation filtering, no genes remain in universe (RNA grid).")
  }

  # Helper: for one (r_rna, p_rna) pair
  get_sets_one_cut <- function(r_cut, p_cut) {
    if (!length(genes_univ)) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    g_sub <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                           , drop = FALSE]

    pass_rows <- !is.na(g_sub$r_rna) &
      !is.na(g_sub$p_adj_rna) &
      g_sub$p_adj_rna < p_cut

    # r_sign == "pos"
    pass_rows <- pass_rows & g_sub$r_rna > r_cut

    genes_pass <- sort(g_sub$gene_key[pass_rows])
    genes_all  <- genes_univ
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else character(0)

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  grid <- expand.grid(
    r_rna_cut = r_rna_cuts,
    p_rna_cut = p_rna_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_rna_cut[i]
    p_cut <- grid$p_rna_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Predicted_regulated",
        gene      = ko_pass$gene,
        log2FC    = ko_pass$log2FC
      ),
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Predicted_nonregulated",
        gene      = ko_non$gene,
        log2FC    = ko_non$log2FC
      ),
      tibble::tibble(
        r_rna_cut = r_cut,
        p_rna_cut = p_cut,
        group     = "Random_background",
        gene      = ko_rand$gene,
        log2FC    = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  r_rna_levels <- r_rna_cuts
  p_rna_levels <- sort(unique(p_rna_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = vapply(r_rna_levels, function(rc) {
          paste0("r_rna > ", rc)
        }, character(1))
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
      )
    )

  grid_type_tag <- "RNA_corr"

  file_stub <- sprintf(
    "%s_%s_%s_%s_%s_%s_%s_%s_%s_summary_violin_plots",
    tf_cell,
    tf_line,
    tf,
    ko_type,
    mode,
    grid_type_tag,
    r_sign_tag,
    fp_filter_tag,
    rna_filter_tag
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = "RNA TF–gene correlated (r_rna>0)",
    plot_df      = plot_df,
    p_xlab       = "p_adj_rna cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}


plot_fp_rna_ko_grid <- function(tf,
                                ko_tbl,
                                tf_gene_links,
                                strict_rna,
                                tf_tfbs,
                                mode = c("canonical", "all"),
                                r_cuts = c(0, 0.1, 0.3, 0.5, 0.7),
                                p_cuts = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
                                cor_method = c("pearson", "spearman"),
                                seed = 1L,
                                out_dir = ".",
                                test = c("wilcox", "t"),
                                alternative = c("two.sided", "greater", "less"),
                                r_sign = "pos",      # we force 'pos'
                                ko_type = "KO",
                                tf_links_scope = c("mode", "all"),
                                filter_by_fp_sig = FALSE,
                                filter_by_rna_sig = FALSE,
                                verbose = TRUE) {

  mode           <- match.arg(mode)
  cor_method     <- match.arg(cor_method)
  test           <- match.arg(test)
  alternative    <- match.arg(alternative)
  tf_links_scope <- match.arg(tf_links_scope)
  r_sign_tag     <- make_r_sign_tag(r_sign)
  r_sign         <- "pos"

  tf_meta  <- get_tf_cell_lineage(tf)
  tf_cell  <- tf_meta$cell
  tf_line  <- tf_meta$lineage

  fp_filter_tag  <- make_fp_filter_tag(tf_line, filter_by_fp_sig)
  rna_filter_tag <- make_rna_filter_tag(tf_line, filter_by_rna_sig)

  # Map gene_keys -> strict_rna row indices
  map_gene_keys_to_rna <- function(keys, strict_rna) {
    idx_hgnc <- match(keys, strict_rna$HGNC)
    idx_ens  <- match(keys, strict_rna$ensembl_gene_id)
    row_idx  <- ifelse(!is.na(idx_hgnc), idx_hgnc, idx_ens)
    keep     <- !is.na(row_idx)
    tibble::tibble(
      gene_key = keys[keep],
      row_idx  = row_idx[keep]
    )
  }

  # TF links: canonical vs all
  tf_links_all <- tf_gene_links
  if (mode == "canonical") {
    if (!"tfs" %in% names(tf_gene_links)) {
      stop("mode = 'canonical' requires a 'tfs' column in tf_gene_links (FP+RNA grid).")
    }
    tf_links_all <- tf_links_all[tf_links_all$tfs == tf, , drop = FALSE]
  }
  if (tf_links_scope == "all") {
    tf_links_all <- tf_gene_links
  }

  tf_links0 <- tf_links_all[tf_links_all$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (FP+RNA grid).")
  }

  # Optional FP sig filter early
  if (filter_by_fp_sig) {
    fp_sig <- get_fp_sig_for_lineage(tf_line)
    if (is.null(fp_sig)) {
      warning("filter_by_fp_sig = TRUE but no fp_sig_* object found for lineage = ", tf_line)
    } else {
      n_before_fp <- nrow(tf_links0)
      tf_links0 <- tf_links0[tf_links0$fp_peak %in% fp_sig, , drop = FALSE]
      if (verbose) {
        message("FP+RNA grid (", tf, "): FP sig filter kept ",
                nrow(tf_links0), " / ", n_before_fp, " TF links.")
      }
      if (nrow(tf_links0) == 0L) {
        stop("FP sig filter removed all TF links for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
    }
  }

  genes_tf <- sort(unique(tf_links0$gene_key))

  # Background genes (unbound TFBS)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes    <- sort(unique(tf_links_bg$gene_key))
  bg_genes    <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  # Universe before RNA sig filter
  genes_univ0_raw <- intersect(genes_tf, ko_tbl$gene)
  if (!length(genes_univ0_raw)) {
    stop("No overlapping genes between KO table and tf_gene_links (FP+RNA grid).")
  }

  genes_univ0 <- genes_univ0_raw

  # RNA sig filter
  if (filter_by_rna_sig) {
    rna_sig <- get_rna_sig_for_lineage(tf_line)
    if (is.null(rna_sig)) {
      warning("filter_by_rna_sig = TRUE but no rna_sig_* object found for lineage = ", tf_line)
    } else {
      genes_univ0 <- intersect(genes_univ0, rna_sig)
      if (verbose) {
        message("FP+RNA grid (", tf, "): RNA sig filter kept ",
                length(genes_univ0), " / ", length(genes_univ0_raw),
                " KO×TF-linked genes (lineage = ", tf_line, ").")
      }
      if (!length(genes_univ0)) {
        stop("RNA sig filter removed all genes for TF = ", tf,
             " (lineage = ", tf_line, ").")
      }
    }
  } else if (verbose) {
    message("FP+RNA grid (", tf, "): no RNA sig filter; universe size = ",
            length(genes_univ0_raw), ".")
  }

  if (!is.null(seed)) set.seed(seed)

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

  # Map universe genes to RNA rows
  map_univ <- map_gene_keys_to_rna(genes_univ0, strict_rna)
  if (!nrow(map_univ)) {
    stop("No universe genes found in strict_rna (FP+RNA grid).")
  }

  genes_univ <- map_univ$gene_key
  ko_univ    <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]

  # Compute TF–gene RNA correlations
  corr_list <- lapply(seq_len(nrow(map_univ)), function(i) {
    g      <- map_univ$gene_key[i]
    row_i  <- map_univ$row_idx[i]
    g_expr <- as.numeric(rna_expr_mat[row_i, ])

    ok <- is.finite(tf_expr) & is.finite(g_expr)
    n  <- sum(ok)
    if (n < 3L) {
      return(list(gene_key = g,
                  n_rna    = n,
                  r_rna    = NA_real_,
                  p_rna    = NA_real_))
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

  rna_corr_full <- dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )

  genes_univ <- intersect(genes_univ, rna_corr_full$gene_key)
  ko_univ    <- ko_univ[ko_univ$gene %in% genes_univ, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                                 , drop = FALSE]

  if (!length(genes_univ)) {
    stop("After RNA correlation filtering, no genes remain in universe (FP+RNA grid).")
  }

  # FP mode-specific links (after FP sig filter; only genes in genes_univ)
  tf_links_mode <- tf_links0[tf_links0$gene_key %in% genes_univ, , drop = FALSE]

  genes_mode_all <- intersect(
    sort(unique(tf_links_mode$gene_key)),
    genes_univ
  )

  # Helper: for one (r, p) pair, get FP+RNA passing / non / random
  get_sets_one_cut <- function(r_cut, p_cut) {
    if (!length(genes_mode_all)) {
      return(list(
        passing     = character(0),
        non_passing = character(0),
        random      = character(0)
      ))
    }

    # RNA filter within genes_mode_all
    g_sub <- rna_corr_full[rna_corr_full$gene_key %in% genes_mode_all,
                           , drop = FALSE]
    pass_rna <- !is.na(g_sub$r_rna) &
      !is.na(g_sub$p_adj_rna) &
      g_sub$p_adj_rna < p_cut

    # FP filter within same genes (via tf_links_mode)
    genes_pass_fp <- character(0)
    if (nrow(tf_links_mode) > 0L) {
      pass_fp <- !is.na(tf_links_mode$r_fp) &
        !is.na(tf_links_mode$p_adj_fp) &
        tf_links_mode$p_adj_fp < p_cut

      # r_sign == "pos"
      pass_rna <- pass_rna & g_sub$r_rna > r_cut
      pass_fp  <- pass_fp  & tf_links_mode$r_fp > r_cut

      genes_pass_fp <- sort(unique(tf_links_mode$gene_key[pass_fp]))
    } else {
      pass_rna <- pass_rna & g_sub$r_rna > r_cut
    }

    genes_pass_rna <- sort(g_sub$gene_key[pass_rna])

    # Must pass BOTH FP + RNA corr thresholds
    genes_pass <- intersect(genes_pass_rna, genes_pass_fp)

    genes_all  <- genes_mode_all
    genes_fail <- setdiff(genes_all, genes_pass)

    n_draw <- min(length(genes_pass), length(bg_genes))
    genes_rand <- if (n_draw > 0L) {
      sort(sample(bg_genes, n_draw, replace = FALSE))
    } else {
      character(0)
    }

    list(
      passing     = genes_pass,
      non_passing = genes_fail,
      random      = genes_rand
    )
  }

  grid <- expand.grid(
    r_cut = r_cuts,
    p_cut = p_cuts,
    stringsAsFactors = FALSE
  )

  plot_list <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    r_cut <- grid$r_cut[i]
    p_cut <- grid$p_cut[i]

    sets_i <- get_sets_one_cut(r_cut, p_cut)

    pass_i <- intersect(sets_i$passing, genes_univ)
    non_i  <- intersect(sets_i$non_passing, genes_univ)
    rand_i <- intersect(sets_i$random, genes_univ)

    ko_pass <- ko_univ[ko_univ$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_univ[ko_univ$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_univ[ko_univ$gene %in% rand_i, , drop = FALSE]

    plot_list[[i]] <- dplyr::bind_rows(
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Predicted_regulated",
        gene    = ko_pass$gene,
        log2FC  = ko_pass$log2FC
      ),
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Predicted_nonregulated",
        gene    = ko_non$gene,
        log2FC  = ko_non$log2FC
      ),
      tibble::tibble(
        r_cut   = r_cut,
        p_cut   = p_cut,
        group   = "Random_background",
        gene    = ko_rand$gene,
        log2FC  = ko_rand$log2FC
      )
    )
  }

  plot_df <- dplyr::bind_rows(plot_list)

  r_levels <- r_cuts
  p_levels <- sort(unique(p_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_cut,
        levels = r_levels,
        labels = vapply(r_levels, function(rc) {
          paste0("r_fp > ", rc, " & r_rna > ", rc)
        }, character(1))
      ),
      p_lab = factor(
        p_cut,
        levels = p_levels,
        labels = scales::label_scientific(digits = 1)(p_levels)
      ),
      group = factor(
        group,
        levels = c(
          "Predicted_regulated",
          "Random_background",
          "Predicted_nonregulated"
        )
      )
    )

  grid_type_tag <- "FP_RNA_corr"

  file_stub <- sprintf(
    "%s_%s_%s_%s_%s_%s_%s_%s_%s_summary_violin_plots",
    tf_cell,
    tf_line,
    tf,
    ko_type,
    mode,
    grid_type_tag,
    r_sign_tag,
    fp_filter_tag,
    rna_filter_tag
  )

  plot_ko_grid_core(
    tf           = tf,
    mode         = mode,
    method_label = "FP & RNA correlated (r_fp>0 & r_rna>0)",
    plot_df      = plot_df,
    p_xlab       = "shared p_adj cutoff",
    file_stub    = file_stub,
    out_dir      = out_dir,
    test         = test,
    alternative  = alternative,
    verbose      = verbose
  )
}



ko_tbl_HNF1A <- HNF1A_KO %>%
  dplyr::transmute(
    gene   = gene_symbol,
    log2FC = log2fc
  ) %>%
  dplyr::filter(!is.na(gene), !is.na(log2FC))


run_all_ko_grids_for_tf <- function(tf,
                                    ko_tbl,
                                    tf_gene_links,
                                    atac_gene_links,
                                    strict_rna,
                                    tf_tfbs,
                                    ko_type        = "KO",
                                    tf_links_scope = "mode",
                                    mode           = "canonical",
                                    r_sign         = "pos",
                                    out_dir        = ".") {

  # Common cutoffs
  r_cuts <- c(0, 0.1, 0.3, 0.5, 0.7)
  p_cuts <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)

  message("=== Running FP-only grids for ", tf, " ===")
  # FP-only grids, but explore all fp_sig / rna_sig combinations
  for (fp_filter in c(FALSE, TRUE)) {
    for (rna_filter in c(FALSE, TRUE)) {
      plot_fp_ko_grid(
        tf                = tf,
        ko_tbl            = ko_tbl,
        tf_gene_links     = tf_gene_links,
        atac_gene_links   = atac_gene_links,
        tf_tfbs           = tf_tfbs,
        mode              = mode,
        r_fp_cuts         = r_cuts,
        p_fp_cuts         = p_cuts,
        r_sign            = r_sign,
        ko_type           = ko_type,
        tf_links_scope    = tf_links_scope,
        filter_by_fp_sig  = fp_filter,
        filter_by_rna_sig = rna_filter,
        out_dir           = out_dir
      )
    }
  }

  message("=== Running RNA-only grids for ", tf, " ===")
  # RNA-only grids, also all fp_sig / rna_sig combinations
  for (rna_filter in c(FALSE, TRUE)) {
    for (fp_filter in c(FALSE)) {
      plot_rna_ko_grid(
        tf                = tf,
        ko_tbl            = ko_tbl,
        tf_gene_links     = tf_gene_links,
        strict_rna        = strict_rna,
        tf_tfbs           = tf_tfbs,
        mode              = mode,
        r_rna_cuts        = r_cuts,
        p_rna_cuts        = p_cuts,
        cor_method        = "pearson",
        r_sign            = r_sign,
        ko_type           = ko_type,
        tf_links_scope    = tf_links_scope,
        filter_by_fp_sig  = fp_filter,
        filter_by_rna_sig = rna_filter,
        out_dir           = out_dir
      )
    }
  }

  message("=== Running FP+RNA combined grids for ", tf, " ===")
  # FP + RNA: all 4 combinations of fp_sig / rna_sig (unchanged)
  for (fp_filter in c(FALSE, TRUE)) {
    for (rna_filter in c(FALSE, TRUE)) {
      plot_fp_rna_ko_grid(
        tf                = tf,
        ko_tbl            = ko_tbl,
        tf_gene_links     = tf_gene_links,
        strict_rna        = strict_rna,
        tf_tfbs           = tf_tfbs,
        mode              = mode,
        r_cuts            = r_cuts,
        p_cuts            = p_cuts,
        cor_method        = "pearson",
        r_sign            = r_sign,
        ko_type           = ko_type,
        tf_links_scope    = tf_links_scope,
        filter_by_fp_sig  = fp_filter,
        filter_by_rna_sig = rna_filter,
        out_dir           = out_dir
      )
    }
  }

  invisible(NULL)
}


run_all_ko_grids_for_tf(
  tf              = "HNF1A",
  ko_tbl          = ko_tbl_HNF1A,
  tf_gene_links   = tf_gene_links_gh,
  atac_gene_links = atac_gene_links_gh,
  strict_rna      = strict_rna,
  tf_tfbs         = tf_tfbs_HNF1A,
  ko_type         = "KO",
  tf_links_scope  = "mode",
  mode            = "all",
  r_sign          = "pos",
  out_dir         = ko_dir
)

run_all_ko_grids_for_tf(
  tf              = "HNF1A",
  ko_tbl          = ko_tbl_HNF1A,
  tf_gene_links   = tf_gene_links_gh,
  atac_gene_links = atac_gene_links_gh,
  strict_rna      = strict_rna,
  tf_tfbs         = tf_tfbs_HNF1A,
  ko_type         = "KO",
  tf_links_scope  = "mode",
  mode            = "canonical",
  r_sign          = "pos",
  out_dir         = ko_dir
)

run_all_ko_grids_for_tf(
  tf              = "HNF4A",
  ko_tbl          = ko_tbl_HNF4A,
  tf_gene_links   = tf_gene_links_gh,
  atac_gene_links = atac_gene_links_gh,
  strict_rna      = strict_rna,
  tf_tfbs         = tf_tfbs_HNF4A,
  ko_type         = "KO",
  tf_links_scope  = "mode",
  mode            = "all",
  r_sign          = "pos",
  out_dir         = ko_dir
)

run_all_ko_grids_for_tf(
  tf              = "HNF4A",
  ko_tbl          = ko_tbl_HNF4A,
  tf_gene_links   = tf_gene_links_gh,
  atac_gene_links = atac_gene_links_gh,
  strict_rna      = strict_rna,
  tf_tfbs         = tf_tfbs_HNF4A,
  ko_type         = "KO",
  tf_links_scope  = "mode",
  mode            = "canonical",
  r_sign          = "pos",
  out_dir         = ko_dir
)


