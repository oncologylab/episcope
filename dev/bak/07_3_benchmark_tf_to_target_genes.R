# -------------------------------------------------------------------
# New: RNA-correlation-based KO split-violin grid (parallel to FP/ATAC)
# -------------------------------------------------------------------
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
                             verbose = TRUE) {

  mode        <- match.arg(mode)
  cor_method  <- match.arg(cor_method)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)

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

  # small helper to run the chosen test on one facet
  run_test <- function(df_sub, test, alternative) {
    df_sub <- df_sub[!is.na(df_sub$log2FC), , drop = FALSE]
    if (nrow(df_sub) == 0L) {
      return(NA_real_)
    }

    # Drop unused levels to avoid 3-level issues
    gtab <- table(droplevels(df_sub$group))
    if (length(gtab) != 2L || any(gtab == 0L)) {
      return(NA_real_)
    }

    if (test == "wilcox") {
      stats::wilcox.test(
        log2FC ~ group,
        data        = df_sub,
        alternative = alternative
      )$p.value
    } else { # test == "t"
      stats::t.test(
        log2FC ~ group,
        data        = df_sub,
        alternative = alternative
      )$p.value
    }
  }

  # ---------------- Extract TF-linked genes & background --------------------
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (RNA grid).")
  }

  # genes linked via TFBS peaks
  genes_tf <- sort(unique(tf_links0$gene_key))

  # background: peaks NOT in tf_tfbs, same logic as FP/ATAC
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  # Universe: genes with KO log2FC and TF links
  genes_univ0 <- intersect(genes_tf, ko_tbl$gene)

  if (!length(genes_univ0)) {
    stop("No overlapping genes between KO table and tf_gene_links (RNA grid).")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---------------- Prepare RNA expression matrix & TF expression ----------
  sample_cols <- setdiff(names(strict_rna),
                         c("ensembl_gene_id", "HGNC"))

  if (!length(sample_cols)) {
    stop("strict_rna has no sample columns.")
  }

  rna_expr_mat <- as.matrix(strict_rna[, sample_cols, drop = FALSE])

  # TF expression: prefer HGNC == tf, else ensembl_gene_id == tf
  idx_tf <- which(strict_rna$HGNC == tf)
  if (!length(idx_tf)) {
    idx_tf <- which(strict_rna$ensembl_gene_id == tf)
  }
  if (!length(idx_tf)) {
    stop(sprintf("TF '%s' not found in strict_rna$HGNC or ensembl_gene_id.", tf))
  }
  idx_tf <- idx_tf[1L]
  tf_expr <- as.numeric(rna_expr_mat[idx_tf, ])

  # ---------------- Map universe genes to RNA rows -------------------------
  map_univ <- map_gene_keys_to_rna(genes_univ0, strict_rna)
  if (!nrow(map_univ)) {
    stop("No universe genes found in strict_rna (RNA grid).")
  }

  genes_univ <- map_univ$gene_key
  ko_tbl     <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  n_univ     <- length(genes_univ)

  if (n_univ == 0L) {
    stop("After matching KO + RNA, no genes remain in universe (RNA grid).")
  }

  # ---------------- Compute TF–gene RNA correlations -----------------------
  n_samples <- length(sample_cols)

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

  rna_corr_full <- dplyr::bind_rows(corr_list)

  rna_corr_full <- rna_corr_full %>%
    dplyr::mutate(
      p_adj_rna = stats::p.adjust(p_rna, method = "BH")
    )

  # Final universe = genes with KO + RNA corr
  genes_univ <- intersect(genes_univ, rna_corr_full$gene_key)
  ko_tbl     <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  rna_corr_full <- rna_corr_full[rna_corr_full$gene_key %in% genes_univ,
                                 , drop = FALSE]
  n_univ     <- length(genes_univ)

  if (n_univ == 0L) {
    stop("After RNA correlation filtering, no genes remain in universe (RNA grid).")
  }

  # ---------------- Precompute mode-specific gene universe -----------------
  if (mode == "canonical") {
    motifs_vec <- toupper(as.character(tf_links0$motifs))
    tf_upper   <- toupper(tf)
    keep_motif <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)
    tf_links_mode <- tf_links0[keep_motif, , drop = FALSE]
  } else {
    tf_links_mode <- tf_links0
  }

  genes_mode_all <- intersect(
    sort(unique(tf_links_mode$gene_key)),
    genes_univ
  )

  # Helper: for one (r_rna, p_adj_rna) pair, get passing / non-passing / random genes
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
      abs(g_sub$r_rna) > r_cut &
      g_sub$p_adj_rna   < p_cut

    genes_pass <- sort(g_sub$gene_key[pass_rows])
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

  # ---------------- Build grid over (r_rna, p_rna) -------------------------
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

    ko_pass <- ko_tbl[ko_tbl$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_tbl[ko_tbl$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_tbl[ko_tbl$gene %in% rand_i, , drop = FALSE]

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

  # Factor labels for facets and x-axis (now r_rna / p_adj_rna)
  r_rna_levels <- r_rna_cuts
  p_rna_levels <- sort(unique(p_rna_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = paste0("|r_rna| > ", r_rna_levels)
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
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

  # -------------------------------------------------------------------
  # p-values for row 1 and row 2
  # -------------------------------------------------------------------
  global_range <- range(plot_df$log2FC, na.rm = TRUE)
  offset       <- 0.05 * diff(global_range)

  y_pos <- plot_df %>%
    dplyr::group_by(r_rna_cut, p_rna_cut) %>%
    dplyr::summarise(
      y = max(log2FC, na.rm = TRUE),
      .groups = "drop"
    )

  # vs Random
  df_rand_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Random_background"),
                         , drop = FALSE]

  pval_rand <- df_rand_all %>%
    dplyr::group_by(r_rna_cut, p_rna_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_rand <- dplyr::left_join(pval_rand, y_pos,
                                 by = c("r_rna_cut", "p_rna_cut")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(p_val)  ~ "N.S.",
        p_val >= 0.05 ~ "N.S.",
        TRUE          ~ paste0(
          "p = ",
          scales::label_scientific(digits = 2)(p_val)
        )
      ),
      y     = y + offset,
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = paste0("|r_rna| > ", r_rna_levels)
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
      )
    )

  # vs Non-regulated
  df_non_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Predicted_nonregulated"),
                        , drop = FALSE]

  pval_non <- df_non_all %>%
    dplyr::group_by(r_rna_cut, p_rna_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_non <- dplyr::left_join(pval_non, y_pos,
                                by = c("r_rna_cut", "p_rna_cut")) %>%
    dplyr::mutate(
      label = dplyr::case_when(
        is.na(p_val)  ~ "N.S.",
        p_val >= 0.05 ~ "N.S.",
        TRUE          ~ paste0(
          "p = ",
          scales::label_scientific(digits = 2)(p_val)
        )
      ),
      y     = y + offset,
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = paste0("|r_rna| > ", r_rna_levels)
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
      )
    )

  # -------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------
  if (verbose) {
    message("=== RNA diagnostics for TF = ", tf,
            " (mode = ", mode, ", tf_filtered = TRUE) ===")
    message("Universe size (n_univ): ", n_univ)
    message("Correlation method: ", cor_method)
    message("Test: ", test, ", alternative = ", alternative)

    # group counts per facet (first few rows)
    message("Random vs predicted (RNA): group counts per facet (head):")
    counts_rand <- df_rand_all %>%
      dplyr::group_by(r_rna_cut, p_rna_cut, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_rand, 10))

    message("Predicted vs non-regulated (RNA): group counts per facet (head):")
    counts_non <- df_non_all %>%
      dplyr::group_by(r_rna_cut, p_rna_cut, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_non, 10))

    if (all(is.na(pval_rand$p_val))) {
      message("Random vs predicted (RNA): all p-values are NA.")
    } else {
      message("Random vs predicted (RNA): ",
              "n facets with data = ", sum(!is.na(pval_rand$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_rand$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_rand$p_val, na.rm = TRUE), 3))
      ord_r <- order(pval_rand$p_val, na.last = NA)
      message("Top 5 smallest p (RNA, random vs predicted):")
      print(utils::head(pval_rand[ord_r, , drop = FALSE], 5))
    }

    if (all(is.na(pval_non$p_val))) {
      message("Predicted vs non-regulated (RNA): all p-values are NA.")
    } else {
      message("Predicted vs non-regulated (RNA): ",
              "n facets with data = ", sum(!is.na(pval_non$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_non$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_non$p_val, na.rm = TRUE), 3))
      ord_n <- order(pval_non$p_val, na.last = NA)
      message("Top 5 smallest p (RNA, predicted vs non-regulated):")
      print(utils::head(pval_non[ord_n, , drop = FALSE], 5))
    }
  }

  # -------------------------------------------------------------------
  # Shared palette (same as FP/ATAC)
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
      width = 0.2,
      alpha = 0.6,
      show.legend = FALSE
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
      values = fill_vals,
      name   = "Group",
      labels = c("Predicted regulated", "Random background")
    ) +
    scale_x_discrete(name = "p_adj_rna cutoff") +
    scale_y_continuous(name = sprintf("log2FC (%s KO vs Ctrl)", tf)) +
    ggtitle(
      sprintf(
        "%s knockout: RNA TF–gene correlated predicted regulated vs random background genes\n(mode = %s; n = %d genes in universe)",
        tf,
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
      width = 0.2,
      alpha = 0.6,
      show.legend = FALSE
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
      values = fill_vals,
      name   = "Group",
      labels = c("Predicted regulated", "Predicted non-regulated")
    ) +
    scale_x_discrete(name = "p_adj_rna cutoff") +
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
  # Row 3: counts (predicted vs non-regulated)
  # -------------------------------------------------------------------
  counts_df <- df_non %>%
    dplyr::group_by(r_rna_cut, p_rna_cut, group) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r_lab = factor(
        r_rna_cut,
        levels = r_rna_levels,
        labels = paste0("|r_rna| > ", r_rna_levels)
      ),
      p_lab = factor(
        p_rna_cut,
        levels = p_rna_levels,
        labels = scales::label_scientific(digits = 1)(p_rna_levels)
      )
    )

  p_counts <- ggplot(
    counts_df,
    aes(x = p_lab, y = n_genes, fill = group)
  ) +
    geom_col(
      position  = position_dodge(width = 0.7),
      width     = 0.6,
      colour    = "black",
      linewidth = 0.2
    ) +
    facet_wrap(~ r_lab, nrow = 1) +
    scale_fill_manual(
      values = fill_vals,
      name   = "Group",
      labels = c("Predicted regulated", "Predicted non-regulated")
    ) +
    scale_x_discrete(name = "p_adj_rna cutoff") +
    scale_y_continuous(name = "Number of unique genes") +
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
  # Combine rows and save
  # -------------------------------------------------------------------
  combined <- p_rand / p_non / p_counts + patchwork::plot_layout(heights = c(3, 3, 1))

  outfile <- file.path(
    out_dir,
    sprintf(
      "%s_KO_rna_corr_tfFiltered_%s_grid_split_violin_vsRand_vsNon.pdf",
      tf,
      mode
    )
  )

  ggplot2::ggsave(
    filename = outfile,
    plot     = combined,
    width    = 18,
    height   = 10,
    dpi      = 600
  )

  invisible(outfile)
}

# HNF1A RNA-corr KO grid
plot_rna_ko_grid(
  tf            = "HNF1A",
  ko_tbl        = ko_tbl,         # HNF1A_KO table you already built
  tf_gene_links = tf_gene_links,  # tf_gene_links_atac_filtered
  strict_rna    = strict_rna,
  tf_tfbs       = tf_tfbs_HNF1A,
  mode          = "all",
  r_rna_cuts    = c(0, 0.1, 0.3, 0.5, 0.7),
  p_rna_cuts    = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
  cor_method    = "pearson",
  test          = "wilcox",
  alternative   = "two.sided",
  verbose       = TRUE,
  out_dir       = ko_dir
)

# HNF4A RNA-corr KO grid
plot_rna_ko_grid(
  tf            = "HNF4A",
  ko_tbl        = ko_tbl_HNF4A,
  tf_gene_links = tf_gene_links,
  strict_rna    = strict_rna,
  tf_tfbs       = tf_tfbs_HNF4A,
  mode          = "all",
  r_rna_cuts    = c(0, 0.1, 0.3, 0.5, 0.7),
  p_rna_cuts    = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
  cor_method    = "pearson",
  test          = "wilcox",
  alternative   = "two.sided",
  verbose       = TRUE,
  out_dir       = ko_dir
)

