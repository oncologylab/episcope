# -------------------------------------------------------------------
# New: ATAC-correlation-based KO split-violin grid (parallel to FP)
# -------------------------------------------------------------------
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
                              verbose = TRUE) {

  mode        <- match.arg(mode)
  test        <- match.arg(test)
  alternative <- match.arg(alternative)

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

  # ---------------- TF links for this TF + its TFBS ------------------
  tf_links0 <- tf_gene_links[tf_gene_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  if (nrow(tf_links0) == 0L) {
    stop("No tf_gene_links rows found for supplied tf_tfbs (ATAC grid).")
  }

  # Universe of genes (as in FP grid): TF-linked genes × KO table
  genes_univ <- intersect(unique(tf_links0$gene_key), ko_tbl$gene)
  ko_tbl <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  n_univ <- length(genes_univ)
  if (n_univ == 0L) {
    stop("No overlapping genes between KO table and tf_gene_links (ATAC grid).")
  }

  # Background gene pool for random (TF-unbound peaks), unchanged
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---------------- Left join ATAC correlations ----------------------
  # This effectively "filters atac_gene_links in the same manner"
  # by restricting to (atac_peak, gene_key) pairs present in tf_gene_links.
  tf_links_atac <- dplyr::left_join(
    tf_links0,
    atac_gene_links,
    by = c("atac_peak", "gene_key")
  )
  # Expect columns: n_atac, r_atac, p_atac, p_adj_atac from atac_gene_links.

  # Helper: for one (r_atac, p_atac) pair, get passing / non-passing / random genes
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
      abs(tf_links_curr$r_atac) > r_cut &
      tf_links_curr$p_adj_atac   < p_cut

    genes_all  <- sort(unique(tf_links_curr$gene_key))
    genes_pass <- sort(unique(tf_links_curr$gene_key[pass_rows]))
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

  # ---------------- Build grid over (r_atac, p_atac) -----------------
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

    ko_pass <- ko_tbl[ko_tbl$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_tbl[ko_tbl$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_tbl[ko_tbl$gene %in% rand_i, , drop = FALSE]

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

  # Factor labels for facets and x-axis (now r_atac / p_atac)
  r_atac_levels <- r_atac_cuts
  p_atac_levels <- sort(unique(p_atac_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_atac_cut,
        levels = r_atac_levels,
        labels = paste0("|r_atac| > ", r_atac_levels)
      ),
      p_lab = factor(
        p_atac_cut,
        levels = p_atac_levels,
        labels = scales::label_scientific(digits = 1)(p_atac_levels)
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
    dplyr::group_by(r_atac_cut, p_atac_cut) %>%
    dplyr::summarise(
      y = max(log2FC, na.rm = TRUE),
      .groups = "drop"
    )

  # vs Random
  df_rand_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Random_background"),
                         , drop = FALSE]

  pval_rand <- df_rand_all %>%
    dplyr::group_by(r_atac_cut, p_atac_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_rand <- dplyr::left_join(pval_rand, y_pos,
                                 by = c("r_atac_cut", "p_atac_cut")) %>%
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
        r_atac_cut,
        levels = r_atac_levels,
        labels = paste0("|r_atac| > ", r_atac_levels)
      ),
      p_lab = factor(
        p_atac_cut,
        levels = p_atac_levels,
        labels = scales::label_scientific(digits = 1)(p_atac_levels)
      )
    )

  # vs Non-regulated
  df_non_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Predicted_nonregulated"),
                        , drop = FALSE]

  pval_non <- df_non_all %>%
    dplyr::group_by(r_atac_cut, p_atac_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_non <- dplyr::left_join(pval_non, y_pos,
                                by = c("r_atac_cut", "p_atac_cut")) %>%
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
        r_atac_cut,
        levels = r_atac_levels,
        labels = paste0("|r_atac| > ", r_atac_levels)
      ),
      p_lab = factor(
        p_atac_cut,
        levels = p_atac_levels,
        labels = scales::label_scientific(digits = 1)(p_atac_levels)
      )
    )

  # -------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------
  if (verbose) {
    message("=== ATAC diagnostics for TF = ", tf,
            " (mode = ", mode, ", tf_filtered = TRUE) ===")
    message("Universe size (n_univ): ", n_univ)
    message("Test: ", test, ", alternative = ", alternative)

    # group counts per facet (first few rows)
    message("Random vs predicted (ATAC): group counts per facet (head):")
    counts_rand <- df_rand_all %>%
      dplyr::group_by(r_atac_cut, p_atac_cut, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_rand, 10))

    message("Predicted vs non-regulated (ATAC): group counts per facet (head):")
    counts_non <- df_non_all %>%
      dplyr::group_by(r_atac_cut, p_atac_cut, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_non, 10))

    if (all(is.na(pval_rand$p_val))) {
      message("Random vs predicted (ATAC): all p-values are NA.")
    } else {
      message("Random vs predicted (ATAC): ",
              "n facets with data = ", sum(!is.na(pval_rand$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_rand$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_rand$p_val, na.rm = TRUE), 3))
      ord_r <- order(pval_rand$p_val, na.last = NA)
      message("Top 5 smallest p (ATAC, random vs predicted):")
      print(utils::head(pval_rand[ord_r, , drop = FALSE], 5))
    }

    if (all(is.na(pval_non$p_val))) {
      message("Predicted vs non-regulated (ATAC): all p-values are NA.")
    } else {
      message("Predicted vs non-regulated (ATAC): ",
              "n facets with data = ", sum(!is.na(pval_non$p_val)),
              ", n facets with p < 0.05 = ",
              sum(pval_non$p_val < 0.05, na.rm = TRUE),
              ", min p = ",
              signif(min(pval_non$p_val, na.rm = TRUE), 3))
      ord_n <- order(pval_non$p_val, na.last = NA)
      message("Top 5 smallest p (ATAC, predicted vs non-regulated):")
      print(utils::head(pval_non[ord_n, , drop = FALSE], 5))
    }
  }

  # -------------------------------------------------------------------
  # Shared palette (same as FP)
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
    scale_x_discrete(name = "p_adj_atac cutoff") +
    scale_y_continuous(name = sprintf("log2FC (%s KO vs Ctrl)", tf)) +
    ggtitle(
      sprintf(
        "%s knockout: ATAC-correlated predicted regulated vs random background genes\n(mode = %s; n = %d genes in universe)",
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
    scale_x_discrete(name = "p_adj_atac cutoff") +
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
    dplyr::group_by(r_atac_cut, p_atac_cut, group) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r_lab = factor(
        r_atac_cut,
        levels = r_atac_levels,
        labels = paste0("|r_atac| > ", r_atac_levels)
      ),
      p_lab = factor(
        p_atac_cut,
        levels = p_atac_levels,
        labels = scales::label_scientific(digits = 1)(p_atac_levels)
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
    scale_x_discrete(name = "p_adj_atac cutoff") +
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
      "%s_KO_ATACcorr_tfFiltered_%s_grid_split_violin_vsRand_vsNon.pdf",
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


# HNF1A
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
  verbose         = TRUE,
  out_dir         = ko_dir
)

# HNF4A
plot_atac_ko_grid(
  tf              = "HNF4A",
  ko_tbl          = ko_tbl_HNF4A,
  tf_gene_links   = tf_gene_links,
  atac_gene_links = atac_gene_links,
  tf_tfbs         = tf_tfbs_HNF4A,
  mode            = "all",
  r_atac_cuts     = c(0, 0.1, 0.3, 0.5, 0.7),
  p_atac_cuts     = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
  test            = "wilcox",
  alternative     = "two.sided",
  verbose         = TRUE,
  out_dir         = ko_dir
)

