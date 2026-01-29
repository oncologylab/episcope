library(magrittr)
# -------------------------------------------------------------------
# Build KO ground truth (Down / Unchanged / Up) and plot
# -------------------------------------------------------------------
build_ko_truth_and_plot <- function(ko_tbl,
                                    tf,
                                    gene_cols  = c("gene_symbol", "symbol", "gene"),
                                    log2fc_cols = c("log2fc", "log2FoldChange", "log2fc_tf"),
                                    padj_cols   = c("p_adj", "padj", "padj_tf"),
                                    lfc_down    = -1,
                                    lfc_up      = 1,
                                    lfc_unch    = 0.25,
                                    padj_sig    = 0.05,
                                    padj_ns     = 0.5,
                                    out_dir     = ".") {

  # --- pick columns heuristically -----------------------------------
  gene_col <- gene_cols[gene_cols %in% names(ko_tbl)][1]
  if (is.na(gene_col)) {
    stop("No gene column found in KO table among: ",
         paste(gene_cols, collapse = ", "))
  }

  log2fc_col <- log2fc_cols[log2fc_cols %in% names(ko_tbl)][1]
  if (is.na(log2fc_col)) {
    stop("No log2FC column found in KO table among: ",
         paste(log2fc_cols, collapse = ", "))
  }

  padj_col <- padj_cols[padj_cols %in% names(ko_tbl)][1]
  if (is.na(padj_col)) {
    stop("No padj column found in KO table among: ",
         paste(padj_cols, collapse = ", "))
  }

  df <- ko_tbl %>%
    dplyr::transmute(
      gene   = .data[[gene_col]],
      log2fc = as.numeric(.data[[log2fc_col]]),
      padj   = as.numeric(.data[[padj_col]])
    ) %>%
    dplyr::filter(!is.na(gene), nzchar(gene), !is.na(log2fc), !is.na(padj))

  if (!nrow(df)) {
    stop("After filtering missing values, KO table is empty.")
  }

  # --- per-row group ------------------------------------------------
  df <- df %>%
    dplyr::mutate(
      ko_group_row = dplyr::case_when(
        log2fc <= lfc_down & padj <= padj_sig ~ "Down",
        log2fc >= lfc_up   & padj <= padj_sig ~ "Up",
        abs(log2fc) < lfc_unch & padj >= padj_ns ~ "Unchanged",
        TRUE ~ "Other"
      )
    )

  # --- collapse to per-gene truth -----------------------------------
  truth <- df %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      ko_group = dplyr::case_when(
        any(ko_group_row == "Down")      ~ "Down",
        any(ko_group_row == "Up")        ~ "Up",
        any(ko_group_row == "Unchanged") ~ "Unchanged",
        TRUE                             ~ "Other"
      ),
      log2fc = stats::median(log2fc, na.rm = TRUE),
      padj   = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      ko_group = factor(ko_group, levels = c("Down", "Unchanged", "Up", "Other"))
    )

  # --- diagnostics ---------------------------------------------------
  cat("=== KO ground truth for TF =", tf, "===\n")
  print(table(truth$ko_group, useNA = "ifany"))

  # --- plot (Down / Unchanged / Up only) -----------------------------
  plot_df <- truth %>%
    dplyr::filter(ko_group %in% c("Down", "Unchanged", "Up"))

  if (!nrow(plot_df)) {
    warning("No genes in Down/Unchanged/Up groups; skipping KO truth plot.")
    return(list(truth_tbl = truth, outfile = NA_character_))
  }

  p_ko <- ggplot(plot_df, aes(x = ko_group, y = log2fc, fill = ko_group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.4, colour = "black") +
    ggplot2::geom_boxplot(width = 0.2, alpha = 0.7, outlier.size = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::scale_fill_manual(
      values = c(
        Down      = "#d73027",
        Unchanged = "#969696",
        Up        = "#1a9850"
      )
    ) +
    ggplot2::labs(
      title = sprintf("%s KO: ground truth groups (log2FC vs Ctrl)", tf),
      x     = "KO ground truth group",
      y     = "log2FC (KO vs Ctrl)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      axis.title.x    = ggplot2::element_text(face = "bold"),
      axis.title.y    = ggplot2::element_text(face = "bold"),
      axis.text.x     = ggplot2::element_text(face = "bold"),
      axis.text.y     = ggplot2::element_text(face = "bold")
    )

  outfile <- file.path(
    out_dir,
    sprintf("%s_KO_ground_truth_groups_violin.pdf", tf)
  )

  ggplot2::ggsave(
    filename = outfile,
    plot     = p_ko,
    width    = 6,
    height   = 4,
    dpi      = 600
  )

  cat("KO truth plot saved to:", outfile, "\n")

  invisible(list(
    truth_tbl = truth,
    outfile   = outfile
  ))
}

# HNF4A KO ground truth
hnf4a_gt <- build_ko_truth_and_plot(
  ko_tbl = HNF4A_KO,
  tf     = "HNF4A",
  out_dir = ko_dir
)
ko_truth_HNF4A <- hnf4a_gt$truth_tbl

# HNF1A KO ground truth (you may want to subset to a specific dataset first)
# e.g., one KO comparison:
# HNF1A_KO_single <- dplyr::filter(HNF1A_KO, dataset_id == "GSEXXXXX", tf_parsed == "HNF1A")
hnf1a_gt <- build_ko_truth_and_plot(
  ko_tbl = HNF1A_KO,
  tf     = "HNF1A",
  out_dir = ko_dir
)
ko_truth_HNF1A <- hnf1a_gt$truth_tbl

