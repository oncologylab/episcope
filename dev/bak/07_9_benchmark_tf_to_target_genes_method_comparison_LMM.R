suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(plyr)
  library(cli)
})

base_dir_window <- "/data/homes/yl814/episcope_test/nutrient_stress/connect_tfs_to_target_genes"
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

## GeneHancer-based (still loaded, but not used in the new KO pipeline)
tf_gene_links_gh   <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_jaspar2024.csv"))
atac_gene_links_gh <- readr::read_csv(file.path(base_dir_window, "atac_gene_corr_full_jaspar2024.csv"))

## Cell-line specific GeneHancer-based GRN sets
## NO LONGER USED in the KO plotting pipeline
# tf_gene_links_gh_Panc1  <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_grn_set_Panc1.csv"))
# tf_gene_links_gh_AsPC1  <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_grn_set_AsPC1.csv"))
# tf_gene_links_gh_HPAFII <- readr::read_csv(file.path(base_dir_window, "fp_gene_corr_full_grn_set_HPAFII.csv"))

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
## (old RNA helpers left here for compatibility if needed elsewhere)
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
## 4) Attach FP (TFBS)–RNA correlation and TF–RNA correlation to pvals
## ==================================================================

attach_fp_r_corr_to_pvals <- function(tf,
                                      pval_tbl,
                                      tf_gene_links_gh,
                                      r_col_source = "r_fp",
                                      r_col_out = "r_fp",
                                      verbose = TRUE) {
  if (!all(c("gene", "fp_peak") %in% names(pval_tbl))) {
    cli::cli_abort("pval_tbl for TF {.val {tf}} must contain 'gene' and 'fp_peak'.")
  }
  if (!all(c("fp_peak", "gene_key", r_col_source, "tfs") %in% names(tf_gene_links_gh))) {
    cli::cli_abort("tf_gene_links_gh must contain 'fp_peak', 'gene_key', {.val {r_col_source}}, and 'tfs'.")
  }

  tf_links <- tf_gene_links_gh[tf_gene_links_gh$tfs == tf, , drop = FALSE]

  if (!nrow(tf_links)) {
    if (verbose) {
      cli::cli_warn("No tf_gene_links_gh rows found for TF {.val {tf}}; {.val {r_col_out}} will be NA.")
    }
    pval_tbl[[r_col_out]] <- NA_real_
    return(pval_tbl)
  }

  tf_links <- tf_links[!is.na(tf_links[[r_col_source]]), , drop = FALSE]

  if (!nrow(tf_links)) {
    if (verbose) {
      cli::cli_warn("TF {.val {tf}} has no non-NA {.val {r_col_source}} in tf_gene_links_gh; {.val {r_col_out}} will be NA.")
    }
    pval_tbl[[r_col_out]] <- NA_real_
    return(pval_tbl)
  }

  # Collapse to one r value per (fp_peak, gene_key) using max |r|
  tf_links_collapsed <- tf_links %>%
    dplyr::group_by(fp_peak, gene_key) %>%
    dplyr::summarise(
      idx          = which.max(abs(.data[[r_col_source]])),
      !!r_col_source := .data[[r_col_source]][idx],
      p_fp         = .data[["p_fp"]][idx],
      .groups      = "drop"
    )


  out <- pval_tbl %>%
    dplyr::left_join(
      tf_links_collapsed,
      by = c("fp_peak" = "fp_peak", "gene" = "gene_key")
    )

  if (r_col_source != r_col_out) {
    out <- dplyr::rename(out, !!r_col_out := .data[[r_col_source]])
  }

  if (verbose) {
    n_match <- sum(!is.na(out[[r_col_out]]))
    cli::cli_inform(
      "Attached {.val {r_col_out}} from tf_gene_links_gh to {.val {n_match}} / {.val {nrow(out)}} rows for TF {.val {tf}}."
    )
  }

  out
}

attach_tf_rna_corr_to_pvals <- function(tf,
                                        pval_tbl,
                                        strict_rna,
                                        r_col_out = "r_rna",
                                        cor_method = c("pearson", "spearman"),
                                        verbose = TRUE) {
  cor_method <- match.arg(cor_method)

  if (!"gene" %in% names(pval_tbl)) {
    cli::cli_abort("pval_tbl for TF {.val {tf}} must contain 'gene'.")
  }

  genes <- sort(unique(pval_tbl$gene))

  if (!length(genes)) {
    pval_tbl[[r_col_out]] <- NA_real_
    return(pval_tbl)
  }

  corr_tbl <- compute_rna_corr_for_tf(
    tf         = tf,
    genes      = genes,
    strict_rna = strict_rna,
    cor_method = cor_method
  )

  if (!nrow(corr_tbl)) {
    if (verbose) {
      cli::cli_warn("No TF–RNA correlation rows computed for TF {.val {tf}}; {.val {r_col_out}} will be NA.")
    }
    pval_tbl[[r_col_out]] <- NA_real_
    return(pval_tbl)
  }

  ## Keep r_rna, p_rna, p_adj_rna
  corr_tbl2 <- corr_tbl %>%
    dplyr::select(
      gene = gene_key,
      r_rna,
      p_rna,
      p_adj_rna
    )

  out <- pval_tbl %>%
    dplyr::left_join(corr_tbl2, by = "gene")

  if (r_col_out != "r_rna") {
    out <- dplyr::rename(out, !!r_col_out := .data[["r_rna"]])
  }

  if (verbose) {
    n_match <- sum(!is.na(out[[r_col_out]]))
    cli::cli_inform(
      "Attached TF–RNA correlation {.val {r_col_out}} to {.val {n_match}} / {.val {nrow(out)}} rows for TF {.val {tf}}."
    )
  }

  out
}

## ==================================================================
## 5) Build KO summary table for combined filters (LMM / FIXED / FIXED_INT)
##      – using pre-collapsed p-values + r_fp & r_rna filters
## ==================================================================

build_ko_summary_table_for_tf <- function(tf,
                                          ko_tbl,        # gene + log2FC
                                          pval_tbl,      # pre-collapsed: gene + ko_group + p_* + r_fp + r_rna + p_rna
                                          panel_type = c("lmm", "fixed", "fixed_int"),
                                          p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                                          verbose = TRUE,
                                          # correlation filters for Predicted group
                                          r_fp_col = "r_fp",
                                          r_rna_col = "r_rna",
                                          rna_p_col = "p_rna",
                                          fp_p_col = "p_fp",
                                          min_abs_fp_cor = 0.3,
                                          min_abs_rna_cor = 0.3,
                                          max_fp_p = 0.01,
                                          max_rna_p = 0.01,
                                          require_fp_cor = TRUE,
                                          require_rna_cor = TRUE,
                                          require_fp_p = TRUE,
                                          require_rna_p = TRUE) {

  panel_type <- match.arg(panel_type)

  if (verbose) {
    message("Building KO summary table for TF = ", tf,
            " (panel = ", panel_type, ", canonical only)")
  }

  tf_meta <- get_tf_cell_lineage(tf)
  tf_cell <- tf_meta$cell
  tf_line <- tf_meta$lineage

  # keep only well-defined KO groups
  if (!all(c("gene", "ko_group") %in% names(pval_tbl))) {
    stop("pval_tbl for TF ", tf,
         " must contain columns 'gene' and 'ko_group'.")
  }

  pdat <- pval_tbl %>%
    dplyr::filter(
      !is.na(ko_group),
      ko_group %in% c("Down", "Unchanged", "Up")
    )

  # basic gene universe: overlap with KO table
  genes_all <- sort(intersect(pdat$gene, ko_tbl$gene))
  if (!length(genes_all)) {
    stop("No overlapping genes between KO table and p-value table for TF ", tf, ".")
  }

  ## ---------------------------------------------------------------
  ## Correlation filters for Predicted group
  ##   - require |r_fp| >= min_abs_fp_cor
  ##   - require |r_rna| >= min_abs_rna_cor
  ##   - require p_rna < max_rna_p
  ## Only applied to Predicted genes; Non-predicted are the complement.
  ## ---------------------------------------------------------------

  genes_good_fp  <- genes_all
  genes_good_rna <- genes_all

  if (require_fp_cor && r_fp_col %in% names(pdat)) {
    r_vec <- pdat[[r_fp_col]]

    if (require_fp_p && fp_p_col %in% names(pdat)) {
      p_vec <- pdat[[fp_p_col]]
      genes_good_fp <- sort(unique(
        pdat$gene[
          !is.na(r_vec) & abs(r_vec) >= min_abs_fp_cor &
            !is.na(p_vec) & p_vec < max_fp_p
        ]
      ))
    } else {
      genes_good_fp <- sort(unique(
        pdat$gene[
          !is.na(r_vec) & abs(r_vec) >= min_abs_fp_cor
        ]
      ))
    }

    if (!length(genes_good_fp) && verbose) {
      if (require_fp_p && fp_p_col %in% names(pdat)) {
        message(
          "  [", tf, "] No genes pass |", r_fp_col,
          "| >= ", min_abs_fp_cor, " and ", fp_p_col,
          " < ", max_fp_p, " (FP correlation). Predicted set will be empty."
        )
      } else {
        message(
          "  [", tf, "] No genes pass |", r_fp_col,
          "| >= ", min_abs_fp_cor, " (FP correlation). Predicted set will be empty."
        )
      }
    }
  } else if (require_fp_cor && !(r_fp_col %in% names(pdat)) && verbose) {
    message(
      "  [", tf, "] Column '", r_fp_col,
      "' not found in pval_tbl; skipping FP-correlation filter."
    )
  }

  if (require_rna_cor && r_rna_col %in% names(pdat)) {
    r_vec2 <- pdat[[r_rna_col]]

    if (require_rna_p && rna_p_col %in% names(pdat)) {
      p_vec <- pdat[[rna_p_col]]
      genes_good_rna <- sort(unique(
        pdat$gene[
          !is.na(r_vec2) & abs(r_vec2) >= min_abs_rna_cor &
            !is.na(p_vec) & p_vec < max_rna_p
        ]
      ))
    } else {
      genes_good_rna <- sort(unique(
        pdat$gene[
          !is.na(r_vec2) & abs(r_vec2) >= min_abs_rna_cor
        ]
      ))
    }

    if (!length(genes_good_rna) && verbose) {
      if (require_rna_p && rna_p_col %in% names(pdat)) {
        message(
          "  [", tf, "] No genes pass |", r_rna_col,
          "| >= ", min_abs_rna_cor, " and ", rna_p_col,
          " < ", max_rna_p, " (TF–RNA correlation). Predicted set will be empty."
        )
      } else {
        message(
          "  [", tf, "] No genes pass |", r_rna_col,
          "| >= ", min_abs_rna_cor, " (TF–RNA correlation). Predicted set will be empty."
        )
      }
    }
  } else if (require_rna_cor && !(r_rna_col %in% names(pdat)) && verbose) {
    message(
      "  [", tf, "] Column '", r_rna_col,
      "' not found in pval_tbl; skipping TF–RNA correlation filter."
    )
  }

  # helper: per-method gene sets
  get_p_sets <- function(tbl, p_col, p_cut) {
    if (!p_col %in% names(tbl)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    # base universe for this method = genes with non-NA p_col
    base_genes <- sort(intersect(genes_all, tbl$gene[!is.na(tbl[[p_col]])]))
    if (!length(base_genes)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    sub <- tbl[tbl$gene %in% base_genes, , drop = FALSE]
    if (!nrow(sub)) {
      return(list(passing = character(0), non_passing = character(0)))
    }

    pass_rows  <- !is.na(sub[[p_col]]) & sub[[p_col]] < p_cut
    genes_pass <- sort(unique(sub$gene[pass_rows]))
    genes_all_m <- sort(unique(sub$gene))

    ## KEY: only keep genes with sufficiently strong FP & RNA correlation
    genes_pass <- intersect(genes_pass, genes_good_fp)
    genes_pass <- intersect(genes_pass, genes_good_rna)

    genes_fail <- setdiff(genes_all_m, genes_pass)

    list(passing = genes_pass, non_passing = genes_fail)
  }

  ## Method definitions per panel -------------------------------------

  if (panel_type == "lmm") {
    method_defs <- tibble::tibble(
      method_id  = c("LMM_FP",
                     "LMM_RNA",
                     "LMM_BOTH_FP",
                     "LMM_BOTH_RNA",
                     "LMM_BOTH_OVERALL"),
      p_col      = c("p_lmm_fp",
                     "p_lmm_rna",
                     "p_lmm_both_fp",
                     "p_lmm_both_rna",
                     "p_lmm_both_overall"),
      method_lab = factor(
        c(
          "LMM – FP only",
          "LMM – RNA only",
          "LMM – both (FP side)",
          "LMM – both (RNA side)",
          "LMM – both (overall)"
        ),
        levels = c(
          "LMM – FP only",
          "LMM – RNA only",
          "LMM – both (FP side)",
          "LMM – both (RNA side)",
          "LMM – both (overall)"
        )
      )
    )
  } else if (panel_type == "fixed") {
    method_defs <- tibble::tibble(
      method_id  = c("FIXED_FP",
                     "FIXED_RNA",
                     "FIXED_BOTH_FP",
                     "FIXED_BOTH_RNA",
                     "FIXED_BOTH_OVERALL"),
      p_col      = c("p_fixed_fp",
                     "p_fixed_rna",
                     "p_fixed_both_fp",
                     "p_fixed_both_rna",
                     "p_fixed_both_overall"),
      method_lab = factor(
        c(
          "Fixed – FP only",
          "Fixed – RNA only",
          "Fixed – both (FP side)",
          "Fixed – both (RNA side)",
          "Fixed – both (overall)"
        ),
        levels = c(
          "Fixed – FP only",
          "Fixed – RNA only",
          "Fixed – both (FP side)",
          "Fixed – both (RNA side)",
          "Fixed – both (overall)"
        )
      )
    )
  } else { # panel_type == "fixed_int"
    method_defs <- tibble::tibble(
      method_id  = c("FIXED_INT_FP",
                     "FIXED_INT_RNA",
                     "FIXED_INT_BOTH_FP",
                     "FIXED_INT_BOTH_RNA",
                     "FIXED_INT_BOTH_OVERALL"),
      p_col      = c("p_fixed_fp_int",
                     "p_fixed_rna_int",
                     "p_fixed_both_fp_int",
                     "p_fixed_both_rna_int",
                     # overall: still use non-int combined (no *_overall_int column)
                     "p_fixed_both_overall"),
      method_lab = factor(
        c(
          "Fixed+Int – FP only",
          "Fixed+Int – RNA only",
          "Fixed+Int – both (FP side)",
          "Fixed+Int – both (RNA side)",
          "Fixed+Int – both (overall)"
        ),
        levels = c(
          "Fixed+Int – FP only",
          "Fixed+Int – RNA only",
          "Fixed+Int – both (FP side)",
          "Fixed+Int – both (RNA side)",
          "Fixed+Int – both (overall)"
        )
      )
    )
  }

  ## Build long-format table ------------------------------------------

  res_list <- list()
  idx      <- 1L

  for (i in seq_len(nrow(method_defs))) {
    mid   <- method_defs$method_id[i]
    mlab  <- as.character(method_defs$method_lab[i])
    p_col <- method_defs$p_col[i]

    # base universe for this method (genes with any p defined)
    base_genes <- sort(intersect(genes_all, pdat$gene[!is.na(pdat[[p_col]])]))

    if (!length(base_genes)) {
      if (verbose) message("Method ", mlab, ": base universe empty (skipping).")
      next
    }

    for (pc in p_cuts) {
      sets <- get_p_sets(pdat, p_col = p_col, p_cut = pc)

      genes_pred <- intersect(sets$passing, ko_tbl$gene)
      genes_non  <- intersect(sets$non_passing, ko_tbl$gene)

      if (!length(genes_pred) && !length(genes_non)) next

      ko_pred <- ko_tbl[ko_tbl$gene %in% genes_pred, , drop = FALSE]
      ko_non  <- ko_tbl[ko_tbl$gene %in% genes_non,  , drop = FALSE]

      if (!nrow(ko_pred) && !nrow(ko_non)) next

      res_list[[idx]] <- dplyr::bind_rows(
        tibble::tibble(
          tf         = tf,
          cell       = tf_cell,
          lineage    = tf_line,
          mode       = "canonical",
          panel      = panel_type,
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
          mode       = "canonical",
          panel      = panel_type,
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
    stop("No data generated for TF ", tf,
         " (panel ", panel_type, ", canonical).")
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
## 6) Plotter: 3-row layout (violin + n_pred + % bins) for each PANEL
## ==================================================================

plot_ko_summary_violin_for_tf <- function(tf,
                                          ko_tbl,
                                          pval_tbl,
                                          panel_type = c("lmm", "fixed", "fixed_int"),
                                          p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                                          ko_type = "KO",
                                          out_dir = ".",
                                          verbose = TRUE,
                                          filename_suffix = "") {

  panel_type <- match.arg(panel_type)

  df <- build_ko_summary_table_for_tf(
    tf         = tf,
    ko_tbl     = ko_tbl,
    pval_tbl   = pval_tbl,
    panel_type = panel_type,
    p_cuts     = p_cuts,
    verbose    = verbose
  )

  tf_meta <- get_tf_cell_lineage(tf)
  tf_cell <- tf_meta$cell
  tf_line <- tf_meta$lineage

  panel_label <- switch(
    panel_type,
    lmm        = "LMM",
    fixed      = "Fixed",
    fixed_int  = "Fixed+Interaction"
  )

  ## 6A) Violin + p-values --------------------------------------------

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
      "%s %s in %s (%s) - canonical\n%s panel: p-value cutoffs",
      tf, ko_type, tf_cell, tf_line, panel_label
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

  ## 6C) Bar: number of predicted genes -------------------------------

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

  ## 6D) Bar: percent of log2FC bins per group ------------------------

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

  ## 6E) Combine rows and save ----------------------------------------

  combined <- p_violin / p_bar_n / p_bar_pct +
    patchwork::plot_layout(heights = c(3, 1.2, 2))

  outfile <- file.path(
    out_dir,
    sprintf(
      "%s_%s_%s_%s_canonical_%s_panel_pval_summary_violin%s.pdf",
      tf,
      ko_type,
      tf_cell,
      tf_line,
      toupper(panel_type),
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
## 7) KO-table helper
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
## 8) Batch runner for KO violin plots (using pre-collapsed lists)
## ==================================================================

run_all_ko_violin <- function(tfs,
                              ko_truth_list,
                              tf_pvals_lmm_collapsed,
                              tf_pvals_fixed_collapsed,
                              tf_pvals_fixed_int_collapsed,
                              strict_rna,
                              ko_type = "KO",
                              out_dir = ".",
                              p_cuts = c(5e-5, 5e-4, 5e-3, 5e-2, 5e-1),
                              verbose = TRUE,
                              filename_suffix = "") {

  # keep only TFs that are present in KO truth and all three pval lists
  tfs_use <- tfs[
    tfs %in% names(ko_truth_list) &
      tfs %in% names(tf_pvals_lmm_collapsed) &
      tfs %in% names(tf_pvals_fixed_collapsed) &
      tfs %in% names(tf_pvals_fixed_int_collapsed)
  ]

  if (!length(tfs_use)) {
    stop("run_all_ko_violin(): no TFs present in all required lists.")
  }

  if (verbose) {
    message("Running KO violin summaries for TFs (canonical, LMM/FIXED/FIXED_INT): ",
            paste(tfs_use, collapse = ", "))
  }

  res <- lapply(tfs_use, function(tf_symbol) {
    if (verbose) {
      message("  -> TF = ", tf_symbol)
    }

    ko_truth_tbl <- ko_truth_list[[tf_symbol]]
    if (!is.data.frame(ko_truth_tbl) || !nrow(ko_truth_tbl)) {
      message("    Skipping TF ", tf_symbol, ": KO truth table empty.")
      return(NULL)
    }

    # compact KO table (gene, log2FC)
    ko_tbl <- make_ko_tbl(ko_truth_tbl, tf_label = tf_symbol)

    # pre-collapsed p-values + attach FP–RNA correlation (r_fp)
    # and TF–RNA correlation (r_rna, p_rna, p_adj_rna)
    pdat_lmm <- tf_pvals_lmm_collapsed[[tf_symbol]] %>%
      attach_fp_r_corr_to_pvals(
        tf               = tf_symbol,
        pval_tbl         = .,
        tf_gene_links_gh = tf_gene_links_gh,
        r_col_source     = "r_fp",
        r_col_out        = "r_fp",
        verbose          = verbose
      ) %>%
      attach_tf_rna_corr_to_pvals(
        tf         = tf_symbol,
        pval_tbl   = .,
        strict_rna = strict_rna,
        r_col_out  = "r_rna",
        verbose    = verbose
      )

    pdat_fixed <- tf_pvals_fixed_collapsed[[tf_symbol]] %>%
      attach_fp_r_corr_to_pvals(
        tf               = tf_symbol,
        pval_tbl         = .,
        tf_gene_links_gh = tf_gene_links_gh,
        r_col_source     = "r_fp",
        r_col_out        = "r_fp",
        verbose          = verbose
      ) %>%
      attach_tf_rna_corr_to_pvals(
        tf         = tf_symbol,
        pval_tbl   = .,
        strict_rna = strict_rna,
        r_col_out  = "r_rna",
        verbose    = verbose
      )

    pdat_fixed_int <- tf_pvals_fixed_int_collapsed[[tf_symbol]] %>%
      attach_fp_r_corr_to_pvals(
        tf               = tf_symbol,
        pval_tbl         = .,
        tf_gene_links_gh = tf_gene_links_gh,
        r_col_source     = "r_fp",
        r_col_out        = "r_fp",
        verbose          = verbose
      ) %>%
      attach_tf_rna_corr_to_pvals(
        tf         = tf_symbol,
        pval_tbl   = .,
        strict_rna = strict_rna,
        r_col_out  = "r_rna",
        verbose    = verbose
      )

    out_files <- list(
      lmm = plot_ko_summary_violin_for_tf(
        tf             = tf_symbol,
        ko_tbl         = ko_tbl,
        pval_tbl       = pdat_lmm,
        panel_type     = "lmm",
        p_cuts         = p_cuts,
        ko_type        = ko_type,
        out_dir        = out_dir,
        verbose        = verbose,
        filename_suffix = paste0(filename_suffix, "_LMM")
      ),
      fixed = plot_ko_summary_violin_for_tf(
        tf             = tf_symbol,
        ko_tbl         = ko_tbl,
        pval_tbl       = pdat_fixed,
        panel_type     = "fixed",
        p_cuts         = p_cuts,
        ko_type        = ko_type,
        out_dir        = out_dir,
        verbose        = verbose,
        filename_suffix = paste0(filename_suffix, "_FIXED")
      ),
      fixed_int = plot_ko_summary_violin_for_tf(
        tf             = tf_symbol,
        ko_tbl         = ko_tbl,
        pval_tbl       = pdat_fixed_int,
        panel_type     = "fixed_int",
        p_cuts         = p_cuts,
        ko_type        = ko_type,
        out_dir        = out_dir,
        verbose        = verbose,
        filename_suffix = paste0(filename_suffix, "_FIXED_INT")
      )
    )

    out_files
  })

  names(res) <- tfs_use
  invisible(res)
}

ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

## ==================================================================
## 9) TF perturbation DB helpers and TFBS loading (TFBS not used in new KO plots)
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
## 10) KO truth tables from DB – ensure HNF1A is included
## ==================================================================

lfc_strong <- 1
padj_sig   <- 0.05
lfc_unch   <- 0.25
padj_unch  <- 0.5

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
        TRUE ~ NA_character_
      )
    )
}

## Rebuild KO truth tables (now including HNF1A) ----------------------

ko_truth_IRF1  <- make_ko_truth_from_db("IRF1",  tf_perturb_db)
ko_truth_RARG  <- make_ko_truth_from_db("RARG",  tf_perturb_db)
ko_truth_SOX9  <- make_ko_truth_from_db("SOX9",  tf_perturb_db)
ko_truth_KLF5  <- make_ko_truth_from_db("KLF5",  tf_perturb_db)
ko_truth_FOXA2 <- make_ko_truth_from_db(
  tf_symbol      = "FOXA2",
  tf_perturb_db  = tf_perturb_db,
  subtype_filter = "PANC1"
)
ko_truth_HNF1A <- make_ko_truth_from_db("HNF1A", tf_perturb_db)

ko_truth_list <- list(
  FOXA2 = ko_truth_FOXA2,
  IRF1  = ko_truth_IRF1,
  KLF5  = ko_truth_KLF5,
  RARG  = ko_truth_RARG,
  HNF1A = ko_truth_HNF1A,
  SOX9  = ko_truth_SOX9
)


ko_truth_FOXA2 |>
  # dplyr::filter(ko_group == "Down") |> # Up Down Unchanged
  ggplot(aes(x = log2fc)) +
  geom_histogram(bins = 100, color = "black", fill = "steelblue") +
  theme_minimal(base_size = 12) +
  labs(x = "log2FC", y = "Count", title = "Distribution of log2FC")

## TFs to plot -------------------------

tfs_all <- c("FOXA2", "IRF1", "KLF5", "RARG", "HNF1A", "SOX9")

## ==================================================================
## 11) Main call example
## ==================================================================

# usage (assuming tf_pvals_*_collapsed and strict_rna are in the workspace):
# strict_rna should have columns: HGNC, ensembl_gene_id, and sample columns.
run_all_ko_violin(
  tfs                           = tfs_all,
  ko_truth_list                 = ko_truth_list,
  tf_pvals_lmm_collapsed        = tf_pvals_lmm_collapsed,
  tf_pvals_fixed_collapsed      = tf_pvals_fixed_collapsed,
  tf_pvals_fixed_int_collapsed  = tf_pvals_fixed_int_collapsed,
  strict_rna                    = strict_rna,
  ko_type                       = "KO",
  out_dir                       = ko_dir,
  p_cuts                        = c(5e-11, 5e-10, 5e-9, 5e-8, 5e-7, 5e-6, 5e-5),
  verbose                       = TRUE,
  filename_suffix               = ""
)


