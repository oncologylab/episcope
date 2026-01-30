library(episcope)
progressr::handlers(global = TRUE)

in_dir <- "/data/homes/yl814/episcope_test/nutrient_stress/connect_tfs_to_target_genes"
tf_perturb_db <- "/data/homes/yl814/episcope/tf_perturb.db"
predicted_tfbs_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

tf_gene_links_gh   <- readr::read_csv(file.path(in_dir, "fp_gene_corr_full_jaspar2024.csv"))
# atac_gene_links_gh <- readr::read_csv(file.path(in_dir, "atac_gene_corr_full_jaspar2024.csv"))
TFLink <- tflink_load("hs")
TFLink <- tflink_load("hs")

unique(tf_gene_links_gh$tfs)

motif_db <- readr::read_tsv(system.file("extdata", "genome/JASPAR2024.txt", package = "episcope"))

# HNF4A is not in the "tf_perturb.db", will load separately
# ---- make_ko_tbl() ---------------------------------------------------------
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
HNF4A_KO <- readr::read_csv("/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/Mayo 5289 siHNF4A RNA-seq.csv")

# TF perturbation DB helpers and TFBS loading
# ---- %||%() ---------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Parse TF from dataset_id like "..._<TF>_(KD|KO|OE|CRISPRa|CRISPRi)"
# ---- .extract_tf_from_dataset_id() -----------------------------------------
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
# ---- get_tf_perturbation_tbl() ---------------------------------------------
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

# ---- get_tf_tfbs() ---------------------------------------------------------
get_tf_tfbs <- function(tf,
                        tfbs_dir      = predicted_tfbs_dir,
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

# Build KO truth tables from DB
lfc_strong <- 1        # |log2FC| > 1 for regulated
padj_sig   <- 0.05     # significant
lfc_unch   <- 0.25     # |log2FC| < 0.25 for unchanged
padj_unch  <- 0.5      # padj > 0.5 for unchanged

# ---- make_ko_truth_from_db() -----------------------------------------------
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

# ---- make_ko_truth_from_deseq2() -------------------------------------------
make_ko_truth_from_deseq2 <- function(tbl,
                                      lfc_strong = 1,
                                      padj_sig   = 0.05,
                                      lfc_unch   = 0.25,
                                      padj_unch  = 0.5) {

  tbl |>
    dplyr::filter(!is.na(symbol), symbol != "") |>
    dplyr::transmute(
      gene     = symbol,
      log2fc   = log2FoldChange,
      ko_group = dplyr::case_when(
        !is.na(padj) & padj <= padj_sig & log2FoldChange <= -lfc_strong ~ "Down",
        !is.na(padj) & padj <= padj_sig & log2FoldChange >=  lfc_strong ~ "Up",
        !is.na(padj) & padj >= padj_unch & abs(log2FoldChange) < lfc_unch ~ "Unchanged",
        TRUE ~ NA_character_
      )
    )
}
## TFBS per TF
tf_tfbs_HNF1A <- get_tf_tfbs("HNF1A", tfbs_dir = predicted_tfbs_dir)
tf_tfbs_HNF4A <- get_tf_tfbs("HNF4A", tfbs_dir = predicted_tfbs_dir)
tf_tfbs_IRF1  <- get_tf_tfbs("IRF1",  tfbs_dir = predicted_tfbs_dir)
tf_tfbs_RARG  <- get_tf_tfbs("RARG",  tfbs_dir = predicted_tfbs_dir)
tf_tfbs_SOX9  <- get_tf_tfbs("SOX9",  tfbs_dir = predicted_tfbs_dir)
tf_tfbs_KLF5  <- get_tf_tfbs("KLF5",  tfbs_dir = predicted_tfbs_dir)
tf_tfbs_FOXA2 <- get_tf_tfbs("FOXA2", tfbs_dir = predicted_tfbs_dir)

## KO truth tables
ko_truth_HNF1A <- make_ko_truth_from_db("HNF1A", tf_perturb_db)
ko_truth_HNF4A <- make_ko_truth_from_deseq2(HNF4A_KO)
ko_truth_IRF1  <- make_ko_truth_from_db("IRF1",  tf_perturb_db)
ko_truth_RARG  <- make_ko_truth_from_db("RARG",  tf_perturb_db)
ko_truth_SOX9  <- make_ko_truth_from_db("SOX9",  tf_perturb_db)
ko_truth_KLF5  <- make_ko_truth_from_db("KLF5",  tf_perturb_db)
ko_truth_FOXA2 <- make_ko_truth_from_db("FOXA2", tf_perturb_db, subtype_filter = "PANC1") # FOXA2: only PANC1 sub type matching with our cell lines

## Lists used by run_all_ko_violin
ko_truth_list <- list(
  HNF1A = ko_truth_HNF1A,
  HNF4A = ko_truth_HNF4A,
  IRF1  = ko_truth_IRF1,
  RARG  = ko_truth_RARG,
  SOX9  = ko_truth_SOX9,
  KLF5  = ko_truth_KLF5,
  FOXA2 = ko_truth_FOXA2
)

tf_tfbs_list <- list(
  HNF1A = tf_tfbs_HNF1A,
  HNF4A = tf_tfbs_HNF4A,
  IRF1  = tf_tfbs_IRF1,
  RARG  = tf_tfbs_RARG,
  SOX9  = tf_tfbs_SOX9,
  KLF5  = tf_tfbs_KLF5,
  FOXA2 = tf_tfbs_FOXA2
)

## TFs to plot
tfs_all <- c("HNF1A", "SOX9")


ko_truth_HNF1A |> dplyr::filter(ko_group == "Down")

ko_truth_HNF1A |>
  dplyr::arrange(abs(log2fc)) |>
  dplyr::slice_head(n = 1000)

# ---- Candidate TF->target links for small GAM test -------------------------
# Strategy:
# - keep ALL "Down" genes from KO truth
# - keep top N "Unchanged" genes closest to 0 (smallest |log2fc|)
# - filter `tf_gene_links_gh` to TF motif hits whose `gene_key` is in that set
.ko_gene_subset <- function(ko_truth_tbl, n_unchanged = 1000) {
  if (!is.data.frame(ko_truth_tbl) ||
      !all(c("gene", "log2fc", "ko_group") %in% names(ko_truth_tbl))) {
    stop("ko_truth_tbl must have columns: gene, log2fc, ko_group")
  }

  genes_down <- ko_truth_tbl |>
    dplyr::filter(.data$ko_group == "Down") |>
    dplyr::pull(.data$gene)

  genes_unch <- ko_truth_tbl |>
    dplyr::filter(.data$ko_group == "Unchanged") |>
    dplyr::arrange(abs(.data$log2fc)) |>
    dplyr::slice_head(n = n_unchanged) |>
    dplyr::pull(.data$gene)

  genes_use <- unique(c(genes_down, genes_unch))
  genes_use[!is.na(genes_use) & genes_use != ""]
}

.filter_tf_gene_links_for_ko_subset <- function(tf,
                                               ko_truth_tbl,
                                               tf_gene_links_gh,
                                               motif_db,
                                               n_unchanged = 1000) {
  tf <- toupper(as.character(tf))

  need_links <- c("fp_peak", "gene_key")
  if (!is.data.frame(tf_gene_links_gh) || !all(need_links %in% names(tf_gene_links_gh))) {
    stop("tf_gene_links_gh must contain: ", paste(need_links, collapse = ", "))
  }

  genes_use <- .ko_gene_subset(ko_truth_tbl, n_unchanged = n_unchanged)
  has_tfs   <- "tfs" %in% names(tf_gene_links_gh)
  has_motifs <- "motifs" %in% names(tf_gene_links_gh)
  use_motifs <- !has_tfs && has_motifs

  if (use_motifs) {
    if (!is.data.frame(motif_db) || !all(c("motif", "HGNC") %in% names(motif_db))) {
      stop("motif_db must contain columns: motif, HGNC (required when tf_gene_links_gh has no 'tfs' column).")
    }
    motifs_tf <- unique(motif_db$motif[toupper(motif_db$HGNC) == tf])
    motifs_tf <- motifs_tf[!is.na(motifs_tf)]
  } else {
    motifs_tf <- character(0)
  }

  if (!length(genes_use) || (use_motifs && !length(motifs_tf))) {
    out <- tf_gene_links_gh[0, , drop = FALSE]
    out$tf <- character(0)
    return(out)
  }

  ko_annot <- ko_truth_tbl |>
    dplyr::select(.data$gene, .data$ko_group, .data$log2fc) |>
    dplyr::distinct(.data$gene, .keep_all = TRUE)

  links_sub <- if (use_motifs) {
    tf_gene_links_gh |>
      dplyr::filter(.data$motifs %in% motifs_tf)
  } else {
    tf_pat <- paste0("(^|[,;\\s])", tf, "($|[,;\\s])")
    tf_gene_links_gh |>
      dplyr::filter(!is.na(.data$tfs), grepl(tf_pat, toupper(.data$tfs)))
  }

  out <- links_sub |>
    dplyr::filter(.data$gene_key %in% genes_use) |>
    dplyr::left_join(ko_annot, by = c("gene_key" = "gene")) |>
    dplyr::mutate(tf = tf) |>
    dplyr::relocate(.data$tf, .before = 1)

  if (use_motifs) {
    out <- dplyr::distinct(out, .data$tf, .data$fp_peak, .data$gene_key, .data$motifs, .keep_all = TRUE)
  } else {
    out <- dplyr::distinct(out, .data$tf, .data$fp_peak, .data$gene_key, .data$tfs, .keep_all = TRUE)
  }

  out
}

candidate_links_HNF1A <- .filter_tf_gene_links_for_ko_subset(
  tf           = "HNF1A",
  ko_truth_tbl = ko_truth_HNF1A,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db     = motif_db,
  n_unchanged  = 1000
)

candidate_links_SOX9 <- .filter_tf_gene_links_for_ko_subset(
  tf           = "SOX9",
  ko_truth_tbl = ko_truth_SOX9,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db     = motif_db,
  n_unchanged  = 1000
)

candidate_links_small_gam <- dplyr::bind_rows(candidate_links_HNF1A, candidate_links_SOX9)

candidate_links_small_gam |>
  dplyr::count(.data$tf, .data$ko_group, name = "n_links") |>
  dplyr::arrange(.data$tf, dplyr::desc(.data$n_links)) |>
  print(n = Inf)

# ---- Small GAM grid test (HNF1A + SOX9) -------------------------------------
# Guarded: won't run on source() unless enabled.
# Enable explicitly via one of:
#   options(episcope.run_small_gam_grid = TRUE)
#   Sys.setenv(EPISCOPE_RUN_SMALL_GAM_GRID = "1")
if (isTRUE(getOption("episcope.run_small_gam_grid")) ||
    identical(Sys.getenv("EPISCOPE_RUN_SMALL_GAM_GRID"), "1")) {

  # after `grn_set` is loaded/created
  source("dev/07_12_benchmark_tf_to_target_genes_method_gam.R")

  # `max_pairs` limits how many (fp_peak, gene_key) pairs are fit (after ordering).
  # Set `max_pairs = Inf` to run all pairs, but note runtime grows quickly:
  #   n_pairs * length(methods) * length(bs_values) * length(k_values)
  methods_use <- c("REML", "ML", "GCV.Cp")
  bs_values_use <- c("tp", "cr", "ps")
  k_values_use <- 3:10
  max_pairs_use <- 100L
  n_workers_use <- 36L

  # ---- Model performance: Down vs Unchanged --------------------------------
  # Goal: p(smooth) < alpha predicts "Down", else predicts "Unchanged".
  alpha_model <- 0.05

  .eval_binary_tbl <- function(tbl, alpha = 0.05) {
    if (!all(c("ko_group", "method", "bs", "smooth_p") %in% names(tbl))) {
      stop("Need columns: ko_group, method, bs, smooth_p")
    }
    p_col <- if ("smooth_p_cap" %in% names(tbl)) "smooth_p_cap" else "smooth_p"

    dplyr::as_tibble(tbl) |>
      dplyr::filter(.data$ok, .data$ko_group %in% c("Down", "Unchanged")) |>
      dplyr::mutate(
        p_model = suppressWarnings(as.numeric(.data[[p_col]])),
        truth_down = (.data$ko_group == "Down"),
        pred_down  = is.finite(.data$p_model) & (.data$p_model < alpha),
        outcome = dplyr::case_when(
          truth_down & pred_down  ~ "TP",
          truth_down & !pred_down ~ "FN",
          !truth_down & pred_down ~ "FP",
          !truth_down & !pred_down ~ "TN",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(.data$outcome))
  }

  .summarise_binary_perf <- function(eval_tbl, group_cols = c("method", "bs")) {
    eval_tbl |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        n = dplyr::n(),
        TP = sum(.data$outcome == "TP"),
        TN = sum(.data$outcome == "TN"),
        FP = sum(.data$outcome == "FP"),
        FN = sum(.data$outcome == "FN"),
        accuracy = (TP + TN) / n,
        tpr = dplyr::if_else((TP + FN) > 0, TP / (TP + FN), NA_real_),  # sensitivity / recall
        tnr = dplyr::if_else((TN + FP) > 0, TN / (TN + FP), NA_real_),  # specificity
        precision = dplyr::if_else((TP + FP) > 0, TP / (TP + FP), NA_real_),
        f1 = dplyr::if_else(
          is.finite(precision) & is.finite(tpr) & (precision + tpr) > 0,
          2 * precision * tpr / (precision + tpr),
          NA_real_
        ),
        mcc = dplyr::if_else(
          (TP + FP) > 0 & (TP + FN) > 0 & (TN + FP) > 0 & (TN + FN) > 0,
          (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
          NA_real_
        ),
        balanced_accuracy = (tpr + tnr) / 2,
        .groups = "drop"
      )
  }

  .pick_one_peak_per_gene <- function(eval_tbl, alpha = 0.05,
                                      group_cols = c("method", "bs", "k_used", "k_selected_reason")) {
    if (!all(c("gene_key", "fp_peak", "ko_group", "p_model") %in% names(eval_tbl))) {
      stop("Need columns: gene_key, fp_peak, ko_group, p_model")
    }
    eval_tbl |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "gene_key")))) |>
      dplyr::mutate(
        margin = dplyr::case_when(
          .data$ko_group == "Down"      ~ alpha - .data$p_model,
          .data$ko_group == "Unchanged" ~ .data$p_model - alpha,
          TRUE ~ NA_real_
        )
      ) |>
      dplyr::arrange(dplyr::desc(.data$margin), .data$p_model) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
  }

  .run_one_tf <- function(tf_label, candidate_links_tbl) {
    pair_tbl <- candidate_links_tbl |>
      dplyr::distinct(.data$fp_peak, .data$gene_key, .keep_all = TRUE)

    gam_grid <- run_fp_gam_grid_for_pairs(
      fp_score = grn_set$fp_score,
      rna_tbl  = grn_set$rna,
      sample_md = grn_set$sample_metadata_used,
      pair_tbl = pair_tbl,
      max_pairs = max_pairs_use,
      order_by = "p_fp",
      methods = methods_use,
      k_values = k_values_use,
      bs_values = bs_values_use,
      run_gam_check = TRUE,
      keep_gamcheck_output = FALSE,
      n_workers = n_workers_use
    )

    res_anno <- annotate_fp_gam_grid_results(gam_grid$results)

    eval_all_pairs <- .eval_binary_tbl(res_anno, alpha = alpha_model)
    eval_all_genes <- .pick_one_peak_per_gene(
      eval_all_pairs,
      alpha = alpha_model,
      group_cols = c("method", "bs", "k_used")
    )
    perf_allk <- .summarise_binary_perf(
      eval_all_genes,
      group_cols = c("method", "bs", "k_used")
    )

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      message("ggplot2 not installed; skipping plot for ", tf_label)
      return(invisible(list(results = res_anno, perf_allk = perf_allk)))
    }

    # ---- "Best k per link" (fp_peak|gene_key) based on k-check --------------
    # For each (pair_key, method, bs): pick the smallest k_used passing k-check;
    # if none pass, pick the largest k_used.
    k_opt_by_link <- pick_optimal_k_by_gamcheck(
      res_anno,
      score = "AIC",
      group_cols = c("pair_key", "method", "bs"),
      k_col = "k_used",
      k_index_min = 0.9,
      k_p_min = 0.05
    )

    k_opt_dist <- k_opt_by_link |>
      dplyr::count(.data$method, .data$bs, .data$k_used, .data$k_selected_reason, name = "n") |>
      dplyr::group_by(.data$method, .data$bs) |>
      dplyr::mutate(frac = .data$n / sum(.data$n)) |>
      dplyr::ungroup() |>
      dplyr::mutate(method_bs = paste(.data$method, .data$bs, sep = " | "))

    p_k <- ggplot2::ggplot(
      k_opt_dist,
      ggplot2::aes(
        x = factor(.data$k_used),
        y = .data$method_bs,
        fill = .data$frac
      )
    ) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = .data$n), size = 3, fontface = "bold") +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
      ggplot2::labs(
        x = "Selected k_used",
        y = NULL,
        fill = "Fraction",
        title = sprintf("%s optimal k by link (fp_peak|gene): k-check driven", tf_label),
        subtitle = "Pick smallest k passing k-check; else largest k"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        legend.position = "top",
        text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(face = "bold")
      )

    if (exists("ko_dir") && is.character(ko_dir) && nzchar(ko_dir)) {
      out_pdf_k <- file.path(
        ko_dir,
        sprintf(
          "%s_small_GAM_grid_optimal_k_by_link_k_%d_%d_alpha_%s.pdf",
          tf_label,
          min(k_values_use),
          max(k_values_use),
          format(alpha_model, scientific = TRUE)
        )
      )
      ggplot2::ggsave(filename = out_pdf_k, plot = p_k, width = 12, height = 6, units = "in")
      message("Saved plot: ", out_pdf_k)
    } else {
      message("ko_dir not set; skipping PDF save.")
    }

    # Mark best k per (method, bs) by balanced accuracy (tie-break: accuracy).
    best_k_tbl <- perf_allk |>
      dplyr::group_by(.data$method, .data$bs) |>
      dplyr::arrange(dplyr::desc(.data$balanced_accuracy), dplyr::desc(.data$accuracy)) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::select(.data$method, .data$bs, best_k = .data$k_used)

    perf_plot_tbl <- perf_allk |>
      dplyr::left_join(best_k_tbl, by = c("method", "bs")) |>
      dplyr::mutate(
        is_best_k = !is.na(.data$best_k) & .data$k_used == .data$best_k,
        model_label = paste(.data$method, .data$bs, paste0("k=", .data$k_used), sep = " | ")
      )

    if (requireNamespace("tidyr", quietly = TRUE)) {
      perf_plot_long <- perf_plot_tbl |>
        tidyr::pivot_longer(
          cols = c("precision", "accuracy", "f1"),
          names_to = "metric",
          values_to = "value"
        )
    } else {
      perf_plot_long <- perf_plot_tbl
      perf_plot_long$metric <- "accuracy"
      perf_plot_long$value <- perf_plot_long$accuracy
    }

    perf_plot_long <- perf_plot_long |>
      dplyr::mutate(
        model_label_metric = paste(.data$metric, as.character(.data$model_label), sep = "___"),
        star = ifelse(.data$is_best_k, "*", "")
      )

    perf_plot_long$metric <- factor(perf_plot_long$metric, levels = c("precision", "accuracy", "f1"))

    p1 <- ggplot2::ggplot(
      perf_plot_long,
      ggplot2::aes(
        x = stats::reorder(model_label_metric, value),
        y = .data$value,
        fill = .data$method
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::geom_text(
        ggplot2::aes(y = -0.02, label = .data$star),
        inherit.aes = TRUE,
        color = "red",
        fontface = "bold",
        size = 3
      ) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
      ggplot2::labs(
        x = NULL,
        y = "Score",
        title = sprintf("%s GAM grid (per gene): p(smooth)<%.2g predicts Down", tf_label, alpha_model),
        subtitle = "* marks best k per (method, bs) by balanced accuracy\nPrecision=TP/(TP+FP); Accuracy=(TP+TN)/N; F1=2PR/(P+R)"
      ) +
      ggplot2::facet_wrap(~metric, ncol = 4, scales = "free_y") +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        legend.position = "top",
        plot.margin = ggplot2::margin(10, 30, 10, 10),
        text = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.text = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(face = "bold")
      )

    if (exists("ko_dir") && is.character(ko_dir) && nzchar(ko_dir)) {
      out_pdf <- file.path(
        ko_dir,
        sprintf(
          "%s_small_GAM_grid_per_gene_metrics_k_%d_%d_alpha_%s.pdf",
          tf_label,
          min(k_values_use),
          max(k_values_use),
          format(alpha_model, scientific = TRUE)
        )
      )
      ggplot2::ggsave(filename = out_pdf, plot = p1, width = 24, height = 18, units = "in")
      message("Saved plot: ", out_pdf)
    } else {
      message("ko_dir not set; skipping PDF save.")
    }

      print(p1)

    invisible(list(results = res_anno, perf_allk = perf_allk, k_opt_by_link = k_opt_by_link))
  }

  gam_grid_runs <- list(
    HNF1A = .run_one_tf("HNF1A", candidate_links_HNF1A),
    SOX9  = .run_one_tf("SOX9",  candidate_links_SOX9)
  )
}




# Fast refactor (from line 270): RNA/FP models

# ---- p-value helpers --------------------------------------------------------
.clip_logp <- function(logp, floor = 1e-300) {
  ifelse(is.na(logp), NA_real_, pmax(logp, log(floor)))
}

.p_from_logp <- function(logp, floor = 1e-300) {
  logp2 <- .clip_logp(logp, floor = floor)
  list(
    p = exp(logp2),
    mlog10p = -logp2 / log(10)
  )
}

.p_from_t <- function(tval, df, floor = 1e-300) {
  logp <- log(2) + stats::pt(-abs(tval), df = df, log.p = TRUE)
  .p_from_logp(logp, floor = floor)
}

.p_from_z <- function(z, floor = 1e-300) {
  logp <- log(2) + stats::pnorm(-abs(z), log.p = TRUE)
  .p_from_logp(logp, floor = floor)
}

.p_from_f <- function(fstat, df1, df2, floor = 1e-300) {
  logp <- stats::pf(fstat, df1 = df1, df2 = df2, lower.tail = FALSE, log.p = TRUE)
  .p_from_logp(logp, floor = floor)
}

# p-value floor + -log10(p) (handles mgcv returning 0)
.cap_p <- function(p, floor = 1e-300) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  p[p < 0] <- NA_real_
  p[p == 0] <- floor
  p[p < floor] <- floor
  p[p > 1] <- 1
  list(p = p, mlog10p = -log10(p))
}

.safe_cor <- function(x, y, floor = 1e-300) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3L) return(list(r = NA_real_, p = NA_real_, mlog10p = NA_real_, n = n))
  x <- x[ok]; y <- y[ok]
  if (stats::sd(x) == 0 || stats::sd(y) == 0) return(list(r = NA_real_, p = NA_real_, mlog10p = NA_real_, n = n))
  r <- stats::cor(x, y)
  r2 <- pmin(pmax(r, -0.999999999999), 0.999999999999)
  tval <- r2 * sqrt((n - 2) / (1 - r2^2))
  pv <- .p_from_t(tval, df = n - 2, floor = floor)
  list(r = as.numeric(r), p = pv$p, mlog10p = pv$mlog10p, n = n)
}

.mgcv_smooth_p <- function(fit, floor = 1e-300) {
  if (is.null(fit)) return(list(p = NA_real_, mlog10p = NA_real_, edf = NA_real_))
  st <- tryCatch(summary(fit)$s.table, error = function(e) NULL)
  if (is.null(st) || !nrow(st)) return(list(p = NA_real_, mlog10p = NA_real_, edf = NA_real_))
  p_raw <- st[1, ncol(st)]
  edf <- st[1, "edf"]
  pv <- .cap_p(p_raw, floor = floor)
  list(p = pv$p, mlog10p = pv$mlog10p, edf = as.numeric(edf))
}

# ---- data helpers -----------------------------------------------------------
.tf_peak_gene_map <- function(tf, tf_gene_links_gh, motif_db) {
  need_links <- c("fp_peak", "gene_key", "motifs")
  if (!all(need_links %in% names(tf_gene_links_gh))) {
    stop("tf_gene_links_gh must contain: ", paste(need_links, collapse = ", "))
  }
  if (!all(c("motif", "HGNC") %in% names(motif_db))) {
    stop("motif_db must contain columns: motif, HGNC")
  }

  motifs_tf <- unique(motif_db$motif[motif_db$HGNC == tf])
  motifs_tf <- motifs_tf[!is.na(motifs_tf)]
  if (!length(motifs_tf)) {
    return(data.frame(fp_peak = character(0), gene = character(0), stringsAsFactors = FALSE))
  }

  map <- tf_gene_links_gh[tf_gene_links_gh$motifs %in% motifs_tf, c("fp_peak", "gene_key"), drop = FALSE]
  map <- map[!is.na(map$fp_peak) & !is.na(map$gene_key) & map$gene_key != "", , drop = FALSE]
  names(map) <- c("fp_peak", "gene")
  unique(map)
}

.extract_tf_expr <- function(rna_tbl, tf, sample_ids) .gene_expr_vec(rna_tbl, tf, sample_ids)

.gene_expr_vec <- function(rna_tbl, gene, sample_ids) {
  idx <- which(rna_tbl$HGNC == gene)
  if (!length(idx)) return(setNames(rep(NA_real_, length(sample_ids)), sample_ids))
  mat <- as.matrix(rna_tbl[idx, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA_real_
  out <- colMeans(mat, na.rm = TRUE)
  setNames(as.numeric(out), sample_ids)
}

.chunk_vec <- function(x, n_chunks) {
  x <- as.list(x)
  n_chunks <- as.integer(n_chunks)
  if (!length(x)) return(vector("list", n_chunks))
  if (!is.finite(n_chunks) || is.na(n_chunks) || n_chunks < 1L) n_chunks <- 1L
  n_chunks <- min(n_chunks, length(x))
  split(unlist(x, use.names = FALSE), rep(seq_len(n_chunks), length.out = length(x)))
}

.chunk_by_weight <- function(ids, weights, n_chunks) {
  ids <- as.character(ids)
  weights <- as.numeric(weights)
  n_chunks <- as.integer(n_chunks)
  if (!length(ids)) return(vector("list", n_chunks))
  if (!is.finite(n_chunks) || is.na(n_chunks) || n_chunks < 1L) n_chunks <- 1L
  n_chunks <- min(n_chunks, length(ids))

  weights[!is.finite(weights) | is.na(weights) | weights < 0] <- 1
  ord <- order(weights, decreasing = TRUE, na.last = TRUE)
  ids <- ids[ord]
  weights <- weights[ord]

  chunks <- vector("list", n_chunks)
  load <- rep(0, n_chunks)
  for (i in seq_along(ids)) {
    j <- which.min(load)
    chunks[[j]] <- c(chunks[[j]], ids[i])
    load[j] <- load[j] + weights[i]
  }
  chunks
}

.build_expr_mat <- function(rna_tbl, genes, sample_ids) {
  if (!"HGNC" %in% names(rna_tbl)) stop("grn_set$rna must contain column 'HGNC'.")
  sub <- rna_tbl[rna_tbl$HGNC %in% genes, c("HGNC", sample_ids), drop = FALSE]
  if (!nrow(sub)) return(matrix(numeric(0), 0, length(sample_ids), dimnames = list(character(0), sample_ids)))

  mat <- as.matrix(sub[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- NA_real_
  grp <- as.character(sub$HGNC)

  sum_mat <- rowsum(mat, group = grp, reorder = FALSE, na.rm = TRUE)
  cnt_mat <- rowsum((!is.na(mat)) * 1, group = grp, reorder = FALSE)
  sum_mat / pmax(cnt_mat, 1)
}

.build_fp_mat_mean <- function(fp_score, map, sample_ids) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need package 'Matrix' for fast footprint aggregation.")
  if (!"peak_ID" %in% names(fp_score)) stop("grn_set$fp_score must contain column 'peak_ID'.")

  peaks <- unique(map$fp_peak)
  fp_sub <- fp_score[fp_score$peak_ID %in% peaks, c("peak_ID", sample_ids), drop = FALSE]
  if (!nrow(fp_sub)) return(matrix(numeric(0), 0, length(sample_ids), dimnames = list(character(0), sample_ids)))

  peak_ids <- as.character(fp_sub$peak_ID)
  fp_mat <- as.matrix(fp_sub[, sample_ids, drop = FALSE])
  storage.mode(fp_mat) <- "numeric"
  rownames(fp_mat) <- peak_ids

  map2 <- unique(map[map$fp_peak %in% peak_ids, , drop = FALSE])
  genes <- sort(unique(map2$gene))
  i <- match(map2$gene, genes)
  j <- match(map2$fp_peak, peak_ids)

  A <- Matrix::sparseMatrix(
    i = i, j = j, x = 1,
    dims = c(length(genes), length(peak_ids)),
    dimnames = list(genes, peak_ids)
  )

  ok <- is.finite(fp_mat)
  fp0 <- fp_mat
  fp0[!ok] <- 0

  sum_mat <- A %*% fp0
  n_mat <- A %*% (ok + 0)
  out <- sum_mat / pmax(n_mat, 1)

  as.matrix(out)
}

# ---- model runners ----------------------------------------------------------
.corr_gene_vs_vec <- function(expr_mat, x, floor = 1e-300) {
  ok_x <- is.finite(x)
  if (sum(ok_x) < 3L || stats::sd(x[ok_x], na.rm = TRUE) == 0) {
    return(list(r = rep(NA_real_, nrow(expr_mat)),
                p = rep(NA_real_, nrow(expr_mat)),
                mlog10p = rep(NA_real_, nrow(expr_mat)),
                n = rep(0L, nrow(expr_mat))))
  }

  em <- expr_mat[, ok_x, drop = FALSE]
  xv <- x[ok_x]
  r <- as.numeric(stats::cor(t(em), xv, use = "pairwise.complete.obs"))
  n <- rowSums(is.finite(em))

  r2 <- pmin(pmax(r, -0.999999999999), 0.999999999999)
  tval <- r2 * sqrt((n - 2) / (1 - r2^2))
  tval[n < 3] <- NA_real_
  pv <- .p_from_t(tval, df = n - 2, floor = floor)

  list(r = r, p = pv$p, mlog10p = pv$mlog10p, n = n)
}

.corr_gene_vs_gene <- function(expr_mat, fp_mat, floor = 1e-300) {
  x <- expr_mat
  y <- fp_mat
  ok <- is.finite(x) & is.finite(y)
  n <- rowSums(ok)

  x0 <- x; x0[!ok] <- 0
  y0 <- y; y0[!ok] <- 0

  sumx <- rowSums(x0)
  sumy <- rowSums(y0)
  sumx2 <- rowSums(x0^2)
  sumy2 <- rowSums(y0^2)
  sumxy <- rowSums(x0 * y0)

  num <- sumxy - (sumx * sumy) / n
  den <- sqrt((sumx2 - sumx^2 / n) * (sumy2 - sumy^2 / n))
  r <- num / den

  bad <- (n < 3L) | !is.finite(r) | den == 0
  r[bad] <- NA_real_

  r2 <- pmin(pmax(r, -0.999999999999), 0.999999999999)
  tval <- r2 * sqrt((n - 2) / (1 - r2^2))
  tval[n < 3] <- NA_real_
  pv <- .p_from_t(tval, df = n - 2, floor = floor)

  list(r = r, p = pv$p, mlog10p = pv$mlog10p, n = n)
}

.residualize <- function(Z, Y) {
  stats::lm.fit(Z, Y)$residuals
}

.lm_effect_x_vec <- function(expr_mat, x, Z, floor = 1e-300) {
  y <- t(expr_mat)
  n <- nrow(y)
  df_cov <- qr(Z)$rank
  df_res <- n - df_cov - 1L
  if (df_res <= 0) stop("Not enough samples for LM with covariates.")

  y_res <- .residualize(Z, y)
  x_res <- as.numeric(.residualize(Z, x))
  sxx <- sum(x_res^2)
  if (!is.finite(sxx) || sxx == 0) {
    return(list(beta = rep(NA_real_, ncol(y)), p = rep(NA_real_, ncol(y)), mlog10p = rep(NA_real_, ncol(y))))
  }

  beta <- colSums(y_res * x_res) / sxx
  resid <- y_res - tcrossprod(x_res, beta)
  sigma2 <- colSums(resid^2) / df_res
  se <- sqrt(sigma2 / sxx)
  tval <- beta / se
  pv <- .p_from_t(tval, df = df_res, floor = floor)

  list(beta = beta, p = pv$p, mlog10p = pv$mlog10p)
}

.lm_effect_x_mat <- function(expr_mat, x_mat, Z, floor = 1e-300) {
  y <- t(expr_mat)
  x <- t(x_mat)
  n <- nrow(y)
  df_cov <- qr(Z)$rank
  df_res <- n - df_cov - 1L
  if (df_res <= 0) stop("Not enough samples for LM with covariates.")

  y_res <- .residualize(Z, y)
  x_res <- .residualize(Z, x)

  sxx <- colSums(x_res^2)
  beta <- colSums(y_res * x_res) / sxx
  beta[!is.finite(beta) | sxx == 0] <- NA_real_

  resid <- y_res - sweep(x_res, 2, beta, `*`)
  sigma2 <- colSums(resid^2) / df_res
  se <- sqrt(sigma2 / sxx)
  tval <- beta / se
  pv <- .p_from_t(tval, df = df_res, floor = floor)

  list(beta = beta, p = pv$p, mlog10p = pv$mlog10p)
}

.spline_gam_x_vec <- function(expr_mat, x, Z, df_spline = 4L, floor = 1e-300) {
  y <- t(expr_mat)
  n <- nrow(y)
  df_cov <- qr(Z)$rank
  df_res0 <- n - df_cov

  B <- splines::ns(x, df = df_spline)
  X <- cbind(Z, B)
  fit1 <- stats::lm.fit(X, y)
  df_res1 <- fit1$df.residual
  if (df_res1 <= 0) {
    return(list(p = rep(NA_real_, ncol(y)), mlog10p = rep(NA_real_, ncol(y))))
  }

  # reduced SSE: covariates only
  sse0 <- colSums(.residualize(Z, y)^2)
  sse1 <- colSums(fit1$residuals^2)
  df1 <- df_res0 - df_res1

  fstat <- ((sse0 - sse1) / df1) / (sse1 / df_res1)
  fstat[!is.finite(fstat) | df1 <= 0 | sse0 < sse1] <- NA_real_
  pv <- .p_from_f(fstat, df1 = df1, df2 = df_res1, floor = floor)

  list(p = pv$p, mlog10p = pv$mlog10p)
}

.spline_gam_x_mat <- function(expr_mat, x_mat, Z, df_spline = 4L, floor = 1e-300) {
  y <- t(expr_mat)
  x <- t(x_mat)
  n <- nrow(y)
  df_cov <- qr(Z)$rank
  df_res0 <- n - df_cov

  y_res0 <- .residualize(Z, y)
  sse0 <- colSums(y_res0^2)

  p_out <- rep(NA_real_, ncol(y))
  mlog10p_out <- rep(NA_real_, ncol(y))

  for (i in seq_len(ncol(y))) {
    xi <- x[, i]
    yi <- y[, i]
    ok <- is.finite(xi) & is.finite(yi)
    if (sum(ok) < (df_cov + df_spline + 2L)) next

    bi <- tryCatch(splines::ns(xi[ok], df = df_spline), error = function(e) NULL)
    if (is.null(bi)) next

    Zi <- Z[ok, , drop = FALSE]
    Xi <- cbind(Zi, bi)
    fit1 <- tryCatch(stats::lm.fit(Xi, yi[ok]), error = function(e) NULL)
    if (is.null(fit1) || fit1$df.residual <= 0) next

    df1 <- (nrow(Zi) - fit1$df.residual) - df_cov
    if (!is.finite(df1) || df1 <= 0) next

    sse1 <- sum(fit1$residuals^2)
    sse0i <- sse0[i]
    if (!is.finite(sse0i) || !is.finite(sse1) || sse0i < sse1) next

    fstat <- ((sse0i - sse1) / df1) / (sse1 / fit1$df.residual)
    pv <- .p_from_f(fstat, df1 = df1, df2 = fit1$df.residual, floor = floor)

    p_out[i] <- pv$p
    mlog10p_out[i] <- pv$mlog10p
  }

  list(p = p_out, mlog10p = mlog10p_out)
}

.sem_linear_mediation <- function(expr_mat, x_rna, x_fp_mat, Z, floor = 1e-300) {
  y <- t(expr_mat)
  x1 <- as.numeric(x_rna)
  x2 <- t(x_fp_mat)

  n <- nrow(y)
  df_cov <- qr(Z)$rank
  df_a <- n - df_cov - 1L
  df_bc <- n - df_cov - 2L
  if (df_a <= 0 || df_bc <= 0) {
    return(list(a = rep(NA_real_, ncol(y)), p_a = rep(NA_real_, ncol(y)), mlog10p_a = rep(NA_real_, ncol(y)),
                b = rep(NA_real_, ncol(y)), p_b = rep(NA_real_, ncol(y)), mlog10p_b = rep(NA_real_, ncol(y)),
                c = rep(NA_real_, ncol(y)), p_c = rep(NA_real_, ncol(y)), mlog10p_c = rep(NA_real_, ncol(y)),
                ab = rep(NA_real_, ncol(y)), p_ab = rep(NA_real_, ncol(y)), mlog10p_ab = rep(NA_real_, ncol(y)),
                total = rep(NA_real_, ncol(y)), p_total = rep(NA_real_, ncol(y)), mlog10p_total = rep(NA_real_, ncol(y))))
  }

  y_res <- .residualize(Z, y)
  x1_res <- as.numeric(.residualize(Z, x1))
  x2_res <- .residualize(Z, x2)

  S11 <- sum(x1_res^2)
  S22 <- colSums(x2_res^2)
  S12 <- colSums(x1_res * x2_res)

  # a: x2 ~ x1
  a <- colSums(x2_res * x1_res) / S11
  resid_a <- x2_res - tcrossprod(x1_res, a)
  sigma2_a <- colSums(resid_a^2) / df_a
  se_a <- sqrt(sigma2_a / S11)
  t_a <- a / se_a
  pva <- .p_from_t(t_a, df = df_a, floor = floor)

  # b,c: y ~ x1 + x2
  Sy1 <- colSums(y_res * x1_res)
  Sy2 <- colSums(y_res * x2_res)
  det <- S11 * S22 - S12^2

  c <- (Sy1 * S22 - Sy2 * S12) / det
  b <- (Sy2 * S11 - Sy1 * S12) / det
  c[!is.finite(c) | det == 0] <- NA_real_
  b[!is.finite(b) | det == 0] <- NA_real_

  y_hat <- tcrossprod(x1_res, c) + sweep(x2_res, 2, b, `*`)
  resid_bc <- y_res - y_hat
  sigma2_bc <- colSums(resid_bc^2) / df_bc

  inv11 <- S22 / det
  inv22 <- S11 / det
  se_c <- sqrt(sigma2_bc * inv11)
  se_b <- sqrt(sigma2_bc * inv22)

  t_c <- c / se_c
  t_b <- b / se_b
  pvc <- .p_from_t(t_c, df = df_bc, floor = floor)
  pvb <- .p_from_t(t_b, df = df_bc, floor = floor)

  ab <- a * b
  se_ab <- sqrt((b^2) * (se_a^2) + (a^2) * (se_b^2))
  z_ab <- ab / se_ab
  pvab <- .p_from_z(z_ab, floor = floor)

  total <- c + ab
  se_total <- sqrt(se_c^2 + se_ab^2)
  z_total <- total / se_total
  pvt <- .p_from_z(z_total, floor = floor)

  list(
    a = a, p_a = pva$p, mlog10p_a = pva$mlog10p,
    b = b, p_b = pvb$p, mlog10p_b = pvb$mlog10p,
    c = c, p_c = pvc$p, mlog10p_c = pvc$mlog10p,
    ab = ab, p_ab = pvab$p, mlog10p_ab = pvab$mlog10p,
    total = total, p_total = pvt$p, mlog10p_total = pvt$mlog10p
  )
}

# ---- main runner (mgcv + lavaan; FP per peak) -------------------------------
.build_tf_long_peak_tbl <- function(tf,
                                    grn_set,
                                    tf_gene_links_gh,
                                    motif_db,
                                    min_samples = 6L) {
  fp_score  <- grn_set$fp_score
  rna_tbl   <- grn_set$rna
  sample_md <- grn_set$sample_metadata_used

  stress_col <- if ("stress_type" %in% names(sample_md)) "stress_type" else if ("stress" %in% names(sample_md)) "stress" else NA_character_
  if (!all(c("id", "cell") %in% names(sample_md)) || is.na(stress_col)) {
    stop("sample_metadata_used must contain: id, cell, stress_type (or stress).")
  }
  md <- sample_md[, c("id", "cell", stress_col), drop = FALSE]
  names(md)[names(md) == stress_col] <- "stress"
  md$cell <- factor(md$cell)
  md$stress <- factor(md$stress)
  sample_ids <- as.character(md$id)

  if (!all(sample_ids %in% names(fp_score))) stop("fp_score missing some sample columns.")
  if (!all(sample_ids %in% names(rna_tbl))) stop("rna_tbl missing some sample columns.")

  map <- .tf_peak_gene_map(tf, tf_gene_links_gh, motif_db)
  if (!nrow(map)) return(tibble::tibble())

  peaks_use <- unique(map$fp_peak)
  fp_sub <- fp_score[fp_score$peak_ID %in% peaks_use, c("peak_ID", sample_ids), drop = FALSE]
  if (!nrow(fp_sub)) return(tibble::tibble())

  fp_long <- tidyr::pivot_longer(fp_sub,
    cols = dplyr::all_of(sample_ids),
    names_to = "sample_id",
    values_to = "tf_fp"
  )
  fp_long <- dplyr::left_join(fp_long, map, by = c("peak_ID" = "fp_peak"))
  names(fp_long)[names(fp_long) == "peak_ID"] <- "fp_peak"

  genes_use <- unique(fp_long$gene)
  rna_sub <- rna_tbl[rna_tbl$HGNC %in% genes_use, c("HGNC", sample_ids), drop = FALSE]
  if (!nrow(rna_sub)) return(tibble::tibble())
  rna_long <- tidyr::pivot_longer(rna_sub,
    cols = dplyr::all_of(sample_ids),
    names_to = "sample_id",
    values_to = "expr"
  )
  names(rna_long)[names(rna_long) == "HGNC"] <- "gene"

  out <- dplyr::inner_join(fp_long, rna_long, by = c("gene", "sample_id"))
  tf_expr <- .extract_tf_expr(rna_tbl, tf, sample_ids)
  out$tf_expr <- tf_expr[out$sample_id]
  out <- dplyr::left_join(out, md, by = c("sample_id" = "id"))
  out <- out |>
    dplyr::filter(is.finite(expr), is.finite(tf_fp), is.finite(tf_expr)) |>
    dplyr::group_by(gene, fp_peak) |>
    dplyr::filter(dplyr::n() >= min_samples) |>
    dplyr::ungroup()

  tibble::as_tibble(out)
}

.fit_rna_models_one_gene <- function(tf, gene, rna_tbl, sample_ids, md,
                                     min_samples = 6L,
                                     p_floor = 1e-300, k = 5L, bs = "tp") {
  t0 <- proc.time()[[3]]
  y <- .gene_expr_vec(rna_tbl, gene, sample_ids)
  tf_expr <- .extract_tf_expr(rna_tbl, tf, sample_ids)
  df <- data.frame(expr = y, tf_expr = tf_expr, md, stringsAsFactors = FALSE)
  df <- df[is.finite(df$expr) & is.finite(df$tf_expr), , drop = FALSE]
  if (nrow(df) < min_samples) {
    return(tibble::tibble(
      tf = tf, gene = gene, n_samples = nrow(df),
      r_rna = NA_real_, p_rna = NA_real_, mlog10p_rna = NA_real_,
      beta_rna_lm = NA_real_, p_rna_lm = NA_real_, mlog10p_rna_lm = NA_real_,
      p_rna_gam = NA_real_, mlog10p_rna_gam = NA_real_, edf_rna_gam = NA_real_,
      time_rna_cor_s = NA_real_, time_rna_lm_s = NA_real_, time_rna_gam_s = NA_real_,
      time_rna_total_s = as.numeric(proc.time()[[3]] - t0)
    ))
  }

  t_cor0 <- proc.time()[[3]]
  corv <- .safe_cor(df$tf_expr, df$expr, floor = p_floor)
  t_cor <- as.numeric(proc.time()[[3]] - t_cor0)

  t_lm0 <- proc.time()[[3]]
  fit_lm <- tryCatch(stats::lm(expr ~ tf_expr + cell + stress, data = df), error = function(e) NULL)
  beta <- NA_real_; p_lm <- NA_real_; mlog10p_lm <- NA_real_
  if (!is.null(fit_lm)) {
    co <- summary(fit_lm)$coefficients
    if ("tf_expr" %in% rownames(co)) {
      beta <- co["tf_expr", "Estimate"]
      pv <- .cap_p(co["tf_expr", "Pr(>|t|)"], floor = p_floor)
      p_lm <- pv$p; mlog10p_lm <- pv$mlog10p
    }
  }
  t_lm <- as.numeric(proc.time()[[3]] - t_lm0)

  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Need package 'mgcv'.")
  s <- mgcv::s
  t_gam0 <- proc.time()[[3]]
  fit_gam <- tryCatch(
    mgcv::gam(expr ~ s(tf_expr, k = k, bs = bs) + cell + stress, data = df, method = "REML"),
    error = function(e) NULL
  )
  pg <- .mgcv_smooth_p(fit_gam, floor = p_floor)
  t_gam <- as.numeric(proc.time()[[3]] - t_gam0)
  t_tot <- as.numeric(proc.time()[[3]] - t0)

  tibble::tibble(
    tf = tf, gene = gene, n_samples = nrow(df),
    r_rna = corv$r, p_rna = corv$p, mlog10p_rna = corv$mlog10p,
    beta_rna_lm = beta, p_rna_lm = p_lm, mlog10p_rna_lm = mlog10p_lm,
    p_rna_gam = pg$p, mlog10p_rna_gam = pg$mlog10p, edf_rna_gam = pg$edf,
    time_rna_cor_s = t_cor,
    time_rna_lm_s = t_lm,
    time_rna_gam_s = t_gam,
    time_rna_total_s = t_tot
  )
}

.fit_fp_models_one_unit <- function(tf, gene, fp_peak, df_unit, p_floor = 1e-300, k = 5L, bs = "tp") {
  t0 <- proc.time()[[3]]
  t_cor0 <- proc.time()[[3]]
  corv <- .safe_cor(df_unit$tf_fp, df_unit$expr, floor = p_floor)
  t_cor <- as.numeric(proc.time()[[3]] - t_cor0)

  t_lm0 <- proc.time()[[3]]
  fit_lm <- tryCatch(stats::lm(expr ~ tf_fp + cell + stress, data = df_unit), error = function(e) NULL)
  beta <- NA_real_; p_lm <- NA_real_; mlog10p_lm <- NA_real_
  if (!is.null(fit_lm)) {
    co <- summary(fit_lm)$coefficients
    if ("tf_fp" %in% rownames(co)) {
      beta <- co["tf_fp", "Estimate"]
      pv <- .cap_p(co["tf_fp", "Pr(>|t|)"], floor = p_floor)
      p_lm <- pv$p; mlog10p_lm <- pv$mlog10p
    }
  }
  t_lm <- as.numeric(proc.time()[[3]] - t_lm0)

  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Need package 'mgcv'.")
  s <- mgcv::s
  t_gam0 <- proc.time()[[3]]
  fit_gam <- tryCatch(
    mgcv::gam(expr ~ s(tf_fp, k = k, bs = bs) + cell + stress, data = df_unit, method = "REML"),
    error = function(e) NULL
  )
  pg <- .mgcv_smooth_p(fit_gam, floor = p_floor)
  t_gam <- as.numeric(proc.time()[[3]] - t_gam0)
  t_tot <- as.numeric(proc.time()[[3]] - t0)

  tibble::tibble(
    tf = tf, gene = gene, fp_peak = fp_peak, n_samples = nrow(df_unit),
    r_fp = corv$r, p_fp = corv$p, mlog10p_fp = corv$mlog10p,
    beta_fp_lm = beta, p_fp_lm = p_lm, mlog10p_fp_lm = mlog10p_lm,
    p_fp_gam = pg$p, mlog10p_fp_gam = pg$mlog10p, edf_fp_gam = pg$edf,
    time_fp_cor_s = t_cor,
    time_fp_lm_s = t_lm,
    time_fp_gam_s = t_gam,
    time_fp_total_s = t_tot
  )
}

.fit_sem_one_unit_lavaan <- function(tf, gene, fp_peak, df_unit, p_floor = 1e-300) {
  t0 <- proc.time()[[3]]
  if (!requireNamespace("lavaan", quietly = TRUE)) stop("Need package 'lavaan'.")
  # Residualize out covariates first (faster + more robust than adding many dummies to SEM)
  Z <- stats::model.matrix(~ cell + stress, data = df_unit)
  dat <- data.frame(
    expr = as.numeric(.residualize(Z, df_unit$expr)),
    tf_expr = as.numeric(.residualize(Z, df_unit$tf_expr)),
    tf_fp = as.numeric(.residualize(Z, df_unit$tf_fp))
  )
  model_txt <- "tf_fp ~ a*tf_expr\nexpr  ~ b*tf_fp + c*tf_expr\nab := a*b\ntotal := c + (a*b)\n"

  fit <- tryCatch(
    lavaan::sem(
      model_txt,
      data = dat,
      meanstructure = FALSE,
      se = "robust.huber.white",
      test = "satorra.bentler"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble::tibble(
      tf = tf, gene = gene, fp_peak = fp_peak, n_samples = nrow(df_unit),
      sem_a = NA_real_, sem_p_a = NA_real_, sem_mlog10p_a = NA_real_,
      sem_b = NA_real_, sem_p_b = NA_real_, sem_mlog10p_b = NA_real_,
      sem_c = NA_real_, sem_p_c = NA_real_, sem_mlog10p_c = NA_real_,
      sem_ab = NA_real_, sem_p_ab = NA_real_, sem_mlog10p_ab = NA_real_,
      sem_total = NA_real_, sem_p_total = NA_real_, sem_mlog10p_total = NA_real_,
      time_sem_s = as.numeric(proc.time()[[3]] - t0)
    ))
  }

  pe <- lavaan::parameterEstimates(fit)
  get_lbl <- function(lbl) {
    ii <- which(pe$label == lbl & pe$op == "~")
    if (!length(ii)) return(list(est = NA_real_, p = NA_real_, mlog10p = NA_real_))
    pv <- .cap_p(pe$pvalue[ii[1]], floor = p_floor)
    list(est = pe$est[ii[1]], p = pv$p, mlog10p = pv$mlog10p)
  }

  a <- get_lbl("a"); b <- get_lbl("b"); c <- get_lbl("c")
  ii_ab <- which(pe$lhs == "ab" & pe$op == ":=")
  pv_ab <- .cap_p(if (length(ii_ab)) pe$pvalue[ii_ab[1]] else NA_real_, floor = p_floor)
  est_ab <- if (length(ii_ab)) pe$est[ii_ab[1]] else NA_real_
  ii_t <- which(pe$lhs == "total" & pe$op == ":=")
  pv_t <- .cap_p(if (length(ii_t)) pe$pvalue[ii_t[1]] else NA_real_, floor = p_floor)
  est_t <- if (length(ii_t)) pe$est[ii_t[1]] else NA_real_

  tibble::tibble(
    tf = tf, gene = gene, fp_peak = fp_peak, n_samples = nrow(df_unit),
    sem_a = a$est, sem_p_a = a$p, sem_mlog10p_a = a$mlog10p,
    sem_b = b$est, sem_p_b = b$p, sem_mlog10p_b = b$mlog10p,
    sem_c = c$est, sem_p_c = c$p, sem_mlog10p_c = c$mlog10p,
    sem_ab = est_ab, sem_p_ab = pv_ab$p, sem_mlog10p_ab = pv_ab$mlog10p,
    sem_total = est_t, sem_p_total = pv_t$p, sem_mlog10p_total = pv_t$mlog10p,
    time_sem_s = as.numeric(proc.time()[[3]] - t0)
  )
}

run_tf_benchmark_fast <- function(tf,
                                  grn_set,
                                  tf_gene_links_gh,
                                  motif_db,
                                  ko_truth_tbl = NULL,
                                  min_samples = 6L,
                                  k = 5L,
                                  bs = "tp",
                                  p_floor = 1e-300,
                                  run_sem = TRUE,
                                  mc_cores = 1L,
                                  cache_dir = NULL,
                                  use_cache = TRUE,
                                  cache_tag = "v1",
                                  cache_level = c("tf", "gene", "none"),
                                  verbose = TRUE) {
  cache_level <- match.arg(cache_level)
  rna_tbl <- grn_set$rna
  md <- grn_set$sample_metadata_used
  stress_col <- if ("stress_type" %in% names(md)) "stress_type" else if ("stress" %in% names(md)) "stress" else NA_character_
  md_use <- md[, c("id", "cell", stress_col), drop = FALSE]
  names(md_use)[names(md_use) == stress_col] <- "stress"
  md_use$cell <- factor(md_use$cell)
  md_use$stress <- factor(md_use$stress)
  sample_ids <- as.character(md_use$id)

  map <- .tf_peak_gene_map(tf, tf_gene_links_gh, motif_db)
  genes <- sort(unique(map$gene))
  genes <- intersect(genes, rna_tbl$HGNC)

  cache_key <- paste0(
    "tag=", cache_tag,
    "_min=", as.integer(min_samples),
    "_k=", as.integer(k),
    "_bs=", bs,
    "_pfloor=", format(p_floor, scientific = TRUE),
    "_sem=", as.integer(isTRUE(run_sem)),
    "_ko=", as.integer(!is.null(ko_truth_tbl))
  )
  cache_base <- if (!is.null(cache_dir)) file.path(cache_dir, tf, cache_key) else NULL
  cache_rna_dir <- if (!is.null(cache_base)) file.path(cache_base, "rna") else NULL
  cache_gene_dir <- if (!is.null(cache_base)) file.path(cache_base, if (isTRUE(run_sem)) "fp_sem" else "fp") else NULL
  cache_tf_path <- if (!is.null(cache_base)) file.path(cache_base, "results_tf.rds") else NULL
  timing_path <- if (!is.null(cache_base)) file.path(cache_base, "timing_summary.csv") else NULL

  if (isTRUE(use_cache) && cache_level == "tf" && !is.null(cache_tf_path) && file.exists(cache_tf_path)) {
    hit <- tryCatch(readRDS(cache_tf_path), error = function(e) NULL)
    if (!is.null(hit)) return(hit)
  }

  .cache_read <- function(path) {
    if (isTRUE(use_cache) && !is.null(path) && file.exists(path)) {
      return(tryCatch(readRDS(path), error = function(e) NULL))
    }
    NULL
  }
  .cache_write <- function(path, obj) {
    if (is.null(path)) return(invisible(FALSE))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tryCatch({ saveRDS(obj, path); TRUE }, error = function(e) FALSE)
  }

  if (isTRUE(verbose)) {
    message(sprintf("[tf=%s] genes_rna=%d genes_fp=? mc_cores=%d cache_level=%s", tf, length(genes), as.integer(mc_cores), cache_level))
  }
  t_tf0 <- proc.time()[[3]]

  # RNA models per TF-gene
  rna_fun <- function(g) {
    cache_path <- if (cache_level == "gene" && !is.null(cache_rna_dir)) file.path(cache_rna_dir, paste0(g, ".rds")) else NULL
    hit <- .cache_read(cache_path)
    if (!is.null(hit)) { hit$from_cache <- TRUE; return(hit) }
    res <- .fit_rna_models_one_gene(tf, g, rna_tbl, sample_ids, md_use[, c("cell", "stress")],
                                    min_samples = min_samples, p_floor = p_floor, k = k, bs = bs)
    res$from_cache <- FALSE
    if (cache_level == "gene") .cache_write(cache_path, res)
    res
  }
  if (.Platform$OS.type != "windows" && mc_cores > 1L) {
    chunks <- .chunk_vec(genes, n_chunks = mc_cores)
    rna_parts <- parallel::mclapply(chunks, function(xx) dplyr::bind_rows(lapply(xx, rna_fun)),
                                    mc.cores = mc_cores, mc.preschedule = FALSE)
    rna_res <- dplyr::bind_rows(rna_parts)
  } else {
    rna_res <- dplyr::bind_rows(lapply(genes, rna_fun))
  }

  # FP (+ optional SEM) per TF-gene (each worker handles all peaks for a target gene)
  tf_long <- .build_tf_long_peak_tbl(tf, grn_set, tf_gene_links_gh, motif_db, min_samples = min_samples)
  if (!nrow(tf_long)) return(list(rna = rna_res, fp = tibble::tibble(), sem = tibble::tibble(), all = tibble::tibble()))

  genes_fp <- sort(unique(tf_long$gene))
  if (isTRUE(verbose)) message(sprintf("[tf=%s] genes_fp=%d (after fp join + min_samples filter)", tf, length(genes_fp)))
  gene_fun <- function(g) {
    cache_path <- if (cache_level == "gene" && !is.null(cache_gene_dir)) file.path(cache_gene_dir, paste0(g, ".rds")) else NULL
    hit <- .cache_read(cache_path)
    if (!is.null(hit)) return(hit)

    df_g <- tf_long[tf_long$gene == g, , drop = FALSE]
    by_peak <- split(df_g, df_g$fp_peak)
    fp_rows <- lapply(names(by_peak), function(pk) {
      .fit_fp_models_one_unit(tf, g, pk, by_peak[[pk]], p_floor = p_floor, k = k, bs = bs)
    })
    fp_tbl <- dplyr::bind_rows(fp_rows)

    sem_tbl <- if (isTRUE(run_sem)) {
      sem_rows <- lapply(names(by_peak), function(pk) {
        .fit_sem_one_unit_lavaan(tf, g, pk, by_peak[[pk]], p_floor = p_floor)
      })
      dplyr::bind_rows(sem_rows)
    } else {
      tibble::tibble()
    }
    out <- list(fp = fp_tbl, sem = sem_tbl)
    if (cache_level == "gene") .cache_write(cache_path, out)
    out
  }

  if (.Platform$OS.type != "windows" && mc_cores > 1L) {
    peaks_per_gene <- tapply(tf_long$fp_peak, tf_long$gene, function(x) length(unique(x)))
    w <- as.numeric(peaks_per_gene[genes_fp])
    chunks_fp <- .chunk_by_weight(genes_fp, weights = w, n_chunks = mc_cores)
    if (isTRUE(verbose)) {
      sizes <- vapply(chunks_fp, length, integer(1))
      message(sprintf("[tf=%s] fp parallel: chunks=%d (min/med/max genes per chunk = %d/%d/%d)",
                      tf, length(chunks_fp),
                      min(sizes), stats::median(sizes), max(sizes)))
    }
    res_gene <- parallel::mclapply(chunks_fp, function(xx) lapply(xx, gene_fun),
                                   mc.cores = mc_cores, mc.preschedule = FALSE)
    res_gene <- unlist(res_gene, recursive = FALSE)
  } else {
    res_gene <- lapply(genes_fp, gene_fun)
  }

  fp_res <- dplyr::bind_rows(lapply(res_gene, `[[`, "fp"))
  sem_res <- if (isTRUE(run_sem)) dplyr::bind_rows(lapply(res_gene, `[[`, "sem")) else tibble::tibble()

  # Attach KO truth (gene-level) to all outputs
  if (!is.null(ko_truth_tbl) && is.data.frame(ko_truth_tbl) && all(c("gene", "log2fc", "ko_group") %in% names(ko_truth_tbl))) {
    ii <- match(rna_res$gene, ko_truth_tbl$gene)
    rna_res$log2fc <- ko_truth_tbl$log2fc[ii]
    rna_res$ko_group <- ko_truth_tbl$ko_group[ii]

    jj <- match(fp_res$gene, ko_truth_tbl$gene)
    fp_res$log2fc <- ko_truth_tbl$log2fc[jj]
    fp_res$ko_group <- ko_truth_tbl$ko_group[jj]

    if (nrow(sem_res)) {
      kk <- match(sem_res$gene, ko_truth_tbl$gene)
      sem_res$log2fc <- ko_truth_tbl$log2fc[kk]
      sem_res$ko_group <- ko_truth_tbl$ko_group[kk]
    }
  }

  all_res <- dplyr::left_join(fp_res, rna_res, by = c("tf", "gene"), suffix = c("", "_rna")) |>
    dplyr::left_join(sem_res, by = c("tf", "gene", "fp_peak"), suffix = c("", "_sem"))

  t_tf <- as.numeric(proc.time()[[3]] - t_tf0)
  timing_summary <- tibble::tibble(
    tf = tf,
    seconds_total = t_tf,
    seconds_rna_total = sum(rna_res$time_rna_total_s, na.rm = TRUE),
    seconds_rna_cor = sum(rna_res$time_rna_cor_s, na.rm = TRUE),
    seconds_rna_lm = sum(rna_res$time_rna_lm_s, na.rm = TRUE),
    seconds_rna_gam = sum(rna_res$time_rna_gam_s, na.rm = TRUE),
    seconds_fp_total = sum(fp_res$time_fp_total_s, na.rm = TRUE),
    seconds_fp_cor = sum(fp_res$time_fp_cor_s, na.rm = TRUE),
    seconds_fp_lm = sum(fp_res$time_fp_lm_s, na.rm = TRUE),
    seconds_fp_gam = sum(fp_res$time_fp_gam_s, na.rm = TRUE),
    seconds_sem_total = if (nrow(sem_res)) sum(sem_res$time_sem_s, na.rm = TRUE) else 0
  )
  if (isTRUE(verbose)) {
    message(sprintf("[tf=%s] timing(s): total=%.1f rna=%.1f fp=%.1f sem=%.1f",
                    tf, timing_summary$seconds_total, timing_summary$seconds_rna_total,
                    timing_summary$seconds_fp_total, timing_summary$seconds_sem_total))
  }
  if (!is.null(timing_path)) {
    dir.create(dirname(timing_path), recursive = TRUE, showWarnings = FALSE)
    tryCatch(utils::write.csv(timing_summary, timing_path, row.names = FALSE), error = function(e) NULL)
  }

  out_all <- list(rna = rna_res, fp = fp_res, sem = sem_res, all = all_res, timing = timing_summary)
  if (isTRUE(use_cache) && cache_level == "tf" && !is.null(cache_tf_path)) .cache_write(cache_tf_path, out_all)
  out_all
}

# ---- Run (HNF1A + SOX9 only) -----------------------------------------------
if (!exists("grn_set")) stop("grn_set not found. Load it first: list(fp_score, rna, sample_metadata_used).")
tfs_run <- intersect(tfs_all, c("HNF1A", "SOX9"))

mc_cores_use <- getOption("mc.cores", NA_integer_)
if (!is.finite(mc_cores_use) || is.na(mc_cores_use)) mc_cores_use <- parallel::detectCores()
mc_cores_use <- max(1L, as.integer(mc_cores_use) - 1L)
if (.Platform$OS.type == "windows") mc_cores_use <- 1L

# Strategy: run TFs sequentially, parallelize within each TF across genes/peaks using mc_cores_use.
results_fast <- setNames(lapply(tfs_run, function(tf) {
  run_tf_benchmark_fast(
    tf = tf,
    grn_set = grn_set,
    tf_gene_links_gh = tf_gene_links_gh,
    motif_db = motif_db,
    ko_truth_tbl = ko_truth_list[[tf]] %||% NULL,
    min_samples = 6L,
    k = 5L,
    bs = "tp",
    p_floor = 1e-300,
    run_sem = TRUE,
    mc_cores = mc_cores_use,
    cache_dir = file.path(ko_dir, "cache_tf_benchmark_fast"),
    use_cache = TRUE,
    cache_tag = "v1",
    cache_level = "tf",
    verbose = TRUE
  )
}), tfs_run)

results_rna_all <- dplyr::bind_rows(lapply(results_fast, `[[`, "rna"))
results_fp_all <- dplyr::bind_rows(lapply(results_fast, `[[`, "fp"))
results_sem_all <- dplyr::bind_rows(lapply(results_fast, `[[`, "sem"))
results_fast_all <- dplyr::bind_rows(lapply(results_fast, `[[`, "all"))
timing_all <- dplyr::bind_rows(lapply(results_fast, `[[`, "timing"))

# ---- P-value grid + plotting (lightweight) ----------------------------------
make_p_grid <- function(p_vals = c(1e-2, 1e-3, 1e-4, 1e-5)) {
  tibble::tibble(
    p_cut_rna   = p_vals,
    p_cut_fp    = p_vals,
    p_cut_model = p_vals,
    cutoff_label = paste0("1e-", -log10(p_vals))
  )
}

.sign_from_group <- function(ko_group) {
  dplyr::case_when(
    ko_group == "Up" ~ 1L,
    ko_group == "Down" ~ -1L,
    TRUE ~ 0L
  )
}

.collapse_min_p <- function(df, p_col, sign_col, extra_cols = character(0)) {
  df <- df[is.finite(df[[p_col]]), , drop = FALSE]
  if (!nrow(df)) return(tibble::tibble())
  df |>
    dplyr::group_by(tf, gene) |>
    dplyr::slice_min(order_by = .data[[p_col]], n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      tf, gene,
      p_value = .data[[p_col]],
      sign = sign(.data[[sign_col]]),
      dplyr::across(dplyr::all_of(extra_cols))
    )
}

build_tf_model_prediction_grid <- function(results_rna,
                                          results_fp,
                                          results_sem,
                                          p_grid = make_p_grid(),
                                          r_cut_rna = 0.3,
                                          r_cut_fp = 0.3) {
  need_truth <- c("tf", "gene", "ko_group", "log2fc")
  if (!all(need_truth %in% names(results_rna))) {
    stop("results_rna must include: tf, gene, ko_group, log2fc (attach KO truth before calling).")
  }

  truth <- results_rna |>
    dplyr::select(tf, gene, ko_group, log2fc) |>
    dplyr::distinct()

  rna_corr <- results_rna |>
    dplyr::transmute(tf, gene, model_label = "RNA: corr",
                     p_value = p_rna, r_value = r_rna, sign = sign(r_rna))
  rna_lm <- results_rna |>
    dplyr::transmute(tf, gene, model_label = "RNA: linear",
                     p_value = p_rna_lm, r_value = NA_real_, sign = sign(beta_rna_lm))
  rna_gam <- results_rna |>
    dplyr::transmute(tf, gene, model_label = "RNA: GAM",
                     p_value = p_rna_gam, r_value = NA_real_, sign = sign(beta_rna_lm))

  fp_corr0 <- results_fp |>
    dplyr::filter(is.finite(r_fp) & abs(r_fp) >= r_cut_fp)
  fp_corr <- .collapse_min_p(fp_corr0, p_col = "p_fp", sign_col = "r_fp") |>
    dplyr::mutate(model_label = "FP: corr", r_value = NA_real_)
  fp_lm <- .collapse_min_p(results_fp, p_col = "p_fp_lm", sign_col = "beta_fp_lm") |>
    dplyr::mutate(model_label = "FP: linear", r_value = NA_real_)
  fp_gam <- .collapse_min_p(results_fp, p_col = "p_fp_gam", sign_col = "beta_fp_lm") |>
    dplyr::mutate(model_label = "FP: GAM", r_value = NA_real_)

  sem_total <- if (nrow(results_sem)) {
    .collapse_min_p(results_sem, p_col = "sem_p_total", sign_col = "sem_total") |>
      dplyr::mutate(model_label = "RNA+FP: SEM total", r_value = NA_real_)
  } else {
    tibble::tibble()
  }

  models <- dplyr::bind_rows(
    rna_corr, rna_lm, rna_gam,
    fp_corr, fp_lm, fp_gam,
    sem_total
  ) |>
    dplyr::left_join(truth, by = c("tf", "gene")) |>
    dplyr::mutate(
      truth_sign = .sign_from_group(ko_group),
      sign_ok = (truth_sign == 0L) | (sign == truth_sign)
    )

  out <- dplyr::bind_rows(lapply(seq_len(nrow(p_grid)), function(i) {
    pc_rna <- p_grid$p_cut_rna[i]
    pc_fp <- p_grid$p_cut_fp[i]
    pc_m <- p_grid$p_cut_model[i]
    lab <- as.character(p_grid$cutoff_label[i])

    models |>
      dplyr::mutate(
        cutoff_label = lab,
        p_cut = dplyr::case_when(
          grepl("^RNA:", model_label) ~ pc_rna,
          grepl("^FP:", model_label) ~ pc_fp,
          TRUE ~ pc_m
        ),
        pass_p = is.finite(p_value) & p_value <= p_cut,
        pass_r = dplyr::case_when(
          model_label == "RNA: corr" ~ is.finite(r_value) & abs(r_value) >= r_cut_rna,
          TRUE ~ TRUE
        ),
        predicted = pass_p & pass_r & sign_ok
      )
  }))
  out$cutoff_label <- factor(out$cutoff_label, levels = as.character(p_grid$cutoff_label))
  desired_models <- c("FP: corr", "FP: GAM", "FP: linear",
                      "RNA: corr", "RNA: GAM", "RNA: linear",
                      "RNA+FP: SEM total")
  out$model_label <- factor(out$model_label, levels = intersect(desired_models, unique(as.character(out$model_label))))
  out
}

filter_tf_pred_grid_by_TFLink <- function(pred_grid_tbl,
                                         tf_name,
                                         tf_link_tbl,
                                         semi = FALSE,
                                         max_non_tf_link = 500L) {
  need <- c("Name.TF", "Name.Target")
  if (!all(need %in% names(tf_link_tbl))) stop("tf_link_tbl must contain: Name.TF, Name.Target")
  tf_targets <- unique(tf_link_tbl$Name.Target[tf_link_tbl$Name.TF == tf_name])
  tf_targets <- tf_targets[!is.na(tf_targets) & tf_targets != ""]
  if (!length(tf_targets)) return(pred_grid_tbl[0, , drop = FALSE])

  in_link <- pred_grid_tbl$gene %in% tf_targets
  if (!isTRUE(semi)) return(pred_grid_tbl[in_link, , drop = FALSE])

  non <- pred_grid_tbl[!in_link, , drop = FALSE]
  if (nrow(non) > max_non_tf_link) non <- non[sample.int(nrow(non), max_non_tf_link), , drop = FALSE]
  dplyr::bind_rows(pred_grid_tbl[in_link, , drop = FALSE], non)
}

plot_tf_model_comparison_grid <- function(df_model,
                                         tf_label = "HNF1A",
                                         out_file = NULL,
                                         lfc_breaks = c(-Inf, -1, -0.5, 0, Inf),
                                         lfc_labels = c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need package 'ggplot2'.")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Need package 'patchwork'.")
  stopifnot(all(c("predicted", "log2fc", "cutoff_label") %in% names(df_model)))

  df_use <- df_model[is.finite(df_model$log2fc), , drop = FALSE]
  if (!nrow(df_use)) return(invisible(NULL))
  df_use$predicted[is.na(df_use$predicted)] <- FALSE
  df_use <- df_use[!is.na(df_use$cutoff_label) & df_use$cutoff_label != "", , drop = FALSE]
  if (!nrow(df_use)) return(invisible(NULL))

  # Normalize expected columns from old pipeline
  if (!("model" %in% names(df_use))) {
    if ("model_label" %in% names(df_use)) {
      df_use$model <- df_use$model_label
    } else {
      stop("df_model must include either 'model' or 'model_label'.")
    }
  }
  if (!("gene" %in% names(df_use))) stop("df_model must include 'gene'.")

  # Preserve order from input factors where possible (matches old plotting behaviour)
  model_levels <- if (is.factor(df_use$model)) levels(df_use$model) else unique(df_use$model)
  df_use$model <- factor(df_use$model, levels = model_levels)
  cutoff_levels <- if (is.factor(df_use$cutoff_label)) levels(df_use$cutoff_label) else unique(df_use$cutoff_label)
  df_use$cutoff_label <- factor(df_use$cutoff_label, levels = cutoff_levels)

  df_use$group <- ifelse(df_use$predicted, "Predicted", "Non-predicted")
  df_use$group <- factor(df_use$group, levels = c("Predicted", "Non-predicted"))

  x_levels <- as.character(unlist(lapply(cutoff_levels, function(cl) paste(cl, levels(df_use$group), sep = " | "))))
  df_use$x_group <- factor(paste(df_use$cutoff_label, df_use$group, sep = " | "), levels = x_levels)

  df_use$log2fc_bin <- cut(
    df_use$log2fc,
    breaks = lfc_breaks,
    labels = lfc_labels,
    right = TRUE,
    include.lowest = TRUE
  )

  lfc_cols <- c(
    "<= -1" = "#d73027",
    "(-1,-0.5]" = "#fc8d59",
    "(-0.5,0]" = "#fee090",
    ">= 0" = "#d9d9d9"
  )

  x_lab_fun <- function(x) {
    parts <- strsplit(x, " | ", fixed = TRUE)
    vapply(parts, function(xx) {
      grp <- xx[2]
      cut <- xx[1]
      grp_short <- if (grp == "Predicted") "Pred" else "Non"
      paste(grp_short, cut, sep = "\n")
    }, character(1L))
  }

  p_violin <- ggplot2::ggplot(df_use, ggplot2::aes(x = x_group, y = log2fc, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.4, colour = "black") +
    ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.2, outlier.alpha = 0.4) +
    ggplot2::facet_grid(rows = ggplot2::vars(1), cols = ggplot2::vars(model)) +
    ggplot2::scale_fill_manual(values = c("Predicted" = "#66c2a5", "Non-predicted" = "#fc8d62"), name = "Group") +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::labs(
      title = sprintf("%s KO - Model comparison across p-value cutoffs", tf_label),
      x = NULL,
      y = "log2FC (KO vs Ctrl)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  df_counts <- df_use[df_use$predicted, , drop = FALSE] |>
    dplyr::distinct(model, gene, cutoff_label) |>
    dplyr::count(model, cutoff_label, name = "n_predicted")
  df_counts <- tidyr::complete(df_counts, model = levels(df_use$model), cutoff_label = levels(df_use$cutoff_label),
                               fill = list(n_predicted = 0))

  p_counts <- ggplot2::ggplot(df_counts, ggplot2::aes(x = cutoff_label, y = n_predicted)) +
    ggplot2::geom_col(fill = "grey70") +
    ggplot2::geom_text(ggplot2::aes(label = n_predicted), vjust = -0.3, size = 2.5) +
    ggplot2::facet_grid(rows = ggplot2::vars(1), cols = ggplot2::vars(model)) +
    ggplot2::labs(x = NULL, y = "Number of predicted target genes") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5)
    )

  df_percent <- df_use[!is.na(df_use$log2fc_bin), , drop = FALSE] |>
    dplyr::distinct(model, cutoff_label, gene, group, log2fc_bin) |>
    dplyr::count(model, cutoff_label, group, log2fc_bin, name = "n") |>
    dplyr::group_by(model, cutoff_label, group) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::ungroup()
  df_percent$x_group <- factor(paste(df_percent$cutoff_label, df_percent$group, sep = " | "), levels = x_levels)

  p_percent <- ggplot2::ggplot(df_percent, ggplot2::aes(x = x_group, y = pct, fill = log2fc_bin)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(rows = ggplot2::vars(1), cols = ggplot2::vars(model)) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_x_discrete(labels = x_lab_fun) +
    ggplot2::scale_fill_manual(values = lfc_cols, name = "log2FC bin") +
    ggplot2::labs(x = "p-value cutoff x Group", y = "Percent of genes within group") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      plot.margin = ggplot2::margin(t = 0, r = 5, b = 5, l = 5),
      legend.position = "right"
    )

  p_all <- p_violin / p_counts / p_percent +
    patchwork::plot_layout(heights = c(3, 1, 1))

  if (!is.null(out_file)) {
    ggplot2::ggsave(filename = out_file, plot = p_all, width = 16, height = 8.5, units = "in", dpi = 300)
  }
  invisible(p_all)
}

# Example: rebuild p-grid plots like before
p_vals <- c(1e-2, 1e-3, 1e-4, 1e-5)
p_grid_all <- make_p_grid(p_vals)

for (tf in tfs_run) {
  res <- results_fast[[tf]]
  if (!nrow(res$rna) || !nrow(res$fp)) next

  tf_pred_grid <- build_tf_model_prediction_grid(
    results_rna = res$rna,
    results_fp = res$fp,
    results_sem = res$sem,
    p_grid = p_grid_all,
    r_cut_rna = 0.3,
    r_cut_fp = 0.3
  )

  tf_pred_grid_tfLink_semi <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl = tf_pred_grid,
    tf_name = tf,
    tf_link_tbl = TFLink,
    semi = TRUE
  )

  plot_tf_model_comparison_grid(
    df_model = tf_pred_grid_tfLink_semi,
    tf_label = tf,
    out_file = file.path(ko_dir, sprintf("%s_model_comparison_p_grid.pdf", tf))
  )

  tf_pred_grid_tfLink <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl = tf_pred_grid,
    tf_name = tf,
    tf_link_tbl = TFLink,
    semi = FALSE
  )

  plot_tf_model_comparison_grid(
    df_model = tf_pred_grid_tfLink,
    tf_label = paste0(tf, " (TFLink only)"),
    out_file = file.path(ko_dir, sprintf("%s_model_comparison_p_grid_TFLink.pdf", tf))
  )
}


# ---- New idea: benchmark using edges_filtered_tidy gene sets -----------------
# Goal: for each TF, take the unique predicted gene list from edges_filtered_tidy
# (tf/gene_key only), then compare KO log2fc distributions:
# - predicted genes (N = n_distinct gene_key for that TF)
# - random genes (N sampled from ko_group == "Unchanged" in that TF's KO table)
#
# This intentionally does NOT use any other columns from edges_filtered_tidy.
options(episcope.save_pdfs = FALSE)
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
                                               return_details = FALSE,
                                               do_plot = TRUE) {
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

  if (!isTRUE(return_details)) {
    if (isTRUE(do_plot)) {
      p <- .plot_predicted_vs_random_3row(
        ko_tbl = ko_tbl,
        tf_label = tf0,
        predicted_genes = pred_in_ko[seq_len(n_use)],
        random_genes = rand_genes,
        out_file = out_file,
        random_label = "Random",
        max_n_counts = max_n_counts
      )
      return(invisible(p))
    }
    return(invisible(NULL))
  }

  .group_stats <- function(ko_tbl, genes, group) {
    ko_map <- ko_tbl |>
      dplyr::transmute(gene_norm = .norm_chr(.data$gene), log2fc = .data$log2fc) |>
      dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
      dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
    genes <- unique(.norm_chr(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    lfc <- ko_map$log2fc[match(genes, ko_map$gene_norm)]
    lfc <- lfc[is.finite(lfc)]
    n_tot <- length(lfc)

    # Match the plotting bins used in row 3:
    # breaks = c(-Inf, -1, -0.5, 0, Inf), right=TRUE, include.lowest=TRUE
    n_leq_m1 <- sum(lfc <= -1)
    n_m1_to_m0.5 <- sum(lfc > -1 & lfc <= -0.5)
    n_m0.5_to_0 <- sum(lfc > -0.5 & lfc <= 0)

    pct_leq_m1 <- if (n_tot > 0) 100 * n_leq_m1 / n_tot else NA_real_
    pct_m1_to_m0.5 <- if (n_tot > 0) 100 * n_m1_to_m0.5 / n_tot else NA_real_
    pct_m0.5_to_0 <- if (n_tot > 0) 100 * n_m0.5_to_0 / n_tot else NA_real_
    pct_nonpos <- if (n_tot > 0) pct_leq_m1 + pct_m1_to_m0.5 + pct_m0.5_to_0 else NA_real_

    tibble::tibble(
      group = group,
      n_genes = n_tot,
      # Percent of genes in the "non-positive" bins, matching the requested figure-2 style.
      pct_of_ko = pct_nonpos,
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

  if (!isTRUE(do_plot)) {
    return(invisible(list(plot = NULL, stats = stats_tbl)))
  }

  p <- .plot_predicted_vs_random_3row(
    ko_tbl = ko_tbl,
    tf_label = tf0,
    predicted_genes = pred_in_ko[seq_len(n_use)],
    random_genes = rand_genes,
    out_file = out_file,
    random_label = "Random",
    max_n_counts = max_n_counts
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
                                                      return_details = FALSE,
                                                      do_plot = TRUE) {
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
    return_details = return_details,
    do_plot = do_plot
  )
}


  # ---- Run across multiple lighting_* folders ------------------------------
  lighting_folder0 <- if (exists("lighting_folder", inherits = TRUE)) {
    get("lighting_folder", inherits = TRUE)
  } else {
    Sys.getenv("EPISCOPE_LIGHTING_FOLDER")
  }
  if (!is.character(lighting_folder0) || !nzchar(lighting_folder0)) {
    stop("Need `lighting_folder` or env var `EPISCOPE_LIGHTING_FOLDER` to locate the lighting_* folders.")
  }

  lighting_parent <- dirname(lighting_folder0)
  lighting_key <- "lighting_fp_tf_corr_FDR_0.05_genehancer_jaspar2024_regulated_genes_1.5_delta_link_1"
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

    do_plots <- isTRUE(getOption("episcope.save_pdfs", TRUE)) &&
      !identical(Sys.getenv("EPISCOPE_SKIP_PDFS"), "1")
    message("[run] do_plots=", do_plots, " (set option episcope.save_pdfs or env EPISCOPE_SKIP_PDFS=1)")

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
        return_details = TRUE,
        do_plot = do_plots
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

    if (isTRUE(do_plots) && length(tf_plots)) {
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
    } else if (isTRUE(do_plots)) {
      message("No TF plots were produced; nothing to save.")
    }

    # ---- TFLink strict + semi (add back up to 200) -------------------------
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
        return_details = TRUE,
        do_plot = do_plots
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

    if (isTRUE(do_plots) && length(tf_plots_tfl)) {
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
    } else if (isTRUE(do_plots)) {
      message("No TFLink-filtered TF plots were produced; nothing to save.")
    }

    # TFLink semi: add back up to 200 smallest p-value filtered-out genes
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
        return_details = TRUE,
        do_plot = do_plots
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

    if (isTRUE(do_plots) && length(tf_plots_tfl_semi)) {
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
    } else if (isTRUE(do_plots)) {
      message("No TFLink semi-filtered TF plots were produced; nothing to save.")
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
