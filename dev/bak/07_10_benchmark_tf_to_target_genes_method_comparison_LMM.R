progressr::handlers(global = TRUE)  # text bar by default; can change handlers later
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"
motif_db <- readr::read_tsv("/data/homes/yl814/package/episcope/inst/extdata/genome/JASPAR2024.txt")

# since the HNF4A is not in the "tf_perturb.db", will load it separately
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

TFLink <- readr::read_tsv("/data/homes/cy232/genome/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv")
tf_perturb_db <- "/data/homes/yl814/episcope/tf_perturb.db"
predicted_tfbs_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
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
tfs_all <- c( "HNF1A", "HNF4A", "IRF1", "RARG", "SOX9", "KLF5", "FOXA2")





# build TF-gene-peak table (no aggregation over peaks)

#' Build long-format FP + RNA table for one TF (GeneHancer links only, per peak)
#'
#' For the given TF, we:
#'   - use motif_db to find its motifs
#'   - keep only tf_gene_links_gh rows with those motifs
#'   - join FP scores and RNA expression without aggregating across peaks
#'
#' @param tf Character scalar, TF symbol (e.g. "HNF1A", "SOX9").
#' @param grn_set List containing fp_score, rna, sample_metadata_used.
#' @param tf_gene_links_gh Tibble/data.frame with GeneHancer TF-gene links.
#'        Must have columns: fp_peak, gene_key, motifs.
#' @param motif_db Tibble with columns motif, cluster, HGNC.
#' @param min_peaks_per_gene Integer, minimum FP peaks per TF-gene to keep.
#'
#' @return tibble with columns:
#'   gene, fp_peak, sample_id, expr, fp_score,
#'   cell, stress_type
build_tf_gene_fp_expr_tbl <- function(tf,
                                      grn_set,
                                      tf_gene_links_gh,
                                      motif_db,
                                      min_peaks_per_gene = 1L) {
  fp_score  <- grn_set$fp_score
  rna_tbl   <- grn_set$rna
  sample_md <- grn_set$sample_metadata_used

  sample_ids <- sample_md$id

  # sanity checks
  missing_fp <- setdiff(sample_ids, names(fp_score))
  missing_rn <- setdiff(sample_ids, names(rna_tbl))
  if (length(missing_fp)) {
    cli::cli_abort("Missing FP columns in fp_score: {.val {missing_fp}}")
  }
  if (length(missing_rn)) {
    cli::cli_abort("Missing RNA columns in rna: {.val {missing_rn}}")
  }
  if (!all(c("fp_peak", "gene_key", "motifs") %in% names(tf_gene_links_gh))) {
    cli::cli_abort("tf_gene_links_gh must have columns 'fp_peak', 'gene_key', 'motifs'.")
  }
  if (!all(c("motif", "HGNC") %in% names(motif_db))) {
    cli::cli_abort("motif_db must have columns 'motif' and 'HGNC'.")
  }
  if (!"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("grn_set$rna must have column 'HGNC'.")
  }
  if (!"peak_ID" %in% names(fp_score)) {
    cli::cli_abort("grn_set$fp_score must have column 'peak_ID'.")
  }

  # ----------------------------
  # 1) TF-specific GeneHancer map via motif_db
  # ----------------------------
  motifs_tf <- motif_db$motif[motif_db$HGNC == tf]
  motifs_tf <- unique(motifs_tf[!is.na(motifs_tf)])

  if (!length(motifs_tf)) {
    cli::cli_warn("No motifs found in motif_db for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  tf_links <- tf_gene_links_gh[
    tf_gene_links_gh$motifs %in% motifs_tf &
      !is.na(tf_gene_links_gh$gene_key) &
      tf_gene_links_gh$gene_key != "",
    c("fp_peak", "gene_key", "motifs"),
    drop = FALSE
  ]

  if (!nrow(tf_links)) {
    cli::cli_warn("No GeneHancer links found for TF {.val {tf}} after motif filter.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  genes_rna  <- rna_tbl$HGNC
  genes_tf   <- unique(tf_links$gene_key)
  genes_keep <- intersect(genes_tf, genes_rna)

  tf_links <- tf_links[tf_links$gene_key %in% genes_keep, , drop = FALSE]
  if (!nrow(tf_links)) {
    cli::cli_warn("No overlapping genes between TF links and RNA for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  # ----------------------------
  # 2) FP per peak (no aggregation yet)
  # ----------------------------
  peaks_use <- unique(tf_links$fp_peak)
  fp_sub    <- fp_score[fp_score$peak_ID %in% peaks_use, c("peak_ID", sample_ids), drop = FALSE]

  if (!nrow(fp_sub)) {
    cli::cli_warn("No FP peaks from fp_score matched GeneHancer links for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  fp_long <- tidyr::pivot_longer(
    data      = fp_sub,
    cols      = dplyr::all_of(sample_ids),
    names_to  = "sample_id",
    values_to = "fp_score"
  )

  fp_long <- dplyr::left_join(
    fp_long,
    tf_links,
    by = c("peak_ID" = "fp_peak")
  )

  fp_long <- fp_long[
    !is.na(fp_long$gene_key) & fp_long$gene_key %in% genes_keep,
    ,
    drop = FALSE
  ]

  if (!nrow(fp_long)) {
    cli::cli_warn("After joining, no FP entries remain for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  # require at least min_peaks_per_gene
  peaks_per_gene <- fp_long |>
    dplyr::distinct(gene_key, peak_ID) |>
    dplyr::count(gene_key, name = "n_peaks")

  keep_genes2 <- peaks_per_gene$gene_key[peaks_per_gene$n_peaks >= min_peaks_per_gene]
  fp_long     <- fp_long[fp_long$gene_key %in% keep_genes2, , drop = FALSE]

  if (!nrow(fp_long)) {
    cli::cli_warn("No TF-gene pairs passed min_peaks_per_gene for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  # ----------------------------
  # 3) RNA and sample metadata
  # ----------------------------
  rna_sub <- rna_tbl[rna_tbl$HGNC %in% genes_keep, c("HGNC", sample_ids), drop = FALSE]

  rna_long <- tidyr::pivot_longer(
    data      = rna_sub,
    cols      = dplyr::all_of(sample_ids),
    names_to  = "sample_id",
    values_to = "expr"
  )
  names(rna_long)[names(rna_long) == "HGNC"] <- "gene"

  tf_long <- dplyr::inner_join(
    fp_long,
    rna_long,
    by = c("gene_key" = "gene", "sample_id" = "sample_id")
  )

  if (!nrow(tf_long)) {
    cli::cli_warn("No overlapping FP+RNA rows after join for TF {.val {tf}}.")
    return(tibble::tibble(
      gene        = character(0),
      fp_peak     = character(0),
      sample_id   = character(0),
      expr        = numeric(0),
      fp_score    = numeric(0),
      cell        = character(0),
      stress_type = character(0)
    ))
  }

  sample_md_use <- sample_md[, c("id", "cell", "stress_type"), drop = FALSE]
  names(sample_md_use)[names(sample_md_use) == "id"] <- "sample_id"

  tf_long <- dplyr::left_join(
    tf_long,
    sample_md_use,
    by = "sample_id"
  )

  # final tidy tibble
  tibble::tibble(
    gene        = tf_long$gene_key,
    fp_peak     = tf_long$peak_ID,
    sample_id   = tf_long$sample_id,
    expr        = tf_long$expr,
    fp_score    = tf_long$fp_score,
    cell        = tf_long$cell,
    stress_type = tf_long$stress_type
  )
}

#' Compute per-gene/per-peak correlation + LMM p-values (ENET disabled)
#'
#' For the given TF, we:
#'   0) Precompute per-gene best correlation p-value from tf_gene_links_gh$p_fp
#'      restricted to this TF via motif_db (p_corr is gene-level).
#'   1) Build per-peak FP+RNA long table via build_tf_gene_fp_expr_tbl().
#'   2) Precompute TF expression per sample and attach as tf_expr.
#'   3) For each gene, split rows by fp_peak and fit LMMs:
#'        (a) expr ~ fp_score + (1 | cell) + (1 | stress_type)
#'        (b) expr ~ tf_expr  + (1 | cell) + (1 | stress_type)
#'        (c) expr ~ fp_score + tf_expr + (1 | cell) + (1 | stress_type)
#'        (d) expr ~ 1 + (1 | cell) + (1 | stress_type)  [null model for LRT]
#'      returning one row per (gene, fp_peak) with:
#'        p_lmm_fp          : p for fp_score in (a)
#'        p_lmm_rna         : p for tf_expr in (b)
#'        p_lmm_both_fp     : p for fp_score in (c)
#'        p_lmm_both_rna    : p for tf_expr in (c)
#'        p_lmm_both_overall: LRT p-value comparing (c) vs (d)
#'   4) (ENET p-values are currently disabled and returned as NA).
#'
#' Parallelization:
#'   - If use_future = TRUE, uses parallel::mclapply() across genes.
#'   - Otherwise, runs a base for-loop with a text progress bar.
#'
#' @param tf Character scalar, TF symbol.
#' @param grn_set List with fp_score, rna, sample_metadata_used.
#' @param tf_gene_links_gh Tibble with columns at least:
#'        fp_peak, gene_key, motifs, p_fp.
#' @param motif_db Tibble with columns motif, HGNC.
#' @param min_samples Integer, minimum sample size per peak to fit LMM.
#' @param use_future Logical; if TRUE, use parallel::mclapply().
#' @param future_chunk_size Ignored (kept for API compatibility).
#' @param cache_dir Optional directory to write chunked CSV results.
#' @param cache_chunk_size Integer number of rows per CSV chunk when caching.
#' @param mc_cores Integer number of cores for mclapply (NULL = detectCores() - 1).
#'
#' @return tibble with columns:
#'   gene, fp_peak, p_corr,
#'   p_lmm_fp, p_lmm_rna, p_lmm_both_fp, p_lmm_both_rna, p_lmm_both_overall, p_enet.
compute_tf_pvals_lmm_enet <- function(tf,
                                      grn_set,
                                      tf_gene_links_gh,
                                      motif_db,
                                      min_samples       = 6L,
                                      use_future        = FALSE,
                                      future_chunk_size = NULL,  # kept for compatibility
                                      cache_dir         = NULL,
                                      cache_chunk_size  = 500L,
                                      mc_cores          = NULL) {

  # ------------------------------------------------------------------
  # 0) Precompute per-gene correlation p from tf_gene_links_gh$p_fp
  #     restricted to this TF via motif_db
  # ------------------------------------------------------------------
  if (!all(c("fp_peak", "gene_key", "motifs", "p_fp") %in% names(tf_gene_links_gh))) {
    cli::cli_abort("tf_gene_links_gh must have columns 'fp_peak', 'gene_key', 'motifs', 'p_fp'.")
  }
  if (!all(c("motif", "HGNC") %in% names(motif_db))) {
    cli::cli_abort("motif_db must have columns 'motif' and 'HGNC'.")
  }

  rna_tbl <- grn_set$rna
  if (!"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("grn_set$rna must have column 'HGNC'.")
  }

  motifs_tf <- unique(motif_db$motif[motif_db$HGNC == tf])
  motifs_tf <- motifs_tf[!is.na(motifs_tf)]
  if (!length(motifs_tf)) {
    cli::cli_warn("No motifs found in motif_db for TF {.val {tf}}.")
    return(tibble::tibble(
      gene               = character(0),
      fp_peak            = character(0),
      p_corr             = numeric(0),
      p_lmm_fp           = numeric(0),
      p_lmm_rna          = numeric(0),
      p_lmm_both_fp      = numeric(0),
      p_lmm_both_rna     = numeric(0),
      p_lmm_both_overall = numeric(0),
      p_enet             = numeric(0)
    ))
  }

  genes_rna <- unique(rna_tbl$HGNC)

  tf_links_corr <- tf_gene_links_gh[
    tf_gene_links_gh$motifs %in% motifs_tf &
      !is.na(tf_gene_links_gh$gene_key) &
      tf_gene_links_gh$gene_key != "" &
      is.finite(tf_gene_links_gh$p_fp),
    c("fp_peak", "gene_key", "p_fp"),
    drop = FALSE
  ]
  tf_links_corr <- tf_links_corr[tf_links_corr$gene_key %in% genes_rna, , drop = FALSE]

  if (!nrow(tf_links_corr)) {
    cli::cli_warn("No valid p_fp entries for TF {.val {tf}} after motif/RNA filter.")
    return(tibble::tibble(
      gene               = character(0),
      fp_peak            = character(0),
      p_corr             = numeric(0),
      p_lmm_fp           = numeric(0),
      p_lmm_rna          = numeric(0),
      p_lmm_both_fp      = numeric(0),
      p_lmm_both_rna     = numeric(0),
      p_lmm_both_overall = numeric(0),
      p_enet             = numeric(0)
    ))
  }

  corr_by_gene <- tf_links_corr |>
    dplyr::group_by(gene_key) |>
    dplyr::summarise(
      p_corr = min(p_fp, na.rm = TRUE),
      .groups = "drop"
    )
  p_corr_map <- stats::setNames(corr_by_gene$p_corr, corr_by_gene$gene_key)

  # ------------------------------------------------------------------
  # 1) Build per-peak FP+RNA long table (for LMM)
  # ------------------------------------------------------------------
  tf_long <- build_tf_gene_fp_expr_tbl(
    tf                 = tf,
    grn_set            = grn_set,
    tf_gene_links_gh   = tf_gene_links_gh,
    motif_db           = motif_db,
    min_peaks_per_gene = 1L
  )

  if (!nrow(tf_long)) {
    cli::cli_warn("No data to compute p-values for TF {.val {tf}}.")
    return(tibble::tibble(
      gene               = character(0),
      fp_peak            = character(0),
      p_corr             = numeric(0),
      p_lmm_fp           = numeric(0),
      p_lmm_rna          = numeric(0),
      p_lmm_both_fp      = numeric(0),
      p_lmm_both_rna     = numeric(0),
      p_lmm_both_overall = numeric(0),
      p_enet             = numeric(0)
    ))
  }

  # restrict to genes that have p_corr defined
  tf_long <- tf_long[tf_long$gene %in% names(p_corr_map), , drop = FALSE]
  if (!nrow(tf_long)) {
    cli::cli_warn("After intersecting with corr_by_gene, no rows remain for TF {.val {tf}}.")
    return(tibble::tibble(
      gene               = character(0),
      fp_peak            = character(0),
      p_corr             = numeric(0),
      p_lmm_fp           = numeric(0),
      p_lmm_rna          = numeric(0),
      p_lmm_both_fp      = numeric(0),
      p_lmm_both_rna     = numeric(0),
      p_lmm_both_overall = numeric(0),
      p_enet             = numeric(0)
    ))
  }

  tf_long$cell        <- factor(tf_long$cell)
  tf_long$stress_type <- factor(tf_long$stress_type)

  # ------------------------------------------------------------------
  # 1b) Attach TF expression per sample as tf_expr
  # ------------------------------------------------------------------
  # We take TF expression from grn_set$rna, averaging across rows if multiple entries.
  sample_cols_rna <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  tf_rows <- which(rna_tbl$HGNC == tf)

  tf_expr_available <- length(tf_rows) > 0L && length(sample_cols_rna) > 0L

  if (tf_expr_available) {
    tf_mat <- as.matrix(rna_tbl[tf_rows, sample_cols_rna, drop = FALSE])
    tf_expr_vec <- colMeans(tf_mat, na.rm = TRUE)

    tf_expr_df <- tibble::tibble(
      sample_id = sample_cols_rna,
      tf_expr   = as.numeric(tf_expr_vec)
    )

    tf_long <- dplyr::left_join(
      tf_long,
      tf_expr_df,
      by = "sample_id"
    )
  } else {
    cli::cli_warn("TF {.val {tf}} not found (or no sample columns) in rna_tbl; TF-expression LMMs will be NA.")
    tf_long$tf_expr <- NA_real_
  }

  # Pre-split by gene once to avoid repeated subsetting
  gene_splits <- split(tf_long, tf_long$gene)
  genes       <- names(gene_splits)

  # ------------------------------------------------------------------
  # 2) Per-gene LMM, per-peak p-values
  # ------------------------------------------------------------------
  compute_for_gene_df <- function(df_gene) {
    gene_name <- df_gene$gene[1L]

    sub_gene <- df_gene[
      is.finite(df_gene$expr) & is.finite(df_gene$fp_score),
      ,
      drop = FALSE
    ]

    if (!nrow(sub_gene)) {
      return(tibble::tibble(
        gene               = character(0),
        fp_peak            = character(0),
        p_corr             = numeric(0),
        p_lmm_fp           = numeric(0),
        p_lmm_rna          = numeric(0),
        p_lmm_both_fp      = numeric(0),
        p_lmm_both_rna     = numeric(0),
        p_lmm_both_overall = numeric(0),
        p_enet             = numeric(0)
      ))
    }

    p_corr_val <- if (gene_name %in% names(p_corr_map)) p_corr_map[[gene_name]] else NA_real_

    peaks <- unique(sub_gene$fp_peak)
    out   <- vector("list", length(peaks))

    for (j in seq_along(peaks)) {
      sub_peak <- sub_gene[sub_gene$fp_peak == peaks[j], , drop = FALSE]
      sub_peak <- sub_peak[
        is.finite(sub_peak$expr) & is.finite(sub_peak$fp_score),
        ,
        drop = FALSE
      ]

      # also require finite tf_expr for TF-expression models
      has_tf_expr <- tf_expr_available && any(is.finite(sub_peak$tf_expr))

      p_lmm_fp_val           <- NA_real_
      p_lmm_rna_val          <- NA_real_
      p_lmm_both_fp_val      <- NA_real_
      p_lmm_both_rna_val     <- NA_real_
      p_lmm_both_overall_val <- NA_real_

      ok_n <- nrow(sub_peak) >= min_samples &&
        length(unique(sub_peak$cell))        >= 2L &&
        length(unique(sub_peak$stress_type)) >= 2L

      if (ok_n) {
        # 2a) FP-only model: expr ~ fp_score + (1|cell) + (1|stress_type)
        lmm_fp <- tryCatch(
          lmerTest::lmer(
            expr ~ fp_score + (1 | cell) + (1 | stress_type),
            data = sub_peak,
            REML = FALSE
          ),
          error = function(e) NULL
        )
        if (!is.null(lmm_fp)) {
          coef_tab <- summary(lmm_fp)$coefficients
          if ("fp_score" %in% rownames(coef_tab)) {
            p_lmm_fp_val <- coef_tab["fp_score", "Pr(>|t|)"]
          }
        }

        # 2b) TF-expression only model: expr ~ tf_expr + (1|cell) + (1|stress_type)
        if (has_tf_expr) {
          sub_peak_tf <- sub_peak[is.finite(sub_peak$tf_expr), , drop = FALSE]
          if (nrow(sub_peak_tf) >= min_samples) {
            lmm_rna <- tryCatch(
              lmerTest::lmer(
                expr ~ tf_expr + (1 | cell) + (1 | stress_type),
                data = sub_peak_tf,
                REML = FALSE
              ),
              error = function(e) NULL
            )
            if (!is.null(lmm_rna)) {
              coef_tab_rna <- summary(lmm_rna)$coefficients
              if ("tf_expr" %in% rownames(coef_tab_rna)) {
                p_lmm_rna_val <- coef_tab_rna["tf_expr", "Pr(>|t|)"]
              }
            }
          }
        }

        # 2c) Combined model: expr ~ fp_score + tf_expr + (1|cell) + (1|stress_type)
        # 2d) Null model for LRT: expr ~ 1 + (1|cell) + (1|stress_type)
        if (has_tf_expr) {
          sub_peak_both <- sub_peak[
            is.finite(sub_peak$tf_expr),
            ,
            drop = FALSE
          ]
          if (nrow(sub_peak_both) >= min_samples) {
            # Fit null model (intercept only)
            lmm_null <- tryCatch(
              lmerTest::lmer(
                expr ~ 1 + (1 | cell) + (1 | stress_type),
                data = sub_peak_both,
                REML = FALSE
              ),
              error = function(e) NULL
            )

            # Fit full model (both predictors)
            lmm_both <- tryCatch(
              lmerTest::lmer(
                expr ~ fp_score + tf_expr + (1 | cell) + (1 | stress_type),
                data = sub_peak_both,
                REML = FALSE
              ),
              error = function(e) NULL
            )

            if (!is.null(lmm_both)) {
              coef_tab_both <- summary(lmm_both)$coefficients
              if ("fp_score" %in% rownames(coef_tab_both)) {
                p_lmm_both_fp_val <- coef_tab_both["fp_score", "Pr(>|t|)"]
              }
              if ("tf_expr" %in% rownames(coef_tab_both)) {
                p_lmm_both_rna_val <- coef_tab_both["tf_expr", "Pr(>|t|)"]
              }

              # Likelihood ratio test for overall model significance
              if (!is.null(lmm_null)) {
                lrt <- tryCatch(
                  anova(lmm_null, lmm_both),
                  error = function(e) NULL
                )
                if (!is.null(lrt) && nrow(lrt) >= 2) {
                  p_lmm_both_overall_val <- lrt$`Pr(>Chisq)`[2]
                }
              }
            }
          }
        }
      }

      out[[j]] <- tibble::tibble(
        gene               = gene_name,
        fp_peak            = peaks[j],
        p_corr             = p_corr_val,
        p_lmm_fp           = p_lmm_fp_val,
        p_lmm_rna          = p_lmm_rna_val,
        p_lmm_both_fp      = p_lmm_both_fp_val,
        p_lmm_both_rna     = p_lmm_both_rna_val,
        p_lmm_both_overall = p_lmm_both_overall_val,
        p_enet             = NA_real_  # ENET disabled for now
      )
    }

    dplyr::bind_rows(out)
  }

  # ------------------------------------------------------------------
  # 3) Loop over genes: sequential or mclapply
  # ------------------------------------------------------------------
  if (use_future) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      cli::cli_abort("use_future=TRUE but package 'parallel' is not available.")
    }
    if (is.null(mc_cores)) {
      mc_cores <- max(1L, parallel::detectCores() - 1L)
    }

    cli::cli_inform("Running per-gene LMM with parallel::mclapply() on {mc_cores} cores.")

    res_list <- parallel::mclapply(
      X            = gene_splits,
      FUN          = compute_for_gene_df,
      mc.cores     = mc_cores,
      mc.preschedule = TRUE
    )

  } else {
    pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
    on.exit(close(pb), add = TRUE)

    res_list <- vector("list", length(genes))
    for (i in seq_along(genes)) {
      utils::setTxtProgressBar(pb, i)
      res_list[[i]] <- compute_for_gene_df(gene_splits[[i]])
    }
  }

  res_df <- dplyr::bind_rows(res_list)

  # ------------------------------------------------------------------
  # 4) Optional caching to CSV in chunks
  # ------------------------------------------------------------------
  if (!is.null(cache_dir)) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    n <- nrow(res_df)
    if (n > 0L) {
      if (is.null(cache_chunk_size) || cache_chunk_size <= 0L) {
        cache_chunk_size <- n
      }
      n_chunks  <- ceiling(n / cache_chunk_size)
      base_name <- sprintf("TF_%s_LMM_pvals", tf)

      for (k in seq_len(n_chunks)) {
        idx_start <- (k - 1L) * cache_chunk_size + 1L
        idx_end   <- min(k * cache_chunk_size, n)
        chunk     <- res_df[idx_start:idx_end, , drop = FALSE]

        outfile_chunk <- file.path(
          cache_dir,
          sprintf("%s_chunk%03d.csv", base_name, k)
        )
        readr::write_csv(chunk, outfile_chunk)
      }

      cli::cli_inform(
        "Wrote {n_chunks} cached LMM chunk(s) for TF {.val {tf}} to {.path {cache_dir}}."
      )
    }
  }

  res_df
}


#' Compute per-gene/per-peak correlation + Fixed Effects Model p-values
#'
#' For the given TF, we:
#'   0) Precompute per-gene best correlation p-value from tf_gene_links_gh$p_fp
#'      restricted to this TF via motif_db (p_corr is gene-level).
#'   1) Build per-peak FP+RNA long table via build_tf_gene_fp_expr_tbl().
#'   2) Precompute TF expression per sample and attach as tf_expr.
#'   3) For each gene, split rows by fp_peak and fit linear models with fixed effects:
#'        (a) expr ~ fp_score + cell + stress_type
#'        (b) expr ~ tf_expr  + cell + stress_type
#'        (c) expr ~ fp_score + tf_expr + cell + stress_type
#'        (d) expr ~ cell + stress_type  [null model for F-test]
#'
#'      Optionally with interactions:
#'        (a_int) expr ~ fp_score * cell + fp_score * stress_type
#'        (b_int) expr ~ tf_expr  * cell + tf_expr  * stress_type
#'        (c_int) expr ~ fp_score * cell + fp_score * stress_type +
#'                       tf_expr  * cell + tf_expr  * stress_type
#'
#'      returning one row per (gene, fp_peak) with:
#'        p_fixed_fp          : p for fp_score in (a)
#'        p_fixed_rna         : p for tf_expr in (b)
#'        p_fixed_both_fp     : p for fp_score in (c)
#'        p_fixed_both_rna    : p for tf_expr in (c)
#'        p_fixed_both_overall: F-test p-value comparing (c) vs (d)
#'
#'        If use_interactions = TRUE, also returns:
#'        p_fixed_fp_int      : F-test for all fp_score terms in (a_int)
#'        p_fixed_rna_int     : F-test for all tf_expr terms in (b_int)
#'        p_fixed_both_fp_int : F-test for all fp_score terms in (c_int)
#'        p_fixed_both_rna_int: F-test for all tf_expr terms in (c_int)
#'
#' Parallelization:
#'   - If use_future = TRUE, uses parallel::mclapply() across genes.
#'   - Otherwise, runs a base for-loop with a text progress bar.
#'
#' @param tf Character scalar, TF symbol.
#' @param grn_set List with fp_score, rna, sample_metadata_used.
#' @param tf_gene_links_gh Tibble with columns at least:
#'        fp_peak, gene_key, motifs, p_fp.
#' @param motif_db Tibble with columns motif, HGNC.
#' @param min_samples Integer, minimum sample size per peak to fit model.
#' @param use_interactions Logical; if TRUE, include interaction terms.
#' @param use_future Logical; if TRUE, use parallel::mclapply().
#' @param cache_dir Optional directory to write chunked CSV results.
#' @param cache_chunk_size Integer number of rows per CSV chunk when caching.
#' @param mc_cores Integer number of cores for mclapply (NULL = detectCores() - 1).
#'
#' @return tibble with columns:
#'   gene, fp_peak, p_corr,
#'   p_fixed_fp, p_fixed_rna, p_fixed_both_fp, p_fixed_both_rna, p_fixed_both_overall
#'   If use_interactions=TRUE, also: p_fixed_fp_int, p_fixed_rna_int,
#'   p_fixed_both_fp_int, p_fixed_both_rna_int
compute_tf_pvals_fixed_effects <- function(tf,
                                           grn_set,
                                           tf_gene_links_gh,
                                           motif_db,
                                           min_samples       = 6L,
                                           use_interactions  = FALSE,
                                           use_future        = FALSE,
                                           cache_dir         = NULL,
                                           cache_chunk_size  = 500L,
                                           mc_cores          = NULL) {

  # ------------------------------------------------------------------
  # 0) Precompute per-gene correlation p from tf_gene_links_gh$p_fp
  # ------------------------------------------------------------------
  if (!all(c("fp_peak", "gene_key", "motifs", "p_fp") %in% names(tf_gene_links_gh))) {
    cli::cli_abort("tf_gene_links_gh must have columns 'fp_peak', 'gene_key', 'motifs', 'p_fp'.")
  }
  if (!all(c("motif", "HGNC") %in% names(motif_db))) {
    cli::cli_abort("motif_db must have columns 'motif' and 'HGNC'.")
  }

  rna_tbl <- grn_set$rna
  if (!"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("grn_set$rna must have column 'HGNC'.")
  }

  motifs_tf <- unique(motif_db$motif[motif_db$HGNC == tf])
  motifs_tf <- motifs_tf[!is.na(motifs_tf)]

  # Build empty return tibble structure
  empty_cols <- c(
    "gene", "fp_peak", "p_corr",
    "p_fixed_fp", "p_fixed_rna", "p_fixed_both_fp",
    "p_fixed_both_rna", "p_fixed_both_overall"
  )
  if (use_interactions) {
    empty_cols <- c(
      empty_cols,
      "p_fixed_fp_int", "p_fixed_rna_int",
      "p_fixed_both_fp_int", "p_fixed_both_rna_int"
    )
  }

  empty_return <- tibble::tibble(
    gene = character(0),
    fp_peak = character(0),
    p_corr = numeric(0),
    p_fixed_fp = numeric(0),
    p_fixed_rna = numeric(0),
    p_fixed_both_fp = numeric(0),
    p_fixed_both_rna = numeric(0),
    p_fixed_both_overall = numeric(0)
  )
  if (use_interactions) {
    empty_return$p_fixed_fp_int <- numeric(0)
    empty_return$p_fixed_rna_int <- numeric(0)
    empty_return$p_fixed_both_fp_int <- numeric(0)
    empty_return$p_fixed_both_rna_int <- numeric(0)
  }

  if (!length(motifs_tf)) {
    cli::cli_warn("No motifs found in motif_db for TF {.val {tf}}.")
    return(empty_return)
  }

  genes_rna <- unique(rna_tbl$HGNC)

  tf_links_corr <- tf_gene_links_gh[
    tf_gene_links_gh$motifs %in% motifs_tf &
      !is.na(tf_gene_links_gh$gene_key) &
      tf_gene_links_gh$gene_key != "" &
      is.finite(tf_gene_links_gh$p_fp),
    c("fp_peak", "gene_key", "p_fp"),
    drop = FALSE
  ]
  tf_links_corr <- tf_links_corr[tf_links_corr$gene_key %in% genes_rna, , drop = FALSE]

  if (!nrow(tf_links_corr)) {
    cli::cli_warn("No valid p_fp entries for TF {.val {tf}} after motif/RNA filter.")
    return(empty_return)
  }

  corr_by_gene <- tf_links_corr |>
    dplyr::group_by(gene_key) |>
    dplyr::summarise(
      p_corr = min(p_fp, na.rm = TRUE),
      .groups = "drop"
    )
  p_corr_map <- stats::setNames(corr_by_gene$p_corr, corr_by_gene$gene_key)

  # ------------------------------------------------------------------
  # 1) Build per-peak FP+RNA long table
  # ------------------------------------------------------------------
  tf_long <- build_tf_gene_fp_expr_tbl(
    tf                 = tf,
    grn_set            = grn_set,
    tf_gene_links_gh   = tf_gene_links_gh,
    motif_db           = motif_db,
    min_peaks_per_gene = 1L
  )

  if (!nrow(tf_long)) {
    cli::cli_warn("No data to compute p-values for TF {.val {tf}}.")
    return(empty_return)
  }

  # restrict to genes that have p_corr defined
  tf_long <- tf_long[tf_long$gene %in% names(p_corr_map), , drop = FALSE]
  if (!nrow(tf_long)) {
    cli::cli_warn("After intersecting with corr_by_gene, no rows remain for TF {.val {tf}}.")
    return(empty_return)
  }

  tf_long$cell        <- factor(tf_long$cell)
  tf_long$stress_type <- factor(tf_long$stress_type)

  # ------------------------------------------------------------------
  # 2) Attach TF expression per sample as tf_expr
  # ------------------------------------------------------------------
  sample_cols_rna <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  tf_rows <- which(rna_tbl$HGNC == tf)

  tf_expr_available <- length(tf_rows) > 0L && length(sample_cols_rna) > 0L

  if (tf_expr_available) {
    tf_mat <- as.matrix(rna_tbl[tf_rows, sample_cols_rna, drop = FALSE])
    tf_expr_vec <- colMeans(tf_mat, na.rm = TRUE)

    tf_expr_df <- tibble::tibble(
      sample_id = sample_cols_rna,
      tf_expr   = as.numeric(tf_expr_vec)
    )

    tf_long <- dplyr::left_join(
      tf_long,
      tf_expr_df,
      by = "sample_id"
    )
  } else {
    cli::cli_warn("TF {.val {tf}} not found (or no sample columns) in rna_tbl; TF-expression models will be NA.")
    tf_long$tf_expr <- NA_real_
  }

  # Pre-split by gene once
  gene_splits <- split(tf_long, tf_long$gene)
  genes       <- names(gene_splits)

  # ------------------------------------------------------------------
  # 3) Per-gene fixed effects models, per-peak p-values
  # ------------------------------------------------------------------
  compute_for_gene_df <- function(df_gene) {
    gene_name <- df_gene$gene[1L]

    sub_gene <- df_gene[
      is.finite(df_gene$expr) & is.finite(df_gene$fp_score),
      ,
      drop = FALSE
    ]

    if (!nrow(sub_gene)) {
      result <- tibble::tibble(
        gene = character(0),
        fp_peak = character(0),
        p_corr = numeric(0),
        p_fixed_fp = numeric(0),
        p_fixed_rna = numeric(0),
        p_fixed_both_fp = numeric(0),
        p_fixed_both_rna = numeric(0),
        p_fixed_both_overall = numeric(0)
      )
      if (use_interactions) {
        result$p_fixed_fp_int <- numeric(0)
        result$p_fixed_rna_int <- numeric(0)
        result$p_fixed_both_fp_int <- numeric(0)
        result$p_fixed_both_rna_int <- numeric(0)
      }
      return(result)
    }

    p_corr_val <- if (gene_name %in% names(p_corr_map)) p_corr_map[[gene_name]] else NA_real_

    peaks <- unique(sub_gene$fp_peak)
    out   <- vector("list", length(peaks))

    for (j in seq_along(peaks)) {
      sub_peak <- sub_gene[sub_gene$fp_peak == peaks[j], , drop = FALSE]
      sub_peak <- sub_peak[
        is.finite(sub_peak$expr) & is.finite(sub_peak$fp_score),
        ,
        drop = FALSE
      ]

      has_tf_expr <- tf_expr_available && any(is.finite(sub_peak$tf_expr))

      p_fixed_fp_val           <- NA_real_
      p_fixed_rna_val          <- NA_real_
      p_fixed_both_fp_val      <- NA_real_
      p_fixed_both_rna_val     <- NA_real_
      p_fixed_both_overall_val <- NA_real_

      p_fixed_fp_int_val       <- NA_real_
      p_fixed_rna_int_val      <- NA_real_
      p_fixed_both_fp_int_val  <- NA_real_
      p_fixed_both_rna_int_val <- NA_real_

      ok_n <- nrow(sub_peak) >= min_samples &&
        length(unique(sub_peak$cell))        >= 2L &&
        length(unique(sub_peak$stress_type)) >= 2L

      if (ok_n) {
        # Model (a): expr ~ fp_score + cell + stress_type
        lm_fp <- tryCatch(
          lm(expr ~ fp_score + cell + stress_type, data = sub_peak),
          error = function(e) NULL
        )
        if (!is.null(lm_fp)) {
          coef_tab <- summary(lm_fp)$coefficients
          if ("fp_score" %in% rownames(coef_tab)) {
            p_fixed_fp_val <- coef_tab["fp_score", "Pr(>|t|)"]
          }
        }

        # Model (a_int): expr ~ fp_score * cell + fp_score * stress_type
        if (use_interactions && !is.null(lm_fp)) {
          lm_fp_int <- tryCatch(
            lm(expr ~ fp_score * cell + fp_score * stress_type, data = sub_peak),
            error = function(e) NULL
          )
          if (!is.null(lm_fp_int)) {
            # F-test for all fp_score terms (main + interactions)
            lm_fp_null <- tryCatch(
              lm(expr ~ cell + stress_type, data = sub_peak),
              error = function(e) NULL
            )
            if (!is.null(lm_fp_null)) {
              f_test <- tryCatch(
                anova(lm_fp_null, lm_fp_int, test = "F"),
                error = function(e) NULL
              )
              if (!is.null(f_test) && nrow(f_test) >= 2) {
                p_fixed_fp_int_val <- f_test$`Pr(>F)`[2]
              }
            }
          }
        }

        # Model (b): expr ~ tf_expr + cell + stress_type
        if (has_tf_expr) {
          sub_peak_rna <- sub_peak[is.finite(sub_peak$tf_expr), , drop = FALSE]
          if (nrow(sub_peak_rna) >= min_samples) {
            lm_rna <- tryCatch(
              lm(expr ~ tf_expr + cell + stress_type, data = sub_peak_rna),
              error = function(e) NULL
            )
            if (!is.null(lm_rna)) {
              coef_tab_rna <- summary(lm_rna)$coefficients
              if ("tf_expr" %in% rownames(coef_tab_rna)) {
                p_fixed_rna_val <- coef_tab_rna["tf_expr", "Pr(>|t|)"]
              }
            }

            # Model (b_int): expr ~ tf_expr * cell + tf_expr * stress_type
            if (use_interactions && !is.null(lm_rna)) {
              lm_rna_int <- tryCatch(
                lm(expr ~ tf_expr * cell + tf_expr * stress_type, data = sub_peak_rna),
                error = function(e) NULL
              )
              if (!is.null(lm_rna_int)) {
                lm_rna_null <- tryCatch(
                  lm(expr ~ cell + stress_type, data = sub_peak_rna),
                  error = function(e) NULL
                )
                if (!is.null(lm_rna_null)) {
                  f_test <- tryCatch(
                    anova(lm_rna_null, lm_rna_int, test = "F"),
                    error = function(e) NULL
                  )
                  if (!is.null(f_test) && nrow(f_test) >= 2) {
                    p_fixed_rna_int_val <- f_test$`Pr(>F)`[2]
                  }
                }
              }
            }
          }
        }

        # Model (c): expr ~ fp_score + tf_expr + cell + stress_type
        # Model (d): expr ~ cell + stress_type (null)
        if (has_tf_expr) {
          sub_peak_both <- sub_peak[is.finite(sub_peak$tf_expr), , drop = FALSE]
          if (nrow(sub_peak_both) >= min_samples) {

            lm_null <- tryCatch(
              lm(expr ~ cell + stress_type, data = sub_peak_both),
              error = function(e) NULL
            )

            lm_both <- tryCatch(
              lm(expr ~ fp_score + tf_expr + cell + stress_type, data = sub_peak_both),
              error = function(e) NULL
            )

            if (!is.null(lm_both)) {
              coef_tab_both <- summary(lm_both)$coefficients
              if ("fp_score" %in% rownames(coef_tab_both)) {
                p_fixed_both_fp_val <- coef_tab_both["fp_score", "Pr(>|t|)"]
              }
              if ("tf_expr" %in% rownames(coef_tab_both)) {
                p_fixed_both_rna_val <- coef_tab_both["tf_expr", "Pr(>|t|)"]
              }

              # F-test for overall model
              if (!is.null(lm_null)) {
                f_test <- tryCatch(
                  anova(lm_null, lm_both, test = "F"),
                  error = function(e) NULL
                )
                if (!is.null(f_test) && nrow(f_test) >= 2) {
                  p_fixed_both_overall_val <- f_test$`Pr(>F)`[2]
                }
              }
            }

            # Model (c_int): expr ~ fp_score*cell + fp_score*stress_type +
            #                       tf_expr*cell + tf_expr*stress_type
            if (use_interactions && !is.null(lm_both)) {
              lm_both_int <- tryCatch(
                lm(expr ~ fp_score * cell + fp_score * stress_type +
                     tf_expr * cell + tf_expr * stress_type,
                   data = sub_peak_both),
                error = function(e) NULL
              )

              if (!is.null(lm_both_int)) {
                # F-test for all fp_score terms in combined model
                lm_both_no_fp <- tryCatch(
                  lm(expr ~ tf_expr * cell + tf_expr * stress_type + cell + stress_type,
                     data = sub_peak_both),
                  error = function(e) NULL
                )
                if (!is.null(lm_both_no_fp)) {
                  f_test_fp <- tryCatch(
                    anova(lm_both_no_fp, lm_both_int, test = "F"),
                    error = function(e) NULL
                  )
                  if (!is.null(f_test_fp) && nrow(f_test_fp) >= 2) {
                    p_fixed_both_fp_int_val <- f_test_fp$`Pr(>F)`[2]
                  }
                }

                # F-test for all tf_expr terms in combined model
                lm_both_no_rna <- tryCatch(
                  lm(expr ~ fp_score * cell + fp_score * stress_type + cell + stress_type,
                     data = sub_peak_both),
                  error = function(e) NULL
                )
                if (!is.null(lm_both_no_rna)) {
                  f_test_rna <- tryCatch(
                    anova(lm_both_no_rna, lm_both_int, test = "F"),
                    error = function(e) NULL
                  )
                  if (!is.null(f_test_rna) && nrow(f_test_rna) >= 2) {
                    p_fixed_both_rna_int_val <- f_test_rna$`Pr(>F)`[2]
                  }
                }
              }
            }
          }
        }
      }

      result_row <- tibble::tibble(
        gene                 = gene_name,
        fp_peak              = peaks[j],
        p_corr               = p_corr_val,
        p_fixed_fp           = p_fixed_fp_val,
        p_fixed_rna          = p_fixed_rna_val,
        p_fixed_both_fp      = p_fixed_both_fp_val,
        p_fixed_both_rna     = p_fixed_both_rna_val,
        p_fixed_both_overall = p_fixed_both_overall_val
      )

      if (use_interactions) {
        result_row$p_fixed_fp_int        <- p_fixed_fp_int_val
        result_row$p_fixed_rna_int       <- p_fixed_rna_int_val
        result_row$p_fixed_both_fp_int   <- p_fixed_both_fp_int_val
        result_row$p_fixed_both_rna_int  <- p_fixed_both_rna_int_val
      }

      out[[j]] <- result_row
    }

    dplyr::bind_rows(out)
  }

  # ------------------------------------------------------------------
  # 4) Loop over genes: sequential or mclapply
  # ------------------------------------------------------------------
  if (use_future) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      cli::cli_abort("use_future=TRUE but package 'parallel' is not available.")
    }
    if (is.null(mc_cores)) {
      mc_cores <- max(1L, parallel::detectCores() - 1L)
    }

    cli::cli_inform("Running per-gene fixed effects models with parallel::mclapply() on {mc_cores} cores.")

    res_list <- parallel::mclapply(
      X            = gene_splits,
      FUN          = compute_for_gene_df,
      mc.cores     = mc_cores,
      mc.preschedule = TRUE
    )

  } else {
    pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
    on.exit(close(pb), add = TRUE)

    res_list <- vector("list", length(genes))
    for (i in seq_along(genes)) {
      utils::setTxtProgressBar(pb, i)
      res_list[[i]] <- compute_for_gene_df(gene_splits[[i]])
    }
  }

  res_df <- dplyr::bind_rows(res_list)

  # ------------------------------------------------------------------
  # 5) Optional caching to CSV in chunks
  # ------------------------------------------------------------------
  if (!is.null(cache_dir)) {
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }

    n <- nrow(res_df)
    if (n > 0L) {
      if (is.null(cache_chunk_size) || cache_chunk_size <= 0L) {
        cache_chunk_size <- n
      }
      n_chunks  <- ceiling(n / cache_chunk_size)
      base_name <- sprintf("TF_%s_FixedEffects_pvals", tf)

      for (k in seq_len(n_chunks)) {
        idx_start <- (k - 1L) * cache_chunk_size + 1L
        idx_end   <- min(k * cache_chunk_size, n)
        chunk     <- res_df[idx_start:idx_end, , drop = FALSE]

        outfile_chunk <- file.path(
          cache_dir,
          sprintf("%s_chunk%03d.csv", base_name, k)
        )
        readr::write_csv(chunk, outfile_chunk)
      }

      cli::cli_inform(
        "Wrote {n_chunks} cached fixed effects chunk(s) for TF {.val {tf}} to {.path {cache_dir}}."
      )
    }
  }

  res_df
}

# RNA-only TF vs gene correlation p-values (unchanged from old logic)
compute_tf_rna_corr_simple <- function(tf,
                                       rna_tbl,
                                       method = c("pearson", "spearman")) {
  method <- match.arg(method)

  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl))) {
    cli::cli_abort("rna_tbl must contain 'ensembl_gene_id' and 'HGNC'.")
  }

  sample_cols <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  if (!length(sample_cols)) {
    cli::cli_abort("No sample columns found in rna_tbl.")
  }

  rna_mat <- as.matrix(rna_tbl[, sample_cols, drop = FALSE])

  tf_rows <- which(rna_tbl$HGNC == tf)
  if (!length(tf_rows)) {
    cli::cli_warn("TF {.val {tf}} not found in rna_tbl$HGNC.")
    return(tibble::tibble(
      gene  = character(0),
      p_rna = numeric(0)
    ))
  }

  tf_expr <- colMeans(rna_mat[tf_rows, , drop = FALSE])
  tf_expr <- suppressWarnings(as.numeric(tf_expr))

  if (sum(is.finite(tf_expr)) < 3L || stats::sd(tf_expr, na.rm = TRUE) == 0) {
    cli::cli_warn("Insufficient variation in TF expression for {.val {tf}}.")
    return(tibble::tibble(
      gene  = character(0),
      p_rna = numeric(0)
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
      suppressWarnings(stats::cor(
        x_num[ok],
        tf_expr[ok],
        use    = "complete.obs",
        method = method
      ))
    }
  )

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

  tibble::tibble(
    gene  = rna_tbl$HGNC,
    p_rna = as.numeric(p_vec)
  ) |>
    dplyr::filter(!is.na(gene), gene != "")
}

# Collapse per-gene/per-peak p-values to one row per gene
# We always use the *overall* combined-model p-value to pick the "best" peak:
#   - If p_lmm_both_overall exists, use that (mixed model).
#   - Else if p_fixed_both_overall exists, use that (fixed effects).
collapse_tf_pvals_by_peak <- function(pval_tbl, ko_truth_tbl) {
  # Basic columns always required
  basic_required <- c("gene", "fp_peak", "p_corr")
  if (!all(basic_required %in% names(pval_tbl))) {
    cli::cli_abort("pval_tbl must contain columns {.val {basic_required}}.")
  }
  if (!all(c("gene", "ko_group") %in% names(ko_truth_tbl))) {
    cli::cli_abort("ko_truth_tbl must contain 'gene' and 'ko_group'.")
  }

  # Decide which overall p-value column to use
  overall_col <- NULL
  if ("p_lmm_both_overall" %in% names(pval_tbl)) {
    overall_col <- "p_lmm_both_overall"
  } else if ("p_fixed_both_overall" %in% names(pval_tbl)) {
    overall_col <- "p_fixed_both_overall"
  }

  if (is.null(overall_col)) {
    cli::cli_abort(
      "pval_tbl must contain either 'p_lmm_both_overall' or 'p_fixed_both_overall'."
    )
  }

  ko_tbl <- ko_truth_tbl |>
    dplyr::filter(
      !is.na(ko_group),
      ko_group %in% c("Down", "Unchanged", "Up")
    ) |>
    dplyr::select(gene, ko_group)

  merged <- pval_tbl |>
    dplyr::inner_join(ko_tbl, by = "gene")

  if (!nrow(merged)) {
    cli::cli_abort("No overlap between pval_tbl and ko_truth_tbl genes.")
  }

  # Collapse to one peak per gene × KO group using overall_col
  collapsed <- merged |>
    dplyr::group_by(gene, ko_group) |>
    dplyr::group_modify(function(df, key) {
      p_overall <- df[[overall_col]]

      if (all(is.na(p_overall))) {
        df[1L, , drop = FALSE]
      } else {
        best <- min(p_overall, na.rm = TRUE)
        df_sub <- df[p_overall == best, , drop = FALSE]
        df_sub[1L, , drop = FALSE]
      }
    }) |>
    dplyr::ungroup()

  collapsed
}

collapse_tf_pvals_by_peak <- function(pval_tbl, ko_truth_tbl) {
  # Basic columns always required
  basic_required <- c("gene", "fp_peak", "p_corr")
  if (!all(basic_required %in% names(pval_tbl))) {
    cli::cli_abort("pval_tbl must contain columns {.val {basic_required}}.")
  }
  if (!all(c("gene", "ko_group") %in% names(ko_truth_tbl))) {
    cli::cli_abort("ko_truth_tbl must contain 'gene' and 'ko_group'.")
  }

  # Decide which overall p-value column to use
  overall_col <- NULL
  if ("p_lmm_both_overall" %in% names(pval_tbl)) {
    overall_col <- "p_lmm_both_overall"
  } else if ("p_fixed_both_overall" %in% names(pval_tbl)) {
    overall_col <- "p_fixed_both_overall"
  }

  if (is.null(overall_col)) {
    cli::cli_abort(
      "pval_tbl must contain either 'p_lmm_both_overall' or 'p_fixed_both_overall'."
    )
  }

  ko_tbl <- ko_truth_tbl |>
    dplyr::filter(
      !is.na(ko_group),
      ko_group %in% c("Down", "Unchanged", "Up")
    ) |>
    dplyr::select(gene, ko_group)

  merged <- pval_tbl |>
    dplyr::inner_join(ko_tbl, by = "gene")

  if (!nrow(merged)) {
    cli::cli_abort("No overlap between pval_tbl and ko_truth_tbl genes.")
  }

  # Collapse to one peak per gene × KO group using overall_col
  collapsed <- merged |>
    dplyr::group_by(gene, ko_group) |>
    dplyr::group_modify(function(df, key) {
      ko_status <- key$ko_group

      # Special case: For "Unchanged" genes with combined models, use first peak
      if (ko_status == "Unchanged") {
        return(df[1L, , drop = FALSE])
      }

      # For "Down" and "Up" genes: select peak with best (minimum) overall p-value
      p_overall <- df[[overall_col]]
      if (all(is.na(p_overall))) {
        df[1L, , drop = FALSE]
      } else {
        best <- min(p_overall, na.rm = TRUE)
        df_sub <- df[p_overall == best, , drop = FALSE]
        df_sub[1L, , drop = FALSE]
      }
    }) |>
    dplyr::ungroup()

  collapsed
}

# Build boxplot tibble from per-peak pvals (after collapsing per gene)
# model_type:
#   "lmm"       -> use p_lmm_* columns         (current behaviour)
#   "fixed"     -> use p_fixed_* columns       (no interactions)
#   "fixed_int" -> use p_fixed_*_int columns   (interaction effects)
build_tf_boxplot_df_from_peaks <- function(tf,
                                           pval_tbl,
                                           rna_tbl,
                                           ko_truth_tbl,
                                           model_type = c("lmm", "fixed", "fixed_int")) {
  model_type <- match.arg(model_type)

  if (!nrow(pval_tbl)) {
    cli::cli_abort("pval_tbl is empty for TF {.val {tf}}.")
  }

  # 1) Collapse per-gene/per-peak rows -> one row per gene
  collapsed <- collapse_tf_pvals_by_peak(
    pval_tbl     = pval_tbl,
    ko_truth_tbl = ko_truth_tbl
  )

  # 2) RNA-only TF vs gene p-values
  rna_corr <- compute_tf_rna_corr_simple(tf = tf, rna_tbl = rna_tbl)
  if (!nrow(rna_corr)) {
    cli::cli_warn("No RNA-only correlations for TF {.val {tf}}.")
  }

  merged <- collapsed |>
    dplyr::inner_join(rna_corr, by = "gene")

  if (!nrow(merged)) {
    cli::cli_abort(
      "No overlap among collapsed FP pvals, RNA pvals, and KO truth for {.val {tf}}."
    )
  }

  safe_log10 <- function(p) {
    p[!is.finite(p) | p <= 0] <- NA_real_
    -log10(p)
  }

  # Decide which columns to use based on model_type
  if (model_type == "lmm") {
    needed <- c("p_lmm_fp", "p_lmm_rna",
                "p_lmm_both_fp", "p_lmm_both_rna", "p_lmm_both_overall")
    if (!all(needed %in% names(merged))) {
      cli::cli_abort("For model_type = 'lmm', merged must contain {.val {needed}}.")
    }
    single_label   <- "LMM (single)"
    combined_label <- "LMM (combined)"

    col_fp_single   <- "p_lmm_fp"
    col_rna_single  <- "p_lmm_rna"
    col_fp_combined <- "p_lmm_both_fp"
    col_rna_combined<- "p_lmm_both_rna"
    col_overall     <- "p_lmm_both_overall"

  } else if (model_type == "fixed") {
    needed <- c("p_fixed_fp", "p_fixed_rna",
                "p_fixed_both_fp", "p_fixed_both_rna", "p_fixed_both_overall")
    if (!all(needed %in% names(merged))) {
      cli::cli_abort("For model_type = 'fixed', merged must contain {.val {needed}}.")
    }
    single_label   <- "Fixed (single)"
    combined_label <- "Fixed (combined)"

    col_fp_single   <- "p_fixed_fp"
    col_rna_single  <- "p_fixed_rna"
    col_fp_combined <- "p_fixed_both_fp"
    col_rna_combined<- "p_fixed_both_rna"
    col_overall     <- "p_fixed_both_overall"

  } else { # model_type == "fixed_int"
    needed_base <- c("p_fixed_both_overall")
    needed_int  <- c("p_fixed_fp_int", "p_fixed_rna_int",
                     "p_fixed_both_fp_int", "p_fixed_both_rna_int")
    if (!all(needed_base %in% names(merged)) ||
        !all(needed_int %in% names(merged))) {
      cli::cli_abort(
        "For model_type = 'fixed_int', merged must contain {.val {c(needed_base, needed_int)}}."
      )
    }
    single_label   <- "Fixed+int (single)"
    combined_label <- "Fixed+int (combined)"

    col_fp_single   <- "p_fixed_fp_int"
    col_rna_single  <- "p_fixed_rna_int"
    col_fp_combined <- "p_fixed_both_fp_int"
    col_rna_combined<- "p_fixed_both_rna_int"
    # overall still uses p_fixed_both_overall (full model vs null)
    col_overall     <- "p_fixed_both_overall"
  }

  # Build long-format boxplot DF
  box_df <- dplyr::bind_rows(
    # --------------------------------------------------------------
    # 1) Correlation: FP vs RNA (same for all model_type)
    # --------------------------------------------------------------
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = "Correlation",
      group      = "FP",
      neg_log10p = safe_log10(merged$p_corr)
    ),
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = "Correlation",
      group      = "RNA",
      neg_log10p = safe_log10(merged$p_rna)
    ),
    # --------------------------------------------------------------
    # 2) Single-predictor models: fp-only vs tf-RNA-only
    # --------------------------------------------------------------
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = single_label,
      group      = "FP",
      neg_log10p = safe_log10(merged[[col_fp_single]])
    ),
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = single_label,
      group      = "RNA",
      neg_log10p = safe_log10(merged[[col_rna_single]])
    ),
    # --------------------------------------------------------------
    # 3) Combined model: fp, tf-RNA, overall LRT
    # --------------------------------------------------------------
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = combined_label,
      group      = "FP",
      neg_log10p = safe_log10(merged[[col_fp_combined]])
    ),
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = combined_label,
      group      = "RNA",
      neg_log10p = safe_log10(merged[[col_rna_combined]])
    ),
    tibble::tibble(
      gene       = merged$gene,
      ko_group   = merged$ko_group,
      method     = combined_label,
      group      = "Overall",
      neg_log10p = safe_log10(merged[[col_overall]])
    )
  )

  box_df <- box_df[is.finite(box_df$neg_log10p), , drop = FALSE]
  if (!nrow(box_df)) {
    cli::cli_abort("All p-values are non-finite for TF {.val {tf}} after log10 transform.")
  }

  # Factor ordering
  box_df$method <- factor(
    box_df$method,
    levels = c("Correlation", single_label, combined_label)
  )

  box_df$group <- factor(
    box_df$group,
    levels = c("FP", "RNA", "Overall")
  )

  box_df$ko_group <- factor(
    box_df$ko_group,
    levels = c("Down", "Unchanged", "Up")
  )

  box_df
}

# Plot boxplots of -log10(p) for FP vs RNA across methods (2 facets)
plot_tf_boxplots_fp_vs_rna <- function(tf,
                                       box_df,
                                       out_dir = ".") {
  if (!nrow(box_df)) {
    cli::cli_abort("box_df is empty for TF {.val {tf}}.")
  }

  # Winsorise per method (cap at 99th percentile for plotting)
  plot_df <- box_df |>
    dplyr::group_by(method) |>
    dplyr::mutate(
      cap_method      = stats::quantile(neg_log10p, probs = 0.99, na.rm = TRUE),
      neg_log10p_plot = pmin(neg_log10p, cap_method)
    ) |>
    dplyr::ungroup()

  plot_df <- plot_df[is.finite(plot_df$neg_log10p_plot), , drop = FALSE]
  if (!nrow(plot_df)) {
    cli::cli_abort("All values became non-finite after winsorisation for TF {.val {tf}}.")
  }

  # KO-group pairwise p-values per (method, group) -> one multi-line label
  wilcox_xy <- function(df, col, g1, g2) {
    x <- df[[col]][df$ko_group == g1]
    y <- df[[col]][df$ko_group == g2]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
  }

  stats_du <- plot_df |>
    dplyr::group_by(method, group) |>
    dplyr::summarise(
      p_DU  = wilcox_xy(dplyr::cur_data_all(), "neg_log10p_plot", "Down", "Unchanged"),
      p_DU2 = wilcox_xy(dplyr::cur_data_all(), "neg_log10p_plot", "Down", "Up"),
      p_UU  = wilcox_xy(dplyr::cur_data_all(), "neg_log10p_plot", "Up",   "Unchanged"),
      y_pos = max(neg_log10p_plot, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      label = paste0(
        "Down vs Non: p=",
        ifelse(is.na(p_DU),  "NA", signif(p_DU,  2)), "\n",
        "Down vs Up: p=",
        ifelse(is.na(p_DU2), "NA", signif(p_DU2, 2)), "\n",
        "Up vs Non: p=",
        ifelse(is.na(p_UU),  "NA", signif(p_UU,  2))
      ),
      y_text = y_pos * 1.15
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = group, y = neg_log10p_plot, fill = ko_group)
  ) +
    ggplot2::geom_boxplot(
      position      = ggplot2::position_dodge(width = 0.7),
      width         = 0.6,
      coef          = 3,
      outlier.size  = 0.2,
      outlier.alpha = 0.4,
      alpha         = 0.9
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(method)
    ) +
    ggplot2::scale_fill_manual(
      values = c(Down = "#4daf4a", Unchanged = "grey60", Up = "#e41a1c"),
      name   = "KO group"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.35))
    ) +
    ggplot2::labs(
      title   = sprintf("%s - GeneHancer FP vs RNA: -log10(p) by method", tf),
      x       = "Predictor / model summary",
      y       = expression(-log[10](p)),
      caption = "For plotting only, -log10(p) values above the 99th percentile are capped per method; boxplot whiskers use coef = 3."
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5)
    ) +
    ggplot2::geom_text(
      data = stats_du,
      ggplot2::aes(
        x     = group,
        y     = y_text,
        label = label
      ),
      inherit.aes = FALSE,
      size  = 2.4,
      hjust = 0.5,
      vjust = 0
    )

  outfile <- file.path(
    out_dir,
    sprintf("TF_%s_GeneHancer_FP_vs_RNA_boxplots_LMM_vs_corr.pdf", tf)
  )

  ggplot2::ggsave(
    filename = outfile,
    plot     = p,
    width    = 12,
    height   = 8,
    units    = "in",
    dpi      = 300
  )

  message("Saved boxplot PDF for ", tf, " to: ", outfile)

  invisible(list(
    box_df  = box_df,
    plot    = p,
    outfile = outfile
  ))
}


# call

future::plan(future::multisession, workers = 30)

pvals_hnf1a <- compute_tf_pvals_lmm_enet(
  tf                = "HNF1A",
  grn_set           = grn_set,
  tf_gene_links_gh  = tf_gene_links_gh,
  motif_db          = motif_db,
  min_samples       = 6L,
  use_future        = TRUE,              # use mclapply
  mc_cores          = 30L,               #
  cache_dir         = file.path(ko_dir, "LMM_pvals_HNF1A_cache"),
  cache_chunk_size  = 500L
)

# Basic fixed effects (no interactions)
pvals_hnf1a_fixed <- compute_tf_pvals_fixed_effects(
  tf = "HNF1A",
  grn_set = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db = motif_db,
  min_samples = 6,
  use_interactions = FALSE,
  use_future = TRUE,
  mc_cores = 30
)

# With interactions (tests if effects vary by condition)
pvals_hnf1a_fixed_int <- compute_tf_pvals_fixed_effects(
  tf = "HNF1A",
  grn_set = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db = motif_db,
  min_samples = 6,
  use_interactions = TRUE,  # Enable interaction terms
  use_future = TRUE,
  mc_cores = 30
)

pvals_sox9 <- compute_tf_pvals_lmm_enet(
  tf               = "SOX9",
  grn_set          = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db         = motif_db,
  min_samples      = 6L,
  use_future        = TRUE,              # use mclapply
  mc_cores          = 30L,
  cache_dir         = file.path(ko_dir, "LMM_pvals_SOX9_cache"),
  cache_chunk_size  = 500L
)

# Basic fixed effects (no interactions)
pvals_sox9_fixed <- compute_tf_pvals_fixed_effects(
  tf = "SOX9",
  grn_set = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db = motif_db,
  min_samples = 6,
  use_interactions = FALSE,
  use_future = TRUE,
  mc_cores = 30
)

# With interactions (tests if effects vary by condition)
pvals_sox9_fixed_int <- compute_tf_pvals_fixed_effects(
  tf = "SOX9",
  grn_set = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db = motif_db,
  min_samples = 6,
  use_interactions = TRUE,  # Enable interaction terms
  use_future = TRUE,
  mc_cores = 30
)

ko_truth_FOXA2
ko_truth_HNF1A
ko_truth_IRF1
ko_truth_KLF5
ko_truth_RARG
ko_truth_SOX9

# usage for HNF1A and SOX9 (same pipeline for both TFs)

out_dir <- ko_dir  # or another directory

# Mixed LMM
box_hnf1a_lmm <- build_tf_boxplot_df_from_peaks(
  tf           = "HNF1A",
  pval_tbl     = pvals_hnf1a,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_HNF1A,
  model_type   = "lmm"
)
plot_tf_boxplots_fp_vs_rna("HNF1A_LMM", box_hnf1a_lmm, out_dir = ko_dir)

# Fixed (no interactions)
box_hnf1a_fixed <- build_tf_boxplot_df_from_peaks(
  tf           = "HNF1A",
  pval_tbl     = pvals_hnf1a_fixed,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_HNF1A,
  model_type   = "fixed"
)
plot_tf_boxplots_fp_vs_rna("HNF1A_fixed", box_hnf1a_fixed, out_dir = ko_dir)

# Fixed with interactions
box_hnf1a_fixed_int <- build_tf_boxplot_df_from_peaks(
  tf           = "HNF1A",
  pval_tbl     = pvals_hnf1a_fixed_int,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_HNF1A,
  model_type   = "fixed_int"
)



plot_tf_boxplots_fp_vs_rna("HNF1A_fixed_int", box_hnf1a_fixed_int, out_dir = ko_dir)



# Mixed LMM
box_sox9_lmm <- build_tf_boxplot_df_from_peaks(
  tf           = "SOX9",
  pval_tbl     = pvals_sox9,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_SOX9,
  model_type   = "lmm"
)
plot_tf_boxplots_fp_vs_rna("SOX9_LMM", box_sox9_lmm, out_dir = ko_dir)

# Fixed (no interactions)
box_sox9_fixed <- build_tf_boxplot_df_from_peaks(
  tf           = "SOX9",
  pval_tbl     = pvals_sox9_fixed,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_SOX9,
  model_type   = "fixed"
)
plot_tf_boxplots_fp_vs_rna("SOX9_fixed", box_sox9_fixed, out_dir = ko_dir)

# Fixed with interactions
box_sox9_fixed_int <- build_tf_boxplot_df_from_peaks(
  tf           = "SOX9",
  pval_tbl     = pvals_sox9_fixed_int,
  rna_tbl      = grn_set$rna,
  ko_truth_tbl = ko_truth_SOX9,
  model_type   = "fixed_int"
)
plot_tf_boxplots_fp_vs_rna("SOX9_fixed_int", box_sox9_fixed_int, out_dir = ko_dir)

#













# Minimal unified pipeline for: LMM → fixed → fixed_int → boxplots

future::plan(future::multisession, workers = 30)

tfs <- c("FOXA2", "IRF1", "KLF5", "RARG") # "HNF1A","SOX9"

run_all_for_tf <- function(tf) {
  message("=== Processing ", tf, " ===")

  ko_name <- paste0("ko_truth_", tf)
  if (!exists(ko_name)) {
    cli::cli_warn("Skipping {.val {tf}} — missing {.val {ko_name}}")
    return(invisible(NULL))
  }
  ko_tbl <- get(ko_name)

  # --- 1) Compute p-values ------------------------------------------------
  p_lmm <- compute_tf_pvals_lmm_enet(
    tf               = tf,
    grn_set          = grn_set,
    tf_gene_links_gh = tf_gene_links_gh,
    motif_db         = motif_db,
    min_samples      = 6L,
    use_future       = TRUE,
    mc_cores         = 30,
    cache_dir        = file.path(ko_dir, paste0("LMM_pvals_", tf, "_cache")),
    cache_chunk_size = 500L
  )

  p_fixed <- compute_tf_pvals_fixed_effects(
    tf               = tf,
    grn_set          = grn_set,
    tf_gene_links_gh = tf_gene_links_gh,
    motif_db         = motif_db,
    min_samples      = 6,
    use_interactions = FALSE,
    use_future       = TRUE,
    mc_cores         = 30
  )

  p_fixed_int <- compute_tf_pvals_fixed_effects(
    tf               = tf,
    grn_set          = grn_set,
    tf_gene_links_gh = tf_gene_links_gh,
    motif_db         = motif_db,
    min_samples      = 6,
    use_interactions = TRUE,
    use_future       = TRUE,
    mc_cores         = 30
  )

  # Save p-values to workspace
  assign(paste0("pvals_", tolower(tf)),            p_lmm,       envir = .GlobalEnv)
  assign(paste0("pvals_", tolower(tf), "_fixed"),  p_fixed,     envir = .GlobalEnv)
  assign(paste0("pvals_", tolower(tf), "_fixed_int"), p_fixed_int, envir = .GlobalEnv)

  # --- 2) Build boxplot tables ------------------------------------------
  box_lmm <-   build_tf_boxplot_df_from_peaks(tf, p_lmm,        grn_set$rna, ko_tbl, model_type = "lmm")
  box_fixed <- build_tf_boxplot_df_from_peaks(tf, p_fixed,      grn_set$rna, ko_tbl, model_type = "fixed")
  box_fixed_int <- build_tf_boxplot_df_from_peaks(tf, p_fixed_int, grn_set$rna, ko_tbl, model_type = "fixed_int")

  # --- 3) Plot PDFs -------------------------------------------------------
  plot_tf_boxplots_fp_vs_rna(paste0(tf, "_LMM"),       box_lmm,       out_dir = ko_dir)
  plot_tf_boxplots_fp_vs_rna(paste0(tf, "_fixed"),     box_fixed,     out_dir = ko_dir)
  plot_tf_boxplots_fp_vs_rna(paste0(tf, "_fixed_int"), box_fixed_int, out_dir = ko_dir)

  message("=== Done ", tf, " ===")
  invisible(NULL)
}

# Run for all TFs
lapply(tfs, run_all_for_tf)




# Process specific TFs first, then all remaining TFs

#' Extract all unique TFs from motif_db, handling dimers
#'
#' @param motif_db Tibble with columns motif, HGNC
#' @return Character vector of unique TF symbols
extract_all_tfs_from_motif_db <- function(motif_db) {
  all_tfs <- motif_db$HGNC %>%
    # Split dimers like "NR1H2::RXRA" into separate TFs
    strsplit("::") %>%
    unlist() %>%
    unique() %>%
    sort()

  # Remove NAs and empty strings
  all_tfs <- all_tfs[!is.na(all_tfs) & all_tfs != ""]

  cli::cli_inform("Found {length(all_tfs)} unique TFs in motif_db")
  all_tfs
}

#' Check if TF is expressed in RNA data
#'
#' @param tf Character scalar, TF symbol
#' @param rna_tbl Tibble with HGNC column
#' @return Logical
is_tf_expressed <- function(tf, rna_tbl) {
  tf %in% rna_tbl$HGNC
}

#' Process one TF: compute all p-values and save to CSV
#'
#' @param tf Character scalar, TF symbol
#' @param grn_set List with fp_score, rna, sample_metadata_used
#' @param tf_gene_links_gh GeneHancer links tibble
#' @param motif_db Motif database tibble
#' @param output_dir Directory to save CSV files
#' @param min_samples Integer, minimum samples for model fitting
#' @param mc_cores Integer, number of cores for parallel processing
#'
#' @return Invisible list with results, or NULL if TF skipped
process_one_tf_all_models <- function(tf,
                                      grn_set,
                                      tf_gene_links_gh,
                                      motif_db,
                                      output_dir,
                                      min_samples = 6L,
                                      mc_cores = 30L) {

  cli::cli_alert_info("Processing TF: {.val {tf}}")

  # Check if TF is expressed
  if (!is_tf_expressed(tf, grn_set$rna)) {
    cli::cli_alert_warning("TF {.val {tf}} not found in RNA data - skipping")
    return(invisible(NULL))
  }

  # --- 1) Compute LMM p-values ------------------------------------------------
  cli::cli_alert("Computing LMM p-values for {.val {tf}}")
  p_lmm <- tryCatch(
    compute_tf_pvals_lmm_enet(
      tf               = tf,
      grn_set          = grn_set,
      tf_gene_links_gh = tf_gene_links_gh,
      motif_db         = motif_db,
      min_samples      = min_samples,
      use_future       = TRUE,
      mc_cores         = mc_cores,
      cache_dir        = NULL,  # Don't cache intermediate results
      cache_chunk_size = 500L
    ),
    error = function(e) {
      cli::cli_alert_danger("LMM failed for {.val {tf}}: {e$message}")
      return(NULL)
    }
  )

  if (is.null(p_lmm) || !nrow(p_lmm)) {
    cli::cli_alert_warning("No LMM results for {.val {tf}} - skipping")
    return(invisible(NULL))
  }

  # --- 2) Compute Fixed effects p-values (no interactions) -------------------
  cli::cli_alert("Computing Fixed effects p-values for {.val {tf}}")
  p_fixed <- tryCatch(
    compute_tf_pvals_fixed_effects(
      tf               = tf,
      grn_set          = grn_set,
      tf_gene_links_gh = tf_gene_links_gh,
      motif_db         = motif_db,
      min_samples      = min_samples,
      use_interactions = FALSE,
      use_future       = TRUE,
      mc_cores         = mc_cores
    ),
    error = function(e) {
      cli::cli_alert_danger("Fixed effects failed for {.val {tf}}: {e$message}")
      return(NULL)
    }
  )

  # --- 3) Compute Fixed effects p-values (with interactions) -----------------
  cli::cli_alert("Computing Fixed effects (interactions) p-values for {.val {tf}}")
  p_fixed_int <- tryCatch(
    compute_tf_pvals_fixed_effects(
      tf               = tf,
      grn_set          = grn_set,
      tf_gene_links_gh = tf_gene_links_gh,
      motif_db         = motif_db,
      min_samples      = min_samples,
      use_interactions = TRUE,
      use_future       = TRUE,
      mc_cores         = mc_cores
    ),
    error = function(e) {
      cli::cli_alert_danger("Fixed effects (int) failed for {.val {tf}}: {e$message}")
      return(NULL)
    }
  )

  # --- 4) Combine all results into one tibble --------------------------------
  result <- p_lmm %>%
    dplyr::select(gene, fp_peak, p_corr,
                  p_lmm_fp, p_lmm_rna, p_lmm_both_fp, p_lmm_both_rna, p_lmm_both_overall)

  # Add fixed effects columns if available
  if (!is.null(p_fixed) && nrow(p_fixed)) {
    result <- result %>%
      dplyr::left_join(
        p_fixed %>%
          dplyr::select(gene, fp_peak,
                        p_fixed_fp, p_fixed_rna, p_fixed_both_fp,
                        p_fixed_both_rna, p_fixed_both_overall),
        by = c("gene", "fp_peak")
      )
  }

  # Add fixed effects with interactions columns if available
  if (!is.null(p_fixed_int) && nrow(p_fixed_int)) {
    result <- result %>%
      dplyr::left_join(
        p_fixed_int %>%
          dplyr::select(gene, fp_peak,
                        p_fixed_fp_int, p_fixed_rna_int,
                        p_fixed_both_fp_int, p_fixed_both_rna_int),
        by = c("gene", "fp_peak")
      )
  }

  # Add TF column at the beginning
  result <- result %>%
    dplyr::mutate(tf = tf, .before = 1)

  # --- 5) Save to CSV ---------------------------------------------------------
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  output_file <- file.path(output_dir, sprintf("TF_%s_all_pvals.csv", tf))
  readr::write_csv(result, output_file)

  cli::cli_alert_success("Saved results for {.val {tf}} to {.path {output_file}}")
  cli::cli_alert_info("  {nrow(result)} gene-peak pairs")

  invisible(list(
    tf = tf,
    results = result,
    output_file = output_file
  ))
}

#' Process TFs in batches: priority TFs first, then remaining
#'
#' @param grn_set List with fp_score, rna, sample_metadata_used
#' @param tf_gene_links_gh GeneHancer links tibble
#' @param motif_db Motif database tibble
#' @param output_dir Directory to save CSV files
#' @param priority_tfs Character vector of TFs to process first
#' @param min_samples Integer, minimum samples for model fitting
#' @param mc_cores Integer, number of cores for parallel processing
#' @param process_remaining Logical, whether to process remaining TFs after priority ones
#'
#' @return Tibble summarizing results for all TFs
process_tfs_priority_batch <- function(grn_set,
                                       tf_gene_links_gh,
                                       motif_db,
                                       output_dir,
                                       priority_tfs = c("FOXA2", "HNF1A", "IRF1", "KLF5", "RARG", "SOX9"),
                                       min_samples = 6L,
                                       mc_cores = 30L,
                                       process_remaining = TRUE) {

  # Extract all unique TFs from motif_db
  all_tfs <- extract_all_tfs_from_motif_db(motif_db)

  # Filter to only expressed TFs
  expressed_tfs <- all_tfs[sapply(all_tfs, is_tf_expressed, rna_tbl = grn_set$rna)]

  cli::cli_alert_info("{length(expressed_tfs)} of {length(all_tfs)} TFs are expressed in RNA data")

  # Split into priority and remaining
  priority_expressed <- intersect(priority_tfs, expressed_tfs)
  remaining_tfs <- setdiff(expressed_tfs, priority_tfs)

  cli::cli_rule("Priority TFs")
  cli::cli_alert_info("Processing {length(priority_expressed)} priority TFs first:")
  cli::cli_alert("{paste(priority_expressed, collapse = ', ')}")

  if (process_remaining) {
    cli::cli_alert_info("Then processing {length(remaining_tfs)} remaining TFs")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Combine TF lists: priority first, then remaining
  if (process_remaining) {
    tfs_to_process <- c(priority_expressed, remaining_tfs)
  } else {
    tfs_to_process <- priority_expressed
  }

  # Process each TF
  results_summary <- tibble::tibble()

  for (i in seq_along(tfs_to_process)) {
    tf <- tfs_to_process[i]

    # Mark if this is a priority TF
    is_priority <- tf %in% priority_tfs
    priority_label <- if (is_priority) " [PRIORITY]" else ""

    cli::cli_rule(sprintf("TF %d/%d: %s%s", i, length(tfs_to_process), tf, priority_label))

    result <- process_one_tf_all_models(
      tf               = tf,
      grn_set          = grn_set,
      tf_gene_links_gh = tf_gene_links_gh,
      motif_db         = motif_db,
      output_dir       = output_dir,
      min_samples      = min_samples,
      mc_cores         = mc_cores
    )

    # Record summary
    if (!is.null(result)) {
      results_summary <- dplyr::bind_rows(
        results_summary,
        tibble::tibble(
          tf = tf,
          priority = is_priority,
          n_gene_peak_pairs = nrow(result$results),
          n_genes = length(unique(result$results$gene)),
          n_peaks = length(unique(result$results$fp_peak)),
          output_file = result$output_file,
          status = "success"
        )
      )
    } else {
      results_summary <- dplyr::bind_rows(
        results_summary,
        tibble::tibble(
          tf = tf,
          priority = is_priority,
          n_gene_peak_pairs = 0L,
          n_genes = 0L,
          n_peaks = 0L,
          output_file = NA_character_,
          status = "failed"
        )
      )
    }
  }

  # Save summary
  summary_file <- file.path(output_dir, "processing_summary.csv")
  readr::write_csv(results_summary, summary_file)
  cli::cli_alert_success("Saved processing summary to {.path {summary_file}}")

  # Print summary statistics
  cli::cli_rule("Processing Summary")
  cli::cli_alert_info("Total TFs processed: {nrow(results_summary)}")

  if (any(results_summary$priority)) {
    priority_summary <- results_summary %>% dplyr::filter(priority)
    cli::cli_alert_info("Priority TFs: {nrow(priority_summary)}")
    cli::cli_alert_success("  Successful: {sum(priority_summary$status == 'success')}")
    cli::cli_alert_danger("  Failed: {sum(priority_summary$status == 'failed')}")
  }

  if (any(!results_summary$priority)) {
    remaining_summary <- results_summary %>% dplyr::filter(!priority)
    cli::cli_alert_info("Remaining TFs: {nrow(remaining_summary)}")
    cli::cli_alert_success("  Successful: {sum(remaining_summary$status == 'success')}")
    cli::cli_alert_danger("  Failed: {sum(remaining_summary$status == 'failed')}")
  }

  return(results_summary)
}

# USAGE

# Set up parallel processing
future::plan(future::multisession, workers = 30)

# Define output directory
output_dir <- file.path(ko_dir, "LMM_pvals_cache")

# Process priority TFs first: FOXA2, IRF1, KLF5, RARG, then HNF1A, SOX9, then rest
results_summary <- process_tfs_priority_batch(
  grn_set          = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db         = motif_db,
  output_dir       = output_dir,
  priority_tfs     = c("FOXA2", "IRF1", "KLF5", "RARG", "HNF1A", "SOX9"),
  min_samples      = 6L,
  mc_cores         = 30L,
  process_remaining = TRUE  # Set to FALSE to only process priority TFs
)

# View summary
print(results_summary)

# View priority TFs results
priority_results <- results_summary %>%
  dplyr::filter(priority) %>%
  dplyr::arrange(tf)
print(priority_results)

# Check which TFs failed
failed_tfs <- results_summary %>%
  dplyr::filter(status == "failed") %>%
  dplyr::pull(tf)

if (length(failed_tfs) > 0) {
  cli::cli_alert_warning("Failed TFs: {paste(failed_tfs, collapse = ', ')}")
}

# Alternative: Process ONLY the 4 specific TFs first

# If you want to run just the 4 TFs first, then manually run the rest later:
results_batch1 <- process_tfs_priority_batch(
  grn_set          = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db         = motif_db,
  output_dir       = output_dir,
  priority_tfs     = c("FOXA2", "IRF1", "KLF5", "RARG"),
  min_samples      = 6L,
  mc_cores         = 30L,
  process_remaining = FALSE  # Only process these 4
)

# Then run HNF1A and SOX9
results_batch2 <- process_tfs_priority_batch(
  grn_set          = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db         = motif_db,
  output_dir       = output_dir,
  priority_tfs     = c("HNF4A"), # "HNF1A", "SOX9",
  min_samples      = 6L,
  mc_cores         = 30L,
  process_remaining = FALSE  # Only process these 2
)

# Then run all remaining TFs
all_tfs <- extract_all_tfs_from_motif_db(motif_db)
expressed_tfs <- all_tfs[sapply(all_tfs, is_tf_expressed, rna_tbl = grn_set$rna)]
already_processed <- c("FOXA2", "IRF1", "KLF5", "RARG", "HNF1A", "SOX9")
remaining <- setdiff(expressed_tfs, already_processed)

results_batch3 <- process_tfs_priority_batch(
  grn_set          = grn_set,
  tf_gene_links_gh = tf_gene_links_gh,
  motif_db         = motif_db,
  output_dir       = output_dir,
  priority_tfs     = remaining,
  min_samples      = 6L,
  mc_cores         = 30L,
  process_remaining = FALSE
)




# Load saved TF pvals → collapse per gene → save results as lists

tfs_all <- c("IRF1", "RARG", "SOX9", "KLF5", "FOXA2", "HNF1A")

# Map TF → its KO truth object name
ko_truth_map <- list(
  IRF1  = ko_truth_IRF1,
  RARG  = ko_truth_RARG,
  SOX9  = ko_truth_SOX9,
  KLF5  = ko_truth_KLF5,
  FOXA2 = ko_truth_FOXA2,
  HNF1A = ko_truth_HNF1A
  # HNF4A = ko_truth_HNF4A
)

# Storage lists
tf_pvals_list      <- list()
# tf_pvals_collapsed <- list()

for (tf in tfs_all) {

  message("Loading TF: ", tf)

  # 1. Load saved CSV
  infile <- file.path(ko_dir, "LMM_pvals_cache", sprintf("TF_%s_all_pvals.csv", tf))
  pval_tbl <- readr::read_csv(infile, show_col_types = FALSE)

  tf_pvals_list[[tf]] <- pval_tbl

  # 2. Get KO truth table
  ko_tbl <- ko_truth_map[[tf]]
  if (is.null(ko_tbl)) {
    cli::cli_warn("KO truth table missing for TF {.val {tf}} — skipping collapse.")
    next
  }

  # # 3. Collapse per-gene/per-peak to one row per gene
  # collapsed_tbl <- collapse_tf_pvals_by_peak(
  #   pval_tbl     = pval_tbl,
  #   ko_truth_tbl = ko_tbl
  # )
  #
  # tf_pvals_collapsed[[tf]] <- collapsed_tbl
}

# Helper: collapse per-gene/per-peak for a specific model family

collapse_tf_pvals_by_peak_model <- function(pval_tbl,
                                            ko_truth_tbl,
                                            overall_col) {
  # Basic columns always required
  basic_required <- c("gene", "fp_peak", "p_corr")
  if (!all(basic_required %in% names(pval_tbl))) {
    cli::cli_abort("pval_tbl must contain columns {.val {basic_required}}.")
  }
  if (!all(c("gene", "ko_group") %in% names(ko_truth_tbl))) {
    cli::cli_abort("ko_truth_tbl must contain 'gene' and 'ko_group'.")
  }
  if (!overall_col %in% names(pval_tbl)) {
    cli::cli_abort("overall_col {.val {overall_col}} not found in pval_tbl.")
  }

  ko_tbl <- ko_truth_tbl |>
    dplyr::filter(
      !is.na(ko_group),
      ko_group %in% c("Down", "Unchanged", "Up")
    ) |>
    dplyr::select(gene, ko_group)

  merged <- pval_tbl |>
    dplyr::inner_join(ko_tbl, by = "gene")

  if (!nrow(merged)) {
    cli::cli_abort("No overlap between pval_tbl and ko_truth_tbl genes.")
  }

  # Collapse to one peak per gene × KO group using chosen overall_col
  collapsed <- merged |>
    dplyr::group_by(gene, ko_group) |>
    dplyr::group_modify(function(df, key) {
      ko_status <- key$ko_group

      # Special case: For "Unchanged" genes, just use the first peak
      if (ko_status == "Unchanged") {
        return(df[1L, , drop = FALSE])
      }

      # For "Down" and "Up" genes: select peak with best (minimum) overall p-value
      p_overall <- df[[overall_col]]

      if (all(is.na(p_overall))) {
        df[1L, , drop = FALSE]
      } else {
        best <- min(p_overall, na.rm = TRUE)
        df_sub <- df[p_overall == best, , drop = FALSE]
        df_sub[1L, , drop = FALSE]
      }
    }) |>
    dplyr::ungroup()

  collapsed
}

# Build model-specific collapsed lists:
#   - LMM (p_lmm_both_overall)
#   - Fixed (p_fixed_both_overall)
#   - Fixed+int (also keyed by p_fixed_both_overall)

tf_pvals_lmm_collapsed       <- list()
tf_pvals_fixed_collapsed     <- list()
tf_pvals_fixed_int_collapsed <- list()

for (tf in tfs_all) {
  message("Collapsing per model for TF: ", tf)

  pval_tbl <- tf_pvals_list[[tf]]
  ko_tbl   <- ko_truth_map[[tf]]

  if (is.null(pval_tbl) || is.null(ko_tbl)) {
    cli::cli_warn("Skipping {.val {tf}} — missing pval_tbl or ko_truth_tbl.")
    next
  }

  # 1) LMM family: use p_lmm_both_overall
  if ("p_lmm_both_overall" %in% names(pval_tbl)) {
    tf_pvals_lmm_collapsed[[tf]] <- collapse_tf_pvals_by_peak_model(
      pval_tbl     = pval_tbl,
      ko_truth_tbl = ko_tbl,
      overall_col  = "p_lmm_both_overall"
    )
  }

  # 2) Fixed-effects family (no interactions): use p_fixed_both_overall
  if ("p_fixed_both_overall" %in% names(pval_tbl)) {
    tf_pvals_fixed_collapsed[[tf]] <- collapse_tf_pvals_by_peak_model(
      pval_tbl     = pval_tbl,
      ko_truth_tbl = ko_tbl,
      overall_col  = "p_fixed_both_overall"
    )
  }

  # 3) Fixed-effects with interactions:
  #    still select peak based on p_fixed_both_overall (full vs null),
  #    but the chosen row carries the *_int p-values.
  if ("p_fixed_both_overall" %in% names(pval_tbl) &&
      any(grepl("^p_fixed_.*_int$", names(pval_tbl)))) {

    tf_pvals_fixed_int_collapsed[[tf]] <- collapse_tf_pvals_by_peak_model(
      pval_tbl     = pval_tbl,
      ko_truth_tbl = ko_tbl,
      overall_col  = "p_fixed_both_overall"
    )
  }
}

# Quick sanity check
str(tf_pvals_lmm_collapsed)
str(tf_pvals_fixed_collapsed,     max.level = 1)
str(tf_pvals_fixed_int_collapsed, max.level = 1)

# So the new method contains useful information and needs to be updated
