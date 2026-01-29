library(episcope)
progressr::handlers(global = TRUE)

in_dir <- "/data/homes/yl814/episcope_test/nutrient_stress/connect_tfs_to_target_genes"
tf_perturb_db <- "/data/homes/yl814/episcope/tf_perturb.db"
predicted_tfbs_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs"
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

tf_gene_links_gh   <- readr::read_csv(file.path(in_dir, "fp_gene_corr_full_jaspar2024.csv"))
# atac_gene_links_gh <- readr::read_csv(file.path(in_dir, "atac_gene_corr_full_jaspar2024.csv"))
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
tfs_all <- c("HNF1A", "HNF4A", "IRF1", "RARG", "SOX9", "KLF5", "FOXA2")

# ---- build_tf_gene_fp_expr_tbl() -------------------------------------------
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
    by           = c("peak_ID" = "fp_peak"),
    relationship = "many-to-many"
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
    by           = c("gene_key" = "gene", "sample_id" = "sample_id"),
    relationship = "many-to-many"
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

  # final tidy tibble: one row per gene × fp_peak × sample
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

# ---- utils_tf_multimodel: Multi-model TF → target gene framework ----------
#
# Simplified multi-model framework for TF → target gene prediction
# using TF RNA, TF footprint (TOBIAS), or both.
#
# Models implemented (per gene × fp_peak unit):
#   - LMM (random intercept)
#   - LMM with random slopes
#   - GAM (tp splines)
#   - GP-GAM (Gaussian process splines)
#   - GAM with interaction smooth ti(tf_expr, tf_fp)
#   - LMM+GAM hybrid (gamm4)
#   - SEM mediation (lavaan)
#
# Caching:
#   - For each TF, all per-unit model results are cached in a single
#     RDS file:  "<cache_dir>/TF_<TF>_all_models.rds"
#   - For each TF, a flattened summary tibble is written in CSV blocks:
#       "<cache_dir>/TF_<TF>_block_XXX.csv"

# ---- 0. Generic model summariser -------------------------------------------

# ---- model_summary() -------------------------------------------------------
#' Summarise a fitted model into coefficients + model-level stats
#'
#' @param model A fitted model object (lm, lmerMod, gam, etc.).
#'
#' @return A list with:
#'   - model_class: character scalar
#'   - coef: tibble of term-level estimates (if available)
#'   - glance: tibble of model-level statistics (if available)
#'   - fitted: numeric vector of fitted values (if available)
#'   - resid: numeric vector of residuals (if available)
model_summary <- function(model) {
  out <- list()
  out$model_class <- class(model)[1]

  # Coefficients
  if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")) {
    if (requireNamespace("broom.mixed", quietly = TRUE)) {
      out$coef <- broom.mixed::tidy(model, effects = "fixed")
    } else {
      out$coef <- tibble::as_tibble(summary(model)$coefficients, rownames = "term")
    }
  } else if (inherits(model, "gam")) {
    if (requireNamespace("broom", quietly = TRUE)) {
      out$coef <- broom::tidy(model)
    } else {
      out$coef <- tibble::as_tibble(summary(model)$p.table, rownames = "term")
    }
  } else if (inherits(model, "lm")) {
    if (requireNamespace("broom", quietly = TRUE)) {
      out$coef <- broom::tidy(model)
    } else {
      out$coef <- tibble::as_tibble(summary(model)$coefficients, rownames = "term")
    }
  } else {
    out$coef <- NULL
  }

  # Model-level statistics
  if (requireNamespace("broom", quietly = TRUE)) {
    out$glance <- tryCatch(
      broom::glance(model),
      error = function(e) NULL
    )
  } else {
    out$glance <- NULL
  }

  # Fitted and residuals (best-effort)
  out$fitted <- tryCatch(stats::fitted(model), error = function(e) NULL)
  out$resid  <- tryCatch(stats::residuals(model), error = function(e) NULL)

  out
}

# ---- 1. Per-unit model families (LMM, GAM, GP, SEM) ------------------------
#    All functions below assume df has columns:
#       expr    : target gene expression
#       tf_expr : TF RNA expression
#       tf_fp   : TF footprint score
#       cell    : factor
#       stress  : factor

# ---- 1a. LMM: random intercepts --------------------------------------------

# ---- tfs_to_target_genes_rna_lmm() -----------------------------------------
tfs_to_target_genes_rna_lmm <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_expr + (1 | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_fp_tobias_lmm() -----------------------------------
tfs_to_target_genes_fp_tobias_lmm <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_fp + (1 | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_rna_fp_tobias_lmm() -------------------------------
tfs_to_target_genes_rna_fp_tobias_lmm <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_expr + tf_fp + (1 | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- 1b. LMM: random slopes ------------------------------------------------

# ---- tfs_to_target_genes_rna_lmm_rand_slope() ------------------------------
tfs_to_target_genes_rna_lmm_rand_slope <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_expr + (tf_expr | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_fp_tobias_lmm_rand_slope() ------------------------
tfs_to_target_genes_fp_tobias_lmm_rand_slope <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_fp + (tf_fp | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_rna_fp_tobias_lmm_rand_slope() --------------------
tfs_to_target_genes_rna_fp_tobias_lmm_rand_slope <- function(df) {
  model <- lmerTest::lmer(
    expr ~ tf_expr + tf_fp + (tf_expr + tf_fp | cell) + (1 | stress),
    data = df,
    REML = FALSE
  )
  model_summary(model)
}

# ---- 1c. GAM models (thin-plate splines, mgcv) -----------------------------

# ---- tfs_to_target_genes_rna_gam() -----------------------------------------
tfs_to_target_genes_rna_gam <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_expr, k = 5) + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_fp_tobias_gam() -----------------------------------
tfs_to_target_genes_fp_tobias_gam <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_fp, k = 5) + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_rna_fp_tobias_gam() -------------------------------
tfs_to_target_genes_rna_fp_tobias_gam <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_expr, k = 5) + s(tf_fp, k = 5) + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- 1d. GP-GAM models (Gaussian process spline basis) --------------------

# ---- tfs_to_target_genes_rna_gam_gp() --------------------------------------
tfs_to_target_genes_rna_gam_gp <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_expr, k = 5, bs = "gp") + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_fp_tobias_gam_gp() --------------------------------
tfs_to_target_genes_fp_tobias_gam_gp <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_fp, k = 5, bs = "gp") + cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- tfs_to_target_genes_rna_fp_tobias_gam_gp() ----------------------------
tfs_to_target_genes_rna_fp_tobias_gam_gp <- function(df) {
  s <- mgcv::s
  model <- mgcv::gam(
    expr ~ s(tf_expr, k = 5, bs = "gp") +
      s(tf_fp, k = 5, bs = "gp") +
      cell + stress,
    data = df
  )
  model_summary(model)
}

# ---- 1e. GAM interaction smooth: RNA-FP synergy ---------------------------

# ---- tfs_to_target_genes_rna_fp_tobias_gam_ti() ----------------------------
tfs_to_target_genes_rna_fp_tobias_gam_ti <- function(df) {
  s  <- mgcv::s
  ti <- mgcv::ti

  model <- mgcv::gam(
    expr ~
      s(tf_expr, k = 5) +
      s(tf_fp,   k = 5) +
      ti(tf_expr, tf_fp, k = 5) +
      cell + stress,
    data = df
  )
  model_summary(model)
}


# ---- 1f. LMM + GAM hybrid (gamm4) ------------------------------------------

# ---- tfs_to_target_genes_rna_lmm_gam() -------------------------------------
tfs_to_target_genes_rna_lmm_gam <- function(df) {
  # Use local alias s <- mgcv::s so that the formula can use s().
  s <- mgcv::s

  fit <- gamm4::gamm4(
    expr ~ s(tf_expr, k = 5),
    random = ~(1 | cell) + (1 | stress),
    data   = df
  )
  model_summary(fit$gam)
}

# ---- tfs_to_target_genes_fp_tobias_lmm_gam() -------------------------------
tfs_to_target_genes_fp_tobias_lmm_gam <- function(df) {
  s <- mgcv::s

  fit <- gamm4::gamm4(
    expr ~ s(tf_fp, k = 5),
    random = ~(1 | cell) + (1 | stress),
    data   = df
  )
  model_summary(fit$gam)
}

# ---- tfs_to_target_genes_rna_fp_tobias_lmm_gam() ---------------------------
tfs_to_target_genes_rna_fp_tobias_lmm_gam <- function(df) {
  s <- mgcv::s

  fit <- gamm4::gamm4(
    expr ~ s(tf_expr, k = 5) + s(tf_fp, k = 5),
    random = ~(1 | cell) + (1 | stress),
    data   = df
  )
  model_summary(fit$gam)
}

# ---- 1g. SEM: mediation model (lavaan) -------------------------------------

# ---- tfs_to_target_genes_rna_fp_tobias_sem() -------------------------------
tfs_to_target_genes_rna_fp_tobias_sem <- function(df) {
  if (all(is.na(df$tf_expr)) || all(is.na(df$tf_fp))) {
    stop("SEM: tf_expr or tf_fp all NA for this unit.")
  }
  if (stats::sd(df$expr, na.rm = TRUE) == 0) {
    stop("SEM: expr has no variability for this unit.")
  }

  model_txt <- "
    # Structural paths
    tf_fp  ~ a * tf_expr
    expr   ~ b * tf_fp
    expr   ~ c * tf_expr

    # Indirect and total effects
    ab     := a * b
    total  := c + (a * b)
  "

  fit <- lavaan::sem(model_txt, data = df, fixed.x = FALSE)

  list(
    model_class  = "lavaan",
    parameters   = lavaan::parameterEstimates(fit),
    fit_measures = lavaan::fitMeasures(fit)
  )
}

# ---- 2. Build per-TF long-format data (gene × fp_peak × sample) ------------

# ---- prepare_tf_long_for_models() ------------------------------------------
#' Prepare TF-specific long-format table for modeling
#'
#' Uses build_tf_gene_fp_expr_tbl() to generate FP+RNA joined table,
#' then adds TF RNA expression (tf_expr). No aggregation across fp_peak;
#' each (gene, fp_peak) is treated as a separate unit.
prepare_tf_long_for_models <- function(tf,
                                       grn_set,
                                       tf_gene_links_gh,
                                       motif_db,
                                       min_peaks_per_gene = 1L,
                                       min_samples        = 6L) {

  tf_long <- build_tf_gene_fp_expr_tbl(
    tf                 = tf,
    grn_set            = grn_set,
    tf_gene_links_gh   = tf_gene_links_gh,
    motif_db           = motif_db,
    min_peaks_per_gene = min_peaks_per_gene
  )

  if (!nrow(tf_long)) {
    cli::cli_warn("No FP+RNA rows returned for TF {.val {tf}}.")
    return(tibble::tibble(
      gene      = character(0),
      fp_peak   = character(0),
      sample_id = character(0),
      expr      = numeric(0),
      tf_expr   = numeric(0),
      tf_fp     = numeric(0),
      cell      = character(0),
      stress    = character(0)
    ))
  }

  rna_tbl <- grn_set$rna
  if (!all(c("HGNC") %in% names(rna_tbl))) {
    cli::cli_abort("grn_set$rna must contain 'HGNC'.")
  }

  sample_cols_rna <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  tf_rows         <- which(rna_tbl$HGNC == tf)
  tf_expr_available <- length(tf_rows) > 0L && length(sample_cols_rna) > 0L

  if (!tf_expr_available) {
    cli::cli_warn(
      "TF {.val {tf}} not found (or no sample columns) in rna_tbl; setting tf_expr = NA."
    )
    tf_expr_df <- tibble::tibble(
      sample_id = sample_cols_rna,
      tf_expr   = NA_real_
    )
  } else {
    tf_mat      <- as.matrix(rna_tbl[tf_rows, sample_cols_rna, drop = FALSE])
    tf_expr_vec <- colMeans(tf_mat, na.rm = TRUE)
    tf_expr_df  <- tibble::tibble(
      sample_id = sample_cols_rna,
      tf_expr   = as.numeric(tf_expr_vec)
    )
  }

  tf_long2 <- dplyr::left_join(
    tf_long,
    tf_expr_df,
    by = "sample_id"
  )

  # Rename FP column to tf_fp; copy stress_type → stress
  tf_long2$tf_fp   <- tf_long2$fp_score
  tf_long2$stress  <- tf_long2$stress_type
  tf_long2$stress_type <- NULL

  # Filter units with at least min_samples samples
  tf_long2 <- tf_long2 |>
    dplyr::group_by(gene, fp_peak) |>
    dplyr::filter(dplyr::n() >= min_samples) |>
    dplyr::ungroup()

  tf_long2
}

# ---- 3. Per-unit model runner (GP + interaction GAMs + SEM) ----------------

# ---- run_all_models_for_gene() ---------------------------------------------
run_all_models_for_gene <- function(df_gene,
                                    min_samples = 6L,
                                    standardize = c("all", "predictors", "none")) {
  standardize <- match.arg(standardize)

  if (nrow(df_gene) < min_samples) {
    return(NULL)
  }

  # Basic sanity: need some variation in expr *before* any scaling
  if (stats::sd(df_gene$expr, na.rm = TRUE) == 0) {
    return(NULL)
  }

  # ------------------------------------------------------------------
  # Optional z-scaling (per unit) for stability
  # ------------------------------------------------------------------
  if (standardize != "none") {
    # Scale predictors (tf_expr, tf_fp)
    if (standardize %in% c("predictors", "all")) {
      if (!all(is.na(df_gene$tf_expr))) {
        sd_tf_expr <- stats::sd(df_gene$tf_expr, na.rm = TRUE)
        if (!is.na(sd_tf_expr) && sd_tf_expr > 0) {
          df_gene$tf_expr <- as.numeric(scale(df_gene$tf_expr))
        }
      }
      if (!all(is.na(df_gene$tf_fp))) {
        sd_tf_fp <- stats::sd(df_gene$tf_fp, na.rm = TRUE)
        if (!is.na(sd_tf_fp) && sd_tf_fp > 0) {
          df_gene$tf_fp <- as.numeric(scale(df_gene$tf_fp))
        }
      }
    }

    # Optionally scale outcome as well
    if (standardize == "all") {
      sd_expr <- stats::sd(df_gene$expr, na.rm = TRUE)
      if (!is.na(sd_expr) && sd_expr > 0) {
        df_gene$expr <- as.numeric(scale(df_gene$expr))
      }
    }
  }

  # Ensure factors (after scaling)
  df_gene$cell   <- factor(df_gene$cell)
  df_gene$stress <- factor(df_gene$stress)

  # Helper: safe execution that returns the error object on failure
  safe_run <- function(fun) {
    tryCatch(fun(df_gene), error = function(e) e)
  }

  list(
    # LMMs (random intercepts)
    tfs_to_target_genes_rna_lmm                   = safe_run(tfs_to_target_genes_rna_lmm),
    tfs_to_target_genes_fp_tobias_lmm             = safe_run(tfs_to_target_genes_fp_tobias_lmm),
    tfs_to_target_genes_rna_fp_tobias_lmm         = safe_run(tfs_to_target_genes_rna_fp_tobias_lmm),

    # LMMs (random slopes)
    tfs_to_target_genes_rna_lmm_rand_slope        = safe_run(tfs_to_target_genes_rna_lmm_rand_slope),
    tfs_to_target_genes_fp_tobias_lmm_rand_slope  = safe_run(tfs_to_target_genes_fp_tobias_lmm_rand_slope),
    tfs_to_target_genes_rna_fp_tobias_lmm_rand_slope =
      safe_run(tfs_to_target_genes_rna_fp_tobias_lmm_rand_slope),

    # GAMs (tp splines)
    tfs_to_target_genes_rna_gam                   = safe_run(tfs_to_target_genes_rna_gam),
    tfs_to_target_genes_fp_tobias_gam             = safe_run(tfs_to_target_genes_fp_tobias_gam),
    tfs_to_target_genes_rna_fp_tobias_gam         = safe_run(tfs_to_target_genes_rna_fp_tobias_gam),

    # GP-GAMs (Gaussian process)
    tfs_to_target_genes_rna_gam_gp                = safe_run(tfs_to_target_genes_rna_gam_gp),
    tfs_to_target_genes_fp_tobias_gam_gp          = safe_run(tfs_to_target_genes_fp_tobias_gam_gp),
    tfs_to_target_genes_rna_fp_tobias_gam_gp      = safe_run(tfs_to_target_genes_rna_fp_tobias_gam_gp),

    # Interaction smooth for RNA-FP synergy
    tfs_to_target_genes_rna_fp_tobias_gam_ti      = safe_run(tfs_to_target_genes_rna_fp_tobias_gam_ti),

    # LMM + GAM hybrids
    tfs_to_target_genes_rna_lmm_gam               = safe_run(tfs_to_target_genes_rna_lmm_gam),
    tfs_to_target_genes_fp_tobias_lmm_gam         = safe_run(tfs_to_target_genes_fp_tobias_lmm_gam),
    tfs_to_target_genes_rna_fp_tobias_lmm_gam     = safe_run(tfs_to_target_genes_rna_fp_tobias_lmm_gam),

    # SEM mediation
    tfs_to_target_genes_rna_fp_tobias_sem         = safe_run(tfs_to_target_genes_rna_fp_tobias_sem)
  )
}

# ---- 4. Helper: flatten per-gene results to tibble for CSV blocks ---------

# ---- flatten_gene_results() ------------------------------------------------
flatten_gene_results <- function(tf, gene, res_gene) {
  # res_gene is now a nested list: names = fp_peak, each element = list of models
  if (is.null(res_gene)) {
    return(NULL)
  }

  fp_peaks <- names(res_gene)
  if (is.null(fp_peaks)) {
    fp_peaks <- rep(NA_character_, length(res_gene))
  }

  rows <- list()
  idx  <- 0L

  for (j in seq_along(res_gene)) {
    fp_peak   <- fp_peaks[j]
    res_unit  <- res_gene[[j]]

    if (is.null(res_unit)) {
      next
    }

    model_names <- names(res_unit)

    for (mn in model_names) {
      mr <- res_unit[[mn]]

      if (inherits(mr, "error")) {
        idx <- idx + 1L
        rows[[idx]] <- tibble::tibble(
          tf          = tf,
          gene        = gene,
          fp_peak     = fp_peak,
          model       = mn,
          ok          = FALSE,
          model_class = NA_character_,
          error       = conditionMessage(mr)
        )
      } else if (is.list(mr) && !is.null(mr$model_class)) {
        idx <- idx + 1L
        rows[[idx]] <- tibble::tibble(
          tf          = tf,
          gene        = gene,
          fp_peak     = fp_peak,
          model       = mn,
          ok          = TRUE,
          model_class = mr$model_class,
          error       = NA_character_
        )
      } else {
        idx <- idx + 1L
        rows[[idx]] <- tibble::tibble(
          tf          = tf,
          gene        = gene,
          fp_peak     = fp_peak,
          model       = mn,
          ok          = TRUE,
          model_class = paste(class(mr), collapse = ";"),
          error       = NA_character_
        )
      }
    }
  }

  if (!length(rows)) {
    return(NULL)
  }

  dplyr::bind_rows(rows)
}

#
# 4a. Helpers for direction + p-value extraction (long-format)
#

# internal: map term to logical effect label
# ---- .infer_effect_from_term() ---------------------------------------------
.infer_effect_from_term <- function(term) {
  if (term %in% c("tf_expr", "s(tf_expr)")) {
    "tf_expr"
  } else if (term %in% c("tf_fp", "s(tf_fp)")) {
    "tf_fp"
  } else if (grepl("tf_expr", term) && grepl("tf_fp", term)) {
    "interaction"
  } else {
    NA_character_
  }
}

# internal: numeric sign
# ---- .sign_from_estimate() -------------------------------------------------
.sign_from_estimate <- function(est) {
  ifelse(
    is.na(est), NA_integer_,
    ifelse(est > 0, 1L, ifelse(est < 0, -1L, 0L))
  )
}

# internal: human-readable direction
# ---- .direction_from_sign() ------------------------------------------------
.direction_from_sign <- function(s) {
  ifelse(
    is.na(s), "unknown",
    ifelse(s > 0, "positive",
           ifelse(s < 0, "negative", "zero"))
  )
}
#
# NEW: full coefficient extraction (no filtering)
#

# ---- .full_coef_empty_tbl() ------------------------------------------------
.full_coef_empty_tbl <- function() {
  tibble::tibble(
    tf          = character(0),
    gene        = character(0),
    fp_peak     = character(0),
    model       = character(0),
    model_class = character(0),
    term        = character(0),
    predictor   = character(0),
    component   = character(0),
    estimate    = numeric(0),
    std_error   = numeric(0),
    statistic   = numeric(0),
    p_value     = numeric(0),
    sign        = integer(0),
    direction   = character(0),
    significant = logical(0)
  )
}

# ---- extract_model_terms_for_gene_full() -----------------------------------
#' Extract ALL terms for a single gene (no filtering)
#'
#' This keeps every row of the coefficient tables for all models
#' (and SEM parameters), and only *adds* annotations (predictor,
#' component, significant) without dropping anything.
#'
#' @param tf Character, TF symbol.
#' @param gene Character, target gene symbol.
#' @param res_gene Nested list of model results for this gene:
#'   names = fp_peak, each element = list(models).
#' @param alpha Numeric, significance threshold for p-values.
#'
#' @return Tibble with one row per coefficient/parameter:
#'   tf, gene, fp_peak, model, model_class, term, predictor, component,
#'   estimate, std_error, statistic, p_value, sign, direction, significant.
extract_model_terms_for_gene_full <- function(tf,
                                              gene,
                                              res_gene,
                                              alpha = 0.05) {
  if (is.null(res_gene) || !length(res_gene)) {
    return(.full_coef_empty_tbl())
  }

  out_list <- list()
  idx <- 0L

  fp_peaks <- names(res_gene)
  if (is.null(fp_peaks)) {
    fp_peaks <- rep(NA_character_, length(res_gene))
  }

  for (j in seq_along(res_gene)) {
    fp_peak  <- fp_peaks[j]
    models_j <- res_gene[[j]]

    if (is.null(models_j) || !length(models_j)) {
      next
    }

    for (mn in names(models_j)) {
      mr <- models_j[[mn]]
      if (inherits(mr, "error") || is.null(mr)) {
        next
      }

      model_class <- if (!is.null(mr$model_class)) mr$model_class else NA_character_

      ## --------------------------------------------------------------
      ## SEM (lavaan): use parameterEstimates
      ## --------------------------------------------------------------
      if (identical(model_class, "lavaan") && !is.null(mr$parameters)) {
        par_tbl <- mr$parameters
        if (!nrow(par_tbl)) next

        term <- paste(par_tbl$lhs, par_tbl$op, par_tbl$rhs)
        est  <- par_tbl$est
        se   <- if ("se" %in% names(par_tbl)) par_tbl$se else NA_real_
        stat <- if ("z" %in% names(par_tbl)) par_tbl$z else NA_real_
        pval <- if ("pvalue" %in% names(par_tbl)) par_tbl$pvalue else NA_real_

        pred <- vapply(
          term,
          .infer_effect_from_term,
          FUN.VALUE = character(1L)
        )

        comp <- "sem_parameter"
        sgn  <- .sign_from_estimate(est)

        idx <- idx + 1L
        out_list[[idx]] <- tibble::tibble(
          tf          = tf,
          gene        = gene,
          fp_peak     = fp_peak,
          model       = mn,
          model_class = model_class,
          term        = term,
          predictor   = pred,
          component   = comp,
          estimate    = est,
          std_error   = se,
          statistic   = stat,
          p_value     = pval,
          sign        = sgn,
          direction   = .direction_from_sign(sgn),
          significant = !is.na(pval) & pval < alpha
        )

        next
      }

      ## --------------------------------------------------------------
      ## Non-SEM models: coefficient tables
      ## --------------------------------------------------------------
      coef_tbl <- mr$coef
      if (is.null(coef_tbl) || !nrow(coef_tbl)) {
        next
      }

      # Normalise key columns
      if (!"term" %in% names(coef_tbl) && "Term" %in% names(coef_tbl)) {
        coef_tbl$term <- coef_tbl[["Term"]]
      }
      if (!"estimate" %in% names(coef_tbl) && "Estimate" %in% names(coef_tbl)) {
        coef_tbl$estimate <- coef_tbl[["Estimate"]]
      }
      if (!"p.value" %in% names(coef_tbl) && "Pr(>|t|)" %in% names(coef_tbl)) {
        coef_tbl$p.value <- coef_tbl[["Pr(>|t|)"]]
      }

      n <- nrow(coef_tbl)

      term <- as.character(coef_tbl$term)

      est <- if ("estimate" %in% names(coef_tbl)) {
        coef_tbl$estimate
      } else {
        rep(NA_real_, n)
      }

      se <- if ("std.error" %in% names(coef_tbl)) {
        coef_tbl$std.error
      } else if ("Std. Error" %in% names(coef_tbl)) {
        coef_tbl[["Std. Error"]]
      } else if ("Std.Error" %in% names(coef_tbl)) {
        coef_tbl[["Std.Error"]]
      } else {
        rep(NA_real_, n)
      }

      stat <- if ("statistic" %in% names(coef_tbl)) {
        coef_tbl$statistic
      } else if ("t.value" %in% names(coef_tbl)) {
        coef_tbl[["t.value"]]
      } else if ("z.value" %in% names(coef_tbl)) {
        coef_tbl[["z.value"]]
      } else {
        rep(NA_real_, n)
      }

      pval <- if ("p.value" %in% names(coef_tbl)) {
        coef_tbl$p.value
      } else {
        rep(NA_real_, n)
      }

      pred <- vapply(
        term,
        .infer_effect_from_term,
        FUN.VALUE = character(1L)
      )

      is_smooth <- grepl("^s\\(", term) | grepl("ti\\(", term)
      component <- ifelse(is_smooth, "smooth", "parametric")

      sgn <- .sign_from_estimate(est)

      idx <- idx + 1L
      out_list[[idx]] <- tibble::tibble(
        tf          = tf,
        gene        = gene,
        fp_peak     = fp_peak,
        model       = mn,
        model_class = model_class,
        term        = term,
        predictor   = pred,
        component   = component,
        estimate    = est,
        std_error   = se,
        statistic   = stat,
        p_value     = pval,
        sign        = sgn,
        direction   = .direction_from_sign(sgn),
        significant = !is.na(pval) & pval < alpha
      )
    }
  }

  if (!idx) {
    return(.full_coef_empty_tbl())
  }

  dplyr::bind_rows(out_list)
}



# ---- collect_model_terms_for_tf_full() -------------------------------------
#' Collect ALL model terms for one TF (no filtering) — memory-safe, data.table-backed
#'
#' This version is designed to avoid crashes on large inputs by *not* keeping the full
#' per-gene results list in memory. Instead, it processes genes in chunks and streams
#' each chunk to a temporary on-disk TSV (or to `out_file` if provided), then reads it
#' back once at the end.
#'
#' @param tf_result Single TF result list as returned by your
#'        non-streaming fit_all_models_for_tf() (must have $tf and $results).
#' @param alpha Numeric p-value threshold (used only to flag `significant`).
#' @param genes Optional character vector of gene symbols to subset.
#' @param mc_cores Optional integer > 1: parallel extraction over genes
#'        using parallel::mclapply (Unix/macOS only).
#' @param chunk_size Integer > 0. Number of genes processed per chunk. Default 250.
#' @param out_file Optional path. If provided, results are written there (TSV).
#'        If NULL, a temp file is used and removed on exit.
#'
#' @return Tibble with all coefficients / SEM parameters for this TF.
collect_model_terms_for_tf_full <- function(tf_result,
                                            alpha      = 0.05,
                                            genes      = NULL,
                                            mc_cores   = NULL,
                                            chunk_size = 250L,
                                            out_file   = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    cli::cli_abort("{.pkg data.table} is required for collect_model_terms_for_tf_full() backend.")
  }

  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) {
    cli::cli_abort("{.arg chunk_size} must be a positive integer.")
  }

  tf       <- tf_result$tf
  res_list <- tf_result$results

  if (!is.null(genes)) {
    keep <- base::intersect(base::names(res_list), genes)
    res_list <- res_list[keep]
  }

  gene_names <- base::names(res_list)
  if (is.null(gene_names) || !base::length(gene_names)) {
    return(.full_coef_empty_tbl())
  }

  worker <- function(g) {
    extract_model_terms_for_gene_full(
      tf       = tf,
      gene     = g,
      res_gene = res_list[[g]],
      alpha    = alpha
    )
  }

  use_parallel <- !is.null(mc_cores) &&
    is.numeric(mc_cores) &&
    as.integer(mc_cores) > 1L &&
    requireNamespace("parallel", quietly = TRUE)

  # Choose output path
  tmp <- FALSE
  if (is.null(out_file)) {
    out_file <- base::tempfile(pattern = "episcope_tf_terms_", fileext = ".tsv")
    tmp <- TRUE
  }
  if (tmp) {
    base::on.exit({
      if (base::file.exists(out_file)) base::unlink(out_file)
    }, add = TRUE)
  }

  wrote_any <- FALSE

  n <- base::length(gene_names)
  idx_starts <- seq.int(1L, n, by = chunk_size)

  for (s in idx_starts) {
    e <- base::min(s + chunk_size - 1L, n)
    genes_chunk <- gene_names[s:e]

    if (use_parallel) {
      mc_cores_i <- as.integer(mc_cores)
      chunk_list <- parallel::mclapply(
        X              = genes_chunk,
        FUN            = worker,
        mc.cores       = mc_cores_i,
        mc.preschedule = TRUE
      )
    } else {
      chunk_list <- base::lapply(genes_chunk, worker)
    }

    chunk_list <- base::Filter(function(x) !is.null(x) && base::nrow(x) > 0L, chunk_list)
    if (!base::length(chunk_list)) next

    chunk_dt <- data.table::rbindlist(
      l         = chunk_list,
      use.names = TRUE,
      fill      = TRUE
    )

    # Stream to disk (append) to keep RAM bounded
    if (!wrote_any) {
      data.table::fwrite(
        chunk_dt,
        file      = out_file,
        sep       = "\t",
        quote     = TRUE,
        na        = "NA",
        col.names = TRUE,
        append    = FALSE
      )
      wrote_any <- TRUE
    } else {
      data.table::fwrite(
        chunk_dt,
        file      = out_file,
        sep       = "\t",
        quote     = TRUE,
        na        = "NA",
        col.names = FALSE,
        append    = TRUE
      )
    }

    # Encourage memory release between chunks
    rm(chunk_dt, chunk_list)
    gc(FALSE)
  }

  if (!wrote_any) {
    return(.full_coef_empty_tbl())
  }

  out_dt <- data.table::fread(
    file      = out_file,
    sep       = "\t",
    na.strings = "NA",
    showProgress = FALSE
  )

  tibble::as_tibble(out_dt)
}

# ---- extract_model_effects_for_gene() --------------------------------------
#' Extract per-model effects for a single gene (long format)
#'
#' @param tf Character, TF symbol.
#' @param gene Character, target gene symbol.
#' @param res_gene Nested list of model results for this gene:
#'   names = fp_peak, each element = list(models).
#' @param alpha Numeric, significance threshold for p-values.
#'
#' @return Tibble with one row per (fp_peak, model, predictor/effect), columns:
#'   tf, gene, fp_peak, model, model_class, component, predictor, term,
#'   estimate, p_value, sign, direction, significant.
extract_model_effects_for_gene <- function(tf,
                                           gene,
                                           res_gene,
                                           alpha = 0.05) {
  if (is.null(res_gene) || !length(res_gene)) {
    return(tibble::tibble(
      tf          = character(0),
      gene        = character(0),
      fp_peak     = character(0),
      model       = character(0),
      model_class = character(0),
      component   = character(0),
      predictor   = character(0),
      term        = character(0),
      estimate    = numeric(0),
      p_value     = numeric(0),
      sign        = integer(0),
      direction   = character(0),
      significant = logical(0)
    ))
  }

  out_list <- list()
  idx <- 0L

  fp_peaks <- names(res_gene)
  if (is.null(fp_peaks)) {
    fp_peaks <- rep(NA_character_, length(res_gene))
  }

  for (j in seq_along(res_gene)) {
    fp_peak  <- fp_peaks[j]
    models_j <- res_gene[[j]]

    if (is.null(models_j) || !length(models_j)) {
      next
    }

    for (mn in names(models_j)) {
      mr <- models_j[[mn]]
      if (inherits(mr, "error") || is.null(mr)) {
        next
      }

      model_class <- if (!is.null(mr$model_class)) mr$model_class else NA_character_

      # SEM / lavaan handled separately
      if (identical(model_class, "lavaan") && !is.null(mr$parameters)) {
        par_tbl <- mr$parameters

        # direct tf_expr → expr
        row_c <- par_tbl[par_tbl$lhs == "expr" & par_tbl$rhs == "tf_expr" & par_tbl$op == "~", , drop = FALSE]
        if (nrow(row_c)) {
          est <- row_c$est[1]
          p   <- row_c$pvalue[1]
          sgn <- .sign_from_estimate(est)
          idx <- idx + 1L
          out_list[[idx]] <- tibble::tibble(
            tf          = tf,
            gene        = gene,
            fp_peak     = fp_peak,
            model       = mn,
            model_class = model_class,
            component   = "sem_direct",
            predictor   = "tf_expr",
            term        = "c (expr ~ tf_expr)",
            estimate    = est,
            p_value     = p,
            sign        = sgn,
            direction   = .direction_from_sign(sgn),
            significant = !is.na(p) && p < alpha
          )
        }

        # tf_fp → expr
        row_b <- par_tbl[par_tbl$lhs == "expr" & par_tbl$rhs == "tf_fp" & par_tbl$op == "~", , drop = FALSE]
        if (nrow(row_b)) {
          est <- row_b$est[1]
          p   <- row_b$pvalue[1]
          sgn <- .sign_from_estimate(est)
          idx <- idx + 1L
          out_list[[idx]] <- tibble::tibble(
            tf          = tf,
            gene        = gene,
            fp_peak     = fp_peak,
            model       = mn,
            model_class = model_class,
            component   = "sem_mediator",
            predictor   = "tf_fp",
            term        = "b (expr ~ tf_fp)",
            estimate    = est,
            p_value     = p,
            sign        = sgn,
            direction   = .direction_from_sign(sgn),
            significant = !is.na(p) && p < alpha
          )
        }

        # indirect (ab)
        row_ab <- par_tbl[par_tbl$label == "ab", , drop = FALSE]
        if (nrow(row_ab)) {
          est <- row_ab$est[1]
          p   <- row_ab$pvalue[1]
          sgn <- .sign_from_estimate(est)
          idx <- idx + 1L
          out_list[[idx]] <- tibble::tibble(
            tf          = tf,
            gene        = gene,
            fp_peak     = fp_peak,
            model       = mn,
            model_class = model_class,
            component   = "sem_indirect",
            predictor   = "tf_expr_indirect",
            term        = "ab (indirect)",
            estimate    = est,
            p_value     = p,
            sign        = sgn,
            direction   = .direction_from_sign(sgn),
            significant = !is.na(p) && p < alpha
          )
        }

        # total
        row_total <- par_tbl[par_tbl$label == "total", , drop = FALSE]
        if (nrow(row_total)) {
          est <- row_total$est[1]
          p   <- row_total$pvalue[1]
          sgn <- .sign_from_estimate(est)
          idx <- idx + 1L
          out_list[[idx]] <- tibble::tibble(
            tf          = tf,
            gene        = gene,
            fp_peak     = fp_peak,
            model       = mn,
            model_class = model_class,
            component   = "sem_total",
            predictor   = "tf_expr_total",
            term        = "total",
            estimate    = est,
            p_value     = p,
            sign        = sgn,
            direction   = .direction_from_sign(sgn),
            significant = !is.na(p) && p < alpha
          )
        }

        next
      }

      # Non-SEM models: use coef table
      coef_tbl <- mr$coef
      if (is.null(coef_tbl) || !nrow(coef_tbl)) {
        next
      }

      # Normalise column names
      if (!"term" %in% names(coef_tbl) && "Term" %in% names(coef_tbl)) {
        coef_tbl$term <- coef_tbl[["Term"]]
      }
      if (!"estimate" %in% names(coef_tbl) && "Estimate" %in% names(coef_tbl)) {
        coef_tbl$estimate <- coef_tbl[["Estimate"]]
      }
      if (!"p.value" %in% names(coef_tbl) && "Pr(>|t|)" %in% names(coef_tbl)) {
        coef_tbl$p.value <- coef_tbl[["Pr(>|t|)"]]
      }

      if (!"term" %in% names(coef_tbl) || !"p.value" %in% names(coef_tbl)) {
        next
      }

      # Keep only rows relevant to tf_expr/tf_fp/interaction
      for (i in seq_len(nrow(coef_tbl))) {
        term_i   <- as.character(coef_tbl$term[i])
        effect_i <- .infer_effect_from_term(term_i)
        if (is.na(effect_i)) {
          next
        }

        est <- if ("estimate" %in% names(coef_tbl)) coef_tbl$estimate[i] else NA_real_
        p   <- coef_tbl$p.value[i]

        # linear vs smooth/interaction components
        comp <- if (grepl("^s\\(", term_i) || grepl("ti\\(", term_i)) {
          if (effect_i == "interaction") "smooth_interaction" else "smooth"
        } else {
          "linear"
        }

        sgn <- if (comp == "linear") .sign_from_estimate(est) else NA_integer_

        idx <- idx + 1L
        out_list[[idx]] <- tibble::tibble(
          tf          = tf,
          gene        = gene,
          fp_peak     = fp_peak,
          model       = mn,
          model_class = model_class,
          component   = comp,
          predictor   = effect_i,
          term        = term_i,
          estimate    = est,
          p_value     = p,
          sign        = sgn,
          direction   = .direction_from_sign(sgn),
          significant = !is.na(p) && p < alpha
        )
      }
    }
  }

  if (!idx) {
    return(tibble::tibble(
      tf          = character(0),
      gene        = character(0),
      fp_peak     = character(0),
      model       = character(0),
      model_class = character(0),
      component   = character(0),
      predictor   = character(0),
      term        = character(0),
      estimate    = numeric(0),
      p_value     = numeric(0),
      sign        = integer(0),
      direction   = character(0),
      significant = logical(0)
    ))
  }

  dplyr::bind_rows(out_list)
}

# ---- collect_model_effects_for_tf() ----------------------------------------
#' Collect effects for all genes of one TF (wrapper)
#'
#' @param tf_result A single TF result list as returned by fit_all_models_for_tf().
#' @param alpha Numeric p-value threshold.
#'
#' @return Long-format tibble for all genes of this TF.
collect_model_effects_for_tf <- function(tf_result, alpha = 0.05) {
  tf <- tf_result$tf
  res_list <- tf_result$results

  out <- vector("list", length(res_list))
  idx <- 0L
  for (g in names(res_list)) {
    eff <- extract_model_effects_for_gene(
      tf       = tf,
      gene     = g,
      res_gene = res_list[[g]],
      alpha    = alpha
    )
    if (!nrow(eff)) next
    idx <- idx + 1L
    out[[idx]] <- eff
  }
  if (!idx) {
    return(tibble::tibble(
      tf          = character(0),
      gene        = character(0),
      fp_peak     = character(0),
      model       = character(0),
      model_class = character(0),
      component   = character(0),
      predictor   = character(0),
      term        = character(0),
      estimate    = numeric(0),
      p_value     = numeric(0),
      sign        = integer(0),
      direction   = character(0),
      significant = logical(0)
    ))
  }
  out <- out[seq_len(idx)]
  dplyr::bind_rows(out)
}

# internal: quick check if any effect is significant for a gene
# ---- .gene_has_significant_effect() ----------------------------------------
.gene_has_significant_effect <- function(tf,
                                         gene,
                                         res_gene,
                                         alpha = 0.05) {
  eff <- extract_model_effects_for_gene(
    tf       = tf,
    gene     = gene,
    res_gene = res_gene,
    alpha    = alpha
  )
  any(eff$significant, na.rm = TRUE)
}

#
# 4b. Helper: per-TF/gene debug plots (scatter + LM + GP-GAM)
#

# ---- plot_tf_gene_relationships() ------------------------------------------
plot_tf_gene_relationships <- function(tf,
                                       gene,
                                       df_gene,
                                       out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # alias for mgcv::s so it can be used in formulas
  s <- mgcv::s

  fp_peaks <- unique(df_gene$fp_peak)
  fp_peaks <- fp_peaks[order(fp_peaks)]

  if (!length(fp_peaks)) {
    cli::cli_inform("No fp_peak entries for TF {.val {tf}}, gene {.val {gene}}; skipping plots.")
    return(invisible(NULL))
  }

  ## ------------------------------------------------------------------
  ## 1) Global RNA-only plot (same for all peaks)
  ## ------------------------------------------------------------------
  d_rna <- df_gene[!is.na(df_gene$tf_expr), , drop = FALSE]

  plot_rna <- NULL
  if (nrow(d_rna) >= 3L && stats::sd(d_rna$tf_expr, na.rm = TRUE) > 0) {
    # Linear model: expr ~ tf_expr
    lm_fit_rna <- stats::lm(expr ~ tf_expr, data = d_rna)
    lm_sum_rna <- summary(lm_fit_rna)
    r2_lm_rna  <- lm_sum_rna$r.squared
    p_lm_rna   <- NA_real_
    if ("tf_expr" %in% rownames(stats::coef(lm_sum_rna))) {
      p_lm_rna <- stats::coef(lm_sum_rna)["tf_expr", "Pr(>|t|)"]
    }

    # GP-GAM: expr ~ s(tf_expr)
    gam_fit_rna <- tryCatch({
      mgcv::gam(expr ~ s(tf_expr, bs = "gp", k = 5), data = d_rna)
    }, error = function(e) NULL)

    if (!is.null(gam_fit_rna)) {
      gam_sum_rna <- summary(gam_fit_rna)
      if (!is.null(gam_sum_rna$s.table) && nrow(gam_sum_rna$s.table) >= 1L) {
        edf_gam_rna <- gam_sum_rna$s.table[1, "edf"]
        p_gam_rna   <- gam_sum_rna$s.table[1, "p-value"]
      } else {
        edf_gam_rna <- NA_real_
        p_gam_rna   <- NA_real_
      }
      dev_gam_rna <- if (!is.null(gam_sum_rna$dev.expl)) gam_sum_rna$dev.expl * 100 else NA_real_
    } else {
      edf_gam_rna <- NA_real_
      p_gam_rna   <- NA_real_
      dev_gam_rna <- NA_real_
    }

    label_rna <- if (!is.null(gam_fit_rna)) {
      sprintf(
        "LM: R^2 = %.2f, p = %.2g\nGP-GAM: edf = %.2f, dev = %.1f%%, p = %.2g",
        r2_lm_rna, p_lm_rna, edf_gam_rna, dev_gam_rna, p_gam_rna
      )
    } else {
      sprintf(
        "LM only: R^2 = %.2f, p = %.2g\n(GP-GAM fit failed)",
        r2_lm_rna, p_lm_rna
      )
    }

    plot_rna <- ggplot2::ggplot(d_rna, ggplot2::aes(x = tf_expr, y = expr)) +
      ggplot2::geom_point(alpha = 0.7, size = 1) +
      ggplot2::geom_smooth(
        method = "lm",
        se     = FALSE,
        colour = "red"
      ) +
      ggplot2::annotate(
        "text",
        x      = Inf,
        y      = Inf,
        label  = label_rna,
        hjust  = 1.02,
        vjust  = 1.02,
        size   = 3
      ) +
      ggplot2::labs(
        title = sprintf("%s \u2192 %s; RNA (tf_expr)", tf, gene),
        x     = "TF expression (tf_expr)",
        y     = "Target gene expression (expr)"
      ) +
      ggplot2::theme_minimal(base_size = 10)

    # Add GP-GAM smooth for RNA using the same data and mapping
    if (!is.null(gam_fit_rna)) {
      plot_rna <- plot_rna +
        ggplot2::geom_smooth(
          method  = "gam",
          formula = y ~ s(x, bs = "gp", k = 5),
          se      = FALSE,
          colour  = "blue"
        )
    }
  }

  ## ------------------------------------------------------------------
  ## 2) Peak-specific FP plots (tf_fp vs expr)
  ## ------------------------------------------------------------------
  pdf_file <- file.path(out_dir, sprintf("%s__%s_scatter.pdf", tf, gene))
  grDevices::pdf(pdf_file, width = 12, height = 4)
  any_page <- FALSE

  for (pk in fp_peaks) {
    d_unit <- df_gene[df_gene$fp_peak == pk, , drop = FALSE]
    d_fp   <- d_unit[!is.na(d_unit$tf_fp), , drop = FALSE]

    if (!nrow(d_fp) || stats::sd(d_fp$tf_fp, na.rm = TRUE) == 0) {
      next
    }

    # Linear model: expr ~ tf_fp
    lm_fit_fp <- stats::lm(expr ~ tf_fp, data = d_fp)
    lm_sum_fp <- summary(lm_fit_fp)
    r2_lm_fp  <- lm_sum_fp$r.squared
    p_lm_fp   <- NA_real_
    if ("tf_fp" %in% rownames(stats::coef(lm_sum_fp))) {
      p_lm_fp <- stats::coef(lm_sum_fp)["tf_fp", "Pr(>|t|)"]
    }

    # GP-GAM: expr ~ s(tf_fp)
    gam_fit_fp <- tryCatch({
      mgcv::gam(expr ~ s(tf_fp, bs = "gp", k = 5), data = d_fp)
    }, error = function(e) NULL)

    if (!is.null(gam_fit_fp)) {
      gam_sum_fp <- summary(gam_fit_fp)
      if (!is.null(gam_sum_fp$s.table) && nrow(gam_sum_fp$s.table) >= 1L) {
        edf_gam_fp <- gam_sum_fp$s.table[1, "edf"]
        p_gam_fp   <- gam_sum_fp$s.table[1, "p-value"]
      } else {
        edf_gam_fp <- NA_real_
        p_gam_fp   <- NA_real_
      }
      dev_gam_fp <- if (!is.null(gam_sum_fp$dev.expl)) gam_sum_fp$dev.expl * 100 else NA_real_
    } else {
      edf_gam_fp <- NA_real_
      p_gam_fp   <- NA_real_
      dev_gam_fp <- NA_real_
    }

    label_fp <- if (!is.null(gam_fit_fp)) {
      sprintf(
        "LM: R^2 = %.2f, p = %.2g\nGP-GAM: edf = %.2f, dev = %.1f%%, p = %.2g",
        r2_lm_fp, p_lm_fp, edf_gam_fp, dev_gam_fp, p_gam_fp
      )
    } else {
      sprintf(
        "LM only: R^2 = %.2f, p = %.2g\n(GP-GAM fit failed)",
        r2_lm_fp, p_lm_fp
      )
    }

    plot_fp <- ggplot2::ggplot(d_fp, ggplot2::aes(x = tf_fp, y = expr)) +
      ggplot2::geom_point(alpha = 0.7, size = 1) +
      ggplot2::geom_smooth(
        method = "lm",
        se     = FALSE,
        colour = "red"
      ) +
      ggplot2::annotate(
        "text",
        x      = Inf,
        y      = Inf,
        label  = label_fp,
        hjust  = 1.02,
        vjust  = 1.02,
        size   = 3
      ) +
      ggplot2::labs(
        title = sprintf("%s \u2192 %s (peak %s; FP)", tf, gene, pk),
        x     = "TF footprint (tf_fp)",
        y     = "Target gene expression (expr)"
      ) +
      ggplot2::theme_minimal(base_size = 10)

    # Add GP-GAM smooth for FP using the same data and mapping
    if (!is.null(gam_fit_fp)) {
      plot_fp <- plot_fp +
        ggplot2::geom_smooth(
          method  = "gam",
          formula = y ~ s(x, bs = "gp", k = 5),
          se      = FALSE,
          colour  = "blue"
        )
    }

    plots <- list()
    if (!is.null(plot_rna)) plots[[length(plots) + 1L]] <- plot_rna
    plots[[length(plots) + 1L]] <- plot_fp

    combined <- patchwork::wrap_plots(plots, nrow = 1)
    print(combined)
    any_page <- TRUE
  }

  grDevices::dev.off()

  if (!any_page) {
    cli::cli_inform(
      "No valid panels produced for TF {.val {tf}}, gene {.val {gene}}; PDF {.path {pdf_file}} has no plots."
    )
  } else {
    cli::cli_inform(
      "Wrote debug scatter PDF for TF {.val {tf}}, gene {.val {gene}} to {.path {pdf_file}}"
    )
  }

  invisible(pdf_file)
}



# ---- 5. Fit all models for one TF (caching + blocking + plots) -------------

# ---- fit_all_models_for_tf() -----------------------------------------------
fit_all_models_for_tf <- function(tf,
                                  grn_set,
                                  tf_gene_links_gh,
                                  motif_db,
                                  cache_dir,
                                  min_peaks_per_gene = 1L,
                                  min_samples        = 6L,
                                  overwrite          = FALSE,
                                  mc_cores           = NULL,
                                  rows_per_block     = 100L,  # here: ~genes per block
                                  max_genes          = NULL,
                                  make_plots         = FALSE,
                                  plot_dir           = NULL,
                                  max_plots_per_tf   = 30L,
                                  plot_only_sig      = FALSE,
                                  plot_alpha         = 0.05,
                                  standardize        = c("all", "predictors", "none")) {

  standardize <- match.arg(standardize)

  if (missing(cache_dir) || is.null(cache_dir)) {
    cli::cli_abort("cache_dir must be provided.")
  }
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ## -----------------------------------------------------------------
  ## Prepare TF-specific long-format data
  ## -----------------------------------------------------------------
  cli::cli_inform("Preparing long-format data for TF {.val {tf}}")

  tf_df <- prepare_tf_long_for_models(
    tf                 = tf,
    grn_set            = grn_set,
    tf_gene_links_gh   = tf_gene_links_gh,
    motif_db           = motif_db,
    min_peaks_per_gene = min_peaks_per_gene,
    min_samples        = min_samples
  )

  if (!nrow(tf_df)) {
    cli::cli_warn("No usable rows for TF {.val {tf}} after preparation; nothing to fit.")
    return(list(
      tf          = tf,
      block_files_csv  = character(0L),
      block_files_rds  = character(0L)
    ))
  }

  gene_splits <- split(tf_df, tf_df$gene)
  genes_all   <- names(gene_splits)

  if (!is.null(max_genes) && max_genes > 0L && length(genes_all) > max_genes) {
    genes_all   <- genes_all[seq_len(max_genes)]
    gene_splits <- gene_splits[genes_all]
    cli::cli_inform(
      "Debug mode: restricting to first {length(genes_all)} genes for TF {.val {tf}}."
    )
  }

  ## -----------------------------------------------------------------
  ## Determine already-processed genes from existing CSV blocks
  ## -----------------------------------------------------------------
  block_pattern_csv <- sprintf("^TF_%s_block_\\d+\\.csv$", tf)
  existing_blocks_csv <- list.files(
    cache_dir,
    pattern    = block_pattern_csv,
    full.names = TRUE
  )

  if (overwrite && length(existing_blocks_csv)) {
    cli::cli_inform(
      "overwrite = TRUE for TF {.val {tf}}: removing {length(existing_blocks_csv)} existing block CSV(s) in {.path {cache_dir}}."
    )
    unlink(existing_blocks_csv)
    existing_blocks_csv <- character(0L)
  }

  genes_done <- character(0L)
  if (!overwrite && length(existing_blocks_csv)) {
    cli::cli_inform(
      "Resuming TF {.val {tf}}: scanning {length(existing_blocks_csv)} existing block CSV(s) for completed genes."
    )

    for (bf in existing_blocks_csv) {
      gvec <- tryCatch({
        tmp <- readr::read_csv(
          bf,
          col_types      = readr::cols_only(gene = readr::col_character()),
          show_col_types = FALSE
        )
        unique(tmp$gene)
      }, error = function(e) character(0L))
      genes_done <- unique(c(genes_done, gvec))
    }
  }

  genes_todo <- setdiff(genes_all, genes_done)

  if (!length(genes_todo)) {
    cli::cli_inform(
      "All {length(genes_all)} genes for TF {.val {tf}} already processed; nothing to do."
    )
    ## Also list any existing RDS blocks
    block_pattern_rds <- sprintf("^TF_%s_block_\\d+_models\\.rds$", tf)
    existing_blocks_rds <- list.files(
      cache_dir,
      pattern    = block_pattern_rds,
      full.names = TRUE
    )
    return(list(
      tf              = tf,
      block_files_csv = sort(existing_blocks_csv),
      block_files_rds = sort(existing_blocks_rds)
    ))
  }

  cli::cli_inform(
    "Fitting models for {length(genes_todo)} remaining genes (of {length(genes_all)} total) for TF {.val {tf}}"
  )

  ## -----------------------------------------------------------------
  ## Parallelisation settings
  ## -----------------------------------------------------------------
  use_parallel <-
    !is.null(mc_cores) &&
    is.numeric(mc_cores) &&
    as.integer(mc_cores) > 1L &&
    requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    mc_cores <- as.integer(mc_cores)
    cli::cli_inform(
      "Using parallel::mclapply() over {length(genes_todo)} genes on {mc_cores} cores for TF {.val {tf}}."
    )
  } else {
    cli::cli_inform(
      "Using sequential processing over {length(genes_todo)} genes for TF {.val {tf}}."
    )
  }

  ## Decide which genes to consider for plotting (subset of all genes)
  genes_for_plots <- character(0L)
  if (make_plots && max_plots_per_tf > 0L) {
    n_plot <- min(max_plots_per_tf, length(genes_all))
    genes_for_plots <- genes_all[seq_len(n_plot)]
  }

  if (is.null(plot_dir)) {
    plot_dir <- file.path(cache_dir, "plots", tf)
  } else {
    plot_dir <- file.path(plot_dir, tf)
  }
  if (make_plots && !dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ## -----------------------------------------------------------------
  ## Streaming over blocks of genes
  ## -----------------------------------------------------------------

  if (is.null(rows_per_block) || rows_per_block <= 0L) {
    rows_per_block <- 100L
  }
  genes_per_block <- as.integer(rows_per_block)

  block_files_csv <- sort(existing_blocks_csv)
  block_files_rds <- character(0L)
  block_index_offset <- length(existing_blocks_csv)

  ## Helper: fit all models for one gene, optionally plot, and return
  ## both the nested model list and the flattened status tibble.
  fit_and_flatten_gene <- function(g) {
    dfg      <- gene_splits[[g]]
    fp_peaks <- unique(dfg$fp_peak)

    res_gene <- vector("list", length(fp_peaks))
    names(res_gene) <- fp_peaks

    for (j in seq_along(fp_peaks)) {
      pk     <- fp_peaks[j]
      dfg_pk <- dfg[dfg$fp_peak == pk, , drop = FALSE]

      res_unit <- run_all_models_for_gene(
        df_gene     = dfg_pk,
        min_samples = min_samples,
        standardize = standardize
      )
      res_gene[[j]] <- res_unit
    }

    ## Optional plotting
    if (make_plots &&
        g %in% genes_for_plots) {

      if (!plot_only_sig ||
          .gene_has_significant_effect(tf, g, res_gene, alpha = plot_alpha)) {
        plot_tf_gene_relationships(
          tf      = tf,
          gene    = g,
          df_gene = dfg,
          out_dir = plot_dir
        )
      }
    }

    flat <- flatten_gene_results(tf = tf, gene = g, res_gene = res_gene)

    list(
      gene     = g,
      res_gene = res_gene,
      flat     = flat
    )
  }

  ## Iterate over gene blocks
  block_counter <- 0L

  for (start_idx in seq(1L, length(genes_todo), by = genes_per_block)) {
    end_idx <- min(start_idx + genes_per_block - 1L, length(genes_todo))
    genes_block <- genes_todo[start_idx:end_idx]

    cli::cli_inform(
      "TF {.val {tf}}: processing genes {start_idx}-{end_idx} of {length(genes_todo)}."
    )

    if (use_parallel) {
      block_list <- parallel::mclapply(
        X            = genes_block,
        FUN          = fit_and_flatten_gene,
        mc.cores     = mc_cores,
        mc.preschedule = TRUE
      )
    } else {
      block_list <- vector("list", length(genes_block))
      for (i in seq_along(genes_block)) {
        block_list[[i]] <- fit_and_flatten_gene(genes_block[i])
      }
    }

    ## Keep only entries with a non-empty 'flat' tibble
    block_list <- Filter(
      function(x) !is.null(x$flat) && nrow(x$flat) > 0L,
      block_list
    )

    if (!length(block_list)) {
      cli::cli_inform(
        "TF {.val {tf}}: no non-empty model summaries in this block; skipping CSV/RDS write."
      )
      next
    }

    ## Combine flattened status rows
    block_df <- dplyr::bind_rows(
      lapply(block_list, function(x) x$flat)
    )

    ## Build nested results list for this block: named by gene
    gene_names_block <- vapply(block_list, function(x) x$gene, character(1L))
    res_block <- setNames(
      lapply(block_list, function(x) x$res_gene),
      gene_names_block
    )

    ## ---- Write CSV block (status) ----
    block_counter <- block_counter + 1L
    block_id      <- block_index_offset + block_counter

    block_file_csv <- file.path(
      cache_dir,
      sprintf("TF_%s_block_%03d.csv", tf, block_id)
    )
    readr::write_csv(block_df, block_file_csv)
    block_files_csv <- c(block_files_csv, block_file_csv)

    ## ---- Write RDS block (full models) ----
    block_file_rds <- file.path(
      cache_dir,
      sprintf("TF_%s_block_%03d_models.rds", tf, block_id)
    )
    saveRDS(
      list(
        tf      = tf,
        results = res_block
      ),
      file = block_file_rds
    )
    block_files_rds <- c(block_files_rds, block_file_rds)

    cli::cli_inform(
      "TF {.val {tf}}: wrote block {block_id} with {length(res_block)} gene(s) to {.path {block_file_csv}} and {.path {block_file_rds}}."
    )
  }

  cli::cli_inform(
    "Finished TF {.val {tf}}: total block CSVs = {length(block_files_csv)}, total block RDS = {length(block_files_rds)}."
  )

  list(
    tf              = tf,
    block_files_csv = sort(unique(block_files_csv)),
    block_files_rds = sort(unique(block_files_rds))
  )
}



# ---- load_tf_models_from_blocks() ------------------------------------------
#' Load all model results for one TF from block RDS files
#'
#' @param cache_dir Directory used for TF_multimodel_full_cache.
#' @param tf        TF symbol, e.g. "HNF1A".
#'
#' @return A list with elements:
#'   - tf: character scalar
#'   - results: named list by gene (each is the nested fp_peak → models list)
load_tf_models_from_blocks <- function(cache_dir, tf) {
  pattern_rds <- sprintf("^TF_%s_block_\\d+_models\\.rds$", tf)
  files_rds <- list.files(
    cache_dir,
    pattern    = pattern_rds,
    full.names = TRUE
  )

  if (!length(files_rds)) {
    cli::cli_abort(
      "No model RDS files for TF {.val {tf}} found in {.path {cache_dir}}."
    )
  }

  results_all <- list()

  for (f in files_rds) {
    blk <- readRDS(f)
    if (!identical(blk$tf, tf)) {
      next
    }
    # blk$results is a named list by gene; genes are unique per block
    for (g in names(blk$results)) {
      if (g %in% names(results_all)) {
        cli::cli_warn(
          "Duplicate gene {.val {g}} for TF {.val {tf}} in file {.path {f}}; overwriting previous entry."
        )
      }
      results_all[[g]] <- blk$results[[g]]
    }
  }

  list(
    tf      = tf,
    results = results_all
  )
}

# ---- 6. Run for many TFs (sequential, with caching and debug) --------------

# ---- run_all_tfs_all_models() ----------------------------------------------
run_all_tfs_all_models <- function(tfs,
                                   grn_set,
                                   tf_gene_links_gh,
                                   motif_db,
                                   cache_dir,
                                   min_peaks_per_gene = 1L,
                                   min_samples        = 6L,
                                   overwrite          = FALSE,
                                   mc_cores           = NULL,
                                   rows_per_block     = 100L,
                                   max_genes_per_tf   = NULL,
                                   make_plots         = FALSE,
                                   plot_dir           = NULL,
                                   max_plots_per_tf   = 30L,
                                   plot_only_sig      = FALSE,
                                   plot_alpha         = 0.05,
                                   standardize        = c("all", "predictors", "none")) {

  standardize <- match.arg(standardize)

  res_list <- vector("list", length(tfs))
  names(res_list) <- tfs

  for (i in seq_along(tfs)) {
    tf <- tfs[i]
    cli::cli_rule(paste0("TF ", i, "/", length(tfs), ": ", tf))

    res_list[[tf]] <- fit_all_models_for_tf(
      tf                 = tf,
      grn_set            = grn_set,
      tf_gene_links_gh   = tf_gene_links_gh,
      motif_db           = motif_db,
      cache_dir          = cache_dir,
      min_peaks_per_gene = min_peaks_per_gene,
      min_samples        = min_samples,
      overwrite          = overwrite,
      mc_cores           = mc_cores,
      rows_per_block     = rows_per_block,
      max_genes          = max_genes_per_tf,
      make_plots         = make_plots,
      plot_dir           = plot_dir,
      max_plots_per_tf   = max_plots_per_tf,
      plot_only_sig      = plot_only_sig,
      plot_alpha         = plot_alpha,
      standardize        = standardize
    )
  }

  res_list
}


# Usage

# small test run to check
# tf <- "HNF1A"
# tf_df <- prepare_tf_long_for_models(
#   tf                 = tf,
#   grn_set            = grn_set,
#   tf_gene_links_gh   = tf_gene_links_gh,
#   motif_db           = motif_db,
#   min_peaks_per_gene = 1L,
#   min_samples        = 6L
# )

# some_gene <- unique(tf_df$gene)[1]
# df_gene   <- tf_df[tf_df$gene == some_gene, , drop = FALSE]

# plot_tf_gene_relationships(
#   tf      = tf,
#   gene    = some_gene,
#   df_gene = df_gene,
#   out_dir = file.path(cache_dir_models, "debug_plots", tf)
# )

# cache_dir_models <- file.path(ko_dir, "TF_multimodel_cache")
# tfs_interest <- c("HNF1A", "HNF4A", "IRF1", "RARG", "SOX9", "KLF5", "FOXA2")

# Debug run with full model set and optional plots (parallel OK)
# debug_results <- run_all_tfs_all_models(
#   tfs                = tfs_interest,
#   grn_set            = grn_set,
#   tf_gene_links_gh   = tf_gene_links_gh,
#   motif_db           = motif_db,
#   cache_dir          = cache_dir_models,
#   min_peaks_per_gene = 1L,
#   min_samples        = 6L,
#   overwrite          = TRUE,
#   mc_cores           = 36L,
#   rows_per_block     = 1000L,
#   # max_genes_per_tf   = 500L,
#   make_plots         = FALSE,
#   plot_dir           = file.path(cache_dir_models, "debug_plots"),
#   # max_plots_per_tf   = 500L,
#   plot_only_sig      = TRUE,
#   plot_alpha         = 0.05,
#   standardize        = "all"
# )

# extract long-format effects for one TF
# effects_HNF1A <- collect_model_effects_for_tf(debug_results[["HNF1A"]], alpha = 0.05)



# Optional: debug run with per TF-gene plots (two panels per PDF, LM + GP-GAM)
# debug_results_plots <- run_all_tfs_all_models(
#   tfs                = c("HNF1A"),
#   grn_set            = grn_set,
#   tf_gene_links_gh   = tf_gene_links_gh,
#   motif_db           = motif_db,
#   cache_dir          = cache_dir_models,
#   min_peaks_per_gene = 1L,
#   min_samples        = 6L,
#   overwrite          = TRUE,
#   mc_cores           = 36L,  # plotting is supported in parallel
#   rows_per_block     = 1000L,
#   max_genes_per_tf   = 50L,
#   make_plots         = TRUE,
#   plot_dir           = file.path(cache_dir_models, "debug_plots"),
#   max_plots_per_tf   = 30L,
#   standardize        = "all"
# )


cache_dir_models <- file.path(ko_dir, "TF_multimodel_full_cache")
tfs_interest <- c("KLF5", "FOXA2", "HNF4A") # "HNF1A", "IRF1", "RARG", "SOX9",
# results_full <- run_all_tfs_all_models(
#   tfs                = tfs_interest,
#   grn_set            = grn_set, # strict_nutrient_stress grn
#   tf_gene_links_gh   = tf_gene_links_gh,
#   motif_db           = motif_db,
#   cache_dir          = cache_dir_models,
#   min_peaks_per_gene = 1L,
#   min_samples        = 6L,
#   overwrite          = FALSE,       # <-- use FALSE to allow resume
#   mc_cores           = 30L,
#   rows_per_block     = 100L,        # ≈ genes per CSV block
#   max_genes_per_tf   = NULL,
#   make_plots         = FALSE,       # no plots here
#   plot_dir           = file.path(cache_dir_models, "debug_plots"),
#   max_plots_per_tf   = 0L,
#   plot_only_sig      = TRUE,
#   plot_alpha         = 0.05,
#   standardize        = "all"
# )









# tf_result_HNF1A <- load_tf_models_from_blocks(cache_dir_models, "HNF1A")
# effects_HNF1A <- collect_model_effects_for_tf(tf_result_HNF1A, alpha = 0.05)

# library(dplyr)
# blocks_HNF1A <- list.files(
#   cache_dir_models,
#   pattern    = "^TF_HNF1A_block_\\d+\\.csv$",
#   full.names = TRUE
# )

# effects_HNF1A <- dplyr::bind_rows(
#   lapply(blocks_HNF1A, readr::read_csv)
# )
# terms_HNF1A <- collect_model_terms_for_tf_full(tf_result_HNF1A, 0.05, mc_cores = 36L)
# terms_HNF1A


# utils_ko_peak_summaries.R
# Simple per-peak TF-gene summary tables for KO benchmarking
# Author: Yaoxiang Li (episcope)

# Small helpers

# ---- .sign_and_direction() -------------------------------------------------
#' Internal: convert numeric estimate -> sign (+1/0/-1) and direction string
#'
#' @param x Numeric estimate.
#' @return List with integer `sign` and character `direction`.
.sign_and_direction <- function(x) {
  if (!is.finite(x) || x == 0) {
    return(list(sign = 0L, direction = NA_character_))
  }
  s <- if (x > 0) 1L else -1L
  d <- if (s > 0) "positive" else "negative"
  list(sign = s, direction = d)
}

# TF-gene RNA correlations (for a given TF)
# ---- compute_tf_rna_corr_for_genes() ---------------------------------------
#' Compute TF-gene RNA correlations for a single TF
#'
#' Correlates TF RNA expression with all genes in `rna_tbl` using Pearson
#' (or Spearman) correlation across samples.
#'
#' @param tf Character scalar, TF symbol (e.g. "HNF1A").
#' @param rna_tbl Tibble/data.frame with columns:
#'   - `HGNC`: gene symbol
#'   - sample columns: one column per sample (numeric)
#' @param genes Optional character vector of gene symbols to keep.
#' @param method Correlation method passed to [stats::cor()].
#'
#' @return Tibble with columns:
#'   - gene  : HGNC symbol
#'   - r_rna : correlation coefficient
#'   - p_rna : two-sided correlation p-value
compute_tf_rna_corr_for_genes <- function(tf,
                                          rna_tbl,
                                          genes  = NULL,
                                          method = c("pearson", "spearman")) {
  method <- match.arg(method)

  if (!"HGNC" %in% names(rna_tbl)) {
    cli::cli_abort("rna_tbl must contain column 'HGNC'.")
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
      r_rna = numeric(0),
      p_rna = numeric(0)
    ))
  }

  tf_expr <- colMeans(rna_mat[tf_rows, , drop = FALSE], na.rm = TRUE)
  tf_expr <- suppressWarnings(as.numeric(tf_expr))

  if (sum(is.finite(tf_expr)) < 3L || stats::sd(tf_expr, na.rm = TRUE) == 0) {
    cli::cli_warn("Insufficient variation in TF expression for {.val {tf}}.")
    return(tibble::tibble(
      gene  = character(0),
      r_rna = numeric(0),
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

  out <- tibble::tibble(
    gene  = rna_tbl$HGNC,
    r_rna = as.numeric(cor_vec),
    p_rna = as.numeric(p_vec)
  )

  if (!is.null(genes)) {
    out <- out[out$gene %in% genes, , drop = FALSE]
  }

  out[!is.na(out$gene) & out$gene != "", , drop = FALSE]
}

# Summaries of model terms per (TF, gene, fp_peak)
# ---- summarise_linear_effect_per_peak() ------------------------------------
#' Internal: summarise best linear effect (tf_expr or tf_fp) per peak
#'
#' This helper:
#'   - keeps only `component == "parametric"` rows,
#'   - keeps only specified `predictor` (e.g. "tf_expr" or "tf_fp"),
#'   - restricts to linear models (`model_class` in c("lmerModLmerTest","lm")),
#'   - within each (tf, gene, fp_peak) chooses the row with smallest p-value.
#'
#' @param terms_tf Tibble from `collect_model_terms_for_tf_full()`.
#' @param predictor Character scalar: "tf_expr" or "tf_fp".
#' @param prefix Prefix for output columns (e.g. "lin_tf_expr").
#'
#' @return Tibble with columns:
#'   tf, gene, fp_peak,
#'   `<prefix>_estimate`, `<prefix>_p_value`,
#'   `<prefix>_sign`, `<prefix>_direction`.
summarise_linear_effect_per_peak <- function(terms_tf,
                                             predictor,
                                             prefix) {
  keep <- terms_tf$component == "parametric" &
    terms_tf$predictor == predictor &
    terms_tf$model_class %in% c("lmerModLmerTest", "lm") &
    is.finite(terms_tf$p_value)

  sub <- terms_tf[keep, , drop = FALSE]

  if (!nrow(sub)) {
    return(tibble::tibble(
      tf                              = character(0),
      gene                            = character(0),
      fp_peak                         = character(0),
      !!paste0(prefix, "_estimate")   := numeric(0),
      !!paste0(prefix, "_p_value")    := numeric(0),
      !!paste0(prefix, "_sign")       := integer(0),
      !!paste0(prefix, "_direction")  := character(0)
    ))
  }

  out <- sub |>
    dplyr::group_by(tf, gene, fp_peak) |>
    dplyr::summarise(
      .groups = "drop",
      estimate = {
        idx <- which.min(p_value)
        estimate[idx]
      },
      p_value = min(p_value, na.rm = TRUE)
    )

  # Add sign and direction
  sign_vec  <- integer(nrow(out))
  dir_vec   <- character(nrow(out))

  for (i in seq_len(nrow(out))) {
    sd <- .sign_and_direction(out$estimate[i])
    sign_vec[i] <- sd$sign
    dir_vec[i]  <- sd$direction
  }

  out[[paste0(prefix, "_estimate")]]  <- out$estimate
  out[[paste0(prefix, "_p_value")]]   <- out$p_value
  out[[paste0(prefix, "_sign")]]      <- sign_vec
  out[[paste0(prefix, "_direction")]] <- dir_vec

  out <- out[, c("tf", "gene", "fp_peak",
                 paste0(prefix, "_estimate"),
                 paste0(prefix, "_p_value"),
                 paste0(prefix, "_sign"),
                 paste0(prefix, "_direction")),
             drop = FALSE]

  out
}

# ---- summarise_sem_effects_per_peak() --------------------------------------
#' Internal: summarise SEM effects (lavaan) per peak
#'
#' Extracts selected lavaan parameters:
#'   - direct TF→gene path via expression (term "expr ~ tf_expr")
#'   - direct TF→gene path via footprint (term "expr ~ tf_fp")
#'   - indirect (mediated) effect ("ab :=")
#'   - total effect ("total :=")
#'
#' For each (tf, gene, fp_peak, term-pattern) we pick the row with smallest
#' p-value (there is usually only one).
#'
#' @param terms_tf Tibble from `collect_model_terms_for_tf_full()`.
#'
#' @return Tibble with columns:
#'   tf, gene, fp_peak,
#'   sem_direct_tf_expr_est,  sem_direct_tf_expr_p,
#'   sem_direct_tf_expr_sign, sem_direct_tf_expr_direction,
#'   sem_direct_tf_fp_est,    sem_direct_tf_fp_p,
#'   sem_direct_tf_fp_sign,   sem_direct_tf_fp_direction,
#'   sem_indirect_ab_est,     sem_indirect_ab_p,
#'   sem_indirect_ab_sign,    sem_indirect_ab_direction,
#'   sem_total_est,           sem_total_p,
#'   sem_total_sign,          sem_total_direction.
summarise_sem_effects_per_peak <- function(terms_tf) {
  is_sem <- terms_tf$model_class == "lavaan" & is.finite(terms_tf$p_value)
  sem_tbl <- terms_tf[is_sem, , drop = FALSE]

  if (!nrow(sem_tbl)) {
    return(tibble::tibble(
      tf                             = character(0),
      gene                           = character(0),
      fp_peak                        = character(0),
      sem_direct_tf_expr_est         = numeric(0),
      sem_direct_tf_expr_p           = numeric(0),
      sem_direct_tf_expr_sign        = integer(0),
      sem_direct_tf_expr_direction   = character(0),
      sem_direct_tf_fp_est           = numeric(0),
      sem_direct_tf_fp_p             = numeric(0),
      sem_direct_tf_fp_sign          = integer(0),
      sem_direct_tf_fp_direction     = character(0),
      sem_indirect_ab_est            = numeric(0),
      sem_indirect_ab_p              = numeric(0),
      sem_indirect_ab_sign           = integer(0),
      sem_indirect_ab_direction      = character(0),
      sem_total_est                  = numeric(0),
      sem_total_p                    = numeric(0),
      sem_total_sign                 = integer(0),
      sem_total_direction            = character(0)
    ))
  }

  # Helper to summarise one pattern
  summarise_pattern <- function(tbl, pattern, prefix) {
    sub <- tbl[grepl(pattern, tbl$term), , drop = FALSE]
    if (!nrow(sub)) {
      return(tibble::tibble(
        tf                      = character(0),
        gene                    = character(0),
        fp_peak                 = character(0),
        !!paste0(prefix, "_est")       := numeric(0),
        !!paste0(prefix, "_p")         := numeric(0),
        !!paste0(prefix, "_sign")      := integer(0),
        !!paste0(prefix, "_direction") := character(0)
      ))
    }

    out <- sub |>
      dplyr::group_by(tf, gene, fp_peak) |>
      dplyr::summarise(
        .groups = "drop",
        est = {
          idx <- which.min(p_value)
          estimate[idx]
        },
        p   = min(p_value, na.rm = TRUE)
      )

    sign_vec <- integer(nrow(out))
    dir_vec  <- character(nrow(out))
    for (i in seq_len(nrow(out))) {
      sd <- .sign_and_direction(out$est[i])
      sign_vec[i] <- sd$sign
      dir_vec[i]  <- sd$direction
    }

    out[[paste0(prefix, "_est")]]       <- out$est
    out[[paste0(prefix, "_p")]]         <- out$p
    out[[paste0(prefix, "_sign")]]      <- sign_vec
    out[[paste0(prefix, "_direction")]] <- dir_vec

    out[, c("tf", "gene", "fp_peak",
            paste0(prefix, "_est"),
            paste0(prefix, "_p"),
            paste0(prefix, "_sign"),
            paste0(prefix, "_direction")),
        drop = FALSE]
  }

  sem_direct_expr <- summarise_pattern(
    sem_tbl,
    pattern = "^expr ~ tf_expr$",
    prefix  = "sem_direct_tf_expr"
  )

  sem_direct_fp <- summarise_pattern(
    sem_tbl,
    pattern = "^expr ~ tf_fp$",
    prefix  = "sem_direct_tf_fp"
  )

  sem_indirect <- summarise_pattern(
    sem_tbl,
    pattern = "^ab :=",
    prefix  = "sem_indirect_ab"
  )

  sem_total <- summarise_pattern(
    sem_tbl,
    pattern = "^total :=",
    prefix  = "sem_total"
  )

  # Join all SEM summaries
  base <- sem_tbl[, c("tf", "gene", "fp_peak"), drop = FALSE] |>
    dplyr::distinct()

  base |>
    dplyr::left_join(sem_direct_expr, by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(sem_direct_fp,   by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(sem_indirect,    by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(sem_total,       by = c("tf", "gene", "fp_peak"))
}

# ---- summarise_gam_effects_per_peak() --------------------------------------
#' Internal: summarise GAM smooth effects per peak
#'
#' Aggregates all `model_class == "gam"` terms into three predictor-level
#' summaries per (tf, gene, fp_peak):
#'   - `gam_tf_expr_*`      for smooths with predictor == "tf_expr"
#'   - `gam_tf_fp_*`        for smooths with predictor == "tf_fp"
#'   - `gam_tf_expr_fp_*`   for interaction smooths with predictor == "interaction"
#'
#' For each predictor family, we:
#'   - keep `component == "smooth"` rows,
#'   - require finite `p_value`,
#'   - within each (tf, gene, fp_peak) take the row with smallest p-value,
#'   - store its `statistic` as `*_stat` and that p-value as `*_p`.
#'
#' Note: GAM smooths do not yield a simple signed slope; therefore we do not
#' attempt to create sign/direction columns here.
summarise_gam_effects_per_peak <- function(terms_tf) {
  is_gam <- terms_tf$model_class == "gam" & is.finite(terms_tf$p_value)
  gam_tbl <- terms_tf[is_gam, , drop = FALSE]

  if (!nrow(gam_tbl)) {
    return(tibble::tibble(
      tf                      = character(0),
      gene                    = character(0),
      fp_peak                 = character(0),
      gam_tf_expr_stat        = numeric(0),
      gam_tf_expr_p           = numeric(0),
      gam_tf_fp_stat          = numeric(0),
      gam_tf_fp_p             = numeric(0),
      gam_tf_expr_fp_stat     = numeric(0),
      gam_tf_expr_fp_p        = numeric(0)
    ))
  }

  summarise_gam_family <- function(tbl, predictor_value, prefix) {
    sub <- tbl[
      tbl$component == "smooth" &
        tbl$predictor == predictor_value &
        is.finite(tbl$p_value),
      ,
      drop = FALSE
    ]

    if (!nrow(sub)) {
      return(tibble::tibble(
        tf                      = character(0),
        gene                    = character(0),
        fp_peak                 = character(0),
        !!paste0(prefix, "_stat") := numeric(0),
        !!paste0(prefix, "_p")    := numeric(0)
      ))
    }

    out <- sub |>
      dplyr::group_by(tf, gene, fp_peak) |>
      dplyr::summarise(
        .groups = "drop",
        stat = {
          idx <- which.min(p_value)
          statistic[idx]
        },
        p   = min(p_value, na.rm = TRUE)
      )

    out[[paste0(prefix, "_stat")]] <- out$stat
    out[[paste0(prefix, "_p")]]    <- out$p

    out[, c("tf", "gene", "fp_peak",
            paste0(prefix, "_stat"),
            paste0(prefix, "_p")),
        drop = FALSE]
  }

  gam_expr <- summarise_gam_family(
    gam_tbl,
    predictor_value = "tf_expr",
    prefix          = "gam_tf_expr"
  )

  gam_fp <- summarise_gam_family(
    gam_tbl,
    predictor_value = "tf_fp",
    prefix          = "gam_tf_fp"
  )

  gam_int <- summarise_gam_family(
    gam_tbl,
    predictor_value = "interaction",
    prefix          = "gam_tf_expr_fp"
  )

  base <- gam_tbl[, c("tf", "gene", "fp_peak"), drop = FALSE] |>
    dplyr::distinct()

  base |>
    dplyr::left_join(gam_expr, by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(gam_fp,   by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(gam_int,  by = c("tf", "gene", "fp_peak"))
}

# Master per-peak summary for one TF
# ---- build_peak_gene_summary_for_tf() --------------------------------------
#' Build per-peak per-gene summary table for KO benchmarking (one TF)
#'
#' This function combines:
#'   - footprint vs RNA correlations (`r_fp`, `p_fp`, `p_adj_fp`)
#'   - TF vs gene RNA correlations (`r_rna`, `p_rna`)
#'   - best linear model effects for TF RNA and TF footprint
#'   - selected SEM effects (direct, indirect, total)
#'   - aggregated GAM smooth evidence (F/chi-square statistic + p-value)
#'
#' The input `terms_tf` is the full term table for a single TF, typically
#' produced by [collect_model_terms_for_tf_full()].
#'
#' The output has one row per `(tf, gene, fp_peak)` and wide columns for each
#' model family, ready to be joined with KO truth tables and/or plotted.
#'
#' @param tf Character scalar, TF symbol (must match `terms_tf$tf`).
#' @param terms_tf Tibble with columns at least:
#'   tf, gene, fp_peak, model, model_class, term, predictor, component,
#'   estimate, p_value, statistic.
#' @param tf_gene_links_gh Tibble with GeneHancer-based correlations for all TFs.
#'   Must contain columns:
#'   - `tfs`      : TF symbol
#'   - `gene_key` : target gene symbol
#'   - `fp_peak`  : footprint peak ID
#'   - `r_fp`, `p_fp`, `p_adj_fp`
#' @param rna_tbl RNA tibble (same as `grn_set$rna`), with columns `HGNC` and
#'   per-sample columns.
#'
#' @return Tibble with columns:
#'   - tf, gene, fp_peak
#'   - r_fp, p_fp, p_adj_fp
#'   - r_rna, p_rna
#'   - lin_tf_expr_estimate,  lin_tf_expr_p_value,  lin_tf_expr_sign,
#'     lin_tf_expr_direction
#'   - lin_tf_fp_estimate,    lin_tf_fp_p_value,    lin_tf_fp_sign,
#'     lin_tf_fp_direction
#'   - sem_direct_tf_expr_*, sem_direct_tf_fp_*,
#'     sem_indirect_ab_*, sem_total_*
#'   - gam_tf_expr_stat, gam_tf_expr_p,
#'     gam_tf_fp_stat,   gam_tf_fp_p,
#'     gam_tf_expr_fp_stat, gam_tf_expr_fp_p.
build_peak_gene_summary_for_tf <- function(tf,
                                           terms_tf,
                                           tf_gene_links_gh,
                                           rna_tbl) {
  if (!nrow(terms_tf)) {
    cli::cli_abort("terms_tf is empty for TF {.val {tf}}.")
  }

  # Clean up NA keys early
  terms_tf <- terms_tf[stats::complete.cases(terms_tf[, c("tf", "gene", "fp_peak")]), , drop = FALSE]

  base <- terms_tf[, c("tf", "gene", "fp_peak"), drop = FALSE] |>
    dplyr::distinct()

  # ---------- FP correlations ----------
  fp_corr <- tf_gene_links_gh |>
    dplyr::filter(tfs == tf, !is.na(gene_key), gene_key != "") |>
    dplyr::select(fp_peak, gene = gene_key, r_fp, p_fp, p_adj_fp) |>
    dplyr::distinct(gene, fp_peak, .keep_all = TRUE)

  # ---------- RNA correlations ----------
  genes_use <- unique(base$gene)
  rna_corr  <- compute_tf_rna_corr_for_genes(
    tf      = tf,
    rna_tbl = rna_tbl,
    genes   = genes_use,
    method  = "pearson"
  )

  # ---------- Model summaries ----------
  lin_expr <- summarise_linear_effect_per_peak(terms_tf, "tf_expr", "lin_tf_expr")
  lin_fp   <- summarise_linear_effect_per_peak(terms_tf, "tf_fp",   "lin_tf_fp")
  sem_sum  <- summarise_sem_effects_per_peak(terms_tf)
  gam_sum  <- summarise_gam_effects_per_peak(terms_tf)

  # ---------- Join all summaries ----------
  out <- base |>
    dplyr::left_join(fp_corr,  by = c("gene", "fp_peak")) |>
    dplyr::left_join(rna_corr, by = "gene") |>
    dplyr::left_join(lin_expr, by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(lin_fp,   by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(sem_sum,  by = c("tf", "gene", "fp_peak")) |>
    dplyr::left_join(gam_sum,  by = c("tf", "gene", "fp_peak")) |>
    dplyr::mutate(
      across(dplyr::where(is.numeric),
             ~ dplyr::if_else(is.infinite(.x), NA_real_, .x))
    )

  out
}


# Attach KO truth (per TF) — optional step
# ---- attach_ko_truth_to_peak_summary() -------------------------------------
#' Join per-peak TF-gene summary with KO truth labels
#'
#' @param peak_summary Output of [build_peak_gene_summary_for_tf()].
#' @param ko_truth_tbl KO truth tibble for the same TF, with columns:
#'   - `gene`
#'   - `log2fc`
#'   - `ko_group` (e.g. "Down", "Unchanged", "Up" or NA).
#' @param keep_intermediate Logical; if FALSE, drops rows with `ko_group` NA.
#'
#' @return Tibble with all columns from `peak_summary` plus `log2fc`, `ko_group`.
attach_ko_truth_to_peak_summary <- function(peak_summary,
                                            ko_truth_tbl,
                                            keep_intermediate = FALSE) {
  if (!all(c("gene", "log2fc", "ko_group") %in% names(ko_truth_tbl))) {
    cli::cli_abort("ko_truth_tbl must contain columns 'gene', 'log2fc', 'ko_group'.")
  }

  out <- peak_summary |>
    dplyr::left_join(
      ko_truth_tbl[, c("gene", "log2fc", "ko_group"), drop = FALSE],
      by = "gene",
      relationship = "many-to-many"
    )

  if (!keep_intermediate) {
    out <- out[!is.na(out$ko_group), , drop = FALSE]
  }

  out
}


# usage HNF1A
# Assuming have:
#   - terms_HNF1A from `collect_model_terms_for_tf_full()`
#   - tf_gene_links_gh (GeneHancer correlations)
#   - grn_set$rna      (RNA matrix tibble)
#   - ko_truth_HNF1A   (KO truth table)

# hnf1a_peak_summary <- build_peak_gene_summary_for_tf(
#   tf               = "HNF1A",
#   terms_tf         = terms_HNF1A,
#   tf_gene_links_gh = tf_gene_links_gh,
#   rna_tbl          = grn_set$rna
# )
# str(hnf1a_peak_summary)
# hnf1a_with_truth <- attach_ko_truth_to_peak_summary(
#   peak_summary      = hnf1a_peak_summary,
#   ko_truth_tbl      = ko_truth_HNF1A,
#   keep_intermediate = TRUE
# )
# str(hnf1a_with_truth)
#
# hnf1a_with_truth[hnf1a_with_truth$gene == "ACSF2", ]



# Model-comparison helpers for KO benchmarks

# Assumes an input tibble like `hnf1a_with_truth` with columns:
#   tf, gene, fp_peak,
#   r_fp, p_fp, p_adj_fp,
#   r_rna, p_rna,
#   lin_tf_expr_estimate, lin_tf_expr_p_value, lin_tf_expr_sign,
#   lin_tf_fp_estimate,   lin_tf_fp_p_value,   lin_tf_fp_sign,
#   sem_direct_tf_expr_est, sem_direct_tf_expr_p, sem_direct_tf_expr_sign,
#   sem_direct_tf_fp_est,   sem_direct_tf_fp_p,   sem_direct_tf_fp_sign,
#   sem_indirect_ab_est, sem_indirect_ab_p, sem_indirect_ab_sign,
#   sem_total_est,       sem_total_p,       sem_total_sign,
#   gam_tf_expr_stat, gam_tf_expr_p,
#   gam_tf_fp_stat,   gam_tf_fp_p,
#   gam_tf_expr_fp_stat, gam_tf_expr_fp_p,
#   log2fc, ko_group
#
# The key idea:
#   * For each model, each TF-gene pair can keep at most ONE fp_peak.
#   * For RNA-only models: pick genes where model direction matches RNA corr.
#   * For FP-only models: FP direction must match RNA corr; require p_fp sig.
#   * For RNA+FP models: RNA + FP directions must agree; p_rna, p_fp,
#                        and model p all must be significant.
#   * For models without an explicit sign, we borrow the FP sign.
#   * Optional |r| filters (r_cut_rna, r_cut_fp) are applied by model type.
#   * After filtering, if multiple fp_peaks per gene remain, keep the one
#     with the best (minimum) model p-value.
#
# This yields, per (gene, model), a single row flagged as Predicted/Non.

# ---- build_default_model_defs() --------------------------------------------
#' Build default model-definition list for a TF summary table
#'
#' @param tf_tbl Tibble like `hnf1a_with_truth`.
#' @return Named list of model definitions.
build_default_model_defs <- function(tf_tbl) {
  defs <- list(
    RNA_corr = list(
      label    = "RNA: corr",
      category = "RNA",
      type     = "RNA",
      p_col    = "p_rna",
      sign_col = NA_character_
    ),
    RNA_lin = list(
      label    = "RNA: linear",
      category = "RNA",
      type     = "RNA",
      p_col    = "lin_tf_expr_p_value",
      sign_col = "lin_tf_expr_sign"
    ),
    RNA_SEM_direct = list(
      label    = "RNA: SEM direct",
      category = "RNA",
      type     = "RNA",
      p_col    = "sem_direct_tf_expr_p",
      sign_col = "sem_direct_tf_expr_sign"
    ),
    RNA_GAM = list(
      label    = "RNA: GAM",
      category = "RNA",
      type     = "RNA",
      p_col    = "gam_tf_expr_p",
      sign_col = NA_character_
    ),
    FP_corr = list(
      label    = "FP: corr",
      category = "FP",
      type     = "FP",
      p_col    = "p_fp",
      sign_col = NA_character_
    ),
    FP_lin = list(
      label    = "FP: linear",
      category = "FP",
      type     = "FP",
      p_col    = "lin_tf_fp_p_value",
      sign_col = "lin_tf_fp_sign"
    ),
    FP_SEM_direct = list(
      label    = "FP: SEM direct",
      category = "FP",
      type     = "FP",
      p_col    = "sem_direct_tf_fp_p",
      sign_col = "sem_direct_tf_fp_sign"
    ),
    FP_GAM = list(
      label    = "FP: GAM",
      category = "FP",
      type     = "FP",
      p_col    = "gam_tf_fp_p",
      sign_col = NA_character_
    ),
    SEM_total = list(
      label    = "RNA+FP: SEM total",
      category = "RNA+FP",
      type     = "RNA+FP",
      p_col    = "sem_total_p",
      sign_col = "sem_total_sign"
    ),
    GAM_int = list(
      label    = "RNA+FP: GAM ti()",
      category = "Other",
      type     = "RNA+FP",
      p_col    = "gam_tf_expr_fp_p",
      sign_col = NA_character_
    )
  )

  # Keep only models whose p_col is actually present and has some finite values
  keep <- vapply(
    defs,
    function(m) {
      m$p_col %in% names(tf_tbl) &&
        any(is.finite(tf_tbl[[m$p_col]]), na.rm = TRUE)
    },
    logical(1L)
  )
  defs[keep]
}

# Per-model definitions directly from fitted model terms
# ---- augment_tf_tbl_with_per_model_terms() -----------------------------------
#' Add per-model p/sign columns from `tf_terms` and build model defs
#'
#' @param tf_tbl   Output of [build_peak_gene_summary_for_tf()].
#' @param terms_tf Full term table from `collect_model_terms_for_tf_full()`.
#'
#' @return List with:
#'   - tf_tbl     : augmented tf_tbl with per-model p/sign columns
#'   - model_defs : list of model definitions (labels = raw model names)
augment_tf_tbl_with_per_model_terms <- function(tf_tbl, terms_tf) {
  if (!nrow(terms_tf)) {
    return(list(tf_tbl = tf_tbl, model_defs = list()))
  }

  sanitize_id <- function(x) {
    out <- gsub("[^A-Za-z0-9]+", "_", x)
    out <- gsub("^_+|_+$", "", out)
    if (identical(out, "") || is.na(out)) "model" else out
  }

  classify_model <- function(model_name) {
    model_lower <- tolower(model_name)
    if (grepl("rna_fp", model_lower, fixed = TRUE)) {
      "RNA+FP"
    } else if (grepl("rna", model_lower, fixed = TRUE)) {
      "RNA"
    } else if (grepl("fp", model_lower, fixed = TRUE)) {
      "FP"
    } else {
      "Other"
    }
  }

  tf_aug      <- tf_tbl
  model_defs  <- list()
  models_use  <- sort(unique(terms_tf$model))

  for (m in models_use) {
    m_id <- sanitize_id(m)

    sub <- terms_tf[
      terms_tf$model == m &
        terms_tf$component == "parametric" &
        terms_tf$predictor %in% c("tf_expr", "tf_fp") &
        is.finite(terms_tf$p_value),
      ,
      drop = FALSE
    ]

    if (!nrow(sub)) next

    # Prefer tf_expr when available; otherwise fall back to tf_fp
    pred_use <- if ("tf_expr" %in% sub$predictor) "tf_expr" else "tf_fp"
    sub <- sub[sub$predictor == pred_use, , drop = FALSE]

    best <- sub |>
      dplyr::group_by(tf, gene, fp_peak) |>
      dplyr::slice_min(order_by = p_value, n = 1L, with_ties = FALSE) |>
      dplyr::ungroup()

    if (!nrow(best)) next

    sign_vec <- vapply(
      best$estimate,
      function(est) .sign_and_direction(est)$sign,
      integer(1L)
    )

    col_p    <- paste0("model_", m_id, "_p")
    col_sign <- paste0("model_", m_id, "_sign")

    best_small <- best[, c("gene", "fp_peak", "p_value"), drop = FALSE] |>
      dplyr::mutate(
        !!col_p    := p_value,
        !!col_sign := sign_vec
      ) |>
      dplyr::select(-p_value) |>
      dplyr::distinct(gene, fp_peak, .keep_all = TRUE)

    tf_aug <- tf_aug |>
      dplyr::left_join(best_small, by = c("gene", "fp_peak"))

    model_defs[[m_id]] <- list(
      label    = m,
      category = classify_model(m),
      type     = classify_model(m),
      p_col    = col_p,
      sign_col = col_sign
    )
  }

  list(tf_tbl = tf_aug, model_defs = model_defs)
}

# Order models by category (RNA → FP → RNA+FP → Other), then label
reorder_model_levels <- function(model_defs) {
  if (!length(model_defs)) return(character())
  labs <- vapply(model_defs, function(m) m$label, character(1L))
  cats <- vapply(model_defs, function(m) m$category, character(1L))
  cat_priority <- match(cats, c("RNA", "FP", "RNA+FP", "Other"))
  cat_priority[is.na(cat_priority)] <- length(c("RNA", "FP", "RNA+FP", "Other")) + 1L
  labs[order(cat_priority, labs)]
}


# ---- build_tf_model_prediction_grid() --------------------------------------
#' Sweep over a grid of p-value cutoffs for all models
#'
#' @param tf_tbl Tibble like `hnf1a_with_truth`.
#' @param model_defs List as from build_default_model_defs().
#' @param p_grid Tibble/data.frame with columns:
#'   p_cut_rna, p_cut_fp, p_cut_model, cutoff_label (character).
#'   Typically the three p_cut_* columns are the same sequence of values,
#'   but you can vary them independently if needed.
#' @param r_cut_rna,r_cut_fp As in build_tf_model_prediction_table().
#'
#' @return Tibble like build_tf_model_prediction_table(), plus
#'   p_cut_rna, p_cut_fp, p_cut_model, cutoff_label.
build_tf_model_prediction_grid <- function(tf_tbl,
                                           model_defs = NULL,
                                           p_grid,
                                           r_cut_rna = 0.3,
                                           r_cut_fp  = 0.3) {
  if (is.null(model_defs)) {
    model_defs <- build_default_model_defs(tf_tbl)
  }

  required_pg <- c("p_cut_rna", "p_cut_fp", "p_cut_model", "cutoff_label")
  missing_pg  <- setdiff(required_pg, names(p_grid))
  if (length(missing_pg)) {
    cli::cli_abort(
      c(
        "p_grid is missing required columns:",
        paste(missing_pg, collapse = ", ")
      )
    )
  }

  grid_tbls <- lapply(seq_len(nrow(p_grid)), function(i) {
    row_i <- p_grid[i, , drop = FALSE]

    pred_i <- build_tf_model_prediction_table(
      tf_tbl      = tf_tbl,
      model_defs  = model_defs,
      p_cut_rna   = row_i$p_cut_rna,
      p_cut_fp    = row_i$p_cut_fp,
      p_cut_model = row_i$p_cut_model,
      r_cut_rna   = r_cut_rna,
      r_cut_fp    = r_cut_fp
    )

    pred_i$p_cut_rna    <- row_i$p_cut_rna
    pred_i$p_cut_fp     <- row_i$p_cut_fp
    pred_i$p_cut_model  <- row_i$p_cut_model
    pred_i$cutoff_label <- as.character(row_i$cutoff_label)

    pred_i
  })

  out <- dplyr::bind_rows(grid_tbls)

  out$cutoff_label <- factor(out$cutoff_label, levels = as.character(p_grid$cutoff_label))
  out
}

# ---- build_tf_model_prediction_table() -------------------------------------
#' Build per-gene prediction flags for each model
#'
#' @param tf_tbl Tibble like `hnf1a_with_truth`.
#' @param model_defs List as returned by build_default_model_defs().
#' @param p_cut_rna Numeric, p-value cutoff for RNA-based filters.
#' @param p_cut_fp Numeric, p-value cutoff for FP-based filters.
#' @param p_cut_model Numeric, p-value cutoff for model-specific p.
#' @param r_cut_rna Numeric, |r_rna| >= r_cut_rna to be eligible (default 0.3).
#' @param r_cut_fp Numeric, |r_fp|  >= r_cut_fp  to be eligible (default 0.3).
#'
#' @return Tibble with one row per (gene, model), columns:
#'   gene, ko_group, log2fc, model, model_label, category,
#'   predicted (logical).
build_tf_model_prediction_table <- function(tf_tbl,
                                            model_defs = NULL,
                                            p_cut_rna   = 0.05,
                                            p_cut_fp    = 0.05,
                                            p_cut_model = 0.05,
                                            r_cut_rna   = 0.3,
                                            r_cut_fp    = 0.3) {
  required_cols <- c("gene", "log2fc", "ko_group", "r_fp", "p_fp", "r_rna", "p_rna")
  missing_cols <- setdiff(required_cols, names(tf_tbl))
  if (length(missing_cols)) {
    cli::cli_abort(
      c(
        "tf_tbl is missing required columns:",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (is.null(model_defs)) {
    model_defs <- build_default_model_defs(tf_tbl)
  }

  base_gene_tbl <- tf_tbl |>
    dplyr::group_by(gene) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(gene, log2fc, ko_group)

  tf_tbl <- tf_tbl |>
    dplyr::mutate(
      rna_sign = dplyr::case_when(
        r_rna > 0  ~  1L,
        r_rna < 0  ~ -1L,
        TRUE       ~  0L
      ),
      fp_sign  = dplyr::case_when(
        r_fp > 0   ~  1L,
        r_fp < 0   ~ -1L,
        TRUE       ~  0L
      ),
      lfc_sign = dplyr::case_when(
        log2fc > 0 ~  1L,
        log2fc < 0 ~ -1L,
        TRUE       ~  0L
      )
    )

  all_genes <- unique(tf_tbl$gene)

  model_tbls <- lapply(names(model_defs), function(model_id) {
    m <- model_defs[[model_id]]

    if (!m$p_col %in% names(tf_tbl)) {
      return(tibble::tibble())
    }

    df <- tf_tbl
    df$model_p <- df[[m$p_col]]

    if (!is.null(m$sign_col) && !is.na(m$sign_col) && m$sign_col %in% names(df)) {
      df$model_sign <- df[[m$sign_col]]
    } else if (identical(m$p_col, "p_rna")) {
      # RNA correlation: sign comes from r_rna itself
      df$model_sign <- df$rna_sign
    } else {
      df$model_sign <- if (identical(m$type, "RNA")) df$rna_sign else df$fp_sign
    }

    # eligibility: finite p_model, non-zero RNA sign
    df <- df[is.finite(df$model_p) & df$rna_sign != 0L, , drop = FALSE]
    if (!nrow(df)) {
      return(
        base_gene_tbl |>
          dplyr::mutate(
            model_id    = model_id,
            model_label = m$label,
            category    = m$category,
            predicted   = FALSE
          )
      )
    }

    df$ok_r_rna <- abs(df$r_rna) >= r_cut_rna
    df$ok_r_fp  <- abs(df$r_fp)  >= r_cut_fp

    keep <- switch(
      m$type,
      "RNA" = {
        pcut <- if (identical(m$p_col, "p_rna")) p_cut_rna else p_cut_model

        df$model_sign == df$rna_sign &
          df$model_p   <= pcut &
          df$ok_r_rna
      },
      "FP" = {
        df$fp_sign     == df$rna_sign &
          df$model_sign == df$rna_sign &
          df$p_fp       <= p_cut_fp &
          df$model_p    <= p_cut_model &
          df$ok_r_fp &
          df$ok_r_rna
      },
      "RNA+FP" = {
        df$fp_sign     == df$rna_sign &
          df$model_sign == df$rna_sign &
          df$p_fp       <= p_cut_fp &
          df$p_rna      <= p_cut_rna &
          df$model_p    <= p_cut_model &
          df$ok_r_fp &
          df$ok_r_rna
      },
      {
        df$fp_sign     == df$rna_sign &
          df$model_sign == df$rna_sign &
          df$model_p    <= p_cut_model &
          df$ok_r_fp &
          df$ok_r_rna
      }
    )

    df <- df[keep, , drop = FALSE]

    if (!nrow(df)) {
      return(
        base_gene_tbl |>
          dplyr::mutate(
            model_id    = model_id,
            model_label = m$label,
            category    = m$category,
            predicted   = FALSE
          )
      )
    }

    if (model_id %in% c("RNA_corr", "FP_corr")) {
      df_best <- df |>
        dplyr::group_by(gene) |>
        dplyr::slice_min(order_by = model_p, n = 1L, with_ties = FALSE) |>
        dplyr::ungroup()
    } else {
      # RNA models must NOT use fp_sign
      df_best <- df |>
        dplyr::mutate(
          dir_mismatch = ifelse(
            if (identical(m$type, "RNA")) {
              lfc_sign != 0L & model_sign != 0L & lfc_sign != model_sign
            } else {
              lfc_sign != 0L & fp_sign != 0L & lfc_sign != fp_sign
            },
            1L, 0L
          )
        ) |>
        dplyr::group_by(gene) |>
        dplyr::arrange(dir_mismatch, model_p, .by_group = TRUE) |>
        dplyr::slice(1L) |>
        dplyr::ungroup()
    }

    predicted_genes <- df_best$gene
    non_genes       <- setdiff(all_genes, predicted_genes)

    pred_tbl <- df_best |>
      dplyr::select(gene, model_p) |>
      dplyr::left_join(base_gene_tbl, by = "gene") |>
      dplyr::mutate(
        model_id    = model_id,
        model_label = m$label,
        category    = m$category,
        predicted   = TRUE
      )

    if (length(non_genes)) {
      non_tbl <- base_gene_tbl[base_gene_tbl$gene %in% non_genes, , drop = FALSE] |>
        dplyr::mutate(
          model_p     = NA_real_,
          model_id    = model_id,
          model_label = m$label,
          category    = m$category,
          predicted   = FALSE
        )
      dplyr::bind_rows(pred_tbl, non_tbl)
    } else {
      pred_tbl
    }
  })

  out <- dplyr::bind_rows(model_tbls)

  model_levels <- vapply(model_defs, function(m) m$label, character(1L))
  out$model    <- factor(out$model_label, levels = model_levels)
  out$category <- factor(out$category, levels = c("RNA", "FP", "RNA+FP", "Other"))
  out$pred_group <- factor(
    ifelse(out$predicted, "Predicted", "Non-predicted"),
    levels = c("Predicted", "Non-predicted")
  )

  out
}



# Single-cutoff composite plot
# ---- plot_tf_model_comparison() --------------------------------------------
plot_tf_model_comparison <- function(df_model,
                                     tf_label   = "HNF1A",
                                     out_file   = NULL,
                                     lfc_breaks = c(-Inf, -1, -0.5, 0, Inf),
                                     lfc_labels = c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")) {
  stopifnot(all(c("model", "predicted", "log2fc") %in% names(df_model)))

  df_use <- df_model[is.finite(df_model$log2fc), , drop = FALSE]

  ## Preserve original model ordering from df_model$model
  model_levels <- if (is.factor(df_model$model)) {
    levels(df_model$model)
  } else {
    unique(df_model$model)
  }
  df_use$model <- factor(df_use$model, levels = model_levels)

  df_use$group <- ifelse(df_use$predicted, "Predicted", "Non-predicted")
  df_use$group <- factor(df_use$group, levels = c("Predicted", "Non-predicted"))

  # Log2FC bin (for bottom panel)
  df_use$log2fc_bin <- cut(
    df_use$log2fc,
    breaks = lfc_breaks,
    labels = lfc_labels,
    right  = TRUE,
    include.lowest = TRUE
  )

  lfc_cols <- c(
    "<= -1"     = "#d73027",
    "(-1,-0.5]" = "#fc8d59",
    "(-0.5,0]"  = "#fee090",
    ">= 0"      = "#d9d9d9"
  )

  # 1) Violin panel
  p_violin <- ggplot2::ggplot(
    df_use,
    ggplot2::aes(x = group, y = log2fc, fill = group)
  ) +
    ggplot2::geom_violin(
      trim   = FALSE,
      alpha  = 0.4,
      colour = "black"
    ) +
    ggplot2::geom_boxplot(
      width         = 0.15,
      outlier.size  = 0.2,
      outlier.alpha = 0.4
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::scale_fill_manual(
      values = c("Predicted" = "#66c2a5", "Non-predicted" = "#fc8d62"),
      name   = "Group"
    ) +
    ggplot2::labs(
      title = sprintf("%s KO - Model comparison", tf_label),
      x     = NULL,
      y     = "log2FC (KO vs Ctrl)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # 2) Count panel: number of predicted target genes per model
  df_counts <- df_use[df_use$predicted, , drop = FALSE] |>
    dplyr::distinct(model, gene) |>
    dplyr::count(model, name = "n_predicted") |>
    dplyr::mutate(x_dummy = "Predicted")

  p_counts <- ggplot2::ggplot(
    df_counts,
    ggplot2::aes(x = x_dummy, y = n_predicted)
  ) +
    ggplot2::geom_col(fill = "grey70") +
    ggplot2::geom_text(
      ggplot2::aes(label = n_predicted),
      vjust = -0.3,
      size  = 3
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of predicted target genes"
    ) +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 0, r = 5, b = 0, l = 5)
    )

  # 3) Percent panel: log2FC-bin composition within each group (Pred/Non)
  df_percent <- df_use[!is.na(df_use$log2fc_bin), , drop = FALSE] |>
    dplyr::distinct(model, gene, group, log2fc_bin) |>
    dplyr::count(model, group, log2fc_bin, name = "n") |>
    dplyr::group_by(model, group) |>
    dplyr::mutate(
      pct = 100 * n / sum(n)
    ) |>
    dplyr::ungroup()

  p_percent <- ggplot2::ggplot(
    df_percent,
    ggplot2::aes(x = group, y = pct, fill = log2fc_bin)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(
      values = lfc_cols,
      name   = "log2FC bin"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Percent of genes within group"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5),
      plot.margin      = ggplot2::margin(t = 0, r = 5, b = 5, l = 5),
      legend.position  = "right"
    )

  # Assemble with patchwork
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort("Package 'patchwork' is required for assembling the plot.")
  }

  p_all <- p_violin / p_counts / p_percent +
    patchwork::plot_layout(heights = c(3, 1, 1), guides = "collect") &
    ggplot2::theme(legend.position = "right")

  if (!is.null(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot     = p_all,
      width    = 24,
      height   = 7,
      units    = "in",
      dpi      = 300
    )
  }

  invisible(p_all)
}

# filter a prediction grid by TFLink prior (Name.TF/Name.Target)
# ---- filter_tf_pred_grid_by_TFLink() ---------------------------------------
filter_tf_pred_grid_by_TFLink <- function(pred_grid_tbl,
                                          tf_name,
                                          tf_link_tbl,
                                          semi            = FALSE,
                                          max_non_tf_link = 1000L) {
  # Basic structure checks
  required_cols <- c("gene")
  miss_pred <- setdiff(required_cols, names(pred_grid_tbl))
  if (length(miss_pred)) {
    cli::cli_abort(
      c(
        "pred_grid_tbl is missing required columns:",
        paste(miss_pred, collapse = ", ")
      )
    )
  }

  if (!all(c("Name.TF", "Name.Target") %in% names(tf_link_tbl))) {
    cli::cli_abort("tf_link_tbl must contain 'Name.TF' and 'Name.Target'.")
  }

  tf_targets <- unique(tf_link_tbl$Name.Target[tf_link_tbl$Name.TF == tf_name])

  # Strict mode: TFLink-only
  if (!semi) {
    out <- pred_grid_tbl[pred_grid_tbl$gene %in% tf_targets, , drop = FALSE]

    if (!nrow(out)) {
      cli::cli_warn(
        "No TFLink matches found for TF {.val {tf_name}} (semi = FALSE). Returning empty tibble."
      )
    }

    return(out)
  }

  # Semi mode:  get TFLink+extra;
  needed_semi <- c("predicted", "model_id", "model_p", "log2fc")
  miss_semi <- setdiff(needed_semi, names(pred_grid_tbl))
  if (length(miss_semi)) {
    cli::cli_abort(
      c(
        "For semi mode, pred_grid_tbl is missing required columns:",
        paste(miss_semi, collapse = ", ")
      )
    )
  }

  is_sem_gam <- grepl("SEM", pred_grid_tbl$model_id) |
    grepl("GAM", pred_grid_tbl$model_id)

  other_models <- pred_grid_tbl[!is_sem_gam, , drop = FALSE]
  sem_gam_tbl <- pred_grid_tbl[is_sem_gam, , drop = FALSE]

  core_tf_sem_gam <- sem_gam_tbl[sem_gam_tbl$gene %in% tf_targets,
                                 , drop = FALSE]
  non_tf_sem_gam <- sem_gam_tbl[!(sem_gam_tbl$gene %in% tf_targets) &
                                  sem_gam_tbl$predicted,
                                , drop = FALSE]

  if (nrow(non_tf_sem_gam) > 0L) {
    non_tf_best <- non_tf_sem_gam |>
      dplyr::group_by(model_id) |>
      dplyr::arrange(
        dplyr::if_else(is.finite(model_p), 0L, 1L),
        model_p,
        dplyr::desc(abs(log2fc)),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(n = max_non_tf_link) |>
      dplyr::ungroup()
  } else {
    non_tf_best <- non_tf_sem_gam
  }

  out <- dplyr::bind_rows(
    other_models,
    core_tf_sem_gam,
    non_tf_best
  )

  if (!nrow(out)) {
    cli::cli_warn(
      "No rows passed the TFLink semi-filter for TF {.val {tf_name}}. Returning empty tibble."
    )
  }

  out
}

# Multi-cutoff composite plot
# Columns = models; within each model, x = cutoff × Group
# ---- plot_tf_model_comparison_grid() ---------------------------------------
plot_tf_model_comparison_grid <- function(df_model,
                                          tf_label   = "HNF1A",
                                          out_file   = NULL,
                                          lfc_breaks = c(-Inf, -1, -0.5, 0, Inf),
                                          lfc_labels = c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")) {
  stopifnot(all(c("model", "predicted", "log2fc", "cutoff_label") %in% names(df_model)))

  df_use <- df_model[is.finite(df_model$log2fc), , drop = FALSE]

  # Preserve original model ordering from df_model$model
  model_levels <- if (is.factor(df_model$model)) {
    levels(df_model$model)
  } else {
    unique(df_model$model)
  }
  df_use$model <- factor(df_use$model, levels = model_levels)

  # Preserve cutoff_label order from df_model
  cutoff_levels_input <- if (is.factor(df_model$cutoff_label)) {
    levels(df_model$cutoff_label)
  } else {
    unique(df_model$cutoff_label)
  }
  df_use$cutoff_label <- factor(df_use$cutoff_label, levels = cutoff_levels_input)

  df_use$group <- ifelse(df_use$predicted, "Predicted", "Non-predicted")
  df_use$group <- factor(df_use$group, levels = c("Predicted", "Non-predicted"))

  # Combined x-level: cutoff × group (so row 3 is split by group)
  cutoff_levels <- levels(df_use$cutoff_label)
  group_levels  <- levels(df_use$group)
  x_levels      <- as.character(
    unlist(lapply(cutoff_levels, function(cl) paste(cl, group_levels, sep = " | ")))
  )
  df_use$x_group <- factor(
    paste(df_use$cutoff_label, df_use$group, sep = " | "),
    levels = x_levels
  )

  # Log2FC bin
  df_use$log2fc_bin <- cut(
    df_use$log2fc,
    breaks = lfc_breaks,
    labels = lfc_labels,
    right  = TRUE,
    include.lowest = TRUE
  )

  # Colours similar to your SOX9 panel
  lfc_cols <- c(
    "<= -1"     = "#d73027",
    "(-1,-0.5]" = "#fc8d59",
    "(-0.5,0]"  = "#fee090",
    ">= 0"      = "#d9d9d9"
  )

  # Shorter, two-line x labels: "Pred"/"Non" on top, p-value below
  x_lab_fun <- function(x) {
    parts <- strsplit(x, " | ", fixed = TRUE)
    vapply(
      parts,
      function(xx) {
        grp <- xx[2]
        cut <- xx[1]
        grp_short <- if (grp == "Predicted") "Pred" else "Non"
        paste(grp_short, cut, sep = "\n")
      },
      character(1L)
    )
  }

  # 1) Violin panel
  p_violin <- ggplot2::ggplot(
    df_use,
    ggplot2::aes(x = x_group, y = log2fc, fill = group)
  ) +
    ggplot2::geom_violin(
      trim   = FALSE,
      alpha  = 0.4,
      colour = "black"
    ) +
    ggplot2::geom_boxplot(
      width         = 0.15,
      outlier.size  = 0.2,
      outlier.alpha = 0.4
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::scale_fill_manual(
      values = c("Predicted" = "#66c2a5", "Non-predicted" = "#fc8d62"),
      name   = "Group"
    ) +
    ggplot2::scale_x_discrete(labels = NULL) +
    ggplot2::labs(
      title = sprintf("%s KO - Model comparison across p-value cutoffs", tf_label),
      x     = NULL,
      y     = "log2FC (KO vs Ctrl)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text       = ggplot2::element_text(face = "bold"),
      axis.text.x      = ggplot2::element_blank(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # 2) Count panel: # predicted genes per model × cutoff
  df_counts <- df_use[df_use$predicted, , drop = FALSE] |>
    dplyr::distinct(model, gene, cutoff_label) |>
    dplyr::count(model, cutoff_label, name = "n_predicted")

  df_counts <- tidyr::complete(
    df_counts,
    model = levels(df_use$model),
    cutoff_label = levels(df_use$cutoff_label),
    fill = list(n_predicted = 0)
  )
  p_counts <- ggplot2::ggplot(
    df_counts,
    ggplot2::aes(x = cutoff_label, y = n_predicted)
  ) +
    ggplot2::geom_col(fill = "grey70") +
    ggplot2::geom_text(
      ggplot2::aes(label = n_predicted),
      vjust = -0.3,
      size  = 2.5
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of predicted target genes"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 0, r = 5, b = 0, l = 5)
    )

  # 3) Percent panel: log2FC-bin composition within group × cutoff
  df_percent <- df_use[!is.na(df_use$log2fc_bin), , drop = FALSE] |>
    dplyr::distinct(model, cutoff_label, gene, group, log2fc_bin) |>
    dplyr::count(model, cutoff_label, group, log2fc_bin, name = "n") |>
    dplyr::group_by(model, cutoff_label, group) |>
    dplyr::mutate(
      pct = 100 * n / sum(n)
    ) |>
    dplyr::ungroup()

  df_percent$x_group <- factor(
    paste(df_percent$cutoff_label, df_percent$group, sep = " | "),
    levels = x_levels
  )

  p_percent <- ggplot2::ggplot(
    df_percent,
    ggplot2::aes(x = x_group, y = pct, fill = log2fc_bin)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(1),
      cols = ggplot2::vars(model)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_x_discrete(labels = x_lab_fun) +
    ggplot2::scale_fill_manual(
      values = lfc_cols,
      name   = "log2FC bin"
    ) +
    ggplot2::labs(
      x = "Adjusted p-value cutoff × Group",
      y = "Percent of genes within group"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size  = 7
      ),
      plot.margin      = ggplot2::margin(t = 0, r = 5, b = 5, l = 5),
      legend.position  = "right"
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort("Package 'patchwork' is required for assembling the plot.")
  }

  p_all <- p_violin / p_counts / p_percent +
    patchwork::plot_layout(heights = c(3, 1, 1), guides = "collect") &
    ggplot2::theme(legend.position = "right")

  if (!is.null(out_file)) {
    ggplot2::ggsave(
      filename = out_file,
      plot     = p_all,
      width    = 16,
      height   = 8.5,
      units    = "in",
      dpi      = 300
    )
  }

  invisible(p_all)
}



# plot order
desired_model_levels <- c(
  "RNA: corr",
  "RNA: linear",
  "RNA: SEM direct",
  "RNA: GAM",
  "FP: corr",
  "FP: linear",
  "FP: SEM direct",
  "FP: GAM",
  "RNA+FP: SEM total",
  "RNA+FP: GAM ti()"
)

# Define model defs and base single-cutoff table
# hnf1a_model_defs <- build_default_model_defs(hnf1a_with_truth)
# hnf1a_pred_tbl <- build_tf_model_prediction_table(
#   tf_tbl      = hnf1a_with_truth,
#   model_defs  = hnf1a_model_defs,
#   p_cut_rna   = 0.05,
#   p_cut_fp    = 0.05,
#   p_cut_model = 0.05,
#   r_cut_rna   = 0.3,
#   r_cut_fp    = 0.3
# )
# hnf1a_pred_tbl$model <- hnf1a_pred_tbl$model_label

# Single-cutoff composite (Figure 1 style)
# p_hnf1a_single <- plot_tf_model_comparison(
#   df_model = hnf1a_pred_tbl,
#   tf_label = "HNF1A",
#   out_file = file.path(ko_dir, "HNF1A_model_comparison.pdf")
# )



# 1) Existing unfiltered grid
# p_vals <- c(1e-2, 1e-3, 1e-4, 1e-5)
# p_grid_all <- tibble::tibble(
#   p_cut_rna    = p_vals,
#   p_cut_fp     = p_vals,
#   p_cut_model  = p_vals,
#   cutoff_label = paste0("1e-", -log10(p_vals))
# )
#
# hnf1a_pred_grid <- build_tf_model_prediction_grid(
#   tf_tbl     = hnf1a_with_truth,
#   model_defs = hnf1a_model_defs,
#   p_grid     = p_grid_all,
#   r_cut_rna  = 0.3,
#   r_cut_fp   = 0.3
# )
# hnf1a_pred_grid$model <- hnf1a_pred_grid$model_label

# p_hnf1a_grid <- plot_tf_model_comparison_grid(
#   df_model = hnf1a_pred_grid,
#   tf_label = "HNF1A",
#   out_file = file.path(ko_dir, "HNF1A_model_comparison_pgrid.pdf")
# )

# 2) TFLink-filtered version and PDF
# hnf1a_pred_grid_tfLink <- filter_tf_pred_grid_by_TFLink(
#   pred_grid_tbl = hnf1a_pred_grid,
#   tf_name       = "HNF1A",
#   tf_link_tbl   = TFLink
# )
# hnf1a_pred_grid_tfLink$model <- hnf1a_pred_grid_tfLink$model_label
#
# p_hnf1a_grid_tfLink <- plot_tf_model_comparison_grid(
#   df_model = hnf1a_pred_grid_tfLink,
#   tf_label = "HNF1A (with TFLink)",
#   out_file = file.path(ko_dir, "HNF1A_model_comparison_pgrid_TFLink.pdf")
# )
#
# hnf1a_pred_grid_tfLink <- filter_tf_pred_grid_by_TFLink(
#   pred_grid_tbl = hnf1a_pred_grid,
#   tf_name       = "HNF1A",
#   tf_link_tbl   = TFLink,
#   semi          = TRUE,
#   max_non_tf_link = 500
# )
# hnf1a_pred_grid_tfLink$model <- hnf1a_pred_grid_tfLink$model_label
# hnf1a_pred_grid_tfLink$model <- factor(hnf1a_pred_grid_tfLink$model, levels = desired_model_levels)
#
# p_hnf1a_grid_tfLink <- plot_tf_model_comparison_grid(
#   df_model = hnf1a_pred_grid_tfLink,
#   tf_label = "HNF1A",
#   out_file = file.path(ko_dir, "HNF1A_model_comparison_pgrid.pdf")
# )


# loop over TFs of interest
tfs_interest <- c("HNF1A", "SOX9") # "IRF1", "FOXA2", "HNF4A", "RARG", "KLF5", "SOX9", "HNF1A"
# Desired facet order
desired_model_levels <- c(
  "RNA: corr",
  "RNA: linear",
  "RNA: SEM direct",
  "RNA: GAM",
  "FP: corr",
  "FP: linear",
  "FP: SEM direct",
  "FP: GAM",
  "RNA+FP: SEM total",
  "RNA+FP: GAM ti()"
)

# P-value grid
p_vals <- c(1e-2, 1e-3, 1e-4, 1e-5)
p_grid_all <- tibble::tibble(
  p_cut_rna    = p_vals,
  p_cut_fp     = p_vals,
  p_cut_model  = p_vals,
  cutoff_label = paste0("1e-", -log10(p_vals))
)

for (tf in tfs_interest) {
  message("Processing TF: ", tf)

  # Load models and collect terms for this TF
  tf_result <- load_tf_models_from_blocks(cache_dir_models, tf)
  tf_terms  <- collect_model_terms_for_tf_full(tf_result, 0.05, mc_cores = 16)

  if (!nrow(tf_terms)) { cli::cli_warn("No terms found for TF {.val {tf}}; skipping.")
    next }

  # Build per-peak summary and attach KO truth
  peak_summary <- build_peak_gene_summary_for_tf(
    tf               = tf,
    terms_tf         = tf_terms,
    tf_gene_links_gh = tf_gene_links_gh,
    rna_tbl          = grn_set$rna
  )

  if (!nrow(peak_summary)) {
    cli::cli_warn("Empty peak summary for TF {.val {tf}}; skipping.")
    next
  }

  if (!tf %in% names(ko_truth_list)) {
    cli::cli_warn("No KO truth found in ko_truth_list for TF {.val {tf}}; skipping.")
    next
  }

  tf_with_truth <- attach_ko_truth_to_peak_summary(
    peak_summary      = peak_summary,
    ko_truth_tbl      = ko_truth_list[[tf]],
    keep_intermediate = TRUE
  )

  if (!nrow(tf_with_truth)) {
    cli::cli_warn("No KO-labelled rows for TF {.val {tf}}; skipping.")
    next
  }

  # Build model defs (core set) and add one entry per fitted model from tf_terms
  tf_model_defs_core <- build_default_model_defs(tf_with_truth)

  aug_models <- augment_tf_tbl_with_per_model_terms(tf_with_truth, tf_terms)
  tf_with_truth_aug <- aug_models$tf_tbl

  tf_model_defs <- if (length(aug_models$model_defs)) {
    c(aug_models$model_defs, tf_model_defs_core)
  } else {
    tf_model_defs_core
  }
  model_levels <- reorder_model_levels(tf_model_defs)

  # Multi-cutoff grid for this TF
  tf_pred_grid <- build_tf_model_prediction_grid(
    tf_tbl     = tf_with_truth_aug,
    model_defs = tf_model_defs,
    p_grid     = p_grid_all,
    r_cut_rna  = 0.3,
    r_cut_fp   = 0.3
  )
  tf_pred_grid$model <- factor(tf_pred_grid$model_label, levels = model_levels)

  # plot Multi-cutoff grid
  tf_pred_grid_tfLink <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl   = tf_pred_grid,
    tf_name         = tf,
    tf_link_tbl     = TFLink,
    semi            = TRUE)

  tf_pred_grid_tfLink$model <- factor(
    tf_pred_grid_tfLink$model_label,
    levels = model_levels
  )

  # Plot p-grid with TFLink filter
  out_grid_tfLink <- file.path(
    ko_dir,
    sprintf("%s_model_comparison_p_grid.pdf", tf)
  )
  plot_tf_model_comparison_grid(
    df_model = tf_pred_grid_tfLink,
    tf_label = tf,
    out_file = out_grid_tfLink
  )

  # TFLink filtered version
  tf_pred_grid_tfLink <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl   = tf_pred_grid,
    tf_name         = tf,
    tf_link_tbl     = TFLink
  )
  tf_pred_grid_tfLink$model <- factor(
    tf_pred_grid_tfLink$model_label,
    levels = model_levels
  )

  # Plot p-grid with TFLink filter
  out_grid_tfLink <- file.path(
    ko_dir,
    sprintf("%s_model_comparison_p_grid_TFLink.pdf", tf)
  )
  plot_tf_model_comparison_grid(
    df_model = tf_pred_grid_tfLink,
    tf_label = tf,
    out_file = out_grid_tfLink
  )
}



# RNA-only: stricter p-value grid
p_vals_rna_strict <- c(1e-16, 1e-17, 1e-18, 1e-19)

p_grid_rna_strict <- tibble::tibble(
  p_cut_rna    = p_vals_rna_strict,
  p_cut_fp     = p_vals_rna_strict,
  p_cut_model  = p_vals_rna_strict,
  cutoff_label = paste0("1e-", -log10(p_vals_rna_strict))
)

for (tf in tfs_interest) {
  message("Processing TF: ", tf)

  # Load models and collect terms for this TF
  tf_result <- load_tf_models_from_blocks(cache_dir_models, tf)
  tf_terms  <- collect_model_terms_for_tf_full(tf_result, 0.05, mc_cores = 16)

  if (!nrow(tf_terms)) { cli::cli_warn("No terms found for TF {.val {tf}}; skipping.")
    next }

  # Build per-peak summary and attach KO truth
  peak_summary <- build_peak_gene_summary_for_tf(
    tf               = tf,
    terms_tf         = tf_terms,
    tf_gene_links_gh = tf_gene_links_gh,
    rna_tbl          = grn_set$rna
  )

  if (!nrow(peak_summary)) {
    cli::cli_warn("Empty peak summary for TF {.val {tf}}; skipping.")
    next
  }

  if (!tf %in% names(ko_truth_list)) {
    cli::cli_warn("No KO truth found in ko_truth_list for TF {.val {tf}}; skipping.")
    next
  }

  tf_with_truth <- attach_ko_truth_to_peak_summary(
    peak_summary      = peak_summary,
    ko_truth_tbl      = ko_truth_list[[tf]],
    keep_intermediate = TRUE
  )

  if (!nrow(tf_with_truth)) {
    cli::cli_warn("No KO-labelled rows for TF {.val {tf}}; skipping.")
    next
  }

  # Build model defs (core) and per-model additions
  tf_model_defs_core <- build_default_model_defs(tf_with_truth)

  aug_models <- augment_tf_tbl_with_per_model_terms(tf_with_truth, tf_terms)
  tf_with_truth_aug <- aug_models$tf_tbl

  tf_model_defs <- if (length(aug_models$model_defs)) {
    c(aug_models$model_defs, tf_model_defs_core)
  } else {
    tf_model_defs_core
  }
  model_levels <- reorder_model_levels(tf_model_defs)
  model_categories <- vapply(tf_model_defs, function(m) m$category, character(1L))

  # Multi-cutoff grid for this TF
  tf_pred_grid <- build_tf_model_prediction_grid(
    tf_tbl     = tf_with_truth_aug,
    model_defs = tf_model_defs,
    p_grid     = p_grid_rna_strict,
    r_cut_rna  = 0.3,
    r_cut_fp   = 0.3
  )
  tf_pred_grid$model <- factor(tf_pred_grid$model_label, levels = model_levels)

  # --- (1) RNA-only plots with stricter p-value grid (1e-16..1e-19) ---
  # Reuse tf_pred_grid (already built with p_grid_rna_strict above)
  tf_pred_grid_rna_strict <- tf_pred_grid

  # Keep RNA models only
  is_rna_model <- as.character(tf_pred_grid_rna_strict$category) == "RNA"

  tf_pred_grid_rna_only <- tf_pred_grid_rna_strict[is_rna_model, , drop = FALSE]

  # Ensure consistent facet order for RNA-only subset
  desired_rna_levels <- model_levels[model_categories == "RNA"]
  tf_pred_grid_rna_only$model <- factor(
    tf_pred_grid_rna_only$model_label,
    levels = desired_rna_levels
  )

  tf_pred_grid_rna_only_semi <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl = tf_pred_grid_rna_only,
    tf_name       = tf,
    tf_link_tbl   = TFLink,
    semi          = TRUE
  )
  tf_pred_grid_rna_only_semi$model <- factor(
    tf_pred_grid_rna_only_semi$model_label,
    levels = desired_rna_levels
  )

  if (nrow(tf_pred_grid_rna_only_semi) == 0L) {
    cli::cli_warn("No RNA-only data after TFLink filtering for TF {.val {tf}}; skipping RNA-only plot.")
  } else {
    out_grid_rna_only <- file.path(
      ko_dir,
      sprintf("%s_model_comparison_pgrid_RNAonly_p1e-16to1e-19.pdf", tf)
    )
    plot_tf_model_comparison_grid(
      df_model = tf_pred_grid_rna_only_semi,
      tf_label = tf,
      out_file = out_grid_rna_only
    )
  }


  # --- (2) RNA-only plot WITH an additional FP evidence filter (per-cutoff) ---
  # Definition (minimal + concrete):
  # keep RNA-model predictions only if the gene has ANY peak with:
  #   p_fp <= cutoff's p_cut_fp
  #   |r_fp| >= r_cut_fp (0.3)
  #   fp_sign == rna_sign (directional agreement)
  tf_sign_tbl <- tf_with_truth_aug |>
    dplyr::mutate(
      rna_sign = dplyr::case_when(
        r_rna > 0  ~  1L,
        r_rna < 0  ~ -1L,
        TRUE       ~  0L
      ),
      fp_sign = dplyr::case_when(
        r_fp > 0   ~  1L,
        r_fp < 0   ~ -1L,
        TRUE       ~  0L
      )
    )
  p_fp_filter_fixed <- 0.01

  fp_ok_tbl <- lapply(seq_len(nrow(p_grid_rna_strict)), function(i) {
    # p_fp_i <- p_grid_rna_strict$p_cut_fp[i] # filter by corresponding fp p-value

    p_fp_i <- p_fp_filter_fixed # filter by a fixed fp p-value

    lab_i  <- as.character(p_grid_rna_strict$cutoff_label[i])

    ok_genes <- unique(tf_sign_tbl$gene[
      is.finite(tf_sign_tbl$p_fp) &
        tf_sign_tbl$p_fp <= p_fp_i &
        is.finite(tf_sign_tbl$r_fp) &
        abs(tf_sign_tbl$r_fp) >= 0.3 &
        tf_sign_tbl$rna_sign != 0L &
        tf_sign_tbl$fp_sign  != 0L &
        tf_sign_tbl$fp_sign == tf_sign_tbl$rna_sign
    ])

    tibble::tibble(gene = ok_genes, cutoff_label = lab_i)
  }) |>
    dplyr::bind_rows() |>
    dplyr::distinct(gene, cutoff_label) |>
    dplyr::mutate(fp_ok = TRUE)

  # Keep cutoff_label factor levels stable for plotting
  fp_ok_tbl$cutoff_label <- factor(
    fp_ok_tbl$cutoff_label,
    levels = levels(tf_pred_grid_rna_only_semi$cutoff_label)
  )

  tf_pred_grid_rna_only_fp <- tf_pred_grid_rna_only |>
    dplyr::left_join(fp_ok_tbl, by = c("gene", "cutoff_label")) |>
    dplyr::mutate(
      fp_ok     = !is.na(fp_ok),
      predicted = predicted & fp_ok
    ) |>
    dplyr::select(-fp_ok)

  tf_pred_grid_rna_only_fp_semi <- filter_tf_pred_grid_by_TFLink(
    pred_grid_tbl = tf_pred_grid_rna_only_fp,
    tf_name       = tf,
    tf_link_tbl   = TFLink,
    semi          = TRUE
  )
  tf_pred_grid_rna_only_fp_semi$model <- factor(
    tf_pred_grid_rna_only_fp_semi$model_label,
    levels = desired_rna_levels
  )

  if (nrow(tf_pred_grid_rna_only_fp_semi) == 0L) {
    cli::cli_warn("No RNA-only+FP data after TFLink filtering for TF {.val {tf}}; skipping RNA-only+FP plot.")
  } else {
    out_grid_rna_only_fp <- file.path(
      ko_dir,
      sprintf("%s_model_comparison_pgrid_RNAonly_p1e-16to1e-19_FPfilter_p0.01.pdf", tf)
    )
    plot_tf_model_comparison_grid(
      df_model = tf_pred_grid_rna_only_fp_semi,
      tf_label = tf,
      out_file = out_grid_rna_only_fp
    )
  }
}

# ---- 7. Minimal p-value histogram loop (HNF1A + SOX9) ----------------------
#
# Assumes these already exist:
# - cache_dir_models, tf_gene_links_gh, grn_set$rna, ko_truth_list
# - ko_dir (output directory)
# - load_tf_models_from_blocks(), collect_model_terms_for_tf_full(),
#   build_peak_gene_summary_for_tf(), attach_ko_truth_to_peak_summary()

# ---- cap_p_for_log10(): Helper for log10 p-values --------------------------
#' Make p-values safe for log10 plotting (avoids log10(0) = -Inf)
cap_p_for_log10 <- function(p, floor = 1e-300) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- NA_real_
  p[p < 0] <- NA_real_
  p[p == 0] <- floor
  p[p > 1] <- 1
  p
}

# ---- Run histogram analysis for TFs of interest ----------------------------
tfs_interest <- c("SOX9") # , "HNF1A"

# ---- TF loop: generate per-TF histograms -----------------------------------
for (tf in tfs_interest) {
  message("Histogram-only: ", tf)

  tf_result <- load_tf_models_from_blocks(cache_dir_models, tf)
  tf_terms  <- collect_model_terms_for_tf_full(tf_result, 0.05, mc_cores = 16)

  if (!nrow(tf_terms)) {
    cli::cli_warn("No terms found for TF {.val {tf}}; skipping.")
    next
  }

  peak_summary <- build_peak_gene_summary_for_tf(
    tf               = tf,
    terms_tf         = tf_terms,
    tf_gene_links_gh = tf_gene_links_gh,
    rna_tbl          = grn_set$rna
  )

  if (!nrow(peak_summary)) {
    cli::cli_warn("Empty peak summary for TF {.val {tf}}; skipping.")
    next
  }

  if (!tf %in% names(ko_truth_list)) {
    cli::cli_warn("No KO truth found in ko_truth_list for TF {.val {tf}}; skipping.")
    next
  }

  tf_with_truth <- attach_ko_truth_to_peak_summary(
    peak_summary      = peak_summary,
    ko_truth_tbl      = ko_truth_list[[tf]],
    keep_intermediate = TRUE
  )

  if (!nrow(tf_with_truth)) {
    cli::cli_warn("No KO-labelled rows for TF {.val {tf}}; skipping.")
    next
  }

  # Per-gene min p-value across peaks (prevents multi-peak genes dominating hist)
  p_tbl <- tf_with_truth |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      p_sem = suppressWarnings(min(sem_direct_tf_expr_p, na.rm = TRUE)),
      p_gam = suppressWarnings(min(gam_tf_expr_p, na.rm = TRUE)),
      .groups = "drop"
    )

  # Some models return extremely small p-values that can underflow to 0;
  # floor them so scale_x_log10() doesn't drop them.
  # For visualization it's usually better to not span 300 orders of magnitude;
  # this is *only* for plotting, not for inference.
  p_floor_hist <- 1e-50
  # Denominator we want to report: N unique TF-gene pairs (here: unique genes for this TF)
  n_pairs <- nrow(p_tbl)
  n0_sem <- sum(is.finite(p_tbl$p_sem) & p_tbl$p_sem == 0, na.rm = TRUE)
  n0_gam <- sum(is.finite(p_tbl$p_gam) & p_tbl$p_gam == 0, na.rm = TRUE)
  if (n0_sem > 0L || n0_gam > 0L) {
    message(
      sprintf(
        "  %s: flooring p==0 for log10 hist (SEM=%d, GAM=%d) to %g",
        tf, n0_sem, n0_gam, p_floor_hist
      )
    )
  }
  # Clean p-values once, then derive plot stats and -log10(p)
  p_sem_clean <- cap_p_for_log10(p_tbl$p_sem, floor = p_floor_hist)
  p_gam_clean <- cap_p_for_log10(p_tbl$p_gam, floor = p_floor_hist)
  n_sem_plot <- sum(is.finite(p_sem_clean), na.rm = TRUE)
  n_gam_plot <- sum(is.finite(p_gam_clean), na.rm = TRUE)
  p_tbl <- p_tbl |>
    dplyr::mutate(
      p_sem = p_sem_clean,
      p_gam = p_gam_clean,
      logp_sem = -log10(p_sem),
      logp_gam = -log10(p_gam)
    )

  # SEM histogram
  x_max_logp <- 50
  p_sem <- ggplot2::ggplot(
    p_tbl[is.finite(p_tbl$logp_sem), , drop = FALSE],
    ggplot2::aes(x = logp_sem)
  ) +
    ggplot2::geom_histogram(bins = 60, colour = "black", fill = "grey80") +
    ggplot2::coord_cartesian(xlim = c(0, x_max_logp)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::annotate(
      geom = "label",
      x = x_max_logp,
      y = Inf,
      hjust = 1.02,
      vjust = 1.2,
      label = sprintf(
        "TF-gene pairs: %d\nplotted: %d\np==0: %d",
        n_pairs, n_sem_plot, n0_sem
      ),
      size = 3
    ) +
    ggplot2::labs(
      title = sprintf("%s RNA SEM direct p-values (per-gene min)", tf),
      subtitle = sprintf("x = -log10(p); p==0 (before flooring): %d; floor used for plotting: %g", n0_sem, p_floor_hist),
      x = expression(-log[10](p)),
      y = "Number of genes"
    ) +
    ggplot2::theme_bw(base_size = 10)

  ggplot2::ggsave(
    filename = file.path(ko_dir, sprintf("%s_RNA_SEM_direct_pvalue_hist.pdf", tf)),
    plot     = p_sem,
    width    = 7,
    height   = 4,
    units    = "in",
    dpi      = 300
  )

  # GAM histogram
  p_gam <- ggplot2::ggplot(
    p_tbl[is.finite(p_tbl$logp_gam), , drop = FALSE],
    ggplot2::aes(x = logp_gam)
  ) +
    ggplot2::geom_histogram(bins = 60, colour = "black", fill = "grey80") +
    ggplot2::coord_cartesian(xlim = c(0, x_max_logp)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::annotate(
      geom = "label",
      x = x_max_logp,
      y = Inf,
      hjust = 1.02,
      vjust = 1.2,
      label = sprintf(
        "TF-gene pairs: %d\nplotted: %d\np==0: %d",
        n_pairs, n_gam_plot, n0_gam
      ),
      size = 3
    ) +
    ggplot2::labs(
      title = sprintf("%s RNA GAM p-values (per-gene min)", tf),
      subtitle = sprintf("x = -log10(p); p==0 (before flooring): %d; floor used for plotting: %g", n0_gam, p_floor_hist),
      x = expression(-log[10](p)),
      y = "Number of genes"
    ) +
    ggplot2::theme_bw(base_size = 10)

  ggplot2::ggsave(
    filename = file.path(ko_dir, sprintf("%s_RNA_GAM_pvalue_hist.pdf", tf)),
    plot     = p_gam,
    width    = 7,
    height   = 4,
    units    = "in",
    dpi      = 300
  )

  # --- FP-only histograms (interactive view; no PDF save) ---
  # Uses fp-related model p-values per gene.
  p_tbl_fp <- tf_with_truth |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      p_sem_fp = suppressWarnings(min(sem_direct_tf_fp_p, na.rm = TRUE)),
      p_gam_fp = suppressWarnings(min(gam_tf_fp_p, na.rm = TRUE)),
      .groups = "drop"
    )

  n_pairs_fp <- nrow(p_tbl_fp)
  n0_sem_fp <- sum(is.finite(p_tbl_fp$p_sem_fp) & p_tbl_fp$p_sem_fp == 0, na.rm = TRUE)
  n0_gam_fp <- sum(is.finite(p_tbl_fp$p_gam_fp) & p_tbl_fp$p_gam_fp == 0, na.rm = TRUE)

  p_sem_fp_clean <- cap_p_for_log10(p_tbl_fp$p_sem_fp, floor = p_floor_hist)
  p_gam_fp_clean <- cap_p_for_log10(p_tbl_fp$p_gam_fp, floor = p_floor_hist)
  n_sem_fp_plot <- sum(is.finite(p_sem_fp_clean), na.rm = TRUE)
  n_gam_fp_plot <- sum(is.finite(p_gam_fp_clean), na.rm = TRUE)

  p_tbl_fp <- p_tbl_fp |>
    dplyr::mutate(
      p_sem_fp = p_sem_fp_clean,
      p_gam_fp = p_gam_fp_clean,
      logp_sem_fp = -log10(p_sem_fp),
      logp_gam_fp = -log10(p_gam_fp)
    )

  p_sem_fp_plot <- ggplot2::ggplot(
    p_tbl_fp[is.finite(p_tbl_fp$logp_sem_fp), , drop = FALSE],
    ggplot2::aes(x = logp_sem_fp)
  ) +
    ggplot2::geom_histogram(bins = 60, colour = "black", fill = "grey80") +
    ggplot2::coord_cartesian(xlim = c(0, x_max_logp)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::annotate(
      geom = "label",
      x = x_max_logp,
      y = Inf,
      hjust = 1.02,
      vjust = 1.2,
      label = sprintf(
        "TF-gene pairs: %d\nplotted: %d\np==0: %d",
        n_pairs_fp, n_sem_fp_plot, n0_sem_fp
      ),
      size = 3
    ) +
    ggplot2::labs(
      title = sprintf("%s FP SEM direct p-values (per-gene min)", tf),
      subtitle = sprintf("x = -log10(p); p==0 (before flooring): %d; floor used for plotting: %g", n0_sem_fp, p_floor_hist),
      x = expression(-log[10](p)),
      y = "Number of genes"
    ) +
    ggplot2::theme_bw(base_size = 10)

  p_gam_fp_plot <- ggplot2::ggplot(
    p_tbl_fp[is.finite(p_tbl_fp$logp_gam_fp), , drop = FALSE],
    ggplot2::aes(x = logp_gam_fp)
  ) +
    ggplot2::geom_histogram(bins = 60, colour = "black", fill = "grey80") +
    ggplot2::coord_cartesian(xlim = c(0, x_max_logp)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::annotate(
      geom = "label",
      x = x_max_logp,
      y = Inf,
      hjust = 1.02,
      vjust = 1.2,
      label = sprintf(
        "TF-gene pairs: %d\nplotted: %d\np==0: %d",
        n_pairs_fp, n_gam_fp_plot, n0_gam_fp
      ),
      size = 3
    ) +
    ggplot2::labs(
      title = sprintf("%s FP GAM p-values (per-gene min)", tf),
      subtitle = sprintf("x = -log10(p); p==0 (before flooring): %d; floor used for plotting: %g", n0_gam_fp, p_floor_hist),
      x = expression(-log[10](p)),
      y = "Number of genes"
    ) +
    ggplot2::theme_bw(base_size = 10)

  # Print to device (RStudio plot pane / interactive session)
  print(p_sem_fp_plot)
  print(p_gam_fp_plot)
}
