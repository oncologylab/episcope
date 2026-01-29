library(DBI)
library(RSQLite)

# prep perturb data -------------------------------------------------------
# helper: safely pick first non-empty value
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

# helper: given TF name, return tibble of DE rows for that TF
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

tf_perturb_db <- "/data/homes/yl814/episcope/tf_perturb.db"
# tf_perturb_db <- "Z:/episcope/tf_perturb.db"

HNF1A_KO <- get_tf_perturbation_tbl(tf_perturb_db, "HNF1A")


# library(episcope)

base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

# base_dir <- "Z:/episcope_test/benchmark_tf_to_target_genes_prediction"
HNF4A_KO <- readr::read_csv("/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/Mayo 5289 siHNF4A RNA-seq.csv")
print(HNF4A_KO)
HNF4A_KO$log2fc  <- HNF4A_KO$log2FoldChange
HNF4A_KO$p_value <- HNF4A_KO$pvalue
HNF4A_KO$p_adj   <- HNF4A_KO$padj
HNF1A_KO
#

source("R/utils_ggvolcano.R")
plot_tf_volcano <- function(df,
                            tf_label,
                            base_dir,
                            symbol_col = c("symbol", "gene_symbol"),
                            logFC_col  = "log2fc",
                            pval_col   = "p_adj",
                            stub_suffix = "KO") {
  symbol_col <- symbol_col[symbol_col %in% names(df)][1]
  if (is.na(symbol_col)) {
    stop("Could not find a symbol column in df (tried: ",
         paste(symbol_col, collapse = ", "), ").")
  }
  if (!logFC_col %in% names(df)) {
    stop("Column '", logFC_col, "' not found in df.")
  }
  if (!pval_col %in% names(df)) {
    stop("Column '", pval_col, "' not found in df.")
  }

  # drop NAs in p_adj
  volcano_data <- df[!is.na(df[[pval_col]]), , drop = FALSE]

  p <- ggvolcano(
    data      = volcano_data,
    labels    = volcano_data[[symbol_col]],  # gene symbols for labeling
    logFC_col = logFC_col,
    pval_col  = pval_col,

    # axis labels & title
    xlab   = bquote(~Log[2]~"fold change ("*.(tf_label)*" KD/KO vs Ctrl)"),
    ylab   = bquote(~-Log[10]~"adjusted"~italic(P)),
    title  = paste0(tf_label, " perturbation RNA-seq: Volcano plot"),
    caption = paste("Total genes:", nrow(volcano_data)),

    # significance cutoffs
    pval_cutoff  = 1e-6,
    logFC_cutoff = 1,

    # point aesthetics
    point_aes = list(
      size  = 0.4,
      shape = c(19, 19, 19, 19),
      color = c("grey70", "#00CD6C", "#009ADE", "#FF1F5B"),
      alpha = 0.7
    ),

    legend_aes = list(
      labels = c(
        "NS",
        expression("|log"[2]*"FC| >= 1"),
        expression("padj <= 1e-6"),
        expression("|log"[2]*"FC| >= 1, padj <= 1e-6")
      ),
      position   = "bottom",
      label_size = 10,
      icon_size  = 4
    ),
    use_significance = TRUE,
    jitter = FALSE
  )

  stub <- paste0(tf_label, "_", stub_suffix)
  outfile <- file.path(
    base_dir,
    paste0(stub, "_RNAseq_volcano_padj.pdf")
  )
  ggplot2::ggsave(outfile, p, width = 8, height = 6, units = "in", dpi = 300)

  message("Saved volcano for ", tf_label, " to:\n  ", outfile)
  invisible(p)
}
p_HNF4A <- plot_tf_volcano(
  df        = HNF4A_KO,
  tf_label  = "HNF4A",
  base_dir  = base_dir,
  symbol_col = c("symbol", "gene_symbol"),
  logFC_col  = "log2fc",
  pval_col   = "p_adj",
  stub_suffix = "KO"
)
p_HNF1A <- plot_tf_volcano(
  df        = HNF1A_KO,
  tf_label  = "HNF1A",
  base_dir  = base_dir,
  symbol_col = c("gene_symbol", "symbol"),
  logFC_col  = "log2fc",
  pval_col   = "p_adj",
  stub_suffix = "KO"
)




db <- "jaspar2024"
base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"
threshold_gene_expr <- 4
threshold_tf_expr   <- 10

# -------------------------------------------------------------------
# Sample metadata and RNA setup
# -------------------------------------------------------------------
sample_metadata <- readxl::read_excel(file.path(base_dir, "sample_metadata.xlsx"), na = "NA")
stric_metadata  <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))

# motif DB + TF list
if (db == "jaspar2024") {
  motif_db <- readr::read_tsv(
    system.file("extdata", "genome", "JASPAR2024.txt", package = "episcope"),
    show_col_types = FALSE
  )
} else if (db == "hocomocov13") {
  motif_db <- readr::read_tsv(
    system.file("extdata", "genome", "HOCOMOCOv13.txt", package = "episcope"),
    show_col_types = FALSE
  )
} else {
  stop("Unsupported db: ", db)
}

tf_list <- motif_db |>
  tidyr::separate_rows(HGNC, sep = "::") |>
  dplyr::filter(!is.na(HGNC), HGNC != "") |>
  dplyr::distinct(HGNC) |>
  dplyr::pull(HGNC)

# RNA
rna <- readr::read_csv(
  file.path(base_dir, "HPAFII_AsPC1_Panc1_combined_smallestGroupSize_3_reads_5_filtered_DESeq2_median_of_ratios_normalized.csv")
)

rna <- episcope::clean_hgnc(rna)
rna <- episcope::filter_rna_expr(
  rna, tf_list,
  hgnc_col   = "HGNC",
  gene_min   = threshold_gene_expr,
  tf_min     = threshold_tf_expr,
  min_samples = 1L
)

# strict_rna: rename columns to ATAC IDs
smap <- dplyr::transmute(stric_metadata, old = strict_match_rna, new = id)
strict_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
nm <- names(strict_rna); nm[match(smap$old, nm)] <- smap$new
strict_rna <- strict_rna |> `names<-`(nm) |> dplyr::as_tibble()

filter_tf_gene_links_by_atac_overlap <- function(tf_gene_links,
                                                 atac_peaks,
                                                 atac_peak_col   = "atac_peak",
                                                 atac_chr_col    = "X1",
                                                 atac_start_col  = "X2",
                                                 atac_end_col    = "X3",
                                                 minoverlap      = 1L) {
  # Requires: GenomicRanges, IRanges, S4Vectors, dplyr, tidyr, rlang

  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    stop(
      "This function requires the Bioconductor packages 'GenomicRanges', 'IRanges', and 'S4Vectors'.",
      call. = FALSE
    )
  }

  # 0) Basic checks
  if (!atac_peak_col %in% colnames(tf_gene_links)) {
    stop(
      sprintf("Column '%s' not found in tf_gene_links.", atac_peak_col),
      call. = FALSE
    )
  }
  if (!all(c(atac_chr_col, atac_start_col, atac_end_col) %in% colnames(atac_peaks))) {
    stop(
      "ATAC peaks must contain columns '", atac_chr_col, "', '",
      atac_start_col, "', and '", atac_end_col, "'.",
      call. = FALSE
    )
  }

  # 1) Unique atac_peak strings from tf_gene_links
  atac_peak_sym <- rlang::sym(atac_peak_col)

  peak_tbl <- tf_gene_links |>
    dplyr::distinct(!!atac_peak_sym) |>
    dplyr::mutate(tf_peak_id = dplyr::row_number())

  if (nrow(peak_tbl) == 0L) {
    return(tf_gene_links[0, , drop = FALSE])
  }

  # 2) Parse "chr:start-end" -> chr, start, end for overlaps
  coord_tbl <- peak_tbl |>
    tidyr::separate(
      col     = !!atac_peak_sym,
      into    = c("chr", "start", "end"),
      sep     = "[:-]",        # split on ":" or "-"
      convert = TRUE
    )

  # 3) Build GRanges for tf_gene_links ATAC peaks
  gr1 <- GenomicRanges::GRanges(
    seqnames = coord_tbl$chr,
    ranges   = IRanges::IRanges(
      start = coord_tbl$start,
      end   = coord_tbl$end
    ),
    tf_peak_id = coord_tbl$tf_peak_id
  )

  # 4) Build GRanges for condition-specific ATAC peaks
  atac_chr   <- atac_peaks[[atac_chr_col]]
  atac_start <- atac_peaks[[atac_start_col]]
  atac_end   <- atac_peaks[[atac_end_col]]

  # drop malformed rows if any
  good <- !(is.na(atac_chr) | is.na(atac_start) | is.na(atac_end))
  atac_chr   <- atac_chr[good]
  atac_start <- atac_start[good]
  atac_end   <- atac_end[good]

  if (length(atac_chr) == 0L) {
    return(tf_gene_links[0, , drop = FALSE])
  }

  gr2 <- GenomicRanges::GRanges(
    seqnames = atac_chr,
    ranges   = IRanges::IRanges(
      start = atac_start,
      end   = atac_end
    )
  )

  # 5) Overlaps
  hits <- GenomicRanges::findOverlaps(
    gr1,
    gr2,
    minoverlap = as.integer(minoverlap)
  )

  if (length(hits) == 0L) {
    return(tf_gene_links[0, , drop = FALSE])
  }

  # 6) Which tf_gene_links peaks overlap?
  keep_ids <- unique(gr1$tf_peak_id[S4Vectors::queryHits(hits)])

  # Map back to original atac_peak strings
  keep_peaks <- peak_tbl |>
    dplyr::filter(tf_peak_id %in% keep_ids) |>
    dplyr::pull(!!atac_peak_sym)

  # 7) Filter tf_gene_links to those peaks
  tf_gene_links[tf_gene_links[[atac_peak_col]] %in% keep_peaks, , drop = FALSE]
}



aspc_10fbs_atac <- readr::read_tsv(file.path(base_dir, "cy414.hg38.rp10m.narrowpeaks.filtered.bed"), col_names = FALSE)
# Shared: read correlation tables once (TF-independent)
tf_gene_links  <- readr::read_csv(file.path(base_dir, "fp_gene_corr_full_jaspar2024.csv"))
atac_gene_links <- readr::read_csv(file.path(base_dir, "atac_gene_corr_full_jaspar2024.csv"))

tf_gene_links_atac_filtered <- filter_tf_gene_links_by_atac_overlap(
  tf_gene_links  = tf_gene_links,
  atac_peaks     = aspc_10fbs_atac,
  atac_peak_col  = "atac_peak",  # column in tf_gene_links
  atac_chr_col   = "X1",         # chrom column in aspc_10fbs_atac
  atac_start_col = "X2",         # start column
  atac_end_col   = "X3",         # end column
  minoverlap     = 1L            # >=1 bp overlap
)

tf_gene_links <- tf_gene_links_atac_filtered

# map motif IDs -> HGNC name(s)
tf_gene_links$motifs <- motif_db$HGNC[match(tf_gene_links$motifs, motif_db$motif)]

filter_tf_gene_links_by_fp_corr <- function(tf_links,
                                            tf_tfbs,
                                            tf,
                                            threshold_r_fp     = 0.3,
                                            threshold_p_adj_fp = 0.001,
                                            seed               = NULL) {
  # Early exit if no rows
  if (nrow(tf_links) == 0L) {
    empty_tbl <- tibble::tibble(gene_key = character(0))
    return(list(
      canonical = list(
        passing     = empty_tbl,
        non_passing = empty_tbl,
        random      = empty_tbl
      ),
      all = list(
        passing     = empty_tbl,
        non_passing = empty_tbl,
        random      = empty_tbl
      )
    ))
  }

  # Background (for random genes): peaks NOT in tf_tfbs, *unfiltered* by motifs
  tf_links_bg <- tf_links[!(tf_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) set.seed(seed)

  # Helper to compute passing / non-passing given a "current" subset
  .split_pass_nonpass <- function(tf_links_current,
                                  threshold_r_fp,
                                  threshold_p_adj_fp) {
    if (nrow(tf_links_current) == 0L) {
      return(list(
        passing     = tibble::tibble(gene_key = character(0)),
        non_passing = tibble::tibble(gene_key = character(0))
      ))
    }

    pass_rows <- !is.na(tf_links_current$r_fp) &
      !is.na(tf_links_current$p_adj_fp) &
      abs(tf_links_current$r_fp) > threshold_r_fp &
      tf_links_current$p_adj_fp   < threshold_p_adj_fp

    genes_all  <- sort(unique(tf_links_current$gene_key))
    genes_pass <- sort(unique(tf_links_current$gene_key[pass_rows]))
    genes_fail <- setdiff(genes_all, genes_pass)

    list(
      passing     = tibble::tibble(gene_key = genes_pass),
      non_passing = tibble::tibble(gene_key = genes_fail)
    )
  }

  # ---------- ALL mode ----------
  tf_links_all <- tf_links[tf_links$fp_peak %in% tf_tfbs, , drop = FALSE]
  split_all <- .split_pass_nonpass(
    tf_links_current   = tf_links_all,
    threshold_r_fp     = threshold_r_fp,
    threshold_p_adj_fp = threshold_p_adj_fp
  )

  n_target_all <- nrow(split_all$passing)
  n_draw_all   <- min(n_target_all, length(bg_genes))
  genes_rand_all <- if (n_draw_all > 0L) {
    sort(sample(bg_genes, n_draw_all, replace = FALSE))
  } else {
    character(0)
  }

  # ---------- CANONICAL mode ----------
  # Keep rows where motifs contains the TF name (case-insensitive)
  tf_upper   <- toupper(tf)
  motifs_vec <- toupper(as.character(tf_links$motifs))
  keep_motif <- !is.na(motifs_vec) & grepl(tf_upper, motifs_vec, fixed = TRUE)

  tf_links_canon <- tf_links[
    tf_links$fp_peak %in% tf_tfbs & keep_motif,
    ,
    drop = FALSE
  ]

  split_canon <- .split_pass_nonpass(
    tf_links_current   = tf_links_canon,
    threshold_r_fp     = threshold_r_fp,
    threshold_p_adj_fp = threshold_p_adj_fp
  )

  n_target_canon <- nrow(split_canon$passing)
  n_draw_canon   <- min(n_target_canon, length(bg_genes))
  genes_rand_canon <- if (n_draw_canon > 0L) {
    sort(sample(bg_genes, n_draw_canon, replace = FALSE))
  } else {
    character(0)
  }

  # Assemble output
  list(
    canonical = list(
      passing     = split_canon$passing,
      non_passing = split_canon$non_passing,
      random      = tibble::tibble(gene_key = genes_rand_canon)
    ),
    all = list(
      passing     = split_all$passing,
      non_passing = split_all$non_passing,
      random      = tibble::tibble(gene_key = genes_rand_all)
    )
  )
}

build_tf_gene_link_sets_with_atac <- function(
    tf,
    tf_tfbs,
    tf_gene_links,
    atac_gene_links,
    threshold_r_fp       = 0.3,
    threshold_p_adj_fp   = 0.01,
    threshold_r_atac     = 0.3,
    threshold_p_adj_atac = 0.05,
    seed                 = 1L
) {
  # --- helper to annotate each tibble with label + meta --------------------
  annotate_tbl <- function(tbl,
                           tf,
                           mode,          # "canonical" or "all"
                           set,           # "passing", "non_passing", "random"
                           atac_filtered, # TRUE/FALSE
                           threshold_r_fp,
                           threshold_p_adj_fp,
                           threshold_r_atac,
                           threshold_p_adj_atac) {

    meta <- list(
      tf                   = tf,
      mode                 = mode,
      set                  = set,
      atac_filtered        = atac_filtered,
      threshold_r_fp       = threshold_r_fp,
      threshold_p_adj_fp   = threshold_p_adj_fp,
      threshold_r_atac     = if (atac_filtered) threshold_r_atac     else NA_real_,
      threshold_p_adj_atac = if (atac_filtered) threshold_p_adj_atac else NA_real_
    )

    lab <- sprintf(
      "TF=%s; mode=%s; set=%s; atac_filtered=%s; |r_fp|>%g; p_adj_fp<%g%s",
      tf,
      mode,
      set,
      atac_filtered,
      threshold_r_fp,
      threshold_p_adj_fp,
      if (atac_filtered) {
        sprintf("; |r_atac|>%g; p_adj_atac<%g",
                threshold_r_atac, threshold_p_adj_atac)
      } else {
        ""
      }
    )

    attr(tbl, "label") <- lab
    attr(tbl, "meta")  <- meta
    tbl
  }

  if (!is.null(seed)) set.seed(seed)

  # ---------------- FP-only (unfiltered by ATAC) ---------------------------
  res_unfiltered <- filter_tf_gene_links_by_fp_corr(
    tf_links           = tf_gene_links,
    tf_tfbs            = tf_tfbs,
    tf                 = tf,
    threshold_r_fp     = threshold_r_fp,
    threshold_p_adj_fp = threshold_p_adj_fp,
    seed               = seed
  )

  # ---------------- Build ATAC-filtered tf_gene_links ----------------------
  atac_pass <- atac_gene_links[
    !is.na(atac_gene_links$r_atac) &
      !is.na(atac_gene_links$p_adj_atac) &
      abs(atac_gene_links$r_atac) > threshold_r_atac &
      atac_gene_links$p_adj_atac < threshold_p_adj_atac,
    ,
    drop = FALSE
  ]

  atac_peaks_pass <- unique(atac_pass$atac_peak)

  tf_gene_links_atac <- tf_gene_links[
    tf_gene_links$atac_peak %in% atac_peaks_pass,
    ,
    drop = FALSE
  ]

  # ---------------- FP + ATAC-filtered version -----------------------------
  res_atac <- filter_tf_gene_links_by_fp_corr(
    tf_links           = tf_gene_links_atac,
    tf_tfbs            = tf_tfbs,
    tf                 = tf,
    threshold_r_fp     = threshold_r_fp,
    threshold_p_adj_fp = threshold_p_adj_fp,
    seed               = seed
  )

  # ---------------- Annotate all 12 tables ---------------------------------
  # Unfiltered
  uf_can_pass     <- annotate_tbl(res_unfiltered$canonical$passing,
                                  tf, "canonical", "passing", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  uf_can_nonpass  <- annotate_tbl(res_unfiltered$canonical$non_passing,
                                  tf, "canonical", "non_passing", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  uf_can_rand     <- annotate_tbl(res_unfiltered$canonical$random,
                                  tf, "canonical", "random", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)

  uf_all_pass     <- annotate_tbl(res_unfiltered$all$passing,
                                  tf, "all", "passing", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  uf_all_nonpass  <- annotate_tbl(res_unfiltered$all$non_passing,
                                  tf, "all", "non_passing", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  uf_all_rand     <- annotate_tbl(res_unfiltered$all$random,
                                  tf, "all", "random", FALSE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)

  # ATAC-filtered
  af_can_pass     <- annotate_tbl(res_atac$canonical$passing,
                                  tf, "canonical", "passing", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  af_can_nonpass  <- annotate_tbl(res_atac$canonical$non_passing,
                                  tf, "canonical", "non_passing", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  af_can_rand     <- annotate_tbl(res_atac$canonical$random,
                                  tf, "canonical", "random", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)

  af_all_pass     <- annotate_tbl(res_atac$all$passing,
                                  tf, "all", "passing", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  af_all_nonpass  <- annotate_tbl(res_atac$all$non_passing,
                                  tf, "all", "non_passing", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)
  af_all_rand     <- annotate_tbl(res_atac$all$random,
                                  tf, "all", "random", TRUE,
                                  threshold_r_fp, threshold_p_adj_fp,
                                  threshold_r_atac, threshold_p_adj_atac)

  list(
    unfiltered = list(
      canonical = list(
        passing     = uf_can_pass,
        non_passing = uf_can_nonpass,
        random      = uf_can_rand
      ),
      all = list(
        passing     = uf_all_pass,
        non_passing = uf_all_nonpass,
        random      = uf_all_rand
      )
    ),
    atac_filtered = list(
      canonical = list(
        passing     = af_can_pass,
        non_passing = af_can_nonpass,
        random      = af_can_rand
      ),
      all = list(
        passing     = af_all_pass,
        non_passing = af_all_nonpass,
        random      = af_all_rand
      )
    )
  )
}

get_tf_tfbs <- function(tf,
                        tfbs_dir      = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs",
                        tfbs_r_cut    = 0.3,
                        tfbs_padj_cut = 0.05) {

  overview_file <- file.path(tfbs_dir, sprintf("%s_overview.txt", tf))
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

make_tf_sets <- function(tf,
                         tfbs_dir             = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs",
                         tfbs_r_cut           = 0.3,
                         tfbs_padj_cut        = 0.05,
                         threshold_r_fp       = 0.3,
                         threshold_p_adj_fp   = 0.01,
                         threshold_r_atac     = 0.3,
                         threshold_p_adj_atac = 0.05,
                         seed                 = 1L) {

  # uses global tf_gene_links and atac_gene_links

  tf_tfbs <- get_tf_tfbs(
    tf           = tf,
    tfbs_dir     = tfbs_dir,
    tfbs_r_cut   = tfbs_r_cut,
    tfbs_padj_cut = tfbs_padj_cut
  )

  build_tf_gene_link_sets_with_atac(
    tf                   = tf,
    tf_tfbs              = tf_tfbs,
    tf_gene_links        = tf_gene_links,
    atac_gene_links      = atac_gene_links,
    threshold_r_fp       = threshold_r_fp,
    threshold_p_adj_fp   = threshold_p_adj_fp,
    threshold_r_atac     = threshold_r_atac,
    threshold_p_adj_atac = threshold_p_adj_atac,
    seed                 = seed
  )
}




sets_HNF1A <- make_tf_sets("HNF1A")
sets_HNF4A <- make_tf_sets("HNF4A")


attr(sets_HNF1A$unfiltered$canonical$passing, "label")
attr(sets_HNF4A$atac_filtered$all$random, "label")

tf_tfbs_HNF1A <- get_tf_tfbs("HNF1A", tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")
tf_tfbs_HNF4A <- get_tf_tfbs("HNF4A", tfbs_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs")


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(plyr)
  library(patchwork)
})
ko_dir  <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

# -------------------------------------------------------------------
# geom_split_violin helper (from introdataviz)
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
# HNF1A KO DEA: prepare KO table
# -------------------------------------------------------------------
ko_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"

# HNF1A_KO is already in memory (17k rows, with log2fc + gene_symbol)
ko_tbl <- HNF1A_KO

if ("gene_symbol" %in% names(ko_tbl)) {
  ko_tbl <- dplyr::transmute(
    ko_tbl,
    gene   = gene_symbol,
    log2FC = log2fc
  )
} else if ("symbol" %in% names(ko_tbl)) {
  ko_tbl <- dplyr::transmute(
    ko_tbl,
    gene   = symbol,
    log2FC = log2fc
  )
} else {
  stop("Cannot find gene_symbol or symbol column in HNF1A_KO.")
}

ko_tbl <- ko_tbl[!is.na(ko_tbl$gene), , drop = FALSE]
tf <- "HNF1A"

# -------------------------------------------------------------------
# Helper: grid plot for one (ATAC filter, mode) combination
# -------------------------------------------------------------------
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

    # *** IMPORTANT FIX: drop unused factor levels ***
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
  ko_tbl <- ko_tbl[ko_tbl$gene %in% genes_univ, , drop = FALSE]
  n_univ <- length(genes_univ)
  if (n_univ == 0L) {
    stop("No overlapping genes between KO table and tf_gene_links.")
  }

  # Background gene pool for random (always unfiltered, excluding TFBS peaks)
  tf_links_bg <- tf_gene_links[!(tf_gene_links$fp_peak %in% tf_tfbs), , drop = FALSE]
  bg_genes <- sort(unique(tf_links_bg$gene_key))
  bg_genes <- bg_genes[!is.na(bg_genes) & nzchar(bg_genes)]

  if (!is.null(seed)) {
    set.seed(seed)
  }

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
      abs(tf_links_curr$r_fp) > r_cut &
      tf_links_curr$p_adj_fp   < p_cut

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

    ko_pass <- ko_tbl[ko_tbl$gene %in% pass_i, , drop = FALSE]
    ko_non  <- ko_tbl[ko_tbl$gene %in% non_i, , drop = FALSE]
    ko_rand <- ko_tbl[ko_tbl$gene %in% rand_i, , drop = FALSE]

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

  # Factor labels for facets and x-axis (now r_fp / p_fp)
  r_fp_levels <- r_fp_cuts
  p_fp_levels <- sort(unique(p_fp_cuts))

  plot_df <- plot_df %>%
    dplyr::mutate(
      r_lab = factor(
        r_fp_cut,
        levels = r_fp_levels,
        labels = paste0("|r_fp| > ", r_fp_levels)
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
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
    dplyr::group_by(r_fp_cut, p_fp_cut) %>%
    dplyr::summarise(
      y = max(log2FC, na.rm = TRUE),
      .groups = "drop"
    )

  # vs Random
  df_rand_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Random_background"),
                         , drop = FALSE]

  pval_rand <- df_rand_all %>%
    dplyr::group_by(r_fp_cut, p_fp_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_rand <- dplyr::left_join(pval_rand, y_pos,
                                 by = c("r_fp_cut", "p_fp_cut")) %>%
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
        r_fp_cut,
        levels = r_fp_levels,
        labels = paste0("|r_fp| > ", r_fp_levels)
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
      )
    )

  # vs Non-regulated
  df_non_all <- plot_df[plot_df$group %in% c("Predicted_regulated", "Predicted_nonregulated"),
                        , drop = FALSE]

  pval_non <- df_non_all %>%
    dplyr::group_by(r_fp_cut, p_fp_cut) %>%
    dplyr::summarise(
      p_val = run_test(dplyr::cur_data_all(), test = test, alternative = alternative),
      .groups = "drop"
    )

  annot_non <- dplyr::left_join(pval_non, y_pos,
                                by = c("r_fp_cut", "p_fp_cut")) %>%
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
        r_fp_cut,
        levels = r_fp_levels,
        labels = paste0("|r_fp| > ", r_fp_levels)
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
      )
    )

  # -------------------------------------------------------------------
  # Diagnostics
  # -------------------------------------------------------------------
  if (verbose) {
    atac_lab <- if (filter_atac) "ATAC-filtered" else "unfiltered"
    message("=== Diagnostics for TF = ", tf,
            " (mode = ", mode, ", ", atac_lab, ") ===")
    message("Universe size (n_univ): ", n_univ)
    message("Test: ", test, ", alternative = ", alternative)

    # group counts per facet (first few rows)
    message("Random vs predicted: group counts per facet (head):")
    counts_rand <- df_rand_all %>%
      dplyr::group_by(r_fp_cut, p_fp_cut, group) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop")
    print(utils::head(counts_rand, 10))

    message("Predicted vs non-regulated: group counts per facet (head):")
    counts_non <- df_non_all %>%
      dplyr::group_by(r_fp_cut, p_fp_cut, group) %>%
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
      ord_r <- order(pval_rand$p_val, na.last = NA)
      message("Top 5 smallest p (random vs predicted):")
      print(utils::head(pval_rand[ord_r, , drop = FALSE], 5))
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
      ord_n <- order(pval_non$p_val, na.last = NA)
      message("Top 5 smallest p (predicted vs non-regulated):")
      print(utils::head(pval_non[ord_n, , drop = FALSE], 5))
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
    scale_x_discrete(name = "p_adj_fp cutoff") +
    scale_y_continuous(name = sprintf("log2FC (%s KO vs Ctrl)", tf)) +
    ggtitle(
      sprintf(
        "%s knockout: predicted regulated vs random background genes\n(mode = %s, ATAC filter = %s; n = %d genes in universe)",
        tf,
        mode,
        if (filter_atac) "yes" else "no",
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
    scale_x_discrete(name = "p_adj_fp cutoff") +
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
    dplyr::group_by(r_fp_cut, p_fp_cut, group) %>%
    dplyr::summarise(
      n_genes = dplyr::n_distinct(gene),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      r_lab = factor(
        r_fp_cut,
        levels = r_fp_levels,
        labels = paste0("|r_fp| > ", r_fp_levels)
      ),
      p_lab = factor(
        p_fp_cut,
        levels = p_fp_levels,
        labels = scales::label_scientific(digits = 1)(p_fp_levels)
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
    scale_x_discrete(name = "p_adj_fp cutoff") +
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

  atac_lab <- if (filter_atac) "atac_filtered" else "unfiltered"

  outfile <- file.path(
    out_dir,
    sprintf(
      "%s_KO_%s_%s_fpCorr_grid_split_violin_vsRand_vsNon.pdf",
      tf,
      atac_lab,
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


plot_fp_ko_grid(
  tf              = "HNF1A",
  ko_tbl          = ko_tbl,
  tf_gene_links   = tf_gene_links,
  atac_gene_links = atac_gene_links,
  tf_tfbs         = tf_tfbs_HNF1A,
  mode            = "all",
  filter_atac     = FALSE,
  test            = "wilcox",
  alternative     = "two.sided",
  verbose         = TRUE,
  out_dir         = ko_dir
)


# -------------------------------------------------------------------
# Run for HNF1A: 4 combinations (unfiltered/ATAC x canonical/all)
# -------------------------------------------------------------------


for (filter_atac in c(FALSE, TRUE)) {
  for (mode in c("all")) {
    plot_fp_ko_grid(
      tf              = "HNF1A",
      ko_tbl          = ko_tbl,           # HNF1A_KO table you built
      tf_gene_links   = tf_gene_links,
      atac_gene_links = atac_gene_links,
      tf_tfbs         = tf_tfbs_HNF1A,    # <---- IMPORTANT
      mode            = mode,
      filter_atac     = filter_atac,
      r_fp_cuts       = c(0, 0.1, 0.3, 0.5, 0.7),
      p_fp_cuts       = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      threshold_r_atac     = 0.3,
      threshold_p_adj_atac = 0.05,
      seed            = 1L,
      out_dir         = ko_dir
    )
  }
}



## ---------------------------------------------------------------
## 1) Prepare KO table for HNF4A
## ---------------------------------------------------------------
# Assuming you already read this earlier:
# HNF4A_KO <- readr::read_csv(file.path(ko_dir, "Mayo 5289 siHNF4A RNA-seq.csv"))

ko_tbl_HNF4A <- HNF4A_KO %>%
  dplyr::transmute(
    gene   = symbol,         # KO gene column
    log2FC = log2FoldChange  # KO log2FC column
  ) %>%
  dplyr::filter(!is.na(gene))


## your HNF4A TFBS vector (analogous to tf_tfbs for HNF1A)
tf_tfbs_HNF4A <- tf_tfbs_HNF4A  # whatever object you already have


## ---------------------------------------------------------------
## 2) Run the same 4 combinations for HNF4A
## ---------------------------------------------------------------
for (filter_atac in c(FALSE, TRUE)) {
  for (mode in c("all")) {
    plot_fp_ko_grid(
      tf              = "HNF4A",
      ko_tbl          = ko_tbl_HNF4A,
      tf_gene_links   = tf_gene_links,
      atac_gene_links = atac_gene_links,
      tf_tfbs         = tf_tfbs_HNF4A,   # <---- IMPORTANT
      mode            = mode,
      filter_atac     = filter_atac,
      r_fp_cuts       = c(0, 0.1, 0.3, 0.5, 0.7),
      p_fp_cuts       = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5),
      threshold_r_atac     = 0.3,
      threshold_p_adj_atac = 0.05,
      seed            = 1L,
      out_dir         = ko_dir
    )
  }
}

