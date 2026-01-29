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
HNF4A_KO <- readr::read_csv(file.path(base_dir, "Mayo 5289 siHNF4A RNA-seq.csv"))
print(HNF4A_KO)
HNF4A_KO$log2fc  <- HNF4A_KO$log2FoldChange
HNF4A_KO$p_value <- HNF4A_KO$pvalue
HNF4A_KO$p_adj   <- HNF4A_KO$padj



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
strict_metadata <- sample_metadata |> dplyr::filter(!is.na(strict_match_rna))

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
smap <- dplyr::transmute(strict_metadata, old = strict_match_rna, new = id)
strict_rna <- rna |> dplyr::select(c("ensembl_gene_id", "HGNC", smap$old))
nm <- names(strict_rna); nm[match(smap$old, nm)] <- smap$new
strict_rna <- strict_rna |> `names<-`(nm) |> dplyr::as_tibble()

# -------------------------------------------------------------------
# Footprint scores + GeneHancer mapping (global)
# -------------------------------------------------------------------
fp_score <- readr::read_csv(
  file.path(base_dir, sprintf("fp_score_strict_tf_filtered_corr_%s.csv", db))
)

gh_std <- episcope::load_genehancer_panc(
  system.file("extdata", "GeneHancer_v5.24_elite_panc.csv", package = "episcope")
)

# keep only GH entries for genes present in strict_rna
gh_std <- gh_std |> dplyr::filter(connected_gene %in% strict_rna$HGNC)

# Cache for footprint↔GeneHancer overlaps
fp_gh_map_file <- file.path(base_dir, sprintf("fp_gh_map_%s.csv", db))

if (file.exists(fp_gh_map_file)) {
  fp_gh_map <- readr::read_csv(fp_gh_map_file, show_col_types = FALSE)
} else {
  # unique footprint BED from peak_ID
  fp_bed <- fp_score |>
    dplyr::distinct(peak_ID) |>
    tidyr::separate(
      peak_ID,
      into = c("chrom", "start", "end"),
      sep  = "[:-]",
      convert = TRUE
    )

  # GH BED with only needed columns
  gh_bed <- gh_std |>
    dplyr::select(
      chrom,
      start,
      end,
      connected_gene
    )

  ov <- valr::bed_intersect(fp_bed, gh_bed)
  ov <- dplyr::filter(ov, .overlap > 0)

  fp_gh_map <- ov |>
    dplyr::transmute(
      fp_chr   = chrom,
      fp_start = start.x,
      fp_end   = end.x,
      gh_chr   = chrom,
      gh_start = start.y,
      gh_end   = end.y,
      connected_gene = connected_gene.y
    ) |>
    dplyr::distinct()

  readr::write_csv(fp_gh_map, fp_gh_map_file)
}

cat("fp_gh_map has", nrow(fp_gh_map), "rows\n")

# -------------------------------------------------------------------
# TF-specific part: choose TF and predicted TFBS
# -------------------------------------------------------------------
tf <- "HNF4A"

predicted_tfbs <- readr::read_tsv(
  "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/predicted_all_tfbs/HNF4A_overview.txt",
  show_col_types = FALSE
)

## Filter predicted_tfbs by corr_fp_tf_r > 0.3 & corr_fp_tf_p_adj < 0.05
predicted_tfbs_sig <- predicted_tfbs[
  !is.na(predicted_tfbs$corr_fp_tf_r) &
    !is.na(predicted_tfbs$corr_fp_tf_p_adj) &
    predicted_tfbs$corr_fp_tf_r > 0.3 &
    predicted_tfbs$corr_fp_tf_p_adj < 0.05,
  ,
  drop = FALSE
]

## Create TFBS ID vector "TFBS_chr:TFBS_start-TFBS_end"
tf_tfbs <- unique(
  paste0(
    predicted_tfbs_sig$TFBS_chr, ":",
    predicted_tfbs_sig$TFBS_start, "-",
    predicted_tfbs_sig$TFBS_end
  )
)

cat("tf_tfbs (", length(tf_tfbs), "entries):\n", sep = "")

# Filter fp_gh_map to current TF peaks → fp_gh_map_current
fp_ids_for_map <- paste0(
  fp_gh_map$fp_chr, ":",
  fp_gh_map$fp_start, "-",
  fp_gh_map$fp_end
)

fp_gh_map_current <- fp_gh_map[
  fp_ids_for_map %in% tf_tfbs,
  ,
  drop = FALSE
]

cat("\nfp_gh_map_current has", nrow(fp_gh_map_current), "rows\n")

# -------------------------------------------------------------------
# Correlate footprint accessibility with connected_gene expression
#   for fp_gh_map_current (with caching + parallel + progress)
# -------------------------------------------------------------------
out_corr_file <- file.path(base_dir, sprintf("fp_gh_map_corr_%s_%s.csv", db, tf))

if (file.exists(out_corr_file)) {
  fp_gh_map_corr <- readr::read_csv(out_corr_file, show_col_types = FALSE)
  if (nrow(fp_gh_map_corr) != nrow(fp_gh_map_current)) {
    warning("Cached ", out_corr_file, " has ", nrow(fp_gh_map_corr),
            " rows but fp_gh_map_current has ", nrow(fp_gh_map_current),
            "; recomputing correlations.")
    fp_gh_map_corr <- NULL
  } else {
    cat("Loaded cached correlations from:\n  ", out_corr_file, "\n", sep = "")
  }
} else {
  fp_gh_map_corr <- NULL
}

if (is.null(fp_gh_map_corr)) {
  # how many workers to use
  workers <- max(1L, parallel::detectCores() - 1L)

  # 1) Build peak IDs and gene keys for fp_gh_map_current
  fp_gh_map_current2 <- fp_gh_map_current |>
    dplyr::mutate(
      peak_ID  = paste0(fp_chr, ":", fp_start, "-", fp_end),
      gene_key = connected_gene
    )

  # 2) Work with unique peak_ID rows in fp_score (there are duplicates)
  fp_score_unique <- fp_score[!duplicated(fp_score$peak_ID), , drop = FALSE]

  # 3) Determine shared sample IDs between fp_score and strict_rna
  sample_ids <- intersect(names(fp_score_unique), names(strict_rna))
  sample_ids <- setdiff(sample_ids, c("peak_ID", "ensembl_gene_id", "HGNC"))

  if (length(sample_ids) == 0L) {
    stop("No overlapping sample IDs between fp_score and strict_rna.")
  }

  # 4) Long format: footprints
  fp_long <- fp_score_unique |>
    dplyr::semi_join(
      fp_gh_map_current2 |> dplyr::distinct(peak_ID),
      by = "peak_ID"
    ) |>
    dplyr::select(peak_ID, dplyr::all_of(sample_ids)) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(sample_ids),
      names_to  = "sample",
      values_to = "acc"
    )

  # 5) Long format: RNA for connected genes
  rna_long <- strict_rna |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::select(
      gene_key = HGNC,
      dplyr::all_of(sample_ids)
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(sample_ids),
      names_to  = "sample",
      values_to = "expr"
    ) |>
    dplyr::distinct(gene_key, sample, .keep_all = TRUE)

  # 6) Unique footprint–gene pairs we need correlations for
  pairs_fg <- fp_gh_map_current2 |>
    dplyr::distinct(peak_ID, gene_key)

  if (!nrow(pairs_fg)) {
    stop("fp_gh_map_current has no (peak_ID, gene_key) pairs.")
  }

  # 7) Split into chunks
  n_pairs   <- nrow(pairs_fg)
  workers   <- max(1L, parallel::detectCores() - 1L)
  n_chunks  <- max(1L, min(workers * 4L, n_pairs))
  chunk_size <- max(1L, ceiling(n_pairs / n_chunks))

  idx_split <- split(
    seq_len(n_pairs),
    ceiling(seq_along(seq_len(n_pairs)) / chunk_size)
  )

  message("Computing correlations for ", n_pairs, " peak–gene pairs in ",
          length(idx_split), " chunks (workers = ", workers, ") ...")

  run_chunk <- function(ii) {
    pp <- pairs_fg[ii, , drop = FALSE]

    ax <- dplyr::inner_join(pp, fp_long, by = "peak_ID")
    if (!nrow(ax)) {
      return(tibble::tibble(
        peak_ID = character(),
        gene_key = character(),
        n_gene = integer(),
        r_gene = double(),
        p_gene = double()
      ))
    }

    axx <- dplyr::inner_join(ax, rna_long, by = c("gene_key", "sample"))
    if (!nrow(axx)) {
      return(tibble::tibble(
        peak_ID = character(),
        gene_key = character(),
        n_gene = integer(),
        r_gene = double(),
        p_gene = double()
      ))
    }

    axx |>
      dplyr::group_by(peak_ID, gene_key) |>
      dplyr::summarise(
        n_gene = sum(is.finite(acc) & is.finite(expr)),
        r_gene = suppressWarnings(stats::cor(acc, expr, use = "complete.obs", method = "pearson")),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        p_gene = dplyr::case_when(
          is.na(r_gene) | n_gene < 3 ~ NA_real_,
          abs(r_gene) >= 1          ~ 0,
          TRUE ~ {
            tt <- r_gene * sqrt((n_gene - 2) / (1 - r_gene^2))
            2 * stats::pt(-abs(tt), df = n_gene - 2)
          }
        )
      )
  }

  # 8) Parallel or serial execution with basic progress
  if (workers > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)

    corr_list <- future.apply::future_lapply(
      seq_along(idx_split),
      function(j) {
        message("  [future] chunk ", j, " / ", length(idx_split))
        run_chunk(idx_split[[j]])
      },
      future.seed = TRUE
    )
  } else {
    if (workers > 1L) {
      message("future.apply not available; falling back to serial execution.")
    }
    pb <- utils::txtProgressBar(min = 0, max = length(idx_split), style = 3)
    corr_list <- vector("list", length(idx_split))
    for (j in seq_along(idx_split)) {
      corr_list[[j]] <- run_chunk(idx_split[[j]])
      utils::setTxtProgressBar(pb, j)
    }
    close(pb)
  }

  corr_tbl <- dplyr::bind_rows(corr_list)

  # 9) Global BH across all pairs
  good_p <- is.finite(corr_tbl$p_gene)
  corr_tbl$p_adj_gene <- NA_real_
  if (any(good_p)) {
    corr_tbl$p_adj_gene[good_p] <- stats::p.adjust(corr_tbl$p_gene[good_p], method = "BH")
  }

  # 10) Join correlations back to fp_gh_map_current (same nrow, extra cols)
  fp_gh_map_current2 <- fp_gh_map_current |>
    dplyr::mutate(
      peak_ID  = paste0(fp_chr, ":", fp_start, "-", fp_end),
      gene_key = connected_gene
    )

  fp_gh_map_corr <- fp_gh_map_current2 |>
    dplyr::left_join(corr_tbl, by = c("peak_ID", "gene_key"))

  stopifnot(nrow(fp_gh_map_corr) == nrow(fp_gh_map_current))

  # 11) Save to CSV (per TF)
  readr::write_csv(fp_gh_map_corr, out_corr_file)
  message("Wrote footprint–GeneHancer correlation table to:\n  ", out_corr_file)
}

# -------------------------------------------------------------------
# TF-specific benchmarking bits (now using fresh correlation file)
# -------------------------------------------------------------------

# Use correlation file as tf_gene_links
tf_gene_links <- fp_gh_map_corr

# Filter to this TF’s TFBS (peak_ID)
tf_gene_links_current <- tf_gene_links[
  tf_gene_links$peak_ID %in% tf_tfbs,
  ,
  drop = FALSE
]
cat("\nFiltered tf_gene_links_current has", nrow(tf_gene_links_current), "rows\n")

# Minimal function to split gene_keys into passing vs non-passing
filter_tf_gene_links_by_gene_corr <- function(tf_links,
                                              threshold_r_fp_gene = 0.3,
                                              threshold_p_adj_fp_gene = 0.001) {
  required_cols <- c("gene_key", "r_gene", "p_adj_gene")
  missing_cols <- setdiff(required_cols, names(tf_links))
  if (length(missing_cols) > 0L) {
    stop(
      "tf_links is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(tf_links) == 0L) {
    return(list(
      passing_gene_keys      = character(0),
      non_passing_gene_keys  = character(0)
    ))
  }

  pass_rows <- !is.na(tf_links$r_gene) &
    !is.na(tf_links$p_adj_gene) &
    abs(tf_links$r_gene) > threshold_r_fp_gene &
    tf_links$p_adj_gene < threshold_p_adj_fp_gene

  # UNIQUE gene-level sets
  genes_all  <- sort(unique(tf_links$gene_key))
  genes_pass <- sort(unique(tf_links$gene_key[pass_rows]))
  genes_fail <- setdiff(genes_all, genes_pass)

  list(
    passing_gene_keys      = genes_pass,
    non_passing_gene_keys  = genes_fail
  )
}


res_genes <- filter_tf_gene_links_by_gene_corr(
  tf_gene_links_current,
  threshold_r_fp_gene     = 0.3,
  threshold_p_adj_fp_gene = 0.01
)

cat("\nPassing gene_keys (", length(res_genes$passing_gene_keys), "):\n", sep = "")
cat("\nNon-passing gene_keys (", length(res_genes$non_passing_gene_keys), "):\n", sep = "")
passing_genes    <- res_genes$passing_gene_keys        # unique gene list
non_passing_genes <- res_genes$non_passing_gene_keys   # unique gene list


## ------------------------------------------------------------------
## Helper: split-violin geom
## ------------------------------------------------------------------
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
        transform(data, x = if (grp %% 2 == 1) xminv else xmaxv),
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

## ------------------------------------------------------------------
## Grid of (r, p) cutoffs and split-violin + counts row with Wilcoxon p
## ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(plyr)
  library(patchwork)
})

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
# HNF4A KO DEA (if not already loaded)
# -------------------------------------------------------------------
ko_dir  <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction"
HNF4A_KO <- readr::read_csv(
  file.path(ko_dir, "Mayo 5289 siHNF4A RNA-seq.csv"),
  show_col_types = FALSE
)

ko_tbl <- HNF4A_KO %>%
  dplyr::select(symbol, log2FoldChange) %>%
  dplyr::filter(!is.na(symbol))

# Gene universe = genes in fp_gh_map_corr that also have KO log2FC
genes_univ <- intersect(
  unique(fp_gh_map_corr$gene_key),
  ko_tbl$symbol
)
n_univ <- length(genes_univ)

# r and p grids (your requested r cutoffs)
r_cuts <- c(0, 0.1, 0.3, 0.5, 0.7)
p_cuts <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)

grid <- expand.grid(
  r_cut = r_cuts,
  p_cut = p_cuts,
  stringsAsFactors = FALSE
)

plot_list <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  r_cut <- grid$r_cut[i]
  p_cut <- grid$p_cut[i]

  pass_rows <- !is.na(fp_gh_map_corr$r_gene) &
    !is.na(fp_gh_map_corr$p_adj_gene) &
    abs(fp_gh_map_corr$r_gene) > r_cut &
    fp_gh_map_corr$p_adj_gene < p_cut

  genes_pass <- intersect(
    genes_univ,
    unique(fp_gh_map_corr$gene_key[pass_rows])
  )
  genes_fail <- setdiff(genes_univ, genes_pass)

  ko_pass <- ko_tbl[ko_tbl$symbol %in% genes_pass, , drop = FALSE]
  ko_fail <- ko_tbl[ko_tbl$symbol %in% genes_fail, , drop = FALSE]

  plot_list[[i]] <- dplyr::bind_rows(
    tibble::tibble(
      r_cut   = r_cut,
      p_cut   = p_cut,
      group   = "Predicted_regulated",
      gene    = ko_pass$symbol,
      log2FC  = ko_pass$log2FoldChange
    ),
    tibble::tibble(
      r_cut   = r_cut,
      p_cut   = p_cut,
      group   = "Predicted_nonregulated",
      gene    = ko_fail$symbol,
      log2FC  = ko_fail$log2FoldChange
    )
  )
}

plot_df <- dplyr::bind_rows(plot_list)

# Factor labels for facets and x-axis
plot_df <- plot_df %>%
  dplyr::mutate(
    r_lab = factor(
      r_cut,
      levels = r_cuts,
      labels = paste0("|r_gene| > ", r_cuts)
    ),
    p_lab = factor(
      p_cut,
      levels = sort(unique(p_cuts)),
      labels = label_scientific(digits = 1)(sort(unique(p_cuts)))
    ),
    group = factor(
      group,
      levels = c("Predicted_regulated", "Predicted_nonregulated")
    )
  )

# -------------------------------------------------------------------
# Compute Wilcoxon p-values per (r_cut, p_cut) – use pick() (no cur_data_all)
# -------------------------------------------------------------------
pval_df <- plot_df %>%
  dplyr::group_by(r_cut, p_cut) %>%
  dplyr::summarise(
    p_wilcox = {
      df_sub <- dplyr::filter(dplyr::pick(dplyr::everything()), !is.na(log2FC))
      gtab   <- table(df_sub$group)
      if (length(gtab) == 2L && all(gtab > 0L)) {
        stats::wilcox.test(log2FC ~ group, data = df_sub)$p.value
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  )

# y-position for labels: slightly above max value in each panel
y_pos_df <- plot_df %>%
  dplyr::group_by(r_cut, p_cut) %>%
  dplyr::summarise(
    y = max(log2FC, na.rm = TRUE),
    .groups = "drop"
  )

global_range <- range(plot_df$log2FC, na.rm = TRUE)
offset       <- 0.05 * diff(global_range)

annot_df <- dplyr::left_join(pval_df, y_pos_df, by = c("r_cut", "p_cut")) %>%
  dplyr::mutate(
    label = dplyr::case_when(
      is.na(p_wilcox)         ~ "N.S.",
      p_wilcox >= 0.05        ~ "N.S.",
      TRUE                    ~ paste0(
        "p = ",
        scales::label_scientific(digits = 2)(p_wilcox)
      )
    ),
    y     = y + offset,
    r_lab = factor(
      r_cut,
      levels = r_cuts,
      labels = paste0("|r_gene| > ", r_cuts)
    ),
    p_lab = factor(
      p_cut,
      levels = sort(unique(p_cuts)),
      labels = label_scientific(digits = 1)(sort(unique(p_cuts)))
    )
  )

# -------------------------------------------------------------------
# Shared fill colors (consistent across both rows)
# -------------------------------------------------------------------
fill_vals <- c(
  Predicted_regulated    = "#1b9e77",
  Predicted_nonregulated = "#d95f02"
)

# -------------------------------------------------------------------
# Row 1: Split-violin grid (no individual dots, just violins + box)
# -------------------------------------------------------------------
p_grid <- ggplot(
  plot_df,
  aes(x = p_lab, y = log2FC, fill = group)
) +
  geom_split_violin(alpha = 0.4, trim = FALSE, colour = "black") +
  geom_boxplot(
    width = 0.2,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  geom_text(
    data = annot_df,
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
  scale_x_discrete(name = "p_adj cutoff") +
  scale_y_continuous(name = "log2FC (HNF4A KO vs Ctrl)") +
  ggtitle(
    sprintf(
      "%s knockout: predicted regulated vs non-regulated genes\n(n = %d genes in universe)",
      tf, n_univ
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
# Row 2: Bar plot of N unique genes per group, aligned by (r_cut, p_cut)
# -------------------------------------------------------------------
counts_df <- plot_df %>%
  dplyr::group_by(r_cut, p_cut, group) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(gene),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    r_lab = factor(
      r_cut,
      levels = r_cuts,
      labels = paste0("|r_gene| > ", r_cuts)
    ),
    p_lab = factor(
      p_cut,
      levels = sort(unique(p_cuts)),
      labels = label_scientific(digits = 1)(sort(unique(p_cuts)))
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
  scale_x_discrete(name = "p_adj cutoff") +
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
combined <- p_grid / p_counts + patchwork::plot_layout(heights = c(3, 1))

ggsave(
  filename = file.path(
    ko_dir,
    sprintf("%s_KO_grid_reg_vs_nonreg_split_violin_by_r_p.pdf", tf)
  ),
  plot   = combined,
  width  = 18,
  height = 8,
  dpi    = 600
)
