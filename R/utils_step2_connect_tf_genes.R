# Build a reference gene table for hg38 (GRCh38) or mm10 (GRCm38)
# Returns tibble: ensembl_gene_id, HGNC, chrom, start, end, strand, tss
#' @export
episcope_build_gene_annot <- function(assembly = c("hg38","mm10"), ensdb = NULL) {
  assembly <- base::match.arg(assembly)

  if (!requireNamespace("ensembldb", quietly = TRUE))
    stop("Please install Bioconductor package 'ensembldb'.")

  # 1) Prefer installed EnsDb.*; else fall back to AnnotationHub
  if (is.null(ensdb)) {
    pkg_candidates <- if (assembly == "hg38") {
      c("EnsDb.Hsapiens.v105", "EnsDb.Hsapiens.v86")   # GRCh38
    } else {
      c("EnsDb.Mmusculus.v86", "EnsDb.Mmusculus.v79")  # GRCm38 (mm10)
    }
    for (p in pkg_candidates) {
      if (requireNamespace(p, quietly = TRUE)) {
        ensdb <- get(p, envir = asNamespace(p))
        break
      }
    }
    if (is.null(ensdb)) {
      if (!requireNamespace("AnnotationHub", quietly = TRUE))
        stop("Install 'AnnotationHub' or one of: ", paste(pkg_candidates, collapse = ", "))
      ah <- AnnotationHub::AnnotationHub()
      q  <- AnnotationHub::query(
        ah,
        c("EnsDb",
          if (assembly=="hg38") "Homo sapiens" else "Mus musculus",
          if (assembly=="hg38") "GRCh38"       else "GRCm38")
      )
      if (!length(q)) stop("No suitable EnsDb found in AnnotationHub for ", assembly)
      ensdb <- q[[length(q)]]  # take the latest hit
    }
  }

  # 2) Pull genes; take coords/strand from the GRanges itself (portable across builds)
  gr <- ensembldb::genes(ensdb, columns = c("gene_id","gene_name","seq_name"))

  # seqnames -> UCSC-style with 'chr' prefix if missing
  sn <- as.character(GenomicRanges::seqnames(gr))
  chrom <- ifelse(grepl("^chr", sn), sn, paste0("chr", sn))

  start_v <- as.integer(GenomicRanges::start(gr))
  end_v   <- as.integer(GenomicRanges::end(gr))
  strand_chr <- as.character(GenomicRanges::strand(gr))
  strand_chr[strand_chr == "*"] <- "."

  tibble::tibble(
    ensembl_gene_id = S4Vectors::mcols(gr)$gene_id,
    HGNC            = S4Vectors::mcols(gr)$gene_name,  # mouse: this is the MGI symbol
    chrom           = chrom,
    start           = start_v,
    end             = end_v,
    strand          = strand_chr,
    tss             = ifelse(strand_chr == "-", end_v, start_v)
  ) |>
    dplyr::filter(!is.na(chrom), !is.na(start), !is.na(end)) |>
    dplyr::distinct()
}

#' Normalize GH columns to a standard set.
.normalize_gh_cols <- function(df) {
  pick <- function(opts) {
    nm <- base::intersect(opts, names(df))
    if (length(nm)) nm[1] else NA_character_
  }
  got <- list(
    chrom          = pick(c("chrom","chr","Chrom","Chromosome","seqnames")),
    start          = pick(c("start","Start","chromStart","txStart","enhancer_start")),
    end            = pick(c("end","End","chromEnd","txEnd","enhancer_end")),
    connected_gene = pick(c("connected_gene","GeneSymbol","gene","symbol","target_gene","gene_name")),
    gh_id          = pick(c("genehancer_id","GHid","gh_id","enhancer_id","id")),
    confidence     = pick(c("score","Score","confidence","GH_Score","GeneHancer_score"))
  )
  req <- c("chrom","start","end","connected_gene")
  if (any(vapply(got[req], is.na, logical(1)))) {
    stop("GeneHancer file is missing required columns (chrom/start/end/connected_gene).")
  }

  rename_if_present <- function(d, old, new) {
    if (!is.na(old) && old %in% names(d) && old != new) {
      dplyr::rename(d, !!new := dplyr::all_of(old))
    } else d
  }

  df2 <- df
  df2 <- rename_if_present(df2, got$chrom,          "chrom")
  df2 <- rename_if_present(df2, got$start,          "start")
  df2 <- rename_if_present(df2, got$end,            "end")
  df2 <- rename_if_present(df2, got$connected_gene, "connected_gene")
  df2 <- rename_if_present(df2, got$gh_id,          "gh_id")
  df2 <- rename_if_present(df2, got$confidence,     "confidence")

  if (!"gh_id" %in% names(df2))       df2$gh_id       <- NA_character_
  if (!"confidence" %in% names(df2))  df2$confidence  <- NA_real_

  df2 |>
    dplyr::mutate(
      chrom      = as.character(.data$chrom),
      start      = as.integer(.data$start),
      end        = as.integer(.data$end),
      gh_id      = as.character(.data$gh_id),
      confidence = suppressWarnings(as.numeric(.data$confidence))
    ) |>
    dplyr::select(.data$chrom, .data$start, .data$end, .data$connected_gene,
                  .data$gh_id, .data$confidence, dplyr::everything())
}
#' Load the pancreas GeneHancer ELITE CSV from inst/extdata and normalize.
#' @export
load_genehancer_panc <- function(path = file.path("inst","extdata","GeneHancer_v5.24_elite_panc.csv")) {
  if (!file.exists(path)) stop("GeneHancer file not found at: ", path)
  gh_raw <- readr::read_csv(path, show_col_types = FALSE)
  .normalize_gh_cols(gh_raw)
}

.parse_peak_id <- function(x) {
  x2 <- base::gsub(":", "-", x, fixed = TRUE)
  sp <- base::strsplit(x2, "-", fixed = TRUE)
  chr   <- vapply(sp, `[[`, character(1), 1)
  start <- suppressWarnings(as.integer(vapply(sp, `[[`, character(1), 2)))
  end   <- suppressWarnings(as.integer(vapply(sp, `[[`, character(1), 3)))
  tibble::tibble(peak_ID = x, chrom = chr, start = start, end = end)
}

#' Build GH-like peak->gene links by distance window (no GeneHancer needed)
#'
#' Create GeneHancer-style rows by linking peaks to nearby genes within a
#' +/-\code{flank_bp} window. Works with either a ready-made gene table or a
#' character vector of gene IDs (HGNC or Ensembl) that are looked up in a
#' provided annotation table.
#'
#' @param peaks Tibble with either \code{peak_ID} or \code{atac_peak}
#'   (both formatted \code{chr:start-end}) or explicit genomic columns
#'   \code{chrom,start,end}.
#' @param genes EITHER:
#'   \itemize{
#'     \item a tibble with \code{chrom}, \code{gene_key}, and either \code{tss}
#'           or \code{start,end,strand}; or
#'     \item a character vector of gene IDs (HGNC symbols or Ensembl IDs) to be
#'           looked up in \code{gene_annot}.
#'   }
#' @param flank_bp Integer. Window size in bp to extend each peak on both sides
#'   (default \code{20000}).
#' @param mode Character. \code{"TSS"} (default) uses transcription start site;
#'   \code{"full"} uses full gene span.
#' @param gene_annot Optional tibble used only when \code{genes} is a character
#'   vector. Must contain \code{chrom,start,end,strand} and an ID column
#'   (e.g., \code{HGNC} or \code{ensembl_gene_id}) matching \code{genes}.
#' @param id_col Optional string naming the ID column in \code{gene_annot}.
#'   If \code{NULL}, tries \code{"HGNC"}, \code{"gene_name"},
#'   \code{"ensembl_gene_id"}, then \code{"gene_id"}.
#'
#' @return Tibble with GH-like columns ready for \code{build_gh_pairs_from_peaks()}:
#'   \code{chrom,start,end,connected_gene,gh_id,confidence,source,feature_name,
#'   strand,frame,is_elite_elem}.
#'
#' @examples
#' # Using a genes tibble that already has TSS:
#' # genes_tbl <- tibble::tibble(chrom="chr1", tss=1e6, gene_key="TP53")
#' # gh_std <- episcope_make_windowed_gh(strict_set$atac_mt, genes_tbl, flank_bp=2e5, mode="TSS")
#'
#' # Using a vector of IDs + annotation look-up:
#' # gene_annot_ref must have: chrom,start,end,strand, and an ID col (e.g. HGNC)
#' # gh_std <- episcope_make_windowed_gh(
#' #   peaks = strict_set$atac_score,   # now supports `atac_peak`
#' #   genes = c("SOX9","HNF1B","ENSG00000141510"),
#' #   flank_bp = 2e5, mode = "TSS",
#' #   gene_annot = gene_annot_ref, id_col = "HGNC"
#' # )
#'
#' @export
episcope_make_windowed_gh <- function(peaks, genes,
                                      flank_bp = 30000L,
                                      mode = c("TSS","full"),
                                      gene_annot = NULL,
                                      id_col = NULL) {
  mode <- base::match.arg(mode)

  # ---- Peaks -> (chrom,start,end)
  # MINIMAL CHANGE: also accept `atac_peak` formatted "chr:start-end"
  if ("peak_ID" %in% names(peaks)) {
    pk <- .parse_peak_id(peaks$peak_ID)
  } else if ("atac_peak" %in% names(peaks)) {
    pk <- .parse_peak_id(peaks$atac_peak)
  } else {
    pk <- tibble::as_tibble(peaks)
  }
  if (!all(c("chrom","start","end") %in% names(pk))) {
    cli::cli_abort("`peaks` must contain either `peak_ID`, `atac_peak`, or columns `chrom,start,end`.")
  }
  pk <- dplyr::transmute(pk,
                         chrom = as.character(.data$chrom),
                         start = as.integer(.data$start),
                         end   = as.integer(.data$end))

  # ---- Genes input normalization -> gn0: (chrom,start,end,gene_key) for chosen mode
  if (is.character(genes)) {
    if (is.null(gene_annot)) {
      cli::cli_abort("When `genes` is a character vector, supply `gene_annot` with coordinates.")
    }
    ga <- tibble::as_tibble(gene_annot)
    if (is.null(id_col)) {
      id_col <- if ("HGNC" %in% names(ga)) "HGNC" else
        if ("gene_name" %in% names(ga)) "gene_name" else
          if ("ensembl_gene_id" %in% names(ga)) "ensembl_gene_id" else
            if ("gene_id" %in% names(ga)) "gene_id" else
              cli::cli_abort("Could not auto-detect an ID column in `gene_annot`. Provide `id_col`.")
    }
    need_cols <- c("chrom","start","end","strand", id_col)
    if (!all(need_cols %in% names(ga))) {
      cli::cli_abort("`gene_annot` must contain columns: {.val {need_cols}}.")
    }
    ga <- dplyr::filter(ga, .data[[id_col]] %in% genes)

    if (mode == "TSS") {
      if ("tss" %in% names(ga)) {
        gn0 <- dplyr::transmute(ga,
                                chrom = .data$chrom,
                                start = as.integer(.data$tss),
                                end   = as.integer(.data$tss) + 1L,
                                gene_key = .data[[id_col]])
      } else {
        gn0 <- dplyr::transmute(ga,
                                chrom = .data$chrom,
                                start = as.integer(ifelse(.data$strand == "-", .data$end, .data$start)),
                                end   = start + 1L,
                                gene_key = .data[[id_col]])
      }
    } else {
      gn0 <- dplyr::transmute(ga,
                              chrom = .data$chrom,
                              start = as.integer(.data$start),
                              end   = as.integer(.data$end),
                              gene_key = .data[[id_col]])
    }

  } else {
    if (!all(c("chrom","gene_key") %in% names(genes))) {
      cli::cli_abort("`genes` tibble must contain `chrom` and `gene_key` plus `tss` or `start,end,strand`.")
    }
    if (mode == "TSS") {
      if ("tss" %in% names(genes)) {
        gn0 <- dplyr::transmute(genes,
                                chrom = .data$chrom,
                                start = as.integer(.data$tss),
                                end   = as.integer(.data$tss) + 1L,
                                gene_key = .data$gene_key)
      } else {
        if (!all(c("start","end","strand") %in% names(genes))) {
          cli::cli_abort("`genes` needs `tss` or the trio `start,end,strand` for mode = 'TSS'.")
        }
        gn0 <- dplyr::transmute(genes,
                                chrom = .data$chrom,
                                start = as.integer(ifelse(.data$strand == "-", .data$end, .data$start)),
                                end   = start + 1L,
                                gene_key = .data$gene_key)
      }
    } else {
      if (!all(c("start","end") %in% names(genes))) {
        cli::cli_abort("`genes` needs `start,end` for mode = 'full'.")
      }
      gn0 <- dplyr::transmute(genes,
                              chrom = .data$chrom,
                              start = as.integer(.data$start),
                              end   = as.integer(.data$end),
                              gene_key = .data$gene_key)
    }
  }

  # ---- Extend peaks by +/-flank and intersect with genes (report original peak span)
  pk_win <- dplyr::transmute(pk,
                             chrom,
                             start = pmax(start - flank_bp, 0L),
                             end   = end + flank_bp,
                             orig_start = start,
                             orig_end   = end)

  ov <- valr::bed_intersect(pk_win, gn0)
  nm <- names(ov)

  chrom_col <- if ("chrom.x" %in% nm) "chrom.x" else if ("chrom" %in% nm) "chrom" else "chrom.y"
  os_col    <- if ("orig_start.x" %in% nm) "orig_start.x" else "orig_start"
  oe_col    <- if ("orig_end.x"   %in% nm) "orig_end.x"   else "orig_end"
  gk_col    <- if ("gene_key.y"   %in% nm) "gene_key.y"   else "gene_key"

  out <- ov |>
    dplyr::transmute(
      chrom          = !!rlang::sym(chrom_col),
      start          = !!rlang::sym(os_col),
      end            = !!rlang::sym(oe_col),
      connected_gene = !!rlang::sym(gk_col)
    ) |>
    dplyr::mutate(
      gh_id        = sprintf("WIN_%s:%d-%d", .data$chrom, .data$start, .data$end),
      confidence   = 1,
      source       = sprintf("%dkb", as.integer(flank_bp/1000)),
      feature_name = "Promoter/Enhancer",
      strand       = ".",
      frame        = ".",
      is_elite_elem = 0L
    ) |>
    dplyr::distinct()

  out
}

#' Build GH overlap pairs from peaks (GenomicRanges).
build_gh_pairs_from_peaks <- function(atac_tbl, gh_tbl_std) {
  if (!"peak_ID" %in% names(atac_tbl)) stop("atac_tbl must contain 'peak_ID'.")
  peaks <- .parse_peak_id(atac_tbl$peak_ID) |>
    dplyr::filter(!is.na(.data$start), !is.na(.data$end))

  peaks_gr <- GenomicRanges::GRanges(seqnames = peaks$chrom,
                                     ranges   = IRanges::IRanges(peaks$start, peaks$end))
  gh_gr    <- GenomicRanges::GRanges(seqnames = gh_tbl_std$chrom,
                                     ranges   = IRanges::IRanges(gh_tbl_std$start, gh_tbl_std$end))

  hits <- GenomicRanges::findOverlaps(peaks_gr, gh_gr, ignore.strand = TRUE)

  tibble::tibble(
    peak_ID        = peaks$peak_ID[as.integer(S4Vectors::queryHits(hits))],
    gene_key       = gh_tbl_std$connected_gene[as.integer(S4Vectors::subjectHits(hits))],
    gh_id          = gh_tbl_std$gh_id[as.integer(S4Vectors::subjectHits(hits))],
    gh_confidence  = gh_tbl_std$confidence[as.integer(S4Vectors::subjectHits(hits))],
    enh_chr        = gh_tbl_std$chrom[as.integer(S4Vectors::subjectHits(hits))],
    enh_start      = gh_tbl_std$start[as.integer(S4Vectors::subjectHits(hits))],
    enh_end        = gh_tbl_std$end[as.integer(S4Vectors::subjectHits(hits))]
  ) |>
    dplyr::distinct()
}

# helpers
.parse_peak_ids <- function(ids) {
  rng <- sub(".*:", "", ids)
  tibble::tibble(
    atac_peak = ids,
    atac_chr  = sub(":.*", "", ids),
    atac_start = suppressWarnings(as.integer(sub("-.*", "", rng))),
    atac_end   = suppressWarnings(as.integer(sub(".*-", "", rng)))
  )
}
.pick_rna_gene_keys <- function(rna_tbl, mode = c("both","ensembl","hgnc")) {
  mode <- match.arg(mode)
  have_hgnc <- "HGNC" %in% names(rna_tbl)
  have_ens  <- "ensembl_gene_id" %in% names(rna_tbl)

  if (mode == "hgnc") {
    if (!have_hgnc) cli::cli_abort("RNA table has no HGNC column for mode='hgnc'.")
    out <- unique(stats::na.omit(rna_tbl$HGNC))
  } else if (mode == "ensembl") {
    if (!have_ens) cli::cli_abort("RNA table has no ensembl_gene_id column for mode='ensembl'.")
    out <- unique(stats::na.omit(rna_tbl$ensembl_gene_id))
  } else {
    keys <- character()
    if (have_hgnc) keys <- c(keys, rna_tbl$HGNC)
    if (have_ens)  keys <- c(keys, rna_tbl$ensembl_gene_id)
    out <- unique(stats::na.omit(keys))
  }
  out[nzchar(out)]
}

# --- main

#' Correlate ATAC peaks to connected genes (generic enhancer/promoter sources) and filter by FDR.
#'
#' @param grn_set   GRN set (from build_grn_set).
#' @param gh_tbl    Tibble with *at least* columns: chrom, start, end, connected_gene
#'                  (names configurable via gh_cols).
#' @param gh_cols   Named list mapping column names:
#'                  list(chrom="chrom", start="start", end="end", gene="connected_gene").
#'                  Optional extras (if present) are kept during overlap: id="gh_id", conf="confidence".
#' @param gene_mode "both", "ensembl", or "hgnc" - which gene identifiers to match to RNA.
#' @param fdr       BH FDR cutoff.
#' @param keep_pos  If TRUE, retain only positive correlations.
#' @param r_abs_min Absolute correlation cutoff (applied to |r_atac|).
#' @param cor_method Correlation method passed to stats::cor: "pearson" or "spearman".
#' @param workers   Parallel workers. Use 1 for serial.
#'
#' Caching options (optional; default = no caching, original behaviour):
#'
#' @param cache_dir Optional directory for on-disk caching. If NULL (default), no
#'                  caching is used and the function behaves as before.
#' @param cache_tag Character label to namespace caches within `cache_dir`. If NULL,
#'                  a tag is derived from `gene_mode` and `cor_method`.
#' @param cache_chunk_size Approximate number of (atac_peak, gene_key) pairs per
#'                  CSV chunk file when caching. More chunks -> more files, finer
#'                  granularity for resume; fewer chunks -> larger files.
#' @param cache_resume Logical; if TRUE (default), re-use any existing cached
#'                  correlations matching the requested (atac_peak, gene_key) pairs.
#'                  If FALSE, ignore old cache contents and overwrite with new results.
#' @param cache_verbose Logical; if TRUE, emit detailed messages about cache/index/chunks.
#'
#' @export
#'
#' @return list(atac_gene_corr_kept, atac_gene_corr_full)
correlate_atac_to_genes <- function(
    grn_set,
    gh_tbl,
    gh_cols   = list(chrom="chrom", start="start", end="end", gene="connected_gene",
                     id=NULL, conf=NULL),
    gene_mode = c("both","ensembl","hgnc"),
    fdr       = 0.05,
    keep_pos  = FALSE,
    r_abs_min = 0.3,   # absolute correlation cutoff
    cor_method = c("pearson", "spearman"),
    workers   = max(1L, parallel::detectCores() - 1L),
    cache_dir        = NULL,
    cache_tag        = NULL,
    cache_chunk_size = 5000L,
    cache_resume     = TRUE,
    cache_verbose    = FALSE
) {
  gene_mode  <- match.arg(gene_mode)
  cor_method <- match.arg(cor_method)

  stopifnot(
    is.data.frame(grn_set$fp_annotation),
    is.data.frame(grn_set$atac_score),
    is.data.frame(grn_set$rna),
    is.data.frame(grn_set$sample_metadata_used)
  )

  # ---------------------------------------------------------------------------
  # samples
  # ---------------------------------------------------------------------------
  ids <- grn_set$sample_metadata_used$id
  ids <- ids[ids %in% intersect(names(grn_set$atac_score), names(grn_set$rna))]
  if (!length(ids)) cli::cli_abort("No overlapping sample IDs between ATAC and RNA.")

  # ATAC peaks from fp_annotation
  ap <- unique(grn_set$fp_annotation$atac_peak)
  ap <- ap[!is.na(ap) & nzchar(ap)]
  if (!length(ap)) cli::cli_abort("No usable atac_peak values in fp_annotation.")

  atac_bed <- .parse_peak_ids(ap) |>
    dplyr::select(atac_chr, atac_start, atac_end, atac_peak)

  # ---------------------------------------------------------------------------
  # GH normalization (generic)
  # ---------------------------------------------------------------------------
  req <- c(gh_cols$chrom, gh_cols$start, gh_cols$end, gh_cols$gene)
  if (!all(req %in% names(gh_tbl))) {
    cli::cli_abort("Enhancer table missing required columns: {.val {req}}.")
  }

  have_id   <- !is.null(gh_cols$id)   && gh_cols$id   %in% names(gh_tbl)
  have_conf <- !is.null(gh_cols$conf) && gh_cols$conf %in% names(gh_tbl)

  gh_bed <- gh_tbl |>
    dplyr::transmute(
      gh_chr   = .data[[gh_cols$chrom]],
      gh_start = as.integer(.data[[gh_cols$start]]),
      gh_end   = as.integer(.data[[gh_cols$end]]),
      gene_key = as.character(.data[[gh_cols$gene]]),
      gh_id    = if (have_id)  as.character(.data[[gh_cols$id]])      else NA_character_,
      gh_conf  = if (have_conf) suppressWarnings(as.numeric(.data[[gh_cols$conf]])) else NA_real_
    ) |>
    dplyr::filter(!is.na(gh_start), !is.na(gh_end), !is.na(gene_key), nzchar(gene_key))

  # ---------------------------------------------------------------------------
  # overlap (any)
  # ---------------------------------------------------------------------------
  ov_raw <- valr::bed_intersect(
    atac_bed |> dplyr::transmute(chrom = atac_chr, start = atac_start, end = atac_end, atac_peak),
    gh_bed   |> dplyr::transmute(chrom = gh_chr,  start = gh_start,  end = gh_end, gh_id, gh_conf, gene_key)
  )

  nm <- names(ov_raw)
  atk_col  <- if ("atac_peak.x" %in% nm) "atac_peak.x" else "atac_peak"
  ghid_col <- if ("gh_id.y"     %in% nm) "gh_id.y"     else "gh_id"
  ghy_col  <- if ("gh_conf.y"   %in% nm) "gh_conf.y"   else "gh_conf"
  ggn_col  <- if ("gene_key.y"  %in% nm) "gene_key.y"  else "gene_key"

  pairs0 <- ov_raw |>
    dplyr::transmute(
      atac_peak = .data[[atk_col]],
      gene_key  = .data[[ggn_col]],
      gh_id     = .data[[ghid_col]],
      gh_conf   = .data[[ghy_col]]
    ) |>
    dplyr::distinct()

  if (!nrow(pairs0)) cli::cli_abort("No ATAC<->enhancer overlaps.")

  # ---------------------------------------------------------------------------
  # keep genes present in RNA under the requested mode
  # ---------------------------------------------------------------------------
  rna_keys <- .pick_rna_gene_keys(grn_set$rna, gene_mode)
  pairs1 <- dplyr::semi_join(pairs0, tibble::tibble(gene_key = rna_keys), by = "gene_key")
  if (!nrow(pairs1)) cli::cli_abort("After gene matching to RNA, no pairs remain.")

  # only ATAC peaks present in atac_score
  pairs1 <- dplyr::semi_join(pairs1, dplyr::distinct(grn_set$atac_score, atac_peak), by = "atac_peak")

  # ---------------------------------------------------------------------------
  # long views
  # ---------------------------------------------------------------------------
  atac_long <- grn_set$atac_score |>
    dplyr::semi_join(dplyr::distinct(pairs1, atac_peak), by = "atac_peak") |>
    dplyr::select(atac_peak, dplyr::all_of(ids)) |>
    tidyr::pivot_longer(dplyr::all_of(ids), names_to = "sample", values_to = "acc")

  rna_long <- dplyr::bind_rows(
    if ("HGNC" %in% names(grn_set$rna)) {
      grn_set$rna |>
        dplyr::filter(!is.na(HGNC), HGNC != "") |>
        dplyr::select(gene_key = HGNC, dplyr::all_of(ids))
    } else NULL,
    if ("ensembl_gene_id" %in% names(grn_set$rna)) {
      grn_set$rna |>
        dplyr::select(gene_key = ensembl_gene_id, dplyr::all_of(ids))
    } else NULL
  ) |>
    dplyr::filter(!is.na(gene_key), nzchar(gene_key)) |>
    tidyr::pivot_longer(dplyr::all_of(ids), names_to = "sample", values_to = "expr") |>
    dplyr::distinct(gene_key, sample, .keep_all = TRUE)

  pairs2 <- dplyr::semi_join(pairs1, dplyr::distinct(rna_long, gene_key), by = "gene_key")
  if (!nrow(pairs2)) cli::cli_abort("After RNA availability filter, no pairs remain.")

  # ---------------------------------------------------------------------------
  # helpers
  # ---------------------------------------------------------------------------
  .empty_corr_tbl <- function() {
    tibble::tibble(
      atac_peak = character(),
      gene_key  = character(),
      n         = integer(),
      r         = double(),
      p         = double()
    )
  }

  # inner correlation engine; can be serial or parallel over index chunks
  run_corr_for_pairs <- function(pairs_df, inner_workers) {
    if (!nrow(pairs_df)) {
      return(.empty_corr_tbl())
    }

    n_pairs <- nrow(pairs_df)
    nchunk  <- max(1L, inner_workers)

    idx_split <- split(
      seq_len(n_pairs),
      ceiling(seq_along(seq_len(n_pairs)) /
                max(1L, ceiling(n_pairs / max(1L, nchunk))))
    )

    run_chunk <- function(ii) {
      pp <- pairs_df[ii, , drop = FALSE]
      ax <- dplyr::inner_join(pp, atac_long, by = "atac_peak", relationship = "many-to-many")
      if (!nrow(ax)) {
        return(.empty_corr_tbl())
      }
      axx <- dplyr::inner_join(
        ax,
        rna_long,
        by = c("gene_key", "sample"),
        relationship = "many-to-many"
      )
      if (!nrow(axx)) {
        return(.empty_corr_tbl())
      }

      axx |>
        dplyr::group_by(atac_peak, gene_key) |>
        dplyr::summarise(
          n = sum(is.finite(acc) & is.finite(expr)),
          r = suppressWarnings(stats::cor(
            acc,
            expr,
            use    = "complete.obs",
            method = cor_method
          )),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          p = dplyr::case_when(
            is.na(r) | n < 3 ~ NA_real_,
            abs(r) >= 1      ~ 0,
            TRUE ~ {
              tt <- r * sqrt((n - 2) / (1 - r^2))
              2 * stats::pt(-abs(tt), df = n - 2)
            }
          )
        )
    }

    if (inner_workers > 1L) {
      out_list <- future.apply::future_lapply(idx_split, run_chunk, future.seed = TRUE)
    } else {
      out_list <- lapply(idx_split, run_chunk)
    }

    dplyr::bind_rows(out_list)
  }

  use_cache <- !is.null(cache_dir)
  cache_root <- NULL

  # ---------------------------------------------------------------------------
  # future plan setup (outer or inner parallel depending on use_cache)
  # ---------------------------------------------------------------------------
  have_parallel <- (workers > 1L)
  if (have_parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      cli::cli_abort("Install future.apply or set workers = 1.")
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
  }

  # ---------------------------------------------------------------------------
  # no-caching path: behaves like original, just with inner parallel
  # ---------------------------------------------------------------------------
  if (!use_cache) {
    inner_workers <- if (have_parallel) workers else 1L
    if (cache_verbose) {
      cli::cli_inform("Running without cache: computing {nrow(pairs2)} ATAC-gene pairs with {inner_workers} workers.")
    }
    corr_all <- run_corr_for_pairs(pairs2, inner_workers)
  } else {
    # -------------------------------------------------------------------------
    # cached path: chunked CSV files + central index
    # -------------------------------------------------------------------------
    cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

    if (is.null(cache_tag) || !nzchar(cache_tag)) {
      cache_tag <- paste0("genes-", gene_mode, "_cor-", cor_method)
    }

    cache_root <- file.path(cache_dir, cache_tag)
    if (!dir.exists(cache_root)) dir.create(cache_root, recursive = TRUE)

    meta_file  <- file.path(cache_root, "meta.rds")
    index_file <- file.path(cache_root, "index.rds")

    meta_new <- list(
      gene_mode   = gene_mode,
      cor_method  = cor_method
    )

    if (file.exists(meta_file)) {
      meta_old <- tryCatch(readRDS(meta_file), error = function(e) NULL)
      if (!is.null(meta_old)) {
        if (!identical(meta_old$gene_mode, gene_mode) ||
            !identical(meta_old$cor_method, cor_method)) {
          cli::cli_abort(
            c(
              "Existing cache at {.path {cache_root}} is not compatible.",
              "i" = "It was created with different gene_mode or cor_method.",
              "i" = "Use a different cache_tag or clean the cache directory."
            )
          )
        }
      }
    } else {
      saveRDS(meta_new, meta_file)
    }

    if (cache_verbose) {
      cli::cli_inform("Using ATAC-gene correlation cache at {.path {cache_root}}.")
    }

    # central index
    if (!cache_resume || !file.exists(index_file)) {
      cache_index <- tibble::tibble(
        atac_peak = character(),
        gene_key  = character(),
        file      = character(),
        chunk_id  = integer()
      )
    } else {
      cache_index <- tryCatch(readRDS(index_file), error = function(e) NULL)
      if (!is.data.frame(cache_index) ||
          !all(c("atac_peak", "gene_key", "file", "chunk_id") %in% names(cache_index))) {
        cli::cli_warn("Cache index at {.path {index_file}} is malformed; restarting index.")
        cache_index <- tibble::tibble(
          atac_peak = character(),
          gene_key  = character(),
          file      = character(),
          chunk_id  = integer()
        )
      }
    }

    run_keys <- pairs2 |>
      dplyr::select(atac_peak, gene_key) |>
      dplyr::distinct()

    total_pairs  <- nrow(run_keys)
    existing_n   <- 0L
    cached_pairs <- 0L
    existing_chunks <- 0L

    if (!nrow(cache_index)) {
      have_idx     <- tibble::tibble(atac_peak = character(), gene_key = character(), file = character(), chunk_id = integer())
      missing_keys <- run_keys
    } else {
      have_idx <- dplyr::inner_join(run_keys, cache_index, by = c("atac_peak", "gene_key"))
      missing_keys <- dplyr::anti_join(run_keys, cache_index, by = c("atac_peak", "gene_key"))
      existing_n <- nrow(have_idx)
      existing_chunks <- length(unique(cache_index$chunk_id))
    }

    missing_n <- nrow(missing_keys)

    # Basic logging on requested size / index contents
    cli::cli_inform("ATAC-gene correlation run: {total_pairs} unique pairs requested.")
    cli::cli_inform("Cache index currently tracks {nrow(cache_index)} pairs across {existing_chunks} chunks.")
    cli::cli_inform("{existing_n} pairs already satisfy this run; {missing_n} pairs remain to compute.")

    corr_cached <- .empty_corr_tbl()

    # ---------------- existing cached correlations ---------------------------
    if (existing_n > 0L) {
      files_used <- unique(have_idx$file)
      if (cache_verbose) {
        cli::cli_inform("Loading cached correlations from {length(files_used)} chunk files.")
      }

      for (ff in files_used) {
        idx_ff <- have_idx[have_idx$file == ff, c("atac_peak", "gene_key"), drop = FALSE]
        if (!file.exists(ff)) {
          cli::cli_warn("Cache chunk file missing: {.path {ff}}; those keys will be recomputed.")
          missing_keys <- dplyr::bind_rows(missing_keys, idx_ff) |>
            dplyr::distinct()
          next
        }
        df_ff <- readr::read_csv(ff, show_col_types = FALSE, progress = FALSE)
        df_sub <- dplyr::inner_join(df_ff, idx_ff, by = c("atac_peak", "gene_key"))
        if (nrow(df_sub)) {
          corr_cached <- dplyr::bind_rows(corr_cached, df_sub)
        }
      }
    }

    # ---------------- new correlations for keys not in cache -----------------
    corr_new <- .empty_corr_tbl()

    if (missing_n > 0L) {
      cache_chunk_size <- as.integer(cache_chunk_size)
      if (is.na(cache_chunk_size) || cache_chunk_size < 1L) cache_chunk_size <- 50000L

      nchunk_new  <- ceiling(missing_n / cache_chunk_size)

      cli::cli_inform(
        "Computing {missing_n} new ATAC-gene pairs in {nchunk_new} chunks (target chunk size ~{cache_chunk_size})."
      )
      if (cache_verbose) {
        cli::cli_inform("Requested workers: {workers}. With {nchunk_new} chunks, max parallelism is min(workers, nchunk_new).")
      }

      # assign chunk IDs
      old_chunks  <- if (nrow(cache_index)) cache_index$chunk_id else integer(0)
      start_id    <- if (length(old_chunks)) max(old_chunks, na.rm = TRUE) + 1L else 1L
      chunk_vec   <- rep(
        seq_len(nchunk_new),
        each = cache_chunk_size,
        length.out = missing_n
      ) + (start_id - 1L)

      missing_keys$chunk_id <- as.integer(chunk_vec)

      # expand back to full pairs2 rows and split by chunk
      missing_full <- dplyr::inner_join(
        pairs2,
        missing_keys,
        by = c("atac_peak", "gene_key")
      )

      split_chunks <- split(missing_full, missing_full$chunk_id)

      # small summary of actual chunk sizes
      chunk_sizes <- vapply(split_chunks, nrow, integer(1L))
      if (cache_verbose) {
        size_summary <- summary(chunk_sizes)
        cli::cli_inform(
          "Chunk sizes (n pairs) summary: Min={min(chunk_sizes)}, Median={median(chunk_sizes)}, Max={max(chunk_sizes)}, Mean={round(mean(chunk_sizes), 1)}."
        )
      }

      # chunk worker: no inner parallel; each chunk is serial but
      # many chunks run in parallel across "workers"
      chunk_fun <- function(df_chunk) {
        chunk_id <- unique(df_chunk$chunk_id)
        if (length(chunk_id) != 1L) {
          cli::cli_warn("Chunk has multiple chunk_id values; using the first.")
          chunk_id <- chunk_id[1L]
        }
        chunk_id <- as.integer(chunk_id[1L])

        chunk_pairs <- df_chunk[, c("atac_peak", "gene_key", "gh_id", "gh_conf"), drop = FALSE]
        chunk_pairs <- dplyr::distinct(chunk_pairs)

        if (cache_verbose) {
          cli::cli_inform("Worker computing chunk {chunk_id} with {nrow(chunk_pairs)} unique pairs.")
        }

        # compute correlations serially inside this worker
        corr_chunk <- run_corr_for_pairs(
          pairs_df      = chunk_pairs,
          inner_workers = 1L
        )

        if (!nrow(corr_chunk)) {
          return(list(
            index = tibble::tibble(
              atac_peak = character(),
              gene_key  = character(),
              file      = character(),
              chunk_id  = integer()
            ),
            corr  = .empty_corr_tbl()
          ))
        }

        chunk_file <- file.path(cache_root, sprintf("chunk_%05d.csv.gz", chunk_id))
        readr::write_csv(corr_chunk, chunk_file)

        index_chunk <- tibble::tibble(
          atac_peak = corr_chunk$atac_peak,
          gene_key  = corr_chunk$gene_key,
          file      = chunk_file,
          chunk_id  = chunk_id
        )

        list(index = index_chunk, corr = corr_chunk)
      }

      if (have_parallel) {
        res_list <- future.apply::future_lapply(split_chunks, chunk_fun, future.seed = TRUE)
      } else {
        res_list <- lapply(split_chunks, chunk_fun)
      }

      index_new <- dplyr::bind_rows(lapply(res_list, function(x) x$index))
      corr_new  <- dplyr::bind_rows(lapply(res_list, function(x) x$corr))

      cache_index <- dplyr::bind_rows(cache_index, index_new) |>
        dplyr::arrange(atac_peak, gene_key) |>
        dplyr::distinct(atac_peak, gene_key, .keep_all = TRUE)

      saveRDS(cache_index, index_file)
      if (cache_verbose) {
        cli::cli_inform("Updated cache index now tracks {nrow(cache_index)} pairs across {length(unique(cache_index$chunk_id))} chunks.")
      }
    } else if (cache_verbose) {
      cli::cli_inform("No new pairs to compute; all requested correlations already present in cache.")
    }

    # union cached + new for this run
    corr_all <- dplyr::bind_rows(corr_cached, corr_new) |>
      dplyr::arrange(atac_peak, gene_key) |>
      dplyr::distinct(atac_peak, gene_key, .keep_all = TRUE)
  }

  # ---------------------------------------------------------------------------
  # finalize results (same semantics as original)
  # ---------------------------------------------------------------------------
  atac_gene_corr_full <- corr_all |>
    dplyr::filter(!is.na(p)) |>
    dplyr::mutate(p_adj = stats::p.adjust(p, method = "BH")) |>
    dplyr::transmute(
      atac_peak,
      gene_key,
      n_atac     = as.integer(n),
      r_atac     = as.numeric(r),
      p_atac     = as.numeric(p),
      p_adj_atac = as.numeric(p_adj)
    )

  atac_gene_corr_kept <- atac_gene_corr_full |>
    dplyr::filter(
      p_adj_atac < fdr,
      abs(r_atac) >= r_abs_min,
      if (keep_pos) r_atac > 0 else TRUE
    ) |>
    dplyr::distinct(atac_peak, gene_key, .keep_all = TRUE)

  list(
    atac_gene_corr_kept = atac_gene_corr_kept,
    atac_gene_corr_full = atac_gene_corr_full
  )
}



#' @export
filter_grn_by_corr <- function(grn_set, corr_kept) {
  stopifnot(is.list(grn_set), is.data.frame(corr_kept))
  if (!"atac_peak" %in% names(corr_kept))
    stop("corr_kept must have an 'atac_peak' column.")

  keep_atac <- unique(stats::na.omit(corr_kept$atac_peak))

  # 1) Filter ATAC matrices by correlated peaks
  grn_set$atac_score <- grn_set$atac_score |>
    dplyr::semi_join(tibble::tibble(atac_peak = keep_atac), by = "atac_peak")

  grn_set$atac_overlap <- grn_set$atac_overlap |>
    dplyr::semi_join(tibble::tibble(atac_peak = keep_atac), by = "atac_peak")

  # 2) Filter fp_annotation by those atac peaks
  grn_set$fp_annotation <- grn_set$fp_annotation |>
    dplyr::semi_join(tibble::tibble(atac_peak = keep_atac), by = "atac_peak")

  # 3) Ensure fp_bound / fp_score match the remaining fp_annotation
  keep_fp <- unique(grn_set$fp_annotation$fp_peak)
  if (!length(keep_fp)) {
    # nothing left; return empty but consistent tables
    grn_set$fp_score <- grn_set$fp_score[0, , drop = FALSE]
    grn_set$fp_bound <- grn_set$fp_bound[0, , drop = FALSE]
    return(grn_set)
  }

  grn_set$fp_score <- grn_set$fp_score |>
    dplyr::semi_join(tibble::tibble(peak_ID = keep_fp), by = "peak_ID")

  grn_set$fp_bound <- grn_set$fp_bound |>
    dplyr::semi_join(tibble::tibble(peak_ID = keep_fp), by = "peak_ID")

  grn_set
}

#' Correlate fp_peak scores to target-gene RNA for ATAC->gene pairs
#'
#' This function correlates footprint peak (fp_peak) scores to target-gene RNA
#' for ATAC->gene pairs, with optional disk caching for the FP->gene correlations.
#'
#' Update (2025-12-16):
#' - Adds optional TF RNA expression -> target-gene RNA expression correlations
#'   (per unique tfs x gene_key) and joins results back onto both outputs as:
#'   r_rna, p_rna, p_adj_rna.
#' - Existing FP->gene correlation logic, filtering, and caching behavior are unchanged.
#'
#' @param grn_set List containing at least `fp_annotation` and (by default)
#'   `fp_score` and `rna` unless overridden by `fp_score_tbl` / `rna_tbl`.
#' @param atac_gene_corr_kept Optional data frame with at least `atac_peak` and
#'   `gene_key`. If NULL, `gh_tbl` is used to derive peak->gene pairs.
#' @param gh_tbl Optional enhancer->gene table (e.g., GeneHancer). If
#'   `atac_gene_corr_kept` is NULL, overlaps are computed between `gh_tbl` and
#'   the ATAC peaks referenced in `grn_set$fp_annotation`.
#' @param gh_cols Named list mapping columns in `gh_tbl`:
#'   list(chrom="chrom", start="start", end="end", gene="connected_gene",
#'   id=NULL, conf=NULL).
#' @param gene_mode Which gene IDs to retain from RNA when building GH pairs:
#'   "both", "ensembl", or "hgnc".
#' @param fp_score_tbl Optional FP score table (peak_ID + samples). Defaults to
#'   `grn_set$fp_score`.
#' @param rna_tbl Optional RNA table (ensembl_gene_id/HGNC + samples). Defaults
#'   to `grn_set$rna`.
#' @param fdr FDR cutoff applied to FP->gene BH-adjusted p-values.
#' @param r_abs_min Minimum absolute correlation for FP->gene.
#' @param keep_pos If TRUE, keep only positive FP->gene correlations.
#' @param method Correlation method.
#' @param workers Number of parallel workers (future).
#' @param cache_dir Optional cache directory for FP->gene correlations.
#' @param cache_tag Cache tag (subfolder name) under cache_dir.
#' @param cache_chunk_size Chunk size for cached FP->gene computations.
#' @param cache_resume Resume from cache index if TRUE.
#' @param cache_verbose Verbosity for cache operations.
#' @param add_tf_rna_corr If TRUE, compute TF RNA -> target gene RNA correlations
#'   and add r_rna, p_rna, p_adj_rna columns to outputs.
#'
#' @return List with:
#'   - fp_gene_corr_kept: filtered FP->gene correlations (+ optional TF->gene RNA corr cols)
#'   - fp_gene_corr_full: all FP->gene correlations (+ optional TF->gene RNA corr cols)
#' @export
correlate_fp_to_genes <- function(
    grn_set,
    atac_gene_corr_kept = NULL,
    gh_tbl = NULL,
    gh_cols   = list(chrom="chrom", start="start", end="end", gene="connected_gene",
                     id=NULL, conf=NULL),
    gene_mode = c("both","ensembl","hgnc"),
    fp_score_tbl = NULL,
    rna_tbl = NULL,
    fdr       = 0.05,
    r_abs_min = 0.3,
    keep_pos  = FALSE,
    method    = c("pearson","spearman","kendall"),
    workers   = max(1L, round(parallel::detectCores(logical = TRUE) / 2)),
    cache_dir        = NULL,
    cache_tag        = NULL,
    cache_chunk_size = 50000L,
    cache_resume     = TRUE,
    cache_verbose    = FALSE,
    add_tf_rna_corr  = FALSE
) {
  method <- match.arg(method)
  gene_mode <- match.arg(gene_mode)

  # ---- sanity checks ----
  stopifnot(is.list(grn_set), "fp_annotation" %in% names(grn_set))

  fp  <- if (!is.null(fp_score_tbl)) fp_score_tbl else grn_set$fp_score
  rna <- if (!is.null(rna_tbl)) rna_tbl else grn_set$rna
  ann <- grn_set$fp_annotation

  stopifnot(is.data.frame(fp), "peak_ID" %in% names(fp))
  stopifnot(is.data.frame(rna), all(c("ensembl_gene_id","HGNC") %in% names(rna)))
  stopifnot(is.data.frame(ann))

  # Normalize fp_annotation schema across canonical/all pipelines.
  fp_col <- if ("fp_peak" %in% names(ann)) {
    "fp_peak"
  } else if ("peak_ID" %in% names(ann)) {
    "peak_ID"
  } else {
    NA_character_
  }
  atac_col <- if ("atac_peak" %in% names(ann)) "atac_peak" else NA_character_
  tf_col <- if ("tfs" %in% names(ann)) {
    "tfs"
  } else if ("TF" %in% names(ann)) {
    "TF"
  } else {
    NA_character_
  }
  motif_col <- if ("motifs" %in% names(ann)) {
    "motifs"
  } else if ("motif" %in% names(ann)) {
    "motif"
  } else {
    NA_character_
  }
  if (is.na(fp_col) || is.na(atac_col) || is.na(tf_col)) {
    cli::cli_abort(
      paste0(
        "`grn_set$fp_annotation` must contain fp/atac/TF columns. ",
        "Accepted names: fp_peak or peak_ID; atac_peak; tfs or TF."
      )
    )
  }
  ann <- ann |>
    dplyr::mutate(
      fp_peak = as.character(.data[[fp_col]]),
      atac_peak = as.character(.data[[atac_col]]),
      tfs = as.character(.data[[tf_col]]),
      motifs = if (is.na(motif_col)) NA_character_ else as.character(.data[[motif_col]])
    )

  if (!is.null(atac_gene_corr_kept)) {
    stopifnot(is.data.frame(atac_gene_corr_kept))
    stopifnot(all(c("atac_peak","gene_key") %in% names(atac_gene_corr_kept)))
  } else if (is.null(gh_tbl)) {
    cli::cli_abort("Provide either `atac_gene_corr_kept` or `gh_tbl`.")
  }

  # samples present in BOTH fp_score and rna
  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id","HGNC"))
  )
  if (length(sample_cols) < 3L) {
    cli::cli_abort("Not enough shared samples between fp_score and rna.")
  }

  # ---- build gene_expr table keyed by gene_key (supports HGNC or ENSG) ----
  rna_hgnc <- rna |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::transmute(gene_key = HGNC, dplyr::across(dplyr::all_of(sample_cols)))

  rna_ensg <- rna |>
    dplyr::filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") |>
    dplyr::transmute(gene_key = ensembl_gene_id, dplyr::across(dplyr::all_of(sample_cols)))

  # prefer HGNC if both exist
  rna_long <- dplyr::bind_rows(
    rna_hgnc,
    dplyr::anti_join(rna_ensg, rna_hgnc, by = "gene_key")
  )

  # ---- derive atac_peak->gene_key pairs from GH if needed ----
  if (is.null(atac_gene_corr_kept)) {
    req <- c(gh_cols$chrom, gh_cols$start, gh_cols$end, gh_cols$gene)
    if (!all(req %in% names(gh_tbl))) {
      cli::cli_abort("Enhancer table missing required columns: {.val {req}}.")
    }

    ap <- unique(ann$atac_peak)
    ap <- ap[!is.na(ap) & nzchar(ap)]
    if (!length(ap)) {
      cli::cli_abort("No usable atac_peak values in fp_annotation.")
    }

    atac_bed <- .parse_peak_ids(ap) |>
      dplyr::select(atac_chr, atac_start, atac_end, atac_peak)

    have_id   <- !is.null(gh_cols$id)   && gh_cols$id   %in% names(gh_tbl)
    have_conf <- !is.null(gh_cols$conf) && gh_cols$conf %in% names(gh_tbl)
    gh_bed <- gh_tbl |>
      dplyr::transmute(
        gh_chr   = .data[[gh_cols$chrom]],
        gh_start = as.integer(.data[[gh_cols$start]]),
        gh_end   = as.integer(.data[[gh_cols$end]]),
        gene_key = as.character(.data[[gh_cols$gene]]),
        gh_id    = if (have_id)  as.character(.data[[gh_cols$id]])      else NA_character_,
        gh_conf  = if (have_conf) suppressWarnings(as.numeric(.data[[gh_cols$conf]])) else NA_real_
      ) |>
      dplyr::filter(!is.na(gh_start), !is.na(gh_end), !is.na(gene_key), nzchar(gene_key))

    ov_raw <- valr::bed_intersect(
      atac_bed |> dplyr::transmute(chrom = atac_chr, start = atac_start, end = atac_end, atac_peak),
      gh_bed   |> dplyr::transmute(chrom = gh_chr,  start = gh_start,  end = gh_end, gh_id, gh_conf, gene_key)
    )

    nm <- names(ov_raw)
    atk_col  <- if ("atac_peak.x" %in% nm) "atac_peak.x" else "atac_peak"
    ghid_col <- if ("gh_id.y"     %in% nm) "gh_id.y"     else "gh_id"
    ghy_col  <- if ("gh_conf.y"   %in% nm) "gh_conf.y"   else "gh_conf"
    ggn_col  <- if ("gene_key.y"  %in% nm) "gene_key.y"  else "gene_key"

    atac_gene_corr_kept <- ov_raw |>
      dplyr::transmute(
        atac_peak = .data[[atk_col]],
        gene_key  = .data[[ggn_col]],
        gh_id     = .data[[ghid_col]],
        gh_conf   = .data[[ghy_col]]
      ) |>
      dplyr::distinct()

    rna_keys <- .pick_rna_gene_keys(rna, gene_mode)
    atac_gene_corr_kept <- dplyr::semi_join(
      atac_gene_corr_kept,
      tibble::tibble(gene_key = rna_keys),
      by = "gene_key"
    )
  }

  genes_keep <- unique(atac_gene_corr_kept$gene_key)
  rna_long   <- rna_long |>
    dplyr::semi_join(tibble::tibble(gene_key = genes_keep), by = "gene_key")

  # ---- map fp_peaks to those kept atac_peaks & genes (explicit many-to-many) ----
  pairs_tbl <- ann |>
    dplyr::select(fp_peak, atac_peak, tfs, motifs) |>
    dplyr::inner_join(
      atac_gene_corr_kept |>
        dplyr::select(atac_peak, gene_key) |>
        dplyr::distinct(),
      by = "atac_peak",
      relationship = "many-to-many"
    ) |>
    dplyr::distinct(fp_peak, atac_peak, gene_key, tfs, motifs)

  if (nrow(pairs_tbl) == 0L) {
    empty <- tibble::tibble()
    return(list(fp_gene_corr_kept = empty, fp_gene_corr_full = empty))
  }

  # ---- prep matrices (rows = fp_peak / gene_key; cols = shared samples) ----
  # FP
  fp_m <- fp |>
    dplyr::semi_join(tibble::tibble(peak_ID = pairs_tbl$fp_peak), by = "peak_ID") |>
    dplyr::distinct(peak_ID, .keep_all = TRUE) |>
    as.data.frame()
  rownames(fp_m) <- fp_m$peak_ID
  fp_m <- as.matrix(fp_m[, sample_cols, drop = FALSE])

  # RNA (targets)
  rna_m <- rna_long |>
    dplyr::semi_join(tibble::tibble(gene_key = pairs_tbl$gene_key), by = "gene_key") |>
    dplyr::distinct(gene_key, .keep_all = TRUE) |>
    as.data.frame()
  rownames(rna_m) <- rna_m$gene_key
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  # keep only computable pairs
  pairs_tbl <- pairs_tbl |>
    dplyr::filter(fp_peak %in% rownames(fp_m), gene_key %in% rownames(rna_m))

  if (nrow(pairs_tbl) == 0L) {
    empty <- tibble::tibble()
    return(list(fp_gene_corr_kept = empty, fp_gene_corr_full = empty))
  }

  # basic per-pair correlation (used by all FP->gene chunks)
  cor_one <- function(fp_peak, gene_key) {
    x <- fp_m[fp_peak, ]
    y <- rna_m[gene_key, ]
    ok <- stats::complete.cases(x, y)
    n  <- sum(ok)
    if (n < 3L) return(list(n = as.integer(n), r = NA_real_, p = NA_real_))

    if (method == "pearson") {
      r <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
      t <- r * sqrt((n - 2) / (1 - r^2))
      p <- 2 * stats::pt(abs(t), df = n - 2, lower.tail = FALSE)
    } else {
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
      r  <- unname(ct$estimate)
      p  <- ct$p.value
    }
    list(n = as.integer(n), r = as.numeric(r), p = as.numeric(p))
  }

  .empty_corr_tbl <- function() {
    tibble::tibble(
      fp_peak  = character(),
      gene_key = character(),
      n        = integer(),
      r        = double(),
      p        = double()
    )
  }

  # inner correlation engine over a set of unique (fp_peak, gene_key)
  run_corr_for_pairs <- function(pairs_keys_df, inner_workers) {
    if (!nrow(pairs_keys_df)) {
      return(.empty_corr_tbl())
    }

    n_pairs <- nrow(pairs_keys_df)
    nchunk  <- max(1L, inner_workers)

    idx_split <- split(
      seq_len(n_pairs),
      ceiling(seq_along(seq_len(n_pairs)) /
                max(1L, ceiling(n_pairs / max(1L, nchunk))))
    )

    run_chunk <- function(ii) {
      sub <- pairs_keys_df[ii, , drop = FALSE]
      res_list <- vector("list", nrow(sub))
      for (i in seq_len(nrow(sub))) {
        fp_peak_i  <- sub$fp_peak[i]
        gene_key_i <- sub$gene_key[i]
        res_list[[i]] <- cor_one(fp_peak_i, gene_key_i)
      }
      tibble::tibble(
        fp_peak  = sub$fp_peak,
        gene_key = sub$gene_key,
        n        = vapply(res_list, `[[`, integer(1), "n"),
        r        = vapply(res_list, `[[`, numeric(1), "r"),
        p        = vapply(res_list, `[[`, numeric(1), "p")
      )
    }

    if (inner_workers > 1L) {
      out_list <- future.apply::future_lapply(idx_split, run_chunk, future.seed = TRUE)
    } else {
      out_list <- lapply(idx_split, run_chunk)
    }

    dplyr::bind_rows(out_list)
  }

  use_cache     <- !is.null(cache_dir)
  have_parallel <- (workers > 1L)

  # ---------------------------------------------------------------------------
  # Set up future plan (outer level only)
  #
  # UPDATE TO PASTE HERE:
  # - Use multicore on Unix/Linux when supported (lower overhead on HPC).
  # - Fall back to multisession otherwise (Windows, unsupported multicore, etc.).
  # ---------------------------------------------------------------------------
  if (have_parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      cli::cli_abort("Install future or set workers = 1.")
    }
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      cli::cli_abort("Install future.apply or set workers = 1.")
    }

    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    use_multicore <- FALSE
    if (isTRUE(future::supportsMulticore())) {
      # supportsMulticore() is already OS-aware (FALSE on Windows)
      use_multicore <- TRUE
    }

    if (use_multicore) {
      future::plan(future::multicore, workers = workers)
      if (cache_verbose) {
        cli::cli_inform("future plan: multicore ({workers} workers).")
      }
    } else {
      future::plan(future::multisession, workers = workers)
      if (cache_verbose) {
        cli::cli_inform("future plan: multisession ({workers} workers).")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # FP->gene correlations (cached or not)
  # ---------------------------------------------------------------------------
  if (!use_cache) {
    pairs_keys <- pairs_tbl |>
      dplyr::select(fp_peak, gene_key) |>
      dplyr::distinct()

    inner_workers <- if (have_parallel) workers else 1L
    if (cache_verbose) {
      cli::cli_inform(
        "Running FP-gene correlations without cache: {nrow(pairs_keys)} unique pairs with {inner_workers} workers."
      )
    }
    corr_pairs <- run_corr_for_pairs(pairs_keys, inner_workers)
  } else {
    cache_dir <- normalizePath(cache_dir, mustWork = FALSE)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

    if (is.null(cache_tag) || !nzchar(cache_tag)) {
      cache_tag <- paste0("fp_genes_", method)
    }

    cache_root <- file.path(cache_dir, cache_tag)
    if (!dir.exists(cache_root)) dir.create(cache_root, recursive = TRUE)

    meta_file  <- file.path(cache_root, "meta.rds")
    index_file <- file.path(cache_root, "index.rds")

    meta_new <- list(method = method)

    if (file.exists(meta_file)) {
      meta_old <- tryCatch(readRDS(meta_file), error = function(e) NULL)
      if (!is.null(meta_old)) {
        if (!identical(meta_old$method, method)) {
          cli::cli_abort(
            c(
              "Existing FP cache at {.path {cache_root}} is not compatible.",
              "i" = "It was created with different correlation method.",
              "i" = "Use a different cache_tag or clean the cache directory."
            )
          )
        }
      }
    } else {
      saveRDS(meta_new, meta_file)
    }

    if (cache_verbose) {
      cli::cli_inform("Using FP-gene correlation cache at {.path {cache_root}}.")
    }

    if (!cache_resume || !file.exists(index_file)) {
      cache_index <- tibble::tibble(
        fp_peak  = character(),
        gene_key = character(),
        file     = character(),
        chunk_id = integer()
      )
    } else {
      cache_index <- tryCatch(readRDS(index_file), error = function(e) NULL)
      if (!is.data.frame(cache_index) ||
          !all(c("fp_peak","gene_key","file","chunk_id") %in% names(cache_index))) {
        cli::cli_warn("FP cache index at {.path {index_file}} is malformed; restarting index.")
        cache_index <- tibble::tibble(
          fp_peak  = character(),
          gene_key = character(),
          file     = character(),
          chunk_id = integer()
        )
      }
    }

    # Recover/augment index from on-disk chunk files (useful after interrupted runs).
    if (isTRUE(cache_resume)) {
      chunk_files <- list.files(
        cache_root,
        pattern = "^chunk_[0-9]+\\.csv\\.gz$",
        full.names = TRUE
      )
      if (length(chunk_files)) {
        files_missing_from_index <- setdiff(chunk_files, unique(cache_index$file))
        if (length(files_missing_from_index)) {
          if (cache_verbose) {
            cli::cli_inform(
              "Recovering FP cache index from {length(files_missing_from_index)} unindexed chunk file(s)."
            )
          }
          recovered_idx <- lapply(files_missing_from_index, function(ff) {
            id_chr <- sub("^chunk_([0-9]+)\\.csv\\.gz$", "\\1", basename(ff))
            chunk_id <- suppressWarnings(as.integer(id_chr))
            df <- tryCatch(
              readr::read_csv(
                ff,
                col_select = c("fp_peak", "gene_key"),
                show_col_types = FALSE,
                progress = FALSE
              ),
              error = function(e) NULL
            )
            if (is.null(df) || !all(c("fp_peak", "gene_key") %in% names(df))) {
              return(tibble::tibble(
                fp_peak = character(),
                gene_key = character(),
                file = character(),
                chunk_id = integer()
              ))
            }
            tibble::tibble(
              fp_peak = as.character(df$fp_peak),
              gene_key = as.character(df$gene_key),
              file = ff,
              chunk_id = chunk_id
            )
          })
          recovered_idx <- dplyr::bind_rows(recovered_idx)
          if (nrow(recovered_idx)) {
            cache_index <- dplyr::bind_rows(cache_index, recovered_idx) |>
              dplyr::arrange(fp_peak, gene_key) |>
              dplyr::distinct(fp_peak, gene_key, .keep_all = TRUE)
            saveRDS(cache_index, index_file)
            if (cache_verbose) {
              cli::cli_inform(
                "Recovered FP cache index now tracks {nrow(cache_index)} pairs across {length(unique(cache_index$chunk_id))} chunks."
              )
            }
          }
        }
      }
    }

    run_keys <- pairs_tbl |>
      dplyr::select(fp_peak, gene_key) |>
      dplyr::distinct()

    total_pairs      <- nrow(run_keys)
    existing_n       <- 0L
    existing_chunks  <- 0L

    if (!nrow(cache_index)) {
      have_idx     <- tibble::tibble(fp_peak = character(), gene_key = character(), file = character(), chunk_id = integer())
      missing_keys <- run_keys
    } else {
      have_idx <- dplyr::inner_join(run_keys, cache_index, by = c("fp_peak","gene_key"))
      missing_keys <- dplyr::anti_join(run_keys, cache_index, by = c("fp_peak","gene_key"))
      existing_n <- nrow(have_idx)
      existing_chunks <- length(unique(cache_index$chunk_id))
    }

    missing_n <- nrow(missing_keys)

    cli::cli_inform("FP-gene correlation run: {total_pairs} unique (fp_peak, gene_key) pairs requested.")
    cli::cli_inform("FP cache index currently tracks {nrow(cache_index)} pairs across {existing_chunks} chunks.")
    cli::cli_inform("{existing_n} pairs already satisfy this run; {missing_n} pairs remain to compute.")

    corr_cached <- .empty_corr_tbl()

    # ---- load cached correlations for pairs we already have ----
    if (existing_n > 0L) {
      files_used <- unique(have_idx$file)
      if (cache_verbose) {
        cli::cli_inform("Loading cached FP correlations from {length(files_used)} chunk files.")
      }
      for (ff in files_used) {
        idx_ff <- have_idx[have_idx$file == ff, c("fp_peak","gene_key"), drop = FALSE]
        if (!file.exists(ff)) {
          cli::cli_warn("FP cache chunk file missing: {.path {ff}}; those keys will be recomputed.")
          missing_keys <- dplyr::bind_rows(missing_keys, idx_ff) |>
            dplyr::distinct()
          next
        }
        df_ff <- readr::read_csv(ff, show_col_types = FALSE, progress = FALSE)
        df_sub <- dplyr::inner_join(df_ff, idx_ff, by = c("fp_peak","gene_key"))
        if (nrow(df_sub)) {
          corr_cached <- dplyr::bind_rows(corr_cached, df_sub)
        }
      }
    }

    # ---- compute correlations for missing pairs, in chunks ----
    corr_new <- .empty_corr_tbl()

    if (missing_n > 0L) {
      cache_chunk_size <- as.integer(cache_chunk_size)
      if (is.na(cache_chunk_size) || cache_chunk_size < 1L) cache_chunk_size <- 50000L

      nchunk_new <- ceiling(missing_n / cache_chunk_size)

      cli::cli_inform(
        "Computing {missing_n} new FP-gene pairs in {nchunk_new} chunks (target chunk size ~{cache_chunk_size})."
      )
      if (cache_verbose) {
        cli::cli_inform("Requested workers: {workers}. With {nchunk_new} chunks, max parallelism is min(workers, nchunk_new).")
      }

      old_chunks <- if (nrow(cache_index)) cache_index$chunk_id else integer(0)
      start_id   <- if (length(old_chunks)) max(old_chunks, na.rm = TRUE) + 1L else 1L
      chunk_vec  <- rep(
        seq_len(nchunk_new),
        each = cache_chunk_size,
        length.out = missing_n
      ) + (start_id - 1L)

      missing_keys$chunk_id <- as.integer(chunk_vec)
      split_chunks <- split(missing_keys, missing_keys$chunk_id)

      chunk_fun <- function(keys_chunk) {
        chunk_id <- unique(keys_chunk$chunk_id)
        if (length(chunk_id) != 1L) {
          cli::cli_warn("FP chunk has multiple chunk_id values; using the first.")
          chunk_id <- chunk_id[1L]
        }
        chunk_id <- as.integer(chunk_id[1L])

        keys_unique <- keys_chunk[, c("fp_peak","gene_key"), drop = FALSE]
        keys_unique <- dplyr::distinct(keys_unique)

        if (cache_verbose) {
          cli::cli_inform("Worker computing FP chunk {chunk_id} with {nrow(keys_unique)} unique pairs.")
        }

        corr_chunk <- run_corr_for_pairs(
          pairs_keys_df = keys_unique,
          inner_workers = 1L
        )

        if (!nrow(corr_chunk)) {
          return(list(
            index = tibble::tibble(
              fp_peak  = character(),
              gene_key = character(),
              file     = character(),
              chunk_id = integer()
            ),
            corr = .empty_corr_tbl()
          ))
        }

        chunk_file <- file.path(cache_root, sprintf("chunk_%05d.csv.gz", chunk_id))
        readr::write_csv(corr_chunk, chunk_file)

        index_chunk <- tibble::tibble(
          fp_peak  = corr_chunk$fp_peak,
          gene_key = corr_chunk$gene_key,
          file     = chunk_file,
          chunk_id = chunk_id
        )

        list(index = index_chunk, corr = corr_chunk)
      }

      if (have_parallel) {
        res_list <- future.apply::future_lapply(split_chunks, chunk_fun, future.seed = TRUE)
      } else {
        res_list <- lapply(split_chunks, chunk_fun)
      }

      index_new <- dplyr::bind_rows(lapply(res_list, function(x) x$index))
      corr_new  <- dplyr::bind_rows(lapply(res_list, function(x) x$corr))

      cache_index <- dplyr::bind_rows(cache_index, index_new) |>
        dplyr::arrange(fp_peak, gene_key) |>
        dplyr::distinct(fp_peak, gene_key, .keep_all = TRUE)

      saveRDS(cache_index, index_file)
      if (cache_verbose) {
        cli::cli_inform(
          "Updated FP cache index now tracks {nrow(cache_index)} pairs across {length(unique(cache_index$chunk_id))} chunks."
        )
      }
    } else if (cache_verbose) {
      cli::cli_inform("No new FP pairs to compute; all requested correlations already present in cache.")
    }

    corr_pairs <- dplyr::bind_rows(corr_cached, corr_new) |>
      dplyr::arrange(fp_peak, gene_key) |>
      dplyr::distinct(fp_peak, gene_key, .keep_all = TRUE)
  }

  # ---------------------------------------------------------------------------
  # Inflate back to per-(fp_peak, atac_peak, gene_key, tfs, motifs) rows
  # ---------------------------------------------------------------------------
  pairs_tbl2 <- pairs_tbl |>
    dplyr::inner_join(corr_pairs, by = c("fp_peak","gene_key"))

  # drop pairs without valid p
  pairs_tbl2 <- pairs_tbl2 |>
    dplyr::filter(!is.na(p))

  if (!nrow(pairs_tbl2)) {
    empty <- tibble::tibble()
    return(list(fp_gene_corr_kept = empty, fp_gene_corr_full = empty))
  }

  # BH adjust across all fp_gene rows (matches original behavior)
  pairs_tbl2$p_adj <- stats::p.adjust(pairs_tbl2$p, method = "BH")

  fp_gene_corr_full <- pairs_tbl2 |>
    dplyr::transmute(
      fp_peak,
      gene_key,
      atac_peak,
      tfs,
      motifs,
      n_fp     = as.integer(n),
      r_fp     = as.numeric(r),
      p_fp     = as.numeric(p),
      p_adj_fp = as.numeric(p_adj)
    )

  # ---------------------------------------------------------------------------
  # NEW: TF RNA -> target gene RNA correlations (per unique tfs x gene_key)
  # ---------------------------------------------------------------------------
  if (isTRUE(add_tf_rna_corr)) {
    tf_gene_keys <- fp_gene_corr_full |>
      dplyr::select(tfs, gene_key) |>
      dplyr::distinct()

    rna_needed <- unique(c(tf_gene_keys$tfs, tf_gene_keys$gene_key))
    rna_needed <- rna_needed[!is.na(rna_needed) & nzchar(rna_needed)]

    rna_all <- grn_set$rna
    rna_hgnc_all <- rna_all |>
      dplyr::filter(!is.na(HGNC), HGNC != "") |>
      dplyr::transmute(gene_key = HGNC, dplyr::across(dplyr::all_of(sample_cols)))
    rna_ensg_all <- rna_all |>
      dplyr::filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") |>
      dplyr::transmute(gene_key = ensembl_gene_id, dplyr::across(dplyr::all_of(sample_cols)))
    rna_long_all <- dplyr::bind_rows(
      rna_hgnc_all,
      dplyr::anti_join(rna_ensg_all, rna_hgnc_all, by = "gene_key")
    )

    rna_long_all <- rna_long_all |>
      dplyr::semi_join(tibble::tibble(gene_key = rna_needed), by = "gene_key") |>
      dplyr::distinct(gene_key, .keep_all = TRUE) |>
      as.data.frame()

    if (nrow(rna_long_all) == 0L) {
      tf_gene_corr <- tibble::tibble(
        tfs       = character(),
        gene_key  = character(),
        r_rna     = double(),
        p_rna     = double(),
        p_adj_rna = double()
      )
    } else {
      rownames(rna_long_all) <- rna_long_all$gene_key
      rna_m_all <- as.matrix(rna_long_all[, sample_cols, drop = FALSE])

      tf_split <- split(tf_gene_keys$gene_key, tf_gene_keys$tfs)

      tf_worker <- function(tf) {
        genes <- unique(tf_split[[tf]])
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (!length(genes)) {
          return(tibble::tibble(tfs = character(), gene_key = character(), r_rna = double(), p_rna = double()))
        }

        if (!(tf %in% rownames(rna_m_all))) {
          return(tibble::tibble(
            tfs      = tf,
            gene_key = genes,
            r_rna    = NA_real_,
            p_rna    = NA_real_
          ))
        }

        genes_ok <- genes[genes %in% rownames(rna_m_all)]
        genes_miss <- setdiff(genes, genes_ok)

        if (!length(genes_ok)) {
          return(tibble::tibble(
            tfs      = tf,
            gene_key = genes,
            r_rna    = NA_real_,
            p_rna    = NA_real_
          ))
        }

        x <- as.numeric(rna_m_all[tf, ])
        Y <- rna_m_all[genes_ok, , drop = FALSE]

        fast_ok <- all(is.finite(x)) && all(is.finite(Y))
        if (identical(method, "pearson") && fast_ok) {
          r <- suppressWarnings(stats::cor(t(Y), x, method = "pearson"))
          r <- as.numeric(r)

          n <- length(x)
          tstat <- r * sqrt((n - 2) / (1 - r^2))
          p <- 2 * stats::pt(abs(tstat), df = n - 2, lower.tail = FALSE)

          out_ok <- tibble::tibble(
            tfs      = tf,
            gene_key = genes_ok,
            r_rna    = as.numeric(r),
            p_rna    = as.numeric(p)
          )
        } else {
          r_vec <- rep(NA_real_, length(genes_ok))
          p_vec <- rep(NA_real_, length(genes_ok))

          for (i in seq_along(genes_ok)) {
            g <- genes_ok[i]
            y <- as.numeric(rna_m_all[g, ])
            ok <- stats::complete.cases(x, y)
            n  <- sum(ok)
            if (n < 3L) next

            if (method == "pearson") {
              rr <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
              tt <- rr * sqrt((n - 2) / (1 - rr^2))
              pp <- 2 * stats::pt(abs(tt), df = n - 2, lower.tail = FALSE)
              r_vec[i] <- as.numeric(rr)
              p_vec[i] <- as.numeric(pp)
            } else {
              ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
              r_vec[i] <- as.numeric(unname(ct$estimate))
              p_vec[i] <- as.numeric(ct$p.value)
            }
          }

          out_ok <- tibble::tibble(
            tfs      = tf,
            gene_key = genes_ok,
            r_rna    = r_vec,
            p_rna    = p_vec
          )
        }

        if (length(genes_miss)) {
          out_miss <- tibble::tibble(
            tfs      = tf,
            gene_key = genes_miss,
            r_rna    = NA_real_,
            p_rna    = NA_real_
          )
          dplyr::bind_rows(out_ok, out_miss)
        } else {
          out_ok
        }
      }

      tf_list <- names(tf_split)

      if (cache_verbose) {
        cli::cli_inform("TF-RNA corr requested: {nrow(tf_gene_keys)} (tfs,gene_key) pairs across {length(tf_list)} TFs.")
        if (have_parallel) cli::cli_inform("TF-RNA corr will use future workers: {future::nbrOfWorkers()}.")
      }

      if (have_parallel && length(tf_list) > 1L) {
        tf_out <- future.apply::future_lapply(tf_list, tf_worker, future.seed = TRUE)
      } else {
        tf_out <- lapply(tf_list, tf_worker)
      }

      tf_gene_corr0 <- dplyr::bind_rows(tf_out) |>
        dplyr::distinct(tfs, gene_key, .keep_all = TRUE)

      tf_gene_corr0$p_adj_rna <- stats::p.adjust(tf_gene_corr0$p_rna, method = "BH")

      tf_gene_corr <- tf_gene_corr0 |>
        dplyr::transmute(
          tfs,
          gene_key,
          r_rna     = as.numeric(r_rna),
          p_rna     = as.numeric(p_rna),
          p_adj_rna = as.numeric(p_adj_rna)
        )
    }

    fp_gene_corr_full <- fp_gene_corr_full |>
      dplyr::left_join(tf_gene_corr, by = c("tfs", "gene_key"))
  }

  fp_gene_corr_kept <- fp_gene_corr_full |>
    dplyr::filter(
      p_adj_fp < fdr,
      abs(r_fp) >= r_abs_min,
      if (keep_pos) r_fp > 0 else TRUE
    )

  list(
    fp_gene_corr_kept = fp_gene_corr_kept,
    fp_gene_corr_full = fp_gene_corr_full
  )
}

#' @export
add_tf_rna_corr_cols <- function(fp_tbl,
                                 rna_tbl,
                                 method = c("pearson", "spearman", "kendall"),
                                 cores  = max(1L, parallel::detectCores(logical = TRUE) - 1L),
                                 verbose = FALSE) {
  method <- match.arg(method)

  stopifnot(is.data.frame(fp_tbl), all(c("tfs", "gene_key") %in% names(fp_tbl)))
  stopifnot(is.data.frame(rna_tbl), all(c("ensembl_gene_id", "HGNC") %in% names(rna_tbl)))

  sample_cols <- setdiff(names(rna_tbl), c("ensembl_gene_id", "HGNC"))
  if (length(sample_cols) < 3L) stop("Not enough RNA sample columns.")
  cores <- as.integer(cores)
  if (is.na(cores) || cores < 1L) cores <- 1L

  # ---- UNIQUE PAIRS ONLY ----
  pairs_df <- unique(fp_tbl[, c("tfs", "gene_key"), drop = FALSE])
  pairs_df <- pairs_df[!is.na(pairs_df$tfs) & nzchar(pairs_df$tfs) &
                         !is.na(pairs_df$gene_key) & nzchar(pairs_df$gene_key), , drop = FALSE]

  if (!nrow(pairs_df)) {
    out <- fp_tbl
    out$r_rna <- NA_real_; out$p_rna <- NA_real_; out$p_adj_rna <- NA_real_
    return(out)
  }

  needed_keys <- unique(c(pairs_df$tfs, pairs_df$gene_key))
  needed_keys <- needed_keys[!is.na(needed_keys) & nzchar(needed_keys)]

  # Build RNA keyed by HGNC first, then ENSG (prefer HGNC)
  hgnc_ok <- !is.na(rna_tbl$HGNC) & nzchar(rna_tbl$HGNC)
  ensg_ok <- !is.na(rna_tbl$ensembl_gene_id) & nzchar(rna_tbl$ensembl_gene_id)

  rna_hgnc <- rna_tbl[hgnc_ok, c("HGNC", sample_cols), drop = FALSE]
  names(rna_hgnc)[1] <- "gene_key"

  rna_ensg <- rna_tbl[ensg_ok, c("ensembl_gene_id", sample_cols), drop = FALSE]
  names(rna_ensg)[1] <- "gene_key"

  rna_long <- rbind(rna_hgnc, rna_ensg[!(rna_ensg$gene_key %in% rna_hgnc$gene_key), , drop = FALSE])
  rna_long <- rna_long[rna_long$gene_key %in% needed_keys, , drop = FALSE]
  rna_long <- rna_long[!duplicated(rna_long$gene_key), , drop = FALSE]

  if (!nrow(rna_long)) {
    out <- fp_tbl
    out$r_rna <- NA_real_; out$p_rna <- NA_real_; out$p_adj_rna <- NA_real_
    return(out)
  }

  rn <- rna_long$gene_key
  rna_m <- as.matrix(rna_long[, sample_cols, drop = FALSE])
  storage.mode(rna_m) <- "double"
  rownames(rna_m) <- rn

  # split UNIQUE pairs by TF (each (tf,gene) appears once here)
  tf_split <- split(pairs_df$gene_key, pairs_df$tfs)
  tf_list  <- names(tf_split)

  tf_worker <- function(tf) {
    genes <- tf_split[[tf]]
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (!length(genes)) return(NULL)

    if (!(tf %in% rownames(rna_m))) {
      return(data.frame(tfs=tf, gene_key=genes, r_rna=NA_real_, p_rna=NA_real_, stringsAsFactors=FALSE))
    }

    x_all <- as.numeric(rna_m[tf, ])
    r_out <- rep(NA_real_, length(genes))
    p_out <- rep(NA_real_, length(genes))

    for (i in seq_along(genes)) {
      g <- genes[i]
      if (!(g %in% rownames(rna_m))) next
      y_all <- as.numeric(rna_m[g, ])
      ok <- stats::complete.cases(x_all, y_all)
      n <- sum(ok)
      if (n < 3L) next

      if (method == "pearson") {
        r <- suppressWarnings(stats::cor(x_all[ok], y_all[ok], method="pearson"))
        tstat <- r * sqrt((n - 2) / (1 - r^2))
        p <- 2 * stats::pt(abs(tstat), df=n - 2, lower.tail=FALSE)
      } else {
        ct <- suppressWarnings(stats::cor.test(x_all[ok], y_all[ok], method=method, exact=FALSE))
        r <- as.numeric(unname(ct$estimate))
        p <- as.numeric(ct$p.value)
      }
      r_out[i] <- as.numeric(r)
      p_out[i] <- as.numeric(p)
    }

    data.frame(tfs=tf, gene_key=genes, r_rna=r_out, p_rna=p_out, stringsAsFactors=FALSE)
  }

  if (isTRUE(verbose)) {
    .log_inform(
      "TF->RNA correlations: {n_pairs} unique pairs; {n_tfs} TFs; cores={cores}.",
      n_pairs = nrow(pairs_df),
      n_tfs = length(tf_list),
      cores = cores
    )
  }

  use_mclapply <- (.Platform$OS.type != "windows") && cores > 1L
  if (use_mclapply) {
    res_list <- parallel::mclapply(tf_list, tf_worker, mc.cores = cores)
  } else if (cores > 1L) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("tf_split", "tf_worker", "rna_m", "method"), envir = environment())
    res_list <- parallel::parLapply(cl, tf_list, tf_worker)
  } else {
    res_list <- lapply(tf_list, tf_worker)
  }

  tf_gene_corr <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (!nrow(tf_gene_corr)) {
    out <- fp_tbl
    out$r_rna <- NA_real_; out$p_rna <- NA_real_; out$p_adj_rna <- NA_real_
    return(out)
  }

  # ensure UNIQUE keys (defensive)
  key_corr <- paste(tf_gene_corr$tfs, tf_gene_corr$gene_key, sep="|")
  keep <- !duplicated(key_corr)
  tf_gene_corr <- tf_gene_corr[keep, , drop = FALSE]

  tf_gene_corr$p_adj_rna <- NA_real_
  okp <- is.finite(tf_gene_corr$p_rna)
  tf_gene_corr$p_adj_rna[okp] <- stats::p.adjust(tf_gene_corr$p_rna[okp], method="BH")

  # ---- attach back to fp_tbl (duplicates inherit) ----
  key_map <- paste(tf_gene_corr$tfs, tf_gene_corr$gene_key, sep="|")
  key_fp  <- paste(fp_tbl$tfs, fp_tbl$gene_key, sep="|")
  m <- match(key_fp, key_map)

  out <- fp_tbl
  out$r_rna <- tf_gene_corr$r_rna[m]
  out$p_rna <- tf_gene_corr$p_rna[m]
  out$p_adj_rna <- tf_gene_corr$p_adj_rna[m]
  out
}

#' Recompute corr_fp_tf_* for fp_annotation using fp_score and RNA
#'
#' @param grn_set List with fp_annotation, fp_score, and rna.
#' @param method Correlation method.
#' @param cores Number of cores for mclapply on Unix.
#' @param chunk_size Pair chunk size.
#' @param min_non_na Minimum paired non-NA values required to test.
#' @return fp_annotation with corr_fp_tf_r/p/p_adj.
#' @export
make_fp_annotation_corr <- function(grn_set,
                                    method = c("pearson", "spearman", "kendall"),
                                    cores = 10L,
                                    chunk_size = 5000L,
                                    min_non_na = 3L) {
  method <- match.arg(method)

  ann <- grn_set$fp_annotation
  fp  <- grn_set$fp_score
  rna <- grn_set$rna

  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id", "HGNC"))
  )
  min_non_na <- as.integer(min_non_na)
  if (length(sample_cols) < min_non_na) {
    stop("Not enough shared samples between fp_score and rna.")
  }

  fp_m <- as.matrix(fp[, sample_cols, drop = FALSE])
  rownames(fp_m) <- fp$peak_ID

  rna_m <- rna |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::select(HGNC, dplyr::all_of(sample_cols)) |>
    dplyr::distinct(HGNC, .keep_all = TRUE) |>
    as.data.frame()
  rownames(rna_m) <- rna_m$HGNC
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  # Optional gating: only correlate on conditions where
  # (i) footprint is bound and (ii) TF is expressed.
  use_gate <- FALSE
  fp_bound_m <- NULL
  rna_expr_m <- NULL
  if (is.data.frame(grn_set$fp_bound_condition) && is.data.frame(grn_set$rna_expressed)) {
    cond_gate <- intersect(
      sample_cols,
      intersect(
        setdiff(names(grn_set$fp_bound_condition), "peak_ID"),
        setdiff(names(grn_set$rna_expressed), c("ensembl_gene_id", "HGNC"))
      )
    )
    if (length(cond_gate)) {
      fpb <- grn_set$fp_bound_condition[, c("peak_ID", cond_gate), drop = FALSE]
      if (anyDuplicated(fpb$peak_ID)) {
        fpb <- fpb[!duplicated(fpb$peak_ID), , drop = FALSE]
      }
      fp_bound_m <- as.matrix(fpb[, cond_gate, drop = FALSE])
      rownames(fp_bound_m) <- fpb$peak_ID

      rxe <- grn_set$rna_expressed |>
        dplyr::filter(!is.na(HGNC), HGNC != "") |>
        dplyr::select(HGNC, dplyr::all_of(cond_gate)) |>
        dplyr::distinct(HGNC, .keep_all = TRUE) |>
        as.data.frame()
      rownames(rxe) <- rxe$HGNC
      rna_expr_m <- as.matrix(rxe[, cond_gate, drop = FALSE])
      use_gate <- TRUE
    }
  }

  fp_id <- match(ann$fp_peak, rownames(fp_m))
  tf_id <- match(ann$tfs, rownames(rna_m))

  pairs_u <- unique(data.frame(fp_id = fp_id, tf_id = tf_id))
  pairs_u <- pairs_u[!is.na(pairs_u$fp_id) & !is.na(pairs_u$tf_id), , drop = FALSE]
  if (!nrow(pairs_u)) {
    out <- ann
    out$corr_fp_tf_r <- NA_real_
    out$corr_fp_tf_p <- NA_real_
    out$corr_fp_tf_p_adj <- NA_real_
    return(out)
  }

  idx <- split(seq_len(nrow(pairs_u)),
               ceiling(seq_len(nrow(pairs_u)) / chunk_size))

  worker <- function(ii) {
    sub <- pairs_u[ii, , drop = FALSE]
    r <- numeric(nrow(sub))
    p <- numeric(nrow(sub))

    for (i in seq_len(nrow(sub))) {
      x <- fp_m[sub$fp_id[i], ]
      y <- rna_m[sub$tf_id[i], ]
      ok <- is.finite(x) & is.finite(y)
      if (isTRUE(use_gate)) {
        fp_key <- rownames(fp_m)[sub$fp_id[i]]
        tf_key <- rownames(rna_m)[sub$tf_id[i]]
        fp_gate <- rep(FALSE, ncol(fp_m))
        tf_gate <- rep(FALSE, ncol(fp_m))
        fp_gate_idx <- match(colnames(fp_m), colnames(fp_bound_m))
        tf_gate_idx <- match(colnames(fp_m), colnames(rna_expr_m))
        fp_row <- match(fp_key, rownames(fp_bound_m))
        tf_row <- match(tf_key, rownames(rna_expr_m))
        if (is.finite(fp_row) && !is.na(fp_row) && any(is.finite(fp_gate_idx))) {
          idx <- which(is.finite(fp_gate_idx))
          fp_gate[idx] <- as.numeric(fp_bound_m[fp_row, fp_gate_idx[idx]]) > 0
        }
        if (is.finite(tf_row) && !is.na(tf_row) && any(is.finite(tf_gate_idx))) {
          idx <- which(is.finite(tf_gate_idx))
          tf_gate[idx] <- as.numeric(rna_expr_m[tf_row, tf_gate_idx[idx]]) > 0
        }
        ok <- ok & fp_gate & tf_gate
      }
      n  <- sum(ok)
      if (n < min_non_na) {
        r[i] <- NA_real_
        p[i] <- NA_real_
        next
      }
      ct <- if (method == "pearson") {
        suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
      } else {
        suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
      }
      r[i] <- as.numeric(unname(ct$estimate))
      p[i] <- as.numeric(ct$p.value)
    }
    data.frame(fp_id = sub$fp_id, tf_id = sub$tf_id, r = r, p = p)
  }

  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)
  corr_tbl <- do.call(rbind, res_list)

  key_u <- paste(corr_tbl$fp_id, corr_tbl$tf_id, sep = "_")
  key_a <- paste(fp_id, tf_id, sep = "_")
  m <- match(key_a, key_u)

  out <- ann
  out$corr_fp_tf_r <- corr_tbl$r[m]
  out$corr_fp_tf_p <- corr_tbl$p[m]
  out$corr_fp_tf_p_adj <- stats::p.adjust(out$corr_fp_tf_p, method = "BH")
  out
}

.row_stats_mean_var <- function(m) {
  m <- as.matrix(m)
  ok <- is.finite(m)
  m[!ok] <- 0
  n <- rowSums(ok)
  sum_x <- rowSums(m)
  sum_x2 <- rowSums(m^2)
  mean_x <- sum_x / n
  var_x <- (sum_x2 - (sum_x^2) / n) / pmax(1, n - 1)
  var_x[n < 2L] <- NA_real_
  list(mean = mean_x, var = var_x)
}

.row_stats_mean_var_chunked <- function(m, cores = 1L, chunk_size = 50000L) {
  m <- as.matrix(m)
  if (cores <= 1L || nrow(m) <= chunk_size) return(.row_stats_mean_var(m))

  idx <- split(seq_len(nrow(m)), ceiling(seq_len(nrow(m)) / chunk_size))
  worker <- function(ii) .row_stats_mean_var(m[ii, , drop = FALSE])
  use_mc <- (.Platform$OS.type != "windows") && cores > 1L
  res_list <- if (use_mc) parallel::mclapply(idx, worker, mc.cores = cores) else lapply(idx, worker)

  mean_x <- numeric(nrow(m))
  var_x <- numeric(nrow(m))
  for (k in seq_along(idx)) {
    i <- idx[[k]]
    mean_x[i] <- res_list[[k]]$mean
    var_x[i] <- res_list[[k]]$var
  }
  list(mean = mean_x, var = var_x)
}

compute_hv_variance_tbl <- function(
    tbl,
    id_cols,
    sample_cols = NULL,
    cores = 1L,
    chunk_size = 50000L
) {
  stopifnot(is.data.frame(tbl))
  if (is.null(sample_cols)) {
    drop_cols <- unique(c(id_cols, "HGNC", "ensembl_gene_id"))
    sample_cols <- setdiff(names(tbl), drop_cols)
  }
  if (!length(sample_cols)) stop("No sample columns found for variance computation.")

  m <- as.matrix(tbl[, sample_cols, drop = FALSE])
  stats_raw <- .row_stats_mean_var_chunked(m, cores = cores, chunk_size = chunk_size)

  out <- tibble::tibble(
    !!!rlang::set_names(lapply(id_cols, function(nm) tbl[[nm]]), id_cols),
    mean_raw = stats_raw$mean,
    var_raw = stats_raw$var
  )

  rsd <- sqrt(stats_raw$var) / stats_raw$mean
  rsd[!is.finite(rsd)] <- NA_real_
  out$rsd <- rsd

  out
}

#' Precompute footprint and gene variance
#'
#' @param grn_set List containing \code{$fp_score} and \code{$rna}.
#' @param cores Number of CPU cores to use.
#' @param chunk_size Row chunk size for parallel processing.
#' @return A list with \code{fp_variance} and \code{rna_variance} tables.
#' @export
precompute_hvf_hvg_variance <- function(
    grn_set,
    cores = 1L,
    chunk_size = 50000L
) {
  if (!is.list(grn_set) || !all(c("fp_score", "rna") %in% names(grn_set))) {
    stop("grn_set must contain $fp_score and $rna")
  }
  list(
    fp_variance = compute_hv_variance_tbl(
      tbl = grn_set$fp_score,
      id_cols = "peak_ID",
      cores = cores,
      chunk_size = chunk_size
    ),
    rna_variance = compute_hv_variance_tbl(
      tbl = grn_set$rna,
      id_cols = c("ensembl_gene_id", "HGNC"),
      cores = cores,
      chunk_size = chunk_size
    )
  )
}

#' Filter FP->gene correlations using TF annotation correlation thresholds
#'
#' @param fp_gene_corr_full fp_gene_corr_full table with fp_peak and tfs columns.
#' @param fp_annotation fp_annotation table with corr_fp_tf_r and corr_fp_tf_p_adj.
#' @param r_thr Minimum corr_fp_tf_r to keep.
#' @param p_adj_thr Maximum corr_fp_tf_p_adj to keep.
#' @return Filtered fp_gene_corr_full table.
#' @export
filter_fp_gene_corr_by_tf_annotation <- function(
    fp_gene_corr_full,
    fp_annotation,
    r_thr = 0.3,
    p_adj_thr = 0.05
) {
  stopifnot(is.data.frame(fp_gene_corr_full), is.data.frame(fp_annotation))
  need_fp <- c("fp_peak", "tfs")
  need_ann <- c("fp_peak", "tfs", "corr_fp_tf_r", "corr_fp_tf_p_adj")
  if (!all(need_fp %in% names(fp_gene_corr_full))) {
    stop("fp_gene_corr_full is missing required columns: fp_peak, tfs.")
  }
  if (!all(need_ann %in% names(fp_annotation))) {
    stop("fp_annotation is missing required columns for filtering.")
  }

  ann_keep <- fp_annotation |>
    dplyr::filter(.data$corr_fp_tf_r > r_thr, .data$corr_fp_tf_p_adj < p_adj_thr) |>
    dplyr::distinct(.data$fp_peak, .data$tfs)

  if (!nrow(ann_keep)) {
    return(fp_gene_corr_full[0, , drop = FALSE])
  }

  fp_gene_corr_full |>
    dplyr::semi_join(ann_keep, by = c("fp_peak", "tfs"))
}

#' Merge FP->gene correlation tables across methods and attach TF->gene RNA corr
#'
#' @param fp_pearson fp_gene_corr_full (pearson).
#' @param fp_spearman fp_gene_corr_full (spearman).
#' @param rna_tbl RNA table for TF->gene correlation.
#' @param rna_method Correlation method for TF->gene RNA (default pearson).
#' @param rna_cores Cores for TF->gene RNA correlation.
#' @return Tibble with TF/gene/peak keys and method-specific columns.
#' @export
combine_fp_gene_corr_methods <- function(
    fp_pearson,
    fp_spearman,
    rna_tbl,
    rna_method = c("pearson", "spearman", "kendall"),
    rna_cores = max(1L, parallel::detectCores(logical = TRUE) - 1L)
) {
  rna_method <- match.arg(rna_method)
  stopifnot(is.data.frame(fp_pearson), is.data.frame(fp_spearman))
  need_cols <- c("fp_peak", "gene_key", "atac_peak", "tfs", "r_fp", "p_fp", "p_adj_fp")
  if (!all(need_cols %in% names(fp_pearson))) {
    stop("fp_pearson is missing required columns.")
  }
  if (!all(need_cols %in% names(fp_spearman))) {
    stop("fp_spearman is missing required columns.")
  }

  base <- fp_pearson |>
    dplyr::transmute(
      peak_ID = fp_peak,
      gene_key,
      atac_peak,
      TF = tfs,
      fp_pearson_r = r_fp,
      fp_pearson_p = p_fp,
      fp_pearson_p_adj = p_adj_fp
    )

  sp <- fp_spearman |>
    dplyr::transmute(
      peak_ID = fp_peak,
      gene_key,
      atac_peak,
      TF = tfs,
      fp_spearman_r = r_fp,
      fp_spearman_p = p_fp,
      fp_spearman_p_adj = p_adj_fp
    )

  out <- base |>
    dplyr::left_join(sp, by = c("peak_ID", "gene_key", "atac_peak", "TF"))

  if (!is.null(rna_tbl)) {
    tmp <- out |>
      dplyr::transmute(tfs = TF, gene_key) |>
      dplyr::distinct()
    tmp <- add_tf_rna_corr_cols(tmp, rna_tbl = rna_tbl, method = rna_method, cores = rna_cores)
    tmp <- tmp |>
      dplyr::transmute(
        TF = tfs,
        gene_key,
        rna_pearson_r = r_rna,
        rna_pearson_p = p_rna,
        rna_pearson_p_adj = p_adj_rna
      )
    out <- out |>
      dplyr::left_join(tmp, by = c("TF", "gene_key"))
  } else {
    out$rna_pearson_r <- NA_real_
    out$rna_pearson_p <- NA_real_
    out$rna_pearson_p_adj <- NA_real_
  }

  out
}

#' Filter TF->gene links using FP and RNA correlation criteria
#'
#' @param tbl Output from combine_fp_gene_corr_methods().
#' @param fp_p_adj_thr FDR threshold for FP correlations.
#' @param fp_r_thr Minimum FP correlation (positive).
#' @param require_pos_rna Logical; require rna_pearson_r >= 0.
#' @return List with filtered tibble and keep indices.
#' @export
filter_links_by_fp_rna_criteria <- function(
    tbl,
    fp_p_adj_thr = 0.01,
    fp_r_thr = 0.3,
    require_pos_rna = TRUE
) {
  stopifnot(is.data.frame(tbl))
  need <- c("fp_pearson_p_adj", "fp_pearson_r", "fp_spearman_p_adj", "fp_spearman_r", "rna_pearson_r")
  if (!all(need %in% names(tbl))) {
    stop("Input table is missing required columns for filtering.")
  }

  idx_keep <- union(
    which(tbl$fp_pearson_p_adj < fp_p_adj_thr & tbl$fp_pearson_r > fp_r_thr),
    which(tbl$fp_spearman_p_adj < fp_p_adj_thr & tbl$fp_pearson_r > fp_r_thr)
  )
  idx_drop <- which(
    tbl$fp_pearson_r < fp_r_thr |
      tbl$fp_spearman_r < fp_r_thr |
      (require_pos_rna & (tbl$rna_pearson_r < 0))
  )
  idx <- setdiff(idx_keep, idx_drop)

  list(
    links = tbl[idx, , drop = FALSE],
    keep_idx = idx
  )
}

#' Build TF->gene link status matrix (per condition) using FP bound and RNA flags
#'
#' @param links Tibble with TF, gene_key, peak_ID columns.
#' @param fp_bound Tibble with peak_ID + condition columns (0/1).
#' @param rna_expressed Tibble with ensembl_gene_id, HGNC + condition columns (0/1).
#' @param tf_col Column name for TF.
#' @param gene_col Column name for gene key.
#' @param peak_col Column name for peak ID.
#' @param fp_score_tbl Optional footprint score table (peak_ID + condition columns).
#' @param fp_score_threshold Numeric threshold on fp_score.
#' @param atac_score_tbl Optional ATAC score table (atac_peak + condition columns).
#' @param atac_score_threshold Numeric threshold on atac_score.
#' @param atac_peak_col Column in `links` containing atac_peak (default "atac_peak").
#' @param require_fp_bound If TRUE, require fp_bound > 0.
#' @param require_gene_expr If TRUE, require TF and gene expression flags > 0.
#' @param require_fp_score If TRUE, require fp_score >= fp_score_threshold.
#' @param require_atac_score If TRUE, require atac_score >= atac_score_threshold.
#' @param return_tbl If TRUE, return the status table in memory.
#' @param out_file Optional CSV path to write the matrix (wide).
#' @param chunk_size Row chunk size for streaming.
#' @param return_keep Logical; if TRUE return per-link keep flags (any condition).
#' @param verbose Verbosity.
#' @return List with out_file and keep flags (if requested).
#' @export
build_link_status_matrix <- function(
    links,
    fp_bound,
    rna_expressed,
    tf_col = "TF",
    gene_col = "gene_key",
    peak_col = "peak_ID",
    out_file = NULL,
    chunk_size = 50000L,
    return_keep = FALSE,
    filter_any = FALSE,
    verbose = TRUE,
    fp_score_tbl = NULL,
    fp_score_threshold = NULL,
    atac_score_tbl = NULL,
    atac_score_threshold = NULL,
    atac_peak_col = "atac_peak",
    require_fp_bound = TRUE,
    require_gene_expr = TRUE,
    require_fp_score = FALSE,
    require_atac_score = FALSE,
    return_tbl = FALSE
) {
  stopifnot(is.data.frame(links), is.data.frame(fp_bound), is.data.frame(rna_expressed))
  stopifnot(all(c(tf_col, gene_col, peak_col) %in% names(links)))
  if (!"peak_ID" %in% names(fp_bound)) stop("`fp_bound` must include peak_ID.")
  if (!all(c("ensembl_gene_id", "HGNC") %in% names(rna_expressed))) {
    stop("`rna_expressed` must include ensembl_gene_id and HGNC.")
  }
  if (isTRUE(require_fp_score) && !is.data.frame(fp_score_tbl)) {
    stop("require_fp_score=TRUE but `fp_score_tbl` is missing or invalid.")
  }
  if (isTRUE(require_atac_score) && !is.data.frame(atac_score_tbl)) {
    stop("require_atac_score=TRUE but `atac_score_tbl` is missing or invalid.")
  }
  if (isTRUE(require_fp_score) && (is.null(fp_score_threshold) || !is.finite(fp_score_threshold))) {
    stop("require_fp_score=TRUE but `fp_score_threshold` is missing or invalid.")
  }
  if (isTRUE(require_atac_score) && (is.null(atac_score_threshold) || !is.finite(atac_score_threshold))) {
    stop("require_atac_score=TRUE but `atac_score_threshold` is missing or invalid.")
  }
  if (isTRUE(require_atac_score) && !atac_peak_col %in% names(links)) {
    stop("require_atac_score=TRUE but links missing atac_peak column.")
  }

  cond_cols <- intersect(
    setdiff(names(fp_bound), "peak_ID"),
    setdiff(names(rna_expressed), c("ensembl_gene_id", "HGNC"))
  )
  if (isTRUE(require_fp_score) && is.data.frame(fp_score_tbl)) {
    cond_cols <- intersect(cond_cols, setdiff(names(fp_score_tbl), "peak_ID"))
  }
  if (isTRUE(require_atac_score) && is.data.frame(atac_score_tbl)) {
    atac_peak_name <- if ("atac_peak" %in% names(atac_score_tbl)) "atac_peak" else "peak_ID"
    cond_cols <- intersect(cond_cols, setdiff(names(atac_score_tbl), atac_peak_name))
  }
  if (!length(cond_cols)) stop("No shared condition columns between fp_bound and rna_expressed.")

  expr_hgnc <- rna_expressed[!is.na(rna_expressed$HGNC) & rna_expressed$HGNC != "",
                             c("HGNC", cond_cols), drop = FALSE]
  names(expr_hgnc)[1] <- "gene_key"
  expr_ensg <- rna_expressed[!is.na(rna_expressed$ensembl_gene_id) & rna_expressed$ensembl_gene_id != "",
                             c("ensembl_gene_id", cond_cols), drop = FALSE]
  names(expr_ensg)[1] <- "gene_key"

  idx_fp <- match(links[[peak_col]], fp_bound$peak_ID)
  idx_tf <- match(links[[tf_col]], expr_hgnc$gene_key)
  idx_gene_hgnc <- match(links[[gene_col]], expr_hgnc$gene_key)
  idx_gene_ensg <- match(links[[gene_col]], expr_ensg$gene_key)
  idx_fp_score <- NULL
  if (isTRUE(require_fp_score) && is.data.frame(fp_score_tbl)) {
    idx_fp_score <- match(links[[peak_col]], fp_score_tbl$peak_ID)
  }
  idx_atac_score <- NULL
  atac_peak_name <- NULL
  if (isTRUE(require_atac_score) && is.data.frame(atac_score_tbl)) {
    atac_peak_name <- if ("atac_peak" %in% names(atac_score_tbl)) "atac_peak" else "peak_ID"
    idx_atac_score <- match(links[[atac_peak_col]], atac_score_tbl[[atac_peak_name]])
  }

  n <- nrow(links)
  keep_any <- if (return_keep) logical(n) else NULL

  wrote_header <- FALSE
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(out_file)) file.remove(out_file)
  }
  idx_chunks <- split(seq_len(n), ceiling(seq_len(n) / as.integer(chunk_size)))
  status_chunks <- if (isTRUE(return_tbl)) vector("list", length = length(idx_chunks)) else NULL
  for (k in seq_along(idx_chunks)) {
    ii <- idx_chunks[[k]]

    fp_mat <- matrix(0L, nrow = length(ii), ncol = length(cond_cols))
    ok_fp <- !is.na(idx_fp[ii])
    if (any(ok_fp)) {
      fp_mat[ok_fp, ] <- as.matrix(fp_bound[idx_fp[ii][ok_fp], cond_cols, drop = FALSE])
    }
    fp_mat[is.na(fp_mat)] <- 0L

    tf_mat <- matrix(0L, nrow = length(ii), ncol = length(cond_cols))
    ok_tf <- !is.na(idx_tf[ii])
    if (any(ok_tf)) {
      tf_mat[ok_tf, ] <- as.matrix(expr_hgnc[idx_tf[ii][ok_tf], cond_cols, drop = FALSE])
    }
    tf_mat[is.na(tf_mat)] <- 0L

    gene_mat <- matrix(0L, nrow = length(ii), ncol = length(cond_cols))
    ok_gene_h <- !is.na(idx_gene_hgnc[ii])
    if (any(ok_gene_h)) {
      gene_mat[ok_gene_h, ] <- as.matrix(expr_hgnc[idx_gene_hgnc[ii][ok_gene_h], cond_cols, drop = FALSE])
    }
    ok_gene_e <- !ok_gene_h & !is.na(idx_gene_ensg[ii])
    if (any(ok_gene_e)) {
      gene_mat[ok_gene_e, ] <- as.matrix(expr_ensg[idx_gene_ensg[ii][ok_gene_e], cond_cols, drop = FALSE])
    }
    gene_mat[is.na(gene_mat)] <- 0L

    fp_score_mat <- NULL
    if (isTRUE(require_fp_score)) {
      fp_score_mat <- matrix(0, nrow = length(ii), ncol = length(cond_cols))
      if (!is.null(idx_fp_score)) {
        ok_fp_score <- !is.na(idx_fp_score[ii])
        if (any(ok_fp_score)) {
          fp_score_mat[ok_fp_score, ] <- as.matrix(fp_score_tbl[idx_fp_score[ii][ok_fp_score], cond_cols, drop = FALSE])
        }
      }
      fp_score_mat[is.na(fp_score_mat)] <- 0
    }

    atac_score_mat <- NULL
    if (isTRUE(require_atac_score)) {
      atac_score_mat <- matrix(0, nrow = length(ii), ncol = length(cond_cols))
      if (!is.null(idx_atac_score)) {
        ok_atac_score <- !is.na(idx_atac_score[ii])
        if (any(ok_atac_score)) {
          atac_score_mat[ok_atac_score, ] <- as.matrix(atac_score_tbl[idx_atac_score[ii][ok_atac_score], cond_cols, drop = FALSE])
        }
      }
      atac_score_mat[is.na(atac_score_mat)] <- 0
    }

    active <- matrix(TRUE, nrow = length(ii), ncol = length(cond_cols))
    if (isTRUE(require_fp_bound)) {
      active <- active & (fp_mat > 0L)
    }
    if (isTRUE(require_gene_expr)) {
      active <- active & (tf_mat > 0L) & (gene_mat > 0L)
    }
    if (isTRUE(require_fp_score)) {
      active <- active & (fp_score_mat >= fp_score_threshold)
    }
    if (isTRUE(require_atac_score)) {
      active <- active & (atac_score_mat >= atac_score_threshold)
    }
    active_int <- matrix(as.integer(active), nrow = nrow(active), ncol = ncol(active))
    colnames(active_int) <- cond_cols

    keep_chunk <- rowSums(active_int) > 0L
    if (return_keep) {
      keep_any[ii] <- keep_chunk
    }

    out_chunk <- cbind(
      links[ii, c(tf_col, gene_col, peak_col), drop = FALSE],
      as.data.frame(active_int, stringsAsFactors = FALSE)
    )
    names(out_chunk)[1:3] <- c("TF", "gene_key", "peak_ID")
    if (isTRUE(filter_any)) {
      out_chunk <- out_chunk[keep_chunk, , drop = FALSE]
    }

    if (!is.null(out_file) && nrow(out_chunk)) {
      data.table::fwrite(
        out_chunk,
        out_file,
        sep = ",",
        append = wrote_header,
        col.names = !wrote_header
      )
      wrote_header <- TRUE
    }
    if (isTRUE(return_tbl)) {
      status_chunks[[k]] <- out_chunk
    }
  }

  status_tbl <- NULL
  if (isTRUE(return_tbl)) {
    status_tbl <- dplyr::bind_rows(status_chunks)
  }

  if (isTRUE(verbose)) {
    .log_inform(sprintf(
      "Built link status matrix for %s links across %s condition(s).",
      format(n, big.mark = ","),
      length(cond_cols)
    ))
    if (!is.null(out_file)) {
      .log_inform(sprintf("Wrote link status matrix: %s", out_file))
    }
  }

  invisible(list(out_file = out_file, keep = keep_any, status_tbl = status_tbl))
}


# Basal links
#' @export
make_basal_links <- function(fp_gene_corr_kept, fp_annotation,
                             out_dir, prefix = "lighting",
                             rna_tbl = NULL,
                             rna_method = c("pearson", "spearman", "kendall"),
                             rna_cores = max(1L, parallel::detectCores(logical = TRUE) - 2L),
                             fp_variance = NULL,
                             rna_variance = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rna_method <- match.arg(rna_method)

  fp_var_tbl <- NULL
  if (is.data.frame(fp_variance) &&
      "peak_ID" %in% names(fp_variance) &&
      all(c("mean_raw", "var_raw", "rsd") %in% names(fp_variance))) {
    fp_var_tbl <- fp_variance |>
      dplyr::transmute(
        peak_ID,
        fp_mean_raw = mean_raw,
        fp_var_raw  = var_raw,
        fp_rsd      = rsd
      ) |>
      dplyr::distinct()
  }

  gene_var_tbl <- NULL
  if (is.data.frame(rna_variance) &&
      all(c("ensembl_gene_id", "HGNC", "mean_raw", "var_raw", "rsd") %in% names(rna_variance))) {
    gene_var_h <- rna_variance[
      !is.na(rna_variance$HGNC) & rna_variance$HGNC != "",
      c("HGNC", "mean_raw", "var_raw", "rsd"), drop = FALSE
    ]
    names(gene_var_h) <- c("gene_key", "gene_mean_raw", "gene_var_raw", "gene_rsd")

    gene_var_e <- rna_variance[
      !is.na(rna_variance$ensembl_gene_id) & rna_variance$ensembl_gene_id != "",
      c("ensembl_gene_id", "mean_raw", "var_raw", "rsd"), drop = FALSE
    ]
    names(gene_var_e) <- c("gene_key", "gene_mean_raw", "gene_var_raw", "gene_rsd")

    gene_var_tbl <- dplyr::bind_rows(
      gene_var_h,
      gene_var_e[!(gene_var_e$gene_key %in% gene_var_h$gene_key), , drop = FALSE]
    ) |>
      dplyr::distinct(gene_key, .keep_all = TRUE)
  }

  fp_gene_corr_use <- fp_gene_corr_kept
  if (!is.null(rna_tbl)) {
    fp_gene_corr_use <- add_tf_rna_corr_cols(
      fp_tbl = fp_gene_corr_use,
      rna_tbl = rna_tbl,
      method = rna_method,
      cores = rna_cores
    )
  } else {
    fp_gene_corr_use$r_rna <- NA_real_
    fp_gene_corr_use$p_rna <- NA_real_
    fp_gene_corr_use$p_adj_rna <- NA_real_
  }

  # from fp_gene_corr_kept -> main columns
  has_atac_peak <- "atac_peak" %in% names(fp_gene_corr_use)
  core <- fp_gene_corr_use |>
    dplyr::transmute(
      TF         = tfs,
      gene_key,
      peak_ID    = fp_peak,
      atac_peak  = if (has_atac_peak) .data$atac_peak else NA_character_,
      edge_weight = 1,                # constant 1 as requested
      r_gene      = r_fp,
      p_gene      = p_fp,
      p_adj_gene  = p_adj_fp,
      n_used_tf   = n_fp,
      r_rna_gene     = dplyr::if_else(is.finite(.data$r_rna), .data$r_rna, NA_real_),
      p_rna_gene     = dplyr::if_else(is.finite(.data$p_rna), .data$p_rna, NA_real_),
      p_rna_adj_gene = dplyr::if_else(is.finite(.data$p_adj_rna), .data$p_adj_rna, NA_real_)
    ) |>
    dplyr::distinct()

  # from fp_annotation -> r_tf / p_tf / p_adj_tf + motif
  tfcols <- fp_annotation |>
    dplyr::transmute(
      peak_ID   = fp_peak,
      TF        = tfs,
      motif     = motifs,
      r_tf      = corr_fp_tf_r,
      p_tf      = corr_fp_tf_p,
      p_adj_tf  = corr_fp_tf_p_adj
    ) |>
    dplyr::distinct()

  basal <- core |>
    dplyr::left_join(tfcols, by = c("peak_ID","TF"))

  if (!is.null(fp_var_tbl)) {
    basal <- basal |> dplyr::left_join(fp_var_tbl, by = "peak_ID")
  }
  if (!is.null(gene_var_tbl)) {
    basal <- basal |> dplyr::left_join(gene_var_tbl, by = "gene_key")
  }

  basal <- basal |>
    dplyr::select(
      TF, gene_key, peak_ID, atac_peak, edge_weight,
      r_gene, p_gene, p_adj_gene,
      n_used_tf, r_tf, p_tf, p_adj_tf, motif,
      r_rna_gene, p_rna_gene, p_rna_adj_gene,
      dplyr::any_of(c("fp_mean_raw", "fp_var_raw", "fp_rsd",
                      "gene_mean_raw", "gene_var_raw", "gene_rsd"))
    )

  readr::write_csv(
    basal, file.path(out_dir, sprintf("%s_overall_tf_gene_links.csv", prefix))
  )
  basal
}

# Lighting by condition
#' Build per-condition TF->gene link tables ("lighting") using a chosen sample label column.
#'
#' @param ds A list containing at least: fp_score (peak_ID + samples), rna (HGNC/ensembl_gene_id + samples),
#'   and sample_metadata_used (must include 'id' and optionally the chosen label column).
#' @param basal_links Basal TF->enhancer->gene link table.
#' @param out_dir Output directory.
#' @param prefix Output file prefix.
#' @param label_col Character scalar **or** character vector of candidates. The first that exists in
#'   `ds$sample_metadata_used` will be used; else falls back to `"id"`. Defaults to
#'   c("strict_match_rna","cell_stress_type") for backward-compatibility.
#' @param link_score_threshold Numeric threshold on |link_score|.
#' @param fp_score_threshold Numeric threshold on fp_score.
#' @param tf_expr_threshold Numeric threshold on TF expression.
#' @param fp_bound_tbl Optional fp_bound matrix (peak_ID + condition columns; 0/1).
#' @param rna_expressed_tbl Optional gene expression flag matrix (HGNC/ensembl + condition columns; 0/1).
#' @param atac_score_tbl Optional ATAC score table (atac_peak + condition columns).
#' @param atac_score_threshold Numeric threshold on atac_score.
#' @param require_atac_score If TRUE, require atac_score >= atac_score_threshold to mark links active.
#' @param fp_annotation_tbl Optional fp_annotation table to map peak_ID -> atac_peak when needed.
#' @param require_fp_bound If TRUE, require fp_bound > 0 to mark links active.
#' @param require_gene_expr If TRUE, require TF and gene expression flags to mark links active.
#' @param gene_expr_threshold Threshold for expression flags (default 1).
#' @param filter_active If TRUE, write only active links per condition.
#' @param fp_variance_tbl Optional footprint variance table (peak_ID + mean_raw/var_raw/rsd).
#' @param rna_variance_tbl Optional RNA variance table (ensembl_gene_id/HGNC + mean_raw/var_raw/rsd).
#' @param use_parallel Whether to parallelize across conditions.
#' @param workers Number of worker processes when `use_parallel = TRUE`.
#' @param reuse_existing If TRUE, skip conditions with existing output files.
#' @param verbose Verbose messages.
#'
#' @export
light_by_condition <- function(ds, basal_links,
                               out_dir, prefix = "lighting",
                               label_col = c("strict_match_rna","cell_stress_type"),
                               link_score_threshold = 0,
                               fp_score_threshold   = 1,
                               tf_expr_threshold    = 10,
                               fp_bound_tbl         = NULL,
                               rna_expressed_tbl    = NULL,
                               atac_score_tbl       = NULL,
                               atac_score_threshold = 0,
                               require_atac_score   = FALSE,
                               fp_annotation_tbl    = NULL,
                               require_fp_bound     = FALSE,
                               require_gene_expr    = FALSE,
                               gene_expr_threshold  = 1L,
                               filter_active        = FALSE,
                               use_parallel = TRUE,
                               workers = max(1L, round(parallel::detectCores(logical = TRUE)/2)),
                               reuse_existing = TRUE,
                               verbose = TRUE,
                               fp_variance_tbl = NULL,
                               rna_variance_tbl = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  matrices_dir <- file.path(out_dir, "per_condition_link_matrices")
  dir.create(matrices_dir, recursive = TRUE, showWarnings = FALSE)
  .log  <- function(...) if (isTRUE(verbose)) message("[light_by_condition] ", sprintf(...))
  .safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
  .cond_file <- function(lab) {
    s_lab <- .safe(lab)
    if (is.character(prefix) && nzchar(prefix)) {
      sprintf("%s_%s_tf_gene_links.csv", prefix, s_lab)
    } else {
      sprintf("%s_tf_gene_links.csv", s_lab)
    }
  }
  old_pat <- if (is.character(prefix) && nzchar(prefix)) {
    sprintf("^%s_cond-(.*)_tf_gene_links\\.csv$", prefix)
  } else {
    "^cond-(.*)_tf_gene_links\\.csv$"
  }
  old_files <- unique(c(
    list.files(out_dir, pattern = old_pat, full.names = TRUE),
    list.files(matrices_dir, pattern = old_pat, full.names = TRUE)
  ))
  if (length(old_files)) {
    for (p in old_files) {
      b <- basename(p)
      lbl <- sub(old_pat, "\\1", b)
      p_new <- file.path(matrices_dir, .cond_file(lbl))
      if (!identical(p, p_new) && !file.exists(p_new)) {
        ok <- file.rename(p, p_new)
        if (isTRUE(ok) && isTRUE(verbose)) {
          .log_inform("[light_by_condition] Renamed legacy file {b} -> {basename(p_new)}.")
        }
      }
    }
  }

  # ------- sanity & shared sample ids -------
  stopifnot(is.data.frame(ds$fp_score), "peak_ID" %in% names(ds$fp_score))
  stopifnot(is.data.frame(ds$rna), all(c("HGNC","ensembl_gene_id") %in% names(ds$rna)))

  samp_fp  <- setdiff(names(ds$fp_score), "peak_ID")
  samp_rna <- setdiff(names(ds$rna),      c("HGNC","ensembl_gene_id"))
  common_ids <- intersect(samp_fp, samp_rna)
  if (length(common_ids) == 0) stop("No shared samples between fp_score and rna.")

  # --- minimal de-dup: ensure fp_score has 1 row per peak_ID to avoid join multiplication ---
  fp_u <- ds$fp_score[, c("peak_ID", common_ids), drop = FALSE]
  if (anyDuplicated(fp_u$peak_ID)) {
    .log(
      "ds$fp_score has duplicated peak_ID (%s rows; %s unique). Using distinct(peak_ID, .keep_all=TRUE).",
      format(nrow(fp_u), big.mark = ","),
      format(dplyr::n_distinct(fp_u$peak_ID), big.mark = ",")
    )
    fp_u <- dplyr::distinct(fp_u, peak_ID, .keep_all = TRUE)
  }

  # ------- metadata (prefer provided label_col when present; otherwise fallback to id) -------
  meta_src <- ds$sample_metadata_used
  stopifnot(is.data.frame(meta_src), "id" %in% names(meta_src))

  # resolve label_col: allow any name; if a vector, take the first that exists; else use "id"
  if (length(label_col) > 1L) {
    pick <- label_col[label_col %in% names(meta_src)]
    label_col <- if (length(pick)) pick[1L] else "id"
  } else {
    if (!is.character(label_col) || length(label_col) != 1L || !(label_col %in% names(meta_src))) {
      label_col <- "id"
    }
  }

  # Match metadata rows to sample columns (id or label_col can map to columns)
  cond_id <- rep(NA_character_, nrow(meta_src))
  if (label_col %in% names(meta_src)) {
    label_vals <- meta_src[[label_col]]
    hit_label <- label_vals %in% common_ids
    cond_id[hit_label] <- label_vals[hit_label]
  }
  hit_id <- meta_src$id %in% common_ids
  cond_id[is.na(cond_id) & hit_id] <- meta_src$id[is.na(cond_id) & hit_id]

  keep <- meta_src[!is.na(cond_id), , drop = FALSE]
  if (!nrow(keep)) {
    stop(
      "No sample metadata ids or labels match fp_score/rna column names. ",
      "Check that sample_metadata_used$id or '", label_col, "' match the column names."
    )
  }

  keep$.__cond_id <- cond_id[!is.na(cond_id)]
  if (isTRUE(verbose)) {
    n_label <- if (exists("hit_label", inherits = FALSE)) sum(hit_label, na.rm = TRUE) else 0L
    n_id <- sum(hit_id, na.rm = TRUE)
    .log_inform(
      "[light_by_condition] Matched {nrow(keep)} sample{?s} (label_col: {n_label}; id: {n_id})."
    )
  }
  base_label <- if (label_col %in% names(keep)) keep[[label_col]] else keep$id
  base_label[is.na(base_label) | base_label == ""] <- keep$id

  # Disambiguate duplicate labels: append "_<id>" for duplicates
  dup_idx <- ave(seq_len(nrow(keep)), base_label, FUN = seq_along)
  label <- ifelse(dup_idx > 1L, paste0(base_label, "_", keep$id), base_label)
  meta <- data.frame(id = keep$.__cond_id, label = label, stringsAsFactors = FALSE)

  have_fp_bound <- is.data.frame(fp_bound_tbl) && "peak_ID" %in% names(fp_bound_tbl)
  have_rna_flags <- is.data.frame(rna_expressed_tbl) &&
    all(c("ensembl_gene_id", "HGNC") %in% names(rna_expressed_tbl))
  if (is.null(atac_score_tbl) && is.list(ds) && is.data.frame(ds$atac_score)) {
    atac_score_tbl <- ds$atac_score
  }
  have_atac_score <- is.data.frame(atac_score_tbl)
  if (isTRUE(require_fp_bound) && !have_fp_bound) {
    stop("require_fp_bound=TRUE but fp_bound_tbl is missing or invalid.")
  }
  if (isTRUE(require_gene_expr) && !have_rna_flags) {
    stop("require_gene_expr=TRUE but rna_expressed_tbl is missing or invalid.")
  }
  if (isTRUE(require_atac_score) && !have_atac_score) {
    stop("require_atac_score=TRUE but atac_score_tbl is missing or invalid.")
  }

  if (!"atac_peak" %in% names(basal_links) || all(is.na(basal_links$atac_peak))) {
    if (is.null(fp_annotation_tbl) && is.list(ds) && is.data.frame(ds$fp_annotation)) {
      fp_annotation_tbl <- ds$fp_annotation
    }
    if (is.data.frame(fp_annotation_tbl) &&
        all(c("fp_peak", "atac_peak", "tfs") %in% names(fp_annotation_tbl))) {
      map_atac <- fp_annotation_tbl |>
        dplyr::transmute(peak_ID = fp_peak, TF = tfs, atac_peak) |>
        dplyr::distinct()
      basal_links <- basal_links |> dplyr::left_join(map_atac, by = c("peak_ID", "TF"))
    }
  }

  fp_bound_u <- fp_bound_tbl
  if (isTRUE(have_fp_bound)) {
    if (anyDuplicated(fp_bound_u$peak_ID)) {
      .log(
        "fp_bound_tbl has duplicated peak_ID (%s rows; %s unique). Using distinct(peak_ID, .keep_all=TRUE).",
        format(nrow(fp_bound_u), big.mark = ","),
        format(dplyr::n_distinct(fp_bound_u$peak_ID), big.mark = ",")
      )
      fp_bound_u <- dplyr::distinct(fp_bound_u, peak_ID, .keep_all = TRUE)
    }
  }

  if (is.null(fp_variance_tbl) && is.list(ds) && is.data.frame(ds$fp_variance)) {
    fp_variance_tbl <- ds$fp_variance
  }
  if (is.null(rna_variance_tbl) && is.list(ds) && is.data.frame(ds$rna_variance)) {
    rna_variance_tbl <- ds$rna_variance
  }

  pick_cond_col <- function(tbl, cond_id, cond_label) {
    if (cond_id %in% names(tbl)) return(cond_id)
    if (cond_label %in% names(tbl)) return(cond_label)
    alt <- sub(paste0("_", cond_id, "$"), "", cond_label)
    if (alt %in% names(tbl)) return(alt)
    NA_character_
  }

  # ------- RNA table keyed by gene_key (prefer HGNC, else ENSG) -------
  sample_cols_rna <- samp_rna
  rna_hgnc <- ds$rna[!is.na(ds$rna$HGNC) & ds$rna$HGNC != "", c("HGNC", sample_cols_rna), drop = FALSE]
  names(rna_hgnc)[1] <- "gene_key"
  if (anyDuplicated(rna_hgnc$gene_key)) {
    .log(
      "rna has duplicated HGNC (%s rows; %s unique). Using distinct(HGNC, .keep_all=TRUE).",
      format(nrow(rna_hgnc), big.mark = ","),
      format(dplyr::n_distinct(rna_hgnc$gene_key), big.mark = ",")
    )
    rna_hgnc <- dplyr::distinct(rna_hgnc, gene_key, .keep_all = TRUE)
  }
  rna_ensg <- ds$rna[!is.na(ds$rna$ensembl_gene_id) & ds$rna$ensembl_gene_id != "", c("ensembl_gene_id", sample_cols_rna), drop = FALSE]
  names(rna_ensg)[1] <- "gene_key"
  if (anyDuplicated(rna_ensg$gene_key)) {
    .log(
      "rna has duplicated ensembl_gene_id (%s rows; %s unique). Using distinct(ensembl_gene_id, .keep_all=TRUE).",
      format(nrow(rna_ensg), big.mark = ","),
      format(dplyr::n_distinct(rna_ensg$gene_key), big.mark = ",")
    )
    rna_ensg <- dplyr::distinct(rna_ensg, gene_key, .keep_all = TRUE)
  }
  # bind_rows preferring HGNC keys; anti-join equivalent:
  rna_hgnc_keys <- rna_hgnc$gene_key
  rna_ensg_only <- rna_ensg[!(rna_ensg$gene_key %in% rna_hgnc_keys), , drop = FALSE]
  rna_gk <- dplyr::bind_rows(rna_hgnc, rna_ensg_only)
  rna_gk <- dplyr::distinct(rna_gk, gene_key, .keep_all = TRUE)

  fp_var_tbl <- NULL
  if (is.data.frame(fp_variance_tbl) &&
      "peak_ID" %in% names(fp_variance_tbl) &&
      all(c("mean_raw", "var_raw", "rsd") %in% names(fp_variance_tbl))) {
    fp_var_tbl <- fp_variance_tbl |>
      dplyr::transmute(
        peak_ID,
        fp_mean_raw = mean_raw,
        fp_var_raw  = var_raw,
        fp_rsd      = rsd
      ) |>
      dplyr::distinct()
  }

  gene_var_tbl <- NULL
  if (is.data.frame(rna_variance_tbl) &&
      all(c("ensembl_gene_id", "HGNC", "mean_raw", "var_raw", "rsd") %in% names(rna_variance_tbl))) {
    gene_var_h <- rna_variance_tbl[
      !is.na(rna_variance_tbl$HGNC) & rna_variance_tbl$HGNC != "",
      c("HGNC", "mean_raw", "var_raw", "rsd"), drop = FALSE
    ]
    names(gene_var_h) <- c("gene_key", "gene_mean_raw", "gene_var_raw", "gene_rsd")

    gene_var_e <- rna_variance_tbl[
      !is.na(rna_variance_tbl$ensembl_gene_id) & rna_variance_tbl$ensembl_gene_id != "",
      c("ensembl_gene_id", "mean_raw", "var_raw", "rsd"), drop = FALSE
    ]
    names(gene_var_e) <- c("gene_key", "gene_mean_raw", "gene_var_raw", "gene_rsd")

    gene_var_tbl <- dplyr::bind_rows(
      gene_var_h,
      gene_var_e[!(gene_var_e$gene_key %in% gene_var_h$gene_key), , drop = FALSE]
    ) |>
      dplyr::distinct(gene_key, .keep_all = TRUE)
  }

  # ------- helper: normalize edge_weight within (peak_ID, gene_key) -------
  normalize_edge_wt <- function(df, floor = 0) {
    df |>
      dplyr::mutate(edge_weight_raw = pmax(tf_expr, 0)) |>
      dplyr::group_by(peak_ID, gene_key) |>
      dplyr::mutate(
        .tot = sum(edge_weight_raw, na.rm = TRUE),
        edge_weight = dplyr::if_else(.tot > 0, edge_weight_raw / .tot, floor),
        .tot = NULL
      ) |>
      dplyr::ungroup()
  }

  # ------- one condition -------
  build_one <- function(cond_id, cond_label) {
    out_file <- file.path(matrices_dir, .cond_file(cond_label))
    if (isTRUE(reuse_existing) && file.exists(out_file)) {
      n_existing <- tryCatch(
        nrow(readr::read_csv(out_file, show_col_types = FALSE)),
        error = function(e) NA_integer_
      )
      if (isTRUE(verbose)) {
        .log_inform("[light_by_condition] Reusing existing file: {basename(out_file)}")
      }
      return(tibble::tibble(label = cond_label, n_links_rows = n_existing))
    }

    # dynamic pulls via base subsetting to avoid NSE/rlang
    fp_score_one <- fp_u[, c("peak_ID", cond_id), drop = FALSE]
    names(fp_score_one) <- c("peak_ID", "fp_score")

    tf_expr_tbl <- rna_hgnc[, c("gene_key", cond_id), drop = FALSE]
    names(tf_expr_tbl) <- c("HGNC", "tf_expr")

    gene_expr_tbl <- rna_gk[, c("gene_key", cond_id), drop = FALSE]
    names(gene_expr_tbl) <- c("gene_key", "gene_expr")

    fp_bound_one <- NULL
    if (have_fp_bound) {
      col_bound <- pick_cond_col(fp_bound_u, cond_id, cond_label)
      if (!is.na(col_bound)) {
        fp_bound_one <- fp_bound_u[, c("peak_ID", col_bound), drop = FALSE]
        names(fp_bound_one) <- c("peak_ID", "fp_bound")
      }
    }

    tf_flag_tbl <- NULL
    gene_flag_tbl <- NULL
    if (have_rna_flags) {
      col_flag <- pick_cond_col(rna_expressed_tbl, cond_id, cond_label)
      if (!is.na(col_flag)) {
        tf_flag_tbl <- rna_expressed_tbl[, c("HGNC", col_flag), drop = FALSE]
        names(tf_flag_tbl) <- c("HGNC", "tf_expr_flag")
        if (anyDuplicated(tf_flag_tbl$HGNC)) {
          tf_flag_tbl <- dplyr::distinct(tf_flag_tbl, HGNC, .keep_all = TRUE)
        }

        gene_flag_h <- rna_expressed_tbl[
          !is.na(rna_expressed_tbl$HGNC) & rna_expressed_tbl$HGNC != "",
          c("HGNC", col_flag), drop = FALSE
        ]
        names(gene_flag_h) <- c("gene_key", "gene_expr_flag")

        gene_flag_e <- rna_expressed_tbl[
          !is.na(rna_expressed_tbl$ensembl_gene_id) & rna_expressed_tbl$ensembl_gene_id != "",
          c("ensembl_gene_id", col_flag), drop = FALSE
        ]
        names(gene_flag_e) <- c("gene_key", "gene_expr_flag")

        gene_flag_tbl <- dplyr::bind_rows(
          gene_flag_h,
          gene_flag_e[!(gene_flag_e$gene_key %in% gene_flag_h$gene_key), , drop = FALSE]
        ) |>
          dplyr::distinct(gene_key, .keep_all = TRUE)
      }
    }

    atac_score_one <- NULL
    if (have_atac_score) {
      atac_peak_col <- if ("atac_peak" %in% names(atac_score_tbl)) "atac_peak" else "peak_ID"
      col_atac <- pick_cond_col(atac_score_tbl, cond_id, cond_label)
      if (!is.na(col_atac)) {
        atac_score_one <- atac_score_tbl[, c(atac_peak_col, col_atac), drop = FALSE]
        names(atac_score_one) <- c("atac_peak", "atac_score")
      }
    }

    links <- basal_links |>
      dplyr::left_join(fp_score_one,  by = "peak_ID") |>
      dplyr::left_join(tf_expr_tbl,   by = c("TF" = "HGNC")) |>
      dplyr::left_join(gene_expr_tbl, by = "gene_key")

    if (!is.null(fp_bound_one)) {
      links <- links |> dplyr::left_join(fp_bound_one, by = "peak_ID")
    } else {
      links <- links |> dplyr::mutate(fp_bound = NA_integer_)
    }

    if (!is.null(atac_score_one)) {
      links <- links |> dplyr::left_join(atac_score_one, by = "atac_peak")
    } else {
      links <- links |> dplyr::mutate(atac_score = NA_real_)
    }

    if (!is.null(tf_flag_tbl)) {
      links <- links |> dplyr::left_join(tf_flag_tbl, by = c("TF" = "HGNC"))
    } else {
      links <- links |> dplyr::mutate(tf_expr_flag = NA_integer_)
    }

    if (!is.null(gene_flag_tbl)) {
      links <- links |> dplyr::left_join(gene_flag_tbl, by = "gene_key")
    } else {
      links <- links |> dplyr::mutate(gene_expr_flag = NA_integer_)
    }

    if (!is.null(fp_var_tbl)) {
      miss_fp <- setdiff(c("fp_mean_raw", "fp_var_raw", "fp_rsd"), names(links))
      if (length(miss_fp)) {
        links <- links |> dplyr::left_join(
          fp_var_tbl[, c("peak_ID", miss_fp), drop = FALSE],
          by = "peak_ID"
        )
      }
    }
    if (!is.null(gene_var_tbl)) {
      miss_gene <- setdiff(c("gene_mean_raw", "gene_var_raw", "gene_rsd"), names(links))
      if (length(miss_gene)) {
        links <- links |> dplyr::left_join(
          gene_var_tbl[, c("gene_key", miss_gene), drop = FALSE],
          by = "gene_key"
        )
      }
    }

    links <- links |>
      dplyr::mutate(
        fp_score  = tidyr::replace_na(fp_score,  0),
        tf_expr   = tidyr::replace_na(tf_expr,   0),
        gene_expr = tidyr::replace_na(gene_expr, 0),
        fp_bound = tidyr::replace_na(fp_bound, 0L),
        tf_expr_flag = tidyr::replace_na(tf_expr_flag, 0L),
        gene_expr_flag = tidyr::replace_na(gene_expr_flag, 0L),
        atac_score = tidyr::replace_na(atac_score, 0)
      ) |>
      normalize_edge_wt(floor = 0) |>
      dplyr::mutate(
        link_score  = dplyr::coalesce(r_gene, 0) * dplyr::coalesce(fp_score, 0),
        active_link = is.finite(link_score) &
          fp_score >= fp_score_threshold &
          (if (isTRUE(require_atac_score)) atac_score >= atac_score_threshold else TRUE) &
          abs(link_score) >= link_score_threshold &
          tf_expr >= tf_expr_threshold &
          (if (isTRUE(require_fp_bound)) fp_bound > 0L else TRUE) &
          (if (isTRUE(require_gene_expr)) (tf_expr_flag >= gene_expr_threshold & gene_expr_flag >= gene_expr_threshold) else TRUE)
      ) |>
      dplyr::select(
        TF, gene_key, peak_ID, atac_peak,
        edge_weight,
        r_gene, p_gene, p_adj_gene,
        n_used_tf, r_tf, p_tf, p_adj_tf, motif,
        dplyr::any_of(c("r_rna_gene", "p_rna_gene", "p_rna_adj_gene")),
        dplyr::any_of(c("fp_mean_raw", "fp_var_raw", "fp_rsd",
                        "gene_mean_raw", "gene_var_raw", "gene_rsd")),
        fp_score, atac_score, tf_expr, gene_expr,
        fp_bound, tf_expr_flag, gene_expr_flag,
        active_link, link_score
      ) |>
      dplyr::arrange(dplyr::desc(abs(link_score)), TF, gene_key, peak_ID)

    if (isTRUE(filter_active)) {
      links <- links |>
        dplyr::filter(.data$active_link)
    }

    readr::write_csv(links, out_file)
    tibble::tibble(label = cond_label, n_links_rows = nrow(links))
  }

  # ------- run all conditions -------
  ids  <- meta$id
  labs <- meta$label

  if (use_parallel && length(ids) > 1L && workers > 1L) {
    oplan <- future::plan(); on.exit(future::plan(oplan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    idx_list <- future.apply::future_mapply(build_one, ids, labs, SIMPLIFY = FALSE)
    index <- dplyr::bind_rows(idx_list)
  } else {
    index <- purrr::map2_dfr(ids, labs, build_one)
  }

  idx_file <- if (is.character(prefix) && nzchar(prefix)) {
    sprintf("%s_per_condition_index.csv", prefix)
  } else {
    "per_condition_index.csv"
  }
  readr::write_csv(index, file.path(matrices_dir, idx_file))
  .log("Wrote per-condition index: %s rows.", format(nrow(index), big.mark = ","))
  invisible(index)
}
