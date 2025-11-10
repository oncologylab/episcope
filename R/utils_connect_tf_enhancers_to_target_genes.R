# Build a reference gene table for hg38 (GRCh38) or mm10 (GRCm38)
# Returns tibble: ensembl_gene_id, HGNC, chrom, start, end, strand, tss
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

  # seqnames → UCSC-style with 'chr' prefix if missing
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

#' Build GH-like peak→gene links by distance window (no GeneHancer needed)
#'
#' Create GeneHancer-style rows by linking peaks to nearby genes within a
#' ±\code{flank_bp} window. Works with either a ready-made gene table or a
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

  # ---- Peaks → (chrom,start,end)
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

  # ---- Genes input normalization → gn0: (chrom,start,end,gene_key) for chosen mode
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

  # ---- Extend peaks by ±flank and intersect with genes (report original peak span)
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
#' @param gh_tbl    Tibble with *at least* columns: chrom, start, end, connected_gene (names configurable via gh_cols).
#' @param gh_cols   Named list mapping column names: list(chrom="chrom", start="start", end="end", gene="connected_gene").
#'                  Optional extras (if present) are kept during overlap: id="gh_id", conf="confidence".
#' @param gene_mode "both", "ensembl", or "hgnc" — which gene identifiers to match to RNA.
#' @param fdr       BH FDR cutoff.
#' @param keep_pos  If TRUE, retain only positive correlations.
#' @param workers   Parallel workers (future.apply). Use 1 for serial.
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
    r_abs_min = 0.3,   # <- NEW: absolute correlation cutoff
    workers   = max(1L, parallel::detectCores() - 1L)
) {
  gene_mode <- match.arg(gene_mode)

  stopifnot(is.data.frame(grn_set$fp_annotation),
            is.data.frame(grn_set$atac_score),
            is.data.frame(grn_set$rna),
            is.data.frame(grn_set$sample_metadata_used))

  # samples
  ids <- grn_set$sample_metadata_used$id
  ids <- ids[ids %in% intersect(names(grn_set$atac_score), names(grn_set$rna))]
  if (!length(ids)) cli::cli_abort("No overlapping sample IDs between ATAC and RNA.")

  # ATAC peaks from fp_annotation
  ap <- unique(grn_set$fp_annotation$atac_peak)
  ap <- ap[!is.na(ap) & nzchar(ap)]
  if (!length(ap)) cli::cli_abort("No usable atac_peak values in fp_annotation.")

  atac_bed <- .parse_peak_ids(ap) |>
    dplyr::select(atac_chr, atac_start, atac_end, atac_peak)

  # --- GH normalization (generic) --------------------------------------------
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

  # --- overlap (any) ----------------------------------------------------------
  ov_raw <- valr::bed_intersect(
    atac_bed |> dplyr::transmute(chrom = atac_chr, start = atac_start, end = atac_end, atac_peak),
    gh_bed   |> dplyr::transmute(chrom = gh_chr,  start = gh_start,  end = gh_end, gh_id, gh_conf, gene_key)
  )

  nm <- names(ov_raw)
  atk_col <- if ("atac_peak.x" %in% nm) "atac_peak.x" else "atac_peak"
  ghid_col <- if ("gh_id.y" %in% nm) "gh_id.y" else "gh_id"
  ghy_col  <- if ("gh_conf.y" %in% nm) "gh_conf.y" else "gh_conf"
  ggn_col  <- if ("gene_key.y" %in% nm) "gene_key.y" else "gene_key"

  pairs0 <- ov_raw |>
    dplyr::transmute(
      atac_peak = .data[[atk_col]],
      gene_key  = .data[[ggn_col]],
      gh_id     = .data[[ghid_col]],
      gh_conf   = .data[[ghy_col]]
    ) |>
    dplyr::distinct()

  if (!nrow(pairs0)) cli::cli_abort("No ATAC↔enhancer overlaps.")

  # --- keep genes present in RNA under the requested mode ---------------------
  rna_keys <- .pick_rna_gene_keys(grn_set$rna, gene_mode)
  pairs1 <- dplyr::semi_join(pairs0, tibble::tibble(gene_key = rna_keys), by = "gene_key")
  if (!nrow(pairs1)) cli::cli_abort("After gene matching to RNA, no pairs remain.")

  # only ATAC peaks present in atac_score
  pairs1 <- dplyr::semi_join(pairs1, dplyr::distinct(grn_set$atac_score, atac_peak), by = "atac_peak")

  # --- long views -------------------------------------------------------------
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

  # --- correlation per (atac_peak, gene_key) ---------------------------------
  idx_split <- split(seq_len(nrow(pairs2)),
                     ceiling(seq_along(seq_len(nrow(pairs2))) /
                               max(1L, ceiling(nrow(pairs2) / max(1L, workers)))))

  run_chunk <- function(ii) {
    pp <- pairs2[ii, , drop = FALSE]
    # many-to-many is expected (each peak × multi-sample matrix)
    ax <- dplyr::inner_join(pp, atac_long, by = "atac_peak",
                            relationship = "many-to-many")
    if (!nrow(ax)) return(tibble::tibble(atac_peak=character(), gene_key=character(), n=integer(), r=double(), p=double()))
    # many-to-many is also expected here (gene × samples)
    axx <- dplyr::inner_join(ax, rna_long, by = c("gene_key","sample"),
                             relationship = "many-to-many")
    if (!nrow(axx)) return(tibble::tibble(atac_peak=character(), gene_key=character(), n=integer(), r=double(), p=double()))

    axx |>
      dplyr::group_by(atac_peak, gene_key) |>
      dplyr::summarise(
        n = sum(is.finite(acc) & is.finite(expr)),
        r = suppressWarnings(stats::cor(acc, expr, use = "complete.obs", method = "pearson")),
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

  if (workers > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      cli::cli_abort("Install future.apply or set workers=1.")
    }
    old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    corr_list <- future.apply::future_lapply(idx_split, run_chunk, future.seed = TRUE)
  } else {
    corr_list <- lapply(idx_split, run_chunk)
  }

  atac_gene_corr_full <- dplyr::bind_rows(corr_list) |>
    dplyr::filter(!is.na(p)) |>
    dplyr::mutate(p_adj = stats::p.adjust(p, method = "BH")) |>
    dplyr::transmute(
      atac_peak,
      gene_key,
      n_atac    = as.integer(n),
      r_atac    = as.numeric(r),
      p_atac    = as.numeric(p),
      p_adj_atac= as.numeric(p_adj)
    )

  atac_gene_corr_kept <- atac_gene_corr_full |>
    dplyr::filter(p_adj_atac < fdr, abs(r_atac) >= r_abs_min, if (keep_pos) r_atac > 0 else TRUE) |> # <- NEW: apply |r| cutoff
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

# Correlate fp_peak scores to target-gene RNA for ATAC→gene pairs that passed your filter
#' @export
correlate_fp_to_genes <- function(
    grn_set,
    atac_gene_corr_kept,
    fdr       = 0.05,
    r_abs_min = 0.3,
    keep_pos  = FALSE,
    method    = "pearson",
    workers = max(1L, round(parallel::detectCores(logical = TRUE)/2))
) {
  method <- match.arg(method)

  # ---- sanity checks ----
  req_slots <- c("fp_score","rna","fp_annotation")
  stopifnot(is.list(grn_set), all(req_slots %in% names(grn_set)))
  stopifnot(is.data.frame(atac_gene_corr_kept))
  stopifnot(all(c("atac_peak","gene_key") %in% names(atac_gene_corr_kept)))
  stopifnot("peak_ID" %in% names(grn_set$fp_score))
  stopifnot(all(c("ensembl_gene_id","HGNC") %in% names(grn_set$rna)))

  fp  <- grn_set$fp_score
  rna <- grn_set$rna
  ann <- grn_set$fp_annotation

  # samples present in BOTH fp_score and rna
  sample_cols <- intersect(
    setdiff(names(fp), "peak_ID"),
    setdiff(names(rna), c("ensembl_gene_id","HGNC"))
  )
  if (length(sample_cols) < 3L)
    stop("Not enough shared samples between fp_score and rna.")

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

  genes_keep <- unique(atac_gene_corr_kept$gene_key)
  rna_long   <- rna_long |> dplyr::semi_join(tibble::tibble(gene_key = genes_keep), by = "gene_key")

  # ---- map fp_peaks to those kept atac_peaks & genes (silence many-to-many warning) ----
  pairs_tbl <- ann |>
    dplyr::select(fp_peak, atac_peak, tfs, motifs) |>
    dplyr::inner_join(
      atac_gene_corr_kept |> dplyr::select(atac_peak, gene_key) |> dplyr::distinct(),
      by = "atac_peak",
      relationship = "many-to-many"  # <- explicitly allow many-to-many
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
    as.data.frame()                      # <- make it a data.frame
  rownames(fp_m) <- fp_m$peak_ID         # <- now rownames actually stick
  fp_m <- as.matrix(fp_m[, sample_cols, drop = FALSE])

  # RNA
  rna_m <- rna_long |>
    dplyr::semi_join(tibble::tibble(gene_key = pairs_tbl$gene_key), by = "gene_key") |>
    dplyr::distinct(gene_key, .keep_all = TRUE) |>
    as.data.frame()                      # <- make it a data.frame
  rownames(rna_m) <- rna_m$gene_key      # <- set rownames on data.frame
  rna_m <- as.matrix(rna_m[, sample_cols, drop = FALSE])

  # keep only computable pairs
  pairs_tbl <- pairs_tbl |>
    dplyr::filter(fp_peak %in% rownames(fp_m), gene_key %in% rownames(rna_m))


  if (nrow(pairs_tbl) == 0L) {
    empty <- tibble::tibble()
    return(list(fp_gene_corr_kept = empty, fp_gene_corr_full = empty))
  }

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
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      r  <- unname(ct$estimate)
      p  <- ct$p.value
    }
    list(n = as.integer(n), r = as.numeric(r), p = as.numeric(p))
  }


  # parallel if requested
  if (workers > 1L) {
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    out <- furrr::future_pmap(pairs_tbl[, c("fp_peak","gene_key")], cor_one, .progress = FALSE)
  } else {
    out <- purrr::pmap(pairs_tbl[, c("fp_peak","gene_key")], cor_one)
  }

  pairs_tbl$n_fp <- vapply(out, `[[`, integer(1), "n")
  pairs_tbl$r_fp <- vapply(out, `[[`, numeric(1), "r")
  pairs_tbl$p_fp <- vapply(out, `[[`, numeric(1), "p")
  pairs_tbl <- pairs_tbl |> dplyr::filter(!is.na(p_fp))
  pairs_tbl$p_adj_fp <- stats::p.adjust(pairs_tbl$p_fp, method = "BH")

  fp_gene_corr_full <- pairs_tbl |>
    dplyr::transmute(
      fp_peak, gene_key, atac_peak, tfs, motifs,
      n_fp = as.integer(n_fp),
      r_fp = as.numeric(r_fp),
      p_fp = as.numeric(p_fp),
      p_adj_fp = as.numeric(p_adj_fp)
    )

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



# Basal links
#' @export
make_basal_links <- function(fp_gene_corr_kept, fp_annotation,
                             out_dir, prefix = "lighting") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # from fp_gene_corr_kept -> main columns
  core <- fp_gene_corr_kept |>
    dplyr::transmute(
      TF         = tfs,
      gene_key,
      peak_ID    = fp_peak,
      edge_weight = 1,                # constant 1 as requested
      r_gene      = r_fp,
      p_gene      = p_fp,
      p_adj_gene  = p_adj_fp,
      n_used_tf   = n_fp
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
    dplyr::left_join(tfcols, by = c("peak_ID","TF")) |>
    dplyr::select(TF, gene_key, peak_ID, edge_weight,
                  r_gene, p_gene, p_adj_gene,
                  n_used_tf, r_tf, p_tf, p_adj_tf, motif)

  readr::write_csv(
    basal, file.path(out_dir, sprintf("%s_overall_tf_gene_links.csv", prefix))
  )
  basal
}

# Lighting by condition
#' @export
light_by_condition <- function(ds, basal_links,
                               out_dir, prefix = "lighting",
                               label_col = c("strict_match_rna","cell_stress_type"),
                               link_score_threshold = 0,
                               fp_score_threshold   = 1,
                               tf_expr_threshold    = 10,
                               use_parallel = TRUE,
                               workers = max(1L, round(parallel::detectCores(logical = TRUE)/2)),
                               verbose = TRUE) {
  label_col <- match.arg(label_col, c("strict_match_rna","cell_stress_type"))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .log  <- function(...) if (isTRUE(verbose)) message("[light_by_condition] ", sprintf(...))
  .safe <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  # ------- sanity & shared sample ids -------
  stopifnot(is.data.frame(ds$fp_score), "peak_ID" %in% names(ds$fp_score))
  stopifnot(is.data.frame(ds$rna), all(c("HGNC","ensembl_gene_id") %in% names(ds$rna)))

  samp_fp  <- setdiff(names(ds$fp_score), "peak_ID")
  samp_rna <- setdiff(names(ds$rna),      c("HGNC","ensembl_gene_id"))
  common_ids <- intersect(samp_fp, samp_rna)
  if (length(common_ids) == 0) stop("No shared samples between fp_score and rna.")

  # ------- metadata (use sample_metadata_used; no mapping needed) -------
  meta_src <- ds$sample_metadata_used
  stopifnot(is.data.frame(meta_src), "id" %in% names(meta_src))
  have_label <- label_col %in% names(meta_src)

  meta <- meta_src |>
    dplyr::filter(.data$id %in% common_ids) |>
    dplyr::mutate(
      .base = if (have_label) .data[[label_col]] else .data$id,
      .base = dplyr::coalesce(.base, .data$id)
    ) |>
    dplyr::group_by(.base) |>
    dplyr::mutate(.dup = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(label = ifelse(.dup > 1, paste0(.base, "_", id), .base)) |>
    dplyr::select(id, label)

  # ------- RNA table keyed by gene_key (prefer HGNC, else ENSG) -------
  sample_cols_rna <- samp_rna
  rna_hgnc <- ds$rna |>
    dplyr::filter(!is.na(HGNC), HGNC != "") |>
    dplyr::transmute(gene_key = HGNC, dplyr::across(dplyr::all_of(sample_cols_rna)))
  rna_ensg <- ds$rna |>
    dplyr::filter(!is.na(ensembl_gene_id), ensembl_gene_id != "") |>
    dplyr::transmute(gene_key = ensembl_gene_id, dplyr::across(dplyr::all_of(sample_cols_rna)))
  rna_gk <- dplyr::bind_rows(rna_hgnc, dplyr::anti_join(rna_ensg, rna_hgnc, by = "gene_key")) |>
    dplyr::distinct(gene_key, .keep_all = TRUE)

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
    atac_one <- ds$fp_score |>
      dplyr::select(peak_ID, fp_score = dplyr::all_of(cond_id)) |>
      dplyr::distinct(peak_ID, .keep_all = TRUE)

    tf_expr_tbl <- ds$rna |>
      dplyr::select(HGNC, tf_expr = dplyr::all_of(cond_id)) |>
      dplyr::distinct(HGNC, .keep_all = TRUE)

    gene_expr_tbl <- rna_gk |>
      dplyr::select(gene_key, gene_expr = dplyr::all_of(cond_id))

    links <- basal_links |>
      dplyr::left_join(atac_one,      by = "peak_ID") |>
      dplyr::left_join(tf_expr_tbl,   by = c("TF" = "HGNC")) |>
      dplyr::left_join(gene_expr_tbl, by = "gene_key") |>
      dplyr::mutate(
        fp_score  = tidyr::replace_na(fp_score,  0),
        tf_expr   = tidyr::replace_na(tf_expr,   0),
        gene_expr = tidyr::replace_na(gene_expr, 0)
      ) |>
      normalize_edge_wt(floor = 0) |>
      dplyr::mutate(
        link_score  = dplyr::coalesce(r_gene, 0) *
          dplyr::coalesce(fp_score,   0),
        active_link = is.finite(link_score) &
          fp_score >= fp_score_threshold &
          abs(link_score) >= link_score_threshold &
          tf_expr >= tf_expr_threshold
      ) |>
      dplyr::select(
        TF, gene_key, peak_ID,
        edge_weight,
        r_gene, p_gene, p_adj_gene,
        n_used_tf, r_tf, p_tf, p_adj_tf, motif,
        fp_score, tf_expr, gene_expr,
        active_link, link_score
      ) |>
      dplyr::arrange(dplyr::desc(abs(link_score)), TF, gene_key, peak_ID)

    # links |> dplyr::filter(link_score != 0)
    # basal_links |>
    #   dplyr::left_join(atac_one,      by = "peak_ID") |>
    #   dplyr::filter(!is.na(fp_score))

    readr::write_csv(
      links,
      file.path(out_dir, sprintf("%s_cond-%s_tf_gene_links.csv", prefix, .safe(cond_label)))
    )
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

  readr::write_csv(index, file.path(out_dir, sprintf("%s_per_condition_index.csv", prefix)))
  .log("Wrote per-condition index: %s rows.", format(nrow(index), big.mark = ","))
  invisible(index)
}
