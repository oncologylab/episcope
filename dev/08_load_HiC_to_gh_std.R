## ==============================================================
## 1) Load all Juicer HiCCUPS merged_loops.bedpe files
## ==============================================================

load_hiccups_loops <- function(base_dir) {
  # 1) Unzip all *.zip into subdirs named after the zip (without .zip)
  zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
  if (length(zip_files)) {
    message("Found ", length(zip_files), " zip files. Unzipping (if needed)...")
  }

  for (zf in zip_files) {
    subdir_name <- sub("\\.zip$", "", basename(zf))
    out_dir <- file.path(base_dir, subdir_name)

    if (!dir.exists(out_dir) || !length(list.files(out_dir))) {
      message("  Unzipping: ", basename(zf), " -> ", subdir_name)
      utils::unzip(zf, exdir = out_dir)
    } else {
      message("  Skipping unzip for ", basename(zf),
              " (directory already populated)")
    }
  }

  # 2) Collect all hiccups directories (top-level dirs ending with "hiccups")
  hiccups_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  hiccups_dirs <- hiccups_dirs[grepl("hiccups$", basename(hiccups_dirs))]
  if (!length(hiccups_dirs)) {
    stop("No directories ending with 'hiccups' found under: ", base_dir)
  }

  message("Found ", length(hiccups_dirs), " hiccups directories:")
  print(basename(hiccups_dirs))

  # 3) Quick peek: list files in the first few directories
  message("\n=== File listing preview (first 3 directories) ===")
  for (d in head(hiccups_dirs, 3)) {
    message("Directory: ", basename(d))
    print(list.files(d))
    message("----")
  }

  # 4) Read any files with "loops" in the name and combine
  all_loops <- list()
  idx <- 1L

  for (d in hiccups_dirs) {
    sample_id <- basename(d)
    files_in_dir <- list.files(d, full.names = TRUE)

    # Use merged_loops.bedpe (or anything with "loops" in the name)
    loop_files <- files_in_dir[
      grepl("loops", basename(files_in_dir), ignore.case = TRUE)
    ]

    if (!length(loop_files)) {
      message("No 'loops' files found in: ", sample_id, " (skipping).")
      next
    }

    for (lf in loop_files) {
      message("Reading loops file: ", basename(lf),
              " (sample: ", sample_id, ")")

      # Force all columns to character to avoid bind_rows type conflicts
      tbl <- tryCatch(
        {
          readr::read_tsv(
            lf,
            col_types = readr::cols(.default = readr::col_character()),
            progress = FALSE,
            show_col_types = FALSE
          ) |>
            tibble::as_tibble()
        },
        error = function(e) {
          message("  Failed to read with readr::read_tsv: ", e$message,
                  " — trying read.delim as character...")
          df <- utils::read.delim(
            lf,
            header = TRUE,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          df[] <- lapply(df, as.character)
          tibble::as_tibble(df)
        }
      )

      tbl$sample_id <- sample_id
      tbl$file_name <- basename(lf)

      all_loops[[idx]] <- tbl
      idx <- idx + 1L
    }
  }

  if (!length(all_loops)) {
    warning("No loops files were successfully read.")
    return(tibble::tibble())
  }

  loops_df <- dplyr::bind_rows(all_loops)

  message("\n=== Preview of combined raw loops tibble (first 20 rows) ===")
  print(utils::head(loops_df, 20))

  loops_df
}

## ==============================================================
## 2) Clean HiCCUPS loops:
##    - drop juicer header/comment row
##    - coerce coords & numeric fields
##    - add chr prefix
##    - parse cell_line (including GSM6086808 -> PDAC1)
## ==============================================================

clean_hiccups_loops <- function(hiccups_loops) {
  # Expect Juicer's '#chr1' column
  if (!("#chr1" %in% colnames(hiccups_loops))) {
    stop("Expected column '#chr1' not found in hiccups_loops.")
  }

  # 1) Drop juicer header / comment lines (they start with "#")
  df <- hiccups_loops[
    !is.na(hiccups_loops$`#chr1`) &
      !startsWith(hiccups_loops$`#chr1`, "#"),
    ,
    drop = FALSE
  ]

  # 2) Coerce key columns to numeric; standardize chr; assign cell_line
  df <- df |>
    dplyr::mutate(
      # standardize chr names
      chr1 = dplyr::if_else(
        startsWith(`#chr1`, "chr"),
        `#chr1`,
        paste0("chr", `#chr1`)
      ),
      chr2 = dplyr::if_else(
        startsWith(chr2, "chr"),
        chr2,
        paste0("chr", chr2)
      ),

      x1 = suppressWarnings(as.integer(x1)),
      x2 = suppressWarnings(as.integer(x2)),
      y1 = suppressWarnings(as.integer(y1)),
      y2 = suppressWarnings(as.integer(y2)),

      observed      = suppressWarnings(as.numeric(observed)),
      expectedBL    = suppressWarnings(as.numeric(expectedBL)),
      expectedDonut = suppressWarnings(as.numeric(expectedDonut)),
      expectedH     = suppressWarnings(as.numeric(expectedH)),
      expectedV     = suppressWarnings(as.numeric(expectedV)),

      fdrBL    = suppressWarnings(as.numeric(fdrBL)),
      fdrDonut = suppressWarnings(as.numeric(fdrDonut)),
      fdrH     = suppressWarnings(as.numeric(fdrH)),
      fdrV     = suppressWarnings(as.numeric(fdrV)),

      numCollapsed = suppressWarnings(as.integer(numCollapsed)),
      centroid1    = suppressWarnings(as.integer(centroid1)),
      centroid2    = suppressWarnings(as.integer(centroid2)),
      radius       = suppressWarnings(as.integer(radius))
    )

  # Base cell line name from sample_id
  base_cell <- sub("\\..*$", "", df$sample_id)

  # 3) Map to a nicer cell_line label
  df$cell_line <- dplyr::case_when(
    # PDAC1 organoid, empty vector Hi-C (GSM6086808) from Kato et al. 2023
    grepl("^GSM6086808", df$sample_id) ~ "PDAC1",
    TRUE ~ base_cell
  )

  # 4) Drop rows with unparsed coords (parsing glitches)
  df <- df |>
    dplyr::filter(
      is.finite(x1), is.finite(x2),
      is.finite(y1), is.finite(y2)
    )

  # 5) Reorder columns to BEDPE-like + metadata
  keep_cols <- c(
    "cell_line", "sample_id", "file_name",
    "chr1", "x1", "x2",
    "chr2", "y1", "y2",
    "name", "score", "strand1", "strand2", "color",
    "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
    "fdrBL", "fdrDonut", "fdrH", "fdrV",
    "numCollapsed", "centroid1", "centroid2", "radius"
  )
  keep_cols <- intersect(keep_cols, colnames(df))

  df[, keep_cols, drop = FALSE]
}

## ==============================================================
## 3) Load GSE149103 loops from FILE PATH (HPNE)
##    and map into the same schema as clean_hiccups_loops()
## ==============================================================

load_gse149103_loops_file <- function(loop_file) {
  if (!file.exists(loop_file)) {
    cli::cli_abort("Could not find file: {.val {loop_file}}")
  }

  message("Reading external loops file: ", basename(loop_file))

  # First try as tab-delimited
  gse_raw <- tryCatch(
    {
      readr::read_delim(
        loop_file,
        delim = "\t",
        col_names = TRUE,
        col_types = readr::cols(.default = readr::col_character()),
        trim_ws = TRUE,
        progress = FALSE,
        show_col_types = FALSE
      )
    },
    error = function(e) {
      message("  read_delim with tab failed (", e$message, "); trying comma...")
      NULL
    }
  )

  # If that failed or produced a single column, retry as comma
  if (is.null(gse_raw) || ncol(gse_raw) == 1L) {
    message("  Detected single-column or failed parse with tab, retrying as CSV/comma...")
    gse_raw <- readr::read_delim(
      loop_file,
      delim = ",",
      col_names = TRUE,
      col_types = readr::cols(.default = readr::col_character()),
      trim_ws = TRUE,
      progress = FALSE,
      show_col_types = FALSE
    )
  }

  if (ncol(gse_raw) < 8L) {
    cli::cli_abort(
      "Expected at least 8 columns (Lchr, Lstart, Lend, Rchr, Rstart, Rend, loop_size, Cell type), got {.val {ncol(gse_raw)}}."
    )
  }

  # Map by position, not by exact header name
  gse_loops <- tibble::tibble(
    Lchr      = gse_raw[[1]],
    Lstart    = suppressWarnings(as.integer(gse_raw[[2]])),
    Lend      = suppressWarnings(as.integer(gse_raw[[3]])),
    Rchr      = gse_raw[[4]],
    Rstart    = suppressWarnings(as.integer(gse_raw[[5]])),
    Rend      = suppressWarnings(as.integer(gse_raw[[6]])),
    loop_size = suppressWarnings(as.integer(gse_raw[[7]])),
    cell_line = as.character(gse_raw[[8]])  # e.g. "HPNE"
  )

  # Map into the same schema as hiccups_loops_clean
  gse_loops <- gse_loops |>
    dplyr::mutate(
      sample_id = paste0("GSE149103_", cell_line),
      file_name = basename(loop_file),

      chr1 = dplyr::if_else(
        startsWith(Lchr, "chr"),
        Lchr,
        paste0("chr", Lchr)
      ),
      chr2 = dplyr::if_else(
        startsWith(Rchr, "chr"),
        Rchr,
        paste0("chr", Rchr)
      ),

      x1 = Lstart,
      x2 = Lend,
      y1 = Rstart,
      y2 = Rend,

      name    = ".",
      score   = ".",
      strand1 = ".",
      strand2 = ".",
      color   = "0,0,0",

      observed      = NA_real_,
      expectedBL    = NA_real_,
      expectedDonut = NA_real_,
      expectedH     = NA_real_,
      expectedV     = NA_real_,

      fdrBL    = NA_real_,
      fdrDonut = NA_real_,
      fdrH     = NA_real_,
      fdrV     = NA_real_,

      numCollapsed = NA_integer_,
      centroid1    = as.integer((Lstart + Lend) / 2),
      centroid2    = as.integer((Rstart + Rend) / 2),
      radius       = NA_integer_
    )

  gse_loops <- gse_loops[, c(
    "cell_line", "sample_id", "file_name",
    "chr1", "x1", "x2",
    "chr2", "y1", "y2",
    "name", "score", "strand1", "strand2", "color",
    "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
    "fdrBL", "fdrDonut", "fdrH", "fdrV",
    "numCollapsed", "centroid1", "centroid2", "radius"
  )]

  message("GSE149103 loops loaded: ", nrow(gse_loops), " rows.")
  gse_loops
}

## ==============================================================
## 4) Run everything & produce a unified loops table
## ==============================================================

hiccups_base_dir <- "/data/homes/yl814/episcope_test/HiC/hiccups_results"

# Juicer HiCCUPS loops: raw + cleaned
hiccups_loops_raw   <- load_hiccups_loops(hiccups_base_dir)
hiccups_loops_clean <- clean_hiccups_loops(hiccups_loops_raw)

# External HPNE loops from GSE149103
gse_file         <- file.path(hiccups_base_dir, "GSE149103_hic_loops_hg38.csv")
gse149103_loops  <- load_gse149103_loops_file(gse_file)

# Unified table: Juicer + GSE149103, with nice cell_line labels (including PDAC1)
hiccups_loops_all <- dplyr::bind_rows(
  hiccups_loops_clean,
  gse149103_loops
)

# sanity checks
dplyr::count(hiccups_loops_all, cell_line)
head(hiccups_loops_all, 10)
unique(hiccups_loops_all$cell_line)

## ==============================================================
## 5) Helper: parse "chr:start-end" peaks to GRanges
## ==============================================================

peak_to_granges <- function(peak_vec) {
  peak_vec <- unique(peak_vec)
  peak_vec <- peak_vec[!is.na(peak_vec)]
  if (!length(peak_vec)) {
    return(GenomicRanges::GRanges())
  }

  parts  <- strsplit(peak_vec, ":", fixed = TRUE)
  chrom  <- vapply(parts, `[`, character(1), 1L)
  range  <- vapply(parts, `[`, character(1), 2L)
  se     <- strsplit(range, "-", fixed = TRUE)
  start  <- as.integer(vapply(se, `[`, character(1), 1L))
  end    <- as.integer(vapply(se, `[`, character(1), 2L))

  ok <- is.finite(start) & is.finite(end) & !is.na(chrom)
  GenomicRanges::GRanges(
    seqnames = chrom[ok],
    ranges   = IRanges::IRanges(start = start[ok], end = end[ok]),
    peak     = peak_vec[ok]
  )
}

## ==============================================================
## 6) Build Hi-C regulatory windows (PP / PE / EE) in gh_std format
##    - TF-agnostic, only uses:
##        * hiccups_loops_all
##        * gene_annot_ref_hg38 (HGNC, chrom, tss)
##        * atac_score_tbl$atac_peak (chr:start-end)
##    - Output columns exactly match gh_std:
##        chrom, start, end, connected_gene,
##        gh_id, confidence, source,
##        feature_name, strand, frame, is_elite_elem
## ==============================================================

build_hic_gh_std <- function(hiccups_loops_all,
                             gene_annot_ref_hg38,
                             atac_score_tbl,
                             promoter_half_width = 2000L,
                             max_enh_gene_dist = 30000L,
                             source_mode = c("hic_cell", "hic"),
                             verbose = TRUE) {
  source_mode <- match.arg(source_mode)

  if (!nrow(hiccups_loops_all)) {
    cli::cli_abort("hiccups_loops_all is empty.")
  }
  required_hic <- c("cell_line", "chr1", "x1", "x2", "chr2", "y1", "y2")
  if (!all(required_hic %in% colnames(hiccups_loops_all))) {
    cli::cli_abort("hiccups_loops_all is missing columns: {.val {setdiff(required_hic, colnames(hiccups_loops_all))}}")
  }

  if (!all(c("HGNC", "ensembl_gene_id", "chrom", "tss") %in% colnames(gene_annot_ref_hg38))) {
    cli::cli_abort("gene_annot_ref_hg38 must contain HGNC, ensembl_gene_id, chrom, tss.")
  }

  if (!("atac_peak" %in% colnames(atac_score_tbl))) {
    cli::cli_abort("atac_score_tbl must contain column 'atac_peak'.")
  }

  ## Promoter GRanges (TSS +/- promoter_half_width)
  gene_tss_ref <- gene_annot_ref_hg38 |>
    dplyr::mutate(is_ensg = startsWith(ensembl_gene_id, "ENSG")) |>
    dplyr::group_by(HGNC) |>
    dplyr::arrange(dplyr::desc(is_ensg), ensembl_gene_id, .by_group = TRUE) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(HGNC, chrom, tss)

  gene_tss_ref <- gene_tss_ref[
    !is.na(gene_tss_ref$HGNC) &
      !is.na(gene_tss_ref$chrom) &
      !is.na(gene_tss_ref$tss),
    ,
    drop = FALSE
  ]

  if (!nrow(gene_tss_ref)) {
    cli::cli_abort("No valid gene TSS rows in gene_annot_ref_hg38 after filtering.")
  }

  promoters_gr <- GenomicRanges::GRanges(
    seqnames = gene_tss_ref$chrom,
    ranges   = IRanges::IRanges(
      start = gene_tss_ref$tss - promoter_half_width,
      end   = gene_tss_ref$tss + promoter_half_width
    ),
    gene_key = gene_tss_ref$HGNC
  )

  ## TSS points for nearest-gene assignment (enhancer-only loops)
  gene_tss_gr <- GenomicRanges::GRanges(
    seqnames = gene_tss_ref$chrom,
    ranges   = IRanges::IRanges(
      start = gene_tss_ref$tss,
      end   = gene_tss_ref$tss
    ),
    gene_key = gene_tss_ref$HGNC
  )

  ## Enhancer GRanges from all ATAC peaks
  peaks_gr <- peak_to_granges(atac_score_tbl$atac_peak)
  if (length(peaks_gr) == 0L) {
    cli::cli_abort("No valid ATAC peaks found in atac_score_tbl$atac_peak.")
  }

  clines <- sort(unique(hiccups_loops_all$cell_line))
  all_out <- list()
  idx <- 1L

  for (cl in clines) {
    loops_sub <- hiccups_loops_all[hiccups_loops_all$cell_line == cl, , drop = FALSE]
    if (!nrow(loops_sub)) next

    if (verbose) {
      message("Building Hi-C gh_std windows for cell_line ", cl,
              " (", nrow(loops_sub), " loops)...")
    }

    loop_ids <- seq_len(nrow(loops_sub))
    loops_sub$loop_id <- loop_ids

    anchor1 <- GenomicRanges::GRanges(
      seqnames = loops_sub$chr1,
      ranges   = IRanges::IRanges(start = loops_sub$x1, end = loops_sub$x2),
      loop_id  = loop_ids
    )

    anchor2 <- GenomicRanges::GRanges(
      seqnames = loops_sub$chr2,
      ranges   = IRanges::IRanges(start = loops_sub$y1, end = loops_sub$y2),
      loop_id  = loop_ids
    )

    ## Overlaps: anchors vs promoters (P)
    ov_g1 <- GenomicRanges::findOverlaps(anchor1, promoters_gr, ignore.strand = TRUE)
    ov_g2 <- GenomicRanges::findOverlaps(anchor2, promoters_gr, ignore.strand = TRUE)

    if (length(ov_g1)) {
      df_g1 <- tibble::tibble(
        anchor_side = "anchor1",
        loop_id     = S4Vectors::mcols(anchor1)$loop_id[S4Vectors::queryHits(ov_g1)],
        gene_key    = S4Vectors::mcols(promoters_gr)$gene_key[S4Vectors::subjectHits(ov_g1)]
      )
    } else {
      df_g1 <- tibble::tibble(
        anchor_side = character(0),
        loop_id     = integer(0),
        gene_key    = character(0)
      )
    }

    if (length(ov_g2)) {
      df_g2 <- tibble::tibble(
        anchor_side = "anchor2",
        loop_id     = S4Vectors::mcols(anchor2)$loop_id[S4Vectors::queryHits(ov_g2)],
        gene_key    = S4Vectors::mcols(promoters_gr)$gene_key[S4Vectors::subjectHits(ov_g2)]
      )
    } else {
      df_g2 <- tibble::tibble(
        anchor_side = character(0),
        loop_id     = integer(0),
        gene_key    = character(0)
      )
    }

    df_gene <- dplyr::bind_rows(df_g1, df_g2)

    ## Overlaps: anchors vs enhancers (E)
    ov_p1 <- GenomicRanges::findOverlaps(anchor1, peaks_gr, ignore.strand = TRUE)
    ov_p2 <- GenomicRanges::findOverlaps(anchor2, peaks_gr, ignore.strand = TRUE)

    if (length(ov_p1)) {
      df_p1 <- tibble::tibble(
        anchor_side = "anchor1",
        loop_id     = S4Vectors::mcols(anchor1)$loop_id[S4Vectors::queryHits(ov_p1)],
        peak        = S4Vectors::mcols(peaks_gr)$peak[S4Vectors::subjectHits(ov_p1)]
      )
    } else {
      df_p1 <- tibble::tibble(
        anchor_side = character(0),
        loop_id     = integer(0),
        peak        = character(0)
      )
    }

    if (length(ov_p2)) {
      df_p2 <- tibble::tibble(
        anchor_side = "anchor2",
        loop_id     = S4Vectors::mcols(anchor2)$loop_id[S4Vectors::queryHits(ov_p2)],
        peak        = S4Vectors::mcols(peaks_gr)$peak[S4Vectors::subjectHits(ov_p2)]
      )
    } else {
      df_p2 <- tibble::tibble(
        anchor_side = character(0),
        loop_id     = integer(0),
        peak        = character(0)
      )
    }

    df_peak <- dplyr::bind_rows(df_p1, df_p2)

    if (!nrow(df_gene) && !nrow(df_peak)) next

    ## Attach genomic coordinates (using anchor coords)
    coord_df <- loops_sub[, c("loop_id", "chr1", "x1", "x2", "chr2", "y1", "y2"),
                          drop = FALSE]

    if (nrow(df_gene)) {
      df_gene <- df_gene |>
        dplyr::left_join(coord_df, by = "loop_id") |>
        dplyr::mutate(
          chrom = dplyr::if_else(anchor_side == "anchor1", chr1, chr2),
          start = dplyr::if_else(anchor_side == "anchor1", x1, y1),
          end   = dplyr::if_else(anchor_side == "anchor1", x2, y2)
        )
    }

    if (nrow(df_peak)) {
      df_peak <- df_peak |>
        dplyr::left_join(coord_df, by = "loop_id") |>
        dplyr::mutate(
          chrom = dplyr::if_else(anchor_side == "anchor1", chr1, chr2),
          start = dplyr::if_else(anchor_side == "anchor1", x1, y1),
          end   = dplyr::if_else(anchor_side == "anchor1", x2, y2)
        )
    }

    ## Map enhancer peaks to loops with promoters (PE / mixed) -----------
    gene_by_loop <- if (nrow(df_gene)) {
      df_gene |>
        dplyr::distinct(loop_id, gene_key)
    } else {
      tibble::tibble(loop_id = integer(0), gene_key = character(0))
    }

    ## Promoter elements: one row per (loop, anchor, gene)
    gh_prom <- tibble::tibble()
    if (nrow(df_gene)) {
      gh_prom <- df_gene |>
        dplyr::transmute(
          chrom          = chrom,
          start          = start,
          end            = end,
          connected_gene = gene_key,
          gh_id          = sprintf("hic_%s:%d-%d", chrom, start, end),
          confidence     = 1,
          source         = if (source_mode == "hic_cell") {
            paste0("hic_", cl)
          } else {
            "hic"
          },
          feature_name  = "Promoter",
          strand        = ".",
          frame         = ".",
          is_elite_elem = 0L
        )
    }

    ## Enhancer elements part 1: loops that also have promoters ----------
    gh_enh_pe <- tibble::tibble()
    if (nrow(df_peak) && nrow(gene_by_loop)) {
      peak_with_gene <- df_peak |>
        dplyr::inner_join(gene_by_loop, by = "loop_id")

      if (nrow(peak_with_gene)) {
        gh_enh_pe <- peak_with_gene |>
          dplyr::transmute(
            chrom          = chrom,
            start          = start,
            end            = end,
            connected_gene = gene_key,
            gh_id          = sprintf("hic_%s:%d-%d", chrom, start, end),
            confidence     = 1,
            source         = if (source_mode == "hic_cell") {
              paste0("hic_", cl)
            } else {
              "hic"
            },
            feature_name  = "Enhancer",
            strand        = ".",
            frame         = ".",
            is_elite_elem = 0L
          )
      }
    }

    ## Enhancer elements part 2: enhancer-enhancer loops (no promoters) --
    gh_enh_ee <- tibble::tibble()
    if (nrow(df_peak)) {
      peak_no_gene <- df_peak |>
        dplyr::anti_join(gene_by_loop, by = "loop_id")

      if (nrow(peak_no_gene)) {
        ## Assign nearest gene TSS to each enhancer anchor center
        ee_centers_gr <- GenomicRanges::GRanges(
          seqnames = peak_no_gene$chrom,
          ranges   = IRanges::IRanges(
            start = as.integer((peak_no_gene$start + peak_no_gene$end) / 2),
            end   = as.integer((peak_no_gene$start + peak_no_gene$end) / 2)
          )
        )

        nearest_idx <- GenomicRanges::nearest(ee_centers_gr, gene_tss_gr, ignore.strand = TRUE)
        nearest_gene <- rep(NA_character_, length(nearest_idx))
        nearest_dist <- rep(NA_integer_, length(nearest_idx))

        has_hit <- !is.na(nearest_idx)
        if (any(has_hit)) {
          nearest_gene[has_hit] <- S4Vectors::mcols(gene_tss_gr)$gene_key[nearest_idx[has_hit]]
          nearest_dist[has_hit] <- GenomicRanges::distance(
            ee_centers_gr[has_hit],
            gene_tss_gr[nearest_idx[has_hit]]
          )
        }

        peak_no_gene$nearest_gene <- nearest_gene
        peak_no_gene$nearest_dist <- nearest_dist

        if (!is.null(max_enh_gene_dist) && is.finite(max_enh_gene_dist)) {
          peak_no_gene <- peak_no_gene[
            !is.na(peak_no_gene$nearest_gene) &
              !is.na(peak_no_gene$nearest_dist) &
              peak_no_gene$nearest_dist <= max_enh_gene_dist,
            ,
            drop = FALSE
          ]
        } else {
          peak_no_gene <- peak_no_gene[
            !is.na(peak_no_gene$nearest_gene),
            ,
            drop = FALSE
          ]
        }

        if (nrow(peak_no_gene)) {
          gh_enh_ee <- peak_no_gene |>
            dplyr::transmute(
              chrom          = chrom,
              start          = start,
              end            = end,
              connected_gene = nearest_gene,
              gh_id          = sprintf("hic_%s:%d-%d", chrom, start, end),
              confidence     = 1,
              source         = if (source_mode == "hic_cell") {
                paste0("hic_", cl)
              } else {
                "hic"
              },
              feature_name  = "Enhancer",
              strand        = ".",
              frame         = ".",
              is_elite_elem = 0L
            )
        }
      }
    }

    gh_cl <- dplyr::bind_rows(gh_prom, gh_enh_pe, gh_enh_ee)

    if (!nrow(gh_cl)) next

    ## Deduplicate per cell line
    gh_cl <- gh_cl |>
      dplyr::distinct(
        chrom, start, end, connected_gene,
        gh_id, confidence, source,
        feature_name, strand, frame, is_elite_elem
      )

    all_out[[idx]] <- gh_cl
    idx <- idx + 1L
  }

  if (!length(all_out)) {
    cli::cli_warn("build_hic_gh_std: no Hi-C windows found.")
    return(
      tibble::tibble(
        chrom          = character(0),
        start          = integer(0),
        end            = integer(0),
        connected_gene = character(0),
        gh_id          = character(0),
        confidence     = numeric(0),
        source         = character(0),
        feature_name   = character(0),
        strand         = character(0),
        frame          = character(0),
        is_elite_elem  = integer(0)
      )
    )
  }

  gh_std_hic <- dplyr::bind_rows(all_out)

  if (verbose) {
    message("Hi-C gh_std windows built: ", nrow(gh_std_hic), " rows.")
    message("Preview (first 10 rows):")
    print(utils::head(gh_std_hic, 10))
  }

  gh_std_hic
}

## ==============================================================
## 7) usage
## ==============================================================

## input
##   - gene_annot_ref_hg38  (HGNC / chrom / tss)
##   - grn_set$atac_score   (with column 'atac_peak')

gh_std_hic <- build_hic_gh_std(
  hiccups_loops_all   = hiccups_loops_all,
  gene_annot_ref_hg38 = gene_annot_ref_hg38,
  atac_score_tbl      = grn_set$atac_score,
  promoter_half_width = 2000L,
  max_enh_gene_dist   = Inf,
  source_mode         = "hic_cell",
  verbose             = TRUE
)

unique(gh_std_hic$source)

src_vec <- unique(gh_std_hic$source)
src_short <- sub("^hic_", "", src_vec)

gh_std_hic_list <- setNames(
  lapply(src_vec, function(s) {
    gh_std_hic[gh_std_hic$source == s, , drop = FALSE]
  }),
  paste0("gh_std_hic_", src_short)
)


list2env(gh_std_hic_list, envir = .GlobalEnv)

gh_std_hic_BxPC3
`gh_std_hic_Capan-1`
gh_std_hic_HPDE6C7
gh_std_hic_HPNE
`gh_std_hic_PANC-1`
gh_std_hic_PDAC1


