library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(valr)
library(data.table)
# ---- helpers --------------------------------------------------------------
.stop_if_missing <- function(p) if (!base::file.exists(p)) cli::cli_abort("File not found: {p}")
.read_tsv_header <- function(path) {
  .stop_if_missing(path)
  df <- try(
    readr::read_tsv(
      path,
      comment   = "#",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      progress  = FALSE
    ),
    silent = TRUE
  )
  if (inherits(df, "try-error")) {
    cli::cli_abort("Could not read as headered TSV: {path}. This function expects an overview table with column names.")
  }
  df
}
.as_bed3 <- function(df, c_chrom, c_start, c_end) {
  bed <- tibble::tibble(
    chrom = base::as.character(df[[c_chrom]]),
    start = base::as.integer(df[[c_start]]),
    end   = base::as.integer(df[[c_end]])
  )
  class(bed) <- c("tbl_df", "tbl", "data.frame", "bed_frame")
  bed
}
.collapse_unique_sites <- function(bed) bed[!base::duplicated(bed[c("chrom", "start", "end")]), , drop = FALSE]
.metrics_basic <- function(tp, fp, fn) {
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1 <- if (!base::is.na(precision) && !base::is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  tibble::tibble(precision = precision, recall = recall, F1 = f1)
}
.metrics_tn <- function(tp, fp, fn, tn) {
  tp <- base::as.numeric(tp); fp <- base::as.numeric(fp)
  fn <- base::as.numeric(fn); tn <- base::as.numeric(tn)
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  denom_mcc <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
  num_mcc   <- (tp * tn) - (fp * fn)
  mcc <- if (base::is.finite(denom_mcc) && denom_mcc > 0) num_mcc / base::sqrt(denom_mcc) else NA_real_
  tibble::tibble(specificity = specificity, MCC = mcc)
}
get_srr <- function(path) {
  fname <- basename(path)
  sub("\\..*$", "", fname)
}

# helper: from a set of intervals (bed_x), return those that have ≥50% of
# THEIR length overlapped by any ChIP peak (used for ATAC peaks -> TP).
.truth_from_overlap_half_atac <- function(bed_x, chip_peaks) {
  ov <- valr::bed_intersect(bed_x, chip_peaks)
  if (!base::nrow(ov)) {
    return(bed_x[0, , drop = FALSE])
  }
  start_x <- base::as.integer(ov[["start.x"]])
  end_x   <- base::as.integer(ov[["end.x"]])
  start_y <- base::as.integer(ov[["start.y"]])
  end_y   <- base::as.integer(ov[["end.y"]])
  len_x   <- base::pmax(1L, end_x - start_x)
  inter_w <- base::pmax(0L, base::pmin(end_x, end_y) - base::pmax(start_x, start_y))
  frac    <- inter_w / len_x
  ov      <- ov[frac >= 0.5, , drop = FALSE]
  if (!base::nrow(ov)) {
    return(bed_x[0, , drop = FALSE])
  }
  tibble::tibble(
    chrom = ov[["chrom"]],
    start = base::as.integer(ov[["start.x"]]),
    end   = base::as.integer(ov[["end.x"]])
  ) |>
    .collapse_unique_sites()
}
.read_chip_bed3 <- function(path) {
  df <- .read_tsv_header(path)
  if (base::ncol(df) < 3L) cli::cli_abort("ChIP file must have at least 3 columns: {path}")
  base::names(df)[1:3] <- c("chrom", "start", "end")
  df$start <- suppressWarnings(base::as.integer(df$start))
  df$end   <- suppressWarnings(base::as.integer(df$end))
  if (base::any(!is.finite(df$start) | !is.finite(df$end))) {
    cli::cli_abort("Non-integer ChIP start/end in: {path}")
  }
  bad <- base::which(df$end < df$start)
  if (base::length(bad)) {
    s <- df$start[bad]; e <- df$end[bad]
    df$start[bad] <- base::pmin(s, e); df$end[bad] <- base::pmax(s, e)
  }
  bed <- tibble::as_tibble(df[, c("chrom", "start", "end"), drop = FALSE])
  class(bed) <- c("tbl_df", "tbl", "data.frame", "bed_frame")
  bed
}

# Annotate ALL ATAC peaks with overlap stats + chip-bound flags
annotate_atac_with_chip50 <- function(pred_bed,
                                      chip_bed,
                                      out_csv = NULL,
                                      overlap_mode  = c("fraction", "bp"),
                                      overlap_frac  = 0.5,
                                      overlap_bp    = 1L) {
  # 1) Read overview (needs peak_chr / peak_start / peak_end)
  pred_raw <- readr::read_tsv(pred_bed, show_col_types = FALSE)

  required_cols <- c("peak_chr", "peak_start", "peak_end")
  if (!all(required_cols %in% names(pred_raw))) {
    missing <- setdiff(required_cols, names(pred_raw))
    cli::cli_abort("Missing required columns in {pred_bed}: {paste(missing, collapse = ', ')}")
  }

  # Unique ATAC peaks
  atac_all <- pred_raw[, required_cols, drop = FALSE] |>
    dplyr::distinct()

  # 2) Read ChIP BED (first 3 columns: chrom/start/end)
  chip_df <- .read_chip_bed3(chip_bed)

  if (ncol(chip_df) < 3L) {
    cli::cli_abort("ChIP file must have at least 3 columns: {chip_bed}")
  }
  colnames(chip_df)[1:3] <- c("chrom", "start", "end")
  chip_df$start <- as.integer(chip_df$start)
  chip_df$end   <- as.integer(chip_df$end)

  chip_peaks <- tibble::as_tibble(chip_df[, c("chrom","start","end"), drop = FALSE])
  class(chip_peaks) <- c("tbl_df","tbl","data.frame","bed_frame")

  # 3) Convert ATAC peaks to bed_frame
  atac_bed <- tibble::tibble(
    chrom = as.character(atac_all$peak_chr),
    start = as.integer(atac_all$peak_start),
    end   = as.integer(atac_all$peak_end)
  )
  class(atac_bed) <- c("tbl_df","tbl","data.frame","bed_frame")

  # 4) Intersect ATAC peaks with ChIP peaks
  ov <- valr::bed_intersect(atac_bed, chip_peaks)

  overlap_mode <- match.arg(overlap_mode)

  if (nrow(ov) > 0L) {

    ov2 <- ov |>
      dplyr::mutate(
        atac_len  = pmax(1L, end.x - start.x),
        frac_atac = .overlap / atac_len) 

    if (overlap_mode == "fraction") {
      ov2 <- ov2 |> dplyr::mutate(chip_hit = frac_atac >= overlap_frac)
    } else if (overlap_mode == "bp") {
      ov2 <- ov2 |> dplyr::mutate(chip_hit = .overlap  >= overlap_bp)
    }
    ov2 <- ov2 |> 
      dplyr::group_by(chrom, start.x, end.x) |>
      dplyr::slice_max(frac_atac, with_ties = FALSE) |>
      dplyr::ungroup()
    
    annot <- ov2 |>
      dplyr::transmute(
        peak_chr   = chrom,
        peak_start = start.x,
        peak_end   = end.x,
        chip_overlap_bp = as.integer(.overlap),
        chip_overlap_frac = frac_atac,
        chip_bound = ifelse(chip_hit, 1L, 0L),
        chip_chr   = dplyr::if_else(chip_hit, chrom,   NA_character_),
        chip_start = dplyr::if_else(chip_hit, start.y, NA_integer_),
        chip_end   = dplyr::if_else(chip_hit, end.y,   NA_integer_)
      )
  } else {
    
    # No overlaps at all
    annot <- atac_all |>
      dplyr::mutate(
        chip_overlap_bp = 0L,
        chip_overlap_frac = 0,
        chip_bound = 0L,
        chip_chr   = NA_character_,
        chip_start = NA_integer_,
        chip_end   = NA_integer_
      )
  }

  # 5) Left-join back onto all ATAC peaks, fill non-overlapped as chip_bound = 0
  atac_annot <- atac_all |>
    dplyr::left_join(
      annot,
      by = c("peak_chr","peak_start","peak_end")
    ) |>
    dplyr::mutate(
      chip_overlap_bp = dplyr::coalesce(chip_overlap_bp, 0L),
      chip_overlap_frac = dplyr::coalesce(chip_overlap_frac, 0),
      chip_bound = dplyr::coalesce(chip_bound, 0L),
      chip_chr   = dplyr::coalesce(chip_chr, NA_character_),
      chip_start = dplyr::coalesce(chip_start, NA_integer_),
      chip_end   = dplyr::coalesce(chip_end,   NA_integer_),
      chip_peak  = dplyr::if_else(
        chip_bound == 1L,
        sprintf("%s:%d-%d", chip_chr, chip_start, chip_end),
        NA_character_
      )
    ) |>
    dplyr::mutate(
      chip_bound_TP = if_else(chip_overlap_bp >= 1L, 1L, 0L),
      chip_bound_FN = if_else(chip_overlap_frac >= overlap_frac, 1L, 0L),
      chip_bound_FP = dplyr::case_when(
        chip_overlap_frac >= overlap_frac ~ 1L,
        chip_overlap_bp == 0L ~ 0L,
        TRUE ~ NA_integer_
      ),
      chip_bound_TN = dplyr::case_when(
        chip_overlap_frac >= overlap_frac ~ 1L,
        chip_overlap_bp == 0L ~ 0L,
        TRUE ~ NA_integer_
      )
    ) |>
    # keep only join keys + chip annotation to avoid duplicate columns on join
    dplyr::select(
      peak_chr, peak_start, peak_end,
      chip_overlap_bp, chip_overlap_frac,
      chip_bound, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN,
      chip_chr, chip_start, chip_end, chip_peak
    )

  # 6) Join annotation back to *every row* of pred_raw
  pred_with_chip <- pred_raw |>
    dplyr::left_join(
      atac_annot,
      by = c("peak_chr","peak_start","peak_end")
    )

  # 7) Optionally write CSV with all original columns + chip_* columns
  if (!is.null(out_csv)) {
    readr::write_csv(pred_with_chip, out_csv)
  }

  pred_with_chip
}

# x <- batch[[1]]
# pred_bed = x$pred
# chip_bed = x$chip
# out_csv  = out_csv
#
# readr::write_tsv(
#   atac_bed,
#   "atac.bed",
#   col_names = FALSE
# )
# readr::write_tsv(
#   chip_peaks,
#   "chip.bed",
#   col_names = FALSE
# )


# process all TFs
base_pred_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/pearson_vs_spearman"
base_exp_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/pearson_vs_spearman"


# TFs and their corresponding CUT&Tag files
chip_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/cutntag"

# batch  <- list(
#
#   list(tf="FOXA2", pred=file.path(base_pred_dir,"predicted_all_tfbs","FOXA2_overview.txt"),
#        chip=file.path(chip_dir,"SRR7826319.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="SOX9", pred=file.path(base_pred_dir,"predicted_all_tfbs","SOX9_overview.txt"),
#        chip=file.path(chip_dir,"SRR13790201.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="CTCF", pred=file.path(base_pred_dir,"predicted_all_tfbs","CTCF_overview.txt"),
#        chip=file.path(chip_dir,"SRR6212679.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="CTCF", pred=file.path(base_pred_dir,"predicted_all_tfbs","CTCF_overview.txt"),
#        chip=file.path(chip_dir,"SRR6212718.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="CTCF", pred=file.path(base_pred_dir,"predicted_all_tfbs","CTCF_overview.txt"),
#        chip=file.path(chip_dir,"SRR6212719.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="CTCF", pred=file.path(base_pred_dir,"predicted_all_tfbs","CTCF_overview.txt"),
#        chip=file.path(chip_dir,"SRR6212680.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="KLF6", pred=file.path(base_pred_dir,"predicted_all_tfbs","KLF6_overview.txt"),
#        chip=file.path(chip_dir,"SRR2834280.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="KLF6", pred=file.path(base_pred_dir,"predicted_all_tfbs","KLF6_overview.txt"),
#        chip=file.path(chip_dir,"SRR2834279.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="FOSL1", pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL1_overview.txt"),
#        chip=file.path(chip_dir,"SRR30594323.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="FOSL1", pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL1_overview.txt"),
#        chip=file.path(chip_dir,"SRR30594324.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="ELF3", pred=file.path(base_pred_dir,"predicted_all_tfbs","ELF3_overview.txt"),
#        chip=file.path(chip_dir,"SRR25419825.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="ELF3", pred=file.path(base_pred_dir,"predicted_all_tfbs","ELF3_overview.txt"),
#        chip=file.path(chip_dir,"SRR25419826.hg38.macs3_peaks.filt.narrowPeak")),
#
#   list(tf="MYB", pred=file.path(base_pred_dir,"predicted_all_tfbs","MYB_overview.txt"),
#        chip=file.path(chip_dir,"SRR32816018.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="MYB", pred=file.path(base_pred_dir,"predicted_all_tfbs","MYB_overview.txt"),
#        chip=file.path(chip_dir,"SRR32816019.hg38.macs3_peaks.filt.narrowPeak")),
#   list(tf="MYB", pred=file.path(base_pred_dir,"predicted_all_tfbs","MYB_overview.txt"),
#        chip=file.path(chip_dir,"SRR32816020.hg38.macs3_peaks.filt.narrowPeak"))
# )
#


# batch <- list(
#
#   ## FOSL2 — PRJNA987138
#   list(tf="FOSL2",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL2_overview.txt"),
#        chip=file.path(chip_dir,"SRR26038772.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="FOSL2",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL2_overview.txt"),
#        chip=file.path(chip_dir,"SRR26038778.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="FOSL2",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL2_overview.txt"),
#        chip=file.path(chip_dir,"SRR26038783.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="FOSL2",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL2_overview.txt"),
#        chip=file.path(chip_dir,"SRR26038789.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="FOSL2",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","FOSL2_overview.txt"),
#        chip=file.path(chip_dir,"SRR26038819.hg38.rp10m.narrowpeaks.bed")),
#
#
#   ## JUN — GSE264148
#   list(tf="JUN",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","JUN_overview.txt"),
#        chip=file.path(chip_dir,"SRR28709070.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="JUN",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","JUN_overview.txt"),
#        chip=file.path(chip_dir,"SRR28709071.hg38.rp10m.narrowpeaks.bed")),
#
#   list(tf="JUN",
#        pred=file.path(base_pred_dir,"predicted_all_tfbs","JUN_overview.txt"),
#        chip=file.path(chip_dir,"SRR28709072.hg38.rp10m.narrowpeaks.bed"))
# )

# base_dir <- base_pred_dir
# # TFs and their corresponding ChIP files
# batch <- list(
#   list(
#     tf   = "HNF4A",
#     pred = file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview.txt"),
#     chip = file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")
#   ),
#   list(
#     tf   = "MAFF",
#     pred = file.path(base_dir, "predicted_all_tfbs", "MAFF_overview.txt"),
#     chip = file.path(base_dir, "cutntag", "cy84.hg38.rp10m.narrowpeaks.bed")
#   ),
#   list(
#     tf   = "ZEB1",
#     pred = file.path(base_dir, "predicted_all_tfbs", "ZEB1_overview.txt"),
#     chip = file.path(base_dir, "cutntag", "cy76.hg38.rp10m.narrowpeaks.bed")
#   )
# )


# out_csv <- file.path(
#   base_exp_dir,
#   paste0(x$tf, get_srr(x$chip), "_ATAC_chip50_annotation_fullrows.csv")
# )


batch <- list(
  list(
    tf = "HNF4A",
    method = "pearson",
    pred_all = file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
    pred_canonical = file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
    chip = file.path(chip_dir, "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "HNF4A",
    method = "spearman",
    pred_all = file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt"),
    pred_canonical = file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt"),
    chip = file.path(chip_dir, "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "HNF4A",
    method = "pearson_or_spearman",
    pred_all = c(
      file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
      file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt")
    ),
    pred_canonical = c(
      file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
      file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt")
    ),
    chip = file.path(chip_dir, "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "HNF4A",
    method = "pearson_and_spearman",
    pred_all = c(
      file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
      file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt")
    ),
    pred_canonical = c(
      file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_HNF4A_overview.txt"),
      file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_HNF4A_overview.txt")
    ),
    chip = file.path(chip_dir, "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "MAFF",
    method = "pearson",
    pred_all = file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
    pred_canonical = file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
    chip = file.path(chip_dir, "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "MAFF",
    method = "spearman",
    pred_all = file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt"),
    pred_canonical = file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt"),
    chip = file.path(chip_dir, "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "MAFF",
    method = "pearson_or_spearman",
    pred_all = c(
      file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
      file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt")
    ),
    pred_canonical = c(
      file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
      file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt")
    ),
    chip = file.path(chip_dir, "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf = "MAFF",
    method = "pearson_and_spearman",
    pred_all = c(
      file.path(base_pred_dir, "predicted_all_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
      file.path(base_pred_dir, "predicted_all_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt")
    ),
    pred_canonical = c(
      file.path(base_pred_dir, "predicted_canonical_tfbs_pearson", "HPAFII_10_FBS_MAFF_overview.txt"),
      file.path(base_pred_dir, "predicted_canonical_tfbs_spearman", "HPAFII_10_FBS_MAFF_overview.txt")
    ),
    chip = file.path(chip_dir, "cy84.hg38.rp10m.narrowpeaks.bed")
  )
)


atac_chip_annot_list <- lapply(batch, function(x) {
  message("Annotating: ", x$tf, " ", x$method)

  out_csv <- file.path(
    base_exp_dir,
    paste0(x$tf, "_", x$method, "_ATAC_chip50_annotation_fullrows.csv")
  )

  if (is.null(x$chip)) {
    pred_raw <- readr::read_tsv(x$pred_all, show_col_types = FALSE)
    if (!"_bound" %in% names(pred_raw)) {
      cli::cli_abort("Missing '_bound' column in {x$pred_all}; provide a ChIP file or add '_bound'.")
    }
    pred_raw$chip_bound <- base::as.integer(pred_raw[["_bound"]])
    pred_raw$chip_bound_TP <- pred_raw$chip_bound
    pred_raw$chip_bound_FN <- pred_raw$chip_bound
    pred_raw$chip_bound_FP <- pred_raw$chip_bound
    pred_raw$chip_bound_TN <- pred_raw$chip_bound
    if (!is.null(out_csv)) {
      readr::write_csv(pred_raw, out_csv)
    }
    pred_raw
  } else if (!x$method %in% c("pearson_or_spearman", "pearson_and_spearman")) {
    annotate_atac_with_chip50(
      pred_bed = x$pred_all,
      chip_bed = x$chip,
      out_csv  = out_csv,
      overlap_mode = "bp",
      overlap_bp    = 1L
    )
  } else {
    pred_raw <- readr::read_tsv(x$pred_all[[1]], show_col_types = FALSE)
    pred_raw <- annotate_atac_with_chip50(
      pred_bed = x$pred_all[[1]],
      chip_bed = x$chip,
      out_csv  = out_csv,
      overlap_mode = "bp",
      overlap_bp    = 1L
    )
    pred_raw
  }
})

names(atac_chip_annot_list) <- vapply(batch, `[[`, character(1), "tf")

atac_overlap <- fread(
  "/data/homes/yl814/episcope_test/nutrient_stress/All_ATAC.master_table.narrowpeaks.10mil.txt",
  select = c("Chr", "Start", "End", "Overlap_cy423")
)

filtered_atac_chip_annot_list <- list()

for (i in seq_along(atac_chip_annot_list)) {

  tf_df <- atac_chip_annot_list[[i]]

  # Only apply ATAC overlap filtering for PANC1
    # LEFT JOIN ATAC overlap
    tf_df_joined <- tf_df %>%
      left_join(
        atac_overlap,
        by = c(
          "peak_chr"   = "Chr",
          "peak_start" = "Start",
          "peak_end"   = "End"
        )
      )

    # Keep only overlapping ATAC peaks
    tf_df_filtered <- tf_df_joined %>%
      filter(Overlap_cy423 == 1)

    filtered_atac_chip_annot_list[[i]] <- tf_df_filtered
}



# New metrics calculation -------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

r_grid  = c(-1, 0, 0.3, 0.5, 0.7)
r_grid  = c(-1, 0, 0.3, 0.5, 0.7)
p_grid  = 10^seq(-4, -1, length.out = 20)

modes <- c("AllFP", "CanonicalFP")
eval_modes <- c("legacy", "fixed_binary")


stack_range <- function(df, comp_levels) {
  if (nrow(df) == 0L) return(c(0, 0))
  tmp <- df %>%
    mutate(component = factor(component, levels = comp_levels)) %>%
    arrange(r, p, component) %>%
    group_by(r, p) %>%
    summarise(
      ymin = min(c(0, cumsum(count_signed))),
      ymax = max(c(0, cumsum(count_signed))),
      .groups = "drop"
    )
  range(c(tmp$ymin, tmp$ymax), na.rm = TRUE)
}

.count_confusion <- function(predicted_bound, predicted_unbound_eval) {
  TP <- predicted_bound %>%
    filter(chip_bound_TP == 1) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  FP <- predicted_bound %>%
    filter(!is.na(chip_bound_FP) & chip_bound_FP == 0) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  TN <- predicted_unbound_eval %>%
    filter(!is.na(chip_bound_TN) & chip_bound_TN == 0) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  FN <- predicted_unbound_eval %>%
    filter(chip_bound_FN == 1) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  c(TP = TP, FP = FP, TN = TN, FN = FN)
}

.count_confusion_fixed_binary <- function(predicted_bound, predicted_unbound_eval) {
  TP <- predicted_bound %>%
    filter(chip_bound_FN == 1) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  FP <- predicted_bound %>%
    filter(chip_bound_FN == 0) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  TN <- predicted_unbound_eval %>%
    filter(chip_bound_FN == 0) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  FN <- predicted_unbound_eval %>%
    filter(chip_bound_FN == 1) %>%
    select(peak_chr, peak_start, peak_end) %>%
    distinct() %>%
    nrow()
  c(TP = TP, FP = FP, TN = TN, FN = FN)
}

for (i in seq_along(batch)) {

  tf     <- batch[[i]]$tf
  method <- batch[[i]]$method
  chip   <- batch[[i]]$chip
  pred_all <- batch[[i]]$pred_all
  pred_canonical <- batch[[i]]$pred_canonical

  message("Running cutoff grid for: ", tf)

  df <- filtered_atac_chip_annot_list[[i]]
  # ALL FP (all annotated peaks)
  df1 <- df %>% distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE)

  if (method %in% c("pearson_or_spearman", "pearson_and_spearman")) {
    df_all_p <- readr::read_tsv(pred_all[[1]], show_col_types = FALSE)
    df_all_s <- readr::read_tsv(pred_all[[2]], show_col_types = FALSE)
    df1_all <- bind_rows(df_all_p, df_all_s) %>%
      inner_join(
        df1 %>% select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN),
        by = c("peak_chr", "peak_start", "peak_end")
      ) %>%
      mutate(
        chip_bound_TP = dplyr::coalesce(chip_bound_TP, 0L),
        chip_bound_FN = dplyr::coalesce(chip_bound_FN, 0L)
      ) %>%
      distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE)

    df_can_p <- readr::read_tsv(pred_canonical[[1]], show_col_types = FALSE)
    df_can_s <- readr::read_tsv(pred_canonical[[2]], show_col_types = FALSE)
    df2 <- bind_rows(df_can_p, df_can_s) %>%
      inner_join(
        df1_all %>% select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN),
        by = c("peak_chr", "peak_start", "peak_end")
      ) %>%
      mutate(
        chip_bound_TP = dplyr::coalesce(chip_bound_TP, 0L),
        chip_bound_FN = dplyr::coalesce(chip_bound_FN, 0L)
      ) %>%
      distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE)
  } else {
    # CANONICAL FP (canonical overview table + chip_bound from all peaks)
    df2 <- readr::read_tsv(pred_canonical, show_col_types = FALSE) %>%
      inner_join(
        df1 %>% select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN),
        by = c("peak_chr", "peak_start", "peak_end")
      ) %>%
      mutate(
        chip_bound_TP = dplyr::coalesce(chip_bound_TP, 0L),
        chip_bound_FN = dplyr::coalesce(chip_bound_FN, 0L)
      ) %>%
      distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE)
  }

  for (eval_mode in eval_modes) {
    message("Running evaluation mode: ", eval_mode)
    grid <- data.frame()
    k <- 1
    for (r in r_grid) {
      for(p in p_grid) {
        print(paste(r,p))

        grid[k, "r"] <- r
        grid[k, "p"] <- p

        for (mode in modes) {

        if (mode==modes[1]) {
          if (method %in% c("pearson_or_spearman", "pearson_and_spearman")) {
            pred_p <- df_all_p %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>%
              select(peak_chr, peak_start, peak_end)
            pred_s <- df_all_s %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>%
              select(peak_chr, peak_start, peak_end)
            predicted_bound <- if (method == "pearson_and_spearman") {
              inner_join(pred_p, pred_s, by = c("peak_chr", "peak_start", "peak_end"))
            } else {
              bind_rows(pred_p, pred_s) %>% distinct
            }
            predicted_bound <- predicted_bound %>%
              distinct %>%
              left_join(
                df1_all %>% select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN),
                by = c("peak_chr", "peak_start", "peak_end")
              ) %>%
              mutate(
                chip_bound_TP = dplyr::coalesce(chip_bound_TP, 0L),
                chip_bound_FN = dplyr::coalesce(chip_bound_FN, 0L)
              )
            predicted_unbound_fn <- anti_join(df1_all, predicted_bound,
                                              by=c("peak_chr", "peak_start", "peak_end")) %>% distinct
          } else {
            predicted_bound <- df1 %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>% distinct
            predicted_unbound_fn <- anti_join(df1,predicted_bound, by=c("peak_chr", "peak_start", "peak_end")) %>% distinct
          }
        } else {
          if (method %in% c("pearson_or_spearman", "pearson_and_spearman")) {
            pred_p <- df_can_p %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>%
              select(peak_chr, peak_start, peak_end)
            pred_s <- df_can_s %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>%
              select(peak_chr, peak_start, peak_end)
            predicted_bound <- if (method == "pearson_and_spearman") {
              inner_join(pred_p, pred_s, by = c("peak_chr", "peak_start", "peak_end"))
            } else {
              bind_rows(pred_p, pred_s) %>% distinct
            }
            predicted_bound <- predicted_bound %>%
              distinct %>%
              left_join(
                df2 %>% select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN),
                by = c("peak_chr", "peak_start", "peak_end")
              ) %>%
              mutate(
                chip_bound_TP = dplyr::coalesce(chip_bound_TP, 0L),
                chip_bound_FN = dplyr::coalesce(chip_bound_FN, 0L)
              )
            predicted_unbound_fn <- anti_join(df2, predicted_bound,
                                              by=c("peak_chr", "peak_start", "peak_end")) %>% distinct
          } else {
            predicted_bound <- df2 %>% filter(corr_fp_tf_r>r & corr_fp_tf_p_adj<p) %>% distinct
            predicted_unbound_fn <- anti_join(df2,predicted_bound, by=c("peak_chr", "peak_start", "peak_end")) %>% distinct
          }
        }

          predicted_unbound_eval <- predicted_unbound_fn %>%
            distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE) %>%
            filter(is.finite(corr_fp_tf_r) & is.finite(corr_fp_tf_p_adj))

          counts <- if (eval_mode == "fixed_binary") {
            .count_confusion_fixed_binary(predicted_bound, predicted_unbound_eval)
          } else {
            .count_confusion(predicted_bound, predicted_unbound_eval)
          }
          TP <- counts[["TP"]]; FP <- counts[["FP"]]
          TN <- counts[["TN"]]; FN <- counts[["FN"]]

          grid[k, paste0(mode, "_TP")] <- TP
          grid[k, paste0(mode, "_FP")] <- FP
          grid[k, paste0(mode, "_TN")] <- TN
          grid[k, paste0(mode, "_FN")] <- FN
          grid[k, paste0(mode, "_accuracy")] <- (TP+TN)/(FP+TN+TP+FN)
          grid[k, paste0(mode, "_precision")] <- TP/(TP+FP)
          grid[k, paste0(mode, "_recall")] <- TP/(TP+FN)
        }
        k=k+1
      }
    }

    # save the raw grid table
    grid_file <- file.path(
      base_exp_dir,
      if (eval_mode == "legacy") {
        paste0(tf, "_", method, "_cutoff_grid_fp.csv")
      } else {
        paste0(tf, "_", method, "_cutoff_grid_fp_", eval_mode, ".csv")
      }
    )
    readr::write_csv(grid, grid_file)
  }
}




#. grid <- read_csv(grid_file)

r_vals <- sort(unique(grid$r))
r_labs <- paste0("r > ", r_vals)

for (i in seq_along(batch)) {

  tf     <- batch[[i]]$tf
  method <- batch[[i]]$method
  chip   <- batch[[i]]$chip
  pred_all <- batch[[i]]$pred_all
  pred_canonical <- batch[[i]]$pred_canonical

  for (eval_mode in eval_modes) {
    metrics_df <- readr::read_csv(
      file.path(
        base_exp_dir,
        if (eval_mode == "legacy") {
          paste0(tf, "_", method, "_cutoff_grid_fp.csv")
        } else {
          paste0(tf, "_", method, "_cutoff_grid_fp_", eval_mode, ".csv")
        }
      ),
      show_col_types = FALSE
    )

  r_vals <- sort(unique(metrics_df$r))
  r_labs <- paste0("r > ", r_vals)
  x_breaks <- scales::trans_breaks("log10", function(x) 10^x)
  x_labels <- scales::trans_format("log10", scales::math_format(bold(10^.x)))

  plot_df <- metrics_df %>%
    select(
      r, p,
      AllFP_accuracy, AllFP_precision, AllFP_recall,
      CanonicalFP_accuracy, CanonicalFP_precision, CanonicalFP_recall
    ) %>%
    pivot_longer(
      cols = -c(r, p),
      names_to = c("line_group", "metric"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    mutate(
      r_lab = factor(r, levels = r_vals, labels = r_labs),
      line_group = recode(line_group, AllFP = "All_FP", CanonicalFP = "Canonical_FP"),
      line_group = factor(line_group, levels = c("All_FP", "Canonical_FP")),
      metric = factor(metric, levels = c("accuracy", "precision", "recall"))
    )

  p_metrics <- ggplot(plot_df, aes(x = p, y = value, color = metric, linetype = line_group)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    labs(
      title = paste0(tf, " cutoff optimization metrics"),
      x = "p-value cutoff",
      y = "Metric value",
      color = "color_group",
      linetype = "line_group"
    ) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_color_manual(values = c(accuracy = "green3", precision = "red2", recall = "blue3")) +
    scale_linetype_manual(values = c(All_FP = "solid", Canonical_FP = "dashed")) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    ) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

  counts_df <- metrics_df %>%
    select(
      r, p,
      AllFP_TP, AllFP_FP, AllFP_TN, AllFP_FN,
      CanonicalFP_TP, CanonicalFP_FP, CanonicalFP_TN, CanonicalFP_FN
    ) %>%
    pivot_longer(
      cols = -c(r, p),
      names_to = c("line_group", "component"),
      names_sep = "_",
      values_to = "count"
    ) %>%
    mutate(
      r_lab = factor(r, levels = r_vals, labels = r_labs),
      line_group = recode(line_group, AllFP = "All_FP", CanonicalFP = "Canonical_FP"),
      component = factor(component, levels = c("FP", "TN", "FN", "TP")),
      count_signed = ifelse(component %in% c("TP", "FN"), count, -count)
    )

  conf_fill <- c(TP = "#31a354", FP = "#de2d26", TN = "#3182bd", FN = "#ff9800")

  counts_all <- counts_df %>% filter(line_group == "All_FP")
  counts_can <- counts_df %>% filter(line_group == "Canonical_FP")

  all_neg <- counts_all %>% filter(component %in% c("FP", "TN")) %>% mutate(component = factor(component, levels = c("FP", "TN")))
  all_pos <- counts_all %>% filter(component %in% c("FN", "TP")) %>% mutate(component = factor(component, levels = c("FN", "TP")))

  can_neg <- counts_can %>% filter(component %in% c("FP", "TN")) %>% mutate(component = factor(component, levels = c("FP", "TN")))
  can_pos <- counts_can %>% filter(component %in% c("FN", "TP")) %>% mutate(component = factor(component, levels = c("FN", "TP")))

  y_lim_all <- range(
    stack_range(all_neg, c("TN", "FP")),
    stack_range(all_pos, c("FN", "TP")),
    0
  )
  y_lim_can <- range(
    stack_range(can_neg, c("TN", "FP")),
    stack_range(can_pos, c("FN", "TP")),
    0
  )
  y_lim_shared <- range(c(y_lim_all, y_lim_can), na.rm = TRUE)

  p_all_counts_fixed <- ggplot(counts_all, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = all_neg) +
    geom_col(data = all_pos) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = y_lim_shared) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p-value cutoff", y = "All (|y| = count)") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

  p_can_counts_fixed <- ggplot(counts_can, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = can_neg) +
    geom_col(data = can_pos) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = y_lim_shared) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p-value cutoff", y = "Canonical (|y| = count)") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )

  p_all_counts_auto <- ggplot(counts_all, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = all_neg) +
    geom_col(data = all_pos) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p-value cutoff", y = "All (|y| = count)") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

  p_can_counts_auto <- ggplot(counts_can, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = can_neg) +
    geom_col(data = can_pos) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p-value cutoff", y = "Canonical (|y| = count)") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )

  combined_fixed <- p_metrics / p_all_counts_fixed / p_can_counts_fixed +
    plot_layout(heights = c(1, 1, 1))

  combined_auto <- p_metrics / p_all_counts_auto / p_can_counts_auto +
    plot_layout(heights = c(1, 1, 1))

    ggsave(
      file.path(
        base_exp_dir,
        if (eval_mode == "legacy") {
          paste0(tf, "_", method, "_cutoff_metrics_and_counts_by_r_fp.pdf")
        } else {
          paste0(tf, "_", method, "_cutoff_metrics_and_counts_by_r_fp_", eval_mode, ".pdf")
        }
      ),
      combined_auto, width = 12, height = 10, dpi = 600
    )
  }
}
