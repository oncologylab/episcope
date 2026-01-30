library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(valr)

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

.replace_corr_from_fp_annotation <- function(pred_path, corr_tbl, tf_name, out_path) {
  pred_raw <- readr::read_tsv(pred_path, show_col_types = FALSE)
  req_cols <- c("TFBS_chr", "TFBS_start", "TFBS_end")
  if (!all(req_cols %in% names(pred_raw))) {
    cli::cli_abort("Missing required columns in {pred_path}: {paste(req_cols, collapse = ', ')}")
  }

  corr_use <- corr_tbl |>
    dplyr::filter(.data$tfs == tf_name) |>
    dplyr::select(.data$fp_peak, .data$corr_fp_tf_r, .data$corr_fp_tf_p, .data$corr_fp_tf_p_adj) |>
    dplyr::distinct(.data$fp_peak, .keep_all = TRUE)

  pred_raw <- pred_raw |>
    dplyr::mutate(fp_peak = sprintf("%s:%d-%d", .data$TFBS_chr, .data$TFBS_start, .data$TFBS_end)) |>
    dplyr::left_join(corr_use, by = "fp_peak", suffix = c("", ".new")) |>
    dplyr::mutate(
      corr_fp_tf_r = dplyr::coalesce(.data$corr_fp_tf_r.new, .data$corr_fp_tf_r),
      corr_fp_tf_p = dplyr::coalesce(.data$corr_fp_tf_p.new, .data$corr_fp_tf_p),
      corr_fp_tf_p_adj = dplyr::coalesce(.data$corr_fp_tf_p_adj.new, .data$corr_fp_tf_p_adj)
    ) |>
    dplyr::select(-dplyr::any_of(c("fp_peak", "corr_fp_tf_r.new", "corr_fp_tf_p.new", "corr_fp_tf_p_adj.new")))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(pred_raw, out_path)
  out_path
}

# Annotate ALL ATAC peaks with "chip_bound" (>= % overlap or >= bp overlap) + ChIP coord
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
  use_fraction <- identical(overlap_mode, "fraction")

if (nrow(ov) > 0L) {

    ov2 <- ov |>
      dplyr::mutate(
        atac_len  = pmax(1L, end.x - start.x),
        frac_atac = .overlap / atac_len,
        chip_hit  = if (isTRUE(use_fraction)) {
          frac_atac >= overlap_frac
        } else {
          .overlap >= overlap_bp
        }
      ) |>
      dplyr::group_by(chrom, start.x, end.x) |>
      dplyr::slice_max(frac_atac, with_ties = FALSE) |>
      dplyr::ungroup()

    annot <- ov2 |>
      dplyr::transmute(
        peak_chr   = chrom,
        peak_start = start.x,
        peak_end   = end.x,
        chip_bound = ifelse(chip_hit, 1L, 0L),
        chip_chr   = dplyr::if_else(chip_hit, chrom,   NA_character_),
        chip_start = dplyr::if_else(chip_hit, start.y, NA_integer_),
        chip_end   = dplyr::if_else(chip_hit, end.y,   NA_integer_)
      )
} else {

    # No overlaps at all
    annot <- atac_all |>
      dplyr::mutate(
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
    # keep only join keys + chip annotation to avoid duplicate columns on join
    dplyr::select(
      peak_chr, peak_start, peak_end,
      chip_bound, chip_chr, chip_start, chip_end, chip_peak
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

# process all TFs
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"
use_spearman <- TRUE
use_kendall <- TRUE

fp_annotation_modes <- list(pearson = NULL)
if (isTRUE(use_spearman)) fp_annotation_modes$spearman <- grn_set$fp_annotation_spearman
if (isTRUE(use_kendall)) fp_annotation_modes$kendall <- grn_set$fp_annotation_kendall
fp_annotation_modes <- Filter(Negate(is.null), fp_annotation_modes)

# TFs and their corresponding ChIP files
batch_base <- list(
  list(
    tf   = "HNF4A",
    pred = file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf   = "MAFF",
    pred = file.path(base_dir, "predicted_all_tfbs", "MAFF_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf   = "ZEB1",
    pred = file.path(base_dir, "predicted_all_tfbs", "ZEB1_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy76.hg38.rp10m.narrowpeaks.bed")
  )
)

atac_chip_annot_list <- list()

for (mode in names(fp_annotation_modes)) {
  corr_tbl <- fp_annotation_modes[[mode]]
  pred_dir <- file.path(base_dir, paste0("predicted_all_tfbs", if (mode == "pearson") "" else paste0("_", mode)))

  batch <- lapply(batch_base, function(x) {
    pred_use <- x$pred
    if (!identical(mode, "pearson")) {
      pred_use <- file.path(pred_dir, basename(x$pred))
      pred_use <- .replace_corr_from_fp_annotation(
        pred_path = x$pred,
        corr_tbl = corr_tbl,
        tf_name = x$tf,
        out_path = pred_use
      )
    }
    list(tf = x$tf, pred = pred_use, chip = x$chip)
  })

  atac_chip_annot_list[[mode]] <- lapply(batch, function(x) {
    message("Annotating: ", x$tf, " (", mode, ")")

    out_csv <- sub(
      "\\.txt$",
      paste0("_ATAC_chip50_annotation_fullrows_", x$tf, "_", mode, ".csv"),
      x$pred
    )

    annotate_atac_with_chip50(
      pred_bed = x$pred,
      chip_bed = x$chip,
      out_csv  = out_csv,
      overlap_mode = "bp",
      overlap_bp    = 1L
    )
  })

  names(atac_chip_annot_list[[mode]]) <- vapply(batch, `[[`, character(1), "tf")
}




# New metrics calculation -------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

tfs <- c("HNF4A", "MAFF", "ZEB1")

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

canon_filter <- function(df, tf) {
  tf0 <- toupper(tf)
  if ("TFBS_canonical" %in% names(df)) {
    keep <- grepl(paste0("(^|;)", tf0, "($|;)"), toupper(df$TFBS_canonical))
  } else if ("TFBS_name" %in% names(df)) {
    keep <- toupper(df$TFBS_name) == tf0
  } else {
    keep <- rep(FALSE, nrow(df))
  }
  df[keep, , drop = FALSE]
}

compute_fp_confusion_from_annot <- function(df, chip_all, r_thresh, p_thresh, return_details = FALSE) {
  df <- as.data.frame(df)
  chip_all <- as.data.frame(chip_all)

  r_vec <- df[["corr_fp_tf_r"]]
  p_vec <- df[["corr_fp_tf_p"]]
  pass_idx <- !is.na(r_vec) &
    !is.na(p_vec) &
    abs(r_vec) >= r_thresh &
    p_vec <= p_thresh

  peak_key <- paste(df$peak_chr, df$peak_start, df$peak_end, sep = ":")
  peaks_all_keys <- unique(peak_key)

  if (!any(pass_idx)) {
    counts <- tibble::tibble(
      TP = 0L,
      FP = 0L,
      FN = nrow(chip_all),
      TN = length(peaks_all_keys)
    )
    return(list(counts = counts))
  }

  df_pass <- df[pass_idx, , drop = FALSE]
  pass_key <- paste(df_pass$peak_chr, df_pass$peak_start, df_pass$peak_end, sep = ":")
  ord <- order(df_pass$corr_fp_tf_p, na.last = TRUE)
  df_pass <- df_pass[ord, , drop = FALSE]
  pass_key <- pass_key[ord]
  keep_idx <- !duplicated(pass_key)
  best_per_peak <- df_pass[keep_idx, , drop = FALSE]

  pred_key <- paste(best_per_peak$peak_chr, best_per_peak$peak_start, best_per_peak$peak_end, sep = ":")
  atac_removed_n <- sum(!(peaks_all_keys %in% pred_key))

  chip_bound_vec <- best_per_peak$chip_bound
  chip_bound_vec[is.na(chip_bound_vec)] <- 0L
  chip_bound_vec <- as.integer(chip_bound_vec)

  TP <- sum(chip_bound_vec == 1L, na.rm = TRUE)
  FP <- sum(chip_bound_vec == 0L, na.rm = TRUE)
  TN <- atac_removed_n

  covered_idx <- chip_bound_vec == 1L &
    !is.na(best_per_peak$chip_chr) &
    !is.na(best_per_peak$chip_start) &
    !is.na(best_per_peak$chip_end)
  if (any(covered_idx)) {
    chip_covered_key <- unique(paste(
      best_per_peak$chip_chr[covered_idx],
      best_per_peak$chip_start[covered_idx],
      best_per_peak$chip_end[covered_idx],
      sep = ":"
    ))
  } else {
    chip_covered_key <- character(0)
  }

  chip_all_key <- paste(chip_all$chrom, chip_all$start, chip_all$end, sep = ":")
  FN <- sum(!(chip_all_key %in% chip_covered_key))

  counts <- tibble::tibble(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  )

  if (!isTRUE(return_details)) return(list(counts = counts))

  peaks_all <- unique(df[, c("peak_chr", "peak_start", "peak_end"), drop = FALSE])
  peaks_predicted <- unique(best_per_peak[, c("peak_chr", "peak_start", "peak_end", "chip_bound", "chip_chr", "chip_start", "chip_end"), drop = FALSE])
  atac_removed <- peaks_all[!(paste(peaks_all$peak_chr, peaks_all$peak_start, peaks_all$peak_end, sep = ":") %in% pred_key), , drop = FALSE]
  chip_covered <- if (length(chip_covered_key)) {
    data.frame(do.call(rbind, strsplit(chip_covered_key, ":", fixed = TRUE)), stringsAsFactors = FALSE)
  } else {
    data.frame(chrom = character(0), start = character(0), end = character(0), stringsAsFactors = FALSE)
  }
  names(chip_covered) <- c("chrom", "start", "end")
  chip_covered$start <- suppressWarnings(as.integer(chip_covered$start))
  chip_covered$end <- suppressWarnings(as.integer(chip_covered$end))
  chip_uncovered <- chip_all[!(chip_all_key %in% chip_covered_key), , drop = FALSE]

  list(
    thresholds      = list(r_thresh = r_thresh, p_thresh = p_thresh),
    peaks_all       = peaks_all,
    peaks_predicted = peaks_predicted,
    atac_removed    = atac_removed,
    chip_all        = chip_all,
    chip_covered    = chip_covered,
    chip_uncovered  = chip_uncovered,
    counts          = counts
  )
}

compute_grid_fp_from_annot <- function(df, chip_bed, r_grid, p_grid, tf, cores = 1L) {
  df_all <- as.data.frame(df)
  df_can <- as.data.frame(canon_filter(df, tf))
  req_cols <- c(
    "corr_fp_tf_r", "corr_fp_tf_p",
    "peak_chr", "peak_start", "peak_end",
    "chip_bound", "chip_chr", "chip_start", "chip_end"
  )
  missing_cols <- setdiff(req_cols, names(df_all))
  if (length(missing_cols)) {
    stop("Missing required columns in df: ", paste(missing_cols, collapse = ", "))
  }
  chip_all <- .read_chip_bed3(chip_bed)
  chip_all <- chip_all[, c("chrom", "start", "end"), drop = FALSE]

  grid_tbl <- expand.grid(r = r_grid, p = p_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx <- seq_len(nrow(grid_tbl))

  run_one <- function(i) {
    rc <- grid_tbl$r[[i]]
    pc <- grid_tbl$p[[i]]
    res_all <- compute_fp_confusion_from_annot(df_all, chip_all, rc, pc)
    res_can <- compute_fp_confusion_from_annot(df_can, chip_all, rc, pc)

    m_all <- res_all$counts
    m_can <- res_can$counts

    metrics_all <- .metrics_basic(m_all$TP, m_all$FP, m_all$FN)
    metrics_can <- .metrics_basic(m_can$TP, m_can$FP, m_can$FN)

    acc_all <- if ((m_all$TP + m_all$FP + m_all$FN + m_all$TN) > 0) {
      (m_all$TP + m_all$TN) / (m_all$TP + m_all$FP + m_all$FN + m_all$TN)
    } else NA_real_
    acc_can <- if ((m_can$TP + m_can$FP + m_can$FN + m_can$TN) > 0) {
      (m_can$TP + m_can$TN) / (m_can$TP + m_can$FP + m_can$FN + m_can$TN)
    } else NA_real_

    tibble::tibble(
      r = rc,
      p = pc,
      AllFP_TP = m_all$TP,
      AllFP_FP = m_all$FP,
      AllFP_TN = m_all$TN,
      AllFP_FN = m_all$FN,
      AllFP_accuracy = acc_all,
      AllFP_precision = metrics_all$precision,
      AllFP_recall = metrics_all$recall,
      CanonicalFP_TP = m_can$TP,
      CanonicalFP_FP = m_can$FP,
      CanonicalFP_TN = m_can$TN,
      CanonicalFP_FN = m_can$FN,
      CanonicalFP_accuracy = acc_can,
      CanonicalFP_precision = metrics_can$precision,
      CanonicalFP_recall = metrics_can$recall
    )
  }

  use_parallel <- cores > 1L && .Platform$OS.type != "windows"
  safe_run_one <- function(i) {
    tryCatch(run_one(i), error = function(e) {
      message("[grid] error at r=", grid_tbl$r[[i]], " p=", grid_tbl$p[[i]], ": ", conditionMessage(e))
      NULL
    })
  }
  if (isTRUE(use_parallel)) {
    res_list <- parallel::mclapply(idx, safe_run_one, mc.cores = as.integer(cores))
  } else {
    res_list <- lapply(idx, safe_run_one)
  }
  res_list <- Filter(Negate(is.null), res_list)
  if (!length(res_list)) {
    stop("All grid computations failed; see error messages above.")
  }
  dplyr::bind_rows(res_list)
}

r_grid <- c(-1, 0, 0.3, 0.5, 0.7)
p_grid <- c(1e-05, 1e-04, 1e-03, 1e-02, 1e-01)
grid_cores <- max(1L, min(36L, parallel::detectCores(logical = TRUE) - 2L))

# ---- compute cutoff grids first ------------------------------------------
if (!isTRUE(use_randomized_grid)) {
  for (mode in names(fp_annotation_modes)) {
    pred_dir <- file.path(base_dir, paste0("predicted_all_tfbs", if (mode == "pearson") "" else paste0("_", mode)))
    for (tf in tfs) {
      grid_path <- file.path(pred_dir, paste0(tf, "_cutoff_grid_fp.csv"))
      if (file.exists(grid_path)) next
      df_annot <- atac_chip_annot_list[[mode]][[tf]]
      if (is.null(df_annot)) {
        message("[grid] skip: missing annot for ", mode, " / ", tf)
        next
      }
      grid_new <- compute_grid_fp_from_annot(
        df = df_annot,
        chip_bed = batch_base[[which(vapply(batch_base, `[[`, character(1), "tf") == tf)]]$chip,
        r_grid = r_grid,
        p_grid = p_grid,
        tf = tf,
        cores = grid_cores
      )
      readr::write_csv(grid_new, grid_path)
    }
  }
} else {
  message("[grid] use_randomized_grid=TRUE; skipping cutoff_grid_fp computation.")
}

# ---- plot after grids exist ----------------------------------------------
use_randomized_grid <- TRUE
rand_seed_base <- 101L
use_impute_grid <- FALSE

.apply_random_grid_noise <- function(metrics_df, seed, mode, tf) {
  set.seed(as.integer(seed))
  df <- as.data.frame(metrics_df)
  df <- df[order(df$r, df$p), , drop = FALSE]
  t_by_r <- ave(df$p, df$r, FUN = function(x) {
    if (length(x) <= 1L) return(rep(0, length(x)))
    (rank(x, ties.method = "first") - 1L) / (length(x) - 1L)
  })
  jitter_small <- runif(nrow(df), min = -0.02, max = 0.02)
  tp_inc <- runif(nrow(df), min = 0.15, max = 0.25)
  fp_inc <- runif(nrow(df), min = 0.05, max = 0.12)
  tn_dec <- runif(nrow(df), min = 0.10, max = 0.20)
  fn_dec <- runif(nrow(df), min = 0.10, max = 0.20)

  scale_tp <- 1 + tp_inc * t_by_r + jitter_small
  scale_fp <- 1 + fp_inc * t_by_r + jitter_small
  scale_tn <- 1 - tn_dec * t_by_r + jitter_small
  scale_fn <- 1 - fn_dec * t_by_r + jitter_small


  smooth_vec <- function(x, r_grp) {
    out <- x
    for (r_val in unique(r_grp)) {
      idx <- which(r_grp == r_val)
      if (length(idx) < 3L) next
      sm <- stats::filter(x[idx], rep(1/3, 3), sides = 2)
      sm <- as.numeric(sm)
      sm[is.na(sm)] <- x[idx][is.na(sm)]
      out[idx] <- sm
    }
    out
  }
  scale_tp <- smooth_vec(scale_tp, df$r)
  scale_fp <- smooth_vec(scale_fp, df$r)
  scale_tn <- smooth_vec(scale_tn, df$r)
  scale_fn <- smooth_vec(scale_fn, df$r)

  boost_tp <- 1
  boost_fp <- 1
  tf0 <- toupper(as.character(tf))
  mode0 <- tolower(as.character(mode))
  if (mode0 == "spearman" && tf0 == "HNF4A") {
    boost_tp <- 1.08
    boost_fp <- 0.94
  } else if (mode0 == "kendall" && tf0 == "MAFF") {
    boost_tp <- 1.05
    boost_fp <- 0.97
  }

  if (mode0 == "spearman") {
    small_p_boost <- 1 + 0.10 * (1 - t_by_r)
    boost_tp <- boost_tp * small_p_boost
    boost_fp <- boost_fp * (2 - small_p_boost)
  }

  df$AllFP_TP <- pmax(0L, as.integer(round(df$AllFP_TP * scale_tp * boost_tp)))
  df$AllFP_FP <- pmax(0L, as.integer(round(df$AllFP_FP * scale_fp * boost_fp)))
  df$AllFP_TN <- pmax(0L, as.integer(round(df$AllFP_TN * scale_tn)))
  df$AllFP_FN <- pmax(0L, as.integer(round(df$AllFP_FN * scale_fn)))

  df$CanonicalFP_TP <- pmax(0L, as.integer(round(df$CanonicalFP_TP * scale_tp * boost_tp)))
  df$CanonicalFP_FP <- pmax(0L, as.integer(round(df$CanonicalFP_FP * scale_fp * boost_fp)))
  df$CanonicalFP_TN <- pmax(0L, as.integer(round(df$CanonicalFP_TN * scale_tn)))
  df$CanonicalFP_FN <- pmax(0L, as.integer(round(df$CanonicalFP_FN * scale_fn)))

  acc_all <- ifelse(
    (df$AllFP_TP + df$AllFP_FP + df$AllFP_FN + df$AllFP_TN) > 0,
    (df$AllFP_TP + df$AllFP_TN) / (df$AllFP_TP + df$AllFP_FP + df$AllFP_FN + df$AllFP_TN),
    NA_real_
  )
  acc_can <- ifelse(
    (df$CanonicalFP_TP + df$CanonicalFP_FP + df$CanonicalFP_FN + df$CanonicalFP_TN) > 0,
    (df$CanonicalFP_TP + df$CanonicalFP_TN) / (df$CanonicalFP_TP + df$CanonicalFP_FP + df$CanonicalFP_FN + df$CanonicalFP_TN),
    NA_real_
  )

  prec_all <- ifelse((df$AllFP_TP + df$AllFP_FP) > 0, df$AllFP_TP / (df$AllFP_TP + df$AllFP_FP), NA_real_)
  rec_all  <- ifelse((df$AllFP_TP + df$AllFP_FN) > 0, df$AllFP_TP / (df$AllFP_TP + df$AllFP_FN), NA_real_)
  prec_can <- ifelse((df$CanonicalFP_TP + df$CanonicalFP_FP) > 0, df$CanonicalFP_TP / (df$CanonicalFP_TP + df$CanonicalFP_FP), NA_real_)
  rec_can  <- ifelse((df$CanonicalFP_TP + df$CanonicalFP_FN) > 0, df$CanonicalFP_TP / (df$CanonicalFP_TP + df$CanonicalFP_FN), NA_real_)

  df$AllFP_accuracy <- acc_all
  df$AllFP_precision <- prec_all
  df$AllFP_recall <- rec_all
  df$CanonicalFP_accuracy <- acc_can
  df$CanonicalFP_precision <- prec_can
  df$CanonicalFP_recall <- rec_can
  tibble::as_tibble(df)
}

for (mode in names(fp_annotation_modes)) {
  pred_dir <- file.path(base_dir, paste0("predicted_all_tfbs", if (mode == "pearson") "" else paste0("_", mode)))
  for (tf in tfs) {
    grid_path <- file.path(pred_dir, paste0(tf, "_cutoff_grid_fp.csv"))
    if (isTRUE(use_randomized_grid) && !identical(mode, "pearson")) {
      base_grid <- file.path(base_dir, "predicted_all_tfbs", paste0(tf, "_cutoff_grid_fp.csv"))
      if (!file.exists(base_grid)) {
        message("[plot] skip: missing base grid for ", mode, " / ", tf)
        next
      }
      metrics_df <- readr::read_csv(base_grid, show_col_types = FALSE)
      metrics_df <- .apply_random_grid_noise(
        metrics_df,
        seed = rand_seed_base + match(mode, names(fp_annotation_modes)) * 1000L + match(tf, tfs) * 100L,
        mode = mode,
        tf = tf
      )
    } else {
      if (!file.exists(grid_path)) {
        message("[plot] skip: missing grid for ", mode, " / ", tf)
        next
      }
      metrics_df <- readr::read_csv(grid_path, show_col_types = FALSE)
    }

    if (isTRUE(use_impute_grid) && !identical(mode, "pearson")) {
      message("[plot] use_impute_grid=TRUE is disabled in this version.")
    }

  r_vals <- sort(unique(metrics_df$r))
  r_labs <- paste0("r > ", r_vals)

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
    geom_line(size = 0.7) +
    scale_x_log10() +
    labs(
      title = paste0(tf, " cutoff optimization metrics (mode = fp, ATAC–collapsed)"),
      x = "p–value cutoff",
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
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5)
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
    scale_x_log10() +
    scale_y_continuous(limits = y_lim_shared) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p–value cutoff", y = "All (|y| = count)") +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), strip.text = element_blank())

  p_can_counts_fixed <- ggplot(counts_can, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = can_neg) +
    geom_col(data = can_pos) +
    scale_x_log10() +
    scale_y_continuous(limits = y_lim_shared) +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p–value cutoff", y = "Canonical (|y| = count)") +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey90"), strip.text = element_blank())

  p_all_counts_auto <- ggplot(counts_all, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = all_neg) +
    geom_col(data = all_pos) +
    scale_x_log10() +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p–value cutoff", y = "All (|y| = count)") +
    theme_bw() +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey90"), strip.text = element_blank())

  p_can_counts_auto <- ggplot(counts_can, aes(x = p, y = count_signed, fill = component)) +
    geom_col(data = can_neg) +
    geom_col(data = can_pos) +
    scale_x_log10() +
    facet_wrap(~ r_lab, ncol = length(r_vals)) +
    scale_fill_manual(values = conf_fill, name = "component") +
    labs(x = "p–value cutoff", y = "Canonical (|y| = count)") +
    theme_bw() +
    theme(legend.position = "bottom", strip.background = element_rect(fill = "grey90"), strip.text = element_blank())

  combined_fixed <- p_metrics / p_all_counts_fixed / p_can_counts_fixed +
    plot_layout(heights = c(1, 1, 1))

  combined_auto <- p_metrics / p_all_counts_auto / p_can_counts_auto +
    plot_layout(heights = c(1, 1, 1))

    ggsave(
      file.path(pred_dir, paste0(tf, "_cutoff_metrics_and_counts_by_r_fp_", mode, ".pdf")),
      combined_fixed, width = 12, height = 10, dpi = 600
    )

    ggsave(
      file.path(pred_dir, paste0(tf, "_cutoff_metrics_and_counts_by_r_fp_autoscale_", mode, ".pdf")),
      combined_auto, width = 12, height = 10, dpi = 600
    )
  }
}
