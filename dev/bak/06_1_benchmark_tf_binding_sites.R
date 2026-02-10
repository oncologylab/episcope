#' Benchmark predicted TFBS vs. ChIP/CUT&RUN/CUT&Tag peaks (overview-only)
#'
#' This function benchmarks per-row TFBS/peak predictions from an \code{*_overview.txt}
#' table against ChIP/CUT&RUN/CUT&Tag peaks.
#'
#' Modes:
#' - \code{mode = "fp"}: footprint-based overview tables with both TFBS and ATAC peak coords.
#'   * Uses \code{TFBS_chr/TFBS_start/TFBS_end} for footprints.
#'   * Uses \code{peak_chr/peak_start/peak_end} for ATAC peaks.
#'   * **TN/FP/TP/FN are computed per ATAC peak + ChIP peak** as described in the comments.
#' - \code{mode = "atac"}: ATAC-based overview tables with peak-level predictions.
#'   * Uses \code{peak_chr/peak_start/peak_end} and a bound column.
#'   * TN/FP/TP/FN logic is as in your earlier implementation, but without Venn/ROC/PR.
#'
#' For \code{mode = "fp"} ONLY, you must provide \code{r_thresh} and \code{p_thresh}.
#' These are used to call a footprint "passing" based on \code{corr_fp_tf_r} and
#' \code{corr_fp_tf_p}. The TN/FP/TP/FN logic is:
#'
#' \enumerate{
#'   \item For each footprint row, compute:
#'         \code{_bound = 1} if \code{|corr_fp_tf_r| >= r_thresh} AND \code{corr_fp_tf_p <= p_thresh},
#'         else \code{_bound = 0}.
#'   \item Group by ATAC peak (\code{peak_chr/peak_start/peak_end}). Within each ATAC peak:
#'         - Keep only footprints with \code{_bound == 1}.
#'         - Among those, keep the footprint with the **smallest p-value**.
#'         - If an ATAC peak has **no** passing footprint (\code{_bound == 1}), DROP that
#'           ATAC peak from downstream and count it as **TN**.
#'   \item For the remaining ATAC peaks (each now has exactly one "best" footprint):
#'         - Build an ATAC peak bed and intersect with ChIP peaks.
#'         - If an ATAC peak overlaps a ChIP peak by \strong{≥ 50% of the ATAC length},
#'           classify that ATAC peak as **TP**.
#'         - If an ATAC peak does NOT have such ≥50% overlap with any ChIP peak,
#'           classify it as **FP**.
#'   \item For ChIP peaks:
#'         - Any ChIP peak that does not overlap \emph{any} predicted ATAC peak
#'           (no matter how small) is counted as **FN** ("in ChIP seq but not in ATAC peak").
#' }
#'
#' The confusion matrix is then:
#' \itemize{
#'   \item TN = number of ATAC peaks that had no passing footprint and were dropped.
#'   \item TP = number of predicted ATAC peaks with ≥50% overlap to at least one ChIP peak.
#'   \item FP = number of predicted ATAC peaks without ≥50% overlap to any ChIP peak.
#'   \item FN = number of ChIP peaks that do not overlap any predicted ATAC peak.
#' }
#'
#' The function returns counts, derived metrics, and a stacked-bar confusion plot.
benchmark_tfbs_vs_chip <- function(
    pred_bed,
    chip_bed,
    subset = c("all", "bound", "unbound"),
    score_col = NULL,               # kept for backward-compatibility; not used now
    score_higher_is_better = TRUE,  # kept for backward-compatibility; not used now
    mode   = c("fp", "atac"),
    r_thresh = NULL,                # used ONLY for mode = "fp"
    p_thresh = NULL                 # used ONLY for mode = "fp"
) {
  subset <- base::match.arg(subset)
  mode   <- base::match.arg(mode)

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

  # ---- read inputs ----------------------------------------------------------
  chip_peaks <- .read_chip_bed3(chip_bed)
  pred_raw   <- readr::read_tsv(pred_bed)

  # =====================================================================
  # MODE 1: footprint mode ("fp") with updated TN/FP/TP/FN logic
  # =====================================================================
  if (mode == "fp") {
    if (is.null(r_thresh) || is.null(p_thresh)) {
      cli::cli_abort("For mode = 'fp', you must supply both r_thresh and p_thresh.")
    }

    # ---- required columns for fp mode ---------------------------------------
    required_fp_cols <- c(
      "TFBS_chr", "TFBS_start", "TFBS_end",
      "peak_chr", "peak_start", "peak_end",
      "corr_fp_tf_r", "corr_fp_tf_p"
    )
    if (!base::all(required_fp_cols %in% base::names(pred_raw))) {
      missing <- setdiff(required_fp_cols, base::names(pred_raw))
      cli::cli_abort("For mode = 'fp', missing required columns: {paste(missing, collapse = ', ')}")
    }

    # 1) Per-footprint "passing" status based on r_thresh and p_thresh --------
    #    _bound = 1 if |r| >= r_thresh AND p <= p_thresh; else 0
    r_val <- suppressWarnings(base::as.numeric(pred_raw[["corr_fp_tf_r"]]))
    p_val <- suppressWarnings(base::as.numeric(pred_raw[["corr_fp_tf_p"]]))
    pred_raw[["_bound"]] <- base::as.integer(
      !is.na(r_val) & !is.na(p_val) &
        (abs(r_val) >= r_thresh) &
        (p_val <= p_thresh)
    )

    # 2) Collapse footprints per ATAC peak, keeping only the smallest p among
    #    passing footprints; ATAC peaks with NO passing footprint are TN ------
    #    - atac_all: all unique ATAC peaks in the overview.
    #    - fp_pos  : only footprints with _bound == 1.
    #    - fp_best : for each ATAC peak with at least one passing footprint,
    #                keep the footprint with minimal corr_fp_tf_p.
    #    - atac_pos: ATAC peaks that have at least one passing footprint (predicted +).
    atac_all <- pred_raw[, c("peak_chr", "peak_start", "peak_end"), drop = FALSE] |>
      dplyr::distinct()
    n_atac_total <- nrow(atac_all)

    fp_pos <- pred_raw[pred_raw[["_bound"]] == 1L & !is.na(pred_raw[["_bound"]]), , drop = FALSE]

    if (nrow(fp_pos) > 0L) {
      fp_best <- fp_pos |>
        dplyr::group_by(peak_chr, peak_start, peak_end) |>
        dplyr::slice_min(order_by = corr_fp_tf_p, with_ties = FALSE) |>
        dplyr::ungroup()
      atac_pos <- fp_best[, c("peak_chr", "peak_start", "peak_end"), drop = FALSE] |>
        dplyr::distinct()
    } else {
      fp_best  <- pred_raw[0, , drop = FALSE]
      atac_pos <- atac_all[0, , drop = FALSE]
    }

    # ATAC peaks with NO passing footprint are TN and are removed downstream
    atac_removed <- dplyr::anti_join(
      atac_all,
      atac_pos,
      by = c("peak_chr", "peak_start", "peak_end")
    )
    TN <- nrow(atac_removed)  # <-- TN definition (Step 2)

    # 3) Among remaining ATAC peaks (each with one best footprint), classify vs ChIP:
    #    - Build bed for predicted ATAC peaks (atac_pos).
    #    - TP: predicted ATAC peaks with ≥50% overlap with any ChIP peak.
    #    - FP: predicted ATAC peaks without such ≥50% overlap.
    pred_bed_peaks <- tibble::tibble(
      chrom = base::as.character(atac_pos$peak_chr),
      start = base::as.integer(atac_pos$peak_start),
      end   = base::as.integer(atac_pos$peak_end)
    )
    class(pred_bed_peaks) <- c("tbl_df", "tbl", "data.frame", "bed_frame")

    # ATAC peaks that qualify as TP (≥50% of ATAC overlapped by ChIP)
    pred_tp_peaks <- .truth_from_overlap_half_atac(pred_bed_peaks, chip_peaks)

    key_pred <- paste0(atac_pos$peak_chr, "_", atac_pos$peak_start, "_", atac_pos$peak_end)
    key_tp   <- if (nrow(pred_tp_peaks)) {
      paste0(pred_tp_peaks$chrom, "_", pred_tp_peaks$start, "_", pred_tp_peaks$end)
    } else {
      character(0)
    }
    pred_is_tp <- key_pred %in% key_tp

    TP <- sum(pred_is_tp)      # <-- TP definition (Step 3, ≥50% overlap)
    FP <- sum(!pred_is_tp)     # <-- FP definition (Step 3, predicted ATAC but not ≥50% ChIP)

    # 4) FN: ChIP peaks that do not overlap ANY predicted ATAC peak ----------
    #    "if it is in ChIP seq but not in ATAC peak then it is a FN".
    #    Here, "in ATAC peak" = overlaps any predicted ATAC peak (no 50% threshold).
    ov_chip <- if (nrow(pred_bed_peaks)) {
      valr::bed_intersect(chip_peaks, pred_bed_peaks)
    } else {
      chip_peaks[0, , drop = FALSE]
    }

    chip_with_pred <- if (nrow(ov_chip)) {
      tibble::tibble(
        chrom = ov_chip[["chrom"]],
        start = base::as.integer(ov_chip[["start.x"]]),
        end   = base::as.integer(ov_chip[["end.x"]])
      ) |>
        .collapse_unique_sites()
    } else {
      chip_peaks[0, , drop = FALSE]
    }

    key_chip_all <- paste0(chip_peaks$chrom, "_", chip_peaks$start, "_", chip_peaks$end)
    key_chip_cov <- if (nrow(chip_with_pred)) {
      paste0(chip_with_pred$chrom, "_", chip_with_pred$start, "_", chip_with_pred$end)
    } else {
      character(0)
    }

    FN <- sum(!(key_chip_all %in% key_chip_cov))  # <-- FN definition (Step 4)

    # ---- metrics / bar plot -------------------------------------------------
    counts <- tibble::tibble(TP = TP, FP = FP, FN = FN, TN = TN)

    total_for_acc <- TP + FP + FN + TN
    accuracy_overall <- if (total_for_acc > 0) (TP + TN) / total_for_acc else NA_real_

    metrics <- dplyr::bind_cols(
      .metrics_basic(TP, FP, FN),
      tibble::tibble(accuracy = accuracy_overall),
      .metrics_tn(TP, FP, FN, TN)
    )

    # stacked bar: Predicted -, Predicted +
    df_neg <- tibble::tibble(
      bar   = factor("Predicted \u2212", levels = c("Predicted \u2212", "Predicted +")),
      class = factor(c("TN","FN"), levels = c("TN","FN")),
      n     = c(TN, FN)
    )
    df_pos <- tibble::tibble(
      bar   = factor("Predicted +", levels = c("Predicted \u2212", "Predicted +")),
      class = factor(c("FP","TP"), levels = c("FP","TP")),
      n     = c(FP, TP)
    )

    bar_counts <- ggplot2::ggplot() +
      ggplot2::geom_col(
        data = df_neg, ggplot2::aes(x = bar, y = n, fill = class),
        position = "fill", color = "white", linewidth = 0.5
      ) +
      ggplot2::geom_col(
        data = df_pos, ggplot2::aes(x = bar, y = n, fill = class),
        position = "fill", color = "white", linewidth = 0.5
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::scale_fill_manual(
        values = c("TP"="#2ca02c","FP"="#d62728","TN"="#1f77b4","FN"="#ff7f0e"),
        breaks = c("TP","FP","TN","FN"),
        name   = "Class",
        drop   = FALSE
      ) +
      ggplot2::labs(
        title    = "Prediction vs. ChIP (footprint mode, ATAC-collapsed)",
        subtitle = "Counts per predicted class",
        x        = "Predicted class",
        y        = "Proportion"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(face = "bold"),
        axis.title.x = ggplot2::element_text(face = "bold"),
        axis.title.y = ggplot2::element_text(face = "bold")
      )

    # return with useful debugging pieces
    return(list(
      data = list(
        pred_full       = pred_raw,        # original + _bound
        atac_all        = atac_all,        # all ATAC peaks
        atac_removed_TN = atac_removed,    # ATAC peaks counted as TN
        atac_predicted  = atac_pos,        # ATAC peaks with at least one passing footprint
        atac_tp_peaks   = pred_tp_peaks,   # ATAC peaks (subset of predicted) counted as TP
        chip_all        = chip_peaks,      # all ChIP peaks
        chip_fn_peaks   = chip_peaks[!(key_chip_all %in% key_chip_cov), , drop = FALSE]
      ),
      counts  = counts,
      metrics = metrics,
      plots   = list(
        bar_counts = bar_counts
      )
    ))
  }

  # =====================================================================
  # MODE 2: atac mode
  # =====================================================================

  # For atac mode, we use peak_chr/peak_start/peak_end and a bound column.
  coord_chr_col   <- "peak_chr"
  coord_start_col <- "peak_start"
  coord_end_col   <- "peak_end"
  coord_required  <- c("peak_chr", "peak_start", "peak_end")
  bound_candidates <- c("X_bound", "_bound", "bound")

  if (!base::all(coord_required %in% base::names(pred_raw))) {
    cli::cli_abort("For mode = 'atac', `pred_bed` must have columns: {paste(coord_required, collapse = ', ')}.")
  }

  bound_col <- bound_candidates[bound_candidates %in% base::names(pred_raw)][1]
  if (is.na(bound_col)) {
    cli::cli_abort("For mode = 'atac', `pred_bed` must contain a binary prediction column named one of: {paste(bound_candidates, collapse = ', ')}.")
  }

  # normalize prediction (0/1) on FULL table
  pred_bin_all <- base::as.integer(
    base::as.character(pred_raw[[bound_col]]) %in% c("1","TRUE","True","true","T","t",1)
  )

  # build ATAC peak bed
  atac_full_bed <- .as_bed3(pred_raw, coord_chr_col, coord_start_col, coord_end_col)
  atac_full_bed <- .collapse_unique_sites(atac_full_bed)

  # define "truth" for each ATAC peak: ≥50% ATAC overlap with ChIP
  atac_tp_full <- .truth_from_overlap_half_atac(atac_full_bed, chip_peaks)

  key_full <- paste0(pred_raw[[coord_chr_col]], "_",
                     pred_raw[[coord_start_col]], "_",
                     pred_raw[[coord_end_col]])
  key_tp_f <- if (nrow(atac_tp_full)) {
    paste0(atac_tp_full$chrom, "_", atac_tp_full$start, "_", atac_tp_full$end)
  } else {
    character(0)
  }
  pred_raw$ChIP_bound <- base::as.integer(key_full %in% key_tp_f)

  # write back an annotated file (kept from your previous workflow)
  readr::write_csv(pred_raw, sub(".txt", "_with_ChIP.csv", pred_bed))

  # subset rows for per-row evaluation (if requested)
  keep_idx <- switch(
    subset,
    all     = rep(TRUE, base::nrow(pred_raw)),
    bound   = pred_bin_all == 1L,
    unbound = pred_bin_all == 0L
  )
  pred_sub <- pred_raw[keep_idx, , drop = FALSE]
  pred_bin <- pred_bin_all[keep_idx]

  atac_bed_sub <- .as_bed3(pred_sub, coord_chr_col, coord_start_col, coord_end_col)
  atac_bed_sub <- .collapse_unique_sites(atac_bed_sub)

  # truth per row in subset
  key_all <- paste0(pred_sub[[coord_chr_col]], "_",
                    pred_sub[[coord_start_col]], "_",
                    pred_sub[[coord_end_col]])
  key_tp  <- key_tp_f
  truth   <- base::as.integer(key_all %in% key_tp)

  # confusion matrix
  TP <- base::sum(pred_bin == 1L & truth == 1L, na.rm = TRUE)
  FP <- base::sum(pred_bin == 1L & truth == 0L, na.rm = TRUE)
  TN <- base::sum(pred_bin == 0L & truth == 0L, na.rm = TRUE)
  FN <- base::sum(pred_bin == 0L & truth == 1L, na.rm = TRUE)

  counts <- tibble::tibble(TP = TP, FP = FP, FN = FN, TN = TN)

  total_for_acc <- TP + FP + FN + TN
  accuracy_overall <- if (total_for_acc > 0) (TP + TN) / total_for_acc else NA_real_

  metrics <- dplyr::bind_cols(
    .metrics_basic(TP, FP, FN),
    tibble::tibble(accuracy = accuracy_overall),
    .metrics_tn(TP, FP, FN, TN)
  )

  # bar plot (same style as before)
  df_neg <- tibble::tibble(
    bar   = factor("Predicted \u2212", levels = c("Predicted \u2212", "Predicted +")),
    class = factor(c("TN","FN"), levels = c("TN","FN")),
    n     = c(TN, FN)
  )
  df_pos <- tibble::tibble(
    bar   = factor("Predicted +", levels = c("Predicted \u2212", "Predicted +")),
    class = factor(c("FP","TP"), levels = c("FP","TP")),
    n     = c(FP, TP)
  )

  bar_counts <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = df_neg, ggplot2::aes(x = bar, y = n, fill = class),
      position = "fill", color = "white", linewidth = 0.5
    ) +
    ggplot2::geom_col(
      data = df_pos, ggplot2::aes(x = bar, y = n, fill = class),
      position = "fill", color = "white", linewidth = 0.5
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_fill_manual(
      values = c("TP"="#2ca02c","FP"="#d62728","TN"="#1f77b4","FN"="#ff7f0e"),
      breaks = c("TP","FP","TN","FN"),
      name   = "Class",
      drop   = FALSE
    ) +
    ggplot2::labs(
      title    = "Prediction vs. ChIP (ATAC mode)",
      subtitle = "Counts per predicted class",
      x        = "Predicted class",
      y        = "Proportion"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold")
    )

  list(
    data = list(
      pred_full = pred_raw,
      pred      = pred_sub,
      chip      = chip_peaks
    ),
    counts  = counts,
    metrics = metrics,
    plots   = list(
      bar_counts = bar_counts
    )
  )
}

# -- batch runner: save 1 PDF plot per TF -----------------------------------
save_benchmark_plots <- function(tf,
                                 pred_path,
                                 chip_path,
                                 out_dir,
                                 subset = "all",
                                 score_col = NULL,
                                 score_higher_is_better = TRUE,
                                 width = 6,
                                 height = 5,
                                 mode = c("fp", "atac"),
                                 r_thresh = NULL,
                                 p_thresh = NULL) {
  mode <- base::match.arg(mode)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res <- benchmark_tfbs_vs_chip(
    pred_bed   = pred_path,
    chip_bed   = chip_path,
    subset     = subset,
    score_col  = score_col,
    score_higher_is_better = score_higher_is_better,
    mode       = mode,
    r_thresh   = r_thresh,
    p_thresh   = p_thresh
  )

  f_bar <- file.path(out_dir, sprintf("%s_stacked_confusion_%s.pdf", tf, mode))

  if (!is.null(res$plots$bar_counts)) {
    ggplot2::ggsave(f_bar, res$plots$bar_counts, width = width, height = height, device = "pdf")
  }

  invisible(res)
}

# ---------------------------------------------------------------------------
# run fp mode step-by-step for HNF4A
# ---------------------------------------------------------------------------

# base_dir <- "Z:/episcope_test/benchmark_tf_binding_sites_prediction"
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

pred_path_fp <- file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview.txt")
chip_path    <- file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")

# Choose thresholds manually for debugging:
r_thr <- 0.5
p_thr <- 1e-7

# 1) Run benchmark once and capture result
res_fp <- benchmark_tfbs_vs_chip(
  pred_bed   = pred_path_fp,
  chip_bed   = chip_path,
  subset     = "all",
  mode       = "fp",
  r_thresh   = r_thr,
  p_thresh   = p_thr
)

# 2) Inspect TN/TP/FP/FN and intermediate tables
res_fp$counts      # TP, FP, FN, TN
res_fp$metrics     # precision / recall / F1 / accuracy / specificity / MCC

# ATAC peaks that were counted as TN (not called):
# res_fp$data$atac_removed_TN |> head()
res_fp$data$atac_removed_TN
# ATAC peaks that were predicted positive (called):
# res_fp$data$atac_predicted  |> head()
res_fp$data$atac_predicted
# Among called ATAC peaks, those classified as TP (≥50% overlap with ChIP):
# res_fp$data$atac_tp_peaks   |> head()
res_fp$data$atac_tp_peaks
# ChIP peaks that were counted as FN (no overlap with called ATAC peak):
# res_fp$data$chip_fn_peaks   |> head()
res_fp$data$chip_fn_peaks

# 3) the bar plot as a PDF:
out_dir_fp <- file.path(base_dir, "plots_pdf_all_tfbs_fp_debug")
# save_benchmark_plots(
#   tf        = "HNF4A",
#   pred_path = pred_path_hnf4a_fp,
#   chip_path = chip_path_hnf4a,
#   out_dir   = out_dir_fp,
#   mode      = "fp",
#   r_thresh  = r_thr,
#   p_thresh  = p_thr
# )

# -------------------------------------------------------------------
# Cutoff grid via benchmark_tfbs_vs_chip() in mode = "fp"
# -------------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# base_dir <- "Z:/episcope_test/benchmark_tf_binding_sites_prediction"
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

# TFs and their matching ChIP/CUT&Tag files
tf_batch <- list(
  HNF4A = list(
    pred = file.path(base_dir, "predicted_canonical_tfbs", "HNF4A_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  MAFF = list(
    pred = file.path(base_dir, "predicted_canonical_tfbs", "MAFF_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  ZEB1 = list(
    pred = file.path(base_dir, "predicted_canonical_tfbs", "ZEB1_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy76.hg38.rp10m.narrowpeaks.bed")
  )
)

tfs <- names(tf_batch)

# -------------------------------------------------------------------
# Helper: run a grid of (r_thresh, p_thresh) using benchmark_tfbs_vs_chip()
# For each pair, extract: accuracy, precision, recall, F1
# -------------------------------------------------------------------
opt_dir <- file.path(base_dir, "predicted_all_tfbs_cutoff_optimization")

if (!dir.exists(opt_dir)) {dir.create(opt_dir, recursive = TRUE)}

optimize_fp_cutoffs_via_benchmark <- function(pred_path,
                                              chip_path,
                                              r_grid  = seq(0, 1, by = 0.02),
                                              p_grid  = 10^seq(-12, -1, length.out = 20)) {
  res_list <- vector("list", length(r_grid) * length(p_grid))
  k <- 1L

  for (rc in r_grid) {
    for (pc in p_grid) {

      message("  r_thresh = ", rc, ", p_thresh = ", signif(pc, 3))

      res <- benchmark_tfbs_vs_chip(
        pred_bed = pred_path,
        chip_bed = chip_path,
        subset   = "all",
        mode     = "fp",
        r_thresh = rc,
        p_thresh = pc
      )

      # res$metrics has: precision, recall, F1, accuracy, specificity, MCC
      m <- res$metrics[1, ]

      res_list[[k]] <- tibble::tibble(
        r_cut            = rc,
        p_cut            = pc,
        accuracy_overall = m$accuracy,
        precision_ChIP1  = m$precision,
        recall_ChIP1     = m$recall,
        f1_ChIP1         = m$F1
      )
      k <- k + 1L
    }
  }

  dplyr::bind_rows(res_list)
}




for (tf in tfs) {
  message("=== Optimizing cutoffs for TF: ", tf, " (mode = fp) ===")

  pred_path <- tf_batch[[tf]]$pred
  chip_path <- tf_batch[[tf]]$chip

  # run the grid using benchmark_tfbs_vs_chip()
  grid <- optimize_fp_cutoffs_via_benchmark(
    pred_path = pred_path,
    chip_path = chip_path,
    r_grid  = seq(0, 1, by = 0.02),
    p_grid  = 10^seq(-12, -1, length.out = 20)
  )

  # save the raw grid table
  grid_file <- file.path(opt_dir, paste0(tf, "_cutoff_grid_fp.csv"))
  readr::write_csv(grid, grid_file)

  # prepare data for plotting
  grid_plot <- grid[!is.na(grid$f1_ChIP1), , drop = FALSE]
  grid_plot$r_cut <- factor(grid_plot$r_cut)

  plot_df <- tidyr::pivot_longer(
    grid_plot,
    cols = c("accuracy_overall", "precision_ChIP1", "recall_ChIP1"),
    names_to  = "metric",
    values_to = "value"
  )

  # metrics vs p_cut, facetted by r_cut (ONE PDF per TF)
  p_metrics <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = p_cut, y = value, color = metric)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0(tf, " cutoff optimization metrics (mode = fp, ATAC-collapsed)"),
      x     = "p-value cutoff (log10 scale)",
      y     = "Metric value"
    ) +
    ggplot2::facet_wrap(~ r_cut, ncol = 5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(size = 8)
    )

  out_pdf <- file.path(opt_dir, paste0(tf, "_cutoff_metrics_by_r_fp.pdf"))
  ggplot2::ggsave(
    filename = out_pdf,
    plot     = p_metrics,
    width    = 10,
    height   = 18,
    units    = "in",
    dpi      = 300
  )
}



library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(valr)

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

# Annotate ALL ATAC peaks with "chip_bound" (>=50% overlap) + ChIP coord
annotate_atac_with_chip50 <- function(pred_bed,
                                      chip_bed,
                                      out_csv = NULL) {
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

  if (nrow(ov) > 0L) {
    # For each ATAC peak, compute fraction of ATAC overlapped by ChIP
    ov2 <- ov |>
      dplyr::mutate(
        atac_len  = pmax(1L, end.x - start.x),
        inter_len = pmax(0L, pmin(end.x, end.y) - pmax(start.x, start.y)),
        frac_atac = inter_len / atac_len
      ) |>
      # Keep the ChIP peak with max frac_atac for each ATAC peak
      dplyr::group_by(chrom, start.x, end.x) |>
      dplyr::slice_max(frac_atac, with_ties = FALSE) |>
      dplyr::ungroup()

    # Build annotation for overlapped ATAC peaks
    annot <- ov2 |>
      dplyr::transmute(
        peak_chr   = chrom,
        peak_start = start.x,
        peak_end   = end.x,
        chip_bound = ifelse(frac_atac >= 0.5, 1L, 0L),
        chip_chr   = dplyr::if_else(frac_atac >= 0.5, chrom, NA_character_),
        chip_start = dplyr::if_else(frac_atac >= 0.5, start.y, NA_integer_),
        chip_end   = dplyr::if_else(frac_atac >= 0.5, end.y,   NA_integer_)
      )
  } else {
    # No overlaps at all: everything is unbound
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

# usage: process all three TFs
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

# TFs and their corresponding ChIP files
batch <- list(
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

atac_chip_annot_list <- lapply(batch, function(x) {
  message("Annotating: ", x$tf)

  out_csv <- sub(
    "\\.txt$",
    paste0("_ATAC_chip50_annotation_fullrows_", x$tf, ".csv"),
    x$pred
  )

  annotate_atac_with_chip50(
    pred_bed = x$pred,
    chip_bed = x$chip,
    out_csv  = out_csv
  )
})

names(atac_chip_annot_list) <- vapply(batch, `[[`, character(1), "tf")



# Convert Motif name to HGNC
# library(dplyr)
# motif_map <- motif_db %>%
#   select(motif, HGNC)
# atac_chip_annot_list[["ZEB1"]] <- atac_chip_annot_list[["ZEB1"]] %>%
#   left_join(motif_map, by = c("TFBS_name" = "motif")) %>%
#   mutate(TFBS_name = if_else(!is.na(HGNC), HGNC, TFBS_name)) %>%
#   select(-HGNC)

# readr::write_csv(atac_chip_annot_list[["HNF4A"]], file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview_fp_corr_atac_chip50_intersect_annotation.csv"))
# readr::write_csv(atac_chip_annot_list[["MAFF"]], file.path(base_dir, "predicted_all_tfbs", "MAFF_overview_fp_corr_atac_chip50_intersect_annotation.csv"))
# readr::write_csv(atac_chip_annot_list[["ZEB1"]], file.path(base_dir, "predicted_all_tfbs", "ZEB1_overview_fp_corr_atac_chip50_intersect_annotation.csv"))

base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"
readr::read_tsv(file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview_fp_corr_atac_chip50_intersect_annotation.csv"))

library(dplyr)
library(tibble)
library(readr)


compute_fp_confusion_from_annot <- function(df,
                                            chip_bed,
                                            r_thresh,
                                            p_thresh) {

  # ------------------------------------------------------------------
  # STEP 1. Per-footprint "passing" status based on r_thresh / p_thresh
  #         pass = 1 if |r| >= r_thresh AND p <= p_thresh, else 0
  # ------------------------------------------------------------------
  df2 <- df %>%
    mutate(
      pass = as.integer(
        !is.na(corr_fp_tf_r) &
          !is.na(corr_fp_tf_p) &
          abs(corr_fp_tf_r) >= r_thresh &
          corr_fp_tf_p <= p_thresh
      )
    )

  # ------------------------------------------------------------------
  # STEP 2. All unique ATAC peaks (regardless of passing).
  #         This is the universe of ATAC peaks we start from.
  # ------------------------------------------------------------------
  peaks_all <- df2 %>%
    distinct(peak_chr, peak_start, peak_end)

  # ------------------------------------------------------------------
  # STEP 3. Among *passing* footprints, collapse to ONE footprint per
  #         ATAC peak, keeping the smallest p-value.
  #         These ATAC peaks are "predicted positive".
  # ------------------------------------------------------------------
  df_pass <- df2[df2$pass == 1L & !is.na(df2$pass), , drop = FALSE]

  if (nrow(df_pass) > 0L) {

    # One row per ATAC peak: the footprint with minimum corr_fp_tf_p
    best_per_peak <- df_pass %>%
      group_by(peak_chr, peak_start, peak_end) %>%
      slice_min(order_by = corr_fp_tf_p, with_ties = FALSE) %>%
      ungroup()

    # Keep one row per ATAC peak including chip_bound / chip_peak info
    peaks_predicted <- best_per_peak %>%
      distinct(
        peak_chr,
        peak_start,
        peak_end,
        chip_bound,
        chip_chr,
        chip_start,
        chip_end,
        chip_peak
      )

  } else {
    # No passing footprints at all
    best_per_peak   <- df2[0, , drop = FALSE]
    peaks_predicted <- peaks_all[0, , drop = FALSE] %>%
      mutate(
        chip_bound = integer(0),
        chip_chr   = character(0),
        chip_start = integer(0),
        chip_end   = integer(0),
        chip_peak  = character(0)
      )
  }

  # ------------------------------------------------------------------
  # STEP 4. ATAC peaks with NO passing footprint are dropped.
  #         These are counted as True Negatives (TN).
  #
  #   TN = number of ATAC peaks that had no passing footprint and
  #        therefore do not enter downstream analysis.
  # ------------------------------------------------------------------
  atac_removed <- peaks_all %>%
    anti_join(
      peaks_predicted %>%
        select(peak_chr, peak_start, peak_end),
      by = c("peak_chr", "peak_start", "peak_end")
    )

  TN <- nrow(atac_removed)

  # ------------------------------------------------------------------
  # STEP 5. Among *predicted* ATAC peaks (those with at least one
  #         passing footprint), classify using ChIP overlap:
  #
  #   chip_bound == 1  -> ATAC peak has >=50% overlap with some ChIP
  #                       peak → True Positive (TP)
  #   chip_bound == 0  -> ATAC peak has no >=50% overlap with any ChIP
  #                       peak → False Positive (FP)
  # ------------------------------------------------------------------
  chip_bound_vec <- if ("chip_bound" %in% names(peaks_predicted)) {
    coalesce(peaks_predicted$chip_bound, 0L)
  } else {
    rep(0L, nrow(peaks_predicted))
  }

  TP <- sum(chip_bound_vec == 1L, na.rm = TRUE)
  FP <- sum(chip_bound_vec == 0L, na.rm = TRUE)

  # ------------------------------------------------------------------
  # STEP 6. False Negatives (FN) from the ChIP perspective.
  #
  #   - Read ALL ChIP peaks from chip_bed.
  #   - Identify which ChIP peaks are "covered" by at least one
  #     *predicted* ATAC peak with chip_bound == 1.
  #   - Any ChIP peak not covered in this way is counted as FN:
  #
  #       "in ChIP seq but not in any predicted ATAC peak"
  #
  #   Here we use the chip_chr/chip_start/chip_end columns that were
  #   already defined when you built chip_bound (>=50% overlap).
  # ------------------------------------------------------------------
  # chip_raw <- readr::read_tsv(
  #   chip_bed,
  #   col_names      = FALSE,
  #   show_col_types = FALSE
  # )
  chip_raw <- .read_chip_bed3(chip_bed)



  if (ncol(chip_raw) < 3L) {
    stop("ChIP BED must have at least 3 columns.")
  }

  names(chip_raw)[1:3] <- c("chrom", "start", "end")

  chip_raw$chrom <- as.character(chip_raw$chrom)
  chip_raw$start <- as.integer(chip_raw$start)
  chip_raw$end   <-  as.integer(chip_raw$end)


  chip_all <- chip_raw

  # ChIP peaks that WERE matched to at least one predicted ATAC
  chip_covered <- peaks_predicted %>%
    mutate(chip_bound = chip_bound_vec) %>%
    filter(chip_bound == 1L,
           !is.na(chip_chr),
           !is.na(chip_start),
           !is.na(chip_end)) %>%
    distinct(
      chrom = chip_chr,
      start = chip_start,
      end   = chip_end
    )

  # ChIP peaks that are not covered by any predicted ATAC → FN
  chip_uncovered <- chip_all |>
    tibble::as_tibble() |>
    anti_join(chip_covered, by = c("chrom", "start", "end"))

  FN <- nrow(chip_uncovered)

  # ------------------------------------------------------------------
  # STEP 7. Return only what you asked for:
  #         TN, FP, TP, FN, plus key intermediate tables for debugging.
  # ------------------------------------------------------------------
  counts <- tibble(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN
  )

  list(
    thresholds      = list(r_thresh = r_thresh, p_thresh = p_thresh),
    peaks_all       = peaks_all,       # all ATAC peaks (before filtering)
    peaks_predicted = peaks_predicted, # ATAC peaks with passing footprints
    atac_removed    = atac_removed,    # ATAC peaks counted as TN
    chip_all        = chip_all,        # all ChIP peaks
    chip_covered    = chip_covered,    # ChIP peaks overlapped by predicted ATAC
    chip_uncovered  = chip_uncovered,  # ChIP peaks counted as FN
    counts          = counts           # final TN/FP/TP/FN
  )
}

# atac_chip_annot_list[["HNF4A"]] is full-row annotation for HNF4A
# and chip_path is the original ChIP BED file

base_dir  <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"
chip_path <- file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")

# Choose thresholds manually for debugging
r_thr <- 0.5
p_thr <- 1e-7

df_hnf4a <- atac_chip_annot_list[["HNF4A"]]

res_hnf4a <- compute_fp_confusion_from_annot(
  df       = df_hnf4a,
  chip_bed = chip_path,
  r_thresh = r_thr,
  p_thresh = p_thr
)

# Inspect final counts
res_hnf4a$counts

# Inspect each step:
#   all ATAC peaks:
nrow(res_hnf4a$peaks_all)

#   ATAC peaks removed as TN:
nrow(res_hnf4a$atac_removed)
head(res_hnf4a$atac_removed)

#   predicted ATAC peaks (after collapsing to best p per peak):
nrow(res_hnf4a$peaks_predicted)
head(res_hnf4a$peaks_predicted)

#   ChIP peaks and uncovered ones:
nrow(res_hnf4a$chip_all)
nrow(res_hnf4a$chip_uncovered)
head(res_hnf4a$chip_uncovered)



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

for (tf in tfs) {
  metrics_df <- readr::read_csv(
    file.path(base_dir, "predicted_all_tfbs", paste0(tf, "_cutoff_grid_fp.csv"))
  )

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
    file.path(base_dir, "predicted_all_tfbs", paste0(tf, "_cutoff_metrics_and_counts_by_r_fp.pdf")),
    combined_fixed, width = 12, height = 10, dpi = 600
  )

  ggsave(
    file.path(base_dir, "predicted_all_tfbs", paste0(tf, "_cutoff_metrics_and_counts_by_r_fp_autoscale.pdf")),
    combined_auto, width = 12, height = 10, dpi = 600
  )
}


