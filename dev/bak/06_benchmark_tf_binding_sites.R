#' Benchmark predicted TFBS vs. ChIP/CUT&RUN/CUT&Tag peaks (overview-only)
#'
#' This function benchmarks per-row TFBS predictions from an \code{*_overview.txt}
#' table against ChIP/CUT&RUN/CUT&Tag peaks.
#'
#' Input is REQUIRED to be an “overview” table containing:
#' \itemize{
#'   \item \code{TFBS_chr}, \code{TFBS_start}, \code{TFBS_end}
#'   \item a binary prediction column named \code{_bound} or \code{bound}
#' }
#'
#' Each overview row is treated as one candidate TFBS. "Truth" is defined by
#' overlap of the row's TFBS interval with any ChIP/CUT&RUN/CUT&Tag peak.
#' The confusion matrix (TP/FP/TN/FN) is computed strictly over the (optionally
#' subsetted) overview rows. No genome background is used in this mode.
#'
#' If \code{score_col} is provided (e.g., \code{"corr_fp_tf_r"}), precision–recall
#' and ROC curves plus AUCs are computed by thresholding that column across the
#' (subsetted) overview rows.
#'
#' Additionally, the function produces a **3-way Euler diagram using TFBS row counts**
#' (not base pairs) for the FULL overview table: \emph{Predicted bound} vs
#' \emph{ChIP-bound (overview-derived)} vs \emph{Predicted unbound}. The plot is
#' returned as \code{plots$venn3}, and the underlying counts as \code{data$venn3_counts}.
benchmark_tfbs_vs_chip <- function(
    pred_bed,
    chip_bed,
    subset = c("all", "bound", "unbound"),
    score_col = NULL,
    score_higher_is_better = TRUE,
    mode = c("fp", "atac")  # BEGIN EDIT: new arg
) {
  subset <- base::match.arg(subset)
  mode   <- base::match.arg(mode)  # END EDIT

  # ---- helpers --------------------------------------------------------------
  .stop_if_missing <- function(p) if (!base::file.exists(p)) cli::cli_abort("File not found: {p}")

  .read_tsv_header <- function(path) {
    .stop_if_missing(path)
    df <- try(
      readr::read_tsv(
        path,
        comment   = "#",
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
    if (base::any(!base::is.finite(df$start) | !base::is.finite(df$end))) cli::cli_abort("Non-integer ChIP start/end in: {path}")
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
    f1 <- if (!base::is.na(precision) && !base::is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
    tibble::tibble(precision = precision, recall = recall, F1 = f1)
  }
  .metrics_tn <- function(tp, fp, fn, tn) {
    tp <- base::as.numeric(tp); fp <- base::as.numeric(fp); fn <- base::as.numeric(fn); tn <- base::as.numeric(tn)
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    denom_mcc <- (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
    num_mcc   <- (tp * tn) - (fp * fn)
    mcc <- if (base::is.finite(denom_mcc) && denom_mcc > 0) num_mcc / base::sqrt(denom_mcc) else NA_real_
    tibble::tibble(specificity = specificity, MCC = mcc)
  }
  .auc_trapz <- function(x, y) {
    if (base::length(x) < 2L) return(NA_real_)
    dx <- base::diff(x); my <- (y[-1] + y[-base::length(y)]) / 2
    base::sum(dx * my, na.rm = TRUE)
  }

  .coords_left <- function(df) {
    nm <- base::names(df)
    if (base::all(c("chrom", "start.x", "end.x") %in% nm)) {
      out <- df[, c("chrom", "start.x", "end.x")]; base::names(out) <- c("chrom", "start", "end"); return(out)
    }
    if (base::all(c("chrom.x", "start.x", "end.x") %in% nm)) {
      out <- df[, c("chrom.x", "start.x", "end.x")]; base::names(out) <- c("chrom", "start", "end"); return(out)
    }
    if (base::all(c("chrom", "start", "end", "chrom2", "start2", "end2") %in% nm)) return(df[, c("chrom", "start", "end")])
    if (base::all(c("chrom1", "start1", "end1") %in% nm)) {
      out <- df[, c("chrom1", "start1", "end1")]; base::names(out) <- c("chrom", "start", "end"); return(out)
    }
    cli::cli_abort("Unexpected columns from bed_intersect(): {paste(nm, collapse = ', ')}")
  }

  .label_for_subset <- function(x) switch(x, all = "All rows", bound = "Subset: bound=1", unbound = "Subset: bound=0")

  # BEGIN EDIT: helper to get "truth" intervals with optional ≥50% overlap rule
  .truth_from_overlap <- function(bed_x, chip_peaks, mode) {
    ov <- valr::bed_intersect(bed_x, chip_peaks)
    if (!base::nrow(ov)) {
      return(bed_x[0, , drop = FALSE])
    }
    if (mode == "atac") {
      # fraction of ATAC peak overlapped by ChIP (using x as ATAC)
      start_x <- base::as.integer(ov[["start.x"]])
      end_x   <- base::as.integer(ov[["end.x"]])
      start_y <- base::as.integer(ov[["start.y"]])
      end_y   <- base::as.integer(ov[["end.y"]])
      len_x   <- base::pmax(1L, end_x - start_x)
      inter_w <- base::pmax(0L, base::pmin(end_x, end_y) - base::pmax(start_x, start_y))
      frac    <- inter_w / len_x
      ov <- ov[frac >= 0.5, , drop = FALSE]
      if (!base::nrow(ov)) {
        return(bed_x[0, , drop = FALSE])
      }
    }
    .coords_left(ov) |> .collapse_unique_sites()
  }
  # END EDIT

  # ---- read inputs ----------------------------------------------------------
  chip_peaks <- .read_chip_bed3(chip_bed)
  pred_raw   <- .read_tsv_header(pred_bed)

  # BEGIN EDIT: choose coordinate and bound columns by mode
  if (mode == "fp") {
    coord_chr_col   <- "TFBS_chr"
    coord_start_col <- "TFBS_start"
    coord_end_col   <- "TFBS_end"
    coord_required  <- c("TFBS_chr", "TFBS_start", "TFBS_end")
    bound_candidates <- c("_bound", "bound")
  } else {
    coord_chr_col   <- "peak_chr"
    coord_start_col <- "peak_start"
    coord_end_col   <- "peak_end"
    coord_required  <- c("peak_chr", "peak_start", "peak_end")
    bound_candidates <- c("X_bound", "_bound", "bound")
  }

  if (!base::all(coord_required %in% base::names(pred_raw))) {
    cli::cli_abort("`pred_bed` must contain coordinate columns: {paste(coord_required, collapse = ', ')}.")
  }

  bound_col <- bound_candidates[bound_candidates %in% base::names(pred_raw)][1]
  if (is.na(bound_col)) {
    cli::cli_abort("`pred_bed` must contain a binary prediction column named one of: {paste(bound_candidates, collapse = ', ')}.")
  }
  # END EDIT

  # ---- normalize prediction (0/1) on FULL table -----------------------------
  pred_bin_all <- base::as.integer(base::as.character(pred_raw[[bound_col]]) %in% c("1","TRUE","True","true","T","t",1))

  # ---- compute ChIP_bound on FULL table -------------------------------------
  tfbs_full_bed <- .as_bed3(pred_raw, coord_chr_col, coord_start_col, coord_end_col)
  tfbs_full_bed <- .collapse_unique_sites(tfbs_full_bed)

  # BEGIN EDIT: use helper with optional ≥50% rule
  tfbs_tp_full <- .truth_from_overlap(tfbs_full_bed, chip_peaks, mode)
  # END EDIT

  key_full <- paste0(pred_raw[[coord_chr_col]], "_",
                     pred_raw[[coord_start_col]], "_",
                     pred_raw[[coord_end_col]])
  key_tp_f <- if (base::nrow(tfbs_tp_full)) paste0(tfbs_tp_full$chrom, "_", tfbs_tp_full$start, "_", tfbs_tp_full$end) else character(0)
  pred_raw$ChIP_bound <- base::as.integer(key_full %in% key_tp_f)

  readr::write_csv(pred_raw, sub(".txt", "_with_ChIP.csv", pred_bed))

  # ---- subset rows for per-row evaluation -----------------------------------
  keep_idx <- switch(subset,
                     all     = rep(TRUE, base::nrow(pred_raw)),
                     bound   = pred_bin_all == 1L,
                     unbound = pred_bin_all == 0L
  )
  pred_sub <- pred_raw[keep_idx, , drop = FALSE]
  pred_bin <- pred_bin_all[keep_idx]

  # build TFBS/peak bed from subset
  tfbs_bed <- .as_bed3(pred_sub, coord_chr_col, coord_start_col, coord_end_col)
  tfbs_bed <- .collapse_unique_sites(tfbs_bed)

  # truth by overlap (subset)
  # BEGIN EDIT: use same helper for subset
  tfbs_tp <- .truth_from_overlap(tfbs_bed, chip_peaks, mode)
  # END EDIT

  # map truth back to each ROW in subset
  key_all <- paste0(pred_sub[[coord_chr_col]], "_",
                    pred_sub[[coord_start_col]], "_",
                    pred_sub[[coord_end_col]])
  key_tp  <- if (base::nrow(tfbs_tp)) paste0(tfbs_tp$chrom, "_", tfbs_tp$start, "_", tfbs_tp$end) else character(0)
  truth   <- base::as.integer(key_all %in% key_tp)

  # confusion matrix over subset
  TP <- base::sum(pred_bin == 1L & truth == 1L, na.rm = TRUE)
  FP <- base::sum(pred_bin == 1L & truth == 0L, na.rm = TRUE)
  TN <- base::sum(pred_bin == 0L & truth == 0L, na.rm = TRUE)
  FN <- base::sum(pred_bin == 0L & truth == 1L, na.rm = TRUE)

  # unique coords by class
  pred_coords <- tibble::tibble(
    chrom = pred_sub[[coord_chr_col]],
    start = base::as.integer(pred_sub[[coord_start_col]]),
    end   = base::as.integer(pred_sub[[coord_end_col]]),
    pred  = pred_bin,
    truth = truth
  )
  tp_sites <- .collapse_unique_sites(pred_coords[pred_coords$pred == 1L & pred_coords$truth == 1L, c("chrom", "start", "end"), drop = FALSE])
  fp_sites <- .collapse_unique_sites(pred_coords[pred_coords$pred == 1L & pred_coords$truth == 0L, c("chrom", "start", "end"), drop = FALSE])
  fn_sites <- .collapse_unique_sites(pred_coords[pred_coords$pred == 0L & pred_coords$truth == 1L, c("chrom", "start", "end"), drop = FALSE])

  # Jaccard over predicted-positive TFBS/peaks in subset
  pred_pos_bed <- .collapse_unique_sites(.as_bed3(pred_coords[pred_coords$pred == 1L, ], "chrom", "start", "end"))
  jac <- valr::bed_jaccard(pred_pos_bed, chip_peaks)
  jaccard_bp <- if (base::nrow(jac)) jac$jaccard[1] else NA_real_

  # optional PR/ROC
  pr_curve <- roc_curve <- NULL
  pr_plot  <- roc_plot  <- NULL
  pr_auc   <- NA_real_
  roc_auc  <- NA_real_

  if (!is.null(score_col) && score_col %in% base::names(pred_sub)) {
    sc <- pred_sub[[score_col]]
    if (!base::is.numeric(sc)) {
      sc2 <- suppressWarnings(base::as.numeric(sc))
      sc  <- if (base::all(base::is.finite(sc2))) sc2 else NULL
    }
    if (!is.null(sc)) {
      y_true <- truth
      s      <- if (isTRUE(score_higher_is_better)) sc else -sc
      ord    <- base::order(s, decreasing = TRUE, na.last = NA)
      y_ord  <- y_true[ord]
      tp_cum <- base::cumsum(y_ord == 1L)
      fp_cum <- base::cumsum(y_ord == 0L)
      tot_pos <- base::sum(y_true == 1L)
      tot_neg <- base::sum(y_true == 0L)
      recall_vec <- if (tot_pos > 0) tp_cum / tot_pos else rep(NA_real_, base::length(tp_cum))
      prec_vec   <- ifelse((tp_cum + fp_cum) > 0, tp_cum / (tp_cum + fp_cum), NA_real_)
      fpr_vec    <- if (tot_neg > 0) fp_cum / tot_neg else rep(NA_real_, base::length(fp_cum))

      pr_curve <- tibble::tibble(recall = recall_vec, precision = prec_vec)
      roc_curve <- tibble::tibble(fpr = fpr_vec, tpr = recall_vec)

      pr_auc  <- .auc_trapz(pr_curve$recall[base::order(pr_curve$recall)], pr_curve$precision[base::order(pr_curve$recall)])
      roc_auc <- .auc_trapz(roc_curve$fpr[base::order(roc_curve$fpr)], roc_curve$tpr[base::order(roc_curve$tpr)])

      # --- inside the PR section where pr_plot is created ---
      pr_plot <- ggplot2::ggplot(pr_curve, ggplot2::aes(x = recall, y = precision)) +
        ggplot2::geom_line() +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        # BEGIN EDIT: clearer axis formulas
        ggplot2::labs(
          title    = "Precision–Recall curve",
          subtitle = sprintf("AUC = %.3f", pr_auc),
          x        = "Recall = TP / (TP + FN)",
          y        = "Precision = TP / (TP + FP)"
        ) +
        # END EDIT
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title  = ggplot2::element_text(face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold")
        )

      if (base::nrow(roc_curve) >= 2L) {
        # --- inside the ROC section where roc_plot is created ---
        roc_plot <- ggplot2::ggplot(roc_curve, ggplot2::aes(x = fpr, y = tpr)) +
          ggplot2::geom_line() +
          ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
          ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
          # BEGIN EDIT: clearer axis formulas
          ggplot2::labs(
            title    = "ROC curve",
            subtitle = sprintf("AUC = %.3f", roc_auc),
            x        = "False Positive Rate (FPR) = FP / (FP + TN)",
            y        = "True Positive Rate (TPR) = TP / (TP + FN)"
          ) +
          # END EDIT
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title  = ggplot2::element_text(face = "bold"),
            axis.title.x = ggplot2::element_text(face = "bold"),
            axis.title.y = ggplot2::element_text(face = "bold")
          )
      }
    }
  }

  counts <- tibble::tibble(TP = TP, FP = FP, FN = FN, TN = TN)
  metrics <- dplyr::bind_cols(
    .metrics_basic(TP, FP, FN),
    tibble::tibble(jaccard_bp = jaccard_bp),
    .metrics_tn(TP, FP, FN, TN)
  )

  # ---- bar plot -------------------------------------------------------------
  # ---- stacked bars (proportions), with per-bar stack order -------------------
  df_neg <- tibble::tibble(
    bar   = factor("Predicted \u2212", levels = c("Predicted \u2212", "Predicted +")),
    class = factor(c("TN","FN"), levels = c("TN","FN")),   # TN bottom, FN top
    n     = c(TN, FN)
  )

  df_pos <- tibble::tibble(
    bar   = factor("Predicted +", levels = c("Predicted \u2212", "Predicted +")),
    class = factor(c("FP","TP"), levels = c("FP","TP")),   # FP bottom, TP top
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
    # NOTE: no y-limits here (we'll clamp with coord_cartesian so counts aren't dropped)
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_fill_manual(
      values = c("TP"="#2ca02c","FP"="#d62728","TN"="#1f77b4","FN"="#ff7f0e"),
      breaks = c("TP","FP","TN","FN"),
      name = "Class",
      drop = FALSE
    ) +
    ggplot2::labs(
      title    = "Prediction vs. ChIP (stacked by predicted class)",
      subtitle = paste0("proportions within each bar"),
      x        = "Predicted class", y = "Proportion"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold")
    )
  # ---- 3-way Chow–Ruskey Venn (TFBS row counts) on FULL overview ------------
  venn3_plot   <- NULL
  venn3_counts <- NULL
  if (!requireNamespace("Vennerable", quietly = TRUE)) {
    cli::cli_warn("Skipping Chow–Ruskey Venn: please install {Vennerable}.")
  } else {
    lab_A <- "Pred. bound"
    lab_B <- "ChIP (overview)"
    lab_C <- "Pred. unbound"

    A  <- base::sum(pred_bin_all == 1L)
    B  <- base::sum(pred_raw$ChIP_bound == 1L)
    C  <- base::sum(pred_bin_all == 0L)
    AB <- base::sum(pred_bin_all == 1L & pred_raw$ChIP_bound == 1L)  # TP
    AC <- 0L
    BC <- base::sum(pred_bin_all == 0L & pred_raw$ChIP_bound == 1L)  # FN
    ABC <- 0L

    venn3_counts <- c(
      stats::setNames(A,  lab_A),
      stats::setNames(B,  lab_B),
      stats::setNames(C,  lab_C),
      stats::setNames(AB, paste0(lab_A, "&", lab_B)),
      stats::setNames(AC, paste0(lab_A, "&", lab_C)),
      stats::setNames(BC, paste0(lab_B, "&", lab_C)),
      stats::setNames(ABC, paste0(lab_A, "&", lab_B, "&", lab_C))
    )

    # Convert to region weights expected by Vennerable
    w111 <- ABC
    w110 <- AB - ABC
    w101 <- AC - ABC
    w011 <- BC - ABC
    w100 <- A  - AB - AC + ABC
    w010 <- B  - AB - BC + ABC
    w001 <- C  - AC - BC + ABC

    regs <- c(`100`=w100, `010`=w010, `001`=w001,
              `110`=w110, `101`=w101, `011`=w011, `111`=w111)
    regs[regs < 0] <- 0L

    # ---- important: epsilon for zero regions to avoid geometry crash ----
    eps <- .Machine$double.eps
    regs_plot <- regs
    regs_plot[regs_plot == 0] <- eps

    # make sure S4 methods are available
    if (!("methods" %in% base::loadedNamespaces())) base::requireNamespace("methods")

    venn3_plot <- Vennerable::Venn(
      SetNames = c(lab_A, lab_B, lab_C),
      Weight   = regs_plot
    )
    # keep the exact integers for later reference
    attr(venn3_plot, "weights_exact_integers") <- regs
  }

  list(
    data = list(
      pred_full = pred_raw,  # includes ChIP_bound
      pred      = pred_sub,
      chip      = chip_peaks,
      tp        = tp_sites,
      fp        = fp_sites,
      fn        = fn_sites,
      subset    = subset,
      venn3_counts = venn3_counts
    ),
    counts  = counts,
    metrics = dplyr::mutate(
      metrics,
      pr_auc  = ifelse(base::is.na(pr_auc), NA_real_, pr_auc),
      roc_auc = ifelse(base::is.na(roc_auc), NA_real_, roc_auc)
    ),
    pr_curve = pr_curve,
    roc_curve = roc_curve,
    plots = list(
      bar_counts = bar_counts,
      pr         = pr_plot,
      roc        = roc_plot,
      venn3      = venn3_plot
    )
  )
}

# -- batch runner: save 4 PDF plots per TF -----------------------------------
save_benchmark_plots <- function(tf, pred_path, chip_path, out_dir,
                                 subset = "all",
                                 score_col = "corr_fp_tf_r",
                                 score_higher_is_better = TRUE,
                                 width = 6, height = 5,
                                 mode = c("fp", "atac")  # BEGIN EDIT: new arg
) {
  mode <- base::match.arg(mode)  # END EDIT
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res <- benchmark_tfbs_vs_chip(
    pred_bed   = pred_path,
    chip_bed   = chip_path,
    subset     = subset,
    score_col  = score_col,
    score_higher_is_better = score_higher_is_better,
    mode       = mode    # BEGIN EDIT: pass through
  )

  # filenames (PDF, include TF)
  f_bar  <- file.path(out_dir, sprintf("%s_stacked_confusion.pdf", tf))
  f_pr   <- file.path(out_dir, sprintf("%s_pr_curve.pdf", tf))
  f_roc  <- file.path(out_dir, sprintf("%s_roc_curve.pdf", tf))
  f_venn <- file.path(out_dir, sprintf("%s_venn_chowruskey.pdf", tf))

  # 1) stacked confusion
  if (!is.null(res$plots$bar_counts)) {
    ggplot2::ggsave(f_bar, res$plots$bar_counts, width = width, height = height, device = "pdf")
  }

  # 2) PR curve
  if (!is.null(res$plots$pr)) {
    ggplot2::ggsave(f_pr, res$plots$pr, width = width, height = height, device = "pdf")
  }

  # 3) ROC curve
  if (!is.null(res$plots$roc)) {
    ggplot2::ggsave(f_roc, res$plots$roc, width = width, height = height, device = "pdf")
  }

  # 4) Chow–Ruskey Venn (already PDF-capable via out_file)
  if (!is.null(res$data$venn3_counts)) {
    plot_venn_chowruskey_counts(res$data$venn3_counts, out_file = f_venn, width = width, height = height)
  }

  invisible(res)
}

library(Vennerable)
# Chow–Ruskey only, using plot(Venn, ...)
plot_venn_chowruskey_counts <- function(counts_named,
                                        out_file = NULL, width = 6, height = 5, dpi = 300) {
  if (!requireNamespace("Vennerable", quietly = TRUE))
    cli::cli_abort("Please install {Vennerable} (e.g., install.packages('Vennerable')).")
  if (!("methods" %in% loadedNamespaces())) requireNamespace("methods")

  nm <- names(counts_named)
  stopifnot(length(counts_named) == 7L, is.numeric(counts_named), !is.null(nm))

  # ---- infer base labels and convert to region weights ----
  .bases <- function(nms) {
    uniq <- unique(unlist(strsplit(nms, "&", fixed = TRUE)))
    uniq[nchar(uniq) > 0]
  }

  if (all(sort(nm) == sort(c("100","010","001","110","101","011","111")))) {
    regs <- counts_named
    base_labs <- attr(counts_named, "base_labels")
    if (is.null(base_labs)) base_labs <- c("Set A","Set B","Set C")
  } else {
    bases <- .bases(grep("^[^&]+$", nm, value = TRUE))
    stopifnot(length(bases) == 3L)
    A <- bases[1]; B <- bases[2]; C <- bases[3]
    G <- function(k) if (k %in% nm) as.numeric(counts_named[[k]]) else 0
    A1  <- G(A); B1  <- G(B); C1  <- G(C)
    AB  <- G(paste0(A,"&",B))
    AC  <- G(paste0(A,"&",C))
    BC  <- G(paste0(B,"&",C))
    ABC <- G(paste0(A,"&",B,"&",C))
    regs <- c(`100`=A1-AB-AC+ABC, `010`=B1-AB-BC+ABC, `001`=C1-AC-BC+ABC,
              `110`=AB-ABC,       `101`=AC-ABC,       `011`=BC-ABC,       `111`=ABC)
    regs[regs < 0] <- 0
    base_labs <- c(A, B, C)
  }

  # ---- adaptive epsilon for zero regions (avoid geometry crashes) ----
  nz <- regs[regs > 0]
  eps <- if (length(nz)) max(1, min(nz) * 1e-3) else 1
  regs_plot <- regs
  regs_plot[regs_plot == 0] <- eps

  # ---- build Venn via the same API your example uses ----
  venn <- Vennerable::Venn(n = 3)
  # rename the base set columns to your labels
  colnames(venn@IndicatorWeight)[1:3] <- base_labs

  # Weights must include the 8th region "000" (set to 0)
  w <- rep(0, 8)
  names(w) <- c("000","100","010","001","110","101","011","111")
  w[names(regs_plot)] <- regs_plot
  Vennerable::Weights(venn) <- w

  # ---- open device if requested ----
  open_dev <- function(path) {
    if (is.null(path)) return(NULL)
    ext <- tolower(tools::file_ext(path))
    if (ext == "pdf") grDevices::pdf(path, width = width, height = height)
    else if (ext == "png") grDevices::png(path, width = width*dpi, height = height*dpi, res = dpi)
    else if (ext == "svg") grDevices::svg(path, width = width, height = height)
    else cli::cli_abort("Unsupported out_file extension: {ext}. Use pdf/png/svg.")
    function() grDevices::dev.off()
  }
  closer <- open_dev(out_file); on.exit(if (!is.null(closer)) closer(), add = TRUE)

  # ---- plot Chow–Ruskey with weights ----
  # (No plotVenn: use the exported S4 plot() method)
  plot(venn, type = "ChowRuskey", doWeights = TRUE, show = list(Quantities = TRUE, SetLabels = TRUE))

  invisible(venn)
}


# -- batch runner: save 4 PDF plots per TF -----------------------------------
# save_benchmark_plots <- function(tf, pred_path, chip_path, out_dir,
#                                  subset = "all",
#                                  score_col = "corr_fp_tf_r",
#                                  score_higher_is_better = TRUE,
#                                  width = 6, height = 5) {
#   if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#
#   res <- benchmark_tfbs_vs_chip(
#     pred_bed   = pred_path,
#     chip_bed   = chip_path,
#     subset     = subset,
#     score_col  = score_col,
#     score_higher_is_better = score_higher_is_better
#   )
#
#   # filenames (PDF, include TF)
#   f_bar  <- file.path(out_dir, sprintf("%s_stacked_confusion.pdf", tf))
#   f_pr   <- file.path(out_dir, sprintf("%s_pr_curve.pdf", tf))
#   f_roc  <- file.path(out_dir, sprintf("%s_roc_curve.pdf", tf))
#   f_venn <- file.path(out_dir, sprintf("%s_venn_chowruskey.pdf", tf))
#
#   # 1) stacked confusion
#   if (!is.null(res$plots$bar_counts)) {
#     ggplot2::ggsave(f_bar, res$plots$bar_counts, width = width, height = height, device = "pdf")
#   }
#
#   # 2) PR curve
#   if (!is.null(res$plots$pr)) {
#     ggplot2::ggsave(f_pr, res$plots$pr, width = width, height = height, device = "pdf")
#   }
#
#   # 3) ROC curve
#   if (!is.null(res$plots$roc)) {
#     ggplot2::ggsave(f_roc, res$plots$roc, width = width, height = height, device = "pdf")
#   }
#
#   # 4) Chow–Ruskey Venn (already PDF-capable via out_file)
#   if (!is.null(res$data$venn3_counts)) {
#     plot_venn_chowruskey_counts(res$data$venn3_counts, out_file = f_venn, width = width, height = height)
#   }
#
#   invisible(res)
# }

# ---------------------------------------------------------------------------
# three TFs

# base_dir <- "Z:/episcope_test/benchmark_tf_binding_sites_prediction"
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

batch <- list(
  list(tf = "HNF4A",
       pred = file.path(base_dir, "predicted_canonical_tfbs", "HNF4A_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy83.hg38.rp10m.narrowpeaks.bed")),
  list(tf = "MAFF",
       pred = file.path(base_dir, "predicted_canonical_tfbs", "MAFF_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy84.hg38.rp10m.narrowpeaks.bed")),
  list(tf = "ZEB1",
       pred = file.path(base_dir, "predicted_canonical_tfbs", "ZEB1_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy76.hg38.rp10m.narrowpeaks.bed"))
)
out_dir <- file.path(base_dir, "plots_pdf_canonical_tfbs")
purrr::walk(batch, function(x) {
  save_benchmark_plots(tf = x$tf, pred_path = x$pred, chip_path = x$chip, out_dir = out_dir)
})


# base_dir <- "Z:/episcope_test/benchmark_tf_binding_sites_prediction"
base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

batch <- list(
  list(tf = "HNF4A",
       pred = file.path(base_dir, "predicted_all_tfbs", "HNF4A_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy83.hg38.rp10m.narrowpeaks.bed")),
  list(tf = "MAFF",
       pred = file.path(base_dir, "predicted_all_tfbs", "MAFF_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy84.hg38.rp10m.narrowpeaks.bed")),
  list(tf = "ZEB1",
       pred = file.path(base_dir, "predicted_all_tfbs", "ZEB1_overview.txt"),
       chip = file.path(base_dir, "cutntag",  "cy76.hg38.rp10m.narrowpeaks.bed"))
)

out_dir <- file.path(base_dir, "plots_pdf_all_tfbs")
purrr::walk(batch, function(x) {
  save_benchmark_plots(tf = x$tf, pred_path = x$pred, chip_path = x$chip, out_dir = out_dir)
})


library(readr)
library(ggplot2)

source("R/utils_ggvolcano.R")

tfs <- c("HNF4A", "MAFF", "ZEB1")

for (tf in tfs) {

  message("Processing: ", tf)

  # ---- Load data ----
  infile <- file.path(base_dir, "predicted_canonical_tfbs", paste0(tf, "_overview_with_ChIP.csv"))
  volcano_data <- readr::read_csv(infile)

  # ---- Make base volcano ----
  p <- ggvolcano(
    volcano_data,
    logFC_col  = "corr_fp_tf_r",
    pval_col   = "corr_fp_tf_p",
    use_significance = FALSE,
    highlight_col    = "ChIP_bound",
    highlight_colors = c("0" = "grey80", "1" = "red"),
    highlight_labels = c("Not ChIP-bound", "ChIP-bound"),
    xlab     = bquote(~"Correlation"~italic(r)),
    x_limits = c(-1, 1),
    sample_frac        = 0.2,
    sample_keep_levels = "1",
    point_aes = list(
      size  = 0.25,
      shape = c(16, 16, 16, 16),
      color = c("grey80", "red"),
      alpha = 0.6
    ),
    title   = paste0(tf, " ChIP-bound TFBS"),
    caption = paste("Total variables:", nrow(volcano_data)),
    legend_aes = list(
      labels     = c("Not ChIP-bound", "ChIP-bound"),
      position   = "bottom",
      label_size = 10,
      icon_size  = 4
    )
  )

  # ---- Make jittered volcano ----
  p_jitter <- ggvolcano(
    volcano_data,
    logFC_col  = "corr_fp_tf_r",
    pval_col   = "corr_fp_tf_p",
    use_significance = FALSE,
    highlight_col    = "ChIP_bound",
    highlight_colors = c("0" = "grey80", "1" = "red"),
    highlight_labels = c("Not ChIP-bound", "ChIP-bound"),
    xlab     = bquote(~"Correlation"~italic(r)),
    x_limits = c(-1, 1),
    sample_frac        = 0.2,
    sample_keep_levels = "1",
    point_aes = list(
      size  = 0.25,
      shape = c(16, 16, 16, 16),
      color = c("grey80", "red"),
      alpha = 0.6
    ),
    title   = paste0(tf, " ChIP-bound TFBS"),
    caption = paste("Total variables:", nrow(volcano_data)),
    legend_aes = list(
      labels     = c("Not ChIP-bound", "ChIP-bound"),
      position   = "bottom",
      label_size = 10,
      icon_size  = 4
    ),
    jitter = TRUE
  )

  # ---- Save Figures ----
  outfile1 <- file.path(base_dir, "predicted_canonical_tfbs", paste0(tf, "_volcano_with_ChIP.pdf"))
  outfile2 <- file.path(base_dir, "predicted_canonical_tfbs", paste0(tf, "_volcano_with_ChIP_jitter.pdf"))

  ggsave(outfile1, p, width = 10, height = 8, units = "in", dpi = 300)
  ggsave(outfile2, p_jitter, width = 10, height = 8, units = "in", dpi = 300)
}


# volcano_data have:
# - corr_fp_tf_r   (correlation)
# - corr_fp_tf_p   (p-value)
# - ChIP_bound     (0/1 ground truth)

optimize_corr_p_cutoffs <- function(dat,
                                    r_col   = "corr_fp_tf_r",
                                    p_col   = "corr_fp_tf_p",
                                    y_col   = "ChIP_bound",
                                    r_grid  = seq(0, 1, by = 0.02),
                                    p_grid  = 10^seq(-12, -1, length.out = 20),
                                    metric  = c("accuracy_overall",
                                                "precision_ChIP1",
                                                "recall_ChIP1",
                                                "f1_ChIP1")) {
  metric <- match.arg(metric)

  # Precompute
  abs_r   <- abs(dat[[r_col]])
  p_val   <- dat[[p_col]]
  y_true  <- as.integer(dat[[y_col]] == 1L)  # 1 = ChIP_bound, 0 = not

  res_list <- vector("list", length(r_grid) * length(p_grid))
  k <- 1L

  for (rc in r_grid) {
    for (pc in p_grid) {

      # classification rule:
      # predicted bound if |r| >= r_cut AND p <= p_cut
      pred <- as.integer(abs_r >= rc & p_val <= pc)

      tp <- sum(pred == 1L & y_true == 1L)
      fp <- sum(pred == 1L & y_true == 0L)
      fn <- sum(pred == 0L & y_true == 1L)
      tn <- sum(pred == 0L & y_true == 0L)

      accuracy_overall <- (tp + tn) / (tp + fp + fn + tn)
      precision_ChIP1  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      recall_ChIP1     <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      f1_ChIP1 <- if (is.na(precision_ChIP1) ||
                      is.na(recall_ChIP1) ||
                      (precision_ChIP1 + recall_ChIP1) == 0) {
        NA_real_
      } else {
        2 * precision_ChIP1 * recall_ChIP1 / (precision_ChIP1 + recall_ChIP1)
      }

      res_list[[k]] <- data.frame(
        r_cut            = rc,
        p_cut            = pc,
        accuracy_overall = accuracy_overall,
        precision_ChIP1  = precision_ChIP1,
        recall_ChIP1     = recall_ChIP1,
        f1_ChIP1         = f1_ChIP1
      )
      k <- k + 1L
    }
  }

  grid_results <- do.call(rbind, res_list)

  # pick the best according to the chosen metric
  best_idx <- switch(
    metric,
    "accuracy_overall" = which.max(grid_results$accuracy_overall),
    "precision_ChIP1"  = which.max(grid_results$precision_ChIP1),
    "recall_ChIP1"     = which.max(grid_results$recall_ChIP1),
    "f1_ChIP1"         = which.max(grid_results$f1_ChIP1)
  )

  best_row <- grid_results[best_idx, , drop = FALSE]

  list(
    grid_results = grid_results,
    best = best_row
  )
}

library(readr)

# TFs want to process
tfs <- c("HNF4A", "MAFF", "ZEB1")

opt_dir <- file.path(base_dir, "predicted_canonical_tfbs_cutoff_optimization")
if (!dir.exists(opt_dir)) {
  dir.create(opt_dir, recursive = TRUE)
}

target_metrics <- c("accuracy_overall",
                    "precision_ChIP1",
                    "recall_ChIP1",
                    "f1_ChIP1")

for (tf in tfs) {
  message("Optimizing cutoffs for TF: ", tf)

  # ----- load data for this TF -----
  infile <- file.path(
    base_dir, "predicted_canonical_tfbs",
    paste0(tf, "_overview_with_ChIP.csv")
  )
  volcano_data <- readr::read_csv(infile, show_col_types = FALSE)

  # ----- loop over target metrics -----
  for (metric in target_metrics) {
    message("  Metric: ", metric)

    res_opt <- optimize_corr_p_cutoffs(
      dat   = volcano_data,
      r_col = "corr_fp_tf_r",
      p_col = "corr_fp_tf_p",
      y_col = "ChIP_bound",
      metric = metric
    )

    grid_file <- file.path(
      opt_dir,
      paste0(tf, "_cutoff_grid_", metric, ".csv")
    )
    best_file <- file.path(
      opt_dir,
      paste0(tf, "_cutoff_best_", metric, ".csv")
    )

    readr::write_csv(res_opt$grid_results, grid_file)
    readr::write_csv(res_opt$best,         best_file)
  }
}

# assumes optimize_corr_p_cutoffs() is already defined as in previous step

tfs <- c("HNF4A", "MAFF", "ZEB1")

opt_dir <- file.path(base_dir, "predicted_canonical_tfbs_cutoff_optimization")
if (!dir.exists(opt_dir)) {
  dir.create(opt_dir, recursive = TRUE)
}

for (tf in tfs) {
  message("Optimizing + plotting metrics for TF: ", tf)

  # ----- load data -----
  infile <- file.path(
    base_dir, "predicted_canonical_tfbs",
    paste0(tf, "_overview_with_ChIP.csv")
  )
  volcano_data <- readr::read_csv(infile, show_col_types = FALSE)

  # ----- run optimization (as before) -----
  res_opt <- optimize_corr_p_cutoffs(
    dat   = volcano_data,
    r_col = "corr_fp_tf_r",
    p_col = "corr_fp_tf_p",
    y_col = "ChIP_bound",
    metric = "f1_ChIP1"
  )

  grid <- res_opt$grid_results

  # drop rows with NA F1 (usually extreme cutoffs with no positives)
  grid <- grid[!is.na(grid$f1_ChIP1), , drop = FALSE]

  # make r_cut a factor so facets are ordered nicely
  grid$r_cut <- factor(grid$r_cut)

  # ----- long format for 3 metrics -----
  plot_df <- tidyr::pivot_longer(
    grid,
    cols = c("accuracy_overall", "precision_ChIP1", "recall_ChIP1"),
    names_to  = "metric",
    values_to = "value"
  )

  # ----- scatter: x = p_cut, y = metric value, color = metric, facet by r_cut -----
  p_metrics <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = p_cut, y = value, color = metric)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0(tf, " cutoff optimization metrics"),
      x     = "p-value cutoff (log10 scale)",
      y     = "Metric value"
    ) +
    ggplot2::facet_wrap(~ r_cut, ncol = 5) +  # facet by |r| cutoff
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(size = 8)
    )

  out_pdf <- file.path(
    opt_dir,
    paste0(tf, "_cutoff_metrics_by_r.pdf")
  )
  ggplot2::ggsave(
    filename = out_pdf,
    plot     = p_metrics,
    width    = 10,
    height   = 18,
    units    = "in",
    dpi      = 300
  )
}


# ATAC peak correlation with TF -------------------------------------------

build_atac_rna_cor_overviews <- function(strict_build_args,
                                         base_dir,
                                         tf_names = strict_build_args$tf_list,
                                         r_cut = 0.3,
                                         p_cut = 0.01) {
  atac_score <- strict_build_args$atac_score
  rna        <- strict_build_args$rna

  # ---------------------------------------------------------------------------
  # 1) Match samples between ATAC and RNA
  # ---------------------------------------------------------------------------
  atac_cols <- colnames(atac_score)
  rna_cols  <- colnames(rna)

  # Numeric sample IDs are those that appear in both ATAC and RNA tables
  sample_ids <- intersect(atac_cols, rna_cols)
  sample_ids <- sample_ids[sample_ids %in% strict_build_args$metadata$id]

  if (length(sample_ids) == 0L) {
    cli::cli_abort("No shared sample IDs between ATAC and RNA matrices.")
  }

  if (!is.null(strict_build_args$expected_n)) {
    if (length(sample_ids) != strict_build_args$expected_n) {
      cli::cli_warn(
        c(
          "!" = "Number of matched samples ({length(sample_ids)})",
          " " = "does not equal expected_n ({strict_build_args$expected_n})."
        )
      )
    }
  }

  # Keep only matched samples
  atac_mat <- as.matrix(atac_score[, sample_ids, drop = FALSE])
  rna_mat  <- as.data.frame(rna[, c("ensembl_gene_id", "HGNC", sample_ids),
                                drop = FALSE])

  # ---------------------------------------------------------------------------
  # 2) Parse peak coordinates once
  # ---------------------------------------------------------------------------
  peak_ids <- atac_score[["atac_peak"]]

  peak_chr <- sub(":.*$", "", peak_ids)
  coord    <- sub("^.*:", "", peak_ids)
  peak_start <- as.integer(sub("-.*$", "", coord))
  peak_end   <- as.integer(sub("^.*-", "", coord))

  # ---------------------------------------------------------------------------
  # 3) Helper: compute correlations for one TF
  # ---------------------------------------------------------------------------
  cor_one_tf <- function(tf_symbol) {
    idx <- which(rna_mat[["HGNC"]] == tf_symbol)

    if (length(idx) == 0L) {
      cli::cli_warn("No RNA row found for TF {.val {tf_symbol}}; skipping.")
      return(NULL)
    }
    if (length(idx) > 1L) {
      cli::cli_warn(
        "Multiple RNA rows for TF {.val {tf_symbol}}; using the first match."
      )
      idx <- idx[1L]
    }

    tf_expr <- as.numeric(rna_mat[idx, sample_ids, drop = TRUE])

    # Center TF expression
    n <- length(tf_expr)
    tf_mean <- mean(tf_expr)
    tf_center <- tf_expr - tf_mean
    tf_ss <- sum(tf_center^2)

    # Center ATAC rows
    atac_means  <- rowMeans(atac_mat)
    atac_center <- atac_mat - atac_means

    # Numerator and denominator for Pearson r (vectorised across peaks)
    num <- atac_center %*% tf_center
    den <- sqrt(rowSums(atac_center^2) * tf_ss)

    r_vals <- as.numeric(num / den)
    r_vals[den == 0] <- NA_real_

    # p-values via t-statistic
    df <- n - 2
    # guard against r ~ +/-1
    denom_p <- pmax(1e-12, 1 - r_vals^2)
    t_stat <- r_vals * sqrt(df / denom_p)
    p_vals <- 2 * stats::pt(-abs(t_stat), df = df)
    p_vals[is.na(r_vals)] <- NA_real_

    p_adj <- stats::p.adjust(p_vals, method = "BH")

    # Bound flag
    bound_flag <- ifelse(
      !is.na(r_vals) & (abs(r_vals) >= r_cut) & (p_vals <= p_cut),
      1L,
      0L
    )

    out <- data.frame(
      peak_chr          = peak_chr,
      peak_start        = peak_start,
      peak_end          = peak_end,
      TF                = tf_symbol,
      corr_atac_tf_r    = r_vals,
      corr_atac_tf_p    = p_vals,
      corr_atac_tf_p_adj = p_adj,
      `_bound`          = bound_flag,
      stringsAsFactors  = FALSE
    )

    out
  }

  # ---------------------------------------------------------------------------
  # 4) Loop over TFs, save results
  # ---------------------------------------------------------------------------
  out_dir <- file.path(base_dir, "predicted_all_atac")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  res_list <- vector("list", length(tf_names))
  names(res_list) <- tf_names

  for (i in seq_along(tf_names)) {
    tf <- tf_names[[i]]
    cli::cli_inform("Processing TF {.val {tf}} ({i}/{length(tf_names)})")

    df_tf <- cor_one_tf(tf)
    if (is.null(df_tf)) {
      next
    }

    res_list[[i]] <- df_tf

    out_file <- file.path(out_dir, paste0(tf, "_overview.txt"))
    readr::write_tsv(df_tf, out_file)
  }

  invisible(res_list)
}

base_dir <- "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction"

res_atac <- build_atac_rna_cor_overviews(
  strict_build_args = strict_build_args,
  base_dir          = base_dir,
  tf_names          = c("HNF4A", "MAFF", "ZEB1"),
  r_cut             = 0.3,
  p_cut             = 0.01
)


batch <- list(
  list(
    tf   = "HNF4A",
    pred = file.path(base_dir, "predicted_all_atac", "HNF4A_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy83.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf   = "MAFF",
    pred = file.path(base_dir, "predicted_all_atac", "MAFF_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy84.hg38.rp10m.narrowpeaks.bed")
  ),
  list(
    tf   = "ZEB1",
    pred = file.path(base_dir, "predicted_all_atac", "ZEB1_overview.txt"),
    chip = file.path(base_dir, "cutntag", "cy76.hg38.rp10m.narrowpeaks.bed")
  )
)

out_dir <- file.path(base_dir, "plots_pdf_all_atac")

purrr::walk(batch, function(x) {
  save_benchmark_plots(
    tf        = x$tf,
    pred_path = x$pred,
    chip_path = x$chip,
    out_dir   = out_dir,
    score_col = "corr_atac_tf_r",
    mode      = "atac"
  )
})




tfs <- c("HNF4A", "MAFF", "ZEB1")

opt_dir <- file.path(base_dir, "predicted_canonical_tfbs_cutoff_optimization")
if (!dir.exists(opt_dir)) {
  dir.create(opt_dir, recursive = TRUE)
}

for (tf in tfs) {
  message("Optimizing + plotting metrics for TF: ", tf)

  # ----- load data -----
  infile <- file.path(
    base_dir, "predicted_canonical_tfbs",
    paste0(tf, "_overview_with_ChIP.csv")
  )
  volcano_data <- readr::read_csv(infile, show_col_types = FALSE)

  # ----- run optimization (as before) -----
  res_opt <- optimize_corr_p_cutoffs(
    dat   = volcano_data,
    r_col = "corr_fp_tf_r",
    p_col = "corr_fp_tf_p",
    y_col = "ChIP_bound",
    metric = "f1_ChIP1"
  )

  grid <- res_opt$grid_results

  # drop rows with NA F1 (usually extreme cutoffs with no positives)
  grid <- grid[!is.na(grid$f1_ChIP1), , drop = FALSE]

  # make r_cut a factor so facets are ordered nicely
  grid$r_cut <- factor(grid$r_cut)

  # ----- long format for 3 metrics -----
  plot_df <- tidyr::pivot_longer(
    grid,
    cols = c("accuracy_overall", "precision_ChIP1", "recall_ChIP1"),
    names_to  = "metric",
    values_to = "value"
  )

  # ----- scatter: x = p_cut, y = metric value, color = metric, facet by r_cut -----
  p_metrics <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = p_cut, y = value, color = metric)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = paste0(tf, " cutoff optimization metrics"),
      x     = "p-value cutoff (log10 scale)",
      y     = "Metric value"
    ) +
    ggplot2::facet_wrap(~ r_cut, ncol = 5) +  # facet by |r| cutoff
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(size = 8)
    )

  out_pdf <- file.path(
    opt_dir,
    paste0(tf, "_cutoff_metrics_by_r.pdf")
  )
  ggplot2::ggsave(
    filename = out_pdf,
    plot     = p_metrics,
    width    = 10,
    height   = 18,
    units    = "in",
    dpi      = 300
  )
}


optimize_corr_p_cutoffs <- function(dat,
                                    r_col   = "corr_fp_tf_r",
                                    p_col   = "corr_fp_tf_p",
                                    y_col   = "ChIP_bound",
                                    r_grid  = seq(0, 1, by = 0.02),
                                    p_grid  = 10^seq(-12, -1, length.out = 20),
                                    metric  = c("accuracy_overall",
                                                "precision_ChIP1",
                                                "recall_ChIP1",
                                                "f1_ChIP1")) {
  metric <- match.arg(metric)

  # Precompute
  abs_r  <- abs(dat[[r_col]])
  p_val  <- dat[[p_col]]
  y_true <- as.integer(dat[[y_col]] == 1L)  # 1 = ChIP_bound, 0 = not
  # y_true can contain NA if ChIP_bound is NA

  res_list <- vector("list", length(r_grid) * length(p_grid))
  k <- 1L

  for (rc in r_grid) {
    for (pc in p_grid) {

      # predicted bound if |r| >= r_cut AND p <= p_cut
      pred <- as.integer(abs_r >= rc & p_val <= pc)
      # pred can be 0/1/NA when abs_r or p_val is NA

      tp <- sum(pred == 1L & y_true == 1L, na.rm = TRUE)
      fp <- sum(pred == 1L & y_true == 0L, na.rm = TRUE)
      fn <- sum(pred == 0L & y_true == 1L, na.rm = TRUE)
      tn <- sum(pred == 0L & y_true == 0L, na.rm = TRUE)

      denom_all <- tp + fp + fn + tn
      accuracy_overall <- if (denom_all > 0) (tp + tn) / denom_all else NA_real_

      precision_ChIP1  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      recall_ChIP1     <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

      f1_ChIP1 <- if (is.na(precision_ChIP1) ||
                      is.na(recall_ChIP1) ||
                      (precision_ChIP1 + recall_ChIP1) == 0) {
        NA_real_
      } else {
        2 * precision_ChIP1 * recall_ChIP1 / (precision_ChIP1 + recall_ChIP1)
      }

      res_list[[k]] <- data.frame(
        r_cut            = rc,
        p_cut            = pc,
        accuracy_overall = accuracy_overall,
        precision_ChIP1  = precision_ChIP1,
        recall_ChIP1     = recall_ChIP1,
        f1_ChIP1         = f1_ChIP1
      )
      k <- k + 1L
    }
  }

  grid_results <- do.call(rbind, res_list)

  # pick the best according to the chosen metric
  best_idx <- switch(
    metric,
    "accuracy_overall" = which.max(grid_results$accuracy_overall),
    "precision_ChIP1"  = which.max(grid_results$precision_ChIP1),
    "recall_ChIP1"     = which.max(grid_results$recall_ChIP1),
    "f1_ChIP1"         = which.max(grid_results$f1_ChIP1)
  )

  best_row <- grid_results[best_idx, , drop = FALSE]

  list(
    grid_results = grid_results,
    best = best_row
  )
}


library(readr)

# TFs to process
tfs <- c("HNF4A", "MAFF", "ZEB1")

# -------------------------------------------------------------------
# 1) GRID SEARCH + SAVE CSVs for both FP and ATAC
# -------------------------------------------------------------------

target_metrics <- c("accuracy_overall",
                    "precision_ChIP1",
                    "recall_ChIP1",
                    "f1_ChIP1")

for (mode in c("atac")) { # c("fp", "atac")
  message("=== Mode: ", mode, " ===")

  if (mode == "fp") {
    pred_subdir <- "predicted_canonical_tfbs"
    opt_dir     <- file.path(base_dir, "predicted_canonical_tfbs_cutoff_optimization")
    r_col       <- "corr_fp_tf_r"
    p_col       <- "corr_fp_tf_p"
  } else if (mode == "atac") {
    # ATAC mode: use the ATAC overview-with-ChIP tables and ATAC columns
    pred_subdir <- "predicted_all_atac"
    opt_dir     <- file.path(base_dir, "predicted_all_atac_cutoff_optimization")
    r_col       <- "corr_atac_tf_r"
    p_col       <- "corr_atac_tf_p"
  }

  if (!dir.exists(opt_dir)) {
    dir.create(opt_dir, recursive = TRUE)
  }

  for (tf in tfs) {
    message("Optimizing cutoffs for TF: ", tf, " (mode = ", mode, ")")

    infile <- file.path(
      base_dir, pred_subdir,
      paste0(tf, "_overview_with_ChIP.csv")
    )
    volcano_data <- readr::read_csv(infile, show_col_types = FALSE)

    # loop over target metrics
    for (metric in target_metrics) {
      message("  Metric: ", metric)

      res_opt <- optimize_corr_p_cutoffs(
        dat    = volcano_data,
        r_col  = r_col,
        p_col  = p_col,
        y_col  = "ChIP_bound",
        metric = metric
      )

      grid_file <- file.path(
        opt_dir,
        paste0(tf, "_cutoff_grid_", metric, ".csv")
      )
      best_file <- file.path(
        opt_dir,
        paste0(tf, "_cutoff_best_", metric, ".csv")
      )

      readr::write_csv(res_opt$grid_results, grid_file)
      readr::write_csv(res_opt$best,         best_file)
    }
  }
}

# -------------------------------------------------------------------
# 2) PLOTTING METRICS (accuracy/precision/recall) FACETTED BY r_cut
#    for both FP and ATAC
# -------------------------------------------------------------------

for (mode in c("atac")) { # c("fp", "atac")
  message("=== Mode (plots): ", mode, " ===")

  if (mode == "fp") {
    pred_subdir <- "predicted_canonical_tfbs"
    opt_dir     <- file.path(base_dir, "predicted_canonical_tfbs_cutoff_optimization")
    r_col       <- "corr_fp_tf_r"
    p_col       <- "corr_fp_tf_p"
  } else if (mode == "atac") {
    pred_subdir <- "predicted_all_atac"
    opt_dir     <- file.path(base_dir, "predicted_all_atac_cutoff_optimization")
    r_col       <- "corr_atac_tf_r"
    p_col       <- "corr_atac_tf_p"
  }

  if (!dir.exists(opt_dir)) {
    dir.create(opt_dir, recursive = TRUE)
  }

  for (tf in tfs) {
    message("Optimizing + plotting metrics for TF: ", tf, " (mode = ", mode, ")")

    infile <- file.path(
      base_dir, pred_subdir,
      paste0(tf, "_overview_with_ChIP.csv")
    )
    volcano_data <- readr::read_csv(infile, show_col_types = FALSE)

    # run optimization once for plotting (e.g. target F1)
    res_opt <- optimize_corr_p_cutoffs(
      dat    = volcano_data,
      r_col  = r_col,
      p_col  = p_col,
      y_col  = "ChIP_bound",
      metric = "f1_ChIP1"
    )

    grid <- res_opt$grid_results
    grid <- grid[!is.na(grid$f1_ChIP1), , drop = FALSE]
    grid$r_cut <- factor(grid$r_cut)

    # long format for accuracy / precision / recall
    plot_df <- tidyr::pivot_longer(
      grid,
      cols = c("accuracy_overall", "precision_ChIP1", "recall_ChIP1"),
      names_to  = "metric",
      values_to = "value"
    )

    p_metrics <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = p_cut, y = value, color = metric)
    ) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        title = paste0(tf, " cutoff optimization metrics (mode = ", mode, ")"),
        x     = "p-value cutoff (log10 scale)",
        y     = "Metric value"
      ) +
      ggplot2::facet_wrap(~ r_cut, ncol = 5) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position   = "bottom",
        strip.background  = ggplot2::element_rect(fill = "grey90"),
        strip.text        = ggplot2::element_text(size = 8)
      )

    out_pdf <- file.path(
      opt_dir,
      paste0(tf, "_cutoff_metrics_by_r_", mode, ".pdf")
    )
    ggplot2::ggsave(
      filename = out_pdf,
      plot     = p_metrics,
      width    = 10,
      height   = 18,
      units    = "in",
      dpi      = 300
    )
  }
}


