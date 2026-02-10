# ATAC-signal benchmark for TFBS prediction (canonical mode)
# -----------------------------------------------------------------------------
# DATA CONTEXT (expected when running this script)
# - Dataset: nutrient_stress_strict, grouped by strict RNA-matched conditions.
# - Grouping label: `label_col = "strict_match_rna"` (default).
# - Expected number of condition groups after aggregation: 23.
# - Typical source object: `omics_data` produced in:
#   dev/pipeline/11_pipeline_nutrient_strict.R
#   via correlate_tf_to_fp(..., mode = "canonical", label_col = "strict_match_rna").
#
# REQUIRED `omics_data` SLOTS USED HERE
# 1) omics_data$atac_score
#    - columns: `atac_peak` + sample id columns (e.g., cy410, cy411, ...)
#    - role: ATAC signal matrix used for ATAC-vs-TF expression correlations.
# 2) omics_data$rna
#    - columns: `HGNC` + sample id columns matching metadata ids.
#    - role: TF expression matrix; target TF rows (e.g., HNF4A, MAFF) are used.
# 3) omics_data$fp_annotation
#    - columns: `fp_peak`, `atac_peak`, `motifs` (plus optional extras).
#    - role: motif-to-peak annotation; links canonical motif hits to ATAC peaks.
# 4) omics_data$motif_db
#    - columns: `motif`, `gene_symbol`
#    - role: maps motif IDs to TF symbols; supports dimer entries (e.g., AHR::ARNT).
# 5) omics_data$sample_metadata_used
#    - columns: `id`, `strict_match_rna` (or selected `label_col`)
#    - role: maps sample ids to condition labels for per-condition aggregation.
#
# OPTIONAL INPUTS
# - `atac_overlap_file` + `atac_overlap_col`: optional external ATAC filter table
#   (disabled by default in this script).
# - `chip_files`: CUT&Tag/ChIP BED paths for benchmark truth labels per TF.
#
# OUTPUT
# - Writes ATAC benchmark CSVs and PDF plots to:
#   /data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/atac_vs_fp
# -----------------------------------------------------------------------------
# Self-sufficient: input is `omics_data` and output is written under `out_dir`.

benchmark_tf_binding_sites_atac <- function(
    omics_data,
    label_col = "strict_match_rna",
    tf_targets = c("HNF4A", "MAFF"),
    methods = c("pearson", "spearman", "pearson_or_spearman", "pearson_and_spearman"),
    chip_files = c(
      HNF4A = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/cutntag/cy83.hg38.rp10m.narrowpeaks.bed",
      MAFF  = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/cutntag/cy84.hg38.rp10m.narrowpeaks.bed"
    ),
    atac_overlap_file = NULL,
    atac_overlap_col = NULL,
    out_dir = "/data/homes/yl814/episcope_test/benchmark_tf_binding_sites_prediction/atac_vs_fp",
    r_grid = c(-1, 0, 0.3, 0.5, 0.7),
    p_grid = 10^seq(-4, -1, length.out = 20),
    min_non_na = 5L,
    mode_label = "canonical",
    debug_print_grid = FALSE,
    verbose = TRUE
) {
  .log_inform <- function(...) {
    if (isTRUE(verbose)) {
      cli::cli_inform(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", sprintf(...)))
    }
  }
  .stop <- function(...) cli::cli_abort(sprintf(...))

  .require_cols <- function(df, cols, where) {
    miss <- setdiff(cols, names(df))
    if (length(miss)) .stop("Missing required columns in %s: %s", where, paste(miss, collapse = ", "))
  }

  .split_gene_symbols <- function(x) {
    x <- as.character(x)
    out <- strsplit(x, "::", fixed = TRUE)
    lapply(out, function(v) unique(trimws(v[nzchar(v)])))
  }

  .parse_peak <- function(x, prefix) {
    m <- stringr::str_match(x, "^([^:]+):(\\d+)-(\\d+)$")
    bad <- which(!is.na(x) & (is.na(m[, 2]) | is.na(m[, 3]) | is.na(m[, 4])))
    if (length(bad)) {
      .stop("Could not parse %s coordinates for %d rows (example: %s)", prefix, length(bad), x[bad[[1]]])
    }
    tibble::tibble(
      chrom = m[, 2],
      start = as.integer(m[, 3]),
      end = as.integer(m[, 4])
    )
  }

  .aggregate_by_label <- function(tbl, id_col, metadata, label_col, value_name) {
    .require_cols(tbl, id_col, sprintf("%s table", value_name))
    .require_cols(metadata, c("id", label_col), "sample_metadata_used")

    sample_cols <- setdiff(names(tbl), id_col)
    meta_use <- metadata[, c("id", label_col), drop = FALSE]
    meta_use <- meta_use[!is.na(meta_use[[label_col]]) & nzchar(meta_use[[label_col]]), , drop = FALSE]
    keep_ids <- intersect(meta_use$id, sample_cols)
    if (!length(keep_ids)) {
      .stop("No overlapping sample ids between %s and sample_metadata_used", value_name)
    }

    meta_use <- meta_use[match(keep_ids, meta_use$id), , drop = FALSE]
    dt <- data.table::as.data.table(tbl[, c(id_col, keep_ids), drop = FALSE])
    long <- data.table::melt(
      dt,
      id.vars = id_col,
      variable.name = "id",
      value.name = "value"
    )
    meta_dt <- data.table::as.data.table(meta_use)
    data.table::setnames(meta_dt, old = label_col, new = "label")
    long <- meta_dt[long, on = "id", nomatch = 0L]

    agg <- long[, .(value = mean(as.numeric(value), na.rm = TRUE)), by = c(id_col, "label")]
    wide <- data.table::dcast(agg, formula = stats::as.formula(sprintf("%s ~ label", id_col)), value.var = "value")
    tibble::as_tibble(wide)
  }

  .as_bed3 <- function(df, c_chrom, c_start, c_end) {
    bed <- tibble::tibble(
      chrom = as.character(df[[c_chrom]]),
      start = as.integer(df[[c_start]]),
      end = as.integer(df[[c_end]])
    )
    class(bed) <- c("tbl_df", "tbl", "data.frame", "bed_frame")
    bed
  }

  .collapse_unique_sites <- function(bed) bed[!duplicated(bed[c("chrom", "start", "end")]), , drop = FALSE]

  .read_chip_bed3 <- function(path) {
    if (!file.exists(path)) .stop("ChIP file not found: %s", path)
    df <- readr::read_tsv(
      path,
      comment = "#",
      col_names = FALSE,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )
    if (ncol(df) < 3L) .stop("ChIP file must have at least 3 columns: %s", path)
    names(df)[1:3] <- c("chrom", "start", "end")
    df$start <- suppressWarnings(as.integer(df$start))
    df$end <- suppressWarnings(as.integer(df$end))
    if (any(!is.finite(df$start) | !is.finite(df$end))) {
      .stop("Non-integer ChIP start/end in: %s", path)
    }
    bad <- which(df$end < df$start)
    if (length(bad)) {
      s <- df$start[bad]
      e <- df$end[bad]
      df$start[bad] <- pmin(s, e)
      df$end[bad] <- pmax(s, e)
    }
    bed <- tibble::as_tibble(df[, c("chrom", "start", "end"), drop = FALSE])
    class(bed) <- c("tbl_df", "tbl", "data.frame", "bed_frame")
    bed
  }

  .annotate_atac_with_chip <- function(
      pred_df,
      chip_bed,
      overlap_frac_atac = 0.5,
      overlap_frac_chip = 0.5
  ) {
    required <- c("peak_chr", "peak_start", "peak_end")
    .require_cols(pred_df, required, "prediction table")

    atac_all <- pred_df[, required, drop = FALSE] |>
      dplyr::distinct()

    chip_peaks <- .read_chip_bed3(chip_bed)
    atac_bed <- .as_bed3(atac_all, "peak_chr", "peak_start", "peak_end")
    ov <- valr::bed_intersect(atac_bed, chip_peaks)

    if (nrow(ov) > 0L) {
      ov2 <- ov |>
        dplyr::mutate(
          atac_len = pmax(1L, end.x - start.x),
          chip_len = pmax(1L, end.y - start.y),
          frac_atac = .overlap / atac_len,
          frac_chip = .overlap / chip_len,
          chip_hit = (frac_atac >= overlap_frac_atac) & (frac_chip >= overlap_frac_chip)
        ) |>
        dplyr::group_by(chrom, start.x, end.x) |>
        dplyr::slice_max(frac_atac, with_ties = FALSE) |>
        dplyr::ungroup()

      annot <- ov2 |>
        dplyr::transmute(
          peak_chr = chrom,
          peak_start = start.x,
          peak_end = end.x,
          chip_overlap_bp = as.integer(.overlap),
          chip_overlap_frac = frac_atac,
          chip_overlap_frac_chip = frac_chip,
          chip_bound = ifelse(chip_hit, 1L, 0L),
          chip_chr = dplyr::if_else(chip_hit, chrom, NA_character_),
          chip_start = dplyr::if_else(chip_hit, start.y, NA_integer_),
          chip_end = dplyr::if_else(chip_hit, end.y, NA_integer_)
        )
    } else {
      annot <- atac_all |>
        dplyr::mutate(
          chip_overlap_bp = 0L,
          chip_overlap_frac = 0,
          chip_overlap_frac_chip = 0,
          chip_bound = 0L,
          chip_chr = NA_character_,
          chip_start = NA_integer_,
          chip_end = NA_integer_
        )
    }

    atac_annot <- atac_all |>
      dplyr::left_join(annot, by = c("peak_chr", "peak_start", "peak_end")) |>
      dplyr::mutate(
        chip_overlap_bp = dplyr::coalesce(chip_overlap_bp, 0L),
        chip_overlap_frac = dplyr::coalesce(chip_overlap_frac, 0),
        chip_overlap_frac_chip = dplyr::coalesce(chip_overlap_frac_chip, 0),
        chip_bound = dplyr::coalesce(chip_bound, 0L),
        chip_chr = dplyr::coalesce(chip_chr, NA_character_),
        chip_start = dplyr::coalesce(chip_start, NA_integer_),
        chip_end = dplyr::coalesce(chip_end, NA_integer_),
        chip_peak = dplyr::if_else(
          chip_bound == 1L,
          sprintf("%s:%d-%d", chip_chr, chip_start, chip_end),
          NA_character_
        ),
        chip_bound_TP = dplyr::if_else((chip_overlap_frac >= overlap_frac_atac) & (chip_overlap_frac_chip >= overlap_frac_chip), 1L, 0L),
        chip_bound_FN = dplyr::if_else((chip_overlap_frac >= overlap_frac_atac) & (chip_overlap_frac_chip >= overlap_frac_chip), 1L, 0L),
        chip_bound_FP = dplyr::if_else((chip_overlap_frac >= overlap_frac_atac) & (chip_overlap_frac_chip >= overlap_frac_chip), 1L, 0L),
        chip_bound_TN = dplyr::if_else((chip_overlap_frac >= overlap_frac_atac) & (chip_overlap_frac_chip >= overlap_frac_chip), 1L, 0L)
      )

    pred_df |>
      dplyr::left_join(
        atac_annot |>
          dplyr::select(
            peak_chr, peak_start, peak_end,
            chip_overlap_bp, chip_overlap_frac,
            chip_overlap_frac_chip,
            chip_bound, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN,
            chip_chr, chip_start, chip_end, chip_peak
          ),
        by = c("peak_chr", "peak_start", "peak_end")
      )
  }

  .count_confusion <- function(predicted_bound, predicted_unbound_eval) {
    tp <- predicted_bound |>
      dplyr::filter(chip_bound_TP == 1) |>
      dplyr::distinct(peak_chr, peak_start, peak_end) |>
      nrow()
    fp <- predicted_bound |>
      dplyr::filter(!is.na(chip_bound_FP) & chip_bound_FP == 0) |>
      dplyr::distinct(peak_chr, peak_start, peak_end) |>
      nrow()
    tn <- predicted_unbound_eval |>
      dplyr::filter(!is.na(chip_bound_TN) & chip_bound_TN == 0) |>
      dplyr::distinct(peak_chr, peak_start, peak_end) |>
      nrow()
    fn <- predicted_unbound_eval |>
      dplyr::filter(chip_bound_FN == 1) |>
      dplyr::distinct(peak_chr, peak_start, peak_end) |>
      nrow()
    c(TP = tp, FP = fp, TN = tn, FN = fn)
  }

  .corr_rows <- function(mat, y, method = c("pearson", "spearman"), min_non_na = 5L) {
    method <- match.arg(method)
    y <- as.numeric(y)
    idx_y <- is.finite(y)
    mat <- mat[, idx_y, drop = FALSE]
    y <- y[idx_y]

    res <- apply(mat, 1L, function(x) {
      keep <- is.finite(x) & is.finite(y)
      n <- sum(keep)
      if (n < min_non_na) return(c(r = NA_real_, n = n))
      r <- suppressWarnings(stats::cor(x[keep], y[keep], method = method))
      c(r = r, n = n)
    })
    r <- as.numeric(res["r", ])
    n <- as.numeric(res["n", ])
    r_clip <- pmin(pmax(r, -0.999999), 0.999999)
    tval <- r_clip * sqrt((n - 2) / pmax(1e-12, 1 - r_clip^2))
    pval <- rep(NA_real_, length(r))
    ok <- is.finite(r) & is.finite(n) & (n > 2)
    pval[ok] <- 2 * stats::pt(-abs(tval[ok]), df = n[ok] - 2)
    padj <- rep(NA_real_, length(pval))
    good <- is.finite(pval)
    if (any(good)) padj[good] <- stats::p.adjust(pval[good], method = "BH")

    tibble::tibble(corr_fp_tf_r = r, corr_fp_tf_p = pval, corr_fp_tf_p_adj = padj)
  }

  if (!is.list(omics_data)) .stop("`omics_data` must be a list.")
  .require_cols(as.data.frame(omics_data$fp_annotation), c("fp_peak", "atac_peak", "motifs"), "omics_data$fp_annotation")
  .require_cols(as.data.frame(omics_data$motif_db), c("motif", "gene_symbol"), "omics_data$motif_db")
  .require_cols(as.data.frame(omics_data$atac_score), c("atac_peak"), "omics_data$atac_score")
  .require_cols(as.data.frame(omics_data$rna), c("HGNC"), "omics_data$rna")
  .require_cols(as.data.frame(omics_data$sample_metadata_used), c("id", label_col), "omics_data$sample_metadata_used")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  .log_inform("Preparing ATAC and RNA matrices grouped by '%s'.", label_col)
  atac_cond <- .aggregate_by_label(
    tbl = as.data.frame(omics_data$atac_score),
    id_col = "atac_peak",
    metadata = as.data.frame(omics_data$sample_metadata_used),
    label_col = label_col,
    value_name = "atac_score"
  )

  rna_cond <- .aggregate_by_label(
    tbl = as.data.frame(omics_data$rna),
    id_col = "HGNC",
    metadata = as.data.frame(omics_data$sample_metadata_used),
    label_col = label_col,
    value_name = "rna"
  ) |>
    dplyr::distinct(HGNC, .keep_all = TRUE)

  cond_cols <- intersect(setdiff(names(atac_cond), "atac_peak"), setdiff(names(rna_cond), "HGNC"))
  if (length(cond_cols) < min_non_na) {
    .stop("Insufficient shared condition columns: %d found, min_non_na=%d", length(cond_cols), min_non_na)
  }

  atac_mat <- as.matrix(atac_cond[, cond_cols, drop = FALSE])
  storage.mode(atac_mat) <- "double"

  motif_map <- as.data.frame(omics_data$motif_db) |>
    dplyr::transmute(motif = as.character(motif), gene_symbol = as.character(gene_symbol))
  motif_expanded <- motif_map |>
    dplyr::mutate(tfs = .split_gene_symbols(gene_symbol)) |>
    tidyr::unnest(tfs) |>
    dplyr::mutate(tfs = toupper(tfs))

  fp_ann <- as.data.frame(omics_data$fp_annotation) |>
    dplyr::transmute(
      fp_peak = as.character(fp_peak),
      atac_peak = as.character(atac_peak),
      motifs = as.character(motifs)
    )

  use_atac_overlap_filter <- !is.null(atac_overlap_file) && !is.null(atac_overlap_col) &&
    nzchar(atac_overlap_file) && nzchar(atac_overlap_col)
  atac_overlap <- NULL
  if (isTRUE(use_atac_overlap_filter)) {
    if (!file.exists(atac_overlap_file)) .stop("ATAC overlap file not found: %s", atac_overlap_file)
    atac_overlap <- data.table::fread(atac_overlap_file, select = c("Chr", "Start", "End", atac_overlap_col))
    if (!atac_overlap_col %in% names(atac_overlap)) {
      .stop("Column '%s' not found in %s", atac_overlap_col, atac_overlap_file)
    }
    .log_inform("Applying ATAC overlap filter using '%s'.", atac_overlap_col)
  } else {
    .log_inform("ATAC overlap filter disabled.")
  }

  all_grids <- list()

  for (tf in tf_targets) {
    tf_u <- toupper(tf)
    if (!tf_u %in% names(chip_files)) .stop("Missing chip file path for TF: %s", tf_u)
    chip_file <- unname(chip_files[[tf_u]])

    .log_inform("Building ATAC-correlation predictions for %s.", tf_u)

    tf_motifs <- motif_expanded |>
      dplyr::filter(tfs == tf_u) |>
      dplyr::distinct(motif) |>
      dplyr::pull(motif)
    if (!length(tf_motifs)) {
      .stop("No motif mapping found for TF '%s' in omics_data$motif_db", tf_u)
    }

    tf_rows <- fp_ann |>
      dplyr::filter(motifs %in% tf_motifs) |>
      dplyr::distinct(fp_peak, atac_peak, motifs)
    if (!nrow(tf_rows)) {
      .stop("No canonical fp_annotation rows for TF '%s'", tf_u)
    }

    fp_coord <- .parse_peak(tf_rows$fp_peak, "fp_peak")
    atac_coord <- .parse_peak(tf_rows$atac_peak, "atac_peak")
    tf_rows <- dplyr::bind_cols(
      tf_rows,
      fp_coord |>
        dplyr::rename(TFBS_chr = chrom, TFBS_start = start, TFBS_end = end),
      atac_coord |>
        dplyr::rename(peak_chr = chrom, peak_start = start, peak_end = end)
    )

    y <- rna_cond |>
      dplyr::filter(HGNC == tf_u) |>
      dplyr::slice(1) |>
      dplyr::select(dplyr::all_of(cond_cols))
    if (!nrow(y)) .stop("TF '%s' not found in RNA table (HGNC)", tf_u)
    y_vec <- as.numeric(y[1, ])

    corr_by_method <- list()
    pred_by_method <- list()

    for (m in c("pearson", "spearman")) {
      .log_inform("Computing %s ATAC-vs-%s correlations.", m, tf_u)
      corr_tbl <- .corr_rows(atac_mat, y_vec, method = m, min_non_na = min_non_na)
      corr_tbl <- tibble::tibble(atac_peak = atac_cond$atac_peak) |>
        dplyr::bind_cols(corr_tbl)
      corr_by_method[[m]] <- corr_tbl

      pred_tbl <- tf_rows |>
        dplyr::left_join(corr_tbl, by = "atac_peak") |>
        dplyr::mutate(
          tf = tf_u,
          method = m,
          mode = mode_label
        )

      pred_tbl <- .annotate_atac_with_chip(
        pred_tbl,
        chip_bed = chip_file,
        overlap_frac_atac = 0.5,
        overlap_frac_chip = 0.5
      )

      if (isTRUE(use_atac_overlap_filter)) {
        pred_tbl <- pred_tbl |>
          dplyr::left_join(
            tibble::as_tibble(atac_overlap),
            by = c("peak_chr" = "Chr", "peak_start" = "Start", "peak_end" = "End")
          )
        pred_tbl <- pred_tbl[!is.na(pred_tbl[[atac_overlap_col]]) & pred_tbl[[atac_overlap_col]] == 1, , drop = FALSE]
      }

      pred_by_method[[m]] <- pred_tbl

      out_pred <- file.path(out_dir, sprintf("%s_%s_%s_ATAC_chip50_annotation_fullrows.csv", tf_u, m, mode_label))
      readr::write_csv(pred_tbl, out_pred)
    }

    method_work <- methods[methods %in% c("pearson", "spearman", "pearson_or_spearman", "pearson_and_spearman")]

    for (m in method_work) {
      .log_inform("Running cutoff grid for %s (%s).", tf_u, m)

      df_p <- pred_by_method[["pearson"]]
      df_s <- pred_by_method[["spearman"]]

      chip_lookup <- df_p |>
        dplyr::select(peak_chr, peak_start, peak_end, chip_bound_TP, chip_bound_FN, chip_bound_FP, chip_bound_TN) |>
        dplyr::distinct()

      grid <- vector("list", length(r_grid) * length(p_grid))
      k <- 1L

      for (r in r_grid) {
        for (p in p_grid) {
          if (m == "pearson") {
            predicted_bound <- df_p |>
              dplyr::filter(corr_fp_tf_r > r & corr_fp_tf_p_adj < p) |>
              dplyr::distinct(peak_chr, peak_start, peak_end)
          } else if (m == "spearman") {
            predicted_bound <- df_s |>
              dplyr::filter(corr_fp_tf_r > r & corr_fp_tf_p_adj < p) |>
              dplyr::distinct(peak_chr, peak_start, peak_end)
          } else {
            pred_p <- df_p |>
              dplyr::filter(corr_fp_tf_r > r & corr_fp_tf_p_adj < p) |>
              dplyr::distinct(peak_chr, peak_start, peak_end)
            pred_s <- df_s |>
              dplyr::filter(corr_fp_tf_r > r & corr_fp_tf_p_adj < p) |>
              dplyr::distinct(peak_chr, peak_start, peak_end)
            predicted_bound <- if (m == "pearson_and_spearman") {
              dplyr::inner_join(pred_p, pred_s, by = c("peak_chr", "peak_start", "peak_end"))
            } else {
              dplyr::bind_rows(pred_p, pred_s) |>
                dplyr::distinct(peak_chr, peak_start, peak_end)
            }
          }

          predicted_bound <- predicted_bound |>
            dplyr::left_join(chip_lookup, by = c("peak_chr", "peak_start", "peak_end"))

          if (m == "pearson") {
            universe <- df_p
          } else if (m == "spearman") {
            universe <- df_s
          } else {
            universe <- dplyr::bind_rows(df_p, df_s) |>
              dplyr::distinct(peak_chr, peak_start, peak_end, .keep_all = TRUE)
          }

          predicted_unbound_eval <- dplyr::anti_join(
            universe,
            predicted_bound,
            by = c("peak_chr", "peak_start", "peak_end")
          ) |>
            dplyr::filter(is.finite(corr_fp_tf_r) & is.finite(corr_fp_tf_p_adj))

          cnt <- .count_confusion(predicted_bound, predicted_unbound_eval)
          tp <- cnt[["TP"]]
          fp <- cnt[["FP"]]
          tn <- cnt[["TN"]]
          fn <- cnt[["FN"]]
          denom <- tp + fp + tn + fn

          grid[[k]] <- tibble::tibble(
            tf = tf_u,
            method = m,
            mode = mode_label,
            r = r,
            p = p,
            ATAC_TP = tp,
            ATAC_FP = fp,
            ATAC_TN = tn,
            ATAC_FN = fn,
            ATAC_accuracy = if (denom > 0) (tp + tn) / denom else NA_real_,
            ATAC_precision = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
            ATAC_recall = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
          )
          k <- k + 1L
        }
      }

      grid_df <- dplyr::bind_rows(grid)
      all_grids[[paste(tf_u, m, sep = "::")]] <- grid_df

      grid_file <- file.path(out_dir, sprintf("%s_%s_%s_cutoff_grid_atac.csv", tf_u, m, mode_label))
      readr::write_csv(grid_df, grid_file)

      if (isTRUE(debug_print_grid)) {
        .log_inform("Debug grid: %s (%s)", tf_u, m)
        print(
          grid_df |>
            dplyr::select(
              tf, method, mode, r, p,
              ATAC_accuracy, ATAC_precision, ATAC_recall,
              ATAC_FP, ATAC_TN, ATAC_FN, ATAC_TP
            ) |>
            dplyr::arrange(r, p)
        )
      }

      r_vals <- sort(unique(grid_df$r))
      r_labs <- paste0("r > ", r_vals)
      metric_pal <- c(accuracy = "green3", precision = "red2", recall = "blue3")
      conf_pal <- c(TP = "#31a354", FP = "#de2d26", TN = "#3182bd", FN = "#ff9800")
      x_breaks <- scales::trans_breaks("log10", function(x) 10^x)
      x_labels <- scales::trans_format("log10", scales::math_format(bold(10^.x)))

      plot_df <- grid_df |>
        dplyr::select(r, p, ATAC_accuracy, ATAC_precision, ATAC_recall) |>
        tidyr::pivot_longer(
          cols = -c(r, p),
          names_to = "metric",
          values_to = "value"
        ) |>
        dplyr::mutate(
          metric = dplyr::recode(
            metric,
            ATAC_accuracy = "accuracy",
            ATAC_precision = "precision",
            ATAC_recall = "recall"
          ),
          metric = factor(metric, levels = c("accuracy", "precision", "recall")),
          r_lab = factor(r, levels = r_vals, labels = r_labs),
          line_group = factor("Canonical_ATAC", levels = "Canonical_ATAC")
        )

      p_metrics <- ggplot2::ggplot(plot_df, ggplot2::aes(x = p, y = value, color = metric, linetype = line_group)) +
        ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
        ggplot2::scale_x_log10(breaks = x_breaks, labels = x_labels) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::facet_wrap(~ r_lab, ncol = length(r_vals)) +
        ggplot2::scale_color_manual(values = metric_pal) +
        ggplot2::scale_linetype_manual(values = c(Canonical_ATAC = "dashed")) +
        ggplot2::labs(
          title = sprintf("%s cutoff optimization metrics", tf_u),
          x = "p-value cutoff",
          y = "Metric value",
          color = "color_group",
          linetype = "line_group"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "bottom",
          legend.box = "vertical",
          strip.background = ggplot2::element_rect(fill = "grey90"),
          strip.text = ggplot2::element_text(size = 9, face = "bold"),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(face = "bold"),
          legend.title = ggplot2::element_text(face = "bold"),
          legend.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(order = 1), linetype = ggplot2::guide_legend(order = 2))

      counts_df <- grid_df |>
        dplyr::select(r, p, ATAC_TP, ATAC_FP, ATAC_TN, ATAC_FN) |>
        tidyr::pivot_longer(
          cols = -c(r, p),
          names_to = "component",
          values_to = "count"
        ) |>
        dplyr::mutate(
          component = dplyr::recode(
            component,
            ATAC_TP = "TP",
            ATAC_FP = "FP",
            ATAC_TN = "TN",
            ATAC_FN = "FN"
          ),
          component = factor(component, levels = c("FP", "TN", "FN", "TP")),
          count_signed = ifelse(component %in% c("TP", "FN"), count, -count),
          r_lab = factor(r, levels = r_vals, labels = r_labs)
        )

      p_all_counts_auto <- ggplot2::ggplot(counts_df, ggplot2::aes(x = p, y = count_signed, fill = component)) +
        ggplot2::geom_col() +
        ggplot2::scale_x_log10(breaks = x_breaks, labels = x_labels) +
        ggplot2::facet_wrap(~ r_lab, ncol = length(r_vals)) +
        ggplot2::scale_fill_manual(values = conf_pal) +
        ggplot2::labs(x = "p-value cutoff", y = "All (|y| = count)", fill = "component") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "none",
          strip.background = ggplot2::element_rect(fill = "grey90"),
          strip.text = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(face = "bold")
        )

      p_can_counts_auto <- ggplot2::ggplot(counts_df, ggplot2::aes(x = p, y = count_signed, fill = component)) +
        ggplot2::geom_col() +
        ggplot2::scale_x_log10(breaks = x_breaks, labels = x_labels) +
        ggplot2::facet_wrap(~ r_lab, ncol = length(r_vals)) +
        ggplot2::scale_fill_manual(values = conf_pal) +
        ggplot2::labs(x = "p-value cutoff", y = "Canonical (|y| = count)", fill = "component") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          legend.position = "bottom",
          strip.background = ggplot2::element_rect(fill = "grey90"),
          strip.text = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(face = "bold"),
          axis.title.y = ggplot2::element_text(face = "bold"),
          axis.text.x = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(face = "bold"),
          legend.title = ggplot2::element_text(face = "bold"),
          legend.text = ggplot2::element_text(face = "bold")
        )

      combined <- p_metrics / p_all_counts_auto / p_can_counts_auto +
        patchwork::plot_layout(heights = c(1, 1, 1))

      ggplot2::ggsave(
        filename = file.path(out_dir, sprintf("%s_%s_%s_cutoff_metrics_and_counts_by_r_atac.pdf", tf_u, m, mode_label)),
        plot = combined,
        width = 12,
        height = 10,
        dpi = 600
      )
    }
  }

  invisible(list(cutoff_grids = all_grids, out_dir = out_dir))
}

res <- benchmark_tf_binding_sites_atac(
  omics_data = omics_data,
  label_col = "strict_match_rna",
  debug_print_grid = TRUE)
