# 5-method TF->target KO benchmark plotting for HNF1A and SOX9 only.
# Methods: ATAC Pearson, RNA Pearson, FP Pearson, FP Spearman, FP/RNA hybrid.
# Reads isolated method tables from lighting_method_comparison_atac_benchmark.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# Standalone defaults (edit here if needed)
workspace_rdata <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/workspace.RData"
method_dir <- "/data/homes/yl814/episcope_test/nutrient_stress/lighting_method_comparison_atac_benchmark"

required_syms <- c("grn_set", "ko_truth_list", "base_dir")
missing_syms <- required_syms[!vapply(required_syms, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
if (length(missing_syms)) {
  if (!file.exists(workspace_rdata)) {
    stop(
      "Missing required objects: ", paste(missing_syms, collapse = ", "),
      ". workspace_rdata not found: ", workspace_rdata
    )
  }
  message("Missing objects (", paste(missing_syms, collapse = ", "), "). Loading workspace: ", workspace_rdata)
  ws_env <- new.env(parent = emptyenv())
  load(workspace_rdata, envir = ws_env)
  for (sym in missing_syms) {
    if (exists(sym, envir = ws_env, inherits = FALSE)) {
      assign(sym, get(sym, envir = ws_env, inherits = FALSE), envir = .GlobalEnv)
    }
  }
  missing_syms <- required_syms[!vapply(required_syms, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
  if (length(missing_syms)) {
    stop("After workspace load, still missing required objects: ", paste(missing_syms, collapse = ", "))
  }
}

gr_n <- get("grn_set", inherits = TRUE)
ko_truth_list <- get("ko_truth_list", inherits = TRUE)
base_dir <- get("base_dir", inherits = TRUE)

if (!dir.exists(method_dir)) {
  stop("Method directory not found: ", method_dir)
}
plot_out_dir <- file.path(method_dir, "ko_plots_5methods")
dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)

fp_r_vals <- c(0.3, 0.5)
fp_p_vals <- c(0.1, 0.05, 0.01, 0.001)
param_grid <- expand.grid(fp_r = fp_r_vals, fp_p = fp_p_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

methods <- c("atac_pearson", "rna_pearson", "fp_pearson", "fp_spearman", "fp_rna_hybrid")
method_to_file <- c(
  atac_pearson = "atac_pearson",
  rna_pearson = "rna_pearson",
  fp_pearson = "pearson",
  fp_spearman = "spearman",
  fp_rna_hybrid = "fp_rna_hybrid"
)
method_labels <- c(
  atac_pearson = "ATAC Pearson",
  rna_pearson = "RNA Pearson",
  fp_pearson = "FP Pearson",
  fp_spearman = "FP Spearman",
  fp_rna_hybrid = "FP/RNA Hybrid"
)

.norm_chr <- function(x) toupper(trimws(as.character(x)))

group_fill <- c(Predicted = "#0F4C81", Random = "#C7CED6")
group_edge <- c(Predicted = "#0B355A", Random = "#7C8794")

.theme_pub <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(face = "bold", color = "#111111"),
      axis.text = ggplot2::element_text(face = "bold", color = "#111111"),
      axis.title = ggplot2::element_text(face = "bold", color = "#111111"),
      plot.title = ggplot2::element_text(face = "bold", color = "#111111", size = base_size + 1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E7EBF0", linewidth = 0.3),
      panel.border = ggplot2::element_rect(color = "#1A1A1A", fill = NA, linewidth = 0.4),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

connect_dir <- file.path(dirname(method_dir), "connect_tfs_to_target_genes")
fp_full_path <- file.path(connect_dir, "fp_gene_corr_full_jaspar2024.csv")
fp_gene_universe <- NULL
if (file.exists(fp_full_path)) {
  fp_gene_universe <- readr::read_csv(
    fp_full_path,
    col_select = c("tfs", "gene_key"),
    show_col_types = FALSE
  ) |>
    dplyr::mutate(tfs = .norm_chr(.data$tfs), gene_key = .norm_chr(.data$gene_key)) |>
    dplyr::filter(!is.na(.data$tfs), .data$tfs != "", !is.na(.data$gene_key), .data$gene_key != "") |>
    dplyr::distinct()
} else {
  warning("GeneHancer universe file not found: ", fp_full_path, ". vs_GeneHancer pool may be empty.")
}

.get_genehancer_pool <- function(tf) {
  if (is.null(fp_gene_universe) || !nrow(fp_gene_universe)) return(character(0))
  tf_norm <- .norm_chr(tf)
  unique(fp_gene_universe$gene_key[fp_gene_universe$tfs == tf_norm])
}

# KO TF -> cell line + baseline sample id mapping
# (restrict to requested TFs)
tf_ko_map <- tibble::tribble(
  ~tf,     ~cell,    ~sample_id,
  "HNF1A", "AsPC1", "cy414",
  "SOX9",  "Panc1", "cy423"
)

.detect_p_col <- function(ko_tbl) {
  for (nm in c("p_adj", "padj", "pval", "p_value")) {
    if (nm %in% names(ko_tbl)) return(nm)
  }
  NA_character_
}

.detect_basal_cols <- function(tbl) {
  tf_col <- intersect(c("TF", "tfs", "tf"), names(tbl))[1]
  gene_col <- intersect(c("gene_key", "gene", "target", "Name.Target"), names(tbl))[1]
  peak_col <- intersect(c("peak_ID", "fp_peak", "peak_id", "peak"), names(tbl))[1]
  if (anyNA(c(tf_col, gene_col, peak_col))) {
    stop("Basal table missing required columns. Found: ", paste(names(tbl), collapse = ", "))
  }
  list(tf = tf_col, gene = gene_col, peak = peak_col)
}

.detect_stat_cols <- function(tbl, method_key) {
  if (identical(method_key, "rna_pearson")) {
    r_col <- if ("r_rna_gene" %in% names(tbl)) "r_rna_gene" else NA_character_
    p_col <- if ("p_rna_gene" %in% names(tbl)) "p_rna_gene" else NA_character_
    p_adj_col <- if ("p_rna_adj_gene" %in% names(tbl)) "p_rna_adj_gene" else NA_character_
    return(list(p = p_col, r = r_col, p_adj = p_adj_col))
  }
  r_col <- if ("r_gene" %in% names(tbl)) "r_gene" else NA_character_
  p_col <- if ("p_gene" %in% names(tbl)) "p_gene" else NA_character_
  p_adj_col <- if ("p_adj_gene" %in% names(tbl)) "p_adj_gene" else NA_character_
  list(p = p_col, r = r_col, p_adj = p_adj_col)
}

.rna_expr_maps <- function(rna_expr_tbl, sample_id) {
  if (!sample_id %in% names(rna_expr_tbl)) stop("Sample id not found in rna_expressed: ", sample_id)
  expr_vec <- rna_expr_tbl[[sample_id]]
  hgnc_map <- expr_vec; names(hgnc_map) <- rna_expr_tbl$HGNC
  ens_map <- expr_vec; names(ens_map) <- rna_expr_tbl$ensembl_gene_id
  list(hgnc = hgnc_map, ens = ens_map)
}

.get_expr_flag <- function(gene_keys, expr_maps) {
  gene_norm <- .norm_chr(gene_keys)
  hgnc_norm <- .norm_chr(names(expr_maps$hgnc))
  ens_norm <- .norm_chr(names(expr_maps$ens))
  hgnc_idx <- match(gene_norm, hgnc_norm)
  out <- expr_maps$hgnc[hgnc_idx]
  miss <- is.na(out)
  if (any(miss)) {
    ens_idx <- match(gene_norm[miss], ens_norm)
    out[miss] <- expr_maps$ens[ens_idx]
  }
  ifelse(is.na(out), 0L, as.integer(out))
}

.pick_fp_annotation <- function(grn_set) {
  if (is.data.frame(grn_set$fp_annotation_pearson)) return(grn_set$fp_annotation_pearson)
  if (is.data.frame(grn_set$fp_annotation)) return(grn_set$fp_annotation)
  stop("No fp_annotation table available in grn_set.")
}

.gate_basal_for_tf <- function(basal_tbl, tf, sample_id, grn_set, fp_annotation_tbl, require_fp_bound = TRUE) {
  cols <- .detect_basal_cols(basal_tbl)
  basal_tbl <- basal_tbl |> dplyr::filter(.data[[cols$tf]] == tf)
  if (!nrow(basal_tbl)) return(basal_tbl)

  tf_norm <- .norm_chr(tf)
  ann_tf <- fp_annotation_tbl |>
    dplyr::mutate(tf_norm = .norm_chr(.data$tfs)) |>
    dplyr::filter(.data$tf_norm == tf_norm) |>
    dplyr::distinct(.data$fp_peak, .data$atac_peak)
  if (!nrow(ann_tf)) return(basal_tbl[0, , drop = FALSE])

  basal_tbl <- basal_tbl |>
    dplyr::mutate(fp_peak_norm = .norm_chr(.data[[cols$peak]])) |>
    dplyr::filter(.data$fp_peak_norm %in% .norm_chr(ann_tf$fp_peak))
  if (!nrow(basal_tbl)) return(basal_tbl[0, , drop = FALSE])

  expr_maps <- .rna_expr_maps(grn_set$rna_expressed, sample_id)
  tf_expr <- .get_expr_flag(tf, expr_maps)
  if (!isTRUE(tf_expr == 1L)) return(basal_tbl[0, , drop = FALSE])

  gene_expr <- .get_expr_flag(basal_tbl[[cols$gene]], expr_maps)

  if (isTRUE(require_fp_bound)) {
    fp_bound_tbl <- grn_set$fp_bound
    if (!sample_id %in% names(fp_bound_tbl)) stop("Sample id not found in fp_bound: ", sample_id)
    fp_idx <- match(basal_tbl[[cols$peak]], fp_bound_tbl[[1]])
    fp_bound_vec <- fp_bound_tbl[[sample_id]][fp_idx]
  } else {
    fp_bound_vec <- rep(1L, nrow(basal_tbl))
  }

  atac_tbl <- grn_set$atac_overlap
  if (!sample_id %in% names(atac_tbl)) stop("Sample id not found in atac_overlap: ", sample_id)
  fp_atac_idx <- match(basal_tbl[[cols$peak]], ann_tf$fp_peak)
  atac_peaks <- ann_tf$atac_peak[fp_atac_idx]
  atac_idx <- match(atac_peaks, atac_tbl[[1]])
  atac_vec <- atac_tbl[[sample_id]][atac_idx]

  keep <- (gene_expr == 1L) & (fp_bound_vec == 1L) & (atac_vec == 1L)
  basal_tbl[which(keep), , drop = FALSE]
}

.pick_best_peak_by_ko <- function(basal_tbl, ko_tbl, method_key) {
  if (!nrow(basal_tbl)) return(basal_tbl)
  cols <- .detect_basal_cols(basal_tbl)
  stat_cols <- .detect_stat_cols(basal_tbl, method_key)
  r_col <- stat_cols$r
  p_col <- stat_cols$p
  if (is.na(r_col) || !r_col %in% names(basal_tbl)) return(basal_tbl)
  if (is.na(p_col) || !p_col %in% names(basal_tbl)) return(basal_tbl)

  ko_map <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::group_by(.data$gene_norm) |>
    dplyr::summarise(ko_log2fc = mean(as.numeric(.data$log2fc), na.rm = TRUE), .groups = "drop")
  if (!nrow(ko_map)) return(basal_tbl)

  if (identical(method_key, "fp_rna_hybrid")) {
    # For hybrid only: choose peak per TF-gene by KO direction.
    # KO log2fc < 0 -> prefer higher r (activation-like), then smaller p.
    # KO log2fc > 0 -> prefer lower r (repression-like), then smaller p.
    return(
      basal_tbl |>
        dplyr::mutate(gene_norm = .norm_chr(.data[[cols$gene]])) |>
        dplyr::left_join(ko_map, by = "gene_norm") |>
        dplyr::mutate(
          r_val = as.numeric(.data[[r_col]]),
          p_val = as.numeric(.data[[p_col]]),
          ko_neg = is.finite(.data$ko_log2fc) & .data$ko_log2fc < 0,
          ko_pos = is.finite(.data$ko_log2fc) & .data$ko_log2fc > 0,
          ord1 = dplyr::case_when(
            ko_neg ~ -r_val,
            ko_pos ~ r_val,
            TRUE ~ abs(r_val)
          ),
          ord2 = p_val
        ) |>
        dplyr::mutate(
          ord1 = ifelse(is.finite(.data$ord1), .data$ord1, Inf),
          ord2 = ifelse(is.finite(.data$ord2), .data$ord2, Inf)
        ) |>
        dplyr::group_by(.data$gene_norm) |>
        dplyr::arrange(.data$ord1, .data$ord2, .by_group = TRUE) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::select(-dplyr::any_of(c("gene_norm", "ko_log2fc", "r_val", "p_val", "ko_neg", "ko_pos", "ord1", "ord2")))
    )
  }

  basal_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data[[cols$gene]])) |>
    dplyr::left_join(ko_map, by = "gene_norm") |>
    dplyr::mutate(
      r_val = as.numeric(.data[[r_col]]),
      p_val = as.numeric(.data[[p_col]]),
      ko_neg = is.finite(.data$ko_log2fc) & .data$ko_log2fc < 0,
      ko_pos = is.finite(.data$ko_log2fc) & .data$ko_log2fc > 0,
      ord1 = dplyr::case_when(
        ko_neg ~ -r_val,
        ko_pos ~ abs(r_val),
        TRUE ~ abs(r_val)
      ),
      ord2 = dplyr::case_when(
        ko_neg ~ p_val,
        ko_pos ~ -p_val,
        TRUE ~ p_val
      )
    ) |>
    dplyr::group_by(.data$gene_norm) |>
    dplyr::arrange(.data$ord1, .data$ord2, .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-dplyr::any_of(c("gene_norm", "ko_log2fc", "r_val", "p_val", "ko_neg", "ko_pos", "ord1", "ord2")))
}

.filter_fp_by_positive_rna <- function(basal_tbl, method_key) {
  if (!identical(method_key, "fp_pearson") && !identical(method_key, "fp_spearman")) {
    return(basal_tbl)
  }
  if (!"r_rna_gene" %in% names(basal_tbl)) {
    warning("r_rna_gene column not found for ", method_key, "; skipping RNA>0 filter.")
    return(basal_tbl)
  }
  basal_tbl |>
    dplyr::filter(is.finite(.data$r_rna_gene), .data$r_rna_gene > 0)
}

.pval_stars <- function(p) {
  if (!is.finite(p)) return("n.s.")
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  "n.s."
}

.make_normal_trimmed_reference <- function(
    ko_tbl,
    n_max,
    seed = 1L,
    clip_min = -3,
    clip_max = 3
) {
  if (!nrow(ko_tbl) || n_max <= 0L) return(character(0))
  pool_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene), log2fc = as.numeric(.data$log2fc)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::filter(.data$log2fc >= clip_min, .data$log2fc <= clip_max) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
  if (!nrow(pool_tbl)) return(character(0))

  n_pick <- min(as.integer(n_max), nrow(pool_tbl))
  if (n_pick <= 0L) return(character(0))

  set.seed(as.integer(seed))
  # deterministic tie-breaking
  pool_tbl <- pool_tbl |>
    dplyr::slice(sample.int(dplyr::n()))

  target_vals <- stats::qnorm(stats::ppoints(n_pick))
  target_vals <- pmax(clip_min, pmin(clip_max, target_vals))

  vals <- as.numeric(pool_tbl$log2fc)
  used <- rep(FALSE, length(vals))
  out_idx <- integer(n_pick)
  for (k in seq_len(n_pick)) {
    d <- abs(vals - target_vals[[k]])
    d[used] <- Inf
    take <- which.min(d)
    used[take] <- TRUE
    out_idx[[k]] <- take
  }
  pool_tbl$gene_norm[out_idx]
}

.make_stable_random_reference <- function(
    ko_tbl,
    baseline_mode = c("random", "genehancer"),
    genehancer_pool = NULL,
    n_max = 0L,
    seed = 1L,
    random_from_unchanged = TRUE
) {
  baseline_mode <- match.arg(baseline_mode)
  if (!nrow(ko_tbl) || n_max <= 0L) return(character(0))

  pool_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene), log2fc = as.numeric(.data$log2fc)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  pool_genes <- if (identical(baseline_mode, "genehancer")) {
    intersect(.norm_chr(genehancer_pool), pool_tbl$gene_norm)
  } else {
    if (isTRUE(random_from_unchanged) && "ko_group" %in% names(pool_tbl)) {
      pool_tbl$gene_norm[pool_tbl$ko_group == "Unchanged"]
    } else {
      pool_tbl$gene_norm
    }
  }

  pool_tbl <- pool_tbl |>
    dplyr::filter(.data$gene_norm %in% pool_genes) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
  if (!nrow(pool_tbl)) return(character(0))

  n_pick <- min(as.integer(n_max), nrow(pool_tbl))
  if (n_pick <= 0L) return(character(0))

  # Deterministic sign-balanced selection:
  # 1) prioritize pos/neg symmetry (mean near 0),
  # 2) preserve variance by spreading over abs(log2fc) quantiles,
  # 3) keep stable under fixed seed.
  set.seed(as.integer(seed))
  pool_tbl <- pool_tbl |>
    dplyr::mutate(
      row_id = dplyr::row_number(),
      abs_fc = abs(.data$log2fc)
    ) |>
    dplyr::slice(sample.int(dplyr::n()))  # deterministic tie-breaking

  pos_tbl <- pool_tbl |>
    dplyr::filter(.data$log2fc > 0) |>
    dplyr::arrange(.data$abs_fc)
  neg_tbl <- pool_tbl |>
    dplyr::filter(.data$log2fc < 0) |>
    dplyr::arrange(.data$abs_fc)
  zero_tbl <- pool_tbl |>
    dplyr::filter(.data$log2fc == 0) |>
    dplyr::arrange(.data$abs_fc)

  n_pos <- min(floor(n_pick / 2), nrow(pos_tbl))
  n_neg <- min(floor(n_pick / 2), nrow(neg_tbl))
  n_zero <- min(n_pick - n_pos - n_neg, nrow(zero_tbl))
  rem <- n_pick - n_pos - n_neg - n_zero
  if (rem > 0L) {
    free_pos <- nrow(pos_tbl) - n_pos
    free_neg <- nrow(neg_tbl) - n_neg
    add_pos <- min(rem, max(0L, free_pos))
    n_pos <- n_pos + add_pos
    rem <- rem - add_pos
    add_neg <- min(rem, max(0L, free_neg))
    n_neg <- n_neg + add_neg
    rem <- rem - add_neg
    if (rem > 0L) {
      n_zero <- n_zero + min(rem, nrow(zero_tbl) - n_zero)
    }
  }

  .pick_spread_ids <- function(tbl, n_take) {
    if (!nrow(tbl) || n_take <= 0L) return(integer(0))
    if (n_take >= nrow(tbl)) return(tbl$row_id)
    idx <- unique(pmax(1L, pmin(nrow(tbl), round(seq(1, nrow(tbl), length.out = n_take)))))
    # Ensure exact size with deterministic fill if rounding collapsed.
    if (length(idx) < n_take) {
      extra <- setdiff(seq_len(nrow(tbl)), idx)
      idx <- c(idx, head(extra, n_take - length(idx)))
    }
    tbl$row_id[idx]
  }

  take_ids <- c(
    .pick_spread_ids(pos_tbl, n_pos),
    .pick_spread_ids(neg_tbl, n_neg),
    .pick_spread_ids(zero_tbl, n_zero)
  )

  # Interleave signs to avoid ordering artifacts.
  out_tbl <- pool_tbl |>
    dplyr::filter(.data$row_id %in% take_ids) |>
    dplyr::mutate(
      sign_group = dplyr::case_when(
        .data$log2fc > 0 ~ "pos",
        .data$log2fc < 0 ~ "neg",
        TRUE ~ "zero"
      ),
      ord_in_sign = dplyr::row_number()
    ) |>
    dplyr::group_by(.data$sign_group) |>
    dplyr::arrange(.data$abs_fc, .by_group = TRUE) |>
    dplyr::mutate(ord_in_sign = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$ord_in_sign, .data$sign_group)

  out_tbl$gene_norm[seq_len(min(n_pick, nrow(out_tbl)))]
}

.plot_pred_vs_random <- function(
    ko_tbl,
    pred_genes,
    title = NULL,
    min_n = 5L,
    baseline_mode = c("random", "genehancer"),
    genehancer_pool = NULL,
    random_reference = NULL,
    random_from_unchanged = TRUE
) {
  baseline_mode <- match.arg(baseline_mode)
  if (!nrow(ko_tbl)) return(list(box = ggplot() + theme_void(), pct = ggplot() + theme_void()))
  baseline_label <- if (identical(baseline_mode, "genehancer")) "GeneHancer" else "Random"

  ko_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  pred_genes <- unique(.norm_chr(pred_genes))
  pred_genes <- pred_genes[pred_genes != "" & !is.na(pred_genes)]
  pred_in_ko <- intersect(pred_genes, ko_tbl$gene_norm)
  if (!length(pred_in_ko)) {
    empty_plot <- ggplot() + theme_void() + ggtitle(if (is.null(title)) "" else title)
    return(list(box = empty_plot, pct = empty_plot))
  }

  pool_tbl <- ko_tbl |>
    dplyr::mutate(log2fc = as.numeric(.data$log2fc))
  pool_genes <- if (identical(baseline_mode, "genehancer")) {
    intersect(.norm_chr(genehancer_pool), pool_tbl$gene_norm)
  } else {
    if (isTRUE(random_from_unchanged) && "ko_group" %in% names(pool_tbl)) {
      pool_tbl$gene_norm[pool_tbl$ko_group == "Unchanged"]
    } else {
      pool_tbl$gene_norm
    }
  }

  pool_tbl <- pool_tbl |>
    dplyr::filter(.data$gene_norm %in% pool_genes) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  # Cap random-selection candidate pool before final sampling.
  if (identical(baseline_mode, "random")) {
    pool_tbl <- pool_tbl |>
      dplyr::filter(is.finite(.data$log2fc), .data$log2fc >= -3, .data$log2fc <= 3)
  } else if (identical(baseline_mode, "genehancer")) {
    pred_abs_max <- suppressWarnings(max(abs(pool_tbl$log2fc[pool_tbl$gene_norm %in% pred_in_ko]), na.rm = TRUE))
    if (!is.finite(pred_abs_max) || pred_abs_max <= 0) pred_abs_max <- 3
    pool_tbl <- pool_tbl |>
      dplyr::filter(is.finite(.data$log2fc), .data$log2fc >= -pred_abs_max, .data$log2fc <= pred_abs_max)
  }

  n_pred <- length(pred_in_ko)
  if (length(random_reference) > 0L) {
    rand_genes <- unique(.norm_chr(random_reference))
    rand_genes <- rand_genes[rand_genes %in% pool_tbl$gene_norm]
    rand_genes <- head(rand_genes, n_pred)
    if (length(rand_genes) < n_pred) {
      # Deterministic back-fill from remaining pool so N_random tracks N_pred when pool permits.
      fill_genes <- setdiff(pool_tbl$gene_norm, rand_genes)
      rand_genes <- c(rand_genes, head(fill_genes, n_pred - length(rand_genes)))
    }
    rand_tbl <- pool_tbl |>
      dplyr::filter(.data$gene_norm %in% rand_genes) |>
      dplyr::mutate(.ord = match(.data$gene_norm, rand_genes)) |>
      dplyr::arrange(.data$.ord) |>
      dplyr::select(-".ord")
  } else {
    pool_tbl2 <- pool_tbl |>
      dplyr::filter(!.data$gene_norm %in% pred_in_ko)
    n_rand <- min(n_pred, nrow(pool_tbl2))
    set.seed(sum(utf8ToInt(paste(pred_in_ko, collapse = ""))) %% 2147483647L)
    rand_tbl <- if (n_rand > 0) dplyr::sample_n(pool_tbl2, n_rand) else pool_tbl2[0, , drop = FALSE]
  }
  rand_tbl <- dplyr::distinct(rand_tbl, .data$gene_norm, .keep_all = TRUE)

  df <- dplyr::bind_rows(
    tibble::tibble(gene_norm = pred_in_ko, group = "Predicted"),
    rand_tbl |> dplyr::transmute(gene_norm = .data$gene_norm, group = baseline_label)
  ) |>
    dplyr::left_join(ko_tbl, by = "gene_norm") |>
    dplyr::filter(is.finite(.data$log2fc))

  p_val <- NA_real_
  if (sum(df$group == "Predicted") > 1 && sum(df$group == baseline_label) > 1) {
    p_val <- suppressWarnings(stats::wilcox.test(log2fc ~ group, data = df)$p.value)
  }

  star <- .pval_stars(p_val)
  n_pred2 <- sum(df$group == "Predicted")
  n_rand2 <- sum(df$group == baseline_label)
  p_txt <- if (n_pred2 >= min_n && n_rand2 >= min_n) sprintf("p=%.2g %s", p_val, star) else "p=N/A"

  y_max <- max(df$log2fc, na.rm = TRUE)
  y_min <- min(df$log2fc, na.rm = TRUE)
  y_rng <- max(1, abs(y_max - y_min))
  y_annot <- y_max + 0.1 * y_rng
  y_n <- y_min + 0.02 * y_rng

  labs_bin <- c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")
  binned <- df |>
    dplyr::mutate(
      bin = dplyr::case_when(
        .data$log2fc <= -1 ~ "<= -1",
        .data$log2fc > -1 & .data$log2fc <= -0.5 ~ "(-1,-0.5]",
        .data$log2fc > -0.5 & .data$log2fc <= 0 ~ "(-0.5,0]",
        .data$log2fc > 0 ~ ">= 0",
        TRUE ~ NA_character_
      ),
      bin = factor(.data$bin, levels = labs_bin)
    ) |>
    dplyr::filter(!is.na(.data$bin)) |>
    dplyr::count(.data$group, .data$bin, name = "n") |>
    tidyr::complete(.data$group, .data$bin, fill = list(n = 0)) |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(frac = .data$n / sum(.data$n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      frac_signed = dplyr::case_when(
        .data$bin %in% c("<= -1", "(-1,-0.5]", "(-0.5,0]") ~ -.data$frac,
        TRUE ~ .data$frac
      ),
      bin = factor(.data$bin, levels = c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0"))
    )

  list(
    df = df,
    binned = binned,
    baseline_label = baseline_label,
    p_txt = p_txt,
    n_txt = sprintf("N=%d/%d", n_pred2, n_rand2),
    y_annot = y_annot,
    y_n = y_n,
    title = title
  )
}

.plot_from_panel_data <- function(panel_data) {
  if (is.null(panel_data) || !is.list(panel_data) || !is.data.frame(panel_data$df) || !nrow(panel_data$df)) {
    blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    return(list(box = blank_plot, pct = blank_plot))
  }

  df <- panel_data$df
  binned <- panel_data$binned
  baseline_label <- if (!is.null(panel_data$baseline_label)) panel_data$baseline_label else "Random"
  group_order <- c("Predicted", baseline_label)
  group_fill_vals <- c("Predicted" = unname(group_fill[["Predicted"]]), "Baseline" = unname(group_fill[["Random"]]))
  names(group_fill_vals)[2] <- baseline_label
  group_edge_vals <- c("Predicted" = unname(group_edge[["Predicted"]]), "Baseline" = unname(group_edge[["Random"]]))
  names(group_edge_vals)[2] <- baseline_label

  p_box <- ggplot(df, aes(x = .data$group, y = .data$log2fc, fill = .data$group)) +
    geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9, color = "#1A1A1A", linewidth = 0.35) +
    geom_hline(yintercept = 0, linetype = 2, color = "#5E6670", linewidth = 0.4) +
    geom_jitter(aes(color = .data$group), width = 0.11, height = 0, alpha = 0.18, size = 0.22, show.legend = FALSE) +
    annotate("text", x = 1.5, y = panel_data$y_annot, label = panel_data$p_txt, fontface = "bold", size = 3) +
    annotate("text", x = 1.5, y = panel_data$y_n, label = panel_data$n_txt, size = 2.8) +
    scale_fill_manual(name = "Group", values = group_fill_vals) +
    scale_color_manual(values = group_edge_vals) +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.14))) +
    labs(x = NULL, y = "KO log2FC", title = panel_data$title) +
    .theme_pub() +
    theme(legend.position = "right")

  p_pct <- ggplot(binned, aes(x = .data$group, y = .data$frac, fill = .data$bin)) +
    geom_col(width = 0.7, color = "#1A1A1A", linewidth = 0.2) +
    scale_fill_manual(
      name = "log2FC Bin",
      values = c("<= -1" = "#8E1B1B", "(-1,-0.5]" = "#C55A11", "(-0.5,0]" = "#E3A857", ">= 0" = "#D3D3D3")
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_x_discrete(limits = group_order) +
    labs(x = NULL, y = "Fraction", title = NULL) +
    .theme_pub() +
    theme(legend.position = "right")

  list(box = p_box, pct = p_pct)
}

.method_file_path <- function(method_key, fp_r, fp_p, method_dir) {
  file_stub <- if (identical(method_key, "rna_pearson")) {
    sprintf("basal_filtered_method_rna_pearson_rna_r%.1f_p%.2g.csv", fp_r, fp_p)
  } else if (identical(method_key, "atac_pearson")) {
    sprintf("basal_filtered_method_atac_pearson_r%.1f_p%.2g.csv", fp_r, fp_p)
  } else {
    sprintf("basal_filtered_method_%s_fp_r%.1f_p%.2g.csv", method_to_file[[method_key]], fp_r, fp_p)
  }
  file.path(method_dir, file_stub)
}

.load_method_tbl <- function(method_key, fp_r, fp_p, method_dir) {
  f <- .method_file_path(method_key, fp_r, fp_p, method_dir)
  if (!file.exists(f)) return(NULL)
  readr::read_csv(f, show_col_types = FALSE)
}

fp_annot_use <- .pick_fp_annotation(gr_n)
cache_dir <- file.path(plot_out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
use_panel_cache <- TRUE
force_rebuild_cache <- FALSE
stable_random_baseline <- TRUE
random_seed_base <- 20260210L
RUN_MAIN_5METHODS_BLOCK <- TRUE
RUN_ONLY_GENEHANCER_PDFS <- TRUE

if (isTRUE(RUN_MAIN_5METHODS_BLOCK)) {
for (i in seq_len(nrow(tf_ko_map))) {
  tf <- tf_ko_map$tf[[i]]
  cell <- tf_ko_map$cell[[i]]
  sample_id <- tf_ko_map$sample_id[[i]]

  ko_tbl <- ko_truth_list[[tf]]
  if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) {
    message("[skip] KO table missing/empty for TF: ", tf)
    next
  }

  param_labels <- sprintf("r=%.1f\np=%.2g", param_grid$fp_r, param_grid$fp_p)
  genehancer_pool_tf <- .get_genehancer_pool(tf)

  baseline_modes <- if (isTRUE(RUN_ONLY_GENEHANCER_PDFS)) c("genehancer") else c("random", "genehancer")
  for (baseline_mode in baseline_modes) {
    cache_file <- file.path(cache_dir, sprintf("%s_%s_5methods_panel_cache.rds", tf, baseline_mode))
    panel_cache <- NULL

    method_files <- unlist(lapply(methods, function(m) {
      vapply(seq_len(nrow(param_grid)), function(j) {
        .method_file_path(m, param_grid$fp_r[[j]], param_grid$fp_p[[j]], method_dir)
      }, character(1))
    }), use.names = FALSE)
    method_files <- unique(method_files[file.exists(method_files)])
    input_mtimes <- stats::setNames(as.character(file.info(method_files)$mtime), method_files)
    cache_signature <- list(
      input_mtimes = input_mtimes,
      stable_random_baseline = stable_random_baseline,
      random_seed_base = random_seed_base,
      stable_random_algorithm = "mean0_sd_total_symmetry_v5_main_random_normal_trimmed"
    )

    if (isTRUE(use_panel_cache) && file.exists(cache_file) && !isTRUE(force_rebuild_cache)) {
      tmp <- readRDS(cache_file)
      if (is.list(tmp) && is.list(tmp$meta) && is.list(tmp$panels)) {
        if (identical(tmp$meta$cache_signature, cache_signature)) {
          panel_cache <- tmp$panels
          message("Using plot cache: ", cache_file)
        } else {
          message("Input changed; rebuilding cache: ", cache_file)
        }
      }
    }

    if (is.null(panel_cache)) {
      panel_cache <- list()
      panel_inputs <- list()
      ko_tbl_norm <- ko_tbl |>
        dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
        dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
        dplyr::distinct(.data$gene_norm, .keep_all = TRUE)
      ko_gene_set <- ko_tbl_norm$gene_norm
      max_n_pred <- 0L
      n_pred_needed <- integer(0)

      for (method_key in methods) {
        panel_inputs[[method_key]] <- vector("list", nrow(param_grid))
        for (j in seq_len(nrow(param_grid))) {
          tbl <- .load_method_tbl(method_key, param_grid$fp_r[[j]], param_grid$fp_p[[j]], method_dir)
          if (is.null(tbl) || !nrow(tbl)) {
            panel_inputs[[method_key]][[j]] <- NULL
            next
          }

          req_fp <- !identical(method_key, "atac_pearson")
          tf_tbl <- .gate_basal_for_tf(tbl, tf, sample_id, gr_n, fp_annot_use, require_fp_bound = req_fp)
          tf_tbl <- .filter_fp_by_positive_rna(tf_tbl, method_key)
          tf_tbl <- .pick_best_peak_by_ko(tf_tbl, ko_tbl, method_key)
          if (!nrow(tf_tbl)) {
            panel_inputs[[method_key]][[j]] <- NULL
            next
          }

          cols <- .detect_basal_cols(tf_tbl)
          pred_genes <- unique(tf_tbl[[cols$gene]])
          pred_in_ko <- intersect(unique(.norm_chr(pred_genes)), ko_gene_set)
          max_n_pred <- max(max_n_pred, length(pred_in_ko))
          n_pred_needed <- c(n_pred_needed, length(pred_in_ko))

          panel_inputs[[method_key]][[j]] <- list(
            pred_genes = pred_genes,
            title = if (method_key == methods[[1]]) param_labels[[j]] else NULL,
            n_pred = length(pred_in_ko)
          )
        }
      }

      stable_seed <- as.integer(random_seed_base + i * 97L + if (identical(baseline_mode, "genehancer")) 13L else 7L)
      random_reference_map <- list()
      if (isTRUE(stable_random_baseline)) {
        n_pred_needed <- sort(unique(n_pred_needed[n_pred_needed > 0L]))
        for (n_req in n_pred_needed) {
          if (identical(baseline_mode, "random")) {
            random_reference_map[[as.character(n_req)]] <- .make_normal_trimmed_reference(
              ko_tbl = ko_tbl,
              n_max = n_req,
              seed = as.integer(stable_seed + n_req * 131L),
              clip_min = -3,
              clip_max = 3
            )
          } else {
            random_reference_map[[as.character(n_req)]] <- .make_stable_random_reference(
              ko_tbl = ko_tbl,
              baseline_mode = baseline_mode,
              genehancer_pool = genehancer_pool_tf,
              n_max = n_req,
              seed = as.integer(stable_seed + n_req * 131L)
            )
          }
        }
      }
      message(
        "Stable random baseline (", baseline_mode, "): max_n=", max_n_pred,
        ", n_levels=", length(random_reference_map),
        ", seed=", stable_seed
      )

      for (method_key in methods) {
        panel_cache[[method_key]] <- vector("list", nrow(param_grid))
        for (j in seq_len(nrow(param_grid))) {
          panel_in <- panel_inputs[[method_key]][[j]]
          if (is.null(panel_in)) {
            panel_cache[[method_key]][[j]] <- NULL
            next
          }

          panel_cache[[method_key]][[j]] <- .plot_pred_vs_random(
            ko_tbl = ko_tbl,
            pred_genes = panel_in$pred_genes,
            title = panel_in$title,
            min_n = 5L,
            baseline_mode = baseline_mode,
            genehancer_pool = genehancer_pool_tf,
            random_reference = random_reference_map[[as.character(panel_in$n_pred)]],
            random_from_unchanged = !identical(baseline_mode, "random")
          )
        }
      }
      saveRDS(list(meta = list(cache_signature = cache_signature), panels = panel_cache), cache_file)
      message("Saved plot cache: ", cache_file)
    }

    row_plots <- list()
    blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
    for (method_key in methods) {
      row1 <- vector("list", nrow(param_grid))
      row2 <- vector("list", nrow(param_grid))
      for (j in seq_len(nrow(param_grid))) {
        panel_data <- panel_cache[[method_key]][[j]]
        if (is.null(panel_data)) {
          row1[[j]] <- blank_plot
          row2[[j]] <- blank_plot
        } else {
          g <- .plot_from_panel_data(panel_data)
          row1[[j]] <- g$box
          row2[[j]] <- g$pct
        }
      }

      label_base <- method_labels[[method_key]]
      if (length(row1)) {
        row1[[1]] <- row1[[1]] + ggplot2::labs(y = paste(label_base, "violin"))
        row1[[1]] <- row1[[1]] + ggplot2::theme(axis.title.y = ggplot2::element_text(face = "bold", size = 9))
      }
      if (length(row2)) {
        row2[[1]] <- row2[[1]] + ggplot2::labs(y = paste(label_base, "bins"))
        row2[[1]] <- row2[[1]] + ggplot2::theme(axis.title.y = ggplot2::element_text(face = "bold", size = 9))
      }
      row_plots[[length(row_plots) + 1L]] <- patchwork::wrap_plots(row1, nrow = 1, ncol = nrow(param_grid))
      row_plots[[length(row_plots) + 1L]] <- patchwork::wrap_plots(row2, nrow = 1, ncol = nrow(param_grid))
    }

    plot_stack <- row_plots[[1]]
    for (k in 2:length(row_plots)) plot_stack <- plot_stack / row_plots[[k]]

    comp_tag <- if (identical(baseline_mode, "genehancer")) "GeneHancer" else "random"
    plot_stack <- plot_stack +
      patchwork::plot_layout(guides = "collect", heights = rep(c(1.05, 1.0), times = length(methods))) +
      patchwork::plot_annotation(
        title = sprintf("%s (%s) | 5-method KO benchmark vs %s", tf, cell, comp_tag),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15, color = "#111111"))
      )

    out_pdf <- if (identical(baseline_mode, "genehancer")) {
      file.path(plot_out_dir, sprintf("%s_%s_5methods_vs_GeneHancer.pdf", tf, cell))
    } else {
      file.path(plot_out_dir, sprintf("%s_%s_5methods_vs_random.pdf", tf, cell))
    }
    ggplot2::ggsave(out_pdf, plot_stack, width = 20, height = 32, units = "in", dpi = 320, bg = "white")
    message("Saved: ", out_pdf)
  }
}
}

# Optional standalone custom gene-set plotting block ---------------------------
# This block does not affect the default pipeline unless explicitly enabled.
RUN_CUSTOM_GENESET_BLOCK <- FALSE

if (isTRUE(RUN_CUSTOM_GENESET_BLOCK)) {
  custom_out_dir <- file.path(plot_out_dir, "custom_gene_sets")
  dir.create(custom_out_dir, recursive = TRUE, showWarnings = FALSE)
  custom_random_mode <- "normal_trimmed"  # options: normal_trimmed, stable

  .plot_custom_one <- function(tf, cell, method_name, genes) {
    ko_tbl <- ko_truth_list[[tf]]
    if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) {
      message("[custom][skip] KO table missing/empty for TF: ", tf)
      return(invisible(NULL))
    }
    if (is.null(genes) || !length(genes)) {
      message("[custom][skip] Empty gene list for TF: ", tf, " method: ", method_name)
      return(invisible(NULL))
    }

    genes <- unique(.norm_chr(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    if (!length(genes)) {
      message("[custom][skip] No valid genes after normalization for TF: ", tf, " method: ", method_name)
      return(invisible(NULL))
    }

    ko_tbl_norm <- ko_tbl |>
      dplyr::mutate(
        gene_norm = .norm_chr(.data$gene),
        ens_norm = if ("ensembl_gene_id" %in% names(ko_tbl)) .norm_chr(.data$ensembl_gene_id) else NA_character_
      )
    keep_idx <- ko_tbl_norm$gene_norm %in% genes | ko_tbl_norm$ens_norm %in% genes
    pred_genes_plot <- unique(ko_tbl_norm$gene_norm[keep_idx])
    pred_genes_plot <- pred_genes_plot[!is.na(pred_genes_plot) & pred_genes_plot != ""]
    pred_in_ko <- pred_genes_plot
    if (!length(pred_in_ko)) {
      message("[custom][skip] No overlap with KO genes for TF: ", tf, " method: ", method_name)
      return(invisible(NULL))
    }

    stable_seed <- as.integer(random_seed_base + sum(utf8ToInt(tf)) + sum(utf8ToInt(method_name)) + 17L)
    random_reference <- if (isTRUE(stable_random_baseline) && identical(custom_random_mode, "normal_trimmed")) {
      .make_normal_trimmed_reference(
        ko_tbl = ko_tbl,
        n_max = length(pred_in_ko),
        seed = stable_seed,
        clip_min = -3,
        clip_max = 3
      )
    } else if (isTRUE(stable_random_baseline)) {
      .make_stable_random_reference(
        ko_tbl = ko_tbl,
        baseline_mode = "random",
        genehancer_pool = NULL,
        n_max = length(pred_in_ko),
        seed = stable_seed,
        random_from_unchanged = FALSE
      )
    } else {
      character(0)
    }
    message(
      "[custom] random_mode=", custom_random_mode,
      " | n_pred=", length(pred_in_ko),
      " | n_rand_ref=", length(random_reference)
    )

    panel_data <- .plot_pred_vs_random(
      ko_tbl = ko_tbl,
      pred_genes = pred_genes_plot,
      title = paste0(method_name, "\n", tf, " (", cell, ")"),
      min_n = 5L,
      baseline_mode = "random",
      genehancer_pool = NULL,
      random_reference = random_reference,
      random_from_unchanged = FALSE
    )
    g <- .plot_from_panel_data(panel_data)
    combo <- g$box / g$pct +
      patchwork::plot_layout(heights = c(1.2, 1.0), guides = "collect") +
      patchwork::plot_annotation(
        title = sprintf("%s (%s) | %s vs random", tf, cell, method_name),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, color = "#111111"),
          legend.position = "right"
        )
      )

    safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", method_name)
    out_pdf <- file.path(custom_out_dir, sprintf("%s_%s_%s_vs_random.pdf", tf, cell, safe_name))
    ggplot2::ggsave(out_pdf, combo, width = 8.5, height = 10, units = "in", dpi = 320, bg = "white")
    message("[custom] Saved: ", out_pdf)
    invisible(out_pdf)
  }

  granie_csv <- file.path(plot_out_dir, "granie_conns_full.csv")
  if (!file.exists(granie_csv)) {
    stop("Custom mode requires file not found: ", granie_csv)
  }

  granie_tbl <- readr::read_csv(granie_csv, show_col_types = FALSE)
  names(granie_tbl) <- trimws(names(granie_tbl))
  required_cols <- c("TF.ID", "gene.ENSEMBL", "gene.name")
  miss_cols <- setdiff(required_cols, names(granie_tbl))
  if (length(miss_cols)) {
    stop("Missing required columns in granie_conns_full.csv: ", paste(miss_cols, collapse = ", "))
  }

  granie_tf_id <- "HNF1A.H12INVIVO.0.PS.A"
  granie_sub <- granie_tbl |>
    dplyr::filter(.data$`TF.ID` == granie_tf_id)
  if (!nrow(granie_sub)) {
    stop("No rows found for TF.ID = ", granie_tf_id, " in ", granie_csv)
  }

  granie_genes <- unique(c(granie_sub$`gene.ENSEMBL`, granie_sub$`gene.name`))
  granie_genes <- granie_genes[!is.na(granie_genes) & trimws(granie_genes) != ""]
  .plot_custom_one(
    tf = "HNF1A",
    cell = "AsPC1",
    method_name = "granie_method",
    genes = granie_genes
  )
}
