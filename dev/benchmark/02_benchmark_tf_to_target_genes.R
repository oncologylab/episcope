# Nutrient method-comparison plotting (refactored)
# Assumes these objects exist in the environment:
#   - grn_set (with rna_expressed, fp_bound, atac_overlap, fp_annotation or fp_annotation_* )
#   - ko_truth_list (named list of KO tables per TF; columns: gene, log2fc, ko_group, optional pval)
#   - base_dir
# If any required object is missing, this script attempts to auto-load a workspace
# from EPISCOPE_TF_TARGET_WORKSPACE_RDATA, defaulting to:
# /data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/workspace.RData
#
# Standalone usage:
# - Set env flags as needed (EPISCOPE_SAVE_MAIN_PDF, EPISCOPE_MAIN_PDF_ONLY,
#   EPISCOPE_SAVE_UNION_CSV, EPISCOPE_SAVE_VENN_PDF).
# - Venn-specific block is controlled by RUN_VENN_BLOCK below.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

required_syms <- c("grn_set", "ko_truth_list", "base_dir")
missing_syms <- required_syms[!vapply(required_syms, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
workspace_rdata <- Sys.getenv(
  "EPISCOPE_TF_TARGET_WORKSPACE_RDATA",
  unset = "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/workspace.RData"
)
if (file.exists(workspace_rdata)) {
  workspace_env <- new.env(parent = emptyenv())
  load(workspace_rdata, envir = workspace_env)

  # Force-refresh base_dir from workspace regardless of whether it already exists.
  if ("base_dir" %in% ls(workspace_env)) {
    assign("base_dir", get("base_dir", envir = workspace_env), envir = .GlobalEnv)
    message("Force-loaded base_dir from workspace: ", workspace_rdata)
  } else {
    warning("Workspace loaded but has no 'base_dir': ", workspace_rdata)
  }

  if (length(missing_syms)) {
    message(
      "Missing objects (", paste(missing_syms, collapse = ", "),
      "). Loading from workspace: ", workspace_rdata
    )
    for (sym in missing_syms) {
      if (sym %in% ls(workspace_env)) {
        assign(sym, get(sym, envir = workspace_env), envir = .GlobalEnv)
      }
    }
  }
} else if (length(missing_syms)) {
  stop(
    "Missing required objects in environment: ", paste(missing_syms, collapse = ", "),
    ". Also could not find workspace file: ", workspace_rdata
  )
}

missing_syms <- required_syms[!vapply(required_syms, exists, logical(1), envir = .GlobalEnv, inherits = TRUE)]
if (length(missing_syms)) {
  stop(
    "After workspace load, still missing required objects: ",
    paste(missing_syms, collapse = ", ")
  )
}

gr_n <- get("grn_set", inherits = TRUE)
ko_truth_list <- get("ko_truth_list", inherits = TRUE)
base_dir <- get("base_dir", inherits = TRUE)

RUN_VENN_BLOCK <- FALSE

plot_out_dir <- if (exists("ko_dir", inherits = TRUE)) {
  get("ko_dir", inherits = TRUE)
} else {
  file.path(base_dir, "lighting_method_comparison", "ko_plots")
}
dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)

save_main_pdf <- tolower(Sys.getenv("EPISCOPE_SAVE_MAIN_PDF", unset = "false")) %in% c("1", "true", "t", "yes", "y")
save_venn_pdf <- tolower(Sys.getenv("EPISCOPE_SAVE_VENN_PDF", unset = "false")) %in% c("1", "true", "t", "yes", "y")
use_existing_venn_csv <- tolower(Sys.getenv("EPISCOPE_USE_EXISTING_VENN_CSV", unset = "true")) %in% c("1", "true", "t", "yes", "y")
main_pdf_only <- tolower(Sys.getenv("EPISCOPE_MAIN_PDF_ONLY", unset = "false")) %in% c("1", "true", "t", "yes", "y")
include_kendall <- tolower(Sys.getenv("EPISCOPE_INCLUDE_KENDALL", unset = "false")) %in% c("1", "true", "t", "yes", "y")
save_union_csv <- tolower(Sys.getenv("EPISCOPE_SAVE_UNION_CSV", unset = "true")) %in% c("1", "true", "t", "yes", "y")
if (isTRUE(main_pdf_only)) {
  save_main_pdf <- TRUE
  save_venn_pdf <- FALSE
  use_existing_venn_csv <- FALSE
}
if (isTRUE(save_union_csv)) {
  save_main_pdf <- FALSE
  save_venn_pdf <- FALSE
  use_existing_venn_csv <- FALSE
  main_pdf_only <- TRUE
}

method_dir <- file.path(base_dir, "lighting_method_comparison")
if (!dir.exists(method_dir)) {
  stop("Missing method directory: ", method_dir)
}

fp_r_vals <- c(0.3, 0.5)
fp_p_vals <- c(0.1, 0.05, 0.01, 0.001)
param_grid <- expand.grid(fp_r = fp_r_vals, fp_p = fp_p_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

methods <- c("fp_pearson", "fp_spearman", "rna_pearson")
if (isTRUE(include_kendall)) methods <- c("fp_pearson", "fp_spearman", "fp_kendall", "rna_pearson")
method_to_file <- c(
  fp_pearson = "pearson",
  fp_spearman = "spearman",
  rna_pearson = "rna_pearson"
)
if (isTRUE(include_kendall)) method_to_file <- c(method_to_file, fp_kendall = "kendall")
method_labels <- c(
  fp_pearson = "FP Pearson",
  fp_spearman = "FP Spearman",
  rna_pearson = "RNA Pearson"
)
if (isTRUE(include_kendall)) method_labels <- c(method_labels, fp_kendall = "FP Kendall")

methods <- methods[
  vapply(
    methods,
    function(m) {
      any(vapply(seq_len(nrow(param_grid)), function(j) {
        f <- if (identical(m, "rna_pearson")) {
          file.path(
            method_dir,
            sprintf(
              "basal_filtered_method_rna_pearson_rna_r%.1f_p%.2g.csv",
              param_grid$fp_r[[j]],
              param_grid$fp_p[[j]]
            )
          )
        } else {
          file.path(
            method_dir,
            sprintf(
              "basal_filtered_method_%s_fp_r%.1f_p%.2g.csv",
              method_to_file[[m]],
              param_grid$fp_r[[j]],
              param_grid$fp_p[[j]]
            )
          )
        }
        file.exists(f)
      }, logical(1)))
    },
    logical(1)
  )
]
if (!length(methods)) stop("No basal_filtered_method_* files found in ", method_dir)

# KO TF -> cell line + baseline sample id mapping
# (cell name for filename; sample id for gating)
tf_ko_map <- tibble::tribble(
  ~tf,     ~cell,      ~sample_id,
  "HNF1A", "AsPC1",    "cy414",
  "SOX9",  "Panc1",    "cy423",
  "HNF4A", "Mayo5289", "cy333",
  "IRF1",  "CFPAC1",   "cy333",
  "RARG",  "Panc1",    "cy423",
  "KLF5",  "AsPC1",    "cy414",
  "FOXA2", "Panc1",    "cy423"
)

.norm_chr <- function(x) toupper(trimws(as.character(x)))
.or_empty <- function(x) if (is.null(x)) character(0) else x
.parse_genes <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  out <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  out <- out[!is.na(out) & out != ""]
  unique(out)
}

.pval_stars <- function(p) {
  if (!is.finite(p)) return("n.s.")
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  "n.s."
}

.detect_p_col <- function(ko_tbl) {
  for (nm in c("p_adj", "padj", "pval", "p_value")) {
    if (nm %in% names(ko_tbl)) return(nm)
  }
  NA_character_
}

.venn_weights_4 <- function(sets, set_names) {
  stopifnot(length(sets) == 4, length(set_names) == 4)
  names(sets) <- set_names

  region_weight <- function(in_sets) {
    inter <- Reduce(intersect, sets[in_sets])
    out_sets <- setdiff(set_names, in_sets)
    if (length(out_sets)) {
      inter <- setdiff(inter, unique(unlist(sets[out_sets])))
    }
    length(unique(inter))
  }

  order_keys <- list(
    "1000" = set_names[1],
    "0100" = set_names[2],
    "1100" = set_names[1:2],
    "0010" = set_names[3],
    "1010" = set_names[c(1, 3)],
    "0110" = set_names[c(2, 3)],
    "1110" = set_names[c(1, 2, 3)],
    "0001" = set_names[4],
    "1001" = set_names[c(1, 4)],
    "0101" = set_names[c(2, 4)],
    "1101" = set_names[c(1, 2, 4)],
    "0011" = set_names[c(3, 4)],
    "1011" = set_names[c(1, 3, 4)],
    "0111" = set_names[c(2, 3, 4)],
    "1111" = set_names
  )

  vapply(order_keys, region_weight, integer(1))
}

.simple_venn_plot <- function(sets, title = NULL) {
  a <- unique(.norm_chr(sets[[1]]))
  b <- unique(.norm_chr(sets[[2]]))
  c <- unique(.norm_chr(sets[[3]]))
  ab <- intersect(a, b)
  ac <- intersect(a, c)
  bc <- intersect(b, c)
  abc <- Reduce(intersect, list(a, b, c))
  a_only <- setdiff(a, union(b, c))
  b_only <- setdiff(b, union(a, c))
  c_only <- setdiff(c, union(a, b))
  ab_only <- setdiff(ab, c)
  ac_only <- setdiff(ac, b)
  bc_only <- setdiff(bc, a)

  circle_df <- function(cx, cy, r, n = 200L, id = 1L) {
    t <- seq(0, 2 * pi, length.out = n)
    data.frame(
      x = cx + r * cos(t),
      y = cy + r * sin(t),
      id = id
    )
  }
  circles <- rbind(
    circle_df(0, 0, 1, id = 1L),
    circle_df(1.2, 0, 1, id = 2L),
    circle_df(0.6, 1.0, 1, id = 3L)
  )
  circles$fill <- factor(circles$id, levels = c(1, 2, 3),
                         labels = c("fp_pearson", "fp_spearman", "fp_kendall"))

  p <- ggplot2::ggplot(circles, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id, fill = .data$fill)) +
    ggplot2::geom_polygon(color = "black", alpha = 0.35) +
    ggplot2::annotate("text", x = -0.35, y = 0, label = length(a_only), size = 3) +
    ggplot2::annotate("text", x = 1.55, y = 0, label = length(b_only), size = 3) +
    ggplot2::annotate("text", x = 0.6, y = 1.55, label = length(c_only), size = 3) +
    ggplot2::annotate("text", x = 0.6, y = -0.15, label = length(ab_only), size = 3) +
    ggplot2::annotate("text", x = 0.1, y = 0.7, label = length(ac_only), size = 3) +
    ggplot2::annotate("text", x = 1.1, y = 0.7, label = length(bc_only), size = 3) +
    ggplot2::annotate("text", x = 0.6, y = 0.35, label = length(abc), size = 3, fontface = "bold") +
    ggplot2::annotate("text", x = -0.65, y = -0.85, label = "fp_pearson", size = 2.7) +
    ggplot2::annotate("text", x = 1.85, y = -0.85, label = "fp_spearman", size = 2.7) +
    ggplot2::annotate("text", x = 0.6, y = 2.05, label = "fp_kendall", size = 2.7) +
    ggplot2::scale_fill_manual(values = c(
      fp_pearson = "#4C78A8",
      fp_spearman = "#F58518",
      fp_kendall = "#54A24B"
    )) +
    ggplot2::coord_fixed(xlim = c(-1.2, 2.4), ylim = c(-1.2, 2.4)) +
    ggplot2::theme_void()

  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

.pick_fp_annotation <- function(grn_set, method_key = "pearson") {
  if (identical(method_key, "spearman") && is.data.frame(grn_set$fp_annotation_spearman)) {
    return(grn_set$fp_annotation_spearman)
  }
  if (identical(method_key, "kendall") && is.data.frame(grn_set$fp_annotation_kendall)) {
    return(grn_set$fp_annotation_kendall)
  }
  if (is.data.frame(grn_set$fp_annotation_pearson)) return(grn_set$fp_annotation_pearson)
  if (is.data.frame(grn_set$fp_annotation)) return(grn_set$fp_annotation)
  stop("No fp_annotation table available in grn_set.")
}

fp_annot_by_method <- list(
  fp_pearson  = .pick_fp_annotation(gr_n, "pearson"),
  fp_spearman = .pick_fp_annotation(gr_n, "spearman"),
  rna_pearson = .pick_fp_annotation(gr_n, "pearson")
)
if (isTRUE(include_kendall)) {
  fp_annot_by_method$fp_kendall <- .pick_fp_annotation(gr_n, "kendall")
}

.rna_expr_maps <- function(rna_expr_tbl, sample_id) {
  if (!sample_id %in% names(rna_expr_tbl)) {
    stop("Sample id not found in rna_expressed: ", sample_id)
  }
  expr_vec <- rna_expr_tbl[[sample_id]]
  hgnc_map <- expr_vec
  names(hgnc_map) <- rna_expr_tbl$HGNC
  ens_map <- expr_vec
  names(ens_map) <- rna_expr_tbl$ensembl_gene_id
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
  out <- ifelse(is.na(out), 0L, as.integer(out))
  out
}

.load_basal <- function(method_key, fp_r, fp_p) {
  file_stub <- if (identical(method_key, "rna_pearson")) {
    sprintf("basal_filtered_method_rna_pearson_rna_r%.1f_p%.2g.csv", fp_r, fp_p)
  } else {
    sprintf("basal_filtered_method_%s_fp_r%.1f_p%.2g.csv", method_to_file[[method_key]], fp_r, fp_p)
  }
  f <- file.path(method_dir, file_stub)
  if (!file.exists(f)) {
    message("[skip] missing file: ", f)
    return(NULL)
  }
  readr::read_csv(f, show_col_types = FALSE)
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

.gate_basal_for_tf <- function(basal_tbl, tf, sample_id, grn_set, fp_annotation_tbl) {
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

  gene_keys <- basal_tbl[[cols$gene]]
  gene_expr <- .get_expr_flag(gene_keys, expr_maps)

  fp_bound_tbl <- grn_set$fp_bound
  if (!sample_id %in% names(fp_bound_tbl)) stop("Sample id not found in fp_bound: ", sample_id)
  fp_idx <- match(basal_tbl[[cols$peak]], fp_bound_tbl[[1]])
  fp_bound_vec <- fp_bound_tbl[[sample_id]][fp_idx]

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
    dplyr::select(-c(.data$gene_norm, .data$ko_log2fc, .data$r_val, .data$p_val, .data$ko_neg, .data$ko_pos, .data$ord1, .data$ord2))
}

.plot_pred_vs_random <- function(ko_tbl, pred_genes, random_pool = NULL, title = NULL, min_n = 5L) {
  if (!nrow(ko_tbl)) return(list(box = ggplot() + theme_void(), pct = ggplot() + theme_void(), stats = NULL))

  ko_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  pred_genes <- unique(.norm_chr(pred_genes))
  pred_genes <- pred_genes[pred_genes != "" & !is.na(pred_genes)]

  pred_in_ko <- intersect(pred_genes, ko_tbl$gene_norm)
  if (!length(pred_in_ko)) {
    empty_plot <- ggplot() + theme_void() + ggtitle(if (is.null(title)) "" else title)
    return(list(box = empty_plot, pct = empty_plot, stats = NULL))
  }

  pool_tbl <- ko_tbl |>
    dplyr::filter(is.finite(.data$log2fc)) |>
    dplyr::mutate(log2fc = as.numeric(.data$log2fc))
  pool_genes <- if ("ko_group" %in% names(pool_tbl)) {
    pool_tbl$gene_norm[pool_tbl$ko_group == "Unchanged"]
  } else {
    pool_tbl$gene_norm
  }
  if (!is.null(random_pool)) {
    pool_genes <- intersect(.norm_chr(random_pool), pool_genes)
  }
  pool_tbl <- pool_tbl |>
    dplyr::filter(.data$gene_norm %in% pool_genes)
  pool_tbl <- pool_tbl |>
    dplyr::filter(.data$gene_norm %in% setdiff(.data$gene_norm, pred_in_ko))

  n_pred <- length(pred_in_ko)
  n_rand <- min(n_pred, nrow(pool_tbl))
  set.seed(sum(utf8ToInt(paste(pred_in_ko, collapse = ""))) %% 2147483647L)

  pool_pos <- pool_tbl |> dplyr::filter(.data$log2fc >= 0)
  pool_neg <- pool_tbl |> dplyr::filter(.data$log2fc < 0)
  n_pos <- min(nrow(pool_pos), ceiling(n_rand / 2))
  n_neg <- min(nrow(pool_neg), n_rand - n_pos)

  pick_pos <- if (n_pos > 0) pool_pos |> dplyr::sample_n(n_pos) else pool_pos[0, , drop = FALSE]
  pick_neg <- if (n_neg > 0) pool_neg |> dplyr::sample_n(n_neg) else pool_neg[0, , drop = FALSE]
  rand_tbl <- dplyr::bind_rows(pick_pos, pick_neg)
  if (nrow(rand_tbl) < n_rand && nrow(pool_tbl) > nrow(rand_tbl)) {
    need <- n_rand - nrow(rand_tbl)
    rand_tbl <- dplyr::bind_rows(rand_tbl, pool_tbl |> dplyr::sample_n(need))
  }
  rand_tbl <- rand_tbl |> dplyr::distinct(.data$gene_norm, .keep_all = TRUE)

  df <- dplyr::bind_rows(
    tibble::tibble(gene_norm = pred_in_ko, group = "Predicted"),
    rand_tbl |> dplyr::transmute(gene_norm = .data$gene_norm, group = "Random")
  ) |>
    dplyr::left_join(ko_tbl, by = "gene_norm") |>
    dplyr::filter(is.finite(.data$log2fc))

  p_val <- NA_real_
  if (sum(df$group == "Predicted") > 1 && sum(df$group == "Random") > 1) {
    p_val <- suppressWarnings(stats::wilcox.test(log2fc ~ group, data = df)$p.value)
  }
  star <- .pval_stars(p_val)
  n_pred <- sum(df$group == "Predicted")
  n_rand <- sum(df$group == "Random")
  p_txt <- if (n_pred >= min_n && n_rand >= min_n) sprintf("p=%.2g %s", p_val, star) else "p=N/A"
  n_txt <- sprintf("N=%d/%d", n_pred, n_rand)

  y_max <- max(df$log2fc, na.rm = TRUE)
  y_min <- min(df$log2fc, na.rm = TRUE)
  y_rng <- max(1, abs(y_max - y_min))
  y_annot <- y_max + 0.1 * y_rng
  y_n <- y_min + 0.02 * y_rng

  p_box <- ggplot(df, aes(x = .data$group, y = .data$log2fc, fill = .data$group)) +
    geom_violin(trim = FALSE, alpha = 0.5, color = "black") +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "black") +
    scale_fill_manual(values = c(Predicted = "#4C78A8", Random = "#A0A0A0")) +
    labs(title = title, x = NULL, y = "log2FC") +
    annotate("text", x = 1.5, y = y_annot, label = p_txt, size = 3.0, fontface = "bold") +
    annotate("text", x = 0.8, y = y_n, label = n_txt, size = 3.0, fontface = "bold", hjust = 0) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      axis.text = element_text(size = 8, face = "bold"),
      axis.title = element_text(size = 9, face = "bold"),
      plot.margin = ggplot2::margin(t = 6, r = 4, b = 0, l = 4)
    )

  lfc_breaks <- c(-Inf, -1, -0.5, 0, Inf)
  lfc_labels <- c("<= -1", "(-1,-0.5]", "(-0.5,0]", ">= 0")
  df$log2fc_bin <- cut(df$log2fc, breaks = lfc_breaks, labels = lfc_labels, right = TRUE, include.lowest = TRUE)

  if (n_pred < min_n || n_rand < min_n) {
    p_pct <- ggplot() + theme_void()
  } else {
    df_pct <- df |>
      dplyr::count(.data$group, .data$log2fc_bin, name = "n") |>
      dplyr::group_by(.data$group) |>
      dplyr::mutate(pct = 100 * .data$n / sum(.data$n)) |>
      dplyr::ungroup() |>
      tidyr::complete(group = c("Predicted", "Random"), log2fc_bin = lfc_labels, fill = list(n = 0L, pct = 0))

    p_pct <- ggplot(df_pct, aes(x = .data$group, y = .data$pct, fill = .data$log2fc_bin)) +
      geom_col(width = 0.6, color = "grey70", linewidth = 0.2) +
      scale_fill_manual(values = c("<= -1" = "#d73027", "(-1,-0.5]" = "#fc8d59", "(-0.5,0]" = "#fee090", ">= 0" = "#fbfbfb")) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
      labs(x = NULL, y = "% genes") +
      theme_bw(base_size = 9) +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 7, face = "bold"),
        axis.title = element_text(size = 8, face = "bold"),
        plot.margin = ggplot2::margin(t = 0, r = 4, b = 2, l = 4)
      )
  }

  list(box = p_box, pct = p_pct, stats = tibble::tibble(p_val = p_val, n_pred = n_pred, n_rand = n_rand))
}

.summarize_method_overlap <- function(tf, ko_tbl, callsets_by_method, effect_thr = -0.5, overlap_high = 0.6, overlap_mid = 0.4) {
  tf0 <- .norm_chr(tf)
  ko_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc))

  ko_down <- ko_tbl$gene_norm[ko_tbl$log2fc < effect_thr]
  ko_down <- unique(ko_down)

  sets <- lapply(callsets_by_method, function(x) intersect(.norm_chr(x), ko_down))
  fp_sets <- sets[intersect(c("fp_pearson", "fp_spearman", "fp_kendall"), names(sets))]
  fp_union <- unique(unlist(fp_sets))
  fp_inter <- Reduce(intersect, fp_sets)

  fp2_counts <- table(unlist(fp_sets))
  fp_2of3 <- if (length(fp_sets) >= 3) names(fp2_counts)[fp2_counts >= 2] else character(0)

  fp_union_n <- length(fp_union)
  fp_inter_n <- length(fp_inter)
  fp_2of3_n <- length(fp_2of3)

  fp_3_frac <- if (fp_union_n > 0) fp_inter_n / fp_union_n else NA_real_
  fp_2_frac <- if (fp_union_n > 0) fp_2of3_n / fp_union_n else NA_real_

  strategy <- if (is.finite(fp_3_frac) && fp_3_frac >= overlap_high) {
    "fp_3method_consensus"
  } else if (is.finite(fp_2_frac) && fp_2_frac >= overlap_mid) {
    "fp_2of3_consensus"
  } else {
    "effect_size_or_sign_stratified"
  }

  rna_set <- sets[["rna_pearson"]]
  if (is.null(rna_set)) rna_set <- character(0)

  tibble::tibble(
    tf = tf0,
    ko_effect_thr = effect_thr,
    fp_union_n = fp_union_n,
    fp_inter_n = fp_inter_n,
    fp_2of3_n = fp_2of3_n,
    fp_3_frac = fp_3_frac,
    fp_2_frac = fp_2_frac,
    rna_pearson_n = length(rna_set),
    strategy = strategy
  )
}

make_venn_plot <- function(callsets_by_method, ko_tbl, log2fc_thr = -0.5, p_thr = 0.05, title = NULL) {
  p_col <- .detect_p_col(ko_tbl)
  ko_tbl <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc))

  ko_sig <- if (is.na(p_col)) {
    ko_tbl |>
      dplyr::filter(.data$log2fc < log2fc_thr) |>
      dplyr::pull(.data$gene_norm) |>
      unique()
  } else {
    ko_tbl |>
      dplyr::filter(.data$log2fc < log2fc_thr, .data[[p_col]] < p_thr) |>
      dplyr::pull(.data$gene_norm) |>
      unique()
  }

  fp_methods <- intersect(c("fp_pearson", "fp_spearman", "fp_kendall"), names(callsets_by_method))
  sets <- lapply(callsets_by_method[fp_methods], function(x) intersect(.norm_chr(x), ko_sig))

  if (requireNamespace("VennDiagram", quietly = TRUE)) {
    p <- tryCatch({
      grob <- VennDiagram::venn.diagram(
        x = sets,
        filename = NULL,
        category.names = names(sets),
        fill = c("#4C78A8", "#F58518", "#54A24B"),
        alpha = 0.5,
        cex = 0.9,
        cat.cex = 0.9,
        lwd = 1,
        col = "black"
      )
      if (requireNamespace("ggplotify", quietly = TRUE)) {
        ggplotify::as.ggplot(grob)
      } else {
        ggplot2::ggplot() +
          ggplot2::annotation_custom(grob) +
          ggplot2::coord_equal()
      }
    }, error = function(e) NULL)
    if (!is.null(p)) {
      p <- p + ggplot2::theme_void() +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0, 0, 0, 0),
          aspect.ratio = 1
        )
      if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
      return(p)
    }
  }
  if (length(sets) < 3) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  .simple_venn_plot(sets, title = title)
}

basal_cache <- new.env(parent = emptyenv())

get_basal_cached <- function(method_key, fp_r, fp_p) {
  cache_key <- sprintf("%s_r%.1f_p%.2g", method_key, fp_r, fp_p)
  if (exists(cache_key, envir = basal_cache, inherits = FALSE)) {
    return(get(cache_key, envir = basal_cache, inherits = FALSE))
  }
  tbl <- .load_basal(method_key, fp_r, fp_p)
  assign(cache_key, tbl, envir = basal_cache)
  tbl
}

summary_overlap <- tibble::tibble()

for (i in seq_len(nrow(tf_ko_map))) {
  tf <- tf_ko_map$tf[[i]]
  cell <- tf_ko_map$cell[[i]]
  sample_id <- tf_ko_map$sample_id[[i]]

  ko_tbl <- ko_truth_list[[tf]]
  if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) {
    message("[skip] KO table missing/empty for TF: ", tf)
    next
  }

  # Predicted gene sets by method and parameter combo
  pred_by_method <- setNames(vector("list", length(methods)), methods)

  for (method_key in methods) {
    combo_sets <- list()
    for (j in seq_len(nrow(param_grid))) {
      fp_r <- param_grid$fp_r[[j]]
      fp_p <- param_grid$fp_p[[j]]

      basal_tbl <- get_basal_cached(method_key, fp_r, fp_p)
      if (is.null(basal_tbl) || !nrow(basal_tbl)) {
        combo_sets[[j]] <- character(0)
        next
      }

      fp_ann_use <- fp_annot_by_method[[method_key]]
      if (is.null(fp_ann_use)) {
        combo_sets[[j]] <- character(0)
        next
      }
      basal_tf <- .gate_basal_for_tf(basal_tbl, tf, sample_id, gr_n, fp_ann_use)
      if (!nrow(basal_tf)) {
        combo_sets[[j]] <- character(0)
        next
      }

      basal_tf <- .pick_best_peak_by_ko(basal_tf, ko_tbl, method_key)
      gene_col <- .detect_basal_cols(basal_tf)$gene
      combo_sets[[j]] <- unique(basal_tf[[gene_col]])
    }
    pred_by_method[[method_key]] <- combo_sets
  }

  # Union callsets per method for overlap summary + Venn
  callsets_by_method <- lapply(pred_by_method, function(x) unique(unlist(x)))
  summary_overlap <- dplyr::bind_rows(summary_overlap, .summarize_method_overlap(tf, ko_tbl, callsets_by_method))

  # Per-TF overlap table across FP methods for each parameter combo.
  fp_methods <- c("fp_pearson", "fp_spearman")
  if (isTRUE(include_kendall)) fp_methods <- c("fp_pearson", "fp_spearman", "fp_kendall")
  ko_lfc05 <- ko_tbl |>
    dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
    dplyr::filter(.data$log2fc < -0.5) |>
    dplyr::pull(.data$gene_norm) |>
    unique()
  overlap_rows <- vector("list", nrow(param_grid))
  for (j in seq_len(nrow(param_grid))) {
    sets_j <- lapply(fp_methods, function(m) .or_empty(pred_by_method[[m]][[j]]))
    names(sets_j) <- fp_methods
    inter_all <- Reduce(intersect, sets_j)
    union_all <- unique(unlist(sets_j))
    uniq_by_method <- lapply(fp_methods, function(m) {
      other_union <- unique(unlist(sets_j[names(sets_j) != m]))
      setdiff(sets_j[[m]], other_union)
    })

    sets_j_lfc <- lapply(sets_j, function(x) intersect(.norm_chr(x), ko_lfc05))
    inter_lfc <- Reduce(intersect, sets_j_lfc)
    uniq_by_method_lfc <- lapply(fp_methods, function(m) {
      other_union_lfc <- unique(unlist(sets_j_lfc[names(sets_j_lfc) != m]))
      setdiff(sets_j_lfc[[m]], other_union_lfc)
    })

    if (isTRUE(include_kendall)) {
      overlap_rows[[j]] <- tibble::tibble(
        tf = tf,
        cell = cell,
        fp_r = param_grid$fp_r[[j]],
        fp_p = param_grid$fp_p[[j]],
        common_n = length(inter_all),
        common_genes = paste(inter_all, collapse = ";"),
        fp_pearson_unique_n = length(uniq_by_method[[1]]),
        fp_pearson_unique_genes = paste(uniq_by_method[[1]], collapse = ";"),
        fp_spearman_unique_n = length(uniq_by_method[[2]]),
        fp_spearman_unique_genes = paste(uniq_by_method[[2]], collapse = ";"),
        fp_kendall_unique_n = length(uniq_by_method[[3]]),
        fp_kendall_unique_genes = paste(uniq_by_method[[3]], collapse = ";"),
        common_n_log2fc_lt_0.5 = length(inter_lfc),
        common_genes_log2fc_lt_0.5 = paste(inter_lfc, collapse = ";"),
        fp_pearson_unique_n_log2fc_lt_0.5 = length(uniq_by_method_lfc[[1]]),
        fp_pearson_unique_genes_log2fc_lt_0.5 = paste(uniq_by_method_lfc[[1]], collapse = ";"),
        fp_spearman_unique_n_log2fc_lt_0.5 = length(uniq_by_method_lfc[[2]]),
        fp_spearman_unique_genes_log2fc_lt_0.5 = paste(uniq_by_method_lfc[[2]], collapse = ";"),
        fp_kendall_unique_n_log2fc_lt_0.5 = length(uniq_by_method_lfc[[3]]),
        fp_kendall_unique_genes_log2fc_lt_0.5 = paste(uniq_by_method_lfc[[3]], collapse = ";")
      )
    } else {
      overlap_rows[[j]] <- tibble::tibble(
        tf = tf,
        cell = cell,
        fp_r = param_grid$fp_r[[j]],
        fp_p = param_grid$fp_p[[j]],
        common_n = length(inter_all),
        common_genes = paste(inter_all, collapse = ";"),
        fp_pearson_unique_n = length(uniq_by_method[[1]]),
        fp_pearson_unique_genes = paste(uniq_by_method[[1]], collapse = ";"),
        fp_spearman_unique_n = length(uniq_by_method[[2]]),
        fp_spearman_unique_genes = paste(uniq_by_method[[2]], collapse = ";"),
        common_n_log2fc_lt_0.5 = length(inter_lfc),
        common_genes_log2fc_lt_0.5 = paste(inter_lfc, collapse = ";"),
        fp_pearson_unique_n_log2fc_lt_0.5 = length(uniq_by_method_lfc[[1]]),
        fp_pearson_unique_genes_log2fc_lt_0.5 = paste(uniq_by_method_lfc[[1]], collapse = ";"),
        fp_spearman_unique_n_log2fc_lt_0.5 = length(uniq_by_method_lfc[[2]]),
        fp_spearman_unique_genes_log2fc_lt_0.5 = paste(uniq_by_method_lfc[[2]], collapse = ";")
      )
    }
  }
  overlap_tbl <- dplyr::bind_rows(overlap_rows)
  if (!isTRUE(save_union_csv) && !isTRUE(main_pdf_only)) {
    overlap_out <- file.path(plot_out_dir, sprintf("%s_%s_fp_method_overlap.csv", tf, cell))
    readr::write_csv(overlap_tbl, overlap_out)
  }

  # Per-TF vectors for 4 methods (KO log2FC < -0.5) by parameter combo.
  if (isTRUE(RUN_VENN_BLOCK) && !isTRUE(save_union_csv) && !isTRUE(main_pdf_only)) {
    venn_out <- file.path(plot_out_dir, sprintf("%s_%s_venn_vectors_log2fc_lt_0.5.csv", tf, cell))
    if (!file.exists(venn_out)) {
      venn_rows <- vector("list", nrow(param_grid))
      for (j in seq_len(nrow(param_grid))) {
        fp_pearson_genes <- intersect(.norm_chr(.or_empty(pred_by_method[["fp_pearson"]][[j]])), ko_lfc05)
        fp_spearman_genes <- intersect(.norm_chr(.or_empty(pred_by_method[["fp_spearman"]][[j]])), ko_lfc05)
        rna_pearson_genes <- intersect(.norm_chr(.or_empty(pred_by_method[["rna_pearson"]][[j]])), ko_lfc05)
        if (isTRUE(include_kendall)) {
          fp_kendall_genes <- intersect(.norm_chr(.or_empty(pred_by_method[["fp_kendall"]][[j]])), ko_lfc05)
          venn_rows[[j]] <- tibble::tibble(
            tf = tf,
            cell = cell,
            fp_r = param_grid$fp_r[[j]],
            fp_p = param_grid$fp_p[[j]],
            fp_pearson_genes = paste(fp_pearson_genes, collapse = ";"),
            fp_spearman_genes = paste(fp_spearman_genes, collapse = ";"),
            fp_kendall_genes = paste(fp_kendall_genes, collapse = ";"),
            rna_pearson_genes = paste(rna_pearson_genes, collapse = ";")
          )
        } else {
          venn_rows[[j]] <- tibble::tibble(
            tf = tf,
            cell = cell,
            fp_r = param_grid$fp_r[[j]],
            fp_p = param_grid$fp_p[[j]],
            fp_pearson_genes = paste(fp_pearson_genes, collapse = ";"),
            fp_spearman_genes = paste(fp_spearman_genes, collapse = ";"),
            rna_pearson_genes = paste(rna_pearson_genes, collapse = ";")
          )
        }
      }
      venn_tbl <- dplyr::bind_rows(venn_rows)
      readr::write_csv(venn_tbl, venn_out)
      message("Saved per-TF Venn vectors: ", venn_out)
    } else {
      message("Using existing Venn vectors: ", venn_out)
    }
  }

  # Build plot grid
  param_labels <- sprintf("r=%.1f\np=%.2g", param_grid$fp_r, param_grid$fp_p)
  row_plots <- list()
  blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()

  for (method_key in methods) {
    row1 <- vector("list", nrow(param_grid))
    row2 <- vector("list", nrow(param_grid))
    for (j in seq_len(nrow(param_grid))) {
      pred_genes <- pred_by_method[[method_key]][[j]]
      if (is.null(pred_genes) || length(pred_genes) == 0) {
        row1[[j]] <- blank_plot
        row2[[j]] <- blank_plot
        next
      }

      res <- .plot_pred_vs_random(
        ko_tbl = ko_tbl,
        pred_genes = pred_genes,
        random_pool = NULL,
        title = if (method_key == methods[[1]]) param_labels[[j]] else NULL,
        min_n = 5L
      )
      row1[[j]] <- res$box
      row2[[j]] <- res$pct
    }

    label_base <- method_labels[[method_key]]
    label_base <- if (is.na(label_base) || !nzchar(label_base)) method_key else label_base

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

  # Venn plots are generated in a separate PDF using Vennerable.

  plot_stack <- row_plots[[1]]
  for (k in 2:length(row_plots)) plot_stack <- plot_stack / row_plots[[k]]

  plot_stack <- plot_stack +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("%s (%s)", tf, cell),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))
    )
  plot_stack <- plot_stack & ggplot2::theme(legend.position = "right")

  if (isTRUE(save_main_pdf)) {
    out_pdf <- file.path(plot_out_dir, sprintf("%s_%s.pdf", tf, cell))
    ggplot2::ggsave(out_pdf, plot_stack, width = 18, height = 14, units = "in", dpi = 300)
    message("Saved: ", out_pdf)
  } else {
    message("Skip main PDF save (EPISCOPE_SAVE_MAIN_PDF=FALSE).")
  }

  # Venn vectors already handled above.

  # Optional: plot Venns from the saved per-TF CSV (one PDF per parameter combo).
  if (isTRUE(RUN_VENN_BLOCK) && isTRUE(save_venn_pdf) && !isTRUE(main_pdf_only) && !isTRUE(save_union_csv)) {
    if (!requireNamespace("Vennerable", quietly = TRUE)) {
      message("Vennerable not installed; skipping Venn PDFs for ", tf)
    } else {
      suppressPackageStartupMessages(library(Vennerable))
      venn_tbl_in <- readr::read_csv(venn_out, show_col_types = FALSE)
      for (j in seq_len(nrow(venn_tbl_in))) {
        sets_j <- list(
          fp_pearson  = .parse_genes(venn_tbl_in$fp_pearson_genes[[j]]),
          fp_spearman = .parse_genes(venn_tbl_in$fp_spearman_genes[[j]]),
          rna_pearson = .parse_genes(venn_tbl_in$rna_pearson_genes[[j]])
        )
        if (isTRUE(include_kendall) && "fp_kendall_genes" %in% names(venn_tbl_in)) {
          sets_j$fp_kendall <- .parse_genes(venn_tbl_in$fp_kendall_genes[[j]])
        }

        title_j <- sprintf("r=%.1f p=%.3g", venn_tbl_in$fp_r[[j]], venn_tbl_in$fp_p[[j]])
        venn_pdf <- file.path(
          plot_out_dir,
          sprintf("%s_%s_venn_chowruskey_r%.1f_p%.3g.pdf", tf, cell, venn_tbl_in$fp_r[[j]], venn_tbl_in$fp_p[[j]])
        )

        tryCatch({
          grDevices::pdf(venn_pdf, width = 6, height = 6, onefile = FALSE)
          v <- Vennerable::Venn(sets_j)
          all_inter <- Reduce(intersect, sets_j)
          any_empty <- any(vapply(sets_j, function(x) length(x) == 0L, logical(1)))
          use_cr <- length(all_inter) > 0 && !any_empty

          if (use_cr) {
            plot(v, doWeights = TRUE, type = "ChowRuskey", show = list(SetLabels = TRUE, DarkMatter = FALSE))
          } else {
            plot(v, doWeights = FALSE, type = "simple", show = list(SetLabels = TRUE, DarkMatter = FALSE))
          }
          graphics::title(main = title_j)
        }, error = function(e) {
          message("Vennerable Venn PDF failed for ", tf, " r=", venn_tbl_in$fp_r[[j]], " p=", venn_tbl_in$fp_p[[j]], ": ", conditionMessage(e))
        }, finally = {
          try(grDevices::dev.off(), silent = TRUE)
        })
      }
    }
  }
}

if (nrow(summary_overlap)) {
  if (!isTRUE(save_union_csv) && !isTRUE(main_pdf_only)) {
    out_csv <- file.path(plot_out_dir, "method_overlap_summary_all_tfs.csv")
    readr::write_csv(summary_overlap, out_csv)
    message("Saved overlap summary: ", out_csv)
  }
}

# Standalone: union genes across fp_pearson/fp_spearman/rna_pearson for r=0.3 p=0.1,
# with per-method p/r and KO log2FC.

.summarize_gene_stats <- function(tbl, gene_col, r_col, p_col, p_adj_col) {
  tbl <- tbl |>
    dplyr::mutate(
      gene_norm = .norm_chr(.data[[gene_col]]),
      r_val = if (!is.na(r_col)) as.numeric(.data[[r_col]]) else NA_real_,
      p_val = if (!is.na(p_col)) as.numeric(.data[[p_col]]) else NA_real_,
      p_adj_val = if (!is.na(p_adj_col)) as.numeric(.data[[p_adj_col]]) else NA_real_,
      fp_rsd_val = if ("fp_rsd" %in% names(tbl)) as.numeric(.data$fp_rsd) else NA_real_
    ) |>
    dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "")

  if (!nrow(tbl)) {
    return(tibble::tibble(gene_norm = character(0), r = numeric(0), p = numeric(0), p_adj = numeric(0)))
  }

  first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) NA_real_ else x[[1]]
  }

  tbl |>
    dplyr::group_by(.data$gene_norm) |>
    dplyr::summarise(
      r = first_non_na(.data$r_val),
      p = first_non_na(.data$p_val),
      p_adj = first_non_na(.data$p_adj_val),
      fp_rsd = first_non_na(.data$fp_rsd_val),
      .groups = "drop"
    )
}

target_fp_r <- 0.3
target_fp_p <- 0.1
export_methods <- c("fp_pearson", "fp_spearman", "rna_pearson")
basal_full_cache <- new.env(parent = emptyenv())
if (is.null(gr_n$fp_variance) && is.data.frame(gr_n$fp_score)) {
  gr_n$fp_variance <- compute_hv_variance_tbl(gr_n$fp_score, id_cols = "peak_ID")
}
if (is.null(gr_n$rna_variance) && is.data.frame(gr_n$rna)) {
  gr_n$rna_variance <- compute_hv_variance_tbl(gr_n$rna, id_cols = c("ensembl_gene_id", "HGNC"))
}
fp_rsd_map <- gr_n$fp_variance |>
  dplyr::select(peak_ID, fp_rsd = rsd) |>
  dplyr::distinct(.data$peak_ID, .keep_all = TRUE)
rna_rsd_map <- gr_n$rna_variance |>
  dplyr::select(ensembl_gene_id, HGNC, rna_rsd = rsd)

.lookup_rna_rsd <- function(gene_norm, rna_tbl) {
  hgnc_norm <- .norm_chr(rna_tbl$HGNC)
  ens_norm <- .norm_chr(rna_tbl$ensembl_gene_id)
  idx <- match(.norm_chr(gene_norm), hgnc_norm)
  out <- rna_tbl$rna_rsd[idx]
  miss <- is.na(out)
  if (any(miss)) {
    idx2 <- match(.norm_chr(gene_norm[miss]), ens_norm)
    out[miss] <- rna_tbl$rna_rsd[idx2]
  }
  out
}

.get_basal_full <- function(method_key) {
  if (exists(method_key, envir = basal_full_cache, inherits = FALSE)) {
    return(get(method_key, envir = basal_full_cache, inherits = FALSE))
  }


  fp_corr <- switch(
    method_key,
    fp_pearson = if (exists("fp_res_full_pearson", inherits = TRUE)) fp_res_full_pearson$fp_gene_corr_full else NULL,
    fp_spearman = if (exists("fp_res_full_spearman", inherits = TRUE)) fp_res_full_spearman$fp_gene_corr_full else NULL,
    rna_pearson = if (exists("fp_res_full_pearson", inherits = TRUE)) fp_res_full_pearson$fp_gene_corr_full else NULL,
    NULL
  )
  fp_ann <- switch(
    method_key,
    fp_pearson = gr_n$fp_annotation_pearson,
    fp_spearman = gr_n$fp_annotation_spearman,
    rna_pearson = gr_n$fp_annotation_pearson,
    NULL
  )
  if (is.null(fp_corr) || is.null(fp_ann)) {
    assign(method_key, NULL, envir = basal_full_cache)
    return(NULL)
  }

  full_tbl <- make_basal_links(
    fp_gene_corr_kept = fp_corr,
    fp_annotation     = fp_ann,
    out_dir           = file.path(tempdir(), paste0("basal_full_", method_key)),
    prefix            = "lighting",
    rna_tbl           = gr_n$rna,
    rna_method        = "pearson",
    rna_cores         = 10
  )
  assign(method_key, full_tbl, envir = basal_full_cache)
  full_tbl
}

if (isTRUE(save_union_csv)) {
  for (i in seq_len(nrow(tf_ko_map))) {
    tf <- tf_ko_map$tf[[i]]
    cell <- tf_ko_map$cell[[i]]
    sample_id <- tf_ko_map$sample_id[[i]]

    ko_tbl <- ko_truth_list[[tf]]
    if (!is.data.frame(ko_tbl) || !nrow(ko_tbl)) {
      message("[skip] KO table missing/empty for TF: ", tf)
      next
    }

    ko_map <- ko_tbl |>
      dplyr::mutate(gene_norm = .norm_chr(.data$gene)) |>
      dplyr::filter(!is.na(.data$gene_norm), .data$gene_norm != "", is.finite(.data$log2fc)) |>
      dplyr::group_by(.data$gene_norm) |>
      dplyr::summarise(
        gene = .data$gene[[1]],
        log2fc = mean(as.numeric(.data$log2fc), na.rm = TRUE),
        .groups = "drop"
      )

    per_method <- list()
    method_genes <- list()

    for (method_key in export_methods) {
      basal_tbl <- get_basal_cached(method_key, target_fp_r, target_fp_p)
      if (is.null(basal_tbl) || !nrow(basal_tbl)) {
        per_method[[method_key]] <- tibble::tibble(gene_norm = character(0), r = numeric(0), p = numeric(0), p_adj = numeric(0))
        method_genes[[method_key]] <- character(0)
        next
      }

      fp_ann_use <- fp_annot_by_method[[method_key]]
      if (is.null(fp_ann_use)) {
        per_method[[method_key]] <- tibble::tibble(gene_norm = character(0), r = numeric(0), p = numeric(0), p_adj = numeric(0))
        method_genes[[method_key]] <- character(0)
        next
      }

      basal_full <- .get_basal_full(method_key)
      basal_src <- if (is.data.frame(basal_full) && nrow(basal_full)) basal_full else basal_tbl
      basal_tf <- .gate_basal_for_tf(basal_src, tf, sample_id, gr_n, fp_ann_use)
      if (!nrow(basal_tf)) {
        per_method[[method_key]] <- tibble::tibble(gene_norm = character(0), r = numeric(0), p = numeric(0), p_adj = numeric(0))
        method_genes[[method_key]] <- character(0)
        next
      }

      basal_tf <- .pick_best_peak_by_ko(basal_tf, ko_tbl, method_key)
      peak_col <- .detect_basal_cols(basal_tf)$peak
      basal_tf <- basal_tf |>
        dplyr::left_join(fp_rsd_map, by = setNames("peak_ID", peak_col))
      gene_col <- .detect_basal_cols(basal_tf)$gene
      stat_cols <- .detect_stat_cols(basal_tf, method_key)
      stats_tbl <- .summarize_gene_stats(basal_tf, gene_col, stat_cols$r, stat_cols$p, stat_cols$p_adj)
      per_method[[method_key]] <- stats_tbl
      method_genes[[method_key]] <- stats_tbl$gene_norm
    }

    union_genes <- unique(unlist(method_genes))
    if (!length(union_genes)) {
      message("[skip] No union genes for TF: ", tf)
      next
    }

    out_tbl <- tibble::tibble(gene_norm = union_genes) |>
      dplyr::left_join(ko_map, by = "gene_norm") |>
      dplyr::mutate(gene = dplyr::coalesce(.data$gene, .data$gene_norm)) |>
      dplyr::select(.data$gene_norm, .data$gene, .data$log2fc)

    for (method_key in export_methods) {
      method_tbl <- per_method[[method_key]]
      method_tbl <- method_tbl |>
        dplyr::rename(
          !!paste0(method_key, "_r") := .data$r,
          !!paste0(method_key, "_p") := .data$p,
          !!paste0(method_key, "_p_adj") := .data$p_adj,
          !!paste0(method_key, "_fp_rsd") := .data$fp_rsd
        )
      out_tbl <- out_tbl |>
        dplyr::left_join(method_tbl, by = "gene_norm")
    }

    out_tbl <- out_tbl |>
      dplyr::mutate(
        fp_rsd = dplyr::coalesce(.data$fp_pearson_fp_rsd, .data$fp_spearman_fp_rsd, .data$rna_pearson_fp_rsd),
        rna_rsd = .lookup_rna_rsd(.data$gene_norm, rna_rsd_map)
      ) |>
      dplyr::select(
        .data$gene,
        .data$log2fc,
        .data$fp_rsd,
        .data$rna_rsd,
        dplyr::any_of(c(
          "fp_pearson_r", "fp_pearson_p", "fp_pearson_p_adj",
          "fp_spearman_r", "fp_spearman_p", "fp_spearman_p_adj",
          "rna_pearson_r", "rna_pearson_p", "rna_pearson_p_adj"
        ))
      )
    out_tbl <- out_tbl |>
      dplyr::filter(!is.na(.data$log2fc))

    out_csv <- file.path(
      plot_out_dir,
      sprintf("%s_%s_union_fp_r%.1f_p%.2g_methods.csv", tf, cell, target_fp_r, target_fp_p)
    )
    readr::write_csv(out_tbl, out_csv)
    message("Saved union method stats: ", out_csv)
  }
}

if (isTRUE(save_union_csv)) {
  union_files_out <- list.files(
    path = plot_out_dir,
    pattern = "_union_fp_r0\\.3_p0\\.1_methods\\.csv$",
    full.names = TRUE
  )
  if (length(union_files_out)) {
    get_col <- function(tbl, nm) if (nm %in% names(tbl)) tbl[[nm]] else NA_real_
    for (f in union_files_out) {
      tbl <- readr::read_csv(f, show_col_types = FALSE)
      pass <- (get_col(tbl, "fp_pearson_r") > 0.3 & get_col(tbl, "fp_pearson_p_adj") < 0.05) |
        (get_col(tbl, "fp_spearman_r") > 0.3 & get_col(tbl, "fp_spearman_p_adj") < 0.05) |
        (get_col(tbl, "rna_pearson_r") > 0.3 & get_col(tbl, "rna_pearson_p_adj") < 0.05)
      tbl <- tbl[which(pass), , drop = FALSE]
      readr::write_csv(tbl, f)
    }
  }
}
