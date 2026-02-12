# Build isolated method-comparison tables including ATAC Pearson and FP/RNA hybrid.
# This script does not overwrite existing tables in lighting_method_comparison.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})


# Standalone defaults (edit here if needed)
workspace_rdata <- "/data/homes/yl814/episcope_test/benchmark_tf_to_target_genes_prediction/workspace.RData"
source_dir <- "/data/homes/yl814/episcope_test/nutrient_stress/lighting_method_comparison"
out_dir <- "/data/homes/yl814/episcope_test/nutrient_stress/lighting_method_comparison_atac_benchmark"
base_dir <- "/data/homes/yl814/episcope_test/nutrient_stress"

required_syms <- c("grn_set", "base_dir")
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
base_dir <- get("base_dir", inherits = TRUE)

if (!dir.exists(source_dir)) {
  stop("Source method directory not found: ", source_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fp_r_vals <- c(0.3, 0.5)
fp_p_vals <- c(0.1, 0.05, 0.01, 0.001)

threshold_atac_gene_corr_p <- 0.05
threshold_atac_gene_corr_abs_r <- 0.3
basal_p_thr <- 0.05
basal_r_thr <- 0.3

message("Copying existing FP/RNA method tables into isolated folder...")
copy_patterns <- c(
  "^basal_filtered_method_pearson_fp_r",
  "^basal_filtered_method_spearman_fp_r",
  "^basal_filtered_method_rna_pearson_rna_r"
)
copy_files <- list.files(source_dir, full.names = TRUE)
copy_files <- copy_files[vapply(
  basename(copy_files),
  function(x) any(grepl(paste(copy_patterns, collapse = "|"), x)),
  logical(1)
)]
for (f in copy_files) {
  file.copy(f, file.path(out_dir, basename(f)), overwrite = TRUE)
}

message("Computing ATAC Pearson (GeneHancer) correlations...")
gh_path <- system.file("extdata", "GeneHancer_v5.24_elite_panc.csv", package = "episcope")
if (!nzchar(gh_path)) {
  gh_path <- file.path("inst", "extdata", "GeneHancer_v5.24_elite_panc.csv")
}
if (!file.exists(gh_path)) {
  stop("Could not find GeneHancer file: ", gh_path)
}

# GeneHancer candidate links (same source style as previous method-comparison block).
gh_std <- episcope::load_genehancer_panc(gh_path)

atac_res <- episcope::correlate_atac_to_genes(
  grn_set = gr_n,
  gh_tbl = gh_std,
  gene_mode = "both",
  fdr = threshold_atac_gene_corr_p,
  r_abs_min = threshold_atac_gene_corr_abs_r,
  cor_method = "pearson",
  cache_dir = file.path(base_dir, "cache", "atac_gene_corr"),
  cache_tag = "nutrient_genehancer_pearson",
  workers = 30,
  cache_verbose = TRUE
)

readr::write_csv(
  atac_res$atac_gene_corr_full,
  file.path(out_dir, "atac_gene_corr_full_genehancer_pearson.csv")
)

message("Building ATAC Pearson basal links...")
fp_annot <- if (is.data.frame(gr_n$fp_annotation_pearson)) gr_n$fp_annotation_pearson else gr_n$fp_annotation

atac_fp_gene <- dplyr::inner_join(
  fp_annot |>
    dplyr::select(fp_peak, atac_peak, motifs, tfs) |>
    dplyr::distinct(),
  atac_res$atac_gene_corr_full |>
    dplyr::select(atac_peak, gene_key, n_atac, r_atac, p_atac, p_adj_atac) |>
    dplyr::distinct(),
  by = "atac_peak"
) |>
  dplyr::distinct(.data$fp_peak, .data$gene_key, .data$tfs, .keep_all = TRUE) |>
  dplyr::transmute(
    fp_peak = .data$fp_peak,
    gene_key = .data$gene_key,
    atac_peak = .data$atac_peak,
    tfs = .data$tfs,
    motifs = .data$motifs,
    n_fp = .data$n_atac,
    r_fp = .data$r_atac,
    p_fp = .data$p_atac,
    p_adj_fp = .data$p_adj_atac
  )

basal_atac <- episcope::make_basal_links(
  fp_gene_corr_kept = atac_fp_gene,
  fp_annotation = fp_annot,
  out_dir = file.path(out_dir, "basal_tmp_atac_pearson"),
  prefix = "lighting",
  rna_tbl = gr_n$rna,
  rna_method = "pearson",
  rna_cores = 20
)

basal_atac <- basal_atac |>
  dplyr::filter(.data$p_adj_tf < basal_p_thr, .data$r_tf > basal_r_thr)

for (fp_r_thr in fp_r_vals) {
  for (fp_p_thr in fp_p_vals) {
    out_atac <- basal_atac |>
      dplyr::filter(!is.na(.data$r_gene), !is.na(.data$p_adj_gene)) |>
      dplyr::filter(abs(.data$r_gene) >= fp_r_thr, .data$p_adj_gene < fp_p_thr)

    readr::write_csv(
      out_atac,
      file.path(out_dir, sprintf("basal_filtered_method_atac_pearson_r%.1f_p%.2g.csv", fp_r_thr, fp_p_thr))
    )
  }
}

message("Building FP/RNA hybrid method tables...")
for (fp_r_thr in fp_r_vals) {
  for (fp_p_thr in fp_p_vals) {
    in_fp_pearson <- file.path(out_dir, sprintf("basal_filtered_method_pearson_fp_r%.1f_p%.2g.csv", fp_r_thr, fp_p_thr))
    in_fp_spearman <- file.path(out_dir, sprintf("basal_filtered_method_spearman_fp_r%.1f_p%.2g.csv", fp_r_thr, fp_p_thr))
    if (!file.exists(in_fp_pearson) || !file.exists(in_fp_spearman)) {
      warning(
        "Missing source table(s) for hybrid at r=", fp_r_thr, " p=", fp_p_thr,
        " | pearson: ", file.exists(in_fp_pearson),
        " | spearman: ", file.exists(in_fp_spearman)
      )
      next
    }
    tbl_p <- readr::read_csv(in_fp_pearson, show_col_types = FALSE)
    tbl_s <- readr::read_csv(in_fp_spearman, show_col_types = FALSE)

    key_cols <- c("TF", "gene_key", "peak_ID")
    if (!all(key_cols %in% names(tbl_p)) || !all(key_cols %in% names(tbl_s))) {
      stop("Hybrid build requires key columns TF, gene_key, peak_ID in both Pearson and Spearman tables.")
    }

    p_tbl <- tbl_p |>
      dplyr::select(
        dplyr::all_of(key_cols),
        r_gene_p = "r_gene",
        p_adj_gene_p = "p_adj_gene",
        r_rna_p = "r_rna_gene"
      ) |>
      dplyr::distinct()

    s_tbl <- tbl_s |>
      dplyr::select(
        dplyr::all_of(key_cols),
        r_gene_s = "r_gene",
        p_adj_gene_s = "p_adj_gene",
        r_rna_s = "r_rna_gene"
      ) |>
      dplyr::distinct()

    # Hybrid as true union of methods:
    # filter positivity per source table independently, then take union.
    p_pos_tbl <- p_tbl |>
      dplyr::filter(
        !is.na(.data$r_rna_p), !is.na(.data$r_gene_p),
        .data$r_rna_p > 0, .data$r_gene_p > 0
      )
    s_pos_tbl <- s_tbl |>
      dplyr::filter(
        !is.na(.data$r_rna_s), !is.na(.data$r_gene_s),
        .data$r_rna_s > 0, .data$r_gene_s > 0
      )

    # 2) Significance check per method separately.
    p_sig_keys <- p_pos_tbl |>
      dplyr::filter(
        !is.na(.data$p_adj_gene_p),
        abs(.data$r_gene_p) >= fp_r_thr,
        .data$p_adj_gene_p < fp_p_thr
      ) |>
      dplyr::select(dplyr::all_of(key_cols))

    s_sig_keys <- s_pos_tbl |>
      dplyr::filter(
        !is.na(.data$p_adj_gene_s),
        abs(.data$r_gene_s) >= fp_r_thr,
        .data$p_adj_gene_s < fp_p_thr
      ) |>
      dplyr::select(dplyr::all_of(key_cols))

    p_rows <- tbl_p |>
      dplyr::inner_join(p_sig_keys, by = key_cols) |>
      dplyr::mutate(hybrid_source = "pearson")

    s_rows <- tbl_s |>
      dplyr::inner_join(s_sig_keys, by = key_cols) |>
      dplyr::mutate(hybrid_source = "spearman")

    support_tbl <- dplyr::full_join(
      p_sig_keys |> dplyr::mutate(in_pearson = TRUE),
      s_sig_keys |> dplyr::mutate(in_spearman = TRUE),
      by = key_cols
    ) |>
      dplyr::mutate(
        hybrid_support = dplyr::case_when(
          dplyr::coalesce(.data$in_pearson, FALSE) & dplyr::coalesce(.data$in_spearman, FALSE) ~ "both",
          dplyr::coalesce(.data$in_pearson, FALSE) ~ "pearson_only",
          dplyr::coalesce(.data$in_spearman, FALSE) ~ "spearman_only",
          TRUE ~ "none"
        )
      ) |>
      dplyr::select(dplyr::all_of(key_cols), "hybrid_support")

    # 3) Union Pearson-significant and Spearman-significant sets.
    hyb <- dplyr::bind_rows(p_rows, s_rows) |>
      dplyr::arrange(dplyr::desc(.data$hybrid_source == "pearson")) |>
      dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE) |>
      dplyr::left_join(support_tbl, by = key_cols)

    message(
      sprintf(
        "[hybrid r=%.1f p=%.2g] pearson=%d spearman=%d pos_p=%d pos_s=%d p_sig=%d s_sig=%d union=%d",
        fp_r_thr, fp_p_thr, nrow(tbl_p), nrow(tbl_s),
        nrow(p_pos_tbl), nrow(s_pos_tbl), nrow(p_sig_keys), nrow(s_sig_keys), nrow(hyb)
      )
    )

    readr::write_csv(
      hyb,
      file.path(out_dir, sprintf("basal_filtered_method_fp_rna_hybrid_fp_r%.1f_p%.2g.csv", fp_r_thr, fp_p_thr))
    )
  }
}

message("Done. New method tables written to: ", out_dir)
